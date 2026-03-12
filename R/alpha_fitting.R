#' @importFrom dplyr %>%
#' @importFrom rlang %||% !! := sym
#' @importFrom stats quantile
#' @importFrom utils head tail
NULL

utils::globalVariables(c(
  "dbh", "log_x", "x", "trunc_output", "variable", "peaks",
  ".plot_id", "crown_plot_id", ".plot_id_majority", "n"
))


#' Stan model string for estimating alpha with LAI and breakpoint corrections
#'
#' This is the default Stan model used by [fit_alpha_model()]. You can modify it by replacing the `stan_model_code` argument.
#' @format A character string containing the full Stan model code.
#' @export
stan_alpha_model <- "
data {
  int<lower=0> N;
  real<lower=0> x_min;
  real<lower=0> trunc_point;
  real<lower=trunc_point> trunc_upper;
  vector<lower=trunc_point>[N] x;
  real<lower=0, upper=1> LAI_norm;
  real<lower=0, upper=1> breakpoint_norm;
  real prior_mean;
  real prior_sd;
}
parameters {
  real<lower=0> alpha;
}
transformed parameters {
  real adjustment_factor = 1 - sqrt(LAI_norm * breakpoint_norm);
}
model {
  alpha ~ normal(prior_mean, prior_sd) T[0, ];
  real p_trunc = pareto_cdf(trunc_upper | x_min, alpha) - pareto_cdf(trunc_point | x_min, alpha);
  for (n in 1:N) {
    target += adjustment_factor * (pareto_lpdf(x[n] | x_min, alpha) - log(p_trunc));
  }
}

"



#' Identify Potential Breakpoint and Prepare Kernel Density Data
#'
#' This function performs Kernel Density Estimation (KDE) and bootstrapping
#' to identify a potential visibility breakpoint in size-abundance data.
#' It returns the estimated breakpoint and the processed KDE data. The name of the
#' column in the data frame that holds the size information must be called "dbh".
#'
#' @param data A data frame containing individual measurements (e.g., trees). Measurements must be called "dbh".
#' @param n_bootstrap Number of bootstrap replicates for KDE smoothing (default = 1000).
#' @param bandwidth Bandwidth method passed to `density()` (default = "SJ").
#' @param trim_max Upper DBH size limit in cm (default = 50). Trees larger than this value are
#'   excluded before breakpoint detection.
#'
#' @return A list with:
#' \describe{
#'   \item{potential_breakpoint}{A numeric value representing the estimated lower log10 size bound breakpoint.}
#'   \item{bootstrap_kde_log}{A data frame with log-scaled KDE values (`log_x`, `mean_log_density`).}
#'   \item{original_data_trimmed}{The original size vector *after* trimming by `trim_max`. This is crucial for quantile calculation in the first function and for the final filtering in the second.}
#'   \item{original_raw_data_df}{The original input data frame, untouched, so the second function can filter it.}
#'   \item{trim_max_value}{The `trim_max` value used in this function, passed along for consistency.}
#' }
#' @examples
#' \dontrun{
#' data("harv_data_sample")
#' df <- harv_data_sample[harv_data_sample$IDhectbest == 1, "dbh", drop = FALSE]
#' kde_output <- potential_break(df)
#' }
#' @export
potential_break <- function(data,
                            n_bootstrap = 1000,
                            bandwidth = "SJ",
                            trim_max = 50) {

  # Store original raw data for later use in the second function
  original_raw_data_df <- data
  trim_max_value <- trim_max

  # Pull and trim size vector
  size_vector_trimmed <- data$dbh
  size_vector_trimmed <- size_vector_trimmed[!is.na(size_vector_trimmed)]
  size_vector_trimmed <- size_vector_trimmed[size_vector_trimmed <= trim_max]
  if (length(size_vector_trimmed) < 25) {
    stop("Not enough observations below trim_max for breakpoint estimation (min 25 required).")
  }

  # KDE smoothing
  kde <- stats::density(size_vector_trimmed, bw = bandwidth)
  kde_df <- data.frame(x = kde$x, y = kde$y)
  observed_vals <- sort(unique(size_vector_trimmed)) # Used for filtering KDE results

  # Filter KDE results to observed data range
  kde_df <- kde_df %>%
    dplyr::filter(x >= min(size_vector_trimmed), x <= max(size_vector_trimmed))

  # Further filter KDE to points near observed values
  filtered_kde <- kde_df %>%
    dplyr::rowwise() %>%
    dplyr::filter(any(abs(x - observed_vals) <= 0.5)) %>%
    dplyr::ungroup()

  x_values <- filtered_kde$x

  # Bootstrap KDEs
  bootstrap_kdes <- replicate(
    n_bootstrap,
    stats::density(sample(size_vector_trimmed, length(size_vector_trimmed), replace = TRUE), bw = bandwidth),
    simplify = FALSE
  )

  densities_matrix <- sapply(bootstrap_kdes, function(k) {
    stats::approx(k$x, k$y, xout = x_values)$y
  })

  mean_density <- rowMeans(densities_matrix, na.rm = TRUE)

  # Log-scaled kernel density data frame
  bootstrap_kde_log <- data.frame(
    log_x = log10(x_values),
    mean_log_density = log10(mean_density)
  )

  # Identify local peaks using splus2R::peaks()
  peak_df <- splus2R::peaks(bootstrap_kde_log$mean_log_density) %>%
    cbind(bootstrap_kde_log) %>%
    dplyr::filter(. == TRUE)

  if (nrow(peak_df) == 0) {
    stop("No distinct peaks detected in the KDE curve for breakpoint identification.")
  }

  # Calculate the candidate breakpoint
  potential_breakpoint <- peak_df %>%
    dplyr::filter(log_x <= quantile(log10(size_vector_trimmed), 0.75)) %>%
    dplyr::filter(log_x == max(log_x)) %>%
    dplyr::pull(log_x)

  # Handle cases where no breakpoint is found after filtering
  if (length(potential_breakpoint) == 0) {
    warning("No suitable breakpoint found after filtering peaks. Returning NA for breakpoint.")
    potential_breakpoint <- NA
  }

  return(list(
    potential_breakpoint = potential_breakpoint,
    bootstrap_kde_log = bootstrap_kde_log,
    original_data_trimmed = size_vector_trimmed, # Return the trimmed size vector
    original_raw_data_df = original_raw_data_df, # Pass the original DF
    trim_max_value = trim_max_value              # Pass trim_max
  ))
}

#' Determine Truncation Points and Filter Data for Power-Law Modeling
#'
#' This function takes the potential breakpoint and KDE data from a previous step,
#' performs segmented regression, determines the upper truncation point,
#' and filters the original dataset to the identified size range.
#'
#' @param breakpoint_kde_results A list returned by `potential_break`.
#'   It must contain `potential_breakpoint`, `bootstrap_kde_log`,
#'   and `original_data_trimmed`.
#' @param min_size Numeric. Minimum tree size (DBH in cm) to retain in the
#'   filtered data (default = 10). Trees below this threshold are excluded
#'   regardless of the detected breakpoint.
#'
#' @return A list with:
#' \describe{
#'   \item{final_breakpoint}{The final determined lower log10 size bound.}
#'   \item{bayesian_data}{Subset of original data between 10^breakpoint and 10^upper_bound.}
#'   \item{kerneldens_logtransform}{The full log-scaled KDE data from the first function.}
#' }
#' @examples
#' \dontrun{
#' data("harv_data_sample")
#' df <- harv_data_sample[harv_data_sample$IDhectbest == 1, "dbh", drop = FALSE]
#' kde_output   <- potential_break(df)
#' trunc_output <- truncate_filter(kde_output)
#' }
#' @export
truncate_filter <- function(breakpoint_kde_results, min_size = 10) {

  # Extract necessary components from the input list
  potential_breakpoint <- breakpoint_kde_results$potential_breakpoint
  bootstrap_kde_log <- breakpoint_kde_results$bootstrap_kde_log
  original_data_trimmed_df <- data.frame(dbh=breakpoint_kde_results$original_data_trimmed)

  # For robustness, let's ensure potential_breakpoint is not NA
  if (is.na(potential_breakpoint)) {
    stop("Potential breakpoint is NA. Cannot proceed with segmented regression. Check 'potential_break' output.")
  }

  # Filter the data for segmented regression *before* passing it to lm
  bootstrap_kde_log_for_lm <- bootstrap_kde_log[bootstrap_kde_log$log_x >= potential_breakpoint, ]
  # Ensure there's enough data for lm after filtering
  if (nrow(bootstrap_kde_log_for_lm) < 2) { # Need at least 2 points for a line for lm
    stop("Not enough data points after breakpoint filtering for segmented regression.")
  }

  # Segmented regression for upper truncation decision
  lm_fit <- stats::lm(mean_log_density ~ log_x, data = bootstrap_kde_log_for_lm)

  lm_fit$call$data<-bootstrap_kde_log_for_lm

  seg_model <- segmented::selgmented(lm_fit,msg = FALSE)

  breakpoints <- seg_model$psi[, "Est."]

  slopes <- segmented::slope(seg_model)$log_x[, 1]
  segment_boundaries <- c(min(bootstrap_kde_log$log_x), breakpoints, max(bootstrap_kde_log$log_x))

  # Create segment data frame
  segments_df <- data.frame(
    left_x = head(segment_boundaries, -1),     # Left boundary of each segment
    right_x = tail(segment_boundaries, -1),   # Right boundary of each segment
    slope = slopes                    # Slope of each segment
  )

  if(segments_df[nrow(segments_df),]$slope<=-min_size | segments_df[nrow(segments_df),]$slope>0){
    bayesian_data<-original_data_trimmed_df%>%
      dplyr::filter(dbh>=10^potential_breakpoint)%>%
      dplyr::filter(dbh>=min_size)%>%
      dplyr::filter(dbh<=10^segments_df[nrow(segments_df),]$left_x)
    bootstrap_kde_log<-bootstrap_kde_log%>%
      dplyr::filter(log_x>=potential_breakpoint)%>%
      dplyr::filter(log_x>=log10(min_size))%>%
      dplyr::filter(log_x<=segments_df[nrow(segments_df),]$left_x)
  }else{
    bayesian_data<-original_data_trimmed_df%>%
      dplyr::filter(dbh>=10^potential_breakpoint)%>%
      dplyr::filter(dbh>=min_size)

    bootstrap_kde_log<-bootstrap_kde_log%>%
      dplyr::filter(log_x>=potential_breakpoint)%>%
      dplyr::filter(log_x>=log10(min_size))
  }
  return(list(
    bayesian_data = bayesian_data,
    kerneldens_logtransform = bootstrap_kde_log, # Return the full KDE data for potential plotting
    final_breakpoint = potential_breakpoint
  ))
}

#' Map Forest Structural Estimates to a Spatial Grid
#'
#' Joins a data frame of per-plot \code{ForestForTrees} results (alpha, Ntot,
#' or any other plot-level estimates) to a spatial polygon layer, producing an
#' \code{sf} object suitable for direct export as a shapefile or conversion to
#' raster format.
#'
#' Two workflows are supported:
#' \itemize{
#'   \item \strong{Pre-existing plot polygons}: Supply an \code{sf} object via
#'     \code{plots_sf}. Each row is a plot polygon already in the right shape.
#'     The results data frame is joined to it by the shared plot ID column.
#'   \item \strong{No pre-existing polygons}: Leave \code{plots_sf = NULL} and
#'     supply a crown segmentation \code{sf} via \code{crowns_sf}. A regular
#'     square grid of \code{cellsize} metres is constructed over the crown
#'     extent, and plot IDs are assigned to cells by the majority of crown
#'     centroids falling within each cell.
#' }
#'
#' @param results_df A data frame of plot-level estimates (e.g., output from
#'   \code{fit_alpha_model()} or \code{estimate_total_trees()} collated across
#'   plots). Must contain one row per plot and a column identifying each plot.
#' @param results_id_col Character. Name of the plot ID column in
#'   \code{results_df}. Default \code{"plot_id"}.
#' @param plots_sf Optional \code{sf} object of pre-existing plot polygons
#'   (e.g., a ForestGEO shapefile). If supplied, \code{crowns_sf} and
#'   \code{cellsize} are ignored.
#' @param plots_id_col Character. Name of the plot ID column in
#'   \code{plots_sf}. Default \code{"plot_id"}.
#' @param crowns_sf Optional \code{sf} object of individual tree crown polygons,
#'   used to construct a grid when \code{plots_sf = NULL}. Must be in a
#'   projected CRS (metres).
#' @param crowns_id_col Character. Name of an existing plot ID column in
#'   \code{crowns_sf}, if one exists. If \code{NULL} (default), plot IDs are
#'   assigned by spatial intersection with the generated grid.
#' @param cellsize Numeric. Side length of grid cells in CRS units (metres)
#'   when building a grid from \code{crowns_sf}. Default \code{100}
#'   (1-hectare cells).
#'
#' @return An \code{sf} polygon object in WGS84 (EPSG:4326) with all columns
#'   from \code{results_df} joined to the plot geometries. Cells with no
#'   matching results retain their geometry with \code{NA} values.
#'
#' @examples
#' \dontrun{
#' # --- Case 1: pre-existing ForestGEO plot polygons ---
#' library(sf)
#' plots  <- read_sf("HARVplotsGEOJan25.shp") |> st_set_crs(32618)
#' # alpha_df and trees_df are data frames collated from fit_alpha_model()
#' # and estimate_total_trees() across plots
#' alpha_map <- map_grid_estimates(
#'   results_df    = alpha_df,
#'   results_id_col = "IDhectbest",
#'   plots_sf      = plots,
#'   plots_id_col  = "IDhectbest"
#' )
#' write_sf(alpha_map, "harv_alpha.shp")
#'
#' # --- Case 2: build grid from crown shapefile ---
#' crowns <- read_sf("crowns.shp") |> st_set_crs(32618)
#' alpha_map <- map_grid_estimates(
#'   results_df    = alpha_df,
#'   results_id_col = "plot_id",
#'   crowns_sf     = crowns,
#'   cellsize      = 100
#' )
#'
#' # Convert to raster
#' library(terra)
#' r <- rast(vect(alpha_map), resolution = 0.001)
#' r_alpha <- rasterize(vect(alpha_map), r, field = "mean")
#' }
#'
#' @export
map_grid_estimates <- function(results_df,
                               results_id_col = "plot_id",
                               plots_sf       = NULL,
                               plots_id_col   = "plot_id",
                               crowns_sf      = NULL,
                               crowns_id_col  = NULL,
                               cellsize       = 100) {

  # ── Input checks ────────────────────────────────────────────────────────────
  if (!is.data.frame(results_df))
    stop("`results_df` must be a data frame.")
  if (!results_id_col %in% names(results_df))
    stop(paste0("Column '", results_id_col, "' not found in `results_df`."))
  if (is.null(plots_sf) && is.null(crowns_sf))
    stop("Supply either `plots_sf` (pre-existing polygons) or `crowns_sf` (to build a grid).")

  # ── Step 1: Get or build plot polygons ──────────────────────────────────────
  if (!is.null(plots_sf)) {

    # Case 1: use pre-existing plot polygons
    if (!inherits(plots_sf, "sf"))
      stop("`plots_sf` must be an sf object.")
    if (!plots_id_col %in% names(plots_sf))
      stop(paste0("Column '", plots_id_col, "' not found in `plots_sf`."))

    grid_sf <- plots_sf %>%
      dplyr::rename(.plot_id = !!rlang::sym(plots_id_col))

  } else {

    # Case 2: build a regular grid from crown extent
    if (!inherits(crowns_sf, "sf"))
      stop("`crowns_sf` must be an sf object.")
    if (is.na(sf::st_crs(crowns_sf)))
      stop("`crowns_sf` has no CRS. Set one with sf::st_set_crs().")

    grid_raw <- sf::st_make_grid(crowns_sf, cellsize = cellsize, square = TRUE) %>%
      sf::st_sf(.plot_id = seq_along(.), crs = sf::st_crs(crowns_sf))

    if (!is.null(crowns_id_col)) {

      # Crown polygons already know which plot they belong to —
      # assign each cell the majority plot ID of crowns within it
      if (!crowns_id_col %in% names(crowns_sf))
        stop(paste0("Column '", crowns_id_col, "' not found in `crowns_sf`."))

      centroids <- crowns_sf %>%
        sf::st_centroid() %>%
        dplyr::select(crown_plot_id = !!rlang::sym(crowns_id_col))

      cell_majority <- sf::st_join(centroids, grid_raw, join = sf::st_within) %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(!is.na(.plot_id), !is.na(crown_plot_id)) %>%
        dplyr::count(.plot_id, crown_plot_id) %>%
        dplyr::group_by(.plot_id) %>%
        dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(.plot_id, .plot_id_majority = crown_plot_id)

      grid_sf <- grid_raw %>%
        dplyr::left_join(cell_majority, by = ".plot_id") %>%
        dplyr::mutate(.plot_id = .plot_id_majority) %>%
        dplyr::select(-.plot_id_majority) %>%
        dplyr::filter(!is.na(.plot_id))

    } else {

      # No crown plot IDs — use the generated grid cell integer as the plot ID
      grid_sf <- grid_raw

    }
  }

  # ── Step 2: Standardise results ID column name and join ─────────────────────
  results_df <- results_df %>%
    dplyr::rename(.plot_id = !!rlang::sym(results_id_col))

  grid_out <- grid_sf %>%
    dplyr::left_join(results_df, by = ".plot_id") %>%
    dplyr::rename(!!rlang::sym(results_id_col) := .plot_id) %>%
    sf::st_transform(crs = 4326)

  return(grid_out)
}

#' Fit Alpha Using Bayesian Pareto Model with LAI and Breakpoint Corrections
#'
#' This function fits a Pareto model to tree size data using a Bayesian approach with Stan,
#' accounting for reduced visibility of small trees and site-level LAI corrections. Note that running
#' this function will require rstan, which itself needs Rtools. Rtools is not a CRAN package and as far
#' as we know must be installed directly. Try this link (https://cran.r-project.org/bin/windows/Rtools/)
#'
#' @param filtered_data A list returned by \code{truncate_filter()}, containing the filtered
#'   tree size data frame (\code{bayesian_data}), the final breakpoint (\code{final_breakpoint}),
#'   and the log-scaled KDE data frame (\code{kerneldens_logtransform}).
#' @param LAI A numeric value for site-level Leaf Area Index. The function assumes your LAI value is on
#'   a 1-10 scale and will divide by 10 to get LAI between 0 and 1 (e.g., if you have an LAI of 5 on a
#'   scale of 1-10, enter 5, while if you have a 50 on a scale of 1-100, enter 5).
#' @param prior_mean Prior mean for the alpha parameter.
#' @param prior_sd Prior standard deviation for the alpha parameter.
#' @param stan_model_code Character string of a Stan model (default = \code{stan_alpha_model}).
#' @param iter Number of Stan iterations (default = 9000).
#' @param warmup Number of warmup iterations (default = 6000).
#' @param chains Number of chains (default = 4).
#' @param cores Number of CPU cores to use for parallel chains (default = 1).
#' @param refresh Frequency of Stan progress output (default = 0, i.e., silent).
#'
#' @return A list with:
#' \describe{
#'   \item{posterior_summary}{Data frame of posterior summaries for alpha with R2.}
#'   \item{stan_fit}{Stan fit object.}
#'   \item{breakpoint_norm}{Normalized breakpoint value passed to downstream functions.}
#'   \item{LAI_norm}{Normalized LAI value passed to downstream functions.}
#'   \item{bayesian_data}{The filtered data frame used for model fitting.}
#' }
#' @examples
#' \dontrun{
#' data("harv_data_sample")
#' df <- harv_data_sample[harv_data_sample$IDhectbest == 1, "dbh", drop = FALSE]
#' kde_output   <- potential_break(df)
#' trunc_output <- truncate_filter(kde_output)
#' fit <- fit_alpha_model(trunc_output, LAI = 5.4, prior_mean = 1.4, prior_sd = 0.3)
#' }
#' @export
fit_alpha_model <- function(filtered_data,
                            LAI,
                            prior_mean,
                            prior_sd,
                            stan_model_code = stan_alpha_model,
                            iter = 9000,
                            warmup = 6000,
                            chains = 4,
                            cores = 1,
                            refresh = 0) {
  bayesian_data<-filtered_data$bayesian_data
  breakpoint <- filtered_data$final_breakpoint
  bootstrap_kde_log <- filtered_data$kerneldens_logtransform

  stopifnot(is.data.frame(bayesian_data),
            "dbh" %in% names(bayesian_data),
            is.numeric(LAI),
            is.numeric(prior_mean),
            is.numeric(prior_sd))

  x_vals <- bayesian_data[["dbh"]]

  if (length(x_vals) < 25) stop("Not enough data to fit model.")

  trunc_point <- 10^breakpoint
  trunc_upper <- max(x_vals)
  breakpoint_norm <- (trunc_point - 10) / (min(trunc_upper, 50) - 10)
  breakpoint_norm <- ifelse(breakpoint_norm > 0, breakpoint_norm, 0)

  stan_data <- list(
    N = length(x_vals),
    x_min = 10,
    trunc_point = trunc_point,
    trunc_upper = trunc_upper,
    x = x_vals,
    LAI_norm = LAI / 10,
    breakpoint_norm = breakpoint_norm,
    prior_mean = prior_mean,
    prior_sd = prior_sd
  )

  stan_fit <- rstan::sampling(
    object = rstan::stan_model(model_code = stan_model_code),
    data = stan_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    refresh = refresh
  )

  summary_df <- posterior::summarise_draws(stan_fit) %>%
    dplyr::filter(variable == "alpha") %>%
    dplyr::mutate(
      R2_kernel = performance::r2(stats::lm(mean_log_density ~ log_x, bootstrap_kde_log))$R2
    )

  return(list(
    posterior_summary = summary_df,
    stan_fit = stan_fit,
    breakpoint_norm = breakpoint_norm,
    LAI_norm = LAI / 10,
    bayesian_data = bayesian_data
  ))
}
