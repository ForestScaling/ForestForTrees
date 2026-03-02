#' @importFrom dplyr %>%
NULL

utils::globalVariables(c("dbh", "log_x", "x", "trunc_output", "variable", "peaks"))


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
#' @param data A data frame containing individual measurements (e.g., trees).
#' @param n_bootstrap Number of bootstrap replicates for KDE smoothing (default = 1000).
#' @param bandwidth Bandwidth method passed to `density()` (default = "SJ").
#' @param trim_max Upper size limit.
#'
#' @return A list with:
#' \describe{
#'   \item{potential_breakpoint}{A numeric value representing the estimated lower log10 size bound breakpoint.}
#'   \item{bootstrap_kde_log}{A data frame with log-scaled KDE values (`log_x`, `mean_log_density`).}
#'   \item{original_data_trimmed}{The original size vector *after* trimming by `trim_max`. This is crucial for quantile calculation in the first function and for the final filtering in the second.}
#'   \item{original_raw_data_df}{The original input data frame, untouched, so the second function can filter it.}
#'   \item{trim_max_value}{The `trim_max` value used in this function, passed along for consistency.}
#' }
#' @export
potential_break <- function(data,
                                             n_bootstrap = 1000,
                                             bandwidth = "SJ",
                                             trim_max = 50) {
  # Input validation
  # stopifnot(is.numeric(data[["dbh"]]))

  # Store original raw data for later use in the second function
  original_raw_data_df <- data
  trim_max_value <- trim_max # Store this as well

  # Pull and trim size vector
  size_vector_trimmed <- data[is.na(data) == FALSE]#[["dbh"]]
  size_vector_trimmed <- size_vector_trimmed[size_vector_trimmed <= trim_max]
  if (length(size_vector_trimmed) < 25) {
    stop("Not enough observations below trim_max for breakpoint estimation (min 25 required).")
  }

  # KDE smoothing
  kde <- stats::density(size_vector_trimmed, bw = bandwidth)
  kde_df <- data.frame(x = kde$x, y = kde$y)
  observed_vals <- sort(unique(size_vector_trimmed)) # Used for filtering KDE results

  # filter KDE results to observed data range
  kde_df <- kde_df %>%
    dplyr::filter(x >= min(size_vector_trimmed), x <= max(size_vector_trimmed))

  # Further  filter KDE to points near observed values
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
  peak_df <- cbind(
    bootstrap_kde_log,
    peaks = splus2R::peaks(bootstrap_kde_log$mean_log_density))%>%
    dplyr::filter(peaks == TRUE)


  if (nrow(peak_df) == 0) {
    stop("No distinct peaks detected in the KDE curve for breakpoint identification.")
  }

  # Calculate the candidate breakpoint
  potential_breakpoint <- peak_df %>%
    dplyr::filter(log_x <= stats::quantile(log10(size_vector_trimmed), 0.75)) %>%
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

#' Determine Truncation Points and filter Data for Power-Law Modeling
#'
#' This function takes the potential breakpoint and KDE data from a previous step,
#' performs segmented regression, determines the upper truncation point,
#' and filters the original dataset to the identified size range.
#'
#' @param breakpoint_kde_results A list returned by `potential_break`.
#'   It must contain `potential_breakpoint`, `bootstrap_kde_log`,
#'   `original_raw_data_df`, and `trim_max_value`.
#' @param min_size The mimimum tree DBH that will be in the final distribution. This is the lower bound to
#' the Pareto distribution, not the breakpoint for obscured vs unobscured trees. As such, this is usually approx
#' 10 cm.
#'
#' @return A list with:
#' \describe{
#'   \item{final_breakpoint}{The final determined lower log10 size bound.}
#'   \item{bayesian_data}{Subset of original data between 10^breakpoint and 10^upper_bound.}
#'   \item{kerneldens_logtransform}{The full log-scaled KDE data from the first function.}
#' }
#' @export
truncate_filter <- function(breakpoint_kde_results, min_size = 10) {

  # Extract necessary components from the input list
  potential_breakpoint <- breakpoint_kde_results$potential_breakpoint
  bootstrap_kde_log <- breakpoint_kde_results$bootstrap_kde_log
  original_data_trimmed_df <- data.frame(dbh=breakpoint_kde_results$original_data_trimmed)
  trim_max <- breakpoint_kde_results$trim_max_value

  # For robustness, let's ensure potential_breakpoint is not NA
  if (is.na(potential_breakpoint)) {
    stop("Potential breakpoint is NA. Cannot proceed with segmented regression. Check 'potential_break' output.")
  }

  # filter the data for segmented regression *before* passing it to lm
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
    left_x = utils::head(segment_boundaries, -1),     # Left boundary of each segment
    right_x = utils::tail(segment_boundaries, -1),   # Right boundary of each segment
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

#' Fit Alpha Using Bayesian Pareto Model with LAI and Breakpoint Corrections
#'
#' This function fits a Pareto model to tree size data using a Bayesian approach with Stan,
#' accounting for reduced visibility of small trees and site-level LAI corrections. Note that running
#' this function will require rstan, which itself needs Rtools. Rtools is not a CRAN package and as far
#' as we know must be installed directly. Try this link (https://cran.r-project.org/bin/windows/Rtools/)
#'
#' @param bayesian_data A data frame of filtered tree sizes (from `truncate_filter()`).
#' @param breakpoint The final lower log10 size threshold (from KDE peak).
#' @param LAI A numeric value for site-level Leaf Area Index. The function assumes your LAI value is on
#' a 1-10 scale and will divide by 10 to get LAI between 0 and 1 (e.g., if you have an LAI of 5 on a scale of 1-10,
#' enter 5, while if you have a 50 on a scale of 1-100, enter 5)
#' the LAI to be between 0 and 1
#' @param prior_mean Prior mean for the alpha parameter.
#' @param prior_sd Prior standard deviation for the alpha parameter.
#' @param stan_model_code Character string of a Stan model (default = `stan_alpha_model`).
#' @param iter Number of Stan iterations (default = 9000).
#' @param warmup Number of warmup iterations (default = 6000).
#' @param chains Number of chains (default = 4).
#' @param refresh Frequency of Stan progress output (default = 0).
#' @param cores Number of CPU cores to use (default = 1).
#'
#' @return A list with:
#' \describe{
#'   \item{posterior_summary}{Data frame of posterior summaries for alpha with R2.}
#'   \item{stan_fit}{Stan fit object.}
#' }
#' @export
fit_alpha_model <- function(bayesian_data,
                            breakpoint,
                            LAI,
                            prior_mean,
                            prior_sd,
                            stan_model_code = stan_alpha_model,
                            iter = 9000,
                            warmup = 6000,
                            chains = 4,
                            cores = 1,
                            refresh = 0) {

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
    cores = cores
  )

  summary_df <- posterior::summarise_draws(stan_fit) %>%
    dplyr::filter(variable == "alpha") %>%
    dplyr::mutate(
      R2_kernel = performance::r2(stats::lm(mean_log_density ~ log_x, trunc_output$kerneldens_logtransform))$R2
    )

  return(list(
    posterior_summary = summary_df,
    stan_fit = stan_fit,
    breakpoint_norm = breakpoint_norm,
    LAI_norm = LAI / 10,
    bayesian_data = bayesian_data
  ))
}

