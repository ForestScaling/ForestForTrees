#' @importFrom dplyr %>%
#' @importFrom rlang %||%
NULL

utils::globalVariables(c("dbh", "log_x", "x", "trunc_output", "variable"))

#' Stan model string for estimating Ntotal with LAI and breakpoint corrections
#'
#' This is the default Stan model used internally by [estimate_total_trees()].
#' @format A character string containing the full Stan model code.
#' @export
stan_totaltrees_model <- "data {
  int<lower=1> K;                // Number of bins
  real<lower=0> bin_min[K];      // Lower bounds of DBH for bins
  real<lower=0> bin_max[K];      // Upper bounds of DBH for bins
  real<lower=0> N_obs[K];        // Observed counts per bin (now continuous)
  real<lower=0> x_min;           // Minimum DBH threshold
  real<lower=0> x_max;           // Maximum DBH threshold (e.g., 50)
  int<lower=1> n_alpha_samples;  // Number of pre-sampled alphas
  real alpha_samples[n_alpha_samples]; // Pre-sampled alpha values
  real<lower=0> N_tot_prior_mean; // Mean of the prior for N_tot
  real<lower=0> N_tot_prior_sd;   // SD of the prior for N_tot
  real<lower=0, upper=1> breakpoint_norm; // Normalized breakpoint value (higher = less confident)
  real<lower=0, upper=1> LAI_norm;        // Normalized LAI value (higher = less confident)
}
parameters {
  real<lower=0> N_tot;           // Total number of trees with DBH >= x_min
  real<lower=0> sigma;           // Base standard deviation for normal likelihood
}
model {
  vector[K] lambda;               // Expected counts for each bin
  real log_lik[n_alpha_samples];  // Log-likelihoods for each alpha sample
  // Compute the adjustment factor
  real adjustment_factor = 1 - sqrt(LAI_norm * breakpoint_norm);
  real adjusted_N_tot = N_tot / adjustment_factor;
  // Prior for N_tot
  N_tot ~ normal(N_tot_prior_mean, N_tot_prior_sd);
  // Prior for sigma
  sigma ~ normal(0, 1);  // Adjust this based on the scale of your data
  // Marginalize over alpha uncertainty
  for (j in 1:n_alpha_samples) {
    real alpha_sample;
    real trunc_correction;
    alpha_sample = alpha_samples[j];  // Use pre-sampled alpha
    // Compute the normalization constant for the truncated power-law distribution
    trunc_correction = (x_min^(1 - alpha_sample) - x_max^(1 - alpha_sample));
    // Calculate expected counts for each bin for this alpha
    for (k in 1:K) {
      lambda[k] = adjusted_N_tot * (bin_min[k]^(1 - alpha_sample) - bin_max[k]^(1 - alpha_sample)) /
                  trunc_correction;  // Normalized lambda
    }
    // Compute the log-likelihood for this alpha
    log_lik[j] = 0;
    for (k in 1:K) {
      real sigma_k = sigma;  // No crown score adjustment
      log_lik[j] += normal_lpdf(N_obs[k] | lambda[k], sigma_k);
    }
  }
  // Combine likelihoods across alpha samples
  target += log_sum_exp(log_lik) - log(n_alpha_samples);
} "


#' Estimate Total Number of Trees Using Fitted Alpha Model and Stan
#'
#' This function fits a Stan model to estimate the total number of trees in a site,
#' using output from \code{fit_alpha_model()}, which includes posterior samples, LAI and
#' breakpoint normalization, and raw DBH data.
#'
#' @param alpha_model_output A named list returned by \code{fit_alpha_model()}.
#' @param siteID Optional site ID tag for output (default = \code{NULL}).
#' @param DBH_max Optional max DBH value (default = max of site data).
#' @param N_tot_prior_mean Prior mean for N_tot (default = 1250).
#' @param N_tot_prior_sd Prior SD for N_tot (default = 625).
#' @param chains Number of MCMC chains (default = 4).
#' @param warmup Warmup iterations (default = 6000).
#' @param iter Total iterations (default = 9000).
#' @param refresh Refresh rate for Stan progress output (default = 0, i.e., silent).
#' @param num_alpha_samples Number of posterior alpha samples to use (default = 500).
#'   Reducing this can speed up runtime at a small cost to precision.
#'
#' @return A list with:
#' \describe{
#'   \item{posterior_summary}{Data frame of posterior summaries for N_tot.}
#'   \item{stan_fit}{Full Stan fit object.}
#' }
#' @examples
#' \dontrun{
#' data("harv_data_sample")
#' df <- harv_data_sample[harv_data_sample$IDhectbest == 1, "dbh", drop = FALSE]
#' kde_output   <- potential_break(df)
#' trunc_output <- truncate_filter(kde_output)
#' fit          <- fit_alpha_model(trunc_output, LAI = 5.4, prior_mean = 1.4, prior_sd = 0.3)
#' trees        <- estimate_total_trees(fit, N_tot_prior_mean = 1250, N_tot_prior_sd = 625)
#' }
#' @export
estimate_total_trees <- function(alpha_model_output,
                                 siteID           = NULL,
                                 DBH_max          = NULL,
                                 N_tot_prior_mean = 1250,
                                 N_tot_prior_sd   = 625,
                                 chains           = 4,
                                 warmup           = 6000,
                                 iter             = 9000,
                                 refresh          = 0,
                                 num_alpha_samples = 500) {

  # Unpack from alpha_model_output
  stopifnot(all(c("posterior_summary", "stan_fit", "breakpoint_norm",
                  "LAI_norm", "bayesian_data") %in% names(alpha_model_output)))

  site_data <- alpha_model_output$bayesian_data
  LAI_norm <- alpha_model_output$LAI_norm
  breakpoint_norm <- alpha_model_output$breakpoint_norm

  # Extract alpha samples (capped to available posterior draws)
  all_alpha_samples <- rstan::extract(alpha_model_output$stan_fit, "alpha")[[1]]
  alpha_samples <- all_alpha_samples[1:min(num_alpha_samples, length(all_alpha_samples))]

  # Validate
  stopifnot("dbh" %in% names(site_data),
            is.numeric(alpha_samples))

  if (nrow(site_data) < 25) stop("Not enough data to estimate total trees.")

  # Convert normalized breakpoint to DBH scale
  trunc_point <- 10 + breakpoint_norm * (min(DBH_max %||% max(site_data$dbh, na.rm = TRUE), 50) - 10)
  site_data <- dplyr::filter(site_data, dbh >= trunc_point)

  if (nrow(site_data) < 25) stop("Not enough data after truncation.")

  # Compute DBH max
  DBH_max <- if (is.null(DBH_max)) max(site_data$dbh, na.rm = TRUE) else DBH_max

  # Log-bin the data
  binned_data <- logbin(site_data$dbh, n = 8) %>%
    tidyr::drop_na()

  if (nrow(binned_data) < 2) stop("Binning failed or produced too few bins.")

  # Prepare Stan data
  stan_data <- list(
    K = nrow(binned_data),
    bin_min = binned_data$bin_min,
    bin_max = binned_data$bin_max,
    N_obs = binned_data$bin_count,
    x_min = 10,
    x_max = max(site_data$dbh),
    n_alpha_samples = length(alpha_samples),
    alpha_samples = alpha_samples,
    N_tot_prior_mean = N_tot_prior_mean,
    N_tot_prior_sd = N_tot_prior_sd,
    breakpoint_norm = breakpoint_norm,
    LAI_norm = LAI_norm
  )

  # Fit Stan model
  stan_fit <- rstan::sampling(
    object = rstan::stan_model(model_code = stan_totaltrees_model),
    data   = stan_data,
    chains = chains,
    warmup = warmup,
    iter   = iter,
    refresh = refresh
  )

  summary_df <- posterior::summarise_draws(stan_fit) %>%
    dplyr::filter(variable == "N_tot") %>%
    dplyr::mutate(site = siteID %||% NA_character_,
                  breakpoint = trunc_point)

  return(list(
    posterior_summary = summary_df,
    stan_fit = stan_fit
  ))
}
