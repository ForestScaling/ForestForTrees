test_that("truncate_filter runs with valid input and returns expected structure", {
  set.seed(42)
  # Create synthetic data and run the first function
  test_data <- data.frame(dbh = c(
    rlnorm(500, meanlog = log(15), sdlog = 0.3),
    rlnorm(500, meanlog = log(35), sdlog = 0.2)
  ))

  kde_results <- potential_break(
    data = test_data$dbh,
    n_bootstrap = 50,
    bandwidth = "SJ",
    trim_max = 50
  )

  # Run the second function
  truncation_results <- truncate_filter(kde_results)

  # Structure and type checks
  expect_type(truncation_results, "list")
  expect_named(
    truncation_results,
    c("bayesian_data", "kerneldens_logtransform", "final_breakpoint")
  )

  expect_s3_class(truncation_results$bayesian_data, "data.frame")
  expect_s3_class(truncation_results$kerneldens_logtransform, "data.frame")
  expect_type(truncation_results$final_breakpoint, "double")

  # Data sanity checks
  expect_true(all(truncation_results$bayesian_data$dbh >= 10^truncation_results$final_breakpoint))
  expect_true(all(truncation_results$kerneldens_logtransform$log_x >= truncation_results$final_breakpoint))
  expect_true(all(is.finite(truncation_results$kerneldens_logtransform$mean_log_density)))
})

test_that("function stops with NA breakpoint", {
  fake_kde_results <- list(
    potential_breakpoint = NA,
    bootstrap_kde_log = data.frame(log_x = 1:10, mean_log_density = rnorm(10)),
    original_data_trimmed = runif(50, 10, 50),
    trim_max_value = 50
  )

  expect_error(
    truncate_filter(fake_kde_results),
    "Potential breakpoint is NA"
  )
})

test_that("function stops with insufficient data after filtering", {
  fake_kde_results <- list(
    potential_breakpoint = 2, # large value means little remaining data
    bootstrap_kde_log = data.frame(
      log_x = c(1.9, 1.95),  # only two points, right on the edge
      mean_log_density = c(-2, -2.5)
    ),
    original_data_trimmed = runif(100, 10, 50),
    trim_max_value = 50
  )

  expect_error(
    truncate_filter(fake_kde_results),
    "Not enough data points after breakpoint filtering"
  )
})

test_that("filtered data are within size limits and above min_size", {
  set.seed(123)
  test_data <- data.frame(dbh = rlnorm(300, meanlog = log(20), sdlog = 0.25))

  kde_results <- potential_break(
    data = test_data$dbh,
    n_bootstrap = 30,
    bandwidth = "SJ",
    trim_max = 50
  )

  truncation_results <- truncate_filter(kde_results, min_size = 15)

  expect_true(all(truncation_results$bayesian_data$dbh >= 15))
  expect_true(all(truncation_results$kerneldens_logtransform$log_x >= log10(15)))
})
