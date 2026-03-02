test_that("potential_break runs and returns correctly structured output", {
  set.seed(42)
  # Create test data with a realistic DBH distribution
  test_data <- data.frame(dbh = c(
    rlnorm(500, meanlog = log(15), sdlog = 0.3),
    rlnorm(500, meanlog = log(35), sdlog = 0.2)
  ))

  # Run function
  result <- potential_break(
    data = test_data,
    n_bootstrap = 50,  # smaller for faster testing
    bandwidth = "SJ",
    trim_max = 50
  )

  # Structure tests
  expect_type(result, "list")
  expect_named(
    result,
    c("potential_breakpoint", "bootstrap_kde_log", "original_data_trimmed",
      "original_raw_data_df", "trim_max_value")
  )

  # Check elements
  expect_type(result$potential_breakpoint, "double")
  expect_true(is.numeric(result$original_data_trimmed))
  expect_s3_class(result$bootstrap_kde_log, "data.frame")
  expect_s3_class(result$original_raw_data_df, "data.frame")

  # Check numeric behavior
  expect_true(all(result$original_data_trimmed <= 50))
  expect_true(result$trim_max_value == 50)
  expect_true(all(is.finite(result$bootstrap_kde_log$log_x)))
  expect_true(all(is.finite(result$bootstrap_kde_log$mean_log_density)))
})

