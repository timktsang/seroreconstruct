test_that("plot_diagnostics runs without error", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot_diagnostics(test_fit))
})

test_that("plot_diagnostics works with param subset", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot_diagnostics(test_fit, params = c("random_error", "twofold_error")))
})

test_that("plot_diagnostics errors on invalid param names", {
  expect_error(plot_diagnostics(test_fit, params = c("nonexistent")), "Unknown")
})

test_that("plot_diagnostics errors for multi", {
  expect_error(plot_diagnostics(test_fit_multi), "not supported")
})

# ---- plot_trajectory tests ----

test_that("plot_trajectory runs without error", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot_trajectory(test_fit, id = 1))
})

test_that("plot_trajectory works with custom arguments", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot_trajectory(test_fit, id = 50, n_samples = 20,
                                  main = "Custom title", show_legend = FALSE))
})

test_that("plot_trajectory errors on invalid id", {
  expect_error(plot_trajectory(test_fit, id = 0), "must be an integer")
  expect_error(plot_trajectory(test_fit, id = 999999), "must be an integer")
  expect_error(plot_trajectory(test_fit, id = 1.5), "must be a single integer")
})

test_that("plot_trajectory errors for multi", {
  expect_error(plot_trajectory(test_fit_multi), "not supported")
})

# ---- plot_trajectory subject_ids tests ----

test_that("plot_trajectory works with character subject_ids", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  # Character lookup should work
  expect_no_error(plot_trajectory(test_fit_ids, id = "S1"))
  # Row index should still work
  expect_no_error(plot_trajectory(test_fit_ids, id = 1))
})

test_that("plot_trajectory works with numeric subject_ids out of range", {
  n <- nrow(inputdata)
  num_ids <- 10000L + seq_len(n)

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1,
                          subject_ids = num_ids)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  # Numeric ID out of row range should trigger lookup
  expect_no_error(plot_trajectory(fit, id = 10001))
  # Row index should still work
  expect_no_error(plot_trajectory(fit, id = 1))
})

test_that("plot_trajectory works with subjects argument", {
  n <- nrow(inputdata)
  ext_ids <- paste0("EXT", seq_len(n))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  # subjects argument should enable lookup
  expect_no_error(plot_trajectory(test_fit, id = "EXT1", subjects = ext_ids))
})

test_that("plot_trajectory errors on nonexistent subject_id", {
  expect_error(plot_trajectory(test_fit_ids, id = "NONEXISTENT"), "not found")
})

test_that("plot_trajectory errors on out-of-range numeric without subject_ids", {
  # No subject_ids → out of range numeric should error
  expect_error(plot_trajectory(test_fit, id = 99999), "must be an integer")
})

test_that("plot_trajectory works in multi-panel layout", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2))
  expect_no_error({
    plot_trajectory(test_fit, id = 1, show_legend = FALSE)
    plot_trajectory(test_fit, id = 50, show_legend = FALSE)
  })
  # Verify mfrow was preserved (not reset)
  expect_equal(par("mfrow"), c(1, 2))
  par(mfrow = c(1, 1))
})
