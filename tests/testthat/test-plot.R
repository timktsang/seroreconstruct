test_that("plot_diagnostics runs without error", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  expect_no_error(plot_diagnostics(fit))
})

test_that("plot_diagnostics works with param subset", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  expect_no_error(plot_diagnostics(fit, params = c("random_error", "twofold_error")))
})

test_that("plot_diagnostics errors on invalid param names", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  expect_error(plot_diagnostics(fit, params = c("nonexistent")), "Unknown")
})

test_that("plot_diagnostics errors for multi", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1,
                          group_by = ~age_group)

  expect_error(plot_diagnostics(fit), "not supported")
})

# ---- plot_trajectory tests ----

test_that("plot_trajectory runs without error", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  expect_no_error(plot_trajectory(fit, id = 1))
})

test_that("plot_trajectory works with custom arguments", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  expect_no_error(plot_trajectory(fit, id = 50, n_samples = 20,
                                  main = "Custom title", show_legend = FALSE))
})

test_that("plot_trajectory errors on invalid id", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  expect_error(plot_trajectory(fit, id = 0), "must be an integer")
  expect_error(plot_trajectory(fit, id = 999999), "must be an integer")
  expect_error(plot_trajectory(fit, id = 1.5), "must be an integer")
})

test_that("plot_trajectory errors for multi", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1,
                          group_by = ~age_group)

  expect_error(plot_trajectory(fit), "not supported")
})

test_that("plot_trajectory works in multi-panel layout", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  par(mfrow = c(1, 2))
  expect_no_error({
    plot_trajectory(fit, id = 1, show_legend = FALSE)
    plot_trajectory(fit, id = 50, show_legend = FALSE)
  })
  # Verify mfrow was preserved (not reset)
  expect_equal(par("mfrow"), c(1, 2))
  par(mfrow = c(1, 1))
})
