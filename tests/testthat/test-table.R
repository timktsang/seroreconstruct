test_that("table_parameters returns correct structure", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  tbl <- table_parameters(fit)
  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), 10)  # 3-group single-season: 10 active params
  expect_true(all(c("Parameter", "Mean", "Median", "Lower", "Upper") %in% names(tbl)))
})

test_that("table_parameters has correct parameter names for 3-group", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  tbl <- table_parameters(fit)
  expect_true("random_error" %in% tbl$Parameter)
  expect_true("inf_prob_children" %in% tbl$Parameter)
  expect_true("inf_prob_adults" %in% tbl$Parameter)
  expect_true("hai_coef" %in% tbl$Parameter)
})

test_that("table_parameters works with custom probs", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  tbl <- table_parameters(fit, probs = c(0.1, 0.9))
  expect_true(all(tbl$Lower <= tbl$Upper))
})

test_that("table_parameters works for seroreconstruct_multi", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1,
                          group_by = ~age_group)

  tbl <- table_parameters(fit)
  expect_true("Group" %in% names(tbl))
  # 3 groups x some params each
  expect_true(nrow(tbl) > 10)
})

test_that("table_infections returns correct structure", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  tbl <- table_infections(fit)
  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), nrow(inputdata))
  expect_true(all(c("Individual", "Infection_prob", "Infection_time_mean",
                     "Baseline_titer_mean") %in% names(tbl)))
})

test_that("table_infections has valid probability values", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  tbl <- table_infections(fit)
  expect_true(all(tbl$Infection_prob >= 0 & tbl$Infection_prob <= 1))
})

test_that("table_infections errors for multi", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1,
                          group_by = ~age_group)

  expect_error(table_infections(fit), "not supported")
})

test_that(".get_param_names returns correct count", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  names <- seroreconstruct:::.get_param_names(fit)
  expect_equal(length(names), ncol(fit$posterior_model_parameter))
})
