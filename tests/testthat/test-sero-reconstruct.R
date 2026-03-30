test_that("sero_reconstruct returns seroreconstruct_fit", {
  expect_s3_class(test_fit, "seroreconstruct_fit")
  expect_true("posterior_model_parameter" %in% names(test_fit))
  expect_true("posterior_baseline_HAI_titer" %in% names(test_fit))
  expect_true("posterior_inf_status" %in% names(test_fit))
  expect_true("data" %in% names(test_fit))
})

test_that("sero_reconstruct fit has correct attributes", {
  expect_equal(attr(test_fit, "n_groups"), 3L)
  expect_equal(attr(test_fit, "n_individuals"), nrow(inputdata))
  expect_equal(attr(test_fit, "n_seasons"), 1L)
  expect_true(is.numeric(attr(test_fit, "runtime_secs")))
})

test_that("sero_reconstruct posterior has correct dimensions", {
  expected_samples <- (200 - 100) / 1

  post <- test_fit$posterior_model_parameter
  expect_equal(nrow(post), expected_samples)
  # Single season: 10 active params
  expect_equal(ncol(post), 10)
})

test_that("sero_reconstruct works with HAI_titer_3 column name", {
  renamed <- inputdata
  colnames(renamed)[colnames(renamed) == "HAI_titer3"] <- "HAI_titer_3"

  fit <- sero_reconstruct(renamed, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)
  expect_s3_class(fit, "seroreconstruct_fit")
})

test_that("sero_reconstruct multi-season produces correct n_seasons", {
  expect_s3_class(test_fit_ms, "seroreconstruct_fit")
  expect_equal(attr(test_fit_ms, "n_seasons"), 2L)
  # 2 seasons: 6 shared + 3*2 inf_prob + 2 hai_coef = 14 active params
  expect_equal(ncol(test_fit_ms$posterior_model_parameter), 14)
})
