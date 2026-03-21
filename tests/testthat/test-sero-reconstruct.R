test_that("sero_reconstruct returns seroreconstruct_fit", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  expect_s3_class(fit, "seroreconstruct_fit")
  expect_true("posterior_model_parameter" %in% names(fit))
  expect_true("posterior_baseline_HAI_titer" %in% names(fit))
  expect_true("posterior_inf_status" %in% names(fit))
  expect_true("data" %in% names(fit))
})

test_that("sero_reconstruct fit has correct attributes", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  expect_equal(attr(fit, "n_groups"), 3L)
  expect_equal(attr(fit, "n_individuals"), nrow(inputdata))
  expect_equal(attr(fit, "n_seasons"), 1L)
  expect_true(is.numeric(attr(fit, "runtime_secs")))
})

test_that("sero_reconstruct posterior has correct dimensions", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  n_iter <- 200
  burnin <- 100
  thinning <- 1
  expected_samples <- (n_iter - burnin) / thinning

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = n_iter, burnin = burnin,
                          thinning = thinning)

  post <- fit$posterior_model_parameter
  expect_equal(nrow(post), expected_samples)
  # Single season: 10 active params
  expect_equal(ncol(post), 10)
})

test_that("sero_reconstruct works with HAI_titer_3 column name", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  # Rename to underscore variant
  renamed <- inputdata
  colnames(renamed)[colnames(renamed) == "HAI_titer3"] <- "HAI_titer_3"

  fit <- sero_reconstruct(renamed, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)
  expect_s3_class(fit, "seroreconstruct_fit")
})

test_that("sero_reconstruct multi-season produces correct n_seasons", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  d1 <- inputdata[1:100, ]
  d1$season <- 0L
  d2 <- inputdata[101:200, ]
  d2$season <- 1L
  multi <- rbind(d1, d2)

  fit <- sero_reconstruct(multi, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  expect_s3_class(fit, "seroreconstruct_fit")
  expect_equal(attr(fit, "n_seasons"), 2L)
  # 2 seasons: 6 shared + 3*2 inf_prob + 2 hai_coef = 14 active params
  expect_equal(ncol(fit$posterior_model_parameter), 14)
})
