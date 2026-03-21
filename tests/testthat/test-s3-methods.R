test_that("print.seroreconstruct_fit works", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  out <- capture.output(print(fit))
  expect_true(any(grepl("seroreconstruct fit", out)))
  expect_true(any(grepl("Individuals:", out)))
  expect_true(any(grepl("summary\\(\\)", out)))
})

test_that("summary.seroreconstruct_fit returns table", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  s <- summary(fit)
  expect_s3_class(s, "summary.seroreconstruct_fit")
  expect_true("table" %in% names(s))
  expect_s3_class(s$table, "data.frame")

  # 3-group fit: 16 rows
  expect_equal(nrow(s$table), 16)
  expect_equal(ncol(s$table), 4)
})

test_that("summary.seroreconstruct_fit with period argument", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  s <- summary(fit, period = c(800, 1100))
  expect_s3_class(s, "summary.seroreconstruct_fit")
  expect_equal(nrow(s$table), 16)
})

test_that("print.summary.seroreconstruct_fit works", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  s <- summary(fit)
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
  expect_true(any(grepl("Random error", out)))
})

test_that("summary errors for multi-season fits", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  d1 <- inputdata[1:100, ]
  d1$season <- 0L
  d2 <- inputdata[101:200, ]
  d2$season <- 1L
  multi <- rbind(d1, d2)

  fit <- sero_reconstruct(multi, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  expect_error(summary(fit), "not yet implemented")
})

test_that("output_model_estimate gives deprecation warning", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  fit <- sero_reconstruct(inputdata, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1)

  expect_warning(output_model_estimate(fit), "deprecated")
})
