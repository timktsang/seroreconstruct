test_that("print.seroreconstruct_fit works", {
  out <- capture.output(print(test_fit))
  expect_true(any(grepl("seroreconstruct fit", out)))
  expect_true(any(grepl("Individuals:", out)))
  expect_true(any(grepl("summary\\(\\)", out)))
})

test_that("summary.seroreconstruct_fit returns table", {
  s <- summary(test_fit)
  expect_s3_class(s, "summary.seroreconstruct_fit")
  expect_true("table" %in% names(s))
  expect_s3_class(s$table, "data.frame")

  # 3-group fit: 16 rows
  expect_equal(nrow(s$table), 16)
  expect_equal(ncol(s$table), 4)
})

test_that("summary.seroreconstruct_fit with period argument", {
  s <- summary(test_fit, period = c(800, 1100))
  expect_s3_class(s, "summary.seroreconstruct_fit")
  expect_equal(nrow(s$table), 16)
})

test_that("print.summary.seroreconstruct_fit works", {
  s <- summary(test_fit)
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
  expect_true(any(grepl("Random error", out)))
})

test_that("summary errors for multi-season fits", {
  expect_error(summary(test_fit_ms), "not yet implemented")
})

test_that("output_model_estimate gives deprecation warning", {
  expect_warning(output_model_estimate(test_fit), "deprecated")
})
