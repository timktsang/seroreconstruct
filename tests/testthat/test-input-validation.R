test_that("sero_reconstruct rejects non-data.frame input", {
  expect_error(sero_reconstruct("not a df", matrix(1), 100, 50, 1),
               "must be a data frame")
})

test_that("sero_reconstruct rejects missing columns", {
  bad <- data.frame(age_group = 0, start_time = 1)
  expect_error(sero_reconstruct(bad, matrix(1), 100, 50, 1),
               "missing required columns")
})

test_that("sero_reconstruct rejects bad MCMC parameters", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  expect_error(sero_reconstruct(inputdata, flu_activity, -1, 0, 1),
               "positive integer")
  expect_error(sero_reconstruct(inputdata, flu_activity, 100, 200, 1),
               "greater than")
  expect_error(sero_reconstruct(inputdata, flu_activity, 100, 50, -1),
               "positive integer")
})

test_that("sero_reconstruct rejects invalid age_group values", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  bad <- inputdata[1:20, ]
  bad$age_group <- 5L
  expect_error(sero_reconstruct(bad, flu_activity, 100, 50, 1),
               "age_group")
})

test_that("sero_reconstruct rejects non-formula group_by", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  expect_error(sero_reconstruct(inputdata, flu_activity, 100, 50, 1,
                                group_by = "age_group"),
               "must be a formula")
})

test_that("sero_reconstruct rejects group_by with missing variable", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  expect_error(sero_reconstruct(inputdata, flu_activity, 100, 50, 1,
                                group_by = ~nonexistent),
               "not found")
})

test_that("sero_reconstruct rejects ILI shorter than max end_time", {
  data("inputdata", package = "seroreconstruct")

  short_ili <- matrix(0.01, nrow = 10, ncol = 1)
  expect_error(sero_reconstruct(inputdata, short_ili, 100, 50, 1),
               "rows")
})

test_that("simulate_data rejects non-data.frame input", {
  expect_error(simulate_data("not a df", matrix(1), 1:10, 1:20),
               "must be a data frame")
})

test_that("simulate_data rejects wrong para1 length", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")
  data("para2", package = "seroreconstruct")

  expect_error(simulate_data(inputdata, flu_activity, 1:5, para2),
               "length")
})

test_that("simulate_data rejects wrong para2 length", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")
  data("para1", package = "seroreconstruct")

  expect_error(simulate_data(inputdata, flu_activity, para1, 1:5),
               "length")
})

test_that("season column validation works", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")

  bad <- inputdata[1:20, ]
  bad$season <- 1L  # not 0-indexed
  expect_error(sero_reconstruct(bad, flu_activity, 100, 50, 1),
               "0-indexed")

  bad2 <- inputdata[1:20, ]
  bad2$season <- c(rep(0L, 10), rep(2L, 10))  # gap: missing season 1
  expect_error(sero_reconstruct(bad2, flu_activity, 100, 50, 1),
               "contiguous")
})
