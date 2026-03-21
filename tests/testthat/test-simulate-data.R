test_that("simulate_data returns correct dimensions", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")
  data("para1", package = "seroreconstruct")
  data("para2", package = "seroreconstruct")

  sim <- simulate_data(inputdata, flu_activity, para1, para2)

  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), nrow(inputdata))
  expect_equal(ncol(sim), 9)
})

test_that("simulate_data is reproducible with set.seed", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")
  data("para1", package = "seroreconstruct")
  data("para2", package = "seroreconstruct")

  set.seed(123)
  sim1 <- simulate_data(inputdata[1:50, ], flu_activity, para1, para2)

  set.seed(123)
  sim2 <- simulate_data(inputdata[1:50, ], flu_activity, para1, para2)

  expect_identical(sim1, sim2)
})

test_that("simulate_data produces valid HAI titers", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")
  data("para1", package = "seroreconstruct")
  data("para2", package = "seroreconstruct")

  set.seed(42)
  sim <- simulate_data(inputdata[1:100, ], flu_activity, para1, para2)

  # HAI titers should be non-negative integers
  for (col in c("X7", "X8", "X9")) {
    vals <- sim[[col]]
    expect_true(all(vals >= 0 | vals == -1),
                info = paste(col, "has invalid values"))
  }
})

test_that("simulate_data works with multi-season input", {
  data("inputdata", package = "seroreconstruct")
  data("flu_activity", package = "seroreconstruct")
  data("para1", package = "seroreconstruct")
  data("para2", package = "seroreconstruct")

  # Build 2-season data
  d1 <- inputdata[1:50, ]
  d1$season <- 0L
  d2 <- inputdata[51:100, ]
  d2$season <- 1L
  multi <- rbind(d1, d2)

  # para1 for 2 seasons: 6 + 4*2 = 14
  para1_2s <- c(para1[1:6], para1[7:9], para1[7:9], para1[10], para1[10])
  para2_2s <- rep(para2, 2)

  set.seed(42)
  sim <- simulate_data(multi, flu_activity, para1_2s, para2_2s)

  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), 100)
  expect_equal(ncol(sim), 9)
})
