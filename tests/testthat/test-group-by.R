test_that("group_by returns seroreconstruct_multi", {
  expect_s3_class(test_fit_multi, "seroreconstruct_multi")
  expect_equal(length(attr(test_fit_multi, "group_labels")), 3)
})

test_that("group_by individual fits are seroreconstruct_fit", {
  labels <- attr(test_fit_multi, "group_labels")
  for (g in labels) {
    expect_s3_class(test_fit_multi[[g]], "seroreconstruct_fit")
    expect_equal(attr(test_fit_multi[[g]], "n_groups"), 1L)
  }
})

test_that("print.seroreconstruct_multi works", {
  out <- capture.output(print(test_fit_multi))
  expect_true(any(grepl("multi-group", out)))
  expect_true(any(grepl("Groups:", out)))
})

test_that("summary.seroreconstruct_multi returns combined table", {
  s <- summary(test_fit_multi)
  expect_s3_class(s, "summary.seroreconstruct_multi")
  expect_true("table" %in% names(s))

  # 3 groups x 6 rows (single-group summary) = 18 rows
  expect_equal(nrow(s$table), 18)
  expect_true("Group" %in% names(s$table))
})

test_that("[[ accessor works for seroreconstruct_multi", {
  # Access by index
  expect_s3_class(test_fit_multi[[1]], "seroreconstruct_fit")

  # Access by name
  labels <- attr(test_fit_multi, "group_labels")
  expect_s3_class(test_fit_multi[[labels[1]]], "seroreconstruct_fit")

  # Access by bad name errors
  expect_error(test_fit_multi[["nonexistent"]], "not found")
})

test_that("group_by errors on too-small groups", {
  small <- inputdata[1:15, ]
  small$age_group <- c(rep(0L, 5), rep(1L, 10))

  expect_error(
    sero_reconstruct(small, flu_activity,
                     n_iteration = 200, burnin = 100, thinning = 1,
                     group_by = ~age_group),
    "< 10 minimum"
  )
})
