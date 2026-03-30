test_that("shared argument validation works", {
  # shared without group_by is ignored (no error)
  fit <- sero_reconstruct(test_inputdata_vax, flu_activity,
                          n_iteration = 200, burnin = 100, thinning = 1,
                          shared = c("error"))
  expect_s3_class(fit, "seroreconstruct_fit")

  # invalid shared values
  expect_error(
    sero_reconstruct(test_inputdata_vax, flu_activity,
                     n_iteration = 200, burnin = 100, thinning = 1,
                     group_by = ~vaccine, shared = c("invalid")),
    "Invalid"
  )

  # shared must be character
 expect_error(
    sero_reconstruct(test_inputdata_vax, flu_activity,
                     n_iteration = 200, burnin = 100, thinning = 1,
                     group_by = ~vaccine, shared = 123),
    "character"
  )
})

test_that("shared joint model returns seroreconstruct_joint", {
  expect_s3_class(test_fit_joint, "seroreconstruct_joint")
  expect_s3_class(test_fit_joint, "seroreconstruct_fit")
  expect_equal(attr(test_fit_joint, "n_groups_joint"), 2L)
  expect_equal(length(attr(test_fit_joint, "group_labels")), 2)
  expect_equal(attr(test_fit_joint, "shared"), c("error", "boosting_waning"))
})

test_that("joint model posterior has correct dimensions", {
  post <- test_fit_joint$posterior_model_parameter
  # 2 groups, 1 season: 6 shared + 2 inf_prob + 1 hai_coef = 9
  expect_equal(ncol(post), 9)
  expect_equal(nrow(post), 100)  # (200 - 100) / 1
})

test_that("print.seroreconstruct_joint works", {
  out <- capture.output(print(test_fit_joint))
  expect_true(any(grepl("joint fit", out)))
  expect_true(any(grepl("Groups:", out)))
  expect_true(any(grepl("Shared:", out)))
})

test_that("summary.seroreconstruct_joint returns group-aware table", {
  s <- summary(test_fit_joint)
  expect_s3_class(s, "summary.seroreconstruct_joint")
  expect_true("table" %in% names(s))

  # 4 shared rows + 2 group-specific infection prob rows = 6
  expect_equal(nrow(s$table), 6)
  expect_equal(ncol(s$table), 4)

  # Check that both groups appear in the output
  vars <- s$table$Variable
  expect_true(any(grepl("0", vars)))
  expect_true(any(grepl("1", vars)))
  expect_true(any(grepl("Random error", vars)))
  expect_true(any(grepl("Infection probability", vars)))
})

test_that("print.summary.seroreconstruct_joint works", {
  s <- summary(test_fit_joint)
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
  expect_true(any(grepl("Shared:", out)))
  expect_true(any(grepl("Groups:", out)))
})

test_that("group_by without shared still gives independent chains", {
  expect_warning(
    fit <- sero_reconstruct(test_inputdata_vax, flu_activity,
                            n_iteration = 200, burnin = 100, thinning = 1,
                            group_by = ~vaccine),
    "multiple age groups"
  )
  expect_s3_class(fit, "seroreconstruct_multi")
})
