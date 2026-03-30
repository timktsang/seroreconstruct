# Shared fit objects — computed once, reused across all test files.
# Avoids redundant MCMC runs to keep total check time under CRAN's 10-min limit.

data("inputdata", package = "seroreconstruct")
data("flu_activity", package = "seroreconstruct")

# Standard 3-group single-season fit
test_fit <- sero_reconstruct(inputdata, flu_activity,
                             n_iteration = 200, burnin = 100, thinning = 1)

# Multi-group fit (seroreconstruct_multi)
test_fit_multi <- sero_reconstruct(inputdata, flu_activity,
                                   n_iteration = 200, burnin = 100, thinning = 1,
                                   group_by = ~age_group)

# Joint model fit (seroreconstruct_joint)
test_inputdata_vax <- inputdata
test_inputdata_vax$vaccine <- rep(c(0L, 1L), length.out = nrow(inputdata))
test_fit_joint <- sero_reconstruct(test_inputdata_vax, flu_activity,
                                   n_iteration = 200, burnin = 100, thinning = 1,
                                   group_by = ~vaccine,
                                   shared = c("error", "boosting_waning"))

# Multi-season fit
test_d1 <- inputdata[1:100, ]
test_d1$season <- 0L
test_d2 <- inputdata[101:200, ]
test_d2$season <- 1L
test_inputdata_ms <- rbind(test_d1, test_d2)
test_fit_ms <- sero_reconstruct(test_inputdata_ms, flu_activity,
                                n_iteration = 200, burnin = 100, thinning = 1)

# Fit with character subject_ids (for plot_trajectory tests)
test_char_ids <- paste0("S", seq_len(nrow(inputdata)))
test_fit_ids <- sero_reconstruct(inputdata, flu_activity,
                                 n_iteration = 200, burnin = 100, thinning = 1,
                                 subject_ids = test_char_ids)
