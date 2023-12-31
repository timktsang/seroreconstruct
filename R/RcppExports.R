# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

prior_loglik <- function(para) {
    .Call(`_seroreconstruct_prior_loglik`, para)
}

sim_data <- function(data1, ILI, para, para2) {
    .Call(`_seroreconstruct_sim_data`, data1, ILI, para, para2)
}

loglik <- function(data11, data111, data21, ILI, para, para2, level1, level2, level3, season, blankmatrix) {
    .Call(`_seroreconstruct_loglik`, data11, data111, data21, ILI, para, para2, level1, level2, level3, season, blankmatrix)
}

all_update <- function(data11, data111, data21, ILI, para, para2, loglik1, loglik2, loglik3) {
    .Call(`_seroreconstruct_all_update`, data11, data111, data21, ILI, para, para2, loglik1, loglik2, loglik3)
}

add_remove_infection <- function(data11, data111, data21, ILI, para, para2, loglik1, loglik2, loglik3) {
    .Call(`_seroreconstruct_add_remove_infection`, data11, data111, data21, ILI, para, para2, loglik1, loglik2, loglik3)
}

mcmc <- function(input1, input2, input3, ILI, mcmc_n, int_para, int_para2, int_para3, paraseason, move, sigma, sigma3, burnin, thinning) {
    .Call(`_seroreconstruct_mcmc`, input1, input2, input3, ILI, mcmc_n, int_para, int_para2, int_para3, paraseason, move, sigma, sigma3, burnin, thinning)
}

