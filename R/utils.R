# Internal utilities for seroreconstruct
# These functions are not exported.

#' Prepare input data for MCMC or simulation
#'
#' Reorders columns, adds padding columns, adjusts times for boosting delay,
#' and converts to matrix format expected by the C++ backend.
#'
#' @param inputdata Data frame with 9 required columns.
#' @param inputILI Data frame or matrix of influenza activity.
#' @return A list with prepared \code{inputdata} and \code{inputILI} matrices.
#' @keywords internal
.prepare_inputs <- function(inputdata, inputILI) {
  if ("HAI_titer_3" %in% colnames(inputdata) && !"HAI_titer3" %in% colnames(inputdata)) {
    colnames(inputdata)[colnames(inputdata) == "HAI_titer_3"] <- "HAI_titer3"
  }
  inputdata <- inputdata[, c("age_group", "start_time", "end_time", "time1", "time2",
                              "time3", "HAI_titer_1", "HAI_titer_2", "HAI_titer3")]
  inputdata <- cbind(0, 0, inputdata)
  inputdata[, 4:8] <- inputdata[, 4:8] - 14
  inputdata[is.na(inputdata)] <- -1
  inputdata <- cbind(inputdata, 2, 0, 2)
  inputILI[inputILI < 0] <- 1e-11
  inputdata <- as.matrix(inputdata)
  inputILI <- as.matrix(inputILI)
  list(inputdata = inputdata, inputILI = inputILI)
}

#' Compute MCMC summary statistics
#'
#' @param mcmc_matrix Matrix of posterior samples (rows = iterations, cols = parameters).
#' @return A matrix with \code{ncol(mcmc_matrix)} rows and 4 columns:
#'   mean, 2.5\% quantile, 97.5\% quantile, and acceptance rate.
#' @keywords internal
.mcmc_summary <- function(mcmc_matrix) {
  y <- matrix(NA, ncol(mcmc_matrix), 4)
  for (i in seq_len(ncol(mcmc_matrix))) {
    y[i, 1:3] <- quantile(mcmc_matrix[, i], c(0.5, 0.025, 0.975), na.rm = TRUE)
    y[i, 4] <- sum(diff(mcmc_matrix[, i]) != 0) / nrow(mcmc_matrix)
    y[i, 1] <- mean(mcmc_matrix[, i], na.rm = TRUE)
  }
  y
}

#' Plot MCMC trace plots
#'
#' @param mcmc_matrix Matrix of posterior samples.
#' @param nrow Number of rows in the plot layout.
#' @param ncol Number of columns in the plot layout.
#' @keywords internal
.plot_traces <- function(mcmc_matrix, nrow, ncol) {
  layout(matrix(seq_len(nrow * ncol), nrow = nrow, byrow = TRUE))
  par(mar = c(2, 4, 1, 1))
  for (i in seq_len(ncol(mcmc_matrix))) {
    plot(mcmc_matrix[, i], type = "l")
  }
}

#' Validate inputs for sero_reconstruct
#'
#' @param inputdata Data frame of individual-level data.
#' @param inputILI Data frame or matrix of influenza activity.
#' @param n_iteration Number of MCMC iterations.
#' @param burnin Burn-in iterations.
#' @param thinning Thinning interval.
#' @param group_by Optional formula; when non-NULL, age_group value checks are skipped.
#' @return The (possibly column-renamed) \code{inputdata}.
#' @keywords internal
.validate_inputs <- function(inputdata, inputILI, n_iteration, burnin, thinning,
                             group_by = NULL) {
  if (!is.data.frame(inputdata)) {
    stop("'inputdata' must be a data frame.", call. = FALSE)
  }
  if ("HAI_titer_3" %in% colnames(inputdata) && !"HAI_titer3" %in% colnames(inputdata)) {
    colnames(inputdata)[colnames(inputdata) == "HAI_titer_3"] <- "HAI_titer3"
  }
  required_cols <- c("age_group", "start_time", "end_time", "time1", "time2",
                     "time3", "HAI_titer_1", "HAI_titer_2", "HAI_titer3")
  missing_cols <- setdiff(required_cols, colnames(inputdata))
  if (length(missing_cols) > 0) {
    stop("'inputdata' is missing required columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!is.data.frame(inputILI) && !is.matrix(inputILI)) {
    stop("'inputILI' must be a data frame or matrix.", call. = FALSE)
  }
  if (NROW(inputILI) == 0) {
    stop("'inputILI' must have at least one row.", call. = FALSE)
  }
  max_time <- max(inputdata$end_time, na.rm = TRUE)
  if (NROW(inputILI) < max_time) {
    stop("'inputILI' has ", NROW(inputILI), " rows but 'inputdata' requires at least ",
         max_time, " rows (max end_time). Provide influenza activity data covering ",
         "the full follow-up period.", call. = FALSE)
  }
  if (!is.numeric(n_iteration) || length(n_iteration) != 1 || n_iteration < 1) {
    stop("'n_iteration' must be a positive integer.", call. = FALSE)
  }
  if (n_iteration != round(n_iteration)) {
    stop("'n_iteration' must be a whole number, got ", n_iteration, ".", call. = FALSE)
  }
  if (!is.numeric(burnin) || length(burnin) != 1 || burnin < 0) {
    stop("'burnin' must be a non-negative integer.", call. = FALSE)
  }
  if (burnin != round(burnin)) {
    stop("'burnin' must be a whole number, got ", burnin, ".", call. = FALSE)
  }
  if (!is.numeric(thinning) || length(thinning) != 1 || thinning < 1) {
    stop("'thinning' must be a positive integer.", call. = FALSE)
  }
  if (thinning != round(thinning)) {
    stop("'thinning' must be a whole number, got ", thinning, ".", call. = FALSE)
  }
  if (n_iteration <= burnin) {
    stop("'n_iteration' must be greater than 'burnin'.", call. = FALSE)
  }
  if (is.null(group_by)) {
    age_vals <- unique(inputdata$age_group[!is.na(inputdata$age_group)])
    invalid_ages <- setdiff(age_vals, c(0, 1, 2))
    if (length(invalid_ages) > 0) {
      stop("'inputdata$age_group' must contain only values 0, 1, 2. Found: ",
           paste(invalid_ages, collapse = ", "), call. = FALSE)
    }
    for (ag in 0:2) {
      n_ag <- sum(inputdata$age_group == ag, na.rm = TRUE)
      if (n_ag > 0 && n_ag < 30) {
        warning("Age group ", ag, " has only ", n_ag, " individuals (< 30). ",
                "Estimates may be unreliable.", call. = FALSE)
      }
    }
  }
  inputdata
}

#' Validate simulation parameters
#'
#' @param para1 Numeric vector of length 10.
#' @param para2 Numeric vector of length 20.
#' @keywords internal
.validate_simulation_params <- function(para1, para2) {
  if (!is.numeric(para1) || length(para1) != 10) {
    stop("'para1' must be a numeric vector of length 10.", call. = FALSE)
  }
  if (!is.numeric(para2) || length(para2) != 20) {
    stop("'para2' must be a numeric vector of length 20.", call. = FALSE)
  }
  if (any(!is.finite(para1))) {
    stop("'para1' must contain only finite values.", call. = FALSE)
  }
  if (any(!is.finite(para2))) {
    stop("'para2' must contain only finite values.", call. = FALSE)
  }
  invisible(NULL)
}

#' Run MCMC for a single fit
#'
#' Core fitting function used by \code{\link{sero_reconstruct}} for both
#' standard 3-age-group fits and single-group sub-fits.
#'
#' @param inputdata Data frame of individual-level data.
#' @param inputILI Data frame or matrix of influenza activity.
#' @param n_iteration Number of MCMC iterations.
#' @param burnin Burn-in iterations.
#' @param thinning Thinning interval.
#' @param n_groups Number of effective age groups (3 for standard, 1 for single-group).
#' @return A \code{seroreconstruct_fit} object.
#' @keywords internal
.fit_single <- function(inputdata, inputILI, n_iteration, burnin, thinning, n_groups = 3L) {
  t_start <- Sys.time()
  n_individuals <- nrow(inputdata)
  keep_iteration <- burnin + 1:((n_iteration - burnin) / thinning) * thinning

  prepared <- .prepare_inputs(inputdata, inputILI)
  inputdata <- prepared$inputdata
  inputILI <- prepared$inputILI

  int_para <- c(0.005, rep(0.6, 17), rep(c(3.5, 0.5), 12),
                rep(c(0.4, 0.2, 0.2), 7), rep(-0.1, 12))

  int_para2 <- c(0.82, 0.08, 0.065, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
                 0.87, 0.04, 0.04, 0.015, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005,
                 0.82, 0.08, 0.065, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
                 0.87, 0.04, 0.04, 0.015, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005,
                 0.82, 0.08, 0.065, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
                 0.87, 0.04, 0.04, 0.015, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005,
                 0.82, 0.08, 0.065, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
                 0.87, 0.04, 0.04, 0.015, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005,
                 0.82, 0.08, 0.065, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
                 0.87, 0.04, 0.04, 0.015, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005,
                 0.82, 0.08, 0.065, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
                 0.87, 0.04, 0.04, 0.015, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005)

  t <- sim_data(inputdata, inputILI, int_para, int_para2)

  input1 <- t[[1]]
  input3 <- t[[3]]

  sigma <- abs(int_para) / 10
  move <- rep(1, length(int_para))
  int_para3 <- rep(1, 12)
  sigma3 <- int_para3 / 10
  move[c(4:18, 35:42, 64:65, 67, 69, 71, 73, 75)] <- 0
  move[c(3, 19:26, 31:51, 55:67, 69:75)] <- 0
  int_para[64:65] <- 0
  paraseason <- c(rep(0, 42), rep(1, 6), rep(2:6, each = 3), rep(1:6, each = 2))

  tt <- mcmc(input1, inputdata, input3, inputILI, n_iteration,
             int_para, int_para2, int_para3, paraseason,
             move, sigma, sigma3, burnin, thinning)

  t_end <- Sys.time()
  runtime <- as.numeric(difftime(t_end, t_start, units = "secs"))

  output <- list(
    posterior_model_parameter = tt[[1]][keep_iteration, which(move == 1)],
    posterior_baseline_HAI_titer = tt[[2]][keep_iteration, 41:60],
    posterior_inf_status = tt[[13]],
    posterior_inf_time = tt[[14]],
    posterior_waning = tt[[15]],
    posterior_boosting = tt[[16]],
    posterior_baseline_titer = tt[[17]],
    data = inputdata,
    ILI_data = inputILI
  )

  message("MCMC complete in ", round(runtime), " seconds. Use summary() to view estimates.")

  .new_seroreconstruct_fit(output,
                            n_groups = n_groups,
                            n_individuals = n_individuals,
                            n_posterior_samples = length(keep_iteration),
                            runtime_secs = runtime)
}
