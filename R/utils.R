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
  # Capture season column before subsetting to required columns
  if ("season" %in% colnames(inputdata)) {
    season_col <- inputdata$season
  } else {
    season_col <- rep(0L, nrow(inputdata))
  }
  inputdata <- inputdata[, c("age_group", "start_time", "end_time", "time1", "time2",
                              "time3", "HAI_titer_1", "HAI_titer_2", "HAI_titer3")]
  inputdata <- cbind(0, 0, inputdata)
  inputdata[, 4:8] <- inputdata[, 4:8] - 14
  inputdata[is.na(inputdata)] <- -1
  inputdata <- cbind(inputdata, season_col, 0, 2)
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
  # Validate season column if present
  if ("season" %in% colnames(inputdata)) {
    season_vals <- inputdata$season
    if (!is.numeric(season_vals) || any(season_vals != as.integer(season_vals))) {
      stop("'inputdata$season' must contain integer values.", call. = FALSE)
    }
    if (min(season_vals) != 0L) {
      stop("'inputdata$season' must be 0-indexed (minimum value must be 0).",
           call. = FALSE)
    }
    n_seasons <- max(season_vals) + 1L
    expected <- seq(0L, n_seasons - 1L)
    if (!all(expected %in% season_vals)) {
      stop("'inputdata$season' must be contiguous (0 to ", n_seasons - 1L,
           "). Missing seasons: ",
           paste(setdiff(expected, unique(season_vals)), collapse = ", "),
           call. = FALSE)
    }
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
#' @param para1 Numeric vector of active model parameters (length \code{6 + 4 * n_seasons}).
#' @param para2 Numeric vector of baseline HAI titer distribution (length \code{20 * n_seasons}).
#' @param n_seasons Number of seasons. Default 1.
#' @keywords internal
.validate_simulation_params <- function(para1, para2, n_seasons = 1L) {
  expected_para1_len <- 6L + 4L * n_seasons
  if (!is.numeric(para1) || length(para1) != expected_para1_len) {
    stop("'para1' must be a numeric vector of length ", expected_para1_len,
         " (6 + 4 * n_seasons).", call. = FALSE)
  }
  expected_para2_len <- 20L * n_seasons
  if (!is.numeric(para2) || length(para2) != expected_para2_len) {
    stop("'para2' must be a numeric vector of length ", expected_para2_len,
         " (20 * n_seasons).", call. = FALSE)
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

  # Determine n_seasons from the season column (col 12 in R 1-indexed, 0-indexed values)
  n_seasons <- as.integer(max(inputdata[, 12]) + 1L)
  hai_start_c <- 42L + 3L * n_seasons  # C++ 0-indexed
  n_para <- 42L + 5L * n_seasons

  # Build parameter vectors dynamically
  # Block 1: random error (1) + two-fold error (17) + boosting/waning (24)
  # Block 2: infection risk (3 per season)
  # Block 3: HAI protection (2 per season)
  int_para <- c(0.005, rep(0.6, 17), rep(c(3.5, 0.5), 12),
                rep(c(0.4, 0.2, 0.2), n_seasons),
                rep(-0.1, 2L * n_seasons))

  # Baseline HAI titer distribution: 10 levels x 2 age groups x n_seasons
  titer_block <- c(0.82, 0.08, 0.065, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
                   0.87, 0.04, 0.04, 0.015, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005)
  int_para2 <- rep(titer_block, n_seasons)

  # Dirichlet hyperparameters: 2 per season
  int_para3 <- rep(1, 2L * n_seasons)

  t <- sim_data(inputdata, inputILI, int_para, int_para2, hai_start_c)

  input1 <- t[[1]]
  input3 <- t[[3]]

  sigma <- abs(int_para) / 10
  sigma3 <- int_para3 / 10

  # Move vector: which parameters are active in MCMC
  move <- rep(0L, n_para)
  move[1] <- 1L  # random error
  move[2] <- 1L  # two-fold error (first)
  move[27:30] <- 1L  # boosting/waning (4 params for boost_group=2)
  # All infection risk params
  move[43:(42L + 3L * n_seasons)] <- 1L
  # First HAI param per season (the one that moves; second is fixed)
  for (s in seq_len(n_seasons)) {
    move[42L + 3L * n_seasons + 2L * (s - 1L) + 1L] <- 1L
  }

  # Season assignment for each parameter (used for selective likelihood updates)
  paraseason <- c(rep(0L, 42L),
                  rep(seq_len(n_seasons), each = 3L),
                  rep(seq_len(n_seasons), each = 2L))

  tt <- mcmc(input1, inputdata, input3, inputILI, n_iteration,
             int_para, int_para2, int_para3, paraseason,
             move, sigma, sigma3, burnin, thinning, n_seasons)

  t_end <- Sys.time()
  runtime <- as.numeric(difftime(t_end, t_start, units = "secs"))

  # Extract posterior baseline HAI titer for the first season
  # para2 columns: 10 levels x 2 age groups x n_seasons
  # For backward compat, extract columns for the first season + duplicate for 3-group display
  hai_cols_start <- 20L * (n_seasons - 1L) + 1L
  hai_cols_end <- 20L * n_seasons
  posterior_hai <- tt[[2]][keep_iteration, hai_cols_start:hai_cols_end, drop = FALSE]

  output <- list(
    posterior_model_parameter = tt[[1]][keep_iteration, which(move == 1), drop = FALSE],
    posterior_baseline_HAI_titer = posterior_hai,
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
                            n_seasons = n_seasons,
                            runtime_secs = runtime)
}
