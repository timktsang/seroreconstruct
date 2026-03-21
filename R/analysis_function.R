#' Run the MCMC for the Bayesian model
#'
#' The main function to run the MCMC for the Bayesian model, to obtain individual dynamics, model parameters such as infection probability, boosting, waning, and measurement error.
#' @param inputdata The data for running MCMC, in dataframe format. It should be in the same format as the data in the package. It includes: 1) age_group (0: children, 1: adults, 2: older adults), 2) start_time: start of follow-up, 3) end_time: end of follow-up, 4) time1: date for first serum collection, 5) time2: date for second serum collection, 6) time3: date for third serum collection, 7) HAI_titer_1: HAI titer for first serum collection, 8) HAI_titer_2: HAI titer for second serum collection, 9) HAI_titer_3: HAI titer for third serum collection.
#' @param inputILI The data for influenza activity used in the inference. The row number should match with the date in the inputdata.
#' @param n_iteration The number of iterations of the MCMC.
#' @param burnin The iteration for burn-in for MCMC.
#' @param thinning The number of thinning in MCMC.
#' @param group_by Optional formula specifying grouping variables (e.g., \code{~age_group}).
#'   When provided, independent MCMCs are fit for each combination of the grouping
#'   variables. The formula uses interaction semantics: \code{~age + vac} means all
#'   age-by-vac combinations. Returns a \code{seroreconstruct_multi} object.
#' @details
#' \strong{Multi-season support (experimental):} If \code{inputdata} contains an
#' optional integer column named \code{season} (0-indexed, contiguous from 0 to
#' \code{n_seasons - 1}), the model fits season-specific infection risk and HAI
#' protection parameters. When no \code{season} column is present, all individuals
#' are assigned to a single season (\code{n_seasons = 1}) and behavior is identical
#' to previous versions.
#'
#' \strong{Current limitation:} \code{summary()} is not yet implemented for fits
#' with \code{n_seasons > 1}. Multi-season posterior samples are accessible directly
#' from the fit object (e.g., \code{fit$posterior_model_parameter}).
#'
#' @return A \code{seroreconstruct_fit} object (when \code{group_by} is \code{NULL})
#'   or a \code{seroreconstruct_multi} object (when \code{group_by} is provided).
#'   Use \code{summary()} to extract model estimates.
#' @examples
#' \donttest{
#' a1 <- sero_reconstruct(inputdata, flu_activity,
#'                         n_iteration = 2000, burnin = 1000, thinning = 1)
#' summary(a1)
#' }
#' @export
sero_reconstruct <- function(inputdata, inputILI, n_iteration = 2000, burnin = 1000,
                              thinning = 1, group_by = NULL) {
  inputdata <- .validate_inputs(inputdata, inputILI, n_iteration, burnin, thinning, group_by)

  if (!is.null(group_by)) {
    if (!inherits(group_by, "formula")) {
      stop("'group_by' must be a formula (e.g., ~age_group).", call. = FALSE)
    }
    group_vars <- all.vars(group_by)
    missing_vars <- setdiff(group_vars, colnames(inputdata))
    if (length(missing_vars) > 0) {
      stop("Variables not found in 'inputdata': ",
           paste(missing_vars, collapse = ", "), call. = FALSE)
    }

    groups <- interaction(inputdata[, group_vars, drop = FALSE], drop = TRUE)
    group_labels <- levels(groups)
    group_sizes <- vapply(group_labels, function(g) sum(groups == g), integer(1))

    for (g in group_labels) {
      if (group_sizes[g] < 10) {
        stop("Group '", g, "' has only ", group_sizes[g],
             " individuals (< 10 minimum).", call. = FALSE)
      }
      if (group_sizes[g] < 30) {
        warning("Group '", g, "' has only ", group_sizes[g],
                " individuals (< 30). Estimates may be unreliable.", call. = FALSE)
      }
    }

    # Warn if subgroups contain multiple age groups (age heterogeneity is
    # collapsed to a single-group model within each subgroup)
    has_mixed_ages <- any(vapply(group_labels, function(g) {
      length(unique(inputdata$age_group[groups == g])) > 1L
    }, logical(1)))
    if (has_mixed_ages) {
      warning("Some groups contain individuals from multiple age groups. ",
              "Each subgroup is fit as a single-group model (age heterogeneity ",
              "within subgroups is not modeled). To preserve age structure, ",
              "include 'age_group' in the group_by formula.", call. = FALSE)
    }

    fits <- list()
    for (g in group_labels) {
      message("Fitting group: ", g)
      subset_data <- inputdata[groups == g, , drop = FALSE]
      subset_data$age_group <- 0L
      fits[[g]] <- .fit_single(subset_data, inputILI, n_iteration, burnin,
                               thinning, n_groups = 1L)
    }

    return(.new_seroreconstruct_multi(fits, group_labels, group_sizes))
  }

  .fit_single(inputdata, inputILI, n_iteration, burnin, thinning, n_groups = 3L)
}

#################################################
#' Extract the model estimates from the fitted MCMC
#'
#' \code{output_model_estimate} is deprecated; use \code{summary()} instead.
#'
#' @param fitted_MCMC A \code{seroreconstruct_fit} object, or a list returned by
#'   an older version of \code{sero_reconstruct()}.
#' @param period A vector indicating the start and the end of a season to compute
#'   the infection probabilities. If empty, the start and end of the season are
#'   inferred from the data.
#' @return A data frame of model estimates (invisibly).
#' @examples
#' \dontrun{
#' a1 <- sero_reconstruct(inputdata, flu_activity,
#'                         n_iteration = 2000, burnin = 1000, thinning = 1)
#' fitted_result <- output_model_estimate(a1)  # deprecated, use summary(a1)
#' }
#' @export
output_model_estimate <- function(fitted_MCMC, period) {
  .Deprecated("summary", package = "seroreconstruct")

  if (!inherits(fitted_MCMC, "seroreconstruct_fit")) {
    if (is.list(fitted_MCMC) && "posterior_model_parameter" %in% names(fitted_MCMC)) {
      fitted_MCMC <- .new_seroreconstruct_fit(
        fitted_MCMC,
        n_groups = 3L,
        n_individuals = nrow(fitted_MCMC$data),
        n_posterior_samples = nrow(fitted_MCMC$posterior_model_parameter),
        runtime_secs = NULL
      )
    } else {
      stop("'fitted_MCMC' must be a seroreconstruct_fit object or a compatible list.",
           call. = FALSE)
    }
  }

  if (missing(period)) {
    s <- summary(fitted_MCMC)
  } else {
    s <- summary(fitted_MCMC, period = period)
  }

  output_print <- s$table
  output_print[, -1] <- round(output_print[, -1], 2)
  print(output_print)

  return(s$table)
}

#################################################
#' Simulation of the dataset of the Bayesian model
#'
#' The function to simulate the dataset, for validation or other purpose.
#' @param inputdata The data with the same format that for running MCMC, in dataframe format.
#' @param inputILI The data for influenza activity used in the inference. The row number should match with the date in the inputdata.
#' @param para1 The parameter vector for the model parameters, in the following format: 1) the parameter value of the random measurement error, 2) the parameter value of the 2-fold error, 3) the boosting in HAI titer after infection for children (in log2 unit), 4) the waning in HAI titer for children (in log2 unit), 5) the boosting in HAI titer after infection for adults (in log2 unit), 6) the waning in HAI titer for adults (in log2 unit), 7) the scale parameter for children, 8) the scale parameter for adults, 9) the scale parameter for older adults, 10) the log of risk ratio of 2-fold increase in baseline HAI titer
#' @param para2 The parameter vector for the parameter for the HAI titer distribution. The first 10 elements are the probability that the HAI titer is equal to 0-9 for children, and the elements of 11-20 is the probability that the HAI titer is equal to 0-9 for adults.
#' @return A simulated data based on the input parameter vectors, with the format equal to the input data.
#' @examples
#' simulated <- simulate_data(inputdata, flu_activity, para1, para2)
#' @export
simulate_data <- function(inputdata, inputILI, para1, para2) {
  if (!is.data.frame(inputdata)) {
    stop("'inputdata' must be a data frame.", call. = FALSE)
  }

  prepared <- .prepare_inputs(inputdata, inputILI)
  inputdata_mat <- prepared$inputdata
  inputILI_mat <- prepared$inputILI

  # Determine n_seasons from the season column (col 12 in R 1-indexed, 0-indexed values)
  n_seasons <- as.integer(max(inputdata_mat[, 12]) + 1L)
  hai_start_c <- 42L + 3L * n_seasons
  n_para <- 42L + 5L * n_seasons

  .validate_simulation_params(para1, para2, n_seasons)

  # Build full parameter vector and insert user params at active positions
  int_para <- c(0.005, rep(0.6, 17), rep(c(3.5, 0.5), 12),
                rep(c(0.4, 0.2, 0.2), n_seasons),
                rep(-0.1, 2L * n_seasons))

  move <- rep(0L, n_para)
  move[1] <- 1L
  move[2] <- 1L
  move[27:30] <- 1L
  move[43:(42L + 3L * n_seasons)] <- 1L
  for (s in seq_len(n_seasons)) {
    move[42L + 3L * n_seasons + 2L * (s - 1L) + 1L] <- 1L
  }
  int_para[which(move == 1)] <- para1

  int_para2 <- para2

  t <- sim_data(inputdata_mat, inputILI_mat, int_para, int_para2, hai_start_c)

  return(data.frame(t[[2]][, 2 + 1:9]))
}
