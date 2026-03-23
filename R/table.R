# Table functions for seroreconstruct

#' Get human-readable parameter names from a fit object
#'
#' @param fit A \code{seroreconstruct_fit} or derived object.
#' @return Character vector of parameter names.
#' @keywords internal
.get_param_names <- function(fit) {
  n_groups <- attr(fit, "n_groups")
  n_seasons <- attr(fit, "n_seasons")
  if (is.null(n_groups)) n_groups <- 3L
  if (is.null(n_seasons)) n_seasons <- 1L

  # Check for joint fit with custom group labels
  group_labels <- attr(fit, "group_labels")
  if (is.null(group_labels)) {
    if (n_groups == 3L) {
      group_labels <- c("children", "adults", "older_adults")
    } else if (n_groups == 1L) {
      group_labels <- "group"
    } else {
      group_labels <- paste0("group", seq_len(n_groups) - 1L)
    }
  }

  # Boosting/waning labels: use first two group labels (boost_wane_group maps

  # to at most 2 groups: group 0 = first label, group 1 = second label)
  bw_labels <- if (length(group_labels) >= 2L) {
    group_labels[1:2]
  } else {
    c(group_labels[1], group_labels[1])
  }

  shared <- c("random_error", "twofold_error",
              paste0("boosting_", bw_labels[1]),
              paste0("waning_", bw_labels[1]),
              paste0("boosting_", bw_labels[2]),
              paste0("waning_", bw_labels[2]))

  inf_names <- character(0)
  for (s in seq_len(n_seasons)) {
    season_suffix <- if (n_seasons > 1L) paste0("_s", s - 1L) else ""
    for (g in seq_len(n_groups)) {
      inf_names <- c(inf_names,
                     paste0("inf_prob_", group_labels[g], season_suffix))
    }
  }

  hai_names <- character(0)
  for (s in seq_len(n_seasons)) {
    season_suffix <- if (n_seasons > 1L) paste0("_s", s - 1L) else ""
    hai_names <- c(hai_names, paste0("hai_coef", season_suffix))
  }

  c(shared, inf_names, hai_names)
}

#' Compute effective sample size from MCMC chain
#'
#' Uses the initial positive sequence estimator (Geyer 1992).
#'
#' @param x Numeric vector of MCMC samples.
#' @return Effective sample size (numeric scalar).
#' @keywords internal
.effective_sample_size <- function(x) {
  n <- length(x)
  if (n < 4L) return(n)

  acf_vals <- acf(x, lag.max = n - 1L, plot = FALSE)$acf[, 1, 1]

  # Sum pairs of consecutive autocorrelations (Geyer's initial positive sequence)
  # acf_vals[1] = lag 0 (always 1.0); start pairing from lag 1 = index 2
  tau <- 1
  max_lag <- length(acf_vals)
  k <- 2L
  while (k + 1L <= max_lag) {
    pair_sum <- acf_vals[k] + acf_vals[k + 1L]
    if (pair_sum < 0) break
    tau <- tau + 2 * pair_sum
    k <- k + 2L
  }
  max(1, n / tau)
}

#' Summary table of model parameters with credible intervals
#'
#' Extracts posterior summaries (mean, median, credible intervals) for all
#' active model parameters.
#'
#' @param fit A \code{seroreconstruct_fit}, \code{seroreconstruct_joint}, or
#'   \code{seroreconstruct_multi} object.
#' @param probs Numeric vector of length 2 giving the lower and upper quantile
#'   probabilities for the credible interval. Default \code{c(0.025, 0.975)}
#'   for a 95\% interval.
#' @return A data frame with columns: Parameter, Mean, Median, Lower, Upper.
#' @examples
#' \donttest{
#' fit <- sero_reconstruct(inputdata, flu_activity,
#'                         n_iteration = 2000, burnin = 1000, thinning = 1)
#' table_parameters(fit)
#' }
#' @export
table_parameters <- function(fit, probs = c(0.025, 0.975)) {
  if (inherits(fit, "seroreconstruct_multi")) {
    # For multi-group fits, combine parameter tables with a Group column
    group_labels <- attr(fit, "group_labels")
    tables <- vector("list", length(group_labels))
    for (j in seq_along(group_labels)) {
      tbl <- table_parameters(unclass(fit)[[group_labels[j]]], probs = probs)
      tbl <- cbind(Group = group_labels[j], tbl)
      tables[[j]] <- tbl
    }
    result <- do.call(rbind, tables)
    rownames(result) <- NULL
    return(result)
  }

  post <- fit$posterior_model_parameter
  param_names <- .get_param_names(fit)

  result <- data.frame(
    Parameter = param_names,
    Mean = colMeans(post),
    Median = apply(post, 2, stats::median),
    Lower = apply(post, 2, stats::quantile, probs[1]),
    Upper = apply(post, 2, stats::quantile, probs[2]),
    stringsAsFactors = FALSE
  )
  rownames(result) <- NULL
  result
}

#' Per-individual infection estimates
#'
#' Summarizes posterior infection status, timing, and baseline titer for each
#' individual in the dataset.
#'
#' @param fit A \code{seroreconstruct_fit} or \code{seroreconstruct_joint} object.
#' @return A data frame with one row per individual and columns:
#'   \code{Individual} (row index), \code{Infection_prob} (posterior mean
#'   probability of infection), \code{Infection_time_mean} (mean infection time
#'   among infected samples), \code{Baseline_titer_mean} (mean imputed baseline
#'   HAI titer).
#' @examples
#' \donttest{
#' fit <- sero_reconstruct(inputdata, flu_activity,
#'                         n_iteration = 2000, burnin = 1000, thinning = 1)
#' head(table_infections(fit))
#' }
#' @export
table_infections <- function(fit) {
  if (inherits(fit, "seroreconstruct_multi")) {
    stop("table_infections() is not supported for seroreconstruct_multi objects. ",
         "Use it on individual group fits via fit[[\"group_name\"]].", call. = FALSE)
  }

  inf_status <- fit$posterior_inf_status
  inf_time <- fit$posterior_inf_time
  baseline <- fit$posterior_baseline_titer

  n_ind <- ncol(inf_status)
  inf_prob <- colMeans(inf_status)

  # Mean infection time (only among infected samples)
  inf_time_weighted <- colSums(inf_time * inf_status)
  inf_count <- colSums(inf_status)
  inf_time_mean <- ifelse(inf_count > 0,
                          inf_time_weighted / inf_count,
                          NA_real_)

  data.frame(
    Individual = seq_len(n_ind),
    Infection_prob = round(inf_prob, 4),
    Infection_time_mean = round(inf_time_mean, 1),
    Baseline_titer_mean = round(colMeans(baseline), 2),
    stringsAsFactors = FALSE
  )
}
