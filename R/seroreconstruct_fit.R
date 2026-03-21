# S3 classes for seroreconstruct

# ---- seroreconstruct_fit (single fit) ----

#' Constructor for seroreconstruct_fit objects
#' @keywords internal
.new_seroreconstruct_fit <- function(output, n_groups, n_individuals,
                                      n_posterior_samples, n_seasons = 1L,
                                      runtime_secs) {
  attr(output, "n_groups") <- n_groups
  attr(output, "n_individuals") <- n_individuals
  attr(output, "n_posterior_samples") <- n_posterior_samples
  attr(output, "n_seasons") <- n_seasons
  attr(output, "runtime_secs") <- runtime_secs
  class(output) <- "seroreconstruct_fit"
  output
}

#' Print method for seroreconstruct_fit
#'
#' @param x A \code{seroreconstruct_fit} object.
#' @param ... Additional arguments (ignored).
#' @export
print.seroreconstruct_fit <- function(x, ...) {
  cat("seroreconstruct fit\n")
  cat("  Individuals:", attr(x, "n_individuals"), "\n")
  cat("  Age groups:", attr(x, "n_groups"), "\n")
  n_seasons <- attr(x, "n_seasons")
  if (!is.null(n_seasons) && n_seasons > 1L) {
    cat("  Seasons:", n_seasons, "\n")
  }
  cat("  Posterior samples:", attr(x, "n_posterior_samples"), "\n")
  runtime <- attr(x, "runtime_secs")
  if (!is.null(runtime)) {
    cat("  Runtime:", round(runtime), "seconds\n")
  }
  cat("\nUse summary() to extract model estimates.\n")
  invisible(x)
}

#' Summary method for seroreconstruct_fit
#'
#' Computes estimates of infection probabilities, boosting, waning, and
#' measurement error from a fitted MCMC object.
#'
#' @param object A \code{seroreconstruct_fit} object.
#' @param period Optional numeric vector of length 2 specifying the start and end
#'   of a season to compute infection probabilities. If omitted, the full
#'   follow-up period is used.
#' @param ... Additional arguments (ignored).
#' @return A \code{summary.seroreconstruct_fit} object with element \code{$table}.
#' @export
summary.seroreconstruct_fit <- function(object, period, ...) {
  n_seasons <- attr(object, "n_seasons")
  if (is.null(n_seasons)) n_seasons <- 1L
  if (n_seasons > 1L) {
    stop("summary() for multi-season fits (n_seasons = ", n_seasons,
         ") is not yet implemented. Use the raw posterior samples directly.",
         call. = FALSE)
  }

  z1 <- .mcmc_summary(object$posterior_model_parameter)
  z1[1, ] <- z1[1, ] * 10 * 100
  z1[2, ] <- 0.25 / exp(z1[2, ]) * 100

  n_groups <- attr(object, "n_groups")

  mean_ci <- function(v) c(mean(v), quantile(v, c(0.025, 0.975)))

  if (n_groups == 3L) {
    z1[3:6, ] <- 2^(z1[3:6, ])

    ILI <- object$ILI_data
    mcmc1 <- object$posterior_model_parameter
    mcmc2 <- cbind(object$posterior_baseline_HAI_titer,
                   object$posterior_baseline_HAI_titer[, 11:20])

    xvec <- 0:100
    d1list <- vector("list", 3)
    d1 <- matrix(NA, nrow(mcmc1), 101)

    if (missing(period)) {
      indrow <- min(object$data[, "start_time"]):max(object$data[, "end_time"])
    } else {
      indrow <- period[1]:period[2]
    }

    for (u in 1:3) {
      for (i in seq_len(nrow(mcmc1))) {
        d1[i, ] <- 1 - exp(-(mcmc1[i, 6 + u] * sum(ILI[indrow, 1])) *
                              exp(xvec / 10 * mcmc1[i, 10]))
      }
      d1list[[u]] <- d1
    }

    vec <- vector("list", 3)
    for (uu in 1:3) {
      vec[[uu]] <- rowSums(d1list[[uu]][, c(6 + 0:9 * 10)] *
                             mcmc2[, 1:10 + 10 * (uu - 1)])
    }

    pm <- matrix(NA, 4, 9)
    for (v in 0:2) {
      pm[1, 1:3 + 3 * v] <- mean_ci(vec[[v + 1]])
      pm[2, 1:3 + 3 * v] <- mean_ci(vec[[v + 1]] / vec[[2]])
      pm[3, 1:3 + 3 * v] <- mean_ci(d1list[[v + 1]][, 1])
      pm[4, 1:3 + 3 * v] <- mean_ci(d1list[[v + 1]][, 1] / d1list[[2]][, 1])
    }

    output <- data.frame(matrix(NA, 16, 4))
    output[1:6, ] <- z1[1:6, c(4, 1:3)]
    output[7:16, 2:4] <- rbind(pm[, 1:3], pm[c(1, 3), 4:6],
                                pm[, 7:9])[c(1, 5, 7, 3, 6, 9, 2, 4, 8, 10), ]

    names(output) <- c("Variable", "Point estimate", "Lower bound", "Upper bound")
    output[, 1] <- c(
      "Random error (%)",
      "Two-fold error (%)",
      "Fold-increase after infection for children (Boosting)",
      "Fold-decrease after 1 year for children (Waning)",
      "Fold-increase after 1 year for adults (Boosting)",
      "Fold-decrease after 1 year for adults (Waning)",
      "Infection probability for children",
      "Infection probability for adults",
      "Infection probability for older adults",
      "Infection probability for children with pre-epidemic HAI titer < 10",
      "Infection probability for adults with pre-epidemic HAI titer < 10",
      "Infection probability for older adults with pre-epidemic HAI titer < 10",
      "Relative risk for children (Ref: Adults)",
      "Relative risk for older adults (Ref: Adults)",
      "Relative risk for children with pre-epidemic HAI titer < 10 (Ref: Adults with pre-epidemic HAI titer < 10)",
      "Relative risk for older adults with pre-epidemic HAI titer < 10 (Ref: Adults with pre-epidemic HAI titer < 10)"
    )
  } else {
    # Single-group fit (from group_by)
    z1[3:4, ] <- 2^(z1[3:4, ])

    ILI <- object$ILI_data
    mcmc1 <- object$posterior_model_parameter
    mcmc2 <- object$posterior_baseline_HAI_titer[, 1:10]

    xvec <- 0:100
    d1 <- matrix(NA, nrow(mcmc1), 101)

    if (missing(period)) {
      indrow <- min(object$data[, "start_time"]):max(object$data[, "end_time"])
    } else {
      indrow <- period[1]:period[2]
    }

    for (i in seq_len(nrow(mcmc1))) {
      d1[i, ] <- 1 - exp(-(mcmc1[i, 7] * sum(ILI[indrow, 1])) *
                            exp(xvec / 10 * mcmc1[i, 10]))
    }

    vec <- rowSums(d1[, c(6 + 0:9 * 10)] * mcmc2)

    output <- data.frame(matrix(NA, 6, 4))
    output[1:4, ] <- z1[1:4, c(4, 1:3)]
    output[5, 2:4] <- mean_ci(vec)
    output[6, 2:4] <- mean_ci(d1[, 1])

    names(output) <- c("Variable", "Point estimate", "Lower bound", "Upper bound")
    output[, 1] <- c(
      "Random error (%)",
      "Two-fold error (%)",
      "Fold-increase after infection (Boosting)",
      "Fold-decrease after 1 year (Waning)",
      "Infection probability",
      "Infection probability (HAI titer < 10)"
    )
  }

  result <- list(table = output, n_groups = n_groups)
  class(result) <- "summary.seroreconstruct_fit"
  result
}

#' Print method for summary.seroreconstruct_fit
#'
#' @param x A \code{summary.seroreconstruct_fit} object.
#' @param digits Number of decimal places for rounding. Default 2.
#' @param ... Additional arguments (ignored).
#' @method print summary.seroreconstruct_fit
#' @export
print.summary.seroreconstruct_fit <- function(x, digits = 2, ...) {
  tbl <- x$table
  numeric_cols <- which(vapply(tbl, is.numeric, logical(1)))
  tbl[, numeric_cols] <- round(tbl[, numeric_cols], digits)
  print(tbl, row.names = FALSE)
  invisible(x)
}

# ---- seroreconstruct_multi (multi-group fit) ----

#' Constructor for seroreconstruct_multi objects
#' @keywords internal
.new_seroreconstruct_multi <- function(fits, group_labels, group_sizes) {
  attr(fits, "group_labels") <- group_labels
  attr(fits, "group_sizes") <- group_sizes
  class(fits) <- "seroreconstruct_multi"
  fits
}

#' Print method for seroreconstruct_multi
#'
#' @param x A \code{seroreconstruct_multi} object.
#' @param ... Additional arguments (ignored).
#' @export
print.seroreconstruct_multi <- function(x, ...) {
  group_labels <- attr(x, "group_labels")
  group_sizes <- attr(x, "group_sizes")
  cat("seroreconstruct multi-group fit\n")
  cat("  Groups:", length(group_labels), "\n")
  for (i in seq_along(group_labels)) {
    marker <- ""
    if (group_sizes[i] < 30) marker <- " (warning: N < 30)"
    cat("    ", group_labels[i], ": N =", group_sizes[i], marker, "\n")
  }
  cat("\nUse summary() to extract model estimates.\n")
  cat("Use x[[\"group_name\"]] or x[[index]] to access individual fits.\n")
  invisible(x)
}

#' Summary method for seroreconstruct_multi
#'
#' Computes estimates for each group and combines into a single table.
#'
#' @param object A \code{seroreconstruct_multi} object.
#' @param period Optional numeric vector of length 2 specifying the start and end
#'   of a season to compute infection probabilities.
#' @param ... Additional arguments (ignored).
#' @return A \code{summary.seroreconstruct_multi} object with element \code{$table}.
#' @export
summary.seroreconstruct_multi <- function(object, period, ...) {
  group_labels <- attr(object, "group_labels")
  tables <- vector("list", length(group_labels))
  for (j in seq_along(group_labels)) {
    g <- group_labels[j]
    fit_g <- unclass(object)[[g]]
    if (missing(period)) {
      s <- summary(fit_g)
    } else {
      s <- summary(fit_g, period = period)
    }
    tbl <- s$table
    tbl <- cbind(Group = g, tbl)
    tables[[j]] <- tbl
  }
  combined <- do.call(rbind, tables)
  rownames(combined) <- NULL

  result <- list(table = combined, group_labels = group_labels)
  class(result) <- "summary.seroreconstruct_multi"
  result
}

#' Print method for summary.seroreconstruct_multi
#'
#' @param x A \code{summary.seroreconstruct_multi} object.
#' @param digits Number of decimal places for rounding. Default 2.
#' @param ... Additional arguments (ignored).
#' @method print summary.seroreconstruct_multi
#' @export
print.summary.seroreconstruct_multi <- function(x, digits = 2, ...) {
  tbl <- x$table
  numeric_cols <- which(vapply(tbl, is.numeric, logical(1)))
  tbl[, numeric_cols] <- round(tbl[, numeric_cols], digits)
  print(tbl, row.names = FALSE)
  invisible(x)
}

#' Subset a seroreconstruct_multi object
#'
#' Access an individual group fit by name or index.
#'
#' @param x A \code{seroreconstruct_multi} object.
#' @param i Group name (character) or index (integer).
#' @param ... Additional arguments (ignored).
#' @return A \code{seroreconstruct_fit} object for the requested group.
#' @method [[ seroreconstruct_multi
#' @export
`[[.seroreconstruct_multi` <- function(x, i, ...) {
  group_labels <- attr(x, "group_labels")
  if (is.character(i)) {
    if (!(i %in% group_labels)) {
      stop("Group '", i, "' not found. Available groups: ",
           paste(group_labels, collapse = ", "), call. = FALSE)
    }
  }
  unclass(x)[[i]]
}
