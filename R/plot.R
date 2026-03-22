# Plot functions for seroreconstruct

#' MCMC diagnostic plots
#'
#' Produces trace plots and posterior density plots for each model parameter.
#' Trace plots show the MCMC chain with the posterior mean (red dashed line).
#' Density plots show the marginal posterior with 95\% credible interval bounds
#' (blue dashed lines).
#'
#' @param fit A \code{seroreconstruct_fit}, \code{seroreconstruct_joint}, or
#'   \code{seroreconstruct_multi} object.
#' @param params Optional character vector of parameter names to plot. If
#'   \code{NULL} (default), all parameters are plotted. Use
#'   \code{table_parameters(fit)$Parameter} to see available names.
#' @return Invisible \code{NULL}. Called for its side effect of producing plots.
#' @examples
#' \donttest{
#' fit <- sero_reconstruct(inputdata, flu_activity,
#'                         n_iteration = 2000, burnin = 1000, thinning = 1)
#' plot_diagnostics(fit)
#' }
#' @export
plot_diagnostics <- function(fit, params = NULL) {
  if (inherits(fit, "seroreconstruct_multi")) {
    stop("plot_diagnostics() is not supported for seroreconstruct_multi objects. ",
         "Use it on individual group fits via fit[[\"group_name\"]].", call. = FALSE)
  }

  post <- fit$posterior_model_parameter
  param_names <- .get_param_names(fit)

  if (!is.null(params)) {
    idx <- match(params, param_names)
    if (any(is.na(idx))) {
      stop("Unknown parameter names: ",
           paste(params[is.na(idx)], collapse = ", "),
           ". Available: ", paste(param_names, collapse = ", "), call. = FALSE)
    }
    post <- post[, idx, drop = FALSE]
    param_names <- param_names[idx]
  }

  n_params <- ncol(post)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::par(mfrow = c(n_params, 2),
                mar = c(3, 3, 2, 1),
                mgp = c(2, 0.7, 0))

  for (i in seq_len(n_params)) {
    x <- post[, i]

    # Trace plot
    graphics::plot(x, type = "l",
                   main = param_names[i],
                   xlab = "Iteration (post-burnin)",
                   ylab = "Value",
                   col = "grey40")
    graphics::abline(h = mean(x), col = "red", lty = 2)

    # Density plot
    d <- stats::density(x)
    graphics::plot(d, main = "",
                   xlab = param_names[i],
                   ylab = "Density")
    ci <- stats::quantile(x, c(0.025, 0.975))
    graphics::abline(v = ci, col = "blue", lty = 2)
    graphics::abline(v = mean(x), col = "red", lty = 2)
  }

  invisible(NULL)
}

#' Plot antibody trajectory with model fit
#'
#' For a single individual, draws posterior-sampled antibody trajectories
#' overlaid on observed HAI titers. Red lines show trajectories where infection
#' occurred; blue lines show trajectories without infection. Matches the
#' visualization style of Figure 1B in Tsang et al. (2022).
#'
#' @param fit A \code{seroreconstruct_fit} or \code{seroreconstruct_joint}
#'   object.
#' @param id Row index (integer) or subject identifier to plot. If
#'   \code{subject_ids} was provided to \code{sero_reconstruct()}, or if
#'   \code{subjects} is supplied here, the value is matched against those IDs
#'   to find the row. Otherwise treated as a 1-based integer row index.
#' @param subjects Optional vector of subject identifiers aligned with fit rows.
#'   Use this when \code{subject_ids} was not provided at fitting time.
#'   Example: \code{subjects = inputdata$subject_id}.
#' @param n_samples Number of posterior samples to draw. Default 100.
#' @param main Optional plot title. If \code{NULL}, a default title with the
#'   individual index and posterior infection probability is generated.
#' @param col_infected Color for infected trajectories. Default semi-transparent
#'   red.
#' @param col_uninfected Color for uninfected trajectories. Default
#'   semi-transparent blue.
#' @param show_legend Logical; whether to draw a legend. Default \code{TRUE}.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#' @return Invisible \code{NULL}. Called for its side effect of producing a plot.
#' @examples
#' \donttest{
#' fit <- sero_reconstruct(inputdata, flu_activity,
#'                         n_iteration = 2000, burnin = 1000, thinning = 1)
#' plot_trajectory(fit, id = 1)
#' }
#' @export
plot_trajectory <- function(fit, id = 1, subjects = NULL, n_samples = 100,
                            main = NULL, col_infected = NULL,
                            col_uninfected = NULL, show_legend = TRUE, ...) {
  if (inherits(fit, "seroreconstruct_multi")) {
    stop("plot_trajectory() is not supported for seroreconstruct_multi objects. ",
         "Use it on individual group fits via fit[[\"group_name\"]].", call. = FALSE)
  }

  n_ind <- attr(fit, "n_individuals")

  # ---- Resolve id to row index ----
  # Numeric id → always treated as row index (1-based).
  # Non-numeric id → looked up in fit$subject_ids or the supplied subjects vector.
  if (is.numeric(id)) {
    if (length(id) != 1 || id < 1 || id > n_ind || id != round(id)) {
      stop("'id' must be an integer between 1 and ", n_ind, ".", call. = FALSE)
    }
    id_label <- as.character(id)
  } else {
    id_lookup <- if (!is.null(subjects)) subjects else fit$subject_ids
    if (is.null(id_lookup)) {
      stop("'id' is not numeric and no subject IDs are available. ",
           "Supply 'subjects' or pass 'subject_ids' to sero_reconstruct().",
           call. = FALSE)
    }
    row_idx <- match(id, id_lookup)
    if (is.na(row_idx)) {
      stop("Subject '", id, "' not found. ",
           "Available IDs (first 10): ",
           paste(utils::head(as.character(id_lookup), 10), collapse = ", "),
           call. = FALSE)
    }
    id_label <- as.character(id)
    id <- row_idx
  }

  n_post <- nrow(fit$posterior_inf_status)
  if (n_samples > n_post) n_samples <- n_post

  # ---- Extract data for this individual ----
  # fit$data columns (1-indexed): 4=start_time-14, 5=end_time-14,
  #   6=time1-14, 7=time2-14, 8=time3-14, 9=HAI1, 10=HAI2, 11=HAI3
  dat <- fit$data
  time1  <- dat[id, 6]
  end_t  <- dat[id, 5]

  # Blood draw times and observed titers
  obs_times <- dat[id, 6:8]   # adjusted scale
  obs_titers <- dat[id, 9:11]

  # Filter missing observations (encoded as -1)
  valid <- obs_times != -1 & obs_titers != -1
  obs_x <- obs_times[valid] - time1   # days relative to first blood draw
  obs_y <- obs_titers[valid]

  # ---- Posterior samples ----
  sample_idx <- sample(seq_len(n_post), n_samples)

  inf_status <- fit$posterior_inf_status[sample_idx, id]
  inf_time   <- fit$posterior_inf_time[sample_idx, id]
  wn_vec     <- fit$posterior_waning[sample_idx, id]
  bo_vec     <- fit$posterior_boosting[sample_idx, id]
  bl_vec     <- fit$posterior_baseline_titer[sample_idx, id]

  # P(infection) from full posterior (not just the sample)
  p_inf <- mean(fit$posterior_inf_status[, id])

  # ---- Colors ----
  if (is.null(col_infected))   col_infected   <- grDevices::rgb(1, 0, 0, 0.15)
  if (is.null(col_uninfected)) col_uninfected <- grDevices::rgb(0, 0, 1, 0.15)

  # ---- Plot setup ----
  x_max <- end_t - time1
  # Y range: use 95th percentile of trajectory peaks (not max)
  peak_titers <- ifelse(inf_status == 1, bl_vec + bo_vec, bl_vec)
  y_max_est <- max(c(obs_y, stats::quantile(peak_titers, 0.95)), na.rm = TRUE)
  y_max <- min(ceiling(y_max_est) + 1, 7)  # cap at 640
  y_max <- max(y_max, max(obs_y, na.rm = TRUE) + 1)  # ensure observed titers visible
  y_ticks <- 0:y_max
  y_labels <- c("<10", as.character(5 * 2^(seq_len(y_max))))

  if (is.null(main)) {
    main <- paste0("Individual ", id_label,
                   "    P(infection) = ", format(round(p_inf, 2), nsmall = 2))
  }

  # Save and restore only mar, not mfrow (so multi-panel layouts work)
  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar))
  graphics::par(mar = c(4, 4, 3, 1))

  graphics::plot(NA, xlim = c(0, x_max), ylim = c(0, y_max),
                 xlab = "Days since start of follow up",
                 ylab = "Antibody titer",
                 axes = FALSE, main = main, ...)
  graphics::axis(1)
  graphics::axis(2, at = y_ticks, labels = y_labels, las = 1)

  # ---- Draw posterior trajectories ----
  for (k in seq_len(n_samples)) {
    bl <- bl_vec[k]
    wn <- wn_vec[k]
    bo <- bo_vec[k]
    it <- inf_time[k]     # infection time (adjusted scale)
    is_inf <- inf_status[k]

    if (is_inf == 1) {
      # ---- Infected trajectory (red) ----
      # Phase 1: decay from time1 to 14 days before infection
      n1 <- it - time1 - 14

      if (n1 >= 0) {
        x1 <- 0:n1
        y1 <- bl * exp(-wn * x1 / 365)
        graphics::lines(x1, y1, col = col_infected)
        startpt <- y1[length(y1)]
      } else {
        # Boost started before first blood draw
        startpt <- bl * exp(-wn * max(n1, 0) / 365)
      }
      endpt <- startpt + bo

      # Phase 2: 14-day boost ramp
      x_bs <- it - time1 - 14
      x_be <- it - time1
      if (x_be >= 0) {
        graphics::lines(c(max(x_bs, 0), x_be),
                        c(if (x_bs < 0) startpt + bo * (abs(x_bs) / 14)
                          else startpt,
                          endpt),
                        col = col_infected)
      }

      # Phase 3: post-infection decay
      n3 <- end_t - it
      if (n3 > 0) {
        x3 <- 0:n3
        y3 <- endpt * exp(-wn * x3 / 365)
        graphics::lines(x3 + (it - time1), y3, col = col_infected)
      }

    } else {
      # ---- Uninfected trajectory (blue) ----
      n_days <- end_t - time1
      if (n_days > 0) {
        x_all <- 0:n_days
        y_all <- bl * exp(-wn * x_all / 365)
        graphics::lines(x_all, y_all, col = col_uninfected)
      }
    }
  }

  # ---- Observed titers ----
  graphics::points(obs_x, obs_y, pch = 16, col = "black", cex = 1.5)

  # ---- Legend ----
  if (show_legend) {
    graphics::legend("topright",
                     legend = c("Observed titer",
                                "Trajectories with infection",
                                "Trajectories without infection"),
                     col = c("black", "red", "blue"),
                     lty = c(NA, 1, 1), pch = c(16, NA, NA),
                     bty = "n", cex = 0.9)
  }

  invisible(NULL)
}
