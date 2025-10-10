#' Plots histograms of the estimates
#'
#' @description This function plots histograms of the R estimates, highlighting
#'   the true value, and histograms of the standard error estimates
#' @param df_R_hat a data frame with R estimates, after replacing the divergent
#'   runs by NAs, containing columns
#' @param R_true the true value of the effective reproduction number used to
#'   generate the trajectory
#' @param model_colors a named vector specifying the model colors. Must have
#'   names "Poiss", "Q-Poiss", "NegBin-L" and "NegBin-Q"
#' @param limits_x a list of the limits for the x-axis, must have 2 elements
#'   called \code{R_hat} and \code{se_hat}. Both elements must contain a vector
#'   with 2 elements - the axis range.
#' @return a list containing the two histograms named as \code{R_hat} and
#'   \code{se_hat}
plot_hists <- function(df_R_hat, R_true, model_colors, limits_x) {
  # Plot the point estimates of R_eff (R_hat)
  R_hat_hist <- ggplot(
    df_R_hat,
    aes(x = R, fill = model, alpha = window_len_fct)
  ) +
    geom_histogram(color = "black") +
    geom_vline(xintercept = R_true, color = "red") +
    scale_fill_manual(values = model_colors) +
    scale_alpha_manual(values = c("short" = 0.6, "long" = 0.3)) +
    facet_wrap(~model) +
    labs(
      x = expression(hat(R)),
      fill = "Model",
      y = NULL,
      alpha = "Estimation\nwindow length"
    ) +
    xlim(limits_x$R_hat) # A couple of values might get clipped

  # Plot the standard errors of the R_eff estimates
  se_hat_hist <- ggplot(
    df_R_hat,
    aes(x = se, fill = model, alpha = window_len_fct)
  ) +
    geom_histogram(color = "black") +
    scale_fill_manual(values = model_colors) +
    scale_alpha_manual(values = c("short" = 0.6, "long" = 0.3)) +
    facet_wrap(~model) +
    labs(
      x = expression(se(hat(R))),
      fill = "Model",
      y = NULL,
      alpha = "Estimation\nwindow length"
    ) +
    xlim(limits_x$se_hat) # A couple of values might get clipped
  list(R_hat = R_hat_hist, se_hat = se_hat_hist)
}

#' Plots nominal vs empirical coverage
#'
#' @description This function plots nominal vs. empirical coverage levels with
#'   a 45-degree reference line, faceted by model.
#' @param R_eff the true value of the effective reproduction number used to
#'   generate the trajectory
#' @param nb_size positive real value, the size parameter of the negative
#'   binomial distribution for "NegBin-Q" and "NegBin-L"
#' @param df_R_hat a data frame with R estimates, after replacing the divergent
#'   runs by NAs, containing columns
#'   \code{R_hat}, \code{se_hat} and \code{converged}
#' @param X a vector, the incidence, can be dropped when we get rid
#'   of the normal approximation part
#' @param Lambda vector, the single covariate, can be dropped when we get rid
#'   of the normal approximation part
#' @param window the length of the estimation window, can be dropped when we get
#'   rid of the normal approximation part
#' @param nominal_covr a vector of the nominal coverage levels
#' @return a patchwork plot with 3 panels - first 1000 generated trajectories,
#'   nominal vs. empirical coverage for the short estimation window and
#'   nominal vs. empirical coverage for the long estimation window
plot_coverage <- function(
  R_eff,
  nb_size,
  df_R_hat,
  X,
  Lambda,
  nominal_covr,
  model_colors
) {
  df_coverage <- create_coverage_df(
    R_eff,
    nb_size,
    df_R_hat,
    X,
    Lambda,
    nominal_covr
  )

  ggplot() +
    geom_line(
      data = df_coverage$lines,
      mapping = aes(x = covr_nominal, y = covr_empirical, color = model),
      linewidth = 1.5,
      alpha = 0.6
    ) +
    geom_point(
      data = df_coverage$points,
      mapping = aes(x = covr_nominal, y = covr_empirical, shape = type)
    ) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = model_colors) +
    scale_shape_manual(
      values = c("normal approx." = 3, "theoretical Poisson coverage" = 4),
      labels = c(
        "normal approx." = "normal approx.",
        "theoretical Poisson coverage" = "theoretical\nPoisson coverage"
      )
    ) +
    labs(
      x = "Nominal coverage",
      y = "Empirical coverage",
      color = "Model",
      shape = NULL
    ) +
    facet_wrap(~window_len_fct)
}

#' Plots the fan-out trajectories
#'
#' @description This function plots the fan-out of the trajectory from its fixed
#'   initial values. The next step is plotting the trajectory for the longer
#'   estimation windows and making the part corresponding to the shorter window
#'   lighter/darker.
#' @param X a matrix of the simulated counts, one column per simulation run
#' @return a ggplot object
plot_trajectories <- function(X, short_window, n_init) {
  df_trajectories <- reshape2::melt(
    # Display only the first 1000 simulation runs to make the plot readable
    X[, seq_len(min(1000, ncol(X)))],
    varnames = c("day", "trajectory"),
    value.name = "cases"
  ) |>
    mutate(
      segment = dplyr::case_when(
        day < n_init ~ "Initial",
        day >= n_init & day < short_window + n_init ~ "Short window",
        day >= short_window + n_init ~ "Long window"
      )
    )
  ggplot() +
    geom_line(
      data = df_trajectories,
      mapping = aes(
        x = day,
        y = cases,
        group = trajectory,
        color = segment,
        alpha = segment
      )
    ) +
    ggpubr::geom_bracket(
      xmin = n_init + 1,
      xmax = nrow(X),
      y.position = max(df_trajectories$cases) + 10,
      label = "Long estimation window"
    ) +
    ggpubr::geom_bracket(
      xmin = n_init + 1,
      xmax = n_init + short_window,
      y.position = max(df_trajectories$cases),
      label = "Short estimation window"
    ) +
    scale_color_manual(
      values = c(
        "Initial" = "black",
        "Short window" = "black",
        "Long window" = "gray30"
      ),
      guide = "none"
    ) +
    scale_alpha_manual(
      values = c(
        "Initial" = 1,
        "Short window" = 0.1,
        "Long window" = 0.1
      ),
      guide = "none"
    ) +
    labs(x = "Day", y = "Cases")
}
