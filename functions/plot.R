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
  # Plot the standard errors of the R_eff estimates
  R_hat_hist <- ggplot(df_R_hat, aes(x = R, fill = model)) +
    geom_histogram(color = "black") +
    geom_vline(xintercept = R_true, color = "red") +
    scale_fill_manual(values = model_colors) +
    facet_wrap(~model) +
    labs(x = expression(hat(R)), fill = "Model", y = NULL) +
    xlim(limits_x$R_hat)  # A couple of values might get clipped

  # Plot the standard errors of the R_eff estimates
  se_hat_hist <- ggplot(df_R_hat, aes(x = se, fill = model)) +
    geom_histogram(color = "black") +
    scale_fill_manual(values = model_colors) +
    facet_wrap(~model) +
    labs(x = expression(se(hat(R))), fill = "Model", y = NULL) +
    xlim(limits_x$se_hat)  # A couple of values might get clipped
  list(R_hat = R_hat_hist, se_hat = se_hat_hist)
}

#' Plots nominal vs empirical coverage
#'
#' @description This function plots histograms of the R estimates, highlighting
#'   the true value, and histograms of the standard error estimates
#' @param df_coverage the data frame containing the nominal and the empirical
#'   coverage levels, returned by \code{create_coverage_df()}
#' @return a ggplot object, 4 facets displaying the nominal vs. empirical
#'   coverage
plot_coverage <- function(df_coverage) {
  ggplot(df_coverage, aes(x = covr_nominal, y = covr_empirical, color = type)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(
      labels = c(
        "'real'" = "'real'",
        "GLM" = "GLM",
        "normal approx." = "normal\napprox."
      )
    ) +
    labs(x = "nominal coverage", y = "empirical coverage", color = "coverage") +
    facet_wrap(~model)
}

#' Plots the fan-out trajectories
#'
#' @description This function plots the fan-out of the trajectory from its fixed
#'   initial values. The next step is plotting the trajectory for the longer
#'   estimation windows and making the part corresponding to the shorter window
#'   lighter/darker.
#' @param X a matrix of the simulated counts, one column per simulation run
#' @return a list containing the two histograms named as \code{R_hat} and
#'   \code{se_hat}
plot_trajectories <- function(X) {
  df_trajectories <- reshape2::melt(
    X,
    varnames = c("day", "trajectory"),
    value.name = "cases"
  )
  ggplot(df_trajectories, aes(x = day, y = cases, group = trajectory)) +
    geom_line(alpha = 0.05, color = "gray")
}
