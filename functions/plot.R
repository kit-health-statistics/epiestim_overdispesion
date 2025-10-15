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
#' @return a list containing 2 ggplots of coverage
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

  # Split the data frames to create 2 separate plots instead of facets. This
  # way, it's easier to put the final plot together from multiple blocks
  df_split_lines <- split.data.frame(
    df_coverage$lines,
    df_coverage$lines$window_len_fct
  )
  df_split_points <- split.data.frame(
    df_coverage$points,
    df_coverage$points$window_len_fct
  )
  p_coverage <- vector("list", 2)
  names(p_coverage) <- names(df_split_lines)
  for (k in 1:2) {
    p_coverage[[k]] <- ggplot() +
      geom_line(
        data = df_split_lines[[k]],
        mapping = aes(x = covr_nominal, y = covr_empirical, color = model),
        linewidth = 1.5,
        alpha = 0.6
      ) +
      geom_point(
        data = df_split_points[[k]],
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
      )
  }
  p_coverage
}

#' Plots the fan-out trajectories
#'
#' @description This function plots the fan-out of the trajectory from its fixed
#'   initial values. The next step is plotting the trajectory for the longer
#'   estimation windows and making the part corresponding to the shorter window
#'   lighter/darker.
#' @param X a matrix of the simulated counts, one column per simulation run
#' @param short_window an integer, the length of the short estimation window
#' @param n_init an integer, the number of initial values of the trajectory
#' @return a ggplot object
plot_trajectories <- function(X, short_window, n_init) {
  df_trajectories <- reshape2::melt(
    # Display only the first 100 simulation runs to make the plot readable
    X[, seq_len(min(100, ncol(X)))],
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

  # Set the offset for the brackets
  max_cases <- max(df_trajectories$cases)
  bracket_offset <- max_cases * 0.08  # 8% of the max value

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
      y.position = max_cases + 2 * bracket_offset,
      label = "Long window",
      label.size = 3
    ) +
    ggpubr::geom_bracket(
      xmin = n_init + 1,
      xmax = n_init + short_window,
      y.position = max_cases + bracket_offset,
      label = "Short window",
      label.size = 3
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
    coord_cartesian(ylim = c(0, max_cases + 3 * bracket_offset)) +
    labs(x = "Day", y = "Cases")
}

#' Plots the metadata of the simulation scenario
#'
#' @description This function plots the information about the true parameter
#'   values used in the simulation scenario.
#' @param R_eff a positive real, the true value of the effective reproduction
#'   number
#' @param dispersion string values, either "high", or "low"
#' @param magnitude string values, either "high", or "low"
#' @return a ggplot object
plot_metadata <- function(R_eff, dispersion, magnitude) {
  df_text <- data.frame(
    x = rep(1, 3),
    y = c(-1, 0, 1),
    label = c(
      paste0("R = ", R_eff),
      paste0("Dispersion: ", gsub("_.*", "", dispersion)),
      paste0("Magnitude: ", gsub("_.*", "", magnitude))
    )
  )
  ggplot(df_text, aes(x = x, y = y, label = label)) +
    geom_text(hjust = 0) +
    coord_cartesian(ylim = c(-3, 3), xlim = c(0.995, 1.02)) +
    theme_void()
}

#' Compose the patches into the final plot
#'
#' @description This function composes plots from all the scenario runs
#' together, including titles and labels, using the patchwork approach.
#' @param plot_panels a list of lists of patchwork patches. The number of
#'   simulation scenarios is the length of the outer list. Each inner list
#'   contains 4 plots - simulation trajectories, coverage for the short and
#'   long estimation window and a plot with parameter values.
#' @param short_window an integer, the length of the short estimation window
#' @param long_window an integer, the length of the long estimation window
#' @return a patchwork plot
compose_patches <- function(plot_panels, short_window, long_window) {
  # We will have as many rows as the simulation scenarios
  n_rows <- length(plot_panels)

  # Set the first row as a ggplot with text only
  p_trajectories <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "Simulation trajectories"), size = 5) +
    theme_void()
  p_meta <- plot_spacer()
  p_coverage <- ggplot() +
    geom_text(
      aes(
        x = 1,
        y = 1,
        label = paste0("Window width: ", short_window)
      ),
      size = 5
    ) +
    theme_void() +
    ggplot() +
    geom_text(
      aes(
        x = 1,
        y = 1,
        label = paste0("Window width: ", long_window)
      ),
      size = 5
    ) +
    theme_void()

  for (k in seq_len(n_rows)) {
    p_trajectories <- p_trajectories / plot_panels[[k]]$trajectories
    p_coverage <- p_coverage + plot_panels[[k]]$short + plot_panels[[k]]$long
    p_meta <- p_meta / plot_panels[[k]]$meta
  }

  p_trajectories <- p_trajectories +
    plot_layout(heights = c(1, rep(6, n_rows)), axes = "collect")
  p_meta <- p_meta + plot_layout(heights = c(1, rep(6, n_rows)))

  p_coverage <- p_coverage +
    plot_layout(
      nrow = n_rows + 1,
      ncol = 2,
      heights = c(1, rep(6, n_rows)),
      guides = "collect",
      axes = "collect"
    )
  p_final <- (p_meta | p_trajectories | p_coverage) +
    plot_layout(guides = "collect", widths = c(1, 3, 2))
  p_final
}
