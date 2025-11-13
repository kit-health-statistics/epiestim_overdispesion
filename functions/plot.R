#' Returns the x-axis limits for the R estimate distribution plot
#'
#' @description This function returns the x-axis limits for the plot of R
#'   estimate distribution and the plot of its standard errors. The limits are
#'   returned as a list
#' @param R_eff positive real, the true value of the effective reproductive
#'   number. The x-axis is centered around this value.
#' @param distribution the true underlying count distribution
#' @return a list of 2 vectors of 2 elements - the limits for R estimate and
#'   the limits for its standard errors.
get_xlim <- function(R_eff, distribution) {
  se_limits <- if (distribution == "NegBin-L") {
    c(0, 0.6)
  } else if (distribution == "NegBin-Q") {
    c(0, 0.4)
  } else if (distribution == "Poiss") {
    c(0, 0.3)
  }
  list(R_hat = R_eff + c(-0.8, 0.8), se_hat = se_limits)
}

#' Plots density distributions of the estimates
#'
#' @description This function plots density distributions of the R estimates,
#'   highlighting the true value, and density distributions of the standard
#'   error estimates
#' @param df_R_hat a data frame with R estimates, after replacing the divergent
#'   runs by NAs, containing columns
#' @param R_true the true value of the effective reproduction number used to
#'   generate the trajectory
#' @param model_colors a named vector specifying the model colors. Must have
#'   names "Poiss", "Q-Poiss", "NegBin-L" and "NegBin-Q"
#' @param limits_x a list of the limits for the x-axis, must have 2 elements
#'   called \code{R_hat} and \code{se_hat}. Both elements must contain a vector
#'   with 2 elements - the axis range.
#' @return a nested list containing density plots: outer list has elements
#'   \code{R_hat} and \code{se_hat}, each containing a list of two plots
#'   (one per window length)
plot_dens <- function(df_R_hat, R_true, model_colors, limits_x) {
  # Split the data frames to create 2 separate plots for each window width. This
  # way, it's easier to put the final plot together from multiple blocks
  df_R_hat_split <- split.data.frame(
    df_R_hat,
    df_R_hat$window_len_fct
  )

  p_R_hat <- p_se_hat <- vector("list", 2)
  names(p_R_hat) <- names(p_se_hat) <- names(df_R_hat_split)
  for (k in 1:2) {
    # Plot the point estimates of R_eff (R_hat)
    p_R_hat[[k]] <- ggplot(df_R_hat_split[[k]], aes(x = R, color = model)) +
      geom_line(stat = "density", linewidth = 1, alpha = 0.6, na.rm = TRUE) +
      geom_vline(aes(xintercept = R_true, linetype = "R_true"), color = "red") +
      scale_color_manual(values = model_colors) +
      scale_linetype_manual(
        values = c("R_true" = "dashed"),
        labels = c("R_true" = expression(R[t]))
      ) +
      labs(
        x = expression(hat(R)),
        color = "Model",
        linetype = NULL,
        y = "density"
      ) +
      # A couple of values might get clipped.
      coord_cartesian(xlim = limits_x$R_hat) +
      theme(legend.text = element_text(size = 10))

    # Plot the standard errors of the R_eff estimates
    p_se_hat[[k]] <- ggplot(df_R_hat_split[[k]], aes(x = se, color = model)) +
      geom_line(
        stat = "density",
        linewidth = 1,
        alpha = 0.6,
        na.rm = TRUE,
        bounds = c(0, Inf)
      ) +
      scale_color_manual(values = model_colors) +
      labs(
        x = expression(widehat(se)(hat(R))),
        color = "Model",
        y = "density"
      ) +
      # A couple of values might get clipped.
      coord_cartesian(xlim = limits_x$se_hat) +
      theme(legend.text = element_text(size = 10))
  }
  # Return the list
  list(R_hat = p_R_hat, se_hat = p_se_hat)
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
#' @param window the length of the estimation window, can be dropped when we get
#'   rid of the normal approximation part
#' @param nominal_covr a vector of the nominal coverage levels
#' @return a list containing 2 ggplots of coverage
plot_coverage <- function(
  R_eff,
  nb_size,
  df_R_hat,
  nominal_covr,
  distribution,
  weekday_effect,
  model_colors
) {
  df_coverage <- create_coverage_df(
    R_eff,
    nb_size,
    df_R_hat,
    nominal_covr,
    distribution,
    weekday_effect
  )

  # Split the data frames to create 2 separate plots instead of facets. This
  # way, it's easier to put the final plot together from multiple blocks
  df_split <- split.data.frame(
    df_coverage,
    df_coverage$window_len_fct
  )

  # Don't show the legend for the theoretical Poisson coverage for the
  # "NegBin-Q" distribution
  if (distribution == "NegBin-L") {
    linetype_guide <- "legend"
  } else {
    linetype_guide <- "none"
  }

  p_coverage <- vector("list", 2)
  names(p_coverage) <- names(df_split)
  for (k in 1:2) {
    p_coverage[[k]] <- ggplot() +
      geom_abline(intercept = 0, slope = 1, color = "grey40") +
      geom_line(
        data = df_split[[k]],
        mapping = aes(
          x = covr_nominal,
          y = covr_empirical,
          color = model,
          linetype = type,
          linewidth = type,
          alpha = type
        )
      ) +
      scale_color_manual(
        values = c(model_colors, "theoretical Poisson" = "black"),
        breaks = names(model_colors),
        labels = names(model_colors)
      ) +
      scale_linetype_manual(
        values = c(
          "Empirical coverage" = "solid",
          "theoretical Poisson coverage" = "dashed"
        ),
        breaks = "theoretical Poisson coverage",
        labels = c(
          "theoretical Poisson coverage" = "theoretical\nPoisson\ncoverage"
        ),
        guide = linetype_guide
      ) +
      scale_linewidth_manual(
        values = c(
          "Empirical coverage" = 1.5,
          "theoretical Poisson coverage" = 0.5
        ),
        guide = "none"
      ) +
      scale_alpha_manual(
        values = c(
          "Empirical coverage" = 0.6,
          "theoretical Poisson coverage" = 1
        ),
        guide = "none"
      ) +
      scale_x_continuous(
        breaks = seq(0, 1, by = 0.2),
        minor_breaks = seq(0, 1, by = 0.1)
      ) +
      scale_y_continuous(
        breaks = seq(0, 1, by = 0.2),
        minor_breaks = seq(0, 1, by = 0.1)
      ) +
      labs(
        x = "Nominal coverage",
        y = "Empirical coverage",
        color = "Model",
        shape = NULL,
        linetype = NULL
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

  # The long window is the length of the trajectory minus the initialization
  long_window <- nrow(X) - n_init

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
      label = paste0(long_window, "-day window"),
      label.size = 3
    ) +
    ggpubr::geom_bracket(
      xmin = n_init + 1,
      xmax = n_init + short_window,
      y.position = max_cases + 0.5 * bracket_offset,
      label = paste0(short_window, "-day window"),
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
#' @param nb_size a positive real, the true value of the dispersion parameter
#' @param magnitude string values, either "high", or "low"
#' @param distribution string value indicating the count distribution:
#'   "NegBin-L", "NegBin-Q", or "Poiss"
#' @return a ggplot object
plot_metadata <- function(R_eff, nb_size, magnitude, distribution) {
  df_text <- data.frame(
    x = rep(1, 3),
    y = c(1, 0, -1),
    label = c(
      paste0("R[t] == ", R_eff),
      if (distribution == "NegBin-L") {
        paste0("xi == ", (1 + 1 / nb_size))
      } else if (distribution == "NegBin-Q") {
        paste0("psi == ", 1 / nb_size)
      } else {
        NA
      },
      paste0("Magnitude: ", gsub("_.*", "", magnitude))
    )
  )
  # Remove the dispersion parameter, when it's not present for the Poisson
  # distribution.
  if (distribution == "Poiss") {
    df_text <- df_text |> filter(!is.na(label))
  }
  ggplot(df_text, aes(x = x, y = y, label = label)) +
    geom_text(hjust = 0, parse = TRUE) +
    coord_cartesian(ylim = c(-3, 3), xlim = c(0.997, 1.02)) +
    theme_void()
}

#' Compose coverage plot patches into the final plot
#'
#' @description This function composes coverage and trajectory plots from all
#'  the scenario runs together, including titles and labels, using the patchwork
#'  approach.
#' @param plot_panels a list of lists of patchwork patches. The number of
#'   simulation scenarios is the length of the outer list. Each inner list
#'   contains 4 plots - simulation trajectories, coverage for the short and
#'   long estimation window and a plot with parameter values.
#' @param short_window an integer, the length of the short estimation window
#' @param long_window an integer, the length of the long estimation window
#' @param panel_widths vector with 3 elements specifying the width of the
#'   vertical plot panels - parameter values (metadata), trajectories and the
#'   coverage panel
#' @return a patchwork plot
compose_coverage_patches <- function(
  plot_panels,
  short_window,
  long_window,
  panel_widths = c(1, 3, 2)
) {
  # We will have as many rows as the simulation scenarios
  n_rows <- length(plot_panels)

  # Set the first row as a ggplot with text only
  p_trajectories <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "Simulation trajectories"), size = 5) +
    theme_void()
  p_meta <- plot_spacer()

  # Compose the right panel corresponding to the empirical vs. nominal coverage
  p_coverage <- compose_subplot_by_windows(
    map(plot_panels, "coverage"),
    short_window,
    long_window
  )

  # Glue the plots under each other into vertical strips. `p_coverage` is
  # composed of two columns containing the coverage plots for both estimation
  # windows in order to allow or collecting the guides and axes.
  for (k in seq_len(n_rows)) {
    p_trajectories <- p_trajectories / plot_panels[[k]]$trajectories
    p_meta <- p_meta / plot_panels[[k]]$meta
  }

  # Adjust the layout of the vertical strips
  p_trajectories <- p_trajectories +
    plot_layout(heights = c(1, rep(6, n_rows)), axes = "collect")
  p_meta <- p_meta + plot_layout(heights = c(1, rep(6, n_rows)))
  # Collect the guides
  p_coverage <- p_coverage +
    plot_layout(
      nrow = n_rows + 1,
      ncol = 2,
      heights = c(1, rep(6, n_rows)),
      guides = "collect",
      axes = "collect"
    ) &
    theme(
      legend.text = element_text(size = 11)
    )

  # Create the final plot
  p_final <- (p_meta | p_trajectories | p_coverage) +
    plot_layout(guides = "collect", widths = panel_widths)
  p_final
}

#' Compose plot patches into two columns by the window length
#'
#' @description This function composes those plot patches, that have 2 versions
#'   based on the window length. Short window is on the left, long window on the
#'   right.
#' @param subplot_panels a list of lists of patchwork patches. The number of
#'   simulation scenarios is the length of the outer list. Each inner list
#'   contains 2 plots - versions of the same plot, once with the short and once
#'   with the long estimation window. The inner list elements must be
#'   accordingly named as "short" and "long".
#' @param short_window an integer, the length of the short estimation window
#' @param long_window an integer, the length of the long estimation window
#' @return a patchwork plot
compose_subplot_by_windows <- function(
  subplot_panels,
  short_window,
  long_window
) {
  n_rows <- length(subplot_panels)
  # Create header row with two column titles using patchwork's + operator
  p <- ggplot() +
    geom_text(
      aes(x = 1, y = 1, label = paste0(short_window, "-day window")),
      size = 5
    ) +
    theme_void() +
    ggplot() +
    geom_text(
      aes(x = 1, y = 1, label = paste0(long_window, "-day window")),
      size = 5
    ) +
    theme_void()

  # Glue the plots under each other into vertical strips. The patchwork plot is
  # composed of two columns containing the plots of desired quantities for both
  # estimation windows in order to allow or collecting the guides and axes.
  for (k in seq_len(n_rows)) {
    p <- p + subplot_panels[[k]]$short + subplot_panels[[k]]$long
  }
  p <- p +
    # Collect the guides
    plot_layout(
      nrow = n_rows + 1,
      ncol = 2,
      heights = c(1, rep(6, n_rows)),
      guides = "collect",
      axes = "collect"
    ) &
    theme(
      legend.text = element_text(size = 11)
    )
  p
}

#' Compose density plot patches into the final plot
#'
#' @description This function composes density plots from all the scenario runs
#' together, including titles and labels, using the patchwork approach.
#' @param plot_panels a list of lists of patchwork patches. The number of
#'   simulation scenarios is the length of the outer list. Each inner list
#'   contains 2 lists - the density plot of R estimate and the density plot of
#'   its SE. Each list than contains one plot for the short and one plot for the
#'   long estimation window.
#' @param plot_meta_panels a list containing plots of the parameter values used
#'   in the simulation scenario.
#' @param short_window an integer, the length of the short estimation window
#' @param long_window an integer, the length of the long estimation window
#' @return a patchwork plot
compose_dens_patches <- function(
  plot_panels,
  plot_meta_panels,
  short_window,
  long_window
) {
  # We will have as many rows as the simulation scenarios
  n_rows <- length(plot_panels)

  # Compose the density plots of R estimates
  p_R_hat <- compose_subplot_by_windows(
    map(plot_panels, "R_hat"),
    short_window,
    long_window
  )
  # Compose the density plots of SEs of R estimates
  p_se_hat <- compose_subplot_by_windows(
    map(plot_panels, "se_hat"),
    short_window,
    long_window
  )

  # Compose the panel of the metadata about the parameter values
  p_meta <- plot_spacer()
  for (k in seq_len(n_rows)) {
    p_meta <- p_meta / plot_meta_panels[[k]]
  }
  p_meta <- p_meta + plot_layout(heights = c(1, rep(6, n_rows)))

  # Extract the legend and remove it from the plots
  p_legend <- wrap_elements(ggpubr::get_legend(plot_panels[[1]]$R_hat$short))
  p_R_hat <- wrap_elements(
    (p_meta | p_R_hat) +
      plot_layout(widths = c(2, 6)) &
      theme(legend.position = "none")
  )
  p_se_hat <- wrap_elements(
    (plot_spacer() | p_se_hat) +
      plot_layout(widths = c(0.25, 6)) & theme(legend.position = "none")
  )

  # Create the upper titles
  p_headline_R_hat <- plot_spacer() +
    ggplot() +
    geom_text(
      aes(x = 1, y = 1, label = "Distribution~of~hat(R)"),
      size = 5.5,
      parse = TRUE
    ) +
    theme_void() +
    plot_layout(widths = c(2, 6))
  p_headline_se_hat <- plot_spacer() +
    ggplot() +
    geom_text(
      aes(x = 1, y = 1, label = "Distribution~of~widehat(se)(hat(R))"),
      size = 5.5,
      parse = TRUE
    ) +
    theme_void() +
    plot_layout(widths = c(0.25, 6))

  # Combine into 2 big panels
  p_headline <- (p_headline_R_hat | p_headline_se_hat) +
    plot_layout(widths = c(6, 5))
  p_plot <- (p_R_hat | p_se_hat) +
    plot_layout(widths = c(6, 5))

  # Combine everything
  p_main <- wrap_elements(
    (p_headline / p_plot) + plot_layout(heights = c(1, 22))
  ) +
    p_legend +
    plot_layout(widths = c(10, 1))
  p_main
}
