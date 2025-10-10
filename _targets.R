library(targets)
library(tarchetypes)
library(qs2)

# Set targets options
tar_option_set(
  packages = c(
    "ggplot2",
    "patchwork",
    "purrr",
    "here",
    "dplyr",
    "qs2",
    "EpiEstim",
    "gamlss",
    "readr"
  ),
  format = "qs", # Use qs format (qs2 is used via repository option)
  memory = "transient", # Free memory after each target completes
  garbage_collection = TRUE, # Run garbage collection
  repository = "local", # Use qs2 backend for storage
  error = "continue" # Continue pipeline when targets fail
)
tar_source(files = "functions")

# set the ggplot theme
ggplot2::theme_set(ggplot2::theme_bw())

# Set the globals
model_colors <- c(
  "Poiss" = "#009E73",
  "Q-Poiss" = "#F0E442",
  "NegBin-L" = "#56B4E9",
  "NegBin-Q" = "#CC79A7"
)

# Define pipeline
list(
  tar_target(
    params,
    list(
      init_magnitude = 50,
      nb_size = 0.3,
      var_infl_factor_true = 1 + 1 / 0.3,
      R_eff = 1.1,
      n_init = 14,
      short_window = 7,
      long_window = 14,
      n_sim = 1e4,
      mean_si = 3,
      std_si = 1
    )
  ),

  # Serial interval
  tar_target(
    si,
    discr_si(seq_len(params$n_init), params$mean_si, params$std_si)
  ),

  # Initial values
  tar_target(init, {
    set.seed(1999)
    round(params$init_magnitude * runif(params$n_init))
  }),

  # Simulations
  tar_target(
    trajectories,
    with(
      params,
      replicate(
        n_sim,
        simulate_renewal(
          init,
          R_eff,
          si,
          n_init + long_window,
          model = "NegBin-L",
          nb_size = nb_size
        ),
        simplify = FALSE
      )
    )
  ),

  # Extract observations and means
  tar_target(X, do.call(cbind, map(trajectories, "X"))),
  tar_target(Lambda, do.call(cbind, map(trajectories, "Lambda"))),

  # Estimation
  tar_target(
    df_R_hat_raw,
    create_results_df(X, Lambda, params$short_window, params$long_window)
  ),

  # Remove unstable estimates
  tar_target(df_R_hat, remove_divergent(df_R_hat_raw)),

  # Set x-axis limits for the histograms
  tar_target(limits_x, list(R_hat = c(0.5, 1.5), se_hat = c(0, NA))),

  # Plot histograms
  tar_target(
    p_hists,
    plot_hists(
      df_R_hat,
      params$R_eff,
      model_colors = model_colors,
      limits_x = limits_x
    )
  ),

  # Plot the trajectories and coverage
  tar_target(
    p_coverage,
    (
      plot_trajectories(X, params$short_window, params$n_init) +
        plot_coverage(
          params$R_eff,
          params$nb_size,
          df_R_hat,
          X[-seq_len(params$n_init), ],
          Lambda[-seq_len(params$n_init), ],
          seq(0, 1, by = 0.05),
          model_colors
        )
    ) +
      plot_layout(widths = c(1, 1))
  ),

  # Save plots
  tar_target(saved_figures, {
    # Ensure directory exists
    dir.create("figure", showWarnings = FALSE, recursive = TRUE)
    ggsave(
      "figure/R_hat_histogram_test.pdf",
      p_hists$R_hat,
      width = 8,
      height = 6
    )
    ggsave(
      "figure/se_hat_histogram_test.pdf",
      p_hists$se_hat,
      width = 8,
      height = 6
    )
    ggsave("figure/coverage_traj_test.pdf", p_coverage, width = 16, height = 4)
  })
)
