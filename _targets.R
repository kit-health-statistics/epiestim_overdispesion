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

# Set the global objects =======================================================
model_colors <- c(
  "Poiss" = "#009E73",
  "Q-Poiss" = "#F0E442",
  "NegBin-L" = "#56B4E9",
  "NegBin-Q" = "#CC79A7"
)

# Parameters common for all simulation runs
global_params <- list(
  mean_si = 6,
  std_si = 1,
  n_init = 14,
  short_window = 7,
  long_window = 14,
  n_sim = 1000L,
  weekday = c(1.05, 0.24, 1.4, 1.69, 0.17, 1.4, 1.05),
  base_seed = 9786L
)

# Size of the final plot
plot_size <- list(width = 14, height = 14.5)

# Scenario blocks for the static branching:
# NegBin-L
# NegBin-L with weekday effects
# NegBin-Q
outer_scenarios <- data.frame(
  distribution = c("NegBin-L", "NegBin-L", "NegBin-Q"),
  weekday_effect = c("weekday_no", "weekday_yes", "weekday_no")
) |>
  dplyr::mutate(
    scenario_id = paste(distribution, weekday_effect, sep = "_")
  )

# Define pipeline ==============================================================
list(
  # Serial interval, common for all scenarios
  tar_target(
    si,
    with(global_params, discr_si(seq_len(n_init), mean_si, std_si))
  ),
  # Static branching over 3 scenario blocks defined in `outer_scenarios`
  tar_map(
    unlist = FALSE,
    values = outer_scenarios,
    names = scenario_id,
    # Get the data frame with scenario parameters
    tar_target(scenarios, create_scenario_grid(distribution)),
    # Initial values. Sample iid counts using a reasonable data generating
    # mechanism.
    tar_target(
      init,
      {
        # Set seed ensuring identical initialization for scenarios with the same
        # magnitude. Creates redundant copies (nrow(scenarios) instead of 2) but
        # simplifies pipeline structure with negligible computational cost.
        set.seed(global_params$base_seed + scenarios$init_seed)
        pmax(
          0,
          round(
            rnorm(
              global_params$n_init,
              scenarios$init_magnitude,
              scenarios$init_sd
            )
          )
        )
      },
      pattern = map(scenarios),
      iteration = "list"
    ),
    # Simulations
    tar_target(
      trajectories,
      with(
        global_params,
        generate_trajectories(
          n_sim = n_sim,
          init,
          R_eff = scenarios$R_eff,
          si = si,
          lgt = n_init + long_window,
          model = distribution,
          nb_size = scenarios$nb_size,
          seed = global_params$base_seed + scenarios$scenario_number,
          # If the weekday effect is not present, we multiply the simulated mean
          # value by 1 under the hood, so there is indeed no effect
          weekday_effect = if (weekday_effect == "weekday_yes") {
            global_params$weekday
          } else {
            rep(1, global_params$n_init + global_params$long_window)
          }
        )
      ),
      pattern = map(init, scenarios),
      iteration = "list"
    ),
    # Estimation
    tar_target(
      df_R_hat,
      create_results_df(
        trajectories$X,
        trajectories$Lambda,
        global_params$short_window,
        global_params$long_window
      ),
      pattern = map(trajectories),
      iteration = "list"
    ),
    # Plot the trajectories and the coverage
    tar_target(
      plot_panels,
      list(
        trajectories = plot_trajectories(
          trajectories$X,
          global_params$short_window,
          global_params$n_init
        ),
        coverage = plot_coverage(
          scenarios$R_eff,
          scenarios$nb_size,
          df_R_hat,
          trajectories$X[-seq_len(global_params$n_init), ],
          trajectories$Lambda[-seq_len(global_params$n_init), ],
          seq(0, 1, by = 0.05),
          distribution,
          model_colors
        ),
        meta = plot_metadata(
          scenarios$R_eff,
          scenarios$nb_size,
          scenarios$magnitude,
          distribution
        )
      ),
      pattern = map(trajectories, df_R_hat, scenarios),
      iteration = "list"
    ),
    # Plot the density of R estimates and its SEs
    tar_target(
      plot_density_panels,
      plot_dens(
        df_R_hat,
        scenarios$R_eff,
        model_colors,
        get_xlim(scenarios$R_eff)
      ),
      pattern = map(df_R_hat, scenarios),
      iteration = "list"
    ),
    # Save plots
    tar_target(saved_figures, {
      # Save the coverage plot
      p_simulation <- compose_coverage_patches(
        plot_panels,
        global_params$short_window,
        global_params$long_window
      )
      save_plot(
        p_simulation,
        paste(scenario_id, "simulation_coverage", sep = "_"),
        width = plot_size$width,
        height = plot_size$height
      )
      # Save the distribution of the estimates
      p_densities <- compose_dens_patches(
        plot_density_panels,
        map(plot_panels, "meta"),
        global_params$short_window,
        global_params$long_window
      )
      save_plot(
        p_densities,
        paste(scenario_id, "Rhat_density", sep = "_"),
        width = plot_size$width,
        height = plot_size$height
      )
    })
  )
)
