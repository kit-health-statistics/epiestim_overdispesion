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
    "tidyr",
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

# If we plot only half of the scenarios, by how much do we divide the height?
plot_halving_coeff <- 1.98

# Scenario blocks for the static branching:
# NegBin-L
# NegBin-L with weekday effects
# NegBin-Q
outer_scenarios <- data.frame(
  distribution = c("NegBin-L", "NegBin-L", "NegBin-Q", "Poiss"),
  weekday_effect = c("weekday_no", "weekday_yes", "weekday_no", "weekday_no")
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
  # Static branching over 4 scenario blocks defined in `outer_scenarios`
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
      df_R_hat_raw,
      create_results_df(
        trajectories$X,
        trajectories$Lambda,
        global_params$short_window,
        global_params$long_window
      ),
      pattern = map(trajectories),
      iteration = "list"
    ),
    # Replace the divergent NegBin estimates by the Poison ones
    tar_target(
      df_R_hat,
      replace_divergent(df_R_hat_raw),
      pattern = map(df_R_hat_raw),
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
    # How many did converge?
    tar_target(
      convergence_tables,
      summarize_convergence(scenarios$scenario_id, df_R_hat),
      pattern = map(df_R_hat, scenarios),
      iteration = "list"
    ),
    # Save the convergence table
    tar_target(
      saved_tables,
      {
        convergence <- convergence_tables |>
          purrr::map("df_convergence") |>
          dplyr::bind_rows()
        unstable <- convergence_tables |>
          purrr::map("df_unstable") |>
          dplyr::bind_rows()
        write_csv(
          convergence,
          paste0("tables/", scenario_id, "_convergence.csv")
        )
        write_csv(unstable, paste0("tables/", scenario_id, "_unstable.csv"))
      }
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
      if (distribution == "Poiss") {
        # For Poisson, we have only half the scenarios as for the rest, so the
        # height of the resulting plot must be divided by 2. We divide by less
        # than 2 to allow for some space for the title.
        plot_height <- plot_size$height / plot_halving_coeff
      } else {
        plot_height <- plot_size$height
      }

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
        height = plot_height
      )
      # Save the distribution of the estimates
      p_densities <- compose_dens_patches(
        plot_density_panels,
        purrr::map(plot_panels, "meta"),
        global_params$short_window,
        global_params$long_window
      )
      save_plot(
        p_densities,
        paste(scenario_id, "Rhat_density", sep = "_"),
        width = plot_size$width,
        height = plot_height
      )
    }),
    tar_target(
      saved_figures_poster,
      {
        # Create the coverage plot for a subset of scenarios. We will use only
        # NegBin-L and NegBin-Q versions, but it's easier to generate everything
        # in one go.
        p_simulation_poster <- compose_coverage_patches(
          plot_panels[scenarios$R_eff == 1.5 & scenarios$magnitude == "low"],
          global_params$short_window,
          global_params$long_window,
          panel_widths = c(1.3, 2, 3)
        ) &
          theme(plot.margin = unit(c(5.5, 5.5, 0.5, 5.5), "points"))
        save_plot(
          p_simulation_poster,
          paste(scenario_id, "simulation_coverage_poster", sep = "_"),
          # Figure dimensions must be set manually to fit the poster page
          width = plot_size$width / 1.4,
          height = (plot_size$height / 1.4) / 2.35
        )
      }
    )
  )
)
