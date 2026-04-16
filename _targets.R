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
  "Q-Poiss" = "#E69F00",
  "NegBin-L" = "#56B4E9",
  "NegBin-Q" = "#CC79A7"
)

# Parameters common for all simulation runs
global_params <- list(
  short_window = 7,
  long_window = 14,
  n_sim = 1000L,
  base_seed = 9786L
)

# Size of the final plot
plot_size <- c(width = 11, height = 11.7)

# If we plot only half of the scenarios, by how much do we divide the height?
plot_halving_coeff <- 1.75

# Scenario blocks for the static branching:
# NegBin-L
# NegBin-Q
# Poisson
outer_scenarios <- data.frame(
  distribution = c("NegBin-L", "NegBin-Q", "Poiss", "Branching")
) |>
  dplyr::mutate(
    scenario_id = distribution
  )


# Define pipeline ==============================================================
list(
  # Static branching over 4 scenario blocks defined in `outer_scenarios`
  tar_map(
    unlist = FALSE,
    values = outer_scenarios,
    names = scenario_id,
    # Get the data frame with scenario parameters
    tar_target(scenarios, create_scenario_grid(distribution)),
    # Serial interval, depends on the scenario
    tar_target(
      si,
      with(
        global_params,
        discr_si(
          seq_len(scenarios$n_burnin),
          scenarios$mean_si,
          scenarios$std_si
        )
      ),
      pattern = map(scenarios)
    ),
    # Initial values. Sample iid counts using a reasonable data generating
    # mechanism.
    tar_target(
      init,
      initialize_trajectory(
        scenarios$n_burnin,
        scenarios$init_magnitude,
        scenarios$init_sd,
        seed = global_params$base_seed + scenarios$init_seed,
        model = distribution
      ),
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
          n_burnin = scenarios$n_burnin,
          init = init,
          R_eff = scenarios$R_eff,
          si = si,
          lgt = long_window,
          model = distribution,
          nb_size = scenarios$nb_size,
          offspring_disp = scenarios$offspring_disp,
          reporting_prob = scenarios$reporting_prob,
          seed = global_params$base_seed + scenarios$scenario_number
        )
      ),
      pattern = map(init, scenarios, si),
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
    # Replace the divergent NegBin estimates by the Poisson ones
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
          # Length of the initialization. We don't plot have fixed initial
          # values to plot for the branching process. For the renewal equation
          # models, the initialization is always as long as the burn-in period
          # for the sake of simplicity.
          if (distribution == "Branching") 0 else scenarios$n_burnin,
          # A longer burn-in period for the branching process to take off
          if (distribution == "Branching") {
            2 * scenarios$n_burnin
          } else {
            scenarios$n_burnin
          }
        ),
        coverage = plot_coverage(
          scenarios$R_eff,
          scenarios$nb_size,
          df_R_hat,
          seq(0, 1, by = 0.01),
          distribution,
          model_colors
        ),
        meta = plot_metadata(
          scenarios$R_eff,
          scenarios$nb_size,
          scenarios$magnitude,
          scenarios$mean_si,
          scenarios$std_si,
          distribution,
          scenarios$offspring_disp
        )
      ),
      pattern = map(trajectories, df_R_hat, scenarios),
      iteration = "list"
    ),
    # How many did converge?
    tar_target(
      convergence_tables,
      summarize_convergence(
        df_R_hat,
        scenarios$magnitude,
        scenarios$nb_size,
        scenarios$R_eff,
        scenarios$serial_interval,
        # Model name from the outer scenarios data frame
        distribution
      ),
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
        scenarios$nb_size,
        scenarios$magnitude,
        scenarios$serial_interval,
        model_colors,
        distribution
      ),
      pattern = map(df_R_hat, scenarios),
      iteration = "list"
    ),
    # Save plots
    tar_target(saved_figures, {
      compose_and_save_plots(
        distribution,
        plot_panels,
        plot_density_panels,
        unlist(global_params[c("short_window", "long_window")]),
        plot_size,
        scenarios[, "scenario_id"],
        plot_halving_coeff
      )
    }),
    tar_target(
      saved_figures_poster,
      {
        # Create the coverage plot for a subset of scenarios. We will use only
        # NegBin-L and NegBin-Q versions.
        # CURRENTLY, THIS PLOT NEEDS ADJUSTMENTS REGARDING THE TEXT SIZES.
        # NEEDS TO BE RESOLVED BEFORE GENERATING THE PLOTS FOR THE POSTER.
        generate_poster_figure <- distribution == "NegBin-Q" ||
          (distribution == "NegBin-L")
        if (generate_poster_figure) {
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
            width = plot_size["width"] / 1.4,
            height = (plot_size["height"] / 1.4) / 2.35
          )
        }
      }
    )
  )
)
