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
  base_seed = 9786L
)

# Define pipeline
list(
  tar_target(scenarios, create_scenario_grid()),

  # Serial interval, common for all scenarios
  tar_target(
    si,
    with(global_params, discr_si(seq_len(n_init), mean_si, std_si))
  ),

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
        model = scenarios$distribution,
        nb_size = scenarios$nb_size,
        seed = global_params$base_seed + scenarios$scenario_number
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
    c(
      trajectories = plot_trajectories(
        trajectories$X,
        global_params$short_window,
        global_params$n_init
      ),
      plot_coverage(
        scenarios$R_eff,
        scenarios$nb_size,
        df_R_hat,
        trajectories$X[-seq_len(global_params$n_init), ],
        trajectories$Lambda[-seq_len(global_params$n_init), ],
        seq(0, 1, by = 0.05),
        model_colors
      ),
      meta = plot_metadata(
        scenarios$R_eff,
        scenarios$nb_size,
        scenarios$magnitude
      )
    ),
    pattern = map(trajectories, df_R_hat, scenarios),
    iteration = "list"
  ),

  # Save plots
  tar_target(saved_figures, {
    # Ensure directory exists
    dir.create("figure", showWarnings = FALSE, recursive = TRUE)
    p_simulation <- compose_patches(
      plot_panels,
      global_params$short_window,
      global_params$long_window
    )
    # Save as a PDF
    ggsave(
      here("figure", "simulation_coverage.pdf"),
      p_simulation,
      width = 12.5,
      height = 14
    )
    # Save as a PNG to make the comparison in a PR on GitHub easier.
    # Can be deleted later.
    ggsave(
      here("figure", "simulation_coverage.png"),
      p_simulation,
      width = 12.5,
      height = 14,
      dpi = 400
    )
  })
)
