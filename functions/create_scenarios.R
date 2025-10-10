#' Create the scenario grid for the simulation
#'
#' @description Creates the scenario grid of the simulation with
#'  \begin{itemize}
#'    \item 2 orders of magnitude (10 and 500)
#'    \item 2 degrees of dispersion (low and high)
#'    \item 2 estimation window sizes to assess the small sample biases
#'    \item 2 true values of the effective reproduction number R
#'    \item 1 distribution (NegBin-L)
#'    \item 1 serial interval distribution
#'  \end{itemize}
#' @return a data frame with scenario names and parameter values
create_scenario_grid <- function() {
  scenario_grid <- expand.grid(
    distribution = "NegBin-L",
    init_magnitude = c(10, 500),
    dispersion = c("low_disp", "high_disp"),
    window = c("week", "fortnight"),
    R_eff = c(1.5, 2.5),
    mean_si = 6,
    sd_si = 1,
    KEEP.OUT.ATTRS = FALSE
  )

  # Pair the window length name and value
  window <- data.frame(
    window = c("week", "fortnight"),
    window_len = c(7, 14)
  )
  # Pair the dispersion degree name and value
  dispersion <- data.frame(
    dispersion = c("low_disp", "high_disp"),
    nb_size = c(2, 0.2)
  )

  # Join the parameter names and its values
  scenarios <- dplyr::left_join(scenario_grid, window, by = "window") |>
    dplyr::left_join(dispersion, by = "dispersion") |>
    dplyr::mutate(
      scenario_id = paste(
        distribution,
        window,
        init_magnitude,
        dispersion,
        paste0("R", R_eff),
        sep = "_"
      )
    )
  scenarios
}
