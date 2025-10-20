#' Create the scenario grid for the simulation
#'
#' @description Creates the scenario grid of the simulation with
#'   \begin{itemize}
#'     \item 2 orders of magnitude (5 and 100)
#'     \item 2 degrees of dispersion (low and high)
#'     \item 2 true values of the effective reproduction number R
#'     \item 1 serial interval distribution
#'  \end{itemize}
#'  The final number of simulation scenarios is 8.
#' @return a data frame with scenario names and parameter values
create_scenario_grid <- function(distribution) {
  scenario_grid <- expand.grid(
    magnitude = c("low", "high"),
    dispersion = c("low_disp", "high_disp"),
    R_eff = c(1.5, 2.5),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # Pair the dispersion degree name and value
  dispersion <- data.frame(
    dispersion = c("low_disp", "high_disp"),
    nb_size = if (distribution == "NegBin-L") {
      c(2, 0.2)
    } else if (distribution == "NegBin-Q") {
      c(4, 0.4)
    }
  )

  magnitude <- data.frame(
    magnitude = c("low", "high"),
    # Parameters for generating the initial values
    init_magnitude = c(5, 100),
    init_sd = c(2, 10),
    init_seed = c(10L, 100L)
  )

  # Join the parameter names and its values
  scenarios <- scenario_grid |>
    dplyr::left_join(dispersion, by = "dispersion") |>
    dplyr::left_join(magnitude, by = "magnitude") |>
    dplyr::mutate(
      scenario_number = seq_len(dplyr::n()),
      scenario_id = paste(
        "sc",
        stringr::str_pad(init_magnitude, 3, pad = "0"),
        dispersion,
        paste0("R", R_eff),
        sep = "_"
      )
    )
  scenarios
}
