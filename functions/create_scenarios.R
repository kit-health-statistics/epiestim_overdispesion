#' Create the scenario grid for the simulation
#'
#' @description Creates the scenario grid of the simulation with
#'   \itemize{
#'     \item 2 orders of magnitude (5 and 100)
#'     \item 2 degrees of dispersion (low and high), does not apply for the
#'     Poisson distribution
#'     \item 2 true values of the effective reproduction number R
#'     \item 1 serial interval distribution
#'   }
#'  The final number of simulation scenarios is 8 for "NegBin-L", or "NegBin-Q"
#'  and 4 for "Poiss"
#' @param distribution the count distribution used in the 8 scenarios (or 4 for
#'   Poisson). Possible values: "Poiss", "NegBin-L", or "NegBin-Q"
#' @return a data frame with scenario names and parameter values
create_scenario_grid <- function(
  distribution = c("NegBin-L", "NegBin-Q", "Poiss")
) {

  distribution <- match.arg(distribution)

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
      c(20, 1 / 0.9)
    } else {
      # Set the dispersion parameter to NA for Poisson
      c(NA, NA)
    }
  )
  # Pair the magnitude degree name, value and initialization parameters
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
    dplyr::left_join(magnitude, by = "magnitude")

  # If the data generating process is Poisson, we don't have scenarios for
  # different dispersion values. Values of the dispersion parameter shall be
  # all NA, thus the rows will be dropped as duplicates.
  if (distribution == "Poiss") {
    scenarios <- scenarios |>
      dplyr::select(-dispersion) |>
      dplyr::distinct(.keep_all = FALSE) |>
      # Add the string denoting the dispersion back, even though it's not
      # technically needed.
      dplyr::mutate(dispersion = "not_applicable")
  }
  # Add a scenario number and ID
  scenarios |> dplyr::mutate(
    scenario_number = seq_len(dplyr::n()),
    scenario_id = paste(
      "sc",
      stringr::str_pad(init_magnitude, 3, pad = "0"),
      dispersion,
      paste0("R", R_eff),
      sep = "_"
    )
  )
}
