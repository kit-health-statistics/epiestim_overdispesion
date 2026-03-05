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
  distribution = c("NegBin-L", "NegBin-Q", "Poiss", "Branching")
) {

  distribution <- match.arg(distribution)

  scenario_grid <- expand.grid(
    magnitude = c("low", "high"),
    dispersion = c("low", "high"),
    R_eff = c(1.5, 2.5),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # Pair the dispersion degree name and values of the corresponding parameters
  if (distribution == "Branching") {
    dispersion <- data.frame(
      dispersion = c("low", "high"),
      nb_size = c(NA, NA),
      R_sd = c(0.1, 0.3),
      reporting_prob = c(0.8, 0.8)
    )
  } else {
    dispersion <- data.frame(
      dispersion = c("low", "high"),
      nb_size = if (distribution == "NegBin-L") {
        c(2, 0.2)
      } else if (distribution == "NegBin-Q") {
        c(1 / 0.02, 1 / 0.06)
      } else {
        # Set the negative binomial dispersion parameter to NA for Poisson
        c(NA, NA)
      },
      # Set the parameters for the branching process to NA
      R_sd = c(NA, NA),
      reporting_prob = c(NA, NA)
    )
  }

  # Pair the magnitude degree name, value and initialization parameters
  if (distribution == "Branching") {
    magnitude <- data.frame(
      magnitude = "low",
      # Initial value
      init_magnitude = 15,
      # No sd and seed, since for the branching process we initialize by a
      # single number - by the number of initial infectious individuals
      init_sd = NA,
      init_seed = NA
    )
  } else {
    magnitude <- data.frame(
      magnitude = c("low", "high"),
      # Parameters for generating the initial values
      init_magnitude = c(5, 100),
      init_sd = c(2, 10),
      init_seed = c(10L, 100L)
    )
  }

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
  # If the data generating process is a branching process, we don't have
  # scenarios for higher magnitudes.
  if (distribution == "Branching") {
    scenarios <- scenarios |>
      filter(magnitude == "low")
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
