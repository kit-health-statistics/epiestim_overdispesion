#' Create the scenario grid for the simulation
#'
#' @description Creates the scenario grid of the simulation with
#'   \itemize{
#'     \item 2 orders of magnitude (around 5 and 100), does not apply for the
#'     Branching process, where we have only lower magnitude
#'     \item 2 degrees of dispersion (low and high), does not apply for the
#'     Poisson distribution
#'     \item 2 true values of the effective reproduction number R +
#'     time-varying R, the time varying R applies only for selected NegBin-L
#'     scenarios
#'     \item 3 serial interval distributions for selected NegBin-L scenarios,
#'     1 serial interval distribution fo the rest
#'   }
#'  The final number of simulation scenarios is 8 for "NegBin-Q", 4 for "Poiss",
#'  4 for a branching process and 20 "NegBin-L"
#' @param distribution the count distribution used in the scenarios.
#'   Possible values: "Poiss", "NegBin-L", "NegBin-Q", or "Branching"
#' @return a data frame with scenario names and parameter values
create_scenario_grid <- function(
  distribution = c("NegBin-L", "NegBin-Q", "Poiss", "Branching")
) {
  distribution <- match.arg(distribution)

  scenario_grid <- expand.grid(
    dispersion = c("low", "high"),
    R_eff = c(1.5, 2.5, "time_dependent"),
    magnitude = c("low", "high"),
    serial_interval = c("RSV", "measles", "influenza"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # Pair the dispersion degree name and values of the corresponding parameters
  if (distribution == "Branching") {
    dispersion <- data.frame(
      dispersion = c("low", "high"),
      nb_size = c(NA, NA),
      offspring_disp = c(1.5, 3),
      reporting_prob = c(0.5, 0.5)
    )
  } else {
    dispersion <- data.frame(
      dispersion = c("low", "high"),
      nb_size = if (distribution == "NegBin-L") {
        c(2, 0.2)
      } else if (distribution == "NegBin-Q") {
        c(50, 16.6)
      } else {
        # Set the negative binomial dispersion parameter to NA for Poisson
        c(NA, NA)
      },
      # Set the parameters for the branching process to NA
      offspring_disp = c(NA, NA),
      R_sd = c(NA, NA),
      reporting_prob = c(NA, NA)
    )
  }

  # Pair the magnitude degree name, value and initialization parameters
  if (distribution == "Branching") {
    magnitude <- data.frame(
      magnitude = "low",
      # Initial value
      init_magnitude = 4,
      init_sd = 2,
      init_seed = 10L
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

  # Pair the disease name and values of the corresponding serial interval
  # parameters
  serial_interval <- data.frame(
    serial_interval = c("RSV", "influenza", "measles"),
    mean_si = c(7.5, 3.7, 13.7),
    std_si = c(2.1, 1.1, 1.5),
    # The burn-in period must be longer for a longer serial interval to keep R
    # constant long enough.
    # Initialization and burn-in are handled differently in the branching
    # process simulation and plotting. For this reason, the length of the
    # burn-in period depends only on the serial interval. The functions
    # corresponding to the branching process will take care of the
    # peculiarities.
    n_burnin = c(14, 14, 21)
  )

  # Join the parameter names and its values
  scenarios <- scenario_grid |>
    dplyr::left_join(dispersion, by = "dispersion") |>
    dplyr::left_join(magnitude, by = "magnitude") |>
    dplyr::left_join(serial_interval, by = "serial_interval")

  # Filter non-applicable scenarios.
  if (distribution == "Poiss") {
    # If the data generating process is Poisson, we don't have scenarios for
    # different dispersion values. Values of the dispersion parameter shall be
    # all NA, thus the rows will be dropped as duplicates.
    scenarios <- scenarios |>
      dplyr::select(-dispersion) |>
      dplyr::distinct(.keep_all = FALSE) |>
      # Add the string denoting the dispersion back, even though it's not
      # technically needed.
      dplyr::mutate(dispersion = "not_applicable") |>
      # We will generate only simulation scenarios using the serial interval of
      # the RSV and constant true values of R.
      filter(serial_interval == "RSV" & R_eff != "time_dependent")
  } else if (distribution == "Branching") {
    # If the data generating process is a branching process, we don't have
    # scenarios for higher magnitudes and we also rewrite the true value of R to
    # be lower in order to make the trajectories less explosive.
    scenarios <- scenarios |>
      dplyr::filter(magnitude == "low") |>
      mutate(R_eff = dplyr::case_when(R_eff == 1.5 ~ 1.2, R_eff == 2.5 ~ 2)) |>
      # We will generate only simulation scenarios using the serial interval of
      # the RSV and constant true values of R.
      filter(serial_interval == "RSV" & R_eff != "time_dependent")
  } else if (distribution == "NegBin-Q") {
    # We will generate only simulation scenarios using the serial interval of
    # the RSV and constant true values of R.
    scenarios <- scenarios |>
      filter(serial_interval == "RSV" & R_eff != "time_dependent")
  } else if (distribution == "NegBin-L") {
    scenarios <- scenarios |>
      # For low magnitude scenarios, which are shown in the main paper, we will
      # generate scenarios with all 3 generation times. For the high magnitude,
      # we will generate only simulation scenarios using the serial interval of
      # the RSV.
      # In addition, we generate trajectories with the RSV serial interval and
      # time dependent true values of the reproductive number. This will be done
      # for both high and low magnitude.
      filter(
        (magnitude == "low" & serial_interval != "RSV" &
           R_eff != "time_dependent") |
          serial_interval == "RSV"
      )
  }
  # Add a scenario number and ID
  scenarios |> dplyr::mutate(
    scenario_number = seq_len(dplyr::n()),
    scenario_id = paste(
      "sc",
      stringr::str_pad(init_magnitude, 3, pad = "0"),
      dispersion,
      paste0("R", R_eff),
      serial_interval,
      sep = "_"
    )
  )
}

#' Generate a vector of true values of R
#'
#' @description Create a vector of the effective reproductive number values,
#'   which are either constant, or in the form of a cosine wave.
#'
#' @param R_eff a value of R, if passed as a string "time_dependent", the
#'   effective reproductive number will be a sine wave. Otherwise the value of
#'   `R_eff` will be coerced into a numeric value and repeated as a constant
#'   vector.
#' @param lgt an integer, the desired length of the vector of true R values
#' @return vector of true R values
get_true_R <- function(R_eff, lgt) {
  if (R_eff == "time_dependent") {
    R_vec <- 0.11 * cos(2.5 * pi * seq_len(lgt) / lgt) + 1.5
  } else {
    R_vec <- rep(as.numeric(R_eff), lgt)
  }
  R_vec
}
