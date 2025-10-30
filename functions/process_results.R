#' Extract the $\hat{R}_t$ estimate
#'
#' @description This function extracts the estimate of $\hat{R}$ and its
#' standard errors from a fitted model.
#'
#' @param fitted_model a fitted \code{gamlss}, or \code{glm} model with one
#'   covariate and no intercept
#' @return a named list with three elements
#'   \describe{
#'      \item{\code{R_hat}}{point estimate of $R_t$}
#'      \item{\code{se_hat}}{standard error of the estimate}
#'      \item{\code{converged}}{logical flag indicating the convergence of the
#'      fitting algorithm}
#'   }
extract_ests <- function(fitted_model) {
  # Initialize the return list
  res <- list(R_hat = NA, se_hat = NA, converged = FALSE)
  if (!inherits(fitted_model, "try-error")) {
    res$converged <- fitted_model$converged
    # The fitted model objects are different for `glm()` and `gamlss()`
    if (class(fitted_model)[1] == "gamlss") {
      res$R_hat <- unlist(fitted_model$mu.coefficients[1])
      # The covariance matrix might not be available, especially for lower
      # sample sizes
      res$se_hat <- sqrt(
        tryCatch(vcov(fitted_model), error = function(e) matrix(NA))[1, 1]
      )
    } else if (class(fitted_model)[1] == "glm") {
      res$R_hat <- unlist(fitted_model$coefficients[1])
      res$se_hat <- unlist(summary(fitted_model)$coefficients[1, "Std. Error"])
    }
  }
  res
}

bind_ests_to_df <- function(results) {
  R_hat <- lapply(results, map, "R_hat") |>
    lapply(unlist) |>
    lapply(unname)
  se_hat <- lapply(results, map, "se_hat") |> lapply(unlist) |> lapply(unname)
  converged <- lapply(results, map, "converged") |>
    lapply(unlist) |>
    lapply(unname)
  # Arrange into a data frame
  df_R_hat <- data.frame(
    R = unlist(R_hat),
    se = unlist(se_hat),
    converged = unlist(converged),
    model = rep(names(results), each = length(R_hat[[1]]))
  )
  df_R_hat
}

#' Replace the divergent estimates
#'
#' @description First, marks unstable estimates as NA across all models
#'   (extreme R/se or non-finite). Then, for NegBin-L/Q rows with
#'   \code{converged == FALSE}, replaces \code{R} and \code{se} with the
#'   corresponding Poisson estimates from the same iteration and window.
#'
#' @param df_R_hat_raw a data frame with raw R estimates containing columns
#'   \code{R}, \code{se}, \code{converged} and \code{model}
#' @return a data frame with the same columns as \code{df_R_hat}, but with some
#'   values replaced by the Poisson estimates.
replace_divergent <- function(df_R_hat_raw) {
  # Replace the divergent and other unstable estimates by NA
  rows_to_keep <- with(df_R_hat_raw, converged & R < 10 & se < 14 & !is.nan(se))
  df_R_hat <- df_R_hat_raw |> mutate(
    R = ifelse(rows_to_keep, R, NA),
    se = ifelse(rows_to_keep, se, NA)
  )

  # Which rows to replace by the Poisson estimates. Only negative binomial ones
  # are affected.
  rows_to_replace <- !df_R_hat$converged &
    df_R_hat$model %in% c("NegBin-L", "NegBin-Q")
  # Create a dummy data frame by repeating the Poisson part as many times as we
  # have models. This assures, that the replacement values are located on
  # correct positions.
  df_R_hat_pois_copy <- df_R_hat |>
    dplyr::filter(model == "Poiss") |>
    replicate(n = length(unique(df_R_hat$model)), simplify = FALSE) |>
    dplyr::bind_rows()

  df_R_hat[rows_to_replace, c("R", "se")] <-
    df_R_hat_pois_copy[rows_to_replace, c("R", "se")]
  df_R_hat
}

#' Calculate the coverage of the CIs
#'
#' @description This function calculates the coverage of the true R value by the
#'   confidence intervals based on multiple simulation runs.
#'
#' @param est a vector of the R point estimates
#' @param se a vector of the estimated standard errors
#' @param true_par a numerical value of the true parameter used to generate the
#'   simulated trajectory
#' @param level a real value between 0 and 1, the coverage level
#' @return the empirical coverage between 0 and 1
calc_coverage <- function(est, se, true_par, level) {
  lower <- est - qnorm(1 - (1 - level) / 2) * se
  upper <- est + qnorm(1 - (1 - level) / 2) * se
  coverage <- (lower < true_par & upper > true_par)
  mean(coverage, na.rm = TRUE)
}

#' Summarize number of convergent fittings in a contingency table
#'
#' @param df_R_hat a data frame with raw R estimates containing columns
#'   \code{R}, \code{se}, \code{converged} \code{window_len_fct} and
#'   \code{model}
#' @param magnitude string, magnitude of the simulation scenario
#' @param nb_size numeric, dispersion parameter of the negative binomial
#'   distribution used in the simulation scenario
#' @param R_eff numeric, value of the effective reproductive number
#'   distribution used in the simulation scenario
#' @true_model string, count distribution used in the simulation scenario
#' @return a list with two data frames:
#'   \itemize{
#'     \item \code{df_convergence}: counts where \code{converged == TRUE} and
#'       R estimate is not \code{NA}, wide by \code{model}, plus
#'       metadata about the parameter values and window length
#'     \item \code{df_unstable}: counts where \code{converged == TRUE} but
#'       R estimate is \code{NA} after post‑processing, wide by \code{model},
#'       plus metadata about the parameter values and window length
#'   }
summarize_convergence <- function(
  df_R_hat,
  magnitude,
  nb_size,
  R_eff,
  true_model
) {
  df_summarized <- df_R_hat |>
    mutate(
      R_eff = R_eff,
      overdispersion = if (true_model == "NegBin-Q") {
        round(1 / nb_size, digits = 2)
      } else if (true_model == "NegBin-L") {
        round(1 + 1 / nb_size, digits = 2)
      } else {
        NA
      },
      magnitude = magnitude
    ) |>
    group_by(R_eff, overdispersion, magnitude, window_len_fct, model) |>
    summarise(
      converged = sum(converged & !is.na(R), na.rm = TRUE),
      # unstable = convergent but masked/flagged by NA in `R` and `se`columns
      unstable = sum(converged & is.na(R), na.rm = TRUE)
    ) |>
    ungroup()
  # Create the table counting convergent runs
  df_convergence <- df_summarized |>
    dplyr::select(-unstable) |>
    pivot_wider(names_from = "model", values_from = "converged")
  # Create the table counting unstable, but convergent runs
  df_unstable <- df_summarized |>
    dplyr::select(-converged) |>
    pivot_wider(names_from = "model", values_from = "unstable")
  list(df_convergence = df_convergence, df_unstable = df_unstable)
}
