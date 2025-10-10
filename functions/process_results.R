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

#' Replace the divergent estimates by NAs
#'
#' @description This function goes through the data frame of raw R estimates
#'   and replaces the divergent fits by NAs.
#'
#' @param df_R_hat a data frame with raw R estimates containing columns
#'   \code{R}, \code{se}, \code{converged} and \code{model}
#' @return a data frame with the same columns as \code{df_R_hat}, but with some
#'   values replaced by NAs.
remove_divergent <- function(df_R_hat) {
  df_R_hat |> mutate(
    R = ifelse(converged, R, NA),
    se = ifelse(converged, se, NA)
  )
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
