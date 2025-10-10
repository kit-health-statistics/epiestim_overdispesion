#' Fit all models and extract results
#'
#' @description This function fits all 4 renewal equation models and extracts
#'   the results into a data frame.
#' @param X a matrix of the simulated counts, one column per simulation run
#' @param Lambda a matrix of the single covariate, one column per simulation run
#' @param window the length of the estimation window
#' @return a data frame with three columns:
#'   \begin{itemize}
#'     \item \code{R}, estimates of the effective reproduction number
#'     \item \code{se}, estimates of the standard errors
#'     \item \code{converged}, an indicator, whether the fitting algorithm
#'     \item \code{model}, count distribution label
#'     converged
#'   \end{itemize}
create_results_df <- function(X, Lambda, window) {
  pre_vectorized_fitting <- function(ind, model) {
    fit_reg_model(
      X = tail(X[, ind], window),
      Lambda = tail(Lambda[, ind], window),
      model = model
    )
  }
  # Retrieve the number of simulation runs
  n_sim <- ncol(X)

  # Fit all models
  fits <- list(
    Poiss = sapply(
      seq_len(n_sim),
      pre_vectorized_fitting,
      model = "Poiss",
      simplify = FALSE
    ),
    `Q-Poiss` = sapply(
      seq_len(n_sim),
      pre_vectorized_fitting,
      model = "Q-Poiss",
      simplify = FALSE
    ),
    `NegBin-L` = sapply(
      seq_len(n_sim),
      pre_vectorized_fitting,
      model = "NegBin-L",
      simplify = FALSE
    ),
    `NegBin-Q` = sapply(
      seq_len(n_sim),
      pre_vectorized_fitting,
      model = "NegBin-Q",
      simplify = FALSE
    )
  )

  # Extract results
  results <- lapply(fits, function(x) lapply(x, extract_ests))
  df_R_hat_raw <- bind_ests_to_df(results)
  df_R_hat_raw
}

#' Calculate the coverage
#'
#' @description This function calculates empirical coverage at given nominal
#'   levels for each model including the normal approximation coverage and the
#'   "true" coverage of the Poisson model under dispersion misspecification.
#' @param R_eff the true value of the effective reproduction number used to
#'   generate the trajectory
#' @param nb_size positive real value, the size parameter of the negative
#'   binomial distribution for "NegBin-Q" and "NegBin-L"
#' @param df_R_hat a data frame with R estimates, after replacing the divergent
#'   runs by NAs, containing columns
#'   \code{R_hat}, \code{se_hat} and \code{converged}
#' @param X a matrix of incidences, can be dropped when we get rid of the normal
#'   approximation part
#' @param Lambda matrix of the single covariate, can be dropped when we get rid
#'   of the normal approximation part
#' @param window the length of the estimation window, can be dropped when we get
#'   rid of the normal approximation part
#' @param nominal_covr a vector of the nominal coverage levels
#' @return a data frame with four columns:
#'   \describe{
#'     \item{\code{covr_nominal}}{the nominal coverage level}
#'     \item{\code{covr_empirical}}{the empirical coverage level from the
#'       simulation}
#'     \item{\code{model}}{the corresponding fitted model}
#'     \item{\code{type}}{string indicating using which quantities the coverage
#'       was calculated. Currently, this serves as an aesthetic for the ggplot
#'       colors in the coverage plot.}
#'   }
create_coverage_df <- function(
  R_eff,
  nb_size,
  df_R_hat,
  X,
  Lambda,
  window,
  nominal_covr
) {
  # Calculate the coverage of true value of R. We take the point estimates and
  # the SEs from the fitted model.
  df_coverage_model <- df_R_hat |>
    group_by(model) |>
    reframe(
      type = "GLM",
      covr_nominal = nominal_covr,
      covr_empirical = sapply(
        nominal_covr,
        calc_coverage,
        est = R,
        se = se,
        true_par = R_eff
      )
    )

  # What would be the coverage, if we assume the data to be normally distributed
  # (i.e. we approximate the NB by a Gaussian)? Same for all models.
  # For simplicity we assume that X / Lambda = R + iid Gaussian error,
  # so that the sample mean and variance give us the point estimate and the
  # standard error directly. This is very crude, but can give us an idea,
  # whether the counts already converge to a Gaussian, or not.
  # Will be removed for the final plot.
  df_coverage_norm_approx <- df_coverage_model |>
    mutate(
      type = "normal approx.",
      covr_empirical = sapply(
        covr_nominal,
        calc_coverage,
        est = apply(X / Lambda, 2, mean, na.rm = TRUE),
        se = apply(X / Lambda, 2, sd, na.rm = TRUE) *
          (window - 1) /
          (window * sqrt(window)), # From the unbiased estimator to the MLE
        true_par = R_eff
      )
    )

  # What is the actual coverage when the Poisson model is misspecified?
  # If it is, the true variance of the R estimate will be the Poisson variance
  # inflated by a factor depending on the dispersion parameter of the NB
  # distribution.
  var_infl_factor_true <- (1 + 1 / nb_size)  # For NegBin-L
  df_coverage_poiss <- data.frame(
    covr_nominal = nominal_covr,
    model = "Poiss",
    type = "'real'"
  ) |>
    mutate(
      covr_empirical = pnorm(
        qnorm(1 - (1 - covr_nominal) / 2) / sqrt(var_infl_factor_true)
      ) -
        pnorm(-qnorm(1 - (1 - covr_nominal) / 2) / sqrt(var_infl_factor_true))
    )

  # Bind together
  df_coverage <- rbind(
    df_coverage_model,
    df_coverage_norm_approx,
    df_coverage_poiss
  )
  df_coverage
}
