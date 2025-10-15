#' Fit the models and return the results in a data frame
#'
#' @description This function is the outermost wrapper of the model fitting
#'   procedure. It fits all the models for both window lengths
#'   (can be generalized to more windows) and binds the results together in a
#'   data frame. The short estimation window is a subset of the long one.
#' @param X a matrix of the simulated counts, one column per simulation run
#' @param Lambda a matrix of the single covariate, one column per simulation run
#' @param short_window the length of the shorter estimation window
#' @param long_window the length of the longer estimation window
#' @return a data frame with six columns:
#'   \begin{itemize}
#'     \item \code{R}, estimates of the effective reproduction number
#'     \item \code{se}, estimates of the standard errors
#'     \item \code{converged}, an indicator, whether the fitting algorithm
#'     converged
#'     \item \code{model}, count distribution label
#'     \item \code{window_len}, length of the estimation window
#'     \item \code{window_len_fct}, length of the estimation window as a string,
#'     either "short", or "long"
#'   \end{itemize}
create_results_df <- function(X, Lambda, short_window, long_window) {
  # Fitting for the short window
  df_R_hat_raw_short <- fit_all_models(X, Lambda, short_window) |>
    mutate(window_len_fct = "short")
  # Fitting for the long window
  df_R_hat_raw_long <- fit_all_models(X, Lambda, long_window) |>
    mutate(window_len_fct = "long")

  # Bind together and make the window length into a factor
  df_R_hat_raw <- rbind(df_R_hat_raw_short, df_R_hat_raw_long) |>
    mutate(window_len_fct = factor(window_len_fct))

  # Remove unstable estimates
  remove_divergent(df_R_hat_raw)
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
#'   runs by NAs, containing columns \code{R}, \code{se}, \code{converged} and
#'   \code{model}
#' @param X a matrix of incidences, can be dropped when we get rid of the normal
#'   approximation part
#' @param Lambda matrix of the single covariate, can be dropped when we get rid
#'   of the normal approximation part
#' @param nominal_covr a vector of the nominal coverage levels
#' @return a list containing two data frames with common columns
#'   \describe{
#'     \item{\code{covr_nominal}}{the nominal coverage level}
#'     \item{\code{covr_empirical}}{the empirical coverage level from the
#'       simulation}
#'     \item \code{window_len}{the length of the estimation window}
#'     \item \code{window_len_fct}{the name of the length of the estimation
#'       window as a factor}
#'     \item \code{type}{string indicating types of coverage, namely 
#'       "Empirical coverage", “theoretical Poisson coverage” under model
#'       misspecification and the coverage using the OLS estimates
#'       ("normal approx."). The normal approximation is only in the 
#'       \code{points} data frame and will be dropped for the final version of
#'       the paper.}
#'   }
#'   The data frame \code{lines} is intended for plotting the coverage as lines
#'   and it contains additionally one column
#'   \describe{
#'     \item \code{model}{the corresponding distributional model}
#'   }
create_coverage_df <- function(
  R_eff,
  nb_size,
  df_R_hat,
  X,
  Lambda,
  nominal_covr
) {
  # Calculate the coverage of true value of R. We take the point estimates and
  # the SEs from the fitted model.
  df_coverage_model <- df_R_hat |>
    group_by(model, window_len_fct, window_len) |>
    reframe(
      covr_nominal = nominal_covr,
      covr_empirical = sapply(
        nominal_covr,
        calc_coverage,
        est = R,
        se = se,
        true_par = R_eff
      )
    ) |>
    mutate(
      type = "Empirical coverage"
    )

  # What would be the coverage, if we assume the data to be normally distributed
  # (i.e. we approximate the NB by a Gaussian)? Same for all models.
  # For simplicity we assume that X / Lambda = R + iid Gaussian error,
  # so that the sample mean and variance give us the point estimate and the
  # standard error directly. This is very crude, but can give us an idea,
  # whether the counts already converge to a Gaussian, or not.
  # Will be removed for the final plot.
  df_coverage_norm_approx <- df_coverage_model |>
    filter(model == "Poiss") |>
    select(-model) |>
    group_by(window_len_fct, window_len) |>
    mutate(
      type = "normal approx.",
      covr_empirical = sapply(
        covr_nominal,
        calc_coverage,
        est = apply(
          X[seq_len(window_len[1]), ] / Lambda[seq_len(window_len[1]), ],
          2,
          mean,
          na.rm = TRUE
        ),
        se = apply(
          X[seq_len(window_len[1]), ] / Lambda[seq_len(window_len[1]), ],
          2,
          sd,
          na.rm = TRUE
        ) *
          (window_len[1] - 1) /
          # From the unbiased estimator to the MLE
          (window_len[1] * sqrt(window_len[1])),
        true_par = R_eff
      )
    )
  # What is the actual coverage when the Poisson model is misspecified?
  # If it is, the true variance of the R estimate will be the Poisson variance
  # inflated by a factor depending on the dispersion parameter of the NB
  # distribution.
  var_infl_factor_true <- (1 + 1 / nb_size)  # For NegBin-L
  df_coverage_poiss <- df_coverage_norm_approx |>
    mutate(
      covr_empirical = pnorm(
        qnorm(1 - (1 - covr_nominal) / 2) / sqrt(var_infl_factor_true)
      ) -
        pnorm(-qnorm(1 - (1 - covr_nominal) / 2) / sqrt(var_infl_factor_true)),
      type = "theoretical Poisson coverage",
      model = "theoretical Poisson"
    )

  # Return data frames in a list
  list(
    lines = rbind(df_coverage_model, df_coverage_poiss),
    points = df_coverage_norm_approx
  )
}
