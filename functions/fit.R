#' Fit the renewal equation as a GLM model
#'
#' @param X vector, the incidence
#' @param Lambda vector, the single covariate
#' @param model the count distribution, one of "Poiss", "Q-Poiss", "NegBin-Q"
#'   and "NegBin-L"
#' @return a fitted GLM of class \code{glm} for "Poiss" and "Q-Poiss", or
#'   \code{gamlss} for "NegBin-L" and "NegBin-Q". Returns an object of class
#'   \code{try-error} if fitting fails.
fit_reg_model <- function(
  X,
  Lambda,
  model = c("Poiss", "Q-Poiss", "NegBin-Q", "NegBin-L")
) {
  # Select the family expression based on the model
  family <- switch(
    model,
    Poiss = poisson(link = "identity"),
    `Q-Poiss` = quasipoisson(link = "identity"),
    `NegBin-L` = NBII(mu.link = "identity", sigma.link = "log"),
    `NegBin-Q` = NBI(mu.link = "identity", sigma.link = "log")
  )

  if (model %in% c("Poiss", "Q-Poiss")) {
    # Use standard `glm()` function for Poisson and Quasi-Poisson
    model_call <- substitute(
      glm(X ~ Lambda - 1, family = fam),
      list(fam = family)
    )
  } else {
    # Use standard `gamlss()` function for the negative binomial distributions
    model_call <-
      substitute(
        gamlss(
          formula = X ~ Lambda - 1,
          family = fam,
          control = gamlss.control(trace = FALSE)
        ),
        list(fam = family)
      )
  }
  try(eval(model_call), silent = TRUE)
}

#' Fit all models and extract results
#'
#' @description This function fits all 4 renewal equation models and extracts
#'   the results into a data frame.
#' @param X a matrix of the simulated counts, one column per simulation run
#' @param Lambda a matrix of the single covariate, one column per simulation run
#' @param window the length of the estimation window
#' @return a data frame with three columns:
#'   \begin{itemize}
#'     \item \code{R_hat}, estimates of the effective reproduction number
#'     \item \code{se_hat}, estimates of the standard errors
#'     \item \code{overdisp}, estimates of the overdispersion parameter, on
#'     different scales for NegBin-L and NegBin-Q
#'     \item \code{converged}, an indicator, whether the fitting algorithm
#'     converged
#'     \item \code{window_len}, length of the estimation window
#'   \end{itemize}
fit_all_models <- function(X, Lambda, window) {
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
  fits <- vector("list", 4)
  names(fits) <- c("Poiss", "Q-Poiss", "NegBin-L", "NegBin-Q")
  for (k in seq_along(fits)) {
    fits[[names(fits)[k]]] <- sapply(
      seq_len(n_sim),
      pre_vectorized_fitting,
      model = names(fits)[k],
      simplify = FALSE
    )
  }

  # Extract results
  results <- lapply(
    fits,
    function(x) lapply(x, extract_ests)
  )
  df_R_hat_raw <- bind_ests_to_df(results) |>
    mutate(window_len = window)
  df_R_hat_raw
}
