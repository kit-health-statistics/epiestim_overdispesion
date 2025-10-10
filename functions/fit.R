#' Fit the renewal equation as a GLM model
#'
#' @param X vector, the incidence
#' @param Lambda vector, the single covariate
#' @param model the count distribution, one of "Poiss", "Q-Poiss", "NegBin-Q" and
#'   "NegBin-L"
#' @return a fitted GLM of class \code{glm} for "Poiss" and "Q-Poiss", or
#'   \code{gamlss} for "NegBin-L" and "NegBin-Q"
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
      glm(X ~ Lambda - 1, family = fam, data = df),
      list(fam = family, df = data.frame(X = X, Lambda = Lambda))
    )
  } else {
    # Use standard `gamlss()` function for the negative binomial distributions
    model_call <-
      substitute(
        gamlss(
          formula = X ~ Lambda - 1,
          family = fam,
          data = df,
          control = gamlss.control(trace = FALSE)
        ),
        list(fam = family, df = data.frame(X = X, Lambda = Lambda))
      )
  }
  try(eval(model_call), silent = TRUE)
}
