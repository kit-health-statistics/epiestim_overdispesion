#' Generate the incidence counts from the corresponding distribution
#'
#' @param mu mean value of the count distribution
#' @param model the count distribution, one of "Poiss", "NegBin-Q" and
#'   "NegBin-L"
#' @param nb_size the size parameter of the negative binomial distribution
#' @return integer value of the generated incidence
generate_from_obs_mod <- function(
  mu,
  model = c("Poiss", "NegBin-Q", "NegBin-L"),
  nb_size = NULL
) {
  model <- match.arg(model)
  switch(
    model,
    Poiss = rpois(1, mu),
    `NegBin-L` = rnbinom(1, mu = mu, size = nb_size * mu),
    `NegBin-Q` = rnbinom(1, mu = mu, size = nb_size)
  )
}

#' Simulate incidence using the renewal equation
#'
#' @description This function simulates an incidence trajectory based on the
#' renewal equation:
#'
#' \deqn{\mu_t = R\times \Lambda_t, \quad \Lambda_t = \sum_{d = 1}^{D} \omega_d X_{t - d}}  # nolint
#'
#' $X_t$ has mean value $\mu_t$ and follows, Poisson, NegBin-L, or NegBin-Q
#' distribution.
#'
#' @param init a vector, incidence at the beginning of the incidence trajectory
#' @param R positive real value, the value of the effective reproduction number,
#'   same across the whole trajectory
#' @param si a non-negative vector, the distribution of the serial interval,
#'   must not be longer than the \code{init} vector
#' @param lgt positive integer, total length of the trajectory
#' @param model the count distribution, one of "Poiss", "NegBin-Q" and
#'   "NegBin-L"
#' @param nb_size positive real value, the size parameter of the negative
#'   binomial distribution for "NegBin-Q" and "NegBin-L"
#' @return a named list with two elements
#'  \begin{itemize}
#'     \item \code{X}, the simulated incidence,
#'     \item \code{Lambda} value of the sole covariate in the renewal equation.
#'  \end{itemize}
simulate_renewal <- function(
  init,
  R,
  si,
  lgt,
  model = c("Poiss", "NegBin-Q", "NegBin-L"),
  nb_size = NULL
) {
  model <- match.arg(model)
  X <- Lambda <- numeric(lgt)
  n_si <- length(si)
  X[seq_along(init)] <- init
  Lambda[seq_along(init)] <- NA

  for (t in (length(init) + 1):lgt) {
    Lambda[t] <- sum(si * X[t - (1:n_si)])
    X[t] <- generate_from_obs_mod(
      R * Lambda[t],
      model = model,
      nb_size = nb_size
    )
  }
  return(list(X = X, Lambda = Lambda))
}
