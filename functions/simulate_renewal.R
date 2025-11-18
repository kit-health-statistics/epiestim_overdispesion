initialize_trajectory <- function(
  n_init,
  init_magnitude,
  init_sd,
  weekday_effect,
  seed
) {
  # Set seed ensuring identical initialization for scenarios with the same
  # magnitude. Creates redundant copies (nrow(scenarios) instead of 2) but
  # simplifies pipeline structure with negligible computational cost.
  set.seed(seed)
  pmax(0, round(rnorm(n_init, init_magnitude * weekday_effect, init_sd)))
}

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
#' @param weekday_effect a vector, the potential weekday effects. For scenarios
#'   with a weekday effect, the mean of the vector must be 1. For scenarios
#'   without the effect, all its elements shall equal to 1.  The choice of the
#'   default to be of length 2 is arbitrary and has no effect compared to other
#'   vector lengths. If the length of the initialization is not a multiple of
#'   the number of weekday effects, we align the beginnings of both vectors and
#'   recycle the weekday effects. For example, if the initialization has length
#'   8 and the weekday effect 7, the first incidence value after the
#'   initialization will be sampled with the 2nd element of the weekday effect.
#' @return a named list with two elements
#'   \describe{
#'     \item{\code{X}}{integer vector, the simulated incidence}
#'     \item{\code{Lambda}}{numeric vector, value of the sole covariate in the
#'     renewal equation}
#'  }
simulate_renewal <- function(
  init,
  R,
  si,
  lgt,
  model = c("Poiss", "NegBin-Q", "NegBin-L"),
  nb_size = NULL,
  weekday_effect = c(1, 1)
) {
  model <- match.arg(model)
  X <- Lambda <- numeric(lgt)
  n_si <- length(si)
  X[seq_along(init)] <- init
  Lambda[seq_along(init)] <- NA
  len_weekday_eff <- length(weekday_effect)

  for (t in (length(init) + 1):lgt) {
    Lambda[t] <- sum(si * X[t - (1:n_si)]) *
      # Multiply by the weekday effect, `weekday_effect` is all ones, when none
      # is present.
      weekday_effect[(t - 1) %% len_weekday_eff + 1]
    X[t] <- generate_from_obs_mod(
      R * Lambda[t],
      model = model,
      nb_size = nb_size
    )
  }
  list(X = X, Lambda = Lambda)
}

#' Generate multiple incidence trajectories
#'
#' @description This function simulates multiple incidence trajectories based on
#'  the renewal equation.
#'
#' @param n_sim an integer, number of trajectories to generate
#' @param init a vector, incidence at the beginning of the incidence trajectory
#' @param R_eff positive real value, the value of the effective reproduction
#'   number, same across the whole trajectory
#' @param si a non-negative vector, the distribution of the serial interval,
#'   must not be longer than the \code{init} vector
#' @param lgt positive integer, total length of one trajectory. All trajectories
#'   are of equal length
#' @param model the count distribution, one of "Poiss", "NegBin-Q" and
#'   "NegBin-L"
#' @param nb_size positive real value, the size parameter of the negative
#'   binomial distribution for "NegBin-Q" and "NegBin-L"
#' @param seed an integer, seed to be set before the `n_sim` replications of the
#'   trajectory simulation
#' @param weekday_effect a vector, the potential weekday effects. For scenarios
#'   with a weekday effect, the mean of the vector must be 1. For scenarios
#'   without the effect, all its elements shall equal to 1. The choice of the
#'   default to be of length 2 is arbitrary and has no effect compared to other
#'   vector lengths.
#' @return a named list with two elements
#'   \describe{
#'     \item{\code{X}}{integer matrix, the simulated incidence, each column
#'     represents a single trajectory}
#'     \item{\code{Lambda}}{numeric matrix, value of the sole covariate in the
#'     renewal equation, each column represents a single trajectory}
#'  }
generate_trajectories <- function(
  n_sim,
  init,
  R_eff,
  si,
  lgt,
  model,
  nb_size = NULL,
  seed = 432,
  weekday_effect = c(1, 1)
) {
  set.seed(seed)
  trajectories <- replicate(
    n_sim,
    simulate_renewal(init, R_eff, si, lgt, model, nb_size, weekday_effect),
    simplify = FALSE
  )
  X <- do.call(cbind, map(trajectories, "X"))
  Lambda <- do.call(cbind, map(trajectories, "Lambda"))
  list(X = X, Lambda = Lambda)
}
