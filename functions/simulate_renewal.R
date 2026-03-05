initialize_trajectory <- function(
  n_init,
  init_magnitude,
  init_sd,
  weekday_effect,
  seed,
  model = c("Poiss", "NegBin-Q", "NegBin-L", "Branching")
) {
  model <- match.arg(model)
  # Set seed ensuring identical initialization for scenarios with the same
  # magnitude. Creates redundant copies (nrow(scenarios) instead of 2) but
  # simplifies pipeline structure with negligible computational cost. We
  # generate some initial values for the renewal equation models only. For the
  # branching process, we need a single initial value, which is the value of
  # `init_magnitude`.
  if (model == "Branching") {
    init_magnitude
  } else {
    set.seed(seed)
    pmax(0, round(rnorm(n_init, init_magnitude * weekday_effect, init_sd)))
  }
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
  X[seq_along(init)] <- init
  Lambda[seq_along(init)] <- NA
  len_weekday_eff <- length(weekday_effect)

  for (t in (length(init) + 1):(lgt)) {
    Lambda[t] <- sum(si * X[t - seq_along(si)]) *
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

#' Simulate incidence using the branching process
#'
#' @description This function simulates an incidence trajectory based on the
#' branching process. 
#' 
#' \deqn{I_g = \sum_{i = 1}^{I_{g - 1}}Z_{i, g - 1},}
#'
#' where \eqn{I_g} is the size of generation \eqn{g} and \eqn{Z_{i,g}} is the
#' number of offspring of individual \eqn{i} from generation \eqn{g}. The
#' offspring distribution here is defined as the Neyman type A distribution.
#' 
#' We run the branching process for multiple generations and
#' then incorporate the generation time to go from the generation sizes \eqn{I_g} to
#' the incidence \eqn{X_t} in calendar time \eqn{t}.
#' 
#' Apart from the offspring distribution we include additional overdispersion by
#' sampling different reproductive number \eqn{R_{g,i}} per each individual
#' according to the log-normal distribution, that is
#' 
#' \deqn{R_{g, i} \sim Lognormal(R, \sigma_R^2),}
#' 
#' where \eqn{R} is the true value of the reproductive number.
#' 
#' As the final source of overdispersion, we introduce underreporting in the
#' form of binomial thinning, where we report on average only a given proportion
#' of the final incidence.
#'
#' @param init an integer, the initial number of infectious individual, i.e. the
#'   size of the 0-th generation of the branching process
#' @param R positive real value, the true value of the effective reproduction
#'   number, same across the whole trajectory
#' @param si a non-negative vector, the distribution of the serial interval,
#'   must not be of the same length, or longer than \code{lgt}
#' @param lgt positive integer, total length of the trajectory
#' @param offspring_disp a positive real value larger than 1 indicating the
#'   degree of overdispersion in the offspring distribution. The variance of the
#'   offspring distribution is \code{R * offspring_disp}
#' @param R_sd positive real value, the standard deviation of the log-normal
#'   distribution on the individual effective reproductive numbers
#' @param reporting_prob proportion of the incidence that gets reported on
#'   average
#' @return a named list with two elements
#'   \describe{
#'     \item{\code{X}}{integer vector, the simulated incidence}
#'     \item{\code{Lambda}}{numeric vector, value of the sole covariate in the
#'     renewal equation}
#'  }
simulate_branching <- function(
  init,
  R,
  si,
  lgt,
 # offspring_disp = 2,
  R_sd = 0.1,
  reporting_prob = 1
) {
  # Create a vector of discrete generation times, from which we will draw using
  # `sample()`
  gen_time_to_sample <- seq_along(si)
  # We sample the individual R_t values from a log-normal distribution with
  # mean `R` and standard deviation `R_sd`. However the `rlnorm()` function
  # accepts parameters on the scale of the non-transformed normal distribution.
  # We have to therefore reparametrize.
  R_sdlog <- sqrt(log(R_sd^2 / R^2 + 1))
  R_meanlog <- log(R) - R_sdlog^2 / 2

  # We store the branching process as a list. Each list element corresponds to
  # one generation with a data frame containing the information about each
  # individual in the given generation.
  df_individuals <- list(
    data.frame(
      # Number of the individual in this generation
      individual = seq_len(init),
      # Individual-specific R
      R_sampled = rlnorm(init, meanlog = R_meanlog, sdlog = R_sdlog),
      # Time of infection
      t = 0
    )
  )
  # We need to stop the simulation as soon as all infections until time `lgt`
  # (included) are generated. `t_actual` stores the calendar time, until which
  # all branching events have been resolved, i.e. no individual before this time
  # can be infected anymore.
  t_actual <- 0
  # Generation counter to index the list with data frames.
  gen <- 1
  while (t_actual < lgt) {
    # How many offspring the individuals from the previous generation generate
    offspring <- rpois(nrow(df_individuals[[gen]]), df_individuals[[gen]]$R_sampled)
    sum_offspring <- sum(offspring)
    # New generation
    df_individuals[[gen + 1]] <- data.frame(
      # Number of the individual in this generation
      individual = seq_len(sum_offspring),
      # Individual-specific R
      R_sampled = exp(rnorm(sum_offspring, mean = R_meanlog, sd = R_sdlog)),
      # Time of infection is the time of infection of the parent plus a random
      # generation time
      t = rep(df_individuals[[gen]]$t, times = offspring) + 
        sample(
          gen_time_to_sample,
          size = sum_offspring,
          replace = TRUE,
          prob = si
        )
    )
    # We take the minimum infection time of the newly infected individuals as
    # the actual time, since all new infection events will occur only after this
    # time point. 
    t_actual <- min(c(df_individuals[[gen + 1]]$t, lgt))
    # Increase the generation counter
    gen <- gen + 1
  }

  # Aggregate the infections into the incidence curve.
  df_incid <- dplyr::bind_rows(df_individuals) |>
    dplyr::group_by(t) |>
    dplyr::summarise(incidence = dplyr::n())
  # We might have some times, where no infection occurred. To pad the incidence
  # curve by zeros, we will have to do a join with a complete dataframe.
  df_time <- data.frame(t = seq_len(lgt + 1) - 1)
  incid <- dplyr::left_join(df_time, df_incid, by = "t") |>
    dplyr::pull(incidence) |>
    tidyr::replace_na(0)

  # Add underreporting as binomial thinning. We drop the initialization here.
  underrep_incid <- rbinom(lgt, incid[-1], reporting_prob)

  # Pre-calculate the Lambda to match the output from the `simulate_renewal()`
  # function. This Lambda will be then used to fit the renewal model. Note that
  # this quantity is not used in the branching process.
  Lambda <- c(
    # NAs at the beginning
    rep(NA, length(si)),
    sapply(
      seq_len(lgt - length(si)),
      function(ind) sum(underrep_incid[ind - 1 + seq_along(si)] * rev(si))
    )
  )
  list(X = underrep_incid, Lambda = Lambda)
}

#' Generate multiple incidence trajectories
#'
#' @description This function simulates multiple incidence trajectories based on
#'  the renewal equation.
#'
#' @param n_sim an integer, number of trajectories to generate
#' @param init a vector or a scalar. For the renewal equation models "Poiss",
#'   "NegBin-Q" and "NegBin-L", this is a vector of the incidence at the
#'   beginning of the incidence trajectory. For the branching process model,
#'   this is an integer representing the size of the initial branching process
#'   generation.
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
  n_burnin,
  init,
  R_eff,
  si,
  lgt,
  model = c("Poiss", "NegBin-Q", "NegBin-L", "Branching"),
  nb_size = NA,
  weekday_effect = c(1, 1),
  R_sd = NA,
  reporting_prob = NA,
  seed = 432
) {
  model <- match.arg(model)
  if (
    model == "Branching" && (is.na(R_sd) || is.na(reporting_prob))
  ) {
    stop("For the branching process option, 'R_sd' and 'reporting_prob' must be specified") # nolint
  } else if ((model == "NegBin-Q" || model == "NegBin-L") && is.na(nb_size)) {
    stop("For the negative binomial options, 'nb_size' must be specified.")
  }

  set.seed(seed)
  if (model == "Branching") {
    trajectories <- replicate(
      n_sim,
      simulate_branching(
        init,
        R_eff,
        si,
        lgt + n_burnin,
        R_sd,
        reporting_prob
      ),
      simplify = FALSE
    )
  } else {
    trajectories <- replicate(
      n_sim,
      simulate_renewal(init, R_eff, si, lgt + n_burnin + length(init), model, nb_size, weekday_effect),
      simplify = FALSE
    )
  }
  X <- do.call(cbind, map(trajectories, "X"))
  Lambda <- do.call(cbind, map(trajectories, "Lambda"))
  list(X = X, Lambda = Lambda)
}
