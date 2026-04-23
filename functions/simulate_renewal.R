#' Initialize a simulation trajectory
#'
#' @description This function initializes a simulation trajectory using an iid
#' sample of normal random variables that are made non-negative and rounded
#' to count values.
#'
#' @param n_init integer, the number of initial points to generate
#' @param init_magnitude the mean of the normal distribution to sample from
#' @param init_sd the standard deviation of the normal distribution to sample
#'   from
#' @param seed seed value for sampling the initial values
#' @param model the count distribution, one of "Poiss", "NegBin-Q", "NegBin-L"
#'   and "Branching"
#' @return integer vector of initial values of a simulation trajectory
initialize_trajectory <- function(
  n_init,
  init_magnitude,
  init_sd,
  seed,
  model = c("Poiss", "NegBin-Q", "NegBin-L", "Branching")
) {
  model <- match.arg(model)
  # Set seed ensuring identical initialization for scenarios with the same
  # magnitude. Creates redundant copies (nrow(scenarios) instead of 2) but
  # simplifies pipeline structure with negligible computational cost.
  set.seed(seed)
  abs(round(rnorm(n_init, init_magnitude, init_sd)))
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
#' @param R a vector of positive real values, the values of the effective
#'   reproduction number after the initializiation until the end of the
#'   trajectory
#' @param si a non-negative vector, the distribution of the serial interval,
#'   must not be longer than the \code{init} vector
#' @param lgt positive integer, total length of the trajectory
#' @param model the count distribution, one of "Poiss", "NegBin-Q" and
#'   "NegBin-L"
#' @param nb_size positive real value, the size parameter of the negative
#'   binomial distribution for "NegBin-Q" and "NegBin-L"
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
  nb_size = NULL
) {
  model <- match.arg(model)
  X <- Lambda <- numeric(lgt)
  X[seq_along(init)] <- init
  Lambda[seq_along(init)] <- NA

  for (t in (length(init) + 1):(lgt)) {
    Lambda[t] <- sum(si * X[t - seq_along(si)])
    X[t] <- generate_from_obs_mod(
      R[t - length(init)] * Lambda[t],
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
#' offspring distribution here is defined as a "double" Poisson distribution,
#' meaning, that each individual has a Poisson distributed number of offspring
#' clusters and each offspring cluster contains again a Poisson distributed
#' number of individuals. All individuals from the same clusters become
#' infectious at the same time, i.e. they have the same serial interval
#' distribution.
#'
#' We run the branching process for multiple generations and
#' then incorporate the generation time to go from the generation sizes
#' \eqn{I_g} to the incidence \eqn{X_t} in calendar time \eqn{t}.
#'
#' As an additional source of overdispersion, we introduce underreporting in the
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
  offspring_disp = 2,
  reporting_prob = 1
) {
  # Create a vector of discrete serial intervals/generation times, from which we
  # will draw using `sample()`
  si_to_sample <- seq_along(si)
  # Total trajectory length including the initialization
  lgt_total <- lgt + length(init)

  # We store the branching process as a list. Each list element corresponds to
  # one generation with a data frame containing the information about each
  # individual in the given generation.
  # The first generation is special, as the initially infectious individuals are
  # distributed across multiple calendar time points.
  df_individuals <- sapply(
    seq_along(init),
    function(ind) {
      data.frame(
        # Number of the individual in this generation
        individual = seq_len(init[ind]),
        # Time of infection
        t = rep(ind - length(init), init[ind])
      )
    },
    simplify = FALSE
  ) |>
    bind_rows() |>
    list()

  # We need to stop the simulation as soon as all infections until time `lgt`
  # (included) are generated. `t_actual` stores the calendar time, until which
  # all branching events have been resolved, i.e. no individual before this time
  # can be infected anymore.
  t_actual <- 0
  # Generation counter to index the list with data frames.
  gen <- 1
  while (t_actual < lgt) {
    # How many offspring clusters the individuals from the previous generation
    # generate
    clusters <- rpois(
      nrow(df_individuals[[gen]]), R / (offspring_disp - 1)
    )
    # How many individuals are there in the cluster.
    inner_offspring <- sapply(clusters, rpois, lambda = offspring_disp - 1) |>
      unlist()
    # New generation
    df_individuals[[gen + 1]] <- data.frame(
      # Number of individuals in this generation
      individual = seq_len(sum(inner_offspring)),
      # Time of infection is the time of infection of the parent plus a random
      # generation time. The generation time is sampled once for a whole
      # cluster.
      t = rep(
        rep(df_individuals[[gen]]$t, times = clusters) +
          sample(si_to_sample, size = sum(clusters), replace = TRUE, prob = si),
        times = inner_offspring
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
    dplyr::mutate(t = t + length(init)) |>  # Shift the calendar time
    dplyr::group_by(t) |>
    dplyr::summarise(incidence = dplyr::n())
  # We might have some times, where no infection occurred. To pad the incidence
  # curve by zeros, we will have to do a join with a complete data frame.
  df_time <- data.frame(t = seq_len(lgt_total))
  incid <- dplyr::left_join(df_time, df_incid, by = "t") |>
    dplyr::pull(incidence) |>
    tidyr::replace_na(0)

  # Add underreporting as binomial thinning, that is, each case has a fixed
  # probability to be reported.
  underrep_incid <- rbinom(lgt_total, incid, reporting_prob)

  # Pre-calculate the Lambda to match the output from the `simulate_renewal()`
  # function. This Lambda will be then used to fit the renewal model. Note that
  # this quantity is not used to calculate any characteristics of the branching
  # process.
  Lambda <- c(
    # NAs at the beginning
    rep(NA, length(si)),
    sapply(
      seq_len(lgt_total - length(si)),
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
#' @param n_burnin an integer, the length of the burn-in period before the
#'   first estimation window begins.
#' @param init a vector of the incidence at the beginning of the incidence
#'   trajectory.
#' @param R_eff positive real value, the value of the effective reproduction
#'   number, same across the whole trajectory
#' @param si a non-negative vector, the distribution of the serial interval,
#'   must not be longer than the \code{init} vector
#' @param lgt positive integer, total length of one trajectory. All trajectories
#'   are of equal length
#' @param model the count distribution, one of "Poiss", "NegBin-Q", "NegBin-L"
#'   and "Branching"
#' @param nb_size positive real value, the size parameter of the negative
#'   binomial distribution for "NegBin-Q" and "NegBin-L"
#' @param offspring_disp a positive real value larger than 1 indicating the
#'   degree of overdispersion in the offspring distribution. The variance of the
#'   offspring distribution is \code{R_eff * offspring_disp}
#' @param reporting_prob a real value between 0 and 1, the proportion of
#'   reported cases in the branching process simulation
#' @param seed an integer, seed to be set before the `n_sim` replications of the
#'   trajectory simulation
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
  offspring_disp = NA,
  reporting_prob = NA,
  seed = 432
) {
  model <- match.arg(model)
  if (model == "Branching") {
    if (is.na(offspring_disp) || is.na(reporting_prob)) {
      stop("For the branching process option, 'offspring_disp' and 'reporting_prob' must be specified") # nolint
    }
    if (offspring_disp <= 1) {
      stop("For the branching process option, 'offspring_disp' must be larger than 1") # nolint
    }
  } else if ((model == "NegBin-Q" || model == "NegBin-L") && is.na(nb_size)) {
    stop("For the negative binomial options, 'nb_size' must be specified.")
  }

  set.seed(seed)
  if (model == "Branching") {
    trajectories <- replicate(
      n_sim,
      simulate_branching(
        init,
        # When `R_eff` is passed as a vector, the branching process will use
        # only the first element, which will be regarded as a constant value of
        # the reproductive number throughout the whole trajectory.
        R_eff[1],
        si,
        lgt + n_burnin,
        offspring_disp,
        reporting_prob
      ),
      simplify = FALSE
    )
  } else {
    trajectories <- replicate(
      n_sim,
      simulate_renewal(
        init,
        R_eff,
        si,
        lgt + n_burnin + length(init),
        model,
        nb_size
      ),
      simplify = FALSE
    )
  }
  X <- do.call(cbind, map(trajectories, "X"))
  Lambda <- do.call(cbind, map(trajectories, "Lambda"))
  list(X = X, Lambda = Lambda)
}
