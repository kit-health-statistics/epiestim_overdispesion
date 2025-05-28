analyse_Rt <- function(incidence, start_date, end_date, window_width, mean_si, std_si) {
  incidence_subset <- incidence |> filter(Date >= start_date & Date <= end_date)
  n_obs <- nrow(incidence_subset) # How many observations

  # Calculate the starting and end points of the time windows
  t_starts <- 2:(n_obs - window_width + 1)
  t_ends <- t_starts + window_width - 1
  # Discretization taken from EpiEstim
  si_distr <- discr_si(seq(0, nrow(incidence) - 1), mean_si, std_si)
  # Compute lambda
  lambda <- vector(mode = "numeric", length = nrow(incidence))
  lambda[1] <- NA
  for (t in seq(2, nrow(incidence))) {
    lambda[t] <- sum(
      si_distr[seq_len(t)] * incidence$Cases[seq(t, 1)]
    )
  }
  lambda_subset <- lambda[incidence$Date >= start_date & incidence$Date <= end_date]

  pairs <- cbind(lambda_subset, incidence_subset$Cases)
  colnames(pairs) <- c("lambda", "Cases")

  models_pois <- models_qpois <- models_nbin2 <- models_nbin1 <-
    vector(mode = "list", length = length(t_starts))
  models_nbin1_log <- models_nbin2_log <- models_qpois_log <- models_pois
  model_mats <- sapply(
    t_starts,
    function(x) {
      pairs[0:(window_width - 1) + x, ]
    },
    simplify = FALSE
  ) %>%
    lapply(as.data.frame)

  for (k in 1:length(t_starts)) {
    # Poisson model
    models_pois[[k]] <- glm(
      data = model_mats[[k]],
      Cases ~ lambda - 1,
      family = poisson(link = "identity")
    )

    # Quasipoisson model
    models_qpois[[k]] <- glm(
      data = model_mats[[k]],
      Cases ~ lambda - 1,
      family = quasipoisson(link = "identity")
    )

    # NegBin2 model
    models_nbin2[[k]] <- gamlss(
      data = model_mats[[k]],
      formula = Cases ~ lambda - 1,
      family = NBI(mu.link = "identity", sigma.link = "log"),
      control = gamlss.control(trace = FALSE)
    )

    # NegBin1 model
    models_nbin1[[k]] <- gamlss(
      data = model_mats[[k]],
      formula = Cases ~ lambda - 1,
      family = NBII(mu.link = "identity", sigma.link = "log"),
      control = gamlss.control(trace = FALSE)
    )
  }

  # Extract the coefficient
  R_hat <- list()
  R_hat$pois <- models_pois %>%
    lapply(coefficients) %>%
    unlist()
  R_hat$qpois <- models_qpois %>%
    lapply(coefficients) %>%
    unlist()
  R_hat$nbin2 <- models_nbin2 %>%
    lapply(coefficients) %>%
    unlist()
  R_hat$nbin1 <- models_nbin1 %>%
    lapply(coefficients) %>%
    unlist()

  # Extract the standard error
  R_hat_sd <- list()
  R_hat_sd$pois <- models_pois %>%
    lapply(summary) %>%
    lapply(function(x) x$coefficients[2]) %>%
    unlist()
  R_hat_sd$qpois <- models_qpois %>%
    lapply(summary) %>%
    lapply(function(x) x$coefficients[2]) %>%
    unlist()
  R_hat_sd$nbin2 <- models_nbin2 %>%
    lapply(summary) %>%
    lapply(function(x) x[3]) %>%
    unlist()
  R_hat_sd$nbin1 <- models_nbin1 %>%
    lapply(summary) %>%
    lapply(function(x) x[3]) %>%
    unlist()

  # Extract the dispersion parameter estimates
  disp <- list()
  disp$qpois <- models_qpois %>%
    lapply(summary) %>%
    lapply(function(x) x$dispersion) %>%
    unlist()
  disp$nbin2 <- models_nbin2 %>%
    lapply(function(x) {
      exp(x$sigma.coefficients[1])
    }) %>%
    unlist() %>%
    unname()
  disp$nbin1 <- models_nbin1 %>%
    lapply(function(x) {
      exp(x$sigma.coefficients[1])
    }) %>%
    unlist() %>%
    unname()

  # Get the AIC values
  AIC_vals <- list()
  AIC_vals$pois <- models_pois %>%
    lapply(function(x) x$aic) %>%
    unlist()
  AIC_vals$nbin2 <- models_nbin2 %>%
    lapply(function(x) x$aic) %>%
    unlist()
  AIC_vals$nbin1 <- models_nbin1 %>%
    lapply(function(x) x$aic) %>%
    unlist()


  # Create plots ----------------------------------------------------
  df_R_hat <- tibble(
    Date = rep(incidence_subset$Date[t_ends], 4),
    R = c(R_hat$pois, R_hat$qpois, R_hat$nbin2, R_hat$nbin1),
    # Wald CI
    lwr = c(
      R_hat$pois - qnorm(0.975) * R_hat_sd$pois,
      R_hat$qpois - qnorm(0.975) * R_hat_sd$qpois,
      R_hat$nbin2 - qnorm(0.975) * R_hat_sd$nbin2,
      R_hat$nbin1 - qnorm(0.975) * R_hat_sd$nbin1
    ),
    upr = c(
      R_hat$pois + qnorm(0.975) * R_hat_sd$pois,
      R_hat$qpois + qnorm(0.975) * R_hat_sd$qpois,
      R_hat$nbin2 + qnorm(0.975) * R_hat_sd$nbin2,
      R_hat$nbin1 + qnorm(0.975) * R_hat_sd$nbin1
    ),
    Model = factor(rep(
      c("Poiss", "Q-Poiss", "NegBin2", "NegBin1"),
      each = length(models_pois)
    ))
  )

  # Create plots ----------------------------------------------------
  model_colors <- c(
    "Poiss" = "forestgreen", "Q-Poiss" = "gray40",
    "NegBin1" = "dodgerblue", "NegBin2" = "firebrick3"
  )
  
  # 1. Incidence
  p_incidence <- ggplot(
    incidence_subset,
    aes(x = Date, y = Cases)
  ) +
    geom_line(key_glyph = "timeseries", linewidth = 0.3) +
    scale_color_manual(values = "black", name = "") +
    scale_x_date(date_labels = "%b %Y") +
    labs(title = "Incidence") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    theme_bw()

  # 2. Poisson vs. QP
  p_pois_vs_qpois <- ggplot(
    df_R_hat,
    aes(
      x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model,
      alpha = Model
    )
  ) +
    geom_line(linewidth = 0.4) +
    geom_ribbon(color = NA) +
    scale_alpha_manual(
      values = c("Poiss" = 0.4, "Q-Poiss" = 0.4, "NegBin1" = 0.0, "NegBin2" = 0),
      guide = guide_legend(override.aes = list(alpha = 0.4))
    ) +
    scale_color_manual(
      name = "Model",
      values = model_colors
    ) +
    scale_fill_manual(
      name = "Model",
      values = model_colors
    ) +
    labs(
      title = "Poisson vs. Quasi-Poisson",
      y = expression(hat(R)[t])
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_bw()

  # 3. NegBin1 vs. NegBin2
  p_nbin1_vs_nbin2 <- ggplot(
    df_R_hat,
    aes(
      x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model,
      alpha = Model
    )
  ) +
    geom_line(linewidth = 0.4) +
    geom_ribbon(color = NA) +
    scale_alpha_manual(
      values = c("Poiss" = 0, "Q-Poiss" = 0, "NegBin1" = 0.4, "NegBin2" = 0.4),
      guide = guide_legend(override.aes = list(alpha = 0.4))
    ) +
    scale_color_manual(
      name = "Model",
      values = model_colors
    ) +
    scale_fill_manual(
      name = "Model",
      values = model_colors
    ) +
    labs(
      title = "NegBin1 vs. NegBin2",
      y = expression(hat(R)[t])
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_bw()

  # 4. Overdisp Params. Negbin2 plotted on using a second plot axis
  sec_axis_scale <- median(disp$nbin2 / disp$nbin1)
  disp$nbin2_transformed <- disp$nbin2 / sec_axis_scale
  df_disp <- tibble(
    Date = rep(incidence_subset$Date[t_ends], 3),
    Dispersion = c(disp$qpois, disp$nbin1, disp$nbin2_transformed),
    Model = factor(rep(c("Q-Poiss", "NegBin1", "NegBin2"), each = length(disp$qpois)))
  )

  p_disp <- ggplot(
    df_disp,
    aes(x = Date, y = Dispersion, color = Model, group = Model)
  ) +
    geom_line(linewidth = 0.4) +
    scale_y_continuous(sec.axis = sec_axis(~ . * sec_axis_scale, name = "NegBin2")) +
    scale_color_manual(
      name = "Model",
      values = model_colors[-1]
    ) +
    labs(
      title = "Overdispersion Parameters",
      y = "Overdispersion"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    theme_bw() +
    guides(color = "none")

  # 5. Generation Time Distribution (GTD)
  gtd_data <- tibble(
    x = seq(0, length(si_distr) - 1, by = 1),
    y = si_distr
  ) |> filter((x <= 30 & y > 0.004) | x == 0)

  p_gtd <- ggplot(gtd_data, aes(x = x, y = y)) +
    geom_histogram(stat = "identity", fill = "#f8e161b3", color = "black", binwidth = 0.1) +
    labs(
      title = "Generation Time Distribution",
      x = "Days",
      y = "Density"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_bw()
  
  # 6. NegBin1 vs. QP (Extra plot, will probably go to the supplement)
  p_nbin1_vs_qpois <- ggplot(
    df_R_hat,
    aes(
      x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model,
      alpha = Model
    )
  ) +
    geom_line(linewidth = 0.4) +
    geom_ribbon(color = NA) +
    scale_alpha_manual(
      values = c("Poiss" = 0, "Q-Poiss" = 0.4, "NegBin1" = 0.4, "NegBin2" = 0),
      guide = guide_legend(override.aes = list(alpha = 0.4))
    ) +
    scale_color_manual(
      name = "Model",
      values = model_colors
    ) +
    scale_fill_manual(
      name = "Model",
      values = model_colors
    ) +
    labs(
      title = "NegBin1 vs. NegBin2",
      y = expression(hat(R)[t])
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_bw()

  # Combine all plots into a single figure
  p_R_hat_comparison <- (p_incidence / p_pois_vs_qpois / p_nbin1_vs_nbin2 / p_disp / p_gtd) +
    plot_layout(heights = c(1, 1, 1, 1, 1), guides = "collect") & 
    theme(legend.position = "none")
  
  p_legend <- ggpubr::get_legend(p_nbin1_vs_nbin2)

  ret <- list(
    R_hat = R_hat, R_hat_sd = R_hat_sd, disp = disp,
    plt = p_R_hat_comparison, AIC = AIC_vals,
    plt_legend = p_legend,
    plt_nbin1_vs_qpois = p_nbin1_vs_qpois
  )
  return(ret)
}
