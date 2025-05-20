library(tidyverse)
library(ggplot2)
library(EpiEstim)
library(gamlss)
library(patchwork)

# NOTES:
# Window size used is 7 days
# SI distribution is from the paper below (mean 3.6 days, sd 1.6 days)

# Data are from:
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011439#pcbi.1011439.s001
# and are saved to /data/flu
incidence <- readRDS("data/flu/daily_flu.rds") |>
  rename("Cases" = "Incidence") |>
  mutate(
    Date = as.Date(Date, format = "%d/%m/%Y")
  )

# Run the models --------------------------------------------------
# Prep everything - window size, SI distribution, etc.
n_obs <- nrow(incidence)
window_width <- 7
t_starts <- 2:(n_obs - window_width + 1)
t_ends <- t_starts + window_width - 1
a_R_prior <- 1
scale_R_prior <- 5
# Parameters of SI distribution same as Nash et al 2023 paper
std_si <- 1.6
mean_si <- 3.6
si_distr <- discr_si(seq(0, nrow(incidence) - 1), mean_si, std_si)
# Compute lambda
lambda <- vector(mode = "numeric", length = nrow(incidence))
lambda[1] <- NA
for (t in seq(2, nrow(incidence))) {
  lambda[t] <- sum(
    si_distr[seq_len(t)] * incidence$Cases[seq(t, 1)])
}
pairs <- cbind(lambda, incidence$Cases)
colnames(pairs) <- c("lambda", "Cases")

# Calculate Rt for each window
models_pois <- models_qpois <- models_nbin2 <- models_nbin1 <- 
  vector(mode = "list", length = length(t_starts))
model_mats <- sapply(
  t_starts, 
  function(x) {pairs[0:(window_width - 1) + x, ]}, simplify = FALSE
) %>% 
  lapply(as.data.frame)

for (k in 1:length(t_starts)) {
  # Poisson model
  models_pois[[k]] <- glm(
    data = model_mats[[k]],
    Cases ~ 1, 
    offset = log(lambda), 
    family = poisson(link ="log")
  )
  
  # Quasipoisson model
  models_qpois[[k]] <- glm(
    data = model_mats[[k]],
    Cases ~ 1, 
    offset = log(lambda), 
    family = quasipoisson(link ="log")
  )
  
  # NegBin2 model (possible to supply the initial value for the dispersion 
  # parameter as the quasipoisson estimates)
  models_nbin2[[k]] <- gamlss(
    data = model_mats[[k]],
    formula = Cases ~ offset(log(lambda)),
    family = NBI(mu.link = "log"), 
    trace = FALSE
  )
  
  # NegBin1 model
  models_nbin1[[k]] <- gamlss(
    data = model_mats[[k]],
    formula = Cases ~ offset(log(lambda)),
    family = NBII(mu.link = "log"), 
    trace = FALSE
  )
}

# Extract the coefficient
R_hat <- list()
R_hat$pois <- models_pois %>% lapply(coefficients) %>% unlist() %>% exp()
R_hat$qpois <- models_qpois %>% lapply(coefficients) %>% unlist() %>% exp()
R_hat$nbin2 <- models_nbin2 %>% lapply(coefficients) %>% unlist() %>% exp()
R_hat$nbin1 <- models_nbin1 %>% lapply(coefficients) %>% unlist() %>% exp()

# Extract the dispersion parameter estimates
disp <- list()
disp$qpois <- models_qpois %>% lapply(summary) %>% 
  lapply(function (x) x$dispersion) %>% unlist()
disp$nbin2 <- models_nbin2 %>% 
  lapply(function (x) {exp(x$sigma.coefficients[1])}) %>% 
  unlist() %>% unname()
disp$nbin1 <- models_nbin1 %>% 
  lapply(function (x) {exp(x$sigma.coefficients[1])}) %>% 
  unlist() %>% unname()

# Extract the standard error on the log scale
log_R_hat_sd <- list()
log_R_hat_sd$pois <- models_pois %>% lapply(summary) %>% 
  lapply(function (x) x$coefficients[2]) %>% 
  unlist()
log_R_hat_sd$qpois <- models_qpois %>% lapply(summary) %>% 
  lapply(function (x) x$coefficients[2]) %>% 
  unlist()
log_R_hat_sd$nbin2 <- models_nbin2 %>% lapply(summary) %>% 
  lapply(function (x) x[3]) %>% 
  unlist()
log_R_hat_sd$nbin1 <- models_nbin1 %>% lapply(summary) %>% 
  lapply(function (x) x[3]) %>% 
  unlist()

# Create plots ----------------------------------------------------
df_R_hat <- tibble(
  Date = rep(incidence$Date[t_ends], 4),
  R = c(R_hat$pois, R_hat$qpois, R_hat$nbin2, R_hat$nbin1),
  # CI via the endpoint transformation to avoid including negative values
  lwr = exp(c(log(R_hat$pois) - qnorm(0.975) * log_R_hat_sd$pois,
              log(R_hat$qpois) - qnorm(0.975) * log_R_hat_sd$qpois,
              log(R_hat$nbin2) - qnorm(0.975) * log_R_hat_sd$nbin2,
              log(R_hat$nbin1) - qnorm(0.975) * log_R_hat_sd$nbin1)),
  upr = exp(c(log(R_hat$pois) + qnorm(0.975) * log_R_hat_sd$pois,
              log(R_hat$qpois) + qnorm(0.975) * log_R_hat_sd$qpois,
              log(R_hat$nbin2) + qnorm(0.975) * log_R_hat_sd$nbin2,
              log(R_hat$nbin1) + qnorm(0.975) * log_R_hat_sd$nbin1)),
  Model = factor(rep(
    c("Poiss", "Q-Poiss", "NegBin2", "NegBin1"), 
    each = length(models_pois)
  ))
)

# 1. Incidence
p_incidence <- ggplot(
  incidence,
  aes(x = Date, y = Cases, color = "Cases")
) +
  geom_line(key_glyph = "timeseries", linewidth = 0.3) +
  scale_color_manual(values = "black", name = "") +
  scale_x_date(date_labels = "%b %Y") +
  labs(title = "Incidence") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

# 2. Poisson vs. QP
p_pois_vs_qpois <- ggplot(
  df_R_hat, 
  aes(x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model,
      alpha = Model)
) +
  geom_line(linewidth = 0.4) +
  geom_ribbon(color = NA) +
  scale_alpha_manual(
    values = c("Poiss" = 0.4, "Q-Poiss" = 0.4, "NegBin1" = 0.0, "NegBin2" = 0), 
    guide = guide_legend(override.aes = list(alpha = 0.4)) 
  ) +
  scale_color_manual(
    name = "Model",
    values = c("Poiss" = "forestgreen", "Q-Poiss" = "gray40", 
               "NegBin2" = "firebrick3", "NegBin1" = "dodgerblue")
  ) +
  scale_fill_manual(
    name = "Model",
    values = c("Poiss" = "forestgreen", "Q-Poiss" = "gray40", 
               "NegBin2" = "firebrick3", "NegBin1" = "dodgerblue")
  ) +
  labs(
    title = "NegBin1 vs. Quasi-Poisson",
    y = expression(hat(R)[t])
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

# 3. NegBin1 vs. NegBin2
p_nbin1_vs_nbin2 <- ggplot(
  df_R_hat, 
  aes(x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model,
      alpha = Model)
) +
  geom_line(linewidth = 0.4) +
  geom_ribbon(color = NA) +
  scale_alpha_manual(
    values = c("Poiss" = 0, "Q-Poiss" = 0, "NegBin1" = 0.4, "NegBin2" = 0.4), 
    guide = guide_legend(override.aes = list(alpha = 0.4)) 
  ) +
  scale_color_manual(
    name = "Model",
    values = c("Poiss" = "forestgreen", "Q-Poiss" = "gray40", 
               "NegBin2" = "firebrick3", "NegBin1" = "dodgerblue")
  ) +
  scale_fill_manual(
    name = "Model",
    values = c("Poiss" = "forestgreen", "Q-Poiss" = "gray40", 
               "NegBin2" = "firebrick3", "NegBin1" = "dodgerblue")
  ) +
  labs(
    title = "NegBin1 vs. NegBin2",
    #  title = expression(Estimates~of~R[t]),
    y = expression(hat(R)[t])
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

# 4. Overdisp Params (just for QP and NegBin1? NegBin2 is not really comparable)
df_disp <- tibble(
  Date = rep(incidence$Date[t_ends], 3),
  Dispersion = c(disp$qpois, disp$nbin1, disp$nbin2),
  Model = factor(rep(c("Q-Poiss", "NegBin1", "NegBin2"), each = length(disp$qpois)))
)
p_disp <- ggplot(
  df_disp,
  aes(x = Date, y = Dispersion, color = Model, group = Model)
) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(
    name = "Model",
    values = c("Q-Poiss" = "gray40", "NegBin1" = "dodgerblue", "NegBin2" = "firebrick3")
  ) +
  labs(
    title = "Overdispersion Parameters",
    y = "Dispersion Parameter"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  ) +
  theme_bw()

# 5. Generation Time Distribution (GTD)
gtd_data <- tibble(
  x = seq(0, 15, by = 1),
  y = si_distr[1:16]
)

p_gtd <- ggplot(gtd_data, aes(x = x, y = y)) +
  geom_histogram(stat = "identity", fill = "#b4e0f1", color = "black", binwidth = 0.1) +
  labs(
    title = "Generation Time Distribution",
    x = "Days",
    y = "Density"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

# Combine all plots into a single figure
p_R_hat_comparison <- (p_incidence / p_pois_vs_qpois / p_nbin1_vs_nbin2 / p_disp / p_gtd) +
  plot_layout(heights = c(1, 1, 1, 1, 1), guides = "collect") &
  plot_annotation(
    title = "Influenza,\nUSA Active Military Personnel,\n2009-2010",
    theme = theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold", lineheight = 0.9)
    )
  )
print(p_R_hat_comparison)
