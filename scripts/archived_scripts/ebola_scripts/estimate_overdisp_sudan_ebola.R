library(tidyverse)
library(ggplot2)
library(EpiEstim)
library(gamlss)
library(patchwork)

# NOTES:
# Window size used is 7 days
# SI distribution is from the paper (mean 13.8 days, sd 7.6 days)

# Data are transcribed by eye from Figure 1 of:
# https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(24)00555-2/fulltext
# date range is 8/8/22 - 11/19/22 
incidence <- read_csv(here::here("data", "ebola", "ebola_2022_sudan.csv")) |>
  rename("Cases" = "cases") |>
  mutate(
    Date = as.Date(date, format = "%m/%d/%Y")
  ) 
# Run the models --------------------------------------------------
# models struggled to fit on early data
# suspect this has to do with lots of zeros
incidence <- incidence |> 
  filter(Date >= "2022-08-27" & Date <= "2022-11-17")
# Prep everything - window size, SI distribution, etc.
n_obs <- nrow(incidence)
window_width <- 7
t_starts <- 2:(n_obs - window_width + 1)
t_ends <- t_starts + window_width - 1
a_R_prior <- 1
scale_R_prior <- 5
# Parameters of SI distribution same as Nash et al 2023 paper
std_si <- 7.6
mean_si <- 13.8
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
  print(k)
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
  
  # NegBin2 model. We use the log link for the dispersion parameter. 
  # Otherwise the algorithm does not converge.
  # models_nbin2[[k]] <- gamlss(
  #   data = model_mats[[k]],
  #   formula = Cases ~ lambda - 1,
  #   family = NBI(mu.link = "identity", sigma.link = "log"),
  #   trace = FALSE
  # )

  # NegBin1 model. We use the log link for the dispersion parameter to be
  # consistent with the NegBin2 fitting.
  # models_nbin1[[k]] <- gamlss(
  #   data = model_mats[[k]],
  #   formula = Cases ~ lambda - 1,
  #   family = NBII(mu.link = "identity", sigma.link = "log"),
  #   trace = FALSE
  # )
}

R_hat <- list()
R_hat$pois <- models_pois %>% lapply(coefficients) %>% unlist()
R_hat$qpois <- models_qpois %>% lapply(coefficients) %>% unlist()
# R_hat$nbin2 <- models_nbin2 %>% lapply(coefficients) %>% unlist()
# R_hat$nbin1 <- models_nbin1 %>% lapply(coefficients) %>% unlist()

# Extract the dispersion parameter estimates
disp <- list()
disp$qpois <- models_qpois %>% lapply(summary) %>% 
  lapply(function (x) x$dispersion) %>% unlist()
# disp$nbin2 <- models_nbin2 %>% 
#   lapply(function (x) {exp(x$sigma.coefficients[1])}) %>% 
#   unlist() %>% unname()
# disp$nbin1 <- models_nbin1 %>% 
#   lapply(function (x) {exp(x$sigma.coefficients[1])}) %>% 
#   unlist() %>% unname()

# Extract the standard error
R_hat_sd <- list()
R_hat_sd$pois <- models_pois %>% lapply(summary) %>% 
  lapply(function (x) x$coefficients[2]) %>% 
  unlist()
R_hat_sd$qpois <- models_qpois %>% lapply(summary) %>% 
  lapply(function (x) x$coefficients[2]) %>% 
  unlist()
# R_hat_sd$nbin2 <- models_nbin2 %>% lapply(summary) %>% 
#   lapply(function (x) x[3]) %>% 
#   unlist()
# R_hat_sd$nbin1 <- models_nbin1 %>% lapply(summary) %>% 
#   lapply(function (x) x[3]) %>% 
#   unlist()
# Plot the R estimates =========================================================

df_R_hat <- tibble(
  Date = rep(incidence$Date[t_ends], 2),
  R = c(R_hat$pois, R_hat$qpois),
  # CI via the endpoint transformation to avoid including negative values
  lwr = c(R_hat$pois - qnorm(0.975) * R_hat_sd$pois,
          R_hat$qpois - qnorm(0.975) * R_hat_sd$qpois),
  upr = c(R_hat$pois + qnorm(0.975) * R_hat_sd$pois,
          R_hat$qpois + qnorm(0.975) * R_hat_sd$qpois),
  Model = factor(rep(
    c("Poiss", "Q-Poiss"), 
    each = length(models_pois)
  ))
) 

# 1. Incidence
p_incidence <- ggplot(
  incidence,
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
# their plot is only for 9/8/2022 onwards 
p_pois_vs_qpois <- ggplot(
  df_R_hat |> filter(Date >= "2022-09-08"), 
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
    title = "Poisson vs. Quasi-Poisson",
    y = expression(hat(R)[t])
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

df_overdisp <- tibble(
  Date = rep(incidence$Date[t_ends], 1),
  overdisp = disp$qpois
)

overdisp_plot <- ggplot(
  df_overdisp, 
  aes(x = Date, y = overdisp)
) +
  geom_line(linewidth = 0.4) +
  labs(
    title = "Overdispersion",
    y = "Overdispersion"
  ) +
  scale_x_date(date_labels = "%b %Y") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()

summary_plot <- p_incidence + p_pois_vs_qpois + overdisp_plot +
  plot_layout(guides = "collect")

