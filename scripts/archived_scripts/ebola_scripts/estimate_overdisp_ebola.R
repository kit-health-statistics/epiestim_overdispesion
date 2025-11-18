# estimating overdispersion for ebola data from 2014 outbreak
# data originally from https://github.com/fintzij/stemr/tree/master/data

library("readxl")
library("tidyverse")
library("patchwork")
library("EpiEstim")
library("gamlss")  # For the gamlss() function
theme_set(theme_bw())

ebola_df <- read_csv(file = here::here("data", "ebola", "ebola_df.csv"))

# mean and sd of ebola gen time will be 12 and 5.2 days
# midpoint weibull estimate is scale 13.6 shape 2.6 
# taken from Chowell and Nishiura (2014) Transmission dynamics and control of Ebola virus disease (EVD): a review 
# the data are in weeks
# converted weibull should be scale 13.6/7, shape 2.6
set.seed(1)
week_samps <- rweibull(10000, shape = 2.6, scale = 13.6/7)
mean(week_samps)
# 1.73
sd(week_samps)
# 0.72

# now modifying code from GLM_estimation_covid.R
# Subset the data ==============================================================

start_date <- as.Date("2014-03-30")
end_date <- max(ebola_df$week)

incidence <- ebola_df %>%
  dplyr::select(Date = week, Cases = guinea)
incidence_subset <- ebola_df %>% 
  filter(week >= start_date) %>%
  dplyr::select(Date = week, Cases = guinea)

# Plot the subset data =========================================================

p_incidence <- ggplot(
  incidence_subset, 
  aes(x = Date, y = Cases, color = "Cases")
) +
  geom_line(key_glyph = "timeseries", linewidth = 0.3) +
  scale_color_manual(values = "black", name = "") +
  scale_x_date(date_labels = "%b %Y") +
  labs(title = "Weekly Ebola incidence in Guinea") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
print(p_incidence)

# Prepare the model parameters and matrices ====================================

n_obs <- nrow(incidence_subset)  # How many observations
window_width <- 13  # Size of the estimation window

# Calculate the starting and end points of the time windows
t_starts <- 2:(n_obs - window_width + 1)
t_ends <- t_starts + 12

# Prior of Rt - gamma distribution with shape and scale parameters.
# Parameter values taken from (AGES; Richter et al. 2020)
a_R_prior <- 1
scale_R_prior <- 5

# Parameters of the serial interval distribution from (Richter et al., 2022)
std_si <- 0.72
mean_si <- 1.73

# Discretization taken from EpiEstim 
si_distr <- discr_si(seq(0, nrow(incidence) - 1), mean_si, std_si)

# Calculate Lambda for all time points (adapted from EpiEstim)
lambda <- vector(mode = "numeric", length = nrow(incidence))
lambda[1] <- NA
for (t in seq(2, nrow(incidence))) {
  lambda[t] <- sum(
    si_distr[seq_len(t)] * incidence$Cases[seq(t, 1)])
}
lambda_subset <- lambda[incidence$Date >= start_date & incidence$Date <= end_date]

# Create the pairs of the mean values and the observations
pairs <- cbind(lambda_subset, incidence_subset$Cases)
colnames(pairs) <- c("lambda", "Cases")

# Calculate Rt for all time windows using GLM with the identity link ===========

models_pois <- models_qpois <- models_nbin2 <- models_nbin1 <- 
  vector(mode = "list", length = length(t_starts))
models_nbin1_log <- models_nbin2_log <- models_qpois_log <- models_pois
model_mats <- sapply(
  t_starts, 
  function(x) {pairs[0:(window_width - 1) + x, ]}, simplify = FALSE
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
  
  # NegBin2 model. We use the log link for the dispersion parameter. 
  # Otherwise the algorithm does not converge.
  models_nbin2[[k]] <- gamlss(
    data = model_mats[[k]],
    formula = Cases ~ lambda - 1,
    family = NBI(mu.link = "identity", sigma.link = "log"), 
    trace = FALSE
  )
  
  # NegBin1 model. We use the log link for the dispersion parameter to be
  # consistent with the NegBin2 fitting.
  models_nbin1[[k]] <- gamlss(
    data = model_mats[[k]],
    formula = Cases ~ lambda - 1,
    family = NBII(mu.link = "identity", sigma.link = "log"), 
    trace = FALSE
  )
}

# Extract the coefficient
R_hat <- list()
R_hat$pois <- models_pois %>% lapply(coefficients) %>% unlist()
R_hat$qpois <- models_qpois %>% lapply(coefficients) %>% unlist()
R_hat$nbin2 <- models_nbin2 %>% lapply(coefficients) %>% unlist()
R_hat$nbin1 <- models_nbin1 %>% lapply(coefficients) %>% unlist()

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

# Extract the standard error
R_hat_sd <- list()
R_hat_sd$pois <- models_pois %>% lapply(summary) %>% 
  lapply(function (x) x$coefficients[2]) %>% 
  unlist()
R_hat_sd$qpois <- models_qpois %>% lapply(summary) %>% 
  lapply(function (x) x$coefficients[2]) %>% 
  unlist()
R_hat_sd$nbin2 <- models_nbin2 %>% lapply(summary) %>% 
  lapply(function (x) x[3]) %>% 
  unlist()
R_hat_sd$nbin1 <- models_nbin1 %>% lapply(summary) %>% 
  lapply(function (x) x[3]) %>% 
  unlist()


# Plot the R estimates =========================================================

df_R_hat <- tibble(
  Date = rep(incidence_subset$Date[t_ends], 2),
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



pois_vs_qpois <- ggplot(
  df_R_hat, 
  aes(x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model,
      alpha = Model)
) +
  geom_line(linewidth = 0.4) +
  geom_ribbon(color = NA) +
  scale_alpha_manual(
    values = c("Poiss" = 0.4, "Q-Poiss" = 0.4), 
    guide = guide_legend(override.aes = list(alpha = 0.4)) 
  ) +
  scale_color_manual(
    name = "Model",
    values = c("Poiss" = "forestgreen", "Q-Poiss" = "gray40")
  ) +
  scale_fill_manual(
    name = "Model",
    values = c("Poiss" = "forestgreen", "Q-Poiss" = "gray40")
  ) +
  labs(
    title = "Poisson vs. Quasi-Poisson",
    y = expression(hat(R)[t])
  ) +
  scale_x_date(date_labels = "%b %Y") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot the overdisp estimates ==================================================

df_overdisp <- tibble(
  Date = rep(incidence_subset$Date[t_ends], 1),
  overdisp = disp$qpois
)

overdisp_plot <- ggplot(
  df_overdisp, 
  aes(x = Date, y = overdisp)
) +
  geom_line(linewidth = 0.4) +
  labs(
    title = "Estimated Quasi-Poisson Overdispersion",
    y = "Overdispersion"
  ) +
  scale_x_date(date_labels = "%b %Y") +
  theme(plot.title = element_text(hjust = 0.5))


combined_plot <- p_incidence + overdisp_plot + pois_vs_qpois + plot_layout(guides = "collect")
ggsave(
  filename = here::here("scripts", "ebola_scripts", "ebola_overdispersion.pdf"),
  plot = combined_plot,
  width = 12, height = 6)