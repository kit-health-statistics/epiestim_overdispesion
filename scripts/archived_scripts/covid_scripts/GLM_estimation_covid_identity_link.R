library("readxl")
library("tidyverse")
library("patchwork")
library("EpiEstim")
library("gamlss")  # For the gamlss() function
theme_set(theme_bw())

# Script estimating Rt for COVID-19 in Austria using GLMs with the identity link
# (October 2021 - March 2022 inclusive).

# Read the data ================================================================

# Data from the ECDC web:
# https://www.ecdc.europa.eu/en/publications-data/data-daily-new-cases-covid-19-eueea-country
incidence <- read_csv("data/covid/covid_ecdc.csv") %>%
  filter(geoId == "AT") %>%
  dplyr::select(dateRep, cases) %>%
  rename("Date" = "dateRep", "Cases" = "cases") %>%
  mutate(
    Date = as.Date(Date, format = "%d/%m/%Y")
  )
incidence <- incidence[nrow(incidence):1, ]  # Order the data chronologically

# Subset the data ==============================================================

start_date <- as.Date("2021-10-01")
end_date <- as.Date("2022-03-31")
incidence_subset <- incidence %>% 
  filter(Date >= start_date & Date <= end_date)

# Plot the subset data =========================================================

p_incidence <- ggplot(
  incidence_subset, 
  aes(x = Date, y = Cases, color = "Cases")
) +
  geom_line(key_glyph = "timeseries", linewidth = 0.3) +
  scale_color_manual(values = "black", name = "") +
  scale_x_date(date_labels = "%b %Y") +
  labs(title = "COVID-19 incidence") +
  theme(plot.title = element_text(hjust = 0.5))
print(p_incidence)

# Prepare the model parameters and matrices ====================================

n_obs <- nrow(incidence_subset)  # How many observations
window_width <- 13  # Size of the estimation window

# Calculate the starting and end points of the time windows
t_starts <- 2:(n_obs - window_width + 1)
t_ends <- t_starts + window_width - 1

# Prior of Rt - gamma distribution with shape and scale parameters.
# Parameter values taken from (AGES; Richter et al. 2020)
a_R_prior <- 1
scale_R_prior <- 5

# Parameters of the serial interval distribution from (Richter et al., 2022)
std_si <- 1.83
mean_si <- 3.37

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

# Get the AIC values
AIC_vals <- list()
AIC_vals$pois <- models_pois %>% lapply(function (x) x$aic) %>% unlist()
AIC_vals$nbin2 <- models_nbin2 %>% lapply(function (x) x$aic) %>% unlist()
AIC_vals$nbin1 <- models_nbin1 %>% lapply(function (x) x$aic) %>% unlist()


# Plot the boxplots of AIC and AIC differences =================================

df_AIC <- tibble(
  Date = rep(tail(incidence_subset$Date, length(AIC_vals$pois)), times = 3),
  AIC = c(AIC_vals$pois, AIC_vals$nbin2, AIC_vals$nbin1),
  model = rep(c("Poisson", "NegBin2", "NegBin1"), each = length(AIC_vals$pois))
)
df_AIC_diff <- tibble(
  Date = tail(incidence_subset$Date, length(AIC_vals$pois)),
  AIC_diff = AIC_vals$nbin1 - AIC_vals$nbin2,
  models = "NegBin1 - NegBin2"
)

# Inner and outer limits for the plot insets
from <- c(xmin = 0.65, xmax = 2.35, ymin = 150, ymax = 300)
to <- c(xmin = 0.6, xmax = 2.5, ymin = 10000, ymax = 40000)

# Boxplots of AIC values
p_AIC_box <- ggplot(df_AIC, aes(x = model, y = AIC, fill = model)) +
  geom_boxplot(staplewidth = 1, width = 0.6) +
  scale_fill_manual(
    name = "Model",
    values = c("NegBin1" = "dodgerblue", "NegBin2" = "firebrick3", 
               "Poisson" = "forestgreen")
  ) +
  labs(x = "Model") +
  ggmagnify::geom_magnify(from = from, to = to, axes = "y")

# Boxplot of AIC differences between NegBin1 and NegBin2
p_AIC_box_diff_nbin <- ggplot(
  df_AIC_diff,
  aes(x = models, y = AIC_diff, fill = models)
) +
  geom_boxplot(staplewidth = 1, width = 0.7) +
  annotate(
    "segment", x = 0.1, y = 0.1, xend = 0.1, yend = 2.5,
    arrow = arrow(type = "closed", length = unit(0.06, "npc")),
    color = "firebrick3"
  ) +
  annotate(
    "text", x = 0.3, y = 1.3, 
    label = "NegBin2 better",
    color = "firebrick3", 
    angle = 90
  ) +
  annotate(
    "segment", x = 0.1, y = -0.1, xend = 0.1, yend = -2.5,
    arrow = arrow(type = "closed", length = unit(0.06, "npc")),
    color = "dodgerblue"
  ) +
  annotate(
    "text", x = 0.3, y = -1.3, 
    label = "NegBin1 better",
    color = "dodgerblue", 
    angle = 90
  ) +
  scale_fill_manual(
    name = "Model pair",
    values = c("NegBin1 - NegBin2" = "mediumpurple")
  ) +
  scale_y_continuous(breaks = seq(-4, 2, by = 2)) +
  labs(x = "Model pair", y = "AIC difference") +
  coord_cartesian(xlim = c(0.5, 1.4))

# Combine the plots
AIC_boxplots_full_plot <- p_AIC_box + p_AIC_box_diff_nbin + 
  plot_layout(guides = "collect", design = "AAB")
ggsave("figure/AIC_boxplots_all_in_1_identity.pdf", 
       AIC_boxplots_full_plot, width = 7, height = 4)

# Plot the R estimates =========================================================

df_R_hat <- tibble(
  Date = rep(incidence_subset$Date[t_ends], 4),
  R = c(R_hat$pois, R_hat$qpois, R_hat$nbin2, R_hat$nbin1),
  # CI via the endpoint transformation to avoid including negative values
  lwr = c(R_hat$pois - qnorm(0.975) * R_hat_sd$pois,
              R_hat$qpois - qnorm(0.975) * R_hat_sd$qpois,
              R_hat$nbin2 - qnorm(0.975) * R_hat_sd$nbin2,
              R_hat$nbin1 - qnorm(0.975) * R_hat_sd$nbin1),
  upr = c(R_hat$pois + qnorm(0.975) * R_hat_sd$pois,
              R_hat$qpois + qnorm(0.975) * R_hat_sd$qpois,
              R_hat$nbin2 + qnorm(0.975) * R_hat_sd$nbin2,
              R_hat$nbin1 + qnorm(0.975) * R_hat_sd$nbin1),
  Model = factor(rep(
    c("Poiss", "Q-Poiss", "NegBin2", "NegBin1"), 
    each = length(models_pois)
  ))
)

p_nbin1_vs_qpois <- ggplot(
  df_R_hat, 
  aes(x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model,
      alpha = Model)
) +
  geom_line(linewidth = 0.4) +
  geom_ribbon(color = NA) +
  scale_alpha_manual(
    values = c("Poiss" = 0, "Q-Poiss" = 0.4, "NegBin1" = 0.4, "NegBin2" = 0), 
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
  theme(plot.title = element_text(hjust = 0.5))

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
  theme(plot.title = element_text(hjust = 0.5))

p_pois_only <- ggplot(
  df_R_hat, 
  aes(x = Date, y = R, ymin = lwr, ymax = upr, color = Model, fill = Model,
      alpha = Model)
) +
  geom_line(linewidth = 0.4) +
  geom_ribbon(color = NA) +
  scale_alpha_manual(
    values = c("Poiss" = 0.4, "Q-Poiss" = 0, "NegBin1" = 0, "NegBin2" = 0), 
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
    title = "Poisson",
    #  title = expression(Estimates~of~R[t]),
    y = expression(hat(R)[t])
  ) +
  theme(plot.title = element_text(hjust = 0.5))

p_R_hat_comparison <- (p_incidence + p_pois_only + p_nbin1_vs_nbin2 +
                         p_nbin1_vs_qpois) + 
  plot_layout(guides = "collect", axes = "collect")
ggsave("figure/R_hat_comparison_identity.pdf", p_R_hat_comparison, width = 7, height = 4)

# Display the differences between NegBin1 and NegBin2 ==========================

df_nbin_diff_boxplots <- df_R_hat %>% 
  filter(Model %in% c("NegBin1", "NegBin2")) %>% 
  group_by(Date) %>% 
  summarise(
    diffs = list(tibble(R_diff = diff(R), lwr_diff = diff(lwr), 
                        upr_diff = diff(upr)))
  ) %>% unnest(diffs) %>%   # NegBin1 - NegBin2
  pivot_longer(
    R_diff:upr_diff, 
    names_to = "Quantity",
    values_to = "Difference"
  ) %>% 
  mutate(
    Difference = -Difference,
    Quantity = factor(
      Quantity, 
      levels = c("R_diff", "lwr_diff", "upr_diff"),
      labels = c("lwr_diff" = "Lower bound", "upr_diff" = "Upper bound", 
                 "R_diff" = "Point estimate")
    )
  )

# Boxplots of differences between the point estimates and lower and upper bounds
# of the confidence intervals
p_boxplots_values_nbin <- ggplot(
  df_nbin_diff_boxplots,
  aes(x = Quantity, y = Difference, fill = Quantity)
) +
  geom_boxplot(staplewidth = 1, width = 0.6) +
  coord_cartesian(xlim = c(0.5, 3)) +
  annotate(
    "segment", x = 0.1, y = 0.01, xend = 0.1, yend = 0.15,
    arrow = arrow(type = "closed", length = unit(0.02, "npc"))
  ) +
  annotate(
    "text", x = 0.35, y = 0.08,
    label = "NegBin2 higher\nvalues",
    angle = 90
  ) +
  annotate(
    "segment", x = 0.1, y = -0.01, xend = 0.1, yend = -0.15,
    arrow = arrow(type = "closed", length = unit(0.02, "npc"))
  ) +
  annotate(
    "text", x = 0.35, y = -0.08,
    label = "NegBin1 higher\nvalues",
    angle = 90
  )

# Plot of diferences in the point estimates in time
p_diff_values_nbin <- ggplot(
  filter(df_nbin_diff_boxplots, Quantity == "Point estimate"), 
  aes(x = Date, y = Difference, color = Quantity),
) +
  geom_line() +
  annotate(
    "segment", 
    x = as.Date("2021-09-15"), y = 0.01, 
    xend = as.Date("2021-09-15"), yend = 0.15,
    arrow = arrow(type = "closed", length = unit(0.02, "npc"))
  ) +
  annotate(
    "text",
    x = as.Date("2021-09-30"), y = 0.08,
    label = "NegBin2\nhigher",
    angle = 90
  ) +
  annotate(
    "segment", 
    x = as.Date("2021-09-15"), y = -0.01, 
    xend = as.Date("2021-09-15"), yend = -0.15,
    arrow = arrow(type = "closed", length = unit(0.02, "npc"))
  ) +
  annotate(
    "text", x = as.Date("2021-09-30"), y = -0.08,
    label = "NegBin1\nhigher",
    angle = 90
  ) +
  coord_cartesian(xlim = as.Date(c("2021-09-15", "2022-04-01"))) +
  scale_color_manual(
    values = "black",
    name = "", 
    labels = "Difference in\npoint estimates"
  )

p_nbin_diffs <- (p_boxplots_values_nbin / p_diff_values_nbin) + 
  plot_layout(axes = "collect_y", design = "A\nA\nB")
ggsave("figure/Nbin_values_diffs_identity.pdf", p_nbin_diffs, width = 6, height = 7)


# Sanity check: compare estimates from GLM to the official estimates ===========

# Load the Austrian estimates of Rt
R_ests <- read_delim(
  "data/covid/R_eff.csv", 
  delim = ";", 
  locale = locale("de", decimal_mark = ",")
) %>% filter(Datum <= max(incidence$Date)) %>% 
  rename("Date" = "Datum", "R" = "R_eff", "lwr" = "R_eff_lwr", 
         "upr" = "R_eff_upr") %>% 
  mutate(Model = "Official estimate")
R_ests_subset <- R_ests %>% 
  filter(Date >= start_date & Date <= end_date)

df_compare_to_original <- df_R_hat %>% 
  filter(Model %in% c("Poiss", "Official estimate")) %>% 
  rbind(R_ests_subset)

p_compare_to_original <- ggplot(
) +
  scale_color_manual(
    name = "Model",
    values = c("Poiss" = "forestgreen", "Official estimate" = "black"),
    labels = c("Poiss", "Official\nestimate")
  ) +
  geom_line(
    data = df_compare_to_original, 
    aes(x = Date, y = R,color = Model),
    linewidth = 0.4
  )
print(p_compare_to_original)  # Looks OK
# ggsave("figure/R_hat_comparison_to_official.pdf", p_compare_to_original, width = 6, height = 3)
