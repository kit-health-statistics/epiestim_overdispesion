library("tidyverse")
library("patchwork")
library("EpiEstim")
library("gamlss")  # For the gamlss() function
theme_set(theme_bw())

# Script illustrating the contributions to the log-likelihood for the 
# Poisson model vs. Negative binomial models. The identity link is used in the
# estimation.

# Read the data ================================================================

window_width <- 13

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

# Select one time window for the illustration
start_date <- as.Date("2022-01-03")
end_date <- as.Date("2022-01-15")

# Subset the incidence
incidence_subset <- incidence %>%
  filter(Date >= start_date & Date <= end_date) %>% 
  mutate(Day = 1:window_width)


# Calculate the serial interval distribution ===================================

# Parameters of the serial interval distribution from (Richter et al., 2022)
std_si <- 1.83
mean_si <- 3.37

# Discretization taken from EpiEstim 
si_distr <- discr_si(seq(0, nrow(incidence) - 1), mean_si, std_si)

# Calculate Lambdas by weighting the past incidence ============================

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
model_mat <- as.data.frame(pairs)

# Fit the models ===============================================================

# Poisson model
mod_pois <- glm(
  data = model_mat,
  Cases ~ lambda - 1, 
  family = poisson(link = "identity")
)

# NegBin2 model. We use the log link for the dispersion parameter. 
# Otherwise the algorithm does not converge.
mod_nbin2 <- gamlss(
  data = model_mat,
  formula = Cases ~ lambda - 1,
  family = NBI(mu.link = "identity", sigma.link = "log"), 
  trace = FALSE
)

# NegBin1 model. We use the log link for the dispersion parameter to be
# consistent with the NegBin2 fitting.
mod_nbin1 <- gamlss(
  data = model_mat,
  formula = Cases ~ lambda - 1,
  family = NBII(mu.link = "identity", sigma.link = "log"), 
  trace = FALSE
)

# Extract the results ==========================================================

# Extract the coefficients
R_hat <- list(
  pois = mod_pois$coefficients,
  nbin2 = mod_nbin2$mu.coefficients,
  nbin1 = mod_nbin1$mu.coefficients
)

# Extract the dispersion parameters
disp_pars <- list(
  nbin2 = exp(mod_nbin2$sigma.coefficients[1]),
  nbin1 = exp(mod_nbin1$sigma.coefficients[1])
)

# Extract the value of the loglikelihood from the AIC values
ML <- list(
  pois = (2 - mod_pois$aic) / 2,
  nbin2 = (4 - mod_nbin2$aic) / 2,
  nbin1 = (4 - mod_nbin1$aic) / 2
)

# # Sanity check of the loglikelihood values
# ML2 <- list(
#   pois = sum(dpois(
#     incidence_subset$Cases, 
#     lambda = lambda_subset * R_hat$pois, 
#     log = TRUE
#   )),
#   nbin2 = sum(dnbinom(
#     incidence_subset$Cases,
#     mu = lambda_subset * R_hat$nbin2, 
#     size = 1 / disp_pars$nbin2, 
#     log = TRUE
#   )),
#   nbin1 = sum(dnbinom(
#     incidence_subset$Cases, 
#     mu = lambda_subset * R_hat$nbin1, 
#     size = lambda_subset * R_hat$nbin1 / disp_pars$nbin1, 
#     log = TRUE
#   ))
# )
# 
# all(unlist(ML) - unlist(ML2) < 1e-12)  # Fine

# Prepare data frames for plotting =============================================

# Sequences for plotting
R_seq_lgt <- 1000
R_seq <- seq(0.3, 3.6, length = R_seq_lgt)
means <- c(R_seq %*% t(lambda_subset))
obs_long_vec <- rep(incidence_subset$Cases, each = R_seq_lgt)

# Create the data frames with the likelihood curves per observation

# Poisson
df_pois_llik <- tibble(
  R = rep(R_seq, times = window_width),
  llik = dpois(obs_long_vec, means, log = TRUE),
  Day = factor(rep(1:window_width, each = R_seq_lgt)),
  log_likelihood_of = "Individual obs.",
  model = "Poisson"
) %>% group_by(Day) %>% 
  mutate(
    # Find the maximum value of the log-likelihood for each curve to shift its
    # peak to 0.
    max_llik = max(llik),
    rel_llik = llik - max_llik
  ) %>% ungroup()
df_pois_llik_total <- df_pois_llik %>% group_by(R) %>% 
  summarise(
    llik = sum(llik),
    max_llik = ML$pois,
    rel_llik = llik - max_llik,
    log_likelihood_of = "Total",
    Day = "Total",
    model = "Poisson"
  )

# NegBin2
df_nbin2_llik <- tibble(
  R = rep(R_seq, times = window_width),
  llik = dnbinom(
    obs_long_vec, 
    mu = means, 
    size = 1 / disp_pars$nbin2,  # Same for all curves
    log = TRUE
  ),
  Day = factor(rep(1:window_width, each = R_seq_lgt)),
  log_likelihood_of = "Individual obs.",
  model = "NegBin2"
) %>% group_by(Day) %>% 
  mutate(
    # Find the maximum value of the log-likelihood for each curve to shift its
    # peak to 0.
    max_llik = max(llik),
    rel_llik = llik - max_llik
  ) %>% ungroup()
df_nbin2_llik_total <- df_nbin2_llik %>% group_by(R) %>% 
  summarise(
    llik = sum(llik),
    max_llik = ML$nbin2,
    rel_llik = llik - max_llik,
    log_likelihood_of = "Total",
    Day = "Total",
    model = "NegBin2"
  )

# NegBin1
df_nbin1_llik <- tibble(
  R = rep(R_seq, times = window_width),
  llik = dnbinom(
    obs_long_vec, 
    mu = means, 
    size = means / disp_pars$nbin1,  # Dispersion parameter same for all curves
    log = TRUE
  ),
  Day = factor(rep(1:window_width, each = R_seq_lgt)),
  log_likelihood_of = "Individual obs.",
  model = "NegBin1"
) %>% group_by(Day) %>% 
  mutate(
    # Find the maximum value of the log-likelihood for each curve to shift its
    # peak to 0.
    max_llik = max(llik),
    rel_llik = llik - max_llik
  ) %>% ungroup()
df_nbin1_llik_total <- df_nbin1_llik %>% group_by(R) %>% 
  summarise(
    llik = sum(llik),
    max_llik = ML$nbin1,
    rel_llik = llik - max_llik,
    log_likelihood_of = "Total",
    Day = "Total",
    model = "NegBin1"
  )

# Prepare the data frame with the point estimates
df_R_hat <- data.frame(
  R_hat = unlist(R_hat),
  model = factor(
    c("Poisson", "NegBin2", "NegBin1"),
    levels = c("Poisson", "NegBin1", "NegBin2")
  )
)

df_llik <- rbind(
  df_pois_llik, df_nbin1_llik, df_nbin2_llik,
  df_pois_llik_total, df_nbin1_llik_total, df_nbin2_llik_total
) %>% 
  mutate(model = factor(model, levels = c("Poisson", "NegBin1", "NegBin2")))


# Generate the black-and-white illustrative plot ===============================

p_incidence_bw <- ggplot() +
  geom_segment(
    data = incidence_subset, 
    aes(x = Date, xend = Date, y = 0, yend = Cases)
  ) +
  geom_vline(aes(xintercept = start_date, linetype = "Cases"), alpha = 0) +
  geom_line(
    aes(x = incidence_subset$Date, y = lambda_subset, color = "lambda"),
    key_glyph = "timeseries"
  ) +
  scale_color_manual(
    values = "gray", 
    labels = expression(Lambda[t]),
    name = ""
  ) +
  scale_linetype_manual(
    values = "solid", 
    name = "",
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  guides(color = guide_legend(ncol = 2)) +
  labs(y = "Cases", color = "Day") +
  scale_x_date(date_labels = "%d %b ")

p_loglik_individual_bw <- ggplot() +
  geom_line(
    data = df_llik, 
    aes(
      x = R, 
      y = rel_llik, 
      group = Day, 
      alpha = log_likelihood_of, 
      linewidth = log_likelihood_of
    )
  ) +
  geom_vline(
    data = df_R_hat, 
    aes(xintercept = R_hat, linetype = "R_ML"),
    linewidth = 0.3
  ) +
  coord_cartesian(ylim = c(-20, 0), xlim = c(0.4, 2.4)) +
  scale_linetype_manual(
    values = "dashed", 
    labels = expression(hat(R)[ML]),
    name = ""
  ) +
  scale_linewidth_manual(
    values = c("Total" = 0.4, "Individual obs." = 0.3), 
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c("Individual obs." = 0.3, "Total" = 1), 
    name = "Relative\nlog-likelihood",
    labels = c("Individual \nobs.", "Total")
  ) +
  labs(
    y = "Relative log-likelihood",
    x = "R"
  ) +
  facet_wrap(~model, nrow = 3) +
  theme(strip.background = element_blank())

p_incidence_loglik_bw <- p_incidence_bw + p_loglik_individual_bw + 
  plot_layout(design = "A\nB\nB\nB\nB")

ggsave("figure/Incidence_and_llik_bw_identity.pdf", p_incidence_loglik_bw, width = 6, 
       height = 7)

# Generate the illustrative plot in color ======================================

manual_hues <- seq(15, 375, length = window_width)
manual_colors <- c(hcl(h = manual_hues, l = 65, c = 100), "#000000")
names(manual_colors) <- c(1:13, "Total")

p_incidence <- ggplot() +
  geom_segment(
    data = incidence_subset, 
    aes(x = Date, xend = Date, y = 0, yend = Cases, color = factor(Day))
  ) +
  geom_line(
    aes(x = incidence_subset$Date, y = lambda_subset, linetype = "lambda"),
    color = "black"
  ) +
  scale_linetype_manual(
    values = "solid", 
    labels = expression(Lambda[t]),
    name = ""
  ) +
  guides(color = guide_legend(ncol = 2)) +
  labs(y = "Cases", color = "Day") +
  scale_x_date(date_labels = "%d %b ")

p_loglik_individual <- ggplot() +
  geom_line(
    data = df_llik, 
    aes(
      x = R, 
      y = rel_llik, 
      group = Day, 
      color = Day, 
      alpha = log_likelihood_of, 
      linewidth = log_likelihood_of
    )
  ) +
  geom_vline(
    data = df_R_hat, 
    aes(xintercept = R_hat, linetype = "R_ML"),
    linewidth = 0.3
  ) +
  coord_cartesian(ylim = c(-20, 0), xlim = c(0.4, 2.4)) +
  scale_linetype_manual(
    values = "dashed", 
    labels = expression(hat(R)[ML]),
    name = ""
  ) +
  scale_linewidth_manual(
    values = c("Total" = 0.4, "Individual obs." = 0.3), 
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c("Individual obs." = 1, "Total" = 1), 
    name = "Relative\nlog-likelihood",
    labels = c("Individual \nobs.", "Total")
  ) +
  scale_color_manual(values = manual_colors, guide = "none") +
  labs(
    y = "Relative log-likelihood",
    x = "R"
  ) +
  facet_wrap(~model, nrow = 3) +
  theme(strip.background = element_blank())

p_incidence_loglik <- p_incidence + p_loglik_individual + 
  plot_layout(design = "A\nB\nB\nB\nB")

ggsave("figure/Incidence_and_llik_identity.pdf", p_incidence_loglik, width = 6,
       height = 7)

# Extra plot comparing the curvature of the total log-likelihood ===============

# We need to set a shorter range for R to make the plot nicer
R_seq_total <- seq(1, 2, length = R_seq_lgt)
means_total <- c(R_seq_total %*% t(lambda_subset))

# Poisson
df_pois_llik_total2 <- tibble(
  R = rep(R_seq_total, times = window_width),
  llik = dpois(obs_long_vec, means_total, log = TRUE)
) %>% group_by(R) %>% 
  summarise(
    llik = sum(llik),
    rel_llik = llik - ML$pois,
    model = "Poisson"
  )

# NegBin2
df_nbin2_llik_total2 <- tibble(
  R = rep(R_seq_total, times = window_width),
  llik = dnbinom(
    obs_long_vec, 
    mu = means_total, 
    size = 1 / disp_pars$nbin2,
    log = TRUE
  )
) %>% group_by(R) %>% 
  summarise(
    llik = sum(llik),
    rel_llik = llik - ML$nbin2,
    model = "NegBin2"
  )

# NegBin1
df_nbin1_llik_total2 <- tibble(
  R = rep(R_seq_total, times = window_width),
  llik = dnbinom(
    obs_long_vec, 
    mu = means_total, 
    size = means_total / disp_pars$nbin1,
    log = TRUE
  )
) %>% group_by(R) %>% 
  summarise(
    llik = sum(llik),
    rel_llik = llik - ML$nbin1,
    model = "NegBin1"
  )

df_llik_total2 <- rbind(df_pois_llik_total2, df_nbin2_llik_total2, 
                        df_nbin1_llik_total2) %>% 
  mutate(
    model = factor(model, levels = c("Poisson", "NegBin1", "NegBin2"))
  )

p_curvature <- ggplot() +
  geom_line(data = df_llik_total2, aes(x = R, y = rel_llik, color = model)) +
  geom_point(
    data = df_R_hat,
    aes(x = R_hat, y = 0, color = model, shape = model)
  ) +
  labs(
    shape = expression(hat(R)[ML]),
    y = "Relative log-likelihood",
    x = "R"
  ) +
  scale_shape_manual(values = c("NegBin1" = 20, "NegBin2" = 15, "Poisson" =  17)) +
  scale_color_manual(
    name = "Model",
    values = c("NegBin1" = "dodgerblue", "NegBin2" = "firebrick3",
               "Poisson" = "forestgreen"),
    guide = guide_legend(
      override.aes = list(
        shape = c("NegBin1" = 20, "NegBin2" = 15, "Poisson" =  17)
      )
    )
  ) +
  coord_cartesian(ylim = c(-10, 0), xlim = c(1.05, 1.95))

ggsave("figure/rel_llik_curvature_identity.pdf", p_curvature, width = 7, height = 3.5)

