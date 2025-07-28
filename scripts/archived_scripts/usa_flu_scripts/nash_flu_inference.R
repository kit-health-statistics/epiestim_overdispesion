# Replicating the analysis from Nash et al 2023
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011439#pcbi.1011439.s001
# for comparison to the GLM models

# Below code is from their GitHub: https://github.com/rebeccanash/EM_EpiEstim_Nash2023 ----------
## US Influenza analysis

# Packages required
packages <- c("EpiEstim", "scales", "dplyr", "ggplot2",
              "hrbrthemes", "cowplot", "grDevices")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

#############
# Load data #
#############

# Load daily influenza data
daily_inc <- readRDS("data/flu/daily_flu.rds")
ndays <- length(daily_inc$Incidence)
index <- seq(1, ndays)

# Serial interval - Cowling 2009
method <- "parametric_si"
mean_si <- 3.6
sd_si <- 1.6
config <- EpiEstim::make_config(mean_si = mean_si,
                                std_si = sd_si)
t_start <- seq(from = 2, to = ndays - 1, 1)
config_daily <- EpiEstim::make_config(mean_si = mean_si,
                                      std_si = sd_si,
                                      t_start = t_start,
                                      t_end = t_start + 1)

# Estimate R from reported daily incidence (weekly sliding windows)
r_rep_sliding <- EpiEstim::estimate_R(incid = daily_inc$Incidence,
                                      method = method,
                                      config = config)

# Save their R estimates to data -----------------------------------------------
#write.csv(r_rep_sliding$R, "data/flu/nash_r_rep_sliding.csv", row.names = FALSE)