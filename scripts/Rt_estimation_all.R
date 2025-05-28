library(tidyverse)
library(ggplot2)
library(EpiEstim)
library(gamlss)
library(patchwork)
source("scripts/analyse_Rt.R")

# Running the analysis on 3 datasets:
# - flu among USA military personnel
# - COVID-19 in Austria
# - ebola in Sudan

# List for storing the parameters for the 3 datasets
params <- replicate(3, list())
names(params) <- c("flu", "covid", "ebola")

# Read the flu data --------------------------------------------------

# Data are from:
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011439#pcbi.1011439.s001
incidence_flu <- read_csv("data/flu/daily_flu.csv")

# NOTES:
# Window size used is 7 days
# SI distribution is from the paper below (mean 3.6 days, sd 1.6 days)

# We run the analysis on the whole set
params$flu$incidence <- incidence_flu
params$flu$start_date <- incidence_flu$Date[1]
params$flu$end_date <- incidence_flu$Date[nrow(incidence_flu)]
params$flu$window_width <- 7
# Parameters of SI distribution same as Nash et al 2023 paper
params$flu$mean_si <- 3.6
params$flu$std_si <- 1.6

# Read the covid data --------------------------------------------------

# Data from the ECDC web:
# https://www.ecdc.europa.eu/en/publications-data/data-daily-new-cases-covid-19-eueea-country
incidence_covid <- read_csv("data/covid/covid_ecdc.csv") %>%
  filter(geoId == "AT") %>%
  dplyr::select(dateRep, cases) %>%
  rename("Date" = "dateRep", "Cases" = "cases") %>%
  mutate(
    Date = as.Date(Date, format = "%d/%m/%Y")
  )
# Order the data chronologically
incidence_covid <- incidence_covid[nrow(incidence_covid):1, ]

# We run the analysis on a subset just to keep it shorter
params$covid$incidence <- incidence_covid
params$covid$start_date <- as.Date("2021-10-01")
params$covid$end_date <- as.Date("2022-03-31")
params$covid$window_width <- 13
# Parameters of the serial interval distribution from (Richter et al, 2022)
params$covid$mean_si <- 3.37
params$covid$std_si <- 1.83

# Read the ebola data --------------------------------------------------

# Data are transcribed by eye from Figure 1 of:
# https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(24)00555-2/fulltext
# date range is 8/8/22 - 11/19/22 
incidence_ebola <- read_csv(here::here("data", "ebola", "ebola_2022_sudan.csv")) |>
  rename("Cases" = "cases") |>
  mutate(
    Date = as.Date(date, format = "%m/%d/%Y")
  )

# NOTES:
# Window size used is 7 days
# SI distribution is from the paper (mean 13.8 days, sd 7.6 days)

# Parameters of SI distribution same as Nash et al 2023 paper
params$ebola$incidence <- incidence_ebola
params$ebola$start_date <- as.Date("2022-09-03")
params$ebola$end_date <- as.Date("2022-09-19")
params$ebola$window_width <- 7
# Parameters of SI distribution same as de Padua et al 2025 paper
params$ebola$mean_si <- 13.8
params$ebola$std_si <- 7.6

# Run the models --------------------------------------------------

flu_results <- with(
  params$flu, 
  analyse_Rt(incidence, start_date, end_date, window_width, mean_si, std_si)
)
covid_results <- with(
  params$covid, 
  analyse_Rt(incidence, start_date, end_date, window_width, mean_si, std_si)
)
ebola_results <- with(
  params$ebola, 
  analyse_Rt(incidence, start_date, end_date, window_width, mean_si, std_si)
)

# Plot the results --------------------------------------------------

flu_plot <- (flu_results$plt &
  plot_annotation(
    title = "Influenza,\nUSA Active Military Personnel,\n2009-2010",
    theme = theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold", lineheight = 0.9)
    )
  )) |> wrap_elements()
covid_plot <- (covid_results$plt &
  plot_annotation(
    title = "COVID-19,\nAustria,\n2021-2022",
    theme = theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold", lineheight = 0.9)
    )
  )) |> wrap_elements()
ebola_plot <- (ebola_results$plt &
  plot_annotation(
    title = "Ebola,\nSudan,\n2022",
    theme = theme(
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold", lineheight = 0.9)
    )
  )) |> wrap_elements()
composite_plot <- (
  flu_plot | covid_plot | ebola_plot | wrap_elements(flu_results$plt_legend)
) +
  plot_layout(widths = c(3, 3, 3, 1))
ggsave("figure/composite_plot.pdf", composite_plot, width = 14, height = 13)
