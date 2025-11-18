# estimating overdispersion for ebola data from 2014 outbreak
# data originally from https://github.com/fintzij/stemr/tree/master/data

library("readxl")
library("tidyverse")
library("patchwork")
library("EpiEstim")
library("gamlss")  # For the gamlss() function
theme_set(theme_bw())
# this is all from https://royalsocietypublishing.org/doi/10.1098/rsif.2021.0429
# and https://github.com/willgreen236/Heterogeneity_transmission/blob/main/heterogeneity_function.R#L842
ebola_df_daily <- read_csv(file = here::here("data", "ebola", "green2022_data_frame.csv"))
mu <- 15.3
sd <- 9.3
# shape
alpha <- (mu/sd)^2
#rate
beta <- (mu/sd^2)
# divide by 7 means
scaled_beta <- beta / (1/7)
scaled_mean = alpha/scaled_beta
scaled_sd = sqrt(alpha/(scaled_beta^2))

# create weekly data ------------------------------------------------------
ebola_df <- ebola_df_daily %>%
  group_by(week = lubridate::floor_date(date, "week")) %>%
  summarise(Cases = sum(count)) %>%
  ungroup() %>%
  mutate(week = as.Date(week)) %>%
  filter(week > "2014-02-23") %>%
  mutate(Date = week)
# now modifying code from GLM_estimation_covid.R
# Subset the data ==============================================================
# 4 weeks after start
start_date <- as.Date("2014-03-30")
end_date <- max(ebola_df$Date)
incidence <- ebola_df
incidence_subset <- ebola_df %>% 
  filter(Date >= start_date & Date <= end_date)

source(here::here("scripts", "analyse_Rt.R"))

window_width = 4
mean_si = scaled_mean
std_si = scaled_sd
# run the analysis --------------------------------------------------------

ebola_res <- analyse_Rt(incidence, start_date, end_date, window_width, mean_si, std_si) 
ebola_plot <- (ebola_res$plt &
                 plot_annotation(
                   title = "Ebola,\nGuinea,\n2014",
                   theme = theme(
                     plot.title = element_text(size = 16, hjust = 0.5, face = "bold", lineheight = 0.9)
                   )
                 )) |> wrap_elements()
ggsave("figure/green2022_plot.pdf", ebola_plot, width = 4, height = 13)
