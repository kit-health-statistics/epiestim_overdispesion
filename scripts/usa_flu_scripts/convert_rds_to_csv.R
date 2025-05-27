library(tidyverse)
# Data are from:
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011439#pcbi.1011439.s001
flu <- readRDS("data/flu/daily_flu.rds") |>
  rename("Cases" = "Incidence") |>
  mutate(
    Date = as.Date(Date, format = "%d/%m/%Y")
  )
# save df
write_csv(flu, file = "data/flu/daily_flu.csv")