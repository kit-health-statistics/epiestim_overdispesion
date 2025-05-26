library(tidyverse)
# file taken from https://github.com/fintzij/stemr/tree/master/data
load(file=here::here("data", "ebola", "ebola.rda"))
ebola_df <- as.data.frame(ebola)
# convert row names to new column
ebola_df <- ebola_df %>%
  rownames_to_column(var = "week") 
# save df
write_csv(ebola_df, file = here::here("data", "ebola", "ebola_df.csv"))