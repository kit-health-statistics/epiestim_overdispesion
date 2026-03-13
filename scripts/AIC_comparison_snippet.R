# For Ebola, we have trouble fitting the NegBin-L model for
# 19. and 26. April 2015. In these 2 estimation windows, the counts are
# underdispersed. We replace these 2 NegBin-L AIC values by the Poisson AIC
# values. The problematic dates `dates_underdisp` are defined in the main
# script `Rt_estimation_all`.
ind_underdisp <- which(
  seq(params$ebola$start_date, params$ebola$end_date, by = 7) %in%
    dates_underdisp
) - params$ebola$window_width  # Shift by the first estimation window
# Check that we really found the problematic indices. This can be verified by
# looking at the dispersion parameter estimate that is in this case very low.
stopifnot(all(results$ebola$disp$nbin_L[ind_underdisp] < 1e-10))
# Replace the weird AIC values
results$ebola$AIC$nbin_L[ind_underdisp] <- results$ebola$AIC$pois[ind_underdisp]
# Extract the aic values for all estimation windows
df_aic <- vector("list", 3)
names(df_aic) <- names(params)
for (disease in names(results)) {
  df_aic[[disease]] <- data.frame(
    aic = unlist(results[[disease]]$AIC),
    model = rep(
      c("Poiss", "NegBin-Q", "NegBin-L"),
      each = length(results[[disease]]$AIC[[1]])
    ),
    pathogen = disease
  ) |>
    mutate(
      # Does the the negative binomial model have a lower AIC than Poisson?
      lower_aic_than_poiss = rep(results[[disease]]$AIC$pois, 3) > aic,
      # Calculate the likelihood ratio test statistic from the AIC values
      lr_test_stat = rep(results[[disease]]$AIC$pois, 3) - aic + 2
    )
}
df_aic <- bind_rows(df_aic) |>
  mutate(
    pathogen = factor(pathogen, levels = c("flu", "covid", "ebola")),
    model = factor(model, levels = c("Poiss", "NegBin-L", "NegBin-Q")),
    # NA values, when we compare Poisson with Poisson
    lower_aic_than_poiss = ifelse(model == "Poiss", NA, lower_aic_than_poiss),
    lr_test_stat = ifelse(model == "Poiss", NA, lr_test_stat)
  )
# Summarise the AIC values: mean AIC across all windows, percentage of times
# when AIC was lower than in the Poisson case and the number of likelihood ratio
# test rejections on the significance level 5 %. In this settings, when we test
# for the presence of overdispersion, the parameter value lies on the boundary
# of the parameteric space, rendering the likelihood ratio test conservative.
# For this reason, we have to correct the test. If the parameter value under the
# null hypothesis lies inside the parametric space, the test statistic
# asymptotically follows a chi^2(1) distribution. However, since the parameter
# under the null lies on the boundary of the space (zero), the test statistic
# is distributed according to an equal mixture of the chi^2(1) and a point mass
# at zero.
table_aic <- df_aic |>
  group_by(pathogen, model) |>
  summarise(
    mean_aic = mean(aic),
    lower_aic_percent = sum(lower_aic_than_poiss) / n() * 100,
    # For a non-corrected LR test we would calculate the p-value
    # as `1 - pchisq(lr_test_stat, 1)`. However we have to account for the half
    # of the probability mass located at zero.
    lr_test_reject_percent = sum(
      1 - 0.5 * (
        ifelse(lr_test_stat < 0, 0, 1) + pchisq(lr_test_stat, 1)
      ) < 0.05
    ) / n() * 100,
    .groups = "drop"
  ) |>
  # Widen the table, so that we have one row per pathogen and three columns
  # per model.
  # We will have 2 columns (`lower_aic_percent` and `lr_test_reject_percent`)
  # for the Poisson model containing NAs only. We will keep them and use them
  # in the final tableas as dummy columns for the purpose of padding.
  pivot_wider(
    names_from = "model",
    values_from = c("mean_aic", "lower_aic_percent", "lr_test_reject_percent")
  ) |>
  arrange(pathogen) |>
  select(-"pathogen") |>
  as.matrix() |>
  unname()

# Column order should be like this:
# mean AIC Poisson
# dummy column
# mean AIC NegBin-L
# percentage of the NegBin-L AIC lower than in the Poisson case
# percentage of LT test rejections for Poisson vs. NegBin-L
# dummy column
# mean AIC NegBin-Q
# percentage of the NegBin-Q AIC lower than in the Poisson case
# percentage of LT test rejections for Poisson vs. NegBin-Q
column_order <- c(1, 4, 2, 5, 8, 7, 3, 6, 9)
table_aic <- table_aic[, column_order] |>
  round(digits = 2) |>
  format(nsmall = 2)
# Replace the NA strings with empty strings
table_aic[grepl("NA", table_aic)] <- ""
rownames(table_aic) <- c("Pandemic Influenza", "COVID-19", "Ebola")
# Print the table
aic_latex <- xtable(table_aic)

# AIC Boxplots ------------------------------------------

# For Flu and COVID-19, we have to zoom in to see the negative binomial boxlots
# in comparison to Poisson
magn_from <- list(
  flu = c(xmin = 1.7, xmax = 3.3, ymin = 110, ymax = 140),
  covid = c(xmin = 1.7, xmax = 3.3, ymin = 150, ymax = 300),
  ebola = c(xmin = NA, xmax = NA, ymin = NA, ymax = NA)
)
magn_to <- list(
  flu = c(xmin = 1.6, xmax = 3.4, ymin = 2000, ymax = 12000),
  covid = c(xmin = 1.6, xmax = 3.4, ymin = 5000, ymax = 40000),
  ebola = c(xmin = NA, xmax = NA, ymin = NA, ymax = NA)
)

# Create the individual AIC boxplots
p_aic_box <- vector("list", 3)
names(p_aic_box) <- names(params)
for (disease in names(params)) {
  p_aic_box[[disease]] <- ggplot(
    data = filter(df_aic, pathogen == disease),
    aes(x = model, y = aic, fill = model)
  ) +
    geom_boxplot(width = 0.6, staplewidth = 1) +
    scale_fill_manual(values = model_colors[-2]) +
    labs(x = "Model", y = "AIC", fill = "Model") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 15)
    )
    if (disease != "ebola") {
      p_aic_box[[disease]] <- p_aic_box[[disease]] +
        ggmagnify::geom_magnify(
          from = magn_from[[disease]],
          to = magn_to[[disease]],
          axes = "y"
        )
    }
}

# Extract the legend
p_aic_box_legend <- wrap_elements(ggpubr::get_legend(p_aic_box$flu))
# Compose the panels. Here we use `annotations` created for other plots in the
# main script.
aic_boxplot <- with(
  p_aic_box,
  wrap_elements(flu & theme(legend.position = "none") & annotations$flu) +
    wrap_elements(covid & theme(legend.position = "none") & annotations$covid) +
    wrap_elements(ebola & theme(legend.position = "none") & annotations$ebola) +
    p_aic_box_legend
) +
  plot_layout(widths = c(3, 3, 3, 1))
ggsave(
  "figure/aic_boxplot.pdf",
  aic_boxplot,
  width = 14,
  height = 6.5
)
ggsave(
  "figure/aic_boxplot.png",
  aic_boxplot,
  width = 14,
  height = 6.5,
  dpi = 400
)
