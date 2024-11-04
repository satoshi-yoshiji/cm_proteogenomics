library(tidyverse)
library(magrittr)
library(ggplot2)
library(survival)
library(survminer)
library(data.table)
library(RNOmni)

setwd("path-to-wd")
system('mkdir -p output')

################
# load the data
################
# Replace the following line with your actual data loading code
comb_df4 <- read.csv("path-to-data-file.csv")

################
# Linear regression
################
linear_res <- lm(RNT_COL6A3 ~ BMI + age + factor(sex) + factor(centre) + time_to_olink_processing + factor(olink_batch), data = comb_df4)
summary(linear_res)

# Save results with 95% confidence intervals
linear_sum <- summary(linear_res)
linear_sum_df <- data.frame(linear_sum$coefficients) %>%
  rename(Beta = Estimate, SE = Std..Error, T.value = t.value, P = Pr...t..) %>%
  mutate(Beta_LCI = Beta - 1.96 * SE, Beta_UCI = Beta + 1.96 * SE) %>%
  rownames_to_column(var = 'Variable')
write_tsv(linear_sum_df, file = 'output/Linear_OutcomeCOL6A3_BMI.tsv')

# Plotting with regression line
ggplot(data = comb_df4, aes(x = BMI, y = RNT_COL6A3)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(
    title = "Relationship Between BMI and RNT_COL6A3",
    x = "BMI",
    y = "RNT_COL6A3"
  ) +
  theme_minimal()

################
# Cox regression
################
cox_res <- coxph(Surv(time_to_event, case_control) ~ RNT_COL6A3 + age + factor(sex) + factor(centre) + time_to_olink_processing + factor(olink_batch), data = comb_df4)
cox_summary <- summary(cox_res)

# Save Cox regression results
cox_sum_df <- data.frame(cox_summary$coefficients) %>%
  mutate(HR_LCI = exp(coef - 1.96 * se.coef.), HR_UCI = exp(coef + 1.96 * se.coef.)) %>%
  rownames_to_column(var = 'Variable')
write_tsv(cox_sum_df, file = 'output/Cox_OutcomeCAD_COL6A3.tsv')

# Plotting the baseline survival function
surv_plot <- ggsurvplot(
  survfit(cox_res, data = comb_df4),
  palette = "#2E9FDF",
  fun = 'event',
  risk.table = TRUE,
  ggtheme = theme_minimal(),
  title = "Baseline Survival Function"
)
print(surv_plot)

################
# Survival analysis by quantile
################
comb_df4 <- comb_df4 %>%
  mutate(quantile_group = cut(
    col6a3,
    breaks = quantile(col6a3, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
    include.lowest = TRUE,
    labels = c(1, 2, 3, 4)
  ))

# Survival plot by quantiles
fit_quantile <- survfit(Surv(time_to_event, case_control) ~ quantile_group, data = comb_df4)
surv_plot_quantile <- ggsurvplot(
  fit_quantile,
  data = comb_df4,
  fun = "event",
  risk.table = TRUE,
  ggtheme = theme_light(),
  palette = c('#3C5488FF', '#00A087FF', '#F39B7FFF', '#DC0000FF'),
  title = "Cumulative Incidence of CAD by Quantile Groups",
  xlab = 'Years'
)
print(surv_plot_quantile)
