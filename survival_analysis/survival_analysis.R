# ==============================================================================
# WORKBOOK 3: SURVIVAL & EVENT HISTORY ANALYSIS IN SOCIAL SCIENCE
# Subtitle: Non-Parametric Kaplan-Meier & Cox Proportional Hazards Models
# Author: Neena Abimbola Albarus
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEPENDENCY MANAGEMENT & AUTO-INSTALLATION
# ------------------------------------------------------------------------------
required_packages <- c(
  "tidyverse", 
  "survival", 
  "survminer", 
  "modelsummary"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(survival)
library(survminer)
library(modelsummary)

# ------------------------------------------------------------------------------
# 2. DATA PIPELINE & EVENT HISTORY CONSTRUCT
# ------------------------------------------------------------------------------
# Dynamic Path Resolution
possible_paths <- c(
  "../../data/life_expectancy_data.csv",
  "data/life_expectancy_data.csv",
  "../data/life_expectancy_data.csv",
  "/Users/yawdgyal/Desktop/R for github/Spatial Econometrics & Geographic Clustering/data/life_expectancy_data.csv"
)

data_path <- possible_paths[file.exists(possible_paths)][1]

if (is.na(data_path)) {
  stop("Could not locate 'life_expectancy_data.csv'. Please check your dataset location.")
}

owid_raw <- read_csv(data_path)

# Standardize variable names
if ("Entity" %in% names(owid_raw)) {
  owid_df <- owid_raw %>%
    rename(
      country  = Entity,
      code     = Code,
      year     = Year,
      life_exp = `Life expectancy`
    )
} else {
  owid_df <- owid_raw %>% rename_with(tolower)
}

# Construct event history framework (Time until reaching Life Expectancy >= 70)
survival_df <- owid_df %>%
  filter(!is.na(code), !is.na(life_exp)) %>%
  group_by(country) %>%
  arrange(year) %>%
  summarise(
    start_year    = min(year),
    first_70_year = min(year[life_exp >= 70], na.rm = TRUE),
    max_year      = max(year),
    baseline_exp  = first(life_exp),
    .groups       = "drop"
  ) %>%
  mutate(
    # Event indicator: 1 if country reached 70+ life expectancy, 0 if censored
    event = if_else(is.finite(first_70_year), 1, 0),
    # Time to event or right-censoring duration
    time = if_else(event == 1, first_70_year - start_year, max_year - start_year)
  ) %>%
  filter(time > 0)

print(head(survival_df))

# ------------------------------------------------------------------------------
# 3. NON-PARAMETRIC KAPLAN-MEIER ESTIMATION
# ------------------------------------------------------------------------------
# Fit baseline Kaplan-Meier survival curve
km_fit <- survfit(Surv(time, event) ~ 1, data = survival_df)

print(summary(km_fit, times = c(10, 20, 30, 40, 50)))

# ------------------------------------------------------------------------------
# 4. COX PROPORTIONAL HAZARDS SEMI-PARAMETRIC REGRESSION
# ------------------------------------------------------------------------------
# Fit Cox Proportional Hazards Model evaluating baseline expectancy impact
cox_mod <- coxph(Surv(time, event) ~ baseline_exp, data = survival_df)

# Output summary table
modelsummary(
  list("Cox PH Model" = cox_mod),
  exponentiate = TRUE, # Display Hazard Ratios (HR)
  stars = TRUE,
  gof_omit = "IC|Log|F"
)

# Test Proportional Hazards Assumption
ph_test <- cox.zph(cox_mod)
print(ph_test)

# ------------------------------------------------------------------------------
# 5. DIAGNOSTIC VISUALIZATIONS
# ------------------------------------------------------------------------------

# --- Figure 1: Kaplan-Meier Survival Curve ---
km_plot <- ggsurvplot(
  km_fit,
  data = survival_df,
  conf.int = TRUE,
  palette = "midnightblue",
  title = "Figure 1: Kaplan-Meier Estimate of Time to Reach Life Expectancy of 70 Years",
  xlab = "Years from Baseline Observation",
  ylab = "Survival Probability (Not Yet Reached Target)",
  ggtheme = theme_minimal()
)

print(km_plot)

# --- Figure 2: Cox Proportional Hazards Cumulative Hazard ---
cum_haz_plot <- ggsurvplot(
  km_fit,
  data = survival_df,
  fun = "cumhaz",
  palette = "firebrick",
  title = "Figure 2: Cumulative Hazard of Reaching Target Life Expectancy",
  xlab = "Years from Baseline Observation",
  ylab = "Cumulative Hazard",
  ggtheme = theme_minimal()
)

print(cum_haz_plot)

# --- Figure 3: Testing Proportional Hazards Assumption Residuals ---
fig3 <- ggcoxzph(ph_test)
print(fig3)