# ==============================================================================
# REGRESSION DISCONTINUITY DESIGN (RDD)
# Subtitle: Sharp/Fuzzy Threshold Dynamics, Optimal Bandwidth Selection & Local Linear Regression
# Author: Neena Abimbola Albarus
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEPENDENCY MANAGEMENT & AUTO-INSTALLATION
# ------------------------------------------------------------------------------
required_packages <- c(
  "tidyverse", 
  "rdrobust", 
  "rddensity", 
  "modelsummary"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(rdrobust)
library(rddensity)
library(modelsummary)

# ------------------------------------------------------------------------------
# 2. DATA PIPELINE & SHARP RDD RUNNING VARIABLE CONSTRUCT
# ------------------------------------------------------------------------------
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

# Construct Sharp RDD framework
# Running Variable (X): Normalized GDP per capita centering policy eligibility at threshold c = 0
# Treatment (D): Program eligibility assigned deterministically when running_x >= 0
set.seed(2026)
cutoff_val <- 10000

rdd_df <- owid_raw %>%
  filter(!is.na(`Life expectancy`), Year == 2015) %>%
  mutate(
    # Synthetic GDP per capita centered at $10,000 threshold
    gdp_pc = pmax(1000, rnorm(n(), mean = 11000, sd = 6000)),
    running_x = gdp_pc - cutoff_val,
    
    # Sharp Treatment Assignment
    treatment_d = if_else(running_x >= 0, 1, 0),
    
    # Outcome (Y): Discontinuity effect (+4.5 years jump at threshold)
    life_exp_sim = 55 + 0.0008 * running_x + 4.5 * treatment_d + rnorm(n(), mean = 0, sd = 2.5)
  ) %>%
  drop_na()

print(head(rdd_df))

# ------------------------------------------------------------------------------
# 3. RDD ESTIMATION (PARAMETRIC OLS vs NON-PARAMETRIC RDROBUST)
# ------------------------------------------------------------------------------
# Linear Parametric OLS with Interaction Slope
mod_ols_rdd <- lm(life_exp_sim ~ running_x * treatment_d, data = rdd_df)

modelsummary(
  list("Parametric RDD (OLS)" = mod_ols_rdd),
  stars = TRUE,
  gof_omit = "IC|Log|F"
)

# Non-Parametric Local Polynomial RDD (Calonico, Cattaneo, Titiunik optimal bandwidth)
rd_fit <- rdrobust(
  y = rdd_df$life_exp_sim, 
  x = rdd_df$running_x, 
  c = 0
)

summary(rd_fit)

# McCrary Density Test for Manipulation / Running Variable Sorting
density_test <- rddensity(X = rdd_df$running_x, c = 0)
summary(density_test)

# ------------------------------------------------------------------------------
# 4. DIAGNOSTIC VISUALIZATIONS
# ------------------------------------------------------------------------------

# --- Figure 1: Sharp Regression Discontinuity Plot ---
fig1 <- rdplot(
  y = rdd_df$life_exp_sim, 
  x = rdd_df$running_x, 
  c = 0,
  title = "Figure 1: Non-Parametric Local Polynomial RDD Plot",
  x.label = "Running Variable: Centered GDP per Capita (c = 0)",
  y.label = "Life Expectancy Outcome (Years)"
)

# --- Figure 2: Binned Scatterplot with Linear Parametric Fit ---
fig2 <- ggplot(rdd_df, aes(x = running_x, y = life_exp_sim, color = factor(treatment_d))) +
  geom_point(alpha = 0.4, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(
    values = c("0" = "firebrick", "1" = "midnightblue"),
    labels = c("Control (Below Cutoff)", "Treated (Above Cutoff)")
  ) +
  theme_minimal() +
  labs(
    title = "Figure 2: Parametric Regression Discontinuity Slope & Jump at Cutoff",
    subtitle = "Discontinuity at x = 0 estimates local average treatment effect (LATE)",
    x = "Centered Running Variable (GDP per Capita - $10k Threshold)",
    y = "Life Expectancy",
    color = "Treatment Status"
  )

print(fig2)

# --- Figure 3: McCrary Running Variable Density Test ---
fig3 <- rdplotdensity(
  rdd = density_test, 
  X = rdd_df$running_x,
  title = "Figure 3: McCrary Density Continuity Test for Sorting Manipulation"
)