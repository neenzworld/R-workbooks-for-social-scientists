# ==============================================================================
# WORKBOOK 7: INSTRUMENTAL VARIABLES & TWO-STAGE LEAST SQUARES (2SLS)
# Subtitle: Endogeneity, Instrument Validity, Weak Instruments and Hausman Tests
# Author: Neena Abimbola Albarus
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEPENDENCY MANAGEMENT & AUTO-INSTALLATION
# ------------------------------------------------------------------------------
required_packages <- c(
  "tidyverse", 
  "AER", 
  "modelsummary", 
  "ggpubr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(AER)
library(modelsummary)
library(ggpubr)

# ------------------------------------------------------------------------------
# 2. DATA PIPELINE & ENDOGENOUS SYSTEM CONSTRUCT
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

# Construct endogenous system: Health Spending (X) correlated with Unobserved Health Policy (U)
set.seed(2026)
iv_df <- owid_df %>%
  filter(!is.na(code), !is.na(life_exp), year == 2015) %>%
  mutate(
    # Unobserved institutional quality / baseline population health (U)
    unobserved_u = rnorm(n(), mean = 0, sd = 1.5),
    
    # Exogenous Instrument (Z): Geographic distance / historical policy shocks
    instrument_z = rnorm(n(), mean = 50, sd = 15),
    
    # Endogenous Regressor (X): Healthcare spending driven by Z and U
    health_exp = 100 + 1.8 * instrument_z + 2.5 * unobserved_u + rnorm(n(), mean = 0, sd = 5),
    
    # True Outcome (Y): Life expectancy driven by health_exp and unobserved_u
    life_exp_sim = 50 + 0.15 * health_exp + 2.0 * unobserved_u + rnorm(n(), mean = 0, sd = 3)
  ) %>%
  drop_na()

print(head(iv_df))

# ------------------------------------------------------------------------------
# 3. TWO-STAGE LEAST SQUARES (2SLS) ESTIMATION & DIAGNOSTICS
# ------------------------------------------------------------------------------
# Biased Baseline OLS Regression
mod_ols <- lm(life_exp_sim ~ health_exp, data = iv_df)

# Stage 1: Regress Endogenous Variable on Instrument
mod_stage1 <- lm(health_exp ~ instrument_z, data = iv_df)

# Stage 2: Structural IV / 2SLS Regression using ivreg()
mod_2sls <- ivreg(life_exp_sim ~ health_exp | instrument_z, data = iv_df)

# Model Comparison Table
modelsummary(
  list(
    "Biased OLS"      = mod_ols, 
    "First Stage"     = mod_stage1, 
    "2SLS Instrumental" = mod_2sls
  ),
  stars = TRUE,
  gof_omit = "IC|Log|F"
)

# Formal IV Diagnostic Tests (Weak Instruments, Wu-Hausman Endogeneity, Sargan Overid)
iv_diagnostics <- summary(mod_2sls, diagnostics = TRUE)
print(iv_diagnostics$diagnostics)

# ------------------------------------------------------------------------------
# 4. DIAGNOSTIC VISUALIZATIONS
# ------------------------------------------------------------------------------

# --- Figure 1: First-Stage Relevance (Z -> X) ---
fig1 <- ggplot(iv_df, aes(x = instrument_z, y = health_exp)) +
  geom_point(color = "midnightblue", alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "firebrick", se = TRUE) +
  theme_minimal() +
  labs(
    title = "Figure 1: First-Stage Relevance Plot (Instrument Z vs. Endogenous X)",
    subtitle = "Strong linear relationship confirms instrument relevance (First-Stage F > 10)",
    x = "Exogenous Instrument (Z)",
    y = "Endogenous Healthcare Expenditure (X)"
  )

print(fig1)

# --- Figure 2: Comparison of OLS Bias vs. 2SLS Causal Estimates ---
fig2 <- ggplot(iv_df, aes(x = health_exp, y = life_exp_sim)) +
  geom_point(alpha = 0.5, color = "gray40") +
  geom_smooth(method = "lm", aes(color = "Biased OLS Slope"), se = FALSE, size = 1.2) +
  geom_abline(
    aes(
      intercept = coef(mod_2sls)[1], 
      slope = coef(mod_2sls)[2], 
      color = "2SLS Causal Slope"
    ), 
    size = 1.2, 
    linetype = "dashed"
  ) +
  scale_color_manual(
    name = "Estimator", 
    values = c("Biased OLS Slope" = "firebrick", "2SLS Causal Slope" = "midnightblue")
  ) +
  theme_minimal() +
  labs(
    title = "Figure 2: Empirical Comparison of Biased OLS vs. 2SLS Structural Estimates",
    subtitle = "OLS is upward-biased due to unobserved confounders; 2SLS recovers causal effect",
    x = "Endogenous Healthcare Expenditure (X)",
    y = "Simulated Life Expectancy Outcome (Y)"
  )

print(fig2)

# --- Figure 3: Structural Residuals vs Instrument (Exogeneity Check) ---
iv_df$residuals_2sls <- residuals(mod_2sls)

fig3 <- ggplot(iv_df, aes(x = instrument_z, y = residuals_2sls)) +
  geom_point(color = "darkgreen", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "Figure 3: 2SLS Residuals vs. Instrument (Exogeneity Diagnostics)",
    subtitle = "Lack of pattern supports orthgonality between instrument and structural error",
    x = "Exogenous Instrument (Z)",
    y = "2SLS Structural Residuals"
  )

print(fig3)