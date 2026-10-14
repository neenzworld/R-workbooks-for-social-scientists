# ==============================================================================
# CAUSAL MEDIATION & STRUCTURAL EQUATION MODELING (SEM)
# Subtitle: Disentangling Direct, Indirect, and Latent Structural Mechanisms
# Author: Neena Abimbola Albarus
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEPENDENCY MANAGEMENT & AUTO-INSTALLATION
# ------------------------------------------------------------------------------
required_packages <- c(
  "tidyverse", 
  "mediation", 
  "lavaan", 
  "semPlot", 
  "modelsummary"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(mediation)
library(lavaan)
library(semPlot)
library(modelsummary)

# ------------------------------------------------------------------------------
# 2. DATA PIPELINE & SYNTHETIC MECH-DATA CONSTRUCT
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

set.seed(2026)
sem_df <- owid_df %>%
  filter(!is.na(code), !is.na(life_exp), year == 2015) %>%
  mutate(
    # Treatment (X)
    gdp_pc       = pmax(500, rnorm(n(), mean = 12000, sd = 8000)),
    log_gdp      = log(gdp_pc),
    
    # Mediator (M)
    health_exp   = 50 + 0.45 * log_gdp + rnorm(n(), mean = 0, sd = 0.5),
    
    # Outcome (Y)
    life_exp_sim = 45 + 1.2 * log_gdp + 2.8 * health_exp + rnorm(n(), mean = 0, sd = 2),
    
    # Latent indicators for SEM
    gov_eff      = 0.7 * log_gdp + rnorm(n(), 0, 1),
    corrupt_r    = 0.6 * log_gdp + rnorm(n(), 0, 1),
    sanitation   = 0.8 * health_exp + rnorm(n(), 0, 1)
  ) %>%
  drop_na()

print(head(sem_df))

# ------------------------------------------------------------------------------
# 3. CAUSAL MEDIATION ANALYSIS
# ------------------------------------------------------------------------------
mod_m <- lm(health_exp ~ log_gdp, data = sem_df)
mod_y <- lm(life_exp_sim ~ log_gdp + health_exp, data = sem_df)

modelsummary(
  list("Mediator Model (M)" = mod_m, "Outcome Model (Y)" = mod_y),
  stars = TRUE,
  gof_omit = "IC|Log|F"
)

results_med <- mediate(
  mod_m, 
  mod_y, 
  treat = "log_gdp", 
  mediator = "health_exp",
  robustSE = TRUE,
  sims = 500
)

summary(results_med)

# ------------------------------------------------------------------------------
# 4. STRUCTURAL EQUATION MODELING (SEM)
# ------------------------------------------------------------------------------
sem_spec <- '
  # Latent Variable Definitions
  INST =~ gov_eff + corrupt_r
  INFRA =~ health_exp + sanitation

  # Structural Regressions
  INFRA ~ INST + log_gdp
  life_exp_sim ~ INST + INFRA + log_gdp
'

fit_sem <- sem(sem_spec, data = sem_df)
summary(fit_sem, fit.measures = TRUE, standardized = TRUE)

# ------------------------------------------------------------------------------
# 5. DIAGNOSTIC VISUALIZATIONS
# ------------------------------------------------------------------------------
plot(results_med, main = "Figure 1: Average Causal Mediation Effects (ACME vs ADE)")

semPaths(
  fit_sem,
  what = "std",
  layout = "tree",
  edge.label.cex = 1.2,
  style = "lisrel",
  main = "Figure 2: Structural Path Diagram & Latent Measurement Model"
)

fig3 <- sem_df %>%
  ggplot(aes(x = health_exp, y = life_exp_sim, color = log_gdp)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_smooth(method = "lm", color = "firebrick", se = TRUE) +
  scale_color_viridis_c(option = "magma") +
  theme_minimal() +
  labs(
    title = "Figure 3: Intermediary Linkage Between Mediator and Outcome",
    subtitle = "Color gradient maps treatment intensity (Log GDP per capita)",
    x = "Healthcare Expenditure Proxy (Mediator M)",
    y = "Life Expectancy Outcome (Y)",
    color = "Log GDP (X)"
  )

print(fig3)