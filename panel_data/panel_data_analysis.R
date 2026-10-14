# ==============================================================================
# PANEL DATA ANALYSIS & FIXED EFFECTS MODELS
# Subtitle: Controlling for Unobserved Heterogeneity in Longitudinal Social Science Data
# Author: Neena Abimbola Albarus
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEPENDENCY MANAGEMENT & AUTO-INSTALLATION
# ------------------------------------------------------------------------------
required_packages <- c(
  "tidyverse", 
  "fixest", 
  "plm", 
  "modelsummary"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(fixest)
library(plm)
library(modelsummary)

# ------------------------------------------------------------------------------
# 2. DATA PIPELINE & PANEL STRUCTURING
# ------------------------------------------------------------------------------
possible_paths <- c(
  "../../data/life_expectancy_data.csv",
  "data/life_expectancy_data.csv",
  "../data/life_expectancy_data.csv",
  "/Users/yawdgyal/Desktop/R for github/Spatial Econometrics & Geographic Clustering/data/life_expectancy_data.csv"
)

data_path <- possible_paths[file.exists(possible_paths)][1]

if (is.na(data_path)) {
  stop("Could not locate 'life_expectancy_data.csv'. Please verify your project data path.")
}

owid_raw <- read_csv(data_path)

if ("Entity" %in% names(owid_raw)) {
  owid_panel <- owid_raw %>%
    rename(
      country  = Entity,
      code     = Code,
      year     = Year,
      life_exp = `Life expectancy`
    )
} else {
  owid_panel <- owid_raw %>% rename_with(tolower)
}

panel_df <- owid_panel %>%
  filter(!is.na(code), !is.na(life_exp), year >= 1990) %>%
  group_by(country) %>%
  mutate(
    life_exp_lag = lag(life_exp, 1),
    time_trend   = year - min(year)
  ) %>%
  ungroup() %>%
  drop_na(life_exp, life_exp_lag)

print(head(panel_df))

# ------------------------------------------------------------------------------
# 3. MODEL ESTIMATION (Pooled OLS vs. One-Way & Two-Way Fixed Effects)
# ------------------------------------------------------------------------------
mod_pooled     <- feols(life_exp ~ life_exp_lag, data = panel_df)
mod_country_fe <- feols(life_exp ~ life_exp_lag | country, data = panel_df)
mod_twfe       <- feols(life_exp ~ life_exp_lag | country + year, data = panel_df)

modelsummary(
  list("Pooled OLS" = mod_pooled, "Country FE" = mod_country_fe, "Two-Way FE" = mod_twfe),
  vcov = ~country,
  stars = TRUE,
  gof_omit = "IC|Log|F"
)

# ------------------------------------------------------------------------------
# 4. DIAGNOSTIC & METHODOLOGICAL VISUALIZATIONS
# ------------------------------------------------------------------------------

# --- Figure 1: Within-Country Trajectories ---
sample_countries <- c("Ghana", "Brazil", "Vietnam", "Norway", "United States", "Kenya")

fig1 <- panel_df %>%
  filter(country %in% sample_countries) %>%
  ggplot(aes(x = year, y = life_exp, color = country, group = country)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  labs(
    title = "Figure 1: Within-Country Health Trajectories Over Time (1990–Present)",
    subtitle = "Country Fixed Effects Isolate Within-Unit Slopes from Between-Unit Baseline Levels",
    x = "Year",
    y = "Life Expectancy at Birth (Years)",
    color = "Nation"
  )

print(fig1)

# --- Figure 2: Model Estimate Comparisons ---
coef_df <- tibble(
  Model = c("Pooled OLS", "Country FE", "Two-Way FE"),
  Estimate = c(coef(mod_pooled)["life_exp_lag"], coef(mod_country_fe)["life_exp_lag"], coef(mod_twfe)["life_exp_lag"]),
  SE = c(se(mod_pooled)["life_exp_lag"], se(mod_country_fe)["life_exp_lag"], se(mod_twfe)["life_exp_lag"])
) %>%
  mutate(
    ci_low = Estimate - 1.96 * SE,
    ci_high = Estimate + 1.96 * SE
  )

fig2 <- ggplot(coef_df, aes(x = Model, y = Estimate, color = Model)) +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high), size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "Figure 2: Persistence Parameter Estimate Across Panel Model Specifications",
    subtitle = "Clustered Standard Errors (Country Level) with 95% Confidence Intervals",
    x = "Specification Structure",
    y = "Lagged Life Expectancy Coefficient Estimate"
  ) +
  theme(legend.position = "none")

print(fig2)

# --- Figure 3: Country Fixed Effects Distribution ---
fe_vals <- fixef(mod_country_fe)$country

fe_df <- tibble(
  country = names(fe_vals),
  alpha_i = as.numeric(fe_vals)
)

fig3 <- ggplot(fe_df, aes(x = alpha_i)) +
  geom_histogram(bins = 30, fill = "midnightblue", color = "white", alpha = 0.8) +
  geom_density(aes(y = after_stat(count)), color = "firebrick", size = 1) +
  theme_minimal() +
  labs(
    title = "Figure 3: Distribution of Unobserved Country-Specific Intercepts (Alpha_i)",
    subtitle = "Quantifying Baseline Heterogeneity Absorbed by Country Fixed Effects",
    x = "Estimated Country Fixed Effect (Baseline Offset)",
    y = "Count of Countries"
  )

print(fig3)

# --- Figure 4: Year Fixed Effects ---
year_fe_vals <- fixef(mod_twfe)$year

year_df <- tibble(
  year = as.numeric(names(year_fe_vals)),
  gamma_t = as.numeric(year_fe_vals)
)

fig4 <- ggplot(year_df, aes(x = year, y = gamma_t)) +
  geom_line(color = "darkgreen", size = 1.1) +
  geom_point(color = "darkgreen", size = 2.5) +
  theme_minimal() +
  labs(
    title = "Figure 4: Macroeconomic & Global Shocks Absorbed by Year Fixed Effects (Gamma_t)",
    subtitle = "Common Macro Dynamics Filtered Out in Two-Way Fixed Effects (TWFE)",
    x = "Year",
    y = "Year Fixed Effect Coefficient Offset"
  )

print(fig4)

# ------------------------------------------------------------------------------
# 5. LAB EXERCISES
# ------------------------------------------------------------------------------

# --- Exercise 1: High-Dimensional Fixed Effects ---
mod_ex1 <- feols(life_exp ~ life_exp_lag | country + year, data = panel_df)
print(summary(mod_ex1, vcov = ~country))

# --- Exercise 2: Hausman Specification Test ---
pdata <- pdata.frame(panel_df, index = c("country", "year"))

fe_plm <- plm(life_exp ~ life_exp_lag, data = pdata, model = "within")
re_plm <- plm(life_exp ~ life_exp_lag, data = pdata, model = "random")

print(phtest(fe_plm, re_plm))