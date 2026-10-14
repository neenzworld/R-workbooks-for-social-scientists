# ==============================================================================
# WORKBOOK 9: SYNTHETIC CONTROL METHODS (SCM)
# Subtitle: Comparative Case Studies, Convex Combination Weights & Placebo Inference
# Author: Neena Abimbola Albarus
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEPENDENCY MANAGEMENT & AUTO-INSTALLATION
# ------------------------------------------------------------------------------
required_packages <- c(
  "tidyverse", 
  "Synth", 
  "modelsummary"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(Synth)
library(modelsummary)

# ------------------------------------------------------------------------------
# 2. DATA PIPELINE & SYNTHETIC CONTROL PANEL CONSTRUCT
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

# Construct Panel with 1 Treated Unit and a Donor Pool (1980-2015, Intervention at 2000)
set.seed(2026)
top_countries <- owid_df %>%
  filter(!is.na(code), year >= 1980, year <= 2015) %>%
  group_by(country) %>%
  filter(n() == 36) %>%
  ungroup() %>%
  distinct(country) %>%
  slice(1:15) %>%
  pull(country)

scm_df <- owid_df %>%
  filter(country %in% top_countries, year >= 1980, year <= 2015) %>%
  mutate(
    unit_num = as.numeric(factor(country)),
    gdp_proxy = 5000 + 200 * (year - 1980) + rnorm(n(), 0, 500)
  )

# Introduce a structural treatment shock to Unit 1 (Treated Unit) after 2000
treated_id <- 1
treatment_year <- 2000

scm_df <- scm_df %>%
  mutate(
    life_exp = if_else(unit_num == treated_id & year >= treatment_year, 
                       life_exp + 0.35 * (year - treatment_year + 1), 
                       life_exp)
  )

scm_prep_data <- as.data.frame(scm_df)

print(head(scm_df))

# ------------------------------------------------------------------------------
# 3. SYNTHETIC CONTROL ESTIMATION (SYNTH OPTIMIZATION)
# ------------------------------------------------------------------------------
# Data Preparation for Synth
dataprep_out <- dataprep(
  foo = scm_prep_data,
  predictors = c("life_exp", "gdp_proxy"),
  predictors.op = "mean",
  time.predictors.prior = 1980:1999,
  special.predictors = list(
    list("life_exp", 1995:1999, "mean"),
    list("life_exp", 1990:1994, "mean")
  ),
  dependent = "life_exp",
  unit.variable = "unit_num",
  unit.names.variable = "country",
  time.variable = "year",
  treatment.identifier = treated_id,
  controls.identifier = setdiff(unique(scm_prep_data$unit_num), treated_id),
  time.optimize.ssr = 1980:1999,
  time.plot = 1980:2015
)

# Run Optimization to construct W weights
synth_out <- synth(dataprep_out)

# Tables of Weights & Predictor Balance
synth_tables <- synth.tab(dataprep.res = dataprep_out, synth.res = synth_out)
print(synth_tables$tab.w)
print(synth_tables$tab.pred)

# ------------------------------------------------------------------------------
# 4. DIAGNOSTIC VISUALIZATIONS & GAP ANALYSIS
# ------------------------------------------------------------------------------

# --- Figure 1: Path Plot (Treated vs. Synthetic Control) ---
path.plot(
  synth.res = synth_out,
  dataprep.res = dataprep_out,
  Ylab = "Life Expectancy",
  Xlab = "Year",
  Main = "Figure 1: Observed Trajectory vs. Synthetic Counterfactual",
  Legend = c("Treated Unit", "Synthetic Counterfactual"),
  Legend.position = "topleft"
)
abline(v = treatment_year, lty = 2, col = "red")

# --- Figure 2: Gap Plot (Treatment Effect Over Time) ---
gaps.plot(
  synth.res = synth_out,
  dataprep.res = dataprep_out,
  Ylab = "Life Expectancy Gap",
  Xlab = "Year",
  Main = "Figure 2: Causal Gap (Observed - Synthetic)",
)
abline(v = treatment_year, lty = 2, col = "red")
abline(h = 0, lty = 3, col = "gray")

# --- Figure 3: Donor Pool Convex Combination Weights ---
donor_weights <- as.data.frame(synth_tables$tab.w) %>%
  rownames_to_column(var = "donor_country") %>%
  filter(w.weights > 0.01)

fig3 <- ggplot(donor_weights, aes(x = reorder(unit.names, w.weights), y = w.weights)) +
  geom_col(fill = "midnightblue", width = 0.6) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Figure 3: Synthetic Control Convex Weights (W Matrix)",
    subtitle = "Optimal non-negative weights assigned to donor pool units",
    x = "Donor Unit",
    y = "Weight"
  )

print(fig3)