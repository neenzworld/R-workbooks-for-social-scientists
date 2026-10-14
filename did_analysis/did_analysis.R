# ==============================================================================
# DIFFERENCE-IN-DIFFERENCES (DiD) & QUASI-EXPERIMENTS
# Subtitle: Evaluating Policy Impacts with Parallel Trends & Event-Study Designs
# Author: Neena Abimbola Albarus
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEPENDENCY MANAGEMENT & AUTO-INSTALLATION
# ------------------------------------------------------------------------------
required_packages <- c(
  "tidyverse", 
  "fixest", 
  "did", 
  "marginaleffects", 
  "modelsummary"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(fixest)
library(did)
library(marginaleffects)
library(modelsummary)

# ------------------------------------------------------------------------------
# 2. DATA PIPELINE & QUASI-EXPERIMENTAL CONSTRUCTION
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

# Define treatment and control cohorts for quasi-experimental evaluation
set.seed(2026)
treated_countries <- c("Ghana", "Kenya", "Rwanda", "Uganda", "Zambia")
control_countries <- c("Nigeria", "Ethiopia", "Tanzania", "Malawi", "Zimbabwe")

did_df <- owid_df %>%
  filter(country %in% c(treated_countries, control_countries), year >= 1995, year <= 2018) %>%
  mutate(
    treated      = if_else(country %in% treated_countries, 1, 0),
    post         = if_else(year >= 2008, 1, 0),
    # Inject policy effect (+2.5 years post-2008 for treated group)
    life_exp_adj = if_else(treated == 1 & post == 1, life_exp + 2.5, life_exp)
  )

print(head(did_df))

# ------------------------------------------------------------------------------
# 3. MODEL ESTIMATION (Classical 2x2 vs. Two-Way Fixed Effects DiD)
# ------------------------------------------------------------------------------
# 1. Classical 2x2 OLS DiD
did_classical <- lm(life_exp_adj ~ treated * post, data = did_df)

# 2. Two-Way Fixed Effects (TWFE) DiD with Country-Clustered Standard Errors
did_twfe <- feols(life_exp_adj ~ treated:post | country + year, data = did_df)

modelsummary(
  list("Classical 2x2 DiD" = did_classical, "TWFE DiD" = did_twfe),
  vcov = ~country,
  stars = TRUE,
  gof_omit = "IC|Log|F"
)

# ------------------------------------------------------------------------------
# 4. DIAGNOSTIC VISUALIZATIONS & EVENT-STUDY
# ------------------------------------------------------------------------------

# --- Figure 1: Testing Parallel Trends Assumption ---
did_summary <- did_df %>%
  group_by(year, treated) %>%
  summarise(mean_life_exp = mean(life_exp_adj), .groups = "drop") %>%
  mutate(Group = if_else(treated == 1, "Treated Nations", "Control Nations"))

fig1 <- ggplot(did_summary, aes(x = year, y = mean_life_exp, color = Group, group = Group)) +
  geom_vline(xintercept = 2007.5, linetype = "dashed", color = "firebrick", size = 1) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  annotate("text", x = 2004.5, y = max(did_summary$mean_life_exp), label = "Policy Implementation (2008)", color = "firebrick", fontface = "bold") +
  scale_color_manual(values = c("Control Nations" = "steelblue", "Treated Nations" = "darkgreen")) +
  theme_minimal() +
  labs(
    title = "Figure 1: Mean Life Expectancy Trajectories Across Groups (1995–2018)",
    subtitle = "Evaluating Pre-Treatment Parallel Trends and Post-Treatment Divergence",
    x = "Year",
    y = "Mean Life Expectancy (Years)"
  )

print(fig1)

# --- Figure 2: Dynamic Event-Study Specification ---
es_model <- feols(
  life_exp_adj ~ i(year, treated, ref = 2007) | country + year,
  data = did_df
)

iplot(
  es_model,
  main = "Figure 2: Dynamic Event-Study Treatment Effect Estimates",
  xlab = "Year Relative to Intervention (2008)",
  ylab = "Estimated Impact on Life Expectancy (Years)",
  zero.line.par = list(col = "gray40", lty = 2)
)

# --- Figure 3: Marginal Effects Pre- vs Post-Policy ---
mfx <- predictions(did_classical, variables = c("treated", "post"))

fig3 <- ggplot(mfx, aes(x = factor(post, labels = c("Pre-Policy", "Post-Policy")), y = estimate, color = factor(treated, labels = c("Control", "Treated")))) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.3), size = 1) +
  scale_color_manual(values = c("Control" = "steelblue", "Treated" = "darkgreen")) +
  theme_minimal() +
  labs(
    title = "Figure 3: Marginal Adjusted Means (Pre vs. Post Policy)",
    x = "Policy Period",
    y = "Predicted Life Expectancy (Years)",
    color = "Group"
  )

print(fig3)

# ------------------------------------------------------------------------------
# 5. LAB EXERCISES
# ------------------------------------------------------------------------------

# --- Exercise 1: Placebo Interventions ---
placebo_df <- did_df %>%
  filter(year < 2008) %>%
  mutate(placebo_post = if_else(year >= 2002, 1, 0))

placebo_mod <- feols(life_exp ~ treated:placebo_post | country + year, data = placebo_df)
print(summary(placebo_mod, vcov = ~country))

# --- Exercise 2: Callaway & Sant'Anna Staggered DiD ---
did_df_cs <- did_df %>%
  mutate(
    first_treat = if_else(treated == 1, 2008, 0),
    country_id  = as.numeric(as.factor(country))
  )

cs_out <- att_gt(
  yname = "life_exp_adj",
  tname = "year",
  idname = "country_id",
  gname = "first_treat",
  data = did_df_cs
)

print(summary(cs_out))