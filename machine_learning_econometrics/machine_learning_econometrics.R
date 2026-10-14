# ==============================================================================
# WORKBOOK 10: MACHINE LEARNING & HIGH-DIMENSIONAL ECONOMETRICS
# Subtitle: Regularized Regression (Lasso/Ridge) & Double/Debiased Machine Learning (DML)
# Author: Neena Abimbola Albarus
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEPENDENCY MANAGEMENT & AUTO-INSTALLATION
# ------------------------------------------------------------------------------
required_packages <- c(
  "tidyverse", 
  "glmnet", 
  "DoubleML", 
  "mlr3", 
  "mlr3learners", 
  "modelsummary"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(glmnet)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(modelsummary)

# ------------------------------------------------------------------------------
# 2. DATA PIPELINE & HIGH-DIMENSIONAL CONSTRUCT
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

# Construct high-dimensional covariate space (p = 20 controls + treatment)
set.seed(2026)
ml_df <- owid_df %>%
  filter(!is.na(code), !is.na(life_exp), year == 2015) %>%
  mutate(
    # Continuous Treatment (D): Healthcare expenditure per capita proxy
    treatment_d = pmax(100, rnorm(n(), mean = 2500, sd = 1200))
  )

# Simulate 20 high-dimensional confounding covariates (X1 to X20)
p_covariates <- 20
n_obs <- nrow(ml_df)
X_mat <- matrix(rnorm(n_obs * p_covariates), nrow = n_obs, ncol = p_covariates)
colnames(X_mat) <- paste0("X_", 1:p_covariates)

ml_df <- bind_cols(ml_df, as_tibble(X_mat))

# Construct Outcome (Y) driven by sparse signals (X_1, X_2, X_3 confound Y)
ml_df <- ml_df %>%
  mutate(
    life_exp_sim = 50 + 0.003 * treatment_d + 1.5 * X_1 - 2.0 * X_2 + 1.1 * X_3 + rnorm(n_obs, 0, 3)
  ) %>%
  drop_na()

print(head(ml_df))

# ------------------------------------------------------------------------------
# 3. HIGH-DIMENSIONAL REGULARIZATION & DOUBLE MACHINE LEARNING (DML)
# ------------------------------------------------------------------------------
# Part A: LASSO Feature Selection (Explicitly namespace dplyr::select)
X_features_df <- ml_df %>% dplyr::select(starts_with("X_"), treatment_d)
X_features    <- as.matrix(X_features_df)
Y_vector      <- ml_df$life_exp_sim

cv_lasso <- cv.glmnet(X_features, Y_vector, alpha = 1)
print(paste("Optimal Lasso Lambda (min):", cv_lasso$lambda.min))

# Part B: Double/Debiased Machine Learning (Partially Linear Model)
dml_data <- DoubleMLData$new(
  data = as.data.frame(ml_df),
  y_col = "life_exp_sim",
  d_cols = "treatment_d",
  x_cols = paste0("X_", 1:p_covariates)
)

# Initialize Lasso learners for treatment (m(X)) and outcome (g(X)) nuisance models
learner_l <- lrn("regr.glmnet", lambda = cv_lasso$lambda.min)
learner_m <- lrn("regr.glmnet", lambda = cv_lasso$lambda.min)

dml_plr <- DoubleMLPLR$new(
  data = dml_data,
  ml_l = learner_l,
  ml_m = learner_m,
  n_folds = 5
)

# Estimate orthogonalized causal coefficient
dml_plr$fit()
print(dml_plr$summary())

# ------------------------------------------------------------------------------
# 4. DIAGNOSTIC VISUALIZATIONS
# ------------------------------------------------------------------------------

# --- Figure 1: LASSO Coefficient Regularization Path ---
par(mfrow = c(1, 2))
plot(cv_lasso, main = "Figure 1A: CV MSE vs. Lambda Penalty")
plot(cv_lasso$glmnet.fit, xvar = "lambda", main = "Figure 1B: Regularization Path")
par(mfrow = c(1, 1))

# --- Figure 2: Double ML Orthogonal Residuals (Y-tilde vs D-tilde) ---
res_y <- Y_vector - predict(cv_lasso, newx = X_features, s = "lambda.min")

X_only <- as.matrix(ml_df %>% dplyr::select(starts_with("X_")))
cv_lasso_d <- cv.glmnet(X_only, ml_df$treatment_d, alpha = 1)
res_d <- ml_df$treatment_d - predict(cv_lasso_d, newx = X_only, s = "lambda.min")

df_res <- tibble(res_y = as.numeric(res_y), res_d = as.numeric(res_d))

fig2 <- ggplot(df_res, aes(x = res_d, y = res_y)) +
  geom_point(color = "midnightblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "firebrick", size = 1.2) +
  theme_minimal() +
  labs(
    title = "Figure 2: Double ML Orthogonal Residual Plot (Frisch-Waugh-Lovell)",
    subtitle = "Nuisance components partialled out via cross-validated machine learning",
    x = "Residualized Treatment Variable (D - E[D|X])",
    y = "Residualized Outcome Variable (Y - E[Y|X])"
  )

print(fig2)

# --- Figure 3: High-Dimensional Feature Importance (Lasso Coefficients) ---
lasso_coefs <- as.matrix(coef(cv_lasso, s = "lambda.min")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Feature") %>%
  rename(Coefficient = 2) %>%
  filter(Feature != "(Intercept)", Coefficient != 0)

fig3 <- ggplot(lasso_coefs, aes(x = reorder(Feature, abs(Coefficient)), y = Coefficient)) +
  geom_col(fill = "steelblue", width = 0.6) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Figure 3: Non-Zero Variable Coefficients Selected by LASSO",
    subtitle = "Sparse regularization isolates true confounding signals",
    x = "Selected Covariate",
    y = "Estimated Penalty-Adjusted Coefficient"
  )

print(fig3)