# ==============================================================================
# SPATIAL ECONOMETRICS & GEOGRAPHIC CLUSTERING
# Subtitle: Spatial Autocorrelation, Distance Weights Matrices, SAR & SEM Models
# Author: Neena Abimbola Albarus
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEPENDENCY MANAGEMENT & AUTO-INSTALLATION
# ------------------------------------------------------------------------------
required_packages <- c(
  "tidyverse", 
  "sf", 
  "spdep", 
  "spatialreg", 
  "modelsummary"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(sf)
library(spdep)
library(spatialreg)
library(modelsummary)

# ------------------------------------------------------------------------------
# 2. DATA PIPELINE & SYNTHETIC SPATIAL CONSTRUCT
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

# Construct spatial coordinate grid for cross-sectional observational units
set.seed(2026)
spatial_df <- owid_df %>%
  filter(!is.na(code), !is.na(life_exp), year == 2015) %>%
  mutate(
    longitude = runif(n(), min = -120, max = 120),
    latitude  = runif(n(), min = -40, max = 60),
    gdp_pc    = pmax(500, rnorm(n(), mean = 15000, sd = 9000)),
    log_gdp   = log(gdp_pc)
  ) %>%
  drop_na()

# Convert data frame to simple features (sf) object
sf_data <- st_as_sf(spatial_df, coords = c("longitude", "latitude"), crs = 4326)

print(head(sf_data))

# ------------------------------------------------------------------------------
# 3. SPATIAL WEIGHTS MATRIX & MORAN'S I AUTOCORRELATION TEST
# ------------------------------------------------------------------------------
# Build k-nearest neighbors graph (k = 5)
coords <- st_coordinates(sf_data)
knn_nb <- knearneigh(coords, k = 5)
nb     <- knn2nb(knn_nb)

# Row-standardized spatial weights matrix (W)
spatial_weights <- nb2listw(nb, style = "W")

# Global Moran's I Test for Spatial Autocorrelation in Life Expectancy
moran_test <- moran.test(sf_data$life_exp, spatial_weights)
print(moran_test)

# ------------------------------------------------------------------------------
# 4. SPATIAL ECONOMETRIC MODEL ESTIMATION (SAR & SEM)
# ------------------------------------------------------------------------------
# Standard OLS Baseline Model
mod_ols <- lm(life_exp ~ log_gdp, data = sf_data)

# Spatial Autoregressive (SAR) / Spatial Lag Model (rho * W * Y)
mod_sar <- lagsarlm(life_exp ~ log_gdp, data = sf_data, listw = spatial_weights)

# Spatial Error Model (SEM) (lambda * W * u)
mod_sem <- errorsarlm(life_exp ~ log_gdp, data = sf_data, listw = spatial_weights)

# Model Summary Comparison
modelsummary(
  list(
    "Baseline OLS"     = mod_ols, 
    "Spatial Lag (SAR)" = mod_sar, 
    "Spatial Error (SEM)" = mod_sem
  ),
  stars = TRUE,
  gof_omit = "IC|Log|F"
)

# ------------------------------------------------------------------------------
# 5. DIAGNOSTIC VISUALIZATIONS
# ------------------------------------------------------------------------------

# --- Figure 1: Moran Scatterplot ---
moran.plot(
  sf_data$life_exp, 
  spatial_weights, 
  labels = FALSE, 
  xlab = "Life Expectancy", 
  ylab = "Spatially Lagged Life Expectancy (W * Life Exp)",
  main = "Figure 1: Moran Scatterplot of Spatial Autocorrelation"
)

# --- Figure 2: Geographic Distribution Map ---
fig2 <- ggplot(data = sf_data) +
  geom_sf(aes(color = life_exp, size = log_gdp), alpha = 0.8) +
  scale_color_viridis_c(option = "viridis") +
  theme_minimal() +
  labs(
    title = "Figure 2: Geographic Spatial Distribution & Economic Intensity",
    subtitle = "Points mapped by coordinates; color maps Life Expectancy",
    color = "Life Exp",
    size = "Log GDP"
  )

print(fig2)

# --- Figure 3: OLS Residual Spatial Mapping ---
sf_data$residuals_ols <- residuals(mod_ols)

fig3 <- ggplot(data = sf_data) +
  geom_sf(aes(color = residuals_ols), size = 3, alpha = 0.8) +
  scale_color_gradient2(low = "firebrick", mid = "white", high = "midnightblue", midpoint = 0) +
  theme_minimal() +
  labs(
    title = "Figure 3: Spatial Clustering of Standard OLS Residuals",
    subtitle = "Non-random residual patterns indicate unmodeled spatial dependence",
    color = "OLS Residuals"
  )

print(fig3)