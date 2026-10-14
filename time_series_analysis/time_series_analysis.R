# ==============================================================================
# WORKBOOK 4: TIME-SERIES ANALYSIS & LAGGED EFFECTS
# Author: Neena Abimbola Albarus
# ==============================================================================

# 1. DEPENDENCY MANAGEMENT
required_packages <- c("tidyverse", "tseries", "forecast", "modelsummary")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(tidyverse)
library(tseries)
library(forecast)
library(modelsummary)

# 2. DATA PIPELINE
possible_paths <- c(
  "../../data/life_expectancy_data.csv",
  "data/life_expectancy_data.csv",
  "../data/life_expectancy_data.csv",
  "/Users/yawdgyal/Desktop/R for github/Spatial Econometrics & Geographic Clustering/data/life_expectancy_data.csv"
)

data_path <- possible_paths[file.exists(possible_paths)][1]

if (is.na(data_path)) {
  stop("Could not locate 'life_expectancy_data.csv'.")
}

owid_raw <- read_csv(data_path)

if ("Entity" %in% names(owid_raw)) {
  owid_df <- owid_raw %>%
    rename(country = Entity, code = Code, year = Year, life_exp = `Life expectancy`)
} else {
  owid_df <- owid_raw %>% rename_with(tolower)
}

# Plain data frame (strip ts attributes)
ts_df <- owid_df %>%
  filter(!is.na(life_exp), year >= 1950) %>%
  group_by(year) %>%
  summarise(mean_life_exp = mean(life_exp, na.rm = TRUE), .groups = "drop") %>%
  arrange(year) %>%
  mutate(
    mean_life_exp = as.numeric(mean_life_exp),
    diff_life_exp = mean_life_exp - lag(mean_life_exp, 1),
    lag1_life_exp = lag(mean_life_exp, 1),
    lag2_life_exp = lag(mean_life_exp, 2)
  ) %>%
  drop_na()

# Separate ts object strictly for ACF/PACF and forecasting functions
life_exp_ts <- ts(ts_df$mean_life_exp, start = min(ts_df$year), frequency = 1)

# 3. STATIONARITY DIAGNOSTICS & MODEL ESTIMATION
print(adf.test(life_exp_ts))
print(adf.test(diff(life_exp_ts)))

# Fit using plain numeric data frame columns
mod_ardl_lm <- lm(mean_life_exp ~ lag1_life_exp + lag2_life_exp, data = ts_df)

modelsummary(
  list("Autoregressive Lag Model (ARDL)" = mod_ardl_lm),
  stars = TRUE,
  gof_omit = "IC|Log|F"
)

# 4. DIAGNOSTIC VISUALIZATIONS
# Figure 1
fig1 <- ggplot(ts_df, aes(x = year, y = mean_life_exp)) +
  geom_line(color = "midnightblue", size = 1.2) +
  geom_point(color = "midnightblue", size = 2) +
  theme_minimal() +
  labs(
    title = "Figure 1: Global Mean Life Expectancy Time-Series (1950–Present)",
    subtitle = "Evaluating Non-Stationarity and Long-Run Stochastic Secular Trends",
    x = "Year", y = "Global Mean Life Expectancy (Years)"
  )
print(fig1)

# Figure 2
par(mfrow = c(1, 2))
acf(life_exp_ts, main = "Figure 2A: Autocorrelation (ACF)")
pacf(life_exp_ts, main = "Figure 2B: Partial ACF (PACF)")
par(mfrow = c(1, 1))

# Figure 3
fit_arima <- auto.arima(life_exp_ts)
fcast     <- forecast(fit_arima, h = 10)

fig3 <- autoplot(fcast) +
  theme_minimal() +
  labs(
    title = "Figure 3: 10-Year ARIMA Forecast with 80% & 95% Confidence Intervals",
    x = "Year", y = "Projected Life Expectancy"
  )
print(fig3)