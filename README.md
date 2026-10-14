# Applied Econometrics and Causal Inference in R

**A ten part workbook series on causal identification, panel data, spatial econometrics, and machine learning for the social sciences**

Author: **Neena Abimbola Albarus**

---

## Overview

This repository contains ten self contained Quarto workbooks that walk through the modern empirical social science toolkit from first principles. Each workbook combines executable R code, methodological exposition, diagnostic visualizations, and lab exercises on real or synthetic datasets drawn from the Our World in Data life expectancy panel. The series moves from foundational techniques such as time series decomposition and panel fixed effects, through the causal inference identification strategies that dominate contemporary applied work, and finally into machine learning methods that are increasingly central to high dimensional empirical research.

The workbooks are designed for students, researchers, and analysts who want to understand not just how to run a model, but why the model is appropriate for a given research question and what can go wrong when its assumptions are violated. Every workbook includes Common Pitfalls callouts that flag the most frequent mistakes in applied practice, and every figure ends with a Key Takeaway that tells the reader what diagnostic signal to look for.

---

## The Ten Workbooks

| # | Title | Core Methods |
|---|---|---|
| **1** | Panel Data Analysis and Fixed Effects Models | Pooled OLS, country fixed effects, two way fixed effects, Hausman test |
| **2** | Difference in Differences and Quasi Experiments | Classical 2x2 DiD, TWFE DiD, dynamic event studies, Callaway and Sant'Anna |
| **3** | Survival and Event History Analysis | Kaplan Meier curves, Cox proportional hazards, Schoenfeld residuals |
| **4** | Time Series Analysis and Lagged Effects | ADF stationarity tests, autoregressive distributed lag models, ARIMA forecasting |
| **5** | Causal Mediation and Structural Equation Modeling | Baron and Kenny decomposition, ACME and ADE, latent variable SEM with lavaan |
| **6** | Spatial Econometrics and Geographic Clustering | Moran's I, spatial weight matrices, SAR and SEM models, LISA clusters |
| **7** | Instrumental Variables and Two Stage Least Squares | Endogeneity, instrument validity, weak instrument diagnostics, Wu Hausman test |
| **8** | Regression Discontinuity Design | Sharp RDD, optimal bandwidth selection, McCrary density test, local polynomial fits |
| **9** | Synthetic Control Methods | Convex combination weights, donor pool construction, placebo inference |
| **10** | Machine Learning and High Dimensional Econometrics | LASSO and Ridge regularization, Double Debiased Machine Learning |

---

## Prerequisites

Across the series, you are expected to be comfortable with the following before you begin.

- **R fundamentals**: reading data with `read_csv`, manipulating data frames with `dplyr`, and building layered plots with `ggplot2`
- **Linear regression**: interpreting coefficients, standard errors, and R squared, and understanding the difference between correlation and causation
- **Basic statistics**: hypothesis tests, confidence intervals, and the logic of the sampling distribution
- **Familiarity with the tidyverse**: `%>%` pipes, `mutate()`, `filter()`, `group_by()`, and `summarise()`

You do not need prior experience with any of the specialized packages. Each workbook installs its own dependencies at the top.

---

## Repository Structure

```
.
├── README.md
├── data/
│   └── life_expectancy_data.csv
└── workbooks/
    ├── 01_panel_data_analysis.qmd
    ├── 02_difference_in_differences.qmd
    ├── 03_survival_analysis.qmd
    ├── 04_time_series_analysis.qmd
    ├── 05_mediation_sem.qmd
    ├── 06_spatial_econometrics.qmd
    ├── 07_instrumental_variables.qmd
    ├── 08_regression_discontinuity.qmd
    ├── 09_synthetic_control.qmd
    └── 10_high_dimensional_econometrics.qmd
```

Every workbook expects `life_expectancy_data.csv` to be reachable via one of several relative or absolute fallback paths. The workbooks search for the file in `data/`, `../data/`, and `../../data/` relative to their own location, and fall back to the absolute path on the author's machine. If you move the data file, update the `possible_paths` vector at the top of each workbook's data pipeline chunk.

---

## Data Source

| Source | Description | Link |
|---|---|---|
| **Our World in Data (OWID)** | National life expectancy panel covering most countries from 1950 to the present | https://ourworldindata.org/life-expectancy |

The workbooks use OWID as either a real analysis dataset (Workbooks 1, 3, 4) or as the raw skeleton for a synthetic data generating process with a known ground truth (Workbooks 2, 5 through 10). When a synthetic mechanism is used, it is explicitly labeled in the prose so that no reader mistakes a simulated coefficient for a real finding.

---

## Required R Packages

Each workbook installs its own dependencies at the top, but the full set across the series is:

```r
install.packages(c(
  # Core
  "tidyverse", "modelsummary", "scales",
  # Panel and DiD
  "fixest", "plm", "did", "marginaleffects",
  # Survival
  "survival", "survminer",
  # Time series
  "tseries", "forecast",
  # Mediation and SEM
  "mediation", "lavaan", "semPlot",
  # Spatial
  "sf", "spdep", "spatialreg",
  # Instrumental variables
  "AER",
  # Regression discontinuity
  "rdrobust", "rddensity",
  # Synthetic control
  "Synth",
  # Machine learning and DML
  "glmnet", "DoubleML", "mlr3", "mlr3learners",
  # Plotting helpers
  "ggpubr", "RColorBrewer", "viridis", "rnaturalearth"
))
```

Some packages such as `sf` require system level geospatial libraries (GDAL, GEOS, PROJ). On macOS, install them via Homebrew with `brew install gdal geos proj`. On Ubuntu, use `sudo apt install libgdal-dev libgeos-dev libproj-dev`. The `r2u` and `pak` binary mirrors on Linux provide fast precompiled binaries that avoid lengthy source compilation.

---

## How to Render

Each workbook is a self contained Quarto document. To render a single workbook:

```bash
cd workbooks
quarto render 01_panel_data_analysis.qmd
```

To render every workbook in sequence:

```bash
cd workbooks
for f in *.qmd; do quarto render "$f"; done
```

Or open any `.qmd` file in RStudio and click **Render**. The rendered HTML files will appear alongside the source files unless you configure a separate output directory in a `_quarto.yml` file.

---

## Suggested Learning Paths

The workbooks are numbered in a recommended reading order, but not every path through them is linear. The following tracks are organized by research interest.

### Path A: Causal Inference Foundations

Start with **Workbook 1** to build intuition for panel data and unobserved heterogeneity. Move to **Workbook 2** for difference in differences. Then **Workbook 7** for instrumental variables and **Workbook 8** for regression discontinuity. Finish with **Workbook 9** for synthetic control methods. This is the standard empirical microeconomics sequence and covers the four identification strategies that dominate applied causal work.

### Path B: Longitudinal and Time Series Analysis

Start with **Workbook 4** for time series fundamentals. Then **Workbook 1** for panel data and **Workbook 3** for survival analysis. Finish with **Workbook 9** for synthetic control, which bridges panel and time series logic.

### Path C: Structural and Mechanism Based Approaches

Start with **Workbook 5** for mediation and SEM. Then **Workbook 6** for spatial econometrics if your data has a geographic dimension. Finish with **Workbook 10** for high dimensional machine learning, which is increasingly used alongside structural models in applied work.

### Path D: Machine Learning for Social Science

Start with **Workbook 4** for time series intuition. Then **Workbook 10** for regularized regression and double machine learning. Return to **Workbook 2** and **Workbook 7** to see how machine learning is used inside causal identification strategies.

---

## Common Threads Across the Series

Several methodological principles recur across every workbook and are worth internalizing once.

- **Identification before estimation.** Every method in this series rests on assumptions that cannot be tested directly, such as parallel trends in DiD, exclusion restriction in IV, and continuity in RDD. The workbook prose emphasizes the assumptions alongside the mechanics.
- **Clustered standard errors are the default.** Panel, spatial, and time series data all exhibit serial correlation within units. Clustering at the appropriate level is not optional.
- **Diagnostic visualization beats summary statistics.** A coefficient table without the corresponding figure is half a result. Every workbook pairs each model with a plot that shows what the model is actually doing.
- **Synthetic data with a known ground truth is a powerful teaching device.** Injecting a real treatment effect into simulated data lets you see whether an estimator recovers what was built in. When the estimator gets it wrong, you learn something about the method.
- **Every figure tells a diagnostic story.** The Key Takeaway lines are not decoration. They summarize the specific signal a reader should extract from each plot.

---

## Exercises and Extensions

Most workbooks end with a Lab Exercises section that asks you to extend the core analysis. These exercises cover topics such as placebo interventions, staggered treatment adoption, sensitivity analyses, and alternative weight matrices. They are designed to be done in an R session alongside the workbook, not to have pre written solutions. If you want the exercises to have published solutions, open an issue or submit a pull request.

---

## Further Reading

Each workbook ends with a Further Reading section listing the canonical papers and textbooks for its topic. The core references across the whole series are:

- **Angrist, J. and Pischke, J.S. (2009).** *Mostly Harmless Econometrics*. The standard reference for the causal inference identification strategies covered in Workbooks 2, 7, and 8.
- **Angrist, J. and Pischke, J.S. (2014).** *Mastering Metrics*. A gentler, example driven introduction to the same material.
- **Cunningham, S. (2021).** *Causal Inference: The Mixtape*. A modern, code forward treatment of the identification toolkit.
- **Imbens, G. and Rubin, D. (2015).** *Causal Inference for Statistics, Social, and Biomedical Sciences*. The formal theoretical foundation.
- **Wooldridge, J. (2010).** *Econometric Analysis of Cross Section and Panel Data*. The classic reference for panel methods.
- **Kuhn, M. and Silge, J. (2022).** *Tidy Modeling with R*. The reference for the tidymodels ecosystem used in the machine learning workbook.
- **Anselin, L. (1988).** *Spatial Econometrics: Methods and Models*. The foundational reference for the spatial workbook.

---

## License

For educational use. Please cite Our World in Data and Natural Earth as the original data providers where applicable. The synthetic mechanisms constructed inside several workbooks are original to this project.
