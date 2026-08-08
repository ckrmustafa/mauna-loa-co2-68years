# The Pulse and the Acceleration: Mauna Loa CO₂ Record (1958-2026) - Companion R Shiny Application

Single-file R Shiny application reproducing every figure and table of the companion paper: a unified 68-year analysis of trend regimes (PELT changepoint detection on the year-over-year growth rate), seasonal cycle evolution (STL-based amplitude and phase analysis with Newey-West inference) and out-of-sample forecasting (SARIMA, harmonic, XGB-hybrid, ensemble, seasonal naive) of the NOAA GML Mauna Loa monthly CO₂ record, March 1958 to July 2026 (821 monthly means).

## Contents

- `app.R` : complete application (R >= 4.2; packages: shiny, readxl, dplyr, tidyr, ggplot2, forecast, changepoint, xgboost, patchwork, lubridate, zoo, sandwich, lmtest)
- `NOAA_Mauna_Loa_Aylik_CO2_1958_2026.xlsx` : NOAA GML monthly means (Lan et al., 2025; https://doi.org/10.15138/wkgj-f215), column names in Turkish

## Run

Place `app.R` and the xlsx file in the same folder, then:

```r
shiny::runApp("app.R")
```

All figures can be downloaded from within the app as 600 DPI PNG files, together with backtest metrics (CSV, with a choice of test window) and an amplitude + phase analysis summary (TXT).

## Citation

Çakır, M., Ural, G. N., Yılmaz, M., & Oral, O. (2026). Mauna Loa CO₂ 1958-2026 companion Shiny application (v1.0.3). Zenodo. https://doi.org/10.5281/zenodo.21829988

Data: Lan, X., Petron, G., Baugh, K., et al. (2025). Atmospheric Carbon Dioxide Dry Air Mole Fractions from the NOAA GML Global Greenhouse Gas Reference Network, Version 2026-07-17. https://doi.org/10.15138/wkgj-f215 (accessed 6 August 2026)

The xlsx file reproduces the NOAA GML trends file `co2_mm_mlo.txt` as downloaded on 6 August 2026 (snapshot last updated by NOAA on 5 August 2026); July 2026 (429.12 ppm) is the most recent month, and all 821 monthly values are carried over unchanged.
