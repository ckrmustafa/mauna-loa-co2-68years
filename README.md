# The Pulse and the Acceleration: Mauna Loa CO2 Record (1958-2026) - Companion R Shiny Application

Single-file R Shiny application reproducing every figure and table of the companion paper: a unified 68-year analysis of trend regimes (PELT changepoint detection on the year-over-year growth rate), seasonal cycle evolution (harmonic model and STL decomposition) and out-of-sample forecasting (SARIMA, harmonic, XGB-hybrid, ensemble, seasonal naive) of the NOAA GML Mauna Loa monthly COb record, March 1958 to July 2026 (821 monthly means).

## Contents

- `app.R` : complete application (R >= 4.2; packages: shiny, readxl, dplyr, tidyr, ggplot2, forecast, changepoint, xgboost, patchwork, lubridate)
- `NOAA_Mauna_Loa_Aylik_CO2_1958_2026.xlsx` : NOAA GML monthly means (Lan et al., 2025; https://doi.org/10.15138/wkgj-f215), column names in Turkish

## Run

Place `app.R` and the xlsx file in the same folder, then:

```r
shiny::runApp("app.R")
```

All figures can be downloaded from within the app as 600 DPI PNG files.

## Citation

Cakir, M. (2026). Mauna Loa CO2 1958-2026 companion Shiny application (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

Data: Lan, X., Petron, G., Baugh, K., et al. (2025). Atmospheric Carbon Dioxide Dry Air Mole Fractions from the NOAA GML Global Greenhouse Gas Reference Network, Version 2025-08-15. https://doi.org/10.15138/wkgj-f215
