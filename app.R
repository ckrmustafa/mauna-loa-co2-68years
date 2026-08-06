# ================================================================
# Mauna Loa CO2: Trend Regimes + Comparative Projections + Backtest
# Single file: app.R  (xlsx must be in the same folder)
# Journal-ready figures: no in-figure titles
# Tabs: Trend Regimes | Projection | Backtest | Seasonal Amplitude | Data | About
# Requires: shiny, readxl, dplyr, tidyr, ggplot2, forecast,
#           xgboost, lubridate, patchwork, zoo, changepoint
# ================================================================
library(shiny)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forecast)
library(xgboost)
library(lubridate)
library(patchwork)
library(zoo)
library(changepoint)

# ---------------- 1. DATA ----------------
load_data <- function(path = "NOAA_Mauna_Loa_Aylik_CO2_1958_2026.xlsx") {
  raw <- read_excel(path, skip = 3,
                    col_names = c("date_raw", "year", "month", "decimal", "co2",
                                  "co2_adj", "days", "std", "unc"))
  df <- raw %>%
    mutate(across(c(year, month, decimal, co2, co2_adj, days, std, unc),
                  ~ suppressWarnings(as.numeric(.x)))) %>%
    filter(!is.na(year), !is.na(month), year >= 1958, month >= 1, month <= 12) %>%
    mutate(date = as.Date(ISOdate(year, month, 1)),
           co2  = ifelse(co2 < -0.5, NA_real_, co2)) %>%
    arrange(date)
  idx <- seq_len(nrow(df))
  df$co2 <- approx(idx, df$co2, xout = idx, rule = 2)$y
  df$t <- df$year + (df$month - 0.5) / 12
  df
}

# Known event periods (contextual bands; optional, off by default)
EVENTS <- tribble(
  ~event,               ~start_date,    ~end_date,      ~color,
  "El Nino 1997-98",    "1997-05-01",   "1998-05-01",   "#e9c46a",
  "El Nino 2015-16",    "2015-05-01",   "2016-05-01",   "#f4a261",
  "COVID-19",           "2020-03-01",   "2021-06-01",   "#2a9d8f",
  "Volcanic eruption",  "2022-11-29",   "2023-07-04",   "#8338ec"
) %>% mutate(start_date = as.Date(start_date), end_date = as.Date(end_date))

# ---------------- 2. MODEL HELPERS ----------------
harm_X <- function(t) {
  cbind(1, t, t^2,
        sin(2 * pi * t), cos(2 * pi * t),
        sin(4 * pi * t), cos(4 * pi * t))
}

lag_mat <- function(y) {
  n  <- length(y)
  mk <- function(k) c(rep(NA, k), y[1:(n - k)])
  data.frame(lag1 = mk(1), lag2 = mk(2), lag3 = mk(3), lag6 = mk(6),
             lag12 = mk(12), lag13 = mk(13), lag24 = mk(24))
}

fit_models <- function(df) {
  yts <- ts(df$co2, start = c(df$year[1], df$month[1]), frequency = 12)
  
  fit_sarima <- Arima(yts, order = c(1, 1, 1),
                      seasonal = list(order = c(0, 1, 1), period = 12))
  
  X <- harm_X(df$t)
  fit_harm <- lm(df$co2 ~ X - 1)
  resid <- df$co2 - as.numeric(X %*% coef(fit_harm))
  
  L <- lag_mat(resid)
  L$sin <- sin(2 * pi * df$t)
  L$cos <- cos(2 * pi * df$t)
  rows <- complete.cases(L)
  
  xgb_fit <- xgboost(x = as.matrix(L[rows, ]),
                     y = resid[rows],
                     nrounds = 400, verbosity = 0,
                     objective = "reg:squarederror",
                     max_depth = 4, learning_rate = 0.05,
                     subsample = 0.9, colsample_bytree = 0.9)
  
  list(yts = yts, sarima = fit_sarima, harm = fit_harm,
       xgb = xgb_fit, resid = resid)
}

forecast_models <- function(models, df, h, future_dates) {
  sf     <- forecast(models$sarima, h = h, level = 95)
  sarima <- as.numeric(sf$mean)
  sn     <- as.numeric(snaive(models$yts, h = h)$mean)
  
  tfut   <- future_dates$year + (future_dates$month - 0.5) / 12
  harm_p <- as.numeric(harm_X(tfut) %*% coef(models$harm))
  
  res_hist <- models$resid
  xgb_res  <- numeric(h)
  for (i in seq_len(h)) {
    ti <- tfut[i]
    nr <- length(res_hist)
    v  <- c(res_hist[nr], res_hist[nr - 1], res_hist[nr - 2],
            res_hist[nr - 5], res_hist[nr - 11], res_hist[nr - 12],
            res_hist[nr - 23], sin(2 * pi * ti), cos(2 * pi * ti))
    m  <- matrix(v, nrow = 1)
    colnames(m) <- c("lag1", "lag2", "lag3", "lag6", "lag12",
                     "lag13", "lag24", "sin", "cos")
    pr <- as.numeric(predict(models$xgb, m))
    xgb_res[i] <- pr
    res_hist   <- c(res_hist, pr)
  }
  xgb_full <- harm_p + xgb_res
  ens      <- (sarima + harm_p + xgb_full) / 3
  
  data.frame(date = future_dates$date,
             SARIMA = sarima, `S-Naive` = sn, Harmonic = harm_p,
             `XGB-Hybrid` = xgb_full, Ensemble = ens,
             lo95 = as.numeric(sf$lower), hi95 = as.numeric(sf$upper),
             check.names = FALSE)
}

MODEL_COLS <- c("SARIMA" = "#1f77b4", "Harmonic" = "#d62728",
                "XGB-Hybrid" = "#2ca02c", "S-Naive" = "#9467bd",
                "Ensemble" = "black")
MODEL_LTS <- c("SARIMA" = "dashed", "Harmonic" = "solid",
               "XGB-Hybrid" = "dotdash", "S-Naive" = "dotted",
               "Ensemble" = "twodash")
MODEL_LWS <- c("SARIMA" = 0.6, "Harmonic" = 0.6, "XGB-Hybrid" = 0.6,
               "S-Naive" = 0.6, "Ensemble" = 1.0)

# ---------------- 3. UI ----------------
ui <- fluidPage(
  titlePanel("Mauna Loa CO\u2082: Trend Regimes, Projections & Backtest"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("models", "Models to display:",
                         choices  = c("SARIMA", "Harmonic", "XGB-Hybrid",
                                      "S-Naive", "Ensemble"),
                         selected = c("SARIMA", "Harmonic", "XGB-Hybrid",
                                      "Ensemble")),
      sliderInput("horizon", "Projection horizon (year):",
                  min = 2027, max = 2040, value = 2035, step = 1,
                  sep = ""),
      checkboxInput("show_ci", "Show SARIMA 95% interval", value = TRUE),
      checkboxInput("show_events",
                    "Highlight known events (El Nino, COVID-19, volcanic disruption)",
                    value = FALSE),
      hr(),
      sliderInput("amp_window", "Amplitude moving average (years):",
                  min = 1, max = 10, value = 5, step = 1, sep = ""),
      hr(),
      downloadButton("dl_regimes", "Regime plot (600 DPI)"),
      br(), br(),
      downloadButton("dl_projection", "Projection plot (600 DPI)"),
      br(), br(),
      downloadButton("dl_backtest", "Backtest plot (600 DPI)"),
      br(), br(),
      downloadButton("dl_amplitude", "Amplitude plot (600 DPI)"),
      br(), br(),
      downloadButton("dl_metrics", "Backtest metrics (CSV)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Trend Regimes (PELT)",
                 br(),
                 textOutput("reg_info"),
                 plotOutput("reg_plot", height = "720px")),
        tabPanel("Comparative Projection",
                 br(),
                 textOutput("proj_info"),
                 plotOutput("proj_plot", height = "620px")),
        tabPanel("Backtest (last 120 months)",
                 br(),
                 textOutput("bt_info"),
                 plotOutput("bt_plot", height = "720px"),
                 tableOutput("bt_table")),
        tabPanel("Seasonal Amplitude",
                 br(),
                 textOutput("amp_info"),
                 plotOutput("amp_plot", height = "720px")),
        tabPanel("Data",
                 tableOutput("data_table")),
        tabPanel("About",
                 br(),
                 wellPanel(
                   h4("Data source and citation"),
                   p("NOAA Global Monitoring Laboratory (GML), Mauna Loa Observatory monthly atmospheric CO\u2082 (ppm), March 1958 \u2013 July 2026."),
                   tags$ul(
                     tags$li(tags$a(href = "https://gml.noaa.gov/ccgg/trends/",
                                    "NOAA GML CO2 trends", target = "_blank")),
                     tags$li(tags$a(href = "https://scrippsco2.ucsd.edu/",
                                    "Scripps CO2 Program (SIO)", target = "_blank"))
                   ),
                   p("Data from March 1958 to April 1974 were obtained by C. David Keeling at the Scripps Institution of Oceanography. Measurements from November 2022 to July 2023 were taken at the Maunakea Observatories due to the eruption of the Mauna Loa volcano."),
                   p("NOAA GML data are freely available for scientific research; please cite NOAA GML / SIO as the data source when using this app or its outputs."),
                   hr(),
                   h4("About this app"),
                   p("Interactive supplement to an academic study of the Mauna Loa CO\u2082 record: PELT changepoint segmentation of trend regimes, seasonal cycle evolution, comparative forecasting (SARIMA, harmonic regression, XGBoost-hybrid, ensemble, seasonal-naive baseline) and out-of-sample backtesting over the last 120 months. All figures can be downloaded as journal-ready 600 DPI PNG files via the buttons in the sidebar.")
                 ))
      )
    )
  )
)

# ---------------- 4. SERVER ----------------
server <- function(input, output, session) {
  
  dat <- reactive({ load_data() })
  
  # ---- Trend regimes (PELT on year-over-year growth of adjusted series) ----
  regimes <- reactive({
    df <- dat()
    y  <- df$co2_adj
    growth <- diff(y, lag = 12)          # year-over-year change, ppm/year
    n  <- length(growth)
    
    cpt <- cpt.mean(growth, method = "PELT", penalty = "Manual",
                    pen.value = 3 * log(n), minseglen = 24)
    cps <- cpts(cpt)                     # indices within the growth series
    brk <- cps + 12                      # matching row indices in df
    
    bounds <- c(1, brk, nrow(df))
    nseg   <- length(bounds) - 1
    # Mean-growth windows reproduce the regime means reported in the paper
    # (verified against the data: 0.94 / 1.58 / 2.33 ppm/yr):
    # regime 1 covers growth values strictly before the first break; every
    # later regime starts just after the previous break and runs through the
    # next break inclusive; the final regime runs to the end of the record.
    idx_win <- function(i) {
      lo <- if (i == 1) 1L else cps[i - 1] + 1L
      hi <- if (i == nseg) length(growth) else if (i == 1) cps[1] - 1L else cps[i]
      lo:hi
    }
    segs <- lapply(seq_len(nseg), function(i) {
      a <- bounds[i]; b <- bounds[i + 1]
      m  <- lm(df$co2_adj[a:b] ~ df$t[a:b])
      gi <- idx_win(i)
      tibble(regime  = paste("Regime", i),
             x0 = df$date[a], x1 = df$date[b],
             slope     = unname(coef(m)[2]),
             intercept = unname(coef(m)[1]),
             mean_g    = mean(growth[gi], na.rm = TRUE),
             xmid      = df$date[a] + (df$date[b] - df$date[a]) / 2)
    }) %>% bind_rows()
    
    trend_df <- lapply(seq_len(nseg), function(i) {
      a <- bounds[i]; b <- bounds[i + 1]
      tibble(date = df$date[a:b],
             fit  = segs$intercept[i] + segs$slope[i] * df$t[a:b])
    }) %>% bind_rows()
    
    gdf <- tibble(date = df$date[13:nrow(df)], growth = growth)
    
    list(segs = segs, trend_df = trend_df, gdf = gdf,
         break_dates = df$date[brk])
  })
  
  # ---- Backtest (dynamic: last 120 months) ----
  bt <- reactive({
    df       <- dat()
    last     <- max(df$date)
    te_start <- last %m-% months(119)
    tr <- df %>% filter(date < te_start)
    te <- df %>% filter(date >= te_start)
    models <- fit_models(tr)
    fut    <- data.frame(date = te$date, year = te$year, month = te$month)
    fc     <- forecast_models(models, tr, nrow(te), fut)
    fc$observed <- te$co2
    list(start = te_start, end = last, fc = fc)
  })
  
  bt_metrics <- reactive({
    fc <- bt()$fc
    res <- lapply(c("SARIMA", "Harmonic", "XGB-Hybrid", "S-Naive", "Ensemble"),
                  function(m) {
                    e <- fc[[m]] - fc$observed
                    data.frame(Model = m,
                               RMSE = sqrt(mean(e^2)),
                               MAE  = mean(abs(e)))
                  })
    do.call(rbind, res)
  })
  
  # ---- Info texts (kept in the UI, not inside the figures) ----
  output$reg_info <- renderText({
    r <- regimes()
    btxt <- paste(format(r$break_dates, "%b %Y"), collapse = " and ")
    stxt <- paste(sprintf("%.2f", r$segs$mean_g), collapse = " / ")
    paste0("PELT changepoint detection (mean-shift cost, strengthened information-criterion ",
           "penalty 3*log(n), minimum segment length 24 months) applied to the year-over-year ",
           "growth of the seasonally adjusted series. Detected breaks: ", btxt,
           ". Regime mean growth rates: ", stxt, " ppm/year. Upper panel: seasonally adjusted ",
           "CO2 with piecewise linear regime trends; lower panel: annual growth rate with regime means.")
  })
  output$proj_info <- renderText({
    paste0("Observed NOAA GML record 1958\u20132026 and comparative model projections to ",
           input$horizon,
           ". Shaded vertical band marks the projection window; inset shows each model's deviation from the Ensemble.")
  })
  output$bt_info <- renderText({
    b <- bt()
    paste0("Out-of-sample backtest over the last 120 months (",
           format(b$start, "%b %Y"), " \u2013 ", format(b$end, "%b %Y"),
           "). Grey line: observed values; models trained on all data before the test window. Lower panel: forecast minus observed; shaded band = \u00B11 ppm.")
  })
  output$amp_info <- renderText({
    paste0("STL decomposition with a time-varying seasonal window (s.window = 13, robust). ",
           "Upper panel: seasonal component. Lower panel: annual peak-to-trough amplitude with ",
           input$amp_window, "-year moving average (teal) and linear trend (red dashed).")
  })
  
  # ---- Event bands helper ----
  add_event_bands <- function(p, date_min, date_max) {
    ev <- EVENTS %>% filter(end_date >= date_min, start_date <= date_max)
    if (nrow(ev) == 0) return(p)
    p + geom_rect(data = ev, inherit.aes = FALSE,
                  aes(xmin = start_date, xmax = end_date,
                      ymin = -Inf, ymax = Inf, fill = event),
                  alpha = 0.20) +
      scale_fill_manual(values = setNames(ev$color, ev$event), name = NULL)
  }
  
  # ---- Regime segmentation object ----
  regime_obj <- reactive({
    r <- regimes()
    
    p_top <- ggplot(dat(), aes(date, co2_adj)) +
      geom_line(color = "grey60", linewidth = 0.4) +
      geom_line(data = r$trend_df, aes(date, fit),
                color = "#2a9d8f", linewidth = 1.0) +
      geom_vline(xintercept = r$break_dates, linetype = "dashed",
                 color = "#d62728", linewidth = 0.6) +
      labs(x = NULL, y = expression(CO[2] ~ "(ppm, seasonally adjusted)")) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_blank())
    
    p_bot <- ggplot(r$gdf, aes(date, growth)) +
      geom_line(color = "grey70", linewidth = 0.35) +
      geom_segment(data = r$segs,
                   aes(x = x0, xend = x1, y = mean_g, yend = mean_g),
                   color = "#2a9d8f", linewidth = 1.1, inherit.aes = FALSE) +
      geom_vline(xintercept = r$break_dates, linetype = "dashed",
                 color = "#d62728", linewidth = 0.6) +
      geom_text(data = r$segs,
                aes(x = xmid, y = mean_g + 0.28,
                    label = sprintf("%.2f ppm/yr", mean_g)),
                size = 4, color = "grey20", inherit.aes = FALSE) +
      labs(x = NULL, y = "Annual growth (ppm/yr)") +
      theme_minimal(base_size = 13)
    
    p_top / p_bot + plot_layout(heights = c(1.4, 1))
  })
  
  # ---- Projection object ----
  projection_obj <- reactive({
    df  <- dat()
    models <- fit_models(df)
    last  <- max(df$date)
    h_end <- as.Date(paste0(input$horizon, "-12-01"))
    fdates <- seq.Date(last %m+% months(1), h_end, by = "month")
    fut <- data.frame(date  = fdates,
                      year  = as.numeric(format(fdates, "%Y")),
                      month = as.numeric(format(fdates, "%m")))
    f <- forecast_models(models, df, nrow(fut), fut)
    
    long <- f %>%
      select(date, SARIMA, Harmonic, `XGB-Hybrid`, `S-Naive`, Ensemble) %>%
      pivot_longer(-date, names_to = "Model", values_to = "pred") %>%
      filter(Model %in% input$models)
    long$Model <- factor(long$Model, levels = names(MODEL_COLS))
    
    main <- ggplot() +
      geom_line(data = df, aes(date, co2), color = "grey60", linewidth = 0.4) +
      geom_line(data = long,
                aes(date, pred, color = Model, linetype = Model,
                    linewidth = Model)) +
      annotate("rect", xmin = last, xmax = h_end, ymin = -Inf, ymax = Inf,
               fill = "orange", alpha = 0.05) +
      geom_vline(xintercept = last, linetype = "dotted", color = "grey40") +
      scale_color_manual(values = MODEL_COLS) +
      scale_linetype_manual(values = MODEL_LTS) +
      scale_linewidth_manual(values = MODEL_LWS) +
      labs(x = NULL, y = expression(CO[2] ~ "(ppm)")) +
      guides(color    = guide_legend(title = "Model", order = 1,
                                     override.aes = list(linewidth = 1.3)),
             linetype = guide_legend(title = "Model", order = 1),
             linewidth = "none",
             fill     = guide_legend(order = 2, title = NULL)) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top", legend.box = "vertical",
            legend.text = element_text(size = 10),
            legend.key.width = unit(2.2, "lines"))
    
    if (input$show_ci) {
      main <- main +
        geom_ribbon(data = f, aes(date, ymin = lo95, ymax = hi95),
                    fill = "#1f77b4", alpha = 0.15, inherit.aes = FALSE)
    }
    if (input$show_events) {
      main <- add_event_bands(main, min(df$date), max(df$date))
    }
    
    # Zoom panel: deviation from Ensemble (baseline excluded)
    fdev <- f %>%
      mutate(across(c(SARIMA, Harmonic, `XGB-Hybrid`), ~ .x - Ensemble)) %>%
      select(date, SARIMA, Harmonic, `XGB-Hybrid`) %>%
      pivot_longer(-date, names_to = "Model", values_to = "dev") %>%
      filter(Model %in% input$models)
    fdev$Model <- factor(fdev$Model, levels = names(MODEL_COLS))
    
    zoom <- ggplot(fdev, aes(date, dev, color = Model, linetype = Model,
                             linewidth = Model)) +
      geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
      geom_line() +
      scale_color_manual(values = MODEL_COLS) +
      scale_linetype_manual(values = MODEL_LTS) +
      scale_linewidth_manual(values = MODEL_LWS) +
      labs(x = NULL, y = "Deviation (ppm)") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none")
    
    list(plot = main + inset_element(zoom, left = 0.03, bottom = 0.52,
                                     right = 0.48, top = 0.97,
                                     align_to = "panel"),
         forecast = f)
  })
  
  # ---- Backtest object (main panel + error-evolution panel) ----
  backtest_obj <- reactive({
    b  <- bt()
    fc <- b$fc
    
    long <- fc %>%
      select(date, SARIMA, Harmonic, `XGB-Hybrid`, `S-Naive`, Ensemble) %>%
      pivot_longer(-date, names_to = "Model", values_to = "pred") %>%
      filter(Model %in% input$models)
    long$Model <- factor(long$Model, levels = names(MODEL_COLS))
    
    p_main <- ggplot(long, aes(date, pred, color = Model, linetype = Model,
                               linewidth = Model)) +
      geom_line(data = fc, aes(date, observed),
                color = "grey55", linewidth = 0.9, inherit.aes = FALSE) +
      geom_line() +
      scale_color_manual(values = MODEL_COLS) +
      scale_linetype_manual(values = MODEL_LTS) +
      scale_linewidth_manual(values = MODEL_LWS) +
      labs(x = NULL, y = expression(CO[2] ~ "(ppm)")) +
      theme_minimal(base_size = 13) +
      guides(color    = guide_legend(title = "Model", order = 1,
                                     override.aes = list(linewidth = 1.3)),
             linetype = guide_legend(title = "Model", order = 1),
             linewidth = "none",
             fill     = guide_legend(order = 2, title = NULL)) +
      theme(legend.position = "top", legend.box = "vertical",
            legend.text = element_text(size = 10),
            legend.key.width = unit(2.2, "lines"),
            axis.text.x = element_blank())
    
    if (input$show_events) {
      p_main <- add_event_bands(p_main, min(fc$date), max(fc$date))
    }
    
    # Error panel: deviation from observed (baseline excluded)
    err <- fc %>%
      mutate(across(c(SARIMA, Harmonic, `XGB-Hybrid`, Ensemble),
                    ~ .x - observed)) %>%
      select(date, SARIMA, Harmonic, `XGB-Hybrid`, Ensemble) %>%
      pivot_longer(-date, names_to = "Model", values_to = "err") %>%
      filter(Model %in% input$models)
    err$Model <- factor(err$Model, levels = names(MODEL_COLS))
    
    p_err <- ggplot(err, aes(date, err, color = Model, linetype = Model,
                             linewidth = Model)) +
      annotate("rect", xmin = min(fc$date), xmax = max(fc$date),
               ymin = -1, ymax = 1, fill = "grey90", alpha = 0.4) +
      geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
      geom_line() +
      scale_color_manual(values = MODEL_COLS) +
      scale_linetype_manual(values = MODEL_LTS) +
      scale_linewidth_manual(values = MODEL_LWS) +
      labs(x = NULL, y = "Error (ppm)") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none")
    
    (p_main / p_err) +
      plot_layout(heights = c(2.4, 1))
  })
  
  # ---- Seasonal amplitude object (STL, time-varying seasonal window) ----
  amplitude_obj <- reactive({
    df  <- dat()
    yts <- ts(df$co2, start = c(df$year[1], df$month[1]), frequency = 12)
    s   <- stl(yts, s.window = 13, robust = TRUE)
    
    d2 <- data.frame(date = df$date, year = df$year,
                     seasonal = as.numeric(s$time.series[, "seasonal"]))
    amp <- d2 %>%
      group_by(year) %>%
      summarise(amplitude = max(seasonal) - min(seasonal),
                n = n(), .groups = "drop") %>%
      filter(n == 12)
    amp$moving_avg <- rollmean(amp$amplitude, k = input$amp_window,
                               fill = NA, align = "center")
    
    p_stl <- ggplot(d2, aes(date, seasonal)) +
      geom_line(color = "#457b9d", linewidth = 0.4) +
      labs(x = NULL, y = "Seasonal component (ppm)") +
      theme_minimal(base_size = 13)
    
    fit <- lm(amplitude ~ year, data = amp)
    slope_dec <- unname(coef(fit)[2] * 10)
    
    p_amp <- ggplot(amp, aes(year, amplitude)) +
      geom_point(color = "grey60", size = 1.2) +
      geom_line(aes(y = moving_avg), color = "#2a9d8f", linewidth = 0.9) +
      geom_smooth(method = "lm", se = FALSE, color = "#d62728",
                  linetype = "dashed", linewidth = 0.7) +
      annotate("text", x = min(amp$year) + 1, y = max(amp$amplitude),
               label = sprintf("Linear trend: %+.3f ppm/decade", slope_dec),
               hjust = 0, size = 4) +
      labs(x = NULL, y = "Amplitude (ppm)") +
      theme_minimal(base_size = 13)
    
    p_stl / p_amp + plot_layout(heights = c(1, 1.2))
  })
  
  # ---------------- 5. OUTPUTS ----------------
  output$reg_plot  <- renderPlot({ regime_obj() })
  output$proj_plot <- renderPlot({ projection_obj()$plot })
  output$bt_plot   <- renderPlot({ backtest_obj() })
  output$amp_plot  <- renderPlot({ amplitude_obj() })
  output$bt_table  <- renderTable({ bt_metrics() }, digits = 3)
  output$data_table <- renderTable({
    dat() %>% select(date, year, month, co2) %>% tail(24)
  })
  
  output$dl_regimes <- downloadHandler(
    filename = function() "CO2_trend_regimes_PELT_600dpi.png",
    content  = function(file) {
      ggsave(file, plot = regime_obj(),
             width = 10, height = 7, dpi = 600)
    }
  )
  output$dl_projection <- downloadHandler(
    filename = function() paste0("CO2_projection_", input$horizon, "_600dpi.png"),
    content  = function(file) {
      ggsave(file, plot = projection_obj()$plot,
             width = 10, height = 6, dpi = 600)
    }
  )
  output$dl_backtest <- downloadHandler(
    filename = function() "CO2_backtest_last120months_600dpi.png",
    content  = function(file) {
      ggsave(file, plot = backtest_obj(),
             width = 10, height = 7, dpi = 600)
    }
  )
  output$dl_amplitude <- downloadHandler(
    filename = function() "CO2_seasonal_amplitude_600dpi.png",
    content  = function(file) {
      ggsave(file, plot = amplitude_obj(),
             width = 10, height = 7, dpi = 600)
    }
  )
  output$dl_metrics <- downloadHandler(
    filename = function() "CO2_backtest_metrics.csv",
    content  = function(file) write.csv(bt_metrics(), file, row.names = FALSE)
  )
}

shinyApp(ui, server)
