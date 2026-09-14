# ─────────────────────────────────────────────
# Analysis of the scraped petrol prices
# ─────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readabs)
library(knitr)
library(kableExtra)
bstar_colors <- c("#DE8C59", "#337179","#C1B2A2", "#622128","#333F48" ,"#000000", "#A5DFD3")

scraped.fuel.raw <- read.csv("data/petrol_prices.csv", header = TRUE)

cities <- c("Sydney", "Melbourne", "Brisbane", "Adelaide", "Perth")
cpi_weights<-c(0.2974,0.2848,0.1599,0.0753,0.1294) #weights of automotive fuels across Capitals in order
cpi_weights<-cpi_weights/sum(cpi_weights)
aus_cpi_weight<- 0.032 #hardcode for 2025 weighting pattern of automotive fuel
#excludes dates which the ACCC site was down
excluded_dates <- seq.Date(
  from = as.Date("2026-04-20"),
  to   = as.Date("2026-05-16"),
  by   = "day"
)

# ─────────────────────────────────────────────
# STEP 1: Keep latest forecast_date per date
# ─────────────────────────────────────────────

### This section maks the fuel data be based of the median to stop data errors ####
# Compute median per date per city (across all forecast_dates)
city_medians <- scraped.fuel.raw %>%
mutate(
    date          = as.Date(date),
    forecast_date = as.Date(forecast_date)
  ) %>%
filter(!forecast_date %in% excluded_dates) %>% #
  group_by(date) %>%
  summarise(across(all_of(cities), ~ median(.x, na.rm = TRUE)), .groups = "drop")

# Keep latest forecast_date per date
fuel.data <- scraped.fuel.raw %>%
  mutate(
    date          = as.Date(date),
    forecast_date = as.Date(forecast_date)
  ) %>%
  group_by(date) %>%
  filter(!forecast_date %in% excluded_dates) %>%
  slice_max(forecast_date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(date) %>%
  select(forecast_date, date, point, all_of(cities)) %>%
  # Join medians and swap in where date is > 1 month old
  left_join(city_medians, by = "date", suffix = c("", "_med")) %>%
  mutate(
    across(
      all_of(cities),
      ~ if_else(date < Sys.Date() %m-% months(1), get(paste0(cur_column(), "_med")), .x)
    )
  ) %>%
  select(forecast_date, date, point, all_of(cities))
######


# ─────────────────────────────────────────────
# STEP 2: Forecast rest of current month
# ─────────────────────────────────────────────
last_obs_date  <- max(fuel.data$date)
next_q_end     <-  ceiling_date(last_obs_date, "quarter") %m+% months(3) - days(1)
forecast_dates <- seq(last_obs_date + 1, next_q_end, by = "day")

if (length(forecast_dates) > 0) {
  
  last_row <- fuel.data %>% slice_tail(n = 1)
  
  # ── Flat carry-forward for Sydney, Melbourne, Brisbane, Adelaide ──
  flat_cities <- c("Sydney", "Melbourne", "Brisbane", "Adelaide")
  flat_values <- last_row %>% select(all_of(flat_cities))
  
  # ── Perth: use regression of weekday to find seasonality
  
  
  perth_history <- fuel.data %>%
    select(date, Perth) %>%
    mutate(
      day_of_week = wday(date, label = TRUE, abbr = TRUE),
      mod_7       = as.factor(as.numeric(max(fuel.data$date) - date) %% 7)
    ) %>%
    arrange(date)
  
  perth_lm <- lm(Perth ~ mod_7, data = perth_history)
  
  # Build forecast rows
  perth_last <- last_row$Perth
  perth_fcst <- numeric(length(forecast_dates))
  
  for (i in seq_along(forecast_dates)) {
    # mod_7 for each forecast date, relative to the last observed date
    mod_val <- as.character(-i %% 7)
    
    # Extract the coefficient for this mod_7 level (falls back to 0 if reference level)
    coef_name <- paste0("mod_7", mod_val)
    step      <- ifelse(coef_name %in% names(coef(perth_lm)),
                        coef(perth_lm)[coef_name],
                        0)
    
    perth_fcst[i] <- perth_last + step
  }
  
  forecast_rows <- tibble(
    forecast_date    = as.Date(NA),          # no scrape date for estimates
    date             = forecast_dates,
    point            = NA_integer_,
    Actual_estimate  = "E",
    Sydney           = as.numeric(flat_values$Sydney),
    Melbourne        = as.numeric(flat_values$Melbourne),
    Brisbane         = as.numeric(flat_values$Brisbane),
    Adelaide         = as.numeric(flat_values$Adelaide),
    Perth            = perth_fcst
  )
  
  # Tag existing rows as Actuals and bind
  fuel.data <- fuel.data %>%
    mutate(Actual_estimate = "A") %>%
    bind_rows(forecast_rows) %>%
    arrange(date)
  
} else {
  fuel.data <- fuel.data %>% mutate(Actual_estimate = "A")
  message("No forecast dates needed — last observation is already at month end.")
}

#────────────────────────────────────────────
# STEP 3: Australia weighted average series
# ─────────────────────────────────────────────
fuel.data <- fuel.data %>%
  mutate(
    Australia = as.numeric(as.matrix(select(., all_of(cities))) %*% cpi_weights)) %>%
  relocate(Australia, .after = Perth)

# ── Daily data table: last actual month + first estimate month ──
last_A_month <- fuel.data %>%
  filter(Actual_estimate == "A") %>%
  pull(date) %>%
  max() %>%
  floor_date("month")

first_E_month <- fuel.data %>%
  filter(Actual_estimate == "E") %>%
  pull(date) %>%
  min() %>%
  floor_date("month")

daily_detail <- fuel.data %>%
  mutate(year_month = floor_date(date, "month")) %>%
  filter(year_month %in% c(last_A_month, first_E_month)) %>%
  select(date, Actual_estimate, all_of(cities), Australia) %>%
  mutate(date = format(date, "%d-%b-%y"),
         across(c(all_of(cities), "Australia"), ~ round(.x, 1)))


# ─────────────────────────────────────────────
# STEP 4: Monthly average 
# ─────────────────────────────────────────────

monthly_avg <- fuel.data %>%
  mutate(year_month = floor_date(date, "month")) %>%
  group_by(year_month) %>%
  summarise(
    across(c(all_of(cities), "Australia"), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    Actual_estimate = if_else(all(Actual_estimate == "A"), "A", "E"),
    .groups = "drop"
  )

mom_growth <- monthly_avg %>%
  arrange(year_month) %>%
  mutate(
    across(c(all_of(cities), "Australia"), ~ .x / lag(.x) - 1),
    MoM_cpi_bps = round(aus_cpi_weight * Australia * 10^4, digits = 0)
  ) %>%
  na.omit()

# ─────────────────────────────────────────────
# Step 5: Generate charts
# ─────────────────────────────────────────────

plot_city_line <- function(city) {
  
  fuel.data %>%
    filter(date >= max(date) %m-% years(2)) %>%
    mutate(
      year      = as.factor(year(date)),
      day_of_yr = yday(date)
    ) %>%
    ggplot(aes(x = day_of_yr, y = .data[[city]], colour = year)) +
    geom_line(aes(linetype = Actual_estimate)) +
    scale_linetype_manual(values = c("A" = "solid", "E" = "dashed"), guide = "none") +
    scale_x_continuous(
      breaks = yday(as.Date(paste0("2024-", 1:12, "-01"))),
      labels = month.abb
    ) +
    labs(
      title  = paste(city, "— Daily Fuel Prices"),
      x      = NULL,
      y      = "Price (cents per litre)",
      caption = "B* Fixed Income, ACCC",
      colour = "Year"
    ) +
    scale_colour_manual(values = bstar_colors) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor  = element_blank(),
      axis.text.x       = element_text(size = 7, colour = "gray20"),
      axis.text.y       = element_text(size = 7, colour = "gray20"),
      legend.text       = element_text(colour = "gray20"),
      legend.title      = element_text(colour = "gray20")
    )
}

plot_city_line("Sydney")
######### Single line chart - clearer without YoY data ###############################

#Reduce the fuel data to only include the next 1m of E's
fuel.data.reduced <- fuel.data %>%
  filter(
    Actual_estimate == "A" |
      (Actual_estimate == "E" & date <= ceiling_date(max(date[Actual_estimate == "A"]), "month") %m+% months(1) - days(1))
  )
  

plot_city_line_2 <- function(city) {
  
  cutoff <- as.Date("2026-02-01") #max(fuel.data.reduced$date) %m-% months(6)
  max_date <- max(fuel.data.reduced$date)
  
  monthly_dots <- monthly_avg %>%
    filter(year_month >= ceiling_date(cutoff, "month"),
           year_month <= floor_date(max_date, "month")) %>%
    mutate(year_month = year_month + days(14))
  
  fuel.data.reduced %>%
    filter(date >= cutoff) %>%
    ggplot(aes(x = date, y = .data[[city]])) +
    geom_line(aes(linetype = Actual_estimate), colour = bstar_colors[1]) +
    geom_point(
      data = monthly_dots,
      aes(x = year_month, y = .data[[city]]),
      colour = bstar_colors[2], size = 3, shape = 16
    ) +
    scale_linetype_manual(values = c("A" = "solid", "E" = "dashed"), guide = "none") +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%d-%b-%y"
    ) +
    labs(
      title  = paste(city, "— Daily Fuel Prices"),
      x      = NULL,
      y      = "Price (cents per litre)",
      caption = "B* Fixed Income, ACCC"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor  = element_blank(),
      axis.text.x       = element_text(size = 7, colour = "gray20"),
      axis.text.y       = element_text(size = 7, colour = "gray20"),
      legend.text       = element_text(colour = "gray20"),
      legend.title      = element_text(colour = "gray20")
    )
}

#plot_city_line_2("Sydney")

# ─────────────────────────────────────────────
# Comparing to the fuel index
# ─────────────────────────────────────────────
# Read ABS 6401.0 - table 3 - A130398320J automotive fuel index Australia
#--- Read raw data ---
raw.3 <- read_abs("6401.0", tables = 3) # National General Government
auto.aus.data <-raw.3 %>%
  filter(series_id=="A130398320J") %>%
  select(date,value) %>%
  rename(year_month=date)

reg.data <- left_join(monthly_avg, auto.aus.data, by = "year_month") %>%
  na.omit() %>%
  arrange(year_month) %>%
  mutate( Australia_mom = Australia / lag(Australia) - 1,
          value_mom     = value / lag(value) - 1
  ) %>%
  na.omit()  # drop first row where lag is NA

cor(reg.data$Australia_mom, reg.data$value_mom)
mod1 <- lm(value_mom ~ Australia_mom, data = reg.data)
summary(mod1)
#plot(mod1)

reg_chart <- reg.data %>%
  select(year_month, Australia_mom, value_mom) %>%
  pivot_longer(cols = c(Australia_mom, value_mom),
               names_to  = "series",
               values_to = "mom") %>%
  mutate(series = recode(series, "Australia_mom" = "Predicted", "value_mom" = "Actual")) %>%
  ggplot(aes(x = year_month, y = mom, colour = series)) +
  geom_line() +
  scale_colour_manual(values = bstar_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%b-%y") +
  labs(
    title   = "Predicted vs Actual — Australian Automotive Fuel Index",
    x       = NULL,
    y       = "MoM Growth",
    colour  = NULL,
    caption = "B* Fixed Income, ACCC, ABS"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 7, colour = "gray20"),
    axis.text.y      = element_text(size = 7, colour = "gray20"),
    legend.text       = element_text(colour = "gray20")
  )

# ─────────────────────────────────────────────
# Output
# ─────────────────────────────────────────────

cat("\n=== Step 1: Deduplicated data (rows after keeping latest forecast_date) ===\n")
print(fuel.data %>% select(forecast_date, date, Actual_estimate, all_of(cities)))

cat("\n=== Step 3: Monthly averages ===\n")
print(monthly_avg)

# ── Helper: render a ggplot to an inline base64 HTML img tag ──
gg_to_img <- function(p, width = 10, height = 5) {
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, plot = p, width = width, height = height, dpi = 150, bg = "white")
  b64 <- base64enc::base64encode(tmp)
  sprintf('<img src="data:image/png;base64,%s" style="width:100%%;margin:16px 0;">', b64)
}

# ── Format mom_growth: all city/Australia columns as % ──
mom_growth_fmt <- mom_growth %>%
  mutate(
    year_month = format(year_month, "%b-%y"),
    across(c(all_of(cities), "Australia"), ~ scales::percent(.x, accuracy = 0.01))
  )

# ── Format monthly_avg: date as mmm-yy ──
monthly_avg_fmt <- monthly_avg %>%
  mutate(year_month = format(year_month, "%b-%y"))

# ── Find which rows to bold ──
ae_index <- monthly_avg$Actual_estimate
last_A_row  <- max(which(ae_index == "A"))
first_E_row <- min(which(ae_index == "E"))
bold_rows   <- c(last_A_row, first_E_row)

# ── 1) Daily detail: last 20 actuals + next 7 forecasts ──
last_A_date <- max(fuel.data$date[fuel.data$Actual_estimate == "A"])

daily_detail_new <- fuel.data %>%
  filter(
    (Actual_estimate == "A" & date > last_A_date - days(20)) |
      (Actual_estimate == "E" & date <= last_A_date + days(7))
  ) %>%
  arrange(date) %>%
  select(date, Actual_estimate, all_of(cities), Australia) %>%
  mutate(
    date = format(date, "%d-%b-%y"),
    across(c(all_of(cities), "Australia"), ~ round(.x, 1))
  )

# Find the row indices of the last 3 actuals to bold
n_actuals <- sum(fuel.data$Actual_estimate[fuel.data$date > last_A_date - days(20) & fuel.data$date <= last_A_date] == "A")
bold_daily_rows <- (n_actuals - 2):n_actuals

# ── Build HTML ──
html_out <- paste0(
  '<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>TayTay tour Report</title>
  <style>
    body { font-family: Arial, sans-serif; max-width: 1100px; margin: 40px auto; color: #333F48; }
    h2   { color: #337179; border-bottom: 2px solid #DE8C59; padding-bottom: 6px; }
    h3   { color: #622128; }
    table.kable { font-size: 12px; }
  </style>
</head>
<body>
<h2>Fuel Price Report</h2>',
  
  # ── 1) Daily detail table FIRST ──
  "<h3>Daily Prices — Last 20 Actuals & Next 7 Forecast Days</h3>",
  daily_detail_new %>%
    kable("html", digits = 1) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = TRUE, font_size = 14) %>%
    row_spec(bold_daily_rows, bold = TRUE) %>%
    as.character(),
  
  # ── 2) Monthly averages table ──
  "<h3>Monthly Averages</h3>",
  monthly_avg_fmt %>%
    kable("html", digits = 1, format.args = list(big.mark = ",")) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = TRUE, font_size = 14) %>%
    row_spec(bold_rows, bold = TRUE) %>%
    as.character(),
  
  # ── 2) Australia line chart (plot_city_line_2) between monthly avg and MoM ──
  "<h3>Australia — Daily Fuel Prices</h3>",
  gg_to_img(plot_city_line_2("Australia")),
  
  # ── 2) MoM growth table ──
  "<h3>Month-on-Month Growth</h3>",
  mom_growth_fmt %>%
    kable("html") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = TRUE, font_size = 14) %>%
    row_spec(bold_rows[bold_rows <= nrow(mom_growth_fmt)] - 1, bold = TRUE) %>%
    as.character(),
  
  # ── Regression chart ──
  "<h3>Australia MoM vs ABS Automotive Fuel Index</h3>",
  gg_to_img(
    reg_chart +
      scale_colour_manual(values = bstar_colors) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
      scale_x_date(date_breaks = "1 year", date_labels = "%b-%y") +
      labs(x = NULL, y = "MoM Growth", colour = "Series",
           caption = "B* Fixed Income, ACCC, ABS") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(size = 7, colour = "gray20"),
            axis.text.y = element_text(size = 7, colour = "gray20"))
  ),
  
  # ── 3) City-specific charts: plot_city_line_2 (short-range) ──
  "<h3>Capital City Daily Fuel Prices — Recent</h3>",
  paste(sapply(cities, function(city) {
    paste0("<h4>", city, "</h4>", gg_to_img(plot_city_line_2(city)))
  }), collapse = "\n"),
  
  # ── 3) All cities + Australia: plot_city_line (YoY overlay) ──
  "<h3>Capital City Daily Fuel Prices — Year-on-Year</h3>",
  "<h4>Australia</h4>",
  gg_to_img(plot_city_line("Australia")),
  paste(sapply(cities, function(city) {
    paste0("<h4>", city, "</h4>", gg_to_img(plot_city_line(city)))
  }), collapse = "\n"),
  
  "</body></html>"
)

# ── Write to docs folder (served by GitHub Pages) ──
if (!dir.exists("docs")) dir.create("docs", recursive = TRUE)
out_path <- "docs/index.html"
writeLines(html_out, out_path)
cat("Report written to:", out_path, "\n")
