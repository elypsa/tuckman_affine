# ==============================================================================
# 01_data.R - Data Preparation for Gauss+ Model
# ==============================================================================

library(dplyr)
library(readr)

#' Load Gurkaynak-Sack-Wright (GSW) yield curve data
#'
#' Reads the Federal Reserve Board's GSW yield curve dataset and extracts
#' zero-coupon yields (SVENY series). Data are converted from percent to
#' decimal form.
#'
#' @param path Path to feds200628.csv
#' @return data.frame with columns: date, SVENY01, SVENY02, ..., SVENY30
#'
load_gsw <- function(path = "data-raw/feds200628.csv") {
  # Read CSV, skipping the 9 header rows to get to the data
  dat <- read_csv(path, skip = 9, show_col_types = FALSE)

  # Select date column and all SVENY (zero-coupon yield) columns
  dat <- dat %>%
    select(date = Date, starts_with("SVENY")) %>%
    # Convert date to Date type
    mutate(date = as.Date(date)) %>%
    # CRITICAL: GSW data is in percent, convert to decimal
    mutate(across(starts_with("SVENY"), ~ .x / 100)) %>%
    # Remove any rows with missing values
    na.omit()

  return(dat)
}


#' Load Federal Funds Target Rate (short rate proxy)
#'
#' Reads the Fed funds target range (upper and lower bounds) from FRED data
#' and computes the midpoint. This serves as the short rate r_t in the model.
#'
#' @param path Path to fed_funds_target.csv
#' @return data.frame with columns: date, r_t
#'
load_target <- function(path = "data-raw/fed_funds_target.csv") {
  dat <- read_csv(path, show_col_types = FALSE)

  dat <- dat %>%
    # Convert date to Date type
    mutate(date = as.Date(observation_date)) %>%
    # Compute midpoint of target range, converting from percent to decimal
    # Formula: (upper + lower) / 2, then divide by 100
    # Combined: (upper + lower) / 200
    mutate(r_t = (DFEDTARU + DFEDTARL) / 200) %>%
    # Keep only date and r_t
    select(date, r_t) %>%
    # Remove any missing values
    na.omit()

  return(dat)
}


#' Build the estimation panel
#'
#' Combines GSW yield curve data with fed funds target, computes one-year
#' forward rates, and creates exponential time-decay weights for estimation.
#'
#' @param start_date Start of sample (default: "2014-01-05")
#' @param end_date End of sample (default: "2022-01-21")
#' @return list with elements:
#'   - panel: data.frame with yields, forwards, short rate
#'   - weights: numeric vector of exponential weights (last obs = 1)
#'   - maturities: numeric vector of maturities used in estimation
#'
build_panel <- function(start_date = "2014-01-05", end_date = "2022-01-21") {
  # Load both datasets
  gsw <- load_gsw()
  target <- load_target()

  # Inner join: keep only dates present in both datasets
  panel <- inner_join(gsw, target, by = "date")

  # Filter to estimation sample period
  panel <- panel %>%
    filter(date >= as.Date(start_date), date <= as.Date(end_date))

  # Compute one-year forward rates from continuous zero yields
  # Formula from spec §5.1: f(tau, tau'=1) = (tau+1)*y(tau+1) - tau*y(tau)
  # These are continuously compounded, matching the model formulation
  panel <- panel %>%
    mutate(
      fwd_02 = 3 * SVENY03 - 2 * SVENY02,
      fwd_03 = 4 * SVENY04 - 3 * SVENY03,
      fwd_10 = 11 * SVENY11 - 10 * SVENY10,
      fwd_14 = 15 * SVENY15 - 14 * SVENY14,
      fwd_15 = 16 * SVENY16 - 15 * SVENY15
    )

  # Compute exponential weights: 0.8^(age_years)
  # Age is measured backward from the last observation
  # This gives weight = 1 for the most recent observation
  age_years <- as.numeric(max(panel$date) - panel$date) / 365.25
  weights <- 0.8^age_years

  # Define maturities used in estimation
  # Need yields at these points for fitting
  maturities <- c(1:10, 12, 15, 17, 20, 25, 27, 30)

  return(list(
    panel = panel,
    weights = weights,
    maturities = maturities
  ))
}


#' Quick diagnostic check for panel data
#'
#' Prints summary statistics to verify data loading was successful.
#' Use this to catch common errors like forgetting to divide by 100.
#'
#' @param panel_result Output from build_panel()
#'
check_panel <- function(panel_result) {
  cat("=== PANEL DIAGNOSTICS ===\n\n")

  cat("1. Dimensions:\n")
  cat("   Rows:", nrow(panel_result$panel), "\n")
  cat("   Columns:", ncol(panel_result$panel), "\n\n")

  cat("2. Date range:\n")
  cat("   Start:", as.character(min(panel_result$panel$date)), "\n")
  cat("   End:", as.character(max(panel_result$panel$date)), "\n\n")

  cat("3. Sample yields (should be in decimals, e.g., 0.02 for 2%):\n")
  print(summary(panel_result$panel[, c("SVENY02", "SVENY10", "r_t")]))
  cat("\n")

  cat("4. Sample forwards:\n")
  fwd_cols <- grep("^fwd_", names(panel_result$panel), value = TRUE)
  if (length(fwd_cols) > 0) {
    print(summary(panel_result$panel[, fwd_cols]))
  } else {
    cat("   WARNING: No forward columns found!\n")
  }
  cat("\n")

  cat("5. Weights:\n")
  cat("   First (oldest):", round(panel_result$weights[1], 4), "\n")
  cat(
    "   Last (newest):",
    round(panel_result$weights[length(panel_result$weights)], 4),
    "\n"
  )
  # cat("   Should be: ~0.05-0.10 and 1.0000\n\n")

  cat("6. Maturities:\n")
  cat("   ", panel_result$maturities, "\n\n")

  cat("=== END DIAGNOSTICS ===\n")
}
