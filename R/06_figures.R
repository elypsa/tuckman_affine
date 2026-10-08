# ==============================================================================
# 06_figures.R - Figure Generation for Gauss+ Model
# ==============================================================================

library(ggplot2)

#' Figure 9.5: Regression coefficients on 2-year and 10-year benchmarks
#'
#' Compare empirical regression slopes (from weighted OLS) to model-implied
#' slopes from the estimated alpha parameters.
#'
#' @param stage1_result Output from stage1_alpha(), containing slopes_emp and slopes_model
#' @param output_path Optional path to save the figure
#' @return ggplot object
#'
figure_9_5 <- function(stage1_result, output_path = NULL) {
  slopes_emp <- stage1_result$slopes_emp
  slopes_model <- stage1_result$slopes_model
  maturities <- stage1_result$maturities

  # Reshape data for plotting
  df <- data.frame(
    maturity = rep(maturities, 2),
    coef_on_2y = c(slopes_emp[1, ], slopes_model[1, ]),
    coef_on_10y = c(slopes_emp[2, ], slopes_model[2, ]),
    source = rep(c("Empirical", "Model"), each = length(maturities))
  )

  # Reshape to long format for ggplot
  df_long <- data.frame(
    maturity = rep(df$maturity, 2),
    coefficient = c(df$coef_on_2y, df$coef_on_10y),
    benchmark = rep(c("on 2-year", "on 10-year"), each = nrow(df)),
    source = rep(df$source, 2)
  )

  p <- ggplot(df_long, aes(x = factor(maturity), y = coefficient,
                           fill = source)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7) +
    facet_wrap(~ benchmark, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("Empirical" = "#2E86AB", "Model" = "#A23B72")) +
    labs(
      title = "Figure 9.5: Comovement of Yields with Benchmark Yields",
      subtitle = "Regression coefficients: empirical vs. model-implied",
      x = "Maturity (years)",
      y = "Coefficient",
      fill = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      strip.text = element_text(face = "bold")
    )

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 10, height = 7, dpi = 300)
    cat("Figure 9.5 saved to:", output_path, "\n")
  }

  return(p)
}


#' Figure 9.6: Yield volatility term structure
#'
#' Compare empirical yield volatilities (from daily changes) to model-implied
#' volatilities from the estimated sigma parameters.
#'
#' @param stage2_result Output from stage2_sigma(), containing vol_emp_bp and vol_model_bp
#' @param output_path Optional path to save the figure
#' @return ggplot object
#'
figure_9_6 <- function(stage2_result, output_path = NULL) {
  vol_emp <- stage2_result$vol_emp_bp
  vol_model <- stage2_result$vol_model_bp
  maturities <- stage2_result$maturities

  df <- data.frame(
    maturity = maturities,
    empirical = vol_emp,
    model = vol_model
  )

  # Reshape to long format
  df_long <- data.frame(
    maturity = rep(df$maturity, 2),
    volatility = c(df$empirical, df$model),
    source = rep(c("Empirical", "Model"), each = nrow(df))
  )

  p <- ggplot(df_long, aes(x = maturity, y = volatility,
                           color = source, shape = source)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Empirical" = "#2E86AB", "Model" = "#A23B72")) +
    scale_shape_manual(values = c("Empirical" = 16, "Model" = 17)) +
    labs(
      title = "Figure 9.6: Term Structure of Yield Volatility",
      subtitle = "Annualized volatilities in basis points: empirical vs. model-implied",
      x = "Maturity (years)",
      y = "Volatility (bp)",
      color = "",
      shape = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10)
    )

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 10, height = 6, dpi = 300)
    cat("Figure 9.6 saved to:", output_path, "\n")
  }

  return(p)
}


#' Figure 9.7: Instantaneous forward loadings by factor
#'
#' Plot how instantaneous forwards at different maturities load on each
#' of the three factors (r, m, l).
#'
#' @param alpha Vector c(a_r, a_m, a_l)
#' @param tau_grid Vector of maturities to plot (defaults to 0.25 to 30 years)
#' @param output_path Optional path to save the figure
#' @return ggplot object
#'
figure_9_7 <- function(alpha, tau_grid = seq(0.25, 30, by = 0.25),
                       output_path = NULL) {
  source("R/02_model.R")

  # Compute instantaneous forward loadings
  Ups_inst <- Ups_fwd_inst(tau_grid, alpha)

  df <- data.frame(
    maturity = rep(tau_grid, 3),
    loading = c(Ups_inst[, 1], Ups_inst[, 2], Ups_inst[, 3]),
    factor = rep(c("Short (r)", "Medium (m)", "Long (l)"), each = length(tau_grid))
  )

  p <- ggplot(df, aes(x = maturity, y = loading, color = factor)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = c(
      "Short (r)" = "#E63946",
      "Medium (m)" = "#2E86AB",
      "Long (l)" = "#06A77D"
    )) +
    labs(
      title = "Figure 9.7: Instantaneous Forward Rate Loadings",
      subtitle = "How forwards at different maturities load on each factor",
      x = "Maturity (years)",
      y = "Loading",
      color = "Factor"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10)
    )

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 10, height = 6, dpi = 300)
    cat("Figure 9.7 saved to:", output_path, "\n")
  }

  return(p)
}


#' Figure 9.8: Extracted factors and 2-year forward rate
#'
#' Time series plot of the three extracted factors (r, m, l) plus the observed
#' 2-year one-year forward rate. The long factor l is transformed to its
#' 10-year conditional mean L(l) for interpretability.
#'
#' @param panel Data frame with date column and fwd_02
#' @param factors Matrix T x 3 with columns (r, m, l)
#' @param mu Scalar long-run mean
#' @param alpha_l Long factor mean reversion speed
#' @param output_path Optional path to save the figure
#' @return ggplot object
#'
figure_9_8 <- function(panel, factors, mu, alpha_l, output_path = NULL) {
  source("R/04_factors.R")

  # Transform long factor to 10-year conditional mean
  L_l <- L_transform(factors[, "l"], mu, alpha_l)

  # Build data frame for plotting
  df <- data.frame(
    date = panel$date,
    r = factors[, "r"],
    m = factors[, "m"],
    L_l = L_l,
    fwd_2y = panel$fwd_02
  )

  # Reshape to long format
  df_long <- data.frame(
    date = rep(df$date, 4),
    value = c(df$r, df$m, df$L_l, df$fwd_2y),
    series = rep(
      c("Short rate (r)", "Medium factor (m)",
        "Long factor L(l)", "2-year forward"),
      each = nrow(df)
    )
  )

  # Set factor levels for legend order
  df_long$series <- factor(df_long$series, levels = c(
    "Short rate (r)", "Medium factor (m)",
    "Long factor L(l)", "2-year forward"
  ))

  p <- ggplot(df_long, aes(x = date, y = value, color = series)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = c(
      "Short rate (r)" = "#E63946",
      "Medium factor (m)" = "#2E86AB",
      "Long factor L(l)" = "#06A77D",
      "2-year forward" = "#F77F00"
    )) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(
      title = "Figure 9.8: Extracted Factors and 2-Year Forward Rate",
      subtitle = "Long factor shown as 10-year conditional mean L(l)",
      x = NULL,
      y = "Rate",
      color = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      panel.grid.minor = element_blank()
    )

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 12, height = 6, dpi = 300)
    cat("Figure 9.8 saved to:", output_path, "\n")
  }

  return(p)
}


#' Figure 9.10: Risk premium on the 10-year forward rate over time
#'
#' Time series plot showing the 10-year forward rate and its associated
#' risk premium over time.
#'
#' @param panel Data frame with date column and fwd_10
#' @param lambda_t Vector of price-of-risk values over time
#' @param alpha Parameter vector c(a_r, a_m, a_l)
#' @param sigma_l Long factor volatility
#' @param output_path Optional path to save the figure
#' @return ggplot object
#'
figure_9_10 <- function(panel, lambda_t, alpha, sigma_l, output_path = NULL) {
  source("R/05_premia.R")

  # Compute risk premium on 10-year forward
  rp_10y <- lambda_t * rp_unit(10, alpha, c(0, sigma_l, 0))

  # Build data frame
  df <- data.frame(
    date = panel$date,
    fwd_10y = panel$fwd_10,
    risk_premium = rp_10y
  )

  # Remove any NAs
  df <- df[complete.cases(df), ]

  # Reshape to long format for plotting
  df_long <- data.frame(
    date = rep(df$date, 2),
    value = c(df$fwd_10y, df$risk_premium),
    series = rep(c("10-Year Forward", "Risk Premium"), each = nrow(df))
  )

  # Set factor levels for legend order
  df_long$series <- factor(df_long$series, levels = c("10-Year Forward", "Risk Premium"))

  p <- ggplot(df_long, aes(x = date, y = value, color = series)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c(
      "10-Year Forward" = "#2E86AB",
      "Risk Premium" = "#E63946"
    )) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(
      title = "Figure 9.10: Risk Premium on 10-Year Forward Rate",
      subtitle = "10-year forward rate and its risk premium over time",
      x = NULL,
      y = "Rate",
      color = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      panel.grid.minor = element_blank()
    )

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 12, height = 6, dpi = 300)
    cat("Figure 9.10 saved to:", output_path, "\n")
  }

  return(p)
}
