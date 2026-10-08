# ==============================================================================
# run_estimation.R - Full Estimation Pipeline
# ==============================================================================
#
# Run the complete three-stage Gauss+ estimation.
# This script assumes you have downloaded the data files to data-raw/
# ==============================================================================

source("R/01_data.R")
source("R/02_model.R")
source("R/03_estimate.R")
source("R/04_factors.R")
source("R/05_premia.R")
source("R/06_figures.R")

library(dplyr)

cat("=== GAUSS+ MODEL ESTIMATION ===\n\n")

# Step 1: Load and prepare data
cat("Step 1: Loading data...\n")
result <- build_panel()
panel <- result$panel
weights <- result$weights
maturities <- result$maturities

cat(
  "  Sample:",
  nrow(panel),
  "business days from",
  as.character(min(panel$date)),
  "to",
  as.character(max(panel$date)),
  "\n"
)
cat("  Maturities:", maturities, "\n\n")

# Quick diagnostic
check_panel(result)
cat("\n")

# Step 2: Three-stage estimation
cat("Step 2: Running three-stage estimation...\n")

est <- estimate_all(
  panel = panel,
  weights = weights,
  maturities = maturities,
  start_alpha = c(1.0, 0.6, 0.02),
  start_sigma = c(0.01, 0.01, 0.2)
)

# Step 3: Extract factors (already computed in stage 3)
cat("Step 3: Extracting latent factors...\n")
factors <- est$params$mu$factors

cat("  Factors extracted for all", nrow(factors), "dates\n\n")

# Step 4: Extract price of risk lambda
cat("Step 4: Extracting time-varying price of risk...\n")
alpha_vec <- c(est$params$a_r, est$params$a_m, est$params$a_l)
sigma_vec <- c(est$params$s_m, est$params$s_l, est$params$rho)

lambda_t <- solve_lambda(
  fwd_tau = panel$fwd_14,
  fwd_tau_prime = panel$fwd_15,
  tau = 14,
  tau_prime = 15,
  alpha = alpha_vec,
  sigma = sigma_vec
)

cat("  Lambda extracted for all", length(lambda_t), "dates\n")
cat("  Mean lambda:", round(mean(lambda_t, na.rm = TRUE), 3), "\n")
cat("  Median lambda:", round(median(lambda_t, na.rm = TRUE), 3), "\n\n")

# Step 5: Summary
cat("=== FINAL PARAMETER ESTIMATES ===\n\n")

cat("Mean reversion speeds (alpha):\n")
cat(
  "  a_r =",
  round(est$params$a_r, 4),
  "  (half-life:",
  round(log(2) / est$params$a_r, 2),
  "years)\n"
)
cat(
  "  a_m =",
  round(est$params$a_m, 4),
  "  (half-life:",
  round(log(2) / est$params$a_m, 2),
  "years)\n"
)
cat(
  "  a_l =",
  round(est$params$a_l, 4),
  "  (half-life:",
  round(log(2) / est$params$a_l, 1),
  "years)\n\n"
)

cat("Volatilities (sigma):\n")
cat(
  "  s_m =",
  round(est$params$s_m, 5),
  " (",
  round(est$params$s_m * 10000, 1),
  "bp )\n"
)
cat(
  "  s_l =",
  round(est$params$s_l, 5),
  " (",
  round(est$params$s_l * 10000, 1),
  "bp )\n"
)
cat("  rho =", round(est$params$rho, 3), "\n\n")

cat("Long-run mean:\n")
cat(
  "  mu =",
  round(est$params$mu$mu, 5),
  " (",
  round(est$params$mu$mu * 100, 2),
  "% )\n\n"
)

# Compare to reference values from Table 9.1
par_ref <- list(
  a_r = 1.0547,
  a_m = 0.6358,
  a_l = 0.0165,
  s_m = 0.01092,
  s_l = 0.00964,
  rho = 0.212,
  mu = 0.10555
)

cat("=== COMPARISON TO TABLE 9.1 ===\n\n")
cat("Parameter  | Estimated | Reference | Difference\n")
cat("-----------|-----------|-----------|------------\n")
cat(sprintf(
  "a_r        | %8.4f  | %8.4f  | %8.4f\n",
  est$params$a_r,
  par_ref$a_r,
  est$params$a_r - par_ref$a_r
))
cat(sprintf(
  "a_m        | %8.4f  | %8.4f  | %8.4f\n",
  est$params$a_m,
  par_ref$a_m,
  est$params$a_m - par_ref$a_m
))
cat(sprintf(
  "a_l        | %8.4f  | %8.4f  | %8.4f\n",
  est$params$a_l,
  par_ref$a_l,
  est$params$a_l - par_ref$a_l
))
cat(sprintf(
  "s_m        | %8.5f  | %8.5f  | %8.5f\n",
  est$params$s_m,
  par_ref$s_m,
  est$params$s_m - par_ref$s_m
))
cat(sprintf(
  "s_l        | %8.5f  | %8.5f  | %8.5f\n",
  est$params$s_l,
  par_ref$s_l,
  est$params$s_l - par_ref$s_l
))
cat(sprintf(
  "rho        | %8.3f  | %8.3f  | %8.3f\n",
  est$params$rho,
  par_ref$rho,
  est$params$rho - par_ref$rho
))
cat(sprintf(
  "mu         | %8.5f  | %8.5f  | %8.5f\n\n",
  est$params$mu$mu,
  par_ref$mu,
  est$params$mu$mu - par_ref$mu
))

# Generate figures
cat("Generating figures...\n")

cat("  Figure 9.5 (regression coefficients)...\n")
fig_9_5 <- figure_9_5(est$stage1)

cat("  Figure 9.6 (yield volatilities)...\n")
fig_9_6 <- figure_9_6(est$stage2)

cat("  Figure 9.7 (forward loadings)...\n")
alpha_vec <- c(est$params$a_r, est$params$a_m, est$params$a_l)
fig_9_7 <- figure_9_7(alpha_vec)

cat("  Figure 9.8 (factors time series)...\n")
fig_9_8 <- figure_9_8(panel, factors, est$params$mu$mu, est$params$a_l)

cat("  Figure 9.10 (risk premium)...\n")
fig_9_10 <- figure_9_10(panel, lambda_t, alpha_vec, sigma_vec[2])

cat("\n")

# Save results
cat("Saving results to data/estimation_results.rds...\n")
saveRDS(
  list(
    params = est$params,
    factors = factors,
    lambda_t = lambda_t,
    panel = panel,
    weights = weights,
    stage1 = est$stage1,
    stage2 = est$stage2,
    figures = list(
      fig_9_5 = fig_9_5,
      fig_9_6 = fig_9_6,
      fig_9_7 = fig_9_7,
      fig_9_8 = fig_9_8,
      fig_9_10 = fig_9_10
    )
  ),
  "data/estimation_results.rds"
)

cat("\nDone!\n")
cat("\nTo view figures, load the results and print them:\n")
cat("  results <- readRDS('data/estimation_results.rds')\n")
cat("  print(results$figures$fig_9_5)\n")
cat("  print(results$figures$fig_9_6)\n")
cat("  print(results$figures$fig_9_7)\n")
cat("  print(results$figures$fig_9_8)\n")
cat("  print(results$figures$fig_9_10)\n")
