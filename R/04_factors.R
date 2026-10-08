# ==============================================================================
# 04_factors.R - Factor Extraction
# ==============================================================================

source("R/02_model.R")

#' Extract latent factors (m_t, l_t) from observed data
#'
#' With parameters fixed, invert the model to recover the unobserved factors.
#' Use 2-year and 10-year one-year forwards to pin down the two factors exactly.
#'
#' Short rate r_t is observed, so only (m, l) need extraction.
#'
#' This function uses a vectorized implementation that processes all time periods
#' at once, avoiding loops for better performance.
#'
#' @param panel Data frame with columns: fwd_02, fwd_10, r_t
#' @param alpha Parameter vector c(a_r, a_m, a_l)
#' @param sigma Parameter vector c(s_m, s_l, rho)
#' @param mu Scalar long-run mean
#' @return Matrix T x 3 with columns (r, m, l)
#'
extract_factors <- function(panel, alpha, sigma, mu) {
  TT <- nrow(panel)

  # Observed short rate
  r_t <- panel$r_t

  # Benchmark maturities for factor extraction
  bench_tau <- c(2, 10)

  # Forward loadings at tau = 2 and tau = 10 (one-year tenor)
  Ups_f <- Ups_fwd(bench_tau, alpha, tenor = 1) # 2 x 3

  # Convexity terms (constant across time)
  C_f <- C_fwd(bench_tau, alpha, sigma, tenor = 1) # length 2

  # Precompute inverse of the 2 x 2 loading matrix on (m, l)
  Lsolve <- solve(Ups_f[, 2:3, drop = FALSE])

  # Observed forwards: 2 x T matrix
  f_obs <- rbind(panel$fwd_02, panel$fwd_10)

  # Short rate contribution: 2 x T matrix
  rt_term <- outer(Ups_f[, 1], r_t)

  # Constants in forward equation (independent of t)
  const <- mu * (1 - rowSums(Ups_f)) - C_f # length 2

  # Extract (m, l) for all dates at once
  # Solve: Ups_f[, 2:3] %*% (m, l)' = f_obs - const - rt_term
  ml <- t(Lsolve %*% (f_obs - const - rt_term)) # T x 2

  # Combine with observed short rate
  factors <- cbind(r = r_t, ml)
  colnames(factors)[2:3] <- c("m", "l")

  return(factors)
}


#' Transform long factor to 10-year conditional mean
#'
#' For plotting purposes, transform l_t to its 10-year-ahead expected value
#' under the risk-neutral measure. This is just the OU conditional mean.
#'
#' L(l_t) = mu * (1 - exp(-10 * alpha_l)) + l_t * exp(-10 * alpha_l)
#'
#' @param l_t Vector of long factor values
#' @param mu Long-run mean
#' @param alpha_l Long factor mean reversion speed
#' @return Transformed long factor
#'
L_transform <- function(l_t, mu, alpha_l) {
  w <- exp(-10 * alpha_l)
  L_t <- mu * (1 - w) + l_t * w
  return(L_t)
}
