# ==============================================================================
# 03_estimate.R - Three-Stage Estimation for Gauss+ Model
# ==============================================================================

source("R/02_model.R")
library(nloptr)

#' Net out the short rate component from yields
#'
#' The short rate loading is Υ_s(τ, α_r) = B_1(τ) = (1-exp(-α_r τ))/(α_r τ)
#' It depends only on α_r, not α_m or α_l.
#'
#' After netting, yields load only on (m, l), reducing the problem from 3D to 2D.
#'
#' @param yields Matrix T x N of yields at N maturities
#' @param r_t Vector of length T, short rate
#' @param tau Vector of length N, maturities (in years)
#' @param alpha_r Scalar, short rate mean reversion speed
#' @return Matrix T x N of netted yields
#'
net_short_rate <- function(yields, r_t, tau, alpha_r) {
  # Short rate loading: first column of B(tau, alpha) when alpha = (alpha_r, *, *)
  alpha_temp <- c(alpha_r, 1, 1) # Other values don't matter for first column
  B_s <- B(tau, alpha_temp)[, 1] # Length N

  # Subtract r_t * B_s(tau) from each yield
  # r_t is T x 1, B_s is 1 x N, outer product is T x N
  yhat <- yields - outer(r_t, B_s)

  return(yhat)
}

#' Short-rate loading Υ_s(τ, α_r) = B₁(τ)
#'
#' The first column of A(α)⁻¹ is (1,0,0)', so the short-rate loading on a
#' yield is just B₁(τ) and depends on α_r alone.
short_rate_loading <- function(tau, alpha_r) {
  x <- alpha_r * tau
  -expm1(-x) / x # (1 - exp(-x)) / x, stable for small x
}

#' Net the observed short rate out of yields (levels)
net_short_rate <- function(yields, r_t, tau, alpha_r) {
  yields - outer(r_t, short_rate_loading(tau, alpha_r))
}

#' Stage 1: Estimate alpha from regression slopes
#'
#' Work in changes so mu and C drop out. For a candidate alpha_r, net the
#' short rate out of the yield changes, regress all netted changes on the
#' two benchmark netted changes, and compare those empirical slopes to the
#' model-implied slopes Ups %*% solve(Ups_benchmark).
#'
#' Both sides of the comparison depend on alpha, so the regression has to
#' live inside the objective.
#'
#' @param panel Data frame with yield columns SVENY02, SVENY03, ...
#' @param r_t Vector of short rates, same length as nrow(panel)
#' @param weights Vector of time-decay weights, same length as nrow(panel)
#' @param maturities Integer vector of maturities in years; must contain 2 and 10
#' @param start_alpha Initial guess c(alpha_r, alpha_m, alpha_l)
#' @return List with alpha, the fitted and empirical slopes, and diagnostics
stage1_alpha <- function(
  panel,
  r_t,
  weights,
  maturities,
  start_alpha = c(1.0, 0.6, 0.02)
) {
  stopifnot(all(c(2, 10) %in% maturities))
  stopifnot(
    start_alpha[1] > start_alpha[2],
    start_alpha[2] > start_alpha[3],
    start_alpha[3] > 0
  )

  Y <- as.matrix(panel[, paste0("SVENY", sprintf("%02d", maturities))])
  stopifnot(length(r_t) == nrow(Y), length(weights) == nrow(Y))

  # Daily changes; weights indexed by the later date of each change
  dY <- diff(Y)
  dr <- diff(r_t)
  w <- weights[-1]

  # match() rather than which(): guarantees column order is (2y, 10y)
  bench_idx <- match(c(2, 10), maturities)

  unpack <- function(p) {
    a_l <- exp(p[3])
    a_m <- a_l + exp(p[2])
    a_r <- a_m + exp(p[1])
    c(a_r, a_m, a_l)
  }

  # Empirical slopes for a given alpha_r: 2 x N
  emp_slopes <- function(a_r) {
    dYhat <- dY - outer(dr, short_rate_loading(maturities, a_r))
    dYbhat <- dYhat[, bench_idx, drop = FALSE]
    solve(crossprod(dYbhat, w * dYbhat), crossprod(dYbhat, w * dYhat))
  }

  # Model-implied slopes for a given alpha: 2 x N
  mod_slopes <- function(alpha) {
    Ups_full <- Ups_yield(maturities, alpha) # N x 3
    t(
      Ups_full[, 2:3, drop = FALSE] %*%
        solve(Ups_full[bench_idx, 2:3, drop = FALSE])
    ) # 2 x N
  }

  objective <- function(p) {
    alpha <- unpack(p)
    norm(mod_slopes(alpha) - emp_slopes(alpha[1]), "F")
  }

  p_start <- log(c(
    start_alpha[1] - start_alpha[2],
    start_alpha[2] - start_alpha[3],
    start_alpha[3]
  ))

  res <- nloptr(
    x0 = p_start,
    eval_f = objective,
    opts = list(
      algorithm = "NLOPT_LN_NELDERMEAD",
      maxeval = 5000,
      xtol_rel = 1e-6
    )
  )

  if (res$status < 0) {
    warning("stage1_alpha: nloptr status ", res$status, " - ", res$message)
  }

  alpha_hat <- unpack(res$solution)
  names(alpha_hat) <- c("a_r", "a_m", "a_l")

  beta_emp <- emp_slopes(alpha_hat["a_r"])
  beta_mod <- mod_slopes(alpha_hat)
  dimnames(beta_emp) <- list(c("on_2y", "on_10y"), paste0("m", maturities))
  dimnames(beta_mod) <- dimnames(beta_emp)

  list(
    alpha = alpha_hat,
    half_life = log(2) / alpha_hat,
    slopes_emp = beta_emp, # for the comovement chart
    slopes_model = beta_mod,
    maturities = maturities,
    objective = res$objective,
    status = res$status,
    message = res$message
  )
}


#' Stage 2: Estimate sigma from yield volatility
#'
#' With alpha fixed, match model-implied yield variance to empirical variance.
#' Model variance: diag(Ups %*% Omega %*% Omega' %*% Ups')
#' Empirical variance: annualized from daily squared changes
#'
#' @param panel Data frame with yields
#' @param r_t Vector of short rates
#' @param weights Vector of time-decay weights
#' @param maturities Vector of maturities
#' @param alpha Fixed alpha from stage 1
#' @param start_sigma Initial guess c(sigma_m, sigma_l, rho)
#' @param scale Indicates whether variance or volatility is used in the objective function
#' @return List with sigma and diagnostics
#'
stage2_sigma <- function(
  panel,
  r_t,
  weights,
  maturities,
  alpha,
  start_sigma = c(0.01, 0.01, 0.2),
  scale = c("vol", "var")
) {
  scale <- match.arg(scale)
  stopifnot(length(alpha) == 3, all(alpha > 0))
  a_r <- unname(alpha[1])

  Y <- as.matrix(panel[, paste0("SVENY", sprintf("%02d", maturities))])
  stopifnot(length(r_t) == nrow(Y), length(weights) == nrow(Y))

  dY <- diff(Y)
  dr <- diff(r_t)
  w <- weights[-1]

  # Net the observed short rate out of the yield changes. The model gives r
  # no diffusion term, so target moves are not variance that sigma should be
  # asked to explain. Matters at the front end on FOMC days (March 2020 most
  # of all); negligible past ~5y.
  dYhat <- dY - outer(dr, short_rate_loading(maturities, a_r))

  v_emp <- colSums(w * dYhat^2) / sum(w) * 252 # annualised variance

  # alpha is fixed at this stage, so the loadings are constant
  Ups <- Ups_yield(maturities, alpha) # N x 3

  # diag(Ups %*% OmOm' %*% t(Ups)) without building the N x N matrix.
  # Omega %*% t(Omega) in closed form: first row and column are zero.
  model_var <- function(sigma) {
    s_m <- sigma[1]
    s_l <- sigma[2]
    rho <- sigma[3]
    S <- matrix(0, 3, 3)
    S[2, 2] <- s_m^2
    S[2, 3] <- S[3, 2] <- rho * s_m * s_l
    S[3, 3] <- s_l^2
    rowSums((Ups %*% S) * Ups) # length N
  }

  target <- if (scale == "vol") sqrt(v_emp) else v_emp

  objective <- function(p) {
    sigma <- c(exp(p[1]), exp(p[2]), tanh(p[3]))
    vm <- model_var(sigma)
    if (any(!is.finite(vm)) || any(vm < 0)) {
      return(1e6)
    }
    fitted <- if (scale == "vol") sqrt(vm) else vm
    sqrt(sum((fitted - target)^2))
  }

  p_start <- c(log(start_sigma[1]), log(start_sigma[2]), atanh(start_sigma[3]))

  res <- nloptr(
    x0 = p_start,
    eval_f = objective,
    opts = list(
      algorithm = "NLOPT_LN_NELDERMEAD",
      maxeval = 10000,
      xtol_rel = 1e-6
    )
  )

  if (res$status < 0) {
    warning("stage2_sigma: nloptr status ", res$status, " - ", res$message)
  }

  p_opt <- res$solution
  sigma_hat <- c(exp(p_opt[1]), exp(p_opt[2]), tanh(p_opt[3]))
  names(sigma_hat) <- c("s_m", "s_l", "rho")

  v_model <- model_var(sigma_hat)

  nm <- paste0("m", maturities)
  names(v_emp) <- names(v_model) <- nm

  list(
    sigma = sigma_hat,
    vol_emp_bp = sqrt(v_emp) * 1e4, # what the chart plots
    vol_model_bp = sqrt(v_model) * 1e4,
    var_emp = v_emp,
    var_model = v_model,
    maturities = maturities,
    scale = scale,
    objective = res$objective,
    status = res$status,
    message = res$message
  )
}

#' Stage 3: Estimate mu from yield levels
#'
#' With alpha and sigma fixed, only mu remains. For each candidate mu,
#' extract factors and compute model yields, then minimize fit error.
#'
#' @param panel Data frame with yields and forwards
#' @param r_t Vector of short rates
#' @param weights Vector of time-decay weights
#' @param maturities Vector of maturities
#' @param alpha Fixed alpha from stage 1
#' @param sigma Fixed sigma from stage 2
#' @param start_mu Initial guess for mu
#' @param interval Interval over which mu is searched
#' @return Scalar mu
#'
stage3_mu <- function(
  panel,
  r_t,
  weights,
  maturities,
  alpha,
  sigma,
  interval = c(0.01, 0.40)
) {
  Y_data <- as.matrix(panel[, paste0("SVENY", sprintf("%02d", maturities))])
  TT <- nrow(Y_data)
  stopifnot(length(r_t) == TT, length(weights) == TT)
  stopifnot(all(c("fwd_02", "fwd_10") %in% names(panel)))

  bench_tau <- c(2, 10)

  # --- all fixed given (alpha, sigma): no dependence on mu or on t ---
  Ups_f <- Ups_fwd(bench_tau, alpha, tenor = 1) # 2 x 3
  C_f <- C_fwd(bench_tau, alpha, sigma, tenor = 1) # length 2
  Lsolve <- solve(Ups_f[, 2:3, drop = FALSE]) # inverse of the 2 x 2

  Ups_y <- Ups_yield(maturities, alpha) # N x 3
  C_y <- C_conv(maturities, alpha, sigma) # length N
  sum_y <- rowSums(Ups_y) # Ups %*% 1

  f_obs <- rbind(panel$fwd_02, panel$fwd_10) # 2 x T
  rt_term <- outer(Ups_f[, 1], r_t) # 2 x T

  # Factors from the two benchmark forwards, all dates at once
  extract <- function(mu) {
    const <- mu * (1 - rowSums(Ups_f)) - C_f # length 2
    t(Lsolve %*% (f_obs - const - rt_term)) # T x 2 = (m, l)
  }

  # Model yields in cascade coordinates
  model_yields <- function(mu, ml) {
    x <- cbind(r_t, ml) # T x 3 = (r, m, l)
    sweep(x %*% t(Ups_y), 2, mu * (1 - sum_y) - C_y, "+")
  }

  objective <- function(mu) {
    resid <- Y_data - model_yields(mu, extract(mu))
    sum(weights * rowSums(resid^2))
  }

  # Coarse grid first, then refine: optimize() is only reliable on a
  # bracketed minimum, and this objective can be very flat in mu.
  grid <- seq(interval[1], interval[2], length.out = 40)
  obj_grid <- vapply(grid, objective, numeric(1))
  i <- which.min(obj_grid)
  res <- optimize(
    objective,
    interval = c(grid[max(1, i - 1)], grid[min(length(grid), i + 1)]),
    tol = 1e-8
  )

  mu_hat <- res$minimum
  ml_hat <- extract(mu_hat)
  Y_fit <- model_yields(mu_hat, ml_hat)
  resid <- Y_data - Y_fit

  colnames(ml_hat) <- c("m", "l")
  rmse_bp <- sqrt(colSums(weights * resid^2) / sum(weights)) * 1e4
  names(rmse_bp) <- paste0("m", maturities)

  list(
    mu = mu_hat,
    factors = cbind(r = r_t, ml_hat), # T x 3, ready for stage 4
    yields_fit = Y_fit,
    resid = resid,
    rmse_bp = rmse_bp, # fit quality by maturity
    mu_grid = data.frame(mu = grid, obj = obj_grid), # check flatness
    objective = res$objective
  )
}


#' Full three-stage estimation
#'
#' Run all three stages sequentially.
#'
#' @param panel Data frame from build_panel()
#' @param weights Weight vector from build_panel()
#' @param maturities Maturity vector from build_panel()
#' @param start_alpha Initial alpha guess
#' @param start_sigma Initial sigma guess
#' @param start_mu Initial mu guess
#' @return List with estimated parameters
#'
estimate_all <- function(
  panel,
  weights,
  maturities,
  start_alpha = c(1.0, 0.6, 0.02),
  start_sigma = c(0.01, 0.01, 0.2)
) {
  r_t <- panel$r_t

  cat("Stage 1: Estimating alpha from regression slopes...\n")
  res1 <- stage1_alpha(panel, r_t, weights, maturities, start_alpha)
  alpha_hat <- res1$alpha
  cat(
    "  alpha_r =",
    round(alpha_hat[1], 4),
    ", alpha_m =",
    round(alpha_hat[2], 4),
    ", alpha_l =",
    round(alpha_hat[3], 4),
    "\n"
  )
  cat("  Half-lives:", round(log(2) / alpha_hat, 2), "years\n\n")

  cat("Stage 2: Estimating sigma from yield volatility...\n")
  res2 <- stage2_sigma(
    panel,
    r_t,
    weights,
    maturities,
    alpha_hat,
    start_sigma,
    scale = 'vol'
  )
  sigma_hat <- res2$sigma
  cat(
    "  sigma_m =",
    round(sigma_hat[1], 5),
    "(",
    round(sigma_hat[1] * 10000, 1),
    "bp )\n"
  )
  cat(
    "  sigma_l =",
    round(sigma_hat[2], 5),
    "(",
    round(sigma_hat[2] * 10000, 1),
    "bp )\n"
  )
  cat("  rho =", round(sigma_hat[3], 3), "\n\n")

  cat("Stage 3: Estimating mu from yield levels...\n")
  mu_hat <- stage3_mu(
    panel,
    r_t,
    weights,
    maturities,
    alpha_hat,
    sigma_hat
  )
  cat("  mu =", round(mu_hat$mu, 5), "(", round(mu_hat$mu * 100, 2), "% )\n\n")

  # Package results
  params <- list(
    a_r = alpha_hat[1],
    a_m = alpha_hat[2],
    a_l = alpha_hat[3],
    s_m = sigma_hat[1],
    s_l = sigma_hat[2],
    rho = sigma_hat[3],
    mu = mu_hat
  )

  return(list(
    params = params,
    stage1 = res1,
    stage2 = res2
  ))
}
