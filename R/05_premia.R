# ==============================================================================
# 05_premia.R - Risk Premia Extraction
# ==============================================================================

source("R/02_model.R")

#' Unit risk premium loading
#'
#' Risk premium per unit of price of risk lambda.
#' RP_unit(tau) = tau * Υ_3(tau) * sigma_l
#'
#' Only the long factor is priced in this specification.
#'
#' @param tau Maturity (in years)
#' @param alpha Parameter vector c(a_r, a_m, a_l)
#' @param sigma Parameter vector c(s_m, s_l, rho)
#' @return Numeric vector of length(tau)
#'
rp_unit <- function(tau, alpha, sigma) {
  if (is.list(sigma)) {
    s_l <- sigma$s_l
  } else {
    s_l <- sigma[2]
  }

  # Instantaneous forward loading on long factor (column 3)
  Ups_inst <- Ups_fwd_inst(tau, alpha)
  Ups_3 <- Ups_inst[, 3]

  RP <- tau * Ups_3 * s_l

  return(RP)
}


#' Extract price of risk lambda_t
#'
#' Under flat expectations beyond some long maturity, the difference between
#' two long forwards is pure risk premium:
#'
#' lambda_t = (f_t(tau') - f_t(tau)) / (RP_unit(tau') - RP_unit(tau))
#'
#' Default: tau = 14, tau' = 15 (one-year forwards at these start times)
#'
#' @param fwd_tau Forward rate at tau
#' @param fwd_tau_prime Forward rate at tau'
#' @param tau First maturity (default 14)
#' @param tau_prime Second maturity (default 15)
#' @param alpha Parameter vector
#' @param sigma Parameter vector
#' @return Numeric vector of lambda values (one per time observation)
#'
solve_lambda <- function(fwd_tau, fwd_tau_prime, tau = 14, tau_prime = 15,
                         alpha, sigma) {

  RP_tau <- rp_unit(tau, alpha, sigma)
  RP_tau_prime <- rp_unit(tau_prime, alpha, sigma)

  lambda_t <- (fwd_tau_prime - fwd_tau) / (RP_tau_prime - RP_tau)

  return(lambda_t)
}


#' Compute expected future short rate
#'
#' Decompose forward rate into expectation, risk premium, and convexity:
#' E_t[r_tau] = f_t(tau) - lambda_t * RP_unit(tau) + C'(tau)
#'
#' Uses CORRECTED formula from spec (A9.29): convexity term is C', not C'/Δτ
#'
#' @param fwd_t Forward rate at time t, starting at tau
#' @param lambda_t Price of risk at time t
#' @param tau Forward start time
#' @param tenor Forward tenor (default 1 year)
#' @param alpha Parameter vector
#' @param sigma Parameter vector
#' @return Expected short rate tau years forward
#'
expected_rate <- function(fwd_t, lambda_t, tau, tenor = 1, alpha, sigma) {

  RP <- rp_unit(tau, alpha, sigma)
  C_prime <- C_fwd(tau, alpha, sigma, tenor)

  E_r <- fwd_t - lambda_t * RP + C_prime

  return(E_r)
}
