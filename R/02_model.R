# ==============================================================================
# 02_model.R - Gauss+ Model Functions
# ==============================================================================
#
# Closed-form solutions for the three-factor Gauss+ model from
# Tuckman & Serrat Chapter 9, Appendix A9.2
#
# All functions use the corrected formulas
# ==============================================================================

#' Mean reversion matrix K
#'
#' Upper triangular matrix defining the cascade dynamics.
#' Eigenvalues are (alpha_r, alpha_m, alpha_l) on the diagonal.
#'
#' @param alpha Named list or vector c(alpha_r, alpha_m, alpha_l)
#' @return 3x3 matrix K
#'
K_matrix <- function(alpha) {
  if (is.list(alpha)) {
    a_r <- alpha$a_r
    a_m <- alpha$a_m
    a_l <- alpha$a_l
  } else {
    a_r <- alpha[1]
    a_m <- alpha[2]
    a_l <- alpha[3]
  }

  K <- matrix(
    c(
      a_r,
      -a_r,
      0,
      0,
      a_m,
      -a_m,
      0,
      0,
      a_l
    ),
    nrow = 3,
    byrow = TRUE
  )

  return(K)
}


#' Eigenvector matrix A(α)
#'
#' Columns are eigenvectors of K, normalized so first row is (1,1,1).
#' This ensures r_t = μ + sum(X_t).
#'
#' Uses CORRECTED formula from spec §2.1:
#' Entry (3,3) has denominator α_r·α_m, not just α_r as printed in book.
#'
#' @param alpha Named list or vector c(alpha_r, alpha_m, alpha_l)
#' @return 3x3 matrix A
#'
A_mat <- function(alpha) {
  if (is.list(alpha)) {
    a_r <- alpha$a_r
    a_m <- alpha$a_m
    a_l <- alpha$a_l
  } else {
    a_r <- alpha[1]
    a_m <- alpha[2]
    a_l <- alpha[3]
  }

  A <- matrix(
    c(
      1,
      1,
      1,
      0,
      (a_r - a_m) / a_r,
      (a_r - a_l) / a_r,
      0,
      0,
      (a_r - a_l) * (a_m - a_l) / (a_r * a_m)
    ),
    nrow = 3,
    byrow = TRUE
  )

  return(A)
}


#' Inverse of A(α)
#'
#' Uses CORRECTED formula from spec §2.2:
#' Entry (2,3) has a MINUS sign (book omits it).
#'
#' @param alpha Named list or vector c(alpha_r, alpha_m, alpha_l)
#' @return 3x3 matrix A^{-1}
#'
A_inv <- function(alpha) {
  if (is.list(alpha)) {
    a_r <- alpha$a_r
    a_m <- alpha$a_m
    a_l <- alpha$a_l
  } else {
    a_r <- alpha[1]
    a_m <- alpha[2]
    a_l <- alpha[3]
  }

  # For robustness, could also use: solve(A_mat(alpha))
  # Closed form is provided here for transparency and speed

  A_inv <- matrix(
    c(
      1,
      -a_r / (a_r - a_m),
      a_r * a_m / ((a_r - a_m) * (a_r - a_l)),
      0,
      a_r / (a_r - a_m),
      -a_r * a_m / ((a_r - a_m) * (a_m - a_l)),
      0,
      0,
      a_r * a_m / ((a_r - a_l) * (a_m - a_l))
    ),
    nrow = 3,
    byrow = TRUE
  )

  return(A_inv)
}


#' Diffusion matrix Ω(σ)
#'
#' 3x3 matrix (but only 2 independent shocks - third column is zero).
#' Row 1 is zero because short rate r has no diffusion.
#'
#' @param sigma Named list or vector c(sigma_m, sigma_l, rho)
#' @return 3x3 matrix Omega
#'
Omega <- function(sigma) {
  if (is.list(sigma)) {
    s_m <- sigma$s_m
    s_l <- sigma$s_l
    rho <- sigma$rho
  } else {
    s_m <- sigma[1]
    s_l <- sigma[2]
    rho <- sigma[3]
  }

  Om <- matrix(
    c(
      0,
      0,
      0,
      rho * s_m,
      sqrt(1 - rho^2) * s_m,
      0,
      s_l,
      0,
      0
    ),
    nrow = 3,
    byrow = TRUE
  )

  return(Om)
}


#' Loading function B(τ, α)
#'
#' "Averaged" or "normalized" loading: B_i(τ) = (1 - exp(-α_i τ)) / (α_i τ)
#' Used in yield formulas (which are averages of instantaneous forwards).
#'
#' Numerically stable for small α_i τ using expm1.
#'
#' @param tau Numeric vector of maturities (in years)
#' @param alpha Named list or vector c(alpha_r, alpha_m, alpha_l)
#' @return Matrix of size length(tau) x 3, each column is B_i(tau)
#'
B <- function(tau, alpha) {
  if (is.list(alpha)) {
    a_r <- alpha$a_r
    a_m <- alpha$a_m
    a_l <- alpha$a_l
  } else {
    a_r <- alpha[1]
    a_m <- alpha[2]
    a_l <- alpha[3]
  }

  # Vectorize over tau
  # expm1(x) = exp(x) - 1, more accurate for small x
  # (1 - exp(-a*tau)) = -expm1(-a*tau)
  B_r <- -expm1(-a_r * tau) / (a_r * tau)
  B_m <- -expm1(-a_m * tau) / (a_m * tau)
  B_l <- -expm1(-a_l * tau) / (a_l * tau)

  return(cbind(B_r, B_m, B_l))
}


#' Unnormalized loading Btilde(τ, α)
#'
#' Btilde_i(τ) = (1 - exp(-α_i τ)) / α_i = τ · B_i(τ)
#' Used in forward rate formulas.
#'
#' @param tau Numeric vector of maturities (in years)
#' @param alpha Named list or vector c(alpha_r, alpha_m, alpha_l)
#' @return Matrix of size length(tau) x 3, each column is Btilde_i(tau)
#'
Btilde <- function(tau, alpha) {
  if (is.list(alpha)) {
    a_r <- alpha$a_r
    a_m <- alpha$a_m
    a_l <- alpha$a_l
  } else {
    a_r <- alpha[1]
    a_m <- alpha[2]
    a_l <- alpha[3]
  }

  Bt_r <- -expm1(-a_r * tau) / a_r
  Bt_m <- -expm1(-a_m * tau) / a_m
  Bt_l <- -expm1(-a_l * tau) / a_l

  return(cbind(Bt_r, Bt_m, Bt_l))
}


#' Yield loading Υ(τ, α)
#'
#' Row vector that maps state X to yield:
#' y(τ) = μ - C(τ) + Υ(τ) · X
#'
#' Formula: Υ(τ, α) = B(τ, α) %*% A(α)^{-1}
#'
#' @param tau Numeric vector of maturities (in years)
#' @param alpha Named list or vector c(alpha_r, alpha_m, alpha_l)
#' @return Matrix of size length(tau) x 3, each row is Υ(tau)
#'
Ups_yield <- function(tau, alpha) {
  B_mat <- B(tau, alpha) # length(tau) x 3
  A_inv_mat <- A_inv(alpha) # 3 x 3

  # Matrix multiply: each row of B times A_inv
  Ups <- B_mat %*% A_inv_mat # length(tau) x 3

  return(Ups)
}


#' Convexity term C(τ, α, σ)
#'
#' Maturity-dependent constant in yield formula. Depends on parameters but
#' not on state, so it differences away in changes.
#'
#' Uses CORRECTED formula from spec §2.5:
#' - Last term in bracket is ADDED (book prints it subtracted)
#' - σ_ij from S = A^{-1} Ω Ω' (A^{-1})', with transpose on last factor
#'
#' @param tau Numeric vector of maturities (in years)
#' @param alpha Named list or vector c(alpha_r, alpha_m, alpha_l)
#' @param sigma Named list or vector c(sigma_m, sigma_l, rho)
#' @return Numeric vector of length(tau)
#'
C_conv <- function(tau, alpha, sigma) {
  if (is.list(alpha)) {
    a_r <- alpha$a_r
    a_m <- alpha$a_m
    a_l <- alpha$a_l
  } else {
    a_r <- alpha[1]
    a_m <- alpha[2]
    a_l <- alpha[3]
  }

  # Compute S = A^{-1} Ω Ω' (A^{-1})'
  A_inv_mat <- A_inv(alpha)
  Om <- Omega(sigma)
  S <- A_inv_mat %*% Om %*% t(Om) %*% t(A_inv_mat)

  alpha_vec <- c(a_r, a_m, a_l)

  # Vectorize over tau
  C_out <- numeric(length(tau))

  for (k in seq_along(tau)) {
    tau_k <- tau[k]
    B_k <- B(tau_k, alpha) # 1 x 3

    sum_val <- 0
    for (i in 1:3) {
      for (j in 1:3) {
        a_i <- alpha_vec[i]
        a_j <- alpha_vec[j]
        s_ij <- S[i, j]
        B_i <- B_k[i]
        B_j <- B_k[j]

        # CORRECTED: last term is ADDED, not subtracted
        bracket <- 1 -
          B_i -
          B_j +
          (1 - exp(-(a_i + a_j) * tau_k)) / ((a_i + a_j) * tau_k)
        sum_val <- sum_val + s_ij / (2 * a_i * a_j) * bracket
      }
    }
    C_out[k] <- sum_val
  }

  return(C_out)
}


#' Forward rate loading Υ'(τ, α, τ')
#'
#' Row vector that maps state to forward rate starting at τ with tenor τ'.
#'
#' Uses CORRECTED formula from spec §2.6:
#' Υ'(τ, α, τ') = (Btilde(τ+τ') - Btilde(τ)) A^{-1} / τ'
#'
#' Book's (A9.16) omits the τ weights and would give wrong loadings.
#'
#' @param tau Numeric vector of forward start times (in years)
#' @param alpha Named list or vector c(alpha_r, alpha_m, alpha_l)
#' @param tenor Forward tenor (default 1 year)
#' @return Matrix of size length(tau) x 3, each row is Υ'(tau)
#'
Ups_fwd <- function(tau, alpha, tenor = 1) {
  Bt_end <- Btilde(tau + tenor, alpha) # length(tau) x 3
  Bt_start <- Btilde(tau, alpha) # length(tau) x 3

  # Difference and divide by tenor
  Bt_diff <- (Bt_end - Bt_start) / tenor # length(tau) x 3

  # Multiply by A^{-1}
  A_inv_mat <- A_inv(alpha)
  Ups_f <- Bt_diff %*% A_inv_mat # length(tau) x 3

  return(Ups_f)
}


#' Instantaneous forward rate loading Υ'_inst(τ, α)
#'
#' Limiting case as τ' → 0.
#' Used for Figure 9.7 (instantaneous forward loadings vs maturity).
#'
#' Formula: [exp(-α_r τ), exp(-α_m τ), exp(-α_l τ)] A^{-1}
#'
#' @param tau Numeric vector of maturities (in years)
#' @param alpha Named list or vector c(alpha_r, alpha_m, alpha_l)
#' @return Matrix of size length(tau) x 3, each row is Υ'_inst(tau)
#'
Ups_fwd_inst <- function(tau, alpha) {
  if (is.list(alpha)) {
    a_r <- alpha$a_r
    a_m <- alpha$a_m
    a_l <- alpha$a_l
  } else {
    a_r <- alpha[1]
    a_m <- alpha[2]
    a_l <- alpha[3]
  }

  # Exponential decay for each factor
  exp_mat <- cbind(
    exp(-a_r * tau),
    exp(-a_m * tau),
    exp(-a_l * tau)
  )

  # Multiply by A^{-1}
  A_inv_mat <- A_inv(alpha)
  Ups_inst <- exp_mat %*% A_inv_mat

  return(Ups_inst)
}


#' Forward convexity term C'(τ, α, σ, τ')
#'
#' Convexity adjustment for forward rates.
#'
#' Uses CORRECTED formula from spec §2.6:
#' C'(τ, α, σ, τ') = ((τ+τ')·C(τ+τ') - τ·C(τ)) / τ'
#'
#' @param tau Numeric vector of forward start times (in years)
#' @param alpha Named list or vector c(alpha_r, alpha_m, alpha_l)
#' @param sigma Named list or vector c(sigma_m, sigma_l, rho)
#' @param tenor Forward tenor (default 1 year)
#' @return Numeric vector of length(tau)
#'
C_fwd <- function(tau, alpha, sigma, tenor = 1) {
  C_end <- C_conv(tau + tenor, alpha, sigma)
  C_start <- C_conv(tau, alpha, sigma)

  # CORRECTED: weight by maturity before differencing
  C_prime <- ((tau + tenor) * C_end - tau * C_start) / tenor

  return(C_prime)
}
