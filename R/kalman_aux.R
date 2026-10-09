# ==============================================================================
# kalman_aux.R - Kalman Filter Auxiliary Functions
# ==============================================================================
#
# Helper functions for Kalman filtering of the Gauss+ term structure model.
# Includes:
# - Discretization of continuous-time model (Van Loan method)
# - State and measurement equation setup
# - Filter initialization (warm, steady-state, diffuse)
# - Main Kalman filter recursion
# ==============================================================================

#' Van Loan discretization method
#'
#' Discretize continuous-time linear system using the Van Loan matrix
#' exponential method. Converts dx = -A x dt + Sig dW to discrete time
#' x_t = Phi x_{t-1} + nu_t with E[nu nu'] = Q.
#'
#' @param A Matrix n x n, drift matrix (with sign convention -A)
#' @param SigSig Matrix n x n, continuous-time diffusion covariance Sigma * Sigma'
#' @param dt Scalar, time step (in years, typically 1/252 for daily)
#'
#' @return List with components:
#'   \describe{
#'     \item{Phi}{Matrix n x n, discrete-time state transition matrix}
#'     \item{Q}{Matrix n x n, discrete-time innovation covariance}
#'   }
#'
van_loan_Q <- function(A, SigSig, dt) {
  n <- nrow(A)
  Z <- matrix(0, n, n)

  M <- dt * rbind(cbind(-A, SigSig), cbind(Z, t(A)))

  E <- expm::expm(M)

  Phi <- t(E[(n + 1):(2 * n), (n + 1):(2 * n)]) # transpose of bottom-right block
  Q <- Phi %*% E[1:n, (n + 1):(2 * n)] # Phi times top-right block
  Q <- (Q + t(Q)) / 2 # kill tiny asymmetries

  list(Phi = Phi, Q = Q)
}

#' Discretize Gauss+ model for Kalman filter state equation
#'
#' Turn model parameters into the discrete-time transition equation
#' z_t = c + Phi z_{t-1} + nu_t for the (m, l) factors. The short rate r
#' is observed and excluded from the state vector.
#'
#' @param alpha Vector length 3, mean reversion speeds c(a_r, a_m, a_l)
#' @param sigma Vector length 3, volatility parameters c(s_m, s_l, rho)
#' @param mu Scalar, long-run mean
#' @param dt Scalar, time step in years (default 1/252 for daily)
#'
#' @return List with components:
#'   \describe{
#'     \item{Phi}{Matrix 2 x 2, discrete-time state transition matrix}
#'     \item{Q}{Matrix 2 x 2, discrete-time innovation covariance}
#'     \item{c}{Vector length 2, transition intercept}
#'   }
#'
kf_discretise <- function(alpha, sigma, mu, dt = 1 / 252) {
  K2 <- K_matrix(alpha)[2:3, 2:3]
  Sig <- Omega(sigma)[2:3, 1:2]
  vl <- van_loan_Q(-K2, Sig %*% t(Sig), dt)
  c <- (diag(2) - vl$Phi) %*% c(mu, mu)

  return(list(Phi = vl$Phi, Q = vl$Q, c = c))
}

#' Build Kalman filter measurement equation
#'
#' Build the constant and loading matrix for the measurement equation
#' yhat_t = a + B z_t + epsilon_t, where yhat are netted yields (with
#' short rate contribution removed) and z_t = (m_t, l_t).
#'
#' @param maturities Integer vector length N, maturities in years
#' @param alpha Vector length 3, mean reversion speeds c(a_r, a_m, a_l)
#' @param sigma Vector length 3, volatility parameters c(s_m, s_l, rho)
#' @param mu Scalar, long-run mean
#'
#' @return List with components:
#'   \describe{
#'     \item{a}{Vector length N, measurement intercept}
#'     \item{B}{Matrix N x 2, measurement loadings on (m, l)}
#'   }
#'
kf_measurement <- function(maturities, alpha, sigma, mu) {
  Ups <- Ups_yield(maturities, alpha) # N x 3
  a <- mu * (1 - rowSums(Ups)) - C_conv(maturities, alpha, sigma)
  B <- Ups[, 2:3, drop = FALSE]

  return(list(a = a, B = B))
}

#' Prepare netted yield data for Kalman filter
#'
#' Produce the observation matrix the filter consumes: observed yields with
#' the short-rate contribution removed, so the measurement intercept is
#' time-invariant. This is the "netting" step that makes the measurement
#' equation depend only on (m, l).
#'
#' @param panel Data frame with yield columns SVENY01, SVENY02, etc.
#' @param r_t Vector length T, observed short rate
#' @param maturities Integer vector length N, maturities in years
#' @param alpha_r Scalar, short rate mean reversion speed
#'
#' @return Matrix T x N, netted yields Yhat_t = Y_t - r_t * B_1(tau)
#'
kf_data <- function(panel, r_t, maturities, alpha_r) {
  stopifnot(length(r_t) == nrow(panel))
  Y <- as.matrix(panel[, paste0("SVENY", sprintf("%02d", maturities))])
  Yhat <- Y - outer(r_t, short_rate_loading(maturities, alpha_r))
  return(Yhat)
}


#' Initialize Kalman filter state and covariance
#'
#' Provide initial state estimate z0 and covariance P0 for the Kalman filter.
#' Three methods are supported: warm start (use known initial factors), steady
#' state (solve for unconditional covariance), or diffuse (large initial
#' uncertainty).
#'
#' @param Phi Matrix 2 x 2, discrete-time state transition matrix
#' @param Q Matrix 2 x 2, discrete-time innovation covariance
#' @param mu Scalar, long-run mean
#' @param z0 Vector length 2, initial state for warm start (default NULL)
#' @param method Character, one of "warm", "steady", or "diffuse"
#'
#' @return List with components:
#'   \describe{
#'     \item{P0}{Matrix 2 x 2, initial state covariance}
#'     \item{z0}{Vector length 2, initial state estimate}
#'   }
#'
kf_init <- function(
  Phi,
  Q,
  mu,
  z0 = NULL,
  method = c("warm", "steady", "diffuse")
) {
  if (!missing(method) && length(method) != 1) {
    stop(
      "`method` must be exactly one of: 'warm', 'steady', 'diffuse'.",
      call. = FALSE
    )
  }
  method <- match.arg(method)

  if (method == "warm") {
    if (is.null(z0[c('m', 'l')])) {
      stop("`z0` is required when method = 'warm'.", call. = FALSE)
    }
    if (length(z0) != 2) {
      stop("`z0` is required only for m and l.", call. = FALSE)
    }
  }

  P0 <- switch(
    method,
    steady = steady_state_solve(Q, Phi),
    diffuse = diffuse_P0(),
    warm = warm_P0()
  )

  z0 <- switch(
    method,
    steady = c(mu, mu),
    diffuse = c(mu, mu),
    warm = z0
  )

  return(list(P0 = P0, z0 = z0))
}

#' Warm start initial covariance
#'
#' Return a tight initial covariance for warm start initialization.
#' Assumes factors are known with high confidence.
#'
#' @return Matrix 2 x 2, initial covariance (1e-4 * I)
#'
warm_P0 <- function() {
  return((0.01)^2 * diag(2))
}

#' Diffuse initial covariance
#'
#' Return a large initial covariance for diffuse initialization.
#' Represents high initial uncertainty about the state.
#'
#' @param kappa Scalar, scale factor for initial variance (default 10)
#'
#' @return Matrix 2 x 2, initial covariance (kappa * I)
#'
diffuse_P0 <- function(kappa = 10) {
  return(kappa * diag(2))
}

#' Solve for steady-state covariance
#'
#' Solve the discrete Lyapunov equation P = Phi P Phi' + Q iteratively to
#' find the unconditional covariance. Used for steady-state initialization.
#'
#' @param Q Matrix 2 x 2, innovation covariance
#' @param Phi Matrix 2 x 2, state transition matrix
#' @param max_iter Integer, maximum iterations (default 10000)
#' @param tol Scalar, convergence tolerance (default 1e-10)
#'
#' @return Matrix 2 x 2, steady-state covariance P
#'
steady_state_solve <- function(Q, Phi, max_iter = 1e4, tol = 1e-12) {
  # Direct solve of vec(P) = (I - Phi (x) Phi)^-1 vec(Q). The iterative
  # recursion needs ~1e5+ steps when a persistent factor has Phi ~ 1 - 6e-5.
  n <- nrow(Phi)
  P <- tryCatch(
    matrix(solve(diag(n^2) - kronecker(Phi, Phi), as.vector(Q)), n, n),
    error = function(e) NULL
  )

  if (is.null(P) || any(!is.finite(P))) {
    warning("Steady-state solve failed; using warm start P0")
    return(warm_P0())
  }

  (P + t(P)) / 2
}

#' Kalman filter recursion
#'
#' Given everything fixed, produce filtered states, their covariances, and the
#' likelihood. Implements the standard Kalman filter with predict-update steps.
#'
#' @param Yhat Matrix T x N, netted observations (yields with short rate removed)
#' @param Phi Matrix 2 x 2, state transition matrix
#' @param Q Matrix 2 x 2, state innovation covariance
#' @param c Vector length 2, transition intercept
#' @param a Vector length N, measurement intercept
#' @param B Matrix N x 2, measurement loadings
#' @param H Matrix N x N, measurement error covariance
#' @param z0 Vector length 2, initial state estimate
#' @param P0 Matrix 2 x 2, initial state covariance
#' @param store_innov_var Logical, whether to store N x N x T innovation
#'   covariance array (default FALSE to save memory)
#'
#' @return List with components:
#'   \describe{
#'     \item{z_filt}{Matrix T x 2, filtered state estimates E[z_t | Y_1:t]}
#'     \item{P_filt}{Array 2 x 2 x T, filtered state covariances}
#'     \item{z_pred}{Matrix T x 2, one-step-ahead predictions E[z_t | Y_1:(t-1)]}
#'     \item{P_pred}{Array 2 x 2 x T, predicted state covariances}
#'     \item{innov}{Matrix T x N, prediction errors (innovations) v_t}
#'     \item{innov_var}{Array N x N x T, innovation covariances (if requested)}
#'     \item{loglik}{Scalar, total log-likelihood}
#'   }
#'
kalman_filter <- function(Yhat, Phi, Q, c, a, B, H, z0, P0,
                          store_innov_var = FALSE) {
  TT <- nrow(Yhat)
  N <- ncol(Yhat)

  # Storage
  z_filt <- matrix(NA, TT, 2)
  P_filt <- array(NA, dim = c(2, 2, TT))
  z_pred <- matrix(NA, TT, 2)
  P_pred <- array(NA, dim = c(2, 2, TT))
  innov <- matrix(NA, TT, N)
  if (store_innov_var) {
    innov_var <- array(NA, dim = c(N, N, TT))
  }

  loglik <- 0
  z <- z0
  P <- P0

  for (t in 1:TT) {
    # --- Predict step ---
    z_p <- c + Phi %*% z
    P_p <- Phi %*% P %*% t(Phi) + Q

    # Store predictions
    z_pred[t, ] <- z_p
    P_pred[, , t] <- P_p

    # --- Innovation step ---
    v <- Yhat[t, ] - a - B %*% z_p
    F_t <- B %*% P_p %*% t(B) + H

    # Store innovation
    innov[t, ] <- v
    if (store_innov_var) {
      innov_var[, , t] <- F_t
    }

    # --- Likelihood contribution ---
    # Use Cholesky for numerical stability
    L <- chol(F_t)  # upper triangular, F = t(L) %*% L
    logdet <- 2 * sum(log(diag(L)))
    Finv_v <- backsolve(L, forwardsolve(t(L), v))
    loglik <- loglik - 0.5 * (N * log(2 * pi) + logdet + sum(v * Finv_v))

    # --- Kalman gain ---
    K <- P_p %*% t(B) %*% solve(F_t)

    # --- Update step ---
    z <- z_p + K %*% v
    P <- P_p - K %*% B %*% P_p

    # Symmetrise to prevent accumulation of numerical asymmetry
    P <- (P + t(P)) / 2

    # Store filtered estimates
    z_filt[t, ] <- z
    P_filt[, , t] <- P
  }

  result <- list(
    z_filt = z_filt,
    P_filt = P_filt,
    z_pred = z_pred,
    P_pred = P_pred,
    innov = innov,
    loglik = loglik
  )

  if (store_innov_var) {
    result$innov_var <- innov_var
  }

  return(result)
}

#' Scalar log-likelihood as a function of measurement error scale h
#'
#' Wrapper around kalman_filter() with H = h^2 * I, so h can be profiled or
#' optimised. Returns -Inf for invalid h or a non-finite / failed filter run.
#'
#' @param h Scalar, measurement error standard deviation (must be > 0)
#' @param Yhat Matrix T x N, netted observations
#' @param Phi,Q,c Transition matrix, innovation covariance, intercept
#' @param a,B Measurement intercept and loadings
#' @param z0,P0 Initial state and covariance
#'
#' @return Scalar log-likelihood (-Inf if h <= 0 or the filter fails)
#'
kf_loglik <- function(h, Yhat, Phi, Q, c, a, B, z0, P0) {
  if (!is.finite(h) || h <= 0) {
    return(-Inf)
  }
  H <- h^2 * diag(ncol(Yhat))
  ll <- tryCatch(
    kalman_filter(Yhat, Phi, Q, c, a, B, H, z0, P0)$loglik,
    error = function(e) -Inf
  )
  if (!is.finite(ll)) {
    return(-Inf)
  }
  ll
}

#' Profile the log-likelihood over a grid of h
#'
#' Evaluate kf_loglik() across h_grid to see the shape of the likelihood in h
#' rather than trusting a single optimiser result. A clear interior peak means
#' h is identified; a plateau or monotone drift to a boundary means it is not,
#' and h should be fixed by judgement. Expect the fitted h to exceed genuine
#' quote noise: with alpha, sigma, mu fixed it absorbs all model misfit.
#'
#' @param h_grid Numeric vector of h values, e.g. 10^seq(-5, -2, length.out = 40)
#' @param ... Passed to kf_loglik(): Yhat, Phi, Q, c, a, B, z0, P0
#' @param plot Logical, plot loglik against log10(h) (default TRUE)
#'
#' @return Data frame with columns h and loglik
#'
profile_h <- function(h_grid = 10^seq(-5, -2, length.out = 40), ...,
                      plot = TRUE) {
  ll <- vapply(h_grid, kf_loglik, numeric(1), ...)
  out <- data.frame(h = h_grid, loglik = ll)

  if (plot) {
    ok <- is.finite(out$loglik)
    plot(log10(out$h[ok]), out$loglik[ok], type = "b", pch = 16,
         xlab = "log10(h)", ylab = "log-likelihood",
         main = "Profile likelihood in h")
    best <- which.max(out$loglik)
    abline(v = log10(out$h[best]), lty = 2, col = "grey40")
  }

  out
}

#' Extract latent factors (m_t, l_t) with the Kalman filter
#'
#' Alternative to extract_factors() / stage3_mu(): instead of fitting the 2y and
#' 10y forwards exactly, filter all N yields with measurement error
#' H = h^2 * I and return the filtered (m, l). Parameters are held fixed.
#'
#' If h is NULL it is chosen by maximising the profile likelihood over h_grid;
#' a warning is issued if the optimum sits on the edge of the grid (h not
#' identified - fix it by judgement instead).
#'
#' @param panel Data frame with yield columns SVENY.. and r_t
#' @param maturities Integer vector of maturities in years
#' @param alpha Vector c(a_r, a_m, a_l)
#' @param sigma Vector c(s_m, s_l, rho)
#' @param mu Scalar long-run mean
#' @param h Scalar measurement error sd, or NULL to profile
#' @param h_grid Grid used when h is NULL
#' @param init Initialisation method for kf_init(): "steady" or "diffuse"
#' @param dt Time step in years
#'
#' @return List with factors (T x 3 matrix, columns r, m, l), h, loglik,
#'   profile (data frame or NULL) and the full filter output
#'
extract_factors_kalman <- function(
  panel,
  maturities,
  alpha,
  sigma,
  mu,
  h = NULL,
  h_grid = 10^seq(-5, -2, length.out = 40),
  init = c("steady", "diffuse"),
  dt = 1 / 252
) {
  init <- match.arg(init)
  r_t <- panel$r_t

  d <- kf_discretise(alpha, sigma, mu, dt)
  meas <- kf_measurement(maturities, alpha, sigma, mu)
  Yhat <- kf_data(panel, r_t, maturities, alpha[1])
  i0 <- kf_init(d$Phi, d$Q, mu, method = init)

  prof <- NULL
  if (is.null(h)) {
    prof <- profile_h(
      h_grid,
      Yhat = Yhat, Phi = d$Phi, Q = d$Q, c = d$c,
      a = meas$a, B = meas$B, z0 = i0$z0, P0 = i0$P0,
      plot = FALSE
    )
    best <- which.max(prof$loglik)
    h <- prof$h[best]
    if (best == 1 || best == length(h_grid)) {
      warning(
        "extract_factors_kalman: profile likelihood peaks at the edge of ",
        "h_grid (h = ", signif(h, 3), "); h is not identified."
      )
    }
  }

  fit <- kalman_filter(
    Yhat, d$Phi, d$Q, d$c, meas$a, meas$B,
    H = h^2 * diag(length(maturities)), i0$z0, i0$P0
  )

  factors <- cbind(r = r_t, m = fit$z_filt[, 1], l = fit$z_filt[, 2])

  list(factors = factors, h = h, loglik = fit$loglik, profile = prof, filter = fit)
}
