# ============================================================
# Laplace null model for nominal multinomial logit mixed model
# ============================================================
# Model:
#   eta_ik = alpha_k + u_ik, k = 1,...,q, q = C - 1
#   u_k ~ N(0, sigma_k * K), independent across logits
# Baseline category is the last category C.
#
# Input y can be either:
#   1) n x C one-hot matrix, as created by model.matrix(~ nominal - 1)
#   2) length-n category vector/factor
#
# Output:
#   alpha: q intercepts
#   vc: q variance components sigma_k
#   U: n x q matrix of random-effect modes
#   fit: optim object
#   nll: final Laplace negative log-likelihood
# ============================================================

nominal_make_Y <- function(y) {
  if (is.matrix(y) || is.data.frame(y)) {
    y <- as.matrix(y)
    if (!all(y %in% c(0, 1))) {
      stop("If y is a matrix, it should be a one-hot 0/1 category matrix.")
    }
    if (any(rowSums(y) != 1)) {
      stop("Each row of the one-hot y matrix should sum to 1.")
    }
    return(y)
  }
  y_fac <- factor(y)
  model.matrix(~ y_fac - 1)
}

.nominal_laplace_build_H <- function(PI, sigma, Kinv, ridge = 1e-6) {
  n <- nrow(PI)
  q <- ncol(PI)
  H <- matrix(0, n * q, n * q)

  # G^{-1}: category-stacked order, u = c(U[,1], U[,2], ...)
  for (k in seq_len(q)) {
    idx <- ((k - 1) * n + 1):(k * n)
    H[idx, idx] <- H[idx, idx] + Kinv / sigma[k]
  }

  # W_i = diag(pi_i) - pi_i pi_i'
  for (i in seq_len(n)) {
    pi_i <- PI[i, ]
    Wi <- diag(pi_i, q, q) - tcrossprod(pi_i)
    for (k1 in seq_len(q)) {
      for (k2 in seq_len(q)) {
        row_idx <- (k1 - 1) * n + i
        col_idx <- (k2 - 1) * n + i
        H[row_idx, col_idx] <- H[row_idx, col_idx] + Wi[k1, k2]
      }
    }
  }

  H + diag(ridge, n * q)
}

.nominal_laplace_mode_U <- function(Yq, alpha, sigma, Kinv,
                                    inner_maxiter = 30,
                                    inner_tol = 1e-6,
                                    step_factor = 0.5,
                                    ridge = 1e-6,
                                    verbose = FALSE) {
  n <- nrow(Yq)
  q <- ncol(Yq)
  U <- matrix(0, n, q)
  trace <- NULL

  for (iter in seq_len(inner_maxiter)) {
    ETA <- sweep(U, 2, alpha, "+")
    expETA <- exp(ETA)
    denom <- 1 + rowSums(expETA)
    PI <- expETA / denom

    GinvU <- matrix(0, n, q)
    for (k in seq_len(q)) {
      GinvU[, k] <- as.vector(Kinv %*% U[, k]) / sigma[k]
    }

    Score <- Yq - PI - GinvU
    score_vec <- as.vector(Score)  # column/category-stacked

    H <- .nominal_laplace_build_H(PI, sigma, Kinv, ridge = ridge)
    step <- solve(H, score_vec)
    U_new <- matrix(as.vector(U) + step_factor * step, n, q)

    diff <- max(abs(U_new - U))
    trace <- rbind(trace, c(iter = iter, diff = diff))
    if (verbose) cat("inner iter:", iter, "diff:", diff, "\n")

    U <- U_new
    if (diff < inner_tol) break
  }

  list(U = U, trace = as.data.frame(trace))
}

.nominal_laplace_nll <- function(Yq, alpha, sigma, Kinv, logdetK,
                                 inner_maxiter = 30,
                                 inner_tol = 1e-6,
                                 step_factor = 0.5,
                                 ridge = 1e-6,
                                 return_mode = FALSE,
                                 verbose = FALSE) {
  n <- nrow(Yq)
  q <- ncol(Yq)

  mode <- .nominal_laplace_mode_U(
    Yq = Yq,
    alpha = alpha,
    sigma = sigma,
    Kinv = Kinv,
    inner_maxiter = inner_maxiter,
    inner_tol = inner_tol,
    step_factor = step_factor,
    ridge = ridge,
    verbose = FALSE
  )
  U <- mode$U

  ETA <- sweep(U, 2, alpha, "+")
  expETA <- exp(ETA)
  denom <- 1 + rowSums(expETA)
  PI <- expETA / denom

  loglik <- sum(Yq * ETA) - sum(log(denom))

  GinvU <- matrix(0, n, q)
  for (k in seq_len(q)) {
    GinvU[, k] <- as.vector(Kinv %*% U[, k]) / sigma[k]
  }
  quad <- sum(U * GinvU)

  logdetG <- q * logdetK + n * sum(log(sigma))

  H <- .nominal_laplace_build_H(PI, sigma, Kinv, ridge = ridge)
  logdetH <- as.numeric(determinant(H, logarithm = TRUE)$modulus)

  nll <- -loglik + 0.5 * quad + 0.5 * logdetG + 0.5 * logdetH

  if (verbose) {
    cat("alpha =", round(alpha, 3),
        " sigma =", round(sigma, 3),
        " nll =", round(nll, 3), "\n")
  }

  if (return_mode) {
    return(list(nll = as.numeric(nll), U = U, PI = PI, mode_trace = mode$trace))
  }
  as.numeric(nll)
}

fit_nominal_laplace_null <- function(y,
                                     kk,
                                     alpha0 = NULL,
                                     sigma0 = NULL,
                                     optim_method = "L-BFGS-B",
                                     alpha_lower = -5,
                                     alpha_upper = 5,
                                     sigma_lower = 1e-4,
                                     sigma_upper = 5,
                                     inner_maxiter = 30,
                                     inner_tol = 1e-6,
                                     step_factor = 0.5,
                                     ridge = 1e-6,
                                     jitter = 1e-6,
                                     verbose = TRUE) {
  Yfull <- nominal_make_Y(y)
  n <- nrow(Yfull)
  c <- ncol(Yfull)
  q <- c - 1

  if (q < 1) stop("Nominal y should contain at least two categories.")
  if (!all(dim(kk) == c(n, n))) stop("kk should be an n x n kinship matrix matching y.")

  kk <- as.matrix(kk)
  kk <- (kk + t(kk)) / 2
  kk <- kk + diag(jitter, n)

  Kinv <- solve(kk)
  logdetK <- as.numeric(determinant(kk, logarithm = TRUE)$modulus)

  Yq <- Yfull[, seq_len(q), drop = FALSE]

  if (is.null(alpha0)) {
    p_obs <- colMeans(Yfull)
    p_obs <- pmax(p_obs, 1e-6)
    alpha0 <- log(p_obs[seq_len(q)] / p_obs[c])
  }
  if (is.null(sigma0)) sigma0 <- rep(0.2, q)

  if (length(alpha0) != q) stop("alpha0 should have length C - 1.")
  if (length(sigma0) != q) stop("sigma0 should have length C - 1.")
  if (any(sigma0 <= 0)) stop("sigma0 must be positive.")

  par0 <- c(alpha0, log(sigma0))

  obj <- function(par) {
    alpha <- par[seq_len(q)]
    sigma <- exp(par[q + seq_len(q)])

    .nominal_laplace_nll(
      Yq = Yq,
      alpha = alpha,
      sigma = sigma,
      Kinv = Kinv,
      logdetK = logdetK,
      inner_maxiter = inner_maxiter,
      inner_tol = inner_tol,
      step_factor = step_factor,
      ridge = ridge,
      verbose = verbose
    )
  }

  lower <- c(rep(alpha_lower, q), rep(log(sigma_lower), q))
  upper <- c(rep(alpha_upper, q), rep(log(sigma_upper), q))

  start_time <- Sys.time()
  fit <- optim(
    par = par0,
    fn = obj,
    method = optim_method,
    lower = lower,
    upper = upper
  )
  elapsed_time <- Sys.time() - start_time

  alpha_hat <- fit$par[seq_len(q)]
  sigma_hat <- exp(fit$par[q + seq_len(q)])

  final <- .nominal_laplace_nll(
    Yq = Yq,
    alpha = alpha_hat,
    sigma = sigma_hat,
    Kinv = Kinv,
    logdetK = logdetK,
    inner_maxiter = inner_maxiter,
    inner_tol = inner_tol,
    step_factor = step_factor,
    ridge = ridge,
    return_mode = TRUE,
    verbose = FALSE
  )

  list(
    method = "laplace",
    alpha = alpha_hat,
    vc = sigma_hat,
    theta = log(sigma_hat),
    U = final$U,
    PI = final$PI,
    nll = final$nll,
    mode_trace = final$mode_trace,
    fit = fit,
    convergence = fit$convergence,
    message = fit$message,
    elapsed_time = elapsed_time,
    Y = Yfull,
    kk = kk
  )
}


