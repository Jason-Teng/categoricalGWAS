# ============================================================
# Ordinal genome scan methods
# Requires: Matrix package for bdiag()
# Assumes ordinal_make_Y() / fit_ordinal_null() are defined in ordinal_null.R
# ============================================================

# ============================================================
# 1. SCORE METHOD
#    Uses pseudo-response and residual covariance from the null model.
#    Ordinal model has one SNP effect shared across thresholds, so df = 1.
# ============================================================
scan_score_ordinal <- function(ps, rr, par, kk, zz, x0, c, n = NULL, outdir = NULL) {
  if (is.null(n)) n <- nrow(kk)
  m <- nrow(zz)
  J <- matrix(1, c - 1, 1)
  h <- diag(n) %x% J

  vi <- solve(h %*% kk %*% t(h) * par + rr)
  P <- vi - vi %*% x0 %*% solve(t(x0) %*% vi %*% x0) %*% t(x0) %*% vi

  out <- NULL
  for (k in 1:m) {
    z <- as.matrix(zz[k, ]) %x% J
    zPz <- t(z) %*% P %*% z
    zPy <- t(z) %*% P %*% ps

    zPz_inv <- solve(zPz)
    score <- t(zPy) %*% zPz_inv %*% zPy
    effect <- zPz_inv %*% zPy
    stderr <- sqrt(diag(zPz_inv))
    p <- 1 - pchisq(score, 1)

    out <- rbind(out, c(k, as.numeric(effect), as.numeric(stderr), as.numeric(score), as.numeric(p)))
  }

  out <- as.data.frame(out)
  colnames(out) <- c("SNP", "Effect", "StdErr", "Score", "p")

  if (!is.null(outdir)) {
    write.csv(out, file = file.path(outdir, "Ordinal-Score-Test.csv"), row.names = FALSE)
  }

  out
}

# ============================================================
# 2. P3D METHOD
#    Fixes the variance component at the null-model estimate.
#    Refits the ordinal pseudo-response model marker by marker.
# ============================================================
scan_p3d_ordinal <- function(y, kk, zz, x0, par, c, n = NULL,
                             maxiter = 100, minerr = 1e-8, outdir = NULL) {
  if (is.null(n)) n <- nrow(y)
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }

  J <- matrix(1, c - 1, 1)
  h <- diag(n) %x% J

  iblup_k <- function(z_marker) {
    a <- matrix(0, c + 1, 1)
    mu <- matrix(0, 1, c)
    a[1] <- -1e5
    a[c + 1] <- 1e5
    a[2:c] <- 0

    marker_index <- c
    err <- 1e8
    iter <- 0
    b0 <- matrix(0, c, 1)
    g0 <- matrix(0, n, 1)
    gamma <- matrix(0, 1, 1)
    out <- NULL
    x <- cbind(x0, z_marker %x% J)

    while ((iter < maxiter) && (err > minerr)) {
      pe <- NULL
      lp <- NULL
      pp <- vector("list", n)

      for (j in 1:n) {
        eta_j <- g0[j] + z_marker[j] * gamma

        for (kkk in 1:c) {
          mu[kkk] <- pnorm(a[kkk + 1] + eta_j) - pnorm(a[kkk] + eta_j)
          if (mu[kkk] < 1e-4) mu[kkk] <- 1e-4
          if (mu[kkk] > 1 - 1e-4) mu[kkk] <- 1 - 1e-4
        }
        mu <- mu / sum(mu)

        d <- matrix(0, c, c - 1)
        if (c > 2) {
          for (kkk in 2:(c - 1)) {
            d[kkk, kkk - 1] <- -dnorm(a[kkk] + eta_j)
            d[kkk, kkk] <- dnorm(a[kkk + 1] + eta_j)
          }
        }
        d[1, 1] <- dnorm(a[2] + eta_j)
        d[c, c - 1] <- -dnorm(a[c] + eta_j)

        w <- diag(1 / as.numeric(mu))
        dwd <- solve(t(d) %*% w %*% d + diag(c - 1) * 1e-5)

        lp <- rbind(lp, as.matrix(a[2:c]) + g0[j] %x% J + (z_marker[j] %x% J) %*% gamma)
        pe <- rbind(pe, dwd %*% t(d) %*% w %*% as.matrix(as.numeric(y[j, ] - mu)))
        pp[[j]] <- dwd
      }

      ps <- pe + lp
      rr <- as.matrix(Matrix::bdiag(pp))
      vi <- solve(h %*% kk %*% t(h) * par + rr)
      xxi <- solve(t(x) %*% vi %*% x)
      b <- xxi %*% t(x) %*% vi %*% ps
      g <- par * kk %*% t(h) %*% vi %*% (ps - x %*% b)

      a[2:c] <- b[1:(c - 1)]
      gamma <- b[marker_index]
      vb <- xxi[marker_index, marker_index]
      wald <- as.numeric(t(gamma) %*% solve(matrix(vb, 1, 1)) %*% gamma)
      p <- 1 - pchisq(wald, 1)

      iter <- iter + 1
      err <- sum((b - b0)^2)
      b0 <- b
      g0 <- g

      out <- rbind(out, c(iter, err, as.numeric(b[1:(c - 1)]), as.numeric(gamma), sqrt(vb), wald, p))
    }

    out[nrow(out), ]
  }

  res <- NULL
  for (k in 1:nrow(zz)) {
    z_marker <- as.matrix(zz[k, ])
    one <- iblup_k(z_marker)
    res <- rbind(res, c(k, one))
  }

  res <- as.data.frame(res)
  colnames(res) <- c("SNP", "iter", "err", paste0("Intercept", 1:(c - 1)), "Effect", "StdErr", "Wald", "p")

  if (!is.null(outdir)) {
    write.csv(res, file = file.path(outdir, "Ordinal-P3D.csv"), row.names = FALSE)
  }

  res
}

# ============================================================
# 3. PSR METHOD WITHOUT SPECTRAL DECOMPOSITION
#    Re-estimates the variance component marker by marker using ps and rr.
# ============================================================
scan_psr_ordinal <- function(ps, rr, kk, zz, x0, c, n = NULL,
                             theta0 = 0, outdir = NULL) {
  if (is.null(n)) n <- nrow(kk)
  J <- matrix(1, c - 1, 1)
  h <- diag(n) %x% J
  m <- nrow(zz)

  out <- NULL
  for (k in 1:m) {
    z <- as.matrix(zz[k, ])
    x <- cbind(x0, z %x% J)
    marker_index <- c

    obj_fun <- function(theta) {
      v <- h %*% kk %*% t(h) * exp(theta) + rr
      v_jit <- v + diag(n * (c - 1)) * 1e-5
      vi <- solve(v_jit)
      xx <- t(x) %*% vi %*% x
      xy <- t(x) %*% vi %*% ps
      b <- solve(xx) %*% xy
      yy <- t(ps - x %*% b) %*% vi %*% (ps - x %*% b)
      0.5 * determinant(v_jit)[[1]] + 0.5 * determinant(xx)[[1]] + 0.5 * yy
    }

    fixed_fun <- function(par) {
      v <- h %*% kk %*% t(h) * par + rr
      vi <- solve(v + diag(n * (c - 1)) * 1e-5)
      xxi <- solve(t(x) %*% vi %*% x)
      b0 <- xxi %*% t(x) %*% vi %*% ps
      effect <- as.matrix(b0[marker_index])
      vb <- xxi[marker_index, marker_index]
      stderr <- sqrt(vb)
      wald <- as.numeric(effect^2 / vb)
      p <- 1 - pchisq(wald, 1)
      c(as.numeric(effect), as.numeric(stderr), wald, p)
    }

    parm <- optim(par = theta0, fn = obj_fun, method = "L-BFGS-B", lower = -1e8, upper = 1e8)
    theta0 <- parm$par
    par <- exp(theta0)
    tests <- fixed_fun(par)
    out <- rbind(out, c(k, par, tests))
  }

  out <- as.data.frame(out)
  colnames(out) <- c("SNP", "Var", "Effect", "StdErr", "Wald", "p")

  if (!is.null(outdir)) {
    write.csv(out, file = file.path(outdir, "Ordinal-PSR-No-SD.csv"), row.names = FALSE)
  }

  out
}

# ============================================================
# 4. PSRSD METHOD
#    PSR with spectral decomposition after whitening by rr.
# ============================================================
scan_psrsd_ordinal <- function(ps, rr, kk, zz, x0, c, n = NULL,
                               theta0 = 0, outdir = NULL) {
  if (is.null(n)) n <- nrow(kk)
  J <- matrix(1, c - 1, 1)

  r <- chol(solve(rr))
  h0 <- diag(n) %x% J
  h <- t(r) %*% h0
  eig <- eigen(h %*% kk %*% t(h))
  uu <- eig$vectors
  dd <- eig$values

  y_star <- t(uu) %*% t(r) %*% ps
  x0_star <- t(uu) %*% t(r) %*% x0

  reml_fun <- function(theta, x, y) {
    par <- exp(theta)
    wt <- 1 / (dd * par + 1)
    xx <- t(x) %*% (x * wt)
    xy <- t(x) %*% (y * wt)
    yy <- t(y) %*% (y * wt)
    ss <- yy - t(xy) %*% solve(xx) %*% xy
    0.5 * sum(log(dd * par + 1)) + 0.5 * determinant(xx)[[1]] + 0.5 * ss
  }

  fixed_fun <- function(par, x, y) {
    wt <- 1 / (dd * par + 1)
    xxi <- solve(t(x) %*% (x * wt))
    b <- xxi %*% t(x) %*% (y * wt)
    effect <- b[c]
    stderr <- sqrt(xxi[c, c])
    wald <- as.numeric((effect / stderr)^2)
    p <- 1 - pchisq(wald, 1)
    c(as.numeric(effect), as.numeric(stderr), wald, p)
  }

  out <- NULL
  for (k in 1:nrow(zz)) {
    z_star <- t(uu) %*% t(r) %*% as.matrix(zz[k, ] %x% J)
    x <- cbind(x0_star, z_star)
    parm <- optim(par = theta0, fn = reml_fun, x = x, y = y_star,
                  method = "L-BFGS-B", lower = -1e5, upper = 1e5)
    theta0 <- parm$par
    par <- exp(theta0)
    tests <- fixed_fun(par, x, y_star)
    out <- rbind(out, c(k, par, tests))
  }

  out <- as.data.frame(out)
  colnames(out) <- c("SNP", "Var", "Effect", "StdErr", "Wald", "p")

  if (!is.null(outdir)) {
    write.csv(out, file = file.path(outdir, "Ordinal-PSRSD.csv"), row.names = FALSE)
  }

  out
}

# ============================================================
# 5. GLM METHOD
#    Fixed-effect ordinal model without kinship control.
# ============================================================
ordinal_glm_one <- function(z, y, x, c, n,
                            maxiter = 100, minerr = 1e-8) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }

  a <- matrix(0, c + 1, 1)
  mu <- matrix(0, 1, c)
  a[1] <- -1e5
  a[c + 1] <- 1e5
  a[2:c] <- 0

  b0 <- matrix(0, c, 1)
  g0 <- 0
  err <- 1e8
  iter <- 0
  out <- NULL

  while ((iter < maxiter) && (err > minerr)) {
    pe <- NULL
    lp <- NULL
    pp <- vector("list", n)

    for (j in 1:n) {
      eta_j <- z[j] * g0
      for (kkk in 1:c) {
        mu[kkk] <- pnorm(a[kkk + 1] + eta_j) - pnorm(a[kkk] + eta_j)
        if (mu[kkk] < 1e-4) mu[kkk] <- 1e-4
        if (mu[kkk] > 1 - 1e-4) mu[kkk] <- 1 - 1e-4
      }
      mu <- mu / sum(mu)

      d <- matrix(0, c, c - 1)
      if (c > 2) {
        for (kkk in 2:(c - 1)) {
          d[kkk, kkk - 1] <- -dnorm(a[kkk] + eta_j)
          d[kkk, kkk] <- dnorm(a[kkk + 1] + eta_j)
        }
      }
      d[1, 1] <- dnorm(a[2] + eta_j)
      d[c, c - 1] <- -dnorm(a[c] + eta_j)

      w <- diag(1 / as.numeric(mu))
      dwd <- solve(t(d) %*% w %*% d + diag(c - 1) * 1e-5)
      lp <- rbind(lp, as.matrix(a[2:c]) + (z[j] %x% matrix(1, c - 1, 1)) %*% g0)
      pe <- rbind(pe, dwd %*% t(d) %*% w %*% as.matrix(as.numeric(y[j, ] - mu)))
      pp[[j]] <- dwd
    }

    ps <- pe + lp
    rr <- as.matrix(Matrix::bdiag(pp))
    vi <- solve(rr)
    xxi <- solve(t(x) %*% vi %*% x)
    b <- xxi %*% t(x) %*% vi %*% ps

    a[2:c] <- b[1:(c - 1)]
    g <- b[c]
    iter <- iter + 1
    err <- sum((b - b0)^2 + (g - g0)^2)
    b0 <- b
    g0 <- g
    out <- rbind(out, c(iter, err, as.numeric(b)))
  }

  effect <- as.numeric(g0)
  stderr <- sqrt(xxi[c, c])
  wald <- effect^2 / stderr^2
  p <- 1 - pchisq(wald, 1)

  c(iter, err, as.numeric(b[1:(c - 1)]), effect, stderr, wald, p)
}

scan_glm_ordinal <- function(zz, y, x0, c = NULL, n = NULL, outdir = NULL) {
  if (is.null(c)) c <- ncol(y)
  if (is.null(n)) n <- nrow(y)

  out <- NULL
  J <- matrix(1, c - 1, 1)
  for (k in 1:nrow(zz)) {
    z <- as.matrix(zz[k, ])
    x <- cbind(x0, z %x% J)
    one <- ordinal_glm_one(z = z, y = y, x = x, c = c, n = n)
    out <- rbind(out, c(k, one))
  }

  out <- as.data.frame(out)
  colnames(out) <- c("SNP", "Iter", "Err", paste0("Intercept", 1:(c - 1)), "Effect", "StdErr", "Wald", "p")

  if (!is.null(outdir)) {
    write.csv(out, file = file.path(outdir, "Ordinal-GLM.csv"), row.names = FALSE)
  }

  out
}
