# ============================================================
# Ordinal null model utilities
# ============================================================
# This file extracts the polygenic null model from the original
# ordinal scripts into reusable functions.
#
# Main user-facing function:
#   fit_ordinal_null()
#
# Core internal functions:
#   ordinal_make_Y()
#   iblup_ordinal()
#   ordinal_fun()
#   fit_ordinal_outer()
# ============================================================

library(Matrix)

# Convert an ordinal phenotype vector into the category indicator matrix.
# If y is already a matrix, it is returned as a numeric matrix.
ordinal_make_Y <- function(y) {
  if (is.matrix(y) || is.data.frame(y)) {
    return(as.matrix(y))
  }

  y_fac <- as.factor(y)
  model.matrix(~ y_fac - 1)
}

# Convert block diagonal rr into the legacy compact format used by PSRSD.
# For q = c - 1, each individual contributes a q by q block, so the
# compact object has n*q rows and q columns.
ordinal_compact_rr <- function(rr, n, q, ids = NULL) {
  if (is.null(ids)) ids <- seq_len(n)

  rrNoDiag <- NULL
  row_ids <- as.matrix(ids) %x% matrix(1, q, 1)

  for (j in seq_len(n)) {
    idx <- ((j - 1) * q + 1):(j * q)
    rrNoDiag <- rbind(rrNoDiag, rr[idx, idx])
  }

  colnames(rrNoDiag) <- paste0("Cov", seq_len(q))
  rownames(rrNoDiag) <- as.vector(row_ids)
  rrNoDiag
}

# ============================================================
# Inner loop: ordinal IBLUP / pseudo-response generation
# ============================================================
iblup_ordinal <- function(par, y, x, kk, maxiter = 100, minerr = 1e-8) {

  y <- as.matrix(y)
  kk <- as.matrix(kk)

  c <- ncol(y)
  n <- nrow(y)
  q <- c - 1

  if (length(par) != 1) {
    stop("For the current ordinal model, 'par' should be a single polygenic variance component.")
  }

  J <- matrix(1, q, 1)
  h <- diag(n) %x% J

  a <- matrix(0, c + 1, 1)
  mu <- matrix(0, 1, c)
  a[1] <- -1e5
  a[c + 1] <- 1e5
  a[2:c] <- 0

  err <- 1e8
  iter <- 0
  b0 <- a[2:c]
  g0 <- matrix(0, n, 1)
  out <- NULL

  while ((iter < maxiter) & (err > minerr)) {

    pe <- NULL
    lp <- NULL
    pp <- vector("list", n)

    for (j in seq_len(n)) {

      for (k in seq_len(c)) {
        mu[k] <- pnorm(a[k + 1] + g0[j]) - pnorm(a[k] + g0[j])
        if (mu[k] < 1e-4) mu[k] <- 1e-4
        if (mu[k] > (1 - 1e-4)) mu[k] <- 1 - 1e-4
      }

      mu <- mu / sum(mu)

      d <- matrix(0, c, q)

      if (c > 2) {
        for (k in 2:q) {
          d[k, k - 1] <- -dnorm(a[k] + g0[j])
          d[k, k] <- dnorm(a[k + 1] + g0[j])
        }
      }

      d[1, 1] <- dnorm(a[2] + g0[j])
      d[c, q] <- -dnorm(a[c] + g0[j])

      w <- diag(1 / c(mu))
      dwd <- solve(t(d) %*% w %*% d + diag(q) * 1e-5)

      lp_j <- as.matrix(a[2:c]) + g0[j] %x% J
      pe_j <- dwd %*% t(d) %*% w %*% c(y[j, ] - mu)

      lp <- rbind(lp, lp_j)
      pe <- rbind(pe, pe_j)
      pp[[j]] <- dwd
    }

    ps <- pe + lp
    rr <- as.matrix(Matrix::bdiag(pp))

    V <- h %*% kk %*% t(h) * par + rr
    Vi <- solve(V)

    XVX <- t(x) %*% Vi %*% x
    XVy <- t(x) %*% Vi %*% ps
    XVXi <- solve(XVX)
    b <- XVXi %*% XVy

    g <- (kk %*% t(h) * par) %*% Vi %*% (ps - x %*% b)

    a[2:c] <- b[1:q]
    iter <- iter + 1
    err <- sum((b - b0)^2)
    b0 <- b
    g0 <- g

    out <- rbind(out, cbind(iter, err, t(b)))
  }

  list(
    ps = ps,
    rr = rr,
    lp = lp,
    pe = pe,
    out = out,
    beta = b,
    g = g0
  )
}

# ============================================================
# Outer objective: variance component estimation
# ============================================================
ordinal_fun <- function(theta, ps, rr, x, kk) {

  ps <- as.matrix(ps)
  rr <- as.matrix(rr)
  kk <- as.matrix(kk)

  n <- nrow(kk)
  q <- ncol(x)

  J <- matrix(1, q, 1)
  h <- diag(n) %x% J

  s2 <- exp(theta)
  V <- h %*% kk %*% t(h) * s2 + rr

  cholV <- chol(V)
  logdetV <- 2 * sum(log(diag(cholV)))
  Vi <- chol2inv(cholV)

  XVX <- t(x) %*% Vi %*% x
  cholXVX <- chol(XVX)
  logdetXVX <- 2 * sum(log(diag(cholXVX)))

  XVy <- t(x) %*% Vi %*% ps
  b <- solve(XVX, XVy)
  res <- ps - x %*% b
  quad <- t(res) %*% Vi %*% res

  as.numeric(0.5 * logdetV + 0.5 * logdetXVX + 0.5 * quad)
}

# ============================================================
# Main fitting loop for the ordinal polygenic null model
# ============================================================
fit_ordinal_outer <- function(y,
                              x,
                              kk,
                              ids = NULL,
                              theta0 = 0,
                              maxiter = 100,
                              minerr = 1e-8,
                              lower = -1e5,
                              upper = 1e5,
                              outdir = NULL) {

  y <- as.matrix(y)
  kk <- as.matrix(kk)

  c <- ncol(y)
  n <- nrow(y)
  q <- c - 1

  if (is.null(ids)) ids <- seq_len(n)

  par0 <- exp(theta0)
  trace <- NULL
  err <- 1e8
  iter <- 0
  parm <- NULL

  start_time <- Sys.time()

  while ((iter < maxiter) & (err > minerr)) {

    pseudo <- iblup_ordinal(
      par = par0,
      y = y,
      x = x,
      kk = kk,
      maxiter = maxiter,
      minerr = minerr
    )

    ps <- pseudo$ps
    rr <- pseudo$rr

    parm <- optim(
      par = theta0,
      fn = ordinal_fun,
      ps = ps,
      rr = rr,
      x = x,
      kk = kk,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper
    )

    theta <- parm$par
    par <- exp(theta)

    err <- sum((par - par0)^2)
    iter <- iter + 1

    theta0 <- theta
    par0 <- par

    trace <- rbind(trace, c(iter, err, par))
  }

  # ===== final reconstruction =====
  pseudo <- iblup_ordinal(
    par = par0,
    y = y,
    x = x,
    kk = kk,
    maxiter = maxiter,
    minerr = minerr
  )

  elapsed_time <- Sys.time() - start_time

  trace <- as.data.frame(trace)
  colnames(trace) <- c("Iter", "Error", "Var")

  ps <- pseudo$ps
  rr <- pseudo$rr

  rrNoDiag <- ordinal_compact_rr(rr = rr, n = n, q = q, ids = ids)

  if (!is.null(outdir)) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    write.csv(trace, file = file.path(outdir, "ordinal_null_trace.csv"), row.names = FALSE)
    write.csv(ps, file = file.path(outdir, "ordinal_ps.csv"), row.names = FALSE)
    write.csv(rr, file = file.path(outdir, "ordinal_rr.csv"), row.names = FALSE)
    write.csv(rrNoDiag, file = file.path(outdir, "ordinal_rrNoDiag.csv"), row.names = TRUE)
  }

  list(
    method = "pseudo",
    par = par0,
    vc = par0,
    theta = theta0,
    ps = ps,
    rr = rr,
    rrNoDiag = rrNoDiag,
    beta = pseudo$beta,
    g = pseudo$g,
    trace = trace,
    elapsed_time = elapsed_time,
    optim = parm
  )
}

# ============================================================
# User-facing null model wrapper
# ============================================================
fit_ordinal_null <- function(y,
                             x0 = NULL,
                             kk,
                             ids = NULL,
                             null_method = c("pseudo", NULL),
                             vc = NULL,
                             maxiter = 100,
                             minerr = 1e-8,
                             outdir = NULL) {

  null_method <- match.arg(null_method)

  Yfull <- ordinal_make_Y(y)
  n <- nrow(Yfull)
  c <- ncol(Yfull)
  q <- c - 1

  if (is.null(x0)) {
    x0 <- matrix(1, n, 1) %x% diag(q)
  }

  if (!is.null(vc)) {
    message(">>> Start ordinal null linearization with provided VC")
    null_start <- Sys.time()

    pseudo <- iblup_ordinal(
      par = vc,
      y = Yfull,
      x = x0,
      kk = kk,
      maxiter = maxiter,
      minerr = minerr
    )

    message(paste("Null linearization time:", format(Sys.time() - null_start)))

    rrNoDiag <- ordinal_compact_rr(rr = pseudo$rr, n = n, q = q, ids = ids)

    return(list(
      method = "provided_vc",
      par = vc,
      vc = vc,
      theta = log(vc),
      ps = pseudo$ps,
      rr = pseudo$rr,
      rrNoDiag = rrNoDiag,
      beta = pseudo$beta,
      g = pseudo$g,
      trace = pseudo$out,
      elapsed_time = Sys.time() - null_start
    ))
  }

  if (null_method == "pseudo") {
    message(">>> Start fitting ordinal null model (pseudo)")

    fit_ordinal <- fit_ordinal_outer(
      y = Yfull,
      x = x0,
      kk = kk,
      ids = ids,
      maxiter = maxiter,
      minerr = minerr,
      outdir = outdir
    )
    message(paste("Ordinal null model time:", format(fit_ordinal$elapsed_time)))

    return(fit_ordinal)
  }
}
