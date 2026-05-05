library(Matrix)


# ============================================================
# Inner loop: exact refactor of original iblup()
# ============================================================
iblup_nominal <- function(par, y, x, kk, maxiter = 100, minerr = 1e-8) {
  c <- ncol(y)
  n <- nrow(y)
  
  a <- matrix(0, c + 1, 1)
  mu <- matrix(0, 1, c)
  a[1] <- -1e5
  a[c + 1] <- 1e5
  a[2:c] <- 0
  
  err <- 1e8
  iter <- 0
  b0 <- matrix(0, c - 1, 1)
  g <- matrix(0, n * (c - 1), 1)
  poly <- matrix(g, n, c - 1, byrow = TRUE)
  out <- NULL
  
  while ((iter < maxiter) & (err > minerr)) {
    pe <- NULL
    lp <- NULL
    pp <- list()
    
    for (j in 1:n) {
      poly.j <- as.matrix(poly[j, ])
      tau <- 1
      
      for (k in 1:(c - 1)) {
        tau <- tau + exp(a[k + 1] + poly.j[k])
      }
      
      for (k in 1:(c - 1)) {
        mu[k] <- exp(a[k + 1] + poly.j[k]) / tau
      }
      mu[c] <- 1 / tau
      
      d <- matrix(0, c, c - 1)
      for (k in 1:(c - 1)) {
        for (i in 1:(c - 1)) {
          d[k, i] <- -exp(a[k + 1] + a[i + 1] + poly.j[k] + poly.j[i])
          if (k == i) {
            d[k, i] <- d[k, i] + tau * exp(a[k + 1] + poly.j[k])
          }
        }
      }
      
      for (i in 1:(c - 1)) {
        d[c, i] <- -exp(a[i + 1] + poly.j[i])
      }
      
      d <- d / tau^2
      w <- diag(1 / c(mu))
      dwd <- solve(t(d) %*% w %*% d + diag(c - 1) * 1e-5)
      
      lp <- rbind(lp, diag(c - 1) %*% as.matrix(a[2:c]) + poly.j)
      pe <- rbind(pe, dwd %*% t(d) %*% w %*% as.matrix(c(y[j, ] - mu)))
      pp[[j]] <- dwd
    }
    
    ps <- pe + lp
    rr <- as.matrix(bdiag(pp))
    vi <- solve(kk %x% diag(par) + rr)
    xx <- t(x) %*% vi %*% x
    xy <- t(x) %*% vi %*% ps
    xxi <- solve(xx)
    b <- xxi %*% xy
    poly <- (kk %x% diag(par)) %*% vi %*% (ps - x %*% b)
    poly <- matrix(poly, n, c - 1, byrow = TRUE)
    a[2:c] <- b[1:(c - 1)]
    
    iter <- iter + 1
    err <- sum((b - b0)^2)
    b0 <- b
    out <- rbind(out, cbind(iter, err, t(b)))
  }
  
  list(
    # par = par, adjsut this in wrapper function
    ps = ps,
    rr = rr,
    lp = lp,
    pe = pe,
    out = out,
    beta = b
  )
}

# ============================================================
# Outer objective: exact refactor of original fun()
# ============================================================

# using Cholesky
nominal_fun <- function(theta, ps, rr, x, kk) {
  
  s2 <- diag(exp(theta))
  G <- kk %x% s2
  V <- G + rr
  
  # --- Cholesky (more stable) ---
  cholV <- chol(V)
  
  # log |V|
  d1 <- 2 * sum(log(diag(cholV)))
  
  # V^{-1}
  Vi <- chol2inv(cholV)
  
  # X' V^{-1} X
  XVX <- t(x) %*% Vi %*% x
  
  # log |X' V^{-1} X|
  cholXVX <- chol(XVX)
  d2 <- 2 * sum(log(diag(cholXVX)))
  
  # beta
  XVy <- t(x) %*% Vi %*% ps
  b <- solve(XVX, XVy)
  
  # quadratic form
  res <- ps - x %*% b
  yy <- t(res) %*% Vi %*% res
  
  like <- 0.5 * (d1 + d2 + yy)
  
  return(as.numeric(like))
}

# ============================================================
# Main fitting: exact refactor of original outer while loop
# ============================================================
fit_nominal_outer <- function(
    y,
    x,
    kk,
    ids,
    iblup_func,
    fun_func,
    theta0 = NULL,
    maxiter = 100,
    minerr = 1e-8,
    lower = -1e5,
    upper = 1e5,
    outdir = "Output-Nominal") {
  
  c <- ncol(y)
  n <- nrow(y)
  
  if (is.null(theta0)) {
    theta0 <- rep(0, c - 1)
  }
  
  par0 <- exp(theta0)
  
  double.iter <- NULL
  err <- 1e8
  iter <- 0
  
  # start_time <- Sys.time()
  
  while ((iter < maxiter) & (err > minerr)) {
    
    # ===== inner loop =====
    myResult <- iblup_func(par = par0, y = y, x = x, kk = kk)
    
    ps <- myResult[[1]]
    rr <- myResult[[2]]
    
    # ===== variance component optimization =====
    parm <- optim(
      par = theta0,
      fn = fun_func,
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
    
    double.iter <- rbind(double.iter, c(iter, err, par))
  }
  
  # ===== final reconstruction =====
  myResult <- iblup_func(par = par0, y = y, x = x, kk = kk)
  ps <- myResult[[1]]
  rr <- myResult[[2]]
  
  # end_time <- Sys.time()
  # elapsed_time <- end_time - start_time
  
  # print(elapsed_time)
  
  double.iter <- as.data.frame(double.iter)
  
  # naming (keep exactly same structure)
  colnames(double.iter) <- c(
    "Iter", "Error",
    paste0("Var", 1:(c - 1))
  )
  
  # ===== compact rr =====
  newID <- as.matrix(ids) %x% matrix(1, c - 1, 1)
  rrNoDiag <- NULL
  
  for (j in 1:n) {
    indx <- which(newID == j)
    rrNoDiag <- rbind(rrNoDiag, rr[indx, indx])
  }
  
  # ===== export =====
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  write.csv(rrNoDiag, file = file.path(outdir, "rrNoDiag.csv"), row.names = FALSE)
  write.csv(ps,       file = file.path(outdir, "ps.csv"),       row.names = FALSE)
  write.csv(rr,       file = file.path(outdir, "rr.csv"),       row.names = FALSE)
  
  return(list(
    vc = par0,
    theta = theta0,
    trace = double.iter,
    ps = ps,
    rr = rr,
    rrNoDiag = rrNoDiag,
    elapsed_time = elapsed_time,
    optim = parm
  ))
}

