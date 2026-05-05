# ============================================================
# 1. SCORE METHOD
# ============================================================
scan_score_nominal <- function(ps, rr, par, kk, zz, x0, c, outdir = NULL) {
  m <- nrow(zz)
  
  vi <- solve(kk %x% diag(par) + rr)
  P  <- vi - vi %*% x0 %*% solve(t(x0) %*% vi %*% x0) %*% t(x0) %*% vi
  
  # start_time <- Sys.time()
  out <- NULL
  
  for (k in 1:m) {
    z <- as.matrix(zz[k, ] %x% diag(c - 1))
    zPz <- t(z) %*% P %*% z
    zPy <- t(z) %*% P %*% ps
    yPz <- t(zPy)
    
    zPz_inv <- solve(zPz)
    score <- yPz %*% zPz_inv %*% zPy
    effect <- zPz_inv %*% zPy
    stderr <- sqrt(diag(zPz_inv))
    p <- 1 - pchisq(score, c - 1)
    
    out <- rbind(out, c(k, as.numeric(effect), stderr, as.numeric(score), as.numeric(p)))
  }
  
  # end_time <- Sys.time()
  # print(end_time - start_time)
  
  out <- as.data.frame(out)
  
  colnames(out) <- c(
    "SNP",
    paste0("Effect", 1:(c - 1)),
    paste0("StdErr", 1:(c - 1)),
    "Score",
    "p"
  )
  
  if (!is.null(outdir)) {
    write.csv(out, file = file.path(outdir, "Nominal-Score-Test.csv"), row.names = FALSE)
  }
  
  return(out)
}

# ============================================================
# 2. P3D METHOD
#    par is fixed at null-model estimate
# ============================================================
scan_p3d_nominal <- function(y, kk, zz, x0, par, c, n, outdir = NULL) {
  
  iblup_k <- function(par, y, kk, x, z_marker) {
    a <- matrix(0, c + 1, 1)
    mu <- matrix(0, 1, c)
    a[1] <- -1e5
    a[c + 1] <- 1e5
    a[2:c] <- 0
    
    indx <- c:(2 * (c - 1))
    maxiter <- 100
    minerr <- 1e-8
    err <- 1e8
    iter <- 0
    
    b0 <- matrix(0, 2 * (c - 1), 1)
    g <- matrix(0, n * (c - 1), 1)
    gamma <- matrix(0, c - 1, 1)
    poly <- matrix(g, n, c - 1, byrow = TRUE)
    out <- NULL
    
    while ((iter < maxiter) & (err > minerr)) {
      rr <- NULL
      pe <- NULL
      lp <- NULL
      pp <- list()
      
      for (j in 1:n) {
        poly_j <- as.matrix(poly[j, ])
        z_j <- as.matrix(z_marker[j] %x% diag(c - 1)) %*% gamma
        
        tau <- 1
        for (k in 1:(c - 1)) {
          tau <- tau + exp(a[k + 1] + z_j[k] + poly_j[k])
        }
        
        for (k in 1:(c - 1)) {
          mu[k] <- exp(a[k + 1] + z_j[k] + poly_j[k]) / tau
        }
        mu[c] <- 1 / tau
        
        d <- matrix(0, c, c - 1)
        for (k in 1:(c - 1)) {
          for (i in 1:(c - 1)) {
            d[k, i] <- -exp(a[k + 1] + a[i + 1] + z_j[k] + z_j[i] + poly_j[k] + poly_j[i])
            if (k == i) {
              d[k, i] <- d[k, i] + tau * exp(a[k + 1] + z_j[k] + poly_j[k])
            }
          }
        }
        
        for (i in 1:(c - 1)) {
          d[c, i] <- -exp(a[i + 1] + z_j[i] + poly_j[i])
        }
        
        d <- d / tau^2
        w <- diag(1 / c(mu))
        dwd <- solve(t(d) %*% w %*% d + diag(c - 1) * 1e-5)
        
        lp <- rbind(lp, (diag(c - 1) %*% as.matrix(a[2:c]) + z_j + poly_j))
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
      gamma <- b[indx]
      vb <- xxi[indx, indx]
      stderr <- sqrt(diag(vb))
      wald <- t(gamma) %*% solve(vb) %*% gamma
      p <- 1 - pchisq(wald, c - 1)
      
      iter <- iter + 1
      err <- sum((b - b0)^2)
      b0 <- b
      
      out <- rbind(
        out,
        c(
          iter,
          err,
          # as.numeric(b[1:(c - 1)]),     # intercepts
          as.numeric(gamma),            # marker effects
          stderr,
          as.numeric(wald),
          as.numeric(p)
        )
      )
    }
    
    out[nrow(out), ]
  }
  
  m <- nrow(zz)
  # start_time <- Sys.time()
  res <- NULL
  
  for (k in 1:m) {
    z_marker <- as.matrix(zz[k, ])
    x <- cbind(x0, z_marker %x% diag(c - 1))
    one <- iblup_k(par = par, y = y, kk = kk, x = x, z_marker = z_marker)
    res <- rbind(res, c(k, one))
  }
  
  # end_time <- Sys.time()
  # print(end_time - start_time)
  
  res <- as.data.frame(res)
  colnames(res) <- c(
    "SNP",
    "iter",
    "err",
    # paste0("Intercept", 1:(c - 1)),
    paste0("Effect", 1:(c - 1)),
    paste0("StdErr", 1:(c - 1)),
    "Wald",
    "p"
  )
  
  if (!is.null(outdir)) {
    write.csv(res, file = file.path(outdir, "MyScan-P3D.csv"), row.names = FALSE)
  }
  
  return(res)
}

# ============================================================
# 3. PSR METHOD
#    variance components re-estimated at each marker
# ============================================================
scan_psr_nominal <- function(ps, rr, kk, zz, x0, c, n, theta0 = NULL, outdir = NULL) {
  if (is.null(theta0)) theta0 <- rep(0, c - 1)
  
  fixed_fun <- function(par, x) {
    v <- kk %x% diag(par) + rr
    vi <- solve(v + diag(n * (c - 1)) * 1e-5)
    xx <- t(x) %*% vi %*% x
    xy <- t(x) %*% vi %*% ps
    xxi <- solve(xx)
    b0 <- xxi %*% xy
    
    indx <- c:(2 * (c - 1))
    b <- as.matrix(b0[indx])
    vb <- xxi[indx, indx]
    wald <- t(b) %*% solve(vb) %*% b
    p <- 1 - pchisq(wald, c - 1)
    
    c(as.numeric(b), sqrt(diag(vb)), as.numeric(wald), as.numeric(p))
  }
  
  obj_fun <- function(theta, x) {
    v <- kk %x% diag(exp(theta)) + rr
    d1 <- determinant(v + diag(n * (c - 1)) * 1e-5)[[1]]
    vi <- solve(v + diag(n * (c - 1)) * 1e-5)
    xx <- t(x) %*% vi %*% x
    xy <- t(x) %*% vi %*% ps
    d2 <- determinant(xx)[[1]]
    b <- solve(xx) %*% xy
    yy <- t(ps - x %*% b) %*% vi %*% (ps - x %*% b)
    like <- 0.5 * d1 + 0.5 * d2 + 0.5 * yy
    as.numeric(like)
  }
  
  m <- nrow(zz)
  # start_time <- Sys.time()
  out <- NULL
  
  for (k in 1:m) {
    z <- as.matrix(zz[k, ])
    x <- cbind(x0, z %x% diag(c - 1))
    
    parm <- optim(
      par = theta0,
      fn = obj_fun,
      x = x,
      method = "L-BFGS-B",
      lower = -1e8,
      upper =  1e8
    )
    
    theta <- parm$par
    theta0 <- theta
    par <- exp(theta)
    
    psr <- fixed_fun(par, x)
    out <- rbind(out, c(k, par, psr))
  }
  
  # end_time <- Sys.time()
  # print(end_time - start_time)
  
  out <- as.data.frame(out)
  colnames(out) <- c(
    "SNP",
    paste0("Var", 1:(c - 1)),
    paste0("Effect", 1:(c - 1)),
    paste0("StdErr", 1:(c - 1)),
    "Wald",
    "p"
  )
  
  if (!is.null(outdir)) {
    write.csv(out, file = file.path(outdir, "MyScan-PSR-No-Spectral-Decom.csv"), row.names = FALSE)
  }
  
  return(out)
}

# ============================================================
# 4. GLM METHOD
#    No GRM control
# ============================================================
# one SNP
nominal_glm_one <- function(z, y, X, c, n) {
  q <- c - 1
  a <- matrix(0, c + 1, 1)
  mu <- matrix(0, 1, c)
  b <- matrix(0, q, 1)
  a[1] <- -1e3
  a[c + 1] <- 1e3
  a[2:c] <- 0
  g0 <- matrix(0, c - 1 + q, 1)
  maxiter <- 100
  minerr <- 1e-8
  err <- 1e8
  iter <- 0
  out <- NULL
  
  while ((iter < maxiter) & (err > minerr)) {
    rr <- NULL
    ps <- NULL
    eta <- NULL
    pp <- list()
    
    for (i in 1:n) {
      xb <- as.matrix(z[i]) %x% diag(c - 1) %*% b
      tau <- 1
      for (k in 1:(c - 1)) {
        tau <- tau + exp(a[k + 1] + xb[k])
      }
      for (k in 1:(c - 1)) {
        mu[k] <- exp(a[k + 1] + xb[k]) / tau
      }
      mu[c] <- 1 / tau
      mu <- mu / sum(mu)
      
      d <- matrix(0, c, c - 1)
      for (k in 1:(c - 1)) {
        for (j in 1:(c - 1)) {
          d[k, j] <- -exp(a[k + 1] + xb[k] + a[j + 1] + xb[j])
          if (k == j) {
            d[k, j] <- d[k, j] + tau * exp(a[k + 1] + xb[k])
          }
        }
      }
      for (j in 1:(c - 1)) {
        d[c, j] <- -exp(a[j + 1] + xb[j])
      }
      d <- d / tau^2
      
      w <- diag(1 / base::c(mu))
      dwd <- solve(t(d) %*% w %*% d + diag(c - 1) * 1e-5)
      eta <- rbind(eta, (diag(c - 1) %*% as.matrix(a[2:c]) + xb))
      ps <- rbind(ps, dwd %*% t(d) %*% w %*% as.matrix(base::c(y[i, ] - mu)))
      pp[[i]] <- dwd
    }
    
    rr <- as.matrix(bdiag(pp))
    vi <- solve(rr)
    xxi <- solve(t(X) %*% vi %*% X)
    g <- xxi %*% t(X) %*% vi %*% as.matrix(ps + eta)
    g <- as.matrix(g)
    a[2:c] <- g[1:(c - 1)]
    b <- as.matrix(g[c:(c - 1 + q)])
    
    iter <- iter + 1
    err <- sum((g - g0)^2)
    out <- rbind(out, cbind(iter, err, t(g)))
    g0 <- g
  }
  
  effect <- as.matrix(b)
  vb <- xxi[c:(c - 1 + q), c:(c - 1 + q)]
  wald <- t(effect) %*% solve(vb) %*% effect
  p <- 1 - pchisq(wald, q)
  
  list(
    effect = effect,
    vb = vb,
    wald = wald,
    p = p
  )
}

# genome scan
scan_glm_nominal <- function(zz, y, x0) {
  c <- ncol(y)
  n <- nrow(y)
  m <- nrow(zz)
  
  # start_time <- Sys.time()
  www <- NULL
  
  for (k in 1:m) {
    z <- as.matrix(zz[k, ])
    X <- cbind(x0, z %x% diag(c - 1))
    pp <- nominal_glm_one(z = z, y = y, X = X, c = c, n = n)
    
    effect <- pp$effect
    vb <- pp$vb
    wald <- pp$wald
    p <- pp$p
    stderr <- sqrt(diag(vb))
    
    www <- rbind(www, c(k, base::c(effect), base::c(stderr), wald, p))
  }
  
  # end_time <- Sys.time()
  # print(end_time - start_time)
  
  www <- data.frame(www)
  colnames(www) <- c(
    "SNP",
    paste0("Effect", 1:(c - 1)),
    paste0("StdErr", 1:(c - 1)),
    "Wald",
    "p"
  )
  
  www
}

