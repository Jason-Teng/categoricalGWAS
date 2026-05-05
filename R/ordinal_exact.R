# ============================================================
# Ordinal exact genome scan
# Marker-by-marker doubly iterative algorithm.
# This method does NOT require a pre-fit null model.
# ============================================================

scan_exact_ordinal <- function(y, kk, zz, x0, c = NULL, n = NULL,
                               theta0 = 0,
                               maxiter = 100,
                               minerr = 1e-8,
                               lower = -1e5,
                               upper = 1e5,
                               ridge = 1e-5,
                               outdir = NULL) {
  if (is.null(c)) c <- ncol(y)
  if (is.null(n)) n <- nrow(y)
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }

  kk <- as.matrix(kk)
  zz <- as.matrix(zz)
  y <- as.matrix(y)
  x0 <- as.matrix(x0)

  m <- nrow(zz)
  J <- matrix(1, c - 1, 1)
  h <- diag(n) %x% J

  start_time <- Sys.time()
  out_all <- NULL
  error_snps <- integer(0)

  for (marker_id in 1:m) {
    z_raw <- as.numeric(zz[marker_id, ])
    z <- as.matrix(z_raw) %x% J
    x <- cbind(x0, z)

    iblup_k <- function(par) {
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
      iter_hist <- NULL

      while ((iter < maxiter) && (err > minerr)) {
        pe <- NULL
        lp <- NULL
        pp <- vector("list", n)

        for (j in 1:n) {
          # Preserve the legacy exact-method behavior: the original script
          # uses z[j] here, where z is already Kronecker-expanded.
          # This is intentionally not z_raw[j], so results match the draft.
          eta_j <- g0[j] + z[j] * gamma

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
          dwd <- solve(t(d) %*% w %*% d + diag(c - 1) * ridge)

          lp <- rbind(lp, as.matrix(a[2:c]) + g0[j] %x% J + z[j] %x% J %*% gamma)
          pe <- rbind(pe, dwd %*% t(d) %*% w %*% as.matrix(as.numeric(y[j, ] - mu)))
          pp[[j]] <- dwd
        }

        ps <- pe + lp
        rr <- as.matrix(Matrix::bdiag(pp))
        vi <- solve(h %*% kk %*% t(h) * par + rr + diag(n * (c - 1)) * ridge)
        xx <- t(x) %*% vi %*% x
        xy <- t(x) %*% vi %*% ps
        xxi <- solve(xx)
        b <- xxi %*% xy
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
        iter_hist <- rbind(iter_hist, c(iter, err, t(b), wald, p))
      }

      last <- iter_hist[nrow(iter_hist), ]
      list(ps = ps, rr = rr, lp = lp, pe = pe, out = last)
    }

    ordinal_fun_marker <- function(theta, ps, rr) {
      s2 <- exp(theta)
      G <- h %*% kk %*% t(h) * s2
      v <- G + rr
      d1 <- determinant(v)[[1]]
      vi <- solve(v + diag(n * (c - 1)) * ridge)
      xx <- t(x) %*% vi %*% x
      xy <- t(x) %*% vi %*% ps
      d2 <- determinant(xx)[[1]]
      b <- solve(xx) %*% xy
      yy <- t(ps - x %*% b) %*% vi %*% (ps - x %*% b)
      like <- 0.5 * d1 + 0.5 * d2 + 0.5 * yy
      as.numeric(like)
    }

    tryCatch({
      theta_old <- theta0
      par_old <- exp(theta_old)
      double_iter <- NULL
      err_outer <- 1e8
      iter_outer <- 0
      my_result <- NULL

      while ((iter_outer < maxiter) && (err_outer > minerr)) {
        my_result <- iblup_k(par_old)
        opt <- optim(
          par = theta_old,
          fn = ordinal_fun_marker,
          ps = my_result$ps,
          rr = my_result$rr,
          method = "L-BFGS-B",
          lower = lower,
          upper = upper
        )

        theta_new <- opt$par
        par_new <- exp(theta_new)
        err_outer <- sum((par_new - par_old)^2)
        iter_outer <- iter_outer + 1
        theta_old <- theta_new
        par_old <- par_new
        double_iter <- rbind(double_iter, c(iter_outer, err_outer, par_new))
      }

      out_all <- rbind(out_all, c(marker_id, par_old, my_result$out))
    }, error = function(e) {
      error_snps <<- c(error_snps, marker_id)
      out_all <<- rbind(out_all, c(marker_id, rep(NA_real_, 7)))
    })
  }

  out <- as.data.frame(out_all)
  colnames(out) <- c("SNP", "V", "Iter", "Error", paste0("a", 1:(c - 1)), "Effect", "Wald", "p")

  elapsed_time <- Sys.time() - start_time
  attr(out, "elapsed_time") <- elapsed_time
  attr(out, "error_snps") <- error_snps

  if (!is.null(outdir)) {
    write.csv(out, file = file.path(outdir, "Ordinal-Exact-Out.csv"), row.names = FALSE)
  }

  out
}
