iblup.k <- function(par, y, z, x, kk) {   # Take the variances and go to the IBLUP loop to update the pseudo response
  c <- ncol(y)
  n <- nrow(y)
  
  a <- matrix(0, c + 1, 1)
  mu <- matrix(0, 1, c)
  a[1] <- -1e5
  a[c + 1] <- 1e5
  a[2:c] <- 0
  
  indx <- c:(2 * (c - 1))                 # Define the subscript of marker effects from c to 2*(c-1), where 1:c is for the intercepts
  maxiter <- 100
  minerr <- 1e-8
  err <- 1e8
  iter <- 0
  
  b0 <- matrix(0, 2 * (c - 1), 1)         # Initialize the fixed parameters from 1 to 2*(c-1)
  g <- matrix(0, n * (c - 1), 1)          # Initialize the polygenic effects with n*(c-1) elements
  gamma <- matrix(0, c - 1, 1)            # Define the marker effects (subset of b0)
  poly <- matrix(g, n, c - 1, byrow = TRUE) # Reshape the polygenic effects n rows and c-1 columns
  out <- NULL
  
  while ((iter < maxiter) & (err > minerr)) {
    rr <- NULL
    pe <- NULL
    lp <- NULL
    pp <- list()
    
    for (j in 1:n) {
      poly.j <- as.matrix(poly[j, ])                         # poly.j is a c-1 vector of linear predictor for j
      z.j <- as.matrix(z[j] %x% diag(c - 1)) %*% gamma      # The added marker effects (different from the null polygenic model)
      
      tau <- 1
      for (k in 1:(c - 1)) {
        tau <- tau + exp(a[k + 1] + z.j[k] + poly.j[k])     # With extra marker effects
      }
      
      for (k in 1:(c - 1)) {
        mu[k] <- exp(a[k + 1] + z.j[k] + poly.j[k]) / tau   # With extra marker effects
      }
      mu[c] <- 1 / tau                                      # Expectation of the last category
      
      d <- matrix(0, c, c - 1)                              # Define the partial derivative matrix
      for (k in 1:(c - 1)) {
        for (i in 1:(c - 1)) {
          d[k, i] <- -exp(a[k + 1] + a[i + 1] + z.j[k] + z.j[i] + poly.j[k] + poly.j[i])
          if (k == i) {
            d[k, i] <- d[k, i] + tau * exp(a[k + 1] + z.j[k] + poly.j[k])
          }
        }
      }
      
      for (i in 1:(c - 1)) {
        d[c, i] <- -exp(a[i + 1] + z.j[i] + poly.j[i])
      }
      
      d <- d / tau^2
      w <- diag(1 / c(mu))                                  # Generalized inverse of var(y)
      dwd <- solve(t(d) %*% w %*% d + diag(c - 1) * 1e-5)   # Covariance structure of the linear pseudo response
      
      lp <- rbind(lp, (diag(c - 1) %*% as.matrix(a[2:c]) + z.j + poly.j))  # Linear predictor
      pe <- rbind(pe, dwd %*% t(d) %*% w %*% as.matrix(c(y[j, ] - mu)))     # Linear residual
      pp[[j]] <- dwd                                                           # Covariance structure stored in a list
    }
    
    ps <- pe + lp                                        # Pseudo response variable
    rr <- as.matrix(Matrix::bdiag(pp))                   # Diagonalize the covariance structure
    vi <- solve(kk %x% diag(par) + rr)                   # V**-1 = (K*par + rr)**-1
    xx <- t(x) %*% vi %*% x                              # X`V**-1X
    xy <- t(x) %*% vi %*% ps                             # X`V**-1y
    xxi <- solve(xx)                                     # (X`V**-1X)**-1
    b <- xxi %*% xy                                      # All fixed effects, including intercepts and marker effects
    
    poly <- (kk %x% diag(par)) %*% vi %*% (ps - x %*% b) # Predicted polygenic effects
    poly <- matrix(poly, n, c - 1, byrow = TRUE)         # Polygenic effects reshaped into n row and c-1 columns
    
    a[2:c] <- b[1:(c - 1)]                               # Assign the first c-1 values to intercepts
    gamma <- b[indx]                                     # Assign the last c-1 effects into marker effects
    vb <- xxi[indx, indx]                                # Variance matrix of the marker effects
    wald <- t(gamma) %*% solve(vb) %*% gamma             # Wald test statistic
    p <- 1 - pchisq(wald, c - 1)                         # P-value of the test statistic with (c-1) df
    
    iter <- iter + 1
    err <- sum((b - b0)^2)
    b0 <- b
    out <- rbind(out, cbind(iter, err, t(b), wald, p))   # Collect iteration history data
  }
  
  niter <- nrow(out)                                     # Find the last iteration
  iter.out <- out[niter, ]                               # Take the last iteration output
  response <- list(ps, rr, iter.out)                     # Create output list
  return(response)                                       # Return output, the last [[3]] stores the effects, the test and the p value
}

fun <- function(theta, kk, rr, x, ps) {    # Parameter theta is a real number from negative infinity to positive infinity
  s2 <- diag(exp(theta))                   # Exponential of the real parameters gives the positive variance parameters
  G <- kk %x% s2                           # Polygenic covariance matrix
  v <- G + rr                              # Var(pseudo response)
  d1 <- determinant(v)[[1]]                # Log determinant of matrix V
  vi <- solve(v)                           # Inverse of V matrix
  xx <- t(x) %*% vi %*% x                  # X`V**-1X
  xy <- t(x) %*% vi %*% ps                 # X`V**-1y
  d2 <- determinant(xx)[[1]]               # Log determinant of X`V**-1X
  b <- solve(xx) %*% xy                    # Fixed effects
  yy <- t(ps - x %*% b) %*% vi %*% (ps - x %*% b)   # Quadratic form
  like <- 0.5 * d1 + 0.5 * d2 + 0.5 * yy   # Negative log likelihood function for minimization
  return(like)
}

nominal.exact.scan <- function(y, zz, kk, x0, output_file = NULL) {
  # start_time <- Sys.time()                 # Record the starting time
  
  c <- ncol(y)
  m <- nrow(zz)
  
  www <- NULL
  error_snps <- c()                        # Vector to store indices of SNPs that cause errors
  
  for (k in 1:m) {
    # cat("Processing SNP", k, "\n")
    tryCatch({
      z <- as.matrix(zz[k, ]) %x% diag(c - 1)
      x <- cbind(x0, z)
      
      theta0 <- rep(0, c - 1)              # Initial value for the parameters in real values
      par0 <- exp(theta0)                  # Convert the real parameters into positive variance parameters
      
      double.iter <- NULL                  # History of doubly iterative algorithm
      maxiter <- 100                       # Define the maximum number of iterations
      minerr <- 1e-8                       # Define the convergence criterion as an error
      err <- 1e8                           # Starts with a large error to start the iteration
      iter <- 0                            # Record iteration
      
      while ((iter < maxiter) & (err > minerr)) {   # Outer loop
        myResult <- iblup.k(par = par0, y = y, z = z, x = x, kk = kk)
        ps <- myResult[[1]]                # Pseudo response variable
        rr <- myResult[[2]]                # Residual covariance structure
        
        obj <- function(theta) {
          fun(theta = theta, kk = kk, rr = rr, x = x, ps = ps)
        }
        
        parm <- optim(
          par = theta0,
          fn = obj,
          method = "L-BFGS-B",
          lower = -1e5,
          upper = 1e5
        )
        
        theta <- parm$par                  # Returned parameters (-Inf to Inf)
        par <- exp(theta)                  # Transform the real parameters into positive variance parameters
        err <- sum((par - par0)^2)         # Update error
        iter <- iter + 1                   # Record iterations
        theta0 <- theta                    # Change the initial value for optim
        par0 <- par                        # Assign new parameters to old parameters
        
        double.iter <- rbind(double.iter, c(iter, err, par))   # Record iteration history
      }
      
      www <- rbind(www, c(k, par, myResult[[3]]))
      
    }, error = function(e) {
      cat("Error at SNP", k, ":", conditionMessage(e), "\n")
      error_snps <<- c(error_snps, k)
    })
  }
  
  # end_time <- Sys.time()                   # Recording the end time
  # elapsed_time <- end_time - start_time
  # print(elapsed_time)
  
  www <- as.data.frame(www)
  colnames(www) <- c("SNP", "V1", "V2", "V3", "Iter", "Error",
                     "a1", "a2", "a3", "effect1", "effect2", "effect3", "Wald", "p")
  
  if (!is.null(output_file)) {
    write.csv(x = www, file = output_file, row.names = FALSE)
  }
  
  return(list(
    results = www,
    error_snps = error_snps,
    elapsed_time = elapsed_time
  ))
}
