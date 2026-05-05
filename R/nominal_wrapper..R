
fit_nominal_null <- function(y, x0, kk, ids = NULL,
                             null_method = c("pseudo", "laplace", NULL),
                             vc = NULL,
                             maxiter = 100,
                             minerr = 1e-8,
                             outdir = NULL) {

  null_method <- match.arg(null_method)

  Yfull <- nominal_make_Y(y)
  n <- nrow(Yfull)
  c <- ncol(Yfull)
  q <- c - 1

  if (!is.null(vc)) {

    message(paste0(">>> Start running Null Model Linearization with provided VC)"))
    null_start <- Sys.time()

    par_val <- vc
    res <- iblup_nominal(par = par_val, y = Yfull, x = x0, kk = kk)

    message(paste("Null Model Time:", format(Sys.time() - null_start)))

    # Explicitly structure the return for the GWAS scan
    return(list(
      par = par_val,
      ps = res$ps,
      rr = res$rr
    ))
  }

  if (null_method == "pseudo") {

    message(paste0(">>> Start fitting Null Model (pseudo)"))
    null_start <- Sys.time()

    if (is.null(ids)) ids <- seq_len(nrow(Yfull))

    fit <- fit_nominal_outer(
      y = Yfull,
      x = x0,
      kk = kk,
      ids = ids,
      iblup_func = iblup_nominal,
      fun_func = nominal_fun,
      maxiter = maxiter,
      minerr = minerr,
      outdir = ifelse(is.null(outdir), tempdir(), outdir)
    )

    message(paste("Null Model Time (pseudo):", format(Sys.time() - null_start)))

    return(list(
      method = "pseudo",
      par = fit$vc,
      theta = fit$theta,
      ps = fit$ps,
      rr = fit$rr,
      trace = fit$trace,
      elapsed_time = fit$elapsed_time
    ))
  }

  if (null_method == "laplace") {

    message(paste0(">>> Start fitting Null Model (laplace)"))
    null_start <- Sys.time()

    # Call the main fitting function from laplace_functions.R
    lap <- fit_nominal_laplace_null(
      y = Yfull,
      kk = kk,
      verbose = FALSE
    )

    # Use the estimated VC (lap$vc) to generate pseudo-data via IBLUP
    pseudo <- iblup_nominal(
      par = lap$vc,
      y = lap$Y,
      x = x0,
      kk = kk,
      maxiter = maxiter,
      minerr = minerr
    )

    message(paste("Null Model Time (laplace):", format(Sys.time() - null_start)))

    # Return the results structured similarly to the pseudo method
    return(list(
      method       = "laplace",
      par          = lap$vc,           # Variance components
      beta         = pseudo$beta,      # Fixed effects
      theta        = lap$theta,        # Log-variance components
      ps           = pseudo$ps,        # Pseudo-response vector
      rr           = pseudo$rr,        # Working weight matrix
      trace        = lap$mode_trace,   # Saved process (inner loop history)
      elapsed_time = lap$elapsed_time  # Timing
    ))
  }
}


nominal_gwas <- function(y, zz, kk = NULL, x0 = NULL,
                         method = c("score", "p3d", "psr", "exact", "glm"),
                         null_method = c("pseudo", "laplace"),
                         null_fit = NULL, # Option 1: Provide full null object
                         vc = NULL,       # Option 2: Provide just variance components
                         ids = NULL,
                         outdir = NULL,
                         maxiter = 100,
                         minerr = 1e-8) {

  method <- match.arg(method)
  null_method <- match.arg(null_method)

  Yfull <- nominal_make_Y(y)
  n <- nrow(Yfull)
  c <- ncol(Yfull)
  q <- c - 1

  if (is.null(x0)) {
    x0 <- model.matrix(~ as.factor(rep(1:(c - 1), n)) - 1)
  }

  # start_time <- Sys.time()

  if (method %in% c("glm", "exact")) {

    if (!is.null(null_fit) || !is.null(vc)) {
      message("Note: '", method, "' does not use null model inputs. 'null_fit' and 'vc' are ignored.")
    }

  }

  if (method == "glm") {

    start_time <- Sys.time()
    message(paste0(">>> Start running method (glm)"))
    res <- scan_glm_nominal(
      zz = zz,
      y = Yfull,
      x0 = x0
    )

    elapsed_time <- Sys.time() - start_time
    message(paste("GWAS scanning Time (glm):", format(elapsed_time)))

    return(list(
      method = "glm",
      result = res,
      null_fit = NULL,
      gwas_scanning_time = elapsed_time
    ))
  }

  if (method == "exact") {

    start_time <- Sys.time()
    message(paste0(">>> Start running method (exact)"))
    res <- nominal.exact.scan(
      y = Yfull,
      zz = zz,
      kk = kk,
      x0 = x0,
      output_file = NULL
    )

    elapsed_time <- Sys.time() - start_time
    message(paste("GWAS scanning Time (exact):", format(elapsed_time)))

    return(list(
      method = "exact",
      result = res$results,
      error_snps = res$error_snps,
      null_fit = NULL,
      gwas_scanning_time = elapsed_time
    ))
  }

  if (is.null(null_fit)) { # if provided vc, the time would be here
    null_fit <- fit_nominal_null(
      y = Yfull,
      x0 = x0,
      kk = kk,
      ids = ids,
      null_method = null_method,
      vc = vc,
      maxiter = maxiter,
      minerr = minerr,
      outdir = outdir
    )
  }

  start_time <- Sys.time()

  if (method == "score") {

    message(paste0(">>> Start running method (score)"))
    res <- scan_score_nominal(
      ps = null_fit$ps,
      rr = null_fit$rr,
      par = null_fit$par,
      kk = kk,
      zz = zz,
      x0 = x0,
      c = c,
      outdir = outdir
    )
    elapsed_time = Sys.time() - start_time
    message(paste("GWAS scanning Time (exact):", format(elapsed_time)))
  }

  if (method == "p3d") {

    message(paste0(">>> Start running method (p3d)"))
    res <- scan_p3d_nominal(
      y = Yfull,
      kk = kk,
      zz = zz,
      x0 = x0,
      par = null_fit$par,
      c = c,
      n = n,
      outdir = outdir
    )
    elapsed_time = Sys.time() - start_time
    message(paste("GWAS scanning Time (exact):", format(elapsed_time)))
  }

  if (method == "psr") {

    message(paste0(">>> Start running method (psr)"))
    res <- scan_psr_nominal(
      ps = null_fit$ps,
      rr = null_fit$rr,
      kk = kk,
      zz = zz,
      x0 = x0,
      c = c,
      n = n,
      theta0 = log(null_fit$par),
      outdir = outdir
    )
    elapsed_time = Sys.time() - start_time
    message(paste("GWAS scanning Time (exact):", format(elapsed_time)))
  }

  list(
    method = method,
    result = res,
    null_fit = null_fit,
    gwas_scanning_time = elapsed_time
  )

}

