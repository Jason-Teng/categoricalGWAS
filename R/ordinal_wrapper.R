# ============================================================
# Ordinal GWAS wrapper
# ============================================================
# This file follows the same design as nominal_wrapper.R:
#   1. fit_ordinal_null() is the null-model interface.
#   2. ordinal_gwas() is the main GWAS interface.
#
# Required files:
#   ordinal_null.R          -> ordinal_make_Y(), fit_ordinal_null()
#   ordinal_genome_scan.R   -> scan_score_ordinal(), scan_p3d_ordinal(),
#                              scan_psr_ordinal(), scan_psrsd_ordinal(),
#                              scan_glm_ordinal()
#   ordinal_exact.R         -> scan_exact_ordinal()
# ============================================================

ordinal_gwas <- function(y,
                         zz,
                         kk = NULL,
                         x0 = NULL,
                         method = c("score", "p3d", "psr", "psrsd", "exact", "glm"),
                         null_method = c("pseudo"),
                         null_fit = NULL,
                         vc = NULL,
                         ids = NULL,
                         outdir = NULL,
                         maxiter = 100,
                         minerr = 1e-8) {

  allowed_methods <- c("score", "p3d", "psr", "psrsd", "exact", "glm")
  method <- match.arg(method, choices = allowed_methods, several.ok = TRUE)
  null_method <- match.arg(null_method)

  Yfull <- ordinal_make_Y(y)
  zz <- as.matrix(zz)

  n <- nrow(Yfull)
  c <- ncol(Yfull)
  q <- c - 1

  if (is.null(x0)) {
    x0 <- matrix(1, n, 1) %x% diag(q)
  }

  if (is.null(ids)) ids <- seq_len(n)

  methods_need_kk <- c("score", "p3d", "psr", "psrsd", "exact")
  methods_need_null <- c("score", "p3d", "psr", "psrsd")
  methods_no_null <- c("glm", "exact")

  if (any(method %in% methods_need_kk) && is.null(kk)) {
    stop("'kk' is required for ordinal methods: score, p3d, psr, psrsd, and exact.")
  }

  if (any(method %in% methods_no_null) && (!is.null(null_fit) || !is.null(vc))) {
    message("Note: glm and exact do not use null model inputs. 'null_fit' and 'vc' are ignored for those methods.")
  }

  if (any(method %in% methods_need_null) && is.null(null_fit)) {
    null_fit <- fit_ordinal_null(
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

  results <- list()
  times <- list()
  error_snps <- list()

  for (one_method in method) {

    start_time <- Sys.time()
    message(paste0(">>> Start running ordinal method (", one_method, ")"))

    if (one_method == "score") {
      res <- scan_score_ordinal(
        ps = null_fit$ps,
        rr = null_fit$rr,
        par = null_fit$par,
        kk = kk,
        zz = zz,
        x0 = x0,
        c = c,
        n = n,
        outdir = outdir
      )
    }

    if (one_method == "p3d") {
      res <- scan_p3d_ordinal(
        y = Yfull,
        kk = kk,
        zz = zz,
        x0 = x0,
        par = null_fit$par,
        c = c,
        n = n,
        maxiter = maxiter,
        minerr = minerr,
        outdir = outdir
      )
    }

    if (one_method == "psr") {
      res <- scan_psr_ordinal(
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
    }

    if (one_method == "psrsd") {
      res <- scan_psrsd_ordinal(
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
    }

    if (one_method == "glm") {
      res <- scan_glm_ordinal(
        zz = zz,
        y = Yfull,
        x0 = x0,
        c = c,
        n = n,
        outdir = outdir
      )
    }

    if (one_method == "exact") {
      res <- scan_exact_ordinal(
        y = Yfull,
        kk = kk,
        zz = zz,
        x0 = x0,
        c = c,
        n = n,
        theta0 = 0,
        maxiter = maxiter,
        minerr = minerr,
        outdir = outdir
      )
      error_snps[[one_method]] <- attr(res, "error_snps")
    }

    elapsed_time <- Sys.time() - start_time
    message(paste("GWAS scanning time (", one_method, "):", format(elapsed_time)))

    results[[one_method]] <- res
    times[[one_method]] <- elapsed_time
  }

  if (length(method) == 1) {
    return(list(
      trait_type = "ordinal",
      method = method,
      result = results[[method]],
      error_snps = error_snps[[method]],
      null_fit = if (method %in% methods_need_null) null_fit else NULL,
      gwas_scanning_time = times[[method]]
    ))
  }

  list(
    trait_type = "ordinal",
    method = method,
    results = results,
    error_snps = error_snps,
    null_fit = if (any(method %in% methods_need_null)) null_fit else NULL,
    gwas_scanning_time = times
  )
}
