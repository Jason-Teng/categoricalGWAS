#' Run categorical GWAS
#'
#' @param y Phenotype object.
#' @param zz Marker genotype matrix, markers in rows and individuals in columns.
#' @param kk Kinship matrix.
#' @param trait_type Either "nominal" or "ordinal".
#' @param method Character vector of methods.
#' @return A list of GWAS results.
#' @export
categorical_gwas <- function(y,
                             zz,
                             kk = NULL,
                             x0 = NULL,
                             trait_type = c("ordinal", "nominal"),
                             method = NULL,
                             null_method = NULL,
                             null_fit = NULL,
                             vc = NULL,
                             ids = NULL,
                             outdir = NULL,
                             maxiter = 100,
                             minerr = 1e-8) {

  trait_type <- match.arg(trait_type)

  if (trait_type == "ordinal") {
    if (is.null(method)) {
      method <- c("score", "p3d", "psr", "psrsd", "exact", "glm")
    }
    if (is.null(null_method)) {
      null_method <- "pseudo"
    }

    return(ordinal_gwas(
      y = y,
      zz = zz,
      kk = kk,
      x0 = x0,
      method = method,
      null_method = null_method,
      null_fit = null_fit,
      vc = vc,
      ids = ids,
      outdir = outdir,
      maxiter = maxiter,
      minerr = minerr
    ))
  }

  if (trait_type == "nominal") {
    if (is.null(method)) {
      method <- c("score", "p3d", "psr", "exact", "glm")
    }
    if (is.null(null_method)) {
      null_method <- "pseudo"
    }

    # The current nominal_gwas() accepts only one method at a time.
    # This dispatcher allows multiple nominal methods by looping.
    allowed_nominal <- c("score", "p3d", "psr", "exact", "glm")
    method <- match.arg(method, choices = allowed_nominal, several.ok = TRUE)

    nominal_results <- list()
    nominal_times <- list()
    nominal_error_snps <- list()
    shared_null_fit <- null_fit

    for (one_method in method) {
      fit <- nominal_gwas(
        y = y,
        zz = zz,
        kk = kk,
        x0 = x0,
        method = one_method,
        null_method = null_method,
        null_fit = shared_null_fit,
        vc = vc,
        ids = ids,
        outdir = outdir,
        maxiter = maxiter,
        minerr = minerr
      )

      nominal_results[[one_method]] <- fit$result
      nominal_times[[one_method]] <- fit$gwas_scanning_time
      nominal_error_snps[[one_method]] <- fit$error_snps

      if (is.null(shared_null_fit) && !is.null(fit$null_fit)) {
        shared_null_fit <- fit$null_fit
      }
    }

    if (length(method) == 1) {
      return(list(
        trait_type = "nominal",
        method = method,
        result = nominal_results[[method]],
        error_snps = nominal_error_snps[[method]],
        null_fit = if (method %in% c("score", "p3d", "psr")) shared_null_fit else NULL,
        gwas_scanning_time = nominal_times[[method]]
      ))
    }

    return(list(
      trait_type = "nominal",
      method = method,
      results = nominal_results,
      error_snps = nominal_error_snps,
      null_fit = if (any(method %in% c("score", "p3d", "psr"))) shared_null_fit else NULL,
      gwas_scanning_time = nominal_times
    ))
  }
}
