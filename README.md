# categoricalGWAS

**categoricalGWAS** provides a unified framework for genome-wide association studies (GWAS) with **categorical phenotypes**, including both **nominal** and **ordinal** traits.

The package implements multiple methods under a consistent interface, allowing flexible comparison and scalable analysis.

---

## Installation

```r
install.packages("remotes")
remotes::install_github("Jason-Teng/categoricalGWAS")
library(categoricalGWAS)
```

---

## Overview

Supported methods:

* Score test
* P3D
* PSR
* PSRSD (ordinal only)
* GLM (baseline)
* Exact method

---

## Example Workflow

### Load data

```r
data(z)
data(kk)
data(nominal)
data(ordinal)
```

### Subset SNPs (for fast example)

```r
zz <- z[1:5, ]
```

---

## Ordinal GWAS

```r
res_ord <- categorical_gwas(
  y = ordinal,
  zz = zz,
  kk = kk,
  trait_type = "ordinal",
  method = c("score", "psrsd")
)

res_ord$results$score
res_ord$results$psrsd
```

---

## Nominal GWAS

```r
res_nom <- categorical_gwas(
  y = nominal,
  zz = zz,
  kk = kk,
  trait_type = "nominal",
  method = c("score", "psr")
)

res_nom$results$score
```

---

## Input Format

* `z`
  Genotype matrix (**markers × individuals**)

* `kk`
  Kinship matrix (**individuals × individuals**)

* `nominal`
  Factor vector

* `ordinal`
  Factor vector

---

## Raw Data Access

```r
gen_path <- system.file("extdata", "IMF2-Genotypes.csv", package = "categoricalGWAS")
phe_path <- system.file("extdata", "IMF2-Phenotypes.csv", package = "categoricalGWAS")

gen <- read.csv(gen_path)
phe <- read.csv(phe_path)
```

---

## Notes

* Examples use a small subset of SNPs (`z[1:5, ]`) for speed
* Full datasets can be used for large-scale GWAS
* Methods requiring a null model will fit it automatically

---

## Author

Chin-Sheng Teng
Ph.D. Candidate, Applied Statistics
University of California, Riverside

---
