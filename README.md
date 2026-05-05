# categoricalGWAS

**categoricalGWAS** provides a unified framework for genome-wide association studies (GWAS) with **categorical phenotypes**, including both **nominal** and **ordinal** traits.

The package implements multiple methods under a consistent interface, allowing flexible comparison and scalable analysis.

---

## Supported methods

The package currently supports several GWAS testing strategies under a common interface.

| Method | Description | Nominal | Ordinal |
|---|---|---:|---:|
| `score` | Score test using the fitted null model | Yes | Yes |
| `p3d` | Population parameters previously determined | Yes | Yes |
| `psr` | Pseudo-response based marker scan | Yes | Yes |
| `psrsd` | Pseudo-response method with ordinal-specific structure | No | Yes |
| `glm` | Generalized linear model baseline without kinship correction | Yes | Yes |
| `exact` | Marker-specific model fitting | Yes | Yes |

Some methods require a null model, while others can be run directly. When needed, `categoricalGWAS` fits the null model automatically unless a precomputed null model or variance component estimate is supplied.

---

## Installation

```r
install.packages("remotes")
remotes::install_github("Jason-Teng/categoricalGWAS")
library(categoricalGWAS)
```

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

## Output structure

The main function returns a list. The most important component is usually `results`.

```r
names(res_nom)
names(res_nom$results)
```

Each method returns a data frame containing marker-level test results. The exact columns may depend on the selected method, but commonly include:

| Column | Meaning |
|---|---|
| `SNP` | Marker index |
| `Effect` or category-specific effects | Estimated marker effect |
| `StdErr` | Standard error |
| `Wald` or test statistic | Marker-level test statistic |
| `pvalue` | Marker-level p-value |
| `iter` | Number of iterations, when applicable |
| `err` | Convergence error, when applicable |

For example:

```r
head(res_nom$results$score)
```

---

## Citation

If you use `categoricalGWAS`, please cite the related manuscript once available.

For now, you may cite the GitHub repository:

```text
Teng, C.-S. categoricalGWAS: Genome-wide association studies for categorical phenotypes.
GitHub repository: https://github.com/Jason-Teng/categoricalGWAS
```

---

## Author

Chin-Sheng Teng 
University of California, Riverside

---
