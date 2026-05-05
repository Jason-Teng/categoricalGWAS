## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)


library(usethis)

# ---- Locate raw data ----
gen_path <- "inst/extdata/IMF2-Genotypes.csv"
phe_path <- "inst/extdata/IMF2-Phenotypes.csv"
map_path <- "inst/extdata/map.csv"
kk_path  <- "inst/extdata/IMF2-KK.csv"


gen <- read.csv(gen_path)
phe <- read.csv(phe_path)
kk0 <- read.csv(kk_path)
map <- read.csv(map_path)

z <- as.matrix(gen[, -c(1:5)])
kk <- as.matrix(kk0[, -c(1, 2)])

nominal <- as.factor(phe$Nominal)
ordinal <- as.factor(phe$Ordinal)

# save FULL data if you want
usethis::use_data(
  z,
  kk,
  nominal,
  ordinal,
  overwrite = TRUE
)
