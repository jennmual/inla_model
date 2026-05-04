###############################################################################
# 0. Clean environment and load required libraries
###############################################################################

rm(list = ls())

library(dplyr)
library(sf)
library(spdep)
library(INLA)
library(Matrix)

###############################################################################
# 1. Load prepared data (from Script 1)
###############################################################################

datos_modelo <- readRDS("data/datos_modelo.rds")
mun_estudio  <- readRDS("data/mun_estudio.rds")

###############################################################################
# 2. Ensure consistent spatial ordering (CRITICAL STEP)
###############################################################################

# The graph and indices MUST use the same ordering
mun_estudio <- mun_estudio |> arrange(cod_mun)

# Check uniqueness
stopifnot(length(unique(mun_estudio$cod_mun)) == nrow(mun_estudio))
###############################################################################
# 3. Build spatial adjacency graph (FINAL FIX)
###############################################################################

coords <- st_centroid(mun_estudio) |> st_coordinates()

nb <- knearneigh(coords, k = 4) |> knn2nb()

# Diagnostics
cat("Min neighbors:", min(sapply(nb, length)), "\n")
cat("Max neighbors:", max(sapply(nb, length)), "\n")

# 🔥 STEP 1: create adjacency file
nb2INLA("grafo_final.adj", nb)

# 🔥 STEP 2: move to clean path (Windows fix)
dir.create("C:/INLA_temp", showWarnings = FALSE)

file.copy("grafo_final.adj",
          "C:/INLA_temp/grafo_final.adj",
          overwrite = TRUE)

g <- inla.read.graph("C:/INLA_temp/grafo_final.adj")

###############################################################################
# 4. Create INLA indices (CRITICAL STEP)
###############################################################################

# IMPORTANT: must match the SAME ordering used in mun_estudio
municipios_unicos <- mun_estudio$cod_mun

datos_modelo <- datos_modelo |>
  mutate(
    idx_espacio = match(cod_mun, municipios_unicos),
    idx_tiempo  = as.integer(factor(year))
  )

# Sanity checks
stopifnot(!any(is.na(datos_modelo$idx_espacio)))
stopifnot(!any(is.na(datos_modelo$idx_tiempo))
)
stopifnot(max(datos_modelo$idx_espacio) == length(municipios_unicos))

cat("Spatial index range:",
    min(datos_modelo$idx_espacio), "-",
    max(datos_modelo$idx_espacio), "\n")

cat("Temporal index range:",
    min(datos_modelo$idx_tiempo), "-",
    max(datos_modelo$idx_tiempo), "\n")

###############################################################################
# 4. Create INLA indices (aligned with graph)
###############################################################################

# IMPORTANT: use the SAME ordering as mun_estudio
municipios_unicos <- mun_estudio$cod_mun

datos_modelo <- datos_modelo |>
  mutate(
    idx_espacio = match(cod_mun, municipios_unicos),
    idx_tiempo  = as.integer(factor(year))
  )

# Sanity checks
stopifnot(!any(is.na(datos_modelo$idx_espacio)))
stopifnot(!any(is.na(datos_modelo$idx_tiempo)))
stopifnot(max(datos_modelo$idx_espacio) == length(municipios_unicos))

cat("Spatial index range:",
    min(datos_modelo$idx_espacio), "-",
    max(datos_modelo$idx_espacio), "\n")

cat("Temporal index range:",
    min(datos_modelo$idx_tiempo), "-",
    max(datos_modelo$idx_tiempo), "\n")

###############################################################################
# 5. MODEL A: Baseline (no covariates)
###############################################################################

formula_A <- log_FOI_mtv ~ 1 +
  f(idx_espacio,
    model       = "bym2",
    graph       = g,   # ahora ES objeto
    scale.model = TRUE,
    hyper = list(
      phi  = list(prior = "pc", param = c(0.5, 2/3)),
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    )) +
  f(idx_tiempo,
    model = "iid",
    hyper = list(
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    ))

modelo_A <- inla(
  formula_A,
  data              = datos_modelo,
  family            = "gaussian",
  control.predictor = list(compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  control.inla      = list(int.strategy = "adaptive"),
  working.directory = getwd()   # 👈 NO tempdir()
)

summary(modelo_A)

###############################################################################
# 6. MODEL B: Multivariable (paper-level model)
###############################################################################

formula_B <- log_FOI_mtv ~ 1 +
  temp_lag + prec_lag + hum_lag +
  f(idx_espacio,
    model       = "bym2",
    graph       = g,
    scale.model = TRUE,
    hyper = list(
      phi  = list(prior = "pc", param = c(0.5, 2/3)),
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    )) +
  f(idx_tiempo,
    model = "iid",
    hyper = list(
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    ))

modelo_B <- inla(
  formula_B,
  data              = datos_modelo,
  family            = "gaussian",
  control.predictor = list(compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE),
  control.inla      = list(int.strategy = "adaptive")
)

summary(modelo_B)

###############################################################################
# 7. Save results
###############################################################################

saveRDS(modelo_A, "results/modelo_A.rds")
saveRDS(modelo_B, "results/modelo_B.rds")