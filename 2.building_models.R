###############################################################################
# 0. Clean environment and load required libraries
###############################################################################

rm(list = ls())

library(dplyr)
library(INLA)

inla.setOption(inla.mode = "classic")

###############################################################################
# 1. Load prepared data and fix spatial graph
###############################################################################

datos_modelo <- readRDS("data/datos_modelo.rds")
g            <- readRDS("data/grafo_inla.rds")

# Fix graph: enforce symmetry and remove internal duplicates
g_fix <- g

for (i in 1:g_fix$n) {
  g_fix$nbs[[i]] <- sort(unique(g_fix$nbs[[i]]))
  g_fix$nnbs[i]  <- length(g_fix$nbs[[i]])
}

for (i in 1:g_fix$n) {
  for (j in g_fix$nbs[[i]]) {
    if (!(i %in% g_fix$nbs[[j]])) {
      g_fix$nbs[[j]] <- sort(unique(c(g_fix$nbs[[j]], i)))
      g_fix$nnbs[j]  <- length(g_fix$nbs[[j]])
    }
  }
}

# Validate: symmetry
es_simetrico <- all(sapply(1:g_fix$n, function(i) {
  all(sapply(g_fix$nbs[[i]], function(j) i %in% g_fix$nbs[[j]]))
}))
cat("Grafo simétrico:", es_simetrico, "\n")
stopifnot(es_simetrico)

# Validate: no internal duplicates within each neighbor list
duplicados_internos <- any(sapply(g_fix$nbs, function(x) any(duplicated(x))))
cat("Duplicados internos:", duplicados_internos, "\n")
stopifnot(!duplicados_internos)

# Validate: neighbors exist
vecinos_fix <- sapply(1:g_fix$n, function(i) length(g_fix$nbs[[i]]))
cat("Nodos:", g_fix$n, "\n")
cat("Mín vecinos:", min(vecinos_fix), "— Máx:", max(vecinos_fix), "\n")
stopifnot(min(vecinos_fix) > 0)

# inla.read.graph does not work on this system — use g_fix directly
g_clean <- g_fix

###############################################################################
# 2. Create INLA indices (space and time)
###############################################################################

municipios_unicos <- sort(unique(datos_modelo$cod_mun))
años_unicos       <- sort(unique(datos_modelo$year))

datos_modelo <- datos_modelo |>
  mutate(
    idx_espacio = as.integer(factor(cod_mun, levels = municipios_unicos)),
    idx_tiempo  = as.integer(factor(year,    levels = años_unicos))
  )

stopifnot(!any(is.na(datos_modelo$idx_espacio)))
stopifnot(!any(is.na(datos_modelo$idx_tiempo)))
stopifnot(max(datos_modelo$idx_espacio) == g_clean$n)

cat("Índices espacio: 1 a", max(datos_modelo$idx_espacio), "\n")
cat("Índices tiempo:  1 a", max(datos_modelo$idx_tiempo),  "\n")

###############################################################################
# 3. MODEL A: Baseline (no covariates)
# - Response: log FOI time-varying
# - Spatial random effect: BESAG
# - Temporal random effect: IID
###############################################################################

formula_A <- log_FOI_mtv ~ 1 +
  f(idx_espacio,
    model       = "besag",
    graph       = g_clean,
    scale.model = TRUE,
    hyper = list(
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
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
)

summary(modelo_A)

###############################################################################
# 4. MODEL B: Full multivariable model
# - Adds lagged and standardized climate covariates
# - Same spatial + temporal structure as Model A
###############################################################################

formula_B <- log_FOI_mtv ~ 1 +
  temp_lag + prec_lag + hum_lag +
  f(idx_espacio,
    model       = "besag",
    graph       = g_clean,
    scale.model = TRUE,
    hyper = list(
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
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

summary(modelo_B)
summary(modelo_B)$fixed

###############################################################################
# 5. UNIVARIATE MODELS (one covariate at a time)
###############################################################################

vars_clima <- c("temp_lag", "prec_lag", "hum_lag")

make_formula_uni <- function(v) {
  as.formula(paste0(
    "log_FOI_mtv ~ 1 + ", v, " + ",
    "f(idx_espacio, model = 'besag', graph = g_clean, scale.model = TRUE,
       hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01)))) + ",
    "f(idx_tiempo, model = 'iid',
       hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))"
  ))
}

modelos_uni <- lapply(vars_clima, function(v) {
  inla(
    make_formula_uni(v),
    data            = datos_modelo,
    family          = "gaussian",
    control.compute = list(dic = TRUE, waic = TRUE)
  )
})

names(modelos_uni) <- vars_clima

###############################################################################
# 6. MODEL COMPARISON (DIC and WAIC)
###############################################################################

comparacion <- data.frame(
  modelo = c("A_base", "B_multivariable", names(modelos_uni)),
  DIC = c(
    modelo_A$dic$dic,
    modelo_B$dic$dic,
    sapply(modelos_uni, function(m) m$dic$dic)
  ),
  WAIC = c(
    modelo_A$waic$waic,
    modelo_B$waic$waic,
    sapply(modelos_uni, function(m) m$waic$waic)
  )
) |>
  arrange(WAIC)

print(comparacion, row.names = FALSE)

###############################################################################
# 7. SAVE RESULTS
###############################################################################

saveRDS(datos_modelo, "data/datos_modelo_indexed.rds")
saveRDS(g_clean,      "data/grafo_clean.rds")
saveRDS(modelo_A,     "results/modelo_A.rds")
saveRDS(modelo_B,     "results/modelo_B.rds")
saveRDS(modelos_uni,  "results/modelos_univariables.rds")
saveRDS(comparacion,  "results/comparacion_modelos.rds")

###############################################################################
# 8. EXTENSION: Repeat models A and B with constant FOI
###############################################################################

formula_A_kte <- update(formula_A, log_FOI_kte ~ .)
formula_B_kte <- update(formula_B, log_FOI_kte ~ .)

modelo_A_kte <- inla(
  formula_A_kte,
  data            = datos_modelo,
  family          = "gaussian",
  control.compute = list(dic = TRUE, waic = TRUE)
)

modelo_B_kte <- inla(
  formula_B_kte,
  data            = datos_modelo,
  family          = "gaussian",
  control.compute = list(dic = TRUE, waic = TRUE)
)
summary(modelo_B_kte)

saveRDS(modelo_A_kte, "results/modelo_A_kte.rds")
saveRDS(modelo_B_kte, "results/modelo_B_kte.rds")
