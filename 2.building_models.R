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

#------------------------------------------------------------------------------------------

###############################################################################
# 3. CONSTANT FOI MODELS (log_FOI_kte)
###############################################################################
# These models evaluate the association between spatial-temporal structure
# and climate covariates using constant FOI estimates.

#------------------------------------------------------------------------------
# 3.1 Baseline model (space + time only)
#------------------------------------------------------------------------------

formula_A_kte <- log_FOI_kte ~ 1 +
  f(idx_espacio,
    model = "besag",
    graph = g_clean,
    scale.model = TRUE,
    hyper = list(
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    )) +
  f(idx_tiempo,
    model = "iid",
    hyper = list(
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    ))

modelo_A_kte <- inla(
  formula_A_kte,
  data = datos_modelo,
  family = "gaussian",
  control.compute = list(dic = TRUE, waic = TRUE)
)

#------------------------------------------------------------------------------
# 3.2 Full climate model
#------------------------------------------------------------------------------

formula_B_kte <- log_FOI_kte ~ 1 +
  temp_lag + prec_lag + hum_lag +
  f(idx_espacio,
    model = "besag",
    graph = g_clean,
    scale.model = TRUE,
    hyper = list(
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    )) +
  f(idx_tiempo,
    model = "iid",
    hyper = list(
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    ))

modelo_B_kte <- inla(
  formula_B_kte,
  data = datos_modelo,
  family = "gaussian",
  control.compute = list(dic = TRUE, waic = TRUE)
)

#------------------------------------------------------------------------------
# 3.3 Univariate climate models
#------------------------------------------------------------------------------

vars_clima <- c("temp_lag", "prec_lag", "hum_lag")

make_formula_uni_kte <- function(v) {
  as.formula(paste0(
    "log_FOI_kte ~ 1 + ", v, " + ",
    "f(idx_espacio, model='besag', graph=g_clean, scale.model=TRUE,
       hyper=list(prec=list(prior='pc.prec', param=c(1,0.01)))) + ",
    "f(idx_tiempo, model='iid',
       hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))"
  ))
}

modelos_uni_kte <- lapply(vars_clima, function(v) {
  inla(
    make_formula_uni_kte(v),
    data = datos_modelo,
    family = "gaussian",
    control.compute = list(dic = TRUE, waic = TRUE)
  )
})

names(modelos_uni_kte) <- vars_clima


###############################################################################
# 4. TIME-VARYING FOI MODELS (log_FOI_mtv)
###############################################################################
# Same structure, but using time-varying FOI estimates.

#------------------------------------------------------------------------------
# 4.1 Baseline model
#------------------------------------------------------------------------------

formula_A_mtv <- update(formula_A_kte, log_FOI_mtv ~ .)

modelo_A_mtv <- inla(
  formula_A_mtv,
  data = datos_modelo,
  family = "gaussian",
  control.compute = list(dic = TRUE, waic = TRUE)
)

#------------------------------------------------------------------------------
# 4.2 Full climate model
#------------------------------------------------------------------------------

formula_B_mtv <- update(formula_B_kte, log_FOI_mtv ~ .)

modelo_B_mtv <- inla(
  formula_B_mtv,
  data = datos_modelo,
  family = "gaussian",
  control.compute = list(dic = TRUE, waic = TRUE)
)

#------------------------------------------------------------------------------
# 4.3 Univariate climate models
#------------------------------------------------------------------------------

make_formula_uni_mtv <- function(v) {
  update(make_formula_uni_kte(v), log_FOI_mtv ~ .)
}

modelos_uni_mtv <- lapply(vars_clima, function(v) {
  inla(
    make_formula_uni_mtv(v),
    data = datos_modelo,
    family = "gaussian",
    control.compute = list(dic = TRUE, waic = TRUE)
  )
})

names(modelos_uni_mtv) <- vars_clima


###############################################################################
# 5. ELEVATION MODELS
###############################################################################
# Same model structure, replacing temperature with elevation.
# This evaluates whether static topography explains FOI variability.

vars_elev <- c("elevation", "prec_lag", "hum_lag")

#------------------------------------------------------------------------------
# 5.1 Constant FOI with elevation
#------------------------------------------------------------------------------

formula_elev_kte <- log_FOI_kte ~ 1 +
  elevation + prec_lag + hum_lag +
  f(idx_espacio, model = "besag", graph = g_clean, scale.model = TRUE) +
  f(idx_tiempo, model = "iid")

modelo_elev_kte <- inla(
  formula_elev_kte,
  data = datos_modelo,
  family = "gaussian",
  control.compute = list(dic = TRUE, waic = TRUE)
)

#------------------------------------------------------------------------------
# 5.2 Time-varying FOI with elevation
#------------------------------------------------------------------------------

formula_elev_mtv <- update(formula_elev_kte, log_FOI_mtv ~ .)

modelo_elev_mtv <- inla(
  formula_elev_mtv,
  data = datos_modelo,
  family = "gaussian",
  control.compute = list(dic = TRUE, waic = TRUE)
)


###############################################################################
# 6. PSEUDO-R² FUNCTION
###############################################################################

calc_R2 <- function(modelo, y_obs) {
  y_pred <- modelo$summary.fitted.values$mean
  ss_res <- sum((y_obs - y_pred)^2, na.rm = TRUE)
  ss_tot <- sum((y_obs - mean(y_obs, na.rm = TRUE))^2, na.rm = TRUE)
  round(1 - (ss_res / ss_tot), 4)
}


###############################################################################
# 7. MODEL COMPARISON TABLE
###############################################################################

comparacion_final <- data.frame(
  Model = c(
    "Baseline_kte",
    "Climate_kte",
    "Elevation_kte",
    "Baseline_mtv",
    "Climate_mtv",
    "Elevation_mtv"
  ),
  WAIC = c(
    modelo_A_kte$waic$waic,
    modelo_B_kte$waic$waic,
    modelo_elev_kte$waic$waic,
    modelo_A_mtv$waic$waic,
    modelo_B_mtv$waic$waic,
    modelo_elev_mtv$waic$waic
  ),
  R2 = c(
    calc_R2(modelo_A_kte, datos_modelo$log_FOI_kte),
    calc_R2(modelo_B_kte, datos_modelo$log_FOI_kte),
    calc_R2(modelo_elev_kte, datos_modelo$log_FOI_kte),
    calc_R2(modelo_A_mtv, datos_modelo$log_FOI_mtv),
    calc_R2(modelo_B_mtv, datos_modelo$log_FOI_mtv),
    calc_R2(modelo_elev_mtv, datos_modelo$log_FOI_mtv)
  )
)

print(comparacion_final)


