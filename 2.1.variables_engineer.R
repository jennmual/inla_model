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

es_simetrico <- all(sapply(1:g_fix$n, function(i) {
  all(sapply(g_fix$nbs[[i]], function(j) i %in% g_fix$nbs[[j]]))
}))
cat("Grafo simétrico:", es_simetrico, "\n")
stopifnot(es_simetrico)

duplicados_internos <- any(sapply(g_fix$nbs, function(x) any(duplicated(x))))
cat("Duplicados internos:", duplicados_internos, "\n")
stopifnot(!duplicados_internos)

vecinos_fix <- sapply(1:g_fix$n, function(i) length(g_fix$nbs[[i]]))
cat("Nodos:", g_fix$n, "— Mín vecinos:", min(vecinos_fix),
    "— Máx:", max(vecinos_fix), "\n")
stopifnot(min(vecinos_fix) > 0)

g_clean <- g_fix

###############################################################################
# 2. Create INLA indices
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
# 3. VERIFICAR VARIABLES DISPONIBLES
###############################################################################
# Todas las covariables deben venir estandarizadas desde script 1.
# temp_lag, prec_lag, hum_lag → estandarizadas en script 1 (scale_global)
# elevation_z, pop_dens_z    → estandarizadas en script 1 (añadidas al final)

vars_requeridas <- c("temp_lag", "prec_lag", "hum_lag",
                     "elevation_z", "pop_dens_z",
                     "log_FOI_mtv", "log_FOI_kte")

faltantes <- setdiff(vars_requeridas, names(datos_modelo))
if (length(faltantes) > 0) {
  stop("Faltan estas columnas en datos_modelo: ",
       paste(faltantes, collapse = ", "),
       "\nVerifica que el script 1 las genera y guarda correctamente.")
}
cat("Todas las variables requeridas están disponibles.\n")

###############################################################################
# 4. CREAR VARIABLES DE INTERACCIÓN
###############################################################################
# Las covariables base ya están estandarizadas → las interacciones son
# productos de variables en la misma escala → coeficientes comparables.

datos_modelo <- datos_modelo |>
  mutate(
    # Temperatura × Precipitación
    # ¿El efecto del calor depende de la disponibilidad de lluvia?
    temp_x_prec = temp_lag * prec_lag,
    
    # Temperatura × Humedad
    # ¿La humedad potencia el efecto de la temperatura en el vector?
    temp_x_hum  = temp_lag * hum_lag,
    
    # Elevación × Temperatura
    # ¿El efecto de la temperatura difiere entre zonas altas y bajas?
    elev_x_temp = elevation_z * temp_lag,
    
    # Elevación × Precipitación
    # ¿La lluvia tiene distinto efecto según la altitud?
    elev_x_prec = elevation_z * prec_lag,
    
    # Densidad × Temperatura
    # ¿En ciudades más densas y calientes el dengue se amplifica?
    dens_x_temp = pop_dens_z * temp_lag,
    
    # Densidad × Precipitación
    # ¿La lluvia tiene mayor efecto en zonas urbanas densas?
    dens_x_prec = pop_dens_z * prec_lag
  )

cat("Variables de interacción creadas.\n")

###############################################################################
# 5. FUNCIONES AUXILIARES
###############################################################################

# Construir fórmula INLA genérica
make_formula <- function(respuesta, covariables = NULL) {
  cov_str <- if (!is.null(covariables)) {
    paste(covariables, collapse = " + ")
  } else {
    ""
  }
  
  formula_str <- paste0(
    respuesta, " ~ 1",
    if (nchar(cov_str) > 0) paste0(" + ", cov_str),
    " + f(idx_espacio, model='besag', graph=g_clean, scale.model=TRUE,
           hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))",
    " + f(idx_tiempo, model='iid',
           hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))"
  )
  as.formula(formula_str)
}

# Ajustar modelo INLA
ajustar_inla <- function(formula, datos) {
  inla(
    formula,
    data              = datos,
    family            = "gaussian",
    control.predictor = list(compute = TRUE),
    control.compute   = list(dic = TRUE, waic = TRUE)
  )
}

# Extraer tabla de coeficientes formateada
extraer_coefs <- function(modelo, nombre_modelo) {
  modelo$summary.fixed |>
    tibble::rownames_to_column("variable") |>
    dplyr::select(
      variable,
      beta    = mean,
      sd      = sd,
      lower95 = `0.025quant`,
      upper95 = `0.975quant`
    ) |>
    dplyr::mutate(
      modelo        = nombre_modelo,
      significativo = ifelse(lower95 > 0 | upper95 < 0, "Yes", "No"),
      .before       = variable
    )
}

# Pseudo-R²
calc_R2 <- function(modelo, y_obs) {
  y_pred <- modelo$summary.fitted.values$mean
  ss_res <- sum((y_obs - y_pred)^2, na.rm = TRUE)
  ss_tot <- sum((y_obs - mean(y_obs, na.rm = TRUE))^2, na.rm = TRUE)
  round(1 - ss_res / ss_tot, 4)
}

###############################################################################
# 6. DEFINIR ESPECIFICACIONES DE MODELOS
###############################################################################
# Nombres de covariables tal como existen en datos_modelo:
# temp_lag, prec_lag, hum_lag  → climáticas (estandarizadas en script 1)
# elevation_z, pop_dens_z      → geográfica y demográfica (estand. en script 1)
# temp_x_prec, etc.            → interacciones (creadas arriba)

especificaciones <- list(
  
  # ── BASELINE ─────────────────────────────────────────────────────────────────
  "00_baseline"          = NULL,
  
  # ── UNIVARIABLES (todos válidos — no hay colinealidad con una sola variable) ──
  "01_temp"              = c("temp_lag"),
  "02_prec"              = c("prec_lag"),
  "03_hum"               = c("hum_lag"),
  "04_elevation"         = c("elevation_z"),
  "05_pop_density"       = c("pop_dens_z"),
  
  # ── DOS COVARIABLES — SIN pares correlacionados ──────────────────────────────
  "06_temp_prec"         = c("temp_lag", "prec_lag"),
  "07_temp_hum"          = c("temp_lag", "hum_lag"),
  "08_prec_hum"          = c("prec_lag", "hum_lag"),
  # "09_temp_elev" → EXCLUIDO: r(temp, elev) = -0.74
  "10_prec_elev"         = c("prec_lag", "elevation_z"),   # verificar r
  "11_temp_dens"         = c("temp_lag", "pop_dens_z"),
  "12_prec_dens"         = c("prec_lag", "pop_dens_z"),
  "13_elev_dens"         = c("elevation_z", "pop_dens_z"),
  
  # ── MODELO CLIMÁTICO COMPLETO (sin elevación) ────────────────────────────────
  "14_climate_full"      = c("temp_lag", "prec_lag", "hum_lag"),
  
  # ── CON ELEVACIÓN — solo sin temperatura ────────────────────────────────────
  # Elevación reemplaza a temperatura como proxy geográfico
  "15_elev_prec_hum"     = c("elevation_z", "prec_lag", "hum_lag"),
  "16_elev_prec_hum_dens"= c("elevation_z", "prec_lag", "hum_lag", "pop_dens_z"),
  # "18_climate_elev_dens" → EXCLUIDO: incluye temp + elev juntas
  
  # ── CON DENSIDAD POBLACIONAL + TEMPERATURA ───────────────────────────────────
  "17_climate_dens"      = c("temp_lag", "prec_lag", "hum_lag", "pop_dens_z"),
  
  # ── INTERACCIONES — solo entre variables no correlacionadas ─────────────────
  "19_temp_prec_inter"   = c("temp_lag", "prec_lag", "temp_x_prec"),
  "20_temp_hum_inter"    = c("temp_lag", "hum_lag",  "temp_x_hum"),
  # "21_elev_temp_inter" → EXCLUIDO: elev y temp correlacionadas
  "22_elev_prec_inter"   = c("elevation_z", "prec_lag", "elev_x_prec"),
  "23_dens_temp_inter"   = c("pop_dens_z", "temp_lag", "dens_x_temp"),
  "24_dens_prec_inter"   = c("pop_dens_z", "prec_lag", "dens_x_prec"),
  
  # ── MODELOS COMPLETOS CON INTERACCIONES ──────────────────────────────────────
  "25_climate_inter_tp"  = c("temp_lag", "prec_lag", "hum_lag", "temp_x_prec"),
  # "26_elev_inter_et" → EXCLUIDO: elev + temp
  "27_dens_inter_dt"     = c("temp_lag", "prec_lag", "hum_lag",
                             "pop_dens_z", "dens_x_temp"),
  
  # ── MODELOS COMPLETOS SIN PARES CORRELACIONADOS ──────────────────────────────
  # Con temperatura (sin elevación)
  "28_all_temp"          = c("temp_lag", "prec_lag", "hum_lag", "pop_dens_z"),
  
  # Con elevación (sin temperatura) — modelo alternativo
  "29_all_elev"          = c("elevation_z", "prec_lag", "hum_lag", "pop_dens_z"),
  
  # Con temperatura + interacciones (sin elevación)
  "30_all_temp_inter"    = c("temp_lag", "prec_lag", "hum_lag",
                             "pop_dens_z", "temp_x_prec", "dens_x_temp")
)
cat("Especificaciones definidas:", length(especificaciones), "\n")
cat("Total de modelos a ajustar:", length(especificaciones) * 2, "\n")

###############################################################################
# 7. AJUSTAR TODOS LOS MODELOS
###############################################################################

ajustar_lista <- function(respuesta, etiqueta) {
  cat("\n── Ajustando modelos", etiqueta, "──\n")
  lista <- lapply(seq_along(especificaciones), function(i) {
    nombre <- names(especificaciones)[i]
    covs   <- especificaciones[[i]]
    cat("  [", i, "/", length(especificaciones), "]", nombre, "\n")
    tryCatch(
      ajustar_inla(make_formula(respuesta, covs), datos_modelo),
      error = function(e) {
        cat("  ✗ Error:", e$message, "\n")
        NULL
      }
    )
  })
  names(lista) <- names(especificaciones)
  Filter(Negate(is.null), lista)
}

modelos_kte <- ajustar_lista("log_FOI_kte", "FOI constante")
modelos_mtv <- ajustar_lista("log_FOI_mtv", "FOI variante")

cat("\nModelos ajustados — kte:", length(modelos_kte),
    "| mtv:", length(modelos_mtv), "\n")

###############################################################################
# 8. TABLA DE COMPARACIÓN (WAIC + R²)
###############################################################################

construir_tabla_comparacion <- function(lista_modelos, y_obs, etiqueta) {
  do.call(rbind, lapply(names(lista_modelos), function(nm) {
    m <- lista_modelos[[nm]]
    data.frame(
      FOI    = etiqueta,
      Modelo = nm,
      n_cov  = nrow(m$summary.fixed) - 1,
      WAIC   = round(m$waic$waic, 2),
      DIC    = round(m$dic$dic,   2),
      R2     = calc_R2(m, y_obs),
      stringsAsFactors = FALSE
    )
  }))
}

tabla_kte <- construir_tabla_comparacion(
  modelos_kte, datos_modelo$log_FOI_kte, "constante"
)
tabla_mtv <- construir_tabla_comparacion(
  modelos_mtv, datos_modelo$log_FOI_mtv, "variante"
)

tabla_comparacion <- rbind(tabla_kte, tabla_mtv) |> arrange(FOI, WAIC)

cat("\n══ COMPARACIÓN DE MODELOS ══\n")
print(tabla_comparacion, row.names = FALSE)

###############################################################################
# 9. TABLA DE COEFICIENTES
###############################################################################

construir_tabla_coefs <- function(lista_modelos, etiqueta) {
  do.call(rbind, lapply(names(lista_modelos), function(nm) {
    extraer_coefs(lista_modelos[[nm]], nm) |>
      mutate(FOI = etiqueta)
  }))
}

tabla_coefs <- rbind(
  construir_tabla_coefs(modelos_kte, "constante"),
  construir_tabla_coefs(modelos_mtv, "variante")
) |>
  arrange(FOI, modelo, variable)

modelos_clave <- c("04_elevation", "02_prec", "03_hum", "05_pop_density", #UNIVARIADOS
                   "00_baseline", "16_elev_prec_hum_dens")

cat("\n══ COEFICIENTES — MODELOS CLAVE (FOI variante) ══\n")
tabla_coefs |>
  filter(FOI == "variante", modelo %in% modelos_clave,
         variable != "(Intercept)") |>
  select(modelo, variable, beta, lower95, upper95, significativo) |>
  print(row.names = FALSE)

cat("\n══ COEFICIENTES — MODELOS CLAVE (FOI constante) ══\n")
tabla_coefs |>
  filter(FOI == "constante", modelo %in% modelos_clave,
         variable != "(Intercept)") |>
  select(modelo, variable, beta, lower95, upper95, significativo) |>
  print(row.names = FALSE)

###############################################################################
# 10. DELTA WAIC RESPECTO AL BASELINE
###############################################################################

calcular_delta_waic <- function(tabla) {
  waic_base <- tabla |> filter(Modelo == "00_baseline") |> pull(WAIC)
  tabla |>
    mutate(delta_WAIC     = round(WAIC - waic_base, 2),
           mejor_que_base = delta_WAIC < 0) |>
    arrange(WAIC)
}

cat("\n══ Δ WAIC — FOI constante ══\n")
calcular_delta_waic(tabla_kte) |>
  select(Modelo, n_cov, WAIC, delta_WAIC, R2, mejor_que_base) |>
  print(row.names = FALSE)

cat("\n══ Δ WAIC — FOI variante ══\n")
calcular_delta_waic(tabla_mtv) |>
  select(Modelo, n_cov, WAIC, delta_WAIC, R2, mejor_que_base) |>
  print(row.names = FALSE)

###############################################################################
# 12. GUARDAR RESULTADOS
###############################################################################

saveRDS(modelos_kte,       "results/modelos_kte_todos.rds")
saveRDS(modelos_mtv,       "results/modelos_mtv_todos.rds")
saveRDS(tabla_comparacion, "results/tabla_comparacion_modelos.rds")
saveRDS(tabla_coefs,       "results/tabla_coeficientes_todos.rds")

write.csv(tabla_comparacion,
          "results/tabla_comparacion_modelos.csv", row.names = FALSE)
write.csv(tabla_coefs |> filter(variable != "(Intercept)"),
          "results/tabla_coeficientes_todos.csv", row.names = FALSE)

cat("\nResultados guardados en results/\n")

###############################################################################
# 12. SELECCIONAR MODELOS FINALES Y UNIVARIADOS
###############################################################################

modelo_final <- "16_elev_prec_hum_dens"

modelos_univariados <- c(
  "02_prec",
  "03_hum",
  "04_elevation",
  "05_pop_density"
)

# Modelo multivariado final
modelo_final_kte <- modelos_kte[[modelo_final]]
modelo_final_mtv <- modelos_mtv[[modelo_final]]

# Modelos univariados
modelos_uni_kte <- modelos_kte[modelos_univariados]
modelos_uni_mtv <- modelos_mtv[modelos_univariados]

# GUARDAR Modelos finales multivariados
saveRDS(modelo_final_kte, "results/modelo_final_kte.rds")
saveRDS(modelo_final_mtv, "results/modelo_final_mtv.rds")

# GUARDAR Modelos univariados
saveRDS(modelos_uni_kte, "results/modelos_univariados_kte.rds")
saveRDS(modelos_uni_mtv, "results/modelos_univariados_mtv.rds")



#########################################################
# Matriz de correlaciones entre todas las covariables
#########################################################

vars_cov <- c("temp_lag", "prec_lag", "hum_lag", "elevation_z", "pop_dens_z")

cor_matrix <- cor(
  datos_modelo[, vars_cov],
  use = "complete.obs"
)

round(cor_matrix, 3)

# Graphic (excluir pares con |r| > 0.6)
library(corrplot)
corrplot(cor_matrix, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45,
         title = "Correlation between covariates", mar = c(0,0,1,0))
