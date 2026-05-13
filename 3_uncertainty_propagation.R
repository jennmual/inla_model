###############################################################################
# 0. Clean environment and load libraries
###############################################################################

rm(list = ls())

library(dplyr)
library(INLA)

inla.setOption(inla.mode = "classic")

###############################################################################
# 1. Load data and fix graph
###############################################################################


datos_modelo <- readRDS("data/datos_modelo.rds")
g            <- readRDS("data/grafo_inla.rds")

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
g_clean <- g_fix
stopifnot(min(sapply(1:g_clean$n, function(i) length(g_clean$nbs[[i]]))) > 0)
cat("Grafo OK\n")

# ── Recrear índices SIEMPRE en esta sesión ─────────────────────────────────
municipios_unicos <- sort(unique(datos_modelo$cod_mun))
años_unicos       <- sort(unique(datos_modelo$year))

datos_modelo <- datos_modelo |>
  arrange(cod_mun, year) |>
  mutate(
    idx_espacio = as.integer(factor(cod_mun, levels = municipios_unicos)),
    idx_tiempo  = as.integer(factor(year,    levels = años_unicos))
  )

stopifnot(!any(is.na(datos_modelo$idx_espacio)))
stopifnot(!any(is.na(datos_modelo$idx_tiempo)))
stopifnot(max(datos_modelo$idx_espacio) == g_clean$n)

cat("Índices OK — espacio: 1 a", max(datos_modelo$idx_espacio),
    "| tiempo: 1 a", max(datos_modelo$idx_tiempo), "\n")
cat("Filas en datos_modelo:", nrow(datos_modelo), "\n")

###############################################################################
# 2. Define Model B formula
# Use the same covariates as your final selected model
###############################################################################

# Adjust vars_modelo to match your final model
vars_modelo <- c("prec_lag", "hum_lag", "elevation_z", "pop_dens_z")

formula_B <- as.formula(paste0(
  "log_FOI_s ~ 1 + ",
  paste(vars_modelo, collapse = " + "),
  " + f(idx_espacio, model='besag', graph=g_clean, scale.model=TRUE,
         hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))",
  " + f(idx_tiempo, model='iid',
         hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))"
))

###############################################################################
# 3. Align FOI estimates with datos_modelo
###############################################################################

samples_raw <- readRDS("data/lambda_estimates_mtv.RDS")

years_modelo <- sort(unique(datos_modelo$year))

# Alinear por cod_mun × year en el MISMO orden que datos_modelo
samples_aligned <- datos_modelo |>
  select(cod_mun, year) |>
  left_join(
    samples_raw |>
      rename(cod_mun = location) |>
      filter(year %in% years_modelo) |>
      select(cod_mun, year, median, lower, upper),
    by = c("cod_mun", "year")
  )

# Verificar alineación
cat("Filas alineadas:", nrow(samples_aligned), "\n")
cat("Filas datos_modelo:", nrow(datos_modelo), "\n")
stopifnot(nrow(samples_aligned) == nrow(datos_modelo))

# Verificar que cod_mun y year coinciden fila a fila
stopifnot(all(samples_aligned$cod_mun == datos_modelo$cod_mun))
stopifnot(all(samples_aligned$year    == datos_modelo$year))

n_missing <- sum(is.na(samples_aligned$median))
cat("Filas sin FOI después del join:", n_missing, "\n")
if (n_missing > 0) {
  cat("Municipio-años sin match:\n")
  samples_aligned |>
    filter(is.na(median)) |>
    select(cod_mun, year) |>
    print(n = 20)
  stop("Hay filas sin FOI — verifica que todos los cod_mun × year de datos_modelo
        están en samples_raw.")
}

# Parámetros log-normales
samples_aligned <- samples_aligned |>
  mutate(
    mu_log    = log(median),
    sigma_log = pmax((log(upper) - log(lower)) / (2 * 1.96), 1e-4)
  )

cat("Parámetros log-normales OK\n")
print(summary(samples_aligned[, c("mu_log", "sigma_log")]))

# Simular muestras
n_samples <- 200
set.seed(123)

muestras_log <- matrix(
  rnorm(
    n    = n_samples * nrow(datos_modelo),
    mean = rep(samples_aligned$mu_log,    each = n_samples),
    sd   = rep(samples_aligned$sigma_log, each = n_samples)
  ),
  nrow = n_samples,
  ncol = nrow(datos_modelo),
  byrow = FALSE
)

cat("Matriz de muestras:", dim(muestras_log), "\n")


###############################################################################
# 4. Uncertainty propagation loop
###############################################################################

vars_modelo <- c("prec_lag", "hum_lag", "elevation_z", "pop_dens_z")

formula_B <- as.formula(paste0(
  "log_FOI_s ~ 1 + ",
  paste(vars_modelo, collapse = " + "),
  " + f(idx_espacio, model='besag', graph=g_clean, scale.model=TRUE,
         hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))",
  " + f(idx_tiempo, model='iid',
         hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))"
))

# Verificar que todas las variables están en datos_modelo ANTES del loop
vars_requeridas <- c(vars_modelo, "idx_espacio", "idx_tiempo", "log_FOI_s")
# log_FOI_s se añadirá en el loop — el resto debe existir ya
vars_check <- c(vars_modelo, "idx_espacio", "idx_tiempo")
faltantes  <- setdiff(vars_check, names(datos_modelo))
if (length(faltantes) > 0) {
  stop("Faltan columnas en datos_modelo: ", paste(faltantes, collapse = ", "))
}
cat("Todas las variables requeridas están en datos_modelo.\n")

resultados_muestras <- vector("list", n_samples)

for (s in 1:n_samples) {
  
  # Añadir log_FOI_s a datos_modelo — NO crear un objeto nuevo
  # para evitar perder columnas
  datos_modelo$log_FOI_s <- muestras_log[s, ]
  
  if (any(!is.finite(datos_modelo$log_FOI_s))) {
    warning(paste("Sample", s, "has non-finite values — skipping"))
    next
  }
  
  mod_s <- tryCatch({
    inla(
      formula_B,
      data            = datos_modelo,      # usar datos_modelo directamente
      family          = "gaussian",
      control.compute = list(config = FALSE),
      control.inla    = list(int.strategy = "eb")
    )
  }, error = function(e) {
    warning(paste("Sample", s, "failed:", e$message))
    NULL
  })
  
  if (!is.null(mod_s)) {
    resultados_muestras[[s]] <- mod_s$summary.fixed
  }
  
  if (s %% 25 == 0) cat("Sample", s, "of", n_samples, "completed\n")
}

# Limpiar la columna temporal
datos_modelo$log_FOI_s <- NULL

resultados_muestras <- Filter(Negate(is.null), resultados_muestras)
cat("Valid samples used:", length(resultados_muestras), "\n")


###############################################################################
# 5. Aggregate results across samples
###############################################################################

coefs_nombres <- rownames(resultados_muestras[[1]])

coefs_propagados <- data.frame(
  variable = coefs_nombres,
  
  media = sapply(coefs_nombres, function(v)
    mean(sapply(resultados_muestras, function(r) r[v, "mean"]))),
  
  sd = sapply(coefs_nombres, function(v)
    sd(sapply(resultados_muestras, function(r) r[v, "mean"]))),
  
  q025 = sapply(coefs_nombres, function(v)
    quantile(sapply(resultados_muestras,
                    function(r) r[v, "mean"]), 0.025)),
  
  q975 = sapply(coefs_nombres, function(v)
    quantile(sapply(resultados_muestras,
                    function(r) r[v, "mean"]), 0.975))
) |>
  filter(variable != "(Intercept)")

cat("\nCoeficientes con propagación de incertidumbre:\n")
print(coefs_propagados)

###############################################################################
# 6. Compare with point estimate model (no propagation)
###############################################################################

# Load the model fitted on median FOI for comparison
modelo_B_punto <- readRDS("results/modelo_B.rds")

coefs_punto <- modelo_B_punto$summary.fixed |>
  as.data.frame() |>
  tibble::rownames_to_column("variable") |>
  filter(variable != "(Intercept)") |>
  select(variable,
         media   = mean,
         q025    = `0.025quant`,
         q975    = `0.975quant`) |>
  mutate(tipo = "Point estimate (median FOI)")

coefs_propagados_plot <- coefs_propagados |>
  mutate(tipo = "Propagated uncertainty")

cat("\n── Point estimate vs propagated uncertainty ──\n")
cat("\nPoint estimate:\n")
print(coefs_punto[, c("variable", "media", "q025", "q975")])
cat("\nPropagated:\n")
print(coefs_propagados[, c("variable", "media", "q025", "q975")])

###############################################################################
# 7. Save results
###############################################################################

dir.create("results", showWarnings = FALSE)
saveRDS(resultados_muestras, "results/uncertainty_samples.rds")
saveRDS(coefs_propagados,    "results/coefs_propagados.rds")
saveRDS(coefs_punto,         "results/coefs_punto.rds")

write.csv(coefs_propagados,
          "results/coefs_propagados.csv", row.names = FALSE)

cat("\nResultados guardados en results/\n")
