###############################################################################
# 0. Clean environment and load libraries
###############################################################################

rm(list = ls())

library(dplyr)
library(INLA)

inla.setOption(inla.mode = "classic")

###############################################################################
# 1. Load data and fix graph (same procedure as script 2)
###############################################################################

datos_modelo <- readRDS("data/datos_modelo_indexed.rds")
g            <- readRDS("data/grafo_inla.rds")

# Fix graph
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

# Verify
vecinos <- sapply(1:g_clean$n, function(i) length(g_clean$nbs[[i]]))
stopifnot(min(vecinos) > 0)
cat("Grafo OK — mín vecinos:", min(vecinos), "\n")

###############################################################################
# 2. Re-define Model B formula
# IMPORTANT: uses besag (not bym2) — consistent with script 2
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

###############################################################################
# 3. Load posterior samples of FOI from Stan model
###############################################################################

samples_raw <- readRDS("data/lambda_estimates_mtv.RDS")

cat("Filas en datos_modelo:", nrow(datos_modelo), "\n")
cat("Dimensiones samples:", dim(samples_raw), "\n")

# Verify alignment: n_cols should equal nrow(datos_modelo)
stopifnot(ncol(samples_raw) == nrow(datos_modelo))

###############################################################################
# 4. Uncertainty propagation loop
###############################################################################

n_samples <- 200  # increase to 500 for final run

resultados_muestras <- vector("list", n_samples)

for (s in 1:n_samples) {
  
  datos_s <- datos_modelo
  
  # Extract sample s
  foi_s <- samples_raw[s, ]
  
  # Skip if non-positive values
  if (any(foi_s <= 0, na.rm = TRUE)) {
    warning(paste("Sample", s, "has non-positive FOI values — skipping"))
    next
  }
  
  datos_s$log_FOI_s <- log(foi_s)
  
  mod_s <- tryCatch({
    inla(
      update(formula_B, log_FOI_s ~ .),
      data            = datos_s,
      family          = "gaussian",
      control.compute = list(config = FALSE),
      control.inla    = list(int.strategy = "eb")
    )
  }, error = function(e) {
    warning(paste("Sample", s, "failed:", e$message))
    NULL
  })
  
  resultados_muestras[[s]] <- mod_s$summary.fixed
  
  if (s %% 25 == 0) cat("Sample", s, "of", n_samples, "completed\n")
}

# Remove NULLs from skipped or failed samples
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
    quantile(sapply(resultados_muestras, function(r) r[v, "mean"]), 0.025)),
  
  q975 = sapply(coefs_nombres, function(v)
    quantile(sapply(resultados_muestras, function(r) r[v, "mean"]), 0.975))
) |>
  filter(variable != "(Intercept)")

print(coefs_propagados)

###############################################################################
# 6. Save results
###############################################################################

dir.create("results", showWarnings = FALSE)
saveRDS(resultados_muestras, "results/uncertainty_samples.rds")
saveRDS(coefs_propagados,    "results/coefs_propagados.rds")