###############################################################################
# 0. Clean environment and load libraries
###############################################################################

rm(list = ls())

library(dplyr)
library(INLA)

###############################################################################
# 1. Load data and models
###############################################################################

datos_modelo <- readRDS("scripts/INLA_models/data/datos_modelo_indexed.rds")
g            <- readRDS("scripts/INLA_models/data/grafo_inla.rds")

# Load Model B formula (re-define to avoid dependency on previous session)
formula_B <- log_FOI_mtv ~ 1 +
  temp_lag + prec_lag + hum_lag +
  f(idx_espacio,
    model = "bym2", graph = g, scale.model = TRUE,
    hyper = list(
      phi  = list(prior = "pc", param = c(0.5, 2/3)),
      prec = list(prior = "pc.prec", param = c(1, 0.01))
    )) +
  f(idx_tiempo, model = "iid",
    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01))))

###############################################################################
# 2. Load posterior samples of FOI from Stan model
###############################################################################

# Expected structure: list or array where you can extract
# a vector of FOI values (one per municipio-year row in datos_modelo)
# for each posterior sample s

# --- ADJUST THIS BLOCK to match your Stan output structure ---
samples_raw <- readRDS("data/data_tesis_dani/lambda_samples_mtv.RDS")

# Verify dimensions match your data
# Each sample s should give a vector of length nrow(datos_modelo)
cat("Rows in datos_modelo:", nrow(datos_modelo), "\n")
cat("Dimensions of samples:", dim(samples_raw), "\n")  # should be [n_samples x nrow]

###############################################################################
# 3. Uncertainty propagation loop
###############################################################################

n_samples <- 200  # start with 200 to test; increase to 500 for final run

resultados_muestras <- vector("list", n_samples)

for (s in 1:n_samples) {
  
  datos_s <- datos_modelo
  
  # Extract sample s — adjust indexing to your structure:
  foi_s <- samples_raw[s, ]   # if matrix [n_samples × n_rows]
  # foi_s <- samples_raw[[s]] # if list of vectors
  
  # Safety check: skip if any zero or negative values
  if (any(foi_s <= 0, na.rm = TRUE)) {
    warning(paste("Sample", s, "has non-positive FOI values — skipping"))
    next
  }
  
  datos_s$log_FOI_s <- log(foi_s)
  
  mod_s <- inla(
    update(formula_B, log_FOI_s ~ .),
    data   = datos_s,
    family = "gaussian",
    control.compute = list(config = FALSE),
    # Speed up: no CPO/WAIC needed per sample
    control.inla = list(int.strategy = "eb")  # empirical Bayes for speed
  )
  
  resultados_muestras[[s]] <- mod_s$summary.fixed
  
  if (s %% 25 == 0) cat("Sample", s, "of", n_samples, "completed\n")
}

# Remove NULLs from skipped samples
resultados_muestras <- Filter(Negate(is.null), resultados_muestras)
cat("Valid samples used:", length(resultados_muestras), "\n")

###############################################################################
# 4. Aggregate results across samples
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
# 5. Save results
###############################################################################

saveRDS(resultados_muestras, "results/uncertainty_samples.rds")
saveRDS(coefs_propagados,    "results/coefs_propagados.rds")