###############################################################################
# 0. Clean environment and load libraries
###############################################################################

rm(list = ls())

library(dplyr)
library(ggplot2)
library(tibble)
library(purrr)
library(sf)
library(knitr)
library(tidyverse)
library(INLA)

inla.setOption(inla.mode = "classic")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1. Load all results
###############################################################################

datos_modelo  <- readRDS("data/datos_modelo.rds")
mun_estudio   <- readRDS("data/mun_estudio.rds")

modelo_A      <- readRDS("results/modelo_A.rds")
modelo_B      <- readRDS("results/modelo_B.rds")
modelo_A_kte  <- readRDS("results/modelo_A_kte.rds")
modelo_B_kte  <- readRDS("results/modelo_B_kte.rds")

# Univariables FOI variante (script anterior — tiene temp_lag, prec_lag, hum_lag)
modelos_uni     <- readRDS("results/modelos_univariables.rds")

# Univariables FOI constante (script nuevo — tiene 02_prec, 03_hum, etc.)
modelos_uni_kte <- readRDS("results/modelos_univariados_kte.rds")

coefs_propagados <- tryCatch(
  readRDS("results/coefs_propagados.rds"),
  error = function(e) { message("Propagated coefs not found — skipping"); NULL }
)

municipios_unicos <- sort(unique(datos_modelo$cod_mun))

# Variables del modelo final seleccionado
# Ajusta esto si tu mejor modelo tiene otras covariables
vars_modelo <- c("prec_lag", "hum_lag", "elevation_z", "pop_dens_z")

etiquetas_modelo <- c(
  "prec_lag"    = "Precipitation",
  "hum_lag"     = "Relative humidity",
  "elevation_z" = "Elevation",
  "pop_dens_z"  = "Population density"
)

###############################################################################
# 2. Fix spatial graph (needed to fit missing univariable models)
###############################################################################

g     <- readRDS("data/grafo_inla.rds")
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

cat("Grafo OK — nodos:", g_clean$n,
    "| mín vecinos:", min(sapply(1:g_clean$n, function(i) length(g_clean$nbs[[i]]))),
    "\n")


###############################################################################
# 2b. Verificar y reajustar modelos finales si sus variables no coinciden
###############################################################################

vars_en_B <- rownames(modelo_B$summary.fixed)
vars_en_B <- vars_en_B[vars_en_B != "(Intercept)"]
cat("Variables actuales en modelo_B:", vars_en_B, "\n")
cat("Variables requeridas:          ", vars_modelo, "\n")

if (!setequal(vars_en_B, vars_modelo)) {
  cat("modelo_B no coincide con vars_modelo — reajustando...\n")
  
  formula_B_nueva <- as.formula(paste0(
    "log_FOI_mtv ~ 1 + ",
    paste(vars_modelo, collapse = " + "),
    " + f(idx_espacio, model='besag', graph=g_clean, scale.model=TRUE,
           hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))",
    " + f(idx_tiempo, model='iid',
           hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))"
  ))
  
  modelo_B <- inla(
    formula_B_nueva,
    data              = datos_modelo,
    family            = "gaussian",
    control.predictor = list(compute = TRUE),
    control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )
  
  saveRDS(modelo_B, "results/modelo_B.rds")
  cat("modelo_B reajustado y guardado.\n")
}

# Mismo para FOI constante
vars_en_B_kte <- rownames(modelo_B_kte$summary.fixed)
vars_en_B_kte <- vars_en_B_kte[vars_en_B_kte != "(Intercept)"]
cat("Variables actuales en modelo_B_kte:", vars_en_B_kte, "\n")

if (!setequal(vars_en_B_kte, vars_modelo)) {
  cat("modelo_B_kte no coincide con vars_modelo — reajustando...\n")
  
  formula_B_kte_nueva <- as.formula(paste0(
    "log_FOI_kte ~ 1 + ",
    paste(vars_modelo, collapse = " + "),
    " + f(idx_espacio, model='besag', graph=g_clean, scale.model=TRUE,
           hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))",
    " + f(idx_tiempo, model='iid',
           hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))"
  ))
  
  modelo_B_kte <- inla(
    formula_B_kte_nueva,
    data              = datos_modelo,
    family            = "gaussian",
    control.predictor = list(compute = TRUE),
    control.compute   = list(dic = TRUE, waic = TRUE)
  )
  
  saveRDS(modelo_B_kte, "results/modelo_B_kte.rds")
  cat("modelo_B_kte reajustado y guardado.\n")
}

cat("Variables finales modelo_B:    ",
    rownames(modelo_B$summary.fixed), "\n")
cat("Variables finales modelo_B_kte:",
    rownames(modelo_B_kte$summary.fixed), "\n")


###############################################################################
# 3. Create INLA indices (needed if not already in datos_modelo)
###############################################################################

if (!"idx_espacio" %in% names(datos_modelo)) {
  municipios_unicos_idx <- sort(unique(datos_modelo$cod_mun))
  años_unicos_idx       <- sort(unique(datos_modelo$year))
  datos_modelo <- datos_modelo |>
    mutate(
      idx_espacio = as.integer(factor(cod_mun, levels = municipios_unicos_idx)),
      idx_tiempo  = as.integer(factor(year,    levels = años_unicos_idx))
    )
}

###############################################################################
# 4. Complete univariable models for FOI variante
###############################################################################
# modelos_uni (from previous script) may not include elevation_z and pop_dens_z
# Fit the missing ones here

cat("\nChecking univariable models for time-varying FOI...\n")
cat("Available:", names(modelos_uni), "\n")
cat("Needed:   ", vars_modelo, "\n")

vars_faltantes_mtv <- setdiff(vars_modelo, names(modelos_uni))

if (length(vars_faltantes_mtv) > 0) {
  cat("Fitting missing univariable models:", vars_faltantes_mtv, "\n")
  
  modelos_uni_extra <- lapply(vars_faltantes_mtv, function(v) {
    cat("  Fitting:", v, "\n")
    f <- as.formula(paste0(
      "log_FOI_mtv ~ 1 + ", v,
      " + f(idx_espacio, model='besag', graph=g_clean, scale.model=TRUE,
             hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))",
      " + f(idx_tiempo, model='iid',
             hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))"
    ))
    tryCatch(
      inla(f, data = datos_modelo, family = "gaussian",
           control.compute = list(dic = TRUE, waic = TRUE)),
      error = function(e) { cat("  Error:", e$message, "\n"); NULL }
    )
  })
  names(modelos_uni_extra) <- vars_faltantes_mtv
  modelos_uni_extra        <- Filter(Negate(is.null), modelos_uni_extra)
  
  modelos_uni_mtv <- c(modelos_uni, modelos_uni_extra)
} else {
  modelos_uni_mtv <- modelos_uni
}

cat("Univariable models FOI variante:", names(modelos_uni_mtv), "\n")

# Mapeos: variable name → list element name
# FOI variante: names are the variable names directly
mapa_uni_mtv <- setNames(vars_modelo, vars_modelo)

# FOI constante: names from the new engineering script
mapa_uni_kte <- c(
  prec_lag    = "02_prec",
  hum_lag     = "03_hum",
  elevation_z = "04_elevation",
  pop_dens_z  = "05_pop_density"
)

# Final validation
stopifnot(all(mapa_uni_mtv %in% names(modelos_uni_mtv)))
stopifnot(all(mapa_uni_kte %in% names(modelos_uni_kte)))
cat("All univariable models validated.\n")

###############################################################################
# 5. Helper functions
###############################################################################

# Extract a specific variable's fixed effect from an INLA model
extraer_coef_var <- function(modelo, variable, tipo) {
  if (is.null(modelo) || is.null(modelo$summary.fixed)) return(NULL)
  if (!variable %in% rownames(modelo$summary.fixed))    return(NULL)
  
  modelo$summary.fixed |>
    as.data.frame() |>
    rownames_to_column("variable") |>
    filter(variable == !!variable) |>
    transmute(
      variable,
      media = mean,
      q025  = `0.025quant`,
      q975  = `0.975quant`,
      tipo  = tipo
    )
}

# Extract all non-intercept fixed effects from an INLA model
extraer_coefs <- function(modelo, tipo) {
  if (is.null(modelo) || is.null(modelo$summary.fixed)) return(NULL)
  modelo$summary.fixed |>
    as.data.frame() |>
    rownames_to_column("variable") |>
    filter(variable != "(Intercept)") |>
    transmute(
      variable,
      media = mean,
      q025  = `0.025quant`,
      q975  = `0.975quant`,
      tipo  = tipo
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
# 6. FIGURE 1a: Forest plot — Time-varying FOI, univariable vs multivariable
###############################################################################

coefs_uni_mtv <- map_dfr(
  vars_modelo,
  ~ extraer_coef_var(
    modelos_uni_mtv[[ mapa_uni_mtv[.x] ]],
    .x,
    "Univariable"
  )
)

coefs_multi_mtv <- map_dfr(
  vars_modelo,
  ~ extraer_coef_var(modelo_B, .x, "Multivariable")
)

coefs_forest_mtv <- bind_rows(coefs_uni_mtv, coefs_multi_mtv) |>
  mutate(variable = factor(variable, levels = rev(vars_modelo)))

fig1_mtv <- ggplot(coefs_forest_mtv,
                   aes(x = media, y = variable,
                       color = tipo, shape = tipo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbarh(aes(xmin = q025, xmax = q975),
                 position = position_dodge(width = 0.6), height = 0.2) +
  scale_color_manual(values = c("Univariable"   = "#2166ac",
                                "Multivariable" = "#1a9641")) +
  scale_y_discrete(labels = etiquetas_modelo) +
  scale_x_continuous(limits = c(-0.2, 0.2)) +
  labs(x        = "Change in log(FOI)",
       y        = NULL,
       color    = NULL,
       shape    = NULL,
       title    = "Effect of variables on log(FOI)"
       #subtitle = "Time-varying FOI — univariable vs multivariable"
       ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(fig1_mtv)
ggsave("results/figures/fig1_forest_plot_mtv.png",
       fig1_mtv, width = 8, height = 5, dpi = 300)

###############################################################################
# 7. FIGURE 1b: Forest plot — Constant FOI, univariable vs multivariable
###############################################################################

coefs_uni_kte <- map_dfr(
  vars_modelo,
  ~ extraer_coef_var(
    modelos_uni_kte[[ mapa_uni_kte[.x] ]],
    .x,
    "Univariable"
  )
)

coefs_multi_kte <- map_dfr(
  vars_modelo,
  ~ extraer_coef_var(modelo_B_kte, .x, "Multivariable")
)

coefs_forest_kte <- bind_rows(coefs_uni_kte, coefs_multi_kte) |>
  mutate(variable = factor(variable, levels = rev(vars_modelo)))

fig1_kte <- ggplot(coefs_forest_kte,
                   aes(x = media, y = variable,
                       color = tipo, shape = tipo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbarh(aes(xmin = q025, xmax = q975),
                 position = position_dodge(width = 0.6), height = 0.2) +
  scale_color_manual(values = c("Univariable"   = "#2166ac",
                                "Multivariable" = "#1a9641")) +
  scale_y_discrete(labels = etiquetas_modelo) +
  scale_x_continuous(limits = c(-0.2, 0.2)) +
  labs(x        = "Change in log(FOI)",
       y        = NULL,
       color    = NULL,
       shape    = NULL,
       title    = "Effect of variables on log(FOI)"#,
       #subtitle = "Constant FOI — univariable vs multivariable"
       ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(fig1_kte)
ggsave("results/figures/fig1_forest_plot_kte.png",
       fig1_kte, width = 8, height = 5, dpi = 300)

###############################################################################
# 8. FIGURE 2: Forest plot — Time-varying vs Constant FOI (multivariable)
###############################################################################

coefs_comparacion <- bind_rows(
  map_dfr(vars_modelo,
          ~ extraer_coef_var(modelo_B,     .x, "Time-varying FOI")),
  map_dfr(vars_modelo,
          ~ extraer_coef_var(modelo_B_kte, .x, "Constant FOI"))
) |>
  mutate(variable = factor(variable, levels = rev(vars_modelo)))

fig2 <- ggplot(coefs_comparacion,
               aes(x = media, y = variable,
                   color = tipo, shape = tipo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(position = position_dodge(0.5), size = 3) +
  geom_errorbarh(aes(xmin = q025, xmax = q975),
                 position = position_dodge(0.5), height = 0.2) +
  scale_color_manual(values = c("Time-varying FOI" = "#d73027",
                                "Constant FOI"     = "#4575b4")) +
  scale_y_discrete(labels = etiquetas_modelo) +
  scale_x_continuous(limits = c(-0.2, 0.2)) +
  labs(x        = "Change in log(FOI)",
       y        = NULL,
       color    = "FOI type",
       shape    = "FOI type",
       #title    = "Comparison: Time-varying vs Constant FOI",
       subtitle = "Selected multivariable model") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(fig2)
ggsave("results/figures/fig2_forest_kte_vs_mtv.png",
       fig2, width = 8, height = 5, dpi = 300)

###############################################################################
# 9. FIGURE 3: Spatial relative risk map (SRR) — Time-varying FOI
###############################################################################

dpto_sf <- st_read("data/gadm41_COL.gpkg",
                   layer = "ADM_ADM_1", quiet = TRUE) |>
  st_transform(st_crs(mun_estudio))

SRR_df <- modelo_B$summary.random$idx_espacio |>
  mutate(
    cod_mun  = municipios_unicos,
    SRR      = exp(mean),
    SRR_q025 = exp(`0.025quant`),
    SRR_q975 = exp(`0.975quant`)
  )

mun_mapa <- mun_estudio |>
  left_join(SRR_df, by = "cod_mun")

fig3 <- ggplot() +
  geom_sf(data = dpto_sf,
          fill = "grey95", color = "grey60", linewidth = 0.12) +
  geom_sf(data = mun_mapa,
          aes(fill = SRR), color = "grey60", linewidth = 0.12) +
  scale_fill_viridis_c(option = "plasma", direction = -1, name = "SRR") +
  labs(title    = "Spatial Relative Risk of FOI",
       subtitle = "Time-varying FOI — 98 municipalities of Colombia") +
  theme_void(base_size = 12) +
  theme(legend.position = "right")

print(fig3)
ggsave("results/figures/fig3_mapa_SRR_mtv.png",
       fig3, width = 6, height = 8, dpi = 300)

###############################################################################
# 9b. FIGURE 3b: Spatial relative risk map (SRR) — Constant FOI
###############################################################################

SRR_df_kte <- modelo_B_kte$summary.random$idx_espacio |>
  mutate(
    cod_mun  = municipios_unicos,
    SRR      = exp(mean),
    SRR_q025 = exp(`0.025quant`),
    SRR_q975 = exp(`0.975quant`)
  )

mun_mapa_kte <- mun_estudio |>
  left_join(SRR_df_kte, by = "cod_mun")

fig3b <- ggplot() +
  geom_sf(data = dpto_sf,
          fill = "grey95", color = "grey60", linewidth = 0.12) +
  geom_sf(data = mun_mapa_kte,
          aes(fill = SRR), color = "grey60", linewidth = 0.12) +
  scale_fill_viridis_c(option = "plasma", direction = -1, name = "SRR") +
  labs(title    = "Spatial Relative Risk of FOI",
       subtitle = "Constant FOI — 98 municipalities of Colombia") +
  theme_void(base_size = 12) +
  theme(legend.position = "right")

print(fig3b)
ggsave("results/figures/fig3_mapa_SRR_kte.png",
       fig3b, width = 6, height = 8, dpi = 300)

###############################################################################
# 10. FIGURE 4: Observed vs fitted — Time-varying FOI
###############################################################################

datos_modelo$predicho_mtv <- modelo_B$summary.fitted.values$mean

fig4 <- ggplot(datos_modelo,
               aes(x = log_FOI_mtv, y = predicho_mtv)) +
  geom_point(alpha = 0.35, size = 1.5, color = "#2166ac") +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(x        = "Estimated log(FOI)",
       y        = "Predicted log(FOI)",
       title    = "Estimated vs fitted values"
       #subtitle = "Time-varying FOI — selected multivariable model"
       ) +
  theme_minimal(base_size = 13)

print(fig4)
ggsave("results/figures/fig4_observed_vs_fitted_mtv.png",
       fig4, width = 6, height = 5, dpi = 300)

###############################################################################
# 10b. FIGURE 4b: Observed vs fitted — Constant FOI
###############################################################################

datos_modelo$predicho_kte <- modelo_B_kte$summary.fitted.values$mean

fig4b <- ggplot(datos_modelo,
                aes(x = log_FOI_kte, y = predicho_kte)) +
  geom_point(alpha = 0.35, size = 1.5, color = "#4575b4") +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(x        = "Estimated log(FOI)",
       y        = "Predicted log(FOI)",
       title    = "Estimated vs fitted values",
       subtitle = "Constant FOI — selected multivariable model") +
  theme_minimal(base_size = 13)

print(fig4b)
ggsave("results/figures/fig4_observed_vs_fitted_kte.png",
       fig4b, width = 6, height = 5, dpi = 300)

###############################################################################
# 11. TABLE 1: Model comparison — Time-varying FOI
###############################################################################

etiquetas_uni_mtv <- setNames(
  etiquetas_modelo[vars_modelo],
  vars_modelo
)

tabla_mtv <- tibble(
  Model     = c("Baseline",
                "Multivariable",
                etiquetas_uni_mtv[names(modelos_uni_mtv)]),
  DIC       = c(modelo_A$dic$dic,
                modelo_B$dic$dic,
                sapply(modelos_uni_mtv, function(m) m$dic$dic)),
  WAIC      = c(modelo_A$waic$waic,
                modelo_B$waic$waic,
                sapply(modelos_uni_mtv, function(m) m$waic$waic)),
  p_eff     = c(modelo_A$dic$p.eff,
                modelo_B$dic$p.eff,
                sapply(modelos_uni_mtv, function(m) m$dic$p.eff)),
  R2        = c(calc_R2(modelo_A,   datos_modelo$log_FOI_mtv),
                calc_R2(modelo_B,   datos_modelo$log_FOI_mtv),
                sapply(modelos_uni_mtv, calc_R2,
                       y_obs = datos_modelo$log_FOI_mtv))
) |>
  arrange(WAIC)

cat("\n── Model comparison: Time-varying FOI ──\n")
print(kable(tabla_mtv, digits = 2,
            col.names = c("Model", "DIC", "WAIC", "Eff. param.", "R²")))

###############################################################################
# 12. TABLE 2: Model comparison — Constant FOI
###############################################################################

etiquetas_uni_kte <- c(
  "02_prec"        = "Precipitation",
  "03_hum"         = "Relative humidity",
  "04_elevation"   = "Elevation",
  "05_pop_density" = "Population density"
)

tabla_kte <- tibble(
  Model     = c("Baseline",
                "Multivariable",
                etiquetas_uni_kte[names(modelos_uni_kte)]),
  DIC       = c(modelo_A_kte$dic$dic,
                modelo_B_kte$dic$dic,
                sapply(modelos_uni_kte, function(m) m$dic$dic)),
  WAIC      = c(modelo_A_kte$waic$waic,
                modelo_B_kte$waic$waic,
                sapply(modelos_uni_kte, function(m) m$waic$waic)),
  p_eff     = c(modelo_A_kte$dic$p.eff,
                modelo_B_kte$dic$p.eff,
                sapply(modelos_uni_kte, function(m) m$dic$p.eff)),
  R2        = c(calc_R2(modelo_A_kte, datos_modelo$log_FOI_kte),
                calc_R2(modelo_B_kte, datos_modelo$log_FOI_kte),
                sapply(modelos_uni_kte, calc_R2,
                       y_obs = datos_modelo$log_FOI_kte))
) |>
  arrange(WAIC)

cat("\n── Model comparison: Constant FOI ──\n")
print(kable(tabla_kte, digits = 2,
            col.names = c("Model", "DIC", "WAIC", "Eff. param.", "R²")))

saveRDS(tabla_mtv, "results/tabla_comparacion_mtv.rds")
saveRDS(tabla_kte, "results/tabla_comparacion_kte.rds")
write.csv(tabla_mtv, "results/tabla_comparacion_mtv.csv", row.names = FALSE)
write.csv(tabla_kte, "results/tabla_comparacion_kte.csv", row.names = FALSE)

###############################################################################
# 13. TABLE 3: Fixed effect coefficients — both FOI types
###############################################################################

coefs_tabla <- bind_rows(
  extraer_coefs(modelo_B,     "Time-varying FOI") |> mutate(FOI = "variante"),
  extraer_coefs(modelo_B_kte, "Constant FOI")     |> mutate(FOI = "constante")
) |>
  mutate(
    significativo = ifelse(q025 > 0 | q975 < 0, "Yes", "No"),
    variable      = etiquetas_modelo[variable]
  ) |>
  select(FOI, variable, media, q025, q975, significativo)

cat("\n── Fixed effects coefficients ──\n")
print(kable(coefs_tabla, digits = 3,
            col.names = c("FOI", "Variable", "β (mean)",
                          "2.5%", "97.5%", "Significant")))

write.csv(coefs_tabla, "results/tabla_coeficientes_finales.csv",
          row.names = FALSE)

###############################################################################
# SECCIÓN 14 EXPANDIDA: Visualización de propagación de incertidumbre
###############################################################################

library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(purrr)

# Cargar resultados
coefs_propagados    <- readRDS("results/coefs_propagados.rds")
resultados_muestras <- readRDS("results/uncertainty_samples.rds")
modelo_B            <- readRDS("results/modelo_B.rds")

vars_modelo <- c("prec_lag", "hum_lag", "elevation_z", "pop_dens_z")

etiquetas_modelo <- c(
  "prec_lag"    = "Precipitation",
  "hum_lag"     = "Relative humidity",
  "elevation_z" = "Elevation",
  "pop_dens_z"  = "Population density"
)

###############################################################################
# 14.1 Extraer coeficientes del modelo puntual (sin propagación)
###############################################################################

coefs_punto <- modelo_B$summary.fixed |>
  as.data.frame() |>
  rownames_to_column("variable") |>
  filter(variable %in% vars_modelo) |>
  transmute(
    variable,
    media = mean,
    q025  = `0.025quant`,
    q975  = `0.975quant`,
    tipo  = "Point estimate\n(median FOI)"
  )

coefs_prop_plot <- coefs_propagados |>
  filter(variable %in% vars_modelo) |>
  select(variable, media, q025, q975) |>
  mutate(tipo = "Propagated\nuncertainty")

# Combinar
coefs_comp <- bind_rows(coefs_punto, coefs_prop_plot) |>
  mutate(
    variable = factor(variable, levels = rev(vars_modelo)),
    tipo     = factor(tipo, levels = c("Point estimate\n(median FOI)",
                                       "Propagated\nuncertainty"))
  )

###############################################################################
# 14.2 FIGURA 5a: Comparación punto vs propagado
###############################################################################

fig5a <- ggplot(coefs_comp,
                aes(x = media, y = variable,
                    color = tipo, shape = tipo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(aes(xmin = q025, xmax = q975),
                 position = position_dodge(width = 0.5),
                 height = 0.2) +
  scale_color_manual(values = c(
    "Point estimate\n(median FOI)" = "#2166ac",
    "Propagated\nuncertainty"      = "#d73027"
  )) +
  scale_y_discrete(labels = etiquetas_modelo) +
  labs(x        = "Change in log(FOI)",
       y        = NULL,
       color    = NULL,
       shape    = NULL,
       title    = "Model coefficients: point estimate vs propagated uncertainty",
       subtitle = paste0("Propagated: ", length(resultados_muestras),
                         " simulated FOI samples")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(fig5a)
ggsave("results/figures/fig5a_propagated_vs_point.png",
       fig5a, width = 8, height = 5, dpi = 300)

###############################################################################
# 14.3 FIGURA 5b: Distribución posterior de cada coeficiente (violin/density)
# Reconstruye la distribución completa de β a partir de las 200 corridas
###############################################################################

# Extraer los 200 valores de β para cada variable
dist_coefs <- map_dfr(vars_modelo, function(v) {
  vals <- sapply(resultados_muestras, function(r) {
    if (v %in% rownames(r)) r[v, "mean"] else NA
  })
  data.frame(
    variable = v,
    beta     = vals[!is.na(vals)]
  )
}) |>
  mutate(variable = factor(variable, levels = rev(vars_modelo)))

fig5b <- ggplot(dist_coefs, aes(x = beta, y = variable, fill = variable)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  ggridges::geom_density_ridges(
    alpha = 0.7, scale = 0.9,
    quantile_lines = TRUE, quantiles = c(0.025, 0.975)
  ) +
  scale_y_discrete(labels = etiquetas_modelo) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x        = "β — Change in log(FOI)",
       y        = NULL,
       title    = "Posterior distribution of coefficients",
       subtitle = paste0("Across ", length(resultados_muestras),
                         " simulated FOI realizations\nVertical lines = 2.5% and 97.5% quantiles")) +
  theme_minimal(base_size = 13)

print(fig5b)
ggsave("results/figures/fig5b_posterior_distributions.png",
       fig5b, width = 8, height = 5, dpi = 300)

# Si no tienes ggridges instalado:
# install.packages("ggridges")
# O versión alternativa con boxplot:
fig5b_alt <- ggplot(dist_coefs, aes(x = beta, y = variable)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(fill = "#d73027", alpha = 0.4, width = 0.4,
               outlier.shape = 16, outlier.size = 1) +
  geom_point(data = coefs_prop_plot |>
               mutate(variable = factor(variable, levels = rev(vars_modelo))),
             aes(x = media), color = "#d73027", size = 3) +
  scale_y_discrete(labels = etiquetas_modelo) +
  labs(x        = "β — Change in log(FOI)",
       y        = NULL,
       title    = "Distribution of β across 200 FOI simulations",
       subtitle = "Box = IQR | Whiskers = 95% | Point = mean") +
  theme_minimal(base_size = 13)

print(fig5b_alt)
ggsave("results/figures/fig5b_alt_boxplot.png",
       fig5b_alt, width = 7, height = 4, dpi = 300)

###############################################################################
# 14.4 FIGURA 5c: Convergencia — β medio acumulado por número de muestras
# Muestra que 200 muestras son suficientes (la media se estabiliza)
###############################################################################

conv_data <- map_dfr(vars_modelo, function(v) {
  vals <- sapply(resultados_muestras, function(r) {
    if (v %in% rownames(r)) r[v, "mean"] else NA
  })
  vals <- vals[!is.na(vals)]
  data.frame(
    variable   = v,
    n_sample   = seq_along(vals),
    media_acum = cumsum(vals) / seq_along(vals)
  )
}) |>
  mutate(variable = factor(variable, levels = vars_modelo))

fig5c <- ggplot(conv_data,
                aes(x = n_sample, y = media_acum, color = variable)) +
  geom_line(linewidth = 0.8) +
  geom_hline(data = coefs_prop_plot |>
               mutate(variable = factor(variable, levels = vars_modelo)),
             aes(yintercept = media, color = variable),
             linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ variable, scales = "free_y",
             labeller = labeller(variable = etiquetas_modelo)) +
  scale_color_brewer(palette = "Set1", guide = "none") +
  labs(x        = "Number of samples",
       y        = "Cumulative mean β",
       title    = "Convergence of propagated estimates",
       subtitle = "Dashed line = final mean across all samples") +
  theme_minimal(base_size = 12)

print(fig5c)
ggsave("results/figures/fig5c_convergence.png",
       fig5c, width = 9, height = 5, dpi = 300)

###############################################################################
# 14.5 TABLA: Comparación numérica punto vs propagado
###############################################################################

tabla_comparacion_prop <- bind_rows(
  coefs_punto |>
    select(variable, media, q025, q975, tipo),
  coefs_prop_plot |>
    select(variable, media, q025, q975, tipo)
) |>
  mutate(
    IC_width      = round(q975 - q025, 4),
    significativo = ifelse(q025 > 0 | q975 < 0, "Yes", "No"),
    variable      = etiquetas_modelo[as.character(variable)]
  ) |>
  arrange(variable, tipo) |>
  select(tipo, variable, media, q025, q975, IC_width, significativo)

cat("\n── Comparación: punto vs propagado ──\n")
print(knitr::kable(tabla_comparacion_prop, digits = 4,
                   col.names = c("Method", "Variable", "β mean",
                                 "2.5%", "97.5%", "IC width", "Significant")))

write.csv(tabla_comparacion_prop,
          "results/tabla_propagacion_vs_punto.csv", row.names = FALSE)
