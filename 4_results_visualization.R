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

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1. Load all results
###############################################################################

datos_modelo <- readRDS("data/datos_modelo_indexed.rds")
mun_estudio  <- readRDS("data/mun_estudio.rds")

modelo_A     <- readRDS("results/modelo_A.rds")
modelo_B     <- readRDS("results/modelo_B.rds")
modelo_A_kte <- readRDS("results/modelo_A_kte.rds")
modelo_B_kte <- readRDS("results/modelo_B_kte.rds")
modelos_uni  <- readRDS("results/modelos_univariables.rds")
comparacion  <- readRDS("results/comparacion_modelos.rds")

coefs_propagados <- tryCatch(
  readRDS("results/coefs_propagados.rds"),
  error = function(e) { message("Propagated coefs not found — skipping"); NULL }
)

municipios_unicos <- sort(unique(datos_modelo$cod_mun))
vars_clima        <- c("temp_lag", "prec_lag", "hum_lag")
etiquetas_clima   <- c(
  "temp_lag" = "Temperature",
  "prec_lag" = "Precipitation",
  "hum_lag"  = "Relative humidity"
)

###############################################################################
# 2. Helper: extract fixed effects from INLA model
###############################################################################

extraer_coefs <- function(modelo, tipo) {
  modelo$summary.fixed |>
    rownames_to_column("variable") |>
    filter(variable != "(Intercept)") |>
    select(variable,
           media = mean,
           q025  = `0.025quant`,
           q975  = `0.975quant`) |>
    mutate(tipo = tipo)
}

###############################################################################
# 3. FIGURE 1: Forest plot — univariable vs multivariable (FOI variante)
###############################################################################

coefs_forest <- bind_rows(
  map2(modelos_uni, vars_clima, ~ extraer_coefs(.x, "Univariable")),
  extraer_coefs(modelo_B, "Multivariable")
)

fig1 <- ggplot(coefs_forest,
               aes(x = media, y = variable, color = tipo, shape = tipo)) +
  geom_point(position = position_dodge(0.5), size = 3) +
  geom_errorbarh(aes(xmin = q025, xmax = q975),
                 position = position_dodge(0.5), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Univariable"   = "#2166ac",
                                "Multivariable" = "#1a9641")) +
  scale_y_discrete(labels = etiquetas_clima) +
  labs(x     = "Change in log(FOI)",
       y     = "",
       color = "",
       shape = "",
       title    = "Effect of climatic variables on log(FOI)",
       subtitle = "FOI time-variant — univariate vs multivariate") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(fig1)
ggsave("results/figures/fig1_forest_plot_mtv.png",
       fig1, width = 7, height = 4, dpi = 300)


#-------------------------------------------------------------------------------------------
###############################################################################
# 3. FIGURE 1 CORREGIDO: Forest plot — univariable vs multivariable (FOI variante) 
###############################################################################
library(tidyverse)
library(INLA)
library(tibble)

###############################################################################
# 1. Load objects
###############################################################################

datos_modelo <- readRDS("data/datos_modelo_indexed.rds")
mun_estudio  <- readRDS("data/mun_estudio.rds")

modelo_A     <- readRDS("results/modelo_A.rds")
modelo_B     <- readRDS("results/modelo_B.rds")
modelo_A_kte <- readRDS("results/modelo_A_kte.rds")
modelo_B_kte <- readRDS("results/modelo_B_kte.rds")
modelos_uni  <- readRDS("results/modelos_univariables.rds")
comparacion  <- readRDS("results/comparacion_modelos.rds")

coefs_propagados <- tryCatch(
  readRDS("results/coefs_propagados.rds"),
  error = function(e) NULL
)

vars_clima <- c("temp_lag", "prec_lag", "hum_lag")

etiquetas_clima <- c(
  temp_lag = "Temperature",
  prec_lag = "Precipitation",
  hum_lag  = "Relative humidity"
)

###############################################################################
# 2. Helper
###############################################################################

extraer_coef_var <- function(modelo, variable, tipo) {
  
  if (is.null(modelo$summary.fixed)) {
    return(NULL)
  }
  
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

###############################################################################
# 3. Extraer coeficientes univariables
###############################################################################

coefs_uni <- map_dfr(
  vars_clima,
  ~ extraer_coef_var(modelos_uni[[.x]], .x, "Univariable")
)

###############################################################################
# 4. Extraer coeficientes multivariables - FOI constante
###############################################################################

coefs_multi_kte <- map_dfr(
  vars_clima,
  ~ extraer_coef_var(modelo_B_kte, .x, "Multivariable")
)

###############################################################################
# 5. Combinar para FOI constante
###############################################################################

coefs_forest_kte <- bind_rows(coefs_uni, coefs_multi_kte) |>
  mutate(
    variable = factor(variable, levels = rev(vars_clima))
  )

print(coefs_forest_kte)

###############################################################################
# 6. Forest plot - FOI constante
###############################################################################

fig1_kte <- ggplot(
  coefs_forest,
  aes(x = media, y = variable, color = tipo, shape = tipo)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 3
  ) +
  geom_errorbarh(
    aes(xmin = q025, xmax = q975),
    position = position_dodge(width = 0.6),
    height = 0.2
  ) +
  scale_color_manual(
    values = c(
      Univariable   = "#2166ac",
      Multivariable = "#1a9641"
    )
  ) +
  scale_y_discrete(labels = etiquetas_clima) +
  labs(
    x = "Change in log(FOI)",
    y = NULL,
    color = NULL,
    shape = NULL,
    title = "Effect of climatic variables on log(FOI)",
    subtitle = "Univariable Constant FOI vs Multivariable Constant FOI"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom"
  )

print(fig1_kte)

ggsave(
  "results/figures/fig1_forest_plot_kte.png",
  fig1_kte,
  width = 8,
  height = 5,
  dpi = 300
)


###############################################################################
# 4. Extraer coeficientes multivariables para FOI Variante
###############################################################################

coefs_multi <- map_dfr(
  vars_clima,
  ~ extraer_coef_var(modelo_B, .x, "Multivariable")
)


###############################################################################
# 5. Combinar para FOI Variante
###############################################################################

coefs_forest <- bind_rows(coefs_uni, coefs_multi) |>
  mutate(
    variable = factor(variable, levels = rev(vars_clima))
  )

print(coefs_forest)

###############################################################################
# 6. Forest plot - FOI Variante
###############################################################################

fig1 <- ggplot(
  coefs_forest,
  aes(x = media, y = variable, color = tipo, shape = tipo)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 3
  ) +
  geom_errorbarh(
    aes(xmin = q025, xmax = q975),
    position = position_dodge(width = 0.6),
    height = 0.2
  ) +
  scale_color_manual(
    values = c(
      Univariable   = "#2166ac",
      Multivariable = "#1a9641"
    )
  ) +
  scale_y_discrete(labels = etiquetas_clima) +
  labs(
    x = "Change in log(FOI)",
    y = NULL,
    color = NULL,
    shape = NULL,
    title = "Effect of climatic variables on log(FOI)",
    subtitle = "Univariate FOI time varying  vs Multivariate FOI time varying "
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom"
  )

print(fig1)

ggsave(
  "results/figures/fig1_forest_plot_mtv.png",
  fig1,
  width = 8,
  height = 5,
  dpi = 300
)


#-----------------------------------------------------------------------------------

###############################################################################
# 4. FIGURE 2: Forest plot — FOI variante vs FOI constante
###############################################################################

coefs_comparacion <- bind_rows(
  extraer_coefs(modelo_B,     "FOI variante"),
  extraer_coefs(modelo_B_kte, "FOI constante")
)

fig2 <- ggplot(coefs_comparacion,
               aes(x = media, y = variable, color = tipo)) +
  geom_point(position = position_dodge(0.5), size = 3) +
  geom_errorbarh(aes(xmin = q025, xmax = q975),
                 position = position_dodge(0.5), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("FOI variante"  = "#d73027",
                                "FOI constante" = "#4575b4")) +
  scale_y_discrete(labels = etiquetas_clima) +
  labs(x        = "Change in log(FOI)",
       y        = "",
       color    = "Type of FOI",
       title    = "Comparison: Variant FOI vs Constant FOI",
       subtitle = "Multivariable Model B") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(fig2)
ggsave("results/figures/fig2_forest_kte_vs_mtv.png",
       fig2, width = 7, height = 4, dpi = 300)

###############################################################################
# 5. FIGURE 3: Spatial relative risk map (SRR)
# NOTE: model uses BESAG (not BYM2) — only one spatial component
# All n rows of summary.random$idx_espacio are the structured effect
###############################################################################

SRR_df <- modelo_B$summary.random$idx_espacio |>
  mutate(
    cod_mun  = municipios_unicos,
    SRR      = exp(mean),
    SRR_q025 = exp(`0.025quant`),
    SRR_q975 = exp(`0.975quant`)
  )

mun_mapa <- mun_estudio |>
  left_join(SRR_df, by = "cod_mun")

fig3 <- ggplot(mun_mapa) +
  geom_sf(aes(fill = SRR), color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(
    option = "plasma",
    name   = "SRR",
    direction = -1
  ) +
  labs(title    = "Spatial relative risk (SRR) of FOI",
       subtitle = "98 Municipalities of Colombia") +
  theme_void(base_size = 12) +
  theme(legend.position = "right")

print(fig3)
ggsave("results/figures/fig3_mapa_SRR.png",
       fig3, width = 6, height = 8, dpi = 300)


###############################################################################
# 5. FIGURE 3: Spatial relative risk map (SRR)
# WITH DEPARTMENTS FROM COLOMBIA ON BACKGROUND
###############################################################################

library(sf)
library(ggplot2)

dpto_sf <- st_read(
  "data/gadm41_COL.gpkg",
  layer = "ADM_ADM_1"
)

# mismo CRS
dpto_sf <- st_transform(dpto_sf, st_crs(mun_estudio))

fig3 <- ggplot() +

  # Fondo: departamentos
  geom_sf(
    data = dpto_sf,
    fill = "grey95",
    color = "grey60",
    linewidth = 0.12
  ) +

  # Municipios con SRR
  geom_sf(
    data = mun_mapa,
    aes(fill = SRR),
    color = "grey60",
    linewidth = 0.12
  ) +

  scale_fill_viridis_c(
    option    = "plasma",
    direction = -1,
    name      = "SRR"
  ) +

  labs(
    title    = "Spatial relative risk (SRR) of FOI",
    subtitle = "98 Municipalities of Colombia"
  ) +

  theme_void(base_size = 12) +
  theme(
    legend.position = "right"
  )

print(fig3)

ggsave(
  "results/figures/fig3_mapa_SRR.png",
  fig3,
  width = 6,
  height = 8,
  dpi = 300
)






###############################################################################
# 6. FIGURE 4: Observed vs fitted
###############################################################################

datos_modelo$predicho <- modelo_B$summary.fitted.values$mean

fig4 <- ggplot(datos_modelo, aes(x = log_FOI_mtv, y = predicho)) +
  geom_point(alpha = 0.35, size = 1.5, color = "#2166ac") +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(x        = "log(FOI) observed",
       y        = "predicted log(FOI)",
       title    = "Model B Adjustment") +
  theme_minimal(base_size = 13)

print(fig4)
ggsave("results/figures/fig4_observed_vs_fitted.png",
       fig4, width = 6, height = 5, dpi = 300)

###############################################################################
# 7. TABLE 1: Model comparison
# NOTE: FOI constante models not comparable in scale — reported separately
###############################################################################

tabla_mtv <- tibble(
  Modelo     = c("A — sin covariables (FOI variante)",
                 "B — multivariable (FOI variante)",
                 paste("B —", names(modelos_uni), "(univariable)")),
  DIC        = c(modelo_A$dic$dic,
                 modelo_B$dic$dic,
                 sapply(modelos_uni, function(m) m$dic$dic)),
  WAIC       = c(modelo_A$waic$waic,
                 modelo_B$waic$waic,
                 sapply(modelos_uni, function(m) m$waic$waic)),
  p_efectiva = c(modelo_A$dic$p.eff,
                 modelo_B$dic$p.eff,
                 sapply(modelos_uni, function(m) m$dic$p.eff))
) |>
  arrange(WAIC)

cat("\n--- Comparación modelos FOI variante ---\n")
print(knitr::kable(tabla_mtv, digits = 2,
                   col.names = c("Modelo", "DIC", "WAIC", "p efectiva")))

# FOI constante — reported separately (different response scale)
tabla_kte <- tibble(
  Modelo     = c("A — sin covariables (FOI constante)",
                 "B — multivariable (FOI constante)"),
  DIC        = c(modelo_A_kte$dic$dic,  modelo_B_kte$dic$dic),
  WAIC       = c(modelo_A_kte$waic$waic, modelo_B_kte$waic$waic),
  p_efectiva = c(modelo_A_kte$dic$p.eff, modelo_B_kte$dic$p.eff)
) |>
  arrange(WAIC)

cat("\n--- Comparación modelos FOI constante (escala diferente) ---\n")
print(knitr::kable(tabla_kte, digits = 2,
                   col.names = c("Modelo", "DIC", "WAIC", "p efectiva")))

saveRDS(tabla_mtv, "results/tabla_comparacion_mtv.rds")
saveRDS(tabla_kte, "results/tabla_comparacion_kte.rds")

###############################################################################
# 8. OPTIONAL: Forest plot with propagated uncertainty
###############################################################################

if (!is.null(coefs_propagados)) {
  
  fig5 <- ggplot(coefs_propagados, aes(x = media, y = variable)) +
    geom_point(size = 3, color = "#d73027") +
    geom_errorbarh(aes(xmin = q025, xmax = q975),
                   height = 0.2, color = "#d73027") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_y_discrete(labels = etiquetas_clima) +
    labs(x        = "Cambio en log(FOI)",
         y        = "",
         title    = "Efecto climático con propagación de incertidumbre",
         subtitle = paste0(length(coefs_propagados),
                           " muestras posteriores de la FOI")) +
    theme_minimal(base_size = 13)
  
  print(fig5)
  ggsave("results/figures/fig5_propagated_uncertainty.png",
         fig5, width = 7, height = 4, dpi = 300)
}
