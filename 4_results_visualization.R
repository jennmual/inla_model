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

dir.create("scripts/INLA_models/results", recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1. Load all results
###############################################################################

datos_modelo  <- readRDS("scripts/INLA_models/data/datos_modelo_indexed.rds")
mun_estudio   <- readRDS("scripts/INLA_models/data/mun_estudio.rds")

modelo_A      <- readRDS("results/modelo_A.rds")
modelo_B      <- readRDS("results/modelo_B.rds")
modelo_A_kte  <- readRDS("results/modelo_A_kte.rds")
modelo_B_kte  <- readRDS("results/modelo_B_kte.rds")
modelos_uni   <- readRDS("results/modelos_univariables.rds")
comparacion   <- readRDS("results/comparacion_modelos.rds")

# Optional: load propagated coefficients if script 3 was run
coefs_propagados <- tryCatch(
  readRDS("results/coefs_propagados.rds"),
  error = function(e) { message("Propagated coefs not found — skipping"); NULL }
)

municipios_unicos <- sort(unique(datos_modelo$cod_mun))
vars_clima        <- c("temp_lag", "prec_lag", "hum_lag")

###############################################################################
# 2. Helper function: extract fixed effects
###############################################################################

extraer_coefs <- function(modelo, tipo) {
  modelo$summary.fixed |>
    rownames_to_column("variable") |>
    filter(variable != "(Intercept)") |>
    select(variable, media = mean,
           q025 = `0.025quant`, q975 = `0.975quant`) |>
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
  scale_y_discrete(labels = c(
    "temp_lag" = "Temperatura",
    "prec_lag" = "Precipitación",
    "hum_lag"  = "Humedad relativa"
  )) +
  labs(x = "Cambio en log(FOI)", y = "", color = "", shape = "",
       title = "Efecto de variables climáticas sobre log(FOI)",
       subtitle = "FOI tiempo-variante — análisis univariable vs multivariable") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(fig1)
ggsave("results/figures/fig1_forest_plot_mtv.png",
       fig1, width = 7, height = 4, dpi = 300)

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
  scale_y_discrete(labels = c(
    "temp_lag" = "Temperatura",
    "prec_lag" = "Precipitación",
    "hum_lag"  = "Humedad relativa"
  )) +
  labs(x = "Cambio en log(FOI)", y = "", color = "Tipo de FOI",
       title = "Comparación: FOI variante vs FOI constante",
       subtitle = "Modelo B multivariable") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(fig2)
ggsave("results/figures/fig2_forest_kte_vs_mtv.png",
       fig2, width = 7, height = 4, dpi = 300)

###############################################################################
# 5. FIGURE 3: Spatial relative risk map (SRR)
###############################################################################

# BYM2 returns 2*n nodes: first n = structured spatial, second n = unstructured
# Use only first n for the structured component (SRR)
SRR_df <- modelo_B$summary.random$idx_espacio |>
  slice(1:length(municipios_unicos)) |>
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
  scale_fill_viridis_c(option = "plasma", name = "SRR",
                       breaks = c(0.5, 1, 1.5, 2),
                       labels = c("0.5", "1.0", "1.5", "2.0")) +
  labs(title = "Riesgo relativo espacial de FOI — Colombia",
       subtitle = "Efecto espacial estructurado BYM2, Modelo B") +
  theme_void(base_size = 12) +
  theme(legend.position = "right")

print(fig3)
ggsave("results/figures/fig3_mapa_SRR.png",
       fig3, width = 6, height = 8, dpi = 300)

###############################################################################
# 6. FIGURE 4: Observed vs fitted (goodness of fit)
###############################################################################

datos_modelo$predicho <- modelo_B$summary.fitted.values$mean

fig4 <- ggplot(datos_modelo, aes(x = log_FOI_mtv, y = predicho)) +
  geom_point(alpha = 0.35, size = 1.5, color = "#2166ac") +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(x = "log(FOI) observado", y = "log(FOI) predicho",
       title = "Ajuste del Modelo B",
       subtitle = "Línea roja = ajuste perfecto") +
  theme_minimal(base_size = 13)

print(fig4)
ggsave("results/figures/fig4_observed_vs_fitted.png",
       fig4, width = 6, height = 5, dpi = 300)

###############################################################################
# 7. TABLE 1: Model comparison (DIC / WAIC)
###############################################################################

tabla_modelos <- tibble(
  Modelo       = c("A — sin covariables",
                   "B — multivariable",
                   "A — FOI constante",
                   "B — FOI constante"),
  DIC          = c(modelo_A$dic$dic,    modelo_B$dic$dic,
                   modelo_A_kte$dic$dic, modelo_B_kte$dic$dic),
  WAIC         = c(modelo_A$waic$waic,  modelo_B$waic$waic,
                   modelo_A_kte$waic$waic, modelo_B_kte$waic$waic),
  p_efectiva   = c(modelo_A$dic$p.eff,  modelo_B$dic$p.eff,
                   modelo_A_kte$dic$p.eff, modelo_B_kte$dic$p.eff)
) |>
  arrange(WAIC)

print(knitr::kable(tabla_modelos, digits = 2,
                   col.names = c("Modelo", "DIC", "WAIC", "p efectiva")))

saveRDS(tabla_modelos, "results/tabla_comparacion_final.rds")

###############################################################################
# 8. OPTIONAL: Forest plot with propagated uncertainty
###############################################################################

if (!is.null(coefs_propagados)) {
  
  fig5 <- ggplot(coefs_propagados,
                 aes(x = media, y = variable)) +
    geom_point(size = 3, color = "#d73027") +
    geom_errorbarh(aes(xmin = q025, xmax = q975),
                   height = 0.2, color = "#d73027") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_y_discrete(labels = c(
      "temp_lag" = "Temperatura",
      "prec_lag" = "Precipitación",
      "hum_lag"  = "Humedad relativa"
    )) +
    labs(x = "Cambio en log(FOI)", y = "",
         title = "Efecto climático con propagación de incertidumbre",
         subtitle = paste0("Basado en ", nrow(coefs_propagados),
                           " muestras posteriores de la FOI")) +
    theme_minimal(base_size = 13)
  
  print(fig5)
  ggsave("results/figures/fig5_propagated_uncertainty.png",
         fig5, width = 7, height = 4, dpi = 300)
}