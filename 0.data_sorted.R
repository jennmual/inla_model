###############################################################################
# 0. Clean environment
###############################################################################

rm(list = ls())


###############################################################################
# 1. Load required data sources
###############################################################################

# Load municipality metadata (ordering, coordinates, department)
source("scripts/tesis_dani/municipalities_order.R")

# Load FOI estimates
FOI_kte <- readRDS("data/data_tesis_dani/lambda_estimates.RDS") |>
  dplyr::filter(year %in% 2007:2023)

FOI_mtv <- readRDS("data/data_tesis_dani/lambda_estimates_mtv.RDS") |>
  dplyr::filter(year %in% 2007:2023)

# Load climate variables
var_clima_mundani <- readRDS("data/data_2007a2023_mundani/var_clima_mundani.rds")


###############################################################################
# 2. Load libraries
###############################################################################

library(dplyr)
library(tidyr)
library(sf)
library(spdep)
library(INLA)


###############################################################################
# 3. Clean FOI datasets
###############################################################################

# FOI constant (kte): one value per municipality-year
FOI_kte_clean <- FOI_kte |>
  select(
    location,
    year,
    FOI_kte       = median,
    FOI_kte_lower = lower,
    FOI_kte_upper = upper
  ) |>
  rename(cod_mun = location)

# FOI time-varying (mtv): one value per municipality-year
FOI_mtv_clean <- FOI_mtv |>
  select(
    location,
    year,
    FOI_mtv       = median,
    FOI_mtv_lower = lower,
    FOI_mtv_upper = upper
  ) |>
  rename(cod_mun = location)


###############################################################################
# 4. Merge FOI datasets
###############################################################################

# Merge constant and time-varying FOI
datos_base <- FOI_mtv_clean |>
  left_join(FOI_kte_clean, by = c("cod_mun", "year"))

# Check missing FOI_kte values
datos_base |>
  filter(is.na(FOI_kte)) |>
  summarise(
    n_missing = n(),
    n_mun     = n_distinct(cod_mun)
  )


###############################################################################
# 5. Prepare climate data
###############################################################################

# Rename and ensure consistent variable types
var_clima_clean <- var_clima_mundani |>
  rename(year = anio) |>
  mutate(cod_mun = as.character(cod_mun))

# Ensure same type in all datasets
datos_base <- datos_base |> mutate(cod_mun = as.character(cod_mun))
data_order <- data_order |> mutate(cod_mun = as.character(cod_mun))


###############################################################################
# 6. Merge all datasets
###############################################################################

datos_completos <- datos_base |>
  # Add climate variables
  left_join(var_clima_clean, by = c("cod_mun", "year")) |>
  
  # Add geographic information
  left_join(
    data_order |>
      select(cod_mun, dpto, latitud, longitud),
    by = "cod_mun"
  )


###############################################################################
# 7. Data validation and diagnostics
###############################################################################

datos_completos |>
  summarise(
    n_municipalities = n_distinct(cod_mun),
    missing_temp     = sum(is.na(mean_2m_temperature_annual_mean)),
    missing_prec     = sum(is.na(acum_total_precipitation_precip_total_anual)),
    missing_hum      = sum(is.na(mean_relative_humidity_annual_mean)),
    missing_geo      = sum(is.na(latitud))
  )


###############################################################################
# 8. Save final dataset
###############################################################################

saveRDS(
  datos_completos,
  "scripts/INLA_models/data/datos_completos.rds"
)
