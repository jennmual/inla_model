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
# 7. Load and prepare elevation data
###############################################################################

# corregir elevation
elevation <- read_excel("data/elevation.xlsx") |>
  mutate(cod_mun = str_pad(cod_mun, width = 5, pad = "0")) |>
  select(cod_mun, elevation)

# corregir datos_completos y rehacer join
datos_completos <- datos_completos |>
  mutate(cod_mun = str_pad(cod_mun, width = 5, pad = "0")) |>
  select(-any_of("elevation")) |>
  left_join(elevation, by = "cod_mun")

datos_completos |> #To verified
  summarise(
    n_municipios = n_distinct(cod_mun),
    missing_elev = sum(is.na(elevation))
  )



###############################################################################
# 8. Load and prepare population density data
###############################################################################

densidad_mun <- readRDS("data/densidad_mun2018_2035.rds") |>
  filter(!is.na(cod_mun)) |>
  mutate(
    cod_mun = str_pad(cod_mun, width = 5, pad = "0"),
    year = as.numeric(Año),
    densidad = densidad |>
      str_replace_all("\\.", "") |>
      str_replace(",", ".") |>
      as.numeric()
  ) |>
  select(cod_mun, year, densidad)

# extend backwards using 2018 value
densidad_expandida <- densidad_mun |>
  group_by(cod_mun) |>
  complete(year = 2007:2023) |>
  arrange(year) |>
  fill(densidad, .direction = "downup") |>
  ungroup()

datos_completos <- datos_completos |>
  left_join(densidad_expandida, by = c("cod_mun", "year"))


###############################################################################
# 9. Data validation and diagnostics
###############################################################################

datos_completos |>
  summarise(
    n_municipalities = n_distinct(cod_mun),
    missing_temp     = sum(is.na(mean_2m_temperature_annual_mean)),
    missing_prec     = sum(is.na(acum_total_precipitation_precip_total_anual)),
    missing_hum      = sum(is.na(mean_relative_humidity_annual_mean)),
    missing_geo      = sum(is.na(latitud)),
    missing_altu     = sum(is.na(elevation)),
    missing_geo      = sum(is.na(densidad))
  )

###############################################################################
# 10. Save final dataset
###############################################################################

saveRDS(
  datos_completos,
  "data/datos_completos.rds"
)
