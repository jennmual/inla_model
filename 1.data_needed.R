###############################################################################
# 0. Clean environment and load data
###############################################################################

rm(list = ls())

library(dplyr)
library(sf)
library(spdep)
library(INLA)
library(stringi)

# Load processed dataset from Script 0
datos_completos <- readRDS("data/datos_completos.rds")


###############################################################################
# 1. Build lagged climate variables and response variables
###############################################################################

# Create lagged climate variables using (t + t-1)/2
# Apply log transformation to FOI variables
datos_modelo <- datos_completos |>
  arrange(cod_mun, year) |>
  group_by(cod_mun) |>
  mutate(
    temp_lag = (mean_2m_temperature_annual_mean +
                  lag(mean_2m_temperature_annual_mean)) / 2,
    
    prec_lag = (acum_total_precipitation_precip_total_anual +
                  lag(acum_total_precipitation_precip_total_anual)) / 2,
    
    hum_lag  = (mean_relative_humidity_annual_mean +
                  lag(mean_relative_humidity_annual_mean)) / 2,
    
    log_FOI_mtv = log(FOI_mtv),
    log_FOI_kte = log(FOI_kte)
  ) |>
  ungroup() |>
  filter(!is.na(temp_lag))  # Remove first year due to lag construction


###############################################################################
# 2. Standardize covariates (global scaling)
###############################################################################

# Function for global standardization (mean = 0, sd = 1)
scale_global <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Apply scaling to climate covariates
datos_modelo <- datos_modelo |>
  mutate(across(c(temp_lag, prec_lag, hum_lag), scale_global))

# (climate variables already scaled above; add geographic/demographic here)
datos_modelo <- datos_modelo |>
  mutate(
    elevation_z  = scale_global(elevation),    # ajusta nombre si difiere
    pop_dens_z   = scale_global(densidad)   # ajusta nombre si difiere
  )


###############################################################################
# 3. Load and prepare spatial data (municipal shapefile)
###############################################################################

mun_sf <- st_read(
  "data/gadm41_COL.gpkg",
  layer = "ADM_ADM_2"
)

# Normalize municipality and department names (uppercase)
mun_sf <- mun_sf |>
  mutate(
    nombre_mun = toupper(NAME_2),
    dpto       = toupper(NAME_1)
  )

datos_modelo <- datos_modelo |>
  mutate(
    nombre_mun = toupper(mun_name),
    dpto       = toupper(dpto)
  )


###############################################################################
# 4. Normalize and harmonize municipality names
###############################################################################

# Step 1: Normalize strings (remove accents, trim spaces)
normalize_name <- function(x) {
  x |>
    toupper() |>
    stringi::stri_trans_general("Latin-ASCII") |>
    trimws()
}

mun_sf <- mun_sf |>
  mutate(
    nombre_mun = normalize_name(nombre_mun),
    dpto       = normalize_name(dpto)
  )

datos_modelo <- datos_modelo |>
  mutate(
    nombre_mun = normalize_name(nombre_mun),
    dpto       = normalize_name(dpto)
  )

# Step 2: Apply manual corrections AFTER normalization
correcciones <- c(
  "CARTAGENA"         = "CARTAGENA DE INDIAS",
  "CUCUTA"            = "SAN JOSE DE CUCUTA",
  "MARIQUITA"         = "SAN SEBASTIAN DE MARIQUITA",
  "MOCOA"             = "SAN MIGUEL DE MOCOA",
  "S.V. DEL CAGUAN"   = "SAN VICENTE DEL CAGUAN",
  "S.V.CH."           = "SAN VICENTE DE CHUCURI",
  "VISTAHERMOSA"      = "VISTA HERMOSA",
  "LEBRIJA"           = "LEBRIJA",
  "EL ZULIA"          = "EL ZULIA",
  "VALLE DEL GUAMUEZ" = "VALLE DEL GUAMUEZ",
  "CARMEN DE BOLIVAR" = "EL CARMEN DE BOLIVAR",
  "CALI"              = "SANTIAGO DE CALI",
  "BUGA"              = "GUADALAJARA DE BUGA",
  "GUAVIARE"          = "SAN JOSE DEL GUAVIARE",
  "SANTIAGO DE TOLU"  = "TOLU"
)

datos_modelo <- datos_modelo |>
  mutate(
    nombre_mun = ifelse(
      nombre_mun %in% names(correcciones),
      correcciones[nombre_mun],
      nombre_mun
    )
  )


###############################################################################
# 5. Create unique municipality identifier (name + department)
###############################################################################

mun_sf <- mun_sf |>
  mutate(nombre_full = paste0(nombre_mun, "_", dpto))

datos_modelo <- datos_modelo |>
  mutate(nombre_full = paste0(nombre_mun, "_", dpto))


###############################################################################
# 6. Validate matching between datasets
###############################################################################

# Check for unmatched municipalities
sin_match <- setdiff(
  unique(datos_modelo$nombre_full),
  unique(mun_sf$nombre_full)
)

if (length(sin_match) > 0) {
  cat("Unmatched municipalities:\n")
  print(sin_match)
  stop("Fix name inconsistencies before proceeding.")
}


###############################################################################
# 7. Match spatial data with model dataset
###############################################################################

# Extract unique municipality IDs
mun_ids <- datos_modelo |>
  select(cod_mun, nombre_full) |>
  distinct()

# Join shapefile with selected municipalities
mun_estudio <- mun_sf |>
  inner_join(mun_ids, by = "nombre_full")

# Ensure consistent ordering
mun_estudio   <- mun_estudio |> arrange(cod_mun)
datos_modelo  <- datos_modelo |> arrange(cod_mun, year)

# Diagnostic check
cat("Number of municipalities:", nrow(mun_estudio), "\n")


###############################################################################
# 8. Validate geometries
###############################################################################

mun_estudio <- st_make_valid(mun_estudio)


###############################################################################
# 9. Build spatial adjacency structure (distance-based)
###############################################################################

# Compute centroids (used for distance-based neighbors)
coords <- st_centroid(mun_estudio) |> st_coordinates()

# Construct k-nearest neighbors graph (robust to disconnected polygons)
nb <- knearneigh(coords, k = 4) |> knn2nb()

# Convert to INLA graph format
nb2INLA("municipios.adj", nb)
g <- inla.read.graph("municipios.adj")

# Inspect network structure
summary(nb)

# Optional: visualize adjacency network (diagnostic plot)
plot(st_geometry(mun_estudio), border = "gray")
plot(nb, coords, add = TRUE, col = "red", lwd = 0.5)


###############################################################################
# Load department-level shapefile (level 1)
###############################################################################

dpto_sf <- st_read(
  "data/gadm41_COL.gpkg",
  layer = "ADM_ADM_1"
)

# Ensure same CRS
dpto_sf <- st_transform(dpto_sf, st_crs(mun_estudio))

###############################################################################
# Plot: departments (background) + municipalities + adjacency network
###############################################################################

plot(st_geometry(dpto_sf),
     col = "lightyellow",
     border = "gray70",
     lwd = 0.8,
     main = "Epidemiological proximity-98 municipalities")

plot(st_geometry(mun_estudio),
     add = TRUE,
     border = "gray40",
     lwd = 0.5)

plot(nb, coords,
     add = TRUE,
     col = "red",
     lwd = 0.7)


###############################################################################
# 10. Save processed objects
###############################################################################

saveRDS(mun_estudio, "data/mun_estudio.rds")
saveRDS(datos_modelo, "data/datos_modelo.rds")
saveRDS(g, "data/grafo_inla.rds")
