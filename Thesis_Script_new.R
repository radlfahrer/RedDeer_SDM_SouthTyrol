# Habitat Suitability Indices and Landscape Metrics in species distribution models: 
# A case study for Red Deer (Cervus elaphus) in the Central Alps

# Author: Timo Stidl

setwd('C:/Users/timo/Documents/Master/UIBK/Master_Thesis/R/data')

# packages
library(terra)
library(sf)
library(tibble)
library(ggplot2)
library(dplyr)
library(randomForest)
library(landscapemetrics)
library(geodata)
library(factoextra)
library(RStoolbox)
library(tidyverse)
library(USE)
library(modEvA)
library(fastfocal)
library(blockCV)
library(quantreg)
library(patchwork)
library(rfPermute)
library(pdp)
library(usdm)

# Predictor Variables ####

## Land use / Land cover ####
# Source: Rüdisser et al. 2015

LULC <- rast('predictor_variables/LULC/landcover.tif')

study_extent <- ext(LULC)
study_CRS <- crs(LULC)

### LULC reclassification and interpolation ####

# Reclassification
m_LULC <- c(1110, 1121, 1, # Artificial surfaces (without sparse settlement areas)
            1130, 1320, 1,
            1330, 1410, 1,
            1420, 1430, 1,
            1510, 1530, 1,
            
            556, 557, 2, # Sparse settlement areas
            1121, 1130, 2, 
            1320, 1330, 2,
            1410, 1420, 2,
            1430, 1510, 2,
            1530, 1610, 2,
            1620, 1630, 2,
            
            1, 3, 3, # Arable land
            
            409, 410, 4, # Permanent crops
            650, 651, 4, 
            654, 655, 4,
            
            8, 10, 5, # Intensively used grasslands
            17, 19, 5, 
            665, 666, 5,
            
            10, 13, 6, # Extensively used grasslands
            637, 638, 6, 
            658, 659, 6,
            
            653, 654, 7, # Natural grasslands
            769, 770, 7,
            
            3000, 3003, 8, # Coniferous forest
            3004, 3005, 8,
            3011, 3013, 8,
            
            3005, 3008, 9, # Deciduous forest
            3010, 3011, 9,
            3013, 3014, 9,
            
            3008, 3010, 10, # Mixed forest
            3159, 3160, 10,
            
            4099, 4300, 11, # Sparsely vegetated alpine land
            
            5209, 5220, 12, # Water bodies
            
            5099, 5100, 13, # Wetland
            
            3003, 3004, 14, # Krummholz
            3139, 3140, 14,
            
            779, 780, NA, # Not classified
            1619, 1620, NA
)

rclmat_LULC <- matrix(m_LULC, ncol = 3, byrow = TRUE)

LULC_classified <- LULC %>% 
  ifel(is.na(.), 0, .) %>% 
  classify(rclmat_LULC, right = TRUE)

writeRaster(LULC_classified, 'predictor_variables/LULC/LULC_classified.tif', overwrite = TRUE)
LULC_classified <- rast('predictor_variables/LULC/LULC_classified.tif')

# Interpolation

custom_modal <- function(x, ...) {
  x_filtered <- x[x != 0]  # exclude zeros
  if (length(x_filtered) == 0) {
    return(0)
  }
  terra::modal(x_filtered, ...)  # compute modal without zeros
}

interpolation <- function(input_raster) {
  
  # Define 3x3 window for neighbors
  window <- matrix(1, 3, 3)
  
  # Initialize raster to track filling
  r_filled <- input_raster
  
  # Loop until no NA values remain
  while(any(is.na(values(r_filled)))) {
    r_prev <- r_filled
    r_filled <- focal(r_filled, w=window, fun=custom_modal, ties= 'random', na.policy="only", na.rm=TRUE)
  }
  r_filled[r_filled == 0] <- NA
  return(r_filled)
}

path <- 'predictor_variables/LULC/LULC_interpolated.tif'

if (!exists("LULC_interpolated")) {
  if (file.exists(path)) {
    LULC_interpolated <- rast(path)
  } else {
    LULC_interpolated <- interpolation(LULC_classified)
    writeRaster(LULC_interpolated, path)
  }
}


## Habitat Suitability Index ####

path <- 'predictor_variables/LULC/'

HSI_values <- c(0, # Artificial surfaces
                10, # Sparse settlement areas
                40, # Arable land
                25, # Permanent crops
                75, # Intensively used grasslands
                75, # Extensively used grasslands
                75, # Natural grasslands
                95, # Coniferous forest
                95, # Deciduous forest
                100, # Mixed forest
                0, # Sparsely vegetated alpine land
                10, # Water bodies
                70, # Wetland
                30) # Krummholz


m_HSI_values <- as.vector(rbind(c(1:length(HSI_values)), HSI_values))

rclmat_HSI_values <- matrix(m_HSI_values, ncol = 2, byrow = TRUE)

HSI_classified <- classify(LULC_interpolated, rclmat_HSI_values) # on point basis

writeRaster(HSI_classified, 'HSI_classified.tif')
HSI_classified <- rast('HSI_classified.tif')

# HSI for different radii

radii <- c(100, 250, 500)
HSI <- vector("list", length(radii))

for (i in seq_along(radii)) {
  window <- focalMat(HSI_classified, radii[i], type='circle')
  HSI[[i]] <- focal(HSI_classified, w = window, fun = mean, na.rm = TRUE)
  names(HSI)[i] <- paste0("HSI_", radii[i])
}

# Save files named by the vector values
for (i in seq_along(radii)) {
  fname <- paste0(path, "HSI_", radii[i], ".tif")
  writeRaster(HSI[[i]], fname)
}

terraOptions(memfrac = 0.9)

window <- focalMat(HSI_classified, 250, type='circle')
HSI_250 <- focal(
  HSI_classified, w = window,
  fun = mean,
  na.rm = TRUE,
  filename = paste0(path, "HSI_250", ".tif"),
  overwrite = TRUE
)

# Tile approach because R crashes with full raster

r <- HSI_classified

# Define grid for 8 tiles with 1km overlap
n_rows <- 2
n_cols <- 4
overlap_m <- 1000  # 1 km

# Get raster extent and resolution
e <- ext(r)
res_m <- res(r)[1]
overlap_cells <- ceiling(overlap_m / res_m)

# Calculate tile extents with overlap
tile_exts <- list()
tile_idx <- expand.grid(row = 1:n_rows, col = 1:n_cols)

for (i in 1:nrow(tile_idx)) {
  row_i <- tile_idx$row[i]
  col_i <- tile_idx$col[i]
  
  # Tile extent with overlap (shrink by overlap on all sides)
  xmin <- e$xmin + (col_i - 1) * ((e$xmax - e$xmin) / n_cols) - overlap_m
  xmax <- e$xmin + (col_i    ) * ((e$xmax - e$xmin) / n_cols) + overlap_m
  ymin <- e$ymin + (row_i - 1) * ((e$ymax - e$ymin) / n_rows) - overlap_m
  ymax <- e$ymin + (row_i    ) * ((e$ymax - e$ymin) / n_rows) + overlap_m
  
  tile_exts[[i]] <- ext(xmin, xmax, ymin, ymax)
}

# Crop tiles and process
out_files <- paste0("tile500_", 1:8, ".tif")
out_tiles <- list()

for (i in 7:8) {
  
  # Crop tile with overlap
  r_tile <- crop(r, tile_exts[[i]])
  
  # Run function on tile
  window <- focalMat(r_tile, 500, type='circle')
  r_tile_out <- focal(
    r_tile, w = window,
    fun = mean,
    na.rm = TRUE)
  
  # Write to file
  writeRaster(r_tile_out, out_files[i], overwrite = TRUE)
  out_tiles[[i]] <- r_tile_out
}


# Mosaic tiles back together

HSI_500 <- do.call(
  mosaic,
  c(out_tiles, 
    list(fun = "mean"))
)

names(HSI_500) <- 'HSI_500'
writeRaster(HSI_500, paste0(path, "HSI_500", ".tif"))

## Landscape metrics ####

path <- 'predictor_variables/LULC/'

m_landscape_metrics <- c(0, 2, 1, # Artificial
                         2, 5, 2, # Agriculture
                         7, 10, 3, # Forest
                         11, 12, 4, # Water
                         10, 11, 5, # Non-vegetated open land
                         5, 7, 6, # Non-intensively used grassland
                         13, 14, 7, # Krummholz
                         12, 13, 8) # Wetland


rclmat_LM <- matrix(m_landscape_metrics, ncol = 3, byrow = TRUE)

LM_classified <- classify(LULC_interpolated, rclmat_LM, right = TRUE)

writeRaster(LM_classified, 'predictor_variables/LULC/LM_classified.tif')
LM_classified <- rast('predictor_variables/LULC/LM_classified.tif')

# Aggregate to 25m resolution
LM_classified_25m <- terra::aggregate(LM_classified, fact = 5, fun = modal, na.rm = TRUE)

writeRaster(LM_classified_25m, 'LM_classified_25m.tif')
LM_classified_25m <- rast('LM_classified_25m.tif')

# Forest/No Forest layer for Edge Density
m_forest <- c(0, 2, 0,
              2, 3, 1,
              3, 8, 0)

rclmat_forest <- matrix(m_forest, ncol = 3, byrow = TRUE)

LM_classified_forest_25m <- classify(LM_classified_25m, rclmat_forest, right = TRUE)

### Class-level metrics (Edge Density) ####

edge_density_100 <- window_lsm(LM_classified_forest_25m,
                               window = focalMat(LM_classified_forest_25m, 100, type = 'circle'),
                               what = c('lsm_l_ed'),
                               progress = TRUE)

writeRaster(edge_density_100$layer_1$lsm_l_ed, 'edge_density_100.tif')

edge_density_250 <- window_lsm(LM_classified_forest_25m,
                               window = focalMat(LM_classified_forest_25m, 250, type = 'circle'),
                               what = c('lsm_l_ed'),
                               progress = TRUE)

writeRaster(edge_density_250$layer_1$lsm_l_ed, 'edge_density_250.tif')

edge_density_500 <- window_lsm(LM_classified_forest_25m,
                               window = focalMat(LM_classified_forest_25m, 500, type = 'circle'),
                               what = c('lsm_l_ed'),
                               progress = TRUE)

writeRaster(edge_density_500$layer_1$lsm_l_ed, 'edge_density_500.tif')

### Landscape-level metrics ####

landscape_metrics_100 <- window_lsm(LM_classified_25m,
                                    window = focalMat(LM_classified_25m, 100, type = 'circle'),
                                    what = c('lsm_l_shdi', 'lsm_l_contag', 'lsm_l_mesh', 
                                             'lsm_l_pladj', 'lsm_l_pd'),
                                    progress = TRUE)

shannon_100 <- landscape_metrics_100$layer_1$lsm_l_shdi
writeRaster(shannon_100, 'LM_shannon_100.tif')

contag_100 <- landscape_metrics_100$layer_1$lsm_l_contag
contag_100_filled <- terra::mask(subst(contag_100, NA, 100), LM_classified_25m) %>% 
  writeRaster('LM_contag_100.tif')

mesh_100 <- landscape_metrics_100$layer_1$lsm_l_mesh
writeRaster(mesh_100, 'LM_mesh_100.tif')

patch_density_100 <- landscape_metrics_100$layer_1$lsm_l_pd
writeRaster(patch_density_100, 'LM_pd_100.tif')

pladj_100 <- landscape_metrics_100$layer_1$lsm_l_pladj
writeRaster(pladj_100, 'LM_pladj_100.tif')


# 250m radius
landscape_metrics_250 <- window_lsm(LM_classified_25m,
                                    window = focalMat(LM_classified_25m, 250, type = 'circle'),
                                    what = c('lsm_l_shdi', 'lsm_l_contag', 'lsm_l_mesh', 
                                             'lsm_l_pladj', 'lsm_l_pd'),
                                    progress = TRUE)


shannon_250 <- landscape_metrics_250$layer_1$lsm_l_shdi
writeRaster(shannon_250, 'LM_shannon_250.tif')

contag_250 <- landscape_metrics_250$layer_1$lsm_l_contag
contag_250_filled <- terra::mask(subst(contag_250, NA, 100), LM_classified_25m) %>% 
  writeRaster('LM_contag_250.tif')

mesh_250 <- landscape_metrics_250$layer_1$lsm_l_mesh
writeRaster(mesh_250, 'LM_mesh_250.tif')

patch_density_250 <- landscape_metrics_250$layer_1$lsm_l_pd
writeRaster(patch_density_250, 'LM_pd_250.tif')

pladj_250 <- landscape_metrics_250$layer_1$lsm_l_pladj
writeRaster(pladj_250, 'LM_pladj_250.tif')

# 500m radius
landscape_metrics_500 <- window_lsm(LM_classified_25m,
                                    window = focalMat(LM_classified_25m, 500, type = 'circle'),
                                    what = c('lsm_l_shdi', 'lsm_l_contag', 'lsm_l_mesh', 
                                             'lsm_l_pladj', 'lsm_l_pd'),
                                    progress = TRUE)


shannon_500 <- landscape_metrics_500$layer_1$lsm_l_shdi
writeRaster(shannon_500, 'LM_shannon_500.tif')

contag_500 <- landscape_metrics_500$layer_1$lsm_l_contag
contag_500_filled <- terra::mask(subst(contag_500, NA, 100), LM_classified_25m) %>% 
  writeRaster('LM_contag_500.tif')

mesh_500 <- landscape_metrics_500$layer_1$lsm_l_mesh
writeRaster(mesh_500, 'LM_mesh_500.tif')

patch_density_500 <- landscape_metrics_500$layer_1$lsm_l_pd
writeRaster(patch_density_500, 'LM_pd_500.tif')

pladj_500 <- landscape_metrics_500$layer_1$lsm_l_pladj
writeRaster(pladj_500, 'LM_pladj_500.tif')

## Environmental predictors ####

### Topography ####

#### Elevation ####
# Source: Geoportal South Tyrol

if (!exists('dem.tif')){
  if (!file.exists('predictor_variables/environmental_predictors/dem.tif')) {
  
    m_dem <- c(-9999, NA)
    rclmat <- matrix(m_dem, ncol = 2, byrow = TRUE)
    
    dem_list <- list.files('predictor_variables/environmental_predictors/DEM/', 'tif$', full.names = TRUE)
    dem <- sprc(lapply(dem_list, rast)) %>%
      merge() %>% 
      classify(rclmat) %>% 
      project(study_CRS) %>% 
      crop(study_extent) %>% 
      writeRaster('predictor_variables/environmental_predictors/dem.tif')
    
  } else {
    dem <- rast('predictor_variables/environmental_predictors/dem.tif')
  }
}

#### Slope ####

if (!exists('slope.tif')) {
  if (file.exists('predictor_variables/environmental_predictors/slope.tif')) {
    slope <- rast('predictor_variables/environmental_predictors/slope.tif')
  } else {
    slope <- terrain(dem, 'slope')
    writeRaster(slope, 'predictor_variables/environmental_predictors/slope.tif')
  }
}

#### Northness & Eastness ####

if (!exists('northness.tif')) {
  if (file.exists('predictor_variables/environmental_predictors/northness.tif')) {
    northness <- rast('predictor_variables/environmental_predictors/northness.tif')
    eastness <- rast('predictor_variables/environmental_predictors/eastness.tif')
  } else {
    aspect_rad <- terrain(dem, 'aspect') * (pi / 180)
    northness <- cos(aspect_rad)
    eastness <- sin(aspect_rad)
    writeRaster(northness, 'predictor_variables/environmental_predictors/northness.tif')
    writeRaster(eastness, 'predictor_variables/environmental_predictors/eastness.tif')
  }
}

### Climate ####
# Source: WorldClim 2: Fick & Hijmans 2017

#### Temperature ####

if (!exists('temp_avg_year.tif')) {
  if (file.exists('predictor_variables/environmental_predictors/temp_avg_year.tif')) {
    temp_avg_year <- rast('predictor_variables/environmental_predictors/temp_avg_year.tif')
    temp_min_month <- rast('predictor_variables/environmental_predictors/temp_min_month.tif')
  } else {
    temp <- worldclim_country('IT', 'tavg', path = tempdir(), version = '2.1') %>% 
      project(study_CRS) %>%
      crop(study_extent)
    
    temp_avg_year <- app(temp, mean, na.rm = TRUE)
    writeRaster(temp_avg_year, 'predictor_variables/environmental_predictors/temp_avg_year.tif')
    
    temp_min_month <- app(temp, min, na.rm = TRUE)
    writeRaster(temp_min_month, 'predictor_variables/environmental_predictors/temp_min_month.tif')
  }
}

#### Precipitation ####

if (!exists('precip_winter.tif')) {
  if (file.exists('predictor_variables/environmental_predictors/precip_winter.tif')) {
    precip_winter <- rast('predictor_variables/environmental_predictors/precip_winter.tif')
    precip_spring <- rast('predictor_variables/environmental_predictors/precip_spring.tif')
    precip_summer <- rast('predictor_variables/environmental_predictors/precip_summer.tif')
    precip_fall <- rast('predictor_variables/environmental_predictors/precip_fall.tif')
  } else {
    precip <- worldclim_country('IT', 'prec', path = tempdir(), version = '2.1') %>%
      project(study_CRS) %>%
      crop(study_extent)
    
    precip_winter <- sum(precip[[c(1,2,12)]])/3 
    writeRaster(precip_winter, 'predictor_variables/environmental_predictors/precip_winter.tif')
    
    precip_spring <- sum(precip[[3:5]])/3 
    writeRaster(precip_spring, 'predictor_variables/environmental_predictors/precip_spring.tif')
    
    precip_summer <- sum(precip[[6:8]])/3
    writeRaster(precip_summer, 'predictor_variables/environmental_predictors/precip_summer.tif')
    
    precip_fall <- sum(precip[[9:11]])/3
    writeRaster(precip_fall, 'predictor_variables/environmental_predictors/precip_fall.tif')
  }
}

#### Solar radiation ####

if (!exists('solar_radiation_avg.tif')) {
  if (file.exists('predictor_variables/environmental_predictors/solar_radiation_avg.tif')) {
    solar_radiation <- rast('predictor_variables/environmental_predictors/solar_radiation_avg.tif')
  } else {
    solar_radiation <- worldclim_country('IT', 'srad', path = tempdir(), version = '2.1') %>%
      project(study_CRS) %>%
      crop(study_extent)
    solar_radiation_avg <- app(solar_radiation, mean, na.rm = TRUE)
    writeRaster(solar_radiation_avg, 'predictor_variables/environmental_predictors/solar_radiation_avg.tif')
  }
}

## Distance-based predictors ####
path <- 'predictor_variables/disturbance_indicators/'

### Distance to nature ####
# Source: Rüdisser et al. 2012

D2N <- rast('predictor_variables/disturbance_indicators/D2N.tif')

### Distance to asphalted roads ####
# Source: GeoPortal South Tyrol

transport_network <- st_read('predictor_variables/disturbance_indicators/LinksSequences-Queryable_line.shp')

# Extract asphalted roads
unique_class<- distinct(transport_network, transport_network$CLASS, .keep_all = TRUE)

tn_classes <- data.frame(
  CLASS = integer(),
  CLASS_DE = character()
)

for (i in 1:nrow(unique_class)){
  class_codes = as.numeric(unlist(strsplit(unique_class[i,]$CLASS, ",")))
  class_description = unlist(strsplit(unique_class[i,]$CLASS_DE, ","))
  for (j in length(class_codes)) {
    if (!(class_codes[j] %in% tn_classes$CLASS)){
      tn_classes <- rbind(tn_classes, data.frame(CLASS = class_codes[j], CLASS_DE = class_description[j]))
    }
  }
  
}

# Asphalted roads: classes 0-8
road_classes = c(0:8)

asphalted_roads <- transport_network[sapply(strsplit(transport_network$CLASS, ","), function(x) any(as.integer(x) %in% road_classes)),] %>% 
  st_transform(crs = study_CRS)

if (!file.exists('predictor_variables/disturbance_indicators/distance_to_roads.tif')) {
  d_roads <- asphalted_roads %>% 
    rasterize(LULC) %>% 
    distance(gdal = c('BIGTIFF = YES')) %>% 
    writeRaster('predictor_variables/disturbance_indicators/distance_to_roads.tif')
} else {d_roads <- rast('predictor_variables/disturbance_indicators/distance_to_roads.tif')}


### Distance to watercourse ####
# Source: GeoPortal South Tyrol

lakes <- st_read('predictor_variables/disturbance_indicators/lakes_ST/Lakes_polygon.shp') %>%
  st_transform(crs = 32632)

watercourses <- st_read('predictor_variables/disturbance_indicators/watercourses_ST/Watercourses_line.shp') %>% 
  st_transform(crs= 32632)

lakes_rast <- rasterize(lakes, LULC)
d_lakes <- distance(lakes_rast, gdal = c('BIGTIFF = YES'))

watercourse_rast <- rasterize(watercourses, LULC)
d_watercourse <- distance(watercourse_rast, gdal = c('BIGTIFF = YES'))

d_water <- min(d_lakes, d_watercourse)
writeRaster(d_water, 'predictor_variables/disturbance_indicators/distance_to_water.tif', overwrite = TRUE)

### Distance to core forest ####

# minimum patch size 30ha (Urbina et al. 2023)

forest_patches <- patches(LM_classified_forest, directions = 8, zeroAsNA = TRUE)
freq_tbl <- freq(forest_patches)

# identify patches >= 30ha
valid_ids <- freq_tbl$value[freq_tbl$count >= (30*10000/(res(forest_patches)[1]^2))]

forest_large <- forest_patches
forest_large[!forest_patches %in% valid_ids] <- NA

dist_to_forest <- distance(forest_large)

writeRaster(forest_large, 'forest_patches_30ha.tif', overwrite = TRUE)
writeRaster(dist_to_forest, paste0(path, 'predictor_variables/disturbance_indicators/distance_to_forest.tif'))

### Distance to forest edge ####

d_to_nonforest <- distance(LM_classified_forest, target = 1, exclude = NA)

d_to_forest <- distance(LM_classified_forest, target = 0, exclude = NA)

dist_edge <- (d_to_forest + d_to_nonforest) - 5

writeRaster(dist_edge, paste0(path, 'distance_to_forest_edge.tif'), overwrite = TRUE)

### Distance to settlement areas ####
# No roads and trains

#### Dense development ####

dense_development_codes <- c(1111, 1112, 1121, 1200, 1340, 1410, 1430, 1520, 1530)
dense_development <- ifel(LULC %in% dense_development_codes, 1, NA)

dist_to_dense_development <- distance(dense_development, gdal = c('BIGTIFF = YES'))
writeRaster(dist_to_dense_development, paste0(path, 'distance_to_dense_development.tif'))

#### Sparse development ####

sparse_development <- ifel(LULC_interpolated == 2, 1, NA)

dist_to_sparse_development <- distance(sparse_development)
writeRaster(dist_to_sparse_development, paste0(path, 'distance_to_sparse_development.tif'))


# Principal Component Analysis ####

# Load all predictor variables and resample

## Landscape metrics ####
path <- 'predictor_variables/LULC/'

# 100m radius
files_LM100 <- list.files('predictor_variables/LULC/', pattern = "^LM_.*_100\\.tif$", full.names = TRUE)
layer_names_LM100 <- sub("^LM_(.*)_100\\.tif$", "\\1", basename(files_LM100))

ras_list_100 <- lapply(seq_along(files_LM100), function(i) {
  r <- rast(files_LM100[i])
  names(r) <- layer_names_LM100[i]
  r
})

LM_100 <- do.call(c, ras_list_100)

LM100_PCA <- rasterPCA(LM_100, spca = TRUE, maskCheck = TRUE)
saveRDS(LM100_PCA, 'rf_maps/PCA_LM100.rds')

# Keep PCs with cumulative explained variance >= 95%
explained_var_cum <- cumsum(LM100_PCA$model$sdev^2 / sum(LM100_PCA$model$sdev^2))

n_comp <- which(explained_var_cum >= .95)[1]

for (k in 1:n_comp){
  fname <- paste0('predictor_variables/predictor_variables_resampled/', "LM100_PC", k, ".tif")
  writeRaster(LM100_PCA$map[[k]], fname, overwrite = TRUE)
}

# 250m radius
files_LM250 <- list.files('predictor_variables/LULC/', pattern = "^LM_.*_250\\.tif$", full.names = TRUE)
layer_names_LM250 <- sub("^LM_(.*)_250\\.tif$", "\\1", basename(files_LM250))

ras_list_250 <- lapply(seq_along(files_LM250), function(i) {
  r <- rast(files_LM250[i])
  names(r) <- layer_names_LM250[i]
  r
})

LM_250 <- do.call(c, ras_list_250)

LM250_PCA <- rasterPCA(LM_250, spca = TRUE, maskCheck = TRUE)
saveRDS(LM250_PCA, 'rf_maps/PCA_LM250.rds')

explained_var_cum <- cumsum(LM250_PCA$model$sdev^2 / sum(LM250_PCA$model$sdev^2))

n_comp <- which(explained_var_cum >= .95)[1]

for (k in 1:n_comp){
  fname <- paste0('predictor_variables/predictor_variables_resampled/', "LM250_PC", k, ".tif")
  writeRaster(LM250_PCA$map[[k]], fname, overwrite = TRUE)
}

# 500m radius
files_LM500 <- list.files('predictor_variables/LULC/', pattern = "^LM_.*_500\\.tif$", full.names = TRUE)
layer_names_LM500 <- sub("^LM_(.*)_500\\.tif$", "\\1", basename(files_LM500))

ras_list_500 <- lapply(seq_along(files_LM500), function(i) {
  r <- rast(files_LM500[i])
  names(r) <- layer_names_LM500[i]
  r
})

LM_500 <- do.call(c, ras_list_500)

LM500_PCA <- rasterPCA(LM_500, spca = TRUE, maskCheck = TRUE)
saveRDS(LM500_PCA, 'rf_maps/PCA_LM500.rds')

explained_var_cum <- cumsum(LM500_PCA$model$sdev^2 / sum(LM500_PCA$model$sdev^2))

n_comp <- which(explained_var_cum >= .95)[1]

for (k in 1:n_comp){
  fname <- paste0('predictor_variables/predictor_variables_resampled/', "LM500_PC", k, ".tif")
  writeRaster(LM500_PCA$map[[k]], fname, overwrite = TRUE)
}

## Environmental Predictors ####

# topography (dem, slope, northness eastness), climate (precip per season, temp avg year & min month),
# interference factors (distance to forest, settlement, asphalted roads, water, D2N)

folder_path <- 'predictor_variables/environmental_predictors'
EP_list <- list.files(path = folder_path, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)

# Resample
EP_list_resampled <- lapply(EP_list, function(f){
  r <- rast(f)
  resample(r, LM_100, method = 'bilinear')
})

EP_layer_names <- tools::file_path_sans_ext(basename(EP_list))

EP_rasters <- rast(EP_list_resampled)
names(EP_rasters) <- EP_layer_names

writeRaster(EP_rasters, 'predictor_variables/EP_resampled.tif')
EP_rasters <- rast('predictor_variables/EP_resampled.tif')

EP_PCA <- rasterPCA(EP_rasters, spca = TRUE, maskCheck = TRUE)
saveRDS(EP_PCA, 'rf_maps/PCA_EP.rds')

explained_var_cum <- cumsum(EP_PCA$model$sdev^2 / sum(EP_PCA$model$sdev^2))

n_comp <- which(explained_var_cum >= .95)[1]

for (k in 1:n_comp){
  fname <- paste0('predictor_variables/predictor_variables_resampled/', "EP_PC", k, ".tif")
  writeRaster(EP_PCA$map[[k]], fname, overwrite = TRUE)
}

## Distance-based Predictors ####

folder_path <- 'predictor_variables/disturbance_indicators'
IF_list <- list.files(path = folder_path, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)

# Resample
IF_list_resampled <- lapply(IF_list, function(f){
  r <- rast(f)
  resample(r, LM_100, method = 'bilinear')
})

IF_layer_names <- tools::file_path_sans_ext(basename(IF_list))

IF_rasters <- rast(IF_list_resampled)
names(IF_rasters) <- IF_layer_names

writeRaster(IF_rasters, 'predictor_variables/IF_resampled.tif')
IF_rasters <- rast('predictor_variables/IF_resampled.tif')

IF_PCA <- rasterPCA(IF_rasters, spca = TRUE, maskCheck = TRUE)
saveRDS(IF_PCA, 'rf_maps/PCA_IF.rds')

explained_var_cum <- cumsum(IF_PCA$model$sdev^2 / sum(IF_PCA$model$sdev^2))

n_comp <- which(explained_var_cum >= .95)[1]

for (k in 1:n_comp){
  fname <- paste0('predictor_variables/predictor_variables_resampled/', "IF_PC", k, ".tif")
  writeRaster(IF_PCA$map[[k]], fname, overwrite = TRUE)
}

## Variance Inflation Factor ####
# Assess multicollinearity among predictors

# Load all predictors

folder_path <- 'predictor_variables/predictor_variables_resampled'
predictors_list <- list.files(path = folder_path, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)

predictors_layer_names <- tools::file_path_sans_ext(basename(predictors_list))

predictors_rasters <- rast(predictors_list)
names(predictors_rasters) <- predictors_layer_names

predictors_rasters$LULC <- as.factor(predictors_rasters$LULC)

# VIF for each model

# Baseline model
baseline_predictors <- predictors_rasters[[grep("^(EP|IF|LULC)", names(predictors_rasters), value = TRUE)]]
vif_baseline <- usdm::vif(na.omit(baseline_predictors))

# LM model
lm_predictors <- predictors_rasters[[grep("^(EP|IF|LULC|LM250)", names(predictors_rasters), value = TRUE)]]
vif_lm <- usdm::vif(na.omit(lm_predictors))

# HSI integrated model
HSI250_predictors <- predictors_rasters[[grep("^(EP|IF|LULC|LM250|HSI_250)", names(predictors_rasters), value = TRUE)]]
vif_250 <- usdm::vif(na.omit(HSI250_predictors))

# Fine scale model
HSI100_predictors <- predictors_rasters[[grep("^(EP|IF|LULC|LM100|HSI_100)", names(predictors_rasters), value = TRUE)]]
vif_100 <- usdm::vif(na.omit(HSI100_predictors))

# Large scale model
HSI500_predictors <- predictors_rasters[[grep("^(EP|IF|LULC|LM500|HSI_500)", names(predictors_rasters), value = TRUE)]]
vif_500 <- usdm::vif(na.omit(HSI500_predictors))


# Presence data (hunting data) ####
# Source: South Tyrolean Hunting Association

hunting_districts <- read_sf('Jagdreviere_nach_1983/Jagdreviere_nach_1983.shp') %>% 
  st_transform(crs = 32632)

area_ST <- hunting_districts$geometry %>%
  st_union() %>%
  st_as_sf()

presence_data <- st_read('RedDeer.gpkg')

presence_data$datetime <- as.POSIXct(
  paste(presence_data$Datum, presence_data$Zeit),
  format = "%d.%m.%Y %H:%M",
  tz = "Europe/Berlin"
)

# Select only records in winter period within South Tyrol
presence_data_winter <- presence_data %>%
  filter(
    month(datetime) %in% 10:12 | 
      month(datetime) %in% 1:4
  ) %>%
  filter(st_within(geom, area_ST, sparse = FALSE))

st_write(presence_data_winter, 'Red_Deer_Winter.gpkg')
presence_data_winter <- read_sf('Red_Deer_Winter.gpkg')


## Pseudo-absence data ####

# Area of South Tyrol without Stelvio National Park (area for model calculation)

hunting_area <- subset(hunting_districts, JAGDBEZIRK != 'Nationalpark' | is.na(JAGDBEZIRK)) %>%
  st_union()

# Randomly sample 50,000 pseudo-absence points
background_data <- spatSample(vect(hunting_area), size = 50000, method = 'random')
background_data <- st_as_sf(background_data) %>% 
  mutate(occ = 0)

## Combined occurrence data ####

presence_data <- presence_data_winter %>% 
  select(geom) %>% 
  mutate(occ = 1) %>% 
  rename(geometry = geom)

occurrence_data <- bind_rows(presence_data, background_data)

# Block cross-validation

# Determine block size
sac1 <- cv_spatial_autocor(r = predictors_rasters, num_sample = 5000) # Block size ~ 2000m

# Assign folds
sb1 <- cv_spatial(x = occurrence_data,
                  column = "occ",
                  k = 5,
                  size = 2000,
                  selection = "random",
                  iteration = 50)

occurrence_data$folds <- sb1$folds_ids

st_write(occurrence_data, 'occurrence_data.gpkg')
occurrence_data <- read_sf('occurrence_data.gpkg')

# Predictor variable extraction
variable_extract <- terra::extract(x = predictors_rasters[[1:nlyr(predictors_rasters)]],
                            y = occurrence_data,
                            ID = FALSE)

occurrence_data_extract <- cbind(occurrence_data, variable_extract) %>% 
  as.data.frame() %>% 
  select((-geom)) %>% 
  na.omit()

write.table(occurrence_data_extract,'predictor_variables/predictor_variables_resampled/occurence_data_extract.txt')
occurrence_data_extract <- read.table('predictor_variables/predictor_variables_resampled/occurence_data_extract.txt')

training_data <- occurrence_data_extract[occurrence_data_extract$folds < 4,] # folds 1-3 (~ 60%)
training_data$occ <- as.factor(training_data$occ)
training_data$LULC <- as.factor(training_data$LULC)

testing_data <- occurrence_data[occurrence_data$folds >= 4,] # folds 4-5 (~ 40%)
testing_data$occ <- as.factor(testing_data$occ)


# SDM ####

# Baseline model: Environmental predictors, LULC at point-scale

baseline_data <- training_data %>% 
  select(occ, starts_with('EP'), starts_with('IF'), LULC) %>% 
  na.omit()

# Landscape metrics model: Standard model + Landscape Metrics 250m radius

LM250_data <- training_data %>% 
  select(occ, starts_with('EP'), starts_with('IF'), LULC, starts_with(('LM250'))) %>% 
  na.omit()

# HSI-integrated model: Standard model + Landscape Metrics & HSI in 250m radius

LM_HSI250_data <- training_data %>% 
  select(occ, starts_with('EP'), starts_with('IF'), LULC, starts_with(('LM250')), HSI_250) %>% 
  na.omit()

# Fine-scale model 100m: Standard model + Landscape Metrics & HSI in 100m radius

LM_HSI100_data <- training_data %>% 
  select(occ, starts_with('EP'), starts_with('IF'), LULC, starts_with(('LM100')), HSI_100) %>% 
  na.omit()

# Larger-scale model 500m: Standard model + Landscape Metrics & HSI in 500m radius

LM_HSI500_data <- training_data %>% 
  select(occ, starts_with('EP'), starts_with('IF'), LULC, starts_with(('LM500')), HSI_500) %>% 
  na.omit()


# Create all SDM maps in one loop

model_data_list <- list(baseline_data, LM250_data, LM_HSI250_data, LM_HSI100_data, LM_HSI500_data)
model_data_names <- c('Baseline', 'LM_250m', 'LM+HSI_250m', 'LM+HSI_100m', 'LM+HSI_500m')

path = 'rf_maps_new/'

for (model_number in 1:length(model_data_list)) {
  prNum <- as.numeric(table(model_data_list[[model_number]]$occ)['1'])
  smpsize <- c('0' = prNum, '1' = prNum)
  
  rp <- rfPermute(formula = occ~., 
                  data = model_data_list[[model_number]],
                  ntree = 1000,
                  sampsize = smpsize,
                  importance = TRUE,
                  replace = TRUE,
                  num.rep = 100)
  
  saveRDS(rp, file = paste0(path, 'rp_', model_data_names[model_number], '.rds'))
  
  rf <- rp$rf
  
  rf_map <- predict(predictors_rasters[[colnames(model_data_list[[model_number]])[2:ncol(model_data_list[[model_number]])]]],
                    rf,
                    type = 'prob')
  
  writeRaster(rf_map, paste0(path, 'rf_map_', model_data_names[model_number], '.tif'))
}

## Results RF Models ####

# Read model outputs

rp_Baseline <- readRDS('rf_maps_new/rp_Baseline.rds')
rp_LM250 <- readRDS('rf_maps_new/rp_LM_250m.rds')
rp_LM_HSI100 <- readRDS('rf_maps_new/rp_LM+HSI_100m.rds')
rp_LM_HSI250 <- readRDS('rf_maps_new/rp_LM+HSI_250m.rds')
rp_LM_HSI500 <- readRDS('rf_maps_new/rp_LM+HSI_500m.rds')

rp_list <- list(Baseline = rp_Baseline, 
                LM250 = rp_LM250,
                LM_HSI100 = rp_LM_HSI100,
                LM_HSI250 = rp_LM_HSI250,
                LM_HSI500 = rp_LM_HSI500)


acc_df <- do.call(rbind, lapply(names(rp_list), function(m) {
  cm <- confusionMatrix(rp_list[[m]], conf.level = 0.95)  # overall + CI
  
  overall_row <- nrow(cm)
  data.frame(
    Model      = m,
    Accuracy   = cm[overall_row, "pct.correct"],        # %
    LCI_95     = cm[overall_row, "LCI_0.95"],           # lower CI (%)
    UCI_95     = cm[overall_row, "UCI_0.95"]            # upper CI (%)
  )
}))

acc_df$Model <- factor(acc_df$Model, levels = names(rp_list))
write.csv(acc_df, 'plots/oob_results.csv')

# Evaluation ####

## Count Data ####

count_data <- read_sf('CountDataRedDeer/count_data_processed.gpkg')
count_data$mean_value <- rowMeans(st_drop_geometry(count_data)[,4:8], na.rm = TRUE)

count_data <- count_data %>% 
  mutate(
    area_km2 = st_area(geom) / 1e6,
    count_per_area = mean_value / as.numeric(area_km2))

st_write(count_data, 'CountDataRedDeer/count_data_final.gpkg')
count_data <- st_read('CountDataRedDeer/count_data_final.gpkg')

# Load rf_maps

folder_path <- 'rf_maps_new'

rf_maps_list <- list.files(path = folder_path, pattern = "\\.tif$", full.names = TRUE)

rf_maps <- list()

for (file in rf_maps_list) {
  r <- rast(file)
  x1_layer <- r[['X1']]
  rf_maps <- c(rf_maps, x1_layer)
}

rf_maps <- rast(rf_maps)
names(rf_maps) <- c('Baseline', 'LM_250m', 'LM_HSI_100m', 'LM_HSI_250m',  'LM_HSI_500m')

results_df <- data.frame(Region = count_data$Jagdbezirk, count_per_area = count_data$count_per_area)

for (i in 1:nlyr(rf_maps)) {
  mean_values <- terra::extract(rf_maps[[i]], vect(count_data), fun = mean, na.rm = TRUE)
  new_col <- paste0(names(rf_maps[[i]]))
  results_df[[new_col]] <- mean_values[,2]

}

## Spearman's rank correlation coefficient ####

results <- results_df %>%
  pivot_longer(
    cols = 3:7,
    names_to = "model",
    values_to = "pred"
  ) %>%
  group_by(model) %>%
  summarise(
    spearman = cor(count_per_area, pred, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )

## Continuous Boyce Index ####

boyce_values <- sapply(rf_maps, function(map){
  boyce_idx <- Boyce(obs = sf::st_coordinates(testing_data[testing_data$occ ==1,]), pred = map)$Boyce
  return(boyce_idx)
})

results$boyce <- boyce_values

write.csv(results, 'plots/results.csv')
results <- read.csv('plots/results.csv', sep = ';')

## Quantile Regression ####

model_cols <- names(rf_maps)
taus <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)

# Function to compute R1 for a given rq model
compute_R1 <- function(fit, tau) {
  rho <- function(u, tau) abs(u) * (tau - (u < 0))
  resid_fit <- residuals(fit)
  Vhat <- sum(rho(resid_fit, tau))
  null_fit <- rq(count_per_area ~ 1, tau = tau, data = results_df)
  resid_null <- residuals(null_fit)
  V0 <- sum(rho(resid_null, tau))
  R1 <- 1 - Vhat / V0
  return(R1)
}

# Initialize results list
results_list <- list()

# Loop over each model
for (mod_col in model_cols) {
  model_results <- data.frame(tau = taus)
  
  for (tau in taus) {
    # Fit quantile regression
    fit <- rq(reformulate(mod_col, "count_per_area"), tau = tau, data = results_df)
    
    # Extract slope
    slope <- coef(fit)[2]
    
    # Compute R1
    R1 <- compute_R1(fit, tau)
    
    model_results[model_results$tau == tau, paste0(mod_col, "_slope")] <- slope
    model_results[model_results$tau == tau, paste0(mod_col, "_R1")] <- R1
  }
  
  results_list[[mod_col]] <- model_results
}

# Combine all model results into one wide table
results_table <- results_list[[1]]
for (i in 2:length(results_list)) {
  results_table <- full_join(results_table, results_list[[i]], by = "tau")
}

# Reorder columns: tau first, then slope and R1 for each model
slope_cols <- paste0(model_cols, "_slope")
r1_cols <- paste0(model_cols, "_R1")
results_table <- results_table[, c('tau', slope_cols, r1_cols)]

write.csv(results_table, 'plots/quantile_regression_results.csv')
results_table <- read.csv('plots/quantile_regression_results.csv', sep = ',')