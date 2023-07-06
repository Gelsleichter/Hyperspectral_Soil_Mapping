#################################################################################################################
## R script Hyperspectral modeling for Itatiaia National park, Rio de Janeiro, Minas Gerais, Brazil             #
## Version: 1, Sep 2022                                                                                         #
## Author: Dr. Yuri Andrei Gelsleichter, parts adapted from: Dr. Elias Mendes Costa                             #
## License: CC-BY-NC-SA                                                                                         #
#################################################################################################################
#       ______        __  _       _           _   __      __  _                   __    ____             __     #
#      /  _/ /_____ _/ /_(_)___ _(_)___ _    / | / /___ _/ /_(_)___  ____  ____ _/ /   / __ \____ ______/ /__   #
#      / // __/ __ `/ __/ / __ `/ / __ `/   /  |/ / __ `/ __/ / __ \/ __ \/ __ `/ /   / /_/ / __ `/ ___/ //_/   #
#    _/ // /_/ /_/ / /_/ / /_/ / / /_/ /   / /|  / /_/ / /_/ / /_/ / / / / /_/ / /   / ____/ /_/ / /  / ,<      #
#   /___/\__/\__,_/\__/_/\__,_/_/\__,_/   /_/ |_/\__,_/\__/_/\____/_/ /_/\__,_/_/   /_/    \__,_/_/  /_/|_|     #
#                                                                                                               #
#  _  _                                    _            _   ___      _ _   __  __                _              #
# | || |_  _ _ __  ___ _ _ ____ __  ___ __| |_ _ _ __ _| | / __| ___(_) | |  \/  |__ _ _ __ _ __(_)_ _  __ _    #
# | __ | || | '_ \/ -_) '_(_-< '_ \/ -_) _|  _| '_/ _` | | \__ \/ _ \ | | | |\/| / _` | '_ \ '_ \ | ' \/ _` |   #
# |_||_|\_, | .__/\___|_| /__/ .__/\___\__|\__|_| \__,_|_| |___/\___/_|_| |_|  |_\__,_| .__/ .__/_|_||_\__, |   #
#       |__/|_|              |_|                                                      |_|  |_|         |___/    #
# figlet -f                                                                                                     #
#################################################################################################################
### ____        _                                                  _   _               ###
### |  _ \  __ _| |_ __ _      _ __  _ __ ___ _ __   __ _ _ __ __ _| |_(_) ___  _ __   ###
### | | | |/ _` | __/ _` |    | '_ \| '__/ _ \ '_ \ / _` | '__/ _` | __| |/ _ \| '_ \  ###
### | |_| | (_| | || (_| |    | |_) | | |  __/ |_) | (_| | | | (_| | |_| | (_) | | | | ###
### |____/ \__,_|\__\__,_|    | .__/|_|  \___| .__/ \__,_|_|  \__,_|\__|_|\___/|_| |_| ###
###                           |_|            |_|                                       ###
###                                                                                    ### 
########################################################################################## 

#' ### Auto setwd
library(rstudioapi); current_path <- getActiveDocumentContext()$path; setwd(dirname(current_path)) ### https://eranraviv.com/r-tips-and-tricks-working-directory/
getwd()
gc(); rm(list=ls())

### Sentinel notes
# Cloud cover percentage: 0.035851 (0.036%)
# Image Sentinell 2 MSI, processing level 2A with 0.036% of cloud cover in the middle of dry season 
# (older images were offline from the page, impossible to download, and even with more clouds rate, they were offline)

################################################################################
### Load the images
library(raster)
library(rgdal)
library(terra)
path_img_folder <- "../input/raster/S2A_MSIL2A_20220814T131301_N0400_R138_T23KNR_20220814T201157.SAFE/GRANULE/L2A_T23KNR_A037315_20220814T131255/IMG_DATA/R10m/"
Sent_b2_B <- terra::rast(paste0(path_img_folder, "T23KNR_20220814T131301_B02_10m.jp2")) # blue 490 nm
Sent_b3_G <- terra::rast(paste0(path_img_folder, "T23KNR_20220814T131301_B03_10m.jp2")) # green 560 nm
Sent_b4_R <- terra::rast(paste0(path_img_folder, "T23KNR_20220814T131301_B04_10m.jp2")) # red 665 nm
Sent_b8_NIR <- terra::rast(paste0(path_img_folder, "T23KNR_20220814T131301_B08_10m.jp2")) # nir 842 nm
# TCi <- terra::rast(paste0(path_img_folder, "T23KNR_20220814T131301_TCi_10m.jp2")) # true color image

################################################################################
### Load vectors
library(ggplot2)
library(RStoolbox)
library(ggspatial)
crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
pni <- terra::vect("../input/shp/PARNA_Itatiaia_SIRGAS_IMG.shp", crs= crs_ref)
pni_alta <- terra::vect("../input/shp/parte alta.shp", crs= crs_ref)
pnipts <- terra::vect("../input/shp/amostragem_cLHS_100.shp", crs= crs_ref)

pni_bf3km <- buffer(pni, width=3000)
pni_alta_bf1hm <- buffer(pni_alta, width=100)

plot(pni_bf3km)
plot(pni, add=T)
plot(pni_alta, add= T)
plot(pni_alta_bf1hm, add= T)

################################################################################
### Crop rasters img by shapefile
# Sent2 <- c(Sent_b2_B, Sent_b3_G, Sent_b4_R, Sent_b8_NIR)
# Sent2 <- crop(Sent2, pni_alta_bf1hm, mask= F)
## terra::plot(Sent2[[1]])
# plotRGB(Sent2, 3, 2, 1, stretch='lin') # RGB

# Sent_b2_B <- crop(Sent_b2_B, pni_alta_bf1hm, mask= F)
# Sent_b3_G <- crop(Sent_b3_G, pni_alta_bf1hm, mask= F)
# Sent_b4_R <- crop(Sent_b4_R, pni_alta_bf1hm, mask= F)
# Sent_b8_NIR <- crop(Sent_b8_NIR, pni_alta_bf1hm, mask= F)

Sent_b2_B <- terra::crop(Sent_b2_B, pni_alta_bf1hm, snap="out", mask= T)
Sent_b3_G <- terra::crop(Sent_b3_G, pni_alta_bf1hm, snap="out", mask= T)
Sent_b4_R <- terra::crop(Sent_b4_R, pni_alta_bf1hm, snap="out", mask= T)
Sent_b8_NIR <- terra::crop(Sent_b8_NIR, pni_alta_bf1hm, snap="out", mask= T)

## terra::plot(Sent_b2_B)

################################################################################
### Load, crop and adjust DEM
path_dem_folder <- "../input/raster/AP_27029_FBS_F6730_RT1__HI_RES_Terrain_Corrected/"
dem <- terra::rast(paste0(path_dem_folder, "AP_27029_FBS_F6730_RT1.dem.tif"))

### Crop dem by shapefile
# dem_cp <- crop(dem, pni_alta_bf1hm, mask= F)
dem_cp <- terra::crop(dem, pni_alta_bf1hm, snap="out", mask= T)

### Adust dem from 12.5 to 10 m
dem_2.5 <- terra::disagg(dem_cp, fact= 5, method="bilinear")
# method: Either "near" for nearest or "bilinear" for bilinear interpolation

### Adust dem from 2.5 to 10 m
dem_10 <- terra::aggregate(dem_2.5, fact= 4, fun=mean, cores= 4) 
## plot(dem_10)

### Register (resample): certify that all rasters (img and dem) are same center cells and origin
# x <- resample(x, y, method="bilinear")
# x: SpatRaster to be resampled
# y: SpatRaster with the geometry that x should be resampled to
# method: character. Method used for estimating the new cell values. One of:
# near: nearest neighbor. This method is fast, and it can be the preferred method if the cell values represent classes. It is not a good choice for 
# continuous values. This is used by default if the first layer of x is categorical.
# bilinear: bilinear interpolation. This is the default if the first layer of x is numeric (not categorical).
dem <- resample(dem_10, Sent_b2_B, method="bilinear")
# The resample function: 
# It can do either nearest neighbor assignments (for categorical data) or 
# bilinear interpolation (for numerical data).

### Write on disk
terra::writeRaster(dem, filename="../output/raster/dem_10m.tif", overwrite=TRUE)

################################################################################
### Geology and geomorfology
Geology <- terra::rast("../output/raster/geology_bf_25m.tif")
            # Geology <- terra::rast("../input/raster/geologia_25.tif")
crs(Geology)  <- "epsg:32723"
Geomorfology <- terra::rast("../output/raster/geomorfology_bf_25m.tif")
            # Geomorfology <- terra::rast("../input/raster/geomorfolo_25.tif")
crs(Geomorfology)  <- "epsg:32723"
# You can also use PROJ-string notation
# crs(Geomorfology) <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

### Adust dem from 25 to 5 m
Geology5 <- terra::disagg(Geology, fact= 5, method="near")
Geomorfology5 <- terra::disagg(Geomorfology, fact= 5, method="near")

### Adust dem from 5 to 10 m
Geology <- terra::aggregate(Geology5, fact= 2, fun=mean, cores= 4) 
Geomorfology <- terra::aggregate(Geomorfology5, fact= 2, fun=mean, cores= 4) 

### Register (resample)
Geology <- resample(Geology, Sent_b2_B, method="near")
Geomorfology <- resample(Geomorfology, Sent_b2_B, method="near")

### crop and mask
# Geology <- crop(Geology, pni_alta_bf1hm, mask= F)
# Geomorfology <- crop(Geomorfology, pni_alta_bf1hm, mask= F)
Geology <- terra::crop(Geology, pni_alta_bf1hm, snap="out", mask= T)
Geomorfology <- terra::crop(Geomorfology, pni_alta_bf1hm, snap="out", mask= T)

Geology <- as.factor(Geology)
Geomorfology <- as.factor(Geomorfology)
# Geology <- as.int(Geology)
# Geomorfology <- as.int(Geomorfology)

levels(Geology)
levels(Geomorfology)

### Write on disk to make all intervals (numbers) as integer (INT1U)
terra::writeRaster(Geology, filename="../output/raster/geology_10m.tif", 
                   datatype= "INT1U", overwrite=TRUE)
terra::writeRaster(Geomorfology, filename="../output/raster/geomorfology_10m.tif", 
                   datatype= "INT1U", overwrite=TRUE)
Geology <- terra::rast("../output/raster/geology_10m.tif")
Geomorfology <- terra::rast("../output/raster/geomorfology_10m.tif")

# dataType(raster(Geology)) # "dataType" is a raster package function
# dataType(raster(Geomorfology)) # "dataType" is a raster package function

################################################################################
### NDVI, SAVI: https://www.usgs.gov/landsat-missions/landsat-normalized-difference-vegetation-index

# ndvi= (nir-red)/(nir+red)
# savi= ((nir-red)/(nir+red + 0.5))*(1.5)
### Generating SAVI and NDVI (raster)  
# NDVI <- raster::overlay(Sent2[[4]], Sent2[[3]], fun=function(x,y){(x-y)/(x+y)})
# SAVI <- raster::overlay(Sent2[[4]], Sent2[[3]], fun=function(x,y){(1.5)*(x-y)/(x+y+0.5)})
# for SAVI adjust factor de ajuste L = 0.5

### Generating SAVI and NDVI (terra)
# NDVI <- terra::lapp(c(Sent2[[4]], Sent2[[3]]), fun=function(x,y){(x-y)/(x+y)})
# SAVI <- terra::lapp(c(Sent2[[4]], Sent2[[3]]), fun=function(x,y){(1.5)*(x-y)/(x+y+0.5)})
NDVI <- terra::lapp(c(Sent_b4_R, Sent_b8_NIR), fun=function(x,y){(x-y)/(x+y)})
SAVI <- terra::lapp(c(Sent_b4_R, Sent_b8_NIR), fun=function(x,y){(1.5)*(x-y)/(x+y+0.5)})
# plot(NDVI)
# plot(SAVI)
################################################################################

### Exporting Treated rasters
# terra::writeRaster(dem, filename="../raster/output/10m_dem_parte_alta.tif", overwrite=TRUE)

# writeVector(x, filename, filetype=NULL, layer=NULL, insert=FALSE,
#             overwrite=FALSE, options="ENCODING=UTF-8")
# shapefile(pni_bf3km, filename="../shapes/buffer_parte_alta_3km.shp")

gc(); rm(dem_2.5, dem_10, dem_cp, Geology5, Geomorfology5)

ext(Geology)
ext(Geomorfology)
ext(dem)
ext(NDVI)
ext(SAVI)

###################################################################################################
###################################################################################################
###################################################################################################
#   ___                                         _      __                        __           __
#  / _ \_______ ___      _______ _  _____ _____(_)__ _/ /____ ___     ____  ___ / /____ _____/ /__
# / ___/ __/ -_) _ \_   / __/ _ \ |/ / _ `/ __/ / _ `/ __/ -_|_-< __ / __/ (_-</ __/ _ `/ __/  '_/
#/_/  /_/  \__/ .__(_)  \__/\___/___/\_,_/_/ /_/\_,_/\__/\__/___/   /_/    /___/\__/\_,_/\__/_/\_\
#            /_/
# figlet -f smslant Prep. covariates r stack ### run on terminal ctrl+alt+enter
# showfigfonts
###########################################################################
### Covariates preparation 
###########################################################################
#' ### Auto setwd
library(rstudioapi); current_path <- getActiveDocumentContext()$path; setwd(dirname(current_path)) ### https://eranraviv.com/r-tips-and-tricks-working-directory/
getwd()
# gc(); rm(list=ls())

# library(sp)
# library(raster)
# library(rgdal)
# library(maptools)

{
################################################################################################
#### Derive terrain attributes with SAGA-GIS from RSAGA package ####
## para usar as funções desse pacote você precisa ter o saga instalado SAGA GIS (2.3 - 6.4.0)

### Você deve setar o caminho onde está os modulos do saga. Se você usa linux deve ser: 
### Path to SAGA modules: Linux: 
library(RSAGA)
rsaga.env()
env = rsaga.env()
## env = rsaga.env(workspace = "../raster/input")
## env = rsaga.env(workspace = "../PNI_R_spModeling/raster/input")
## env <- RSAGA::rsaga.env(path="/usr/bin")
## env = RSAGA::rsaga.env(workspace = "../raster", modules="/usr/lib/x86_64-linux-gnu/saga") ### Antigo, veio de Elias
## env = rsaga.env(workspace = "../input/raster", path = "/usr/bin", modules = "/usr/lib/x86_64-linux-gnu/saga")
env = rsaga.env(workspace = "../output", path = "/usr/bin", modules = "/usr/lib/x86_64-linux-gnu/saga")

### Se você usa aquele outro sistema (windows) deve ser:   
### Path to SAGA modules: Windows: 
env = RSAGA::rsaga.env(workspace = "../output", ### must be one directory above the getwd() (in my case getwd() = script, always)
                       path="C:/Program Files (x86)/SAGA-GIS",
                       modules="C:/Program Files (x86)/SAGA-GIS/tools")
### RSAGA needs SAGA-GIS 6.3.0 or older

# Test SAGA functions (examples)
# RSAGA::rsaga.get.libraries()
# RSAGA::rsaga.get.libraries(path = env$modules)
# RSAGA::rsaga.get.modules(libs = "ta_morphometry", env = env)
# RSAGA::rsaga.get.modules(libs = "ta_compound", env = env)
# RSAGA::rsaga.get.modules(libs = "ta_lighting", env = env)
# RSAGA::rsaga.get.usage("ta_compound" , module = 0)  ### http://www.saga-gis.org/saga_tool_doc/ ### http://www.saga-gis.org/saga_tool_doc/2.3.0/ta_compound_0.html
## RSAGA::rsaga.get.modules("ta_morphometry")[[1]][, -3]
## RSAGA::rsaga.search.modules("sink")
## RSAGA::rsaga.get.modules()

## Preencher as depressões expúrias do DEM
RSAGA::rsaga.fill.sinks(in.dem="../output/raster/dem_10m.tif", out.dem="raster/terrain/DEM_10m", method = "wang.liu.2006", env=env)
# RSAGA::rsaga.fill.sinks(in.dem="../raster/input/MDE.tif", out.dem="DEM", method = "wang.liu.2006", env=env)

### Derive terrain attributes from DEM
# RSAGA::rsaga.geoprocessor(lib="ta_compound", module = 0, 
#                           param = list(ELEVATION = paste(getwd(),"/raster/terrain/DEM_10m.sgrd", sep = ""),
#                                        SHADE = paste(getwd(),"/raster/terrain/shade.sgrd", sep = ""),
#                                        SLOPE = paste(getwd(),"/raster/terrain/slope.sgrd", sep = ""), # slope default (radians)
#                                        ASPECT = paste(getwd(),"/raster/terrain/aspect.sgrd", sep = ""), # aspect default (radians, north=0)
#                                        HCURV = paste(getwd(),"/raster/terrain/plan_curv.sgrd", sep = ""),
#                                        VCURV = paste(getwd(),"/raster/terrain/prof_curv.sgrd", sep = ""),
#                                        CONVERGENCE = paste(getwd(),"/raster/terrain/convergence.sgrd", sep = ""),
#                                        FLOW = paste(getwd(),"/raster/terrain/cat_area.sgrd", sep = ""),
#                                        WETNESS = paste(getwd(),"/raster/terrain/twi.sgrd", sep = ""),
#                                        LSFACTOR = paste(getwd(),"/raster/terrain/ls_factor.sgrd", sep = ""),
#                                        RSP = paste(getwd(),"/raster/terrain/rsp.sgrd", sep = ""),
#                                        CHANNELS = paste(getwd(),"/raster/terrain/channels.sgrd", sep = ""),
#                                        CHNL_DIST = paste(getwd(),"/raster/terrain/chnd.sgrd", sep = ""),
#                                        CHNL_BASE = paste(getwd(),"/raster/terrain/chnb.sgrd", sep = "")
#                                        ), env = env)


RSAGA::rsaga.geoprocessor(lib="ta_compound", module = 0, 
                          param = list(ELEVATION = "raster/terrain/DEM_10m.sgrd",
                                       SHADE = "raster/terrain/shade.sgrd",
                                       # SLOPE = "raster/terrain/slope.sgrd", # slope default (radians)
                                       # ASPECT = "raster/terrain/aspect.sgrd", # aspect default (radians, north=0)
                                       HCURV = "raster/terrain/plan_curv.sgrd",
                                       VCURV = "raster/terrain/prof_curv.sgrd",
                                       CONVERGENCE = "raster/terrain/convergence.sgrd",
                                       FLOW = "raster/terrain/cat_area.sgrd",
                                       WETNESS = "raster/terrain/twi.sgrd",
                                       LSFACTOR = "raster/terrain/ls_factor.sgrd",
                                       RSP = "raster/terrain/rsp.sgrd",
                                       CHANNELS = "raster/terrain/channels.sgrd",
                                       CHNL_DIST = "raster/terrain/chnd.sgrd",
                                       CHNL_BASE = "raster/terrain/chnb.sgrd"
                          ), env = env)

### Derive slope and aspect in degrees and slope in percentage 
RSAGA::rsaga.slope.asp.curv(in.dem ="../output/raster/terrain/DEM_10m.sgrd", 
                            out.slope = "raster/terrain/slope", 
                            out.aspect = "raster/terrain/aspect", 
                            unit.slope = 2, # 2 "percent", 1 "degrees", 0 "radians"
                            unit.aspect = 1, # 1 "degrees", 0 "radians"
                            env=env)
### Load terrain attributes 
{
  DEM <- terra::rast("../output/raster/terrain/DEM_10m.sdat") # Digital elevation model unit = m
  Slope <- terra::rast("../output/raster/terrain/slope.sdat")   # slope unit = %
  Shade <- terra::rast("../output/raster/terrain/shade.sdat")   # unit = degrees
  Aspect <- terra::rast("../output/raster/terrain/aspect.sdat") # aspect unit = degrees
  Northernness <- (abs(180-Aspect)) # northernness unit = degrees
  Plan_curv <- terra::rast("../output/raster/terrain/plan_curv.sdat") # plan curvature unit = 1/m
  Prof_curv <- terra::rast("../output/raster/terrain/prof_curv.sdat") # profile curvature unit = 1/m
  Convergence <- terra::rast("../output/raster/terrain/convergence.sdat") # Convergence index unit = %
  Cat_area <- terra::rast("../output/raster/terrain/cat_area.sdat") #  catchment area unit = m2
  TWI <- terra::rast("../output/raster/terrain/twi.sdat") # Topographic Wetness Index unit = dimensionless
  LS_factor <- terra::rast("../output/raster/terrain/ls_factor.sdat") # LS-Factor unit = dimensionless
  RSP <- terra::rast("../output/raster/terrain/rsp.sdat") # Relative slope position unit = dimensionless
  ChND <- terra::rast("../output/raster/terrain/chnd.sdat") # Channel network distance unit = m
  ChNB <- terra::rast("../output/raster/terrain/chnb.sdat") # Channel network base level unit = m
}

### Stack all rasters
stk <- c(DEM, 
         Slope, 
         Shade, 
         Aspect, 
         Northernness, 
         Plan_curv, 
         Prof_curv, 
         Convergence, 
         Cat_area, 
         TWI, 
         LS_factor, 
         RSP, 
         ChND, 
         ChNB,
         Sent_b2_B,
         Sent_b3_G,
         Sent_b4_R,
         Sent_b8_NIR,
         NDVI,
         SAVI,
         Geology,
         Geomorfology)

names(stk)
{
names(stk[[1]]) <- "DEM"
names(stk[[2]]) <- "Slope"
names(stk[[3]]) <- "Shade"
names(stk[[4]]) <- "Aspect"
names(stk[[5]]) <- "Northernness"
names(stk[[6]]) <- "Plan_curv"
names(stk[[7]]) <- "Prof_curv"
names(stk[[8]]) <- "Convergence"
names(stk[[9]]) <- "Cat_area"
names(stk[[10]]) <- "TWI"
names(stk[[11]]) <- "LS_factor"
names(stk[[12]]) <- "RSP"
names(stk[[13]]) <- "ChND"
names(stk[[14]]) <- "ChNB"
names(stk[[15]]) <- "Sent_b2_B"
names(stk[[16]]) <- "Sent_b3_G"
names(stk[[17]]) <- "Sent_b4_R"
names(stk[[18]]) <- "Sent_b8_NIR"
names(stk[[19]]) <- "NDVI"
names(stk[[20]]) <- "SAVI"
names(stk[[21]]) <- "Geology"
names(stk[[22]]) <- "Geomorfology"
}

stk_dsm <- stk[[c(1,2,3,4,6,9,13,17,18,19,21,22)]]
# names(stk_dsm)
# a <- c(1,2,3,4,6,9,13,17,18,19,21,22)
# l <- 1:22
# l[!l %in% a] %>% length()
stk_subsuf_img <- stk[[c(5,7,8,10,11,12,14,15,16,20)]]
# names(stk_subsuf_img)

### Write on disk
# terra::writeRaster(stk, filename="../output/raster/stk_10m.tif", overwrite=TRUE)
terra::writeRaster(stk_dsm, filename="../output/raster/stk_dsm_10m.tif", overwrite=TRUE)
terra::writeRaster(stk_subsuf_img, filename="../output/raster/stk_subsuf_img_10m.tif", overwrite=TRUE)
}

################################################################################
### Correlation, to perform the corr. it is necessary to "convert" all rasters 
### in a dataframe, same data is used in prediction (df.pred)
################################################################################

### Extract the spatial location of every pixel in the area of the image
### Load the raster stack
stk_dsm <- terra::rast("../output/raster/stk_dsm_10m.tif")

### With Raster package
# sp_points <- stk_dsm[[1]]
# sp_points[sp_points] <- 1 ### fill all with 1
# sp_points <- setValues(sp_points, 1) ### fill all with 1

# library(raster)
# sp_points <- raster(sp_points)
# names(sp_points) <- "cloud_points"
# upper_points <- rasterToPoints(sp_points, fun=NULL, spatial=T) # convert pixels in points raster
# ### here package terra was slower
# 
# ### Extract values from all covariates using all points from the raster "image"
# library(snow) ### parallel
# stk_dsm_r <- stack(stk_dsm) # convert from terra (SpatRaster) to raster (RasterLayer)
# 
# beginCluster()
# system.time({a = extract(stk_dsm_r,
#               upper_points,
#               method="simple", sp= T)#;
#               # df.pred = as.data.frame(a)
#               })
# endCluster()
# 
# ###Remove points
# df.pred$clous_points <- NULL

###################################################################################
### Since terra package has a "weird" behave when crop (no matter how you set), ###
### the output raster seems smaller, because after run the RSAGA functions some ###
### rasters (like half) get some NaNs. Then I used a bubffer 100 m (parte_alta) ###
### to run all pre processing, but to extract I just used the parte_alta shape. ###
###################################################################################

### With Terra package
# upper_points <- as.points(sp_points, values=F, na.rm=TRUE, na.all=T) # convert pixels in points terra
# writeVector(upper_points, "../output/shp/upper_points.shp")

system.time({df.pred = terra::extract(stk_dsm, 
                                      pni_alta, 
                                      method= "simple", 
                                      xy= T, 
                                      ID= F)})

dim(df.pred)
df.pred
str(df.pred)
# table(is.na(df.pred))

################################################################################
### save the dataset 

# With write.csv(), very slow and heavy (88 seconds and 500 MB)
# system.time({
# write.csv(df.pred, file = "../output/dataset/dataset_22_covariates.csv", row.names = F)
# })

# With write_csv() and zip, fast but still heavy (2 seconds and 500 MB)
# library(readr)
# system.time({
# write_csv(df.pred, file = "../output/dataset/dataset_22_covariates.csv.gz")
# })
# 
# system.time({ # (40 seconds and 200 MB)
# zip(zipfile = "../output/dataset/dataset_22_covariates.zip", 
#     files = "../output/dataset/dataset_22_covariates.csv")
# })
# 
# system.time({ # read from zip: 18 seconds
# df.pred <- read_csv("../output/dataset/dataset_22_covariates.zip")
# })

### Or write and read directly as "csv.gz"
# system.time({ # 36 seconds, 228 MB
#   write_csv(df.pred, file = "../output/dataset/dataset_22_covariates.csv.gz")
# })
# system.time({ # read from csv.gz: 19 seconds
# df.pred <- read_csv("../output/dataset/dataset_22_covariates.csv.gz")
# })

# ### With R RDS fast and light
# system.time({ # 27 seconds and 120 MB and compress= TRUE
# saveRDS(df.pred, "../output/dataset/dataset_22_covariates.rds", compress= TRUE) # lighter than zip
# })
# 
# system.time({ # 2 seconds to read
# df.pred <- readRDS("../output/dataset/dataset_22_covariates.rds") # faster than read_csv
# })

### With "fst" faster and lighter
# install.packages("fst")
library(fst)
system.time({ # 5 seconds to write and 100 MB
write_fst(df.pred, "../output/dataset/dataset_22_covariates_df_pred_all_pixels.fst", compress = 100)
})

system.time({ # .5 second to read
df.pred <- read_fst("../output/dataset/dataset_22_covariates_df_pred_all_pixels.fst")
})

# Source: 
# https://data.nozav.org/post/2019-r-data-frame-benchmark/
# https://edomt.github.io/File-IO-Storage/ # faster and lighter
# https://waterdata.usgs.gov/blog/formats/ # very good
# https://community.rstudio.com/t/compres-data-to-export-it/48000/2

str(df.pred)

################################################################################
### Covariates correlation
library(corrplot)

df.pred$Geology <- as.integer(df.pred$Geology)
df.pred$Geomorfology <- as.integer(df.pred$Geomorfology)

names(df.pred)

df.pred_mt <- as.matrix(df.pred[,1:12])
M = cor(df.pred_mt)

# save correlation
# png(height=1800, width=1800, res= 200, 
    # file= "../output/graphs_plots/INP_covariates_correlation.png")

corrplot(M, order = 'original', addCoef.col = 'black', #tl.pos = 'd', 
         number.cex = .7, tl.cex = .8, number.font = 1, number.digits = 2,
         tl.srt = 45, 
         title = "INP covariates correlation",
         type = "full", cl.pos = 'n', col = COL2('RdYlBu'),
         mar = c(0, 0, 2, 0)
         )
# dev.off()

library(caret)
hc = findCorrelation(M, cutoff= 0.45)
hc = sort(hc)
reduced_Data = df.pred[, -c(hc)]
names(reduced_Data)
# print(reduced_Data)



################################################################################
### Extract data for spectral predition
################################################################################

stk_subsuf_img <- terra::rast("../output/raster/stk_subsuf_img_10m.tif")

### With Terra package
# upper_points <- as.points(sp_points, values=F, na.rm=TRUE, na.all=T) # convert pixels in points terra
# writeVector(upper_points, "../output/shp/upper_points.shp")

system.time({df.pred_spectr_subsuf = terra::extract(stk_subsuf_img,
                                                    pni_alta,
                                                    method= "simple",
                                                    xy= T,
                                                    ID= F)})

dim(df.pred_spectr_subsuf)
df.pred_spectr_subsuf
str(df.pred_spectr_subsuf)
# table(is.na(df.pred_spectr_subsuf))

################################################################################
### save the dataset with "fst"
library(fst)
system.time({ 
  write_fst(df.pred_spectr_subsuf, "../output/dataset/dataset_df.pred_spectr_subsuf.fst", compress = 100)
})

# system.time({ 
#   df.pred_spectr_subsuf <- read_fst("../output/dataset/dataset_df.pred_spectr_subsuf.fst")
# })


################################################################################
### Preparation data to extract points information form the stack rasters (covs)
################################################################################

################################################################################
### Load wet chemistry, spectral and spacial data ###
#' ### Auto setwd
library(rstudioapi); current_path <- getActiveDocumentContext()$path; setwd(dirname(current_path)) ### https://eranraviv.com/r-tips-and-tricks-working-directory/
getwd()
# gc(); rm(list=ls())
ld_spc <- read.csv("../input/dataset/dados_PNI_sc_spctr_data.csv", sep = ";", dec = ".", skip = 0)

### indicação de linhas ld_spc (lab data + spectral)
### de 1 a 335: dados de solo + dados espectrais espectrorrdiômetro ASD
### de 336 a 359: dados Paula Soares, pode ser usado na validação de mapas
### de 360 a 376: Rock outcrops
### de 377 a 431: pontos do google que Elias utilizou para validação de classes de solo
### de 432 a 476: calibração com spectralon (white reference)
### de 477 a 479: calibração absoluta com spectralon da caixa de madeira
### de 480 a 482: espectros PNI lido em vazio (pni00164.asd sem amostra de solo), e pni00167.asd e pni00168.asd não foi feita analise quimica das amostra de solo
### Observação: Elias quando foi devolver (unificar) amostras PNI utilizadas por Yuri para a embalagem principal verificou a possibilidade de ter errado 
### nas etiquetas do P77. Ficando então a sugestão: 		
### pni00216.asd	P77 	O1 ; pni00219.asd	P77 	O2 ; pni00220.asd	P77 	OB ; pni00217.asd	P77 	Bi1 ; pni00218.asd	P77 	Bi2
### sugestão minha, talvez remover o P77 das análises espectrais

### removendo linhas sem espectro ### http://www.datasciencemadesimple.com/drop-variables-columns-r-using-dplyr/
ld_spc <- ld_spc[-c(3, 10, 24, 35, 39, 56, 67, 69, 143, 145, 208, 226, 261, 327), ]

### Ajustando dados de solo e espectros PNI ###
# ld_spc[1:4, 1:45]
# names(ld_spc[1:45])
# names(ld_spc[1200:NCOL(ld_spc)])

dim(ld_spc)
### Selecting the surface horizon only
library(dplyr)
ld_spc <- ld_spc %>% filter(top == "0")

### removendo P77
# ld_spc <- ld_spc %>% filter(id != "P77")

### Look on data
# dim(ld_spc)
# head(names(ld_spc), 50)
# tail(names(ld_spc), 50)
# ld_spc$name
# ld_spc[, 5:6]
# ld_spc[1:6, 1:32]
# ld_spc[, c(1:4, 23)]

################################################################################
### Prepare data to run the model
################################################################################
ld_spc[1:84, c(1:4, 23)] ### Spectral data, model data (Elias, Yuri)
ld_spc[85:90, c(1:4, 23)] ### No spectral data, external validation (Paula Soares)

### Remove high correlated covariables (spectres)
# ld_spc[1:84, 42:43]
# ld_spc[1:84, 2190:ncol(ld_spc)]
# ld_spc[1:84, 42:ncol(ld_spc)]
library(caret)
spctral_data <- ld_spc[1:84, (42+80):(ncol(ld_spc)-80)]
m2 = cor(spctral_data)
hc = findCorrelation(m2, cutoff=0.99998)
hc = sort(hc)
reduced_Data = spctral_data[,-c(hc)]
names(reduced_Data)
print(reduced_Data)
lab_spectral <- cbind(ld_spc[1:84, 1:41], reduced_Data)

### Shuffle data (randomly)
set.seed(1234)
shuf_lab_spectral = lab_spectral[sample(1:nrow(lab_spectral)), ]
shuf_lab_spectral[, c(1:4, 23)]
dim(shuf_lab_spectral)
# write_fst(shuf_lab_spectral, "../output/dataset/shuf_lab_spectral.fst", compress = 100)

### Random sample for external validation 
set.seed(1234)
s1 <- shuf_lab_spectral[, c(1:4, 23)] %>% sample_frac(0.14)

modelData <- shuf_lab_spectral[!shuf_lab_spectral$id.num %in% s1$id.num, ]
dim(modelData)
modelData[, c(1:4, 23)]
# write_fst(modelData, "../output/dataset/modelData.fst", compress = 100)

ext_valData_spectral <- shuf_lab_spectral[shuf_lab_spectral$id.num %in% s1$id.num, ]
dim(ext_valData_spectral)
ext_valData_spectral[, c(1:4, 23)]
# write_fst(ext_valData_spectral, "../output/dataset/ext_valData_spectral.fst", compress = 100)

### Add the random selection with the external data
ext_valData_soil_prop <- rbind(ext_valData_spectral[, 1:30], ld_spc[85:90, 1:30])
dim(ext_valData_soil_prop)
# write_fst(ext_valData_soil_prop, "../output/dataset/ext_valData_soil_prop.fst", compress = 100)

################################################################################
### Extract points information from the stack rasters (covariates) for (df.mod)
################################################################################

################################### Extract
# library(raster)
### Raster package
# ### Load spectral data; ### convert plain data in a spatial object (SpatialPointsDataFrame)
# coordinates(modelData)=~X+Y # tranform data frame to spatialpointdataframe (transformando os dados da tabela em shape de pontos, 
# apos transformar x e y minusculo somem)
# # definindo a projeção    
# # CRS("+init=epsg:4326") # obtendo o codigo WGS84 lat long
# CRS("+init=epsg:32723") # obtendo o codigo WGS84 utm
# wgs84_utm = CRS("+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"); # (colocando sistema de referencia no shape de pontos)
# proj4string(modelData) = wgs84_utm # To define the projection
# # plot(modelData)
# # class(modelData)
# # names(modelData)
# # modelData
# 
# # library(snow)
# # beginCluster() # use all cores
# # extraindo os dados das covariáveis usando os dados de solos
# a = raster::extract(raster(stk_dsm), modelData, method="simple", sp=T); df.mod = as.data.frame(a)
# # endCluster()
# names(df.mod)

library(terra)
### Terra package
################################################################################
### Convert plain data in a spatial object (SpatVector)
### df.model DSM
crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
points <- vect(modelData, geom=c("X", "Y"), crs= crs_ref, keepgeom=FALSE)
plot(points)

### Extract
system.time({df.mod = terra::extract(stk_dsm, 
                                     points, 
                                     method= "simple", 
                                     xy= T, 
                                     ID= F,
                                     bind= T # return a SpatVector object, with all variables
                                     )})

df.mod
dim(df.mod)
df.mod <- as.data.frame(df.mod) # from SpatVector to df
# str(df.mod)

library(fst)
system.time({ #
  write_fst(df.mod, "../output/dataset/dataset_HDSM_covariates_df_mod_point_pixels.fst", compress = 100)
})

# system.time({ #
#   df.mod <- read_fst("../output/dataset/dataset_HDSM_covariates_df_mod_point_pixels.fst")
# })

class(df.mod)
dim(df.mod)
str(df.mod)
names(df.mod)
dput(names(df.mod))
table(is.na(df.mod))
table(is.na(df.mod[40:181]))

############################################################
### df.external.validation
### Convert plain data in a spatial object (SpatVector)
crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
points <- vect(ext_valData_soil_prop, geom=c("X", "Y"), crs= crs_ref, keepgeom=FALSE)
plot(points)

### Extract
system.time({df.ex.val = terra::extract(stk_dsm, 
                                        points, 
                                        method= "simple", 
                                        xy= T, 
                                        ID= F,
                                        bind= T # return a SpatVector object, with all variables
)})

df.ex.val
dim(df.ex.val)
df.ex.val <- as.data.frame(df.ex.val) # from SpatVector to df
str(df.ex.val)

library(fst)
system.time({ #
  write_fst(df.ex.val, "../output/dataset/dataset_Ext_val_covariates_df_external_validation.fst", compress = 100)
})

# system.time({ #
#   df.ex.val <- read_fst("../output/dataset/dataset_Ext_val_covariates_df_external_validation.fst")
# })

################################################################################
### For spectral image
################################################################################
### Extract df.mod 
### Convert plain data in a spatial object (SpatVector)
points <- vect(modelData, geom=c("X", "Y"), crs= crs_ref, keepgeom=FALSE)
plot(points)

system.time({df.mod_subsuf_img = terra::extract(stk_subsuf_img, 
                                                points, 
                                                method= "simple", 
                                                xy= T, 
                                                ID= F,
                                                bind= T # return a SpatVector object, with all variables
)})

df.mod_subsuf_img
dim(df.mod_subsuf_img)
df.mod_subsuf_img <- as.data.frame(df.mod_subsuf_img) # from SpatVector to df
# str(df.mod_subsuf_img)

library(fst)
system.time({ #
  write_fst(df.mod_subsuf_img, "../output/dataset/dataset_spectral_image_covariates_df_mod_point_pixels.fst", compress = 100)
})

# system.time({ #
#   df.mod_subsuf_img <- read_fst("../output/dataset/dataset_spectral_image_covariates_df_mod_point_pixels.fst")
# })

############################################################
### ext_valData_spectral
### Convert plain data in a spatial object (SpatVector)
crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
points <- vect(ext_valData_spectral, geom=c("X", "Y"), crs= crs_ref, keepgeom=FALSE)
plot(points)

### Extract
system.time({ext_valData_spectral = terra::extract(stk_subsuf_img, 
                                                   points, 
                                                   method= "simple", 
                                                   xy= T, 
                                                   ID= F,
                                                   bind= T # return a SpatVector object, with all variables
)})

ext_valData_spectral
dim(ext_valData_spectral)
ext_valData_spectral <- as.data.frame(ext_valData_spectral) # from SpatVector to df
str(ext_valData_spectral)

library(fst)
system.time({ #
  write_fst(ext_valData_spectral, "../output/dataset/dataset_Ext_val_spectral_covariates.fst", compress = 100)
})

# system.time({ #
#   ext_valData_spectral <- read_fst("../output/dataset/dataset_Ext_val_spectral_covariates.fst")
# })


################################################################################
# ### Several trials with grid to find best tune for RF
# library(tidymodels)
# library(ranger)
# # install.packages("ranger")
# # install.packages("themis")
# library(themis)
# 
# # Prepare data
# test_models <- df.mod[, c(21, 170:(191-10))] # spatial covs, DEM, Slope, etc
# test_models <- df.mod[, c(21, 40:169)] # spctral covs
# 
# set.seed(1234)
# # Split the data 70/30 for training/testing
# train_test_split <- initial_split(test_models, 0.7, strata = 'C')
# 
# # Extract training data
# train <- training(train_test_split)
# dim(train)
# 
# # Create 5-fold 3-repeat CV splits
# cv_folds <- vfold_cv(train, v = 5, repeats = 3, strata = 'C')
# 
# # Simple recipe with ROSE
# rec <- recipe(C ~ ., data = train) #%>% step_rose(C)
# 
# # Specify model
# model <- rand_forest(
#   mode = 'regression',
#   # mode = 'classification',
#   engine = 'randomForest',
#   mtry = tune(), # number of predictors that will be randomly sampled at each split (mtry)
#   trees = tune(), # number of trees contained in the ensemble (ntree)
#   min_n = tune() # minimum number of data points in a node that are required for the node to be split further (nodesize)
# )
# 
# # Define workflow
# wflow <- workflow(
#   preprocessor = rec,
#   spec = model
# )
# 
# # Perform CV
# res <-
#   tune_grid(
#     object = wflow,
#     resamples = cv_folds,
#     control = control_grid(save_pred = TRUE)
#   )
# 
# # Find the best hyperparameters for the randomForest model
# best_hyperparams <- select_best(res)
# 
# # Add those hyperparams to the workflow
# final_wflow <- finalize_workflow(wflow, best_hyperparams)
# 
# # Fit the updated workflow to the whole train data
# # Evaluate performance in the held-out test data
# final_fit <- last_fit(final_wflow, train_test_split)

################################################################################
# ??tuneRF
# data(fgl, package="MASS")
# fgl.res <- tuneRF(fgl[,-10], fgl[,10], stepFactor=1.5)
# 
# test_models <- df.mod[, c(21, 170:(191-10))] # spatial covs, DEM, Slope, etc
# fgl.res <- tuneRF(test_models[,-1], test_models[,1], stepFactor=1.5)

################################################################################
"
DDDDDDDDDDDDD           SSSSSSSSSSSSSSS MMMMMMMM               MMMMMMMM
D::::::::::::DDD      SS:::::::::::::::SM:::::::M             M:::::::M
D:::::::::::::::DD   S:::::SSSSSS::::::SM::::::::M           M::::::::M
DDD:::::DDDDD:::::D  S:::::S     SSSSSSSM:::::::::M         M:::::::::M
D:::::D    D:::::D S:::::S            M::::::::::M       M::::::::::M
D:::::D     D:::::DS:::::S            M:::::::::::M     M:::::::::::M
D:::::D     D:::::D S::::SSSS         M:::::::M::::M   M::::M:::::::M
D:::::D     D:::::D  SS::::::SSSSS    M::::::M M::::M M::::M M::::::M
D:::::D     D:::::D    SSS::::::::SS  M::::::M  M::::M::::M  M::::::M
D:::::D     D:::::D       SSSSSS::::S M::::::M   M:::::::M   M::::::M
D:::::D     D:::::D            S:::::SM::::::M    M:::::M    M::::::M
D:::::D    D:::::D             S:::::SM::::::M     MMMMM     M::::::M
DDD:::::DDDDD:::::D  SSSSSSS     S:::::SM::::::M               M::::::M
D:::::::::::::::DD   S::::::SSSSSS:::::SM::::::M               M::::::M
D::::::::::::DDD     S:::::::::::::::SS M::::::M               M::::::M
DDDDDDDDDDDDD         SSSSSSSSSSSSSSS   MMMMMMMM               MMMMMMMM
"
################################################################################
### "Hand made" K-fold cross validation
#' ### Auto setwd
library(rstudioapi); current_path <- getActiveDocumentContext()$path; setwd(dirname(current_path)) ### https://eranraviv.com/r-tips-and-tricks-working-directory/
getwd()
# gc(); rm(list=ls())

{
library(terra)
library(fst)
library(randomForest)
library(readr)
# install.packages("devtools")
# library(devtools)
# install_bitbucket("brendo1001/ithir/pkg") 
library(ithir)
library(hexbin)
library(grid)

### Load data
# shuf_lab_spectral <- read_fst("../output/dataset/shuf_lab_spectral.fst")
# modelData <- read_fst("../output/dataset/modelData.fst")
# ext_valData_spectral <- read_fst("../output/dataset/ext_valData_spectral.fst")
# ext_valData_soil_prop <- read_fst("../output/dataset/ext_valData_soil_prop.fst")

stk_dsm <- terra::rast("../output/raster/stk_dsm_10m.tif")

### old names
# df.mod <- read_fst("../output/dataset/dataset_22_covariates_df_mod_point_pixels.fst")
# df.ex.val <- read_fst("../output/dataset/dataset_22_covariates_df_external_validation.fst")
# df.pred <- read_fst("../output/dataset/dataset_22_covariates_df_pred_all_pixels.fst")

### Adjusted dataset names
df.mod <- read_fst("../output/dataset/dataset_HDSM_covariates_df_mod_point_pixels.fst")
df.ex.val <- read_fst("../output/dataset/dataset_Ext_val_covariates_df_external_validation.fst")
df.pred <- read_fst("../output/dataset/dataset_22_covariates_df_pred_all_pixels.fst")
}

# df.mod$Geology <- as.factor(df.mod$Geology)
# df.mod$Geomorfology <- as.factor(df.mod$Geomorfology)

# str(df.mod[, 160:193-10])
# str(df.ex.val)
# str(df.pred)

### Look on data
# df.mod_spatial_data <- df.mod[, c(21, 170:(191-10), 192:193-10)] # spatial covs, DEM, Slope, etc
# df.mod_spectral_data <- df.mod[, c(21, 40:169-10, 192:193-10)] # spctral covs

### Check with indexes
# df.mod_spatial_data$idnum <- 1:nrow(df.mod_spatial_data)

### Split data into 10 folds
folds <- cut(seq(1,nrow(df.mod)), breaks= 8, labels=FALSE)
nflds <- length(unique(folds))
# table(folds)

### Test to see the behavior of k-fold
# for(i in 1:nflds){ ### This loop is to run and validate with the K-folds
#   ### Segment your data by fold using the which() function 
#   testIndexes <- which(folds== i, arr.ind=TRUE)
#   testData <- df.mod_spatial_data[testIndexes, ]
#   cat("Test data: ")
#   print(testData$idnum)
#   print(" ")
#   
#   trainData <- df.mod_spatial_data[-testIndexes, ]
#   cat("Train data: ")
#   print(trainData$idnum)
#   print(" ")
#   Sys.sleep(1)
# }
### Test OK, remove indexes
# df.mod_spatial_data$idnum <- NULL

### list to store the predicted raster
# predicted_raster <- list()

{
### create null objects (lists) to fill with the results: rf
rf_models <- list() ### results from the models

fitted_value_rf_cv <- list() ### results from the models (for cv)
goof_rf_cv <- list() ### metrics
goof_rf_cv_df <- list() ### metrics
# goof_rf_mean_cv_df <- list() ### metrics

fitted_value_rf_ext_val <- list() ### results from the models (for external validation)
goof_rf_ext_val <- list() ### metrics
goof_rf_ext_val_df <- list() ### metrics
# goof_rf_mean_ext_val_df <- list() ### metrics

fitted_value_rf_calibr <- list()  ### results from the models (for calibration)
goof_rf_calibr <- list()
goof_rf_calibr_df <- list()
# goof_rf_mean_calibr_df <- list()
}

system.time(
for(j in c(21)){ ### C=21
  ### building the formula to use in the algorithm
  # allVars <- colnames(df.mod[, c(j, 170:191)]) # trainData[, c(j, 170:191)] select the terrain and images covariates
  allVars <- colnames(df.mod[, c(j, 170:(191-10))]) # trainData[, c(j, 170:(191-10))] select the terrain and images covariates
  sp_spc <- df.mod[, c(j, 170:(191-10))]
  predictorVars <- allVars[!allVars%in%names(sp_spc[1])]
  predictorVars <- paste(predictorVars, collapse = "+")
  form = as.formula(paste(names(sp_spc[1]), "~", predictorVars, collapse = "+"))
  
  for(i in 1:nflds){ ### This loop is to run and validate with the K-folds
    ### Segment your data by fold using the which() function 
    testIndexes <- which(folds== i, arr.ind=TRUE)
    # testData <- df.mod_spatial_data[testIndexes, ]
    # trainData <- df.mod_spatial_data[-testIndexes, ]
    testData <- df.mod[testIndexes, ]
    trainData <- df.mod[-testIndexes, ]
    
    ######################################
    #### rf_model (trainig data)
    set.seed(1234) # to fix the model randomness
    rf_models[[i]] <- randomForest(form, data= trainData, ntree= 150, mtry= 3, nodesize= 3, importance= T) 
    
    ################ Calibration
    ### predict with the model (test data), and store in "i", with means in the loop of the folders 
    fitted_value_rf_calibr[[i]] <- predict(rf_models[[i]], newdata= trainData[ , 170:(191-10)]) 
    
    ### perform goof for calibr and store as list
    goof_rf_calibr[[i]] <- goof(trainData[,j], fitted_value_rf_calibr[[i]], type='DSM') 
    goof_rf_calibr[[i]]$model_number <- i # paste0("Model ", i)
    
    ### extrated the values from the list inside each property group and store as list object
    goof_rf_calibr_df[[j]] <- (as.data.frame(do.call(rbind, goof_rf_calibr)))
    
    ### compute mean for each metric of all properties across all folders store as list object
    # goof_rf_mean_calibr_df[[j]] <- colMeans(x= goof_rf_calibr_df[[j]])
    
    ### pass to a dataframe
    # df_goof_rf_mean_calibr <- as.data.frame(t(goof_rf_mean_calibr_df[[j]]))
    # df_goof_rf_mean_calibr$model_number <- NULL
    
    ### extrated the values from the list and store as data.frame object
    # goof_rf_mean_calibr_df <- (as.data.frame(do.call(rbind, goof_rf_mean_calibr_df)))
    # goof_rf_mean_calibr_df$propertie <- names(sp_spc[1])
    # goof_rf_mean_calibr_df$model_number <- NULL
    
    ### compute mean for each metric of all properties across all folders store as list object
    ### and bind with folders metrics
    goof_rf_calibr_metrics <- as.data.frame(rbind(goof_rf_calibr_df[[j]], colMeans(x= goof_rf_calibr_df[[j]])))
    
    ### fill "mean" model with zero
    if (!is.na(goof_rf_calibr_metrics$model_number[9])) {
      goof_rf_calibr_metrics$model_number[9] <- 0
    }
    
    ### Write results on disk
    ds <- "DSM"
    # write_csv(round(goof_rf_calibr_df[[j]], 2),
    #           file = paste0("../output/dataset/", ds, "/", ds, "_",
    #                         names(sp_spc[1]), "_goof_rf_folds_calibr.csv"))
    # write_csv(round(df_goof_rf_mean_calibr, 2),
    #           file = paste0("../output/dataset/", ds, "/", ds, "_",
    #                         names(sp_spc[1]), "_goof_rf_mean_calibra.csv"))
    write_csv(round(goof_rf_calibr_metrics, 2),
              file = paste0("../output/dataset/", ds, "/", ds, "_",
                            names(sp_spc[1]), "_goof_rf_calibr_metrics.csv"))
    
    ################ k-fold
    ### predict with the model (test data), and store in "i", with means in the loop of the folders 
    fitted_value_rf_cv[[i]] <- predict(rf_models[[i]], newdata= testData[ , 170:(191-10)]) 
    
    ### perform goof for cv and store as list
    goof_rf_cv[[i]] <- goof(testData[,j], fitted_value_rf_cv[[i]], type='DSM') 
    goof_rf_cv[[i]]$model_number <- i # paste0("Model ", i)
    
    ### extrated the values from the list inside each property group and store as list object
    goof_rf_cv_df[[j]] <- (as.data.frame(do.call(rbind, goof_rf_cv)))
    
    ### compute mean for each metric of all properties across all folders store as list object
    # goof_rf_mean_cv_df[[j]] <- colMeans(x= goof_rf_cv_df[[j]])
    
    ### pass to a dataframe
    # df_goof_rf_mean_cv <- as.data.frame(t(goof_rf_mean_cv_df[[j]]))
    # df_goof_rf_mean_cv$model_number <- NULL
    
    ### extrated the values from the list and store as data.frame object
    # goof_rf_mean_cv_df <- (as.data.frame(do.call(rbind, goof_rf_mean_cv_df)))
    # goof_rf_mean_cv_df$propertie <- names(sp_spc[1])
    # goof_rf_mean_cv_df$model_number <- NULL
    
    ### compute mean for each metric of all properties across all folders store as list object
    ### and bind with folders metrics
    goof_rf_k_cv_metrics <- as.data.frame(rbind(goof_rf_cv_df[[j]], colMeans(x= goof_rf_cv_df[[j]])))
    
    ### fill "mean" model with zero
    if (!is.na(goof_rf_k_cv_metrics$model_number[9])) {
      goof_rf_k_cv_metrics$model_number[9] <- 0
    }
    
    ### Write results on disk
    # write_csv(round(goof_rf_cv_df[[j]], 2),
    #           file = paste0("../output/dataset/", ds, "/", ds, "_",
    #                         names(sp_spc[1]), "_goof_rf_folds_k_fold_cv.csv"))
    # write_csv(round(df_goof_rf_mean_cv, 2),
    #           file = paste0("../output/dataset/", ds, "/", ds, "_",
    #                         names(sp_spc[1]), "_goof_rf_mean_k_fold_cv.csv"))
    write_csv(round(goof_rf_k_cv_metrics, 2),
              file = paste0("../output/dataset/", ds, "/", ds, "_",
                            names(sp_spc[1]), "_goof_rf_k_fold_cv_metrics.csv"))
    
    names(df.ex.val)
    ################ External validation
    ### predict with the model (test data), and store in "i", with means in the loop of the folders 
    fitted_value_rf_ext_val[[i]] <- predict(rf_models[[i]], newdata= df.ex.val[ , 29:(50-10)]) 
    
    ### perform goof for external validation and store as list
    goof_rf_ext_val[[i]] <- goof(df.ex.val[,j], fitted_value_rf_ext_val[[i]], type='DSM') 
    goof_rf_ext_val[[i]]$model_number <- i # paste0("Model ", i)
    
    ### extrated the values from the list inside each property group and store as list object
    goof_rf_ext_val_df[[j]] <- (as.data.frame(do.call(rbind, goof_rf_ext_val)))
    
    ### compute mean for each metric of all properties across all folders store as list object
    # goof_rf_mean_ext_val_df[[j]] <- colMeans(x= goof_rf_ext_val_df[[j]])
    
    ### pass to a dataframe
    # df_goof_rf_mean_ext_val <- as.data.frame(t(goof_rf_mean_ext_val_df[[j]]))
    # df_goof_rf_mean_ext_val$model_number <- NULL
    
    ### extrated the values from the list and store as data.frame object
    # goof_rf_mean_ext_val_df <- (as.data.frame(do.call(rbind, goof_rf_mean_ext_val_df)))
    # goof_rf_mean_ext_val_df$propertie <- names(sp_spc[1])
    # goof_rf_mean_ext_val_df$model_number <- NULL
    
    ### compute mean for each metric of all properties across all folders store as list object
    ### and bind with folders metrics
    goof_rf_ext_val_metrics <- as.data.frame(rbind(goof_rf_ext_val_df[[j]], colMeans(x= goof_rf_ext_val_df[[j]])))
    
    ### fill "mean" model with zero
    if (!is.na(goof_rf_ext_val_metrics$model_number[9])) {
      goof_rf_ext_val_metrics$model_number[9] <- 0
    }
    
    ### Write results on disk
    # write_csv(round(goof_rf_ext_val_df[[j]], 2),
    #           file = paste0("../output/dataset/", ds, "/", ds, "_",
    #                         names(sp_spc[1]), "_goof_rf_folds_ext_val.csv"))
    # write_csv(round(df_goof_rf_mean_ext_val, 2),
    #           file = paste0("../output/dataset/", ds, "/", ds, "_",
    #                         names(sp_spc[1]), "_goof_rf_mean_ext_val.csv"))
    write_csv(round(goof_rf_ext_val_metrics, 2),
              file = paste0("../output/dataset/", ds, "/", ds, "_",
                            names(sp_spc[1]), "_goof_rf_mean_ext_val.csv"))
    ### Save outputs
    w_list <- c("rf_models", 
                ### Calibration
                "fitted_value_rf_calibr",
                "goof_rf_calibr",
                "goof_rf_calibr_df",
                # "goof_rf_mean_calibr_df",
                # "df_goof_rf_mean_calibr",
                "goof_rf_calibr_metrics",
                ### k-fold cv
                "fitted_value_rf_cv",
                "goof_rf_cv",
                "goof_rf_cv_df",
                # "df_goof_rf_mean_cv",
                # "goof_rf_mean_cv_df",
                "goof_rf_k_cv_metrics",
                ### ext val
                "fitted_value_rf_ext_val",
                "goof_rf_ext_val",
                "goof_rf_ext_val_df",
                # "goof_rf_mean_ext_val_df",
                "goof_rf_ext_val_metrics")
    
    # w_list[1] # write lists
    # w_list[19] # write lists
    
    ### Save in a loop
    for (sv in seq_along(w_list)) {
      saveRDS(eval(parse(text= w_list[[sv]])),
              paste0("../output/dataset/", ds, "/", ds, "_",
                     names(sp_spc[1]), "_", w_list[[sv]], ".rds"), compress= TRUE)
    }
    
  }
}
)
### Liberate memory (except)
# rm(list=setdiff(ls(), c("rf_models", "stk_dsm", "df.mod", "df.ex.val", "df.pred", "folds", "nflds")))

# ### Read outputs
# ds <- "DSM"
# propertie <- "C"
# for (sv in seq_along(w_list)) {
#   assign(w_list[[sv]], readRDS(
#     paste0("../output/dataset/", ds, "/", ds, "_", 
#            propertie, "_", w_list[[sv]],".rds"))) 
# }

### Or one by one
# ds <- "DSM"
# propertie <- "C"
# rf_models <- readRDS(paste0("../output/dataset/", ds, "/", ds, "_", propertie, "_rf_models.rds")) 
# fitted_value_rf_cv <- readRDS(paste0("../output/dataset/", ds, "/", ds, "_", propertie, "_fitted_value_rf_cv.rds")) 
# fitted_value_rf_ext_val <- readRDS(paste0("../output/dataset/", ds, "/", ds, "_", propertie, "_fitted_value_rf_ext_val.rds")) 
# goof_rf_cv <- readRDS(paste0("../output/dataset/", ds, "/", ds, "_", propertie, "_goof_rf_cv.rds")) 
# goof_rf_ext_val <- readRDS(paste0("../output/dataset/", ds, "/", ds, "_", propertie, "_goof_rf_ext_val.rds")) 
# goof_rf_cv_df <- readRDS(paste0("../output/dataset/", ds, "/", ds, "_", propertie, "_goof_rf_cv_df.rds")) 
# goof_rf_ext_val_df <- readRDS(paste0("../output/dataset/", ds, "/", ds, "_", propertie, "_goof_rf_ext_val_df.rds")) 


################ Select best model to predict the map
# goof_rf_cv_df[[j]][which(goof_rf_cv_df[[j]] == max(goof_rf_cv_df[[j]]$R2)), "model_number"]

# best_model_kf_cv <- goof_rf_cv_df[[j]][goof_rf_cv_df[[j]]$R2 == max(goof_rf_cv_df[[j]]$R2), "model_number"]
# coef_best_model_kf_cv <- goof_rf_cv_df[[j]][goof_rf_cv_df[[j]]$R2 == max(goof_rf_cv_df[[j]]$R2), ]
# 
# best_model_val_ext <- goof_rf_ext_val_df[[j]][goof_rf_ext_val_df[[j]]$R2 == max(goof_rf_ext_val_df[[j]]$R2), "model_number"]
# coef_best_model_val_ext <- goof_rf_ext_val_df[[j]][goof_rf_ext_val_df[[j]]$R2 == max(goof_rf_ext_val_df[[j]]$R2), ]

# best_model_calibr <- goof_rf_calibr_df[[j]][goof_rf_calibr_df[[j]]$MSE == min(goof_rf_calibr_df[[j]]$MSE), "model_number"]
# coef_best_model_calibr <- goof_rf_calibr_df[[j]][goof_rf_calibr_df[[j]]$MSE == min(goof_rf_calibr_df[[j]]$MSE), ]
# 
# best_model_kf_cv <- goof_rf_cv_df[[j]][goof_rf_cv_df[[j]]$MSE == min(goof_rf_cv_df[[j]]$MSE), "model_number"]
# coef_best_model_kf_cv <- goof_rf_cv_df[[j]][goof_rf_cv_df[[j]]$MSE == min(goof_rf_cv_df[[j]]$MSE), ]
# 
# best_model_val_ext <- goof_rf_ext_val_df[[j]][goof_rf_ext_val_df[[j]]$MSE == min(goof_rf_ext_val_df[[j]]$MSE), "model_number"]
# coef_best_model_val_ext <- goof_rf_ext_val_df[[j]][goof_rf_ext_val_df[[j]]$MSE == min(goof_rf_ext_val_df[[j]]$MSE), ]

fst_qt_cal <- goof_rf_calibr_df[[j]][quantile(goof_rf_calibr_df[[j]]$RMSE, 0.25) > goof_rf_calibr_df[[j]]$RMSE, ]
best_model_calibr <- fst_qt_cal[fst_qt_cal$R2 == max(fst_qt_cal$R2), "model_number"]
coef_best_model_calibr <- fst_qt_cal[fst_qt_cal$R2 == max(fst_qt_cal$R2), ]

fst_qt_kcv <- goof_rf_cv_df[[j]][quantile(goof_rf_cv_df[[j]]$RMSE, 0.25) > goof_rf_cv_df[[j]]$RMSE, ]
best_model_kf_cv <- fst_qt_kcv[fst_qt_kcv$R2 == max(fst_qt_kcv$R2), "model_number"]
coef_best_model_kf_cv <- fst_qt_kcv[fst_qt_kcv$R2 == max(fst_qt_kcv$R2), ]

fst_qt_ex_val <- goof_rf_ext_val_df[[j]][quantile(goof_rf_ext_val_df[[j]]$RMSE, 0.25) > goof_rf_ext_val_df[[j]]$RMSE, ]
best_model_val_ext <- fst_qt_ex_val[fst_qt_ex_val$R2 == max(fst_qt_ex_val$R2), "model_number"]
coef_best_model_val_ext <- fst_qt_ex_val[fst_qt_ex_val$R2 == max(fst_qt_ex_val$R2), ]

### Predict the map pixels with best model
pred_map <- stats::predict(rf_models[[best_model_kf_cv]], df.pred[, 1:(22-10)]) # to map
# pred_map <- stats::predict(rf_models[[best_model_val_ext]], df.pred[, 1:(22-10)]) # to map

### Join predicted with X, Y
Soil_prop_to_map = cbind(pred_map, df.pred[, 23:24-10]) # df.pred[, 23:24] = X, Y

### Convert to points
crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
points <- terra::vect(Soil_prop_to_map, geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)

### Convert points to raster
map <- terra::rasterize(points, stk_dsm[[1]] , "pred_map")
# plot(map)

# Export each raster
terra::writeRaster(map, filename= paste0("../output/raster/predict_rasters_results/", ds, "/", ds, "_predict_", names(sp_spc[1]), ".tif"), overwrite=TRUE)
### limpando memoria durante
# rm(list=setdiff(ls(), c("rf_model", "df.mod", "form", "raster_stack_22_bands", "predicted_raster", "tipo", "sp_spc")))

################################################################################
### Applying the best model to the entire data
### Calibration global
pred_df.mod <- stats::predict(rf_models[[best_model_kf_cv]], df.mod[, 170:(191-10)]) # to map

### perform goof for entire data
### Calibration global
goof_rf_global <- goof(df.mod[, j], pred_df.mod, type='DSM') 
goof_rf_global

### Values from inside the model, do not be used
### goof(rf_models[[i]]$predicted, rf_models[[i]]$y)

################################################################################
### Plots

## ----correlation_plot, message=F-----------------------------------------
library(grid)
library(ithir)
library(RColorBrewer)
library(hexbin)


################################################################################
### Calibration (average across folders)
################################################################################
crp <- colorRampPalette(rev(brewer.pal(9,'YlOrBr')))

### Unidade
### os dados PNI são todo em bquote((cmol[c]~dm^-3)) e NÂO em bquote((cmol[c]~kg^-1))
# unit <- bquote((cmol[c]~dm^-3)) ### Na, K, Ca, Mg, Al, H_AL []subescrito, ^sobreescreito, ~espaço 
unit <- "(%)" ### C, H, N ### quando C convertido para g/kg, unit = bquote((g[]~kg^-1))
# unit <- bquote((g~dm^-3)) ### BD

############### Bring the data from the best model ############################
testIndexes <- which(folds== best_model_calibr, arr.ind=TRUE)
testData <- df.mod[testIndexes, ]
trainData <- df.mod[-testIndexes, ]

### save plot RF calibration
plt_val_rf <- as.data.frame(cbind(trainData[,j], fitted_value_rf_calibr[[best_model_calibr]]))

rf_path <- file.path(paste0("../output/graphs_plots/", ds, "/", ds, "_",  names(sp_spc[1]), "_Scater_plot_Calibration_average_across_folers", ".png"))
png(file= rf_path, height=12, width=12, units="cm", res=120)
### run plot 
### hex <- hexbin(plt_val_rf[,1] ~ plt_val_rf[,2]) ### https://rdrr.io/cran/hexbin/man/hexbinplot.html 
bin <- hexbinplot(plt_val_rf[,1] ~ plt_val_rf[,2], data= plt_val_rf, colramp= crp, ### "plt_val_rf[[j]][,1]"= colums 1 from: element "j" in the list
                  type = c("g", "r"), style = "colorscale", border = T,
                  cex.title= .9, cex.labels = .9, xbins = 30, colorcut = seq(0, 1, length = 8), aspect = "1",
                  # main = c(paste("Rapid Eye + terrain + geographic")),
                  xlab = bquote(paste("Observed", " ", .(names(sp_spc[1])), " ", .(unit))),
                  ylab = bquote(paste("Predicted", " ", .(names(sp_spc[1])), " ", .(unit))),
                  panel = function(x, y, ...) { ### to add 1:1 line # https://stackoverflow.com/a/53515302/14361772
                    panel.hexbinplot(x, y, ...)
                    lattice::panel.abline(a = 0, b = 1, lty= 9, lwd = .5)})
p <- plot(bin) ### https://stackoverflow.com/questions/29958561/how-to-plot-points-on-hexbin-graph-in-r
grid.text(label = c(paste("R²", round((data.frame(coef_best_model_calibr))[,1], 2), "\n",      ## position == $R2
                          "MSE", round((data.frame(coef_best_model_calibr))[,3], 2), "\n",     ## position == $MSE
                          "RMSE", round((data.frame(coef_best_model_calibr))[,4], 2), "\n",    ## position == $RMSE
                          "bias", round((data.frame(coef_best_model_calibr))[,5], 2))),         ## position == $RPD
          x = unit(0.26, "npc"), y = unit(0.745, "npc"))
grid.text(label = c(paste("Calibration DSM \n Sentinel + Terrain + Geographic")), x = unit(0.477, "npc"), y = unit(0.9, "npc"))
dev.off()


################################################################################
### K-fold cross validation average coefficients across folders

############### Bring the data from the best model ############################
testIndexes <- which(folds== best_model_kf_cv, arr.ind=TRUE)
testData <- df.mod[testIndexes, ]
trainData <- df.mod[-testIndexes, ]

plt_val_rf <- as.data.frame(cbind(testData[,j], fitted_value_rf_cv[[best_model_kf_cv]]))

crp <- colorRampPalette(brewer.pal(9,'YlOrBr'))

rf_path <- file.path(paste0("../output/graphs_plots/", ds, "/", ds, "_",  names(sp_spc[1]), "_Scater_plot_K_fold_Validation_average_across_folers", ".png"))
png(file= rf_path, height=12, width=12, units="cm", res=120)
### run plot 
### hex <- hexbin(plt_val_rf[,1] ~ plt_val_rf[,2]) ### https://rdrr.io/cran/hexbin/man/hexbinplot.html 
bin <- hexbinplot(plt_val_rf[,1] ~ plt_val_rf[,2], data= plt_val_rf, colramp= crp, ### "plt_val_rf[[j]][,1]"= colums 1 from: element "j" in the list
                  type = c("g", "r"), style = "colorscale", border = T,
                  cex.title= .9, cex.labels = .9, xbins = 30, colorcut = seq(0, 1, length = 8), aspect = "1",
                  # main = c(paste("Rapid Eye + terrain + geographic")),
                  xlab = bquote(paste("Observed", " ", .(names(sp_spc[1])), " ", .(unit))),
                  ylab = bquote(paste("Predicted", " ", .(names(sp_spc[1])), " ", .(unit))),
                  panel = function(x, y, ...) { ### to add 1:1 line # https://stackoverflow.com/a/53515302/14361772
                    panel.hexbinplot(x, y, ...)
                    lattice::panel.abline(a = 0, b = 1, lty= 9, lwd = .5)})
p <- plot(bin) ### https://stackoverflow.com/questions/29958561/how-to-plot-points-on-hexbin-graph-in-r
grid.text(label = c(paste("R²", round((data.frame(coef_best_model_kf_cv))[,1], 2), "\n",      ## position == $R2
                          "MSE", round((data.frame(coef_best_model_kf_cv))[,3], 2), "\n",     ## position == $MSE
                          "RMSE", round((data.frame(coef_best_model_kf_cv))[,4], 2), "\n",    ## position == $RMSE
                          "bias", round((data.frame(coef_best_model_kf_cv))[,5], 2))),         ## position == $RPD
          x = unit(0.26, "npc"), y = unit(0.745, "npc"))
grid.text(label = c(paste("K-fold cross validation DSM \n Sentinel + Terrain + Geographic")), x = unit(0.477, "npc"), y = unit(0.9, "npc"))
dev.off()

################################################################################
### External validation average coefficients across folders
################################################################################

############### Bring the data from the best model ############################
testIndexes <- which(folds== best_model_val_ext, arr.ind=TRUE)
testData <- df.mod[testIndexes, ]
trainData <- df.mod[-testIndexes, ]

plt_val_rf <- as.data.frame(cbind(df.ex.val[,j], fitted_value_rf_ext_val[[best_model_val_ext]]))

rf_path <- file.path(paste0("../output/graphs_plots/", ds, "/", ds, "_",  names(sp_spc[1]), "_Scater_plot_External_Validation_average_across_folers", ".png"))
png(file= rf_path, height=12, width=12, units="cm", res=120)
### run plot 
### hex <- hexbin(plt_val_rf[,1] ~ plt_val_rf[,2]) ### https://rdrr.io/cran/hexbin/man/hexbinplot.html 
bin <- hexbinplot(plt_val_rf[,1] ~ plt_val_rf[,2], data= plt_val_rf, colramp= crp, ### "plt_val_rf[[j]][,1]"= colums 1 from: element "j" in the list
                  type = c("g", "r"), style = "colorscale", border = T,
                  cex.title= .9, cex.labels = .9, xbins = 30, colorcut = seq(0, 1, length = 8), aspect = "1",
                  # main = c(paste("Rapid Eye + terrain + geographic")),
                  xlab = bquote(paste("Observed", " ", .(names(sp_spc[1])), " ", .(unit))),
                  ylab = bquote(paste("Predicted", " ", .(names(sp_spc[1])), " ", .(unit))),
                  panel = function(x, y, ...) { ### to add 1:1 line # https://stackoverflow.com/a/53515302/14361772
                    panel.hexbinplot(x, y, ...)
                    lattice::panel.abline(a = 0, b = 1, lty= 9, lwd = .5)})
p <- plot(bin) ### https://stackoverflow.com/questions/29958561/how-to-plot-points-on-hexbin-graph-in-r
grid.text(label = c(paste("R²", round((data.frame(coef_best_model_val_ext))[,1], 2), "\n",      ## position == $R2
                          "MSE", round((data.frame(coef_best_model_val_ext))[,3], 2), "\n",     ## position == $MSE
                          "RMSE", round((data.frame(coef_best_model_val_ext))[,4], 2), "\n",    ## position == $RMSE
                          "bias", round((data.frame(coef_best_model_val_ext))[,5], 2))),         ## position == $RPD
          x = unit(0.26, "npc"), y = unit(0.745, "npc"))
grid.text(label = c(paste("K-fold cross validation DSM \n Sentinel + Terrain + Geographic")), x = unit(0.477, "npc"), y = unit(0.9, "npc"))
dev.off()


################################################################################
### plot predicted soil map
################################################################################
# map <- terra::rast(paste0("../output/raster/predict_rasters_results/", ds, "/", ds, "_predict_", names(sp_spc[1]), ".tif"))
plot(map)

library(ggplot2)
library(ggspatial)
crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# pni <- terra::vect("../input/shp/PARNA_Itatiaia_SIRGAS_IMG.shp", crs= crs_ref)
pni_alta <- terra::vect("../input/shp/parte alta.shp", crs= crs_ref)
# pnipts <- terra::vect("../input/shp/amostragem_cLHS_100.shp", crs= crs_ref)

### Convert train and val to points
train_kfold_poins <- terra::vect(df.mod[, 192:193-10], geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)
ext_val_poins <- terra::vect(df.ex.val[, 41:42], geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)

### Plot
ggplot() +
  layer_spatial(map, aes(alpha = stat(band1)), fill = "#8c510a") +
  scale_alpha_continuous(na.value = 0, name = bquote(paste(.(names(sp_spc[1])), " ", .(unit)))) +
  layer_spatial(pni_alta, aes(), color= "black", fill = "transparent", size=0.5) +
  # layer_spatial(pnipts, aes(), color= "black", size=0.5) +
  layer_spatial(train_kfold_poins, aes(), color= "blue", size=0.5) +
  layer_spatial(ext_val_poins, aes(), color= "green", size=0.5) +
  ggtitle("Itatiaia National Park", subtitle = paste0("Predicted ", names(sp_spc[1]), "\n", "DSM model: Sentinel + Terrain + Geographic")) + xlab("") + ylab("") + 
  annotation_scale(location = "tl") + annotation_north_arrow(location = "br", which_north = "true") + 
  annotate("text", x = 525900, y = 7537500, label = c(paste("R²", round((data.frame(coef_best_model_kf_cv))[,1], 2), "\n",
                                                            "MSE", round((data.frame(coef_best_model_kf_cv))[,3], 2), "\n",
                                                            "RMSE", round((data.frame(coef_best_model_kf_cv))[,4], 2), "\n",
                                                            "bias", round((data.frame(coef_best_model_kf_cv))[,5], 2))), parse = F) + 
  theme_light() + coord_sf(datum= 32723)
ggsave(paste0("../output/graphs_plots/", ds, "/", ds, "_", "Map_of_predicted_", names(sp_spc[1]), ".png"), 
       plot = last_plot(), width = 15, height = 15, units = "cm", dpi = 300)
# ggsave(paste0("../output/graphs_plots/", ds, "/", ds, "_", "Map_of_predicted_", names(sp_spc[1]), ".pdf"), 
# plot = last_plot(), width = 15, height = 15, units = "cm", dpi = 300)

### coeficients
# MSE <- round(rf_model$mse[rf_model$ntree],1)  
# R2 <- round(rf_model$rsq[rf_model$ntree],2)  
# RMSE <- round(sqrt(rf_model$mse[rf_model$ntree]),1)  
# result <- data.frame(MSE, R2, RMSE, ds)
# result
# write.csv(result, file = paste0("../data/", tipo, "_result_",length(raster_stack_chris_raw_84@layers),"_bands.csv"))
## ----cumulative_plots, message=F-----------------------------------------
# plot(rf_model$mse, xlab = "trees", ylab = "MSE")
# plot(rf_model$rsq, xlab = "trees", ylab = "Rsquared")

# save.image("../data/results_RF_rfe.RData")
# load("../data/data.RData")



################################################################################
################################################################################

# Graficos mostrando a importancia das covariaveis #
# round(randomForest::importance(rf_models[[i]]),2) # mostra a importancia de todas as covariaveis #
# plot(sort(rf_models[[i]]$importance),axes=T, ylab="Importancia", xlab="covariaveis") # plotar importancia - verificar a figura sem nome das covariaveis #

# for (h in 1:13){ # plotar as 5 primeiras classes
#    plot(sort(rf_models[[i]]$importance[,h], dec=T),
#         type="h", main=paste("Medida", h), ylab= "importancia")}

# varImpPlot(rf_models[[i]], pch = 19)
# varImpPlot(rf_models[[i+1]], pch = 19)
# varImpPlot(rf_models[[i+2]], pch = 19)

# Analise para determinar o n de arvores #
# par(mfrow=c(1,1))
# plot(rf_models[[i]]) 

# rf_models[[i+1]] <- randomForest(form, data= trainData, ntree= 250, mtry= 5, nodesize= 5) 
# rf_models[[i+2]] <- randomForest(form, data= trainData) 

# rf_models[[i]]$mse


# inspect the output
# str(rf_model, max.level=2)

# inspect the top of the first tree in the forest.
# randomForest::getTree(rf_models[[i]])
# head(getTree(rf_models[[i]], k=1, labelVar = TRUE))


###########################################################################################################################
###########################################################################################################################
###################################################################################
#      _/_/_/  _/_/_/    _/_/_/
#   _/        _/    _/  _/    _/
#    _/_/    _/_/_/    _/_/_/
#       _/  _/    _/  _/
#  _/_/_/    _/    _/  _/
# figlet -f lean SRP  ### run on terminal ctrl+alt+enter
#   ____             __           __  ___           __
#  / __/__  ___ ____/ /________ _/ / / _ \___ ____ / /____ ____
# _\ \/ _ \/ -_) __/ __/ __/ _ `/ / / , _/ _ `(_-</ __/ -_) __/
#/___/ .__/\__/\__/\__/_/  \_,_/_/ /_/|_|\_,_/___/\__/\__/_/
#   /_/
#   ___             ___     __  _
#  / _ \_______ ___/ (_)___/ /_(_)__  ___
# / ___/ __/ -_) _  / / __/ __/ / _ \/ _ \
#/_/  /_/  \__/\_,_/_/\__/\__/_/\___/_//_/
# 
# figlet -f smslant Spectral  Raster  Prediction ### run on terminal ctrl+alt+enter
# showfigfonts
###################################################################################

### "Hand made" K-fold cross validation
#' ### Auto setwd
library(rstudioapi); current_path <- getActiveDocumentContext()$path; setwd(dirname(current_path)) ### https://eranraviv.com/r-tips-and-tricks-working-directory/
getwd()
# gc(); rm(list=ls())

{
  library(terra)
  library(fst)
  library(randomForest)
  library(readr)
  library(ithir)
  library(hexbin)
  library(grid)
  
  ### Load data
  stk_subsuf_img <- terra::rast("../output/raster/stk_subsuf_img_10m.tif")
  df.mod_subsuf_img <- read_fst("../output/dataset/dataset_spectral_image_covariates_df_mod_point_pixels.fst")
  ext_valData_spectral <- read_fst("../output/dataset/dataset_Ext_val_spectral_covariates.fst")
  df.pred_spectr_subsuf <- read_fst("../output/dataset/dataset_df.pred_spectr_subsuf.fst")
  
  ### Split data into 10 folds
  folds <- cut(seq(1,nrow(df.mod_subsuf_img)), breaks= 8, labels=FALSE)
  nflds <- length(unique(folds))
  # table(folds)
}

# names(df.mod_subsuf_img)
df.mod_subsuf_img_sbs <- df.mod_subsuf_img[, c(40:179)] ### 40:169 the spectral data
# names(df.mod_subsuf_img_sbs)
ext_valData_spectral_sbs <- ext_valData_spectral[, c(40:179)] ### 40:169 the spectral data
# names(ext_valData_spectral_sbs)

# rm(var_x, var_y, predictorVars, form)

system.time(
  for(j in c(1:130)){ # Across entire Spectral data
  # for(j in c(1:3)){ # Test
    ### building the formula to use in the algorithm
    var_x <- df.mod_subsuf_img_sbs[, 131:140]
    var_y <- df.mod_subsuf_img_sbs[j]
    predictorVars <- paste(names(var_x), collapse = "+")
    form = as.formula(paste0(names(var_y), "~", predictorVars))
    
    ### create null objects (lists) to fill with the results: rf
    rf_models <- list() ### results from the models
    
    fitted_value_rf_cv <- list() ### results from the models (for cv)
    goof_rf_cv <- list() ### metrics
    goof_rf_cv_df <- list() ### metrics
    
    fitted_value_rf_ext_val <- list() ### results from the models (for external validation)
    goof_rf_ext_val <- list() ### metrics
    goof_rf_ext_val_df <- list() ### metrics
    
    fitted_value_rf_calibr <- list()  ### results from the models (for calibration)
    goof_rf_calibr <- list()
    goof_rf_calibr_df <- list()
    
    for(i in 1:nflds){ ### This loop is to run and validate with the K-folds
      ### Segment your data by fold using the which() function 
      testIndexes <- which(folds== i, arr.ind=TRUE)
      testData <- df.mod_subsuf_img_sbs[testIndexes, ]
      trainData <- df.mod_subsuf_img_sbs[-testIndexes, ]
      
      ######################################
      #### rf_model (trainig data)
      ## set.seed(1234) # to fix the model randomness (do not need in these cases)
      rf_models[[i]] <- randomForest(form, data= trainData, ntree= 150, mtry= 3, nodesize= 3, importance= T) 
      
      ################ Calibration
      ### predict with the model (test data), and store in "i", with means in the loop of the folders 
      fitted_value_rf_calibr[[i]] <- predict(rf_models[[i]], newdata= trainData) 
      
      ### perform goof for calibr and store as list
      goof_rf_calibr[[i]] <- goof(trainData[,j], fitted_value_rf_calibr[[i]], type='DSM') 
      goof_rf_calibr[[i]]$model_number <- i # paste0("Model ", i)
      
      ### extrated the values from the list inside each property group and store as list object
      goof_rf_calibr_df[[j]] <- (as.data.frame(do.call(rbind, goof_rf_calibr)))
      
      ### compute mean for each metric of all properties across all folders store as list object
      ### and bind with folders metrics
      goof_rf_calibr_metrics <- as.data.frame(rbind(goof_rf_calibr_df[[j]], colMeans(x= goof_rf_calibr_df[[j]])))
      
      ### fill "mean" model with zero
      # if (!is.na(goof_rf_calibr_metrics$model_number[9])) {
      #   goof_rf_calibr_metrics$model_number[9] <- 0
      # }
      
      ### Write results on disk
      ds <- "Spect_subsuf"
      write_csv(round(goof_rf_calibr_metrics, 2),
                file = paste0("../output/dataset/", ds, "/", ds, "_",
                              names(var_y), "_goof_rf_calibr_metrics.csv"))
      
      ################ k-fold
      ### predict with the model (test data), and store in "i", with means in the loop of the folders 
      fitted_value_rf_cv[[i]] <- predict(rf_models[[i]], newdata= testData) 
      
      ### perform goof for cv and store as list
      goof_rf_cv[[i]] <- goof(testData[,j], fitted_value_rf_cv[[i]], type='DSM') 
      goof_rf_cv[[i]]$model_number <- i # paste0("Model ", i)
      
      ### extrated the values from the list inside each property group and store as list object
      goof_rf_cv_df[[j]] <- (as.data.frame(do.call(rbind, goof_rf_cv)))
      
      ### compute mean for each metric of all properties across all folders store as list object
      ### and bind with folders metrics
      goof_rf_k_cv_metrics <- as.data.frame(rbind(goof_rf_cv_df[[j]], colMeans(x= goof_rf_cv_df[[j]])))
      
      ### fill "mean" model with zero
      # if (!is.na(goof_rf_k_cv_metrics$model_number[9])) {
      #   goof_rf_k_cv_metrics$model_number[9] <- 0
      # }
      
      ### Write results on disk
      write_csv(round(goof_rf_k_cv_metrics, 2),
                file = paste0("../output/dataset/", ds, "/", ds, "_",
                              names(var_y), "_goof_rf_k_fold_cv_metrics.csv"))
      
      # names(ext_valData_spectral_sbs)
      ################ External validation
      ### predict with the model (test data), and store in "i", with means in the loop of the folders 
      fitted_value_rf_ext_val[[i]] <- predict(rf_models[[i]], newdata= ext_valData_spectral_sbs) 
      
      ### perform goof for external validation and store as list
      goof_rf_ext_val[[i]] <- goof(ext_valData_spectral_sbs[,j], fitted_value_rf_ext_val[[i]], type='DSM') 
      goof_rf_ext_val[[i]]$model_number <- i # paste0("Model ", i)
      
      ### extrated the values from the list inside each property group and store as list object
      goof_rf_ext_val_df[[j]] <- (as.data.frame(do.call(rbind, goof_rf_ext_val)))
      
      ### compute mean for each metric of all properties across all folders store as list object
      ### and bind with folders metrics
      goof_rf_ext_val_metrics <- as.data.frame(rbind(goof_rf_ext_val_df[[j]], colMeans(x= goof_rf_ext_val_df[[j]])))
      
      ### fill "mean" model with zero
      # if (!is.na(goof_rf_ext_val_metrics$model_number[9])) {
      #   goof_rf_ext_val_metrics$model_number[9] <- 0
      # }
      
      ### Write results on disk
      write_csv(round(goof_rf_ext_val_metrics, 2),
                file = paste0("../output/dataset/", ds, "/", ds, "_",
                              names(var_y), "_goof_rf_mean_ext_val.csv"))
      ### Save outputs
      w_list <- c("rf_models", 
                  ### Calibration
                  "fitted_value_rf_calibr",
                  "goof_rf_calibr",
                  "goof_rf_calibr_df",
                  "goof_rf_calibr_metrics",
                  ### k-fold cv
                  "fitted_value_rf_cv",
                  "goof_rf_cv",
                  "goof_rf_cv_df",
                  "goof_rf_k_cv_metrics",
                  ### ext val
                  "fitted_value_rf_ext_val",
                  "goof_rf_ext_val",
                  "goof_rf_ext_val_df",
                  "goof_rf_ext_val_metrics")
      
      ### Save in a loop
      for (sv in seq_along(w_list)) {
        saveRDS(eval(parse(text= w_list[[sv]])),
                paste0("../output/dataset/", ds, "/", ds, "_",
                       names(var_y), "_", w_list[[sv]], ".rds"), compress= TRUE)
      }
      
    }

################ Select best model to predict the map
# best_model_calibr <- goof_rf_calibr_df[[j]][goof_rf_calibr_df[[j]]$MSE == min(goof_rf_calibr_df[[j]]$MSE), "model_number"]
# coef_best_model_calibr <- goof_rf_calibr_df[[j]][goof_rf_calibr_df[[j]]$MSE == min(goof_rf_calibr_df[[j]]$MSE), ]

best_model_kf_cv <- goof_rf_cv_df[[j]][goof_rf_cv_df[[j]]$MSE == min(goof_rf_cv_df[[j]]$MSE), "model_number"]
# coef_best_model_kf_cv <- goof_rf_cv_df[[j]][goof_rf_cv_df[[j]]$MSE == min(goof_rf_cv_df[[j]]$MSE), ]

# best_model_val_ext <- goof_rf_ext_val_df[[j]][goof_rf_ext_val_df[[j]]$MSE == min(goof_rf_ext_val_df[[j]]$MSE), "model_number"]
# coef_best_model_val_ext <- goof_rf_ext_val_df[[j]][goof_rf_ext_val_df[[j]]$MSE == min(goof_rf_ext_val_df[[j]]$MSE), ]

### Predict the map pixels with best model
pred_map <- stats::predict(rf_models[[best_model_kf_cv]], df.pred_spectr_subsuf[, 1:10]) # to map

### Join predicted with X, Y
spectr_pred = cbind(pred_map, df.pred_spectr_subsuf[, 11:12])
names(spectr_pred)[1] <- names(var_y)
head(spectr_pred) 

### Export
### As fst
# system.time({ #
#   write_fst(spectr_pred, paste0("../output/raster/predict_rasters_results/", ds, "/", ds, "_predict_", names(var_y), ".fst"), compress = 100)
# })
### As RDS
# system.time({ 
#   saveRDS(spectr_pred, paste0("../output/raster/predict_rasters_results/", ds, "/", ds, "_predict_", names(var_y), ".rds"), compress= TRUE) # lighter than zip
# })
## As tif
system.time({ #
### Convert to points
crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
points <- terra::vect(spectr_pred, geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)

### Convert points to raster
map <- terra::rasterize(points, stk_subsuf_img[[1]] , names(spectr_pred)[1])
# plot(map)

# Export each raster
terra::writeRaster(map, filename= paste0("../output/raster/predict_rasters_results/", ds, "/", ds, "_predict_", names(var_y), ".tif"), overwrite=TRUE)
})

### Liberate memory (except)
rm(list=setdiff(ls(), c("stk_subsuf_img", 
                        "df.mod_subsuf_img", 
                        "ext_valData_spectral", 
                        "df.pred_spectr_subsuf", 
                        "folds", 
                        "nflds", 
                        "df.mod_subsuf_img_sbs", 
                        "ext_valData_spectral_sbs")))

  }
)

### Elapsed time
# 1134/60 # with fst 19 min
# 4165/60 # with writeRaster 70 min

################################################################################
### Load every predicted spectral image (as dataframe in fst format)
################################################################################
# stk_dsm <- terra::rast("../output/raster/stk_dsm_10m.tif")
# df.mod <- read_fst("../output/dataset/dataset_HDSM_covariates_df_mod_point_pixels.fst")
# ### Names of objects (spectral images) to read
# fst_list <- names(df.mod[, 40:169])
# 
# ### Read predicted spectral image
# ds <- "Spect_subsuf"
# system.time(
#   for (fsl in seq_along(fst_list)) {
#     assign(fst_list[[fsl]], read_fst(
#       paste0("../output/raster/predict_rasters_results/", ds, "/", 
#              ds, "_predict_", fst_list[[fsl]], ".fst")))
#   }) # 20 seconds
# 
# ### one predicted image used to get the dimensions of df
# X430 
# # df.pred.spectral <- data.frame(X430)
# df.pred.spectral <- data.frame(a= matrix(NA, ncol=1, nrow= length(X430[[1]])))[-1]
# for (fsl in seq_along(fst_list)) {
#   df.pred.spectral <- cbind(df.pred.spectral, eval(parse(text= fst_list[[fsl]]))[1]) ### [1] selects only the predict, and leave all those x and y
# }
# ### Write on disk
# system.time({ #
#   write_fst(df.pred.spectral, paste0("../output/raster/predict_rasters_results/", ds, "/", ds, "_df.pred.spectral.fst"), compress = 100)
# })
# ### Liberate memory (except)
# rm(list=setdiff(ls(), c("df.pred.spectral")))
# ### Read df.pred.spectral
# ds <- "Spect_subsuf"
# system.time({
#   df.pred.spectral <- read_fst(paste0("../output/raster/predict_rasters_results/", ds, "/", ds, "_df.pred.spectral.fst"))
# }) # 12 seconds 

#################################################################################################################
###  _   _ ____  __  __        _       _                                                  _   _               ### 
### | | | / ___||  \/  |    __| | __ _| |_ __ _      _ __  _ __ ___ _ __   __ _ _ __ __ _| |_(_) ___  _ __    ### 
### | |_| \___ \| |\/| |   / _` |/ _` | __/ _` |    | '_ \| '__/ _ \ '_ \ / _` | '__/ _` | __| |/ _ \| '_ \   ### 
### |  _  |___) | |  | |  | (_| | (_| | || (_| |    | |_) | | |  __/ |_) | (_| | | | (_| | |_| | (_) | | | |  ### 
### |_| |_|____/|_|  |_|   \__,_|\__,_|\__\__,_|    | .__/|_|  \___| .__/ \__,_|_|  \__,_|\__|_|\___/|_| |_|  ### 
###                                                 |_|            |_|                                        ### 
##### HSM data preparation                                                                                    ### 
#################################################################################################################
#' ### Auto setwd
library(rstudioapi); current_path <- getActiveDocumentContext()$path; setwd(dirname(current_path)) ### https://eranraviv.com/r-tips-and-tricks-working-directory/
getwd()
# gc(); rm(list=ls())
{
  library(terra)
  library(fst)
  library(randomForest)
  library(readr)
  library(ithir)
  library(hexbin)
  library(grid)
}

### Load every predicted spectral image (as .tif)
library(terra)
ds <- "Spect_subsuf"
pf <- paste0("../output/raster/predict_rasters_results/", ds, "/")
rastlist <- list.files(path = pf, pattern='.tif$', all.files= T, full.names= T)
spec_stk <- terra::rast(rastlist)
# plot(spec_stk[[120]])

{
  stk_dsm <- terra::rast("../output/raster/stk_dsm_10m.tif")
  df.mod <- read_fst("../output/dataset/dataset_HDSM_covariates_df_mod_point_pixels.fst")
  df.ex.val <- read_fst("../output/dataset/dataset_Ext_val_covariates_df_external_validation.fst")
  df.pred <- read_fst("../output/dataset/dataset_22_covariates_df_pred_all_pixels.fst")
}

# names(df.mod)
# names(df.pred)
# names(stk_dsm)

stk_hsm <- c(spec_stk, stk_dsm)
names(stk_hsm)[1:130] <- names(df.mod[40:(130+39)])

# names(stk_hsm)
# names(stk_hsm)[131:142] == names(df.mod)[131:142+39]
# names(stk_hsm)[1:142] == names(df.mod)[1:142+39]
# table(df.pred["x"] == df.pred.spectral["x"]) ### It was tested and they are all equal
# table(df.pred["y"] == df.pred.spectral["y"]) ### It was tested and they are all equal

################################################################################
### Extract data for spectral external validation
################################################################################
### Get the validation data X and Y from
df.ex.val

### Get points X and Y
crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
points <- vect(df.ex.val, geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)
plot(points)

### Extract all info from stack
system.time({df.ex.val_subsuf_img = terra::extract(stk_hsm, 
                                                   points, 
                                                   method= "simple", 
                                                   xy= T, 
                                                   ID= F,
                                                   bind= T # return a SpatVector object, with all variables
)})

df.ex.val_subsuf_img
dim(df.ex.val_subsuf_img)
df.ex.val_subsuf_img <- as.data.frame(df.ex.val_subsuf_img) # from SpatVector to df
# str(df.ex.val_subsuf_img)

system.time({ #
  write_fst(df.ex.val_subsuf_img, "../output/dataset/df.ex.val_subsuf_img.fst", compress = 100)
})

df.ex.val_subsuf_img <- read_fst("../output/dataset/df.ex.val_subsuf_img.fst")

################################################################################
### Extract data for spectral prediction
################################################################################
### Get points X and Y
pni_alta <- terra::vect("../input/shp/parte alta.shp", crs= crs_ref)
# plot(pni_alta)

### Extract all info from stack
system.time({df.subsuf.pred = terra::extract(stk_hsm, 
                                             pni_alta, 
                                             method= "simple", 
                                             xy= T, 
                                             ID= F,
                                             # bind= T # return a SpatVector object, with all variables
)}) # 132 seconds

df.subsuf.pred
dim(df.subsuf.pred)
# df.subsuf.pred <- as.data.frame(df.subsuf.pred) # from SpatVector to df
# str(df.subsuf.pred)

system.time({ #
  write_fst(df.subsuf.pred, "../output/dataset/df.subsuf.pred.fst", compress = 100)
}) # 31 seconds

system.time({ #
  df.subsuf.pred <- read_fst("../output/dataset/df.subsuf.pred.fst")
}) # 2 seconds


################################################################################
"
HHHHHHHHH     HHHHHHHHH   SSSSSSSSSSSSSSS MMMMMMMM               MMMMMMMM
H:::::::H     H:::::::H SS:::::::::::::::SM:::::::M             M:::::::M
H:::::::H     H:::::::HS:::::SSSSSS::::::SM::::::::M           M::::::::M
HH::::::H     H::::::HHS:::::S     SSSSSSSM:::::::::M         M:::::::::M
  H:::::H     H:::::H  S:::::S            M::::::::::M       M::::::::::M
  H:::::H     H:::::H  S:::::S            M:::::::::::M     M:::::::::::M
  H::::::HHHHH::::::H   S::::SSSS         M:::::::M::::M   M::::M:::::::M
  H:::::::::::::::::H    SS::::::SSSSS    M::::::M M::::M M::::M M::::::M
  H:::::::::::::::::H      SSS::::::::SS  M::::::M  M::::M::::M  M::::::M
  H::::::HHHHH::::::H         SSSSSS::::S M::::::M   M:::::::M   M::::::M
  H:::::H     H:::::H              S:::::SM::::::M    M:::::M    M::::::M
  H:::::H     H:::::H              S:::::SM::::::M     MMMMM     M::::::M
HH::::::H     H::::::HHSSSSSSS     S:::::SM::::::M               M::::::M
H:::::::H     H:::::::HS::::::SSSSSS:::::SM::::::M               M::::::M
H:::::::H     H:::::::HS:::::::::::::::SS M::::::M               M::::::M
HHHHHHHHH     HHHHHHHHH SSSSSSSSSSSSSSS   MMMMMMMM               MMMMMMMM
"
### HSM data loading and run
################################################################################
### "Hand made" K-fold cross validation
#' ### Auto setwd
library(rstudioapi); current_path <- getActiveDocumentContext()$path; setwd(dirname(current_path)) ### https://eranraviv.com/r-tips-and-tricks-working-directory/
getwd()
# gc(); rm(list=ls())
{
  library(terra)
  library(fst)
  library(randomForest)
  library(readr)
  library(ithir)
  library(hexbin)
  library(grid)
  
  stk_dsm <- terra::rast("../output/raster/stk_dsm_10m.tif")
  df.mod <- read_fst("../output/dataset/dataset_HDSM_covariates_df_mod_point_pixels.fst")
  df.ex.val <- read_fst("../output/dataset/dataset_Ext_val_covariates_df_external_validation.fst")
  df.pred <- read_fst("../output/dataset/dataset_22_covariates_df_pred_all_pixels.fst")
  
  ### Load every predicted spectral image (as .tif)
  ds <- "Spect_subsuf"
  pf <- paste0("../output/raster/predict_rasters_results/", ds, "/")
  rastlist <- list.files(path = pf, pattern='.tif$', all.files= T, full.names= T)
  spec_stk <- terra::rast(rastlist)
  # plot(spec_stk[[120]])
  
  stk_hsm <- c(spec_stk, stk_dsm)
  names(stk_hsm)[1:130] <- names(df.mod[40:(130+39)])
  
  df.ex.val_subsuf_img <- read_fst("../output/dataset/df.ex.val_subsuf_img.fst")
  df.subsuf.pred <- read_fst("../output/dataset/df.subsuf.pred.fst")
}

################################################################################
### Start the Random forest model for HSM
################################################################################

### Split data into 10 folds
folds <- cut(seq(1,nrow(df.mod)), breaks= 8, labels=FALSE)
nflds <- length(unique(folds))
# table(folds)

{
  ### create null objects (lists) to fill with the results: rf
  rf_models <- list() ### results from the models
  
  fitted_value_rf_cv <- list() ### results from the models (for cv)
  goof_rf_cv <- list() ### metrics
  goof_rf_cv_df <- list() ### metrics
  # goof_rf_mean_cv_df <- list() ### metrics
  
  fitted_value_rf_ext_val <- list() ### results from the models (for external validation)
  goof_rf_ext_val <- list() ### metrics
  goof_rf_ext_val_df <- list() ### metrics
  # goof_rf_mean_ext_val_df <- list() ### metrics
  
  fitted_value_rf_calibr <- list()  ### results from the models (for calibration)
  goof_rf_calibr <- list()
  goof_rf_calibr_df <- list()
  # goof_rf_mean_calibr_df <- list()
}

system.time(
  for(j in c(21)){ ### C=21
    ### building the formula to use in the algorithm
    allVars <- colnames(df.mod[, c(j, 40:(191-10))]) ## DSM 170:181, HSM 40:181
    sp_spc <- df.mod[, c(j, 40:(191-10))]
    predictorVars <- allVars[!allVars%in%names(sp_spc[1])]
    predictorVars <- paste(predictorVars, collapse = "+")
    form = as.formula(paste(names(sp_spc[1]), "~", predictorVars))
    
    for(i in 1:nflds){ ### This loop is to run and validate with the K-folds
      ### Segment your data by fold using the which() function 
      testIndexes <- which(folds== i, arr.ind=TRUE)
      testData <- df.mod[testIndexes, ]
      trainData <- df.mod[-testIndexes, ]
      
      ######################################
      #### rf_model (trainig data)
      set.seed(1234) # to fix the model randomness
      rf_models[[i]] <- randomForest(form, data= trainData, ntree= 150, mtry= 3, nodesize= 3, importance= T) 
      
      ################ Calibration
      ### predict with the model (test data), and store in "i", with means in the loop of the folders 
      fitted_value_rf_calibr[[i]] <- predict(rf_models[[i]], newdata= trainData[ , 40:(191-10)]) 
      
      ### perform goof for calibr and store as list
      goof_rf_calibr[[i]] <- goof(trainData[,j], fitted_value_rf_calibr[[i]], type='DSM') 
      goof_rf_calibr[[i]]$model_number <- i # paste0("Model ", i)
      
      ### extrated the values from the list inside each property group and store as list object
      goof_rf_calibr_df[[j]] <- (as.data.frame(do.call(rbind, goof_rf_calibr)))
      
      ### compute mean for each metric of all properties across all folders store as list object
      # goof_rf_mean_calibr_df[[j]] <- colMeans(x= goof_rf_calibr_df[[j]])
      
      ### pass to a dataframe
      # df_goof_rf_mean_calibr <- as.data.frame(t(goof_rf_mean_calibr_df[[j]]))
      # df_goof_rf_mean_calibr$model_number <- NULL
      
      ### extrated the values from the list and store as data.frame object
      # goof_rf_mean_calibr_df <- (as.data.frame(do.call(rbind, goof_rf_mean_calibr_df)))
      # goof_rf_mean_calibr_df$propertie <- names(sp_spc[1])
      # goof_rf_mean_calibr_df$model_number <- NULL
      
      ### compute mean for each metric of all properties across all folders store as list object
      ### and bind with folders metrics
      goof_rf_calibr_metrics <- as.data.frame(rbind(goof_rf_calibr_df[[j]], colMeans(x= goof_rf_calibr_df[[j]])))
      
      ### fill "mean" model with zero
      if (!is.na(goof_rf_calibr_metrics$model_number[9])) {
        goof_rf_calibr_metrics$model_number[9] <- 0
      }
      
      ### Write results on disk
      ds <- "HSM"
      # write_csv(round(goof_rf_calibr_df[[j]], 2),
      #           file = paste0("../output/dataset/", ds, "/", ds, "_",
      #                         names(sp_spc[1]), "_goof_rf_folds_calibr.csv"))
      # write_csv(round(df_goof_rf_mean_calibr, 2),
      #           file = paste0("../output/dataset/", ds, "/", ds, "_",
      #                         names(sp_spc[1]), "_goof_rf_mean_calibra.csv"))
      write_csv(round(goof_rf_calibr_metrics, 2),
                file = paste0("../output/dataset/", ds, "/", ds, "_",
                              names(sp_spc[1]), "_goof_rf_calibr_metrics.csv"))
      
      ################ k-fold
      ### predict with the model (test data), and store in "i", with means in the loop of the folders 
      fitted_value_rf_cv[[i]] <- predict(rf_models[[i]], newdata= testData[ , 40:(191-10)]) 
      
      ### perform goof for cv and store as list
      goof_rf_cv[[i]] <- goof(testData[,j], fitted_value_rf_cv[[i]], type='DSM') 
      goof_rf_cv[[i]]$model_number <- i # paste0("Model ", i)
      
      ### extrated the values from the list inside each property group and store as list object
      goof_rf_cv_df[[j]] <- (as.data.frame(do.call(rbind, goof_rf_cv)))
      
      ### compute mean for each metric of all properties across all folders store as list object
      # goof_rf_mean_cv_df[[j]] <- colMeans(x= goof_rf_cv_df[[j]])
      
      ### pass to a dataframe
      # df_goof_rf_mean_cv <- as.data.frame(t(goof_rf_mean_cv_df[[j]]))
      # df_goof_rf_mean_cv$model_number <- NULL
      
      ### extrated the values from the list and store as data.frame object
      # goof_rf_mean_cv_df <- (as.data.frame(do.call(rbind, goof_rf_mean_cv_df)))
      # goof_rf_mean_cv_df$propertie <- names(sp_spc[1])
      # goof_rf_mean_cv_df$model_number <- NULL
      
      ### compute mean for each metric of all properties across all folders store as list object
      ### and bind with folders metrics
      goof_rf_k_cv_metrics <- as.data.frame(rbind(goof_rf_cv_df[[j]], colMeans(x= goof_rf_cv_df[[j]])))
      
      ### fill "mean" model with zero
      if (!is.na(goof_rf_k_cv_metrics$model_number[9])) {
        goof_rf_k_cv_metrics$model_number[9] <- 0
      }
      
      ### Write results on disk
      # write_csv(round(goof_rf_cv_df[[j]], 2),
      #           file = paste0("../output/dataset/", ds, "/", ds, "_",
      #                         names(sp_spc[1]), "_goof_rf_folds_k_fold_cv.csv"))
      # write_csv(round(df_goof_rf_mean_cv, 2),
      #           file = paste0("../output/dataset/", ds, "/", ds, "_",
      #                         names(sp_spc[1]), "_goof_rf_mean_k_fold_cv.csv"))
      write_csv(round(goof_rf_k_cv_metrics, 2),
                file = paste0("../output/dataset/", ds, "/", ds, "_",
                              names(sp_spc[1]), "_goof_rf_k_fold_cv_metrics.csv"))
      
      names(df.ex.val_subsuf_img)
      ################ External validation
      ### predict with the model (test data), and store in "i", with means in the loop of the folders 
      fitted_value_rf_ext_val[[i]] <- predict(rf_models[[i]], newdata= df.ex.val_subsuf_img[ , 41:182])
      
      ### perform goof for external validation and store as list
      goof_rf_ext_val[[i]] <- goof(df.ex.val_subsuf_img[,j], fitted_value_rf_ext_val[[i]], type='DSM') 
      goof_rf_ext_val[[i]]$model_number <- i # paste0("Model ", i)
      
      ### extrated the values from the list inside each property group and store as list object
      goof_rf_ext_val_df[[j]] <- (as.data.frame(do.call(rbind, goof_rf_ext_val)))
      
      ### compute mean for each metric of all properties across all folders store as list object
      # goof_rf_mean_ext_val_df[[j]] <- colMeans(x= goof_rf_ext_val_df[[j]])
      
      ### pass to a dataframe
      # df_goof_rf_mean_ext_val <- as.data.frame(t(goof_rf_mean_ext_val_df[[j]]))
      # df_goof_rf_mean_ext_val$model_number <- NULL
      
      ### extrated the values from the list and store as data.frame object
      # goof_rf_mean_ext_val_df <- (as.data.frame(do.call(rbind, goof_rf_mean_ext_val_df)))
      # goof_rf_mean_ext_val_df$propertie <- names(sp_spc[1])
      # goof_rf_mean_ext_val_df$model_number <- NULL
      
      ### compute mean for each metric of all properties across all folders store as list object
      ### and bind with folders metrics
      goof_rf_ext_val_metrics <- as.data.frame(rbind(goof_rf_ext_val_df[[j]], colMeans(x= goof_rf_ext_val_df[[j]])))
      
      ### fill "mean" model with zero
      if (!is.na(goof_rf_ext_val_metrics$model_number[9])) {
        goof_rf_ext_val_metrics$model_number[9] <- 0
      }
      
      ### Write results on disk
      # write_csv(round(goof_rf_ext_val_df[[j]], 2),
      #           file = paste0("../output/dataset/", ds, "/", ds, "_",
      #                         names(sp_spc[1]), "_goof_rf_folds_ext_val.csv"))
      # write_csv(round(df_goof_rf_mean_ext_val, 2),
      #           file = paste0("../output/dataset/", ds, "/", ds, "_",
      #                         names(sp_spc[1]), "_goof_rf_mean_ext_val.csv"))
      write_csv(round(goof_rf_ext_val_metrics, 2),
                file = paste0("../output/dataset/", ds, "/", ds, "_",
                              names(sp_spc[1]), "_goof_rf_mean_ext_val.csv"))
      ### Save outputs
      w_list <- c("rf_models", 
                  ### Calibration
                  "fitted_value_rf_calibr",
                  "goof_rf_calibr",
                  "goof_rf_calibr_df",
                  # "goof_rf_mean_calibr_df",
                  # "df_goof_rf_mean_calibr",
                  "goof_rf_calibr_metrics",
                  ### k-fold cv
                  "fitted_value_rf_cv",
                  "goof_rf_cv",
                  "goof_rf_cv_df",
                  # "df_goof_rf_mean_cv",
                  # "goof_rf_mean_cv_df",
                  "goof_rf_k_cv_metrics",
                  ### ext val
                  "fitted_value_rf_ext_val",
                  "goof_rf_ext_val",
                  "goof_rf_ext_val_df",
                  # "goof_rf_mean_ext_val_df",
                  "goof_rf_ext_val_metrics")
      
      # w_list[1] # write lists
      # w_list[19] # write lists
      
      ### Save in a loop
      for (sv in seq_along(w_list)) {
        saveRDS(eval(parse(text= w_list[[sv]])),
                paste0("../output/dataset/", ds, "/", ds, "_",
                       names(sp_spc[1]), "_", w_list[[sv]], ".rds"), compress= TRUE)
      }
      
    }
  }
)
### Liberate memory (except)
# rm(list=setdiff(ls(), c("rf_models", "stk_dsm", "df.mod", "df.ex.val", "df.pred", "folds", "nflds")))

### Read outputs
# ds <- "HSM"
# propertie <- "C"
# for (sv in seq_along(w_list)) {
#   assign(w_list[[sv]], readRDS(
#     paste0("../output/dataset/", ds, "/", ds, "_",
#            propertie, "_", w_list[[sv]],".rds")))
# }
# j=21

################ Select best model to predict the map
fst_qt_cal <- goof_rf_calibr_df[[j]][quantile(goof_rf_calibr_df[[j]]$RMSE, 0.25) > goof_rf_calibr_df[[j]]$RMSE, ]
best_model_calibr <- fst_qt_cal[fst_qt_cal$R2 == max(fst_qt_cal$R2), "model_number"]
coef_best_model_calibr <- fst_qt_cal[fst_qt_cal$R2 == max(fst_qt_cal$R2), ]

fst_qt_kcv <- goof_rf_cv_df[[j]][quantile(goof_rf_cv_df[[j]]$RMSE, 0.25) > goof_rf_cv_df[[j]]$RMSE, ]
best_model_kf_cv <- fst_qt_kcv[fst_qt_kcv$R2 == max(fst_qt_kcv$R2), "model_number"]
coef_best_model_kf_cv <- fst_qt_kcv[fst_qt_kcv$R2 == max(fst_qt_kcv$R2), ]

fst_qt_ex_val <- goof_rf_ext_val_df[[j]][quantile(goof_rf_ext_val_df[[j]]$RMSE, 0.25) > goof_rf_ext_val_df[[j]]$RMSE, ]
best_model_val_ext <- fst_qt_ex_val[fst_qt_ex_val$R2 == max(fst_qt_ex_val$R2), "model_number"]
coef_best_model_val_ext <- fst_qt_ex_val[fst_qt_ex_val$R2 == max(fst_qt_ex_val$R2), ]

### Predict the map pixels with best model
# df.subsuf.pred <- read_fst("../output/dataset/df.subsuf.pred.fst")
### Liberate memory (except)
rm(list=setdiff(ls(), c("rf_models", 
                        "best_model_kf_cv",
                        "df.subsuf.pred",
                        "coef_best_model_kf_cv")))
# names(df.subsuf.pred)
system.time(
  pred_map <- stats::predict(rf_models[[best_model_kf_cv]], df.subsuf.pred[, 1:142]) # to map
) # 75 seconds

### Join predicted with X, Y
Soil_prop_to_map = cbind(pred_map, df.subsuf.pred[, 143:144]) # df.pred[, 143:144] = X, Y

### Convert to points
crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
points <- terra::vect(Soil_prop_to_map, geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)

### Convert points to raster
rm(list=setdiff(ls(), c("rf_models", "best_model_kf_cv", "coef_best_model_kf_cv", "points"))) ### Liberate memory (except)
stk_dsm <- terra::rast("../output/raster/stk_dsm_10m.tif")
system.time(
  map <- terra::rasterize(points, stk_dsm[[1]] , "pred_map")
) # 22 seconds
# plot(map)

str(rf_models[[best_model_kf_cv]])
pred.property <- rf_models[[best_model_kf_cv]]$terms[[2]]
sp_spc <- as.character(pred.property)

# Export each raster
ds <- "HSM"
terra::writeRaster(map, filename= paste0("../output/raster/predict_rasters_results/", ds, "/", ds, "_predict_", pred.property, ".tif"), overwrite=TRUE)


##############################################################################################
########..##........#######..########..######..........##.....##....###....########...######  
##.....##.##.......##.....##....##....##....##.........###...###...##.##...##.....##.##....##  
##.....##.##.......##.....##....##....##...............####.####..##...##..##.....##.##.....  
########..##.......##.....##....##.....######..........##.###.##.##.....##.########...######  
##........##.......##.....##....##..........##.........##.....##.#########.##..............##  
##........##.......##.....##....##....##....##.........##.....##.##.....##.##........##....##  
##........########..#######.....##.....######..........##.....##.##.....##.##.........######  
##################### Generate figure (maps)  
##############################################################################################
#' ### Auto setwd
library(rstudioapi); current_path <- getActiveDocumentContext()$path; setwd(dirname(current_path)) ### https://eranraviv.com/r-tips-and-tricks-working-directory/
getwd()
gc(); rm(list=ls())

{
library(terra)
library(fst)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(tidyverse)
library(hrbrthemes) # theme ipsum
library(viridis)
library(reshape) # melt data
library(dplyr)
}

################################################################################
### Reload data
################################################################################
### Read outputs
w_list <- c("rf_models", 
            ### Calibration
            "fitted_value_rf_calibr",
            "goof_rf_calibr",
            "goof_rf_calibr_df",
            # "goof_rf_mean_calibr_df",
            # "df_goof_rf_mean_calibr",
            "goof_rf_calibr_metrics",
            ### k-fold cv
            "fitted_value_rf_cv",
            "goof_rf_cv",
            "goof_rf_cv_df",
            # "df_goof_rf_mean_cv",
            # "goof_rf_mean_cv_df",
            "goof_rf_k_cv_metrics",
            ### ext val
            "fitted_value_rf_ext_val",
            "goof_rf_ext_val",
            "goof_rf_ext_val_df",
            # "goof_rf_mean_ext_val_df",
            "goof_rf_ext_val_metrics")

### Reload outputs
ds <- "DSM" # Manual
ds <- "HSM" # Manual

## In loop
mthd <- c("DSM", "HSM")
system.time(
for (ds in mthd) {
  print(ds)

### Reload outputs
propertie <- "C"
for (sv in seq_along(w_list)) {
  assign(w_list[[sv]], readRDS(
    paste0("../output/dataset/", ds, "/", ds, "_",
           propertie, "_", w_list[[sv]],".rds")))
}

### Assign the predicted property
j=21 # 21 == C

### Calculate the best model and coefficients
fst_qt_kcv <- goof_rf_cv_df[[j]][quantile(goof_rf_cv_df[[j]]$RMSE, 0.25) > goof_rf_cv_df[[j]]$RMSE, ]
best_model_kf_cv <- fst_qt_kcv[fst_qt_kcv$R2 == max(fst_qt_kcv$R2), "model_number"]
coef_best_model_kf_cv <- fst_qt_kcv[fst_qt_kcv$R2 == max(fst_qt_kcv$R2), ]

### Get the predicted property (C)
pred.property <- rf_models[[best_model_kf_cv]]$terms[[2]]
sp_spc <- as.character(pred.property)
names(sp_spc) <- as.character(pred.property)
names(sp_spc[1])

### Load predicted map
map <- terra::rast(paste0("../output/raster/predict_rasters_results/", ds, "/", ds, "_predict_", pred.property, ".tif"))
# plot(map)

crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# pni <- terra::vect("../input/shp/PARNA_Itatiaia_SIRGAS_IMG.shp", crs= crs_ref)
pni_alta <- terra::vect("../input/shp/parte alta.shp", crs= crs_ref)
# pnipts <- terra::vect("../input/shp/amostragem_cLHS_100.shp", crs= crs_ref)

### Convert train and val to points
df.mod <- read_fst("../output/dataset/dataset_HDSM_covariates_df_mod_point_pixels.fst")
df.ex.val <- read_fst("../output/dataset/dataset_Ext_val_covariates_df_external_validation.fst")
train_kfold_poins <- terra::vect(df.mod[, 192:193-10], geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)
ext_val_poins <- terra::vect(df.ex.val[, 41:42], geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)

################################################################################
### Plot TC map
################################################################################
ggplot() +
  layer_spatial(map, aes(alpha = stat(band1)), fill = "#8c510a") +
  scale_alpha_continuous(na.value = 0, name = bquote(paste("T",.(names(sp_spc[1])), " ", .(unit)))) +
  layer_spatial(pni_alta, aes(), color= "black", fill = "transparent", size=0.5) +
  # layer_spatial(pnipts, aes(), color= "black", size=0.5) +
  layer_spatial(train_kfold_poins, aes(), color= "black", size=0.5) + # it was "blue"
  layer_spatial(ext_val_poins, aes(), color= "black", size=0.5) + # it was "green"
  ggtitle("Itatiaia National Park", subtitle = paste0("Predicted TC with ", ds)) + xlab("") + ylab("") + 
  annotation_scale(location = "tl") + annotation_north_arrow(location = "br", which_north = "true") + 
  annotate("text", x = 525900, y = 7537500, label = c(paste("R²", round((data.frame(coef_best_model_kf_cv))[,1], 2), "\n",
                                                            # "MSE", round((data.frame(coef_best_model_kf_cv))[,3], 2), "\n",
                                                            "RMSE", round((data.frame(coef_best_model_kf_cv))[,4], 2), "\n",
                                                            # "bias", round((data.frame(coef_best_model_kf_cv))[,5], 2)
                                                            "Model", round((data.frame(coef_best_model_kf_cv))[,6], 2)
                                                            )), parse = F) + theme_light() + coord_sf(datum= 32723)

ggsave(paste0("../output/graphs_plots/", ds, "/", ds, "_", "Map_of_predicted_", names(sp_spc[1]), ".png"), 
       plot = last_plot(), width = 15, height = 15, units = "cm", dpi = 300)
# ggsave(paste0("../output/graphs_plots/", ds, "/", ds, "_", "Map_of_predicted_", names(sp_spc[1]), ".pdf"), 
# plot = last_plot(), width = 15, height = 15, units = "cm", dpi = 300)

### Export to overleaf folder
# pwd <- "C:/Users/Gel8695/Dropbox/Aplicativos/Overleaf/Gelsleichter_etal_HSM_2022_v2_Geoderma_Regional/figures/"
# ggsave(paste0(pwd, ds, "_", "Map_of_predicted_", names(sp_spc[1]), ".pdf"),
# plot = last_plot(), width = 18, height = 15, units = "cm", dpi = 300)


################################################################################
### Re-plot the IncNodePurity ### https://hackernoon.com/random-forest-regression-in-r-code-and-interpretation

library(randomForest)
### randomForest original plot
### (1=mean decrease in accuracy, 2=mean decrease in node impurity).
# varImpPlot(rf_models[[best_model_kf_cv]], type=1, n.var= 40)
# varImpPlot(rf_models[[best_model_kf_cv]], type=2, n.var= 40)

# ### Plot random forest varImpPlot and error model
# # varImpPlot(rf_models[[i]], main = "RF model - C, INP")
# # varImpPlot(rf_models[[i]], main = "RF model - C, INP", type=1)
# # varImpPlot(rf_models[[i]], main = "RF model - C, INP", n.var=20, type=1)
# varImpPlot(rf_models[[best_model_kf_cv]], type=1)
# varImpPlot(rf_models[[best_model_kf_cv]], type=2)
# 
# rf_path <- file.path(paste0("../output/graphs_plots/", ds, "/", ds, "_",  names(sp_spc[1]), "_Var_importance_", ".png"))
# png(file= rf_path, height=18, width=16, units="cm", res=200)
# varImpPlot(rf_models[[best_model_kf_cv]], main = paste0(ds, "\n", "Sentinel+Terrain+Geographic"), 
#            # n.var=22, 
#            type=2) ### (1=mean decrease in accuracy, 2=mean decrease in node impurity).
# dev.off()
# rf_path <- file.path(paste0("../output/graphs_plots/", ds, "/", ds, "_",  names(sp_spc[1]), "_Model_", ".png"))
# png(file= rf_path, height=22, width=18, units="cm", res=200)
# plot(rf_models[[best_model_kf_cv]], main = paste0("DSM"))
# dev.off()

### Plot with ggplot, much better. ### But always look on the random forest plot to have an idea of correct order, etc, ggplot can revert the order.

ImpData <- as.data.frame(importance(rf_models[[best_model_kf_cv]]))
ImpData$Var.Names <- row.names(ImpData)

ImpData %>%
  arrange(( IncNodePurity)) %>% # First sort by val. This sort the dataframe but NOT the factor levels
  slice_tail(n= 60) %>% ### get only the first 30 rows with higest values of IncNodePurity
  mutate(`%IncMSE`= round(`%IncMSE`, 0)) %>% # make the %IncMSE round to one digit
  # mutate(Var.Names = str_replace_all(Var.Names, "X", "HSI_")) %>%
  mutate(Var.Names = gsub("X", "H.", Var.Names)) %>%
  mutate(Var.Names = paste0(Var.Names, "  r", length(Var.Names):1)) %>%
  mutate(Var.Names= factor(Var.Names, levels= Var.Names)) %>%   # update the factor levels and order the plot

  ggplot(aes(x=Var.Names, y= IncNodePurity)) +
  geom_segment(aes(x=Var.Names, xend=Var.Names, y=0, yend= IncNodePurity), color="skyblue") +
  geom_point(aes(size = `%IncMSE`), color="blue", alpha= 0.7) +
  xlab("Covariates") + 
  ggtitle(paste0(ds)) + # ggtitle(paste0(ds, "\n", "Sentinel+Terrain+Geographic")) + 
  theme_light() + coord_flip() + #scale_x_reverse() +
  theme(legend.position="bottom",
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size= 16))

w <- if (ds == "DSM") {17} else 17 # width of the plot 
h <- if (ds == "DSM") {12} else 44 # height of the plot 
ggsave(paste0("../output/graphs_plots/", ds, "/", ds, "_", "IncNodePuryty_IncMSE_", 
              names(sp_spc[1]), ".png"), plot = last_plot(), width = w, height = h, units = "cm", dpi = 300)
ggsave(paste0("../output/graphs_plots/", ds, "/", ds, "_", "IncNodePuryty_IncMSE_", 
              names(sp_spc[1]), ".pdf"), plot = last_plot(), width = w, height = h, units = "cm", dpi = 300)
### Export to overleaf folder
# pwd <- "C:/Users/Gel8695/Dropbox/Aplicativos/Overleaf/Gelsleichter_etal_HSM_2022_v2_Geoderma_Regional/figures/"
# ggsave(paste0(pwd, ds, "_", "IncNodePuryty_IncMSE_", names(sp_spc[1]), ".pdf"),
#        plot = last_plot(), width = w, height = h, units = "cm", dpi = 300)

}) # close printing loop

################################################################################
### Boxplot with the outputs
################################################################################
################################################################################
### Load data
{
ds <- "DSM"
pf <- paste0("../output/dataset/", ds, "/")
fllist <- list.files(path = pf, pattern='.csv$', all.files= T, full.names= T)

DSM_calib_metrics <- read_csv(fllist[1], n_max= 8)
DSM_kfcv_metrics <- read_csv(fllist[2], n_max= 8)
DSM_ext.val_metrics <- read_csv(fllist[3], n_max= 8)

ds <- "HSM"
pf <- paste0("../output/dataset/", ds, "/")
fllist <- list.files(path = pf, pattern='.csv$', all.files= T, full.names= T)

HSM_calib_metrics <- read_csv(fllist[1], n_max= 8)
HSM_kfcv_metrics <- read_csv(fllist[2], n_max= 8)
HSM_ext.val_metrics <- read_csv(fllist[3], n_max= 8)

DSM_calib_metrics$Method <- "DSM"
DSM_kfcv_metrics$Method <- "DSM"
DSM_ext.val_metrics$Method <- "DSM"
HSM_calib_metrics$Method <- "HSM"
HSM_kfcv_metrics$Method <- "HSM"
HSM_ext.val_metrics$Method <- "HSM"

DSM_calib_metrics$Metric <- "Calibration"
DSM_kfcv_metrics$Metric <- "8-fold CV"
DSM_ext.val_metrics$Metric <- "EV"
HSM_calib_metrics$Metric <- "Calibration"
HSM_kfcv_metrics$Metric <- "8-fold CV"
HSM_ext.val_metrics$Metric <- "EV"

df <- rbind(DSM_calib_metrics,
            DSM_kfcv_metrics,
            DSM_ext.val_metrics,
            HSM_calib_metrics,
            # HSM_kfcv_metrics[, 2:ncol(HSM_kfcv_metrics)], ### strange
            HSM_kfcv_metrics, ### the above one was commented and this one applied
            HSM_ext.val_metrics)

# df$model_number <- NULL
names(df)[6] <- "Fold"
df$Fold <- as.factor(df$Fold)

dff <- as.data.frame(df)
# mdata <- reshape::melt(dff[, c(1,4,6,7)], id= c("Method", "Metric"))
# ggplot(mdata, aes(x = Metric, y= value, fill= Method)) +
#   geom_boxplot(alpha = 0.80) + 
#   facet_grid(~ variable)
  
dff$Metric <- factor(dff$Metric, levels=c("Calibration", "8-fold CV", "EV"))
dff$Metric %>% levels()

# R2 <- reshape::melt(dff[, c(1,6,7)], id= c("Method", "Metric"))
# MSE <- reshape::melt(dff[, c(3,6,7)], id= c("Method", "Metric"))
# RMSE <- reshape::melt(dff[, c(4,6,7)], id= c("Method", "Metric"))

R2 <- reshape::melt(dff[, c(1,6,7,8)], id= c("Method", "Metric", "Fold"))
MSE <- reshape::melt(dff[, c(3,6,7,8)], id= c("Method", "Metric", "Fold"))
RMSE <- reshape::melt(dff[, c(4,6,7,8)], id= c("Method", "Metric", "Fold"))

runlist_m3 <- list(R2, MSE, RMSE)

}

### Select the metric(s) for plot
## as the line bellow, but in loop
## runlist_m3[[1]] <- filter(runlist_m3[[1]], Metric != "Calibration")
runlist <- list()
for (mtr in 1:3) {
runlist[[mtr]] <- filter(runlist_m3[[mtr]], Metric != "Calibration")
# runlist[[mtr]] <- filter(runlist_m3[[mtr]], Metric == "EV")
}


plt <- list()
for (p in 1:3) {
  (
    plt[[p]] <- ggplot(runlist[[p]], aes(x = Metric, y= value, fill= Method)) +
      geom_boxplot(alpha = 0.80) + facet_grid(~ variable) +
      theme_ipsum(axis_text_size  = 12,
                  strip_text_size   = 18) +
      geom_point(aes(fill = Method), size = 1.5, shape = 1, position = position_dodge(1)) +
      scale_fill_viridis(discrete = T, alpha= 0.6) +
      stat_summary(fun= mean, geom="point", shape=23, size= 2, ### place the mean in the boxplot
      position=position_dodge(.75)) +
      theme(#text = element_text(size = 18),
        axis.text.x = element_text(angle = 30, vjust= 1, hjust= 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # panel.grid.minor.x = element_blank(),
        # panel.grid.major.x = element_blank(), 
        # get the legend only in the last one
        legend.position = 
          if (p != 3) {
            paste0("none") # paste "none" when p=1 and p=2
          } else 
            paste0(c(0.94,-0.2)) # paste c(0.94,-0.4) when p=3
      ) + 
      if (p != 1) {
        # scale_y_continuous() + # no action need
      } else 
        scale_y_continuous(breaks=seq(0, 1, 0.1)) 
  )
}
plt[[1]]
plt[[2]]
plt[[3]]

library(gridExtra)
pr <- grid.arrange(plt[[1]],
                   # plt[[2]],
                   plt[[3]],
                   nrow = 1)
pr
ggsave(paste0("../output/graphs_plots/", 
              "Boxplot_metrics", ".png"), 
       plot = pr, width = 22, height = 15, units = "cm", dpi = 300)

plt <- list()
for (p in 1:3) {
  (
    plt[[p]] <- ggplot(runlist[[p]], aes(x = Metric, y= value, fill= Method)) +
      geom_boxplot(alpha = 0.80) + facet_grid(~ variable) +
      theme_ipsum(base_family = "sans",
                  axis_text_size  = 12,
                  strip_text_size   = 14) +
      geom_point(aes(fill = Method), size = 1.5, shape = 1, position = position_dodge(1)) +
      scale_fill_viridis(discrete = T, alpha= 0.4) +
      stat_summary(fun= mean, geom="point", shape= 23, size= 2, ### place the mean in the boxplot
                   position=position_dodge(.75)) + # mean
      # geom_label_repel(aes(label = ifelse(Metric == "8-fold CV", Fold,'')), # ifelse makes the label only in the "8-fold CV"
      #                  size= 3,
      #                  # direction = "x",
      #                  min.segment.length = 0, # Draw all line segments
      #                  box.padding   = 0.9, # Length of line
      #                  label.padding = 0.15, # box size
      #                  point.padding = 0, 
      #                  position = position_dodge(1)#,
      #                  # colour="black", 
      #                  # segment.colour="black"
      #                  ) +
      theme(#text= element_text(size= 16,  family="sans"),
            axis.text.x = 
              if (p != 3) {
                (element_blank())
              } else 
                (element_text(angle = 30, vjust= 1, hjust= 1)),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = # Adjust margins inside ggplot
              if (p != 3) {
                (unit(c(1,1,0.5,1), "cm")) # top, right, bottom, and left
              } else 
                (unit(c(1,1,1,1.35), "cm")),
            # get the legend only in the last one
            legend.position = 
              if (p != 3) {
                paste0("none") # paste "none" when p=1 and p=2
              } else 
                paste0(c(0.94,-0.2)) # paste c(0.94,-0.4) when p=3
      ) + 
      if (p != 1) {
        # scale_y_continuous() + # no action need
      } else 
        scale_y_continuous(breaks= seq(0, 1, 0.1)) 
  )
}
plt[[1]]
plt[[2]]
plt[[3]]

# library(gridExtra)
# pr <- arrangeGrob(plt[[1]],
#                   # plt[[2]],
#                   plt[[3]],
#                   ncol = 1)
# ggsave(paste0("../output/graphs_plots/", 
#               "Boxplot_metrics_2", ".png"), 
#        plot = pr, width = 12, height = 25, units = "cm", dpi = 300)
# 
# ### Export to overleaf folder
# pwd <- "C:/Users/Gel8695/Dropbox/Aplicativos/Overleaf/Gelsleichter_etal_HSM_2022_v2_Geoderma_Regional/figures/"
# ggsave(paste0(pwd, "Boxplot_metrics_2.pdf"), # does not accept pdf
#        plot = pr, width = 12, height = 25, units = "cm", dpi = 300)

############################################## 
### With 8-fold CV and EV in the same scale

plt <- list()
for (p in 1:3) {
  (
    plt[[p]] <- ggplot(runlist[[p]], aes(x = Metric, y= value, fill= Method)) +
      geom_boxplot(alpha = 0.80) + facet_grid(~ variable) + # facet_grid adds the name to the top
      theme_ipsum(base_family = "sans",
                  axis_text_size  = 12,
                  strip_text_size   = 14) +
      geom_point(aes(fill = Method), size = 1.5, shape = 1, position = position_dodge(1), col= "transparent") +
      scale_fill_viridis(discrete = T, alpha= 0.4) +
      stat_summary(fun= mean, geom="point", shape= 23, size= 2, ### place the mean in the boxplot
                   position=position_dodge(.75)) + # mean
      theme(#text= element_text(size= 16,  family="sans"),
        axis.text.x = 
          if (p != 3) {
            (element_blank())
          } else 
            (element_text(angle = 30, vjust= 1, hjust= 1)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = # Adjust margins inside ggplot
          if (p != 3) {
            (unit(c(1,1,0.5,1), "cm")) # top, right, bottom, and left
          } else 
            (unit(c(1,1,1,1.35), "cm")),
        # get the legend only in the last one
        legend.position = 
          if (p != 3) {
            paste0("none") # paste "none" when p=1 and p=2
          } else 
            paste0(c(0.94,-0.2)) # paste c(0.94,-0.4) when p=3
      ) + 
      if (p != 1) {
        # scale_y_continuous() + # no action need
      } else 
        scale_y_continuous(breaks= seq(0, 1, 0.1)) 
  )
}

### Add the labels later (after plot), then labels are over the figure
# geom_label_repel(aes(label = ifelse(Metric == "8-fold CV", Fold,'')), # ifelse makes the label only in the "8-fold CV"
gr <- geom_label_repel(aes(label = Fold), # ifelse makes the label only in the "8-fold CV"
                       size= 3.5,
                       direction = "x",
                       min.segment.length = 0, # Draw all line segments
                       box.padding   = 0.1, # Length of line
                       label.padding = 0.18, # box size
                       point.padding = 0,
                       seed = 123,
                       alpha = 0.85,
                       fontface = 'bold',
                       show.legend = F,
                       # fill= "red", 
                       max.overlaps = 20, 
                       position = position_dodge(1),
                       colour="royalblue2",
                       segment.colour= NA, #"transparent", red, blue
) 
plt[[3]] + gr

for (cl in 1:3) {
  plt[[cl]] <- plt[[cl]] + gr
}

library(gridExtra)
# pr <- grid.arrange(plt[[1]], # to look
#                    # plt[[2]],
#                    plt[[3]],
#                    ncol = 1)
pr <- arrangeGrob(plt[[1]], # to export ligther
                  # plt[[2]],
                  plt[[3]],
                  ncol = 1)

### Adjust margins in the grid.arrange function (smaller, closer plots)
# margin = theme(plot.margin = unit(c(-0.7, 1, 1, 1) , "cm")) #  top, right, bottom, and left 
# pr <- grid.arrange(grobs= lapply(plt, "+", margin), ncol = 1)

# library(grid) # for viewport
# pr <- grid.arrange(grobs= lapply(list(plt[[1]],
#                                       # plt[[2]],
#                                       plt[[3]]), "+", margin),  
#                    vp= viewport(width= 1.0, height=0.9), ncol = 1)
# pr <- arrangeGrob(grobs= lapply(list(plt[[1]],
#                                       # plt[[2]],
#                                       plt[[3]]), "+", margin),  
#                    vp= viewport(width= 1.0, height=0.9), ncol = 1)
ggsave(paste0("../output/graphs_plots/", 
              "Boxplot_metrics_2", ".png"), 
       plot = pr, width = 12, height = 25, units = "cm", dpi = 300)

### Export to overleaf folder
# extrafont::font_import() ### to save in PDF, takes a lot of time ~20 min
# extrafont::loadfonts()
# extrafont::font_import(pattern = 'ARIAL')
# extrafont::loadfonts(device= "win")
# pwd <- "C:/Users/Gel8695/Dropbox/Aplicativos/Overleaf/Gelsleichter_etal_HSM_2022_v2_Geoderma_Regional/figures/"
# ggsave(paste0(pwd, "Boxplot_metrics_2.pdf"), # does not accept pdf
#        plot = pr, width = 12, height = 25, units = "cm", dpi = 300)

############################################## 
### With 8-fold CV and EV in each scale


###############################################
### R2, RMSE for 8-fold CV

### Select the metric(s) for plot
runlist <- list()
for (mtr in 1:3) {
  runlist[[mtr]] <- filter(runlist_m3[[mtr]], Metric == "8-fold CV")
}

plt <- list()
for (p in 1:3) {
  (
    plt[[p]] <- ggplot(runlist[[p]], aes(x = Metric, y= value, fill= Method)) +
      geom_boxplot(alpha = 0.80) + # facet_grid(~ variable) + # facet_grid adds the name to the top
      theme_ipsum(base_family = "sans",
                  axis_text_size  = 12,
                  strip_text_size   = 14) +
      geom_point(aes(fill = Method), size = 1.5, shape = 1, position = position_dodge(1), col= "transparent") +
      scale_fill_viridis(discrete = T, alpha= 0.4) +
      stat_summary(fun= mean, geom="point", shape= 23, size= 2, ### place the mean in the boxplot
                   position=position_dodge(.75)) + # mean
      # ylab(as.character(unique(runlist[[p]]$variable))) +
      ylab(if (as.character(unique(runlist[[p]]$variable)) == "R2") {
        bquote(R^2)
        } else as.character(unique(runlist[[p]]$variable)) ) +
    theme(#text= element_text(size= 16,  family="sans"),
      axis.text.x = 
        if (p != 3) {
          (element_blank())
        } else 
          (element_text(angle = 30, vjust= 1, hjust= 1)),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size= 14),
      # axis.title.y = element_blank(),
      plot.margin = # Adjust margins inside ggplot
        if (p != 3) {
          (unit(c(1,1,0.5,1), "cm")) # top, right, bottom, and left
        } else 
          (unit(c(1,1,1,1.35), "cm")),
      # get the legend only in the last one
      legend.position = 
        # if (p != 3) {
          paste0("none") # paste "none" when p=1 and p=2
        # } else 
          # paste0(c(0.94,-0.2)) # paste c(0.94,-0.4) when p=3
    ) + 
      if (p != 1) {
        # scale_y_continuous() + # no action need
      } else 
        scale_y_continuous(breaks= seq(0, 1, 0.1)) 
  )
}

### Add the labels later (after plot), then labels are over the figure
# geom_label_repel(aes(label = ifelse(Metric == "8-fold CV", Fold,'')), # ifelse makes the label only in the "8-fold CV"
gr <- geom_label_repel(aes(label = Fold), # ifelse makes the label only in the "8-fold CV"
                 size= 3.5,
                 direction = "x",
                 min.segment.length = 0, # Draw all line segments
                 box.padding   = 0.1, # Length of line
                 label.padding = 0.18, # box size
                 point.padding = 0,
                 seed = 123,
                 alpha = 0.85,
                 fontface = 'bold',
                 show.legend = F,
                 # fill= "red", 
                 max.overlaps = 20, 
                 position = position_dodge(1),
                 colour="royalblue2",
                 segment.colour= NA, #"transparent", red, blue
) 
plt[[1]] + gr
plt[[3]] + gr

for (cl in 1:3) {
plt[[cl]] <- plt[[cl]] + gr
}

###############################################
### R2, RMSE for EV

### Select the metric(s) for plot
runlist <- list()
for (mtr in 1:3) {
  runlist[[mtr]] <- filter(runlist_m3[[mtr]], Metric == "EV")
}

plt2 <- list()
for (p in 1:3) {
  (
    plt2[[p]] <- ggplot(runlist[[p]], aes(x = Metric, y= value, fill= Method)) +
      geom_boxplot(alpha = 0.80) + # facet_grid(~ variable) + # facet_grid adds the name to the top
      theme_ipsum(base_family = "sans",
                  axis_text_size  = 12,
                  strip_text_size   = 14) +
      geom_point(aes(fill = Method), size = 1.5, shape = 1, position = position_dodge(1), col= "transparent") +
      scale_fill_viridis(discrete = T, alpha= 0.4) +
      stat_summary(fun= mean, geom="point", shape= 23, size= 2, ### place the mean in the boxplot
                   position=position_dodge(.75)) + # mean
      theme(#text= element_text(size= 16,  family="sans"),
        axis.text.x = 
          if (p != 3) {
            (element_blank())
          } else 
            (element_text(angle = 30, vjust= 1, hjust= 1)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = # Adjust margins inside ggplot
          if (p != 3) {
            (unit(c(1,1,0.5,1), "cm")) # top, right, bottom, and left
          } else 
            (unit(c(1,1,1,1), "cm")), # changed from 1.35 to 1
        # get the legend only in the last one
        legend.position = 
          if (p != 3) {
            paste0("none") # paste "none" when p=1 and p=2
          } else 
            paste0(c(0.94, -0.08)) # paste c(0.94,-0.4) when p=3
      ) + 
      if (p != 1) {
        scale_y_continuous() # no action need
      } else 
        scale_y_continuous(breaks= seq(0, 1, 0.05)) 
  )
}

plt2[[3]] + gr

for (cl in 1:3) {
  plt2[[cl]] <- plt2[[cl]] + gr
}
library(gridExtra)
pr <- grid.arrange(plt[[1]], # to look
                   plt2[[1]],
                   plt[[3]],
                   plt2[[3]],
                   ncol = 2)
pr <- arrangeGrob(plt[[1]], # to export (ligther for ggsave)
                  plt2[[1]],
                  plt[[3]],
                  plt2[[3]],
                  ncol = 2)

### Adjust margins in the grid.arrange function (smaller, closer plots)
ggsave(paste0("../output/graphs_plots/", 
              "Boxplot_metrics_4_4", ".png"), 
       plot = pr, width = 16, height = 26, units = "cm", dpi = 300)

### Export to overleaf folder
# pwd <- "C:/Users/Gel8695/Dropbox/Aplicativos/Overleaf/Gelsleichter_etal_HSM_2022_v2_Geoderma_Regional/figures/"
# ggsave(paste0(pwd, "Boxplot_metrics_4_4.pdf"), # does not accept pdf
#        plot = pr, width = 16, height = 26, units = "cm", dpi = 300)


################################################################################
################################################################################
### plot the folds data in map with C
################################################################################
library(viridis)
library(RColorBrewer)
library(fst)
library(raster)
### Convert train and val to points
stk_dsm <- terra::rast("../output/raster/stk_dsm_10m.tif")
map <- terra::rast(paste0("../output/raster/predict_rasters_results/", ds, "/", ds, "_predict_", pred.property, ".tif"))

crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
pni_alta <- terra::vect("../input/shp/parte alta.shp", crs= crs_ref)

df.mod <- read_fst("../output/dataset/dataset_HDSM_covariates_df_mod_point_pixels.fst")
df.ex.val <- read_fst("../output/dataset/dataset_Ext_val_covariates_df_external_validation.fst")
train_kfold_poins <- terra::vect(df.mod[, 192:193-10], geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)
ext_val_poins <- terra::vect(df.ex.val[, 41:42], geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)

names(df.mod)
tt <- df.mod[, c(21, 182:183)]
### Split data into 10 folds
folds <- cut(seq(1,nrow(tt)), breaks= 8, labels=FALSE)
nflds <- length(unique(folds))
# table(folds)

for(i in 1:nflds){ ### This loop is to run and validate with the K-folds
  ### Segment your data by fold using the which() function 
  testIndexes <- which(folds== i, arr.ind=TRUE)
  testData <- tt[testIndexes, ]
  trainData <- tt[-testIndexes, ]
  test_poins <- terra::vect(testData, geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)
  train_poins <- terra::vect(trainData, geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)
  test_poins$Datatype <- "Test"
  train_poins$Datatype <- "Train"
  ext_val_poins$Datatype <- "Ext. validation"
  trts <- terra::union(test_poins, train_poins)
  shp <- terra::union(trts, ext_val_poins)
  
ggplot() +
  layer_spatial(map, aes(alpha = stat(band1)), fill = "#8c510a") +
  scale_alpha_continuous(na.value = 0, name = bquote(paste(.(names(sp_spc[1])), " ", .(unit)))) +
  layer_spatial(pni_alta, aes(), color= "black", fill = "transparent", size=0.5) +
  layer_spatial(shp, aes(col= Datatype), size=1.2, shape= c(1)) +
  scale_colour_manual(values = c("pink", "green2", "blue4")) + 
  ggtitle("Datasets (C)", subtitle = paste0("Fold: ", i)) + xlab("") + ylab("") + 
  annotation_scale(location = "tl") + annotation_north_arrow(location = "br", which_north = "true") + 
  # annotate("text", x = 525900, y = 7537500, label = c(paste("R²", round((data.frame(coef_best_model_kf_cv))[,1], 2), "\n",
                                                            # "MSE", round((data.frame(coef_best_model_kf_cv))[,3], 2), "\n",
                                                            # "RMSE", round((data.frame(coef_best_model_kf_cv))[,4], 2), "\n",
                                                            # "bias", round((data.frame(coef_best_model_kf_cv))[,5], 2))), parse = F) + 
theme_light() + coord_sf(datum= 32723)
ggsave(paste0("../output/graphs_plots/", 
              "Train_test_data_ext_val_Spatial_distribution/", 
              "C_Train_Test_Val_Datasets_fold_", i, ".png"), 
       plot = last_plot(), width = 15, height = 15, units = "cm", dpi = 300)
}

################################################################################
###############
###############
################################################################################
### Plot folds data train, test and valida with DEM in ggplot
################################################################################
###############
###############
################################################################################
library(raster)
library(viridis) ### to the color scale
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(RStoolbox)
library(ggspatial) ### for north arrow and scale

### Load data and calculate the hill
MDE <- stk_dsm[[1]]
MDE <- terra::crop(MDE, pni_alta, snap="out", mask= T)
MDE <- raster(MDE)
terr <- raster::terrain(MDE, c("slope", "aspect"))
hill <- raster::hillShade(terr[["slope"]], terr[["aspect"]])
ggR(MDE)
ggR(terr)
ggR(hill)
plot(hill)

### Load TC and Coordinates
names(df.mod)
C_sp_df <- df.mod[, c(21,182,183)]
# Make int
C_sp_df$C <- as.integer(C_sp_df$C)
names(C_sp_df) <- c("TC", "X", "Y")
# Change to number of bins you want
# C_sp_df$breaks <- cut(C_sp_df$TC, 5) 

summary(C_sp_df) # min 3 max 27 (limits)

# Convert from SpatVector to SpatialPolygonsDataFrame 
pnisp <- as(pni_alta, "Spatial")
tt_points <- as(shp, "Spatial")
tt_points <- fortify(as.data.frame(tt_points)) 

# Colors
# crp <- colorRampPalette(rev(brewer.pal(9,'YlOrBr')))
myPalette <- colorRampPalette((brewer.pal(9,'YlOrBr')))

for(i in 1:nflds){ ### This loop is to run and validate with the K-folds
  ### Segment your data by fold using the which() function 
  testIndexes <- which(folds== i, arr.ind=TRUE)
  testData <- tt[testIndexes, ]
  trainData <- tt[-testIndexes, ]
  test_poins <- terra::vect(testData, geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)
  train_poins <- terra::vect(trainData, geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)
  test_poins$Datatype <- "Test"
  train_poins$Datatype <- "Train"
  ext_val_poins$Datatype <- "Ext. validation"
  trts <- terra::union(test_poins, train_poins)
  shp <- terra::union(trts, ext_val_poins)
  
  ggR(hill, forceCat = F) +
    geom_raster(data= MDE, aes(x=x, y=y, fill= DEM), alpha=0.6) +
    scale_fill_gradient(low = "green4", high = "firebrick4", na.value = NA, name = "Elevation (m)", guide = "legend") +
    layer_spatial(shp, aes(col = Datatype), size= 1.2, shape= c(1)) +
    scale_color_manual(values = c("pink1", "green", "blue4")) + 
    geom_polygon(data = pnisp, aes(x = long, y = lat, group = group), colour = "black", size = 0.5, fill = NA) + 
    ggtitle("Datasets (C)", subtitle = paste0("Fold: ", i)) + xlab("") + ylab("") + 
    theme_map() +
    annotation_scale(location = "tl") +
    annotation_north_arrow(location = "br", which_north = "true") + theme_classic() + #theme_bw() + 
    theme(legend.position="right", legend.key.width=unit(0.4, "cm"), ### position and thickness of scale bar
          axis.title.x=element_blank(), ### remove the axis tittle 
          axis.title.y=element_blank(),
    ) # + coord_sf(crs = 4326), + coord_sf(datum= 32723) 
  
  ggsave(paste0("../output/graphs_plots/", 
                "Train_test_data_ext_val_Spatial_distribution/", 
                "INP_hillshape_TTV_fold_", i, ".png"), 
         plot = last_plot(), width = 15, height = 15, units = "cm", dpi = 300)
}


################################################################################
################### plot the folds data with DEM
################################################################################

### Split data into 10 folds
folds <- cut(seq(1,nrow(tt)), breaks= 8, labels=FALSE)
nflds <- length(unique(folds))

tt$Datatype <- as.factor(paste0("Fold ", folds))

df.ex.val <- read_fst("../output/dataset/dataset_Ext_val_covariates_df_external_validation.fst")
df.ex.val$id
df.ex <- df.ex.val[1:12, c(21,41:42)]
df.ex_legacy <- df.ex.val[13:18, c(21,41:42)]

df.ex$Datatype <- as.factor(" \nExternal \nvalidation \n ")
df.ex_legacy$Datatype <- as.factor("Legacy \nexternal \nvalidation")

ttdf <- rbind(tt, df.ex, df.ex_legacy)

### Convert to SpatVector
shp <- terra::vect(ttdf, geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)

# Convert from SpatVector to SpatialPolygonsDataFrame 
pnisp <- as(pni_alta, "Spatial")
tt_points <- as(shp, "Spatial")
tt_points <- fortify(as.data.frame(tt_points)) 

ggR(hill, forceCat = F) +
  geom_raster(data= MDE, aes(x=x, y=y, fill= DEM), alpha=0.6) +
  # scale_fill_gradient(low = "darkgreen", high = "darkred", na.value = NA, name = "Elevation (m)", guide = "legend") +
  scale_fill_gradient(low = "blue2", high = "goldenrod1", na.value = NA, name = "Elevation (m)", guide = "legend") +
  layer_spatial(shp, aes(col = Datatype, shape= Datatype, size= Datatype)) +
  scale_color_manual(values = c("red", "green", "cyan",
                                "blue", "yellow", "darkviolet",
                                "magenta", "gray10", "greenyellow", "lawngreen")) +
  scale_shape_manual(values=c(1,1,1,1,1,1,1,1,0,6)) + 
  scale_size_manual(values=c(1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 2.5, 2.5)) + 
  geom_polygon(data = pnisp, aes(x = long, y = lat, group = group), colour = "black", size = 0.5, fill = NA) + 
  ggtitle("Datasets (TC)") + xlab("") + ylab("") + 
  theme_map() +
  annotation_scale(location = "tl") +
  annotation_north_arrow(location = "br", which_north = "true") + theme_classic() + #theme_bw() + 
  theme(legend.position="right", legend.key.width=unit(0.4, "cm"), ### position and thickness of scale bar
        axis.title.x=element_blank(), ### remove the axis tittle 
        axis.title.y=element_blank()) + coord_sf(datum= 32723) # + coord_sf(crs = 4326)

# ggsave(paste0("../output/graphs_plots/", 
#               "Folds_valida_data_DEM.png"), 
#        plot = last_plot(), width = 18, height = 15, units = "cm", dpi = 300)
ggsave(paste0("../output/graphs_plots/", 
              "Folds_valida_data_DEM.pdf"), 
       plot = last_plot(), width = 18, height = 15, units = "cm", dpi = 300)

### Export to overleaf folder
# pwd <- "C:/Users/Gel8695/Dropbox/Aplicativos/Overleaf/Gelsleichter_etal_HSM_2022_v2_Geoderma_Regional/figures/"
# ggsave(paste0(pwd, "Folds_valida_data_DEM.pdf"), 
#        plot = last_plot(), width = 18, height = 15, units = "cm", dpi = 300)

################################################################################
################### plot the folds data as Histograms
names(df.mod)
# tt <- df.mod[, c(170:181)]
tt <- df.mod
### Split data into 10 folds
folds <- cut(seq(1,nrow(tt)), breaks= 8, labels=FALSE)
nflds <- length(unique(folds))
# table(folds)
class(tt[[180]])
tt[[180]] <- as.integer(tt[[180]])
tt[[181]] <- as.integer(tt[[181]])
names(tt)[180]
names(tt)[181]

for (w in c(21, 170:181)) {
# save plots to png
rf_path <- file.path(paste0("../output/graphs_plots/Histograms/", "Histogram_Train_Test_Val_Datasets_propertie_", names(tt)[w], ".png"))
png(file= rf_path, height=45, width=20, units="cm", res=200)
par(mfrow = c(8, 2))
### Segment your data by fold using the which() function 
### This loop is to run and validate with the K-folds
for(i in 1:nflds){ 
  testIndexes <- which(folds== i, arr.ind=TRUE)
  testData <- tt[testIndexes, ]
  trainData <- tt[-testIndexes, ]
  summary(testData[w])
  summary(trainData[w])
  h <- hist(trainData[[w]], main= paste0("Train Data fold: ", i), 
            col= i, xlab = names(trainData[w]))
  x <- trainData[[w]]
  xfit <- seq(min(x), max(x), length= 100)
  yfit <- dnorm(xfit, mean= mean(x), sd=sd(x))
  yfit <- yfit*diff(h$mids[1:2])*length(x)
  lines(xfit, yfit, col="blue", lwd=2)
  
  h2 <- hist(testData[[w]], main= paste0("Test Data fold: ", i), 
             col= i, xlab = names(trainData[w]))
  x2 <- testData[[w]]
  xfit <- seq(min(x2), max(x2), length= 100)
  yfit <- dnorm(xfit, mean= mean(x2), sd=sd(x2))
  yfit <- yfit*diff(h2$mids[1:2])*length(x2)
  lines(xfit, yfit, col="blue", lwd=2)
}
# turn off png plotting
dev.off()
}

################################################################################
### Bar plot 
ttc <- data.frame(df.mod[21])
### Split data into 10 folds
folds <- cut(seq(1,nrow(ttc)), breaks= 8, labels=FALSE)
nflds <- length(unique(folds))

ttc$folds_num <- as.factor(folds)
ttc$ids_num <- 1:nrow(ttc)

head(ttc)

# Plot
ggplot(ttc, aes(x= ids_num, y= C, fill= folds_num)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_fill_viridis(discrete = T, direction = -1) 
  
################################################################################
### Covariates bar plots
library(dplyr)

### Split data into 10 folds
names(df.mod)
ttc <- data.frame(df.mod)
folds <- cut(seq(1,nrow(ttc)), breaks= 8, labels=FALSE)
nflds <- length(unique(folds))

names(df.mod[, c(21,170:180)])
# k=21 # TC
# k=180 # geology
# k=170 # DEM
# summary(df.mod[[k]])

### pl_list <- list() # list to take away and plot together

# for (k in 170:181) {
for (k in c(21,170:180)) {
ttc <- data.frame(df.mod[k])
ttc[[1]] <- as.numeric(ttc[[1]]) # due geology and geomofology
ttc$folds_num <- as.factor(folds)
# ttc$folds_num <- as.character(ttc$folds_num)

ttc <- ttc %>% group_by(ttc[[1]]) %>% arrange(folds_num, ttc[[1]]) 
ttc <- ttc[1:2]
ttc$ids_num <- 1:nrow(ttc)

# Plot
{
p1 <- ggplot(ttc, aes(x= ids_num, y= ttc[[1]], fill= folds_num)) + 
    geom_bar(position="dodge", stat="identity") +
    xlab("Sample points") + ylab(names(df.mod[k])) + labs(fill = "Folds") +
    scale_fill_viridis(discrete = T, direction = -1) + 
    theme_classic() + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) + 
    if (names(ttc[1]) != "DEM") {
      coord_cartesian()
      } else 
        coord_cartesian(ylim=c(1480, 2500)) 
}
p1 <- p1 + if (names(ttc[1]) == "C") {ylab("TC (%)")
  } else 
    if (names(ttc[1]) == "DEM") {ylab("DEM (m)")
      } else 
        ylab(names(ttc[1]))
# p1

### ### take away to plot together
### if (names(ttc[1]) == "C") {pl_list[[k]] <- p1
### } else 
###   if (names(ttc[1]) == "DEM") {pl_list[[k]] <- p1
###   } else p1 <- p1

ggsave(paste0("../output/graphs_plots/", 
              "Barplot_covariates/", 
              "Barplot_covariates_property_", names(df.mod[k]), ".png"), 
       plot = p1, width = 15, height = 15, units = "cm", dpi = 300)

### Export to overleaf folder (but only TC or DEM)
if (any(names(ttc[1]) == c("DEM","C"))) {
pwd <- "C:/Users/Gel8695/Dropbox/Aplicativos/Overleaf/Gelsleichter_etal_HSM_2022_v2_Geoderma_Regional/figures/"

### loop to export into pdf and png
for (filetype in c(".png", ".pdf")) {
ggsave(paste0(pwd, "Barplot_covariates_property_", names(df.mod[k]), filetype), 
       plot = p1, width = 15, height = 15, units = "cm", dpi = 300)
}
}

} # close plot loop

### Ply with boxplot
# class(ttc$ids_num)
# ttc$ids_num <- as.factor(ttc$ids_num)
# ttc$folds_num <- as.factor(ttc$folds_num)
# ggplot(ttc, aes(x= ids_num, y= C, fill= folds_num)) +
#   # geom_bar(position="dodge", stat="identity") +
#   geom_boxplot() + # width = .4 # https://stackoverflow.com/questions/55816951/cants-seem-to-adjust-boxplot-width
#   # stat_boxplot(geom = "errorbar", width = 0.1) +
#   xlab("Sample points") + ylab(names(df.mod[k])) + labs(fill = "Folds") +
#   scale_fill_viridis(discrete = T, direction = -1) +
#   theme_classic() +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank()) 

################################################################################
### Plot spectral predicted
################################################################################
{
  library(terra)
  stk_dsm <- terra::rast("../output/raster/stk_dsm_10m.tif")
  ### Load every predicted spectral image (as .tif)
  ds <- "Spect_subsuf"
  pf <- paste0("../output/raster/predict_rasters_results/", ds, "/")
  rastlist <- list.files(path = pf, pattern='.tif$', all.files= T, full.names= T)
  spec_stk <- terra::rast(rastlist)
  # plot(spec_stk[[120]])
  
  stk_hsm <- c(spec_stk, stk_dsm)
  names(stk_hsm)[1:130] <- names(df.mod[40:(130+39)])
  
  crs_ref <- "+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  pni_alta <- terra::vect("../input/shp/parte alta.shp", crs= crs_ref)
  
  names(stk_hsm)
  # names(stk_hsm[[3]])
  # names(stk_hsm[[6]])
  # names(stk_hsm[[17]])
  # names(stk_hsm[[c(3,6,17)]])
  # cm <- stk_hsm[[c(3,6,17)]]
  # cm <- stk_hsm[[c(3,4,5)]]
  # cm <- stk_hsm[[c(1,4,17)]]
  # cm <- stk_hsm[[c(1,4,20)]]
  # cm <- stk_hsm[[c(1,4,125)]]
  cm <- stk_hsm[[c(17,4,1)]]
  # cm <- stk_hsm[[c(5,4,3)]]
  # cm <- stk_hsm[[c(6,5,4)]]

  ### Convert train and val to points
  df.mod <- read_fst("../output/dataset/dataset_HDSM_covariates_df_mod_point_pixels.fst")
  df.ex.val <- read_fst("../output/dataset/dataset_Ext_val_covariates_df_external_validation.fst")
  train_kfold_poins <- terra::vect(df.mod[, 192:193-10], geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)
  ext_val_poins <- terra::vect(df.ex.val[, 41:42], geom=c("x", "y"), crs= crs_ref, keepgeom=FALSE)
}

names(cm)
# RGB
library(ggplot2); library(RStoolbox); library(ggspatial)
ggplot() +
  # ggRGB(img = cm, r = 1, g = 2, b = 3 , stretch = 'lin', ggLayer = T) + xlim(523643.3, 544173.3) + ylim(7520437, 7540391) +
  ggRGB(img = cm, r = 1, g = 2, b = 3 , stretch = 'hist', ggLayer = T) + xlim(523643.3, 544173.3) + ylim(7520437, 7540391) +
  layer_spatial(pni_alta, aes(), color= "black", fill = "transparent", size=0.5) +
  layer_spatial(train_kfold_poins, aes(), color= "blue", size=0.5) +
  layer_spatial(ext_val_poins, aes(), color= "green", size=0.5) +
  ggtitle("Itatiaia National Park", subtitle = "Composite of Hyperspectral Subsurface Image (bands 1422; 1001; 430 nm)") + 
  annotation_scale(location = "tl") +
  annotation_north_arrow(location = "br", which_north = "true") +
  theme_light() + coord_sf(datum= 32723) # theme_dark() # theme_bw() # + coord_sf(crs = 4326)

# ggsave(paste0("../output/graphs_plots/",
#               "Subsurface_image_1422_1001_430.png"),
#        plot = last_plot(), width = 16, height = 15, units = "cm", dpi = 300)
ggsave(paste0("../output/graphs_plots/", 
              "Subsurface_image_1422_1001_430.pdf"), 
       plot = last_plot(), width = 16, height = 15, units = "cm", dpi = 300)


################################################################################
### Load and export spectral raster prediction
################################################################################
#' ### Auto setwd
library(rstudioapi); current_path <- getActiveDocumentContext()$path; setwd(dirname(current_path)) ### https://eranraviv.com/r-tips-and-tricks-working-directory/
getwd()
gc(); rm(list=ls())

library(data.table)  
library(stringr)

### load files names and path
files <- list.files(path = "../output/dataset/Spect_subsuf/", 
                    pattern = "goof_rf_k_fold_cv_metrics",
                    all.files= T, full.names= T)
### sort files
files <- str_sort(files, numeric = T)

### load files and name them for each wavelength
Sp_kfcv <- list()
for (i in 1:130) {
  Sp_kfcv[[i]] <- read_csv(files[i], n_max= 8) 
  Sp_kfcv[[i]]$Wavelength <- paste0("W", str_extract_all(files[i], "\\d+")[[1]])
}

### bind as table
Sp_kfcv_df <- rbindlist(Sp_kfcv)

### export as csv
library(readr)
write_csv(Sp_kfcv_df, file = "../output/dataset/Metrics_spectral_raster_prediction.csv")

#   #####   #####   #####   #####   #####   #####   #####   #####   #####   ####
#####   #####   #####   #####   #####   #####   #####   #####   #####   #####   
#   #####   #####   #####   #####   #####   #####   #####   #####   #####   ####
### End script
#####   #####   #####   #####   #####   #####   #####   #####   #####   #####   
#   #####   #####   #####   #####   #####   #####   #####   #####   #####   ####
#####   #####   #####   #####   #####   #####   #####   #####   #####   #####   


