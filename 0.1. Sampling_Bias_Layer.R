rm(list=ls())
td<-tempdir()
dir.create(td,recursive = TRUE)
#
############################################################################################
#                     Generates a sampling bias raster layer                               #----
############################################################################################
#
#' This script produce two different sampling-bias layers based on the downloaded Gbif records:
#' The first bias-layer is based no the number of records recorded for each cell
#' The second bias-layer is based on a density 2d kernell based on a sub-sample of the records   
#
# 0. load the needed libraries----
list.of.packages<-c("dplyr","terra","sf","dismo","rJava","doParallel",
                    "foreach","rstudioapi","stringr","prettymapr","MASS")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0.1 Configure the working environment----
library(rstudioapi)    
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored

# 0.2 Basic information----
crs_prop<-"+proj=utm +zone=28 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # set the project CRS
list.dirs(getwd(),recursive=FALSE) # List of folder in our working environment

# a. Load some base information----
base_route<-"./Data/Spatial/Processed_information/Vectorial"
CanPol <- base_route %>% list.files(pattern = "Pro_BaseLayer.shp$",full.names = TRUE) %>% st_read()
CanPol <- CanPol %>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop)

# 1. Produce the bias-layer ----
# a. load the data and transform it into spatial information----
ALL_gbif_records <- "./Data/Records" %>% list.files(pattern = "All_records.csv",full.names = TRUE) %>% read.csv()
dat_gbif <- ALL_gbif_records %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs="EPSG:4326")
dat_gbif <- dat_gbif %>% st_transform(crs=crs_prop)

# b. Rasterize and re-scale points----
dat_terrestrial <- dat_gbif %>% st_intersection(CanPol)
ref.rast<-terra::rasterize(CanPol %>% vect(),rast(extent=CanPol %>% vect() %>% ext(),resolution=c(10000,10000)))

# Rasterize the records counts
count_records<-terra::rasterize(vect(dat_terrestrial),y=ref.rast,fun=length)
count_records %>% plot()

# Scale the raster between 0-1 
r.scale01 <- function(x){
  rasterMinMax <- minmax(x)
  rescaledRaster <- (x - rasterMinMax[1])/(rasterMinMax[2] - rasterMinMax[1])
  return(rescaledRaster)
}

par(mfrow=c(2,1))
records01<-r.scale01(count_records)
records01 %>% mask(CanPol %>% st_geometry() %>% vect()) %>% plot()

ref.rast2 <- terra::rasterize(CanPol %>% vect(),rast(extent=CanPol %>% vect() %>% ext(),resolution=c(100,100)))
records01 <- records01 %>% terra::resample(y=ref.rast2) %>% mask(CanPol %>% vect())

records01 %>% plot()

# c. Calculate a density kernel
# Building a kernel density raster of the presence points used to create the background weighting data:
source("./Code/Functions/d.fromXY.R")
obs_density.kernell<-d.fromXY(x=dat_terrestrial[sample(1:nrow(dat_terrestrial),size=200000),],
                              y.size=dim(ref.rast2)[1:2],mask.p = CanPol)

obs_density.kernell %>% plot()

# 2. Resample and reproject the layers to meet the format of the rest of the project----
proj_rast<-"./Data/Spatial/Processed_information/Raster/Clim100m.tif" %>% rast()

par(mfrow=c(2,1))
records01 <- records01 %>% terra::resample(proj_rast) ; records01 %>% plot()
obs_density.kernell <- obs_density.kernell %>% terra::resample(proj_rast) ; obs_density.kernell %>% plot()

# a. Export the bias-layers---- 
writeRaster(records01,paste("./Data/Spatial/Processed_information/Raster","Bias_records01.tif",sep="/"),overwrite=TRUE)
writeRaster(obs_density.kernell,paste("./Data/Spatial/Processed_information/Raster","Bias_kernel01.tif",sep="/"),overwrite=TRUE)

unlink(td,recursive = TRUE)
#
# End of the script!
#