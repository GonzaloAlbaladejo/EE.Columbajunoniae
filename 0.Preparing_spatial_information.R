rm(list=ls())
#
#
# Spatial information pre-processing and correction
#
# 0. load the needed libraries----
list.of.packages<-c("dplyr","terra","sf","dismo","rJava","readxl","doParallel",
                    "foreach","beepr","rstudioapi","stringr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0.1 Configure the working environment----
library(rstudioapi)    
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored

# 1. Routes to the records and spatial information----
list.dirs(getwd(),recursive=TRUE) # List of folder in our working environment

rast_routes<-"./Data/Spatial/Raster_information"
Vect_routes<-"./Data/Spatial/Vector_information"

results_route<-"./Data/Spatial/Processed_information" ; dir.create(results_route,showWarnings = FALSE)

# 1.b CRS for the project
# crs_prop <- "+init=epsg:4079" # Reg Can 95 (a precise UTM projection for the canary islands) 
epsg <- "+init=epsg:32628"
crs_prop <- "+proj=utm +zone=28 +datum=WGS84 +units=m +no_defs"

# 2. Load the basic vectorial information----
list.dirs(Vect_routes,full.names = "TRUE",recursive=FALSE)

Can_polRef<-list.files("./Data/Spatial/Vector_information/Polygon",pattern = ".shp$",full.names = TRUE) %>% sf::st_read()# Administrative Area
Can_polRef %>% sf::st_crs()

# Change the CRS and projection
# First we need to change from psudo-mercator to mercator using the WGS84 projection and then change this to RegCan95
Can_polRef <- Can_polRef %>% st_transform(crs_prop) %>% st_set_crs(crs_prop) 

# Split the object to have individual islands
Can_polRef$Island<-NA
plot(Can_polRef[1,] %>% st_geometry()) # Gran Canaria, Lanzarote, Fuerteventura & La Graciosa (this one needs to be split)
plot(Can_polRef[2,] %>% st_geometry()) ; Can_polRef$Island[2]<- "Tenerife" # Tenerife
plot(Can_polRef[3,] %>% st_geometry()) ; Can_polRef$Island[3]<- "La_Gomera"# La Gomera
plot(Can_polRef[4,] %>% st_geometry()) ; Can_polRef$Island[4]<- "EL_Hierro"# El Hierro
plot(Can_polRef[5,] %>% st_geometry()) ; Can_polRef$Island[5]<- "La_Palma" # La Palma

# Split the MUNTYPOLYGON into individual islands
Can_polRef[1,1] %>% plot()
Can_to_Split<-Can_polRef[1,] %>% st_geometry() %>% st_cast("POLYGON")

Area_pols<-data.frame(Area_pols=Can_to_Split %>% st_area(by_element=TRUE),
                      index=1:length(Can_to_Split))

Area_pols$Area_pols %>% range()
Area_pols[Area_pols$Area_pols==max(Area_pols$Area_pols),"index"]
head(Area_pols[order(Area_pols$Area_pols,decreasing = TRUE),]) # check the first values

# Check the polygons
plot(Can_to_Split[63] %>% st_geometry()) # Fuerteventura   
plot(Can_to_Split[46] %>% st_geometry()) # GranCanaria
plot(Can_to_Split[37] %>% st_geometry()) # Lanzarote
plot(Can_to_Split[152]%>% st_geometry()) # La Graciosa
plot(Can_to_Split[114]%>% st_geometry()) # Montana clara
plot(Can_to_Split[41] %>% st_geometry())  # Alegranza

# Create a data.frame with the same arguments as Can_polRef
names(Can_polRef)
d <- data.frame(id=c(1:6),
                Island=c("Fuerteventura","Gran_Canaria","Lanzarote",
                         "La_graciosa","Montana_clara","Alegranza")
)

pols <- st_sfc(Can_to_Split[[63]],Can_to_Split[[46]],
               Can_to_Split[[37]],Can_to_Split[[152]],
               Can_to_Split[[114]],Can_to_Split[[41]])

Can_2<-st_sf(data.frame(id=c(1:6),
                        Island=c("Fuerteventura","Gran_Canaria","Lanzarote",
                                 "La_graciosa","Montana_clara","Alegranza"),
                        geom=pols),crs=crs_prop)

Can_polRef2<-rbind(Can_polRef,Can_2) %>% st_cast("MULTIPOLYGON")
Can_polRef2<-Can_polRef2[!is.na(Can_polRef2$Island),]

plot(st_geometry(Can_polRef2),col=2:10) # now we can select each individual island to plot!

# Export the corrected spatial layer
dir.create(paste(results_route,"Vectorial",sep="/"),showWarnings = FALSE)

st_write(Can_polRef2, paste(results_route,"Vectorial","Pro_BaseLayer.shp",sep="/"),append = FALSE)
rm(Can_polRef,pols,d,Can_2,Area_pols,Can_to_Split)

Can_polRef<-Can_polRef2

rm(Can_polRef2)

# 2. Land-use cover for the first models:----
# We are going to compile the land-cover information, transform it into plant-functional types and calculate the mean cover for each pixel
# Load the land-use (cover information)----
LU_route<-"D:/Data/Spatial information/Land-use/ESA LandCoverClimateChangeInitiative/Yearly_tiff" # This is for the raw land-Use information
layers_route<-list.files(LU_route,recursive=TRUE,pattern = ".tif$",full.names = TRUE) # The files are sorted by year from 1992 to 2020

# Load the tiff ESA Lu information and stack the layers
ESA_LU<-rast(layers_route) #%>% terra::project(y=crs_prop)
ESA_LU %>% crs(describe=T)

temp_pol<-Can_polRef

temp_pol %>% st_set_crs("EPSG:4326") # This doesn't change the projection of the layer
temp_pol<-temp_pol %>% st_transform(crs="EPSG:4326") # Reproject the layer

ESA_LU<-terra::crop(ESA_LU,y=terra::ext(temp_pol)) # We don't need the full tif! Cut it to meet the size of the Canary Islands
plot(ESA_LU[[1]]) ; rm(temp_pol)

# 2. Prepare the FPT information for the Cross transformation of raster layers----
# Design a function to work with the stack of layers
source("./Code/Functions/ESA_to_PFT.R")

# 2.1 Cross walk table for PFT-ESA transformation----
# to run the ESA_to_PFT we need to load the cross walk table from Li, W., MacBean et al., 2017
Cross_walkESA<-read.csv("D:/Data/Spatial information/Land-use/ESA LandCoverClimateChangeInitiative/Transform_ESA_into_PFT/PFT_CorrelatesTable.csv")

# 2.1.b Transform the pixel values into vegetation cover metrics----
names(Cross_walkESA)
PFT_types<-names(Cross_walkESA[4:ncol(Cross_walkESA)])

# 2.1.c Create a table to contain the values per PFT and pixel----
PFT<-c("Tree_Broadleaf_Evergreen","Tree_Broadleaf_Deciduous","Tree_Needleleaf_Evergreen", 
       "Tree_Needleleaf_Deciduous","Shrub_Broadleaf_Evergreen","Shrub_Broadleaf_Deciduous", 
       "Shrub_Needleleaf_Evergreen","Shrub_Needleleaf_Deciduous","Natural_Grass","Crop",
       "Bare.Soil","Water","Snow.Ice","Urban","No.data")

# 2.2. Transform the layers into PFT----
# A simple loop to process the information
PFT_route<-paste(tempdir(),"Temp_rast",sep="/") ; dir.create(PFT_route,recursive = TRUE,showWarnings = FALSE)
# PFT_route<-"./Data/Spatial/Raster_information/Land-Use information"
years<-1992:2020

for(k in 1:nlyr(ESA_LU)){  
  # Select the raster
  x.rast<-ESA_LU[[k]]
  
  for(i in 1:length(PFT)){
    # PFT to pixels and values
    pixels<-Cross_walkESA[,colnames(Cross_walkESA)%in%c("Pixel_value",PFT[i])] %>% na.omit()
    
    # Transfer the values to the raster and save the individual PFT in a stack
    PFT_pix<-pixels$Pixel_value
    
    T.matrix<-as.matrix(pixels)
    P.na<-Cross_walkESA[!Cross_walkESA$Pixel_value%in%PFT_pix,"Pixel_value"]
    
    T.na<-matrix(c(P.na,151,rep(NA,times=length(P.na)+1)),ncol=2,byrow = FALSE)
    T.matrix<-rbind(T.matrix,T.na)
    
    # Classify the pixel values into PFT cover values
    t1<-Sys.time()  
    P.1<-classify(x.rast,rcl=T.matrix)
    P.1 <- P.1 %>% terra::project(y=crs_prop)
    
    t2<-Sys.time()
    t1-t2
    
    #
    names(P.1)<-paste(years[k],PFT[i],sep="_")
    route_layer<-paste(PFT_route,gsub(pattern="ESACCI-LC-L4-LCCS-Map-300m-P1Y-",
                                       replacement ="",names(x.rast)) %>% 
                         gsub(pattern = "-v2.0.7cds",replacement = "") %>%
                         gsub(pattern = "-v2.1.1",replacement = ""),sep="/") # create the route of the file
    gc()
    
    dir.create(route_layer,showWarnings = FALSE)
    f<-file.path(route_layer,paste0(names(P.1),".tif"))
    terra::writeRaster(P.1,f,overwrite=TRUE)
    
    gc()
    print(i)
  }
  
  print(k)
  rm(x.rast,T.na,T.matrix,P.1,PFT_pix)
}

rm(f,i,k)

# 2.3 Calculate the mean PFT cover by PFT type---- (the routes and some information processing can be further optimize in this stage)
# Export the PFT crosswalk table
Raster_processed_route<-paste(results_route,"Raster",sep="/") ; Raster_processed_route %>% dir.create(showWarnings = FALSE)
write.csv(Cross_walkESA,paste(Raster_processed_route,"PFT_ESA_CrossWalkTable.csv",sep="/"))

# Lets group some of the PFT 
group_PFT <- c("Tree","Shrub","Crop","Bare","Grass","Urban")
index_dirs <- list.dirs(PFT_route,recursive=FALSE) %>% basename()

list_dirs <- list.dirs(PFT_route,recursive=FALSE)[index_dirs %in% years]

for(w in 1:length(list_dirs)){
  xfiles<-list.files(list_dirs[w],pattern = ".tif$",recursive = FALSE,full.names = TRUE)
  
  for(k in 1:length(group_PFT)){
    r<-xfiles[grep(xfiles,pattern = group_PFT[k])] %>% rast() %>% terra::app(fun=sum,na.rm=TRUE)
    names(r)<-group_PFT[k] ; plot(r)
    
    terra::writeRaster(r,filename = paste(PFT_route,paste0(w,"_sum_",group_PFT[k],".tif"),sep="/"),overwrite=TRUE)
  }

  print(paste(paste(rep("1-UP",times=w),collapse="-"),"!",sep=" "))
  # beepr::beep(sound = 2, expr = NULL)
}

dev.off()
rm(w,r,k,xfiles,list_dirs)

# Calculate the mean of the reduced PFT groups!
list_rasters<-list.files(PFT_route,pattern = "_sum_",recursive = TRUE,full.names = TRUE)

# Another simple loop to process to calculate the mean land_cover
temp_dir2<-paste(PFT_route,"Means",sep="/") ; dir.create(temp_dir2)

for(i in 1:length(group_PFT)){
  
  x<-list_rasters[list_rasters %>% grep(pattern = group_PFT[i])] %>% rast()
  
  x[is.na(x)]<-0 ; x <- x %>% mask(Can_polRef %>% st_geometry() %>% vect()) # convert the NA into 0 and mask using the CanayIsland hihg definition Polygon
  
  mean_cover<-function(w){mean(w,na.rm=TRUE)}
  y<-terra::app(x,fun=mean_cover,cores=detectCores(logical=FALSE)-2)
  
  names(y)<-group_PFT[i]
  plot(y,main=names(y))
  
  terra::writeRaster(y,paste(temp_dir2,paste0("Mean_",group_PFT[i],".tif"),sep="/"),overwrite=TRUE)
  rm(x,y)
  }

# Load the Land-use information and export the rast brick object
Lu_cover<-list.files(temp_dir2,pattern = ".tif$",full.names = "TRUE") %>% rast()
terra::writeRaster(Lu_cover,
                   paste(Raster_processed_route,"MeanLandUseCover.tif",sep="/"),
                   overwrite=TRUE)

# Remove the temporal information
unlink(PFT_route,recursive = TRUE) # Remove the temporal files
dev.off() # clean the plotting display

Lu_cover<-list.files(Raster_processed_route,
                     pattern ="MeanLandUseCover.tif",
                     full.names = "TRUE") %>% rast()

plot(Lu_cover) # check the values of the pixels
names(Lu_cover)

# Check the crs and projection of the project
Lu_cover %>% crs()
Lu_cover <- Lu_cover %>% terra::project(y=crs_prop)

# 3. Topographic information----
# We are going to use the 25m EDM as a template to extract all the topographic information we needed to run the analysis. This 
# would need to be upscalled to meet the spatial resolution of the land-use and climatic informatino. 100 m is going to be our finest resolution for
# the climatic and topographic DSM models. Another model including also the land-use information will be included but not used to make future predictions of 
# species distribution since there is no high-resolution scenarios for land-use change and surely no ones that take into considereations islands
#
# 3.1 Elevation digital data----
# https://www.ign.es/web/ign/portal/inicio - information source MTD25
DEM_route<-"D:/Data/Spatial information/DEM/Canarias_25m"
DEM_list<-list.files(DEM_route,pattern = ".asc$",full.names = TRUE) %>% sprc()

DEM_canarias<-DEM_list %>% terra::merge(na.rm=TRUE)
DEM_canarias<-classify(DEM_canarias,cbind(-Inf,0,NA)) ; plot(DEM_canarias)

# Check the crs and projection of the project
crs(DEM_canarias)<-"+proj=utm +zone=28 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # See https://www.ign.es/web/ign/portal/inicio for the crs of the layer
DEM_canarias <- DEM_canarias %>% terra::project(y=crs_prop)

plot(DEM_canarias)

# Save the DEM information
exit_dem<-"./Data/Spatial/Raster_information/DEM_25" ; dir.create(exit_dem,recursive = TRUE,showWarnings = FALSE)
writeRaster(DEM_canarias,paste(exit_dem,"MTD25_Can_Composite.tiff",sep="/"),overwrite=TRUE)

# 3.2 Claculate topographic variables----
# We are going to use the high definition 25m DTM to make the calculations and then upscale the information
# Create a blank reference raster to coerce all the spatial information to it size and resolution
ref_100_rast<-rast(resolution=c(100,100),crs=crs_prop,extent=ext(Lu_cover),vals=1) # units of the CRS in meters!! Important to fix the resolution
ref_300_rast<-rast(resolution=c(300,300),crs=crs_prop,extent=ext(Lu_cover),vals=1)

# Calculate the distance to the coast from each cell in the raster
ref_100<-ref_100_rast
ref_100<-terra::mask(x=ref_100,mask=vect(st_geometry(Can_polRef)),updatevalue=-999) ; plot(ref_100)

# Calculate the distance of each cell to the ocean----
dist_coast_100m<-terra::costDist(mask(x=ref_100,mask=vect(st_geometry(Can_polRef)),updatevalue=-999),
                                 target=-999,scale=1000) %>% mask(mask=vect(Can_polRef))# scale in meters

# Calculate the distance of each cell to urban areas----
ref_urban<-Lu_cover$Urban
ref_urban[!is.na(ref_urban)]<-0 ; ref_urban[is.na(ref_urban)]<-1

dist_urban_300m<-terra::costDist(ref_urban,target=0,scale=1000) %>% mask(mask=vect(Can_polRef))

plot(dist_coast_100m)
plot(dist_urban_300m)

# 5. Terrain features----
# The DEM is a litle bit to precise for what we need, we should decrease its resolution so it follows the same as the land-cover data
DEM_100<-terra::resample(x=DEM_canarias,y=ref_100_rast,method="bilinear")
DEM_300<-terra::resample(x=DEM_canarias,y=ref_300_rast,method="bilinear")# %>% mask(mask=vect(Can_polRef))

par(mfrow=c(3,1))
plot(DEM_canarias,main="25m") ; plot(DEM_100,main="100m") ; plot(DEM_300,main="300m")

# LEST GET SOME TERRAIN FEATURES! (Using the 200m layer)
# High definition distance to coastline

# Slope and aspect
SlpAs_100<-terra::terrain(DEM_100,v=c("slope","aspect","TRI","flowdir"),neighbors=8)
SlpAs_300<-terra::terrain(DEM_300,v=c("slope","aspect","TRI","flowdir"),neighbors=8)

# Terrain Roughness
f <- matrix(1, nrow=3, ncol=3)
rough_100 <- focal(DEM_100, w=f, fun=function(x, ...) max(x) - min(x), na.rm=TRUE)
names(rough_100)<-"Roughness"

rough_300 <- focal(DEM_300, w=f, fun=function(x, ...) max(x) - min(x), na.rm=TRUE)
names(rough_300)<-"Roughness"

# 6. Get the climatic information----
# Patiño, J. & Collart, F., Vanderpoorten, V., Martin-Esquivel, J.L., Naranjo-Cigala, A., Mirolo, 
# S., Karger, D.N. Spatial resolution determines projected plant responses to climate change on 
# topographically complex islands. Diversitity and Distributions.
# doi: 10.1111/ddi.13757
#
Clim_route <- paste(tempdir(),"Climate",sep="/") ; dir.create(Clim_route,recursive = TRUE,showWarnings = FALSE)

unzip("D:/Data/Spatial information/Climatic information/CanaryClim/22060340/1979-2013.zip", 
      exdir = Clim_route)

clim_layers <- list.files(paste(Clim_route,"combined-100m/bio",sep="/"),pattern=".tif$",full.names = TRUE)
clim_vars <- clim_layers %>% rast() ; plot(clim_vars)

# adapt the variables
clim_vars <- clim_vars %>% project(y = crs_prop) # Set the crs to Eckert 4 equal area-projection
clim_vars <- clim_vars %>% resample(y=rough_100)

clim_vars300 <- clim_vars %>% resample(y=rough_300) ; plot(clim_vars300)

to_rm<-c("CHELSA_CanaryIslands_","_1979_2013")

names_clim<-names(clim_vars) %>% gsub(pattern=to_rm[1],replacement="") ; names_clim<-names_clim %>% gsub(pattern=to_rm[2],replacement="") 
names(clim_vars) <- names_clim
names(clim_vars300) <- names_clim

# Export the present information (land-use,climate,)
LCover<-Lu_cover %>% resample(rough_300) ; writeRaster(LCover,filename=paste(Raster_processed_route,"land-use300m.tif",sep="/"),overwrite=TRUE)
LCover100<-Lu_cover %>% resample(rough_100,method="bilinear") ; writeRaster(LCover100,filename=paste(Raster_processed_route,"land-use100m.tif",sep="/"),overwrite=TRUE)

Top100<-c(SlpAs_100,rough_100,DEM_100) ; writeRaster(Top100,filename=paste(Raster_processed_route,"Topography100m.tif",sep="/"),overwrite=TRUE)
Top300<-c(SlpAs_300,rough_300,DEM_300) ; writeRaster(Top300,filename=paste(Raster_processed_route,"Topography300m.tif",sep="/"),overwrite=TRUE)

writeRaster(clim_vars,filename=paste(Raster_processed_route,"Clim100m.tif",sep="/"),overwrite=TRUE)
writeRaster(clim_vars300,filename=paste(Raster_processed_route,"Clim300m.tif",sep="/"),overwrite=TRUE)

# 6.b Get the different climatic scenarios----
Clim_route_2 <- paste(tempdir(),"Scenarios",sep="/") ; dir.create(Clim_route,recursive = TRUE,showWarnings = FALSE)
unzip("D:/Data/Spatial information/Climatic information/CanaryClim/22060340/2071-2100.zip", 
      exdir = Clim_route_2)

Scenario_dirs <- paste(Clim_route_2,"2071-2100",sep="/") %>% list.dirs(recursive = FALSE)

# For each atmospheric circulation model there are 3 socio economic  pathways (ssp) scenarios for 2071
# ssp126,360 and 585. Create a raster stack for each climatic and ssp scenario.
Scenario_dirs

for(i in 1:length(Scenario_dirs)){
  scenario<-basename(Scenario_dirs[i])
  ssp<-list.dirs(Scenario_dirs[i],full.names = TRUE,recursive = FALSE)
  
  for(k in 1:length(ssp)){
    clim_f<- paste(ssp[k],"bio",sep="/") %>% list.files(pattern = ".tif$",full.names = TRUE) %>% rast()
    clim_f100 <- clim_f %>% resample(rough_100)  
    clim_f300 <- clim_f %>% resample(rough_300,method="bilinear")  
    
    to_rm <- c("CHELSA_CanaryIslands_GFDL-ESM4_r1i1p1f1_ssp126_","_2071_2100")
    
    nClim<-names(clim_f100) %>% gsub(pattern=to_rm[1],replacement="") ; nClim <- nClim %>% gsub(pattern=to_rm[2],replacement="")
    names(clim_f100)<-nClim
    names(clim_f300)<-nClim
    
    route_ex <- paste(Raster_processed_route,ssp[k] %>% basename(),sep="/") 
    route_ex %>% dir.create(recursive = TRUE,showWarnings = FALSE) 
    
    writeRaster(clim_f100,paste(route_ex,paste0(scenario,"2071-2100_100m.tif"),sep="/"),overwrite=TRUE)
    writeRaster(clim_f300,paste(route_ex,paste0(scenario,"2071-2100_300m.tif"),sep="/"),overwrite=TRUE)
    
    paste(scenario,ssp[k]%>%basename,sep="/")
    }
print(scenario)
  }

# 8. Dowload the point information from GBIF----
# 8.1 All C. junoniae observations
source("./Code/Functions/GBIF_function.R")
Dowload_gbif(sp_list="Columba junoniae",
                      initial_date=1970,
                      exit_route = "./Data/Records")

points_palomas<-list.files("./Data/Records/",pattern = "Columba junoniae.csv",full.names = TRUE) %>% read.csv()

# Check the coordinates
palomas_limpios<-CoordinateCleaner::clean_coordinates(points_palomas,
                                              lat="decimalLatitude",
                                              lon="decimalLongitude")

palomas_limpios %>% write.csv(paste("./Data/Records","C.junoniaeGBIF_flagged.csv",sep="/"))
palomas_limpios <- palomas_limpios %>% filter(.summary==TRUE)

# Remove the points that fall into the sea
points_palomas<-palomas_limpios %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs="EPSG: 4326") %>% st_transform(crs=crs_prop) # Create the spatial object
points_palomas<-sf::st_intersection(points_palomas,Can_polRef)

names_vars <- c("key","scientificName","geometry","issues","phylum","order","family",                          
                "genus","species","genericName","specificEpithet","taxonRank","taxonomicStatus",                 
                "iucnRedListCategory","continent","stateProvince","year","month","day",                             
                "eventDate","lastInterpreted","license","identifier","facts","relations",                       
                "institutionKey","isInCluster","recordedBy","geodeticDatum","class","countryCode")

sf::st_write(points_palomas[,names_vars],dsn=paste(getwd(),"Data/Records",paste0("C.junoniaeGBIF_points",".csv"),sep="/"),
             layer_options = "GEOMETRY=AS_XY",append=TRUE)

# 8.2 All the records in Gbif for the Canary islands----
wtk_can<-st_as_sfc(st_bbox(Can_polRef %>% st_transform(crs="EPSG: 4326"))) %>% st_as_text()
wtk_can

Dowload_gbif(sp_list=NULL,
             area=wtk_can,
             initial_date=1970,
             exit_route = "./Data/Records")

points_GBIF<-read.csv(list.files("./Data/Records",pattern = "All_records.csv",full.names = TRUE)) #%>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs="EPSG: 4326")
x_flags<-CoordinateCleaner::clean_coordinates(points_GBIF,
                                     lat="decimalLatitude",
                                     lon="decimalLongitude")

x_flags %>% write.csv(paste("./Data/Records","FlaggedGBIF_All.csv",sep="/"))
x_flags <- x_flags %>% filter(.summary==TRUE)

x_flags <- x_flags %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs="EPSG: 4326") %>% st_transform(crs=crs_prop)
x_flags<-sf::st_intersection(x_flags,Can_polRef)

names(x_flags)
x_flags<-x_flags %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs="EPSG: 4326") %>% st_transform(crs=crs_prop) # Create the spatial object

sf::st_write(x_flags[,names_vars],dsn=paste(getwd(),"Data/Records",paste0("GBIF_All_points",".csv"),sep="/"),
             layer_options = "GEOMETRY=AS_XY",append=TRUE,overwrite=TRUE)

# clean the temoporal folder
unlink(tempdir(),recursive = TRUE)

#
# End of the code
#