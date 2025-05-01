rm(list=ls())
td<-tempdir()
dir.create(td,recursive = TRUE)
#
############################################################################################
#                     Prepare the records for C. junoniae                                  #----
############################################################################################
#
# This script automatically run the whole analysis to produce the different 
# MAXENT species distribution models as well as the output sensitivity metrics
#
#
# 0. load the needed libraries----
list.of.packages<-c("dplyr","terra","sf","dismo","rJava","readxl","doParallel","data.table",
                    "foreach","beepr","rstudioapi","stringr","rgbif","prettymapr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0.1 Configure the working environment----
library(rstudioapi)    
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored

# load some functions
# 0.1 Load the functions to run the analysis----
functions<-"./Code/Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 1. Routes to the records and spatial information----
crs_prop<-"+proj=utm +zone=28 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # set the project CRS
list.dirs(getwd(),recursive=FALSE) # List of folder in our working environment

# a. Load the GBIF points records----
points_routes<-"./Data/Records"
list.files(points_routes,pattern = ".csv$")

ColumbaGbif<- paste(points_routes,"C.junoniaeGBIF_points.csv",sep="/") %>% st_read(options=c("X_POSSIBLE_NAMES=X","Y_POSSIBLE_NAMES=Y"))
ColumbaGbif <- ColumbaGbif %>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop) # Set the CRS

# b. Load the Columba junonia monitoring records from the LIFE+ Rabiche----
ColumbaLIFE<-paste(points_routes,"DatosSEG_Life_Gonzalo.xlsx",sep="/") %>% read_xlsx(sheet="DATOSTOTALES")
ColumbaLIFE <- ColumbaLIFE %>% filter(!X_UTM %in% c("-"," ","",NA)) # Remove the records with no coordinates
ColumbaLIFE$X_UTM <- ColumbaLIFE$X_UTM %>% gsub(pattern=",",replacement=".") ; ColumbaLIFE$Y_UTM <- ColumbaLIFE$Y_UTM %>% gsub(pattern=",",replacement=".")

ColumbaLIFE <- ColumbaLIFE %>% st_as_sf(coords=c("X_UTM","Y_UTM"), crs=crs_prop) %>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop)

# c. Load the Reference Polygons
Vect_routes<-"./Data/Spatial/Processed_information/Vectorial"
list.files(Vect_routes,pattern = ".shp$")

CanPol<-Vect_routes %>% list.files(pattern = "Pro_BaseLayer.shp$",full.names = TRUE) %>% st_read() %>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop)# Load the polygons of the Canary Islands
AtlasCj<-Vect_routes %>% list.files(pattern = "UTM_AtlasColumbaJunoniaeSEOBIRDLIFE.shp$",full.names = TRUE) %>% st_read() %>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop)

# 2. Clean the records ----
# a. GBIF data ----
# Remove the records with known issues
head(ColumbaGbif)
issues_data<-ColumbaGbif$issues %>% unique() %>% str_split(pattern = ",") %>% unlist() %>% unique()

Gbif_codes<-rgbif::gbif_issues() ; head(Gbif_codes)
Problems_GBIF<- Gbif_codes %>% filter(code %in% issues_data)

write.csv(Problems_GBIF,paste("./Results","Gbif_issues.csv",sep="/"))

# Remove the record with incomis; inmafu;and bri
issues_g<-c("incomis","inmafu","bri","cdround")

detect_issue<-function(x,y_list=issues_g){
  
  y <- x %>% str_split(pattern = ",") %>% unlist()
  y.t <- TRUE %in% c(y_list %in% y)
  
  return(TRUE %in% y.t)
  
}

index_y<-lapply(ColumbaGbif$issues,detect_issue) %>% unlist()
summary(index_y)

xtabs(~ColumbaGbif[index_y,"issues"]$issues)

C_Gbif_f<-ColumbaGbif[!index_y,]

# Remove the records that doesn't align with SEObirdlife atlas
C_Gbif_f %>% st_crs() == AtlasCj %>% st_crs() # set the CRS and projections
C_Gbif_f <- C_Gbif_f %>% sf::st_intersection(AtlasCj %>% filter(P_Palomas==1) %>% st_geometry())

# Remove duplicated geometries
C_Gbif_f <- C_Gbif_f %>% dplyr::distinct(geometry, .keep_all = TRUE)

# Plot the information
fig_route<-"./Results/Figures" ; dir.create(fig_route,recursive = TRUE,showWarnings = FALSE)

png(paste(fig_route,"Figura.1.png",sep="/"),width = 21,height = 18,units="cm",res=600)
  
  CanPol %>% filter(Island %in% c("Gran_Canaria","Tenerife","EL_Hierro","La_Gomera","La_Palma")) %>% st_geometry() %>% plot()
  ColumbaGbif %>% st_geometry() %>% plot(col="lightgreen" %>% adjustcolor(alpha.f = 0.5),add=TRUE,pch=19)
  C_Gbif_f %>% st_geometry() %>% plot(col="darkgreen" %>% adjustcolor(alpha.f = 0.5),add=TRUE,pch=19)
  legend("bottomleft",horiz=FALSE,legend=c("Localizaciones iniciales","Localizaciones filtradas"),
         col=c("lightgreen","darkgreen"),
         pch=19,pt.cex = 1.2,bty="n")
  
  prettymapr::addnortharrow(pos="bottomright",scale=1)
  
  par(xpd=TRUE)
  prettymapr::addscalebar(pos="bottomright",padin=c(0,-0.1))
  
  mtext(side=3,line=0,font=2,
        text="Figura 1.",cex=1.5,
        adj=0)
  
  mtext(side=3,line=-3,font=3,
        text="Distribucion de los puntos de presencia georeferenciados para C. junoniae obtenidos\nde la plataforma de datos de biodiversidad GBif",
        adj=0)

dev.off()


# b. Clean the records from the LIFE+ Rabiche----
# we have a lot of records from observations points, these are likely to overwheight the values of those 
# pixels on the model, therefore we are going to remove the records with shared coordinates
C_Life_f <- ColumbaLIFE %>% filter(T.Grupo>0 & !T.Grupo %in% NA) # Remove the points/instances in which pigeons were not detected
C_Life_f <- C_Life_f %>% dplyr::distinct(geometry, .keep_all = TRUE) # Remove the duplicated points

C_Life_f <- C_Life_f %>% sf::st_intersection(CanPol %>% filter(Island=="Gran_Canaria") %>% st_geometry())

# Check the records by monitoring type
xtabs(~C_Life_f$TipoMuestreo) # GPS information account for more than half of the individual monitoring points, but they only account for 4 animals

CanPol %>% filter(Island=="Gran_Canaria") %>% st_geometry() %>% plot()

col_gps<-colorRampPalette(c("skyblue4","purple3","green4"))
col_gps<-col_gps(9) ; col_gps <- col_gps %>% adjustcolor(alpha.f = 0.5)

to_thin<-c(paste0("GPS",1:4),"GPS01","GPS02","RadioTracking","TRANS","FM")

for(k in 1:length(to_thin)){
  x <- C_Life_f %>% filter(TipoMuestreo==to_thin[k])
  x %>% st_geometry() %>% plot(add=TRUE,col=col_gps[k],pch=19)

  print(nrow(x))
  
  }

rm(k)

# The bulk of the remote sensing data account for a low diversity of habitats and its highly concentrated, which may cause spatial autocorrelaiton problems
# We need to trim these points a bit
# Filter the points based on distance - for each GPS and radiotracking dataset
  trim_records <- list()
  
  for(i in 1:length(to_thin)){
  set.seed(200)
  trim_records[[i]]<-sp_thining(x=C_Life_f %>% filter(TipoMuestreo==to_thin[i]),
                        dist.t=200, # The hihg resolution raster has a resolution of 100 meters, lets take two pixels of mininum separation between points
                        units.d="m",
                        reps=1 # This is a first thinning so we are goin to repeat the process just one time
                        )[[1]] # just one dataset created
  names(trim_records)[i] <- to_thin[i]
  
  }

Records_trim<-trim_records %>% rbindlist()# %>% st_sf(geom=geometry)
Records_trim <- Records_trim %>% st_sf(crs=crs_prop)

print("observations trimmed!")

# Remove the trimmed data from the initial dataset and merge the filtered records
C_Life_f <- C_Life_f %>% filter(!TipoMuestreo %in% to_thin)
Life_filter <- rbind(C_Life_f,Records_trim)

# Display the data
# png(paste(fig_route,"Figura.2.png",sep="/"),width = 21,height = 18,units="cm",res=600)

  CanPol %>% filter(Island %in% c("Gran_Canaria","Tenerife","EL_Hierro","La_Gomera","La_Palma")) %>% st_geometry() %>% plot()
  ColumbaLIFE %>% st_geometry() %>% plot(col="#ffd1faff" %>% adjustcolor(alpha.f = 0.5),add=TRUE,pch=19)
  Life_filter %>% st_geometry() %>% plot(col="#9005b5ff" %>% adjustcolor(alpha.f = 0.5),add=TRUE,pch=19)
  legend("bottomleft",horiz=FALSE,legend=c("Localizaciones iniciales","Localizaciones filtradas"),
         col=c("#dd97fcff","#9005b5ff"),
         pch=19,pt.cex = 1.2,bty="n")
  
  prettymapr::addnortharrow(pos="bottomright",scale=1)
  
  par(xpd=TRUE)
  prettymapr::addscalebar(pos="bottomright",padin=c(0,-0.1))
  
  mtext(side=3,line=0,font=2,
        text="Figura 2.",cex=1.5,
        adj=0)
  
  mtext(side=3,line=-3.5,font=3,
        text="Distribucion de los puntos de presencia georeferenciados para C. junoniae en la\nisla de Gran Canaria obtenidos de censos oficilaes, observaciones ocasionales y acciones\nde seguimiento incluidas en el Life+ Rabiche",
        adj=0)

# dev.off()
# 3. Export the data for the analysis----

# a. Combine the two datasets to create a unique spatial object for the analysis----
C_Gbif_f <- C_Gbif_f %>% mutate(Id_rec=c(1:nrow(C_Gbif_f)),type="GBIF" %>% rep(times=nrow(C_Gbif_f)))
Life_filter <- Life_filter %>% mutate(Id_rec=c(1:nrow(Life_filter)),type="LIFE" %>% rep(times=nrow(Life_filter)))

Rec_junoniae <- rbind(C_Gbif_f[,colnames(C_Gbif_f) %in% c("Id_rec","type","geometry")],
                      Life_filter[,colnames(Life_filter) %in% c("Id_rec","type","geometry")])
     

# b. Export the filtered records----
exit_route<-"./Data/Records/Processed" ; exit_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)

write_sf(C_Gbif_f,paste(exit_route,"Gbif_ColumbaJunoniaeClean.shp",sep="/"))
write_sf(Life_filter,paste(exit_route,"LIFe_ColumbaJunoniaeClean.shp",sep="/")) 
write_sf(Rec_junoniae,paste(exit_route,"Records_ColumbaJunoniaeAnalysis.shp",sep="/")) 

write_sf(ColumbaGbif,paste(exit_route,"raw_gbif_ColumbaJunoniae.shp",sep="/"))
write_sf(ColumbaLIFE,paste(exit_route,"raw_Life_ColumbaJunoniae.shp",sep="/"))

unlink(td,recursive = TRUE)
#
# End of the script
#