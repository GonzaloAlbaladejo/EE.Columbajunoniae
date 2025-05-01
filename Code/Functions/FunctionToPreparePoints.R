#
###//\/\/\/\/\/\/\////\/\/\/\/\/\/\//\/\/\/\//\////\/\/\///////////////////\\\\\\##-#
#                 Cleaning the points for the maxent analysis
###///\/\/\/\/\/\/\////\/\/\/\/\/\/\\/\////\/\/\/\/\/\/\/\/\/\\\\\/\/\/\/\/\/\///##-#
#
#
Prepare_points<-function(points_sp, # Spatial points data.frame, spatial points or matrix containing the x y coordinates and the name of the species of the points
                         range_sp=NULL,
                         xy.c=c("decimalLongitude","decimalLatitude"), # Variables describing the points coordinates (check GBIF variable naming to correctly assign thsi)
                         b.width=0, # if range_sp is != NULL, a buffer to apply to the distribution data in order to include or exclude points
                         crs.r="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", # CRS for the spatial information, default = wgs84
                         wrld_pol=NULL # optional we can add a polygon to remove the observations that fall into the ocean or water bodies
                         ){
   # 0. Load the required packages if there are not loaded ----
   list.of.packages<-c("dplyr","rgeos","sp","sf","data.table","CoordinateCleaner","raster")
   
   new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
   if(length(new.packages)) install.packages(new.packages)
   
   lapply(list.of.packages,require,character.only=TRUE)
   rm(list.of.packages,new.packages)
   
   T1<-Sys.time()
   
   # a. Clean the GBIF points using Coordinate cleaner----
   # Does the data came directly from GBIF?
   if(class(points_sp)=="gbif_data"){
         points_sp <-points_sp$data %>% as.data.frame()
      }
      
      if(class(points_sp) %in% c("SpatialPoints","SpatialPointsDataFrame")){
         points_sp<-cbind(coordinates(points_sp),points_sp@data)
      }
   
   # Process the spatial information   
   if(is.data.frame(points_sp)){
         
      y<-CoordinateCleaner::clean_coordinates(points_sp,
                                              lon = xy.c[1],
                                              lat = xy.c[2],
                                              tests = c("capitals","centroids","equal","gbif",
                                                           "institutions","outliers","zeros"))$.summary
      points_sp<-points_sp[y,] %>% as.data.frame()
    
     }
   
   # b.1 Create a buffer around the range data and remove the observations that fall outside----
   if(!is.null(range_sp)){
      
      # configure the range data
      range_sp<-sp::spTransform(range_sp,crs.r)
      r.buff<-rgeos::gBuffer(range_sp,byid=FALSE,width=b.width) %>% rgeos::gUnaryUnion()
      
      if(!is.null(wrld_pol)){
         # Configure the wrld_pol data
         wrld_pol<-sp::spTransform(wrld_pol,crs.r)
         r.buff<-rgeos::gIntersection(wrld_pol,r.buff)# %>% plot()
      }
         # Configure the point data      
         points_sp$ID_p<-1:nrow(points_sp)
      
         coordinates(points_sp)<-paste0("~",paste0(xy.c,collapse = "+")) %>% as.formula()
         raster::crs(points_sp)<-crs.r
         points_sp<-sp::spTransform(points_sp,crs.r)
      
         # Check the overlapping of the information
         points_sp<-raster::intersect(points_sp,r.buff)
         points_sp<-cbind(coordinates(points_sp),points_sp@data)
   }

   return(points_sp)

}

# End of the function