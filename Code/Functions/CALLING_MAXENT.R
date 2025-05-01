#
# _____/------\    ______  /------\\\  \\
#/             \__|      \/        / \\\
#                                      //|/\/\/\/_____---
#     Function to call Maxent       \/  __  /\\\\\
#        Analysis and testing         /\  ___  /\\\\ \\\
#                                  _      \/   |/\\ \  \ \
#        ____           / |     |\  \  \\ \\\\ 
#_______/    \_________/  |____/  \\\\\  \
#
#
do_MAX<-function(spp_points, # Distribution points for the species, either a spatialpoints/dataframe, or route to a dataframe with the coordinates into columns
                  sp.coord=NULL, # if the spp_points is not a SpatialPoints object, provide the coordinates fields e.g x,y for ~x+y, otherwise it assumes that the coordiantes are of the GBIF form
                  n=1, # number of times the model is going to be repeated to average the results (default 1)
                  predictors, # The predictors or environmental variables
                  background_n=10000, # number of background points
                  weights_bk="Random",# [Character]["BwData","Random"] Method to draw the background points: BwData, background weighted data, a density kernell is first calculated and 
                                   # background data is created based on the probabilities of this kernell; Random, points are randomly generated within 
                                   # the limits of the specified range (default) 
                  TrainTest=0.8, # proportion of data used for the testing and training of the models, default is 80% train 20% testing
                  SPP.range=NULL, # Route to the range information of species
                  buffer.dist=1, # Buffer to apply to the range of species, for the delimitation of the background sampling area, default is 1 degree (latlon) or 100km (utm). If no range data is supplied the MCP is calculated and used
                  world_smpl=NULL, # Route to the world_profile to use to trim the buffered areas 
                  crs.r="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", # CRS for the spatial information, default = wgs84
                  results.maxent=paste(getwd(),"spp_Maxent",sep="/"), # folder to store MAXENT results
                  mod.name,
                  bias.lyr=NULL, # Is a bias layer available for the are where the species was sampled? default is to not included this into the model
                  replicates=1, # Number of times we want to run the model 
                  report=TRUE # An rmarkdown report summarising some model parameters
                  ){
   
   # 0. load the needed libraries----
   list.of.packages<-c("doParallel","dplyr","foreach","terra","sf","dismo","rJava","raster","predicts")
   
   new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
   if(length(new.packages)) install.packages(new.packages)
   
   lapply(list.of.packages,require,character.only=TRUE)
   rm(list.of.packages,new.packages)
   
   # b. Creates the Temporal file folder and the folder to allocated the results----
   if(!dir.exists(results.maxent)){dir.create(results.maxent)}
   
   # 1. load the needed functions (testing functions)----
   source("./Code/Functions/BackgroundPOINTS.R") # Prepare background points
   
   # 1. b Some basic information for the configuration of points and other layers----
   res.r<-terra::res(predictors)
   
   # 2. Prepare the presence data----
   # a. Check the format of the point data
   if(!"sf" %in% class(spp_points)){
      if(is.character(spp_points)){
         spp_points<-fread(spp_points) %>% as.data.frame()
      }
      if(TRUE %in% c(class(spp_points) %in% c("data.frame","data.table"))){
        spp_points <- spp_points %>% st_as_sf(coords=sp.coord,crs=crs.r) # takes the default crs in this case wgs 84
         }
   }
   
   # b. Split the data into training and test data----
   # add and option here to split the data using a field specified by the used
   if(is.numeric(TrainTest)){
      r.samp<-sample(1:nrow(spp_points),size=c(nrow(spp_points)*TrainTest) %>% round(digits=0))
   
      train.p<-spp_points[r.samp,] %>% st_geometry()
      test.p<-spp_points[!1:nrow(spp_points) %in% r.samp,] %>% st_geometry()
   }else{
      train.p<-spp_points %>% filter(type==TrainTest) %>% st_geometry()
      test.p<-spp_points %>% filter(type!=TrainTest) %>% st_geometry()
   }
   
   # 3. Create the background information
   # a. Delimit the area to draw the psudo-absences information----
   if(is.null(SPP.range)){
      bk_points<-backgroundPOINTS(presence=spp_points, 
                                  background_n=10000, 
                                  TrainTest=0.8, 
                                  buffer.dist,
                                  range_samp=SPP.range,
                                  cut_area=world_smpl,
                                  weights=weights_bk)
   }else{
      bk_points<-backgroundPOINTS(presence=spp_points, 
                                  background_n=10000, 
                                  TrainTest=0.8, 
                                  range_samp=SPP.range,
                                  cut_area=world_smpl,
                                  buffer.dist,
                                  weights=weights_bk)
   }
   
   # 4. Add some parameters for the model----
   # 4.1.a If the model is run multiple times, change some of the arguments so the models behave differently----
   if(n>1){
   for(w in 1:n){   
      # set.seed(185)
      
      supress_features<-c("nolinear","noquadratic","noproduct","nothreshold")
      supress_features<-sample(supress_features,size=sample(1:3,1),replace=F)
   
   argsMX=c(
      "autofeature=false", # Prevent the creation of automatic features or relatiohips between predictors provided to the model and the response variable
      supress_features, # Fearures to be included into the model creation
      "defaultprevalence=1.00",
      paste0("betamultiplier=",sample(1:10,1)), # random beta multiplier set between 1-10
      "pictures=true",
      "hinge=true",
      "writeplotdata=true", 
      paste0('threads=',ifelse((detectCores(logical=FALSE)-2)>1,detectCores(logical=FALSE)-2,1)), # Numeric. The number of processor threads to use. 
      # Matching this number to the number of cores on your
      # computer speeds up some operations, especially variable jackknifing.
      "responsecurves=false",
      "jackknife=false", #Logical if TRUE measures the importance of each environmental variable by training 
      # with each environmental variable first omitted, then used in isolation
      "askoverwrite=false"
      )
   
   # 4.1.b Run the model ----
   # 4.1.b.1 Extract the variables values and configure the data.frames and vectors to run MAXENT----
      
   x.rast<-as(predictors,"Raster")
   Train.m <- dismo::maxent(x=x.rast, 
                            removeDuplicates=T, 
                            p=train.p %>% as_Spatial(),
                            a=bk_points$Train %>% as_Spatial(),
                            args=argsMX)
   
   # 5. Test the Train model against the test data
   Test.m <- dismo::evaluate(x=x.rast,
                             p=test.p %>% as_Spatial(), 
                             a=bk_points$Test %>% as_Spatial(), 
                             model=Train.m)
   
   # 5.2 Run some predictions on the data
   pred.m<-terra::predict(Train.m,predictors) #%>% crop(vect(SPP.range)) %>% mask(vect(SPP.range))
   pred.m %>% mask(world_smpl %>% vect()) %>% plot()
   
   plot(pred.m) ; plot(st_geometry(bk_points$Samp_area),add=TRUE)      
     
   # 6. Store the results in the local folder
   mx_mod<-list(Train_model=Train.m,
                Test=Test.m,
                prediction=pred.m %>% raster::raster(),
                Area_mod=bk_points$Samp_area,
                MAXENT_args=argsMX)
   
   results.sp<-paste(results.maxent,mod.name,sep="/") ; dir.create(results.sp,recursive = TRUE,showWarnings = FALSE)
   save(mx_mod, file=paste(results.sp,paste0("Maxent_",mod.name,w,".Rdata"),sep="/"))
   terra::writeRaster(pred.m,paste(results.sp,paste0("Maxent_",mod.name,w,".tif"),sep="/"),overwrite=TRUE)
   
   }
   }else{
   # List of model arguments
   argsMX=c(
      "autofeature=true", # Prevent the creation of automatic features or relatiohips between predictors provided to the model and the response variable
      "defaultprevalence=1.00",
      paste0("betamultiplier=",sample(1:10,1)), # random beta multiplier set between 1-10
      "pictures=true",
      "hinge=true",
      "writeplotdata=true", 
      paste0('threads=',ifelse((detectCores(logical=FALSE)-2)>1,detectCores(logical=FALSE)-2,1)), # Numeric. The number of processor threads to use. 
                                                                                                  # Matching this number to the number of cores on your
                                                                                                  # computer speeds up some operations, especially variable jackknifing.
      "responsecurves=false",
      "jackknife=false",        #Logical if TRUE measures the importance of each environmental variable by training 
                                # with each environmental variable first omitted, then used in isolation
      "askoverwrite=false")
   
   # Add the sampling bias layer
   if(!is.null(bias.lyr)){
      # Include the file.name into the argsMX list
      bias_lyr<-raster(bias.lyr)

      # Need to write the layer to a temporal file
      temp_dir<-paste(getwd(),"TEMP_MAXENT",sep="/") ; dir.create(temp_dir,showWarnings = FALSE)
      raster::writeRaster(bias.lyr, paste(temp_dir,"bias_lyr.tif",sep="/"),overwrite=TRUE)
      argsMX<-c(argsMX,"biasfile=bias_lyr")
   }

   # 4.1.b Run the model ----
   # 4.1.b.1 Extract the variables values and configure the data.frames and vectors to run MAXENT----
   x.rast<-as(predictors,"Raster")
   
 Train.m <- dismo::maxent(x=x.rast, 
                           removeDuplicates=T, 
                           p=train.p %>% as_Spatial(),
                           a=bk_points$Train %>% as_Spatial(),
                           args=argsMX)

   # 5. Test the Train model against the test data
         Test.m <- dismo::evaluate(x=rast_vars,
                                   p=test.p %>% as_Spatial(), 
                                   a=bk_points$Test %>% as_Spatial(), 
                                   model=Train.m)
   
   # 5.2 Run some predictions on the data
         pred.m<-terra::predict(Train.m,predictors) %>% crop(vect(SPP.range)) %>% mask(vect(SPP.range))
         plot(pred.m)      
         
   # 6. Store the results in the local folder
         mx_mod<-list(Train_model=Train.m,
                      Test=Test.m,
                      prediction=pred.m %>% raster::raster(),
                      Area_mod=bk_points$Samp_area,
                      MAXENT_args=argsMX)
         
         save(mx_mod, file=paste(results.maxent,paste0("Maxent_",mod.name,".Rdata"),sep="/"))
         terra::writeRaster(pred.m,paste(results.maxent,paste0("Maxent_",mod.name,".tif"),sep="/"),overwrite=TRUE)
   }
} 

# End of The function
