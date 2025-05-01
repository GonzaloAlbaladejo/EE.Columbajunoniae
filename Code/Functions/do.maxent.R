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
do.maxent<-function(spp_points, # Distribution points for the species, either a spatialpoints/dataframe, or route to a dataframe with the coordinates into columns
                    sp.coord=NULL, # if the spp_points is not a SpatialPoints object, provide the coordinates fields e.g x,y for ~x+y, otherwise it assumes that the coordiantes are of the GBIF form
                    n=1, # number of times the model is going to be repeated to average the results (default 1)
                    predictors, # The predictors or environmental variables
                    Test_data, # proportion of data used for the testing and training of the models, default is 80% train 20% testing
                    results.maxent=paste(getwd(),"spp_Maxent",sep="/"), # folder to store MAXENT results
                    mod.name,
                    bias.lyr=NULL # Is a bias layer available for the are where the species was sampled? default is to not included this into the model
                    #report=TRUE # An rmarkdown report summarising some model parameters
){
  
  # 0. load the needed libraries----
  list.of.packages<-c("doParallel","dplyr","foreach","terra","sf","dismo","rJava","raster","predicts")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # a. Creates the Temporal file folder and the folder to allocated the results----
  if(!dir.exists(results.maxent)){dir.create(results.maxent)}
  
  # 1. Prepare the presence data----
  # a. Check the format of the point data
  if(!"sf" %in% class(spp_points)){
    if(is.character(spp_points)){
      spp_points<-fread(spp_points) %>% as.data.frame()
    }
    if(TRUE %in% c(class(spp_points) %in% c("data.frame","data.table"))){
      spp_points <- spp_points %>% st_as_sf(coords=sp.coord,crs=crs.r) # takes the default crs in this case wgs 84
    }
  }
  
  # 2. Add some parameters for the model----
  # 2.1.a If the model is run multiple times, change some of the arguments so the models behave differently----
  if(n>1){
    for(w in 1:n){   
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
      
      Train.m <- dismo::maxent(x=data_records, 
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
      
      # 7. Return the results         
      plot(pred.m) ; plot(st_geometry(bk_points$Samp_area),add=TRUE)
      
      # 7.b Create a model report
      # if(report==TRUE){
      #    
      #    names_mod<-paste(mod.name,w,sep="_")
      #    
      #    date<-Sys.Date()
      #    rmarkdown::render(input = "./Base_code/Scripts/SDM_portion/ShortReport.Rmd",clean=TRUE, 
      #                      output_file = paste(results.maxent,mod.name,sep="/"))
      # }
      
      #return(mx_mod)
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
    if(length(train.p)>250){
      Train.m <- dismo::maxent(x=rast_vars, 
                               removeDuplicates=T, 
                               p=train.p %>% as_Spatial(),
                               a=bk_points$Train %>% as_Spatial(),
                               args=argsMX)
    }else{
      Train.m <- dismo::maxent(x=rast_vars, 
                               removeDuplicates=F, 
                               p=train.p%>% as_Spatial(),
                               a=bk_points$Train%>% as_Spatial(),
                               args=argsMX)
    }
    
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
    
    # 7. Return the results
    plot(pred.m);plot(st_geometry(bk_points$Samp_area),add=TRUE)
    
    # 7.b Create a model report
    if(report==TRUE){
      
      names_mod<-mod.name
      
      date<-Sys.Date()
      rmarkdown::render(input = "./Base_code/Scripts/SDM_portion/ShortReport.Rmd",clean=TRUE, 
                        output_file = paste(results.maxent,mod.name,sep="/"))
    }      
    
    
    #     return(mx_mod)
  }
} 

# End of The function
