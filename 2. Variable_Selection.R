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
list.of.packages<-c("dplyr","terra","sf","dismo","rJava","doParallel",
                    "foreach","rstudioapi","stringr","prettymapr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0.1 Configure the working environment----
library(rstudioapi)    
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored

# 0.1.a load some functions----
functions<-"./Code/Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 1. Routes to the records and needed spatial information----
crs_prop<-"+proj=utm +zone=28 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # set the project CRS
list.dirs(getwd(),recursive=FALSE) # List of folder in our working environment

# a. Load some base information
base_route<-"./Data/Spatial/Processed_information/Vectorial"
CanPol <- base_route %>% list.files(pattern = "Pro_BaseLayer.shp$",full.names = TRUE) %>% st_read()
CanPol <- CanPol%>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop)

# b. Load the points records----
points_routes<-"./Data/Records/Processed"
list.files(points_routes,pattern = ".shp$")

points<-st_read(paste(points_routes,"Records_ColumbaJunoniaeAnalysis.shp",sep="/")) 
points <- points %>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop)

points$Presence<-rep(1,nrow(points))

# c. load the variables----
var_route<-"./Data/Spatial/Processed_information/Raster"
list.files(var_route)

m_vars_100 <- list.files(var_route,full.names = TRUE)[list.files(var_route) %>% grepl(pattern="100")] %>% rast()
m_vars_100 # lets change some variables names

names(m_vars_100)[names(m_vars_100) %in% c("TRI","PNOA_MDT25_REGCAN95_HU28_1079_LID")]<-c("T.rough_index","Elevation")

m_vars_300 <- list.files(var_route,full.names = TRUE)[list.files(var_route) %>% grepl(pattern="300")] %>% rast()
m_vars_300

names(m_vars_300)[names(m_vars_300) %in% c("TRI","PNOA_MDT25_REGCAN95_HU28_1079_LID")]<-c("T.rough_index","Elevation")

# c.2 Future scenarios----
list.dirs(var_route,recursive=FALSE)
ssp126_01<-list.dirs(var_route,recursive=FALSE)[1] %>% list.files(full.names = TRUE)

ssp_01<-ssp126_01[1] %>% rast()
names(ssp_01) ; plot(ssp_01)

# c.3 Exmaple Background information----
# Topographic variables and climatic variables can be use for all models, land-use can only 
# be use for present distribution models.
# We have more than 600 point on our "clean" presence points, to explore how variables are distributed and how the model migh look at these
# variables we are going to generate 10000 background points much like MAXENT does.
bk_points<-st_sample(CanPol,size=10000,type="random",by_polygon=FALSE)

names(points)
dat_bk<-data.frame(Id_rec=paste("b",1:length(bk_points),sep="_"),
                   type="background",Presence=0)

points0<-st_sf(dat_bk, geometry = bk_points)
points0<-points0 %>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop)

par(mfrow=c(1,1))
CanPol %>% st_geometry() %>% plot()
bk_points %>% st_geometry() %>% plot(add=TRUE,col="tomato",cex=0.1,pch=19) # quite a dense distribution of points

# 2. Variable selection and spatial autocorrelation----
# 2.1 Extract the variables values from the high definition rasters----
all_points<-rbind(points,points0)

values_100m <- terra::extract(x=m_vars_100,y=all_points %>% st_geometry() %>% vect())
values_100m <- cbind(st_coordinates(all_points),values_100m)

values_100m<-cbind(all_points$Presence,values_100m)
names(values_100m)[1]<-"Presence_abs"

values_100m <- values_100m[values_100m %>% complete.cases(),]

# 2.2 Check the distribution of the variables----
exit_route<-"./Results/Variable_selection" ; exit_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)

# List of variables types
all_vars<-names(m_vars_100) ; names(values_100m)

clim_vars <- paste("bio",1:19,sep="_")
land_vars <- c("Bare","Crop","Grass","Shrub","Tree","Urban")
top_vars <- c("slope","aspect","T.rough_index","flowdir","Roughness","Elevation")

var_list <-list(all_vars,clim_vars,land_vars,top_vars) ; names(var_list) <- c("All","Clim","Land","Top")


png(paste(exit_route,"Variable_distribution.png",sep="/"),width=25,height = 65,units="cm",res=600)
par(mfrow=c(16,6))

for (i in 1:length(var_list$All)){
  par(cex=0.5)
  values_100m[,var_list$All[i]] %>%
  hist(xlab=var_list$All[i],col="grey88",border=NA,main="",las=2,axes=F)
  axis(1) ; axis(2,las=2)
  mtext(side=3,adj=0,text=paste0("Raw Hist ",var_list$All[i]),font=2,cex=0.5)
  
  values_100m[,var_list$All[i]] %>% log1p() %>%
    hist(xlab=var_list$All[i],col="cyan2",border=NA,main="",las=2,axes=F)
  axis(1) ; axis(2,las=2)
  mtext(side=3,adj=0,text=paste0("Log Hist ",var_list$All[i]),font=2,cex=0.5)
  
  values_100m[,var_list$All[i]] %>% scale(center = TRUE,scale = TRUE) %>%
    hist(xlab=var_list$All[i],col="tomato2",border=NA,main="",las=2,axes=F)
  axis(1) ; axis(2,las=2)
  mtext(side=3,adj=0,text=paste0("Scale Hist ",var_list$All[i]),font=2,cex=0.5)
  }

dev.off()


# 2.3 Check the VIF and correlation (0.7 threshold) values of the variables----
var_list
var_selection <- list()
var_cor <- list()

m.route <- paste(exit_route,"VIF_matrix",sep="/") ; dir.create(m.route,showWarnings = FALSE,recursive = TRUE)

for(k in 1:length(var_list)){
  # a. RAW Values----
  xvif<-VIF_vars(yp=values_100m,vars=var_list[[k]])
  write.csv(xvif$matrix,paste(m.route,paste0("RAW_VIF_matrix_",names(var_list[k]),".csv"),sep="/"))
  x.index<-apply(xvif$matrix[,-1]<5,2,function(y) (y[!is.na(y)] %>% length())==sum(y[!is.na(y)])) # check whihc colum as all true values
  
  vif.name <-x.index[x.index==TRUE][1] %>% names() ; vif.index <- xvif$matrix[,c("Variables",vif.name)]
  first_selection<- vif.index[!is.na(vif.index[,2]),] ; names(first_selection)[2]<-"raw.dat"
  first_selection %>% write.csv(paste(m.route,paste0("RAW_VIF_selection_",names(var_list[k]),".csv"),sep="/"))
  
  # RAW values correlation
  cor_raw<-cor(values_100m[,colnames(values_100m) %in% var_list[[k]]],method=c("pearson")) ; cor_raw[upper.tri(cor_raw,diag=TRUE)]<-NA
  cor_raw_selected <-cor(values_100m[,colnames(values_100m) %in% first_selection$Variables],method=c("pearson")) ; cor_raw_selected[upper.tri(cor_raw_selected,diag=TRUE)]<-NA
  x.cor1<-list(cor_raw,cor_raw_selected) ; names(x.cor1) <- c("All","Selected")
  
  # b. LOG transformation----
  xvif<-VIF_vars(yp=values_100m %>% log1p() ,vars=var_list[[k]])
  write.csv(xvif$matrix,paste(m.route,paste0("Log1p_VIF_matrix_",names(var_list[k]),".csv"),sep="/"))
  x.index<-apply(xvif$matrix[,-1]<5,2,function(y) (y[!is.na(y)] %>% length())==sum(y[!is.na(y)])) # check whihc colum as all true values
  
  vif.name <-x.index[x.index==TRUE][1] %>% names() ; vif.index <- xvif$matrix[,c("Variables",vif.name)]
  first_selection_log<- vif.index[!is.na(vif.index[,2]),] ; names(first_selection_log)[2]<-"log.transform"
  first_selection_log %>% write.csv(paste(m.route,paste0("Log1p_VIF_selection_",names(var_list[k]),".csv"),sep="/"))
  
  # Log transform correlation
  y <- values_100m[,colnames(values_100m) %in% var_list[[k]]] ; y <- y %>% log1p()
  cor_raw<-cor(y,method=c("pearson")) ; cor_raw[upper.tri(cor_raw,diag=TRUE)]<-NA
  cor_raw_selected <-cor(y[,colnames(y) %in% first_selection_log$Variables],method=c("pearson")) ; cor_raw_selected[upper.tri(cor_raw_selected,diag=TRUE)]<-NA
  x.cor2<-list(cor_raw,cor_raw_selected) ; names(x.cor2) <- c("All","Selected")
  
  # c. NORMALIZED transformation----
  xvif<-VIF_vars(yp=apply(values_100m,2,function(y) scale(y,center = TRUE,scale = TRUE)) %>% as.data.frame(),vars=var_list[[k]])
  write.csv(xvif$matrix,paste(m.route,paste0("Stand_VIF_matrix_",names(var_list[k]),".csv"),sep="/"))
  x.index<-apply(xvif$matrix[,-1]<5,2,function(y) (y[!is.na(y)] %>% length())==sum(y[!is.na(y)])) # check whihc colum as all true values
  
  vif.name <-x.index[x.index==TRUE][1] %>% names() ; vif.index <- xvif$matrix[,c("Variables",vif.name)]
  first_selection_Z<- vif.index[!is.na(vif.index[,2]),] ; names(first_selection_Z)[2]<-"Standardize"
  first_selection_Z %>% write.csv(paste(m.route,paste0("Stand_VIF_selection_",names(var_list[k]),".csv"),sep="/"))
 
  # NORMALIZED transform correlation
  yy<-values_100m[,colnames(values_100m) %in% var_list[[k]]]
  
  y <- apply(yy,2,function(y) scale(y,center = TRUE,scale = TRUE)) %>% as.data.frame()
  cor_raw<-cor(y,method=c("pearson")) ; cor_raw[upper.tri(cor_raw,diag=TRUE)]<-NA
  cor_raw_selected <-cor(y[,colnames(y) %in% first_selection_Z$Variables],method=c("pearson")) ; cor_raw_selected[upper.tri(cor_raw_selected,diag=TRUE)]<-NA
  x.cor3<-list(cor_raw,cor_raw_selected) ; names(x.cor3) <- c("All","Selected")
  
  # Combine the information
  var_selection[[k]]<-full_join(first_selection_Z,first_selection_log,by='Variables') %>% full_join(first_selection,by='Variables')
  names(var_selection)[length(var_selection)] <- names(var_list)[k]
  
  x.cor <- list(x.cor1,x.cor2,x.cor3) ; names(x.cor) <- c("raw","log1","Z_st")
  var_cor[[k]] <-x.cor ; names(var_cor)[k] <- names(var_list)[k]
  }

rm(x.cor,x.cor1,x.cor2,x.cor3,y,cor_raw,yy,first_selection,first_selection_log,
   first_selection_Z,xvif,x.index,k,i)

var_selection %>% saveRDS(paste(exit_route,paste0("VIF_selection",".rds"),sep="/"))
var_cor %>% saveRDS(paste(exit_route,paste0("Cor_variables",".rds"),sep="/"))

# 3. Check the spatial autocorrelation of the variables and records----
# This will help us to stablish if further spatial thinning is needed for the analysis
type_obs <- list(presence=1,absence=0,all=c(1,0))

for(k in type_obs){ # Evaluate the type of records we have (presence,peudo-absence,all)
  dat_x <-values_100m %>% filter(Presence_abs %in% c(k %>% unlist()))
  
  svg(paste(exit_route,paste0("Mantel_correlogram",k %>% paste(collapse = "_"),".svg"),sep="/"),
       width = 18,height = 20,pointsize = 18)
   
  par(mfrow=c(4,3))
  for(i in 1:length(var_list)){
    
    mantel_correlogram(x=dat_x,colxy = c("X","Y"),n=500, env.vars = var_list[[i]],
                       distances = seq(100,1500,by=100))
    
    mtext(side=3,adj=0,text=paste0("Mantel correlogram for ",
           ifelse(1 %in% k,"presence data\n","absence data\n"),"and ",names(var_list)[i]," variables"),cex=0.5)
    
    mantel_correlogram(x=dat_x %>% log1p(),colxy = c("X","Y"),n=500, env.vars = var_list[[i]],
                       distances = seq(100,1500,by=100))
    
    mtext(side=3,adj=0,text=paste0("Mantel correlogram for ",
          ifelse(1 %in% k,"presence data\n","absence data\n"),"and ",names(var_list)[i]," Log1 variables"),cex=0.5)
    
    mantel_correlogram(x=dat_x %>% scale(center=TRUE,scale=TRUE),colxy = c("X","Y"),n=500, env.vars = var_list[[i]],
                       distances = seq(100,1500,by=100))
    
    mtext(side=3,adj=0,text=paste0("Mantel correlogram for ",
          ifelse(1 %in% k,"presence data\n","absence data\n"),"and ",names(var_list)[i]," Z.std variables"),cex=0.5)
  }
  
dev.off()
rm(dat_x)
}  

# 3.1 Observation thinning----
# The previous test shows that the spatial autocorrelation of our presence data is low, however we are going to remove the points that are
# really close to each other (less than a pixel of separation aprox 50m)
#
obs_points <- values_100m %>% filter(Presence_abs==1) %>% st_as_sf(coords = c("X","Y"))
points_sp_cor <- sp_thining(obs_points,dist.t=50,units.d="m",reps=1)
points_sp_cor <- points_sp_cor[[1]]

points_cor<-cbind(points_sp_cor,st_coordinates(points_sp_cor)) %>% as.data.frame()

mantel_correlogram(x=points_cor,colxy = c("X","Y"),n=500, env.vars = var_list[[1]],
                   distances = seq(50,1500,by=100))

# 3.2 Export the final list of species records----
results_route <- "./Data/Data_analysis" ; results_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)
write.csv(values_100m %>% filter(Presence_abs==1),paste(results_route,"var_values.csv",sep="/")) # just the presence data
obs_points %>% st_write(paste("./Data/Data_analysis","Presence_points.shp",sep="/"),append=FALSE)

#
# End of the script
#

unlink(td,recursive = TRUE)