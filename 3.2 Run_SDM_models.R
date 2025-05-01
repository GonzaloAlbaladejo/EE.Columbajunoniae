rm(list=ls())
td<-tempdir()
dir.create(td,recursive = TRUE)
#
############################################################################################
#                             Final models and results                                  #----
############################################################################################
#
# This script automatically run the whole analysis to produce the different 
# MAXENT species distribution models as well as the output sensitivity metrics
#
#
# 0. load the needed libraries----
list.of.packages<-c("dplyr","terra","sf","dismo","rJava","doParallel","mapsf",
                    "foreach","rstudioapi","stringr","prettymapr","data.table")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0.1 Configure the working environment----
library(rstudioapi)    
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored

fig_results_route <- "./Results/Figures" ; fig_results_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)
results_route <-  "./Results" ; results_route %>% dir.create(recursive = TRUE,showWarnings = FALSE)

# 0.1.a load some functions----
functions<-"./Code/Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 1. Routes to the models and needed information----
crs_prop<-"+proj=utm +zone=28 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs" # set the project CRS
list.dirs(getwd(),recursive=FALSE) # List of folder in our working environment

# a. Load some base information----
base_route<-"./Data/Spatial/Processed_information/Vectorial"
CanPol <- base_route %>% list.files(pattern = "Pro_BaseLayer.shp$",full.names = TRUE) %>% st_read()
CanPol <- CanPol %>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop)

# a.1 Load the Tenoya ravine index layer
Tenoya<- st_read("./Data/Spatial/Vector_information/TenoyaRavine.shp") %>% st_transform(crs=crs_prop)

# b. load the selected models----
load(paste("./Data/Data_analysis","MAXENT_selected_models.RData",sep="/"))

# b.2 In order to map the results and make predictions we need also to load the original variables----
# b. load the variables----
var_route<-"./Data/Spatial/Processed_information/Raster"
list.files(var_route)

m_vars_100 <- list.files(var_route,full.names = TRUE)[list.files(var_route) %>% grepl(pattern="100")] %>% rast()
m_vars_100 # lets change some variables names

names(m_vars_100)[names(m_vars_100) %in% c("TRI","PNOA_MDT25_REGCAN95_HU28_1079_LID")]<- c("T.rough_index","Elevation")

# c. load the data for the analysis----
points_route <- "./Data/Records/Processed"
points_cjunoniae <- points_route %>% list.files(pattern = "LIFe_ColumbaJunoniaeClean.shp",full.names = TRUE) %>% st_read()

# 2. Summarizing final models----
# Create some basic plots to summarize the final selected models
# Get the mean AUC values and extract the AUC at .9 likelihood threshold
mean_auc<-function(y,test="auc",return="mean"){
  x <- y$Test[rownames(y$Test) %in% test,]
  x[x==1]<-NA ; mean_x <- x %>% mean(na.rm=TRUE)
  
  if(return=="mean"){
    return(mean_x)
  }
  if(return=="values"){
    return(x %>% t() %>% as.data.frame())
  }
}

# a. Model performance----
auc<-lapply(y.MAXENT,mean_auc)
values_auc <- lapply(y.MAXENT,function(x) mean_auc(x,return="values")) %>% rbindlist()
values_errrorI <- lapply(y.MAXENT,function(x) mean_auc(x,test="typeI.error",return="values")) %>% rbindlist()
values_errorII <- lapply(y.MAXENT,function(x) mean_auc(x,test="typeII.error",return="values")) %>% rbindlist()
values_PPP <- lapply(y.MAXENT,function(x) mean_auc(x,test="PCC",return="values")) %>% rbindlist()

# a.2 Variable contribution----
  vars_contribution<-lapply(y.MAXENT,function(x) x$Var_contrib$Var_contribution[1,] %>% t() %>% as.data.frame() ) %>% rbindlist(use.names=TRUE)
  vars_Permutation<-lapply(y.MAXENT,function(x) x$Var_contrib$Var_contribution[2,] %>% t() %>% as.data.frame() ) %>% rbindlist(use.names=TRUE)
  
  mean_contri <- apply(vars_contribution,MARGIN = 2,mean)
  sd_contri <- apply(vars_contribution,MARGIN = 2,sd)
  
# a.3 plot the model performance
png(paste(fig_results_route,"Final_models_contrib.png",sep="/"),width=22,height=16,res = 600,units="cm")
  par(bg="grey95",col.axis = "black", col.lab = "black",col="black",col.main="black")
  par(mar=c(6,4,4,1))
  par(mfrow=c(1,2)) # c(bottom, left, top, right)
  
  plot(c(1,1),xlim=c(0.05,1),ylim=c(0,1),type="n",axes=F,
       ylab="AUC",xlab="Probability threshold")
  
  for(k in 1:nrow(values_auc)){
    
    lines(x=names(values_auc)%>%as.numeric(),
          y=values_auc[k,],col="skyblue4" %>% adjustcolor(alpha.f = 0.25),type="l",pch=20)
    
    lines(x=names(values_errrorI)%>%as.numeric(),
          y=values_errrorI[k,],col="tomato" %>% adjustcolor(alpha.f = 0.25),type="l",pch=20) 
    
    lines(x=names(values_errorII)%>%as.numeric(),
          y=values_errorII[k,],col="gold4"%>% adjustcolor(alpha.f = 0.25),type="l",pch=20)
    
    lines(x=names(values_PPP)%>%as.numeric(),
          y=values_PPP[k,]/100,col="green4"%>% adjustcolor(alpha.f = 0.25),type="l",pch=20) 
    
  }
  
  # plot the means
  lines(x=names(values_errorII)%>%as.numeric(),
        y=apply(values_errorII,2,function(x) mean(x,na.rm=TRUE)),col="gold4",type="b",pch=20,lwd=2,cex=1.2) 
  
  lines(x=names(values_errrorI)%>%as.numeric(),
        y=apply(values_errrorI,2,function(x) mean(x,na.rm=TRUE)),col="tomato",type="b",pch=20,lwd=2,cex=1.2) 
  
  lines(x=names(values_auc)%>%as.numeric(),
        y=apply(values_auc,2,function(x) mean(x,na.rm=TRUE)),col="skyblue4",type="b",pch=20,lwd=2,cex=1.2)
  
  lines(x=names(values_PPP)%>%as.numeric(),
        y=apply(values_PPP,2,function(x) mean(x,na.rm=TRUE))/100,col="green4",type="b",pch=20,lwd=2,cex=1.2) 
  
  
  mtext(side=3,adj=-0.2,text="a)",line=2,cex=1)
  
  axis(1) ; axis(2,las=2) ; #abline(h=0.7,col="tomato",lwd=2,lty=3,xpd=F)
  legend("right",legend=c("AUC","PPP","Type I error","Type II error"),pt.cex=1.5,cex=0.8,xpd=TRUE,ncol=1,
         lty=1,pch=20,lwd=2,col=c("skyblue4","green4","tomato","gold4"),bty="n",bg=NA,horiz=FALSE)

# Display the variable contribution
  par(mar=c(6,6,4,1))# c(bottom, left, top, right)
  color_ramp<-colorRampPalette(c("grey99","skyblue4"))
  color_ramp<-color_ramp(length(mean_contri))
  
  kp<-barplot(sort(mean_contri,decreasing = TRUE) %>% log(),las=2,
              xlab="Environmental variables",cex.names=0.75,
              axes=F,ylab="Mean contribution to model\n (log)",col=rev(color_ramp))
  axis(2,las=2)
  mtext(side=3,adj=-0.3,text="b)",line=2,cex=1)
  
dev.off()

# b. Model results and summary----
png(paste(fig_results_route,"Best_present_models.png",sep="/"),
    width=23,height=20,res=600,units="cm")

lt<-layout(matrix(c(1:5,rep(11,5),
                    rep(11,5),
                    rep(11,5),
                    rep(11,5),6:10),
                  ncol=5,nrow =6,
                  byrow=TRUE),respect = TRUE) # this is cool but can be improved

layout.show(lt)
rast_pred<-list()

slope <- terrain(m_vars_100$Elevation, "slope", unit="radians")
aspect <- terrain(m_vars_100$Elevation, "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270) ; hill_gc <- hill %>% crop(CanPol[6,] %>% vect())

par(bg="grey95")
for(i in 1:length(y.MAXENT)){
  # Run the predictions
  x<- predict(y.MAXENT[[i]]$Train_model,m_vars_100) %>% crop(CanPol[6,] %>% vect())

  # Plot the maps
  x %>% plot.r(color.r=c("skyblue","cyan4","#54bb64ff","gold"),range.r = c(0,1),legend.r = FALSE)
  plot(hill_gc,legend.r=F,color.r=c("black",NA,"white"),bg="grey95", col=grey(0:100/100), legend=FALSE,alpha=0.15,add=TRUE)

  mtext(side=3,adj=0,text=paste0(letters[i],")"))

rast_pred[[i]]<-x
}

rast_pred2 <- rast_pred %>% rast()
mean_p <- rast_pred2 %>% app(mean,na.rm=TRUE)
sd_p <- rast_pred2 %>% app(sd,na.rm=TRUE)

plot.r(mean_p,color.r=c("skyblue","cyan4","#54bb64ff","gold"),range.r = c(0,1),legend.r = TRUE)
plot(hill_gc,col=grey(0:100/100), legend=FALSE,alpha=0.15,add=TRUE, bg="grey95")
plot(Tenoya %>% st_geometry(),add=T,col="tomato")

mf_scale(pos="bottomleft",cex=1.5,lwd=2)

mtext(side=3,adj=0,text=paste0(letters[i+1],")"),cex=2)
legend("left",legend="Tenoya\nRavine",lty=1,lwd=2,col="tomato",bty="n")
mtext(side=4,adj=0.5,text="Mean probability\n of occurrence",cex=1,line=-4)

dev.off()

# 3. Reclassify models and the areas marked by the models ----
# (extract the areas for the individual models and combine them into a unified map/polygon)
threshold_r<-function(y, # Object containing the Threshold_kappa of the model
                      y.r # rast with the model predictions
                      ){
  # get the treshold of the model        
  t.r <- y$Threshold_kappa
  
  # Classify the raster using the threshold
  rcl_m<-matrix(c(-Inf,t.r,NA,t.r,Inf,1),ncol=3,nrow=2,byrow = TRUE)
  y.c<-terra::classify(y.r,rcl=rcl_m)
  
  y.size<- cellSize(y.c,mask=TRUE,unit="km")
  area_selected <- y.size %>% values() %>% sum(na.rm=TRUE) # Area in km^2
  
  # Get the important areas for the species
  pol.r <- terra::as.polygons(y.c) %>% sf::st_as_sf()
  pol.r <- sf::st_as_sf(pol.r)
  
  return(list(area_km2=area_selected,select=y.c,pol_area=pol.r))
  }

# Run a loop to process the information----
results_present<-list()
  for(i in 1:length(y.MAXENT)){
    results_present[[i]]<-threshold_r(y=y.MAXENT[[i]],y.r=rast_pred2[[i]])
  }

selected_areas_r <- lapply(results_present,function(x) x$select) %>% rast()
selected_total <- selected_areas_r %>% app(sum,na.rm=TRUE) ; select_pol_present <- selected_total %>% terra::classify(cbind(0,Inf,1)) %>% as.polygons()

# Total area
Total_present <- selected_total %>% cellSize(mask=TRUE,unit="km") %>% values() %>% sum(na.rm=TRUE) 
dat_areas <- data.frame(Scenario="Present",model="Present",area=Total_present)

# Get and export the general polygons for the different models----
polygons_preds <- lapply(results_present, function(x) x$pol_area)
Polygons_models <- do.call("rbind",polygons_preds)

Polygons_models$maxent<- 1:10
Polygons_models %>% st_geometry() %>% plot(col="tomato3"%>%adjustcolor(alpha.f = 0.15),border=NA)
Polygons_models %>% st_union() %>% plot(border="black",add=T)

total_pol <- Polygons_models %>% st_union()

# 3.2 Count the cases within the AMPO and the different probabilty areas----
ColumbaLIFE<-paste("./Data/Records","DatosSEG_Life_Gonzalo.xlsx",sep="/") %>% read_xlsx(sheet="DATOSTOTALES")
ColumbaLIFE <- ColumbaLIFE %>% filter(!X_UTM %in% c("-"," ","",NA)) # Remove the records with no coordinates
ColumbaLIFE$X_UTM <- ColumbaLIFE$X_UTM %>% gsub(pattern=",",replacement=".") ; ColumbaLIFE$Y_UTM <- ColumbaLIFE$Y_UTM %>% gsub(pattern=",",replacement=".")

ColumbaLIFE <- ColumbaLIFE %>% st_as_sf(coords=c("X_UTM","Y_UTM"), crs=crs_prop) %>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop)
values_p <- mean_p %>% extract(ColumbaLIFE %>% vect())

Percentage_overlap <- (table(cut(values_p$mean,10))/length(values_p$ID)) *100
percentage_overlap_AMPO <- sum(ColumbaLIFE %>% st_intersects(total_pol,sparse=F))/nrow(ColumbaLIFE)

Percentage_overlap %>% write.csv("./Results/P_overlap.csv")

# 3.3 Export the polygon data----
"./Results/Polygon_maps" %>% dir.create(showWarnings = F,recursive = T)

st_write(Polygons_models,paste("./Results/Polygon_maps","Mod_polygons.shp",sep="/"))
st_write(total_pol,paste("./Results/Polygon_maps","Full_pol.shp",sep="/"))

# Area returned by the models
dims_crop<-crop(hill_gc,select_pol_present) %>% dim()

png(paste(fig_results_route,"present_models_predictions_areas_highest_probability.png",sep="/"),
    width=12000,height=10000,res=600,units="px")

lt<-layout(matrix(c(1:5,rep(11,5),
                    rep(11,5),
                    rep(11,5),
                    rep(11,5),6:10),
                  ncol=5,nrow =6,
                  byrow=TRUE),respect = TRUE) # this is cool but can be improved

layout.show(lt)
rast_pred<-list()

for(i in 1:nlyr(selected_areas_r)){
  par(mar=c(0,0,0,0))
  plot(select_pol_present,axes=F)
  plot.r(hill_gc,legend.r=F,color.r=c("black","white"),add=TRUE)
  
  plot(select_pol_present,axes=F,add=TRUE,lwd=0.1,border="black")
  
  plot(selected_areas_r[[i]],add=TRUE,legend=F,alpha=0.50,col="tomato" %>% adjustcolor(alpha.f = 1))
  mtext(side=3,adj=0.1,line=-4,text=paste0(letters[i],")"),cex=1.5)
}

plot.r(hill_gc,legend.r=F,color.r=c("black","white"))
plot(select_pol_present,col="tomato",add=TRUE,border=NA,alpha=0.7)

mtext(side=1,adj=0,line=-4,text="Area of maximun\npresence probability")
legend("bottomleft",legend="Aggregated selected area",pch=15,col="tomato",bty="n",bg=NA,
       pt.cex=3,xpd=TRUE)

mtext(side=3,adj=0.1,line=-4,text=paste0(letters[i+1],")"),cex=3)


mf_scale(pos="bottomrigth",cex=2.5,lwd=2)
dev.off()

#' 4. Make the predictions using the different scenarios----
#' We only have future data for the variables related with climate (the others are going to remain constant)
vars_100c <- m_vars_100[[-c(1:19)]] 
#plot(vars_100c)

# Load the variables of the different climatic scenarios----
route_scenarios <- "./Data/Spatial/Processed_information/Raster" #"/ssp126"

# 4.1 ssp126 (5 potential scenarios)----
route_scenarios<-list.dirs(route_scenarios,recursive = FALSE)
layers <- route_scenarios %>% list.files(recursive=FALSE,pattern = ".tif$",full.names = TRUE)
layers<-layers[grep(layers,pattern="_100m")]

predictions_future <- list()
areas_predicted <-list()

for(k in 1:length(layers)){
  # a. Get the base information of the layers    
    name<-dirname(layers[k]) %>% basename() # get the name of the scenario
    name_lyr <- layers[k] %>% basename() %>% gsub(pattern=".tif$",replacement="")
  
  # b. load and project the variables----
      rast_future <- layers[k] %>% rast() 
      
      names(rast_future)<-names(rast_future) %>% gsub(pattern = "CHELSA_CanaryIslands_IPSL-CM6A-LR_r1i1p1f1_",replacement = "") %>% gsub(pattern=paste0(name,"_"),replacement="")
      names(rast_future)<-names(rast_future) %>% gsub(pattern = "CHELSA_CanaryIslands_GFDL-ESM4_r1i1p1f1_",replacement = "") %>% gsub(pattern=paste0(name,"_"),replacement="")
      names(rast_future)<-names(rast_future) %>% gsub(pattern = "CHELSA_CanaryIslands_mpi-esm1-2-lr_r1i1p1f1_",replacement = "") %>% gsub(pattern=paste0(name,"_"),replacement="")
      names(rast_future)<-names(rast_future) %>% gsub(pattern = "CHELSA_CanaryIslands_MRI-ESM2-0_r1i1p1f1_",replacement = "") %>% gsub(pattern=paste0(name,"_"),replacement="")
      names(rast_future)<-names(rast_future) %>% gsub(pattern = "CHELSA_CanaryIslands_UKESM1-0-LL_r1i1p1f2_",replacement = "") %>% gsub(pattern=paste0(name,"_"),replacement="")
      names(rast_future)<-names(rast_future) %>% gsub(pattern = "CHELSA_CanaryIslands_GFDL-ESM4_r1i1p1f1_ssp585_",replacement = "") %>% gsub(pattern=paste0(name,"_"),replacement="")      
      names(rast_future)<-names(rast_future) %>% gsub(pattern = "CHELSA_CanaryIslands_mpi-esm1-2-lr_r1i1p1f1_ssp585_",replacement = "") %>% gsub(pattern=paste0(name,"_"),replacement="")      
      
      
      rast_future <- rast_future %>% terra::project(y=hill)
      future_c<-c(rast_future,vars_100c)

  # c. run the predictions----
      fut_pred<-list()
      
    for(l in 1:length(y.MAXENT)){
        # Run the predictions
        fut_pred[[l]]<- predict(y.MAXENT[[l]]$Train_model,future_c) %>% crop(CanPol[6,] %>% vect())
        }    
  
      x <- fut_pred %>% rast()
      
  # d. plot the predictions
      scenario_out <- paste(fig_results_route,"complete_scenarios",name,sep="/")
      
      dir.create(scenario_out,recursive = TRUE)
      png(paste(scenario_out,paste0(paste(name,name_lyr,sep="_"),".png"),sep="/"),width=23,height=23,res = 600,units="cm")
      
      lt<-layout(matrix(c(1:6,rep(11,5),7,
                          rep(11,5),8,
                          rep(11,5),9,
                          rep(11,5),10),
                        ncol=6,nrow =5,
                        byrow=TRUE)) # this is cool but can be improved
      
      # Plot the maps
      for(ll in 1:nlyr(x)){
        x[[ll]] %>% plot.r(color.r=c("#e7f8ceff","#baac53ff","#d19a7fff","#ce4140ff"),range.r = c(0,1),legend.r = F)
        plot(hill_gc,legend.r=F,color.r=c("black",NA,"white"), col=grey(0:100/100), legend=FALSE,alpha=0.15,add=TRUE)

        mtext(side=3,adj=0,text=paste0(letters[ll],")"))
      }
      
      med_pred <- x %>% app(mean,na.rm=TRUE)
      
      plot.r(med_pred,color.r=c("#e7f8ceff","#baac53ff","#d19a7fff","#ce4140ff"),range.r = c(0,1),legend.r = T)
      plot(hill_gc,col=grey(0:100/100), legend=FALSE,alpha=0.15,add=TRUE)
      #plot(st_geometry(CanPol[6,]),add=TRUE)
      mtext(side=3,adj=0,text=paste0(letters[ll+1],")"))
      mtext(side=4,adj=0.5,text="Mean probability\n of occurrence",cex=1,line=-5)
      mtext(side=2,adj=0,text=paste(name,name_lyr %>% gsub(pattern="-2100_100m",replacement="")))
      
      dev.off()
      
  # e. extract the predicted areas
      results_future<-list()
      
      for(i in 1:nlyr(x)){
        results_future[[i]]<-threshold_r(y=y.MAXENT[[i]],y.r=x[[i]])
      }
      
      selected_areas_future <- lapply(results_future,function(x) x$select) %>% rast()
      selected_total <- selected_areas_future %>% app(sum,na.rm=TRUE) ; select_pol_future <- selected_total %>% terra::classify(cbind(0,Inf,1)) %>% as.polygons()
      
      # Total area
      Total_area_future <- selected_total %>% cellSize(mask=TRUE,unit="km") %>% values() %>% sum(na.rm=TRUE) 
      
      # save the results to export
      predictions_future[[k]] <- med_pred ; names(predictions_future)[k] <- paste(name,name_lyr,sep="_")
      areas_predicted[[k]] <- select_pol_future ; names(predictions_future)[k] <- paste(name,name_lyr,sep="_")
      y <- data.frame(Scenario=name,model=name_lyr,area=Total_area_future)
      dat_areas<-rbind(dat_areas,y)
      }

# 4.2 Check and plot the changes in areas across the different scenarios----
# 4.2.a Changes in the most probable ocupation area (MPO)----
dat_areas
dat_areas %>% group_by(Scenario) %>% summarize(mean(area,na.rm=TRUE)) # areas for each ssp scenario
dat_areas %>% group_by(Scenario) %>% summarize(sd(area,na.rm=TRUE)) 

# PRESENT MPO
area_present<- dat_areas %>% filter(Scenario=="Present") %>% dplyr::select(area) %>% unlist() 

diff_area_ssp126 <- dat_areas %>% filter(Scenario=="ssp126") %>% dplyr::select(area) - area_present
diff_area_ssp370 <- dat_areas %>% filter(Scenario=="ssp370") %>% dplyr::select(area) - area_present
diff_area_ssp585 <- dat_areas %>% filter(Scenario=="ssp585") %>% dplyr::select(area) - area_present

# mean and standard deviation
diff_media<-data.frame(Scenario=c("ssp126","ssp370","ssp585"),Area_Change_Mean=c((mean(diff_area_ssp126$area)*100)/area_present,
                        (mean(diff_area_ssp370$area)*100)/area_present,
                          (mean(diff_area_ssp585$area)*100)/area_present),
           
           sd=c((sd(diff_area_ssp126$area)*100)/area_present,
                  (sd(diff_area_ssp370$area)*100)/area_present,
                      (sd(diff_area_ssp585$area)*100)/area_present))

# Representation of differences
dat_areas$coords_plot<-c(0,rep(c(1,2,3),each=5))

svg(paste(fig_results_route,paste0("Areas_future_scenarios",".svg"),sep="/"),width=22,height=60,pointsize = 32)

  lt<-layout(rbind(matrix(1:15,ncol=3,nrow=5,byrow = FALSE),matrix(rep(16,6),ncol=3,nrow=2,byrow = FALSE)),respect = TRUE)
  layout.show(lt)

  color_r<-colorRampPalette(c("purple4","gold","skyblue"))
  color_r<-color_r(length(unique(dat_areas$model)[-1])) %>% rev
  
# Plot the AMPO areas by SPP and Climatic scenario
  ssp128<-areas_predicted[c(1:5)] ; ssp370<-areas_predicted[c(6:10)] ; ssp585<-areas_predicted[c(11:15)]
  
  bounding_box<-select_pol_present %>% sf::st_bbox()

for(i in 1:length(ssp128)){
  select_pol_present %>% plot(border="black",col="grey88" %>% adjustcolor(alpha.f = 0.3),axes=F,ext=A)
  plot(hill_gc,add=TRUE,alpha=0.75,col=grey(0:100/100),legend=F)
  plot(ssp128[[i]],col=color_r[i] %>% adjustcolor(alpha=0.75),add=TRUE,border=color_r[i],lty=3)
  select_pol_present %>% plot(border="black",col=NA,add=TRUE,axes=F,ext=A)
  mtext(text=c(unique(dat_areas$model)[-1] %>% gsub(pattern="2071-2100_100m",replacement=""))[i],side=2,adj=0.5,line=0,font=2)

  if(i==1){
    mtext(side=3,adj=0.5,text="SSP1-2.8",cex=1.5,font=2)
  }
}

for(i in 1:length(ssp370)){
  select_pol_present %>% plot(border="black",col="grey88" %>% adjustcolor(alpha.f = 0.3),axes=F,ext=A)
  plot(hill_gc,add=TRUE,alpha=0.75,col=grey(0:100/100),legend=F)
  plot(ssp370[[i]],col=color_r[i] %>% adjustcolor(alpha=0.75),add=TRUE,border=color_r[i],lty=3)
  select_pol_present %>% plot(border="black",col=NA,add=TRUE,axes=F,ext=A)
  
  if(i==1){
    mtext(side=3,adj=0.5,text="SSP3-7.0",cex=1.5,font=2)
  }
}

for(i in 1:length(ssp585)){
  select_pol_present %>% plot(border="black",col="grey88" %>% adjustcolor(alpha.f = 0.3),axes=F,ext=A)
  plot(hill_gc,add=TRUE,alpha=0.75,col=grey(0:100/100),legend=F)
  plot(ssp585[[i]],col=color_r[i] %>% adjustcolor(alpha=0.75),add=TRUE,border=color_r[i],lty=3)
  select_pol_present %>% plot(border="black",col=NA,add=TRUE,axes=F,ext=A)
  
  if(i==1){
    mtext(side=3,adj=0.5,text="SSP5-8.5",cex=1.5,font=2)
  }
}

# Evolution of areas graph
par(mar=c(8, 6, 1, 8.5)) # c(bottom, left, top, right)
plot(y=dat_areas$area[-1],x=dat_areas$coords_plot[-1],ylab=expression(paste("AMPO in ",km^2)),
     type="p",axes=F,pch=NA,xlab="Shared Socioeconomic Pathways",cex=2)

axis(1,at=c(1,2,3),labels=c("SSP1-2.6","SSP3-7.0","SSP5-8.5"),lwd=1,lwd.ticks = 2,font=2)
axis(2,labels=seq(0,35,by=5),at=seq(0,35,by=5),las=2)


abline(h=dat_areas$area[1],lwd=2,col="tomato3",lty=2)
text(x=3,y=dat_areas$area[1]+1,"Present occupation",adj=1,cex=1.2,col="tomato3",font=2)
text(x=3.3,y=dat_areas$area[1]+1,col="tomato3",labels=dat_areas$area[1] %>% round(digits=2),xpd=TRUE,cex=1.2)

for(i in 1:length(unique(dat_areas$model)[-1])){
  y.lines <- dat_areas %>% filter(model %in% c(unique(dat_areas$model[-1])[i]))
  lines(y=y.lines$area,x=y.lines$coords_plot,type="b",pch=19,lty=1,col=color_r[i],cex=2.1,lwd=1.5)
  
  # add the mean and standard deviation
  mean_mod<-paste(y.lines$area %>% mean() %>% round(digits=2),"+/-",y.lines$area %>% sd() %>% round(digits=2))
  text(x=3.3,y=y.lines %>% filter(Scenario=="ssp585") %>% dplyr::select("area") %>% unlist(),
       labels=mean_mod,xpd=TRUE,cex=1.2)  
}

legend("top",legend=unique(dat_areas$model)[-1] %>% gsub(pattern="2071-2100_100m",replacement=""),
       pch=19,pt.cex = 1.2,col=color_r,lwd=2,
       bty="n",bg=NA,horiz=FALSE,ncol=3,cex=0.8,title="Atmospheric General Circulation Models",
       title.cex = 1,title.font = 2)

dev.off()

#
# Present and future of the White Tailed Laurel Pigeon!! 
# To retain the models with different data!!!!! For a potential paper!!
# update(best_model[[1]],data = rbind(train,test))
#
#
# End of the script
unlink(td,recursive=TRUE)