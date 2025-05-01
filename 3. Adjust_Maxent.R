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
list.of.packages<-c("dplyr","terra","sf","dismo","rJava","doParallel","mapsf",
                    "foreach","rstudioapi","stringr","prettymapr","data.table")

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

# a. Load some base information----
base_route<-"./Data/Spatial/Processed_information/Vectorial"
CanPol <- base_route %>% list.files(pattern = "Pro_BaseLayer.shp$",full.names = TRUE) %>% st_read()
CanPol <- CanPol%>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop)

# b. load the variables----
var_route<-"./Data/Spatial/Processed_information/Raster"
list.files(var_route)

m_vars_100 <- list.files(var_route,full.names = TRUE)[list.files(var_route) %>% grepl(pattern="100")] %>% rast()
m_vars_100 # lets change some variables names

names(m_vars_100)[names(m_vars_100) %in% c("TRI","PNOA_MDT25_REGCAN95_HU28_1079_LID")]<- c("T.rough_index","Elevation")

# c. load the data for the analysis----
results_route <- "./Data/Data_analysis"
points_cjunoniae <- results_route %>% list.files(pattern = "Presence_points.shp",full.names = TRUE) %>% st_read()

points_cjunoniae<-points_cjunoniae %>% st_set_crs(crs_prop) %>% st_transform(crs=crs_prop)

# d. Extra information----
# Sampling bias layers
bias_r01 <- paste("./Data/Spatial/Processed_information/Raster","Bias_records01.tif",sep="/") %>% rast()
bias_dk  <- paste("./Data/Spatial/Processed_information/Raster","Bias_kernel01.tif",sep="/") %>% rast()

par(mfrow=c(2,1))
bias_r01 <- resample(bias_r01,m_vars_100) ; bias_r01 %>% plot()
bias_dk  <- resample(bias_dk,m_vars_100) ; bias_dk %>% plot()

# e. Selected variables for the analysis (VIF-Cor filtering)----
# VIF values
vars.vif<-"./Results/Variable_selection" %>% list.files(pattern = "VIF_selection.rds",full.names = TRUE) %>% readRDS()
select_vars<-list(All=vars.vif$All[vars.vif$All%>%complete.cases(),1],
                  Clim=vars.vif$Clim[vars.vif$All%>%complete.cases(),1],
                  land=vars.vif$Land[vars.vif$All%>%complete.cases(),1],
                  top=vars.vif$Top[vars.vif$All%>%complete.cases(),1]
)
# COR values
# load data
vars.cor<-"./Results/Variable_selection" %>% list.files(pattern = "Cor_variables.rds",full.names = TRUE) %>% readRDS()

extract_cor<-function(x,limit.c=07){
  # Remove correlations higher than the threshold
  x[x>limit.c]<-NA
  
  # Transform the correlation matrix into a 3 entrance data.frame
  y<-as.data.frame(as.table(x)) ; y<-y[complete.cases(y),]
  
  vars.c <-c(unique(y$Var1),unique(y$Var2)) %>% unique() %>% as.character()
  return(vars.c)
}

# Extract the variables below a 0.7 corr limit
all_vars <- vector()
clim_vars <- vector()
land_vars <- vector()
top_vars <- vector()

for(k in 1:3){
  all_vars<-c(all_vars,vars.cor$All[[k]]$Selected %>% extract_cor())
  clim_vars <-c(clim_vars,vars.cor$Clim[[k]]$Selected %>% extract_cor())
  land_vars <-c(land_vars,vars.cor$Land[[k]]$Selected %>% extract_cor())
  top_vars <-c(top_vars,vars.cor$Top[[k]]$Selected %>% extract_cor())
}

# Combine the information
cor_vars<-list( All = all_vars %>% unique(),
                Clim = clim_vars %>% unique(),
                land = land_vars %>% unique(),
                top = top_vars %>% unique())

# e.2 Get the intial set of variables
select.vars1<-list(All=select_vars$All[select_vars$All %in% cor_vars$All],
                   Clim=select_vars$Clim[select_vars$Clim %in% cor_vars$Clim],
                   land=select_vars$land[select_vars$land %in% cor_vars$land],
                   top=select_vars$top[select_vars$top %in% cor_vars$top])

# clean the trash
rm(all_vars,clim_vars,land_vars,top_vars,k)

# 2. Preliminary Maxent model----
#' Create the presence and background points for the models 
#' we are going to test three different types of psudoabsence/backgorund points):
#' a. using a random sampling of points
#' b. using a weighted sampling bias compensation based on presence records
#' c. using a weighted sampling bias compensation based on all the records found on Gbif for the study area
#' 
#' Split the presence data into train and test (usually 0.8-0.2)
set.seed(69)
index_p<-sample(1:nrow(points_cjunoniae),size=(nrow(points_cjunoniae)*0.8) %>% round(digits = 0))

Train_obs <-  points_cjunoniae[index_p,]
Test_obs <-   points_cjunoniae[-c(index_p),]

prop <- (nrow(Train_obs)/(nrow(Train_obs)+(nrow(Test_obs)))) %>% round(digits=1)

set.seed(1188) # set the random parameters
# a. Random distribution----
bk_random<-backgroundPOINTS(presence=points_cjunoniae, background_n=10000, 
                            TrainTest=prop, buffer.dist=1,
                            range_samp=CanPol, cut_area=CanPol,
                            weights="Random")

# b. Presence weighted distribution----
bk_pres_bias <- backgroundPOINTS(presence=points_cjunoniae, background_n=10000, 
                                 TrainTest=prop, buffer.dist=1,
                                 range_samp=CanPol,cut_area=CanPol,
                                 weights="BwData")

# c. Global sampling weighted distribution----
bk_tot_biasRecords <- backgroundPOINTS(presence=points_cjunoniae, background_n=10000, 
                                       TrainTest=prop, buffer.dist=1,
                                       range_samp=NULL,cut_area=CanPol,
                                       weights=NULL,bias.sampling=bias_r01)

bk_tot_biasKernell <- backgroundPOINTS(presence=points_cjunoniae, background_n=10000, 
                                       TrainTest=prop, buffer.dist=1,
                                       range_samp=NULL,cut_area=CanPol,
                                       weights=NULL,bias.sampling=bias_dk)

# 2.2. Rperesent the different spatial distribution of the sampling methods----
mod_figs_route<-"./Results/Figures/Mod_selection" ; mod_figs_route %>% dir.create(showWarnings = FALSE,recursive = TRUE)
mod_route<-"./Results/Mod_selection"; mod_route %>% dir.create(showWarnings = FALSE,recursive = TRUE)

# plot layout and parameters
png(paste(mod_figs_route,"Background_sampling_comparison.png",sep="/"),
    width = 23,height = 23,units="cm",res=600)

ty <- layout(matrix(c(1,2,3,4),ncol=2,nrow=2,byrow=TRUE))
# col.fun <- colorRampPalette(c("skyblue3","cyan4","#54bb64ff","gold"))
# color.r <- col.fun(4) 

par(mar=c(0,1.8,2.5,0),cex=0.5) # c(bottom, left, top, right)
par(bg="grey99",col.axis = 'black', col.lab = 'black',col='black',col.main='black')

CanPol[6,] %>% st_geometry() %>% plot(col="grey85")
bk_random$Train %>% plot(add=TRUE,col="skyblue3" %>% adjustcolor(alpha.f = 0.6),pch=19,csx=0.5)
plot(points_cjunoniae %>% st_geometry(),pch=19,col=c("tomato","skyblue3"),add=TRUE,title="Records")
mtext(side=3,adj=0,text="a)",cex=1.2)
legend("bottomleft",legend=c("Presence","background"),pch=19,col=c("tomato","skyblue3"),bty="n",cex=2)
mapsf::mf_scale(lwd=2,cex=2)

CanPol[6,] %>% st_geometry() %>% plot(col="grey85")
bk_pres_bias$Train %>% plot(add=TRUE,col="skyblue3" %>% adjustcolor(alpha.f = 0.6),pch=19,csx=0.5)
plot(points_cjunoniae %>% st_geometry(),pch=19,col="tomato",add=TRUE)
mtext(side=3,adj=0,text="b)",cex=1.2)

CanPol[6,] %>% st_geometry() %>% plot(col="grey85")
bk_tot_biasKernell$Train %>% plot(add=TRUE,col="skyblue3" %>% adjustcolor(alpha.f = 0.6),pch=19,csx=0.5)
plot(points_cjunoniae %>% st_geometry(),pch=19,col="tomato",add=TRUE)
mtext(side=3,adj=0,text="c)",cex=1.2)

CanPol[6,] %>% st_geometry() %>% plot(col="grey85")
bk_tot_biasRecords$Train %>% plot(add=TRUE,col="skyblue3" %>% adjustcolor(alpha.f = 0.6),pch=19,csx=0.5)
plot(points_cjunoniae %>% st_geometry(),pch=19,col="tomato",add=TRUE)
mtext(side=3,adj=0,text="d)",cex=1.2)

dev.off()

# 3. Call the initial models and check the contributions of the different layers variable contribution----
background_Train <- list(bk_random$Train,bk_pres_bias$Train,bk_tot_biasKernell$Train,bk_tot_biasRecords$Train)
bk_names <- c("random","pres_bias","tot_pres","tot_kernell")

background_Test <- list(bk_random$Test,bk_pres_bias$Test,bk_tot_biasKernell$Test,bk_tot_biasRecords$Test)

# 3.2 Run the models with the different background information and all the available variables----
#' First we are going to run multiple models using the different types of psudo-absences with all the available variables
#' by doing this we can check if the selection of variables we have is relevant for the modelling of the spatial distribution of
#' C. junoniae, this way we might be able to refine our variable selection and shape the final model to run the predictions 
#
# Test_psudo_abs
Test_psudo_all <- list()

for(k in 1:length(background_Test)){
  ty <- layout(matrix(c(1,2,3,4),ncol=2,nrow=2,byrow=TRUE))
  # All variables
  Test_psudo_all[[k]]<-run_maxent(train.p = Train_obs,
                                  test.p = Test_obs,
                                  bk.train = background_Train[[k]],
                                  bk.test= background_Test[[k]],
                                  predictors = m_vars_100[[c(names(m_vars_100) %in% c(select.vars1$All))]],
                                  clip = CanPol,
                                  dp.t = FALSE,
                                  n.m = 50,
                                  jackknife=TRUE,
                                  results.maxent = "./maxent_models/Test_models",
                                  mod.name = bk_names[k],
                                  return = TRUE,
                                  plot.maxent=TRUE)
  }  

# 3.3 Check variable importance and contribution for the different models----
# a. Check the importance of the variables----
for(k in 1:length(Test_psudo_all)){
  names(Test_psudo_all) <- bk_names  
  vars_contribution<-lapply(Test_psudo_all[[k]],function(x) x$Var_contrib$Var_contribution[1,] %>% t() %>% as.data.frame() ) %>% rbindlist(use.names=TRUE)
  vars_Permutation<-lapply(Test_psudo_all[[k]],function(x) x$Var_contrib$Var_contribution[2,] %>% t() %>% as.data.frame() ) %>% rbindlist(use.names=TRUE)

mean_contri <- apply(vars_contribution,MARGIN = 2,mean)
sd_contri <- apply(vars_contribution,MARGIN = 2,sd)

mean_per <- apply(vars_Permutation,MARGIN = 2,mean)
sd_per <- apply(vars_Permutation,MARGIN = 2,mean)

# Combine
if(k==1){
  mean_sd_var_contri <- rbind(mean_contri,sd_contri,mean_per,sd_per) %>% as.data.frame()  
  mean_sd_var_contri <- mean_sd_var_contri %>% mutate(mod=names(Test_psudo_all)[k],metric=row.names(mean_sd_var_contri),.before=1)
} else {
  x <- rbind(mean_contri,sd_contri,mean_per,sd_per) %>% as.data.frame()  
  x <- x %>% mutate(mod=names(Test_psudo_all)[k],metric=row.names(x),.before=1)
    
  mean_sd_var_contri <- rbind(mean_sd_var_contri,x)
  }
}

mean_sd_var_contri %>% write.csv(paste(mod_route,"Mean_sd_variable_contribution.csv",sep="/"))
mean_sd_var_contri %>% print()

# a.2 Plot the information----
svg(paste(mod_figs_route,paste0("Contri_perm_var_test.svg"),sep="/"),
    width = 24,height = 12,pointsize = 18,bg=NA)

ty <- layout(matrix(c(1,2),ncol=2,nrow=1,byrow=TRUE))
col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
color.r <- col.fun(4) 

  par(mar=c(5,7,2.5,2)) # c(bottom, left, top, right)
  par(bg=NA,col.axis = 'black', col.lab = 'black',col='black',col.main='black')
  
  # Prepare the data for the plots
  mean_contri <- mean_sd_var_contri %>% filter(metric=="mean_contri") ; row.names(mean_contri) <- mean_contri$mod ; mean_contri<-mean_contri[,!colnames(mean_contri) %in% c("mod","metric")] %>% as.matrix()
  sd_contri <- mean_sd_var_contri %>% filter(metric=="sd_contri") ; row.names(sd_contri) <- sd_contri$mod ; sd_contri<-sd_contri[,!colnames(sd_contri) %in% c("mod","metric")] %>% as.matrix()
  
  mean_contri <- mean_contri %>% log1p()
  sd_contri <- sd_contri %>% log1p()
  
  # Mean contribution of the variables to the model
  xp<-barplot(mean_contri[,order(colnames(mean_contri),decreasing = TRUE)],beside=TRUE,horiz=TRUE,las=2,cex.names=1,axes=F,
              col=col.fun(nrow(mean_contri)),border=NA,
              xlim=c(min(mean_contri-(sd_contri*-1),na.rm=TRUE),max(mean_contri+sd_contri,na.rm=TRUE)),
              xlab="log(x+1) Contribucion media al modelo")
  mtext(side=3,adj=0,col='black',text=paste("a)"),line=1,cex=1.2)
  error.bar(y=xp,x=mean_contri[,order(colnames(mean_contri),decreasing = TRUE)],
            sd_contri[,order(colnames(mean_contri),decreasing = TRUE)],col='black',length=0.05,angle.x=90,lwd=1)
  abline(v=0,lty=3,col='black',xpd=FALSE)
  axis(1,col='black'); axis(2,lwd=NA,tick = TRUE,lwd.ticks = 1,labels=NA,at=apply(xp,2,mean))
  
  legend("bottomright",legend=rownames(mean_contri),col=col.fun(nrow(mean_contri)),pch=15,bty="n",bg=NA)
  
  # Permutation training gain
  mean_per <- mean_sd_var_contri %>% filter(metric=="mean_per") ; row.names(mean_per) <- mean_per$mod ; mean_per<-mean_per[,!colnames(mean_per) %in% c("mod","metric")] %>% as.matrix()
  sd_per <- mean_sd_var_contri %>% filter(metric=="sd_per") ; row.names(sd_per) <- sd_per$mod ; sd_per<-sd_per[,!colnames(sd_per) %in% c("mod","metric")] %>% as.matrix()
  
  mean_per <- mean_per %>% log1p()
  sd_per <- sd_per %>% log1p()
  
  xp<-barplot(mean_per[,order(colnames(mean_per),decreasing = TRUE)],beside=TRUE,horiz=TRUE,las=2,cex.names=1,axes=F,
              col=col.fun(nrow(mean_per)),border=NA,
              xlim=c(min(mean_per-(sd_per*-1),na.rm=TRUE),max(mean_per+sd_per,na.rm=TRUE)),
              xlab="log(x+1) Importancia durante la permutacion")
  mtext(side=3,adj=0,col='black',text=paste("b)"),line=1)
  error.bar(y=xp,x=mean_per[,order(colnames(mean_per),decreasing = TRUE)],
            sd_per[,order(colnames(mean_per),decreasing = TRUE)],col='black',length=0.05,angle.x=90,lwd=1)
  abline(v=0,lty=3,col='black',xpd=FALSE)
  axis(1,col='black') ; axis(2,lwd=NA,tick = TRUE,lwd.ticks = 1,labels=NA,at=apply(xp,2,mean))

  legend("bottomright",legend=rownames(mean_contri),col=col.fun(nrow(mean_contri)),pch=15,bty="n",bg=NA)
  
dev.off()
  
# b. Check the Performance of the different models----
# AUC (mean and trend for all the models with different background information)----
mod_AUC <- lapply(Test_psudo_all,function(y)
                  lapply(y,function(x) x$Test[1,] %>% t() %>% as.data.frame()) %>% rbindlist())

mod_parameters <- lapply(Test_psudo_all,function(y) 
                  lapply(y,function(x) x$MAXENT_args[x$MAXENT_args %>% grep(pattern="betamultiplier")] %>%
                           gsub(pattern="\\D", replacement="\\1") %>% as.numeric()) %>% unlist())


# b.2 Display the ACU trends of the models----
svg(paste(mod_figs_route,paste0("AUC_models.svg"),sep="/"),width = 24,height = 12,pointsize = 20,bg=NA)

col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
color.r <- col.fun(4) 

par(bg=NA,col.axis = "black", col.lab = "black",col="black",col.main="black")
par(mar=c(5,4,2.5,1))
par(mfrow=c(1,2)) # c(bottom, left, top, right)

plot(c(1,1),xlim=c(0.05,1),ylim=c(0.4,1),type="n",axes=F,
     ylab="AUC",xlab="Umbral de probabilidad")

# for(k in 1:length(mod_AUC)){
#   mod_performance<-mod_AUC[[k]]
#   mod_performance[mod_performance==1]<-NA
# 
#   beta_mod<-mod_parameters[[k]]
# 
# for(i in 1:nrow(mod_performance)){
#   # Remove the values with 1's (these are included to remove NA's values)
#   lines(x=names(mod_performance)%>%as.numeric(),y=mod_performance[i,],col=color.r[k] %>% adjustcolor(alpha.f = 0.1*beta_mod[i]),lwd=1.2)
#   }
# }

# Means only
for(k in 1:length(mod_AUC)){
  mod_performance<-mod_AUC[[k]]
  mod_performance[mod_performance==1]<-NA

  mod_mean<-mod_performance %>% apply(2,function(x) mean(x,na.rm=TRUE))
  mod_sd<-mod_performance %>% apply(2,function(x) sd(x,na.rm=TRUE))

  lines(x=names(mod_performance)%>%as.numeric(),y=mod_mean,col=color.r[k],lwd=2,type="b")
  lines(x=names(mod_performance)%>%as.numeric(),y=mod_mean+mod_sd,lwd=1,lty=3,col=color.r[k],type="l",pch=20)
  lines(x=names(mod_performance)%>%as.numeric(),y=mod_mean-mod_sd,lwd=1,lty=3,col=color.r[k],type="l",pch=20)
}

axis(1) ; axis(2,las=2) ; abline(h=0.7,col="tomato",lwd=2,lty=3,xpd=F)
legend("topleft",title="Pseudo ausencias",legend=names(mod_AUC),
       pch=15,col=color.r,bty="n",bg=NA)

col_legend<-c()
for(i in 1:10){
  col_legend <- cbind(col_legend,adjustcolor("grey",alpha.f = 0.1*i))
  print(i)
}

legend("bottom",title="Parametro beta",legend=c(1:10),pt.cex=1.5,
       pch=15,col=col_legend,bty="n",bg=NA,horiz=TRUE)
mtext(side=3,adj=0,text="a)",cex=1.2)

# b.2 Create a boxplot with the AUC for the different models----
# Create a dataframe to hold all the information
nrow_auc<-mod_AUC[[1]] %>% nrow()
AUC_data.frame <- lapply(mod_AUC,function(x) data.frame(AUC=x %>% as.data.frame() %>% unlist() %>% unname(),threshold=rep(colnames(x),each=nrow_auc)))

AUC_data.frame$random$beta<- rep(mod_parameters$random,each=nrow(AUC_data.frame$random)/nrow_auc) ; AUC_data.frame$random$mod <- "random"
AUC_data.frame$pres_bias$beta<- rep(mod_parameters$pres_bias,each=nrow(AUC_data.frame$pres_bias)/nrow_auc) ; AUC_data.frame$pres_bias$mod <- "pres_bias"
AUC_data.frame$tot_pres$beta<- rep(mod_parameters$tot_pres,each=nrow(AUC_data.frame$tot_pres)/nrow_auc) ; AUC_data.frame$tot_pres$mod <- "tot_pres"
AUC_data.frame$tot_kernell$beta<- rep(mod_parameters$tot_kernell,each=nrow(AUC_data.frame$tot_kernell)/nrow_auc) ; AUC_data.frame$tot_kernell$mod <- "tot_kernell"

AUC_data.frame <- AUC_data.frame %>% rbindlist()
AUC_data.frame[AUC_data.frame$AUC==1,"AUC"]<-NA

# b.2.1 Create the boxplot----
par(bg=NA,col.axis = "black", col.lab = "black",col="black",col.main="black")
par(mar=c(5,4,2.5,1)) # c(bottom, left, top, right)

xb<-boxplot(AUC~mod+threshold,data=AUC_data.frame,axes=F,xlab="Umbral de probabilidad",col=color.r[c(2,1,3,4)],ylim=c(0.4,1))
axis(2,las=2)  
axis(1,las=2,at=seq(4,76,length.out=19),labels=seq(0.05,0.95,by=0.05),lwd=NA,lwd.ticks = 1,las=1,cex=0.5)                                                                                                                        

legend("topleft",legend=names(mod_AUC),pch=15,col=color.r,bty="n",bg=NA)
mtext(side=3,adj=0,text="b)",cex=1.2)

dev.off()

# Is there any difference between the AUC values of the different models?
sink(paste(mod_route,paste0("AUC_models_TEST.txt"),sep="/"))
  # Anova - mod
  print("One way anova test")
  print(aov(AUC~mod,data=AUC_data.frame) %>% summary())
  
  print("Tukey test")
  print(mod_tukey <- aov(AUC~mod,data=AUC_data.frame) %>% TukeyHSD())
  
  # Anova - Beta
  print("One way anova test")
  print(aov(AUC~as.factor(beta),data=AUC_data.frame %>% filter(mod=="random")) %>% summary())
  
  print("Tukey test")
  print(beta_tukey <- aov(AUC~as.factor(beta),data=AUC_data.frame %>% filter(mod=="random")) %>% TukeyHSD())
  
  print("Beta levels")
  print(beta_tukey$`as.factor(beta)` %>% as.data.frame() %>% filter(`p adj`<0.05))

sink()

# 4. Base on the previous data, select the variables and re-run the tests----
# a.1 Which are the best performing models?----
mean_auc<-function(y,test="auc"){
  x <- y$Test[rownames(y$Test) %in% test,]
  if(test %in% c("typeI.error","typeII.error")){
    x[x==1]<-NA ; mean_x <- x %>% min(na.rm=TRUE)
    return(mean_x)  
    
  }else{
  x[x==1]<-NA ; mean_x <- x %>% mean(na.rm=TRUE)
  return(mean_x)
  }
}

mod_filtering<-list(
            mod_random = data.frame(mod_index=1:length(Test_psudo_all$random),auc=lapply(Test_psudo_all$random,mean_auc) %>% unlist(),
                                    typeI.error=lapply(Test_psudo_all$random,function(x) mean_auc(x,test="typeI.error")) %>% unlist()),
            
            mod_presbias = data.frame(mod_index=1:length(Test_psudo_all$pres_bias),auc=lapply(Test_psudo_all$pres_bias,mean_auc) %>% unlist(),
                                      typeI.error=lapply(Test_psudo_all$pres_bias,function(x) mean_auc(x,test="typeI.error")) %>% unlist()),
            
            mod_gbifall = data.frame(mod_index=1:length(Test_psudo_all$tot_pres),auc=lapply(Test_psudo_all$tot_pres,mean_auc) %>% unlist(),
                                     typeI.error=lapply(Test_psudo_all$tot_pres,function(x) mean_auc(x,test="typeI.error")) %>% unlist()),
            
            mod_gbifkernell = data.frame(mod_index=1:length(Test_psudo_all$tot_kernell),auc=lapply(Test_psudo_all$tot_kernell,mean_auc) %>% unlist(),
                                         typeI.error=lapply(Test_psudo_all$tot_kernell,function(x) mean_auc(x,test="typeI.error")) %>% unlist())
      )

# a.1.2 Remove models that have the identical results----
mod_filtering$mod_random<-mod_filtering$mod_random[!paste(mod_filtering$mod_random$auc,mod_filtering$mod_random$kappa) %>% duplicated(),] # some models are identical, we are interested in variation
mod_filtering$mod_presbias<-mod_filtering$mod_presbias[!paste(mod_filtering$mod_presbias$auc,mod_filtering$mod_presbias$kappa) %>% duplicated(),] 
mod_filtering$mod_gbifall<-mod_filtering$mod_gbifall[!paste(mod_filtering$mod_gbifall$auc,mod_filtering$mod_gbifall$kappa) %>% duplicated(),] 
mod_filtering$mod_gbifkernell<-mod_filtering$mod_gbifkernell[!paste(mod_filtering$mod_gbifkernell$auc,mod_filtering$mod_gbifkernell$kappa) %>% duplicated(),] 

# a.1.3 Filter the models base on their mean AUC and TypeI.error values----
mod_select<-lapply(mod_filtering,function(w) w %>% filter(auc>0.7,typeI.error < 0.1)) # Only the models with random background information pass the filtering
best_models<-Test_psudo_all$random[c(mod_select$mod_random$mod_index)]

# 4.1 Variable contribution----
# a. Check the importance of the variables----
  vars_contribution_best<-lapply(best_models,function(x) x$Var_contrib$Var_contribution[1,] %>% t() %>% as.data.frame() ) %>% rbindlist(use.names=TRUE)
  vars_Permutation_best<-lapply(best_models,function(x) x$Var_contrib$Var_contribution[2,] %>% t() %>% as.data.frame() ) %>% rbindlist(use.names=TRUE)
  
  vars_contribution_best[vars_contribution_best==0]<-NA
  vars_select<-apply(vars_contribution_best,2,function(w) is.na(w)%>%sum())
  vars_x<-names(vars_select)[vars_select<c(nrow(vars_contribution_best)-nrow(vars_contribution_best)*0.9)] # select the variables that are important for 80% of the models
  
# 4.2 Repeat the models to check the performance and variable importance again----
  Test_psudo_best <- list()
  
  for(k in 1:length(background_Test)){
    ty <- layout(matrix(c(1,2,3,4),ncol=2,nrow=2,byrow=TRUE))
    # All variables
    Test_psudo_best[[k]]<-run_maxent(train.p = Train_obs,test.p = Test_obs,
                                    bk.train = background_Train[[k]],bk.test= background_Test[[k]],
                                    predictors = m_vars_100[[c(names(m_vars_100) %in% vars_x)]],
                                    clip = CanPol,dp.t = FALSE,n.m = 50,
                                    jackknife=TRUE,results.maxent = "./maxent_models/Test_models",
                                    mod.name = bk_names[k],return = TRUE,
                                    plot.maxent=TRUE)
    }    
  
# 4.3 Repeat the previous visualization of the models----  
  # a. Check the importance of the variables----
  for(k in 1:length(Test_psudo_best)){
    names(Test_psudo_best) <- bk_names  
    vars_contribution<-lapply(Test_psudo_best[[k]],function(x) x$Var_contrib$Var_contribution[1,] %>% t() %>% as.data.frame() ) %>% rbindlist(use.names=TRUE)
    vars_Permutation<-lapply(Test_psudo_best[[k]],function(x) x$Var_contrib$Var_contribution[2,] %>% t() %>% as.data.frame() ) %>% rbindlist(use.names=TRUE)
    
    mean_contri <- apply(vars_contribution,MARGIN = 2,mean)
    sd_contri <- apply(vars_contribution,MARGIN = 2,sd)
    
    mean_per <- apply(vars_Permutation,MARGIN = 2,mean)
    sd_per <- apply(vars_Permutation,MARGIN = 2,mean)
    
    # Combine
    if(k==1){
      mean_sd_var_contri <- rbind(mean_contri,sd_contri,mean_per,sd_per) %>% as.data.frame()  
      mean_sd_var_contri <- mean_sd_var_contri %>% mutate(mod=names(Test_psudo_best)[k],metric=row.names(mean_sd_var_contri),.before=1)
    } else {
      x <- rbind(mean_contri,sd_contri,mean_per,sd_per) %>% as.data.frame()  
      x <- x %>% mutate(mod=names(Test_psudo_best)[k],metric=row.names(x),.before=1)
      
      mean_sd_var_contri <- rbind(mean_sd_var_contri,x)
    }
  }
  
  mean_sd_var_contri %>% write.csv(paste(mod_route,"Mean_sd_variable_contribution_best.csv",sep="/"))
  mean_sd_var_contri %>% print()
  
  
# a.2 Plot the information----
svg(paste(mod_figs_route,paste0("Contri_perm_var_best.svg"),sep="/"),
    width = 24,height = 12,pointsize = 18,bg=NA)

ty <- layout(matrix(c(1,2),ncol=2,nrow=1,byrow=TRUE))
col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
color.r <- col.fun(4) 

par(mar=c(5,7,2.5,2)) # c(bottom, left, top, right)
par(bg=NA,col.axis = 'black', col.lab = 'black',col='black',col.main='black')

# Prepare the data for the plots
mean_contri <- mean_sd_var_contri %>% filter(metric=="mean_contri") ; row.names(mean_contri) <- mean_contri$mod ; mean_contri<-mean_contri[,!colnames(mean_contri) %in% c("mod","metric")] %>% as.matrix()
sd_contri <- mean_sd_var_contri %>% filter(metric=="sd_contri") ; row.names(sd_contri) <- sd_contri$mod ; sd_contri<-sd_contri[,!colnames(sd_contri) %in% c("mod","metric")] %>% as.matrix()

# mean_contri <- mean_contri %>% log1p()
# sd_contri <- sd_contri %>% log1p()

# Mean contribution of the variables to the model
xp<-barplot(mean_contri[,order(colnames(mean_contri),decreasing = TRUE)],beside=TRUE,horiz=TRUE,las=2,cex.names=1,axes=F,
            col=col.fun(nrow(mean_contri)),border=NA,
            xlim=c(min(mean_contri-(sd_contri*-1),na.rm=TRUE),max(mean_contri+sd_contri,na.rm=TRUE)),
            xlab="log(x+1) Contribucion media al modelo")
mtext(side=3,adj=0,col='black',text=paste("a)"),line=1,cex=1.2)
error.bar(y=xp,x=mean_contri[,order(colnames(mean_contri),decreasing = TRUE)],
          sd_contri[,order(colnames(mean_contri),decreasing = TRUE)],col='black',length=0.05,angle.x=90,lwd=1)
abline(v=0,lty=3,col='black',xpd=FALSE)
axis(1,col='black'); axis(2,lwd=NA,tick = TRUE,lwd.ticks = 1,labels=NA,at=apply(xp,2,mean))

legend("bottomright",legend=rownames(mean_contri),col=col.fun(nrow(mean_contri)),pch=15,bty="n",bg=NA)

# Permutation training gain
mean_per <- mean_sd_var_contri %>% filter(metric=="mean_per") ; row.names(mean_per) <- mean_per$mod ; mean_per<-mean_per[,!colnames(mean_per) %in% c("mod","metric")] %>% as.matrix()
sd_per <- mean_sd_var_contri %>% filter(metric=="sd_per") ; row.names(sd_per) <- sd_per$mod ; sd_per<-sd_per[,!colnames(sd_per) %in% c("mod","metric")] %>% as.matrix()

mean_per <- mean_per %>% log1p()
sd_per <- sd_per %>% log1p()

xp<-barplot(mean_per[,order(colnames(mean_per),decreasing = TRUE)],beside=TRUE,horiz=TRUE,las=2,cex.names=1,axes=F,
            col=col.fun(nrow(mean_per)),border=NA,
            xlim=c(min(mean_per-(sd_per*-1),na.rm=TRUE),max(mean_per+sd_per,na.rm=TRUE)),
            xlab="log(x+1) Importancia durante la permutacion")
mtext(side=3,adj=0,col='black',text=paste("b)"),line=1)
error.bar(y=xp,x=mean_per[,order(colnames(mean_per),decreasing = TRUE)],
          sd_per[,order(colnames(mean_per),decreasing = TRUE)],col='black',length=0.05,angle.x=90,lwd=1)
abline(v=0,lty=3,col='black',xpd=FALSE)
axis(1,col='black') ; axis(2,lwd=NA,tick = TRUE,lwd.ticks = 1,labels=NA,at=apply(xp,2,mean))

legend("bottomright",legend=rownames(mean_contri),col=col.fun(nrow(mean_contri)),pch=15,bty="n",bg=NA)

dev.off()


# b. Check the Performance of the different models----
# AUC (mean and trend for all the models with different background information)----
mod_AUC <- lapply(Test_psudo_best,function(y)
  lapply(y,function(x) x$Test[1,] %>% t() %>% as.data.frame()) %>% rbindlist())

mod_parameters <- lapply(Test_psudo_best,function(y) 
  lapply(y,function(x) x$MAXENT_args[x$MAXENT_args %>% grep(pattern="betamultiplier")] %>%
           gsub(pattern="\\D", replacement="\\1") %>% as.numeric()) %>% unlist())


# b.2 Display the ACU trends of the models----
svg(paste(mod_figs_route,paste0("AUC_models_best.svg"),sep="/"),width = 24,height = 12,pointsize = 20,bg=NA)

col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
color.r <- col.fun(4) 

par(bg=NA,col.axis = "black", col.lab = "black",col="black",col.main="black")
par(mar=c(5,4,2.5,1))
par(mfrow=c(1,2)) # c(bottom, left, top, right)

plot(c(1,1),xlim=c(0.05,1),ylim=c(0.4,1),type="n",axes=F,
     ylab="AUC",xlab="Umbral de probabilidad")

# for(k in 1:length(mod_AUC)){
#   mod_performance<-mod_AUC[[k]]
#   mod_performance[mod_performance==1]<-NA
# 
#   beta_mod<-mod_parameters[[k]]
# 
# for(i in 1:nrow(mod_performance)){
#   # Remove the values with 1's (these are included to remove NA's values)
#   lines(x=names(mod_performance)%>%as.numeric(),y=mod_performance[i,],col=color.r[k] %>% adjustcolor(alpha.f = 0.1*beta_mod[i]),lwd=1.2)
#   }
# }

# Means only
for(k in 1:length(mod_AUC)){
  mod_performance<-mod_AUC[[k]]
  mod_performance[mod_performance==1]<-NA
  
  mod_mean<-mod_performance %>% apply(2,function(x) mean(x,na.rm=TRUE))
  mod_sd<-mod_performance %>% apply(2,function(x) sd(x,na.rm=TRUE))
  
  lines(x=names(mod_performance)%>%as.numeric(),y=mod_mean,col=color.r[k],lwd=2,type="b")
  lines(x=names(mod_performance)%>%as.numeric(),y=mod_mean+mod_sd,lwd=1,lty=3,col=color.r[k],type="l",pch=20)
  lines(x=names(mod_performance)%>%as.numeric(),y=mod_mean-mod_sd,lwd=1,lty=3,col=color.r[k],type="l",pch=20)
}

axis(1) ; axis(2,las=2) ; abline(h=0.7,col="tomato",lwd=2,lty=3,xpd=F)
legend("topleft",title="Pseudo ausencias",legend=names(mod_AUC),
       pch=15,col=color.r,bty="n",bg=NA)

col_legend<-c()
for(i in 1:10){
  col_legend <- cbind(col_legend,adjustcolor("grey",alpha.f = 0.1*i))
  print(i)
}

legend("bottom",title="Parametro beta",legend=c(1:10),pt.cex=1.5,
       pch=15,col=col_legend,bty="n",bg=NA,horiz=TRUE)
mtext(side=3,adj=0,text="a)",cex=1.2)

# b.2 Create a boxplot with the AUC for the different models----
# Create a dataframe to hold all the information
nrow_auc<-mod_AUC[[1]] %>% nrow()
AUC_data.frame <- lapply(mod_AUC,function(x) data.frame(AUC=x %>% as.data.frame() %>% unlist() %>% unname(),threshold=rep(colnames(x),each=nrow_auc)))

AUC_data.frame$random$beta<- rep(mod_parameters$random,each=nrow(AUC_data.frame$random)/nrow_auc) ; AUC_data.frame$random$mod <- "random"
AUC_data.frame$pres_bias$beta<- rep(mod_parameters$pres_bias,each=nrow(AUC_data.frame$pres_bias)/nrow_auc) ; AUC_data.frame$pres_bias$mod <- "pres_bias"
AUC_data.frame$tot_pres$beta<- rep(mod_parameters$tot_pres,each=nrow(AUC_data.frame$tot_pres)/nrow_auc) ; AUC_data.frame$tot_pres$mod <- "tot_pres"
AUC_data.frame$tot_kernell$beta<- rep(mod_parameters$tot_kernell,each=nrow(AUC_data.frame$tot_kernell)/nrow_auc) ; AUC_data.frame$tot_kernell$mod <- "tot_kernell"

AUC_data.frame <- AUC_data.frame %>% rbindlist()
AUC_data.frame[AUC_data.frame$AUC==1,"AUC"]<-NA

# b.2.1 Create the boxplot----
par(bg=NA,col.axis = "black", col.lab = "black",col="black",col.main="black")
par(mar=c(5,4,2.5,1)) # c(bottom, left, top, right)

xb<-boxplot(AUC~mod+threshold,data=AUC_data.frame,axes=F,xlab="Umbral de probabilidad",col=color.r[c(2,1,3,4)],ylim=c(0.4,1))
axis(2,las=2)  
axis(1,las=2,at=seq(4,76,length.out=19),labels=seq(0.05,0.95,by=0.05),lwd=NA,lwd.ticks = 1,las=1,cex=0.5)                                                                                                                        

legend("topleft",legend=names(mod_AUC),pch=15,col=color.r,bty="n",bg=NA)
mtext(side=3,adj=0,text="b)",cex=1.2)
abline(h=0.7,lty=3,col="tomato",lwd=2)
dev.off()

#
# Is there any difference between the AUC values of the different models?
sink(paste(mod_route,paste0("AUC_models_BEST.txt"),sep="/"))
# Anova - mod
print("One way anova test")
print(aov(AUC~mod,data=AUC_data.frame) %>% summary())

print("Tukey test")
print(mod_tukey <- aov(AUC~mod,data=AUC_data.frame) %>% TukeyHSD())

# Anova - Beta
print("One way anova test")
print(aov(AUC~as.factor(beta),data=AUC_data.frame %>% filter(mod=="random")) %>% summary())

print("Tukey test")
print(beta_tukey <- aov(AUC~as.factor(beta),data=AUC_data.frame %>% filter(mod=="random")) %>% TukeyHSD())

print("Beta levels")
print(beta_tukey$`as.factor(beta)` %>% as.data.frame() %>% filter(`p adj`<0.05))

sink()

# a.2 Which are the contribution of variables for those models?
sink(paste(mod_route,"Var_contribution_mean_BEST.txt",sep="/"))
  mean_contri<-mean_sd_var_contri %>% filter(metric=="mean_contri") %>% select_if(is.numeric) %>% apply(2,mean)
  mean_contri[order(mean_contri,decreasing = TRUE)]
sink()

# Checking the selection of variables from the second step and the values of its VIF and correlation, which variables are more likely to improve our
# predictions and model performance?
mean_contri[order(mean_contri,decreasing = TRUE)]
names(mean_contri[order(mean_contri,decreasing = TRUE)])
x.var<-c("bio_14","bio_15","Crop","slope","Grass","bio_03","aspect")## Selected variables

# 4. Run the final models with the selected paramters and variables----
x.var<-c(x.var,"Tree","bio_04") # lets add some other variables that can be relevant for the species

# Check VIF values again
dt<-terra::extract(m_vars_100,background_Train[[1]] %>% vect())
VIF_vars(yp=dt,vars=x.var) # bio_12 has a higher VIF value when combined with the other variables


# 4.2. Run a number of maxent models with different parameters----
x.MAXENT_PresBias <- run_maxent(train.p = Train_obs,test.p = Test_obs,
                                bk.train = background_Train[[1]],bk.test= background_Test[[1]],
                                predictors = m_vars_100[[c(names(m_vars_100) %in% x.var)]],
                                clip = CanPol,dp.t = FALSE,n.m =100,jackknife=TRUE,
                                results.maxent = "./maxent_models/Test_models",
                                mod.name = bk_names[k],return = TRUE,plot.maxent=FALSE)

par(mfrow=c(2,5))
for(i in 11:20){x.MAXENT_PresBias[[i]]$prediction %>% crop(CanPol[6,] %>% vect()) %>% plot(axes=F)}

# b. Select a sample of models with high performance (AUC values), we are going to select 20 models----
mod_filtering <- data.frame(mod_index=1:length(x.MAXENT_PresBias),auc=lapply(x.MAXENT_PresBias,mean_auc) %>% unlist(),
                            Error1=lapply(x.MAXENT_PresBias,function(x) mean_auc(x,test="typeI.error")) %>% unlist())

dup_index <- paste(mod_filtering$auc,mod_filtering$Error1) %>% duplicated() # some models are identical, we are interested in variation
mod_filtering <- mod_filtering[!dup_index,]

mods_f <- mod_filtering %>% filter(auc>0.75,Error1 < 0.05) # less than 5% false positive error
mods_f <- mods_f[order(mods_f$auc,decreasing = TRUE),] ; mods_f <- mods_f[c(1:10),]


id<-c(mods_f$mod_index %>% unlist())
y.MAXENT<-x.MAXENT_PresBias[c(id)]

par(mfrow=c(2,5))
for(i in 1:length(y.MAXENT)){y.MAXENT[[i]][["prediction"]] %>% crop(CanPol[6,] %>% vect()) %>% plot(axes=F)}


# Export the results
save(x.MAXENT_PresBias,file = paste(results_route,paste0("MAXENT_models",".RData"),sep="/"))
save(y.MAXENT,file = paste(results_route,paste0("MAXENT_selected_models",".RData"),sep="/"))


#
# End of the script
#

unlink(td,recursive=TRUE)