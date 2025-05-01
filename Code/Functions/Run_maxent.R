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
run_maxent<-function(train.p, # Distribution points for the species, either a spatialpoints/dataframe, or route to a dataframe with the coordinates into columns
                 test.p,
                 bk.train,
                 bk.test,
                 dp.t=TRUE,
                 n.m=1, # number of times the model is going to be repeated to average the results (default 1)
                 predictors, # The predictors or environmental variables
                 clip=NULL,
                 results.maxent=paste(getwd(),"spp_Maxent",sep="/"), # folder to store MAXENT results
                 mod.name,
                 jackknife=NULL,        #Logical if TRUE measures the importance of each environmental variable by training 
                 plot.maxent=TRUE,
                 return=FALSE # should models be returned to the main environment?
){
  
  # 0. load the needed libraries----
  options("rgdal_show_exportToProj4_warnings"="none")
  list.of.packages<-c("dplyr","rfUtilities","terra","sf","dismo","rJava","raster")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # 4. Add some parameters for the model----
  # 4.1.a If the model is run multiple times, change some of the arguments so the models behave differently----
  if(n.m>1){
    models <- list()
    
    for(w in 1:n.m){   
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
      
      if(is.null(jackknife)){
        argsMX <- c(argsMX,"jackknife=false")
      }else{
        argsMX <- c(argsMX,"jackknife=true")
      }
      
      # 4.1.b Run the model ----
      # 4.1.b.1 Extract the variables values and configure the data.frames and vectors to run MAXENT----
      x.rast<-as(predictors,"Raster")
      Train.m <- dismo::maxent(x=x.rast, 
                               removeDuplicates=dp.t, 
                               p=train.p %>% as_Spatial(),
                               a=bk.train %>% as_Spatial(),
                               args=argsMX)
      
      # # 5. Test the Train model against the test data
      Test.m <- dismo::evaluate(x=x.rast,
                                p=test.p %>% as_Spatial(),
                                a=bk.test %>% as_Spatial(),
                                model=Train.m,type="response")
      
      kappa_threshold<-threshold(Test.m)$kappa # NEED TO CHECK why the values of this test and the threshold thing are different!!
       
      # 5. Alternative testing
      dat.new<-rbind(terra::extract(x=predictors,y=vect(test.p)),
                     terra::extract(x=predictors,y=vect(bk.test))
      )
      
      index.values<-dat.new %>% complete.cases()
      y.new<-c(rep(1,nrow(test.p)),rep(0,nrow(bk.test)))[index.values]
      
      # 5.1 Calculate the accuracy and specificity of our model at different probability thresholds
      threshold_seq<-seq(0.05,0.95,by=0.05)
      # accuracy_m<-c("kappa","PCC","auc","typeI.error","typeII.error")
      accuracy_m<-c("PCC","auc","typeI.error","typeII.error")
      
      y.acc<-data.frame(test=accuracy_m)
      
      for(l in 1:length(threshold_seq)){
        
          A<-accuracy(ifelse(predict(Train.m,dat.new,type="response") >= threshold_seq[l], 1, 0),y.new)[accuracy_m]
          b <- A %>% unlist() %>% as.data.frame() 
          b$test<- rownames(b)
           
          y.acc <- merge(y.acc,b,by="test",all=TRUE) ; colnames(y.acc)[ncol(y.acc)] <- threshold_seq[l]
            
           rm(A,b)
          }
      
      y.acc[is.na(y.acc)]<-1 ; rownames(y.acc)<-y.acc$test
      y.acc<-y.acc[,-1] %>% as.matrix()
      
      
      # 5.2 Plot the accuracy metrics for the different tresholds
  if(plot.maxent==TRUE){
      par(mar=c(5,5,3.5,3)) # c(bottom, left, top, right)
      par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
      
      col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
      color.r <- col.fun(nrow(y.acc)) 
      
      plot(y=1,x=1,xlim=c(0,1),ylim=c(0,1.1),
           axes=F,
           ylab="Parameter performance",
           xlab="Probability threshold",
           type="n")
      
      mtext(text = "a) Model Accuracy and performance",adj=0,line=1)
      
      for(ll in 1:nrow(y.acc)){
        if(row.names(y.acc)[ll] %in% "PCC"){
          lines(y=y.acc[ll,]/100,x=colnames(y.acc) %>% as.numeric(),
                col=color.r[ll],type="b",pch=19)  
          
        }else{
          lines(y=y.acc[ll,],x=colnames(y.acc) %>% as.numeric(),
                col=color.r[ll],type="b",pch=19)
        }
      }
      
      abline(v=kappa_threshold,lty=3,col="tomato3",xpd=FALSE)
      
      axis(1,seq(0,1,by=0.05),labels =seq(0,1,by=0.05) %>% round(digits = 1) ,gap.axis = 2,col="white")
      axis(2,seq(0,1,by=0.05),labels =seq(0,1,by=0.05) %>% round(digits = 1) ,gap.axis = 2,las=2,col="white")
      
      legend("topleft",legend=rownames(y.acc),lty=1,pch=19,
             col=color.r,bty="n",horiz=FALSE,ncol = 5,x.intersp = 0.5,cex=0.8)
    }
      # 5.2 Run some predictions on the data
      pred.m<-terra::predict(Train.m,predictors)
      
      if(!is.null(clip)){
        pred.m %>% mask(clip %>% vect())
      }
      
      # 5.3 Check the variable importance and variable contribution
      if(is.null(jackknife)){
        jackniffe_var_contrybution<-jack.maxent(Train.m@results,jack=FALSE,plot.j=plot.maxent)
      }else{
        jackniffe_var_contrybution<-jack.maxent(Train.m@results,jack=TRUE,plot.j=plot.maxent)
      }
      
      # 6. Store the results in the local folder
      mx_mod<-list(Train_model=Train.m,
                   Test=y.acc,
                   Test_evaluate=Test.m,
                   Threshold_kappa=kappa_threshold,
                   Var_contrib=jackniffe_var_contrybution,
                   prediction=pred.m,
                   MAXENT_args=argsMX)
      
      models[[w]]<-mx_mod
      
      # 7. Plot the model
  if(plot.maxent==TRUE){  
        par(mar=c(0,0,0,0)) # c(bottom, left, top, right)
        plot(pred.m,col=col.fun(550),axes=F)
        mtext("b) Model predictions",adj=0,line=1)
      }
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
      "askoverwrite=false")
    
    if(is.null(jackknife)){
      argsMX <- c(argsMX,"jackknife=false")
    }else{
      argsMX <- c(argsMX,"jackknife=true")
    }
    
    # 4.1.b Run the model ----
    # 4.1.b.1 Extract the variables values and configure the data.frames and vectors to run MAXENT----
    x.rast<-as(predictors,"Raster")
    
    Train.m <- dismo::maxent(x=x.rast, 
                             removeDuplicates=dp.t, 
                             p=train.p %>% as_Spatial(),
                             a=bk.train %>% as_Spatial(),
                             args=argsMX)
    
    # # 5. Test the Train model against the test data
    Test.m <- dismo::evaluate(x=x.rast,
                              p=test.p %>% as_Spatial(),
                              a=bk.test %>% as_Spatial(),
                              model=Train.m,type="response")

    kappa_threshold<-threshold(Test.m)$kappa
    
    # 5. Alternative testing
    dat.new<-rbind(terra::extract(x=predictors,y=vect(test.p)),
                   terra::extract(x=predictors,y=vect(bk.test))
                   )
    
    index.values<-dat.new %>% complete.cases()
    y.new<-c(rep(1,nrow(test.p)),rep(0,nrow(bk.test)))[index.values]
    
    # 5.1 Calculate the accuracy and specificity of our model at different probability thresholds
    threshold_seq<-seq(0.05,0.95,by=0.05)
    # accuracy_m<-c("kappa","PCC","auc","typeI.error","typeII.error")
    accuracy_m<-c("PCC","auc","typeI.error","typeII.error")
    
    y.acc<-data.frame(test=accuracy_m)
    
    for(l in 1:length(threshold_seq)){
      
      A<-accuracy(ifelse(predict(Train.m,dat.new,type="response") >= threshold_seq[l], 1, 0),y.new)[accuracy_m]
      b <- A %>% unlist() %>% as.data.frame() 
      b$test<- rownames(b)
      
      y.acc <- merge(y.acc,b,by="test",all=TRUE) ; colnames(y.acc)[ncol(y.acc)] <- threshold_seq[l]
      
      rm(A,b)
    }
    
    y.acc[is.na(y.acc)]<-1 ; rownames(y.acc)<-y.acc$test
    y.acc<-y.acc[,-1] %>% as.matrix()
    
    # 5.2 Plot the accuracy metrics for the different tresholds
    if(plot.maxent==TRUE){
    par(mar=c(5,5,3.5,3)) # c(bottom, left, top, right)
    par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
    
    col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
    color.r <- col.fun(nrow(y.acc)) 
    
    plot(y=1,x=1,xlim=c(0,1),ylim=c(0,1.1),
         axes=F,
         ylab="Parameter performance",
         xlab="Probability threshold",
         type="n")
    
    mtext(text = "a) Model Accuracy and performance",adj=0,line=1)
    
  for(ll in 1:nrow(y.acc)){
    if(row.names(y.acc)[ll] %in% "PCC"){
      lines(y=y.acc[ll,]/100,x=colnames(y.acc) %>% as.numeric(),
              col=color.r[ll],type="b",pch=19)  
        
    }else{
      lines(y=y.acc[ll,],x=colnames(y.acc) %>% as.numeric(),
            col=color.r[ll],type="b",pch=19)
      }
  }
    abline(v=kappa_threshold,lty=3,col="tomato3",xpd=FALSE)
    
    axis(1,seq(0,1,by=0.05),labels =seq(0,1,by=0.05) %>% round(digits = 1) ,gap.axis = 2,col="white")
    axis(2,seq(0,1,by=0.05),labels =seq(0,1,by=0.05) %>% round(digits = 1) ,gap.axis = 2,las=2,col="white")
    
    legend("top",legend=rownames(y.acc),lty=1,pch=19,
           col=color.r,bty="n",horiz=TRUE,x.intersp = 0.8)
    }
    # 5.3 Check the variable importance and variable contribution
    if(is.null(jackknife)){
      jackniffe_var_contrybution<-jack.maxent(Train.m@results,jack=FALSE,plot.j=plot.maxent)
      }else{
      jackniffe_var_contrybution<-jack.maxent(Train.m@results,jack=TRUE,plot.j=plot.maxent)
      }
     
    # 5.2 Run some predictions on the data
    pred.m<-terra::predict(Train.m,predictors)
    
    if(!is.null(clip)){
      pred.m %>% mask(clip %>% vect())
    }

    # 6. Store the results in the local folder
    mx_mod<-list(Train_model=Train.m,
                 Test=y.acc,
                 Test_evaluate=Test.m,
                 Threshold_kappa=kappa_threshold,
                 Var_contrib=jackniffe_var_contrybution,
                 prediction=pred.m,
                 MAXENT_args=argsMX)
    
    models<-mx_mod
    
    # 7. Plot the model
  if(plot.maxent==TRUE){ 
       par(mar=c(0,0,0,0)) # c(bottom, left, top, right)
       plot(pred.m,col=col.fun(550),axes=F)
        mtext("b) Model predictions",adj=0)
      }
  }
  
  if(return==FALSE){
      # Creates the Temporal file folder and the folder to allocated the results----
      if(!dir.exists(results.maxent)){dir.create(results.maxent)}
    
      results.sp<-paste(results.maxent,mod.name,sep="/") ; dir.create(results.sp,recursive = TRUE,showWarnings = FALSE)
      save(models, file=paste(results.sp,paste0("Maxent_",mod.name,".Rdata"),sep="/"))
      return(paste(mod.name,"finish!"))
  }else{
        return(models)
      }
} 

# _____/------\    ______  /------\\\  \\
#/             \__|      \/        / \\\
#                                      //|/\/\/\/_____---
#     Function to plot Maxent       \/  __  /\\\\\
#       Analysis and testing         /\  ___  /\\\\ \\\
#                                  _      \/   |/\\ \  \ \
#        ____           / |     |\  \  \\ \\\\ 
#_______/    \_________/  |____/  \\\\\  \
#
#
jack.maxent<-function(x,jack=TRUE,plot.j=TRUE){
  
  # Process the raw data from maxent
  y.x <- x %>% as.data.frame()
  y.x <- y.x %>% mutate(.before = 1,Term=row.names(y.x))
  
  # A) contribution and permutation performance
  # Contribution
  contribution.x <- y.x[y.x$Term %>% grepl(pattern = ".contribution$"),]
  contribution.x <- contribution.x[order(contribution.x$V1,decreasing = FALSE),]
  contribution.x$Term <- contribution.x$Term %>% gsub(pattern=".contribution$",replacement="") 
  names(contribution.x)[2] <- "Contribution"
  
  # Permutation importance
  permutation.x <- y.x[y.x$Term %>% grepl(pattern = ".permutation.importance$"),]
  permutation.x$Term <- permutation.x$Term %>% gsub(pattern = ".permutation.importance$",replacement="") 
  names(permutation.x)[2] <- "Permutation"
  
  ycontr.x <- merge(contribution.x,permutation.x,by="Term") ; ycontr.x<-ycontr.x[order(ycontr.x$Contribution),]
  rownames(ycontr.x) <- ycontr.x$Term ; ycontr.x <- ycontr.x[,-1] %>% t()
  
  # A) Variable contribution and
  # par(mfrow=c(1,2))
  if(plot.j==TRUE){
  par(mar=c(5,7,3.5,5)) # c(bottom, left, top, right)
  par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
  
  col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
  
  x<-barplot(ycontr.x,horiz=TRUE,border=NA,beside=TRUE,
             axes=F,col=col.fun(nrow(ycontr.x)),las=2,cex.names=0.8)
  
  legend("right",legend=rownames(ycontr.x),
         col=col.fun(nrow(ycontr.x)),pch=15,bty="n")
  
  axis(1,col="white") ; abline(v=5,lty=3,lwd=1,col="tomato3",xpd=F)
  axis(2,at=apply(x,2,mean),col="white",line = NA,tick=TRUE,labels=NA,lwd=0,lwd.ticks=1)
  mtext("a) Analysis of variable Contributions and Permutation Importance",side=3,adj=0.9,las=1,xpd=TRUE)
  }
  if(jack==TRUE){
  
  # Jacknife
  Training.x <- y.x[y.x$Term %>% grepl(pattern = "Training."),]
  jack01 <- Training.x[Training.x$Term %>% grepl(pattern = "only."),] ; names(jack01)[2] <- "Only"
  jack01$Term <- jack01$Term %>% gsub(pattern="Training.gain.with.only.",replacement="")
  
  jack02 <- Training.x[Training.x$Term %>% grepl(pattern = "without."),] ; names(jack02)[2] <- "Without" 
  jack02$Term <- jack02$Term %>% gsub(pattern="Training.gain.without.",replacement="")
  
  jack03 <- merge(jack01,jack02,by="Term") ; rownames(jack03) <- jack03$Term ; jack03 <- jack03[order(jack03$Only,decreasing = FALSE),]
  jack04 <- jack03[,-1] %>% t()
  
  # B) jacknife results
  if(plot.j==TRUE){
  par(mar=c(5,7,3.5,5)) # c(bottom, left, top, right)
  par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
  
  col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
  par(mar=c(5,7,3.5,5)) # c(bottom, left, top, right)
  barplot(jack04,beside=TRUE,horiz=TRUE,
          las=2,col=color.r[c(1,2)],border=NA,axes=F,
          xlab="Regularized Training Gain")
  axis(1,col="white",line=0.5)
  mtext("b) Jackknife or regularized training gain for the MAXENT variables",side=3,adj=0.9,las=1,xpd=TRUE)
  
  legend("right",legend=row.names(jack04),
         col=color.r[c(1,2)],border = F,pch=15,inset = -0.2, xpd=TRUE,cex=0.8,bty="n")
    }
  }
  # Returnt the data
    if(jack==TRUE){
    return(list(jackknife=jack04,Var_contribution=ycontr.x))
  
    }else{
    return(list(Var_contribution=ycontr.x))  
  }
}

#
# End of The function
#