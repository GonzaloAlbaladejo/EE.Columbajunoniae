# 
# 
# Add error bards to a basic barplot
error.bar <- function(x, y, upper, lower=upper, length=0.1,angle.x=90,horiz.p=TRUE,...){
  
  if(horiz.p==TRUE){
    arrows(x - lower, y, x + upper, angle=angle.x, code=3, length=length, ...)
  }else{
    yx<-x
    xy<-y
    arrows(xy,yx+upper, xy, yx-lower, angle=angle.x, code=3, length=length, ...)
  }}

# Aggregate the accuracy_data 
acc_data<-function(x,name.m){ 
            w <- x %>% t() %>% as.data.frame() %>% mutate(type=name.m,treshold=colnames(x),.before=1)
            return(w)
}

# Aggregate the variable data
var_data <- function(x,name.m){
  w <- x %>% as.data.frame() %>% mutate(type=name.m,obs=row.names(x),.before=1)
  return(w)
}

# Plot accuracy results
plot.accuracy<-function(x,name.m){
      par(mar=c(5,5,4,5)) # c(bottom, left, top, right)
      par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
      
      plot(y=1,x=1,xlim=c(0,1),ylim=c(0,1.1),axes=F, ylab="Parameter performance",xlab="Probability threshold",type="n")
      mtext(text = name.m,adj=0,line=1)
      pch_test<-c(21:25)

      col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
      color.r <- col.fun(nrow(x))
  
  for(ll in 1:nrow(x)){
    if(row.names(x)[ll] %in% "PCC"){
      lines(y=x[ll,]/100,x=colnames(x) %>% as.numeric(),
            bg=color.r[ll],col=color.r[ll],type="b",pch=pch_test[k])  
    }else{
      lines(y=x[ll,],x=colnames(x) %>% as.numeric(),
            bg=color.r[ll],col=color.r[ll],type="b",pch=pch_test[k])
    }}
  
axis(1,seq(0,1,by=0.05),labels =seq(0,1,by=0.05) %>% round(digits = 1) ,gap.axis = 2,col="white")
axis(2,seq(0,1,by=0.05),labels =seq(0,1,by=0.05) %>% round(digits = 1) ,gap.axis = 2,las=2,col="white")

legend("top",legend=rownames(x),lty=1,pch=19,col=color.r,bty="n",horiz=FALSE,ncol = 5,x.intersp = 0.5,cex=0.8)
}

# Plot variable contributions
plot.contribution <- function(w,legend.inset=-0.3,name.m){
  
  par(mar=c(5,7,3.5,9)) # c(bottom, left, top, right)
  par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
  col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
  
  x<-barplot(w,horiz=TRUE,border=NA,beside=TRUE,xlab="Regularized Training Gain",
             axes=F,col=col.fun(nrow(w)),las=2,cex.names=0.8)
  
  legend("right",legend=rownames(w),inset=legend.inset,
         col=col.fun(nrow(w)),pch=15,bty="n",xpd=TRUE)
  
  axis(1,col="white") ; abline(v=5,lty=3,lwd=1,col="tomato3",xpd=F)
  axis(2,at=apply(x,2,mean),col="white",line = NA,tick=TRUE,labels=NA,lwd=0,lwd.ticks=1)
  mtext(name.m,side=3,adj=0.9,las=1,xpd=TRUE)
}

# Plot raw jackknife results
plot.jack<-function(x,name.m,xlim.m=c(0,1)){
    par(mar=c(5,7,3.5,5)) # c(bottom, left, top, right)
    par(bg="grey35",col.axis = 'white', col.lab = 'white',col="white",col.main="white")
    
    col.fun <- colorRampPalette(c("skyblue","cyan4","#54bb64ff","gold"))
    par(mar=c(5,7,3.5,5)) # c(bottom, left, top, right)
    barplot(x,beside=TRUE,horiz=TRUE,xlim=xlim.m,
            las=2,col=color.r[c(1,2)],border=NA,axes=F,
            xlab="Regularized Training Gain")
    axis(1,col="white",line=0.5)
    mtext(name.m,side=3,adj=0,las=1,xpd=TRUE)
    
    legend("right",legend=row.names(x),
           col=color.r[c(1,2)],border = F,pch=15,inset = -0.2, xpd=TRUE,cex=0.8,bty="n")
}




