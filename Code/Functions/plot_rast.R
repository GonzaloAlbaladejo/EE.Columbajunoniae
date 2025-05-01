# Plot bii information
plot.r<-function(x, # Rast or stack of rast objects
                 color.r=c("skyblue","tomato"), # Colors to construct the legend and raster levels
                 #lines.r=NULL,
                 #col.lines="white",
                 add.r=FALSE,
                 cero=TRUE, # Does the scale need to be centred at cero?
                 range.r=NULL, # Do you want to specify a range of values
                 legend.r=TRUE # Add a legend to the plot?
                 ){
  # set the range of values
  if(is.null(range.r)){
      x.vals <- values(x[[1]])
      x.r <- range(x.vals,na.rm=TRUE) %>% as.numeric()  
      
      if(cero==TRUE){
        b.1 <- seq(from=min(x.r),to=0,length.out=100) ; b.1 <- b.1[b.1!=0]
        b.2 <- seq(from=0,to=max(x.r),length.out=100) ; b.2 <- b.2[b.2 !=0]
      
        breaks.r <- c(b.1,0,b.2)
        
        }else{
        breaks.r <- seq(from=min(x.r),to=max(x.r),length.out=200)
        }
  }else{
    x.r <- range.r
    
    if(cero==TRUE){
      b.1 <- seq(from=min(x.r),to=0,length.out=200) ; b.1 <- b.1[!b.1 %in% 0]
      b.2 <- seq(from=0,to=max(x.r),length.out=200) ; b.2 <- b.2[!b.2 %in% 0]
      
      breaks.r <- c(b.1,0,b.2)
      
    }else{
      breaks.r <- seq(from=min(x.r),to=max(x.r),length.out=200)
    }
  }
  
  # set the color scale
  col.fun<-colorRampPalette(color.r)
  col.r <- col.fun(length(breaks.r))
  
  # Base raster plot
  plot(x,
       breaks=breaks.r,
       col=col.r,
       legend=F,
       box="n",
       axes=F,
       add=add.r,
       frame.plot = FALSE,
       main="")
  
  # if(!is.null(lines.r)){
  #   lines(lines.r %>% vect(),col=col.lines,lwd=0.5)
  #   }
  
  #opar <- par()
  
  # add the legend
  if(legend.r==TRUE){ # in the future potentially adding a paramter  to move the legend around!
  
    # Function to produce a continuous legend
    
    n <- length(col.r)
    bx <- par("usr") # Get the dimensions of the plotting area
    
    bx_range <- abs(diff(bx[3:4])/2) # adjust the height for the legend
    bx_mid <- mean(bx[3:4])
    
    bx[c(3,4)] <- c(bx_mid-(bx_range*0.25),bx_mid+(bx_range*0.25))
    bx[2]<-bx[2]-diff(bx[1:2])*0.04
    
    box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
                bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
    box.cy <- c(bx[3], bx[3])
    box.sy <- (bx[4] - bx[3]) / n
    
    xx <- rep(box.cx, each = 2)
    
    par(xpd = TRUE, new=TRUE)
    
    for(i in 1:n){
      
      yy <- c(box.cy[1] + (box.sy * (i - 1)),
              box.cy[1] + (box.sy * (i)),
              box.cy[1] + (box.sy * (i)),
              box.cy[1] + (box.sy * (i - 1)))
      polygon(xx, yy, col = col.r[i], border = col.r[i])
      
    }
    axis(side = 4,at=seq(bx[3],bx[4],length.out=10) %>% round(digits = 1),
         labels = seq(min(breaks.r),max(breaks.r),length.out=10) %>% round(digits = 1),
         las = 2, tick = TRUE, line =-1,cex=0.7,lwd=NA,lwd.ticks = 0.8)
    
    par(xpd = FALSE, new=FALSE)
    }
  
  # par() <- opar
}


#
# End of the script
#