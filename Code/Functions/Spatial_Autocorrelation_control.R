# Function to perform a mantel correlogram to check spatial autocorrelation
# adapted from ecodist;
# Goslee, S.C. and Urban, D.L. 2007. The ecodist package for dissimilarity-based analysis of
# ecological data. Journal of Statistical Software 22(7):1-19. DOI:10.18637/jss.v022.i07
#
#
#' @export mantel_correlogram
#' @title Mantel correlogram from random samples
#' 
#' @description \code{mantel_correlogram} applies a mantel test over a random sample of spatial point records and returns a correlogram
#' # for the distances specified
#' 
#' @param x data frame with points coordinates and environmental variables
#' @param colxy names of the columns containing the coordinates
#' @param n number of random occurrences for the test
#' @param env.vars names of the columns containing the environmental variables
#' @param max.dist maximun distance, in coordinates units at which to evaluate the spatial-autocorrelation
#' @param distances a list of distances to evaluate spatial autocorrelation (it supress max.dist)
#' @param n.dist number of distances to evaluate 
#' @param nperm number of permutations to execute the mantel test
#' @return mantel_corr: A list object of containing the results from the mantel test and the informaiton for plotting

mantel_correlogram <- function (
          x, # data.frame with environmental varialbes 
          colxy, # columns with XY coordinates
          n, # Number of random occurrences for the test 
          env.vars, # variables
          max.dist, # maximun distance in coordinates units (UTM='m',LonLat=degrees)
          distances=NULL, # A list of distances to evaluate spatial autocorrelation (it supress max.dist)  
          n.dist=10, # number of distances to evaluate the correlation
          nperm=1000 # permutation number)
          ){
    
    # a. Load and install the required packages  
    list.of.packages<-c("dplyr","ecodist")
  
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)
  
    lapply(list.of.packages,require,character.only=TRUE)
    rm(list.of.packages,new.packages)
  
    # b. Subset the information                 
    y <- x[,colnames(x) %in% env.vars]
    y.s <- x[,colnames(x) %in% colxy]        
    
    # c. Normalize the environmental variables----
    y.m <- apply(y,2,mean) # get the mean
    y.sd <- apply(y,2,sd) # Standard deviation
    
    var_norm <- data.frame(t((t(y)-y.m)/y.sd))
    
    # d. Take a random 'n' sample of observations----
    obs_sample <- sample(1:nrow(x), n, replace = TRUE)
    
    # e. Calculate the distance matrices for the records and variables----
    var_dist <- dist(var_norm[obs_sample, ],diag = FALSE, upper = FALSE)
    rec_dist <- dist(y.s[obs_sample,],diag = FALSE, upper = FALSE)
  
    # f. Get a secuence of distances to evaluate the autocorrelation
    if(is.null(distances)){
    b <- seq(from = min(rec_dist), to = max.dist, length.out = n.dist) %>% round(digits=0)
    }else{
    b <- c(min(rec_dist),distances)  
    }
    
    # g. Calculate the mantel correlogram using the ecodist package
    mantel_corr <- ecodist::mgram(species.d=var_dist, space.d=rec_dist, 
                           breaks = b,nboot = 500,
                           nperm = nperm)
    
    f <- mantel_corr$mgram %>% as.data.frame()
    
    ecodist::plot.mgram(mantel_corr) ; abline(h=0,lty=3,col="tomato")
    abline(h = 0)
    
    return(mantel_corr)
  
}

#
# End of the function
#