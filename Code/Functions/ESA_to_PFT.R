# Function to transform ESA maps to PFT
ESAtoPFT<-function(map,        # ESA derived raster
                   cross_walk, # Table containing the ESA-PFT values
                   pixels,     # Data.frame containing the equivalence between ESA pixels values and PFT values
                   group){     # If PFT should be grouped or not
  
  # a. Loading require packages or install them----
  list.of.packages<-c("tidyr","raster","rgdal","maptools","data.table","dplyr")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages)
  
  # b. Define de functions that are needed to process the data----
  # b.1 Create a list of PFT empty matrixs one for each type of PFT
  matrix_fun<-function(rast_matrix, # raster redived matrix
                       PFT # vector contaning PFT groups names
  ){
    
    PFT_matrix<-list()
    for(l in 1:length(PFT)){
      
      PFT_matrix[[l]]<-matrix(data=NA,nrow=nrow(rast_matrix),ncol=ncol(rast_matrix))
      names(PFT_matrix)[l]<-PFT[l]
      
    }
    return(PFT_matrix)
  } 
  
  
# b.2 Tranform the values of the raster to their corresponding covers  
  Transfer_px<-function(rast_matrix, # matrix derived from the population raster
                        matrix_PFT, # list containing the empty list of PFT matrices
                        pixels, # cross_pixel data.frame containing the $PFT names and the $Pixel_values for each type of PFT (derived from the cross walk table) 
                        values  # The cross walk table containing the $Pixel_value and its corresponding land-cover value for each PFT and ESA class
  ){
    
    px_vals<-unique(as.numeric(rast_matrix)) # unique pixel values in the raster
    px_vals
    
    PFT_covers<-matrix_PFT
    
    for(p in 1:length(px_vals)){
      
      for (w in 1:nrow(pixels)){
        
        x<-strsplit(pixels$Pixel_value[w],split=",") %>% unlist() %>% as.numeric()
        y<-pixels$PFT[w]
        
        if(px_vals[p] %in% x){
          
          z<-values[values$Pixel_value==px_vals[p],colnames(values)%in% y ]
          PFT_covers[[y]][rast_matrix==px_vals[p]]<-z/100
          
        }else{next}
        
      }
    }
    
    # Transform into a raster stack
    rast_01<-lapply(PFT_covers,FUN=raster) %>% stack()
    #return(unique(values(sum(rast_01,na.rm=TRUE))))
    return(rast_01)
  }
  
  # c.1 Run the functions with the parameters----
  area_cells<-raster::area(map)%>%as.matrix()
  
  m<-as.matrix(map)
  m_pft<-matrix_fun(rast_matrix = m,PFT=pixels$PFT)
  
  m_01<-Transfer_px(rast_matrix=m, 
                        matrix_PFT=m_pft, # list containing the empty list of PFT matrices
                        pixels=pixels, # cross_pixel data.frame containing the $PFT names and the $Pixel_values for each type of PFT (derived from the cross walk table) 
                        values=cross_walk)  
  
   # c.2 Prepare the data to export----
  if(missing(group)){
    
   # c.2.1. Grouped the data into 
    
    g_PFT<-stack(
        raster::subset(m_01,grep(names(m_01),pattern = "Tree")) %>% sum(na.rm=TRUE),
        raster::subset(m_01,grep(names(m_01),pattern = "Shrub")) %>% sum(na.rm=TRUE),
        m_01$Natural_Grass,
        m_01$Crop,
        m_01$Bare.Soil,
        m_01$Snow.Ice , 
        m_01$Water,
        m_01$Urban,
        m_01$No.data
      )
    
    groups<-c("Forest_trees","Shrubs","Grass",
              "Crop","Bare_Soil","Snow_ice",
              "Water","Urban","No_data")
    
    for (i in 1:length(groups)){
      
      names(g_PFT[[i]])<-groups[i] # set the names of the groups
      g_PFT[[i]]<-raster(area_cells)*g_PFT[[i]] # Calculate the area in km^2 for the cover in each cell
     
    }
    
    # c.2.2 Calculate the total percentage of cover
    cover<-data.frame(PFT_group=groups,cover_percentage=NA)
    
    for (l in 1:length(groups)){
    
      cover$cover_percentage[l]<-sum(as.matrix(g_PFT[[l]]),na.rm=TRUE)/sum(area_cells)
    
    }
    
    return(cover)
    
  }else{
    
  for (i in 1:nlayers(m_01)){
      
    m_01[[i]]<-raster(area_cells)*m_01[[i]] # Calculate the area in km^2 for the cover in each cell
      
    }  
  
    
  cover<-data.frame(PFT_original=names(m_01),cover_percentage=NA)
  groups<-cover$PFT_original
  
  for (f in 1:nlayers(m_01)){
    
    cover$cover_percentage[l]<-sum(as.matrix(m_01[[l]]),na.rm=TRUE)/sum(area_cells)
    
  }
  
  return(cover)  
    
  }

}
