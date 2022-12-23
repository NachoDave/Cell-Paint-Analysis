require(tidyr)

## Functions related to SOMs and minimal spanning trees

makeMSTFromSOM <- function(SOM_map){
  
  ## Function for making MST froms SOMs
  #browser()
  SOM_Dists <- dist(SOM_map$codes[[1]])
  SOM_dst_mat <- matrix(ncol = 3, nrow = length(SOM_Dists))
  
  rowSOM_N <- nrow(as.matrix(SOM_Dists)) 
  colSOM_N <- ncol(as.matrix(SOM_Dists)) 
  
  rowSOM <- unlist(sapply(2:rowSOM_N, function(x){x:rowSOM_N}))
  colSOM <- unlist(sapply(1:(rowSOM_N - 1), function(x){rep(x,(rowSOM_N - x) ) }))
  
  cnt <- 1
  
  for (dx in 1:length(SOM_Dists)){
    
    SOM_dst_mat[dx, ] <- c(colSOM[dx], rowSOM[dx], SOM_Dists[dx])
    
  }
  
  ## Make spanning matrix
  SOM_G <- graph.data.frame(SOM_dst_mat[,1:2] ,directed = FALSE)
  plot(SOM_G, asp = 1)
  
  E(SOM_G)$weight = SOM_dst_mat[,3]       # specify the actual weight
  E(SOM_G)$label  = E(SOM_G)$weight # labels for the weight
  V(SOM_G)$size = 15            # node size
  plot(SOM_G, asp=0)
  
  SOM_MST <- minimum.spanning.tree(SOM_G)
  #plot(SOM_MST)
  return(SOM_MST)
  
}

## Plotting the numbers of each class in each group 
plotSOM_nodes <- function(SOM_map, dat, grp = c("Treatment", "Dose")){
  
  ## 
  #browser()
  TrtDose <- dat[, grp] %>% unite("TreatDose", grp, sep = "_")
  pltDf <-  data.frame(TreatDose = TrtDose, SOM_Class = SOM_map$unit.classif, stringsAsFactors = T)
  pltDfGrp <- pltDf %>% group_by(TreatDose, SOM_Class) %>% count()# %>% group_by("TreatDose") %>% mutate(Percent = 100*n/sum(n))
  
  ppp <- lapply(unique(pltDfGrp$TreatDose), function(x){100*pltDfGrp[pltDfGrp$TreatDose == x,]$n/sum(pltDfGrp[pltDfGrp$TreatDose == x,]$n)})
  pltDfGrp$Percent <- unlist(ppp)
  
  pltDfGrp
  # plt <- ggplot(pltDfGrp, aes(x = SOM_Class, y = Percent, fill = TreatDose)) + 
  #   geom_col(position = "dodge")
  
  
}

## Create image plots for each SOM node

nucPlot <- function(data, imDr, outDr, prefix = "", fn = "Filename", clses = "P1_CisNucMap$unit.classif", 
                    xLoc = "Location_Center_X_nuc", yLoc = "Location_Center_Y_nuc", imSufix = "_Cis.jpeg", objectNoDx = "ObjectNumber_nuc"){
  #browser()
  # Make image directories
  for (ix in unique(data[, clses])){
    
    dir.create(paste0(outDr, prefix, "_", ix))
    
  }
  # Make image patches
  for (dx in 1:nrow(data)){
    #browser()
    I <- load.image(file.path( imDr, gsub(".tif", imSufix,data$Filename[dx])))
    ## Nuclei images
    I1 <- imsub(I, x > round(max(1, data[dx,xLoc])) - 30 &
                  x < round(min(1024, data[dx,xLoc])) + 30,
                y > round(max(1, data[dx,yLoc])) - 30 &
                  y < round(min(1024, data[dx,yLoc])) + 30
    )
    #browser()
    save.image(I1, paste0(outDr, "/", prefix, "_", data[dx, clses],  
                          "/", gsub(".tif", paste0("NucDx_",
                                                   data[dx, objectNoDx],  imSufix
                          ),data$Filename[dx])
    ))
    
    ## Cell images
    I2 <- imsub(I, x > round(max(1025, 1024 + data[dx,xLoc])) - 50 &
                  x < round(min(2*1024, 1024+data[dx,xLoc])) + 50,
                y > round(max(1, data[dx,yLoc])) - 50 &
                  y < round(min(1024, data[dx,yLoc])) + 50
    )
    
    save.image(I2, paste0(outDr, "/", prefix, "_", data[dx, clses],  
                          "/", gsub(".tif", paste0("Cell_NucDx_",
                                                   data[dx, objectNoDx],  imSufix
                          ),data$Filename[dx])
    ))
    
    
  }
  
}
