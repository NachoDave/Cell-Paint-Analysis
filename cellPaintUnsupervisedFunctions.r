require(tidyr)
require(imager)
require(purrr)

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


### Create images of cells from SOM results

makeImagePatches <- function(metaData, SOM_class, nucLocs, outDr, prefix,
                          xLoc = "Location_Center_X_nuc", yLoc = "Location_Center_Y_nuc",
                          objectNoDx = "ObjectNumber_nuc", smpN = 10, imSegSuffix = ".jpeg", xYBd = 100, setSeed = FALSE){
  
  if (!assertthat::are_equal(nrow(metaData), length(SOM_class))){
    
    print("The ncol of metaData is no equal to the length of the SOM class vector")
    return()
  }
  
  if (!assertthat::are_equal(nrow(metaData), nrow(nucLocs))){
    
    print("The ncol of metaData is no equal to the nrow of the nuclei location data")
    return()
  }  
  if (setSeed) set.seed(25)
  
  SOM_classN <- length(unique(SOM_class))
  
  ## Create directories
  for (ix in 1:SOM_classN){
    if (!dir.exists(paste0(outDr, "/",prefix, "_", ix))){
      dir.create(paste0(outDr, "/",prefix, "_", ix))
    }
  }
  
  #browser()
  
  ## Open the images for each cell
  for (dx in 1:SOM_classN){ # loop through each SOM_class
    smp <- sample(sum(SOM_class == dx), min(smpN, sum(SOM_class == dx)))
    
    tMet <- metaData[SOM_class == dx, ][smp,]
    tLoc <- nucLocs[SOM_class == dx, ][smp,]
    
    for(ix in 1:nrow(tMet)){
      ## Tiff images --------------------------- ##
      
      xMin <- floor(max(1, tLoc[ix,xLoc] - xYBd))
      xMax <- ceiling(min(1024, tLoc[ix,xLoc] + xYBd))
      yMin <- floor(max(1, tLoc[ix,yLoc] - xYBd))
      yMax <- ceiling(min(1024, tLoc[ix,yLoc] + xYBd))
      
      # Nuclei
      I_DNA <- load.image(file.path(tMet$Imagepath[ix], tMet$FilenameDNA[ix]))
      I_DNA_Patch <- imsub(I_DNA, x > xMin &
                                   x < xMax,
                                 y > yMin &
                                   y < yMax
      )
      
      I_DNA_Patch <- imappend(list(I_DNA_Patch, I_DNA_Patch, I_DNA_Patch), "c")
      
      # Cytoplasm
      I_Cyto <- load.image(file.path(tMet$Imagepath[ix], tMet$FilenameCyto[ix]))
      I_Cyto_Patch <- imsub(I_Cyto, x > xMin &
                              x < xMax,
                            y > yMin &
                              y < yMax
      )
      
      I_Cyto_Patch <-  imappend(list(imfill(dim(I_Cyto_Patch)[1],dim(I_Cyto_Patch)[2],1, val = 0), I_Cyto_Patch, imfill(dim(I_Cyto_Patch)[1],dim(I_Cyto_Patch)[2],1, val = 0)), "c") # convert to RGB
      
      I_DNA_Cyto_Patch <- imappend(list(I_DNA_Patch, I_Cyto_Patch*2.5), "x")
      I_DNA_Cyto_Patch <- I_DNA_Cyto_Patch/max(I_DNA_Cyto_Patch)
      
      #browser()
      ## Segmentation Images --------------------------- ##
      
      segPth <- file.path(tMet$Segmentationpath[ix], gsub(".tif", paste0("_", tMet$cell_line[ix], imSegSuffix), tMet$FilenameDNA[ix]))
      
      if (file.exists(segPth)){ ## check that the segmentation images exists
        I_Seg <- load.image(segPth)
        I_Seg_DNA_Patch <- imsub(I_Seg, x > xMin &
                                   x < xMax,
                                 y > yMin &
                                   y < yMax
        )
        I_Seg_Cyto_Patch <- imsub(I_Seg, x > floor(max(1025, 1024 + tLoc[ix,xLoc] - xYBd)) &
                                  x < ceiling(min(2048, 1024 + tLoc[ix,xLoc] + xYBd)),
                                y > floor(max(1, tLoc[ix,yLoc] - xYBd)) &
                                  y < ceiling(min(1024, tLoc[ix,yLoc] + xYBd))
        )
        
        I_Seg_DNA_Cyto_Patch <- imappend(list(I_Seg_DNA_Patch, I_Seg_Cyto_Patch), "x")
        I_DNA_Cyto_Patch <- imappend(list(I_DNA_Cyto_Patch, I_Seg_DNA_Cyto_Patch), "y")

      }
      #browser()
      
      ## Write image to file
      svPth <- file.path(outDr, paste0(prefix, "_", dx), paste0(gsub("\\ |\\(|\\))","" , gsub(" wv UV - DAPI\\)\\.tif", "", tMet$FilenameDNA[ix])),
                         "_Ob", tMet$ObjectNumber_nuc[ix], "_", tMet$cell_line[ix], "_", tMet$plate[ix], "_", tMet$Treatment[ix],  gsub("\\.","p" , tMet$Dose[ix]), "_SOM",dx, ".png") )
      
      save.image(I_DNA_Cyto_Patch, svPth)
      #browser()
    }
    
  }
  
  
}

