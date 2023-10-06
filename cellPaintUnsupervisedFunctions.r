require(tidyr)
require(imager)
require(purrr)

## Functions related to SOMs and minimal spanning trees

makeMSTFromSOM <- function(SOM_map, distMethod = "euclidean"){
  
  ## Function for making MST froms SOMs
  #browser()
  
  if(class(SOM_map) == "kohonen"){
    
    SOM_map <- SOM_map$codes[[1]]
    
  }
  
  SOM_Dists <- dist(SOM_map, method = distMethod)
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
                          objectNoDx = "ObjectNumber_nuc", smpN = 10, imSegSuffix = ".jpeg", xYBd = 100, setSeed = FALSE, seed = 1){
  
  if (!assertthat::are_equal(nrow(metaData), length(SOM_class))){
    
    print("The ncol of metaData is no equal to the length of the SOM class vector")
    return()
  }
  
  if (!assertthat::are_equal(nrow(metaData), nrow(nucLocs))){
    
    print("The ncol of metaData is no equal to the nrow of the nuclei location data")
    return()
  }  
  if (setSeed) set.seed(seed)
  
  SOM_classN <- length(unique(SOM_class))
  
  ## Create directories
  for (ix in 1:SOM_classN){
    if (!dir.exists(paste0(outDr, "/",prefix, "_", ix))){
      dir.create(paste0(outDr, "/",prefix, "_", ix))
    }
  }
  
 # browser()
  
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

### Create summary plots of number of cells and percentage in eacg cluster
pltCellsInClusters <- function(dat, SOMVals, somCls, metCols, savDr){
  
  xx <- hclust(dist(SOMVals, "canberra"))
  yy <- cutree(xx, k = 15)
  
  # Data.frame of SOM/cluster classes
  WT_Cis_FiltParmMapSOM36_Class <- cbind(dat[,metCols], som_class = somCls)
  WT_Cis_FiltParmMapSOM36_Class$Cluster <- 0
  
  ## Get the mapping between SOM classes and clusters
  somClustListFiltParm <- yy
  
  myLocIC50WT <- 3 + 1.43/3
  myLocIC50Cis <- 6 + 8000/(40000-20000)
  
  for (idx in 1:length(somClustListFiltParm)){
    WT_Cis_FiltParmMapSOM36_Class$Cluster[WT_Cis_FiltParmMapSOM36_Class$som_class == as.numeric(gsub("V", "", names(somClustListFiltParm)[idx]))] <- somClustListFiltParm[idx]
  }
  
  ## Get possible data columns and number of cells in each
  cellsByGroup <- WT_Cis_FiltParmMapSOM36_Class[,c("Dose", "Treatment", "plate", "Metadata_Well_nuc", "cell_line") ] %>% count(Dose, Treatment, plate, Metadata_Well_nuc, cell_line, name = "total_cells")# %>% group_by(Dose, cell_line, Treatment) %>% summarise(n = sum(n))   
  
  ## Plot some stuff  
  
  for (px in sort(unique(WT_Cis_FiltParmMapSOM36_Class$Cluster))){
    
    tDf <- WT_Cis_FiltParmMapSOM36_Class[WT_Cis_FiltParmMapSOM36_Class$Cluster == px,]
    
    ## Count cells
    tDsCnts <- tDf %>% count(Dose, Treatment, plate, Metadata_Well_nuc, cell_line)
    
    tCnts <- cellsByGroup
    tCnts <- dplyr::left_join(tCnts, tDsCnts, by = c("Dose", "Treatment", "Metadata_Well_nuc", "plate", "cell_line"))
    tCnts$n[is.na(tCnts$n)] = 0
    
    tCnts$percent <- 100*tCnts$n/cellsByGroup$total_cells
    
    ## Get mean and standard deviations of wells for each treatment
    
    tDsPltCnts <- tCnts %>% group_by(Dose, cell_line, Treatment) %>% summarise(mnN = mean(n), sdN = sd(n), mnPer = mean(percent), sdPer = sd(percent))
    tDsPltCnts$Dose <- factor( tDsPltCnts$Dose, levels = sort(unique(as.numeric(tDsPltCnts$Dose))))  
    
    ## Get percentage of total population for each cell type cisplatin
    doxPerPlt <- ggplot(tDsPltCnts[grepl("Doxo|untreated", tDsPltCnts$Treatment), ], aes(x = (Dose), y = mnPer, fill = cell_line)) + 
      geom_bar(stat = "identity", position= position_dodge(), color = "black") + geom_errorbar(aes(ymin=mnPer-sdPer, ymax=mnPer+sdPer), width=.2,position=position_dodge(.9)) +
      xlab("Dose (nm)") + ylab("Percent of total cells") + theme_bw() + ggtitle("Percentage of cells Doxo treat") +
      scale_fill_discrete(name = "Cell line")  + 
      geom_vline(xintercept = myLocIC50WT, linetype = "dashed", linewidth = 1.2) + geom_vline(xintercept = myLocIC50Cis, linetype = "dashed", linewidth = 1.2) # +
    # annotate("text", x=myLocIC50WT - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="WT Dox IC50", angle=90) + 
    # annotate("text", x=myLocIC50Cis - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="Cis Dox IC50", angle=90)
    
    
    cisPerPlt <- ggplot(tDsPltCnts[grepl("Cis|untreated", tDsPltCnts$Treatment), ], aes(x = Dose, y = mnPer, fill = cell_line)) + 
      geom_bar(stat = "identity", position= position_dodge(), color = "black") + geom_errorbar(aes(ymin=mnPer-sdPer, ymax=mnPer+sdPer), width=.2,position=position_dodge(.9)) +
      xlab("Dose (nm)") + ylab("Percent of total cells") + theme_bw() + ggtitle("Percentage of cells Cis treat") +
      scale_fill_discrete(name = "Cell line")  + 
      geom_vline(xintercept = myLocIC50WT, linetype = "dashed", linewidth = 1.2) + geom_vline(xintercept = myLocIC50Cis, linetype = "dashed", linewidth = 1.2) +
      annotate("text", x=myLocIC50WT - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="WT Cis IC50", angle=90) + 
      annotate("text", x=myLocIC50Cis - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="Cis Cis IC50", angle=90) 
    
    
    ## Get cell counts for each group 
    ### Cisplatin
    cntsWTCis <- ggplot(tDsPltCnts[grepl("Cis|untreated", tDsPltCnts$Treatment) & grepl("WT", tDsPltCnts$cell_line), ], aes(x = Dose, y = mnN)) + 
      geom_bar(stat = "identity", position= position_dodge(), color = "black", fill = "#00bfc4") + geom_errorbar(aes(ymin=mnN-sdN, ymax=mnN+sdN), width=.2,position=position_dodge(.9)) +
      xlab("Dose (nm)") + ylab("Number of cells") + theme_bw() + ggtitle("No of cells Cis treat WT") + 
      scale_fill_discrete(name = "Cell line")  + 
      geom_vline(xintercept = myLocIC50WT, linetype = "dashed", linewidth = 1.2) + geom_vline(xintercept = myLocIC50Cis, linetype = "dashed", linewidth = 1.2) +
      annotate("text", x=myLocIC50WT - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="WT Cis IC50", angle=90) + 
      annotate("text", x=myLocIC50Cis - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="Cis Cis IC50", angle=90)       
    
    cntsCisCis <- ggplot(tDsPltCnts[grepl("Cis|untreated", tDsPltCnts$Treatment) & grepl("Cis", tDsPltCnts$cell_line), ], aes(x = Dose, y = mnN)) + 
      geom_bar(stat = "identity", position= position_dodge(), color = "black", fill = "#f8766d") + geom_errorbar(aes(ymin=mnN-sdN, ymax=mnN+sdN), width=.2,position=position_dodge(.9)) +
      xlab("Dose (nm)") + ylab("Number of cells") + theme_bw() + ggtitle("No of cells Cis treat Cis") + 
      scale_fill_discrete(name = "Cell line")  + 
      geom_vline(xintercept = myLocIC50WT, linetype = "dashed", linewidth = 1.2) + geom_vline(xintercept = myLocIC50Cis, linetype = "dashed", linewidth = 1.2) #+
    # annotate("text", x=myLocIC50WT - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="WT Cis IC50", angle=90) + 
    # annotate("text", x=myLocIC50Cis - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="Cis Cis IC50", angle=90)     
    
    ### Doxorubicin
    cntsWTDox <- ggplot(tDsPltCnts[grepl("Dox|untreated", tDsPltCnts$Treatment) & grepl("WT", tDsPltCnts$cell_line), ], aes(x = Dose, y = mnN)) + 
      geom_bar(stat = "identity", position= position_dodge(), color = "black", fill = "#00bfc4") + geom_errorbar(aes(ymin=mnN-sdN, ymax=mnN+sdN), width=.2,position=position_dodge(.9)) +
      xlab("Dose (nm)") + ylab("Number of cells") + theme_bw() + ggtitle("No cells Doxo treat WT") + 
      scale_fill_discrete(name = "Cell line")  + 
      geom_vline(xintercept = myLocIC50WT, linetype = "dashed", linewidth = 1.2) + geom_vline(xintercept = myLocIC50Cis, linetype = "dashed", linewidth = 1.2) +
      annotate("text", x=myLocIC50WT - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="WT Cis IC50", angle=90) + 
      annotate("text", x=myLocIC50Cis - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="Cis Cis IC50", angle=90)       
    
    cntsCisDox <- ggplot(tDsPltCnts[grepl("Dox|untreated", tDsPltCnts$Treatment) & grepl("Cis", tDsPltCnts$cell_line), ], aes(x = Dose, y = mnN)) + 
      geom_bar(stat = "identity", position= position_dodge(), color = "black", fill = "#f8766d") + geom_errorbar(aes(ymin=mnN-sdN, ymax=mnN+sdN), width=.2,position=position_dodge(.9)) +
      xlab("Dose (nm)") + ylab("Number of cells") + theme_bw() + ggtitle("No of cells Dox treat Cis") +
      scale_fill_discrete(name = "Cell line")  + 
      geom_vline(xintercept = myLocIC50WT, linetype = "dashed", linewidth = 1.2) + geom_vline(xintercept = myLocIC50Cis, linetype = "dashed", linewidth = 1.2) #+
    # annotate("text", x=myLocIC50WT - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="WT Cis IC50", angle=90) + 
    # annotate("text", x=myLocIC50Cis - 0.15, y=2*max(tDsPltCnts$mnPer)/11, label="Cis Cis IC50", angle=90)      
    
    plt <- cisPerPlt + doxPerPlt + cntsWTCis + cntsWTDox + cntsCisCis + cntsCisDox + plot_layout(ncol = 2, guides = "collect")
    ggsave(paste0(savDr, "Cluster_", px, ".png"), plt, width = 8, height = 12)
    
  }  
}

## make composite plots from clusters

compPlotsFromClusters <- function(imagePth, imPrefix, clusters, outPth, numberOfImagesFromEachSom = 2){
  
  n <- max(clusters)
  #browser()
  for (cx in 1:n){
    compI <- NULL
    tClsters <- gsub("V", "", names(clusters[clusters == cx]))

    for (sx in 1:length(tClsters)){

      
    pth <- paste0(imagePth, imPrefix, tClsters[sx]) # get the path for the SOM node in the cluster
    fn <- list.files(pth) # get the file names from the directory   
          
    smp <- sample(length(fn), numberOfImagesFromEachSom)
      for(ix in 1:numberOfImagesFromEachSom){
        ## load image
  
        iPth <- file.path(pth, fn[smp[ix]])
        i <- load.image(iPth)
        ## Add title to image
        i <- draw_text(i, 10, 10, gsub(".png", "", fn[smp[ix]]), "white", fsize = 10)
        
        ## Append image
        if (is.null(compI)){
          compI <- i
        }else{
          
          compI <- imappend(list(compI, i), "y")
          
          }
        } # end loop through number of images for each SOM
    } ## end loop through classes in each cluster
    #browser()
    ## Save Composites images
    save.image(compI, paste0(outPth, imPrefix, "_comp_", cx, "SOM", paste(tClsters, collapse = "_"), ".png"))
    
    }
}


