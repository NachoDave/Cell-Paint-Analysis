require(tidyr)
require(imager)
require(purrr)
require(MOFA2)
library(viridis)
library(grid)

## SOM functions
cell_paint_som <- function(data, datCols, xDim = 6, yDim = 6, seed = NULL){
  
  # General function to perform the SOM analysis

  if (!is.null(seed)){
    
    set.seed(seed)
    
  }
  g <- somgrid(xdim = xDim, ydim = yDim, topo = "rectangular")
  data_som <- som(as.matrix(data[, datCols]), grid = g)

  return(data_som)
}

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

### Create summary plots of number of cells and percentage in eacg cluster
pltCellsInClusters <- function(dat, SOMVals, somCls, metCols, savDr, k = 15){
  
  xx <- hclust(dist(SOMVals, "canberra"))
  yy <- cutree(xx, k = k)
  
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
    
    #plt <- cisPerPlt + doxPerPlt + cntsWTCis + cntsWTDox + cntsCisCis + cntsCisDox + plot_layout(ncol = 2, guides = "collect")
    plt <- cisPerPlt +  cntsWTCis + cntsCisCis + plot_layout(ncol = 1, guides = "collect")
    ggsave(paste0(savDr, "Cluster_", px, ".png"), plt, width = 5, height = 12)
    
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

## Theta mapping

# Function to perfomr theta mapping based on Carragher 2016

thetaMap <- function(){
  
  
  
}

## Get parameter groups
## All features

get_feature_groups <- function(nms, ftGrps, ims, areas, area_shape = "AreaShape"){
  #browser()
  grp <- c()
  feat <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(feat) <- c("group", "feature")
  for (fdx in ftGrps){
    for(idx in ims){
      for(adx in areas){
        tGrp <- paste(fdx, idx, adx, sep = "_")
        tFet <- nms[grepl(paste(fdx, idx, adx, sep = ".*"), nms)]
        ## features
        #browser()
        if (length(tFet) > 0){
          feat <- rbind(feat, data.frame(group = tGrp, feature = tFet )) 
          grp <- c(grp, tGrp)                  
          
        }
      }
    }
  }
  ## Area_shape parameters
  for (adx in areas){
    tGrp <- paste("AreaShape", adx, sep = "_")
    tFet <- nms[grepl(paste("AreaShape", adx, sep = ".*"), nms)]
    
    if (length(tFet) > 0){
      feat <- rbind(feat, data.frame(group = tGrp, feature = tFet )) 
      grp <- c(grp, tGrp)           
      
    }
    
  }
  #browser()
  ## Correlation parameters
  tGrp <- "Correlation"
  tFet <- nms[grepl("^Correlation_(?!Correlation)", nms, perl = TRUE)]
  
  if (length(tFet) > 0){
    feat <- rbind(feat, data.frame(group = tGrp, feature = tFet )) 
    grp <- c(grp, tGrp)             
  }
  
  ## Correlation parameters
  # tGrp <- "Correlation"
  # tFet <- nms[grepl("Correlation_(?!Correlation)", nms, perl = TRUE)]
  # 
  # if (length(tFet) > 0){
  #   feat <- rbind(feat, data.frame(group = tGrp, feature = tFet )) 
  #   grp <- c(grp, tGrp)             
  # }
  
  ## Correlation_correlation parameters
  tGrp <- "Correlation_Correlation"
  tFet <- nms[grepl("Correlation_Correlation", nms)]
  
  if (length(tFet) > 0){
    feat <- rbind(feat, data.frame(group = tGrp, feature = tFet )) 
    grp <- c(grp, tGrp)             
  }
  
  ## Neighbour parameters
  tGrp <- "Neighbors"
  tFet <- nms[grepl("Neighbors", nms)]
  
  if (length(tFet) > 0){
    feat <- rbind(feat, data.frame(group = tGrp, feature = tFet )) 
    grp <- c(grp, tGrp)             
  }
  
  #browser()
  if (!assertthat::are_equal(length(nms), length((feat$feature)))){
    
    print("The number of features in the named features is not equal to the number output")
    
  }
  
  #return(list(grp, feat))    
  rownames(feat) <- feat$feature
    
  return(feat)
}

## MOFA related functions =============================================================================== ##
summarizeLongFeat <- function(x){
  
  rownames(x) <- gsub("Granularity", "Grn", rownames(x))
  rownames(x) <- gsub("RadialDistribution", "RD", rownames(x))
  rownames(x) <- gsub("Texture", "Tx", rownames(x))
  rownames(x) <- gsub("Correlation", "Cr", rownames(x))
  rownames(x) <- gsub("InfoMeas", "InfM", rownames(x))
  rownames(x) <- gsub("Granularity", "Grn", rownames(x))
  rownames(x) <- gsub("CytoSkelCorr", "CySk", rownames(x))
  rownames(x) <- gsub("Intensity", "Inty", rownames(x))
  rownames(x) <- gsub("Difference", "Diff", rownames(x))
  rownames(x) <- gsub("Mitochondria", "Mito", rownames(x))
  rownames(x) <- gsub("AreaShape", "ArShp", rownames(x))
  x
  
}

## MAke cell paint MOFA objects

makeCpMOFAOb <- function(cpDat){
  
  MOFADat <- list()
  
  # DNA parameters
  MOFADat$DNA <- t((cpDat %>% select(matches("DNA_Corr")) %>% select(!matches("^Correlation")))[,-1:-2])
  colnames(MOFADat$DNA) <- paste0("Cell", 1:nrow(cpDat))
  MOFADat$DNA <- summarizeLongFeat(MOFADat$DNA)
  
  # Cytoskeleton parameters
  MOFADat$Cyto <- t((cpDat %>% select(matches("CytoSkelCorr"))%>% select(!matches("^Correlation")))[,-1:-2])
  colnames(MOFADat$Cyto) <- paste0("Cell", 1:nrow(cpDat))
  MOFADat$Cyto <- summarizeLongFeat(MOFADat$Cyto)
  
  #rownames(MOFADat$Cyto) <- myFun(nrow(MOFADat$Cyto))
  
  # ER parameters 
  MOFADat$ER <- t((cpDat %>% select(matches("_ER_"))%>% select(!matches("^Correlation")))[,-1:-2])
  colnames(MOFADat$ER) <- paste0("Cell", 1:nrow(cpDat))
  MOFADat$ER <- summarizeLongFeat(MOFADat$ER)
  
  # Mitochondria
  MOFADat$Mito <- t((cpDat %>% select(matches("Mitochondria"))%>% select(!matches("^Correlation")))[,-1:-2])
  colnames(MOFADat$Mito) <- paste0("Cell", 1:nrow(cpDat))
  MOFADat$Mito <- summarizeLongFeat(MOFADat$Mito)
  
  # AreaShapeNeighbour parameters
  MOFADat$AreaShape <- t((cpDat %>% select(matches("AreaShape"))%>% select(!matches("^Correlation")))[,-1:-2])
  colnames(MOFADat$AreaShape) <- paste0("Cell", 1:nrow(cpDat))
  MOFADat$AreaShape <- summarizeLongFeat(MOFADat$AreaShape)
  
  # Correlation parameters
  MOFADat$Correlation <- t((cpDat %>% select(matches("^Correlation")))[,-1:-2])
  colnames(MOFADat$Correlation) <- paste0("Cell", 1:nrow(cpDat))
  MOFADat$Correlation <- summarizeLongFeat(MOFADat$Correlation)

  grps <- cpDat$cell_line
  
  ## Make the MOFA object...
  
  MOFAobject <- create_mofa(MOFADat, groups = grps)
  
  return(MOFAobject)
  
  
}

## Plot 
plot_data_scatter_dj <- function(MOFA, group = "cell_line", view = "DNA", factor = 1, n = 9, cov = NULL, x = "factor", logX = FALSE, colBrw = FALSE){
  #browser()
  W <- get_weights(MOFA)[[view]][, factor]
  Z <- get_factors(MOFA, factor = factor)$group1
  
  W_sort <- W[order(abs(W), decreasing = T) ][1:n]
  
  features <- t(MOFAOb_WT_Cis_Trained@data[[view]]$group1)
  #browser()
  df <- data.frame( Z, features[,names(W_sort)], MOFA@samples_metadata)
  
  if (!is.null(cov)){
    
    df <- df[df[, group] == cov, ]
    
  }
  
  if (x == "factor"){
    
    x <- paste0("Factor", factor)
    
  } 
  
  if (logX){
    
    df[, x] <- log2(df[, x] + 250)
    
  }
  #browser()
  df_melt <- gather(df, key = "feature", value = "value", colnames(df)[2:(n+1)])
  
  p <- ggplot(data = df_melt, aes_string(x = x, y = "value", color = group)) + 
    geom_point(size = 0.5) +
    facet_wrap(~ feature, scales = "free_y") + stat_smooth(formula = y ~ x, aes_string(color = group), 
                                                           method = "lm", alpha = 0.4) + ggpubr::stat_cor(aes_string(color = group,
                                                                                                                     label = "..r.label.."), method = "pearson", 
                                                                                                          label.sep = "\n", output.type = "latex", size = 3) 
  
  if (colBrw){
    
    p <- p + scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint =  10000)
    
  }
  
  p
  
}

## Image, dendrogram heatmap image function

imageDendroHM <- function(imDir, SOM_map, subDr_prefix, savDr, savFn,  clstCut = 16, imN = 5){
  
  SOM_codes <- SOM_map$codes[[1]]
  # 1 Clustering of SOM nodes
  
  clst <- hclust(dist(SOM_codes, "canberra"))
  clstCut <- cutree(clst, k = cut)
  nodeN <- length(clst$order)
  
  row_phm_cut <- cumsum(sapply(unique(clstCut[rev(clst$order)] ), function(x){sum(clstCut[rev(clst$order)] == x)}))  # for colouring the SOM classes by cluster
  
  # 2 Make and create dendrogram plot
  dend <- as.dendrogram(clst)
  dend_data <- dendro_data(dend, type = "rectangle")
  
  dendCol_df <- cbind(dend_data$labels, data.frame(clst = clstCut[dend_data$labels$label])) # for colouring the SOM classes by cluster
  
  p <- ggplot(dend_data$segments) + 
    geom_segment(aes(x = y, y = x, xend = yend, yend = xend))+
    geom_text(data = dendCol_df, aes(y, x, label = label, color = clst),
              hjust = 0.1, angle = 0, size = 5, vjust = 1.3) + ylim(1, nodeN) + coord_fixed(7*max(dend_data$segments$y)/nodeN) +
    theme_bw() + ylab("") + xlab("") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0,0,-5,-5), "cm")
    ) + scale_color_viridis(discrete = FALSE, option = "E")
  
  # 3 Make the empty plot to display the 
  p2 <- ggplot(dend_data$segments) + 
    #geom_segment(aes(x = y, y = x, xend = yend, yend = xend))+
    xlim(0, imN*0.96) +
    ylim(1, nodeN) + theme_bw() + 
    theme(axis.line = element_line(),
          axis.line.x = element_line(color = "white"),
          axis.title.x = element_text(color = "white"),
          axis.text.x = element_text(color = "white"),
          axis.ticks.x = element_line(color = "white"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,0,-5,-1), "cm"),
          panel.background = element_blank()) +
    coord_fixed()
  
  
  # Get the plot limits of the dendrogram plot
  ylimP <- ylim$layout$panel_params[[1]]$y.range
  
  whiteIm <- as.cimg(rep(1, 3*5*200), x = 5, y = 200, cc = 3) # white image for spacing
  
  # 4 Make the image grobs
  for (idx in 1:nodeN){
    
    cd <- clst$order[idx]
    
    tImDir <- dir(paste0(imDir, subDr_prefix, cd), full.names = TRUE, )
    tImFn <- dir(paste0(imDir, subDr_prefix, cd), full.names = FALSE, )
    
    smp <- sample(length(tImDir), imN)
    
    tIms <- lapply(tImDir[smp], load.image) # load images
    tIms <- lapply(tIms, function(tI){imsub(tI, x > (dim(tI)[1]/2 + 1), y > (dim(tI)[2]/2 + 1))})
    tIms  <- lapply(tIms, function(tI){pad(tI, 200 - dim(tI)[1], pos = 1,  "x")})
    tIms  <- lapply(tIms, function(tI){pad(tI, 200 - dim(tI)[2], pos = 1,  "y")})
    
    for(dx in 1:imN){
      tIms[[dx]] <- draw_text(tIms[[dx]], 10, 10, gsub(".png", "", tImFn[smp[dx]]), "white", fsize =8)
    }
    
    
    
    ttIms <- (vector("list", length = (length(image_list) - 1) * 2 + 1))
    
    for (i in seq_along(tIms)) {
      ttIms[[(i - 1) * 2 + 1]] <- tIms[[i]]
      if (i < length(tIms)) {
        ttIms[[(i - 1) * 2 + 2]] <- whiteIm
      }
    }
    
    
    im <- imappend(ttIms, axis = "x")
    
    # 4 Format image in 
    
    p2 <- p2 + annotation_custom(rasterGrob(im, width = 1, height = 1),
                                 xmin = 0, xmax = imN*0.96,
                                 ymin = idx - 0.78, ymax = idx + 0.18) 
    
  }
  
  # 5 make the heatmap
  
  ### All features
  All_feature_Grps <- get_feature_groups(colnames(SOM_map$data[[1]]), c( "Granularity", "Intensity", "RadialDistribution", "Texture"),
                                         c("DNA", "CytoSkel", "Mitochondria", "ER"), c("nuc", "cell", "cyto"))
  
  All_feature_Grps <- All_feature_Grps[order(All_feature_Grps$group),]
  All_feature_Grps_ColN <- cumsum(table(All_feature_Grps$group))
  
  SOM_Code_HM_plot <- SOM_codes[rev(clst$order) ,All_feature_Grps$feature]
  colnames(SOM_Code_HM_plot) <- paste0("cl", 1:ncol(SOM_Code_HM_plot))
  rownames(All_feature_Grps) <- colnames(SOM_Code_HM_plot)
  
  phm <- pheatmap(SOM_Code_HM_plot, show_colnames =  TRUE, breaks = breaks, color = colors, cluster_rows = F,
                   clustering_distance_rows = "canberra", main = subDr_prefix, cutree_rows = clstCut, cluster_cols = FALSE, annotation_col = All_feature_Grps[, 1, drop = FALSE], 
                   gaps_col = All_feature_Grps_ColN, gaps_row = row_phm_cut
  )
  

  # Join plots together and write to file
  
  pp <- p2 + p + phm[[4]]

  ggsave(pp, filename = paste0(savDr, savFn), 
         height = 10, width = 30)
  
  
}
