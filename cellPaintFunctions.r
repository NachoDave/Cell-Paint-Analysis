## Cell Paint functions
## Libaries 
library(data.table)
library(WGCNA)
library(dplyr)
library(ggplot2)

## Load data from csv ================================================================================= ##

load_CP_csv <- function(pth, prefix = "", nucSuf = "Nuclei", cellSuf = "Cell", cytoSuf = "Cytoplasm"){
  
  dt_Nuc <- as.data.frame(fread(paste0(pth, "/", prefix, nucSuf, ".csv"))) # load nuclei data
  colnames(dt_Nuc) <- paste(colnames(dt_Nuc), "nuc", sep = "_")
  
  dt_Cell <- as.data.frame(fread(paste0(pth, "/", prefix, cellSuf, ".csv"))) # load cell data
  colnames(dt_Cell) <- paste(colnames(dt_Cell), "cell", sep = "_")
  
  if(file.exists(paste0(pth, "/", prefix, cytoSuf, ".csv"))){  dt_Cyto <- as.data.frame(fread(paste0(pth, "/", prefix, cytoSuf, ".csv"))) # load cytoplasm data
    colnames(dt_Cyto) <- paste(colnames(dt_Cyto), "cyto", sep = "_")
  #browser()
    dt <- cbind(dt_Nuc, dt_Cell[, 31:ncol(dt_Cell)], dt_Cyto[, 31:ncol(dt_Cyto)])
  } else {
    
    dt <- cbind(dt_Nuc, dt_Cell[, 31:ncol(dt_Cell)])
    
  }

  dt <- dt[, !duplicated(colnames(dt))] # Remove duplicated columns
    
  dt_Im <- as.data.frame(fread(paste0(pth, "/", prefix, "Image.csv")))
  dt_Im <- dt_Im[, c("ImageNumber", "FileName_DNA")]
  #browser()
  dt <- dplyr::inner_join(dt_Im, dt, by = c("ImageNumber" = "ImageNumber_nuc"))
  dt
}


## Change the metadata
cpMakeMeta <- function(plate, doses, utRow, df, cellline, cols, treatment, metaColN = 20, suffix = "_nuc", reverse = T){
  
  if (reverse){
    
    doses <- rev(doses)
    
  }
  #browser()
  metaDat <- data.frame(ImageNumber = df$ImageNumber, Filename = df$FileName_DNA, ObjectNumber = df[paste0("ObjectNumber", suffix)],
                        Metadata_Well = df[paste0("Metadata_Well", suffix)], Metadata_Row = df[paste0("Metadata_Row", suffix)],
                        Metadata_Column = df[paste0("Metadata_Column", suffix)],
                        Metadata_Field_nuc = df[paste0("Metadata_Site", suffix)],
                         plate = plate, cell_line = cellline)
  
  metaDat$Treatment <- "X"
  metaDat$Dose <- NA

  rowNms <- unique(df[paste0("Metadata_Row", suffix)],)
  
  ## Get the treatments
  for (cdx in 1:length(cols)){
    
    metaDat$Treatment[metaDat$Metadata_Column_nuc == cols[cdx]] <- treatment[cdx]
    metaDat$Treatment[metaDat$Metadata_Column_nuc == cols[cdx] & metaDat$Metadata_Row_nuc == utRow[cdx]] <- "untreated"
    metaDat$Dose[metaDat$Metadata_Column_nuc == cols[cdx] & metaDat$Metadata_Row_nuc == utRow[cdx]] <- 0
    # Set the treatments
    
    tRowNms <- rowNms[!rowNms == utRow[cdx]]
    for (rdx in 1:length(tRowNms)){
      #browser()
      metaDat$Dose[metaDat$Metadata_Column_nuc == cols[cdx] & metaDat$Metadata_Row_nuc == tRowNms[rdx]] <- doses[[cdx]][rdx]
      
    }
    
  }
 #browser()
  df <- df[, (metaColN + 1):ncol(df)]
  df <- cbind(metaDat, df)
  
return(df)
  
}

## Normalize Data ===================================================================================== ##
# Function to perform Median and MAD normalization on datasets and remove MAD = 0
# Returns a list containing the dataframe with the normalized data and removes features with MAD of 0
# and the remaining features
normalizeCpData <- function(data, features, wells, sepGroups = "Metadata_Well_nuc", plateIder = "Metadata_Plate_nuc", cellType = "CellType"){
  #browser()
  
  medianList <- list()
  madList <- list()
  
  mad0 <- c()
  
  ## Loop through plates
  cnt <- 1
  for (pdx in unique(data[, plateIder])){
    #browser()
    print(paste(cnt, "plate of", length(unique(data[, plateIder])))) 
    cnt <- cnt + 1
    
    plateIdx <- which(data[, plateIder] %in% pdx )
    dataPlt <- data[plateIdx, ]
    #dataPlt <- data[data[ , plateIder] == pdx,]# get the data from the plate
    wellIdx <- which(dataPlt[, sepGroups] %in% wells ) # well index
    
    dataPlt <- dataPlt[wellIdx,] # get the wells to normalize
    
    dataPlt_Median <-  dataPlt[wellIdx, features] %>% summarise_all("median") # get median value
    
    medianList[[pdx]] <- dataPlt_Median
    
    dataPlt_MAD <-  dataPlt[wellIdx, features] %>% summarise_all("mad") # get mad value 
    
    madList[[pdx]] <- dataPlt_MAD
    
    data[plateIdx, features] <- sweep(data[plateIdx, features], 2, t(dataPlt_Median))
    data[plateIdx, features] <- sweep(data[plateIdx, features], 2, t(dataPlt_MAD), FUN = "/")
    
    ## MAD zero list
    #browser()
    mad0 <- unique(c(mad0, names(dataPlt_MAD)[dataPlt_MAD == 0] ))
    
    ## Get medians/MADs of different cell types (if available)
    
    if (length(unique(data[, cellType])) > 1){
      
      for (ctdx in unique(dataPlt[, cellType])){
        #browser()
        dataPltCt <- dataPlt[which(dataPlt[, cellType] %in% ctdx ),] # normalization of plates
        dataPltCt_Median <- dataPltCt[, features] %>% summarise_all("median") # get median va
        medianList[[paste(pdx, ctdx, sep = "_")]] <- dataPltCt_Median
        
        dataPltCt_Mad <- dataPltCt[, features] %>% summarise_all("mad") # get median va
        madList[[paste(pdx, ctdx, sep = "_")]] <- dataPltCt_Mad
        
        
      }
      
    }
    
  }
  
  return(list(data, medianList, madList, mad0))
  
  # Get the row indices of the untreated wells used to calculate the median and the MAD
  # wellIdx <- which(data[, sepGroups] %in% wells )
  # 
  # data_Median <-  data[wellIdx, features] %>% summarise_all("median") # get median value
  # data_MAD <-  data[wellIdx, features] %>% summarise_all("mad") # get mad value
  # 
  # featuresMAD_Rm <- features[which(apply(data_MAD, 2,function(x){any(x == 0)}))] # features with MAD = 0
  # 
  # # Remove features that have MAD0 
  # data <- dplyr::select(data, -featuresMAD_Rm) # remove MAD0 features from the dataset
  # features <- features[!(features %in% featuresMAD_Rm)] # Remove the MAD0 features from the features
  # 
  # data_Median <- dplyr::select(data_Median, -featuresMAD_Rm) # remove MAD0 features from the median
  # data_MAD <- dplyr::select(data_MAD, -featuresMAD_Rm) # finalyy remove MAD0 features from the MAD 
  # 
  # #browser()
  # # Perform normalization
  # dataMed <- cbind(data[, !(colnames(data) %in% features)], sweep(data[, features], 2, t(data_Median)))
  # dataMedMAD <- cbind(data[, !(colnames(data) %in% features)],sweep(dataMed[, features], 2, 1.4826*t(data_MAD), FUN = "/"))
  # 
  # # Return the dataframe with the nromalized data
  # #return(list(dataMedMAD, data_Median, data_MAD))
  # return(list(dataMedMAD, features))
}

## Find Correlated Features ============================================================================ ##

findCorrelatedFeatures <- function(data, features, corThres = 0.9){
  #browser()
  
  st <- Sys.time()
  # Get the correlations between all features
  featureCorr <- WGCNA::cor(data[, features])
  
  ed <- Sys.time()
  
  print(ed - st)
  
  #browser()
  # Create list to store features that are correlated
  corFeatures <- vector(mode = "list", length = length(features))
  names(corFeatures) <- features # make the column names the same as the features
  
  # Search the correlation matrix to determine 
  cnt <- 1
  
  for(dx in colnames(featureCorr)){
    
    corFeatures[[cnt]] <- colnames(featureCorr)[abs(featureCorr[, dx]) > corThres]
    cnt <- cnt + 1
    
  }
  
  # Remove the correlated features, leaving 1 representative
  repLst <- list()
  #browser()
  
  for (dx in names(corFeatures)){
    
    if (!(dx %in% unique(unlist(repLst)))){
      lsFt <- corFeatures[[dx]] # get the correlated features
      lsFt2 <- unique(unlist(corFeatures[lsFt]))
      
      repLst[[dx]] <- lsFt2
    }
  }
  
  #browser()
  return(repLst)
}

# Plot morph feature distributions

plotMorphParamaterDist <- function(dt, param, regions, channels = "AreaShapeFeature", groupby = c("Plate_nuc", "Well_nuc"), ttl = ""){
  
  params <- list() # empty list to store parameter names
  #browser()
  cnt <- 1
  for (cdx in channels){
    for (rdx in regions){
      #browser()
      if (cdx == "AreaShapeFeature"){
        tP <- paste(param, rdx, sep = "_")
      } else{
      tP <- paste(param, cdx, rdx, sep = "_")
      }
      
      #GroupCols <- dt[, paste("Metadata", groupby, sep = "_")]
      GroupCols <- dt[, groupby]
      Group <- sapply(1:nrow(GroupCols), function(x) paste(GroupCols[x, ], collapse =   " "))
      
      tDat <- data.frame(Param = dt[, tP], x = GroupCols[,1], Location = Group) 
      #browser()
      yMed <- median(tDat$Param)
      yMad <- mad(tDat$Param)
      
      yMax <- max(0, min(4*yMad + yMed, max(tDat$Param)))
      yMin <- min(0, max(-4*yMad + yMed, min(tDat$Param)))
      
      params[[cnt]] <- ggplot(tDat, aes(x = Location, y = Param, fill = x)) + 
        geom_violin() + 
        geom_boxplot(width = 0.3) +
        ylab(tP) + theme_bw() + labs(fill = groupby[1]) + 
          theme(axis.text.x = element_text(angle = 90)) + ylim(yMin,yMax) + 
        ggtitle(ttl)
      
      cnt <- cnt + 1
    }
  }
  return(params)
}

## Get Cell Cycle data ================================================================================================== ##
makeCellCycleHists <- function(dat, param = "Intensity_IntegratedIntensity_DNA_Corr_nuc", 
                               paramFilt = "AreaShape_EquivalentDiameter_nuc", groupby = c("Plate_nuc", "Well_nuc")){
  
  params <- list()
  cols <- colnames(dat)
  #browser()
  dat <- dat[, c(cols[grepl("Metadata_", cols) ], param, paramFilt)]
  GroupCols <- dat[, paste("Metadata", groupby, sep = "_"), drop = FALSE]
  dat$Group <- sapply(1:nrow(GroupCols), function(x) paste(GroupCols[x, ], collapse =   " "))
  #browser()
  cnt <- 1
  for (dx in unique(dat$Group)){
    
    params[[dx]] <- dat[dat$Group == dx, ]
    cnt <- cnt + 1
  }
  
  return(params)
  
  
}

## Remove columns ================================================================================================== ##

rmColumns <- function(dat, rmCols){
  
  nms <- paste(rmCols, collapse = "|")
  dat <- dat[, !grepl(nms, colnames(dat))]
  
  return(dat)
  
}

## Get nuclei locations ------====================================================================================== ##

nucLoc <- function(dat, xLoc = "Location_Center_X_nuc", yLoc = "Location_Center_Y_nuc"){
  
  meta <- colnames(dat)[grepl("Metadata", colnames(dat))]
  dat <- dat[ , c("ImageNumber", "ObjectNumber_nuc", "FileName_DNA", meta, xLoc, yLoc)]
  
  return(dat)
}

## Remove na columns =============================================================================================== ##
rmNaCols <- function(dat, featureList){
  
  dat <- dat[, featureList]
  
  dat <- dat[ ,!as.logical(colSums(is.na(dat)))]
  
  return(dat)
  
}

## Remove signle value columns ===================================================================================== ##
rmSingleValueCols <- function(dat, features){
  
  dat <- dat[, features]
  l <- nrow(dat)
  #dat <- dat[, apply(dat, 2, function(x) !(sum(unique(x)) == 1))]
  dat <- dat[, apply(dat[,1:100, (l - 99):l], 2, function(x) !(sum(x)/200 == x[1]))]
  
}

## Remove nearzero variance columns 
rmNearZeroVar <- function(dat, features){
  #browser()
  ## This function takes a long time to run, so ranodmly sample 10000 data points
  dat <- dat[, features]
  #browser()
  dat <- dat[, apply(dat[sample(nrow(dat),10000) , ], 2, function(x) !length(nearZeroVar(x)))]
  
}