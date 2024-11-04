## 
# nucLocs contains the x, y locations of all cells size: (407244, 15)
# WT_Cis_CisplatinDoxUntreatedNrm contains the cell paint information size: (407244, 15)
# WT_Cis_CisplatinDoxUntreatedNrm[!(WT_Cis_CisplatinDoxUntreatedNrm$Treatment == "Doxorubicin"), ] is what was used to train the som
# The image patches are in /data/CEAT/ImageDatasets/CellPaint_WTCis_EpiCompoundDoseResponse/R_Results/SOM_Results_011223_CisOnly_2/all_features/25_node/All_features_25_1/
# We should know the which sub folder they are in as we know the SOM node from the diagram

source("/data/CEAT/ImageDatasets/CellPaint_WTCis_EpiCompoundDoseResponse/RScripts_011122/cellPaintUnsupervisedFunctions.r")
# Load the packages
library(tiff)
library(magick)
library(ggplot2)
library(gridExtra)
library(imager)
source("/data/CEAT/ImageDatasets/CellPaint_WTCis_EpiCompoundDoseResponse/RScripts_011122/cellPaintFunctions.r")
source("/data/CEAT/ImageDatasets/CellPaint_WTCis_EpiCompoundDoseResponse/RScripts_011122/cellPaintUnsupervisedFunctions.r")
library("EBImage")

im_adjust <- function(im, min_px = 0.01, max_px = 0.99){
  
  qt=as.numeric(quantile(im, c(min_px,max_px)))
  
  im[im < qt[1]] <- qt[1]
  im[im > qt[2]] <- qt[2]
  
  im <- renorm(im)
  
  return(im)
  
}

# 
WT_Cis_cell_ims <- dplyr::inner_join(WT_Cis_CisplatinDoxUntreatedNrm[, metColsInd], 
                                     nucLocs, by = c("FilenameDNA" = "FileName_DNA", "ObjectNumber_nuc" = "ObjectNumber_nuc"))


plates = c('P1', 'P5', 'P1', 'P9', 'P5', 'P9')
wells = c('E - 07', 'F - 11', 'H - 02', 'F - 02', 'F - 12', 'E - 01')
flds = c('15', '09', '18', '12', '21', '24')
obs = c(24, 104, 244, 160, 188, 15)
clst = c(5, 3, 8, 8, 3, 3)
cyto_fac_ar <- c(0.995, 0.99, 0.98, 0.99 ,0.88, 1)
dna_fac_ar <- c(0.99, 0.98, 1, 0.98 ,0.9, 1)
#fac

for (dx in 1:6){#:6){

  
  ## Construct the images strings
  dnaI = paste0(plates[dx], '_', wells[dx], '(', 
                'fld', ' ', flds[dx], ' wv UV - DAPI).tif')
  cytoI = paste0(plates[dx], '_', wells[dx], '(', 
                'fld', ' ', flds[dx], ' wv Green - dsRed).tif')
  mitoI = paste0(plates[dx], '_', wells[dx], '(', 
                'fld', ' ', flds[dx], ' wv Red - Cy5).tif')
  
  tDt = WT_Cis_cell_ims[WT_Cis_cell_ims$FilenameDNA == dnaI & WT_Cis_cell_ims$ObjectNumber_nuc == obs[dx],]
  
  ## Load the original image channels and make image patches (look in cell paint unsupervised images)
  dnaPth <- paste0(tDt$Imagepath, dnaI)
  I_DNA <- load.image(dnaPth)

  cytoPth <- paste0(tDt$Imagepath, cytoI)
  I_cyto <- load.image(cytoPth)
  
  mitoPth <- paste0(tDt$Imagepath, mitoI)
  I_mito <- load.image(mitoPth)
  ## Cut out the images
  xYBd = 100
  xMin <- floor(max(1, tDt$Location_Center_X_nuc - xYBd))
  xMax <- ceiling(min(1024, tDt$Location_Center_X_nuc + xYBd))
  yMin <- floor(max(1, tDt$Location_Center_Y_nuc - xYBd))
  yMax <- ceiling(min(1024, tDt$Location_Center_Y_nuc + xYBd))
  
  I_DNA_Patch <- imsub(I_DNA, x > xMin &
                         x < xMax,
                       y > yMin &
                         y < yMax
  )

  I_cyto_Patch <- imsub(I_cyto, x > xMin &
                         x < xMax,
                       y > yMin &
                         y < yMax
  )
  
  I_mito_Patch <- imsub(I_mito, x > xMin &
                         x < xMax,
                       y > yMin &
                         y < yMax
  )
  ## Adjust contrast
  # if (tDt$cell_line == 'WT'){
  #   dna_fac = 0.9
  #   mito_fac = 0.85
  #   cyto_fac = 0.85
  # } else {
  #   dna_fac = 0.98
  #   mito_fac = 0.91
  #   cyto_fac = 0.93
  # }
  #dna_fac <- dna_fac_ar[dx]
  cyto_fac <- cyto_fac_ar[dx]
  dna_fac <- dna_fac_ar[dx]
  mito_fac <- 1
  
  I_DNA_Patch = im_adjust(I_DNA_Patch, max_px = dna_fac)
  I_cyto_Patch = im_adjust(I_cyto_Patch, max_px = cyto_fac)
  I_mito_Patch = im_adjust(I_mito_Patch, max_px = mito_fac)
  
  print(mean(I_cyto_Patch[I_cyto_Patch > 0]))
  
  ## Load the segmented image and get the patch
  I_seg_pth <- paste0(tDt$Segmentationpath, gsub('.tif', paste0('_', tDt$cell_line, '.jpeg'), tDt$FilenameDNA))
  I_seg <- load.image(I_seg_pth)
  I_seg <- imsub(I_seg, x > 1024)
  
  I_seg_ptch <- imsub(I_seg, x > xMin &
                        x < xMax,
                      y > yMin &
                        y < yMax
  )
  ## Make the mask (need the function isolate_cell in cellPaintUnsupervisedFunctions - need to adapt to get just mask)
  I_msk = isolate_cell(I_seg_ptch, x = 100, y = 100, msk_ret = 1)
  
  ## Combine Images
  empty_image <- matrix(0, nrow = dim(I_DNA_Patch)[1], ncol = dim(I_DNA_Patch)[2])
  I1 <- (imappend(list(as.cimg(empty_image), as.cimg(empty_image), as.cimg(I_DNA_Patch*I_msk)), "c"))
  I2 <- (imappend(list(as.cimg(empty_image), as.cimg(I_cyto_Patch*I_msk), as.cimg(empty_image)), "c"))
  I3 <- (imappend(list(as.cimg(I_mito_Patch*I_msk), as.cimg(empty_image), as.cimg(empty_image)), "c"))
  
  #I <- I1 + I2# + I3  
  #I <- I1 + I2 + I3
  I <- I1 + I2 + I3
  plot(I/255.01, rescale = FALSE)
  #hist(I2[I2 > 0])
  
}