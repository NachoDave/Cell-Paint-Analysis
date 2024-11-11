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

## Cluster 1
plates = c('P1', 'P5', 'P1', 'P9', 'P5', 'P9')
wells = c('E - 07', 'F - 11', 'H - 02', 'F - 02', 'F - 12', 'E - 01')
flds = c('15', '09', '18', '12', '21', '24')
obs = c(24, 104, 244, 160, 188, 15)
clst = c(5, 3, 8, 8, 3, 3)
cyto_fac_ar <- c(0.995, 0.99, 0.98, 0.99 ,0.88, 1)
dna_fac_ar <- c(0.99, 0.98, 1, 0.98 ,0.9, 1)
mito_fac_ar <- c()

im_ar1 <- make_image_array(plates, wells, flds, clst, cyto_fac_ar, dna_fac_ar)

## Cluster 2
plates = c('P9', 'P5', 'P5', 'P1', 'P9', 'P9')
wells = c('F - 05', 'F - 11', 'G - 07', 'G - 02', 'A - 02', 'E - 08')
flds = c('20', '13', '22', '05', '07', '23')
obs = c(75, 340, 21, 38, 231, 16)
clst = c(1, 10, 9, 9, 1, 2)
cyto_fac_ar <- c(0.995, 0.99, 0.99, 0.99 ,0.99, 0.98)
dna_fac_ar <- c(0.99, 0.98, 0.98, 0.98 ,0.9, 0.97)
mito_fac_ar <- c()

im_ar2 <- make_image_array(plates, wells, flds, clst, cyto_fac_ar, dna_fac_ar)


## Cluster 3
plates = c('P1', 'P9', 'P1', 'P5', 'P1', 'P9')
wells = c('B - 02', 'E - 02', 'D - 02', 'F - 12', 'E - 04', 'A - 01')
flds = c('08', '14', '20', '04', '13', '14')
obs = c(138, 118, 6, 65, 411, 286)
clst = c(6, 17, 17, 6, 7, 13)
cyto_fac_ar <- c(0.95, 0.99, 0.98, 0.99 ,0.88, 1)
dna_fac_ar <- c(0.99, 0.98, 1, 0.98 ,0.9, 1)

im_ar3 <- make_image_array(plates, wells, flds, clst, cyto_fac_ar, dna_fac_ar)

## Cluster 4
plates = c('P1', 'P1', 'P1', 'P1', 'P1', 'P1')
wells = c('G - 07', 'G - 02', 'F - 11', 'H - 08', 'E - 10', 'E - 09')
flds = c('03', '05', '15', '10', '18', '12')
obs = c(9, 54, 15, 82, 47, 111)
clst = c(14, 15, 20, 19, 24, 20)
cyto_fac_ar <- c(0.995, 0.99, 0.98, 0.99 ,0.88, 0.97)
dna_fac_ar <- c(0.99, 0.98, 0.97, 0.98 ,0.9, 0.93)

im_ar4 <- make_image_array(plates, wells, flds, clst, cyto_fac_ar, dna_fac_ar)

## Cluster 5
plates = c('P1', 'P1', 'P9', 'P1', 'P5', 'P1')
wells = c('E - 09', 'F - 05', 'F - 11', 'A - 08', 'F - 02', 'E - 10')
flds = c('12', '25', '16', '11', '18', '07')
obs = c(89, 22, 56, 35, 5, 71)
clst = c(18, 18, 18, 25, 25, 25)
cyto_fac_ar <- c(0.995, 0.99, 0.98, 0.99 ,0.88, 1)
dna_fac_ar <- c(0.99, 0.98, 1, 0.98 ,0.9, 1)

im_ar5 <- make_image_array(plates, wells, flds, clst, cyto_fac_ar, dna_fac_ar)

#make_image_array <- function(plates, wells, flds, clst, cyto_fac_ar, dna_fac_ar){

  white_image_array <- array(255, dim = c(20, 200, 1, 3))
  white_image <- as.cimg(white_image_array)
  
  clst1_list <- list()
  
  white_image_array2 <- array(255, dim = c(1320, 20, 1, 3))
  white_image2 <- as.cimg(white_image_array2)
  
  II_list <- list(im_ar1, white_image2,
                  im_ar2, white_image2,
                  im_ar3, white_image2,
                  im_ar4, white_image2,
                  im_ar5
                  )
  
  II <-imappend(II_list, "y")
  
#   for (dx in 1:length(plates)){#:6){
#   
#     
#     ## Construct the images strings
#     dnaI = paste0(plates[dx], '_', wells[dx], '(', 
#                   'fld', ' ', flds[dx], ' wv UV - DAPI).tif')
#     cytoI = paste0(plates[dx], '_', wells[dx], '(', 
#                   'fld', ' ', flds[dx], ' wv Green - dsRed).tif')
#     mitoI = paste0(plates[dx], '_', wells[dx], '(', 
#                   'fld', ' ', flds[dx], ' wv Red - Cy5).tif')
#     
#     tDt = WT_Cis_cell_ims[WT_Cis_cell_ims$FilenameDNA == dnaI & WT_Cis_cell_ims$ObjectNumber_nuc == obs[dx],]
#     
#     ## Load the original image channels and make image patches (look in cell paint unsupervised images)
#     dnaPth <- paste0(tDt$Imagepath, dnaI)
#     I_DNA <- load.image(dnaPth)
#   
#     cytoPth <- paste0(tDt$Imagepath, cytoI)
#     I_cyto <- load.image(cytoPth)
#     
#     mitoPth <- paste0(tDt$Imagepath, mitoI)
#     I_mito <- load.image(mitoPth)
#     ## Cut out the images
#     xYBd = 100
#     xMin <- floor(max(1, tDt$Location_Center_X_nuc - xYBd))
#     xMax <- ceiling(min(1024, tDt$Location_Center_X_nuc + xYBd))
#     yMin <- floor(max(1, tDt$Location_Center_Y_nuc - xYBd))
#     yMax <- ceiling(min(1024, tDt$Location_Center_Y_nuc + xYBd))
#     
#     I_DNA_Patch <- imsub(I_DNA, x > xMin &
#                            x < xMax,
#                          y > yMin &
#                            y < yMax
#     )
#   
#     I_cyto_Patch <- imsub(I_cyto, x > xMin &
#                            x < xMax,
#                          y > yMin &
#                            y < yMax
#     )
#     
#     I_mito_Patch <- imsub(I_mito, x > xMin &
#                            x < xMax,
#                          y > yMin &
#                            y < yMax
#     )
#     ## Adjust contrast
#     # if (tDt$cell_line == 'WT'){
#     #   dna_fac = 0.9
#     #   mito_fac = 0.85
#     #   cyto_fac = 0.85
#     # } else {
#     #   dna_fac = 0.98
#     #   mito_fac = 0.91
#     #   cyto_fac = 0.93
#     # }
#     #dna_fac <- dna_fac_ar[dx]
#     cyto_fac <- cyto_fac_ar[dx]
#     dna_fac <- dna_fac_ar[dx]
#     mito_fac <- 1
#     
#     I_DNA_Patch = im_adjust(I_DNA_Patch, max_px = dna_fac)
#     I_cyto_Patch = im_adjust(I_cyto_Patch, max_px = cyto_fac)
#     I_mito_Patch = im_adjust(I_mito_Patch, max_px = mito_fac)
#     
#     print(mean(I_cyto_Patch[I_cyto_Patch > 0]))
#     
#     ## Load the segmented image and get the patch
#     I_seg_pth <- paste0(tDt$Segmentationpath, gsub('.tif', paste0('_', tDt$cell_line, '.jpeg'), tDt$FilenameDNA))
#     I_seg <- load.image(I_seg_pth)
#     I_seg <- imsub(I_seg, x > 1024)
#     
#     I_seg_ptch <- imsub(I_seg, x > xMin &
#                           x < xMax,
#                         y > yMin &
#                           y < yMax
#     )
#     ## Make the mask (need the function isolate_cell in cellPaintUnsupervisedFunctions - need to adapt to get just mask)
#     I_msk = isolate_cell(I_seg_ptch, x = 100, y = 100, msk_ret = 1)
#     
#     ## Combine Images
#     empty_image <- matrix(0, nrow = dim(I_DNA_Patch)[1], ncol = dim(I_DNA_Patch)[2])
#     I1 <- (imappend(list(as.cimg(empty_image), as.cimg(empty_image), as.cimg(I_DNA_Patch*I_msk)), "c"))
#     I2 <- (imappend(list(as.cimg(empty_image), as.cimg(I_cyto_Patch*I_msk), as.cimg(empty_image)), "c"))
#     I3 <- (imappend(list(as.cimg(I_mito_Patch*I_msk), as.cimg(empty_image), as.cimg(empty_image)), "c"))
#     
#     #I <- I1 + I2# + I3  
#     #I <- I1 + I2 + I3
#     I <- I1 + I2 + I3
#     
#     
#     # Get the current dimensions of the image
#     width <- dim(I)[1]
#     height <- dim(I)[2]
#     
#     # Calculate padding for each side
#     pad_x <- (200 - width)
#     pad_y <- (200 - height)
#     
#     if (pad_x){
#       I <- pad(I, nPix = pad_x, axes = 'x', pos = 1)
#     }
#     
#     if (pad_y){
#       I <- pad(I, nPix = pad_y, axes = 'y', pos = 1)
#     }
#     
#     # Ensure padding values are non-negative (only pad smaller images)
#   
#     #plot(I/255.01, rescale = FALSE)
#     #hist(I2[I2 > 0])
#     clst1_list[[2*dx - 1]] <- I
#     clst1_list[[2*dx]] <- white_image
#     #clst1_list <- append(clst1_list, white_image)
#   }
#   
#   IIII <- imappend(clst1_list, "x")
# 
#   plot(IIII)
#   
# #   return(II)
# # }
# 
# xx <- make_image_array(plates, wells, flds, clst, cyto_fac_ar, dna_fac_ar)
