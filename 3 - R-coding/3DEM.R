##############################################################
### LOADING PACKAGES AND SETTING UP DIRECTORIES
##############################################################
library(raster)

###############################
####### MAINDIRECTORY:
###############################

# Set the main directory, SDD
mainDir <- "D:/Data"

##############################################################

# SET WORKING DIRECTORY FOR RAW DATA
DataDir.Raw <- paste(mainDir, "/X - Raw.Data", sep="")
dir.create(DataDir.Raw)
setwd(DataDir.Raw)

# SET DIRECTORY FOR CLIMATOLOGIES
Dir.Clima <- paste(mainDir, "/1 - Climatology", sep="")  

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES OF DEM DATA
Dir.Clima.DEM <- paste(Dir.Clima, "/3 - DEM", sep="")
dir.create(Dir.Clima.DEM)

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES OF NDVI
Dir.Clima.NDVI <- paste(Dir.Clima, "/1 - NDVI", sep="")  

# SETTING UP A TEMPORARY DIRECTORY
dirtemporary <- paste(mainDir, "/ZZ - Temporary_Storage", sep="")

##############################################################

# MOVING NDVI RASTER TO TEMPORARY DIRECTORY FOR RESAMPLING
setwd(Dir.Clima.NDVI)
# moving .gri file
fileNames.ndvi <- list.files(path = Dir.Clima.NDVI, pattern = paste("NDVI.Climatology1982-2013",".gri", sep=""))
print(paste("copying", fileNames.ndvi, sep=" ")) # move data of year in question to temporary directory
file.copy(fileNames.ndvi, dirtemporary)
# moving .grd file
fileNames.ndvi <- list.files(path = Dir.Clima.NDVI, pattern = paste("NDVI.Climatology1982-2013",".grd", sep=""))
print(paste("copying", fileNames.ndvi, sep=" ")) # move data of year in question to temporary directory
file.copy(fileNames.ndvi, dirtemporary)

setwd(dirtemporary)
ndvi <- list.files(path= dirtemporary, pattern = "NDVI.Climatology1982-2013.grd")
ndvi <- raster(ndvi)

##############################################################
### HANDLING OF DEM DATA
##############################################################
DEM <- function(what){
  print(paste("Caclualtion for DEM product: ", what,sep=""))
  
  # SETTING UP COLOUR FOR PLOT
  col_elevation <- c("black",colorRampPalette(c("darkgreen", "yellow", "gold3", "darkgoldenrod3", "peru", "chocolate4"))(10000))
  
  # GRABBING THE RAW DATA
  Raw <- paste(mainDir,"/X - Raw.Data/3 - GMTED2010/30arc/",what,"30_grd/",what,"30_grd",sep="")
  setwd(Raw)
  files <- list.files(path = Raw, pattern = "w001001.adf")
  raster <- raster(files)
  area <- extent(-180, 180, -60, 90)
  raster <- crop(raster, area)
  plot(raster, main ="Digital Elevation Model (GMTED2010) [m]", col= col_elevation, breaks=c(minValue(raster),seq(0,maxValue(raster),1)), cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2, axis.args=list(at=c(minValue(raster),seq(0,maxValue(raster),500),maxValue(raster)), labels=c(minValue(raster),seq(0,maxValue(raster),500),maxValue(raster)), cex.axis=0.9))
  
  # RESAMPLING AND MASKING ACCORDING TO LANDSEA
  raster.res.bilin <- resample(raster, ndvi, method="ngb") # resample to ndvi resolution
  raster.res <- raster.res.bilin
  raster.res[raster.res < 0] <- NA
  
  # SAVING THE RASTER AND PLOT
  setwd(Dir.Clima.DEM)
  writeRaster(raster.res, paste("DEM_",what,".grd",sep=""), format = "raster", overwrite = TRUE)
  
  jpeg(file=paste(Dir.Clima.DEM, "/DEM_",what,".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(raster.res, main ="Digital Elevation Model (GMTED2010) [m]", col= col_elevation, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2, axis.args=list(at=seq(0,maxValue(raster.res),500), labels=seq(0,maxValue(raster.res),500), cex.axis=0.9))
  dev.off()
} # end of DEM-function


###############################
####### RUN THE FUNCTIONS:
###############################
DEM(what <- "mn")

print("done")