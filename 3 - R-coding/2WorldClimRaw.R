##############################################################
### LOADING PACKAGES AND SETTING UP DIRECTORIES
##############################################################
library("raster")

###############################
####### MAINDIRECTORY:
###############################

# Set the main directory, SDD
mainDir <- "D:/Data"

##############################################################

# SET WORKING DIRECTORY FOR RAW DATA
DataDir.Raw <- paste(mainDir, "/X - Raw.Data", sep="")
setwd(DataDir.Raw)

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES
Dir.Clima <- paste(mainDir, "/1 - Climatology", sep="")  

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES OF NDVI
Dir.Clima.NDVI <- paste(Dir.Clima, "/1 - NDVI", sep="")  

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES OF WORLDLCIM
Dir.Clima.WC <- paste(Dir.Clima, "/2 - WorldClim", sep="")  
dir.create(Dir.Clima.WC)

# SETTING UP A TEMPORARY DIRECTORY
dirtemporary <- paste(mainDir, "/ZZ - Temporary_Storage", sep="")
dir.create(dirtemporary)

###############################################################


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
### DOWNLOADING, CALCULATING AND SAVING WORLDCLIM-CLIMATOLOGIES
##############################################################
WorldClim <- function(){
  # SET WORKING DIRECTORY FOR RAW WORLDCLIM (WC) DATA
  Raw.WC <- paste(DataDir.Raw, "/2 - WorldClim", sep="")
  dir.create(Raw.WC)
  
  ##########################
  ###### PRECIPITATION DATA
  # DOWNLOAD PRECIPITATION CLIMATOLOGY DATA AT 2.5 MINUTE-RESOLUTION
  setwd(Raw.WC)
  prec <- getData('worldclim', var='prec', res=2.5)
  prec.mean <- calc(prec, mean, na.rm=T) # calculating the annual climatology
  
  # RESAMPLING STEP
  setwd(dirtemporary)
  writeRaster(prec.mean, "prec.Climatology.grd", format = "raster", overwrite = TRUE)
  prec.mean.1 <- list.files(path = dirtemporary , pattern = "prec.Climatology.grd")
  prec.mean.1 <- raster(prec.mean.1)
  prec.mean.2 <- resample(prec.mean.1, ndvi, method="ngb")
  
  # PLOT PRECIPITATION CLIMATOLOGIES
  col.prec <- colorRampPalette(c("yellow","blue", "darkblue", "royalblue"))
  plot(prec.mean.2, col= (col.prec(15000)), main = "Mean annual precipitation [mm]", xlab = "Longitude [°]", ylab = "Latitude[°]")  
  
  # SAVE THE PLOTS AND RASTER OF CLIMATOLOGIES TO CLIMATOLOGY DIRECTORY
  setwd(Dir.Clima.WC)
  writeRaster(prec.mean.2, "prec.Climatology.grd", format = "raster", overwrite = TRUE)
  
  jpeg(file=paste(Dir.Clima.WC, "/", "prec.Climatology", ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(prec.mean.2, main="Mean Annual Precipitation [mm]", col=(col.prec(10000)), xlab = "Longitude [°]", ylab = "Latitude[°]",cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black",legend.width=2, axis.args=list(at=seq(0, 900, 100),labels=seq(0, 900, 100),cex.axis=0.9))
  dev.off()
  
  ##########################
  ###### TEMPERATURE DATA
  # DOWNLOAD TEMPERATURE CLIMATOLOGY DATA AT 2.5 MINUTE-RESOLUTION
  setwd(Raw.WC)
  temp <- getData('worldclim', var='tmean', res=2.5)
  temp.kelvin <- temp/10 + 273 # worldlcim temperature data is saved as 10-times °C, this step converts into kelvin
  temp.mean <- calc(temp.kelvin, mean, na.rm=T)
  
  # RESAMPLING STEP
  setwd(dirtemporary)
  writeRaster(temp.mean, "tmean.Climatology.grd", format = "raster", overwrite = TRUE)
  temp.mean.1 <- list.files(path = dirtemporary , pattern = "tmean.Climatology.grd")
  temp.mean.1 <- raster(temp.mean.1)
  temp.mean.2 <- resample(temp.mean.1, ndvi, method="ngb")
  
  # PLOT TEMPERATURE CLIMATOLOGIES
  col.temp <- colorRampPalette(c("darkblue","royalblue","yellow", "red")) # Setting the colour gradient for temp
  plot(temp.mean.1, col= (col.temp(15000)), main = "Mean annual temperature [K]", xlab = "Longitude [°]", ylab = "Latitude[°]")
  
  # SAVE THE PLOTS AND RASTER OF CLIMATOLOGIES TO CLIMATOLOGY DIRECTORY  
  setwd(Dir.Clima.WC)
  writeRaster(temp.mean.2, "tmean.Climatology.grd", format = "raster", overwrite = TRUE)
  
  temp.mean.2[4444:4445] <- 305
  temp.mean.2[4432:4433] <- 245
  jpeg(file=paste(Dir.Clima.WC, "/", "tmean.Climatology", ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(temp.mean.2, main="Mean Annual Temperature [K]", col=(col.temp(10000)), xlab = "Longitude [°]", ylab = "Latitude[°]",cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black",legend.width=2, axis.args=list(at=seq(245, 305, 10),labels=seq(245, 305, 10),cex.axis=0.9))
  dev.off()
  
  # CLEANING THE TEMPORARY DIRECTORY
  setwd(dirtemporary)
  print("Clearing Temporary Directory")
  fileNames.delete <- list.files(path=dirtemporary)
  do.call(file.remove, list(fileNames.delete))
} # end of WorldClim-function


###############################
####### RUN THE FUNCTIONS:
###############################
WorldClim()

print("done")