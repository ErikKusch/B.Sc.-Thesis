##############################################################
### LOADING PACKAGES AND SETTING UP DIRECTORIES
##############################################################

## !!!! TO RUN THIS CODE PLEASE KEEP IN MIND THAT A CONNECTION TO THE INTERNET HAS TO BE ESTABLISHED THROUGHOUT THE ENTIRE RUNNING TIME !!!!

library(gimms)
library(raster)
library(sp)

###############################
####### MAINDIRECTORY:
###############################

# Set the main directory, SDD
mainDir <- "D:/Data"

##############################################################

# SET WORKING DIRECTORY FOR RAW DATA
DataDir.Raw <- paste(mainDir, "/X - Raw.Data", sep="")
setwd(DataDir.Raw)

# SET WORKING DIRECTORY FOR RAW GIMMS DATA
Dir.Gimms <- paste(DataDir.Raw, "/1 - Gimms", sep="")
dir.create(Dir.Gimms)
setwd(Dir.Gimms)

# UPDATE FILE LIST IN GIMMS DIRECTORY
gimms_files <- updateInventory

# SET DIRECTORY FOR RASTERS OF ANNUAL NDVI-COMPOSITES
Dir.Gimms.Annual <- paste(DataDir.Raw, "/1 - Gimms-Annual.Mean", sep="")
dir.create(Dir.Gimms.Annual)
setwd(Dir.Gimms.Annual)

# SET DIRECTORY FOR RASTERS OF NDVI SEASONALITY COMPOSITES
Dir.Gimms.Seasonality <- paste(DataDir.Raw, "/1 - Gimms-Seasonality", sep="")
dir.create(Dir.Gimms.Seasonality)
setwd(Dir.Gimms.Seasonality)

##############################################################
### SELECTING DATA, EXPORTING PLOTS AND RASTERS
##############################################################
RasterGIMMs <- function(from, to){
  for(year in from:to){
    # PREPARING DATA
    gimms_files <- downloadGimms(x = as.Date(paste(year,"-01-01",sep="")), y = as.Date(paste(year,"-12-31",sep="")), dsn = Dir.Gimms)
    # Rasterize files
    gimms_raster <- rasterizeGimms(x = gimms_files, remove_header = TRUE)
    
    # Create monthly maximum value composites
    indices <- monthlyIndices(gimms_files)
    gimms_raster_mvc <- monthlyComposite(gimms_raster, indices = indices)
    
    gimms_annual <- calc(gimms_raster_mvc, fun=mean, progress='text')
    gimms_annual <- crop(gimms_annual, extent(-180,180,-60,90))
    
    gimms_annual[gimms_annual<0] <- 0 # set threshold for barren land (NDVI<0)
    gimms_annual[gimms_annual>1] <- 1 # set threshold for saturated NDVI (NDVI > 1)
    
    # SETTING UP PLOTS
    col.evi <- colorRampPalette(c("bisque3","yellow","springgreen", "darkgreen"))(10000) # Setting the colour gradient for ndvi
    
    # SAVING DATA
    print(paste("Saving raster of annual means for year", year, sep=" "))
    setwd(Dir.Gimms.Annual)
    writeRaster(gimms_annual, paste("NDVI_",year, ".grd",sep=""), format = "raster", overwrite = TRUE)
    
    print(paste("Saving rasters of monthly means for year", year, sep=" "))
    dir.create(paste(Dir.Gimms.Annual, "/Monthly",sep=""))
    setwd(paste(Dir.Gimms.Annual, "/Monthly",sep=""))
    gimms_raster_mvc <- crop(gimms_raster_mvc, extent(-180,180,-60,90))
    gimms_raster_mvc[gimms_raster_mvc<0] <- 0 # set threshold for barren land (NDVI<0)
    gimms_raster_mvc[gimms_raster_mvc>1] <- 1 # set threshold for saturated NDVI (NDVI > 1)
    for(k in 1:12){
      month <- gimms_raster_mvc[[k]]
      writeRaster(month, paste("NDVI_",year,"_",k, ".grd",sep=""), format = "raster", overwrite = TRUE)
      plot(month, main = paste("Monthly data ", year, " month ", k,sep=""))
    }
    
    plot(gimms_raster_mvc[[1]])
    
    # PLOTTING AND SAVING PLOTS
    # MONTHLY PLOTS
    gimms_raster_mvc[gimms_raster_mvc<0] <-0
    gimms_raster_mvc[gimms_raster_mvc>1] <-1
    
    gimms_raster_mvc[4331] <- 1
    gimms_raster_mvc[4330] <- 0
    
    names(gimms_raster_mvc) <- paste(month.abb, year)
    
    print(paste("Saving plots of annual means for year", year, sep=" "))
    jpeg(file=paste(Dir.Gimms.Annual, "/", "Months_NDVI_",year,".jpg", sep = ""), width = 30, height = 30, units = "cm", quality = 100, res = 1000)
    print(spplot(gimms_raster_mvc, main = paste("Normalized Difference Vegetation Index (NDVI)", year, sep=" "),col.regions=col.evi))
    dev.off()
    
    # ANNUAL PLOTS
    gimms_annual[4330] <- 1
    gimms_annual[4331] <- 0
    
    jpeg(file=paste(Dir.Gimms.Annual, "/", "NDVI_",year, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(gimms_annual, main = paste ("Normalized Difference Vegetation Index (NDVI)", year, sep=" "), col=col.evi, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2, axis.args=list(at=seq(0, 1, 0.1),labels=seq(0, 1, 0.1),cex.axis=0.9))
    dev.off()
  } # end of year-loop
}# end of RasterGIMMS-function


##############################################################
### SELECTING DATA, EXPORTING PLOTS AND RASTERS
##############################################################
GIMMsSeasons <- function(from, to){
  for(year in from:to){
    # PREPARING DATA
    gimms_files <- downloadGimms(x = as.Date(paste(year,"-01-01",sep="")), y = as.Date(paste(year,"-12-31",sep="")), dsn = Dir.Gimms)
    # Rasterize files
    gimms_raster <- rasterizeGimms(x = gimms_files, remove_header = TRUE)
    
    # Create monthly maximum value composites
    indices <- monthlyIndices(gimms_files)
    gimms_raster_mvc <- monthlyComposite(gimms_raster, indices = indices)
    
    maxi <- calc(gimms_raster_mvc, fun=max, progress = 'text')
    mini <- calc(gimms_raster_mvc, fun=min, progress = 'text')
    
    gimms_seasonality <- maxi-mini
    gimms_seasonality[gimms_seasonality>1]<-1
    gimms_seasonality[gimms_seasonality<0]<-0
    gimms_seasonality <- crop(gimms_seasonality, extent(-180,180,-60,90))
    
    # SETTING UP PLOTS
    col.evi <- colorRampPalette(c("bisque3","yellow","springgreen", "darkgreen"))(10000) # Setting the colour gradient for ndvi
    
    # SAVING DATA
    print(paste("Saving raster of ndvi seasonality for year", year, sep=" "))
    setwd(Dir.Gimms.Seasonality)
    writeRaster(gimms_seasonality, paste("NDVI_",year, "_Seasonality.grd",sep=""), format = "raster", overwrite = TRUE)
    
    # PLOTTING AND SAVING PLOTS
    # ANNUAL PLOTS
    gimms_seasonality[4330] <- 1
    gimms_seasonality[4331] <- 0
    # gimms_annual[gimms_annual<0] <- NA # set threshold for barren land (NDVI<0)
    
    print(paste("Saving plots of ndvi seasonality for year", year, sep=" "))
    jpeg(file=paste(Dir.Gimms.Seasonality, "/", "NDVI_",year, "_Seasonality.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(gimms_seasonality, main = paste ("Normalized Difference Vegetation Index (NDVI) Seasonality", year, sep=" "), col=col.evi, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2, axis.args=list(at=seq(0, 1, 0.1),labels=seq(0, 1, 0.1),cex.axis=0.9))
    dev.off()
  } # end of year-loop
}# end of GIMMsSeasons-function


##############################################################
### SAVING NDVI BASED CLIMATOLOGIES
##############################################################
GIMMS.Climat <- function(parameter, timespan){
  year.range <- 1982:2013
  
  title <- paste(parameter,".Climatology", min(timespan),"-", max(timespan), sep="")
  
  col.plot <- colorRampPalette(c("bisque3","yellow","springgreen", "darkgreen"))(10000) # Setting the colour gradient for evi
  
  if(parameter == "NDVI"){
    parameter.long <- "Normalized Difference Vegetation Index (NDVI)"
    Dir.GIMMS.NEW <- paste(DataDir.Raw, "/1 - Gimms-Annual.Mean", sep="")  
  }else{
    Dir.GIMMS.NEW <- paste(DataDir.Raw, "/1 - Gimms-Seasonality", sep="")  
    parameter.long <- "NDVI Seasonality"
  }
  
  setwd(Dir.GIMMS.NEW)
  # select the actual data rasters needed
  fileNames <- list.files(path = Dir.GIMMS.NEW, pattern = ".grd")[match(timespan, year.range)]
  
  # CALCULATING THE CLIMATOLOGY
  Stack <- stack(fileNames)
  Climatology <- calc(Stack, mean, na.rm=T)
  
  Dir.Change <- paste(mainDir,"/1 - Climatology",sep="")
  dir.create(Dir.Change)
  Dir.Change.Climat <- paste(Dir.Change, "/1 - NDVI",sep="")
  dir.create(Dir.Change.Climat)
  
  # SAVING THE RASTERS
  setwd(Dir.Change.Climat)
  writeRaster(Climatology, paste(title,".grd",sep=""), format = "raster", overwrite = TRUE)
  
  Climatology[4330] <- 1
  jpeg(file=paste(Dir.Change.Climat, "/", title,".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(Climatology, col = col.plot, main=  paste(parameter.long," Climatology ", min(timespan)," - ", max(timespan), sep=""), cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  dev.off()
} # end of GIMMS.climat-function


###############################
####### RUN THE FUNCTIONS:
###############################
RasterGIMMs(from=1982, to=2013)
GIMMsSeasons(from=1982, to=2013)

GIMMS.Climat(parameter = "NDVI", timespan = 1982:1986)
GIMMS.Climat(parameter = "NDVI.Seasonality", timespan = 1982:1986)
GIMMS.Climat(parameter = "NDVI", timespan = 2009:2013)
GIMMS.Climat(parameter = "NDVI.Seasonality", timespan = 2009:2013)
GIMMS.Climat(parameter = "NDVI", timespan = 1982:2013)
GIMMS.Climat(parameter = "NDVI.Seasonality", timespan = 1982:2013)

print("done")