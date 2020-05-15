##############################################################
### LOADING PACKAGES AND SETTING UP DIRECTORIES
##############################################################
library("raster")
library("rgdal")

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

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES OF WORLDLCIM
Dir.Clima.WC <- paste(Dir.Clima, "/2 - WorldClim", sep="")  

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES OF NDVI
Dir.Clima.NDVI <- paste(Dir.Clima, "/1 - NDVI", sep="")  

# SET DIRECTORY FOR DEM DATA
Dir.DEM <- paste(Dir.Clima, "/3 - DEM", sep="")

# SET UP DIRECTORY FOR SAVING COLINEARITY PLOTS
Dir.PLOTS.COLIN <- paste(Dir.Clima, "/X - Collinearity",sep="")
dir.create(Dir.PLOTS.COLIN)

Dir.Shapes <- paste(mainDir, "/X - Raw.Data/X - ShapeFiles", sep="")

# set up a temporary directory to which the files are moved during the process of the calculations
dirtemporary <- paste(mainDir, "/ZZ - Temporary_Storage", sep="")

###############################
####### PREPARING THE PROGRAM
###############################
setwd(dirtemporary)
print("Clearing temporary directory to avoid remnants in said directory")
fileNames.delete <- list.files(path= dirtemporary)
do.call(file.remove, list(fileNames.delete))

##############################################################
### LOAD DATA, CROP AND ANALYSE FOR COLINEARITY
##############################################################
Colinearity <- function(region, Where, variables, variables.long){
  print(paste("Colinearity for ", paste(variables, collapse = ', '), " in ", Where, sep=""))
  
  #####################################
  # MOVE CLIMATOOGIES INTO PREDEFINED AND SHARED DIRECTORY
  directories <- c(Dir.Clima.NDVI, Dir.Clima.WC, Dir.DEM)
  file.ending <- c(".grd", ".gri")
  
  # MOVE ALL CLIMATOLOGIES TO TEMPOIRARY STORAGE/DIRECTORY    
  for(i in 1:length(directories)){
    for(g in 1:2){
      setwd(directories[i])
      fileNames <- list.files(path = directories[i], pattern = file.ending[g])
      print(paste("copying", fileNames, sep=" ")) # move data of year in question to temporary directory
      file.copy(fileNames, dirtemporary)
    } # end of g-loop (file selection)
  } # end of i-loop (direcotry selection)
  
  #####################################
  # REGION SELECTION WITH SHAPEFILES
  print("Reading ShapeFiles:")
  setwd(Dir.Shapes)
  Urban <- readOGR('.','ne_10m_urban_areas') # reading the shapefiles for urban areas
  Lakes <- readOGR('.','ne_10m_lakes') # reading the shapefiles for lakes and streams
  
  if(region == "Country"){
    Shapes=readOGR('.','ne_50m_admin_0_countries') # reading the shapefiles
  }else{ if(region == "State"){
    Shapes=readOGR('.','ne_10m_admin_1_states_provinces')
  }else{
    print("You have specified no area which to select. No cropping will be done. If you want global data sets ignore this message.")
  }
  }
  
  #####################################
  # IMPORT FURTHER IMPORTANT RASTERS
  setwd(dirtemporary)
  RasterX <- list.files(path = dirtemporary, pattern = "NDVI.Climatology1982-2013.grd")
  RasterX <- raster(RasterX)
  
  #####################################
  # DATA PREPARATION
  # ESTABLISHING AN EMPTY VECTOR TO BE FILLED WITH INDICES OF SHAPEFILES IN POLYGOMFRAME (SHAPES)
  location <- rep(NA, length(Where)) 
  
  # FILLING THE LOCATIONS VECTOR WITH INDICES
  if(region != "Global"){
    for(i in 1:length(Where)){ 
      location[i] <- which(as.vector(Shapes$name) == Where[i])
    }
  }
  
  # EXCEPTIONS FOR RECTANGULAR CROPPING AROUND SHAPES WHERE THE SHAPEFILE ALONE DOESN#T CROP TIGHTLY
  if(Where == "Alaska"){
    area <- extent(-170,-130,52,72)
  }else{ if(Where == c("United States", "Canada", "Mexico")){
    area <- extent(-170,-50,10,90)
    Where <- "NoAm" # set this here for title management later
  }else{ if(Where == "Global"){
    area <- extent(-180,180,-60,90) # this is global data
    Where <- "Global" # set this here for title management later
  }else{
    area <- extent(Shapes[location,]) # this is for simple state selection
  }
  }
  }
  
  # CROPPING IMPORTANT RASTERS
  RasterX <- crop(RasterX,area)
  if(region != "Global"){
    RasterX <- mask(RasterX, Shapes[location,])
  }
  
  Urban <- crop(Urban,area)
  Lakes <- crop(Lakes, area)
  
  # CREATE RASTERS TO MASK URBAN AND LAKE AREAS
  if(is.null(Urban)){
    Antro <- raster(matrix(rep(NA, length(RasterX)), ncol=RasterX@ncols), 
                    xmn=area@xmin, xmx=area@xmax, ymn=area@ymin, ymx=area@ymax)
  }else{
    Antro <- crop(RasterX, area)
    Antro <- mask(Antro, Urban)
    Antro[!is.na(Antro)] <- -8888
  }
  
  if(is.null(Lakes)){
    Water <- raster(matrix(rep(NA, length(RasterX)), ncol=RasterX@ncols), 
                    xmn=area@xmin, xmx=area@xmax, ymn=area@ymin, ymx=area@ymax)
  }else{
    Water <- crop(RasterX, area)
    Water <- mask(Water, Lakes)
    Water[!is.na(Water)] <- -8888
  }
  
  # make the first parameter column, to determine length which is needed to create empty matrix
  setwd(dirtemporary)
  if(variables[1] == "DEM_mn" | variables[1] == "NDVI.Seasonality.Climatology1982-1986" | variables[1] == "NDVI.Climatology1982-1986" | variables[1] == "NDVI.Seasonality.Climatology1982-2013" | variables[1] == "NDVI.Climatology1982-2013"){
    ras <- list.files(path = dirtemporary, pattern = paste(variables[1], sep=""))[1]
  }else{
    ras <- list.files(path = dirtemporary, pattern = paste(variables[1],".Climatology.grd", sep=""))
  }
  ras <- raster(ras)
  ras <- crop(ras,area)
  if(region != "Global"){
    ras <- mask(ras, Shapes[location,])  
  }
  
  # create empty matrix which is to be filled and make first column into row numbers
  Length <- length(as.vector(ras))
  Matrix <- matrix(-8888, nrow = Length, ncol = length(variables)+1)
  Matrix[,1] <- seq(from=1, to=Length, by=1)
  print(paste("Matrix of column length", Length,"created", sep=" "))
  
  # fill the matrix
  for(i in 1:length(variables)){
    if(variables[i] == "DEM_mn" | variables[i] == "NDVI.Seasonality.Climatology1982-1986" | variables[i] == "NDVI.Climatology1982-1986" | variables[i] == "NDVI.Seasonality.Climatology1982-2013" | variables[i] == "NDVI.Climatology1982-2013"){
      ras <- list.files(path = dirtemporary, pattern = paste(variables[i], sep=""))[1]
    }else{
      ras <- list.files(path = dirtemporary, pattern = paste(variables[i],".Climatology.grd", sep=""))
    }
    ras <- raster(ras)
    ras <- crop(ras,area)
    if(region != "Global"){
      ras <- mask(ras, Shapes[location,])
    }
    ras[Antro == -8888] <- NA # masking for urban areas here
    ras[Water == -8888] <- NA # masking for lake areas here
    ras[is.na(RasterX)] <- NA
    Matrix[,i+1] <- as.vector(ras)
    print(paste(variables[i], "masked and fitted into matrix", sep=" "))
    plot(ras, main = variables[i])
  }
  
  # finishing up data preparation
  print("Matrix Ready!")
  colnames(Matrix) <- c("X",variables.long)
  
  # PREPARING THE DATA
  data.set <- Matrix[,-1]
  data.values <- data.set[!rowSums(!is.finite(data.set)),] # remove all rows with non-finite values
  
  # DEFINING COLOUR FOR PLOTTING SYMBOLS
  circles.col=rgb(0,0,0,alpha=0.008)
  
  # DEFINING FUNCTION FOR PRINTING CORRELATION PARAMETERS  
  panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    test <- cor.test(x,y)
    Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.25, paste("r=",txt), cex = cex.cor * r/2.0)
    text(.5, .75, Signif)
  }
  
  # DEFINING COLOUR FOR CORRELATION PLOTS AND SMOOTHING LINES  
  panel.smooth<-function (x, y, col = circles.col, bg = NA, pch = 1, 
                          cex = 0.5, col.smooth = "red", span = 2/3, iter = 3, ...) {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)
  }
  
  # SETTING UP DIRECTORY FOR INDICIDUAL AREAS
  Dir.PLOTS.COLIN.Ind <- paste(Dir.PLOTS.COLIN, "/", Where,sep="")
  dir.create(Dir.PLOTS.COLIN.Ind)
  
  # SAVING THE PLOT
  jpeg(file=paste(Dir.PLOTS.COLIN.Ind, "/", paste(colnames(data.values[]), collapse = '_'), ".jpg", sep = ""), width = 16, height = 11, units = "cm", quality = 100, res = 1000)
  pairs(data.values, upper.panel=panel.smooth,lower.panel=panel.cor, cex.labels=0.6, cex.axis = 0.8)
  dev.off()
  
  # CLEAR TEMPORARY DIRECTORY
  setwd(dirtemporary)
  print("Clearing Temporary Directory")
  fileNames.delete <- list.files(path= dirtemporary)
  do.call(file.remove, list(fileNames.delete))
} # end of Colinearity-function


###############################
####### RUN THE FUNCTION:
###############################
# used in modelling
Colinearity(region = "State", Where = "Alaska", variables <- c("NDVI.Climatology1982-2013", "NDVI.Seasonality.Climatology1982-2013", "tmean", "prec", "DEM_mn"), variables.long <- c("NDVI 1982-2013", "NDVI Seasonality 1982-2013", "Mean Annual Temperature", "Mean Annual Precipitation", "DEM"))
Colinearity(region = "State", Where = "Minnesota", variables <- c("NDVI.Climatology1982-2013", "NDVI.Seasonality.Climatology1982-2013", "tmean", "prec", "DEM_mn"), variables.long <- c("NDVI 1982-2013", "NDVI Seasonality 1982-2013", "Mean Annual Temperature", "Mean Annual Precipitation", "DEM"))
# used in evaluation custering
Colinearity(region = "State", Where = "Alaska", variables <- c("NDVI.Climatology1982-2013", "NDVI.Seasonality.Climatology1982-2013", "tmean.JJA", "prec"), variables.long <- c("NDVI 1982-2013", "NDVI Seasonality 1982-2013", "Mean Summer Temperature", "Mean Annual Precipitation"))
Colinearity(region = "State", Where = "Minnesota", variables <- c("NDVI.Climatology1982-2013", "NDVI.Seasonality.Climatology1982-2013", "tmean.JJA", "prec"), variables.long <- c("NDVI 1982-2013", "NDVI Seasonality 1982-2013", "Mean Summer Temperature", "Mean Annual Precipitation"))

print("Done")