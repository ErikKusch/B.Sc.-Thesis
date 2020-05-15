##############################################################
### LOADING PACKAGES AND SETTING UP DIRECTORIES
##############################################################
library("raster")
library("mclust")
library("fpc")
library("rgl")
library("rgdal")
library("gtools")
library("mgcv")
library("Epi")
library("xlsx")

###############################
####### MAINDIRECTORY:
###############################

# Set the main directory, SDD
mainDir <- "D:/Data"

##############################################################

# SET WORKING DIRECTORY FOR RAW DATA
DataDir.Raw <- paste(mainDir, "/X - Raw.Data", sep="")

# SET DIRECTORY FOR CHANGE ANALYSIS
Dir.Change.Analysis <- paste(mainDir, "/3 - Change", sep="")  
dir.create(Dir.Change.Analysis)

# SET DIRECTORY FOR SHAPE FILES
Dir.Shapes <- paste(mainDir, "/X - Raw.Data/X - ShapeFiles", sep="")

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES (SET UP HERE FOR FOLLOWING DIRECTORY SPECIFICATION)
Dir.Clima <- paste(mainDir,"/1 - Climatology", sep="")  

# SET DIRECTORY FOR DEM DATA
Dir.DEM <- paste(Dir.Clima, "/3 - DEM", sep="") 

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES OF WORLDLCIM
Dir.Clima.WC <- paste(Dir.Clima, "/2 - WorldClim", sep="")  

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES OF EVI
Dir.Clima.NDVI <- paste(Dir.Clima, "/1 - NDVI", sep="")  

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
### SELECTING DATA, MASKING IT FOR REGION, RUN ANALYSIS
##############################################################
CHANGE.FUN <- function(region, Where, clusters, sample.num, DEM){
  variables <- c("NDVI", "NDVI.Seasonality", "prec", "tmean")
  variables.long <- c("NDVI", "NDVI Seasonality", "Precipitation [mm]", "Temperature [K]")
  
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
  # IMPORT DATA AND FURTHER IMPORTANT RASTERS
  # DEM DATA
  setwd(Dir.DEM)
  files <- list.files(path = Dir.DEM, pattern = paste("DEM_",DEM,".grd",sep=""))
  DEM.ras <- raster(files)
  
  # TEMPORAL AUTOCORRELATION (NEEDED FOR MASKING PROCESSING LATER)
  setwd(Dir.Clima.NDVI)
  RasterX <- list.files(path = Dir.Clima.NDVI, pattern = "NDVI.Climatology1982-2013.grd")
  RasterX <- raster(RasterX)
  
  # CLIMATOLOGIES
  file.ending <- c(".grd", ".gri")
  print("copying data")
  for(i in 1:2){
    setwd(paste(mainDir, "/1 - Climatology/2 - WorldClim",sep=""))
    prec.move <- list.files(path=paste(mainDir, "/1 - Climatology/2 - WorldClim",sep=""), pattern=paste("prec.Climatology", file.ending[i], sep=""))
    file.copy(prec.move, dirtemporary)
    
    tmean.move <- list.files(path=paste(mainDir, "/1 - Climatology/2 - WorldClim",sep=""), pattern=paste("tmean.Climatology", file.ending[i], sep=""))
    file.copy(tmean.move, dirtemporary)
    
    setwd(Dir.Clima.NDVI)
    NDVI.move <- list.files(path=Dir.Clima.NDVI, pattern=paste(file.ending[i], sep=""))
    file.copy(NDVI.move, dirtemporary)
  } # end of moving data
  
  files <- list.files(path = dirtemporary, pattern=".grd")[-1][-2][-2][-3] # data spanning 1982-2013
  files.past <- list.files(path = dirtemporary, pattern=".grd")[-2][-2][-3][-3] # data spanning 1982-1986
  files.present <- list.files(path = dirtemporary, pattern=".grd")[-1][-1][-2][-2] # data spanning 2009-2013
  
  ####################################
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
  
  # CROPPING DEM-DATA
  DEM.ras <- crop(DEM.ras, area)
  if(region != "Global"){
    DEM.ras <- mask(DEM.ras, Shapes[location,])
  }
  
  # CROPPING RASTERX DATA
  RasterX <- crop(RasterX,area)
  if(region != "Global"){
    RasterX <- mask(RasterX, Shapes[location,])
  }
  
  # CROPPING MASKING FILES
  Urban <- crop(Urban,area)
  Lakes <- crop(Lakes, area)
  
  # CREATE RASTERS TO MASK URBAN AND LAKE AREAS
  if(is.null(Urban)){
    Antro <- raster(matrix(rep(NA, length(RasterX)), ncol=RasterX@ncols), xmn=area@xmin, xmx=area@xmax, ymn=area@ymin, ymx=area@ymax)
  }else{
    Antro <- crop(RasterX, area)
    Antro <- mask(Antro, Urban)
    Antro[!is.na(Antro)] <- -8888
  }
  
  if(is.null(Lakes)){
    Water <- raster(matrix(rep(NA, length(RasterX)), ncol=RasterX@ncols), xmn=area@xmin, xmx=area@xmax, ymn=area@ymin, ymx=area@ymax)
  }else{
    Water <- crop(RasterX, area)
    Water <- mask(Water, Lakes)
    Water[!is.na(Water)] <- -8888
  }
  
  # CREATE AN EMPTY MATRIX TO HOLD THE PAST DATA (NDVI; NDVI-SEASONALITY; PREC; TMEAN) AND MAKE FIRST COLUMN INTO ROW NUMBERS
  setwd(dirtemporary)
  ras.ndvi <- files[1]
  ras.ndvi <- raster(ras.ndvi)
  ras.ndvi <- crop(ras.ndvi,area)
  if(region != "Global"){
    ras.ndvi <- mask(ras.ndvi, Shapes[location,])
  }
  
  Length <- length(as.vector(ras.ndvi))
  Matrix <- matrix(-8888, nrow = Length, ncol = length(variables)+2)
  Matrix[,1] <- seq(from=1, to=Length, by=1)
  print(paste("Matrix of column length", Length,"created", sep=" "))
  
  ############################ ENTIRE TIME SERIES DATA
  Matrix.Full <- Matrix
  # FILLING THE MATRIX
  for(i in 1:length(variables)){
    ras <- files[i]
    ras <- raster(ras)
    ras <- crop(ras,area)
    if(region != "Global"){
      ras <- mask(ras, Shapes[location,])
    }
    ras[Antro == -8888] <- NA # masking for urban areas here
    ras[Water == -8888] <- NA # masking for lake areas here
    ras[is.na(ras.ndvi)] <- NA # masking for further lake areas here
    Matrix.Full[,i+1] <- as.vector(ras)
    print(paste(variables[i], "masked and fitted into full matrix", sep=" "))
    plot(ras, main = variables[i])
  }
  Matrix.Full[,length(variables)+2] <- as.vector(DEM.ras)
  
  ############################ PAST DATA
  Matrix.past <- Matrix 
  # FILLING THE MATRIX
  for(i in 1:length(variables)){
    ras <- files.past[i]
    ras <- raster(ras)
    ras <- crop(ras,area)
    if(region != "Global"){
      ras <- mask(ras, Shapes[location,])
    }
    ras[Antro == -8888] <- NA # masking for urban areas here
    ras[Water == -8888] <- NA # masking for lake areas here
    ras[is.na(ras.ndvi)] <- NA # masking for further lake areas here
    Matrix.past[,i+1] <- as.vector(ras)
    print(paste(variables[i], "masked and fitted into past matrix", sep=" "))
    plot(ras, main = variables[i])
  }
  Matrix.past[,length(variables)+2] <- as.vector(DEM.ras)
  
  ############################ PRESENT DATA
  Matrix.present <- Matrix 
  # FILLING THE MATRIX
  for(i in 1:length(variables)){
    ras <- files.present[i]
    ras <- raster(ras)
    ras <- crop(ras,area)
    if(region != "Global"){
      ras <- mask(ras, Shapes[location,])
    }
    ras[Antro == -8888] <- NA # masking for urban areas here
    ras[Water == -8888] <- NA # masking for lake areas here
    ras[is.na(ras.ndvi)] <- NA # masking for further lake areas here
    Matrix.present[,i+1] <- as.vector(ras)
    print(paste(variables[i], "masked and fitted into present matrix", sep=" "))
    plot(ras, main = variables[i])
  }  
  Matrix.present[,length(variables)+2] <- as.vector(DEM.ras)
  
  # FINISHING UP DATA PREPARATION
  colnames(Matrix.Full) <- c("X",variables.long,"DEM")
  colnames(Matrix.past) <- c("X",variables.long,"DEM")
  colnames(Matrix.present) <- c("X",variables.long, "DEM")
  rm(Matrix)
  print("Matrices Ready!")
  
  # SAVE SOME RAM
  rm(RasterX)
  rm(NDVI.move)
  rm(prec.move)
  rm(tmean.move)
  
  #####################################
  # INCLUDING THE ACTUAL DATA
  data.set.Full <- Matrix.Full[,1:3]
  data.set.past <- Matrix.past[,1:3]
  data.set.present <- Matrix.present[,1:3]
  
  # remove all rows which contain NAs because Mclust can't handle these
  data.values.Full <- data.set.Full[,-1][!rowSums(!is.finite(data.set.Full[,-1])),]
  data.values.past <- data.set.past[,-1][!rowSums(!is.finite(data.set.past[,-1])),]
  data.values.present <- data.set.present[,-1][!rowSums(!is.finite(data.set.present[,-1])),]
  
  # needed later on in the script, already processed here to be able to empty data.set since it is taking up huge amounts of RAM
  empty <- as.vector(rep(NA, length(data.set.Full[,1])))
  
  ### figure out which cells correspond to rows without NAs in data-matrix
  fill.cells.Full <- as.vector(data.set.Full[,1][!is.na(rowSums(data.set.Full[,-1]))])
  fill.cells.past <- as.vector(data.set.past[,1][!is.na(rowSums(data.set.past[,-1]))])
  fill.cells.present <- as.vector(data.set.present[,1][!is.na(rowSums(data.set.present[,-1]))])     
  
  # save yourself some RAM
  rm(data.set.Full)
  rm(data.set.past)
  rm(data.set.present) 
  
  #####################################
  # RANDOMLY SAMPLE DATA
  if(sample.num == "All"){
    sample <- data.values.Full
  }else{
    set.seed(42)
    sample <- data.values.Full[sample(x = 1:nrow(data.values.Full), size= sample.num, replace=FALSE),]
  }  
  
  #####################################
  # MCLUST-ANALYSIS
  G <- clusters
  rm(clusters) 
  
  # LETTING MCLUST DETERMINE THE CORRECT AMOUNT OF CLUSTERS; MODEL TO BE USED AND MODEL ITSELF
  print("Calculating BIC for Sample")
  dataBIC <- mclustBIC(sample, G=G)
  print(summary(dataBIC))
  plot(dataBIC)
  
  print("Calculating MODEL according to BIC")
  mod <- mclustModel(sample, dataBIC, G=G)
  mod <-  Mclust(sample, G=G)
  print("Calculated MODEL according to BIC")
  
  NClusters <- mod$G # put optimal number of clusters as defined by mclust into environment for later use
  
  #####################################
  # MCLUST-PREDICTION
  # PREDICTING DATA ASSOCIATION IN LARGE DATA SET USING THE MODEL DEFINED BY SAMPLED DATA
  pred.Full <- predict.Mclust(mod, data.values.Full)
  pred.past <- predict.Mclust(mod, data.values.past)
  pred.present <- predict.Mclust(mod, data.values.present)
  print("Prediction done")
  
  #####################################
  # PREAPRING DATA FOR PLOTTING
  n.col <- ras@ncols
  by.row=TRUE
  y1 <- ras@extent@ymin
  y2 <- ras@extent@ymax
  x1 <- ras@extent@xmin
  x2 <- ras@extent@xmax
  
  # DATA FOR CLASSIFICATION PLOT, turn classification as defined by mclust into raster
  done <- empty
  done[fill.cells.Full] <- as.vector(pred.Full$classification)
  classes <- done
  classes <- matrix(classes, ncol=n.col, byrow=by.row)
  classes <- raster(classes, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  
  done <- empty
  done[fill.cells.past] <- as.vector(pred.past$classification)
  classes.past <- done
  classes.past <- matrix(classes.past, ncol=n.col, byrow=by.row)
  classes.past <- raster(classes.past, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  
  done <- empty
  done[fill.cells.present] <- as.vector(pred.present$classification)
  classes.present <- done
  classes.present <- matrix(classes.present, ncol=n.col, byrow=by.row)
  classes.present <- raster(classes.present, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  
  # DATA FOR PROBABILITY PLOT, turn porbabilities as defined by mclust into raster
  done <- empty
  done[fill.cells.Full] <- apply(pred.Full$z, 1, max)
  probs <- done
  probs <- matrix(probs, ncol=n.col, byrow=by.row)
  probs <- raster(probs, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  
  done <- empty
  done[fill.cells.past] <- apply(pred.past$z, 1, max)
  probs.past <- done
  probs.past <- matrix(probs.past, ncol=n.col, byrow=by.row)
  probs.past <- raster(probs.past, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  
  done <- empty
  done[fill.cells.present] <- apply(pred.present$z, 1, max)
  probs.present <- done
  probs.present <- matrix(probs.present, ncol=n.col, byrow=by.row)
  probs.present <- raster(probs.present, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  
  #####################################
  # SETTING UP PLOTS
  # COLOURING OF CLASSIFICATION PLOT
  my_palette_classes <- colorRampPalette(c("wheat1", "wheat4", "yellow", "yellow4", "lightgoldenrod", "gold", "darkgoldenrod", "darkolivegreen", "darkolivegreen2", "yellowgreen", "mediumspringgreen", "green", "forestgreen", "darkgreen", "indianred1", "indianred4", "red", "red3", "darkred", "orangered3", "chocolate1", "tan1", "lightpink1", "palevioletred1", "hotpink", "violet", "violetred1", "violetred", "magenta", "darkorchid4", "blueviolet", "navyblue"))(n = NClusters)
  col_breaks = c(seq(from=0, to=NClusters, by=1)) # breaks for colouring
  
  # COLOURING OF PROBABILITY PLOT
  my_palette_probs <- colorRampPalette(c("darkred","brown2","chocolate1","darkgoldenrod1","darkkhaki", "darkolivegreen1", "chartreuse"))(n=1000)
  
  # COLOURING FOR ELEVATION PLOT
  col_elevation <- c("grey",colorRampPalette(c("darkgreen", "yellow", "gold3", "darkgoldenrod3", "peru", "chocolate4"))(10000))
  
  # SET WORKING DIRECTORY FOR AREA
  dirNew1 <- paste(Dir.Change.Analysis, "/",Where, sep="")
  dir.create(dirNew1)
  
  # SAVING THE MEANS OF THE CLUSTERS ACCORDING TO THE MCLUST-MODEL
  setwd(dirNew1)
  capture.output(mod$parameters, file = paste("MClustModel_Parameters.R", sep=""), append = FALSE, type="output")
  
  #####################################
  # CLASSIFICATION PLOT
  # SET THE POSITION AND STYLE OF LEGEND CORRESPONDING TO CONTENT OF PLOT
  if(NClusters < 30){
    legendpos <- -0.1
    legendcol <- 1
  }else{
    legendpos <- 0
    legendcol <- 2
  }
  
  if(region != "Global"){
    legendposition <- "topright"
  }else{
    legendposition <- "top"
    legendcol <- NClusters/5
  }
  
  jpeg(file=paste(dirNew1, "/", "1 - Observation_MClust_Full.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(classes, col = my_palette_classes, main="Classifications 1982-2013", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2, legend=FALSE)
  par(xpd = TRUE)
  legend(legendposition, inset=c(legendpos,0), legend = c(seq(1,NClusters,1), "Urban"), fill = c(my_palette_classes,"grey"), ncol=legendcol)
  par(new=TRUE)
  plot(Antro, col= "grey", legend=FALSE, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
  dev.off()
  
  jpeg(file=paste(dirNew1, "/", "3 - Observation_MClust_Past.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(classes.past, col = my_palette_classes, main="Classifications 1982-1986", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2, legend=FALSE)
  par(xpd = TRUE)
  legend(legendposition, inset=c(legendpos,0), legend = c(seq(1,NClusters,1), "Urban"), fill = c(my_palette_classes,"grey"), ncol=legendcol)
  par(new=TRUE)
  plot(Antro, col= "grey", legend=FALSE, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
  dev.off()
  
  jpeg(file=paste(dirNew1, "/", "5 - Observation_MClust_Present.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(classes.present, col = my_palette_classes, main="Classifications 2009-2013", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2, legend=FALSE)
  par(xpd = TRUE)
  legend(legendposition, inset=c(legendpos,0), legend = c(seq(1,NClusters,1), "Urban"), fill = c(my_palette_classes,"grey"), ncol=legendcol)
  par(new=TRUE)
  plot(Antro, col= "grey", legend=FALSE, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
  dev.off()
  
  # SAVING BIOME PROPORTION CHANGES AS EXCEL SHEET
  past <- as.vector(classes.past)
  present <- as.vector(classes.present)
  
  # this matrix will hold the data, rows will show past state, columns will show present state
  changematrix <- matrix(rep(NA, NClusters^2), nrow=NClusters, ncol=NClusters)
  changevec <- rep(NA, NClusters)
  
  for(k in 1:NClusters){
    changerun <- changevec
    changeperc <- changevec
    
    for(m in 1:NClusters){ # fill rates into matrix
      presentcells <- which(present==m) # figure out which cells hold value m
      pastcells <- which(past==k) # figure out which cells hold value k
      
      rate <- length(Reduce(intersect, list(pastcells,presentcells))) # figure out how many of the cell denominators are shared by the two vectors
      changerun[m] <- rate
    }
    changematrix[k,] <- changerun
    
    for(n in 1:NClusters){ # turn rates into percentages
      changeperc[n] <- changematrix[k,n] / sum(changematrix[k,])
    }
    changematrix[k,] <- changeperc
  }
  changematrix <- changematrix*100
  rownames(changematrix) <- seq(1,NClusters, 1)
  colnames(changematrix) <- seq(1,NClusters, 1)
  
  # SETTING UP PERCENTAGE MATRIX FOR STACKED BAR PLOTS OF BIOMES
  percentages <- matrix(rep(NA, NClusters*3), ncol=3, byrow=T)
  percentages[,1] <- prop.table(table(as.vector(classes)))*100
  percentages[,2] <- prop.table(table(as.vector(classes.past)))*100
  percentages[,3] <- prop.table(table(as.vector(classes.present)))*100
  colnames(percentages) <- c("1982-2013", "1982-1986", "2009-2013")
  
  export <- round(cbind(changematrix, rep(NA, NClusters), percentages[,2], percentages[,3]), digits=2)
  colnames(export) <- c(seq(1,NClusters,1), "Empty", "Past Proportions", "Present Porportions")
  
  setwd(dirNew1)
  write.csv(export, "7b - Proportions.csv")
  
  # define plot elements
  legendtext <- c(paste(round(as.vector(percentages[,1]), digits=2), rep("%", NClusters), sep=""), paste(round(as.vector(percentages[,2]), digits=2), rep("%", NClusters), sep= ""), paste(round(as.vector(percentages[,3]), digits=2), rep("%", NClusters), sep= ""))
  
  ypos1 <- rep(NA, NClusters*1)
  perc1 <- as.vector(percentages[,1])
  ypos2 <- rep(NA, NClusters*1)
  perc2 <- as.vector(percentages[,2])
  ypos3 <- rep(NA, NClusters*1)
  perc3 <- as.vector(percentages[,3])
  for(i in 1:NClusters){
    ypos1[i] <- sum(perc1[1:i])-perc1[i]/2
    ypos2[i] <- sum(perc2[1:i])-perc2[i]/2
    ypos3[i] <- sum(perc3[1:i])-perc3[i]/2
  }
  
  jpeg(file=paste(dirNew1, "/", "7 - BiomeProportions.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  barplot <- barplot(percentages, width = c(10, 10, 10), space=0.4, col=my_palette_classes, yaxt='n', main=paste("Biome proportions according to mclust - ", Where, sep=""))
  text(x = c(rep(barplot[1], NClusters), rep(barplot[2], NClusters), rep(barplot[3], NClusters)), y = c(ypos1, ypos2, ypos3),labels = legendtext, cex=0.7)
  par(xpd=NA,mar=par()$mar+c(0,0,0,6))
  legend("topleft", inset=c(0.3,0), legend = rev(seq(1,NClusters,1)), fill = rev(my_palette_classes), ncol=legendcol)
  dev.off()
  
  #####################################
  # CONFIDENCE PLOT
  probs[1] <- 0
  probs[2] <- 1
  
  jpeg(file=paste(dirNew1, "/", "2 - Observation_MClust_Probabilities_Full.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(probs, col=my_palette_probs, main = "Confidence of assignment 1982-2013", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  par(new=TRUE)
  plot(Antro, col= "grey", legend=FALSE, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
  dev.off()
  
  probs.past[1] <- 0
  probs.past[2] <- 1
  
  jpeg(file=paste(dirNew1, "/", "4 - Observation_MClust_Probabilities_Past.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(probs.past, col=my_palette_probs, main = "Confidence of assignment 1982-1986", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  par(new=TRUE)
  plot(Antro, col= "grey", legend=FALSE, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
  dev.off()
  
  probs.present[1] <- 0
  probs.present[2] <- 1
  
  jpeg(file=paste(dirNew1, "/", "6 - Observation_MClust_Probabilities_Present.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(probs.present, col=my_palette_probs, main = "Confidence of assignment 2013", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  par(new=TRUE)
  plot(Antro, col= "grey", legend=FALSE, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
  dev.off()
  
  #####################################
  # CROPPED CLIMATOLOGIES
  Dir.Climats.Cropped <- paste(Dir.Change.Analysis, "/", Where, "/X - Climatologies", sep="")
  dir.create(Dir.Climats.Cropped)
  setwd(Dir.Climats.Cropped)
  
  ### NDVI
  # full
  climats <- matrix(Matrix.Full[,2], ncol=n.col, byrow=by.row)
  climats <- raster(climats, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  col.Climat <- colorRampPalette(c("bisque3","yellow","springgreen", "darkgreen"))(10000)
  climats[1] <- 1
  climats[2] <- 0
  
  jpeg(file=paste(Dir.Climats.Cropped, "/NDVI_1982-2013.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(climats, main = paste (variables.long[1], "1982-2013", sep=" "), col=col.Climat, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  dev.off()
  
  # past
  climats <- matrix(Matrix.past[,2], ncol=n.col, byrow=by.row)
  climats <- raster(climats, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  col.Climat <- colorRampPalette(c("bisque3","yellow","springgreen", "darkgreen"))(10000)
  climats[1] <- 1
  climats[2] <- 0
  
  jpeg(file=paste(Dir.Climats.Cropped, "/NDVI_1982-1986.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(climats, main = paste (variables.long[1], "1982-1986", sep=" "), col=col.Climat, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  dev.off()
  
  # present
  climats <- matrix(Matrix.present[,2], ncol=n.col, byrow=by.row)
  climats <- raster(climats, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  col.Climat <- colorRampPalette(c("bisque3","yellow","springgreen", "darkgreen"))(10000)
  climats[1] <- 1
  climats[2] <- 0
  
  jpeg(file=paste(Dir.Climats.Cropped, "/NDVI_2009-2013.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(climats, main = paste (variables.long[1], "2009-2013", sep=" "), col=col.Climat, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  dev.off()
  
  ### NDVI Seasonality
  # full
  climats <- matrix(Matrix.Full[,3], ncol=n.col, byrow=by.row)
  climats <- raster(climats, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  col.Climat <- colorRampPalette(c("bisque3","yellow","springgreen", "darkgreen"))(10000)
  climats[1] <- 1
  climats[2] <- 0
  
  jpeg(file=paste(Dir.Climats.Cropped, "/NDVI_Seasonality_1982-2013.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(climats, main = paste(variables.long[2], "1982 - 2013", sep=" "), col=col.Climat, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  dev.off()
  
  # past
  climats <- matrix(Matrix.past[,3], ncol=n.col, byrow=by.row)
  climats <- raster(climats, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  col.Climat <- colorRampPalette(c("bisque3","yellow","springgreen", "darkgreen"))(10000)
  climats[1] <- 1
  climats[2] <- 0
  
  jpeg(file=paste(Dir.Climats.Cropped, "/NDVI_Seasonality_1982-1986.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(climats, main = paste (variables.long[2], "1982 - 1986", sep=" "), col=col.Climat, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  dev.off()
  
  # present
  climats <- matrix(Matrix.present[,3], ncol=n.col, byrow=by.row)
  climats <- raster(climats, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  col.Climat <- colorRampPalette(c("bisque3","yellow","springgreen", "darkgreen"))(10000)
  climats[1] <- 1
  climats[2] <- 0
  
  jpeg(file=paste(Dir.Climats.Cropped, "/NDVI_Seasonality_2013.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(climats, main = paste (variables.long[2], "2013", sep=" "), col=col.Climat, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black",
       legend.width=2)
  dev.off()
  
  ### PREC
  climats <- matrix(Matrix.Full[,4], ncol=n.col, byrow=by.row)
  climats <- raster(climats, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  col.Climat <- colorRampPalette(c("yellow","blue", "darkblue", "royalblue"))(10000)
  climats[1] <- 0
  
  jpeg(file=paste(Dir.Climats.Cropped, "/Prec.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(climats, main = paste (variables.long[3], sep=" "), col=col.Climat, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  dev.off()
  
  ### TEMP
  climats <- matrix(Matrix.Full[,5], ncol=n.col, byrow=by.row)
  climats <- raster(climats, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  col.Climat <- colorRampPalette(c("darkblue","royalblue","yellow", "red"))(10000)
  if(Where == "Alaksa"){
    climats[1] <- 250
  }
  
  jpeg(file=paste(Dir.Climats.Cropped, "/Temp.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(climats, main = paste (variables.long[4], sep=" "), col=col.Climat, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  dev.off()
  
  ### DEM
  climats <- matrix(Matrix.Full[,6], ncol=n.col, byrow=by.row)
  climats <- raster(climats, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  col.Climat <- col_elevation
  
  if(minValue(climats)<0){
    jpeg(file=paste(Dir.Climats.Cropped, "/DEM.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(climats, main ="Digital Elevation Model (GMTED2010) [m]",col= col_elevation, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2, breaks=seq(0,maxValue(climats),1), axis.args=list(at=c(minValue(climats),seq(0,maxValue(climats),200),maxValue(climats)), labels=c(minValue(climats),seq(0,maxValue(climats),200),maxValue(climats)),cex.axis=0.9))
    dev.off()
  }else{
    jpeg(file=paste(Dir.Climats.Cropped, "/DEM.jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(climats, main ="Digital Elevation Model (GMTED2010) [m]", col= col_elevation, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
    dev.off()  
  }
  
  #####################################
  # GIMMS TIME SERIES
  print("Loading the NDVI time series")
  Dir.Gimms.Series <- paste(mainDir, "/X - Raw.Data/1 - Gimms-Annual.Mean/Monthly", sep="")
  setwd(Dir.Gimms.Series)
  TS <- list.files(path = Dir.Gimms.Series, pattern = ".grd")
  TS <- mixedsort(TS, decreasing=FALSE) # sort in correct order
  
  TimeSeries <- stack(TS)
  TimeSeries <- crop(TimeSeries, area)
  TimeSeries <- mask(TimeSeries, Shapes[location,])
  
  #####################################
  # CHANGE ANALYSIS
  setwd(Dir.Change.Analysis)
  rm(DEM) # removing DEM to avoid confusion of the term
  
  ROC.Matrix <- matrix(-8888, nrow = length(climats), ncol = NClusters)
  
  for(j in 1:NClusters){
    
    # SET WORKING DIRECTORY FOR AREA
    Dir.Biome.Ind <- paste(dirNew1, "/Biome_",j, sep="")
    dir.create(Dir.Biome.Ind)
    
    binary <- classes
    binary[binary == j] <- 888
    binary[binary != 888] <- 0
    binary[binary == 888] <- 1
    
    col.binary <- colorRampPalette(c("white", my_palette_classes[j]))(2)
    
    jpeg(file=paste(Dir.Biome.Ind, "/", "1 - Distribution_Past_Biome_", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(binary, col = col.binary, main = "Biome Distribution according to mclust in 1982-2013", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend=FALSE)
    par(xpd = TRUE)
    legend("topright", inset=c(0,0), legend = c("Presence", "Absence"), fill = c(my_palette_classes[j],"white"), ncol=legendcol)
    dev.off()
    
    # DATA FRAME FOR DATA ON WHICH MODEL IS TO BE BUILT
    biome.df <- data.frame(biome = as.vector(binary), temp = as.vector(Matrix.Full[,5]), prec = as.vector(Matrix.Full[,4]),DEM = as.vector(Matrix.Full[,6]))
    
    # ESTABLISH MODEL
    print(paste("Caclulating GAM for Biome", j, sep=" "))
    model <- gam(biome ~ s(temp) + s(prec) + s(DEM), data = biome.df, family = binomial(link = logit))
    print("GAM calculated")
    
    setwd(Dir.Biome.Ind)
    capture.output(summary(model), file = paste("Model_Biome",j,".R", sep=""), append = FALSE, type="output")
    
    jpeg(file=paste(Dir.Biome.Ind, "/", "x - ModelStatistics", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    op <- par(mfrow = c(2,2), mar=c(5,4,1,2))
    gam.check(model)
    dev.off()
    
    jpeg(file=paste(Dir.Biome.Ind, "/", "x - ModelSmoothers", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    op <- par(mfrow = c(2,2), mar=c(5,4,1,2))
    plot(model)
    dev.off()
    
    # PREDICTION ON COMPLETE DATA, RETURNS PROBABILITIES OF ASSIGNING PIXEL X TO BIOME Y
    dist <- predict.gam(model, newdata = biome.df, type="response")
    
    dist.curr <- matrix(as.vector(dist), ncol=n.col, byrow=by.row)
    dist.curr <- raster(dist.curr, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
    
    dist.curr[1] <- 1
    dist.curr[2] <- 0
    jpeg(file=paste(Dir.Biome.Ind, "/", "2 - Biome_Probabilites_GAM_", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(dist.curr, col=my_palette_probs, main = "Distribution probabilities according to GAM", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
    dev.off()
    
    #####################################
    # RESPONSE CURVES
    # tmean
    response.df <- biome.df[,2:4]
    response.df[,2] <- rep(mean(response.df[,2], na.rm = TRUE), length(response.df[,2]))
    response.df[,3] <- rep(mean(response.df[,3], na.rm = TRUE), length(response.df[,3]))
    response.df <- response.df[order(response.df[,1]),] # order values to avoid spaghetti plots
    response.df <- na.omit(response.df)
    response.dist <- predict.gam(model, newdata = response.df, type="response")
    
    loess_fit <- loess(response.dist ~ response.df[,1], response.df, span = 0.1)
    loessline <- predict(loess_fit)
    loessline[loessline <0] <- 0
    loessline[loessline >1] <- 1
    
    jpeg(file=paste(Dir.Biome.Ind, "/", "3 - ResponseCurve_Temp_Biome_", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(x = response.df[,1],y = response.dist, main = paste("Response curve of biome assignment to variation in temperature data - Biome", j, sep=" "), xlab = "Temperature [K]", ylab="Probability of correct Biome association", pch=20, type="l", ylim = c(0,1))
    lines(response.df[,1], loessline, col = "green", lwd=3)
    dev.off()
    
    # prec
    response.df <- biome.df[,2:4]
    response.df[,1] <- rep(mean(response.df[,1], na.rm = TRUE), length(response.df[,1]))
    response.df[,3] <- rep(mean(response.df[,3], na.rm = TRUE), length(response.df[,3]))
    response.df <- response.df[order(response.df[,2]),] # order values to avoid spaghetti plots
    response.df <- na.omit(response.df)
    response.dist <- predict.gam(model, newdata = response.df, type="response")
    
    loess_fit <- loess(response.dist ~ response.df[,2], response.df, span = 0.1)
    loessline <- predict(loess_fit)
    loessline[loessline <0] <- 0
    loessline[loessline >1] <- 1
    
    jpeg(file=paste(Dir.Biome.Ind, "/", "5 - ResponseCurve_Prec_Biome_", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(x = response.df[,2],y = response.dist, main = paste("Response curve of biome assignment to variation in precipitation data - Biome", j, sep=" "), xlab = "Precipitation [mm]", ylab="Probability of correct Biome association", pch=20, type="l", ylim = c(0,1))
    lines(response.df[,2], loessline, col = "green", lwd=3)
    dev.off()
    
    # DEM
    response.df <- biome.df[,2:4]
    response.df[,1] <- rep(mean(response.df[,1], na.rm = TRUE), length(response.df[,1]))
    response.df[,2] <- rep(mean(response.df[,2], na.rm = TRUE), length(response.df[,2]))
    response.df <- response.df[order(response.df[,3]),] # order values to avoid spaghetti plots
    response.df <- na.omit(response.df)
    response.dist <- predict.gam(model, newdata = response.df, type="response")
    
    loess_fit <- loess(response.dist ~ response.df[,3], response.df, span = 0.1)
    loessline <- predict(loess_fit)
    loessline[loessline <0] <- 0
    loessline[loessline >1] <- 1
    
    jpeg(file=paste(Dir.Biome.Ind, "/", "7 - ResponseCurve_DEM_Biome_", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(x = response.df[,3],y = response.dist, main = paste("Response curve of biome assignment to variation in elevation data - Biome", j, sep=" "), xlab = "Elevation [m]", ylab="Probability of correct Biome association", pch=20, type="l", ylim = c(0,1))
    lines(response.df[,3], loessline, col = "green", lwd=3)
    dev.off()
    
    #####################################
    # OBSERVED RELATIONS
    # tmean
    jpeg(file=paste(Dir.Biome.Ind, "/", "4 - ObservationCurve_Temp_Biome_", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(x = biome.df[,2],y = as.vector(binary), main = paste("Observed Biome presence in relation to variation in temperature data - Biome", j, sep=" "), xlab = "Temperature [K]", ylab="Biome presence/absence", yaxt='n')
    axis(2, at=c(0,1),labels=c("Absence", "Presence"), col.axis="black", las=3, pch=20, cex=1)
    dev.off()
    # prec
    jpeg(file=paste(Dir.Biome.Ind, "/", "6 - ObservationCurve_Prec_Biome_", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(x = biome.df[,3],y = as.vector(binary), main = paste("Observed Biome presence in relation to variation in precipitation data - Biome", j, sep=" "), xlab = "Precipitation [mm]", ylab="Biome presence/absence", yaxt='n')
    axis(2, at=c(0,1),labels=c("Absence", "Presence"), col.axis="black", las=3, pch=20, cex=1)
    dev.off()
    # DEM
    jpeg(file=paste(Dir.Biome.Ind, "/", "8 - ObservationCurve_DEM_Biome_", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(x = biome.df[,4],y = as.vector(binary), main = paste("Observed Biome presence in relation to variation in elevation data - Biome", j, sep=" "), xlab = "Elevation [m]", ylab="Biome presence/absence", yaxt='n')
    axis(2, at=c(0,1),labels=c("Absence", "Presence"), col.axis="black", las=3, pch=20, cex=1)
    dev.off()
    
    #####################################
    # TIME SERIES
    BiomeTS <- TimeSeries
    BiomeTS[binary == 0] <- NA
    TS.vec <- rep(NA, length(TS))
    for(i in 1: length(TS)){
      TS.vec[i] <- mean(as.vector(BiomeTS[[i]]), na.rm = TRUE)
    }
    
    jpeg(file=paste(Dir.Biome.Ind, "/", "9 - NDVITimeSeries_FULL_Biome", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(TS.vec, col = "green", main = paste("NDVI Time Series - Biome ", j, " (", Where ,")", sep=""), type="line", ylim = c(0,1), lwd = 3, ylab = "NDVI", xlab = "Time since 01/1982 [months]")
    dev.off()
    
    jpeg(file=paste(Dir.Biome.Ind, "/", "10 - NDVITimeSeries_1982-1988_Biome", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(TS.vec[1:84], col = "green", main = paste("NDVI Time Series - Biome ", j, " (", Where ,")", sep=""), type="line", ylim = c(0,1), lwd = 3, ylab = "NDVI", xlab = "Time since 01/1982 [months]")
    abline(v = 12, col = "blue", lwd=2)
    abline(v = 24, col = "blue", lwd=2)
    abline(v = 36, col = "blue", lwd=2)
    abline(v = 48, col = "blue", lwd=2)
    abline(v = 60, col = "blue", lwd=2)
    abline(v = 72, col = "blue", lwd=2)
    dev.off()
    
    #####################################
    # ROC/AUC
    dist.curr[1:2] <- NA
    labels <- factor(as.vector(binary)[!is.na(as.vector(binary))])
    predictions <- as.vector(dist.curr)[!is.na(as.vector(binary))]
    
    jpeg(file=paste(Dir.Biome.Ind, "/", "11 - ROC_Biome", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1) , main = paste("Receiver Operating Characteristic (ROC) - Biome ", j, " (", Where, ")", sep=""))
    par(new = TRUE)
    ROC <- ROC(predictions, labels, plot = "ROC")
    dev.off()
    
    # sensitivity value at maximizing cutoff point (sensitivity + specifcity = MAX)
    opt <- which.max(rowSums(ROC$res[, c("sens", "spec")]))
    # optimal cut-off point
    MaxSens <- ROC$res$predictions[opt]
    
    # re-define my_palette_probs for histograms with colour gradient
    my_palette_probs <- colorRampPalette(c("darkred","brown2","chocolate1","darkgoldenrod1","darkkhaki","darkolivegreen1", "chartreuse"))(n=100)
    dist.curr[1] <- 1
    
    # PLOTTING ASSIGNMENT-PROBABILITIES IN HISTOGRAMM, INCLUDING ROC-THRESHOLD
    jpeg(file=paste(Dir.Biome.Ind, "/", "12 - Histogram_ROC_Biome", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    hist(dist.curr,main=paste("Assignment Confidence Scores - Biome ", j, " (", Where, ")", sep=""),xlab=paste("Confidence", sep=" - "), border="black", col=my_palette_probs, breaks=100)
    abline(v = MaxSens, col = "blue", lwd=2)
    dev.off()
    
    # PLOTTING BIOME DISTRIBUTION WITH ROC-THRESHOLD
    ROCras <- dist.curr
    ROCras[ROCras >= MaxSens] <- 1
    ROCras[ROCras != 1] <- 0
    
    ROCras[1] <- 1
    ROCras[0] <- 0
    
    jpeg(file=paste(Dir.Biome.Ind, "/", "13 - Distribution_ROC_Biome", j,"_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(ROCras, col = col.binary, main = "Biome Distribution according to GAM prediction with ROC threshold", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend=FALSE)
    par(xpd = TRUE)
    legend("topright", inset=c(0,0), legend = c("Presence", "Absence"), fill = c(my_palette_classes[j],"white"), ncol=legendcol)
    dev.off()
    
    ROC.Matrix[,j] <- as.vector(ROCras)
    
  } # loop for each biome
  
  #####################################
  # COMBINED ROC-ANALYSIS
  # build rowsumms of ROC.Matrix and make resulting vector into raster
  ROC.states <- rowSums(ROC.Matrix)
  ROC.states <- matrix(ROC.states, ncol=n.col, byrow=by.row)
  ROC.states <- raster(ROC.states, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  
  jpeg(file=paste(dirNew1, "/", "8 - MultipleStates_ROC_", Where, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(ROC.states, main = "Possible biomes per cell according to GAMs and ROC", col = colorRampPalette(c("grey","cadetblue", "cyan", "yellow", "darkgoldenrod","orange","red"))(NClusters+1), breaks = c(0, seq(0.5,NClusters+0.5,1)), colNA="black", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, legend.width=2, axis.args=list(at=c(0.25, seq(1, NClusters, 1)),labels=c(0, seq(1, NClusters, 1))))
  dev.off()
  
  # MCLUST UNCERTAINTY PLOT
  jpeg(file=paste(dirNew1, "/", "X - Mclust-Clusters", ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(mod, what="uncertainty", col=my_palette_classes, main = "Mclust clustering")
  dev.off()
  
  #####################################
  # CLEARING TEMPORARY DIRECTORY FOR NEXT RUN
  setwd(dirtemporary)
  print("Clearing Temporary Directory")
  fileNames.delete <- list.files(path= dirtemporary)
  do.call(file.remove, list(fileNames.delete))
} # end of function


###############################
####### RUN THE FUNCTION:
###############################

### Minnesota, provinces
CHANGE.FUN(region <- "State", Where <- "Minnesota", clusters <- 4, sample.num <- "All", DEM <-"mn")
### Alaska, Scheffer
CHANGE.FUN(region <- "State", Where <- "Alaska", clusters <- 5, sample.num <- 20000, DEM <-"mn")

print("done")