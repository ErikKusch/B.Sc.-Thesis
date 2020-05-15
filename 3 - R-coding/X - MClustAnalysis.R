##############################################################
### LOADING PACKAGES AND SETTING UP DIRECTORIES
##############################################################
library("raster")
library("mclust")
library("fpc")
library("rgl")
library("rgdal")

###############################
####### MAINDIRECTORY:
###############################

# Set the main directory, SDD
mainDir <- "D:/Data"

##############################################################

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES
Dir.Clima <- paste(mainDir, "/1 - Climatology", sep="")  

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES OF WORLDLCIM
Dir.Clima.WC <- paste(Dir.Clima, "/2 - WorldClim", sep="")  

# SET DIRECTORY FOR SAVING THE CLIMATOLOGIES OF EVI
Dir.Clima.NDVI <- paste(Dir.Clima, "/1 - NDVI", sep="")  

# SET DIRECTORY FOR SAVING THE MCLUST-DATA
Dir.Mclust <- paste(mainDir, "/2 - Mclust", sep="")  
dir.create(Dir.Mclust)

# set up a temporary directory to which the files are moved during the process of the calculations
dirtemporary <- paste(mainDir, "/ZZ - Temporary_Storage", sep="")
dir.create(dirtemporary)

Dir.Shapes <- paste(mainDir, "/X - Raw.Data/X - ShapeFiles", sep="")

###############################
####### PREPARING THE PROGRAM
###############################
setwd(dirtemporary)
print("Clearing temporary directory to avoid remnants in said directory")
fileNames.delete <- list.files(path= dirtemporary)
do.call(file.remove, list(fileNames.delete))

##############################################################
### SELECTING DATA, RUN ANALYSIS, SAVE PLOTS
##############################################################
MCLUST <- function(region, Where, clusters, sample.num){
  
  variables <- c("NDVI.Climatology1982-2013", "NDVI.Seasonality.Climatology1982-2013", "prec", "tmean")
  variables.long <- c("NDVI 1982-2913", "NDVI Seasonality 1982-2013", "Precipitation [mm]", "Temperature [K]")
  col.vec <- c(1,1,3,2)
  
  #####################################
  # MOVE CLIMATOOGIES INTO PREDEFINED AND SHARED DIRECTORY
  directories <- c(Dir.Clima.NDVI, Dir.Clima.WC)
  file.ending <- c(".grd", ".gri")
  
  for(i in 1:length(directories)){
    for(g in 1:2){
      setwd(directories[i])
      fileNames <- list.files(path = directories[i], pattern = file.ending[g])
      print(paste("copying Climatologies from", directories[i], sep=" ")) # move data of year in question to temporary directory
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
  
  #####################################
  # INCLUDING THE ACTUAL DATA
  data.set <- Matrix
  rm(Matrix) # save yourself some RAM, the complete raw data is not needed later on
  
  # remove all rows which contain NAs because Mclust can't handle these
  data.values <- data.set[,-1][!rowSums(!is.finite(data.set[,-1])),] 
  
  # create empty vector to be used later on
  empty <- as.vector(rep(NA, length(data.set[,1]))) 
  # figure out which cells correspond to rows without NAs in data-matrix    
  fill.cells <- as.vector(data.set[,1][!is.na(rowSums(data.set[,-1]))]) 
  rm(data.set) # save yourself some RAM
  
  #####################################
  # RANDOMLY SAMPLE DATA
  if(sample.num == "All"){
    sample <- data.values
  }else{
    set.seed(42)
    sample <- data.values[sample(x = 1:nrow(data.values), size= sample.num, replace=FALSE),]
  }  
  
  #####################################
  # MCLUST-ANALYSIS
  G <- clusters
  rm(clusters) 
  
  ### LETTING MCLUST DETERMINE THE CORRECT AMOUNT OF CLUSTERS; MODEL TO BE USED AND MODEL ITSELF
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
  pred <- predict.Mclust(mod, data.values)
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
  done[fill.cells] <- as.vector(pred$classification)
  classes <- done
  classes <- matrix(classes, ncol=n.col, byrow=by.row)
  classes <- raster(classes, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  
  # DATA FOR PROBABILITY PLOT, turn porbabilities as defined by mclust into raster
  done <- empty
  done[fill.cells] <- apply(pred$z, 1, max)
  probs <- done
  probs <- matrix(probs, ncol=n.col, byrow=by.row)
  probs <- raster(probs, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
  
  #####################################
  # SETTING UP PLOTS
  # COLOURING OF CLASSIFICATION PLOT
  my_palette_classes <- colorRampPalette(c("wheat1", "wheat4", "yellow", "yellow4", "lightgoldenrod", "gold", "darkgoldenrod", "darkolivegreen", "darkolivegreen2", "yellowgreen", "mediumspringgreen", "green", "forestgreen", "darkgreen", "indianred1", "indianred4", "red", "red3", "darkred", "orangered3", "chocolate1", "tan1", "lightpink1", "palevioletred1", "hotpink", "violet", "violetred1", "violetred", "magenta", "darkorchid4", "blueviolet", "navyblue"))(n = NClusters)
  col_breaks = c(seq(from=0, to=NClusters, by=1)) # breaks for colouring
  
  # COLOURING OF DIFFERENCE PLOT
  my_palette_differences <- colorRampPalette(c("darkblue", "blue", "cornflowerblue", "maroon3", "yellow", "orange","darkred"))(n=10000)
  
  # COLOURING OF PROBABILITY PLOT
  my_palette_probs <- colorRampPalette(c("darkred","brown2","chocolate1","darkgoldenrod1","darkkhaki", "darkolivegreen1", "chartreuse"))(n=1000)
  
  # COLOURING FOR MCLUST-PLOT
  col.mclust= c("aquamarine","azure4","bisque2","blue","blueviolet","brown1","burlywood4","cadetblue2", "chartreuse","chartreuse4","chocolate1","chocolate4","coral3","cornflowerblue","cyan4", "darkgoldenrod1","darkgoldenrod4","darkgreen","darkred")
  col <- colorRampPalette(col.mclust)
  
  # COLOURS REFERRED TO BY COLOUR VECTOR IN CALL TO FUNCTION
  col.ndvi <- colorRampPalette(c("bisque3","yellow","springgreen", "darkgreen"))(10000)
  col.temp <- colorRampPalette(c("darkblue","royalblue","yellow", "red"))(10000)
  col.prec <- colorRampPalette(c("yellow","blue", "darkblue", "royalblue"))(10000)
  
  # setting up title for all following elements as basis
  title <- paste(Where,"_Classifications_Clusters", NClusters, "_SampleSize", sample.num, sep="")
  
  # setting up new directory where to save the plots
  dirNew1 <- paste(Dir.Mclust, "/",Where, sep="")
  dir.create(dirNew1)
  dirNew <- paste(dirNew1, "/",paste(variables, collapse='_'), sep="")
  dir.create(dirNew)
  setwd(dirNew)
  
  #####################################
  # SAVING THE PLOTS AND MODELS
  print("Saving the Plots")
  # PLOTTING EXPECTED AND ACTUAL VALUES
  for(j in 1:length(variables)){
    
    # colour selection
    if(col.vec[j] == 1){
      colouring <- col.ndvi
      dig <- 2
    }else{ if(col.vec[j] == 2){
      colouring <- col.temp
      dig <- 2
    }else{ if(col.vec[j] == 3){
      colouring <- col.prec
      dig <- 2
    }
    }
    }
    
    # DATA FOR EXPECTED VALUES
    means <- as.vector(mod$parameters$mean[j,])
    done <- empty
    expected.EVI <- rep(NA, length(pred$classification))
    for(i in 1:length(pred$classification)){
      expected.EVI[i] <-  means[pred$classification[i]]
    }
    done[fill.cells] <- expected.EVI
    ras.exEVI <- done
    ras.exEVI <- matrix(ras.exEVI, ncol=n.col, byrow=by.row)
    ras.exEVI <- raster(ras.exEVI, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
    
    # DATA FOR ACTUAL EVI
    done <- empty
    done[fill.cells] <- data.values[,j]
    ras.acEVI <- done
    ras.acEVI <- matrix(ras.acEVI, ncol=n.col, byrow=by.row)
    ras.acEVI <- raster(ras.acEVI, xmn=x1, xmx=x2, ymn=y1, ymx=y2)
    
    # setting up different directorz for each variable
    dirNew2 <- paste(dirNew, "/", variables[j], sep="")
    dir.create(dirNew2)
    
    ####################################
    # DATA FOR DIFFERENCES IN EXPECTED AND ACTUAL EVI
    ras.EVI.difference <- ras.acEVI - ras.exEVI
    
    if(Where == "NoAm"){
      pos2 <- 12
    }else{
      pos2 <- 2
    }
    
    pos1 <- 1
    
    # do this so the colours scale will always be the same on both sides of the spectrum
    if(abs(maxValue(ras.EVI.difference)) > abs(minValue(ras.EVI.difference))){
      ras.EVI.difference[pos1] <- round(maxValue(ras.EVI.difference), digits = dig)
      ras.EVI.difference[pos2] <- -round(maxValue(ras.EVI.difference), digits = dig)
    }else{
      ras.EVI.difference[pos1] <- round(minValue(ras.EVI.difference), digits = dig)
      ras.EVI.difference[pos2] <- -round(minValue(ras.EVI.difference), digits = dig)
    }
    
    max <- maxValue(ras.EVI.difference)
    if(max > 1000){
      r.range <- c(round(minValue(ras.EVI.difference)), round(maxValue(ras.EVI.difference)))
    }else{ if(max > 100){
      r.range <- c(round(minValue(ras.EVI.difference), digits = dig), round(maxValue(ras.EVI.difference), digits =dig))
    }else{
      r.range <- c(round(minValue(ras.EVI.difference), digits = dig), round(maxValue(ras.EVI.difference), digits =dig))
    }
    }
    
    span <- r.range[2]-r.range[1]
    span <- span/10
    
    ras.EVI.difference[1] <- r.range[1]
    ras.EVI.difference[2] <- r.range[2]
    
    jpeg(file=paste(dirNew2, "/", "3 -", title, "_",variables[j],"-Difference", ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(ras.EVI.difference, main = paste("Difference between actual and expected",variables.long[j], "data",sep=" "), col = my_palette_differences, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", axis.args=list(at=seq(r.range[1], r.range[2], span),labels=round(seq(r.range[1], r.range[2], span), digits = dig),cex.axis=0.9), legend.width=2)
    par(new=TRUE)
    plot(Antro, col= "grey", legend=FALSE, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
    dev.off()
    
    ####################################
    # ACTUAL AND EXPECTED DATA PLOT
    if(maxValue(ras.acEVI) > maxValue(ras.exEVI)){
      ras.exEVI[pos1] <- round(maxValue(ras.acEVI), digits = dig)
      ras.acEVI[pos1] <- round(maxValue(ras.acEVI), digits = dig)
    }else{
      ras.acEVI[pos1] <- round(maxValue(ras.exEVI), digits = dig)
      ras.exEVI[pos1] <- round(maxValue(ras.exEVI), digits = dig)
    }
    if(minValue(ras.acEVI) < minValue(ras.exEVI)){
      ras.exEVI[pos2] <- round(minValue(ras.acEVI), digits = dig)
      ras.acEVI[pos2] <- round(minValue(ras.acEVI), digits = dig)
    }else{
      ras.acEVI[pos2] <- round(minValue(ras.exEVI), digits = dig)
      ras.exEVI[pos2] <- round(minValue(ras.exEVI), digits = dig)
    }
    
    max <- maxValue(ras.acEVI)
    if(max > 1000){
      r.range <- c(round(minValue(ras.acEVI)), round(maxValue(ras.acEVI)))
      span <- r.range[2]-r.range[1]
      span <- floor(span/10)
    }else{ if(max > 100){
      r.range <- c(round(minValue(ras.acEVI), digits = dig), round(maxValue(ras.acEVI), digits =dig))
      span <- r.range[2]-r.range[1]
      span <- floor(span)/10
    }else{
      r.range <- c(round(minValue(ras.acEVI), digits = dig), round(maxValue(ras.acEVI), digits =dig))
      span <- r.range[2]-r.range[1]
      span <- span/10
    }
    }
    
    ras.acEVI[1] <- r.range[1]
    ras.acEVI[2] <- r.range[2]
    
    ras.exEVI[1] <- r.range[1]
    ras.exEVI[2] <- r.range[2]
    
    jpeg(file=paste(dirNew2, "/", "2 -", title, "_Actual", variables[j], ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(ras.acEVI, main = paste("Actual",variables.long[j] , "data",sep=" "), col=colouring, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black",axis.args=list(at=seq(r.range[1], r.range[2], span),labels=round(seq(r.range[1], r.range[2], span), digits = dig),cex.axis=0.9), legend.width=2)
    par(new=TRUE)
    plot(Antro, col= "grey", legend=FALSE, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
    dev.off()
    
    jpeg(file=paste(dirNew2, "/", "1 -", title, "_Expected", variables[j], ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    plot(ras.exEVI, main = paste("Expected",variables.long[j], "data",sep=" "), col=colouring, cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", axis.args=list(at=seq(r.range[1], r.range[2], span),labels=round(seq(r.range[1], r.range[2], span), digits = dig),cex.axis=0.9), legend.width=2)
    par(new=TRUE)
    plot(Antro, col= "grey", legend=FALSE, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
    dev.off()
    
    ####################################
    # HISTOGRAMS
    if(sample.num != "All"){
      jpeg(file=paste(dirNew2, "/", "5 -", title, "_",variables[j],"-Histogram_SampleData", ".jpg", sep = ""), width = 64, height = 44, units = "cm", quality = 100, res = 1000)
      hist(sample[,j],main=paste(Where," ", variables.long[j]," data (sampled data set)", sep=""), xlab=variables[j], border="black", col=colouring[10000], probability = TRUE, breaks = 100, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
      dev.off()
    }
    
    jpeg(file=paste(dirNew2, "/", "5 -", title, "_",variables[j],"-Histogram_FullData", ".jpg", sep = ""), width = 64, height = 44, units = "cm", quality = 100, res = 1000)
    hist(data.values[,j],main=paste(Where," ", variables.long[j]," data (full data set)", sep=""), xlab=variables[j], border="black", col=colouring[10000], probability = TRUE, breaks = 100, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
    dev.off()
    
    jpeg(file=paste(dirNew2, "/", "4 -", title, "_",variables[j],"-DifferencesHistogram", ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
    hist(as.vector(ras.EVI.difference),main=paste("Differences in predicted and observed ", variables.long[j]," values", sep=""), xlab="Differences in predicted and observed data", border="black", col="maroon3", probability = TRUE, breaks = 100, cex.lab=1, cex.axis=0.75, cex.main=2, cex.sub=0.5)
    dev.off()
    
  } # end of differences plotting
  
  ####################################
  # CLASSIFICATION PLOT
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
  
  jpeg(file=paste(dirNew, "/", "1 -", title, ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(classes, col = my_palette_classes, main="Classifications", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2, legend=FALSE)
  par(xpd = TRUE)
  legend(legendposition, inset=c(legendpos,0), legend = c(seq(1,NClusters,1), "Urban"), fill = c(my_palette_classes,"grey"), ncol=legendcol)
  par(new=TRUE)
  plot(Antro, col= "grey", legend=FALSE, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
  dev.off()
  
  ####################################
  # CONFIDENCE PLOT
  probs[1] <- 0
  probs[2] <- 1
  jpeg(file=paste(dirNew, "/", "2 -", title, "_Probabilities", ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(probs, col=my_palette_probs, main = "Confidence of assignment", cex.lab=1, cex.axis=1, cex.main=2, cex.sub=1, legend.shrink=1, colNA="black", legend.width=2)
  par(new=TRUE)
  plot(Antro, col= "grey", legend=FALSE, cex.lab=1, cex.axis=1, cex.main=2.5, cex.sub=1)
  dev.off()
  
  ####################################
  # MCLUST UNCERTAINTY PLOT
  jpeg(file=paste(dirNew, "/", "X -", title, "_Mclust-Clusters", ".jpg", sep = ""), width = 32, height = 22, units = "cm", quality = 100, res = 1000)
  plot(mod, what="uncertainty", col=my_palette_classes, main = "Mclust clustering")
  dev.off()
  
  ####################################
  # SAVING RDATA
  output <- cbind(sample, as.vector(mod$classification), as.vector(mod$uncertainty))
  col.names <- colnames(output[,1:length(variables)])
  colnames(output) <- c(col.names, "Cluster", "Uncertainty")
  mean <- mod$parameters$mean
  variance <- mod$parameters$variance
  Name <- mod$modelName
  BIC <- mod$BIC
  
  setwd(dirNew)
  save(mean, file=paste("3 -" ,title, "_ParameterMeans.RData", sep=""))
  save(variance, file=paste("4 -" ,title, "_ParameterVariance.RData", sep=""))
  save(Name, file=paste("5 -" ,title, "_ModelName.RData", sep=""))
  save(BIC, file=paste("6 -" ,title, "_BICValues.RData", sep=""))
  capture.output(summary(mod), file = paste("7 - Model.R", sep=""), append = FALSE, type="output")
}# end of function


###############################
####### RUN THE FUNCTION:
###############################

### Minnesota, provinces
# MCLUST(region <- "State", Where <- "Minnesota", clusters <- 4, sample.num <- "All")

### Minnesota, Full
MCLUST(region <- "State", Where <- "Minnesota", clusters <- 1:20, sample.num <- "All")

### Alaska, according to 10 cluster classification
MCLUST(region <- "State", Where <- "Alaska", clusters <- 10, sample.num <- 20000)

### Alaska, Scheffer
MCLUST(region <- "State", Where <- "Alaska", clusters <- 5, sample.num <- 20000)

### Alaska, Full
MCLUST(region <- "State", Where <- "Alaska", clusters <- 1:20, sample.num <- 20000)

setwd(dirtemporary)
print("Clearing Temporary Directory")
fileNames.delete <- list.files(path= dirtemporary)
do.call(file.remove, list(fileNames.delete))

print("done")
