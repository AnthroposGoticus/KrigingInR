###################################################################################################
###################################################################################################
####################### Kriging Operation for Partial Shapes of a Modelarea ####################### 
###################################################################################################
###################################################################################################
# Autor: David Hoffmann ###########################################################################
# This code is licensed by GNU General Public License v3.0 ########################################
# more information: https://choosealicense.com/licenses/gpl-3.0/ ##################################
###################################################################################################
#install
install.packages("gstat")
install.packages("sp")
install.packages("rgdal")
install.packages("sf")
install.packages("rgeos")
#Packages for Plots
install.packages("ggplot2")
install.packages("ggpubr")
#Parallelisation
install.packages("foreach")
install.packages("doParallel")
#Packages for anything else
install.packages("magrittr")
install.packages("tidyverse")
install.packages("stars")
install.packages("rlist")

#Packages for spacial calculations
library(gstat)
library(sp)
library(rgdal)
library(sf)
library(rgeos)
library(raster)
#Packages for Plots
library(ggplot2)
#Parallelisation
library(foreach)
library(doParallel)
  #Packages for anything else
library(magrittr)
library(tidyverse)
library(stars)
library(rlist)

# USERINPUT
samplepointshape <- "Allpoints.shp" # put it in "Geodaten/Shapes/SamplepointsKriging"
modelborder <- "Modellgrenze_Polygon.shp" # put it in "Geodaten/Shapes/Modell
datafield <- "TiefeNN" #Field for Interpolation in the samplepointshape
res <- 25 #Grid Resolution
Singlefile <- FALSE # If you have only one shapefile with multiple polygons set TRUE, if you have all polygons in multiple shapefiles named as "POLYGON_n set FALSE
  polygonfile <- ("Polygonfile.shp") #only if Singlefile == TRUE put it in "Geodaten/Shapes/Model"


# Workspace
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd('..')
if(.Platform$OS.type == "unix") {sep <- "/"} else {sep <- "\\"}
setwd(paste(getwd(),"",sep=sep))
iniSamplepoints <- st_read(paste("Geodaten","Shapes","SamplepointsKriging",samplepointshape, sep=sep))
fullModelarea <- st_read(paste("Geodaten","Shapes","Modell",modelborder, sep=sep))
subpolygonFolder <- paste("Geodaten","Shapes","Modell","Polygone",sep=sep)
source(paste("Metadaten","Functions.R", sep=sep))

# extract all polygons in shapefile to a individual Shapefile
if(Singlefile == TRUE) {
  fgdb <- st_read(paste("Geodaten","Shapes","Modell",polygonfile, sep=sep))
  dir.create(paste("Geodaten","Shapes","Modell","Polygone", sep=sep), showWarnings = FALSE)
  for (i in c(1:nrow(fgdb))) {
    stnew <- fgdb[i,]
    stnew <- as_Spatial(stnew)
    writeOGR(stnew, paste("Geodaten","Shapes","Modell","Polygone","",sep=sep), layer=paste("Polygon_",i,sep=""),  driver = "ESRI Shapefile" )
  }
}

# Detect computing cores
numCores <- detectCores()
registerDoParallel(numCores)

# Calculates number of shapefiles in the Folder
fils <- list.files(subpolygonFolder, pattern="shp$", full.names = TRUE, recursive = FALSE)
n <- sum(table(fils))

# Defining the function
calcForIndex <- function(index){
  #Calculations
  modelarea <- paste(subpolygonFolder, sep,"Polygon_",index,".shp",sep="")
  if (file.exists(modelarea)){
    print(paste("File", modelarea, "exists"))
    shape <- intersectPointsToModelarea(iniSamplepoints,modelarea, datafield = datafield)
    if(class(shape) != "sf"){
      spoints <- nrow(shape@data)
      if (spoints > 4) { # under 4 Samplepoints no Kriging calculation is possible
        print("File is valid for Kriging")
        modGrid <- convertShapeToRaster(modelarea,res)
        var <- semivariance(shape,datafield)
        fittedSemivariance <- fitVariogram(var)
        krigRes <- krig(shape, datafield, modGrid, fittedSemivariance)
        
        #Plots
        residuals <- plotResiduals(modelarea, krigRes[["Cross Validation"]])
        map <- krigingMap(krigRes[["Kriging Grid"]])
        modelareaSemivariogram <- plotVariogram(var,fittedSemivariance)
        
        #DataStorage
        result <- list(krigRes, fittedSemivariance, modelareaSemivariogram, residuals, map)
        names(result) <- c("Result Kriging","Variogram Parameters", "Semivariogram", "CV-Map","KrigMap")
        return(result)
      } else {return("File is not valid for Kriging and will be excluded")}
    } else {return("ERROR File has zero samplepoints and is not valid for Kriging and will be excluded")}
  } else {return(paste("ERROR File", modelarea, "does not exist"))}
}

# Running the function for each partial shape
allResult <- foreach(index=1:n) %dopar% calcForIndex(index=index) # if you want to see the response for the singulcar functions replace dopar with do

#Save data permanently
list.serialize(allResult, paste("Metadaten","allResult.rdata", sep=sep))

  # to view the data open ViewResults.R