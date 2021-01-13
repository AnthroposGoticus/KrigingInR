###################################################################################################
###################################################################################################
############################## Result Viewer for Krigfing Operation ###############################  
###################################################################################################
###################################################################################################
# Autor: David Hoffmann ###########################################################################
# This code is licensed by GNU General Public License v3.0 ########################################
# more information: https://choosealicense.com/licenses/gpl-3.0/ ##################################
###################################################################################################

library("ggplot2")
library("ggpubr")
library("sf")
library("raster")


modelborder <- "Modellgrenze_Polygon.shp" # put it in "Geodaten/Shapes/Modell
legendname <- "Kupferschieferbasis \n (mNN)" # for the map
crs <- 5678
exportMapName <- "allResult"

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd('..')
if(.Platform$OS.type == "unix") {sep <- "/"} else {sep <- "\\"}
setwd(paste(getwd(),"",sep=sep))
fullModelarea <- st_read(paste("Geodaten","Shapes","Modell",modelborder, sep=sep))
allResult <- readRDS(paste("Metadaten","allResult.rdata",sep=sep))
source(paste("Metadaten","Functions.R", sep=sep))

#Preprocessing
n=length(allResult)
temp = data.frame()
for(i in 1:n){temp = rbind(temp, as.data.frame(allResult[[i]][["Result Kriging"]][["Kriging Grid"]]))}

#Show Map
viewMap(data = temp, values = "var1.pred", legend = legendname, outline = fullModelarea)
  # Colors: https://bootstrappers.umassmed.edu/bootstrappers-courses/pastCourses/rCourse_2016-04/Additional_Resources/Rcolorstyle.html#creating-vectors-of-contiguous-colors

# Export map as geopackage (pointshape)
temp <- st_as_sf(temp, crs=crs)
st_write(temp, paste("Geodaten","Shapes","Kriging",paste(exportMapName,"gpkg",sep="."),sep=sep))


viewDataquality(4) #pick a partial area to view modelquality
    