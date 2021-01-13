###################################################################################################
###################################################################################################
################ Functions for Kriging Operation for Partial Shapes of a Modelarea ################
###################################################################################################
###################################################################################################
# Autor: David Hoffmann ###########################################################################
# This code is licensed by GNU General Public License v3.0 ########################################
# more information: https://choosealicense.com/licenses/gpl-3.0/ ##################################
###################################################################################################

intersectPointsToModelarea <- function(shape, modShp, datafield){
    modShpst <- st_read (modShp)
    #shape <- point.in.poly(shape, modShpst)
    shape <- st_intersection(shape, modShpst)
    shape <- data.frame(shape)
    shape <- shape[!is.na(shape[datafield]),]
    shape <- st_as_sf(shape)
    if (nrow(shape) != 0) {
      shape <- as_Spatial(shape)
    }
    return(shape)
  } #coord = c("x", "y") #### cuts shape to Modelarea

convertShapeToRaster <- function(modShp, resolution){
  modShpgdal <- readOGR(modShp)
  modShpst <- st_read(modShp)
  Mcoord <- st_crs(modShpst)
  ## Transform the Modelarea to a point grid representing the further kriging raster
  #### see https://stackoverflow.com/questions/14777129/convert-a-shapefile-from-polygons-to-points
  res <- resolution
  BB <- bbox(modShpgdal)
  BB <- res*round(BB/res)
  GT <- GridTopology(cellcentre.offset = BB[,1], 
                     cellsize = c(res, res),
                     cells.dim = (c(diff(BB[1,]), diff(BB[2,]))/res) + 1)
  SP <- SpatialPoints(GT, proj4string = crs(proj4string(modShpgdal)))
  proj4string(modShpgdal) <- proj4string(SP)
  vals <- over(SP, modShpgdal)
  res <- cbind(coordinates(SP), vals)
  ModGrid <- res[!is.na(res$id),]
  ModGrid <- st_as_sf(ModGrid, coords = c("x", "y"), crs = Mcoord)
  ModGrid
} #returns Modelarea grid

bb_grid <- function(shape, cellsize){
  BBGrid <- st_make_grid(st_bbox(shape), cellsize)
}

semivariance <- function(shape, datafield){ #, psill, model, range, kappa
  myVariogram <- variogram(get(datafield)~1, shape) # calculates sample variogram values
  return(myVariogram)
}

fitVariogram <- function(myVariogram,model="Sph"){
  #print(model)
  fittedSemivariance <- fit.variogram(object=myVariogram, model=vgm(model))
  return(fittedSemivariance)
}

krig <- function(shape, datafield, ModGrid, VarFit){
  #Param <- c(log(get(datafield)) ~ 1, shape, ModGrid, model=VarFit)
  #lzn.kriged <<- do.call(krige, Param)
  #CrossVal <<- do.call(krige.cv, Param)
  lzn.kriged <- krige(get(datafield) ~ 1, shape, ModGrid, model=VarFit, nmin= 0,nmax=50)
  CrossVal <- krige.cv(get(datafield) ~ 1, shape, ModGrid, model=VarFit, nmin=0, nmax=50)
  KrigGrid <- lzn.kriged %>%
    mutate(x = unlist(map(lzn.kriged$geometry,1)),
           y = unlist(map(lzn.kriged$geometry,2)))
  coordinates(CrossVal)<-~coords.x1+coords.x2
  result <- list(KrigGrid, CrossVal)
  names(result) <- c("Kriging Grid","Cross Validation")
  return(result)
}
  
plotResiduals <- function(modelarea, CrossVal, colourLow ="yellow", colourHigh="blue"){
  ResMap <- ggplot() +
    geom_sf(data = (st_read(modelarea)), size = 1, color = "black", fill = NA) + 
    coord_sf()+
    geom_point(data = as.data.frame(CrossVal), aes(x = coords.x1, y = coords.x2, colour = residual))+
    scale_color_gradientn(colours = rainbow(5))+
    xlab("x")+ ylab("y")
  return(ResMap)
}

krigingMap <- function(KrigGrid, colourLow ="yellow", colourHigh="red"){
  #ResMap <- bubble(CrossVal, "residual", main = "Residuen")
  Map <-  KrigGrid %>% as.data.frame %>%
    ggplot(aes(x=x,y=y)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
    scale_fill_gradient(low = colourLow, high=colourHigh) +
    theme_bw()
  return(Map)
}


plotVariogram <- function(semivariance, fitDistance) {
  preds = variogramLine(fitDistance, maxdist = max(semivariance$dist))
  varPlot <- ggplot() + 
    geom_point(data = semivariance, aes(x = dist, y = gamma, size=np)) +
    geom_line(data = preds, aes(x = dist, y = gamma))
  return(varPlot)
}

viewSemivariogram <- function(modelareaNumber){allResult[[modelareaNumber]][["Semivariogram"]]}

mapCrossvalResiduals <- function(modelareaNumber){
  allResult[[modelareaNumber]][["CV-Map"]]
}


viewDataquality <- function(viewnumber){
  plotOverview <-ggarrange(viewSemivariogram(viewnumber), mapCrossvalResiduals(viewnumber)+ rremove("x.text"), 
                          labels = c("Semivariogram", "Crossvalidation Map"),
                          ncol = 1, nrow = 2)
  annotate_figure(plotOverview,
                  top = text_grob(paste("Dataquality Kriging of the Partial Area", viewnumber ),face = "bold",size=14))
}

viewMap <- function(data, values, legend, outline, colorPalette = topo.colors(10)){ #terrain.colors(10)
  ggplot() +
    geom_tile(data = data, aes(x = x, y = y, fill = get(values))) +
    scale_fill_gradientn(legend, colours = colorPalette)+
    geom_sf(data = (outline), size = 1, color = "black", fill = NA) + 
    coord_sf()
}

extractXY <-function(Shpsf){Shpsf %>%
    mutate(x = unlist(map(Shpsf$geometry,1)),
           y = unlist(map(Shpsf$geometry,2)))
}