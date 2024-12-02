# cp_270
# Final project for Conservation Planning Practicum (ESM270P)

### Gradient Habitat Suitability ###

library(spocc) 
library(spThin)
library(dismo)
library(sf)
library(ENMeval)
library(wallace)

occs_queryDb

# Query selected database for occurrence records
queryDb_Fo <- occs_queryDb(
  spNames = "Fritillaria ojaiensis", 
  occDb = "gbif", 
  occNum = 1000,
  RmUncertain = FALSE)
occs_Fo <- queryDb_Fo$Fritillaria_ojaiensis$cleaned

# Download environmental data 
envs_Fo <- envs_worldclim(
  bcRes = 2.5, 
  bcSel = c('bio01', 'bio03', 'bio04', 'bio11', 'bio18'), 
  mapCntr = c(-119.826, 34.713), # Mandatory for 30 arcsec resolution   
  doBrick = TRUE)
occs_xy_Fo <- occs_Fo[c('longitude', 'latitude')]
occs_vals_Fo <- as.data.frame(raster::extract(envs_Fo, occs_xy_Fo, cellnumbers = TRUE))

# Remove duplicated same cell values
occs_Fo <- occs_Fo[!duplicated(occs_vals_Fo[, 1]), ]
occs_vals_Fo <- occs_vals_Fo[!duplicated(occs_vals_Fo[, 1]), -1]

# remove occurrence records with NA environmental values
occs_Fo <- occs_Fo[!(rowSums(is.na(occs_vals_Fo)) >= 1), ]

# also remove variable value rows with NA environmental values
occs_vals_Fo <- na.omit(occs_vals_Fo)

# add columns for env variable values for each occurrence record
occs_Fo <- cbind(occs_Fo, occs_vals_Fo)
occs_Fo <- poccs_selectOccs(
 occs = occs_Fo,
 polySelXY = matrix(c(-121.266708, -118.735318, -118.735318, -119.954339, -121.0251, -121.266708, 34.998293, 35.173598, 33.970484, 34.275163, 34.306931, 34.998293), ncol = 2, byrow = FALSE),
 polySelID = 2345)
 
# Create a background extent based on user drawn polygon
bgExt_Fo <- penvs_drawBgExtent(
  polyExtXY = matrix(c(-121.315858, -118.800941, -118.74603, -119.970542, -121.057777, -121.315858, 34.980291, 35.209511, 33.961373, 34.252463, 34.247923, 34.980291),ncol=2,byrow=FALSE), 
  polyExtID = 3441, 
  drawBgBuf = 0, 
  occs = occs_Fo)
  
# Mask environmental data to provided extent
bgMask_Fo <- penvs_bgMask(
  occs = occs_Fo,
  envs = envs_Fo,
  bgExt = bgExt_Fo)
  
# Sample background points from the provided area
bgSample_Fo <- penvs_bgSample(
  occs = occs_Fo,
  bgMask =  bgMask_Fo,
  bgPtsNum = 670)
  
# Extract values of environmental layers for each background point
bgEnvsVals_Fo <- as.data.frame(raster::extract(bgMask_Fo,  bgSample_Fo))
##Add extracted values to background points table
bgEnvsVals_Fo <- cbind(scientific_name = paste0("bg_", "Fritillaria ojaiensis"), bgSample_Fo,
                            occID = NA, year = NA, institution_code = NA, country = NA,
                            state_province = NA, locality = NA, elevation = NA,
                            record_type = NA, bgEnvsVals_Fo)
                            
# R code to get partitioned data
groups_Fo <- part_partitionOccs(
 occs = occs_Fo ,
 bg =  bgSample_Fo, 
 method = "block",
 bgMask = bgMask_Fo,
 aggFact = 2)
 
# Run maxent model for the selected species
model_Fo <- model_maxent(
 occs = occs_Fo,
 bg = bgEnvsVals_Fo,
 user.grp = groups_Fo, 
 bgMsk = bgMask_Fo,
 rms = c(1, 1), 
 rmsStep =  1,
 fcs = 'LQ',
 clampSel = TRUE,
 algMaxent = "maxnet",
 parallel = FALSE,
 numCores = 15)
 
# Select current model and obtain raster prediction
m_Fo <- model_Fo@models[["fc.LQ_rm.1"]]
predSel_Fo <- predictMaxnet(m_Fo, bgMask_Fo,
                                          type = "cloglog", 
                                          clamp = TRUE)
#Get values of prediction
mapPredVals_Fo <- getRasterVals(predSel_Fo, "cloglog")
#Define colors and legend  
rasCols <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
legendPal <- colorNumeric(rev(rasCols), mapPredVals_Fo, na.color = 'transparent')
rasPal <- colorNumeric(rasCols, mapPredVals_Fo, na.color = 'transparent')
#Generate map
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap) 
m  %>%
  leaflet::addLegend("bottomright", pal = legendPal,
            title = "Predicted Suitability<br>(Training)",
            values = mapPredVals_Fo, layerId = "train",
            labFormat = reverseLabel(2, reverse_order = TRUE)) %>% 
  #add occurrence data
  addCircleMarkers(data = occs_Fo, lat = ~latitude, lng = ~longitude,
                   radius = 5, color = 'red', fill = TRUE, fillColor = "red",
                   fillOpacity = 0.2, weight = 2, popup = ~pop) %>% 
  ##Add model prediction
  addRasterImage(predSel_Fo, colors = rasPal, opacity = 0.7,
                 group = 'vis', layerId = 'mapPred', method = "ngb") %>%
 ##add background polygons
  addPolygons(data = bgExt_Fo,fill = FALSE,
              weight = 4, color = "blue", group = 'proj')
              
# Retrieve env variables
n <- mxNonzeroCoefs(model_Fo@models[["fc.LQ_rm.1"]], "maxnet")

# Create response curves
for (i in n) {
maxnet::response.plot(
  model_Fo@models[["fc.LQ_rm.1"]],
  v = i,
  type = "cloglog")
}
Transferring the model to the same modelling area with no threshold
rule. New time based on WorldClim 2.1 variables for 2061-2080 using a
MIROC6 GCM and an SSP of 370.
#Download variables for transferring 
xferTimeEnvs_Fo <- geodata::cmip6_world(
  model = "MIROC6",
  ssp = "370",
  time = "2061-2080",
  var = "bio",
  res = round((raster::res(bgMask_Fo) * 60)[1],1),
  path = tempdir())
names(xferTimeEnvs_Fo) <- paste0('bio', c(paste0('0',1:9), 10:19))

# Select variables for transferring to match variables used for modelling 
xferTimeEnvs_Fo <- xferTimeEnvs_Fo[[names(bgMask_Fo)]]

# Convert to rasterstack
xferTimeEnvs_Fo <- raster::stack(xferTimeEnvs_Fo)

# Generate a transfer of the model to the desired area and time
xfer_time_Fo <-xfer_time(
  evalOut = model_Fo,
  curModel = "fc.LQ_rm.1",
  envs = xferTimeEnvs_Fo,
  xfExt = bgExt_Fo,
  alg = "maxnet",
  outputType = "cloglog",
  clamp = TRUE
  ) 
  
# store the cropped variables of transfer
xferExt_Fo <- xfer_time_Fo$xferExt
###Make map of transfer
bb_Fo <-  bgExt_Fo@bbox
bbZoom <- polyZoom(bb_Fo[1, 1], bb_Fo[2, 1], bb_Fo[1, 2], 
                   bb_Fo[2, 2], fraction = 0.05)
mapXferVals_Fo <- getRasterVals(xfer_time_Fo$xferTime,"cloglog")
rasCols_Fo <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")

# if no threshold specified
legendPal <- colorNumeric(rev(rasCols_Fo), mapXferVals_Fo, na.color = 'transparent')
rasPal_Fo <- colorNumeric(rasCols_Fo, mapXferVals_Fo, na.color = 'transparent')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap) 
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  leaflet::addLegend("bottomright", pal = legendPal,
            title = "Predicted Suitability<br>(Transferred)",
            values = mapXferVals_Fo, layerId = 'xfer',
            labFormat = reverseLabel(2, reverse_order = TRUE)) %>%
            
# map model prediction raster and polygon of transfer
  clearMarkers() %>% clearShapes() %>% removeImage('xferRas') %>%
  addRasterImage(xfer_time_Fo$xferTime, colors = rasPal_Fo, opacity = 0.7,
                 layerId = 'xferRas', group = 'xfer', method = "ngb") %>%
 ##add polygon of transfer (same modeling area)
  addPolygons(data = bgExt_Fo, fill = FALSE,
              weight = 4, color = "blue", group = 'xfer')
