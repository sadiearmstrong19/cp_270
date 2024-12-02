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
  polySelXY = matrix(c(-120.81114, -119.743124, -118.603112, -118.817264, -119.838606, -120.711689, -120.81114, 35.025393, 35.157986, 34.741401, 33.956816, 34.279702, 34.347759, 35.025393), ncol = 2, byrow = FALSE),
 polySelID = 2488)
 
# Create a background extent based on user drawn polygon
bgExt_Fo <- penvs_drawBgExtent(
  polyExtXY = matrix(c(-120.865589, -120.744785, -119.866211, -118.822905, -118.559333, -119.734425, -120.865589, 35.020789, 34.329616, 34.266084, 33.934032, 34.800061, 35.209511, 35.020789),ncol=2,byrow=FALSE), 
  polyExtID = 3362, 
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
                                          
# extract the suitability values for all occurrences
occs_xy_Fo <- occs_Fo[c('longitude', 'latitude')]

# determine the threshold based on the current prediction
occPredVals_Fo <- raster::extract(predSel_Fo, occs_xy_Fo)

# Define probability of quantile based on selected threshold
thresProb_Fo <- switch("p10", 
                              "mtp" = 0, "p10" = 0.1, "qtp" = 0)
                              
# Define threshold value
thres_Fo <- stats::quantile(occPredVals_Fo, probs = thresProb_Fo)

# Applied selected threshold
predSel_Fo <- predSel_Fo > thres_Fo

# Get values of prediction
mapPredVals_Fo <- getRasterVals(predSel_Fo, "cloglog")

# Define colors and legend  
rasCols <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
legendPal <- colorNumeric(rev(rasCols), mapPredVals_Fo, na.color = 'transparent')
rasPal <- c('gray', 'blue')

# Generate map
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap) 
m  %>%
  leaflet::addLegend("bottomright", colors = c('gray', 'blue'),
            title = "Thresholded Suitability<br>(Training)",
            labels = c("predicted absence", "predicted presence"),
            opacity = 1, layerId = "train") %>% 
  #add occurrence data
  addCircleMarkers(data = occs_Fo, lat = ~latitude, lng = ~longitude,
                   radius = 5, color = 'red', fill = TRUE, fillColor = "red",
                   fillOpacity = 0.2, weight = 2, popup = ~pop) %>% 
  ##Add model prediction
  addRasterImage(predSel_Fo, colors = rasPal, opacity = 0.7,
                 group = 'vis', layerId = 'mapPred', method = "ngb") %>%
 ##add background polygons
  addPolygons(data = bgExt_Fo, fill = FALSE,
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

# extract the suitability values for all occurrences
occs_xy_Fo <- occs_Fo[c('longitude', 'latitude')]

# determine the threshold based on the current prediction
occPredVals_Fo <- raster::extract(predSel_Fo, occs_xy_Fo)

# Define probability of quantile based on selected threshold
xfer_thresProb_Fo <- switch("p10", 
                                   "mtp" = 0, "p10" = 0.1, "qtp" = 0)
                                   
# Add threshold if specified 
xfer_time_Fo <- xfer_time_Fo$xferTime > xfer_thresProb_Fo
##Make map
###Make map of transfer
bb_Fo <-  bgExt_Fo@bbox
bbZoom <- polyZoom(bb_Fo[1, 1], bb_Fo[2, 1], bb_Fo[1, 2], 
                   bb_Fo[2, 2], fraction = 0.05)
mapXferVals_Fo <- getRasterVals(xfer_time_Fo,"cloglog")

  # if threshold specified
rasPal_Fo <- c('gray', 'red')
m <- leaflet() %>% addProviderTiles(providers$Esri.WorldTopoMap) 
m %>%
  fitBounds(bbZoom[1], bbZoom[2], bbZoom[3], bbZoom[4]) %>%
  leaflet::addLegend("bottomright", colors = c('gray', 'red'),
            title = "Thresholded Suitability<br>(Transferred)",
            labels = c("predicted absence", "predicted presence"),
            opacity = 1, layerId = 'xfer')%>%
            
# map model prediction raster and polygon of transfer
  clearMarkers() %>% clearShapes() %>% removeImage('xferRas') %>%
  addRasterImage(xfer_time_Fo, colors = rasPal_Fo, opacity = 0.7,
                 layerId = 'xferRas', group = 'xfer', method = "ngb") %>%
 ##add polygon of transfer (same modeling area)
  addPolygons(data = bgExt_Fo, fill = FALSE,
              weight = 4, color = "blue", group = 'xfer')
