######################################
# BIAS-CORRECTION MODELING (Max temperature)
######################################
#####################################
# Libraries
library(sf)
library(tidyr)
library(dplyr)
######################################
# Dataset: HUMID + PWS + covariates
###############
HUMID_PWS_NLCD_NED_Pro_idwNOAA_JOINED <- read.csv("H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/Corrected NOAAdata/HUMID_PWS_NLCD_NED_Pro_idwNOAA_JOINED.csv")

########################################
# At some point replace "atl_pws with obs_ga442"
atl_pws <- HUMID_PWS_NLCD_NED_Pro_idwNOAA_JOINED
class(atl_pws)
names(atl_pws)
crs(atl_pws)
class(atl_pws$Date2)
atl_pws$Date2 <- as.Date(atl_pws$Date2)
atl_pws$time_id <- as.numeric(as.factor(as.numeric(atl_pws$Date2)))
atl_pws$spacetime_id <- 1
atl_pws$x <- atl_pws$long_PWS
atl_pws$y <- atl_pws$lat_PWS
# Re-number space_id/time_id
atl_pws$space_id <- as.numeric(as.factor(atl_pws$x)) # to RENUMBER AFTER CLEANING DATA INSTEAD
atl_pws_wgs <- st_as_sf(atl_pws, coords= c("x", "y"), crs=4326)
atl_pws2 <- atl_pws
names(atl_pws)
class(atl_pws)
atl_pws_utm <- st_transform(atl_pws_wgs, crs = 26917)
atl_pws_utm <- atl_pws_utm %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2])

# CHECK for NAs
########################################

atl_pws_utm2 <- atl_pws_utm %>% drop_na(RH)

######################################## REPLACE type of surface and NLCD variables value with indicator variables
# NLCD
# NEW
atl_pws_utm2$NLCD <- as.numeric(as.factor(atl_pws_utm2$NLCD))
#replace NA values with zero in column named col1

atl_pws_utm3 <- atl_pws_utm3 %>% mutate(typeSurf = ifelse(is.na(typeSurf), 7, typeSurf))

####################################### SUBSET the dataset
atl_pws_utm4 <- atl_pws_utm3[c(2,7,19,21,26,27,30:40)]
atl_pws_utm4 <- atl_pws_utm4 %>% drop_na(RH,ceil_hgt,cl,ws,typeSurf)
names(atl_pws_utm4)
class(atl_pws_utm4)

# SAVE ready dataset
write.csv(atl_pws_utm3, "H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/Corrected NOAAdata/atl_pws_utm3.csv")
write.csv(atl_pws_utm4, "H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/Corrected NOAAdata/HUMID_PWS_ensembleDownscaler_ready.csv")

atl_pws_utm4 <- read.csv("H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/Corrected NOAAdata/HUMID_PWS_ensembleDownscaler_ready.csv")

######################################
# Bayesian Downscaler
######################################
install.packages("devtools")
devtools::install_github("WyattGMadden/ensembleDownscaleR")
library(ensembleDownscaleR)
set.seed(42)

###############
### Stage 1 ###
###############

n.iter <- 25e2
burn <- 5e2
thin <- 20

n.iter.pred <- 1e3

pws_fit <- grm(
  Y = atl_pws_utm4$PWS_tmxC,
  X = atl_pws_utm4$tmax,
  L = atl_pws_utm4[, c("NLCD", "NED")],
  M = atl_pws_utm4[, c("RH", "ceil_hgt", "cl", 
                        "ws")],
  n.iter = n.iter,
  burn = burn,
  thin = thin,
  covariance = "matern",
  matern.nu = 0.5,
  coords = atl_pws_utm4[, c("x", "y")],
  space.id = atl_pws_utm4$space_id,
  time.id = atl_pws_utm4$time_id,
  spacetime.id = atl_pws_utm4$spacetime_id,
  verbose.iter = 10
)


saveRDS(pws_fit, "H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/pws_fit.rds")
pws_fit <- readRDS("H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/pws_fit.rds")


###############
### GETTING READY FOR Stage 2 ###
###############

# UPLOAD the HUMID data with idw NOAA

HUMID_idwNOAA_JOINED_NED_NLCD_sf <- st_read("H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/Corrected NOAAdata/HUMID_idwNOAA_JOINED_NED_NLCD.shp")

names(HUMID_idwNOAA_JOINED_NED_NLCD_sf)
length(unique(HUMID_idwNOAA_JOINED_NED_NLCD_sf$geometry))
length(unique(HUMID_idwNOAA_JOINED_NED_NLCD_sf$FID))


# Assign unique code to unique geometry - HUMID TO CREATE space_id variable
HUMID_idwNOAA_JOINED_NED_NLCD_sf <- HUMID_idwNOAA_JOINED_NED_NLCD_sf %>% group_by(geometry) %>% mutate(space_id=cur_group_id())
length(unique(HUMID_idwNOAA_JOINED_NED_NLCD_sf$space_id))

###
atl_pws_wgs <- st_as_sf(atl_pws, coords= c("x", "y"), crs=4326)
###
atl_pws2 <- atl_pws
names(atl_pws)
class(atl_pws)
#atl_pws_utm <- st_as_sf(atl_pws, coords= c("x", "y"), crs = 26917)
atl_pws_utm <- st_transform(atl_pws_wgs, crs = 26917)
atl_pws_utm <- atl_pws_utm %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2])
####
HUMID_for_predictions <- as.data.frame(HUMID_idwNOAA_JOINED_NED_NLCD_sf)
names(HUMID_for_predictions)
HUMID_for_predictions <- HUMID_for_predictions[c(4,11,5,18,19,12,13,9,7,8,14:17,21,23,26)]
names(HUMID_for_predictions)

# RENAME columns
HUMID_for_predictions <- HUMID_for_predictions %>% 
  rename(
    FID = TARGET_F_1,
    x = long_HUMID,
    y = lat_HUMID,
    NLCD = gridcode ,
    NED = grid_code
  )
names(HUMID_for_predictions)
class(HUMID_for_predictions$date)
HUMID_for_predictions$spacetime_id <- 1
HUMID_for_predictions$time_id <- as.numeric(as.factor(as.numeric(HUMID_for_predictions$date)))

######################################## REPLACE NLCD variables value with indicator variables
# NLCD
HUMID_for_predictions$NLCD <- as.numeric(as.factor(HUMID_for_predictions$NLCD))
class(HUMID_for_predictions)
###############
### Stage 2 ###
###############

pws_pred <- grm_pred(
  grm.fit = pws_fit,
  X = HUMID_for_predictions$tmax,
  L = HUMID_for_predictions[, c("NLCD", "NED")],
  M = HUMID_for_predictions[, c("RH", "ceil_hgt", "cl", 
                                "ws")],
  coords = HUMID_for_predictions[, c("x", "y")],
  space.id = HUMID_for_predictions$space_id,
  time.id = HUMID_for_predictions$time_id,
  spacetime.id = HUMID_for_predictions$spacetime_id,
  n.iter = 100, #n.iter.pred,
  verbose = T
)

saveRDS(pws_pred, "H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/pws_pred.rds")

pws_pred <- readRDS("H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/pws_pred.rds")

#################################################### End Stage 2
##############################################

###############
### Stage 3 ###
###############
names(atl_pws_utm4)


cv_id_pws_ord <- create_cv(
  time.id = atl_pws_utm4$time_id, 
  space.id = atl_pws_utm4$space_id,
  spacetime.id = atl_pws_utm4$spacetime_id,
  num.folds = 10L,
  type = "ordinary"
)

saveRDS(cv_id_pws_ord, "H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/cv_id_pws_ord.rds")
cv_id_pws_ord <- readRDS("H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/cv_id_pws_ord.rds")

unique(cv_id_pws_ord$cv.id)
sum(cv_id_pws_ord$cv.id==0)

pws_fit_cv <- grm_cv(
  Y = atl_pws_utm4$PWS_tmxC,
  X = atl_pws_utm4$tmax,
  cv.object = cv_id_pws_ord,
  L = atl_pws_utm4[, c("NLCD", "NED")],
  M = atl_pws_utm4[, c("RH", "ceil_hgt", "cl", 
                        "ws")],
  n.iter = n.iter,
  burn = burn,
  thin = thin,
  coords = atl_pws_utm4[, c("x", "y")],
  space.id = atl_pws_utm4$space_id,
  time.id = atl_pws_utm4$time_id,
  spacetime.id = atl_pws_utm4$spacetime_id,
  verbose.iter = 10
)

# CHECK for NAs
pws_fit_cv[!complete.cases(pws_fit_cv), ]

saveRDS(pws_fit_cv, "H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/pws_fit_cv.rds")
pws_fit_cv <- readRDS("H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/pws_fit_cv.rds")

# REMOVE NAs
########################################
library(tidyr)
names(pws_fit_cv)
pws_fit_cv2 <- pws_fit_cv %>% drop_na(estimate)

# EVALUATE performance metrics
########################################
library(Metrics)
library(caret)

rmse(pws_fit_cv2$obs, pws_fit_cv2$estimate) 
mae(pws_fit_cv2$obs, pws_fit_cv2$estimate) 
R2(pws_fit_cv2$obs, pws_fit_cv2$estimate) 
# PEARSON CORRELATION TEST
cor.test(pws_fit_cv2$obs, pws_fit_cv2$estimate,
         method = "pearson")

mean(pws_fit_cv2$sd)
##################################################################################### End STAGE 3
#####################################################################################


################################################### 

# RANDOM FOREST
################################################### 
###############
### Stage 1 ###
###############

n.iter <- 25e2
burn <- 5e2
thin <- 20

n.iter.pred <- 1e3

rf_pws_fit <- grm(
  Y = rf_pred$pred_rf,
  X = rf_pred$tmax,
  L = rf_pred[, c("NLCD", "NED")],
  M = rf_pred[, c("RH", "ceil_hgt", "cl", 
                        "ws")],
  n.iter = n.iter,
  burn = burn,
  thin = thin,
  covariance = "matern",
  matern.nu = 0.5,
  coords = rf_pred[, c("x", "y")],
  space.id = rf_pred$space_id,
  time.id = rf_pred$time_id,
  spacetime.id = rf_pred$spacetime_id,
  verbose.iter = 10
)


saveRDS(rf_pws_fit, "H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/rf_pws_fit.rds")
rf_pws_fit <- readRDS("H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/pws_fit.rds")

###############
### Stage 2 ###
###############

#cmaq_for_predictions <- readRDS("data/cmaq_for_predictions.rds")


rf_pws_pred <- grm_pred(
  grm.fit = rf_pws_fit,
  X = HUMID_for_predictions$tmax,
  L = HUMID_for_predictions[, c("NLCD", "NED")],
  M = HUMID_for_predictions[, c("RH", "ceil_hgt", "cl", 
                                "ws")],
  coords = HUMID_for_predictions[, c("x", "y")],
  space.id = HUMID_for_predictions$space_id,
  time.id = HUMID_for_predictions$time_id,
  spacetime.id = HUMID_for_predictions$spacetime_id,
  n.iter = 100, #n.iter.pred,
  verbose = T
)

saveRDS(rf_pws_pred, "H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/rf_pws_pred.rds")

pws_pred <- readRDS("H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/pws_pred.rds")

#################################################### End Stage 2
##############################################

###############
### Stage 3 ###
###############
names(atl_pws_utm4)


cv_id_pws_ord <- create_cv(
  time.id = atl_pws_utm4$time_id, 
  space.id = atl_pws_utm4$space_id,
  spacetime.id = atl_pws_utm4$spacetime_id,
  num.folds = 10L,
  type = "ordinary"
)

saveRDS(cv_id_pws_ord, "H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/cv_id_pws_ord.rds")
cv_id_pws_ord <- readRDS("H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/cv_id_pws_ord.rds")

unique(cv_id_pws_ord$cv.id)
sum(cv_id_pws_ord$cv.id==0)

rf_fit_cv <- grm_cv(
  Y = rf_pred$pred_rf,
  X = rf_pred$tmax,
  cv.object = cv_id_pws_ord,
  L = rf_pred[, c("NLCD", "NED")],
  M = rf_pred[, c("RH", "ceil_hgt", "cl", 
                        "ws")],
  n.iter = n.iter,
  burn = burn,
  thin = thin,
  coords = rf_pred[, c("x", "y")],
  space.id = rf_pred$space_id,
  time.id = rf_pred$time_id,
  spacetime.id = rf_pred$spacetime_id,
  verbose.iter = 10
)

# CHECK for NAs
rf_fit_cv[!complete.cases(rf_fit_cv), ]

saveRDS(rf_fit_cv, "H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/rf_fit_cv.rds")# Data with nas (783 obs.)
pws_fit_cv <- readRDS("H:/HUMID/ALL_DATASOURCE/ArcGIS_Pro/ensembleDownscaleR/pws_fit_cv.rds")

# REMOVE NAs
########################################
library(tidyr)
names(rf_fit_cv)
rf_fit_cv2 <- rf_fit_cv %>% drop_na(estimate)
# extracting positions of NA values
print ("Row and Col positions of NA values")
which(is.na(rf_fit_cv2), arr.ind=TRUE)

# EVALUATE performance metrics
########################################
library(Metrics)
library(caret)

rmse(rf_fit_cv2$obs, rf_fit_cv2$estimate) #0.6606373
mae(rf_fit_cv2$obs, rf_fit_cv2$estimate) # 0.4709259
R2(rf_fit_cv2$obs, rf_fit_cv2$estimate) #0.9411253
# PEARSON CORRELATION TEST
cor.test(rf_fit_cv2$obs, rf_fit_cv2$estimate,
         method = "pearson")
# 95 percent confidence interval:
#   0.9694981 0.9707219
# sample estimates:
#   cor 
# 0.9701161
mean(pws_fit_cv2$sd)
##################################################################################### End STAGE 3
#####################################################################################
################################################### END Bayesian Downscaler
