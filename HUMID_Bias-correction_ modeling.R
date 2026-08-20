######################################### MODELING 1: Bayesian
########################################
#####################################################################################
# Upload libraries
##############################################

library(randomForest)
library(mlbench)
library(caret)
library(dplyr)
library(tidyr)
library(sf)
library(Metrics)


# Upload the dataset

pws_obs_utm <- readRDS("/DATA/Modeling2016-2018/pws_obs_utm.rds")
###################################### APPLY the Bayesian model
install.packages("devtools")
devtools::install_github("WyattGMadden/ensembleDownscaleR")
library(ensembleDownscaleR)
set.seed(42)

names(pws_obs_utm)
###############
### Stage 1 ###
###############

n.iter <- 25e2
burn <- 5e2
thin <- 20

n.iter.pred <- 1e3

pws_fit_max <- grm(
  Y = pws_obs_utm$PWS_tmaxC,
  X = pws_obs_utm$tmax,
  L = pws_obs_utm[, c("NLCD", "NED")],
  M = pws_obs_utm[, c("specific_humidity", "shortwave_radiation", "wind_direction", 
                        "wind_speed")],
  n.iter = n.iter,
  burn = burn,
  thin = thin,
  covariance = "matern",
  matern.nu = 0.5,
  coords = pws_obs_utm[, c("x", "y")],
  space.id = pws_obs_utm$space_id,
  time.id = pws_obs_utm$time_id,
  spacetime.id = pws_obs_utm$spacetime_id,
  verbose.iter = 10
)

saveRDS(pws_fit_max, "/DATA/Modeling2016-2018/pws_fit_max.rds")



###############
### Stage 2 ###
###############

# Uloading the HUMID data

HUMID_for_predictions <- readRDS("/DATA/Modeling2016-2018/HUMID_for_predictions.rds")


pws_humid_pred <- grm_pred(
  grm.fit = pws_fit_max,
  X = HUMID_for_predictions$tmax,
  L = HUMID_for_predictions[, c("NLCD", "NED")],
  M = HUMID_for_predictions[, c("specific_humidity", "shortwave_radiation", "wind_direction", 
                                "wind_speed")],
  coords = HUMID_for_predictions[, c("x", "y")],
  space.id = HUMID_for_predictions$space_id,
  time.id = HUMID_for_predictions$time_id,
  spacetime.id = HUMID_for_predictions$spacetime_id,
  n.iter = 100, #n.iter.pred,
  verbose = T
)

saveRDS(pws_humid_pred, "/DATA/Modeling2016-2018/pws_humid_pred.rds")

#################################################### End Stage 2
##############################################

###############
### Stage 3 ### ordinary
###############

cv_id_pws_ORDINARY <- create_cv(
  time.id = pws_obs_utm$time_id, 
  space.id = pws_obs_utm$space_id,
  spacetime.id = pws_obs_utm$spacetime_id,
  num.folds = 10L,
  type = "ordinary"
)

saveRDS(cv_id_pws_ORDINARY, "/DATA/Modeling2016-2018/cv_id_pws_ORDINARY.rds")

cv_ORDINARY <- grm_cv(
  Y = pws_obs_utm$PWS_tmaxC,
  X = pws_obs_utm$tmax,
  cv.object = cv_id_pws_ORDINARY,
  L = pws_obs_utm[, c("NLCD", "NED")],
  M = pws_obs_utm[, c("specific_humidity", "shortwave_radiation", "wind_direction", 
                      "wind_speed")],
  n.iter = n.iter,
  burn = burn,
  thin = thin,
  coords = pws_obs_utm[, c("x", "y")],
  space.id = pws_obs_utm$space_id,
  time.id = pws_obs_utm$time_id,
  spacetime.id = pws_obs_utm$spacetime_id,
  verbose.iter = 10
)

cv_ORDINARY <- cv_ORDINARY %>% drop_na(estimate)

saveRDS(cv_ORDINARY, "/DATA/Modeling2016-2018/cv_ORDINARY.rds")


# EVALUATE performance metrics
########################################

rmse(cv_ORDINARY$obs, cv_ORDINARY$estimate)
mae(cv_ORDINARY$obs, cv_ORDINARY$estimate) 
R2(cv_ORDINARY$obs, cv_ORDINARY$estimate) 
cv_ORDINARY$coverage <- (cv_ORDINARY$lower.95<= cv_ORDINARY$obs) &
  (cv_ORDINARY$upper.95>= cv_ORDINARY$obs)
sum(cv_ORDINARY$coverage==TRUE)*100/nrow(cv_ORDINARY)

###############
### Stage 3 ### Spatial
###############
names(pws_obs_utm)
n.iter <- 10e2
burn <- 5e2
thin <- 20

n.iter.pred <- 1e3

cv_id_pws_SPATIAL <- create_cv(
  time.id = pws_obs_utm$time_id, 
  space.id = pws_obs_utm$space_id,
  spacetime.id = pws_obs_utm$spacetime_id,
  num.folds = 10L,
  type = "spatial_clustered",
  coords = pws_obs_utm[, c("x", "y")],
)

saveRDS(cv_id_pws_SPATIAL, "/DATA/Modeling2016-2018/cv_id_pws_SPATIAL_clustered.rds")


cv_SPATIAL <- grm_cv(
  Y = pws_obs_utm$PWS_tmaxC,
  X = pws_obs_utm$tmax,
  cv.object = cv_id_pws_SPATIAL,
  L = pws_obs_utm[, c("NLCD", "NED")],
  M = pws_obs_utm[, c("specific_humidity", "shortwave_radiation", "wind_direction", 
                      "wind_speed")],
  n.iter = n.iter,
  burn = burn,
  thin = thin,
  coords = pws_obs_utm[, c("x", "y")],
  space.id = pws_obs_utm$space_id,
  time.id = pws_obs_utm$time_id,
  spacetime.id = pws_obs_utm$spacetime_id,
  verbose.iter = 100
)

cv_SPATIAL <- cv_SPATIAL %>% drop_na(estimate)

saveRDS(cv_SPATIAL, "/DATA/Modeling2016-2018/cv_SPATIAL_clustered.rds")


# EVALUATE performance metrics
########################################
rmse(cv_SPATIAL$obs, cv_SPATIAL$estimate) 
mae(cv_SPATIAL$obs, cv_SPATIAL$estimate) 
R2(cv_SPATIAL$obs, cv_SPATIAL$estimate) 

cv_SPATIAL$coverage <- (cv_SPATIAL$lower.95<= cv_SPATIAL$obs) &
  (cv_SPATIAL$upper.95>= cv_SPATIAL$obs)
sum(cv_SPATIAL$coverage==TRUE)*100/nrow(cv_SPATIAL)
##################################################################################### End STAGE 3
#####################################################################################

###################################################################################################

######################################### MODELING 2: Random Forest
########################################
#####################################################################################


rf_maxTemp <- randomForest(formula= PWS_tmaxC~ tmax + NLCD + NED + specific_humidity+ shortwave_radiation+ 
wind_direction+ wind_speed+ x+ y  , data = pws_obs_utm, ntree=500, mtry=3 )
coef(rf_maxTemp)

saveRDS(rf_maxTemp,  "/DATA/Modeling2016-2018/RF_modeling/rf_maxTemp.rds")

rf_prediction  <- pws_obs_utm
##########################################################
########################## RF Stage 1: MODEL AND PREDICTIONS

# Add predictions
rf_prediction <- rf_prediction%>%
  mutate(pred_rf = predict(rf_maxTemp,.))

# Add residuals
rf_prediction <- rf_prediction%>%
  mutate(residuals_rf =  rf_prediction$PWS_tmaxC - predict(rf_maxTemp,.))


saveRDS(rf_prediction,  "/DATA/Modeling2016-2018/RF_modeling/rf_prediction.rds")

##########################################################
########################## RF Stage 2: PREDICTING HUMID

# Before MERGE original HUMID with GRM HUMID predictions

names(pws_humid_pred)
names(HUMID_for_predictions)

# Group
HUMID_GRM <- inner_join(pws_humid_pred, HUMID_for_predictions, by=c("time.id"="time_id", "space.id"="space_id"))

class(HUMID_GRM)
names(HUMID_GRM)
# RENAME variables
HUMID_GRM <- HUMID_GRM %>% rename(
  grmMax= estimate,
  grmSD= sd
)

# Convert to sf
HUMID_GRM_sf <- st_as_sf (HUMID_GRM, coords = c("x", "y"), crs=4326)

# Get lon and lat coordinates
HUMID_GRM_sf <- HUMID_GRM_sf %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2])
# Predict RF HUMID
HUMID_GRM_RF <- HUMID_GRM%>%
  mutate(rfMax = predict(rf_maxTemp,.))

saveRDS(HUMID_GRM_RF, "/DATA/Modeling2016-2018/RF_modeling/HUMID_GRM_RF.rds")

##########################################################
########################## RF Stage 3: CROSS-VALIDATION

################################################ ORDINARY 

# Define the training control with 10-fold cross-validation
train_control <- trainControl(method = "cv", number = 10)

# Train the model 
cv_RF_ordinary <- train(PWS_tmaxC~ tmax + NLCD + NED + specific_humidity+ shortwave_radiation+ 
                          wind_direction+ wind_speed+ x+ y, data = rf_prediction, method = "rf",
                        trControl = train_control)

cv_RF_ordinary_predictions <- as.data.frame(cv_RF_ordinary$trainingData)
cv_RF_ordinary_predictions$rf_cv <- cv_RF_ordinary$finalModel$predictions
cv_RF_ordinary_predictions <- cv_RF_ordinary_predictions %>% rename(
  PWS_tmaxC= .outcome)
cv_RF_ordinary_predictions <- cv_RF_ordinary_predictions%>%
  mutate(residuals_cv =  cv_RF_ordinary_predictions$PWS_tmaxC - cv_RF_ordinary_predictions$rf_cv)

cv_RF_ordinary_predictions_PWS <- as.data.frame(rf_prediction) %>% left_join(as.data.frame(cv_RF_ordinary_predictions), 
                                                                           by=c("x",
                                                                                "y",
                                                                                "tmax", "PWS_tmaxC",
                                                                                "specific_humidity","shortwave_radiation","wind_direction", "wind_speed",
                                                                                "NLCD", "NED"))

saveRDS(cv_RF_ordinary_predictions_PWS,  "/DATA/Modeling2016-2018/RF_modeling/cv_RF_ordinary_predictions_PWS.rds")

# EVALUATE performance metrics
########################################

rmse(cv_RF_ordinary_predictions$PWS_tmaxC, cv_RF_ordinary_predictions$rf_cv)
mae(cv_RF_ordinary_predictions$PWS_tmaxC, cv_RF_ordinary_predictions$rf_cv) 
R2(cv_RF_ordinary_predictions$PWS_tmaxC, cv_RF_ordinary_predictions$rf_cv) 

########################################End ORDINARY


##############################
################################################ RF SPATIAL CV for RF

cv_RF_spatial <- create_cv(
  time.id = cv_RF_ordinary_predictions_PWS$time_id, 
  space.id = cv_RF_ordinary_predictions_PWS$space_id,
  spacetime.id = cv_RF_ordinary_predictions_PWS$spacetime_id,
  num.folds = 10L,
  type = "spatial_clustered",
  coords = cv_RF_ordinary_predictions_PWS[, c("x", "y")]
)


# Upload spatial folds
df_cv_id_SPATIAL <- readRDS("/df_cv_id_SPATIAL.rds")


for(i in 1:10){
  test_indices <- which(df_cv_id_SPATIAL$cv.id == i, arr.ind = TRUE)
  train_data <- df_cv_id_SPATIAL[-test_indices, ]
  test_data <- df_cv_id_SPATIAL[test_indices, ]
  
  rf_model <- ranger(PWS_tmaxC ~ tmax + NLCD + NED + specific_humidity + 
                       shortwave_radiation + wind_direction + wind_speed + x + y, 
                     data = train_data, 
                     num.trees = 500)
  
  fold_predictions <- predict(rf_model, data = test_data)
  
  predictions[[i]] <- fold_predictions$predictions
  actuals[[i]] <- test_data$PWS_tmaxC 

  cv_id[[i]] <- test_data$cv.id
}

combined_predictions2 <- unlist(predictions)
combined_actuals2 <- unlist(actuals)
combined_cv_id <- unlist(cv_id)

predicted_df2 <- data.frame(
  id = 1:length(combined_predictions2),
  cv_id = combined_cv_id,
  PWS_tmaxC = combined_actuals2,
  rf_cv_max = combined_predictions2
)

saveRDS(predicted_df2, "/DATA/Modeling2016-2018/RF_modeling/predicted_df_cv_RF_spatial_clustered.rds")

##################################################################
# EVALUATE performance metrics
########################################

rmse(predicted_df2$PWS_tmaxC, predicted_df2$rf_cv_max) 
R2(predicted_df2$PWS_tmaxC, predicted_df2$rf_cv_max) 

##################################################################End RF SPATIAL

##################################################################################### End STAGE 3
#####################################################################################
################################################### END RANDOM FOREST

###################################################################################################

######################################### MODELING 3: Hybrid modeling: RANDOM FOREST AND BAYESIAN
########################################
#####################################################################################

install.packages("devtools")
devtools::install_github("WyattGMadden/ensembleDownscaleR")
library(ensembleDownscaleR)
set.seed(42)
# Applying Bayesian modeling
#############################################
###############
### Stage 1 ###
###############

n.iter <- 25e2
burn <- 5e2
thin <- 20

n.iter.pred <- 1e3

names(rf_prediction)

rf_grm_fit <- grm(
  Y = rf_prediction$residuals_rf,
  X = rf_prediction$tmax,
  n.iter = n.iter,
  burn = burn,
  thin = thin,
  covariance = "matern",
  matern.nu = 0.5,
  coords = rf_prediction[, c("x", "y")],
  space.id = rf_prediction$space_id,
  time.id = rf_prediction$time_id,
  spacetime.id = rf_prediction$spacetime_id,
  verbose.iter = 10
)

saveRDS(rf_grm_fit, "/DATA/Modeling2016-2018/RF_GRM_modeling/rf_grm_fit.rds")

###############
### Stage 2 ###
###############

rf_grm_pred <- grm_pred(
  grm.fit = rf_grm_fit,
  X = HUMID_GRM_RF$tmax,
  coords = HUMID_GRM_RF[, c("x", "y")],
  space.id = HUMID_GRM_RF$space.id,
  time.id = HUMID_GRM_RF$time.id,
  spacetime.id = HUMID_GRM_RF$spacetime.id,
  n.iter = 100, #n.iter.pred,
  verbose = T
)

HUMID_GRM_RF_RFGRM <- inner_join(HUMID_GRM_RF, rf_grm_pred, by=c("time.id", "space.id", "spacetime.id"))

HUMID_GRM_RF_RFGRM <- HUMID_GRM_RF_RFGRM %>% 
  rename(
    rfResiduals = estimate,
    rfResiduals_SD = sd)

HUMID_GRM_RF_RFGRM$rfgrmMax <- HUMID_GRM_RF_RFGRM$rfMax + HUMID_GRM_RF_RFGRM$rfResiduals
# Convert to sf
HUMID_GRM_RF_RFGRM_sf <- st_as_sf (HUMID_GRM_RF_RFGRM, coords = c("x", "y"), crs=4326)
HUMID_GRM_RF_RFGRM_sf <- HUMID_GRM_RF_RFGRM_sf %>%
  dplyr::mutate(x = sf::st_coordinates(.)[,1],
                y = sf::st_coordinates(.)[,2])
saveRDS(HUMID_GRM_RF_RFGRM_sf, "/DATA/Modeling2016-2018/RF_GRM_modeling/HUMID_GRM_RF_RFGRM_sf.rds")

#################################################### End Stage 2
##############################################

###############
### Stage 3 ### RF-GRM
###############

#############################
#################################
# ORDINARY

cv_id_RFGRM_ordinary <- create_cv(
  time.id = cv_RF_ordinary_predictions_PWS$time_id, 
  space.id = cv_RF_ordinary_predictions_PWS$space_id,
  spacetime.id = cv_RF_ordinary_predictions_PWS$spacetime_id,
  num.folds = 10L,
  type = "ordinary"
)

cv_RFGRM_ordinary <- grm_cv(
  Y = cv_RF_ordinary_predictions_PWS$residuals_cv,
  X = cv_RF_ordinary_predictions_PWS$tmax,
  cv.object = cv_id_RFGRM_ordinary,
  #  L = rf_pred[, c("NLCD", "NED")],
  #  M = rf_pred[, c("RH", "ceil_hgt", "cl", 
  #                  "ws")],
  n.iter = n.iter,
  burn = burn,
  thin = thin,
  coords = cv_RF_ordinary_predictions_PWS[, c("x", "y")],
  space.id = cv_RF_ordinary_predictions_PWS$space_id,
  time.id = cv_RF_ordinary_predictions_PWS$time_id,
  spacetime.id = cv_RF_ordinary_predictions_PWS$spacetime_id,
  verbose.iter = 10
)

cv_RFGRM_ordinary <- cv_RFGRM_ordinary %>% drop_na(estimate)
saveRDS(cv_RFGRM_ordinary, "/DATA/Modeling2016-2018/RF_GRM_modeling/cv_RFGRM_ordinary.rds")


# Join with rf_pred
cv_RFGRM_ordinary_PWS <- as.data.frame(cv_RF_ordinary_predictions_PWS) %>% left_join(as.data.frame(cv_RFGRM_ordinary), 
                                                                                       by=c("time_id" = "time.id",
                                                                                            "space_id" = "space.id",
                                                                                            "spacetime_id" = "spacetime.id",
                                                                                            "x", "y"))
cv_RFGRM_ordinary_PWS <- cv_RFGRM_ordinary_PWS %>% drop_na(estimate)

# rf_pred$pred_rf with rf_grm_fit_cv_ORD2$estimate
cv_RFGRM_ordinary_PWS$rf_grm_cv <- cv_RFGRM_ordinary_PWS$rf_cv + cv_RFGRM_ordinary_PWS$estimate# 90% RF predictions

# For performance metrics
cv_RFGRM_ordinary_PWS$upper95_rf_grm <- cv_RFGRM_ordinary_PWS$rf_cv + cv_RFGRM_ordinary_PWS$upper.95 # 90% RF predictions
cv_RFGRM_ordinary_PWS$lower95_rf_grm <- cv_RFGRM_ordinary_PWS$rf_cv + cv_RFGRM_ordinary_PWS$lower.95 # 90% RF predictions

# EVALUATE performance metrics
########################################

rmse(cv_RFGRM_ordinary_PWS$PWS_tmaxC, cv_RFGRM_ordinary_PWS$rf_grm_cv) 
R2(cv_RFGRM_ordinary_PWS$PWS_tmaxC, cv_RFGRM_ordinary_PWS$rf_grm_cv) 
# PEARSON CORRELATION TEST
cor.test(cv_RFGRM_ordinary_PWS$PWS_tmaxC, cv_RFGRM_ordinary_PWS$rf_grm_cv,
         method = "pearson")

cv_RFGRM_ordinary_PWS$coverage <- (cv_RFGRM_ordinary_PWS$lower95_rf_grm<= cv_RFGRM_ordinary_PWS$PWS_tmaxC) &
  (cv_RFGRM_ordinary_PWS$upper95_rf_grm>= cv_RFGRM_ordinary_PWS$PWS_tmaxC)

##################################################################################### End STAGE 3

###############################################
# SPATIAL
###############################################

cv_RFGRM_spatial <- grm_cv(
  
  Y = cv_RF_ordinary_predictions_PWS$residuals_cv,
  X = cv_RF_ordinary_predictions_PWS$tmax,
  cv.object = cv_id_pws_SPATIAL,
  #  L = rf_pred[, c("NLCD", "NED")],
  #  M = rf_pred[, c("RH", "ceil_hgt", "cl", 
  #                  "ws")],
  n.iter = n.iter,
  burn = burn,
  thin = thin,
  coords = cv_RF_ordinary_predictions_PWS[, c("x", "y")],
  space.id = cv_RF_ordinary_predictions_PWS$space_id,
  time.id = cv_RF_ordinary_predictions_PWS$time_id,
  spacetime.id = cv_RF_ordinary_predictions_PWS$spacetime_id,
  verbose.iter = 100
)


cv_RFGRM_spatial <- cv_RFGRM_spatial %>% drop_na(estimate)
saveRDS(cv_RFGRM_spatial, "/DATA/Modeling2016-2018/RF_GRM_modeling/cv_RFGRM_spatial_clustered.rds")


# Join with rf_pred
cv_RFGRM_spatial_PWS <- as.data.frame(cv_RF_ordinary_predictions_PWS) %>% left_join(as.data.frame(cv_RFGRM_spatial), 
                                                                                           by=c("time_id" = "time.id",
                                                                                                "space_id" = "space.id",
                                                                                                "spacetime_id" = "spacetime.id",
                                                                                                "x", "y"))

cv_RFGRM_spatial_PWS <- cv_RFGRM_spatial_PWS %>% drop_na(estimate)

# rf_pred$pred_rf with rf_grm_fit_cv_SPATIAL2$estimate
cv_RFGRM_spatial_PWS$rf_grm_cv <- cv_RFGRM_spatial_PWS$rf_cv + cv_RFGRM_spatial_PWS$estimate# 90% RF predictions
# 3- For performance metrics
cv_RFGRM_spatial_PWS$upper95_rf_grm <- cv_RFGRM_spatial_PWS$rf_cv + cv_RFGRM_spatial_PWS$upper.95 # 90% RF predictions
cv_RFGRM_spatial_PWS$lower95_rf_grm <- cv_RFGRM_spatial_PWS$rf_cv + cv_RFGRM_spatial_PWS$lower.95 # 90% RF predictions

# EVALUATE performance metrics
########################################

rmse(cv_RFGRM_spatial_PWS$PWS_tmaxC, cv_RFGRM_spatial_PWS$rf_grm_cv) 
R2(cv_RFGRM_spatial_PWS$PWS_tmaxC, cv_RFGRM_spatial_PWS$rf_grm_cv) 

cv_RFGRM_spatial_PWS$coverage <- (cv_RFGRM_spatial_PWS$lower95_rf_grm<= cv_RFGRM_spatial_PWS$PWS_tmaxC) &
  (cv_RFGRM_spatial_PWS$upper95_rf_grm>= cv_RFGRM_spatial_PWS$PWS_tmaxC)


###############################################END RF GRM
#####################################################################################
################################################### END Hybrid modeling: RANDOM FOREST AND BAYESIAN

saveRDS(HUMID_GRM_RF_RFGRM_sf, "/DATA/Modeling2016-2018/RF_GRM_modeling/HUMID_GRM_RF_RFGRM_sf.rds")
