###################################################################################################
# Code file containing functions to:
#    1. Fit an ensemble prediction model
#    2. Summarize the results of the prediction model in the form of a plot
#    3. Summarize the results of the prediction model in the form of a table
#
#    These functions are utilized within the main R script 'main_counterfactual_drought_script.R'
# 
## Written by Jennifer R Head
## Last updated: August 23, 2022
###################################################################################################

## Function 1: Fit an ensemble model that predicts the monthly # of cases per census tract.

## The input to the function (dataframe) is a dataframe object. The dataset contains monthly timeseries for all census tracts 
#    in the study region, which consists of 14 counties/sub-counties. Variables used in the prediction model are:

##  - the outcome variable, cases per month per census tract, termed N_all
##  - a continuous time indicator, termed OnsetMonth, where 1 = Jan, 2000, 2 = Feb, 2000, 13 = Jan, 2001, etc.
##  - a numerical value, 1-12 for calendar month, with 1 = Jan, etc.
##  - a character variable for county name (county)
##  - variables for the mean rainfall and temperature in a given county and calendar month (termed Raincountymean, Tempcountymean)
##  - timeseries of *observed* monthly values of rainfall (TotalRain), temperature (AvgTmean), and days per month where temps exceed 30C (NHotDays30) per census tract
##  - timeseries of *observed* and *lagged* monthly values of rainfall (lagXmsTRain) and temperature (lagXmsTmean)
##  - the epi_year, which is the year spanning a full season, assumed to start April 1
##  - indicator for census tract (GEOID)
##  - binary variable for one or two year post drought (OneYrPostDrought/TwoYrPostDrought)
##  - population data for model offset (CTpopulation)

fit_ensemble_model <- function(dataframe){

  require(ranger)
  
  ## STEP 1. Set up
  
  # 1A. Set up the code so that the model can be trained using a leave-out-one-year cross validation procedure.
  blocks <- unique(dataframe$epi_yr)
  nblocks <- length(blocks)
  
  # 1B. Initialize the number of models in the ensemble. Here, six
  nModels <- 6
  
  # 1C. Initialize storage matrices 
  
  # storage for model loss (sum squared errors) and r2 values 
  loss <- r2 <- matrix(NA, nrow = nblocks, ncol = nModels + 1)
  
  # storage matrices for model predictions
  ensemble.obs <- ensemble.CF12_15 <- ensemble.CF07_09 <- #ensemble.SD12_15 <- ensemble.SA12_15 <- 
    matrix(NA, nrow = nrow(dataframe), ncol = nblocks)
  
  # 1D. Generate the counterfactual dataframes by calling the functions in counterfactual_drought_functions.R
  # These are the timeseries of temp and precip that might have occurred absent the 2007-09 and 2012-15 drought
  CF12_15 <- makeCF12_15(dataframe)
  CF07_09 <- makeCF07_09(dataframe)
  
  ## STEP 2. Loop through each of the blocks. 
  ##         Fit each of the 6 models to all but one year. 
  ##         Calculate the SSE for each of the 6 models 
  ##         Weight model predictions by how they perform, using sum of squared errors
  ##         Calculate weighted average of model predictions for observed and CF environmental timeseries

  for (i in 1:nblocks){

    ## A. Establish the training and the testing dataset
    train <- dataframe %>% subset(epi_yr != blocks[i]) # Train - all years EXCEPT one
    test <- dataframe %>% subset(epi_yr == blocks[i]) # Test - the left-out year
    
    #######################################
    ###           Model 1               ###
    #######################################
    
    # Fit Model 1: 
    #    GLM with seasonality, long-term trends, env conditons, and interactions between post drought and env conditions
    mod1 <- glm(N_all ~ offset(log(CTpopulation)) + ns(epi_yr) +
                  sand_mean + AvgImSP+ AvgElev + winter + spring + summer + 
                  lag1msTRain +
                  lag3msTRain*OneYrPostDrought +
                  lag3msTRain*TwoYrPostDrought +
                  lag6msTRain*OneYrPostDrought +
                  lag6msTRain*TwoYrPostDrought +
                  lag9msTRain*OneYrPostDrought +
                  lag9msTRain*TwoYrPostDrought +
                  lag12msTRain*OneYrPostDrought +
                  lag12msTRain*TwoYrPostDrought + 
                  lag15msTRain + lag18msTRain + lag21msTRain + lag24msTRain + 
                  lag27msTRain + lag30msTRain + lag33msTRain + lag36msTRain +
                  lag1msTmean + lag3msTmean + lag6msTmean +lag9msTmean + lag12msTmean + 
                  lag15msTmean + lag18msTmean + lag21msTmean + lag24msTmean + 
                  lag27msTmean + lag30msTmean + lag33msTmean + lag36msTmean, 
                data = train, family = "quasipoisson")

    # Generate model predictions for: 
    #   1) training set, 
    #   2) test set, 
    #   3) all observations, 
    #   4) counterfactual data (2007-09 and 2012-15)
    pred1.train <- predict(mod1, newdata = train, type = "response")
    pred1.test <- predict(mod1, newdata = test, type = "response")
    pred1.obs <- predict(mod1, newdata = dataframe, type = "response")
    pred1.CF07_09 <- predict(mod1, newdata = CF07_09, type = "response")
    pred1.CF12_15 <- predict(mod1, newdata = CF12_15, type = "response")

    # Calculate out of sample error and r squared for Model 1
    loss[i,1] <- sum((pred1.obs - dataframe$N_all)^2)
    r2[i,1] <- cor(pred1.obs, dataframe$N_all)
    
    
    #######################################
    ###           Model 2               ###
    #######################################
    
    # Fit Model 2: GLM with Long-term and Seasonal trends only (minimal model)
    mod2 <- glm(N_all ~ offset(log(CTpopulation)) + ns(epi_yr) + 
                  winter + spring + summer,
                data = train, family = "quasipoisson")
    
    # Generate model predictions for: 
    #   1) training set, 
    #   2) test set, 
    #   3) all observations, 
    #   4) counterfactual data (2007-09 and 2012-15)
    pred2.train <- predict(mod2, newdata = train, type = "response")
    pred2.test <- predict(mod2, newdata = test, type = "response")
    pred2.obs <- predict(mod2, newdata = dataframe, type = "response")
    pred2.CF07_09 <- predict(mod2, newdata = CF07_09, type = "response")
    pred2.CF12_15 <- predict(mod2, newdata = CF12_15, type = "response")
    
    # Calculate out of sample error and r squared for Model 2
    loss[i,2] <- sum((pred2.obs - dataframe$N_all)^2)
    r2[i,2] <- cor(pred2.obs, dataframe$N_all)
    
    #######################################
    ###           Model 3               ###
    #######################################
    
    # Fit Model 3: GLM with Long-term and Seasonal trends only plus lagged environmental variables
    mod3 <- glm(N_all ~ offset(log(CTpopulation)) + ns(epi_yr) + 
                  winter + spring + summer + 
                  sand_mean + AvgImSP + AvgElev + 
                  lag1msTRain + lag3msTRain + lag6msTRain +lag9msTRain + lag12msTRain+ 
                  lag15msTRain + lag18msTRain + lag21msTRain + lag24msTRain + 
                  lag27msTRain + lag30msTRain + lag33msTRain + lag36msTRain +
                  lag1msTmean + lag3msTmean + lag6msTmean +lag9msTmean + lag12msTmean + 
                  lag15msTmean + lag18msTmean + lag21msTmean + lag24msTmean + 
                  lag27msTmean + lag30msTmean + lag33msTmean + lag36msTmean, 
                data = train, family = "quasipoisson")
    
    # Generate model predictions for: 
    #   1) training set, 
    #   2) test set, 
    #   3) all observations, 
    #   4) counterfactual data (2007-09 and 2012-15)
    pred3.train <- predict(mod3, newdata = train, type = "response")
    pred3.test <- predict(mod3, newdata = test, type = "response")
    pred3.obs <- predict(mod3, newdata = dataframe, type = "response")
    pred3.CF07_09 <- predict(mod3, newdata = CF07_09, type = "response")
    pred3.CF12_15 <- predict(mod3, newdata = CF12_15, type = "response")
    
    # Calculate out of sample error and r squared for Model 3
    loss[i,3] <- sum((pred3.obs - dataframe$N_all)^2)
    r2[i,3] <- cor(pred3.obs, dataframe$N_all)
    
    
    #######################################
    ###           Model 4               ###
    #######################################
    
    # Fit Model 4: GLM with interactions between seasonal and lagged environmental variables
    mod4 <- glm(N_all ~ offset(log(CTpopulation)) + ns(epi_yr) +
                  sand_mean + AvgImSP + AvgElev +
                  
                  winter*OneYrPostDrought + spring*OneYrPostDrought + summer*OneYrPostDrought +
                  winter*TwoYrPostDrought + spring*TwoYrPostDrought + summer*TwoYrPostDrought +
                  
                  winter*lag1msTRain + spring*lag1msTRain + summer*lag1msTRain +
                  winter*lag3msTRain + spring*lag3msTRain + summer*lag3msTRain +
                  winter*lag6msTRain + spring*lag6msTRain + summer*lag6msTRain +
                  winter*lag9msTRain + spring*lag9msTRain + summer*lag9msTRain +
                  winter*lag12msTRain + spring*lag12msTRain + summer*lag12msTRain +
                  winter*lag15msTRain + spring*lag15msTRain + summer*lag15msTRain +
                  winter*lag18msTRain + spring*lag18msTRain + summer*lag18msTRain +
                  winter*lag21msTRain + spring*lag21msTRain + summer*lag21msTRain +
                  winter*lag24msTRain + spring*lag24msTRain + summer*lag24msTRain +
                  winter*lag27msTRain + spring*lag27msTRain + summer*lag27msTRain +
                  winter*lag30msTRain + spring*lag30msTRain + summer*lag30msTRain +
                  winter*lag33msTRain + spring*lag33msTRain + summer*lag33msTRain +
                  winter*lag36msTRain + spring*lag36msTRain + summer*lag36msTRain +
                  
                  winter*lag1msTmean + spring*lag1msTmean + summer*lag1msTmean +
                  winter*lag3msTmean + spring*lag3msTmean + summer*lag3msTmean +
                  winter*lag6msTmean + spring*lag6msTmean + summer*lag6msTmean +
                  winter*lag9msTmean + spring*lag9msTmean + summer*lag9msTmean +
                  winter*lag12msTmean + spring*lag12msTmean + summer*lag12msTmean +
                  winter*lag15msTmean + spring*lag15msTmean + summer*lag15msTmean +
                  winter*lag18msTmean + spring*lag18msTmean + summer*lag18msTmean +
                  winter*lag21msTmean + spring*lag21msTmean + summer*lag21msTmean +
                  winter*lag24msTmean + spring*lag24msTmean + summer*lag24msTmean +
                  winter*lag27msTmean + spring*lag27msTmean + summer*lag27msTmean +
                  winter*lag30msTmean + spring*lag30msTmean + summer*lag30msTmean +
                  winter*lag33msTmean + spring*lag33msTmean + summer*lag33msTmean +
                  winter*lag36msTmean + spring*lag36msTmean + summer*lag36msTmean, 
                data = train, family = "quasipoisson")
    
    # Generate model predictions for: 
    #   1) training set, 
    #   2) test set, 
    #   3) all observations, 
    #   4) counterfactual data (2007-09 and 2012-15)
    pred4.train <- predict(mod4, newdata = train, type = "response")
    pred4.test <- predict(mod4, newdata = test, type = "response")
    pred4.obs <- predict(mod4, newdata = dataframe, type = "response")
    pred4.CF07_09 <- predict(mod4, newdata = CF07_09, type = "response")
    pred4.CF12_15 <- predict(mod4, newdata = CF12_15, type = "response")
    
    # Calculate out of sample error and r squared for Model 4
    loss[i,4] <- sum((pred4.obs - dataframe$N_all)^2)
    r2[i,4] <- cor(pred4.obs, dataframe$N_all)
    
    #######################################
    ###           Model 5               ###
    #######################################
    
    # Fit Model 5: GLM with interactions between lagged environmental variables
    mod5 <- glm(N_all ~ offset(log(CTpopulation)) + ns(epi_yr) +
                  sand_mean + AvgImSP+ AvgElev + winter + spring + summer + 
                  lag1msTRain +
                  lag3msTRain*I(lag15msRainDev <= -0.5) +
                  lag3msTRain*I(lag27msRainDev <= -0.5) +
                  lag6msTRain*I(lag18msRainDev <= -0.5) + 
                  lag6msTRain*I(lag30msRainDev <= -0.5) +
                  lag9msTRain*I(lag21msRainDev <= -0.5) +
                  lag9msTRain*I(lag33msRainDev <= -0.5) +
                  lag12msTRain*I(lag24msRainDev <= -0.5) +
                  lag12msTRain*I(lag36msRainDev <= -0.5) + 
                  
                  lag1msTmean + lag3msTmean + lag6msTmean + lag9msTmean + lag12msTmean +
                  lag3msTRain*I(lag15msHotDays >= 7) +
                  lag3msTRain*I(lag27msHotDays >= 7) +
                  lag6msTRain*I(lag18msHotDays >= 7) + 
                  lag6msTRain*I(lag30msHotDays >= 7) +
                  lag9msTRain*I(lag21msHotDays >= 7) +
                  lag9msTRain*I(lag33msHotDays >= 7) +
                  lag12msTRain*I(lag24msHotDays >= 7) +
                  lag12msTRain*I(lag36msHotDays >= 7), 
                data = train, family = "quasipoisson")
    
    # Generate model predictions for: 
    #   1) training set, 
    #   2) test set, 
    #   3) all observations, 
    #   4) counterfactual data (2007-09 and 2012-15)
    pred5.train <- predict(mod5, newdata = train, type = "response")
    pred5.test <- predict(mod5, newdata = test, type = "response")
    pred5.obs <- predict(mod5, newdata = dataframe, type = "response")
    pred5.CF07_09 <- predict(mod5, newdata = CF07_09, type = "response")
    pred5.CF12_15 <- predict(mod5, newdata = CF12_15, type = "response")

    # Calculate out of sample error and r squared for Model 5
    loss[i,5] <- sum((pred5.obs - dataframe$N_all)^2)
    r2[i,5] <- cor(pred5.obs, dataframe$N_all)
    
    #######################################
    ###           Model 6               ###
    #######################################
    
    # Fit Model 6: Random Forest Model
    mod6 <- ranger(N_all ~ CTpopulation + epi_yr +
                  winter + spring + summer +
                  sand_mean + AvgImSP + AvgElev +
                  lag1msTRain + lag3msTRain + lag6msTRain +lag9msTRain + lag12msTRain+
                  lag15msTRain + lag18msTRain + lag21msTRain + lag24msTRain +
                  lag27msTRain + lag30msTRain + lag33msTRain + lag36msTRain +
                  lag1msTmean + lag3msTmean + lag6msTmean +lag9msTmean + lag12msTmean +
                  lag15msTmean + lag18msTmean + lag21msTmean + lag24msTmean +
                  lag27msTmean + lag30msTmean + lag33msTmean + lag36msTmean +
                  OneYrPostDrought + TwoYrPostDrought,
                data = train)

    # Generate model predictions for: 
    #   1) training set, 
    #   2) test set, 
    #   3) all observations, 
    #   4) counterfactual data (2007-09 and 2012-15)
    pred6.train <- predict(mod6, data = train)$predictions
    pred6.test <- predict(mod6, data = test)$predictions
    pred6.obs <- predict(mod6, data = dataframe)$predictions
    pred6.CF07_09 <- predict(mod6, data = CF07_09)$predictions
    pred6.CF12_15 <- predict(mod6, data = CF12_15)$predictions
    
    # Calculate out of sample error and r squared for Model 6
    loss[i,6] <- sum((pred6.obs - dataframe$N_all)^2)
    r2[i,6] <- cor(pred6.obs, dataframe$N_all)
    
    #######################################
    ###           ENSEMBLE              ###
    #######################################
    
    ## Calculate weights based on the inverse SSE - models with lower SSE get higher weights
    invLoss <- 1/loss[i,1:nModels]
    invLoss[2] <- 0 # for paper, decided not to use Model 2
    weights <- invLoss/sum(invLoss)

    ## Generate model predictions for: test set, all observations, counterfactual data - weighted average
    
    # Observed timeseries
    ensemble.obs[,i] <- pred1.obs*weights[1] + pred2.obs*weights[2] +
      pred3.obs*weights[3] + pred4.obs*weights[4] +
      pred5.obs*weights[5]  + pred6.obs*weights[6]
    
    # Counterfactual timeseries - no 2007-09 drought
    ensemble.CF07_09[,i] <- pred1.CF07_09*weights[1] + pred2.CF07_09*weights[2] +
      pred3.CF07_09*weights[3] + pred4.CF07_09*weights[4] +
      pred5.CF07_09*weights[5] + pred6.CF07_09*weights[6]
    
    # Counterfactual timeseries - no 2012-15 drought
    ensemble.CF12_15[,i] <- pred1.CF12_15*weights[1] + pred2.CF12_15*weights[2] +
      pred3.CF12_15*weights[3] + pred4.CF12_15*weights[4] + 
      pred5.CF12_15*weights[5] +pred6.CF12_15*weights[6]
    
    ## Calculate out of sample error of the ensemble model
    loss[i,7] <- sum((ensemble.obs[,i] - dataframe$N_all)^2)
    r2[i,7] <- cor(ensemble.obs[,i], dataframe$N_all)
    
  }
  
  #Return a list of the SSE, R2, the observed prediction, and the 2 counterfactual predictions
  return(list("loss" = loss,
              "r2" = r2,
              "observed" = ensemble.obs, 
              "CF07_09" = ensemble.CF07_09,
              "CF12_15" = ensemble.CF12_15))

}

###############################################################################
## Function 2: Summarize the excess or averted cases due to drought in a table
###############################################################################

summarize.table <- function(dataframe, # the dataframe input used in function 1 above
                            ensemble_res){ # the result of the ensemble model produced by function 1 above

  ## STEP 1: extract the observed and counterfactual preductions and take the prediction means across the leave-out-one-year
  obs <- ensemble_res[[3]]
  CF07_09 <- ensemble_res[[4]]
  CF12_15 <- ensemble_res[[5]]

  dataframe$Nhat <- rowMeans(obs)
  dataframe$NhatMin <- rowMins(obs)
  dataframe$NhatMax <- rowMaxs(obs)
  
  dataframe$Nhat.CF07_09 <- rowMeans(CF07_09)
  dataframe$NhatMin.CF07_09 <- rowMins(CF07_09)
  dataframe$NhatMax.CF07_09 <- rowMaxs(CF07_09)
  
  dataframe$Nhat.CF12_15 <- rowMeans(CF12_15)
  dataframe$NhatMin.CF12_15 <- rowMins(CF12_15)
  dataframe$NhatMax.CF12_15 <- rowMaxs(CF12_15)
  
  ## STEP 2: summarize the values across each month, to get the total cases per county
  plot <- dataframe %>% group_by(OnsetMonth) %>% summarize(N = sum(N_all),
                                                           Nhat = sum(Nhat),
                                                           NhatMin = sum(NhatMin),
                                                           NhatMax = sum(NhatMax),
                                                           Nhat.CF07_09 = sum(Nhat.CF07_09),
                                                           NhatMin.CF07_09 = sum(NhatMin.CF07_09),
                                                           NhatMax.CF07_09 = sum(NhatMax.CF07_09),
                                                           Nhat.CF12_15 = sum(Nhat.CF12_15),
                                                           NhatMin.CF12_15 = sum(NhatMin.CF12_15),
                                                           NhatMax.CF12_15 = sum(NhatMax.CF12_15))
  
  #Calculate overall r2
  overall.r2 <- cor(plot$N, plot$Nhat)
  
  ## STEP 3: Calculate cases averted during the 2007-2009 drought
  ## DURING: March 2007 - Nov 2009 (move to Mar 2010) (mo:87- 123)
  ## POST: April 2010 - Mar 2012 (month 124:147)
  
  # overall difference and deviation
  overallDiff.07_09 <- sum(dataframe$Nhat) - sum(dataframe$Nhat.CF07_09)
  
  overallExpected.07_09 <- sum(dataframe[which(dataframe$OnsetMonth %in% c(87:147)),]$Nhat)
  
  overallDeviation.07_09 <- overallDiff.07_09/overallExpected.07_09*100
  
  # difference and deviation during drought
  expectedDuring.07_09 <- sum(dataframe[which(dataframe$OnsetMonth %in% c(87:123)),]$Nhat)
  
  diffDuring.07_09 <- expectedDuring.07_09 - 
    sum(dataframe[which(dataframe$OnsetMonth %in% c(87:123)),]$Nhat.CF07_09)
  
  deviationDuring.07_09 <- diffDuring.07_09/expectedDuring.07_09*100

  # difference and deviation following drought (2 years)
  expectedAfter.07_09 <- sum(dataframe[which(dataframe$OnsetMonth %in% c(124:147)),]$Nhat)
  
  diffAfter.07_09 <- expectedAfter.07_09 - 
    sum(dataframe[which(dataframe$OnsetMonth %in% (124:147)),]$Nhat.CF07_09)
  
  deviationAfter.07_09 <- diffAfter.07_09/expectedAfter.07_09*100
  
  ## STEP 4: Calculate cases averted during the 2012-2015 drought
  ## DURING: May2012 - Oct 2016 (move to Mar 2016) (mo: 149- 195)
  ## POST: April 2016 - March 2018 (month 196 - 219)
  
  # overall difference and deviation
  overallDiff.12_15 <- sum(dataframe$Nhat) - sum(dataframe$Nhat.CF12_15)
  
  overallExpected.12_15 <- sum(dataframe[which(dataframe$OnsetMonth %in% c(149:219)),]$Nhat)
  
  overallDeviation.12_15 <- overallDiff.12_15/overallExpected.12_15*100
  
  # difference and deviation during drought
  expectedDuring.12_15 <- sum(dataframe[which(dataframe$OnsetMonth %in% c(149:195)),]$Nhat)
  
  diffDuring.12_15 <- expectedDuring.12_15 - 
    sum(dataframe[which(dataframe$OnsetMonth %in% c(149:195)),]$Nhat.CF12_15)
  
  deviationDuring.12_15 <- diffDuring.12_15/expectedDuring.12_15*100
  
  # difference and deviation following drought (2 years)
  expectedAfter.12_15 <- sum(dataframe[which(dataframe$OnsetMonth %in% c(196:219)),]$Nhat)
  
  diffAfter.12_15 <- expectedAfter.12_15 - 
    sum(dataframe[which(dataframe$OnsetMonth %in% (196:219)),]$Nhat.CF12_15)
  
  deviationAfter.12_15 <- diffAfter.12_15/expectedAfter.12_15*100
  
  ## STEP 5: combine into a table
  table <- rbind(cbind(overallDiff.07_09,
                       overallExpected.07_09,
                       overallDeviation.07_09,
                       diffDuring.07_09,
                       expectedDuring.07_09,
                       deviationDuring.07_09,
                       diffAfter.07_09,
                       expectedAfter.07_09,
                       deviationAfter.07_09,
                       overall.r2),
                 cbind(overallDiff.12_15,
                       overallExpected.12_15,
                       overallDeviation.12_15,
                       diffDuring.12_15,
                       expectedDuring.12_15,
                       deviationDuring.12_15,
                       diffAfter.12_15,
                       expectedAfter.12_15,
                       deviationAfter.12_15,
                       overall.r2))
  
  colnames(table) <- c("OverallDiff", "OverallExp", "OverallDev",
                       "DuringDiff", "DuringExp", "DuringDev",
                       "AfterDiff", "AfterExp", "AfterDev", "r2")
  rownames(table) <- c("2007-2009", "2012-2015")
  
  return(table)
}

###############################################################################
## Function 3: Summarize the excess or averted cases due to drought in a plot
###############################################################################
summarize.plot <- function(dataframe, # the dataframe input used in function 1 above
                           ensemble_res, # the result of the ensemble model produced by function 1 above
                           name){ #name of county (for ggtitle purposes)
  
  require(ggplot2)
  require(matrixStats)
  
  ## STEP 1: extract the observed and counterfactual preductions and take the prediction means across the leave-out-one-year
  obs <- ensemble_res[[3]]
  CF07_09 <- ensemble_res[[4]]
  CF12_15 <- ensemble_res[[5]]
  
  dataframe$Nhat <- rowMeans(obs)
  dataframe$NhatMin <- rowMins(obs)
  dataframe$NhatMax <- rowMaxs(obs)
  
  dataframe$Nhat.CF07_09 <- rowMeans(CF07_09)
  dataframe$NhatMin.CF07_09 <- rowMins(CF07_09)
  dataframe$NhatMax.CF07_09 <- rowMaxs(CF07_09)
  
  dataframe$Nhat.CF12_15 <- rowMeans(CF12_15)
  dataframe$NhatMin.CF12_15 <- rowMins(CF12_15)
  dataframe$NhatMax.CF12_15 <- rowMaxs(CF12_15)
  
  ## STEP 2: summarize the values across each month, to get the total cases per county
  plot <- dataframe %>% group_by(OnsetMonth) %>% summarize(N = sum(N_all),
                                                           Nhat = sum(Nhat),
                                                           NhatMin = sum(NhatMin),
                                                           NhatMax = sum(NhatMax),
                                                           Nhat.CF07_09 = sum(Nhat.CF07_09),
                                                           NhatMin.CF07_09 = sum(NhatMin.CF07_09),
                                                           NhatMax.CF07_09 = sum(NhatMax.CF07_09),
                                                           Nhat.CF12_15 = sum(Nhat.CF12_15),
                                                           NhatMin.CF12_15 = sum(NhatMin.CF12_15),
                                                           NhatMax.CF12_15 = sum(NhatMax.CF12_15))
  
  ## STEP 3: Make the plot
  plots <- ggplot(plot, aes(x = OnsetMonth/12+2000)) +
              geom_rect(aes(xmin = 2007 + 3/12, xmax = 2010 + 3/12, ymin = 0, ymax = Inf), fill = "gray", alpha = 0.01) +
              geom_rect(aes(xmin = 2012 + 5/12, xmax = 2016 + 3/12, ymin = 0, ymax = Inf), fill = "gray", alpha = 0.01) +
              geom_point(aes(y = N), color = "black") + 
              geom_line(aes(y = Nhat.CF07_09, color = "No 2007-2009 drought"), size = 1) +
              geom_line(aes(y = Nhat.CF12_15, color = "No 2012-2015 drought"), size = 1) +
              geom_line(aes(y = Nhat, color = "Observed"), size = 1) +
              theme_bw() +
              ylab("Monthly case count") +
              scale_fill_manual("", breaks = c("Observed", "No 2007-2009 drought", "No 2012-2015 drought"),
                                values = c("darkcyan", "thistle4", "goldenrod")) + 
              scale_color_manual("", breaks = c("Observed", "No 2007-2009 drought", "No 2012-2015 drought"),
                                values = c("darkcyan", "thistle3", "black")) + 
              ggtitle(name) + xlab("Year") 

  
  print(plots)
}


