###################################################################################################
# Code file containing functions to:
#    Produce a timeseries of counterfactual temperatures and precipitation 
#    that might have been observed if, counter to fact, the 2007-2009 and the 
#    2012-2015 drought had not been observed.
#
#    These functions are utilized within the main R script 'main_counterfactual_drought_script.R'
# 
## Written by Jennifer R Head
## Last updated: August 23, 2022
###################################################################################################

## Function 1: Create a counterfactual dataframe containing temperature and precipitation values that
##              might have been observed absent the drought that occurred 2007-2009

## The input to the function (dataframe) is a dataframe object. The dataset contains monthly timeseries for all census tracts 
#    in the study region, which consists of 14 counties/sub-counties. Variables needed are:

##  - a continuous time indicator, termed OnsetMonth, where 1 = Jan, 2000, 2 = Feb, 2000, 13 = Jan, 2001, etc.
##  - a numerical value, 1-12 for calendar month, with 1 = Jan, etc.
##  - a character variable for county name (county)
##  - variables for the mean rainfall and temperature in a given county and calendar month (termed Raincountymean, Tempcountymean)
##  - timeseries of *observed* monthly values of rainfall (TotalRain), temperature (AvgTmean), and days per month where temps exceed 30C (NHotDays30) per census tract
##  - timeseries of *observed* and *lagged* monthly values of rainfall (lagXmsTRain) and temperature (lagXmsTmean)
##  - the epi_year, which is the year spanning a full season, assumed to start April 1
##  - indicator for census tract (GEOID)
##  - binary variable for one or two year post drought (OneYrPostDrought/TwoYrPostDrought)

## For the 2007-2009 drought, OnsetMonth 87-119 is the within-drought period, 123-146 are the two years after

makeCF07_09 <- function(dataframe){
  
  ## STEP 1. To create the counterfactual data frame, we intervene deterministically on values
  ##          where observed temperature is higher than the typical county mean for that month or where
  ##          observed precipitation is lower than the typical county mean for that month.
  
  # Rename dataframe to be the counterfactual on
  NDM7 <- dataframe
  
  # 1A. Create a dataframe that contains the mean total rain (mm) and average temperature in a county within 1 of 12 calendar months
  RainMonths <- data.frame(cbind("Raincountymean" = dataframe$Raincountymean, 
                                 "Tempcountymean" = dataframe$Tempcountymean, 
                                 "month" = dataframe$month))
  
  RainMonths <- RainMonths %>% group_by(month) %>% summarize(RainMean = mean(Raincountymean),
                                                             TempMean = mean(Tempcountymean))
  NDM7 <- full_join(NDM7, RainMonths, by = "month")
  
  # 1B. During drought months only, replace rainfall less than mean with the mean, and replace temperature higher than the mean
  NDM7 <- mutate(NDM7, TotalRainDryNew = ifelse(OnsetMonth %in% 87:119 & TotalRain < RainMean, RainMean, TotalRain))
  NDM7 <- NDM7 %>% mutate(RainDevNew = (TotalRainDryNew - RainMean)/RainMean) # calculcate new deviation from mean
  
  NDM7 <- mutate(NDM7, AvgTmeanNew = ifelse(OnsetMonth %in% 87:119 & AvgTmean > TempMean, TempMean, AvgTmean))
  
  # 1C. During drought months only, if there are >7 days per month were temps exceed 30, replace with 0
  NDM7 <- mutate(NDM7, NHotDaysNew = ifelse(OnsetMonth %in% 87:119 & NHotDays30 >= 7, 0, NHotDays30)) #some # < 7
  

  ## STEP 2. Recalculate the lags using the new, counterfactual values for temperature, 
  ##         rainfall, rainfall deviaiton, and # days per month >30C
  
  monthsconsidered <- c(1,seq(3,36,3)) # lags involved are months 1-36, including 1 3, and every 3 thereafter.
  
  # 2A. Loop though the lags, and calculate the new smoothed, 3-month lags for total rain and average temperature.
  for (i in monthsconsidered){
    NDM7 <- NDM7 %>% group_by(GEOID) %>% 
      mutate(!!paste("nlag",i,"msTRain", sep = "") := (lag(TotalRainDryNew, i-1) + lag(TotalRainDryNew, i) + lag(TotalRainDryNew, i + 1))/3) %>%
      ungroup()
  }
  for (i in monthsconsidered){
    NDM7 <- NDM7 %>% group_by(GEOID) %>% 
      mutate(!!paste("nlag",i,"msTmean", sep = "") := (lag(AvgTmeanNew, i-1) + lag(AvgTmeanNew, i) + lag(AvgTmeanNew, i + 1))/3) %>%
      ungroup()
  }
  
  # 2B. Replace the obsrved lags (lagXms..) with the counterfactual lags (nlagXms...)
  NDM7$lag1msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag1msTRain, NDM7$lag1msTRain)
  NDM7$lag3msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag3msTRain, NDM7$lag3msTRain)
  NDM7$lag6msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag6msTRain, NDM7$lag6msTRain)
  NDM7$lag9msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag9msTRain, NDM7$lag9msTRain)
  NDM7$lag12msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag12msTRain, NDM7$lag12msTRain)
  NDM7$lag15msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag15msTRain, NDM7$lag15msTRain)
  NDM7$lag18msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag18msTRain, NDM7$lag18msTRain)
  NDM7$lag21msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag21msTRain, NDM7$lag21msTRain)
  NDM7$lag24msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag24msTRain, NDM7$lag24msTRain)
  NDM7$lag27msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag27msTRain, NDM7$lag27msTRain)
  NDM7$lag30msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag30msTRain, NDM7$lag30msTRain)
  NDM7$lag33msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag33msTRain, NDM7$lag33msTRain)
  NDM7$lag36msTRain <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag36msTRain, NDM7$lag36msTRain)
  
  NDM7$lag1msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag1msTmean, NDM7$lag1msTmean)
  NDM7$lag3msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag3msTmean, NDM7$lag3msTmean)
  NDM7$lag6msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag6msTmean, NDM7$lag6msTmean)
  NDM7$lag9msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag9msTmean, NDM7$lag9msTmean)
  NDM7$lag12msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag12msTmean, NDM7$lag12msTmean)
  NDM7$lag15msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag15msTmean, NDM7$lag15msTmean)
  NDM7$lag18msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag18msTmean, NDM7$lag18msTmean)
  NDM7$lag21msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag21msTmean, NDM7$lag21msTmean)
  NDM7$lag24msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag24msTmean, NDM7$lag24msTmean)
  NDM7$lag27msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag27msTmean, NDM7$lag27msTmean)
  NDM7$lag30msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag30msTmean, NDM7$lag30msTmean)
  NDM7$lag33msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag33msTmean, NDM7$lag33msTmean)
  NDM7$lag36msTmean <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag36msTmean, NDM7$lag36msTmean)
  
  # 2C. Recalculate the lags for rainfall deviation
  monthsconsidered2 <- seq(15,36,3)
  
  for (i in monthsconsidered){
    NDM7 <- NDM7 %>% group_by(GEOID) %>% 
      mutate(!!paste("nlag",i,"msRainDev", sep = "") := (lag(RainDevNew, i-1) + lag(RainDevNew, i) + lag(RainDevNew, i + 1))/3) %>%
      ungroup()
  }
  for (i in monthsconsidered){
    NDM7 <- NDM7 %>% group_by(GEOID) %>% 
      mutate(!!paste("nlag",i,"msHotDays", sep = "") := (lag(NHotDaysNew, i-1) + lag(NHotDaysNew, i) + lag(NHotDaysNew, i + 1))/3) %>%
      ungroup()
  }
  
  # 2D. Replace the current, observed lagged values with the counterfactual lagged values
  NDM7$lag15msRainDev <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag15msRainDev, NDM7$lag15msRainDev)
  NDM7$lag18msRainDev <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag18msRainDev, NDM7$lag18msRainDev)
  NDM7$lag21msRainDev <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag21msRainDev, NDM7$lag21msRainDev)
  NDM7$lag24msRainDev <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag24msRainDev, NDM7$lag24msRainDev)
  NDM7$lag27msRainDev <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag27msRainDev, NDM7$lag27msRainDev)
  NDM7$lag30msRainDev <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag30msRainDev, NDM7$lag30msRainDev)
  NDM7$lag33msRainDev <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag33msRainDev, NDM7$lag33msRainDev)
  NDM7$lag36msRainDev <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag36msRainDev, NDM7$lag36msRainDev)
  
  NDM7$lag15msHotDays <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag15msHotDays, NDM7$lag15msHotDays)
  NDM7$lag18msHotDays <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag18msHotDays, NDM7$lag18msHotDays)
  NDM7$lag21msHotDays <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag21msHotDays, NDM7$lag21msHotDays)
  NDM7$lag24msHotDays <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag24msHotDays, NDM7$lag24msHotDays)
  NDM7$lag27msHotDays <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag27msHotDays, NDM7$lag27msHotDays)
  NDM7$lag30msHotDays <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag30msHotDays, NDM7$lag30msHotDays)
  NDM7$lag33msHotDays <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag33msHotDays, NDM7$lag33msHotDays)
  NDM7$lag36msHotDays <- ifelse(NDM7$epi_yr > 2005, NDM7$nlag36msHotDays, NDM7$lag36msHotDays)
  
  ## STEP 3. Recalculate the bindary indicator valariable for post drought (replace 1 with 0)
  NDM7$OneYrPostDrought <- ifelse(NDM7$OnsetMonth %in% 123:134, 0, NDM7$OneYrPostDrought)
  NDM7$TwoYrPostDrought <- ifelse(NDM7$OnsetMonth %in% 135:146, 0, NDM7$TwoYrPostDrought)
  
  # return the counterfactual dataframe
  return(NDM7)
  
}

## Function 2: Create a counterfactual dataframe containing temperature and precipitation values that
##              might have been observed absent the drought that occurred 2012-2015

## The input to the function (dataframe) is a dataframe object. The dataset contains monthly timeseries for all census tracts 
#    in the study region, which consists of 14 counties/sub-counties. Variables needed are:

##  - a continuous time indicator, termed OnsetMonth, where 1 = Jan, 2000, 2 = Feb, 2000, 13 = Jan, 2001, etc.
##  - a numerical value, 1-12 for calendar month, with 1 = Jan, etc.
##  - a character variable for county name (county)
##  - variables for the mean rainfall and temperature in a given county and calendar month (termed Raincountymean, Tempcountymean)
##  - timeseries of *observed* monthly values of rainfall (TotalRain), temperature (AvgTmean), and days per month where temps exceed 30C (NHotDays30) per census tract
##  - timeseries of *observed* and *lagged* monthly values of rainfall (lagXmsTRain) and temperature (lagXmsTmean)
##  - the epi_year, which is the year spanning a full season, assumed to start April 1
##  - indicator for census tract (GEOID)
##  - binary variable for one or two year post drought (OneYrPostDrought/TwoYrPostDrought)

## For the 2012-15 drought, OnsetMonth 149-190 is the within-drought period, 195-218 are the two years after

makeCF12_15 <- function(dataframe){
  
  ## STEP 1. To create the counterfactual data frame, we intervene deterministically on values
  ##          where observed temperature is higher than the typical county mean for that month or where
  ##          observed precipitation is lower than the typical county mean for that month.
  
  # Rename dataframe
  NDM <- dataframe
  
  # 1A. Create a dataframe that contains the mean total rain (mm) and average temperature in a county within 1 of 12 calendar months
  RainMonths <- data.frame(cbind("Raincountymean" = dataframe$Raincountymean, 
                                 "Tempcountymean" = dataframe$Tempcountymean, 
                                 "month" = dataframe$month))
  
  RainMonths <- RainMonths %>% group_by(month) %>% summarize(RainMean = mean(Raincountymean),
                                                             TempMean = mean(Tempcountymean))
  NDM <- full_join(NDM, RainMonths, by = "month")
  
  # 1B. During drought months only, replace rainfall less than mean with the mean, and replace temperature higher than the mean
  NDM <- mutate(NDM, TotalRainDryNew = ifelse(OnsetMonth %in% 149:190 & TotalRain < RainMean, RainMean, TotalRain))
  NDM <- NDM %>% mutate(RainDevNew = (TotalRainDryNew - RainMean)/RainMean) #calculate new rainfall deviation from mean
  
  NDM <- mutate(NDM, AvgTmeanNew = ifelse(OnsetMonth %in% 149:190 & AvgTmean > TempMean, TempMean, AvgTmean))
  
  # 1C. During drought months only, if there are >7 days per month were temps exceed 30, replace with 0
  NDM <- mutate(NDM, NHotDaysNew = ifelse(OnsetMonth %in% 149:190 & NHotDays30 >= 7, 0, NHotDays30)) #some # < 7
  
  ## STEP 2. Recalculate the lags using the new, counterfactual values for temperature, 
  ##         rainfall, rainfall deviaiton, and # days per month >30C
  
  monthsconsidered <- c(1,seq(3,36,3)) # lags involved are months 1-36, including 1 3, and every 3 thereafter.
  
  # 2A. Loop though the lags, and calculate the new smoothed, 3-month lags for total rain and average temperature.
  for (i in monthsconsidered){
    NDM <- NDM %>% group_by(GEOID) %>% 
      mutate(!!paste("nlag",i,"msTRain", sep = "") := (lag(TotalRainDryNew, i-1) + lag(TotalRainDryNew, i) + lag(TotalRainDryNew, i + 1))/3) %>%
      ungroup()
  }
  for (i in monthsconsidered){
    NDM <- NDM %>% group_by(GEOID) %>% 
      mutate(!!paste("nlag",i,"msTmean", sep = "") := (lag(AvgTmeanNew, i-1) + lag(AvgTmeanNew, i) + lag(AvgTmeanNew, i + 1))/3) %>%
      ungroup()
  }
  
  # 2B. Replace the observed lags (lagXms..) with the counterfactual lags (nlagXms...)
  NDM$lag1msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag1msTRain, NDM$lag1msTRain)
  NDM$lag3msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag3msTRain, NDM$lag3msTRain)
  NDM$lag6msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag6msTRain, NDM$lag6msTRain)
  NDM$lag9msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag9msTRain, NDM$lag9msTRain)
  NDM$lag12msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag12msTRain, NDM$lag12msTRain)
  NDM$lag15msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag15msTRain, NDM$lag15msTRain)
  NDM$lag18msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag18msTRain, NDM$lag18msTRain)
  NDM$lag21msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag21msTRain, NDM$lag21msTRain)
  NDM$lag24msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag24msTRain, NDM$lag24msTRain)
  NDM$lag27msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag27msTRain, NDM$lag27msTRain)
  NDM$lag30msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag30msTRain, NDM$lag30msTRain)
  NDM$lag33msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag33msTRain, NDM$lag33msTRain)
  NDM$lag36msTRain <- ifelse(NDM$epi_yr > 2010, NDM$nlag36msTRain, NDM$lag36msTRain)
  
  NDM$lag1msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag1msTmean, NDM$lag1msTmean)
  NDM$lag3msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag3msTmean, NDM$lag3msTmean)
  NDM$lag6msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag6msTmean, NDM$lag6msTmean)
  NDM$lag9msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag9msTmean, NDM$lag9msTmean)
  NDM$lag12msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag12msTmean, NDM$lag12msTmean)
  NDM$lag15msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag15msTmean, NDM$lag15msTmean)
  NDM$lag18msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag18msTmean, NDM$lag18msTmean)
  NDM$lag21msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag21msTmean, NDM$lag21msTmean)
  NDM$lag24msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag24msTmean, NDM$lag24msTmean)
  NDM$lag27msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag27msTmean, NDM$lag27msTmean)
  NDM$lag30msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag30msTmean, NDM$lag30msTmean)
  NDM$lag33msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag33msTmean, NDM$lag33msTmean)
  NDM$lag36msTmean <- ifelse(NDM$epi_yr > 2010, NDM$nlag36msTmean, NDM$lag36msTmean)
  
  # 2C. Recalculate the lags for rainfall deviation
  monthsconsidered2 <- seq(15,36,3)
  
  for (i in monthsconsidered){
    NDM <- NDM %>% group_by(GEOID) %>% 
      mutate(!!paste("nlag",i,"msRainDev", sep = "") := (lag(RainDevNew, i-1) + lag(RainDevNew, i) + lag(RainDevNew, i + 1))/3) %>%
      ungroup()
  }
  for (i in monthsconsidered){
    NDM <- NDM %>% group_by(GEOID) %>% 
      mutate(!!paste("nlag",i,"msHotDays", sep = "") := (lag(NHotDaysNew, i-1) + lag(NHotDaysNew, i) + lag(NHotDaysNew, i + 1))/3) %>%
      ungroup()
  }
  
  # 2D. Replace the current, observed lagged values with the counterfactual lagged values
  NDM$lag15msRainDev <- ifelse(NDM$epi_yr > 2005, NDM$nlag15msRainDev, NDM$lag15msRainDev)
  NDM$lag18msRainDev <- ifelse(NDM$epi_yr > 2005, NDM$nlag18msRainDev, NDM$lag18msRainDev)
  NDM$lag21msRainDev <- ifelse(NDM$epi_yr > 2005, NDM$nlag21msRainDev, NDM$lag21msRainDev)
  NDM$lag24msRainDev <- ifelse(NDM$epi_yr > 2005, NDM$nlag24msRainDev, NDM$lag24msRainDev)
  NDM$lag27msRainDev <- ifelse(NDM$epi_yr > 2005, NDM$nlag27msRainDev, NDM$lag27msRainDev)
  NDM$lag30msRainDev <- ifelse(NDM$epi_yr > 2005, NDM$nlag30msRainDev, NDM$lag30msRainDev)
  NDM$lag33msRainDev <- ifelse(NDM$epi_yr > 2005, NDM$nlag33msRainDev, NDM$lag33msRainDev)
  NDM$lag36msRainDev <- ifelse(NDM$epi_yr > 2005, NDM$nlag36msRainDev, NDM$lag36msRainDev)
  
  NDM$lag15msHotDays <- ifelse(NDM$epi_yr > 2005, NDM$nlag15msHotDays, NDM$lag15msHotDays)
  NDM$lag18msHotDays <- ifelse(NDM$epi_yr > 2005, NDM$nlag18msHotDays, NDM$lag18msHotDays)
  NDM$lag21msHotDays <- ifelse(NDM$epi_yr > 2005, NDM$nlag21msHotDays, NDM$lag21msHotDays)
  NDM$lag24msHotDays <- ifelse(NDM$epi_yr > 2005, NDM$nlag24msHotDays, NDM$lag24msHotDays)
  NDM$lag27msHotDays <- ifelse(NDM$epi_yr > 2005, NDM$nlag27msHotDays, NDM$lag27msHotDays)
  NDM$lag30msHotDays <- ifelse(NDM$epi_yr > 2005, NDM$nlag30msHotDays, NDM$lag30msHotDays)
  NDM$lag33msHotDays <- ifelse(NDM$epi_yr > 2005, NDM$nlag33msHotDays, NDM$lag33msHotDays)
  NDM$lag36msHotDays <- ifelse(NDM$epi_yr > 2005, NDM$nlag36msHotDays, NDM$lag36msHotDays)
  
  ## STEP 3. Recalculate the bindary indicator valariable for post drought (replace 1 with 0)
  NDM$OneYrPostDrought <- ifelse(NDM$OnsetMonth %in% 195:206, 0, NDM$OneYrPostDrought)
  NDM$TwoYrPostDrought <- ifelse(NDM$OnsetMonth %in% 207:218, 0, NDM$TwoYrPostDrought)
  
  # return counterfactual dataframe
  return(NDM)
}
