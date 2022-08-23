###################################################################################################
## Main code script to estimate excess cases of coccidioidomycosis caused by or attrbiuted to drought
##  Cases are estimated in 14 counties/sub-counties within California and for two droughts, 2007-09, 2012-15
# 
## Written by Jennifer R Head
## Last updated: August 23, 2022
###################################################################################################

## STEP 0. Setup the script
# Load required packages
library(readr)
library(dplyr)
library(ggplot2)
library(splines)
library(matrixStats)
library(ranger)

# load the functions needed for the main script
source("ensemble_model_function.R")
source("counterfactual_drought_functions.R")

## STEP 1. Prepare the data

# 1A. Load the data and subset to counties of interest
setwd("C:\\Users\\Jennifer\\Documents\\Research\\Cocci\\CDPH docs")
CasesMo<-readRDS("FullCocciDataset2021.rds") %>% arrange(OnsetMonth)

# 14 counties of interest
counties <- c("WKern", "EKern", "Kings", "San Luis Obispo", "WTulare", "WFresno", "Monterey", "Merced", "WMadera", "Ventura",
                "San Joaquin", "Stanislaus", "Santa Barbara", "NLA")  #Incidence  > 10 per 100K AND >20cases/yr, last decade
# subset data
CasesMo <- CasesMo %>% subset(epi_yr %in% 2000:2019 & county2 %in% counties)

# 1B. Create some new variables for use in the models

# Create variable for Fall season following drought
CasesMo$OneYrPostDrought <- ifelse(CasesMo$OnsetMonth %in% 123:134 | #March 2010 - Feb 2011
                                     CasesMo$OnsetMonth %in% 195:206, 1, 0) #March 2016 - Feb 2017
CasesMo$TwoYrPostDrought <- ifelse(CasesMo$OnsetMonth %in% 135:146 | #March 2011 - Feb 2012
                                     CasesMo$OnsetMonth %in% 207:218, 1, 0) #March 2017 - Feb 2018
                                     
# Correct population for census tracts with population < 100 (only a few)
CasesMo$CTpopulation[CasesMo$CTpopulation < 100] <- 1000

# Calculate means in each county per month - this is used to create the counterfactual series (see counterfactual_drought_function.R)
CasesMo <- CasesMo %>% group_by(month, county2) %>%
  mutate(Tempcountymean = mean(AvgTmean, na.rm = T),
         Raincountymean = mean(TotalRain, na.rm = T)) %>% ungroup()

#get list of complete cases
CasesMo <- CasesMo %>% select(-c(Lag1YearlyRainQ, Lag2YearlyRainQ)) #remove variable that won't be used and has missingness
CasesMo <- CasesMo[complete.cases(CasesMo),]

## STEP 2. Initialize the storage objects

# table storage - store as a dataframe
table_full <- data.frame(cbind("OverallDiff"=as.numeric(), "OverallExp"=as.numeric(), "OverallDev"=as.numeric(),
                                "DuringDiff"=as.numeric(), "DuringExp"=as.numeric(), "DuringDev"=as.numeric(),
                                "AfterDiff"=as.numeric(), "AfterExp"=as.numeric(), "AfterDev"=as.numeric(), 
                                "r2"=as.numeric()))

# list of outputs from the fit_ensemble_model function
res_full <- list()

## STEP 3. Loop through each county and subset the dataset to that county.
##         Fit the ensemble model for each county and store the result in a list
##         Generate and save a plot, and store results in a data frame format

for (i in counties[1:14])  {
  
  # subset the data to a single county
  data <- subset(CasesMo, county2 == i)

  # call the function for fitting the ensemble model (see ensemble_model_function.R for details)
  res <- fit_ensemble_model(dataframe = data)

  # append the res values
  res_full <- append(res_full, res)
  
  # extract some values for inspection - sum of square errors and r2.
  SSE <- colMeans(res[[1]])
  r2 <- colMeans(res[[2]])

  # store the results also as a data frame, using summarize.table function (see ensemble_model_function.R for details)
  table_full <- rbind(table_full,
                        cbind(summarize.table(data, res), "county" = i))

  # plot and save for each county, using the sumamrize.plot function (see ensemble_model_function.R for details)
  plot <- summarize.plot(data, res, name = i)
  print(plot)
  ggsave(plot,
         file = paste0("Manuscript//Figures//EnsembleModel_CFDrought_", i,".jpg"), #update file path as needed
         dpi = 600, width = 10, height = 3)

  # print so we know which iteration we are on
  print(c(i, SSE, r2))
  
}

## STEP 4. Clean and save results for further processing/visualizations.
## Edit table rows and variable types
table_full$row <- c("2007-2009", "2012-2015")
table_full[,1:10] <- apply(table_full[,1:10], 2, FUN = function(x) as.numeric(as.character(x)))

## Save results
save(res_full, file = "ModelResults_Mar13_droughtCFs.RData")
write.csv(table_full, "ModelResults_TableSummary_Mar13_droughtCFs.csv")
