# ValleyFever_and_Drought

This repo contains the code associated with the paper, *Drought in Western United States displaces and amplifies coccidioidomycosis, an emerging fungal disease: a longitudinal surveillance study* by Head, et al., published in 2022 in Lancet Planetary Health.

Human case data are protected health information (PHI) with access restricted to authorized California Department of Public Health (CDPH) staff. This repository thus contains only the code used for analyses, and not any of the datasets involved. For those seeking to conduct a similar multi-location study using distributed lag non-linear models (DLNMs), consider referencing the paper *Multivariate meta-analysis for non-linear and other multi-parameter associations* by Gasparrini, et al. in Stat Med, which contains both code and associated data.

The code is divided into two folders:

1. **Counterfactual_Drought_Analysis** contains the code to estimate the number of cases of coccidioidomycosis averted by or attributed to the 2007-2009 and the 2012-2015 droughts in California. Within this folder, the R scripts are as follows:

 - **main_counterfactual_drought_script.R** is the main script that performs the analysis, and calls the functions located in the second script.

 - **ensemble_model_function.R** contains functions that fit the ensemble model and summarize outcomes
 
 - **counterfactual_drought_functions.R** contains functions that produce a timeseries of counterfactual temperatures and precipiation that might have been observed if, counter to fact, the 2007-2009 and the 2012-2015 drought had not been observed.

2. **Meteorological_DLNM_Analysis** contains the code to estimate delayed and interactive associations between temperature, precipitation and coccidioidomycosis incidence. The three scrips are as follows:

 - **DLNM_MetaAnal_1.36Lag_Rainfall.R** examines the association between precipitation lagged 1-36 months and coccidioidomycosis incidence. This is the best code to use to understand the main analysis that produced Figure 2 in the manuscript. This script can be easily modified to change model formula and specification of the spline structure. This code examines multiple lags, and uses a univariate meta-analysis to obtain pooled IRRs. 

 - **DLNM_MultiVarMetaAnal_WinterRain_SummerTemp.R** examines variables that explain heterogeneity in the exposure-response relationship between precipitation and temperature *at a single lag* with incidence across counties. It is used to produce Figure 3 in the manuscript. In particular, it examines heterogeneity in the relationship between rainfall in the prior winter (lagged by 9 months) and incidence, and temperature in the prior summer (lagged by 3 months) and incidence. This code is similar to that in DLNM_MetaAnal_1.36Lag_Rainfall.R, but focuses on one lag and uses multivariate meta-analysis to understand effect modification by spatial factors.
 
  - **EffectModificationRainfall_AntecedentConditions.R** examines how antecedent conditions (e.g., those experienced >12 months ago) modifies the the exposure-response relationship between precipitation *at a single lag* and incidence. In particular, it examines heterogeneity in the relationship between rainfall in the prior winter (lagged by 9 months) and incidence, when the dataset is stratified by rainfall in the winter 2-3 years prior. This code is similar to that in DLNM_MetaAnal_1.36Lag_Rainfall.R, but focuses on one lag stratifies the data to understand effect modification by temporal factors.


