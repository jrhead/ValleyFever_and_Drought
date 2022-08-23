# ValleyFever_and_Drought

This repo contains the code associated with the paper, *Drought in Western United States displaces and amplifies coccidioidomycosis, an emerging fungal disease: a longitudinal surveillance study* by Head, et al., published in Lancet Planetary Health.

Human case data are protected health information (PHI) with access restricted to authorized California Department of Public Health (CDPH) staff. This repository thus contains only the code used for analyses, and not any of the datasets involved. For those seeking to conduct a similar multi-location study using distributed lag non-linear models (DLNMs), consider referencing the paper *Multivariate meta-analysis for non-linear and other multi-parameter associations* by Gasparrini, et al. in Stat Med, which contains both code and associated data.

The code is divided into two folders:

1. **Counterfactual_Drought_Analysis** contains the code to estimate the number of cases of coccidioidomycosis averted by or attributed to the 2007-2009 and the 2012-2015 droughts in California. Within this folder, the R scripts are as follows:

 - **main_counterfactual_drought_script.R** is the main script that performs the analysis, and calls the functions located in the second script.

 - **ensemble_model_function.R** contains functions that fit the ensemble model and summarize outcomes
 
 - **counterfactual_drought_functions.R** contains functions that produce a timeseries of counterfactual temperatures and precipiation that might have been observed if, counter to fact, the 2007-2009 and the 2012-2015 drought had not been observed.

2. **Meteorological_DLNM_Analysis** contains the code to estimate delayed and interactive associations between temperature, precipiation and coccidioidomycosis incidence.


