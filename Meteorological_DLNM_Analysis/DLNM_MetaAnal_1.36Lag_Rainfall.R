###################################################################################################
## Main code script to estimate associations between precipitation and coccidioidomycosis incidence
##  The method uses a two-stage meta-analytic approach with distributed lag non-linear models
##  The code is adapted from Gasparrini, et al., Stat Med: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3546395/.
# 
## Written by Jennifer R Head
## Last updated: August 23, 2022
###################################################################################################

## STEP 0. LOAD PACKAGES
library(dlnm) ; library(splines) ; library(tsModel)
library(readr); library(dplyr); library(ggplot2);
library(mvmeta)

## STEP 1. LOAD AND PREPARE DATA

# The data here is pre-prepared to have lags already calculated, and to be subset so
# that outcome variable (N_all) is only for months Sept - Nov, and to contain information
# only on the 14 counties/sub-counties of interest

# The data contains a timeseries of month cases for all census tracts in CA, linked
# to data on environmental conditions.

# 1A. Load the data
CasesMo<-readRDS("DatasetForGAMAnalysis.rds") %>% arrange(OnsetMonth) %>% subset(!is.na(N_all))

# adjust population variable
CasesMo$CTpopulation[CasesMo$CTpopulation < 100] <- 1000

# 1B. Define the locations for each of the models as counties
dataFrame <- CasesMo
dataFrame$county <- as.character(dataFrame$county2)
counties <- unique(dataFrame$county)
m <- length(counties)

## STEP 2. INITIALIZE STORAGE VECTORS

# storage vector where we will store the 25th and the 75th percentile of the main 
# exposure -- here 36 lags for precipitation. We pick 25 and 75 so we can calculate
# the IRR associated with an IQR increase (aka from the 25th to the 75th percentile)
per25 <- per75 <- rep(NA, 36)

# storage vector to store the IRR associated with an IQR increase in the exposure
# of interest -- here rainfall, lagged 1-36 months. IRR_75 is the point estimate
# and the low and high stores the upper and lower 95% CI's
IRR_75 <- lowIRR_75 <- highIRR_75<- rep(NA, 36)


## STEP 3. LOOP THROUGH EACH LAGGED VALUE OF INTEREST, 
##         RUN FIRST AND SECOND STAGE MODEL, 
##         EXTRACT IRR

for (index in 1:36){ # 1-36 month lag
  
  # 1. EXTRACT THE EXPOSURE VARIABLE OF INTEREST, THE 25TH and 75TH %ILE, AND OTHER MODEL TERMS
  
  # Extract the exposure variable of interest -- total monthly rain lagged by 1-36 months
  dataFrame$exposure <- pull(dataFrame, paste0("lag", index, "msTRain")) # note: User can replace with other variables, like temp
  
  # Extract the 25th and the 75th percentile of this variable from the data and store it
  iqr25 <- unname(quantile(dataFrame$exposure, prob = c(0.25), na.rm = T))
  per25[index] <- iqr25
  iqr75 <- unname(quantile(dataFrame$exposure, prob = c(0.75), na.rm = T))
  per75[index] <- iqr75
  
  # Define the other variables that we want to control for in the model: here mean monthly temps lagged by 3-33 months, and % sand
  dataFrame$ME1 <- dataFrame$lag3msTmean
  dataFrame$ME2 <- dataFrame$lag9msTmean
  dataFrame$ME3 <- dataFrame$lag15msTmean
  dataFrame$ME4 <- dataFrame$lag21msTmean
  dataFrame$ME5 <- dataFrame$lag27msTmean
  dataFrame$ME6 <- dataFrame$lag33msTmean
  dataFrame$ME7 <- dataFrame$sand_mean
  NumME <- 7 # number of other variables to control for in the model
  
  # 2. DEFINE THE MODEL FORMULA
  
  # N_all is cases of cocci
  # splinevar will be the spline on the exposure variable (defined later)
  # log offset of population allows for interpretation of coefficients as IRRs
  # control for long-tern trends on year
  
  formula <- "N_all ~ splinevar + log(offset(CTpopulation)) + ns(year) + ns(ME1, knots = 1)"
  
  # add the other confounders ME2-7 to the formula
  for (i in 2:NumME){
    formula <- paste(formula, " + ns(ME", i, ", knots=1)", sep = "")
  }
  
  # Add in controls for precipitation at lags other than the lag for the main exposure
  
  # First, define these other lags as every 3rd month, except that closest to the one with the main exposure
  if(index == 1){
    other_vars <- c(3,6,9,12,15,18,21,24,27,30,33)
  } else{
    other_vars <- (c(3,6,9,12,15,18,21,24,27,30,33))[-round(index/3,0)]}

  # Then, add these other lags for rainfall into the model
  for (i in unique(other_vars)){
    formula <- paste(formula, " + ns(lag", i, "msTRain , knots=1)", sep = "")
  }
  
  # Set formula as formula
  formula <- as.formula(formula)
  
  # 3. DEFINE COUNTIES AND CREATE LIST OF DF
  
  rm(dat_list) #remove prior instances of dat_list 
  dat_list<-split(dataFrame, dataFrame$county) # make a list of the dataframe with selected counties
  
  # RANGES (FOR VARIABLE OF INTEREST -- CAN BE CHANGED BY THE USER, here 0 to 90th percentile)
  ranges <- t(sapply(dat_list,function(x) quantile(x$exposure, probs = c(0.0, 0.90),na.rm=T)))
  
  # DEFINE THE AVERAGE RANGE, CENTERING POINT, DEGREE AND TYPE OF THE SPLINE
  # (THESE PARAMETERS CAN BE CHANGED BY THE USER FOR ADDITIONAL ANALYSES)
  cen <- iqr25 # center it at the 25th percentile
  (bound <- c(min(ranges[,1]), max(ranges[,2]))) # ANOTHER WAY TO DEFINE RANGE
  degree <- 2
  type <- "bs"
  
  # 4. DEFINE KNOTS AS ABSOLUTE VALS
  # UPDATE THESE VALS ACCORDING TO VARIABLE CONSIDERED
  allvals <- allknots <- seq(from = min(ranges[,1]) + 1, to = max(ranges[,2])-1, 
                             by = ((max(ranges[,2])-1) - (min(ranges[,1]) + 1))/10)
  
  ####################################################################
  # 5. DETERMINE HOW MANY KNOTS GIVES THE LOWEST QAIC AND WHERE
  
  # COMBINATIONS OF KNOTS FOR INCREASING DF, AT MOST 3 KNOTS
  #ranges #observe if error
  comb <- list(c(1),c(2),c(3),c(4),c(5),c(6),c(7),c(8),c(9),c(10),c(11))
  
  # BUILD THE MATRIX OF Q-AIC VALUES FOR EACH CITY/COMBINATION
  qaicmat <- matrix(0,m,length(comb),dimnames=list(counties,NULL))
  
  # RUN THE MODEL FOR EACH COUNTY, FOR EACH COMBINATION
  
  # LOOP FOR COUNTY
  system.time(
    for(i in seq(m)) {
      
      # LOAD DATA
      data <- dat_list[[i]]
      
      # LOOP FOR KNOTS COMBINATIONS
      for(j in seq(comb)) {
        
        if(sum(data$N_all)>20){
          
          # SELECT THE KNOTS
          vknots <- allvals[comb[[j]]] #allvals or allperc
          # CREATE THE CENTERED SPLINE
          # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
          suppressWarnings(
            splinevar <- onebasis(data$exposure,#type=type,degree=degree, #these two options give error so letting it be default
                                  knots=vknots,bound=bound,cen=cen)
          )
          
          # RUN THE MODEL
          model <- glm(formula, family=quasipoisson(), data)
          
          # COMPUTE AND SAVE THE Q-AIC VALUE
          loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
          phi <- summary(model)$dispersion
          qaicmat[i,j] <- -2*loglik + 2*summary(model)$df[3]*phi
        }else {qaicmat[i,j] <- NA}
      } 
      
    })
  
  ####################################################################
  # RESULTS COMPARING Q_AIC BY KNOT LOCATION AND DF
  
  # SUM THE Q-AIC ACROSS COUNTIES
  qaictot <- colSums(qaicmat, na.rm = T)
  
  # CHOSEN DF, PERCENTILES AND KNOTS AT MINIMUM Q-AIC
  df <- length(comb[[which.min(qaictot)]]) + 1 + (type=="bs")*(degree-1)
  perc <- allvals[comb[[which.min(qaictot)]]] #allvalls OR allperc
  (knots <- allknots[comb[[which.min(qaictot)]]])
  
  ####################################################################
  # 6. RUN THE FIRST STAGE MODEL
  ####################################################################
  
  # BUILT OBJECTS WHERE RESULTS WILL BE STORED:
  #   ymat IS THE MATRIX FOR THE OUTCOME PARAMETERS
  #   Slist IS THE LIST WITH (CO)VARIANCE MATRICES
  ymat <- matrix(NA,m,df-1,dimnames=list(counties,paste("spl",seq(df-1),sep="")))
  Slist <- vector("list",m)
  names(Slist) <- counties
  
  ####################################################################
  # RUN THE FIRST-STAGE ANALYSIS
  
  system.time(
    for(i in seq(m)) {
      
      # LOAD
      data <- dat_list[[i]]
      if(sum(data$N_all)>20){
        # CREATE THE CENTERED SPLINE NB: KNOTS AND BOUNDARIES FIXED AT SAME VALUES
        # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
        suppressWarnings(
          splinevar <- onebasis(data$exposure,#type=type,degree=degree,
                                knots=knots,bound=bound,cen=cen)
        )
        
        # RUN THE MODEL
        model <- glm(formula, family=quasipoisson(),data)
        summary(model)
        # EXTRACT AND SAVE THE RELATED COEF AND VCOV
        # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
        suppressWarnings(
          predsplinevar <- crosspred(splinevar,model)
        )
        ymat[i,] <- predsplinevar$coef
        Slist[[i]] <- predsplinevar$vcov
      }
      else {ymat[i,] <- NA
      Slist[[i]] <- NA}
    })
  
  #
  ####################################################################
  # 7) RUN THE MODELS WITH mvmeta
  #
  # 8) CREATE BASIS VARIABLES USING onebasis, USED FOR PREDICTION
  #
  # 9) PREDICT THE OUTCOME PARAMETERS OVER SPECIFIC VALUES OF STUDY-LEVEL
  #   COVARIATES USING predict (mvmeta),THEN RE-BUILD THE PREDICTED CURVE 
  #   USING FOR THOSE crosspred AGAIN
  #
  # NOTE: THE USE OF dlnm FUNCTIONS FACILITATES PREDICTION AND PLOTTING
  #
  ####################################################################
  
  ####################################################################
  # 7. PERFORM META-ANALYSIS
  ####################################################################
  
  mv <- mvmeta(ymat,Slist,method="reml",bscov = "fixed") #would prefer unstr for bscov...might need fixed??
  summary(mv)
  
  # MULTIVARIATE META-REGRESSION
  # Note: we perform this in another script -- but we could do it here if desired.
  
  ####################################################################
  # 8. CREATE BASIS FOR PREDICTION
  ####################################################################
  
  # BASIS USED TO PREDICT VARIABLE, EQUAL TO THAT USED FOR ESTIMATION
  #   NOTE: INTERNAL AND BOUNDARY KNOTS PLACED AT SAME VALUES AS IN ESTIMATION
  #   NOTE: SPLINE CENTERED ON THE SAME VALUE AS IN ESTIMATION
  predictorvar <- seq(bound[1],bound[2],length=30)
  bpredictorvar <- onebasis(predictorvar,#type=type,degree=degree,
                            knots=knots,
                            bound=bound,cen=cen)
  
  ####################################################################
  # 9. PREDICTION FROM MODELS
  ####################################################################
  
  # USE OF crosspred TO PREDICT THE EFFECTS FOR THE CHOSEN VALUES
  
  # PREDICTION FROM SIMPLE META-ANALYSES WITH NO PREDICTORS
  cp <- crosspred(bpredictorvar,coef=coef(mv),vcov=vcov(mv),model.link="log",by=0.1)
  
  # 10. EXTRACT VALUES OF INTEREST AND STORE THEM
  
  # store vectors of the IRRs, lower and upper 95% CI on the IRR across values of the exposure
  a <- cp$matRRfit # vector of IRRs at various values of the exposure
  TotalRain <- rownames(a) 
  dfa <- data.frame(cbind(a, TotalRain))
  dfa[,1] <- as.numeric(as.character(dfa[,1]))
  dfa[,2] <- as.numeric(as.character(dfa[,2]))
  
  low <- cp$matRRlow # vector of the lower 95% CI's of the IRRs
  dflow <- data.frame(cbind(low, TotalRain))
  high <- cp$matRRhigh # vector of the upper 95% CI's of the IRRs
  dfhigh <- data.frame(cbind(high, TotalRain))
  
  # extract the IRR's, lower and upper 95% CI on the IRR at the 75th percentile. Recall, IRR=1 at 25th percentile due to centering
  IRR_75[index] <- as.numeric(as.character(dfa[which(dfa$TotalRain == round(iqr75,1)),]$lag0))
  lowIRR_75[index] <- as.numeric(as.character(dflow[which(dflow$TotalRain == round(iqr75,1)),]$lag0))
  highIRR_75[index] <- as.numeric(as.character(dfhigh[which(dfhigh$TotalRain == round(iqr75,1)),]$lag0))
  
  
  print(index)
  
  ## OPTIONAL COADE TO PLOT
  # plot(cp,"overall",col=1,lwd=2,ylab="IRR",ylim=c(0,max(dfa[,1])+1),xlim=c(min(ranges[,1]),max(ranges[,2])),
  #      xlab=paste("Total Rain (mm), lagged", index, "months"))
  # points(iqr25,1,pch=19,cex=1)
  # points(round(iqr75,1),cp$allRRfit[as.character(round(iqr75,1))],pch=19,cex=1)
  # abline(h = cp$allRRfit[as.character(round(iqr75,1))], lty = "dotted")
  
  
}

## STEP 4. SAVE THE MODEL OUTPUT IN A .CSV FOR FUTHER MANIUPULATION (e.g, plots, etc.)

sig <- ifelse(lowIRR_75 >= 1 | highIRR_75 <= 1, 1, 0) # is IRR significant at the 95% confidence level?

# gather variables into a dataframe

Univar_Vals <- data.frame(cbind(c(rep("Year0", 12), rep("Year1", 12), rep("Year2", 12)), # year
                                c(1:36), 1, # continuous lag
                                c(1:12,1:12,1:12), # lag within a year
                                c(IRR_75), c(lowIRR_75), c(highIRR_75), # IRR, lower and upper CI
                                c(per25), c(per75), # 25th and 75th percentile of exposure
                                c(sig))) # binary variable for significance
colnames(Univar_Vals) <- c("Name", "ContLag", "one", "lag", 
                           "IRR", "lowerCI", "upperCI", 
                           "Lower25th", "Upper25th",   
                           "sig")
# Ensure variables are numeric
Univar_Vals$lag <- as.numeric(as.character(Univar_Vals$lag))
Univar_Vals$ContLag <- as.numeric(as.character(Univar_Vals$ContLag))
Univar_Vals$IRR <- as.numeric(as.character(Univar_Vals$IRR))
Univar_Vals$Lower25th <- as.numeric(as.character(Univar_Vals$Lower25th))
Univar_Vals$Upper25th <- as.numeric(as.character(Univar_Vals$Upper25th))
Univar_Vals$sig <- as.numeric(as.character(Univar_Vals$sig))

# Save
write.csv(Univar_Vals, "DLNM_MetaAnal_Precip_1.36Lag.csv")

