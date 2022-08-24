###################################################################################################
## Main code script to estimate effect modification of the relationship between precipitation and 
##  coccidioidomycosis incidence by prior conditions
##
##  The method uses a two-stage meta-analytic approach with distributed lag non-linear models
##  and code is similar to the code in this github, "DLNM_MetaAnal_1.36Lag_Rainfall.R", but
##  examines effect modification by stratifying the dataset by prior conditions.
##
##  The code is adapted from Gasparrini, et al., Stat Med: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3546395/.
# 
## Written by Jennifer R Head
## Last updated: August 24, 2022
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

# STEP 2. STRATIFY THE DATASET BY PRIOR CONDITIONS BASED ON RAINFALL DEVIATION FROM MEAN

D1.1 <- dataFrame %>% filter(lag21msRainDev < 0) # 2 winters ago drier than average
D1.2 <- dataFrame %>% filter(lag21msRainDev >= 0) # 2 winters ago wetter than average

D2.1 <- dataFrame %>% filter(lag33msRainDev < 0) # 3 winters ago drier than average
D2.2 <- dataFrame %>% filter(lag33msRainDev >= 0) # 3 winters ago wetter than average

# Note, in analysis, also examined when rainfall deviation from mean was < 0 for winters both 2 AND 3 years ago

data.list <- list(D1.1, D1.2, D2.1, D2.2)

## STEP 3. INITIALIZE STORAGE VECTORS

# storage vector where we will store the 25th and the 75th percentile of the main 
# exposure -- here precipitation at 9 mo, across 4 stratified datasets. 25 and 75 allow 
# calculation of IRR associated with an IQR increase (from the 25th to 75th percentile)
per25 <- per75 <- rep(NA, 4)

# storage vector to store the IRR associated with an IQR increase in the exposure
# of interest -- here rainfall at 9 mos -- across 4 stratified datasets. 
# IRR_75 is the point estimate and the low and high stores the upper and lower 95% CI's
IRR_75 <- lowIRR_75 <- highIRR_75<- rep(NA, 4)

## STEP 3. LOOP THROUGH EACH STRATIFIED DATASET, 
##         RUN FIRST AND SECOND STAGE MODEL, 
##         EXTRACT IRR

index <- 9 # Here, the exposure of interest is just rainfall lagged 9 months -- this correspond to winter

for (z in 1:4){ # loop through each of the 4 stratified dataset;
                # note, for analysis, also examined data for two consecutive drier than average winters
  
  dFrame <- data.list[[z]] # loop through the different stratified frames

  # 1. EXTRACT THE EXPOSURE VARIABLE OF INTEREST, THE 25TH and 75TH %ILE, AND OTHER MODEL TERMS
  
  # Extract the exposure variable of interest -- total monthly rain lagged by 9 months
  exposure <- dFrame[,paste0("lag", index, "msTRain")]
  dFrame <- cbind(dFrame, exposure)
  colnames(dFrame)[ncol(dFrame)] <- "exposure"
  
  # Extract the 25th and the 75th percentile of this variable from the data and store it
  iqr25 <- unname(quantile(CasesMo$lag9msTRain, prob = c(0.25), na.rm = T))
  per25[z] <- iqr25
  iqr75 <- unname(quantile(CasesMo$lag9msTRain, prob = c(0.75), na.rm = T))
  per75[z] <- iqr75
  
  # Define the other variables that we want to control for in the model: here mean monthly temps lagged by 3-33 months, and % sand
  dFrame$ME1 <- dFrame$lag3msTmean
  dFrame$ME2 <- dFrame$lag9msTmean
  dFrame$ME3 <- dFrame$lag15msTmean
  dFrame$ME4 <- dFrame$lag21msTmean
  dFrame$ME5 <- dFrame$lag27msTmean
  dFrame$ME6 <- dFrame$lag33msTmean
  dFrame$ME7 <- dFrame$sand_mean
  NumME <- 7
  
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
  rm(dat_list)
  dat_list<-split(dFrame, dFrame$county) #make a list of the dFrame with selected counties
  
  #  RANGES (FOR VARIABLE OF INTEREST -- CAN BE CHANGED BY THE USER)
  #ranges <- t(sapply(dat_list,function(x) range(x$exposure,na.rm=T)))
  
  # DEFINE THE AVERAGE RANGE, CENTERING POINT, DEGREE AND TYPE OF THE SPLINE
  # (THESE PARAMETERS CAN BE CHANGED BY THE USER FOR ADDITIONAL ANALYSES)
  cen <- iqr25 # centered at the 25th percentile
  #(bound <- c(min(ranges[,1]), max(ranges[,2]))) # colMeans(ranges) ##ANOTHER WAY TO DEFINE RANGE
  bound <- c(0,150) # manually set bound. See above line for other options
  type <- "bs"
  degree <- 2

  # 4. DEFINE KNOTS AS ABSOLUTE VALS
  #UPDATE THESE VALS ACCORDING TO VARIABLE CONSIDERED
  allvals <- allknots <- seq(from = bound[1] + 1, to = bound[2]-1, 
                             by = (bound[2]-1 - (bound[1] + 1))/10)
  
  ####################################################################
  # 5. DETERMINE HOW MANY KNOTS GIVES THE LOWEST QAIC AND WHERE
  #    RUN THE MODEL FOR EACH COUNTY, FOR EACH COMBINATION
  
  # COMBINATIONS OF KNOTS FOR INCREASING DF, AT MOST 3 KNOTS
  #ranges #observe if error
  comb <- list(c(1),c(2),c(3),c(4),c(5),c(6),c(7),c(8),c(9),c(10),c(11))
  
  # BUILT THE MATRIX OF Q-AIC VALUES FOR EACH CITY/COMBINATION
  qaicmat <- matrix(0,m,length(comb),dimnames=list(counties,NULL))
  
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
  # 7. PERFORM MULTIVARIATE META-ANALYSIS
  ####################################################################

  # META-ANALYSIS
  mv <- mvmeta(ymat,Slist,method="reml",bscov = "fixed") #would prefer unstr for bscov...might need fixed??
  summary(mv)
  

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
  a<-cp$matRRfit # vector of IRRs at various values of the exposure
  TotalRain <- rownames(a)
  dfa <- data.frame(cbind(a, TotalRain))
  dfa[,1] <- as.numeric(as.character(dfa[,1]))
  dfa[,2] <- as.numeric(as.character(dfa[,2]))
  
  low <- cp$matRRlow # vector of the lower 95% CI's of the IRRs
  dflow <- data.frame(cbind(low, TotalRain))
  high <- cp$matRRhigh # vector of the upper 95% CI's of the IRRs
  dfhigh <- data.frame(cbind(high, TotalRain))
  
  # extract the IRR's, lower and upper 95% CI on the IRR at the 75th percentile. Recall, IRR=1 at 25th percentile due to centering
  IRR_75[z] <- as.numeric(as.character(dfa[which(dfa$TotalRain == round(iqr75,1)),]$lag0))
  lowIRR_75[z] <- as.numeric(as.character(dflow[which(dflow$TotalRain == round(iqr75,1)),]$lag0))
  highIRR_75[z] <- as.numeric(as.character(dfhigh[which(dfhigh$TotalRain == round(iqr75,1)),]$lag0))
  
  ## OPTIONAL COADE TO PLOT
  plot(cp,"overall",col=1,lwd=2,ylab="IRR",ylim=c(0,max(dfa[,1])+1),xlim=c(bound[1],bound[2]),
       xlab=paste("Total Rain (mm), lagged", index, "months"))
  points(iqr25,1,pch=19,cex=1)
  points(round(iqr75,1),cp$allRRfit[as.character(round(iqr75,1))],pch=19,cex=1)
  abline(h = cp$allRRfit[as.character(round(iqr75,1))], lty = "dotted")
  
}

## STEP 4. SAVE THE MODEL OUTPUT IN A TABLE, PLOT

res.df <- data.frame(cbind("IRR" = IRR_75,
                           "LCI" = lowIRR_75,
                           "UCI" = highIRR_75))

res.df$dev <- c("Q1: 0-25%ile", "Q2: 26-50%ile", "Q3: 51-75%ile", "Q4: 76-100%ile")


ggplot(res.df) +
  geom_point(aes (x = dev, y = IRR, color = dev), position = position_dodge(width = 0.2)) +
  geom_errorbar(aes (x = dev, ymin = LCI, ymax = UCI, color = dev), position = position_dodge(width = 0.2), width = 0.1) +
  scale_color_manual("Antecedent winters \n compared to average",
                    breaks = c("Drier", "Wetter"),
                    values = c("thistle3", "darkcyan"))
                

