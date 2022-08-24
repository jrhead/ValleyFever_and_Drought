###################################################################################################
## Main code script to estimate associations between precipitation, temperature and coccidioidomycosis 
##   incidence, and to understand factors that explain heterogeneity in the relationship
##   It particularly examines heterogeneity for the relationship between winter rainfall and summer temps
##   and incidence.
##
##  The method uses a two-stage meta-analytic approach with distributed lag non-linear models
##  and code is similar to the code in this github, "DLNM_MetaAnal_1.36Lag_Rainfall.R", but
##  uses a multivariate meta-analysis to understand heterogeneity in the exposure-response
##  relationship across counties. This code makes Figure 3 in the paper.
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

###################################################################################################
### FIRST,  WE PREPARE THE DATA
### SECOND, WE RUN THE DLNM AND MULTIVARIATE ANALYSIS FOR RAINFALL
### THIRD,  WE RUN THE NLNM AND MULTIVARIATE ANALYSIS FOR TEMPERATURE
### FOURTH, WE COMPUTE PREDICTION FOR MULTIVARIATE META-REGRESSION MODELS, AND PLOT BOTH
###################################################################################################

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

####################################################################
####  FIRST 1/2 OF SCRIPT:  RECENT WINTER RAINFALL AS EXPOSURE   ###
####################################################################

# STEP 2. DEFINE METAPREDICTORS -- here, they are the median, 25th, 75th %tile of typical winter rain and summer temps by county
submod <- dataFrame %>% group_by(county) %>% summarize(rain25 = quantile(lag9msTRain, probs = 0.25, na.rm = T),
                                                       rain50 = quantile(lag9msTRain, probs = 0.50, na.rm = T),
                                                       rain75 = quantile(lag9msTRain, probs = 0.75, na.rm = T),
                                                       temp25 = quantile(lag3msTmean, probs = 0.25, na.rm = T),
                                                       temp50 = quantile(lag3msTmean, probs = 0.50, na.rm = T),
                                                       temp75 = quantile(lag3msTmean, probs = 0.75, na.rm = T))

# STEP 3. DEFINE THE MODEL FORMULA

# 3A. Define the main exposure of interest - here rainfall lagged 9 months (winter)
index <- 9
dataFrame$exposure <- pull(dataFrame, paste0("lag", index, "msTRain"))

# 3B. Define other confounders in the model -- here, lagged temperatures (lagged by 6 mo intervals)
dataFrame$ME1 <- dataFrame$lag3msTmean
dataFrame$ME2 <- dataFrame$lag9msTmean
dataFrame$ME3 <- dataFrame$lag15msTmean
dataFrame$ME4 <- dataFrame$lag21msTmean
dataFrame$ME5 <- dataFrame$lag27msTmean
dataFrame$ME6 <- dataFrame$lag33msTmean
dataFrame$ME7 <- dataFrame$sand_mean
NumME <- 7

# 3B. Initialize model formula -- add to it the confounders ME1-7 and other lags on rainfall

formula <- "N_all ~ splinevar + log(offset(CTpopulation)) + ns(year) + ns(ME1, knots = 1)"

for (i in 2:NumME){
  formula <- paste(formula, " + ns(ME", i, ", knots=1)", sep = "")
}

if(index == 1){
  other_vars <- c(3,6,9,12,15,18,21,24,27,30,33)
} else{
  other_vars <- (c(3,6,9,12,15,18,21,24,27,30,33))[-round(index/3,0)]}

for (i in unique(other_vars)){
  formula <- paste(formula, " + ns(lag", i, "msTRain , knots=1)", sep = "")
}

formula <- as.formula(formula)

# 3D. define the 25th percentile of the exposure
iqr25 <- unname(quantile(dataFrame$exposure, prob = c(0.25), na.rm = T))

# STEP 4. DEFINE COUNTIES AND CREATE LIST OF DF

dat_list<-split(dataFrame, dataFrame$county) #make a list of the dataframe with selected counties

#  RANGES (FOR VARIABLE OF INTEREST -- CAN BE CHANGED BY THE USER)
ranges <- t(sapply(dat_list,function(x) quantile(x$exposure, probs = c(0.0, 0.90),na.rm=T)))

# DEFINE THE AVERAGE RANGE, CENTERING POINT, DEGREE AND TYPE OF THE SPLINE
# (THESE PARAMETERS CAN BE CHANGED BY THE USER FOR ADDITIONAL ANALYSES)
cen <- iqr25 # ceneter at the 25th percentile
(bound <- c(min(ranges[,1]), max(ranges[,2])))
degree <- 2
type <- "bs"


##  STEP 5. DETERMINE HOW MANY KNOTS GIVES THE LOWEST QAIC AND WHERE

## 5A. Define knots as absolute values
#UPDATE THESE VALS ACCORDING TO VARIABLE CONSIDERED
allvals <- allknots <- seq(from = min(ranges[,1]) + 1, to = max(ranges[,2])-1, 
                           by = ((max(ranges[,2])-1) - (min(ranges[,1]) + 1))/10)

# 5B. COMBINATIONS OF KNOTS FOR INCREASING DF, AT MOST 3 KNOTS
#ranges #observe if error
comb <- list(c(1),c(2),c(3),c(4),c(5),c(6),c(7),c(8),c(9),c(10),c(11))

# 5C. BUILD THE MATRIX OF Q-AIC VALUES FOR EACH CITY/COMBINATION
qaicmat <- matrix(0,m,length(comb),dimnames=list(counties,NULL))


# 5D. RUN THE MODEL FOR EACH COUNTY, FOR EACH COMBINATION

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
# 5E. COMPARING Q_AIC BY KNOT LOCATION AND DF

# SUM THE Q-AIC ACROSS COUNTIES
qaictot <- colSums(qaicmat, na.rm = T)

# 5F. CHOSEN DF, PERCENTILES AND KNOTS AT MINIMUM Q-AIC
df <- length(comb[[which.min(qaictot)]]) + 1 + (type=="bs")*(degree-1)
perc <- allvals[comb[[which.min(qaictot)]]] #allvalls OR allperc
(knots <- allknots[comb[[which.min(qaictot)]]])


####################################################################
# STEP 6. RUN THE FIRST STAGE MODEL
####################################################################

# 6A. BUILD OBJECTS WHERE RESULTS WILL BE STORED:
#   ymat IS THE MATRIX FOR THE OUTCOME PARAMETERS
#   Slist IS THE LIST WITH (CO)VARIANCE MATRICES
ymat <- matrix(NA,m,df-1,dimnames=list(counties,paste("spl",seq(df-1),sep="")))
Slist <- vector("list",m)
names(Slist) <- counties

####################################################################
# 6B. RUN THE FIRST-STAGE ANALYSIS

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

# STEP 7. UNIVARIATE META-ANALYSIS
mv <- mvmeta(ymat,Slist,method="reml",bscov = "fixed") #would prefer unstr for bscov...might need fixed??
summary(mv)


####################################################################
# STEP 8. CREATE BASIS FOR PREDICTION
####################################################################

# BASIS USED TO PREDICT VARIABLE, EQUAL TO THAT USED FOR ESTIMATION
#   NOTE: INTERNAL AND BOUNDARY KNOTS PLACED AT SAME VALUES AS IN ESTIMATION
#   NOTE: SPLINE CENTERED ON THE SAME VALUE AS IN ESTIMATION
predictorvar <- seq(bound[1],bound[2],length=30)
bpredictorvar <- onebasis(predictorvar,#type=type,degree=degree,
                          knots=knots,
                          bound=bound,cen=cen)

####################################################################
# STEP 9. PREDICTION FROM MODELS
####################################################################

# USE OF crosspred TO PREDICT THE EFFECTS FOR THE CHOSEN VALUES

# PREDICTION FROM SIMPLE META-ANALYSES WITH NO PREDICTORS
cp <- crosspred(bpredictorvar,coef=coef(mv),vcov=vcov(mv),model.link="log",by=0.1)
#### PLOT INDIVIDUAL LINES
for(i in seq(m)) {
  # WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
  suppressWarnings(
    lines(crosspred(bpredictorvar,coef=ymat[i,],vcov=Slist[[i]],model.link="log"),
          col=grey(0.8),lty=5)
  )
}
points(round(knots,1),cprel$allRRfit[as.character(knots)],pch=19,cex=0.6)
abline(h=1)

####################################################################
# STEP 10. PERFORM MULTIVARIATE META-ANALYSIS
####################################################################

mvprecip <- mvmeta(ymat~rain50,Slist,submod,method="reml",bscov = "fixed")
summary(mvprecip)


########################################################################
####  SECOND 1/2 OF SCRIPT:  RECENT SUMMER TEMPERATURE AS EXPOSURE   ###
########################################################################

# STEP 3. DEFINE THE MODEL FORMULA

# 3A. Define the main exposure of interest - here mean temperature lagged 3 months (summer)
index <- 3
dataFrame$exposure <- pull(dataFrame, paste0("lag", index, "msTmean"))

# 3B. Define other confounders in the model -- lagged rainfall
dataFrame$ME1 <- dataFrame$lag3msTRain
dataFrame$ME2 <- dataFrame$lag6msTRain
dataFrame$ME3 <- dataFrame$lag9msTRain
dataFrame$ME4 <- dataFrame$lag12msTRain
dataFrame$ME5 <- dataFrame$lag15msTRain
dataFrame$ME6 <- dataFrame$lag18msTRain
dataFrame$ME7 <- dataFrame$lag21msTRain
dataFrame$ME8 <- dataFrame$lag24msTRain
dataFrame$ME9 <- dataFrame$lag27msTRain
dataFrame$ME10 <- dataFrame$lag33msTRain
dataFrame$ME11 <- dataFrame$sand_mean
NumME <- 11

# 3C. Compile the model formula, adding in variables on rainfall and other lags for temperature
formula <- "N_all ~ splinevar + log(offset(CTpopulation)) + ns(year) + ns(ME1, knots = 1)"

for (i in 2:NumME){
  formula <- paste(formula, " + ns(ME", i, ", knots=1)", sep = "")
}

if(index == 1 | index == 2){
  other_vars <- c(3,9,15,21,27,33)
} else{
  other_vars <- (c(3,9,15,21,27,33))[-round(index/6,0)]}

for (i in unique(other_vars)){
  formula <- paste(formula, " + ns(lag", i, "msTmean , knots=1)", sep = "")
}

formula <- as.formula(formula)

# 3D. define the 25th percentile of the exposure
iqr25 <- unname(quantile(dataFrame$exposure, prob = c(0.25), na.rm = T))

# STEP 4: DEFINE COUNTIES AND CREATE LIST OF DF

dat_list<-split(dataFrame, dataFrame$county) #make a list of the dataframe with selected counties

#  RANGES (FOR VARIABLE OF INTEREST -- CAN BE CHANGED BY THE USER)
ranges <- t(sapply(dat_list,function(x) range(x$exposure,na.rm=T)))

# DEFINE THE AVERAGE RANGE, CENTERING POINT, DEGREE AND TYPE OF THE SPLINE
# (THESE PARAMETERS CAN BE CHANGED BY THE USER FOR ADDITIONAL ANALYSES)
cen <- iqr25 # center at the 25th percentile
(bound <- c(min(ranges[,1]), max(ranges[,2])))
degree <- 2
type <- "bs"

##  STEP 5. DETERMINE HOW MANY KNOTS GIVES THE LOWEST QAIC AND WHERE

## 5A. Define knots as absolute values
#UPDATE THESE VALS ACCORDING TO VARIABLE CONSIDERED
allvals <- allknots <- seq(from = min(ranges[,1]) + 1, to = max(ranges[,2])-1, 
                           by = ((max(ranges[,2])-1) - (min(ranges[,1]) + 1))/10)

# 5B. COMBINATIONS OF KNOTS FOR INCREASING DF, AT MOST 3 KNOTS
#ranges #observe if error
comb <- list(c(1),c(2),c(3),c(4),c(5),c(6),c(7),c(8),c(9),c(10),c(11))

# 5C. BUILD THE MATRIX OF Q-AIC VALUES FOR EACH CITY/COMBINATION
qaicmat <- matrix(0,m,length(comb),dimnames=list(counties,NULL))

# 5D. DETERMINE HOW MANY KNOTS GIVES THE LOWEST QAIC AND WHERE
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

# 5E. COMPARE Q_AIC BY KNOT LOCATION AND DF

# SUM THE Q-AIC ACROSS COUNTIES
qaictot <- colSums(qaicmat, na.rm = T)

# 5F. CHOSEN DF, PERCENTILES AND KNOTS AT MINIMUM Q-AIC
df <- length(comb[[which.min(qaictot)]]) + 1 + (type=="bs")*(degree-1)
perc <- allvals[comb[[which.min(qaictot)]]] #allvalls OR allperc
(knots <- allknots[comb[[which.min(qaictot)]]])


####################################################################
# STEP 6. RUN THE FIRST STAGE MODEL
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

# STEP 7. UNIVARIATE META-ANALYSIS
mv <- mvmeta(ymat,Slist,method="reml",bscov = "fixed") #would prefer unstr for bscov...might need fixed??
summary(mv)

####################################################################
# STEP 8. CREATE BASIS FOR PREDICTION
####################################################################

# BASIS USED TO PREDICT VARIABLE, EQUAL TO THAT USED FOR ESTIMATION
#   NOTE: INTERNAL AND BOUNDARY KNOTS PLACED AT SAME VALUES AS IN ESTIMATION
#   NOTE: SPLINE CENTERED ON THE SAME VALUE AS IN ESTIMATION
predictorvar <- seq(bound[1],bound[2],length=30)
bpredictorvart <- onebasis(predictorvar,#type=type,degree=degree,
                          knots=knots,
                          bound=bound,cen=cen)

####################################################################
# STEP 9. PREDICTION FROM MODELS
####################################################################

# USE OF crosspred TO PREDICT THE EFFECTS FOR THE CHOSEN VALUES

# PREDICTION FROM SIMPLE META-ANALYSES WITH NO PREDICTORS
cp <- crosspred(bpredictorvart,coef=coef(mv),vcov=vcov(mv),model.link="log",by=0.1)

####################################################################
# STEP 10. PERFORM MULTIVARIATE META-ANALYSIS
####################################################################

mvtemp <- mvmeta(ymat~temp50,Slist,submod,method="reml",bscov = "fixed")
summary(mvtemp)

##############################################################################################################
#  FINAL COMPONENT OF THE SCRIPT                                  
#
# COMPUTE PREDICTION FOR MULTIVARIATE META-REGRESSION MODELS
#   1ST STEP: PREDICT THE OUTCOME PARAMETERS FOR SPECIFIC VALUES OF META-PREDICTOR
#   2ND STEP: CREATE COLOR PALETTE
#   3RD STEP: PREDICT THE RELATIONSHIP AT CHOSEN VALUES GIVEN THE PARAMETERS, PLOTTING AS YOU GO
##############################################################################################################

# STEP 1: PREDICT THE OUTCOME PARAMETERS FOR SPECIFIC VALUES OF META-PREDICTOR

predrain <- predict(mvprecip,data.frame(rain50=submod[,"rain50"]),vcov=T) # RAINFALL
predtemp <- predict(mvtemp,data.frame(temp50=submod[,"temp50"]),vcov=T) # TEMPERATURE

# STEP 2: CREATE COLOR PALETTES
library(paletteer)
library(gameofthrones)

nColor <- m
colorsr <- paletteer_c("pals::coolwarm", n = nColor, direction = -1) #rain
colorst <- paletteer_c("gameofthrones::lannister", n = nColor, direction = -1) #temp

# Transform the numeric variable in bins - rain
rank <- as.factor(as.numeric( cut(submod$rain50, nColor)))
colorsr <- colorsr[rank]

# Transform the numeric variable in bins - temp
rank <- as.factor(as.numeric( cut(submod$temp50, nColor)))
colorst <- colorst[rank]


# STEP 3. PREDICT THE RELATIONSHIP AT CHOSEN VALUES GIVEN THE PARAMETERS, PLOTTING AS YOU GO

# 3A. Open pdf file for saving

pdf("Manuscript\\Figures\\Fig3_MetaPredictors.pdf", width = 15, height = 6)

par(mfrow=c(1,2), mar=c(4,6,4,3)+.1)

#3B. Make the plot for rain

# plot the first line
s_i <- 1 
plot(crosspred(bpredictorvar,coef=predrain[[s_i]]$fit,
               vcov=predrain[[s_i]]$vcov,model.link="log",by=0.1), 
     lwd = 2,
     col = colorsr[s_i],
     ci = "area",
     ci.arg=list(density = 10,col=colorsr[s_i]),
     ylab = "IRR",
     xlab = "Winter rainfall (mm)",
     cex.lab=2, cex.axis=2,
     ylim = c(0.6,2.5))

abline(h=1, col = "white") #remove straight black line that goes all the way across

# plot the remaining lines
for (i in 2:14){#hot_counties){
  lines(crosspred(bpredictorvar,coef=predrain[[i]]$fit,
                  vcov=predrain[[i]]$vcov,model.link="log",by=0.1),
        lwd = 2,
        col = colorsr[i],
        ci = "area",
        ci.arg=list(density=10,col=colorsr[i]))
}

lines(x = c(0,145), y = c(1,1), col = "black") #remove straight black line that goes all the way across

# add legend
colr <- paletteer_c("pals::coolwarm", n = nColor, direction = -1)
legend.col(col = colr, lev = submod$rain50)

#3C. Make the plot for temp
plot(crosspred(bpredictorvart,coef=predtemp[[s_i]]$fit,
               vcov=predtemp[[s_i]]$vcov,model.link="log",by=0.1), 
     lwd = 2,
     col = colorst[s_i],
     ci = "area",
     ci.arg=list(density = 10,col=colorst[s_i]),
     ylab = "IRR",
     xlab = "Summer temperature (\u00B0C)",
     cex.lab=2, cex.axis=2,
     ylim = c(0,4))

abline(h=1, col = "white") #remove straight black line that goes all the way across

# plot the remaining lines
for (i in 2:14){ #hot_counties){ #2-14
  lines(crosspred(bpredictorvart,coef=predtemp[[i]]$fit,
                  vcov=predtemp[[i]]$vcov,model.link="log",by=0.1),
        lwd = 2,
        col = colorst[i],
        ci = "area",
        ci.arg=list(density=10,col=colorst[i]))
}

lines(x = c(12.43,30.9), y = c(1,1), col = "black") #remove straight black line that goes all the way across

# add legend
colr <- paletteer_c("gameofthrones::lannister", n = nColor, direction = -1)
legend.col(col = colr, lev = submod$temp50)

# close pdf
dev.off()

