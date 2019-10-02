
#########################################################################
## EILecture_POLI592.R
## Created by Yuki Atsusaka
## Since: April 16, 2019
## Last updated: April 18, 2019
## Aim: to build and estimate Binomial-Beta models
## THIS CODE WAS PREPARED FOR "Introduction to Ecological Inference" 
## IN POLI 592 (ADVANVED POLITICAL METHODOLOGY) (GRADUATE COURSE)
#########################################################################


# READING DATA #########################################################################################################
# THIS IS FOR MYSELF
# library(dplyr)
# library(tidyverse)
# library(readstata13)
# library(gmodels)
#  dat <- read.dta13("NC_Voting_history_2016_regs.dta")
#  dat2 <- dat %>% 
#          mutate(Vote = ifelse(vote_decision=="Did Not Vote", 0, 1)) %>%
#          mutate(NonWhite = ifelse(white==0, 1,0)) %>%
#          dplyr::select(Vote, NonWhite, county_fips) %>%
#          drop_na() 
# write.csv(dat2, "NCCounty_ind.csv", row.names=F) # SAVE THE AGGREGATE DATA 
########################################################################################################################

# dat2 <- read.table("NCCounty_ind.csv", sep=",", header=T) # IF YOU WANT TO SEE THIS
# # CHECK THE CONTINGENCY TABLE (GROUND TRUTH)
# CrossTable(dat2$NonWhite, dat2$Vote, prop.r=TRUE, prop.c=FALSE, prop.t=FALSE, prop.chisq=FALSE,
#            format="SPSS", dnn=c("Non White", "Vote"))

# # CONSTRUCT COUNTY LEVEL DATASET
#  dat_agg <- dat2 %>%
#             group_by(county_fips) %>%
#             mutate(N_Vote = sum(Vote)) %>%
#             mutate(N_NonWhite = sum(NonWhite)) %>%
#             mutate(N_total = length(Vote)) %>%
#             mutate(Xi = N_NonWhite/N_total) %>%
#             mutate(Ti = N_Vote/N_total) %>%
#             distinct(county_fips, N_Vote, N_NonWhite, N_total, Xi, Ti)
# 
#  head(dat_agg)
# write.csv(dat_agg, "NCCounty.csv", row.names=F) # SAVE THE AGGREGATE DATA
####################################################################################################

 

####################################################################################################
# LECTURE CODE STARTS HERE
####################################################################################################

rm(list=ls())
set.seed(12262017)
library(pscl)
library(coda)
library(R2jags)
library(mcmcplots)

# HERE IS THE GROUND TRUTH THAT WE WANT TO ESTIMATE BELOW

#              | Vote 
#    Non White |        0  |        1  | Row Total | 
# -------------|-----------|-----------|-----------|
#            0 |  1602236  |  2830021  |  4432257  | 
#              |   36.149% |   63.851% |   67.857% | 
# -------------|-----------|-----------|-----------|
#            1 |   936989  |  1162482  |  2099471  | 
#              |   44.630% |   55.370% |   32.143% | 
# -------------|-----------|-----------|-----------|
# Column Total |  2539225  |  3992503  |  6531728  | 
# -------------|-----------|-----------|-----------|


dat_agg <- read.table("NCCounty.csv", sep=",", header=T)

# SCATTER PLOT AND TOMOGRAPHY PLOT-----------------------------------------------------------------#
cons  = (dat_agg$Ti/(1-dat_agg$Xi))
slope = (dat_agg$Xi/(1-dat_agg$Xi))

j = 4   # SELECTED COUNTY
k = 88  # SELECTED HOMOGENEOUS WHITE COUNTY

   
par(mfrow=c(1,2)) 
 plot(dat_agg$Ti ~ dat_agg$Xi, xlim=c(0,1), ylim=c(0,1), xlab="Xi",ylab="Ti", pch=16, main="Scatter Plot") 
 points(dat_agg$Ti[j] ~ dat_agg$Xi[j], col="firebrick", pch=16, cex=2)
 points(dat_agg$Ti[k] ~ dat_agg$Xi[k], col="green", pch=16, cex=2)
 
 plot(dat_agg$Ti ~ dat_agg$Xi, xlim=c(0,1), ylim=c(0,1), xlab=expression(beta[i]^b),ylab=expression(beta[i]^w),
      col="white", pch=16, main="Tomography Plot") 
 for(i in 1:100){ abline(a=cons[i], b=slope[i]) }  
 abline(a=cons[j], b=slope[j], lwd=8, col="firebrick")  
 abline(a=cons[k], b=slope[k], lwd=8, col="green")  
#--------------------------------------------------------------------------------------------------# 
 
# CREATE VARIABLES FOR JAGS MODELS 
 N_obs = dim(dat_agg)[1]
 N_Vote = dat_agg$N_Vote
 N_total = dat_agg$N_total
 N_b = dat_agg$N_NonWhite
 N_w = N_total - N_b
 N_b_state = sum(N_b)
 N_w_state = sum(N_w)
 X = dat_agg$Xi        # PROP OF BLACKS
 Z = X

##############################################################################
####### WRITE DOWN A BUGS FILE WITHIN R ######################################
##############################################################################
# BINOMIAL-BETA EI MODEL WITH NO COVARIATE
cat("model{
    
# (1) BINOMIAL LIKELIHOOD
    for(i in 1:N_obs){ 
    N_Vote[i] ~ dbinom(phi[i], N_total[i])

    phi[i] <- beta_b[i]*X[i] + beta_w[i]*(1 - X[i])

    temp_b[i] <- N_b[i]*beta_b[i]
    temp_w[i] <- N_w[i]*beta_w[i]
    }

    B_b <- sum(temp_b)/N_b_state  # STATEWIDE-PARAMETER
    B_w <- sum(temp_w)/N_w_state  # STATEWIDE-PARAMETER
    

# (2) BETA PRIOR FOR beta_b AND beta_w (YOU CAN COMBINE THE LOOPS)
    for(i in 1:N_obs){
    beta_b[i] ~ dbeta(c_b, d_b)
    beta_w[i] ~ dbeta(c_w, d_w)
    } 

# beta_w[i] <- (T[i]/(1-X[i])) + (X[i]/(1-X[i]))*beta_b[i]    

# (3) EXPONENTIAL HYPER PRIOR
    c_b ~ dexp(0.5)
    d_b ~ dexp(0.5)
    c_w ~ dexp(0.5)
    d_w ~ dexp(0.5)

    }", fill=TRUE, file="BinomBeta.bug") # OUTPUT BUGS FILE
###############################################################################  
###############################################################################  
# PROVIDE FIXED INFORMATION (i.e., DATA IN JAGS) 
datalist <- list("N_obs", "N_Vote", "N_total", "X",                      # DATA TO USE
                 "N_b", "N_w", "N_b_state", "N_w_state")
inits <- function() { list(beta_b=runif(100,min=0,max=1),                # INITIAL VALUES OF CHAINS
                           beta_w=runif(100,min=0,max=1),                # THAT ARE PLAUSIBLE FROM 
                           c_b=rexp(1,rate=0.5), d_b=rexp(1,rate=0.5),   # THE SPECIFIED PDF/PMF
                           c_w=rexp(1,rate=0.5), d_w=rexp(1,rate=0.5))}  # AS PRIOR/HYPERPRIOR

parameters <- c("beta_b", "beta_w", "B_b", "B_w")                        # PARAMETERS TO SAVE

model.sim<-jags(datalist, inits, parameters, model.file="BinomBeta.bug", # ESTIMATE THE MODEL
                n.chains=3, n.iter=80000, n.burnin=8000, n.thin=1)
############################################################################### 


# NOW, LET'S CHECK THE RESULTS
model.sim                                  # ALL RESULTS
model.sim$BUGSoutput$summary[1:2,]         # ONLY STATE-WIDE RESULTS
denplot(model.sim, parms=c("B_b", "B_w"),  # PLOT THE POSTEIOR DENSITIRES
        main=c("Statewide Prop of Non-Whites Who Turned Out",
               "Statewide Prop of Whites Who Turned Out"))

denplot(model.sim, parm="phi") # WE CAN ALSO SEE OTHER PARAMETERS (IF INCLUDED IN parameters)


# CONVERGENCE DIAGNOSTICS (ALWAYS CHECK CONVERGENCE NUMERICALLY)
gelman.diag(as.mcmc(model.sim))     # IF 1.00: OK
geweke.diag(model.sim)              # IF NO SIGNIFICANT DIFFERENCE: OK
heidel.diag(as.mcmc(model.sim))     # IF ALL PASS: OK 
raftery.diag(as.mcmc(model.sim))    # IF MINIMUM LENGTH: NOT OK
mcmcplot(model.sim)                 # VISUALLY (THIS WILL POP UP AN html FILE)



# PLOT THE MEAN POINTS ON THE TOMORAPHY PLOT--------------------------------------------------------------#
mean_b <- model.sim$BUGSoutput$summary[3:102,1]    # COUNTY LEVEL MEAN
mean_w <- model.sim$BUGSoutput$summary[103:202,1]  # COUNTY LEVEL MEAN 

par(mfrow=c(1,2)) 
plot(dat_agg$Ti ~ dat_agg$Xi, xlim=c(0,1), ylim=c(0,1), xlab="Xi",ylab="Ti", pch=16, main="Scatter Plot") 
points(dat_agg$Ti[j] ~ dat_agg$Xi[j], col="firebrick", pch=16, cex=2)
points(dat_agg$Ti[k] ~ dat_agg$Xi[k], col="green", pch=16, cex=2)

plot(dat_agg$Ti ~ dat_agg$Xi, xlim=c(0,1), ylim=c(0,1), xlab=expression(beta[i]^b),ylab=expression(beta[i]^w),
     col="white", pch=16, main="Tomography Plot") 
for(i in 1:100){ abline(a=cons[i], b=slope[i]) }  
abline(a=cons[j], b=slope[j], lwd=8, col="firebrick")  
abline(a=cons[k], b=slope[k], lwd=8, col="green")  
#mean_w_est  = (dat_agg$Ti/(1-dat_agg$Xi)) + (dat_agg$Xi/(1-dat_agg$Xi))*mean_b
points(mean_w ~ mean_b, col="navy", pch=16, cex=1)  
# DO POSTERIOR MEANS MAP ONTO THE TOMOGRAPHY LINES? WHY SO?
#---------------------------------------------------------------------------------------------------------#


##############################################################################
# LET'S INCLUDE COVARIATES NEXT
##############################################################################
# BINOMIAL-BETA EI MODEL WITH COVARIATES
cat("model{
    
# (1) BINOMIAL LIKELIHOOD (SAME AS ABOVE)
    for(i in 1:N_obs){ 
    N_Vote[i] ~ dbinom(phi[i], N_total[i])

    phi[i] <- beta_b[i]*X[i] + beta_w[i]*(1 - X[i])

    temp_b[i] <- N_b[i]*beta_b[i]
    temp_w[i] <- N_w[i]*beta_w[i]
    
    }
    
    B_b <- sum(temp_b)/N_b_state  # STATEWIDE-PARAMETER
    B_w <- sum(temp_w)/N_w_state  # STATEWIDE-PARAMETER
    
# (2) BETA PRIOR FOR beta_b AND beta_w
    for(i in 1:N_obs){
    beta_b[i] ~ dbeta(d_bX[i], d_b)
    d_bX[i] <- d_b * exp(alpha + beta*Z[i])
  } 

    for(i in 1:N_obs){
    beta_w[i] ~ dbeta(d_wX[i], d_w)
    d_wX[i] <- d_w * exp(gamma + delta*Z[i])
    } 

# (3) EXPONENTIAL AND FLAT HYPER PRIOR
    d_b ~ dexp(0.5)
    d_w ~ dexp(0.5)
    alpha ~ dnorm(0, 1e-3)   
    beta  ~ dnorm(0, 1e-3)   
    gamma ~ dnorm(0, 1e-3)   
    delta ~ dnorm(0, 1e-3)   
   
    }", fill=TRUE, file="BinomBeta_Cov.bug") # OUTPUT BUGS FILE
###############################################################################  
###############################################################################  
# PROVIDE FIXED INFORMATION (i.e., DATA IN JAGS) 
datalist <- list("N_obs", "N_Vote", "N_total", "X", "N_b", "N_w", "N_b_state", "N_w_state", "Z")
inits <- function() { list(beta_b=runif(100,min=0,max=1), 
                           beta_w=runif(100,min=0,max=1), 
                           d_b=rexp(1,rate=0.5), d_w=rexp(1,rate=0.5), 
                           alpha=rnorm(1,mean=0,sd=1), beta=rnorm(1,mean=0,sd=1),
                           gamma=rnorm(1,mean=0,sd=1), delta=rnorm(1,mean=0,sd=1))}
parameters <- c("beta_b", "beta_w", "B_b", "B_w", "beta", "delta")

model.sim2<-jags(datalist, inits, parameters, model.file="BinomBeta_Cov.bug",
                n.chains=3, n.iter=80000, n.burnin=8000, n.thin=1)
############################################################################### 

# CHECK THE RESULTS
model.sim2 # ALL RESULTS
model.sim2$BUGSoutput$summary[1:2,]         # ONLY STATE-WIDE RESULTS

denplot(model.sim2, parms=c("B_b", "B_w"),  # PLOT THE POSTEIOR DENSITIRES
        main=c("Statewide Prop of Non-Whites Who Turned Out",
               "Statewide Prop of Whites Who Turned Out"))

denplot(model.sim2, parms=c("beta", "delta"),  # PLOT THE POSTEIOR DENSITIRES FOR COVARIATE
        main=c("Coefficient on the Covariate (beta)",
               "Coefficient on the Covariate (delta)"))

mcmcplot(model.sim2)

# CONVERGENCE DIAGNOSTICS
gelman.diag(as.mcmc(model.sim2))     # IF 1.00: OK
geweke.diag(model.sim2)              # IF NO SIGNIFICANT DIFFERENCE: OK
heidel.diag(as.mcmc(model.sim2))     # IF ALL PASS: OK 
raftery.diag(as.mcmc(model.sim2))    # IF MINIMUM LENGTH: NOT OK



# PLOT THE MEAN POINTS ON THE TOMORAPHY PLOT
mean_b <- model.sim2$BUGSoutput$summary[3:102,1]  # COUNTY LEVEL MEAN
mean_w <- model.sim2$BUGSoutput$summary[103:202,1]  # COUNTY LEVEL MEAN 

par(mfrow=c(1,2)) 
plot(dat_agg$Ti ~ dat_agg$Xi, xlim=c(0,1), ylim=c(0,1), xlab="Xi",ylab="Ti", pch=16, main="Scatter Plot") 
points(dat_agg$Ti[j] ~ dat_agg$Xi[j], col="firebrick", pch=16, cex=2)
points(dat_agg$Ti[k] ~ dat_agg$Xi[k], col="green", pch=16, cex=2)

plot(dat_agg$Ti ~ dat_agg$Xi, xlim=c(0,1), ylim=c(0,1), xlab=expression(beta[i]^b),ylab=expression(beta[i]^w),
     col="white", pch=16, main="Tomography Plot") 
for(i in 1:100){ abline(a=cons[i], b=slope[i]) }  
abline(a=cons[j], b=slope[j], lwd=8, col="firebrick")  
abline(a=cons[k], b=slope[k], lwd=8, col="green")  
#mean_w_est  = (dat_agg$Ti/(1-dat_agg$Xi)) + (dat_agg$Xi/(1-dat_agg$Xi))*mean_b
points(mean_w ~ mean_b, col="navy", pch=16, cex=1)  





##############################################################################  
# END OF THIS R SOURCE FILE
##############################################################################  