

# Phytoplankton forward PSM driver for use with environmental time series model 
############################################################################################  
# Runs Bayesian inversion on calibration data to determine slopes and intercepts for the:
#   1) mu(i) = f(po4, rm) relationship and 
#   2) lith length and coccosphere diameter to mean cell radius transfer functions
#
# Then runs inversion evaluated against proxy data using prior inputs 
#
# Dustin T. Harper
############################################################################################  


# Load libraries
############################################################################################ 
library(rjags)
library(R2jags)
############################################################################################ 


# Multi linear regression model for prior slopes and intercepts to calculate mu(i) 
############################################################################################    
# Read in instantaneous growth rate (mu,i) culture calibration data from Aloisi et al. (2015)
# with calculated mean radius from measured coccosphere using Henderiks and Pagani (2007) transfer functions
cal.df <- read.csv('data/caldata_culture.csv')

# Generate multiple linear regression model mu,i (as a function of [PO4] and mean radius)
igr.model = lm(formula = cal.df$mui ~ cal.df$po4 + cal.df$r, data = cal.df)
igr.model.sum <- summary(igr.model)

#    Load prior coefficients and SEs for mui = f(po4, rm) from multi linear regression 
po4.co.lr <- igr.model.sum$coefficients[2,1]
po4.se.lr <- igr.model.sum$coefficients[2,2]
r.co.lr <- igr.model.sum$coefficients[3,1]
r.se.lr <- igr.model.sum$coefficients[3,2]
y.int.lr <- igr.model.sum$coefficients[1,1]
y.int.se.lr <- igr.model.sum$coefficients[1,2]
############################################################################################


# Determine size-based transfer functions using linear regression and Henderiks and Pagani 07 data
############################################################################################
hp07 <- read.csv('data/HP07_sizedata.csv')
lm.lith <- lm(hp07$lith.size ~ hp07$cell.r, data = hp07)
lith.sum <- summary(lm.lith)
lm.cell <- lm(hp07$sphere.d ~ hp07$cell.r, data = hp07)
cell.sum <- summary(lm.cell)

lith.r.m <- lith.sum$coefficients[2,1]
lith.r.mu <- lith.sum$coefficients[2,2]
lith.r.b <- lith.sum$coefficients[1,1]
lith.r.bu <- lith.sum$coefficients[1,2]

cell.r.m <- cell.sum$coefficients[2,1]
cell.r.mu <- cell.sum$coefficients[2,2]
cell.r.b <- cell.sum$coefficients[1,1]
cell.r.bu <- cell.sum$coefficients[1,2]
############################################################################################


# Load and groom data for coefficient inversion
############################################################################################  

# remove blanks from individual data columns; clean data 
clean.rad <- cal.df[complete.cases(cal.df$radius), which(colnames(cal.df)=='radius')]
clean.po4 <- cal.df[complete.cases(cal.df$po4), which(colnames(cal.df)=='po4')]
clean.mui <- cal.df[complete.cases(cal.df$mui), which(colnames(cal.df)=='mui')]
clean.diamcc <- cal.df[complete.cases(cal.df$diam.cc), which(colnames(cal.df)=='diam.cc')]

# Cull out of range data sets
cal.df <- subset(cal.df, is.na(cal.df$radius) | is.na(cal.df$po4) | cal.df$radius <= 10 & cal.df$po4 <= 2)
nrow.cal.df <- nrow(cal.df)

# Parameters to save in inversion output
parms1 = c("mui", "diam.cocco", "rm", "po4", "lith.m", "lith.b", "diam.m", "diam.b", "y.int", "co.r", "co.po4") 

# data to pass to jags
data.pass.coeff = list("po4.co.lr" = po4.co.lr, "po4.se.lr" = po4.se.lr, "r.co.lr" = r.co.lr, "r.se.lr" = r.se.lr,
                       "y.int.lr" = y.int.lr, "y.int.se.lr" = y.int.se.lr, "cell.r.m" = cell.r.m, "cell.r.mu"=cell.r.mu, 
                       "cell.r.b"= cell.r.b, "cell.r.bu"= cell.r.bu, "lith.r.m" = lith.r.m, "lith.r.mu"= lith.r.mu, 
                       "lith.r.b"= lith.r.b, "lith.r.bu"= lith.r.bu,"nrow.cal.df" = nrow.cal.df, "clean.rad" = clean.rad, 
                       "clean.po4" = clean.po4, "clean.mui" = clean.mui, "clean.diamcc" = clean.diamcc)
############################################################################################  


# Coefficient model for inversion 
############################################################################################
model.string="model {
# DATA MODEL
rad.p = 1/(0.5^2)
for (i in 1:length(nrow.cal.df)){
  clean.rad[i] ~ dnorm(rm[i], rad.p)
}
po4.p = 1/(0.2^2)
for (i in 1:length(nrow.cal.df)){
  clean.po4[i] ~ dnorm(po4[i], po4.p)
}
mui.p = 1/(0.3*10^-6)^2
for (i in 1:length(nrow.cal.df)){
  clean.mui[i] ~ dnorm(mui[i], mui.p)
}
diamcc.p = 1/1^2
for (i in 1:length(nrow.cal.df)){
  clean.diamcc[i] ~ dnorm(diam.cocco[i], diamcc.p)
}

# PROXY MODEL
# eps.f = max fractionation by RuBisCO at infinite CO2 (generally 25 to 28‰)
eps.f <- 28
# epsilon of diffusive transport of CO2(aq) in water 
eps.d <- 0.7 
  
for (i in 1:nrow.cal.df){
# Calculate instantaneous growth rate (mu,i) from [PO4] and rmean with 
mui[i] <- y.int + co.po4*po4[i] + co.r*rm[i]

# Calculate length of the coccolith and diameter of coccosphere from mean radius (rm)
len.lith[i] <- rm[i]*lith.m + lith.b  # in um
diam.cocco[i] <- rm[i]*diam.m + diam.b # in um
}

# PRIORS
for (i in 1:nrow.cal.df){
# Concentration of phosphate (PO4; umol/kg)
po4[i] ~ dnorm(1.2, 0.5)T(0,2)
# Mean cell radius (m)
rm[i] ~ dnorm(1.5, 0.25)T(0,3) # in um
}

# Slope and intercept for lith length vs. cell radius linear regression 
lith.m ~ dnorm(lith.r.m, 1/(lith.r.mu)^2)
lith.b ~ dnorm(lith.r.b, 1/(lith.r.bu)^2)

# Slope and intercept for coccosphere diameter vs. cell radius linear regression 
diam.m ~ dnorm(cell.r.m, 1/(cell.r.mu)^2)
diam.b ~ dnorm(cell.r.b, 1/(cell.r.bu)^2)

# Coefficient for mu(i) - PO4 multi linear regression 
co.po4 ~ dnorm(po4.co.lr, 1/(po4.se.lr)^2)

# Coefficient for radius - PO4 multi linear regression 
co.r ~ dnorm(r.co.lr, 1/(r.se.lr)^2)

# Coefficient for y intercept for mu(i) multi linear regression 
y.int ~ dnorm(y.int.lr, 1/(y.int.se.lr)^2)
}"

writeLines(model.string, con = "model_out/coeff_model.txt")
coeff.out = jags.parallel(data=data.pass.coeff, model.file = "model_out/coeff_model.txt", parameters.to.save = parms1,
                          inits = NULL, n.chains = 3, n.iter = 8000,
                          n.burnin = 2000, n.thin = 5)

coeff.mat <- cbind(coeff.out$BUGSoutput$sims.list$lith.m, coeff.out$BUGSoutput$sims.list$lith.b, 
                   coeff.out$BUGSoutput$sims.list$diam.m, coeff.out$BUGSoutput$sims.list$diam.b,
                   coeff.out$BUGSoutput$sims.list$co.po4, coeff.out$BUGSoutput$sims.list$co.r,
                   coeff.out$BUGSoutput$sims.list$y.int)
############################################################################################


# Generate look up tables for equilibrium constants 
############################################################################################
# Set upper and lower STP bounds for equil constant array 
tempC.lb = 0
tempC.ub = 65
sal.lb = 15
sal.ub = 60

# Step increments for sal (ppt) temp (degrees C) and press (bar)
t.inc = 0.25
s.inc = 0.25

# Ranges of variables over which to evaluate
tempC.vr = seq(tempC.lb, tempC.ub, by=t.inc)
sal.vr = seq(sal.lb, sal.ub, by=s.inc)

# Initiate arrays 
temp.vr = c(1:length(tempC.vr))
base2Darray = c(1:(length(tempC.vr)*length(sal.vr)))
dim(base2Darray) = c((length(tempC.vr)), (length(sal.vr)))
Ksw_sta = base2Darray
K0a = base2Darray

# Constant (cm^3 bar mol^-1 K^-1)
R <- 83.131 

# Calculate 2D array for K0 and Ksw (temp and sal dependent) 
for (i in 1:length(tempC.vr)){
  for (j in 1:length(sal.vr)){
    temp.vr[i] <- tempC.vr[i]+273.15
    Ksw_sta[i,j] <- exp(148.96502-13847.26/temp.vr[i]-23.6521*(log(temp.vr[i]))+(118.67/temp.vr[i]-5.977+1.0495*(log(temp.vr[i])))*(sal.vr[j]^0.5)-0.01615*sal.vr[j])
    K0a[i,j] <- exp(9345.17/temp.vr[i]-60.2409+23.3585*(log(temp.vr[i]/100))+sal.vr[j]*(0.023517-0.00023656*temp.vr[i]+0.0047036*((temp.vr[i]/100)^2)))
  }
}
############################################################################################


# Load proxy data to evaluate against 
############################################################################################

# Read in proxy time series data
prox.in <- read.csv('data/timeseries925.csv')
prox.in <- prox.in[,c(1:9)]
names(prox.in) <- c("age","d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", "d13Cpf.data.sd", 
                    "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd")

# Setup age range and bins 
ages.prox <- unique(round(prox.in$age, digits=1))
ages.prox.max <- max(ages.prox)
dt <- abs(diff(ages.prox, lag=1))
ages.prox.ai <- seq(1,length(ages.prox), by=1)

# Age index proxy data
prox.in <- transform(prox.in,ai=as.numeric(factor(round(age*-1, digits=1))))

# Groom input data 
clean.d13Cmarker <- prox.in[complete.cases(prox.in$d13Cmarker.data), ]
clean.d13Cpf <- prox.in[complete.cases(prox.in$d13Cpf.data), ]
clean.len.lith <- prox.in[complete.cases(prox.in$len.lith.data), ]
clean.Uk <- prox.in[complete.cases(prox.in$Uk.data), ]

# Vector of age indexes that contain d13Cmarker proxy data (with duplicates)
ai.d13Cmarker <- c(clean.d13Cmarker$ai)    

# Vector of age indexes that contain d13Cpf proxy data (with duplicates)
ai.d13Cpf <- c(clean.d13Cpf$ai)    

# Vector of age indexes that contain len.lith proxy data (with duplicates)
ai.len.lith <- c(clean.len.lith$ai)     

# Vector of age indexes that contain Uk'37 proxy data (with duplicates)
ai.Uk <- c(clean.Uk$ai)

# Vector of age indexes for all data
ai.all <- c(ai.d13Cmarker, ai.d13Cpf, ai.len.lith, ai.Uk)

# Index vector which contains each environmental time step that has one or more proxy data

ai.prox <-  unique(ai.all)     
ai.prox <- sort(ai.prox, decreasing = FALSE) 
ages.prox <- sort(ages.prox, decreasing = TRUE) 
n.steps <- length(ai.prox)

# Prior time bin vectors for which there are proxy data (includes duplicates)
ai.d13Cmarker <- match(ai.d13Cmarker, ai.prox)
ai.d13Cpf <- match(ai.d13Cpf, ai.prox)
ai.len.lith <- match(ai.len.lith, ai.prox)
ai.Uk <- match(ai.Uk, ai.prox)
############################################################################################


# Set prior distributions to send to model
############################################################################################
# Temperature (degrees C)
tempC.m = 25
tempC.p = 1/5^2

# Salintiy (ppt)
sal.m = 35
sal.p = 1/0.5^2

# pCO2 (uatm)
pco2.u = 500
pco2.l = 100

# d13C of aqueous CO2 (per mille)
d13C.co2.m = -8
d13C.co2.p = 1/1^2

# Concentration of phosphate (PO4; umol/kg)
po4.m = 1.2
po4.p = 1/0.2^2

# Mean cell radius (m)
rm.m = 2*10^-6
rm.p = 1/(0.5*10^-6)^2
############################################################################################

# Option to load ice core data-derived mui and size equation parameters 
############################################################################################
coeff.mat.ice <- read.csv("model_out/coeff_mat_ice.csv")
############################################################################################


# Select data to pass to jags 
############################################################################################
data.pass = list("coeff.mat" = coeff.mat,#coeff.mat.ice,
                  "K0a" = K0a,
                  "Ksw_sta" = Ksw_sta,
                  "sal.lb" = sal.lb,
                  "tempC.lb" = tempC.lb,
                  "t.inc" = t.inc,
                  "s.inc" = s.inc,
                  "d13Cmarker.data" = clean.d13Cmarker$d13Cmarker.data,
                  "d13Cmarker.data.sd" = clean.d13Cmarker$d13Cmarker.data.sd,
                  "d13Cpf.data" = clean.d13Cpf$d13Cpf.data,
                  "d13Cpf.data.sd" = clean.d13Cpf$d13Cpf.data.sd,
                  "len.lith.data" = clean.len.lith$len.lith.data,
                  "len.lith.data.sd" = clean.len.lith$len.lith.data.sd,
                  "Uk.data" = clean.Uk$Uk.data,
                  "Uk.data.sd" = clean.Uk$Uk.data.sd,
                  "n.steps" = n.steps,
                  "dt" = dt,
                  "ages.prox" = ages.prox,
                  "ai.prox" = ai.prox, 
                  "ai.d13Cmarker" = ai.d13Cmarker,
                  "ai.d13Cpf" = ai.d13Cpf,
                  "ai.len.lith" = ai.len.lith, 
                  "ai.Uk" = ai.Uk, 
                  "tempC.m" = tempC.m,
                  "tempC.p" = tempC.p,
                  "sal.m" = sal.m,
                  "sal.p" = sal.p,
                  "pco2.u" = pco2.u,
                  "pco2.l" = pco2.l,
                  "d13C.co2.m" = d13C.co2.m,
                  "d13C.co2.p" = d13C.co2.p,
                  "po4.m" = po4.m,
                  "po4.p" = po4.p,
                  "rm.m" = rm.m,
                  "rm.p" = rm.p) 
############################################################################################


# Parameters to save as output 
############################################################################################
parms2 = c("tempC", "sal", "pco2", "d13C.co2", "po4", "rm", "b", "coeff.po4", "coeff.rm", "mui.y.int", 
           "lith.m", "lith.b", "cocco.m", "cocco.b")
############################################################################################


# Run the inversion using jags 
############################################################################################
inv.out = jags.parallel(data = data.pass, model.file = "phytoPSM.R", parameters.to.save = parms2,
                          inits = NULL, n.chains = 3, n.iter = 10000,
                          n.burnin = 3000, n.thin = 5)
############################################################################################
#130k iteration + 30k burn-in takes ~1.5 hours for ~20 samples

