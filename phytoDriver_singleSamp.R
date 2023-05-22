

# Phytoplankton forward PSM driver
############################################################################################  
# Runs Bayesian inversion on calibration data to determine slopes and intercepts for the:
#   1) mu(i) = f(po4, rm) relationship and 
#   2) lith length and coccosphere diameter to mean cell radius transfer functions
#
# Sets up environmental model inversion analysis using PSM (separate file)
#
# Dustin T. Harper
# 2 May 2023
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
cal.df.lr <- read.csv('data/caldata_culture_rad.csv')

# Generate multiple linear regression model mu,i (as a function of [PO4] and mean radius)
igr.model = lm(formula = cal.df.lr$mui ~ cal.df.lr$po4 + cal.df.lr$r, data = cal.df.lr)
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
cal.df <- read.csv('data/caldata_culture.csv')

# remove blanks from individual data columns; clean data 
clean.rad <- cal.df[complete.cases(cal.df$radius), which(colnames(cal.df)=='radius')]
clean.po4 <- cal.df[complete.cases(cal.df$po4), which(colnames(cal.df)=='po4')]
clean.mui <- cal.df[complete.cases(cal.df$mui), which(colnames(cal.df)=='mui')]
clean.diamcc <- cal.df[complete.cases(cal.df$diam.cc), which(colnames(cal.df)=='diam.cc')]

# Cull out of range data sets
cal.df <- subset(cal.df, is.na(cal.df$radius) | is.na(cal.df$po4) | cal.df$radius <= 10 & cal.df$po4 <= 2)
nrow.cal.df <- nrow(cal.df)

# index data rows
rad.index <- which(!is.na(cal.df$radius))
po4.index <- which(!is.na(cal.df$po4))
mui.index <- which(!is.na(cal.df$mui))
diamcc.index <- which(!is.na(cal.df$diam.cc))

# Set priors for cells without data
po4.blank.prior <- 1.2
cal.df$po4[is.na(cal.df$po4)] <- po4.blank.prior
cal.po4 <- cal.df$po4

r.blank.prior <- 2.5
cal.df$radius[is.na(cal.df$radius)] <- r.blank.prior
cal.radius <- cal.df$radius

# Parameters to save in inversion output
parms1 = c("mui", "diam.cocco", "rm", "po4",  "lith.m", "lith.b", "diam.m", "diam.b", "y.int", "co.r", "co.po4") 

# data to pass to jags
data.pass.coeff = list("po4.co.lr" = po4.co.lr, "po4.se.lr" = po4.se.lr, "r.co.lr" = r.co.lr, "r.se.lr" = r.se.lr,
                 "y.int.lr" = y.int.lr, "y.int.se.lr" = y.int.se.lr, "cell.r.m" = cell.r.m, "cell.r.mu"=cell.r.mu, 
                 "cell.r.b"=cell.r.b, "cell.r.bu"=cell.r.bu, "lith.r.m" = cell.r.m, "lith.r.mu"=cell.r.mu, 
                 "lith.r.b"=cell.r.b, "lith.r.bu"=cell.r.bu,"cal.po4" = cal.po4, "cal.radius" = cal.radius, 
                 "nrow.cal.df" = nrow.cal.df, "clean.rad" = clean.rad, "clean.po4" = clean.po4, "clean.mui" = clean.mui, 
                 "clean.diamcc" = clean.diamcc, "rad.index" = rad.index, "po4.index" = po4.index, 
                 "mui.index" = mui.index, "diamcc.index" = diamcc.index)
############################################################################################  


# Coefficient model for inversion 
############################################################################################
model.string="model {
# DATA MODEL
rad.p = 1/(1*10^-6)^2
#for (i in 1:length(rad.index)){
#  clean.rad[i] ~ dnorm(rm[rad.index[i]], rad.p)
#}
po4.p = 1/(0.1*10^-6)^2
for (i in 1:length(po4.index)){
  clean.po4[i] ~ dnorm(po4[po4.index[i]], po4.p)
}
mui.p = 1/(0.3*10^-6)^2
for (i in 1:length(mui.index)){
  clean.mui[i] ~ dnorm(mui[mui.index[i]], mui.p)
}
diamcc.p = 1/0.4^2
for (i in 1:length(diamcc.index)){
  clean.diamcc[i] ~ dnorm(diam.cocco[diamcc.index[i]], diamcc.p)
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
len.lith[i] <- (rm[i]*10^6)*lith.m + lith.b  # in um
diam.cocco[i] <- (rm[i]*10^6)*diam.m + diam.b # in um
}

# PRIORS
for (i in 1:nrow.cal.df){
# Concentration of phosphate (PO4; umol/kg)
po4[i] ~ dnorm(cal.po4[i], po4.p)T(0,10)
# Mean cell radius (m)
rm[i] ~ dnorm(cal.radius[i]*10^-6, rad.p)T(0,10)
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


# Read in proxy data to evaluate against 
############################################################################################
prox.in <- read.csv('data/timeseries_data.csv')
prox.in <- prox.in[,c(1:9)]
names(prox.in) <- c("age","d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", "d13Cpf.data.sd", 
                    "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd")

############################################################################################


# Set prior distributions to send to model
############################################################################################
# Temperature (degrees C)
tempC.m = 25
tempC.sd = 5
# Salintiy (ppt)
sal.m = 35
sal.sd = 0.5
# pCO2 (uatm)
pco2.u = 500
pco2.l = 100
# d13C of aqueous CO2 (per mille)
d13C.co2.m = -8
d13C.co2.sd = 1
# Concentration of phosphate (PO4; umol/kg)
po4.m = 1.2
po4.sd = 0.2
# Mean cell radius (m)
rm.m = 2.5*10^-6
rm.sd = 0.25*10^-6
############################################################################################


# Select data to pass to jags 
############################################################################################
data.pass2ka = list("coeff.mat" = coeff.mat,
                      "len.lith.data" = prox.in$len.lith.data[1],
                      "len.lith.data.sd" = prox.in$len.lith.data.sd[1],
                      "d13C.pf.data" = prox.in$d13Cpf.data[1],
                      "d13C.pf.data.sd" = prox.in$d13Cpf.data.sd[1],
                      "d13C.marker.data" = prox.in$d13Cmarker.data[1],
                      "d13C.marker.data.sd" = prox.in$d13Cmarker.data.sd[1],
                      "Uk.data" = prox.in$Uk.data[1],
                      "Uk.data.sd" = prox.in$Uk.data.sd[1],
                      "tempC.m" = tempC.m,
                      "tempC.sd" = tempC.sd,
                      "sal.m" = sal.m,
                      "sal.sd" = sal.sd,
                      "pco2.u" = pco2.u,
                      "pco2.l" = pco2.l,
                      "d13C.co2.m" = d13C.co2.m,
                      "d13C.co2.sd" = d13C.co2.sd,
                      "po4.m" = po4.m,
                      "po4.sd" = po4.sd,
                      "rm.m" = rm.m,
                      "rm.sd" = rm.sd) 

data.pass40.8ka = list("coeff.mat" = coeff.mat,
                      "len.lith.data" = prox.in$len.lith.data[2],
                      "len.lith.data.sd" = prox.in$len.lith.data.sd[2],
                      "d13C.pf.data" = prox.in$d13Cpf.data[2],
                      "d13C.pf.data.sd" = prox.in$d13Cpf.data.sd[2],
                      "d13C.marker.data" = prox.in$d13Cmarker.data[2],
                      "d13C.marker.data.sd" = prox.in$d13Cmarker.data.sd[2],
                      "Uk.data" = prox.in$Uk.data[2],
                      "Uk.data.sd" = prox.in$Uk.data.sd[2],
                      "tempC.m" = tempC.m,
                      "tempC.sd" = tempC.sd,
                      "sal.m" = sal.m,
                      "sal.sd" = sal.sd,
                      "pco2.u" = pco2.u,
                      "pco2.l" = pco2.l,
                      "d13C.co2.m" = d13C.co2.m,
                      "d13C.co2.sd" = d13C.co2.sd,
                      "po4.m" = po4.m,
                      "po4.sd" = po4.sd,
                      "rm.m" = rm.m,
                      "rm.sd" = rm.sd) 

data.pass348ka = list("coeff.mat" = coeff.mat,
                  "len.lith.data" = prox.in$len.lith.data[3],
                  "len.lith.data.sd" = prox.in$len.lith.data.sd[3],
                  "d13C.pf.data" = prox.in$d13Cpf.data[3],
                  "d13C.pf.data.sd" = prox.in$d13Cpf.data.sd[3],
                  "d13C.marker.data" = prox.in$d13Cmarker.data[3],
                  "d13C.marker.data.sd" = prox.in$d13Cmarker.data.sd[3],
                  "Uk.data" = prox.in$Uk.data[3],
                  "Uk.data.sd" = prox.in$Uk.data.sd[3],
                  "tempC.m" = tempC.m,
                  "tempC.sd" = tempC.sd,
                  "sal.m" = sal.m,
                  "sal.sd" = sal.sd,
                  "pco2.u" = pco2.u,
                  "pco2.l" = pco2.l,
                  "d13C.co2.m" = d13C.co2.m,
                  "d13C.co2.sd" = d13C.co2.sd,
                  "po4.m" = po4.m,
                  "po4.sd" = po4.sd,
                  "rm.m" = rm.m,
                  "rm.sd" = rm.sd) 

data.pass403ka = list("coeff.mat" = coeff.mat,
                    "len.lith.data" = prox.in$len.lith.data[4],
                    "len.lith.data.sd" = prox.in$len.lith.data.sd[4],
                    "d13C.pf.data" = prox.in$d13Cpf.data[4],
                    "d13C.pf.data.sd" = prox.in$d13Cpf.data.sd[4],
                    "d13C.marker.data" = prox.in$d13Cmarker.data[4],
                    "d13C.marker.data.sd" = prox.in$d13Cmarker.data.sd[4],
                    "Uk.data" = prox.in$Uk.data[4],
                    "Uk.data.sd" = prox.in$Uk.data.sd[4],
                    "tempC.m" = tempC.m,
                    "tempC.sd" = tempC.sd,
                    "sal.m" = sal.m,
                    "sal.sd" = sal.sd,
                    "pco2.u" = pco2.u,
                    "pco2.l" = pco2.l,
                    "d13C.co2.m" = d13C.co2.m,
                    "d13C.co2.sd" = d13C.co2.sd,
                    "po4.m" = po4.m,
                    "po4.sd" = po4.sd,
                    "rm.m" = rm.m,
                    "rm.sd" = rm.sd) 
############################################################################################

# Parameters to save as output 
parms2 = c("tempC", "sal", "pco2", "d13C.co2", "po4", "rm", "b")

# Run the inversion using jags - 2 ka ODP Hole 688B
inv.out.2ka = jags.parallel(data = data.pass2ka, model.file = "phytoPSM_singleSamp.R", parameters.to.save = parms2,
                              inits = NULL, n.chains = 3, n.iter = 10000,
                              n.burnin = 2000, n.thin = 10)

# Run the inversion using jags - 40.8 ka ODP Hole 688B
inv.out.40.8ka = jags.parallel(data = data.pass40.8ka, model.file = "phytoPSM_singleSamp.R", parameters.to.save = parms2,
                              inits = NULL, n.chains = 3, n.iter = 10000,
                              n.burnin = 2000, n.thin = 10)

# Run the inversion using jags - 348.4 ka ODP Hole 688B
inv.out.348ka = jags.parallel(data = data.pass348ka, model.file = "phytoPSM_singleSamp.R", parameters.to.save = parms2,
                          inits = NULL, n.chains = 3, n.iter = 10000,
                          n.burnin = 2000, n.thin = 10)

# Run the inversion using jags - 403.1 ka ODP Hole 688B
inv.out.403ka = jags.parallel(data=data.pass403ka, model.file = "phytoPSM_singleSamp.R", parameters.to.save = parms2,
                        inits = NULL, n.chains = 3, n.iter = 10000,
                        n.burnin = 2000, n.thin = 10)
############################################################################################




