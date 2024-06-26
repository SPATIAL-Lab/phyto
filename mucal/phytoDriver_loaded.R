

# Phytoplankton forward PSM driver for paleo reconstruction
# Matrix of mui coefficients and y intercepts ('coeff.mat') is loaded from 
# previous mui calibration inversion (line ~117 below) 
# Dustin T. Harper
############################################################################################  

# Load libraries 
############################################################################################  
library(rjags)
library(R2jags)
############################################################################################ 


# Determine size-based transfer functions using linear regression and Henderiks and Pagani 07 data
############################################################################################
hp07 <- read.csv('data/HP07_sizedata.csv')
lm.lith <- lm(hp07$lith.size ~ hp07$cell.r, data = hp07)
lith.sum <- summary(lm.lith)
lm.cell <- lm(hp07$sphere.d ~ hp07$cell.r, data = hp07)
cell.sum <- summary(lm.cell)

lith.m <- lith.sum$coefficients[2,1]
lith.b <- lith.sum$coefficients[1,1]
cocco.m <- cell.sum$coefficients[2,1]
cocco.b <- cell.sum$coefficients[1,1]
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
prox.in <- read.csv('data/timeseriesGIG.csv')
prox.in <- prox.in[,c(1:11)]
names(prox.in) <- c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", 
                    "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd")
prox.in <- prox.in[complete.cases(prox.in[,c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", 
                                             "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd")]), ]

# Site index proxy data
prox.in <- prox.in[order(prox.in$site),]
prox.in <- transform(prox.in,site.index=as.numeric(factor(site)))
site.index <- c(prox.in$site.index)

# Data ages
ages.prox <- unique(round(prox.in$age, digits=1))
ages.prox <- sort(ages.prox, decreasing = TRUE) 
############################################################################################


# Set prior distributions to send to model
############################################################################################
# Temperature (degrees C)
tempC.m = 22
tempC.p = 1/5^2

# Salintiy (ppt)
sal.m = 35
sal.p = 1/0.5^2

# pCO2 (uatm)
pco2.m = 250
pco2.p = 1/40^2

# d13C of aqueous CO2 (per mille)
d13C.co2.m = -8
d13C.co2.p = 1/0.5^2

# Concentration of phosphate (PO4; umol/kg)
po4.m = unique(prox.in$po4.prior)
po4.p = 1/0.1^2

# Mean cell radius (m)
rm.m = 1.5*10^-6
rm.p = 1/(0.5*10^-6)^2
############################################################################################

# Load mui and size equation parameters 
############################################################################################
coeff.mat <- read.csv("model_out/coeff_mui_noModern.csv")
coeff.mat <- coeff.mat[,2:4]
############################################################################################


# Select data to pass to jags 
############################################################################################
data.pass = list("coeff.mat" = coeff.mat,
                 "lith.m" = lith.m,
                 "lith.b" = lith.b,
                 "cocco.m" = cocco.m,
                 "cocco.b" = cocco.b,
                 "K0a" = K0a,
                 "Ksw_sta" = Ksw_sta,
                 "sal.lb" = sal.lb,
                 "tempC.lb" = tempC.lb,
                 "t.inc" = t.inc,
                 "s.inc" = s.inc,
                 "d13Cmarker.data" = prox.in$d13Cmarker.data,
                 "d13Cmarker.data.sd" = prox.in$d13Cmarker.data.sd,
                 "d13Cpf.data" = prox.in$d13Cpf.data,
                 "d13Cpf.data.sd" = prox.in$d13Cpf.data.sd,
                 "len.lith.data" = prox.in$len.lith.data,
                 "len.lith.data.sd" = prox.in$len.lith.data.sd,
                 "Uk.data" = prox.in$Uk.data,
                 "Uk.data.sd" = prox.in$Uk.data.sd,
                 "ages.prox" = ages.prox,
                 "site.index" = site.index,
                 "tempC.m" = tempC.m,
                 "tempC.p" = tempC.p,
                 "sal.m" = sal.m,
                 "sal.p" = sal.p,
                 "pco2.m" = pco2.m,
                 "pco2.p" = pco2.p,
                 "d13C.co2.m" = d13C.co2.m,
                 "d13C.co2.p" = d13C.co2.p,
                 "po4.m" = po4.m,
                 "po4.p" = po4.p,
                 "rm.m" = rm.m,
                 "rm.p" = rm.p) 
############################################################################################


# Parameters to save as output 
############################################################################################
parms2 = c("tempC", "sal", "pco2", "d13C.co2", "po4", "rm", "b", "coeff.po4", "coeff.rm", "mui.y.int")
############################################################################################


# Run the inversion using jags 
############################################################################################
inv.out = jags.parallel(data = data.pass, model.file = "phytoPSM_loaded.R", parameters.to.save = parms2,
                          inits = NULL, n.chains = 3, n.iter = 5e5,
                          n.burnin = 2e5, n.thin = 100)
############################################################################################



View(inv.out$BUGSoutput$summary)


