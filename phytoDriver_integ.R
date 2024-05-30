

# Phytoplankton forward PSM driver for paleo reconstruction, mui calibration integrated in jags model
# Incorporates both modern observational data and GIG data to calibrate mui = f(radius, po4)
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
cal.df <- subset(cal.df, is.na(cal.df$radius) | is.na(cal.df$po4) | cal.df$radius <= 10 & cal.df$po4 <= 2)
cal.df$radius <- cal.df$radius*1e-6

# Generate multiple linear regression model mu,i (as a function of [PO4] and mean radius)
igr.model = lm(formula = cal.df$mui ~ cal.df$po4 + cal.df$radius, data = cal.df)
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


# Load GIG calibration data to evaluate against 
############################################################################################
prox.in.gig <- read.csv('data/timeseriesGIG.csv')
#prox.in.gig <- prox.in.gig[,c(1:12)]
names(prox.in.gig) <- c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", 
                    "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd", "iceco2.data")
prox.in.gig <- prox.in.gig[complete.cases(prox.in.gig[,c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", 
                                             "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd", "iceco2.data")]), ]

# Site index proxy data
prox.in.gig <- prox.in.gig[order(prox.in.gig$site),]
prox.in.gig <- transform(prox.in.gig,site.index=as.numeric(factor(site)))
site.index.gig <- c(prox.in.gig$site.index)

# Set prior distributions for GIG calibration data 

# Temperature (degrees C)
tempC.m.gig = 22
tempC.p.gig = 1/5^2

# Salintiy (ppt)
sal.m.gig = 35
sal.p.gig = 1/0.5^2

# pCO2 (uatm)
pco2.m.gig = 250
pco2.p.gig = 1/40^2

# d13C of aqueous CO2 (per mille)
d13C.co2.m.gig = -8
d13C.co2.p.gig = 1/1^2

# Concentration of phosphate (PO4; umol/kg)
po4.m.gig = unique(prox.in.gig$po4.prior)
po4.p.gig = 1/0.25^2 # 1/0.1^2

# Mean cell radius (m)
rm.m.gig = 1.5*10^-6
rm.p.gig = 1/(0.5*10^-6)^2
############################################################################################


# Load time series proxy data for reconstruction
############################################################################################
prox.in <- read.csv('data/timeseriesGIG.csv')
#prox.in <- prox.in[,c(1:12)]
names(prox.in) <- c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", 
                    "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd", "iceco2.data")
prox.in <- prox.in[complete.cases(prox.in[,c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", 
                                             "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd", "iceco2.data")]), ]

# Site index proxy data
prox.in <- prox.in[order(prox.in$site),]
prox.in <- transform(prox.in,site.index=as.numeric(factor(site)))
site.index <- c(prox.in$site.index)
ages.prox <- unique(prox.in$age)
n.steps <- length(ages.prox)

# Set prior distributions for time series data 
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
d13C.co2.p = 1/1^2

# Concentration of phosphate (PO4; umol/kg)
po4.m.cd = 1.5
po4.m = unique(prox.in$po4.prior)
po4.p = 1/0.25^2 # 1/0.1^2

# Mean cell radius (m)
rm.m = 1.5*10^-6
rm.p = 1/(0.5*10^-6)^2
############################################################################################


# Select data to pass to jags 
############################################################################################
data.pass = list("po4.co.lr" = po4.co.lr, 
                 "po4.se.lr" = po4.se.lr, 
                 "r.co.lr" = r.co.lr, 
                 "r.se.lr" = r.se.lr,
                 "y.int.lr" = y.int.lr, 
                 "y.int.se.lr" = y.int.se.lr, 
                 "cocco.m" = cocco.m,
                 "cocco.b" = cocco.b,
                 "lith.m" = lith.m,
                 "lith.b" = lith.b,
                 "radius.cd" = cal.df$radius,
                 "po4.cd.data" = cal.df$po4,
                 "mui.cd.data" = cal.df$mui,
                 "K0a" = K0a,
                 "Ksw_sta" = Ksw_sta,
                 "sal.lb" = sal.lb,
                 "tempC.lb" = tempC.lb,
                 "t.inc" = t.inc,
                 "s.inc" = s.inc,
                 "d13Cmarker.data.gig" = prox.in.gig$d13Cmarker.data,
                 "d13Cmarker.data.sd.gig" = prox.in.gig$d13Cmarker.data.sd,
                 "d13Cpf.data.gig" = prox.in.gig$d13Cpf.data,
                 "d13Cpf.data.sd.gig" = prox.in.gig$d13Cpf.data.sd,
                 "len.lith.data.gig" = prox.in.gig$len.lith.data,
                 "len.lith.data.sd.gig" = prox.in.gig$len.lith.data.sd,
                 "Uk.data.gig" = prox.in.gig$Uk.data,
                 "Uk.data.sd.gig" = prox.in.gig$Uk.data.sd,
                 "iceco2.data.gig" = prox.in.gig$iceco2.data,
                 "site.index.gig" = site.index.gig,
                 "tempC.m.gig" = tempC.m.gig,
                 "tempC.p.gig" = tempC.p.gig,
                 "sal.m.gig" = sal.m.gig,
                 "sal.p.gig" = sal.p.gig,
                 "pco2.m.gig" = pco2.m.gig,
                 "pco2.p.gig" = pco2.p.gig,
                 "d13C.co2.m.gig" = d13C.co2.m.gig,
                 "d13C.co2.p.gig" = d13C.co2.p.gig,
                 "po4.m.gig" = po4.m.gig,
                 "po4.p.gig" = po4.p.gig,
                 "rm.m.gig" = rm.m.gig,
                 "rm.p.gig" = rm.p.gig,
                 "d13Cmarker.data" = prox.in$d13Cmarker.data,
                 "d13Cmarker.data.sd" = prox.in$d13Cmarker.data.sd,
                 "d13Cpf.data" = prox.in$d13Cpf.data,
                 "d13Cpf.data.sd" = prox.in$d13Cpf.data.sd,
                 "len.lith.data" = prox.in$len.lith.data,
                 "len.lith.data.sd" = prox.in$len.lith.data.sd,
                 "Uk.data" = prox.in$Uk.data,
                 "Uk.data.sd" = prox.in$Uk.data.sd,
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
                 "po4.m.cd" = po4.m.cd,
                 "po4.p" = po4.p,
                 "rm.m" = rm.m,
                 "rm.p" = rm.p) 
############################################################################################


# Parameters to save as output 
############################################################################################
parms = c("tempC", "sal", "pco2", "pco2.gig", "d13C.co2", "po4", "rm", "b", "coeff.po4", "coeff.rm", "mui.y.int")
############################################################################################


# Run the inversion using jags 
############################################################################################
inv.out = jags.parallel(data = data.pass, model.file = "phytoPSM_integ.R", parameters.to.save = parms,
                          inits = NULL, n.chains = 3, n.iter = 5e5,
                          n.burnin = 2e5, n.thin = 100)

############################################################################################


View(inv.out$BUGSoutput$summary)

