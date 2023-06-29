

# Phytoplankton forward PSM driver without environmental time series model; includes ice core CO2 data
############################################################################################  
#
# Dustin T. Harper
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
prox.in <- prox.in[,c(1:12)]
names(prox.in) <- c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", 
                    "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd", "iceco2.data")
prox.in <- prox.in[complete.cases(prox.in[,c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data", 
                                             "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd", "iceco2.data")]), ]

# Site index proxy data
prox.in <- transform(prox.in,site.index=as.numeric(factor(site)))
site.index.d13Cmarker <- c(prox.in$site.index)
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
pco2.u = 500
pco2.l = 100

# d13C of aqueous CO2 (per mille)
d13C.co2.m = -8
d13C.co2.p = 1/1^2

# Concentration of phosphate (PO4; umol/kg)
po4.m = unique(prox.in$po4.prior)
po4.p = 1/0.1^2

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
                 "diam.cocco.cd" = cal.df$diam.cc,
                 "po4.cd" = cal.df$po4,
                 "mui.cd.data" = cal.df$mui,
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
                 "iceco2.data" = prox.in$iceco2.data,
                 "site.index.d13Cmarker" = site.index.d13Cmarker,
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
parms = c("tempC", "sal", "pco2", "d13C.co2", "po4", "rm", "b", "coeff.po4", "coeff.rm", "mui.y.int")
############################################################################################


# Run the inversion using jags 
############################################################################################
inv.out = jags.parallel(data = data.pass, model.file = "phytoPSM_mucal.R", parameters.to.save = parms,
                          inits = NULL, n.chains = 3, n.iter = 100000,
                          n.burnin = 30000, n.thin = 10)
############################################################################################


# Save model parameters governing mui relationship w/ size and po4
############################################################################################
coeff.mat <- cbind(inv.out[["BUGSoutput"]][["sims.list"]][["coeff.po4"]], 
                       inv.out[["BUGSoutput"]][["sims.list"]][["coeff.rm"]], 
                       inv.out[["BUGSoutput"]][["sims.list"]][["mui.y.int"]])
write.csv(coeff.mat, file = "model_out/coeff_mat_ice_culture.csv")


