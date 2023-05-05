
library(rjags)
library(R2jags)
library(tidyverse)
library(viridis)

# Multi linear regression model for priors
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

# Contour plot of  mu(i), [PO4] and r(mean) multi linear regression model
############################################################################################

fig.font <- "Helvetica"
fontsize.axislabels <- 11
fontsize.scalelabels <- 9

po4.seq <- seq(0,10,0.1)
rm.seq <- seq(0,10,0.1)
mlr.fun <- function(x, y) {final_value = y.int.lr + x*po4.co.lr + y*r.co.lr}

mlr.mod <- expand.grid(X1 = po4.seq, X2 = rm.seq) %>%
  mutate(mui.p = mlr.fun(X1, X2)) %>%
  ggplot(aes(X1, X2, z = mui.p)) +
  geom_contour() +
  geom_contour_filled() +
  labs(fill=expression(paste(mu[i]," (",s^-1,")")), x= expression(paste("[", PO[4], "] ","(", mu,"mol ", kg^-1,")")), y = expression(paste(radius[bar(x)]," (",mu,m,")"))) +
  scale_y_continuous(expand = expansion(mult = 0, add = 0), limits = c(min(rm.seq), max(rm.seq))) +
  scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(min(po4.seq), max(po4.seq))) +
  scale_fill_viridis(discrete = TRUE, option = "D", labels = c("<1.00e-5", "1.00e-5 - 1.05e-5", "1.05e-5 - 1.10e-5", "1.10e-5 - 1.15e-5","1.15e-5 - 1.20e-5", "1.20e-5 - 1.25e-5","1.25e-5 - 1.30e-5", "1.30e-5 - 1.35e-5", "1.35e-5 - 1.40e-5", "1.40e-5 - 1.45e-5", "1.45e-5 - 1.50e-5","1.50e-5 - 1.55e-5", "1.55e-5 - 1.60e-5","1.60e-5 - 1.65e-5", ">1.65e-5")) + 
  theme_linedraw() +
  theme(legend.title = element_text(family = fig.font, size = 14,color = "#000000"),
        legend.text = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        legend.box.background = element_rect(size = 0.5),
        axis.text.x = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_text(family = fig.font, size = fontsize.axislabels,color = "#000000"), 
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000")) 
print(mlr.mod)

############################################################################################

# Determine size-based transfer functions using linear regression and Henderiks and Pagani 07 data

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

# Load and groom data for coefficient inversion
############################################################################################  
cal.df <- read.csv('data/caldata_culture.csv')

# remove blanks from individual data columns; clean data 
clean.rad <- cal.df[complete.cases(cal.df$radius), which(colnames(cal.df)=='radius')]
clean.po4 <- cal.df[complete.cases(cal.df$po4), which(colnames(cal.df)=='po4')]
clean.mui <- cal.df[complete.cases(cal.df$mui), which(colnames(cal.df)=='mui')]
clean.diamcc <- cal.df[complete.cases(cal.df$diam.cc), which(colnames(cal.df)=='diam.cc')]

# Cull out of range data sets (include nutrient replete calibration data only)
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
parms = c("mui", "diam.cocco", "rm", "po4",  "lith.m", "lith.b", "diam.m", "diam.b", "y.int", "co.r", "co.po4") 

# data to pass to jags
data.pass = list("po4.co.lr" = po4.co.lr, "po4.se.lr" = po4.se.lr, "r.co.lr" = r.co.lr, "r.se.lr" = r.se.lr,
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
# Data model
############################################################################################    
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
############################################################################################

# Proxy model
############################################################################################    
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
############################################################################################

# Priors
############################################################################################  
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
############################################################################################
}"
############################################################################################

writeLines(model.string, con = "model.out/coeff.model.txt")
coeff.out = jags.parallel(data=data.pass, model.file = "model.out/coeff.model.txt", parameters.to.save = parms,
                     inits = NULL, n.chains = 3, n.iter = 8000,
                     n.burnin = 2000, n.thin = 5)

coeff.mat <- cbind(coeff.out$BUGSoutput$sims.list$lith.m, coeff.out$BUGSoutput$sims.list$lith.b, 
                   coeff.out$BUGSoutput$sims.list$diam.m, coeff.out$BUGSoutput$sims.list$diam.b,
                   coeff.out$BUGSoutput$sims.list$co.po4, coeff.out$BUGSoutput$sims.list$co.r,
                   coeff.out$BUGSoutput$sims.list$y.int)

# Contour plot of  mu(i), [PO4] and r(mean) Bayesian linear regression model
############################################################################################
fig.font <- "Helvetica"
fontsize.axislabels <- 11
fontsize.scalelabels <- 9

po4.seq <- seq(0,10,0.1)
rm.seq <- seq(0,10,0.1)
mlr.fun <- function(x, y) {final_value = mean(coeff.out$BUGSoutput$sims.list$y.int) + 
                  x*mean(coeff.out$BUGSoutput$sims.list$co.po4) + 
                  y*mean(coeff.out$BUGSoutput$sims.list$co.r)}

mlr.mod <- expand.grid(X1 = po4.seq, X2 = rm.seq) %>%
  mutate(mui.p = mlr.fun(X1, X2)) %>%
  ggplot(aes(X1, X2, z = mui.p)) +
  geom_contour() +
  geom_contour_filled() +
  labs(fill=expression(paste(mu[i]," (",s^-1,")")), x= expression(paste("[", PO[4], "] ","(", mu,"mol ", kg^-1,")")), y = expression(paste(radius[bar(x)]," (",mu,m,")"))) +
  scale_y_continuous(expand = expansion(mult = 0, add = 0), limits = c(min(rm.seq), max(rm.seq))) +
  scale_x_continuous(expand = expansion(mult = 0, add = 0), limits = c(min(po4.seq), max(po4.seq))) +
  scale_fill_viridis(discrete = TRUE, option = "D", labels = c("<1.00e-5", "1.00e-5 - 1.05e-5", "1.05e-5 - 1.10e-5", "1.10e-5 - 1.15e-5","1.15e-5 - 1.20e-5", "1.20e-5 - 1.25e-5","1.25e-5 - 1.30e-5", "1.30e-5 - 1.35e-5", "1.35e-5 - 1.40e-5", "1.40e-5 - 1.45e-5", "1.45e-5 - 1.50e-5","1.50e-5 - 1.55e-5", "1.55e-5 - 1.60e-5","1.60e-5 - 1.65e-5", ">1.65e-5")) + 
  theme_linedraw() +
  theme(legend.title = element_text(family = fig.font, size = 14,color = "#000000"),
        legend.text = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        legend.box.background = element_rect(size = 0.5),
        axis.text.x = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"), 
        axis.text.y = element_text(family = fig.font, size = fontsize.scalelabels,color = "#000000"),
        axis.title.x = element_text(family = fig.font, size = fontsize.axislabels,color = "#000000"), 
        axis.title.y = element_text(family = fig.font, size = fontsize.axislabels, color = "#000000")) 
print(mlr.mod)

############################################################################################

# Monte Carlo forward proxy model to generate distributions of predicted measurements 
############################################################################################
# number of iterations in monte carlo method
iter <- 1000

# Initialize vectors 
############################################################################################
P.c <- vector("numeric", iter)
Uk.sl<- vector("numeric", iter)
Uk.int<- vector("numeric", iter)
Uk<- vector("numeric", iter)
bayco.mui<- matrix(nrow=iter, ncol=7)
bayco.rad<- matrix(nrow=iter, ncol=7)
mui<- vector("numeric", iter)
Qr<- vector("numeric", iter)
Qs<- vector("numeric", iter)
len.lith<- vector("numeric", iter)
diam.cocco<- vector("numeric", iter)
b<- vector("numeric", iter)
eps.p<- vector("numeric", iter)
d13C.pf<- vector("numeric", iter)
d13C.biomass<- vector("numeric", iter)
d13C.marker<- vector("numeric", iter)
sal<- vector("numeric", iter)
tempC<- vector("numeric", iter)
pco2<- vector("numeric", iter)
d13C.co2<- vector("numeric", iter)
po4<- vector("numeric", iter)
rm<- vector("numeric", iter)
temp<- vector("numeric", iter)
K0<- vector("numeric", iter)
fco2<- vector("numeric", iter)
co2<- vector("numeric", iter)
cell.vol<- vector("numeric", iter)
cell.vol.um<- vector("numeric", iter)
gam.c<- vector("numeric", iter)
DTi<- vector("numeric", iter)
DT<- vector("numeric", iter)
k.one<- vector("numeric", iter)
k.two<- vector("numeric", iter)
Ksw.noP<- vector("numeric", iter)
hydrox<- vector("numeric", iter)
k.pr<- vector("numeric", iter)
rk<- vector("numeric", iter)
eps.aog<- vector("numeric", iter)
d13C.co2g<- vector("numeric", iter)
############################################################################################

for (i in 1:iter){
  
# Free Parameters
############################################################################################    
# Free parameters to vary; in the stand alone forward PSM (i.e., no inversion), 
# these are input to produce predicted measured values
# Salinity (ppt)  
sal[i] <- 35 
# sal <- rnorm(1, mean=35, sd=1)
  
# Temp in C
tempC[i] <- 25
# tempC[i] <- rnorm(1, mean=25, sd=1)
  
# pCO2 (uatm)
pco2[i] <- 280 
# pco2[i] <- rnorm(1, mean=280, sd=20)
  
# d13C of co2 (aqueous; per mille)
d13C.co2[i] <- -8 
# d13C.co2[i] <- rnorm(1, mean=-8, sd=1)
  
# Concentration of phosphate (PO4; umol/kg)
po4[i] <- 1.2 
# po4[i] <- rnorm(1, mean=1.2, sd=0.2)
  
# Mean cell radius (m)
rm[i] <- 2.5*10^-6 
# rm[i] <- rnorm(1, mean=2.5*10^-6, sd=0.1*10^-6)
############################################################################################

# Proxy System Model
############################################################################################    
# Constants
# eps.f = max fractionation by RuBisCO at infinite CO2 (generally 25 to 28‰)
eps.f <- 28
# epsilon of diffusive transport of CO2(aq) in water 
eps.d <- 0.7 
# epsilon of biomass/biomarker (C isotope fractionation b/w algal biomass and biomarkers, or bob [b over b])
eps.bob <- 4.80
# cell wall permeability to CO2(aq) in m/s taken from Zhang et al. (2020) 
P.c[i] <- rnorm(1, mean=5.09*10^-5, sd=0.16*10^-5)
############################################################################################    

# Calculate Uk'37 from temperature (Conte et al., 2006; sediment - AnnO linear model)
# STAND DEV NEEDS TO BE UPDATED! Determine slope and intercept via inversion 
Uk.sl[i] <- rnorm(1, mean=0.0322, sd=0.001)
Uk.int[i] <- rnorm(1, mean=0.0709, sd=0.002)
Uk[i] <- tempC[i]*Uk.sl[i] + Uk.int[i]

# Calculate aqueous CO2 from atmospheric CO2 with Henry's law (carb chem) 
temp[i] <- tempC[i] + 273.15
K0[i] <- exp(9345.17/temp[i]-60.2409+23.3585*(log(temp[i]/100))+sal[i]*(0.023517-0.00023656*temp[i]+0.0047036*((temp[i]/100)^2)))
fco2[i] <- pco2[i] * 0.9968
co2[i] <- fco2[i] * K0[i] * 10^-3 # mol/m^3 (uM)

# Calculate instantaneous growth rate (mu,i) from [PO4] and rmean by 
# sampling Bayesian model coefficients and intercepts
bayco.mui[i,] <- coeff.mat[sample(nrow(coeff.mat),size=1,replace=TRUE),] 
mui[i] <- bayco.mui[i,5]*po4[i] + bayco.mui[i,6]*rm[i] + bayco.mui[i,7]

# Calculate Qs (the co2 flux into the cell per unit surface area of the cell membrane)
cell.vol[i] <- 4/3*pi*((rm[i])^3) # in m^3
cell.vol.um[i] <- 4/3*pi*((rm[i]*10^6)^3) # in um^3
gam.c[i] <- 3.154*(10^-14)*cell.vol.um[i]
Qr[i] <- mui[i] * gam.c[i]
Qs[i] <- Qr[i] / (4*pi*(rm[i]^2))

# Calculate length of the coccolith and diameter of coccosphere from mean radius (rm)
bayco.rad[i,] <- coeff.mat[sample(nrow(coeff.mat),size=1,replace=TRUE),] 
len.lith[i] <- (rm[i]*10^6)*bayco.rad[i,1] + bayco.rad[i,2]  # in um
diam.cocco[i] <- (rm[i]*10^6)*bayco.rad[i,3] + bayco.rad[i,4] # in um

# Calculate DT (temperature-sensitive diffusivity of C02(aq)in seawater) from SST using eqn. 8 of Rau et al. (1996)
R.gc <- 8.3143 # gas constant in J / K*mol
DTi[i] <- 5.019*(10^-6)*exp(-(19510/(R.gc*temp[i])))
DT[i] <- DTi[i]*(0.9508 - 7.389*(10^-4)*tempC[i])

# Calculate rK (reacto-diffusive length) from SSS and SST using eqns. 6 and 7 of Rau et al. (1996)
k.one[i] <- 8500 * ((exp(-(62800/(R.gc*temp[i])))) / (exp(-(62800/(R.gc*298.15)))))
k.two[i] <- (3*10^-5) * ((exp(-(62800/(R.gc*temp[i])))) / (exp(-(62800/(R.gc*298.15)))))
Ksw.noP[i] <- exp(148.96502-13847.26/temp[i]-23.6521*(log(temp[i]))+(118.67/temp[i]-5.977+1.0495*(log(temp[i])))*(sal[i]^0.5)-0.01615*sal[i])
pH.const <- 8
hyd.const <- 10^(-pH.const)
hydrox[i] <- Ksw.noP[i] / hyd.const
k.pr[i] <- k.one[i]*hydrox[i] + k.two[i]
rk[i] <- sqrt(DT[i]/k.pr[i])

#	Calculate b in uM from Qs, rmean, rK, DT, Pc, eps.f and eps.d using eqn. 15 of Rau et al (1996) 
b[i] <- ((eps.f - eps.d) * Qs[i] * ((rm[i]/((1 + rm[i]/rk[i])*DT[i])) + 1/P.c[i])) * 10^3

# Calculate eps.p from [CO2](aq), eps.f and b
eps.p[i] <- eps.f - (b[i]/(co2[i]*10^3))

# Calculate foram d13C from env parm d13C.co2 (rearranged eqns 19 and 20) 
eps.aog[i] <- (-373/temp[i]) + 0.19  # fractionation b/w CO2 (aq) and CO2 (g)
d13C.co2g[i] <- ((d13C.co2[i] + 1000) / (eps.aog[i]/1000 +1)) - 1000 # corrected; this appears to be incorrect when rearranged from supplement eqn 19b
d13C.pf[i] <- ((11.98 - 0.12*tempC[i])/1000 + 1)*(d13C.co2g[i] + 1000) - 1000 # corrected; this appears to be incorrect when rearranged from supplement eqn 19b
 
# Calculate d13C.biomass from d13C.co2 and epsilon.p (eqn 2)
d13C.biomass[i] <- ((d13C.co2[i] + 1000) / ((eps.p[i] / 1000) + 1)) - 1000
  
# Calculate d13C.marker from d13C.biomass (eqn 22)
d13C.marker[i] <- ((d13C.biomass[i] + 1000) / ((eps.bob/1000)+1)) - 1000
}
# Predicted parameters 
M.out <- cbind(d13C.pf, d13C.marker, len.lith, diam.cocco, Uk, b, eps.p)

############################################################################################

# Histograms and density plots representing proxy model derived uncertainty in model parameters
############################################################################################    

# Density plots and histograms of params resulting from uncertainty in intermediate PSM parms and relationships (e.g., transfer functions) 
# i.e., no uncertainty in free parameters - they are fixed. 

# b
hist(b, probability = TRUE, xlim=c(50,250), ylim=c(0,0.12))
lines(density(b), col="blue")
abline(v=mean(b), col="red")
abline(v=(2*sd(b)+mean(b)), col="pink")
abline(v=(mean(b)-2*sd(b)), col="pink")

# eps.p
hist(eps.p, probability = TRUE, xlim=c(-5,25), ylim=c(0,1))
lines(density(eps.p), col="blue")
abline(v=mean(eps.p), col="red")
abline(v=(2*sd(eps.p)+mean(eps.p)), col="pink")
abline(v=(mean(eps.p)-2*sd(eps.p)), col="pink")

# Uk
hist(Uk, probability = TRUE, xlim=c(0.6,1.2), ylim=c(0,20))
lines(density(Uk), col="blue")
abline(v=mean(Uk), col="red")
abline(v=(2*sd(Uk)+mean(Uk)), col="pink")
abline(v=(mean(Uk)-2*sd(Uk)), col="pink")

# diam.cocco
hist(diam.cocco, probability = TRUE, xlim=c(0,14), ylim=c(0,1.5))
lines(density(diam.cocco), col="blue")
abline(v=mean(diam.cocco), col="red")
abline(v=(2*sd(diam.cocco)+mean(diam.cocco)), col="pink")
abline(v=(mean(diam.cocco)-2*sd(diam.cocco)), col="pink")

# d13C.marker
hist(d13C.marker, probability = TRUE, xlim=c(-35,-10), ylim=c(0,1))
lines(density(d13C.marker), col="blue")
abline(v=mean(d13C.marker), col="red")
abline(v=(2*sd(d13C.marker)+mean(d13C.marker)), col="pink")
abline(v=(mean(d13C.marker)-2*sd(d13C.marker)), col="pink")

# d13C.pf
hist(d13C.pf, probability = TRUE, xlim=c(-4,8), ylim=c(0,2), breaks=seq(-4,8,0.2))
lines(density(d13C.pf), col="blue")
abline(v=mean(d13C.pf), col="red")
abline(v=(2*sd(d13C.pf)+mean(d13C.pf)), col="pink")
abline(v=(mean(d13C.pf)-2*sd(d13C.pf)), col="pink")

############################################################################################
