
model {

# Likelihood function: this block evaluates the data against the modeled values
############################################################################################ 
len.lith.data ~ dnorm(len.lith, len.lith.p)
len.lith.p =  1/len.lith.data.sd^2
d13C.pf.data ~ dnorm(d13C.pf, d13C.pf.p)
d13C.pf.p = 1/d13C.pf.data.sd^2
d13C.marker.data ~ dnorm(d13C.marker, d13C.pf.p)
d13C.marker.p = 1/d13C.marker.data.sd^2
Uk.data ~ dnorm(Uk, Uk.p)
Uk.p = 1/Uk.data.sd^2
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
P.c ~ dnorm(5.09*10^-5, 1/(0.16*10^-5)^2)

# Calculate Uk'37 from temperature (Conte et al., 2006; sediment - AnnO linear model)
# STAND DEV NEEDS TO BE UPDATED! Determine slope and intercept via inversion 
Uk.sl <- 29.876
Uk.int <- 1.334
Uk <- (tempC + Uk.int)/Uk.sl

# Calculate aqueous CO2 from atmospheric CO2 with Henry's law (carb chem) 
temp <- tempC + 273.15
K0 <- exp(9345.17/temp-60.2409+23.3585*(log(temp/100))+sal*(0.023517-0.00023656*temp+0.0047036*((temp/100)^2)))
fco2 <- pco2*0.9968
co2 <- fco2*K0*10^-3 # mol/m^3 (uM)

# Calculate instantaneous growth rate (mu,i) from [PO4] and rmean by 
# sampling Bayesian model coefficients and intercepts
coeff.ind ~ dcat(1:length(coeff.mat[,1]))
mui <- coeff.mat[coeff.ind,5]*po4 + coeff.mat[coeff.ind,6]*rm + coeff.mat[coeff.ind,7]

# Calculate length of the coccolith and diameter of coccosphere from mean radius (rm)
len.lith <- (rm*10^6)*coeff.mat[coeff.ind,1] + coeff.mat[coeff.ind,2]  # in um
diam.cocco <- (rm*10^6)*coeff.mat[coeff.ind,3] + coeff.mat[coeff.ind,4] # in um

# Calculate Qs (the co2 flux into the cell per unit surface area of the cell membrane)
cell.vol <- 4/3*3.141593*((rm)^3) # in m^3
cell.vol.um <- 4/3*3.141593*((rm*10^6)^3) # in um^3
gam.c <- 3.154*(10^-14)*cell.vol.um
Qr <- mui * gam.c
Qs <- Qr / (4*3.141593*(rm^2))
  
# Calculate DT (temperature-sensitive diffusivity of C02(aq)in seawater) from SST using eqn. 8 of Rau et al. (1996)
R.gc <- 8.3143 # gas constant in J / K*mol
DTi <- 5.019*(10^-6)*exp(-(19510/(R.gc*temp)))
DT <- DTi*(0.9508 - 7.389*(10^-4)*tempC)

# Calculate rK (reacto-diffusive length) from SSS and SST using eqns. 6 and 7 of Rau et al. (1996)
k.one <- 8500 * ((exp(-(62800/(R.gc*temp)))) / (exp(-(62800/(R.gc*298.15)))))
k.two <- (3*10^-5) * ((exp(-(62800/(R.gc*temp)))) / (exp(-(62800/(R.gc*298.15)))))
Ksw.noP <- exp(148.96502-13847.26/temp-23.6521*(log(temp))+(118.67/temp-5.977+1.0495*(log(temp)))*(sal^0.5)-0.01615*sal)
pH.const <- 8
hyd.const <- 10^(-pH.const)
hydrox <- Ksw.noP / hyd.const
k.pr <- k.one*hydrox + k.two
rk <- sqrt(DT/k.pr)

#	Calculate b in uM from Qs, rmean, rK, DT, Pc, eps.f and eps.d using eqn. 15 of Rau et al (1996) 
b <- ((eps.f - eps.d) * Qs * ((rm/((1 + rm/rk)*DT)) + 1/P.c)) * 10^3

# Calculate eps.p from [CO2](aq), eps.f and b
eps.p <- eps.f - (b/(co2*10^3))

# Calculate foram d13C from env parm d13C.co2 (rearranged eqns 19 and 20) 
eps.aog <- (-373/temp) + 0.19  # fractionation b/w CO2 (aq) and CO2 (g)
d13C.co2g <- ((d13C.co2 + 1000) / (eps.aog/1000 +1)) - 1000 # corrected; this appears to be incorrect when rearranged from supplement eqn 19b
d13C.pf <- ((11.98 - 0.12*tempC)/1000 + 1)*(d13C.co2g + 1000) - 1000 # corrected; this appears to be incorrect when rearranged from supplement eqn 19b

# Calculate d13C.biomass from d13C.co2 and epsilon.p (eqn 2)
d13C.biomass <- ((d13C.co2 + 1000) / ((eps.p / 1000) + 1)) - 1000

# Calculate d13C.marker from d13C.biomass (eqn 22)
d13C.marker <- ((d13C.biomass + 1000) / ((eps.bob/1000)+1)) - 1000
############################################################################################  


# Environmental prior distributions 
############################################################################################   
# Temperature (degrees C)
tempC ~ dnorm(tempC.m, tempC.p) 
tempC.p = 1/tempC.sd^2
# Salintiy (ppt)
sal ~ dnorm(sal.m, sal.p) 
sal.p = 1/sal.sd^2
# pCO2 (uatm)
pco2 ~ dunif(pco2.l, pco2.u) 
#pco2.p = 1/pco2.sd^2
# d13C of aqueous CO2 (per mille)
d13C.co2 ~ dnorm(d13C.co2.m, d13C.co2.p) 
d13C.co2.p = 1/d13C.co2.sd^2
# Concentration of phosphate (PO4; umol/kg)
po4 ~ dnorm(po4.m, po4.p) 
po4.p = 1/po4.sd^2
# Mean cell radius (m)
rm ~ dnorm(rm.m, rm.p) 
rm.p = 1/rm.sd^2
############################################################################################   
}
