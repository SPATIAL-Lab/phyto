
model {

# Likelihood function: this block evaluates the data against the modeled values
############################################################################################ 

for (i in 1:length(ai.d13Cmarker)){
  d13Cmarker.data[i] ~ dnorm(d13Cmarker[ai.d13Cmarker[i]], d13Cmarker.p[i])
  d13Cmarker.p[i] = 1/d13Cmarker.data.sd[i]^2 # Gaussian precision for d13Cmarker measurements from sd 
}

for (i in 1:length(ai.d13Cpf)){
  d13Cpf.data[i] ~ dnorm(d13Cpf[ai.d13Cpf[i]], d13Cpf.p[i])
  d13Cpf.p[i] = 1/d13Cpf.data.sd[i]^2 # Gaussian precision for d13Cpf measurements from sd
}

for (i in 1:length(ai.len.lith)){
  len.lith.data[i] ~ dnorm(len.lith[ai.len.lith[i]], len.lith.p[i])
  len.lith.p[i] = 1/(len.lith.data.sd[i])^2  # Gaussian precision for length of coccolith measurements from sd
}

for (i in 1:length(ai.Uk)){
  Uk.data[i] ~ dnorm(Uk[ai.Uk[i]], Uk.p[i])
  Uk.p[i] = 1/(Uk.data.sd[i])^2  # Gaussian precision for Uk'37 measurements from sd
}
############################################################################################  
  
# Constants - time invariant
############################################################################################   
# eps.f = max fractionation by RuBisCO at infinite CO2 (generally 25 to 28‰)
eps.f <- 28

# epsilon of diffusive transport of CO2(aq) in water 
eps.d <- 0.7 

# epsilon of biomass/biomarker (C isotope fractionation b/w algal biomass and biomarkers, or bob [b over b])
eps.bob <- 4.80

# cell wall permeability to CO2(aq) in m/s taken from Zhang et al. (2020) 
P.c ~ dnorm(5.09*10^-5, 1/(0.16*10^-5)^2)

# Uk'37 temperature sensitvity (Conte et al., 2006; sediment - AnnO linear model)
# STAND DEV NEEDS TO BE UPDATED! Determine slope and intercept via inversion 
Uk.sl <- 29.876 #~dnorm(0.0322, 1/(0.001)^2)
Uk.int <- 1.334 #~dnorm(0.0709, 1/(0.002)^2)

# Sample coefficient matrix for mui = f(po4, rm) relationship
coeff.ind ~ dcat(1:length(coeff.mat[,1]))

# gas constant in J / K*mol
R.gc <- 8.3143 

# pH value for calculating rk - held constant here; varying this has almost zero effect on the model 
pH.const <- 8
hyd.const <- 10^(-pH.const)
############################################################################################  

# Proxy System Model
############################################################################################  
for (i in 1:length(ai.prox)){
# Calculate Uk'37 from temperature 
Uk[i] <- (tempC[ai.prox[i]] + Uk.int)/Uk.sl

# Calculate aqueous CO2 from atmospheric CO2 with Henry's law (carb chem) 
temp[i] <- tempC[ai.prox[i]] + 273.15
K0[i] <- exp(9345.17/temp[i]-60.2409+23.3585*(log(temp[i]/100))+sal[ai.prox[i]]*(0.023517-0.00023656*temp[i]+0.0047036*((temp[i]/100)^2)))
fco2[i] <- pco2[ai.prox[i]]*0.9968
co2[i] <- fco2[i]*K0[i]*10^-3 # mol/m^3 (uM)

# Calculate instantaneous growth rate (mu,i) from [PO4] and rmean by 
# sampling Bayesian model coefficients and intercepts
mui[i] <- coeff.mat[coeff.ind,5]*po4[ai.prox[i]] + coeff.mat[coeff.ind,6]*rm[ai.prox[i]] + coeff.mat[coeff.ind,7]

# Calculate length of the coccolith and diameter of coccosphere from mean radius (rm)
len.lith[i] <- (rm[ai.prox[i]]*10^6)*coeff.mat[coeff.ind,1] + coeff.mat[coeff.ind,2]  # in um
diam.cocco[i] <- (rm[ai.prox[i]]*10^6)*coeff.mat[coeff.ind,3] + coeff.mat[coeff.ind,4] # in um

# Calculate Qs (the co2 flux into the cell per unit surface area of the cell membrane)
cell.vol[i] <- 4/3*3.141593*((rm[ai.prox[i]])^3) # in m^3
cell.vol.um[i] <- 4/3*3.141593*((rm[ai.prox[i]]*10^6)^3) # in um^3
gam.c[i] <- 3.154*(10^-14)*cell.vol.um[i]
Qr[i] <- mui[i] * gam.c[i]
Qs[i] <- Qr[i] / (4*3.141593*(rm[ai.prox[i]]^2))
  
# Calculate DT (temperature-sensitive diffusivity of C02(aq)in seawater) from SST using eqn. 8 of Rau et al. (1996)
DTi[i] <- 5.019*(10^-6)*exp(-(19510/(R.gc*temp[i])))
DT[i] <- DTi[i]*(0.9508 - 7.389*(10^-4)*tempC[ai.prox[i]])

# Calculate rK (reacto-diffusive length) from SSS and SST using eqns. 6 and 7 of Rau et al. (1996)
k.one[i] <- 8500 * ((exp(-(62800/(R.gc*temp[i])))) / (exp(-(62800/(R.gc*298.15)))))
k.two[i] <- (3*10^-5) * ((exp(-(62800/(R.gc*temp[i])))) / (exp(-(62800/(R.gc*298.15)))))
Ksw.noP[i] <- exp(148.96502-13847.26/temp[i]-23.6521*(log(temp[i]))+(118.67/temp[i]-5.977+1.0495*(log(temp[i])))*(sal[ai.prox[i]]^0.5)-0.01615*sal[ai.prox[i]])
hydrox[i] <- Ksw.noP[i] / hyd.const
k.pr[i] <- k.one[i]*hydrox[i] + k.two[i]
rk[i] <- sqrt(DT[i]/k.pr[i])

#	Calculate b in uM from Qs, rmean, rK, DT, Pc, eps.f and eps.d using eqn. 15 of Rau et al (1996) 
b[i] <- ((eps.f - eps.d) * Qs[i] * ((rm[ai.prox[i]]/((1 + rm[ai.prox[i]]/rk[i])*DT[i])) + 1/P.c)) * 10^3

# Calculate eps.p from [CO2](aq), eps.f and b
eps.p[i] <- eps.f - (b[i]/(co2[i]*10^3))

# Calculate foram d13C from env parm d13C.co2 (rearranged eqns 19 and 20) 
eps.aog[i] <- (-373/temp[i]) + 0.19  # fractionation b/w CO2 (aq) and CO2 (g)
d13C.co2g[i] <- ((d13C.co2[ai.prox[i]] + 1000) / (eps.aog[i]/1000 +1)) - 1000 # corrected; this appears to be incorrect when rearranged from supplement eqn 19b
d13Cpf[i] <- ((11.98 - 0.12*tempC[ai.prox[i]])/1000 + 1)*(d13C.co2g[i] + 1000) - 1000 # corrected; this appears to be incorrect when rearranged from supplement eqn 19b

# Calculate d13C.biomass from d13C.co2 and epsilon.p (eqn 2)
d13C.biomass[i] <- ((d13C.co2[ai.prox[i]] + 1000) / ((eps.p[i] / 1000) + 1)) - 1000

# Calculate d13C.marker from d13C.biomass (eqn 22)
d13Cmarker[i] <- ((d13C.biomass[i] + 1000) / ((eps.bob/1000)+1)) - 1000
}
############################################################################################  


# Environmental model - time series priors
############################################################################################   
# Environmental time-dependent prior initial conditions 
# .phi = temporal autocorrelation of a parameter
# .eps =  error term 
# .tau = error precision for dt = 1
# .pc = error precision of temporal autocorrelation error term (.eps)

# Temp in C
tempC[1] ~ dnorm(tempC.m, tempC.p) 
tempC.phi ~ dbeta(5,2) 
tempC.eps[1] = 0 
tempC.tau ~ dgamma(1, 100)

# Salinity (ppt)  
sal[1] ~ dnorm(sal.m, sal.p)  
sal.phi ~ dbeta(5,2)         
sal.eps[1] = 0                 
sal.tau ~ dgamma(1e2, 5e-3) 

# pCO2 (uatm)
pco2[1] ~ dunif(pco2.l, pco2.u)   
pco2.phi ~ dbeta(5,2)
pco2.eps[1] = 0 
pco2.tau ~ dgamma(0.1, 1) 

# d13C of aqueous CO2 (per mille)
d13C.co2[1] ~ dnorm(d13C.co2.m, d13C.co2.p)    
d13C.co2.phi ~ dbeta(5,2)  
d13C.co2.eps[1] = 0 
d13C.co2.tau ~ dgamma(100,10)

# Concentration of phosphate (PO4; umol/kg)
po4[1] ~ dnorm(po4.m, po4.p)      
po4.phi ~ dbeta(5,2)  
po4.eps[1] = 0 
po4.tau ~ dgamma(1e3, 1e-5)

# Mean cell radius (m)
rm[1] ~ dnorm(rm.m, rm.p)
rm.phi ~ dbeta(5,2) 
rm.eps[1] = 0 
rm.tau ~ dgamma(1e12, 1e-3)
############################################################################################  

# Environmental time series random walk model
############################################################################################  
for (i in 2:n.steps){
  
  # Temp in C
  tempC.pc[i] <- tempC.tau*((1-tempC.phi^2)/(1-tempC.phi^(2*dt[i-1]))) 
  tempC.eps[i] ~ dnorm(tempC.eps[i-1]*(tempC.phi^dt[i-1]), tempC.pc[i])T(-12, 12)
  tempC[i] <- tempC[i-1] + tempC.eps[i]
  
  # Salinity (ppt)  
  sal.pc[i] <- sal.tau*((1-sal.phi^2)/(1-sal.phi^(2*dt[i-1])))   
  sal.eps[i] ~ dnorm(sal.eps[i-1]*(sal.phi^dt[i-1]), sal.pc[i])T(-0.3, 0.3)
  sal[i] <- sal[i-1] * (1 + sal.eps[i])
  
  # pco2 (uatm)
  pco2.pc[i] <- pco2.tau*((1-pco2.phi^2)/(1-pco2.phi^(2*dt[i-1])))
  pco2.eps[i] ~ dnorm(pco2.eps[i-1]*(pco2.phi^dt[i-1]), pco2.pc[i])
  pco2[i] <- pco2[i-1] + pco2.eps[i]
  
  # d13C of atmospheric co2 (per mille VPDB) 
  d13C.co2.pc[i] <- d13C.co2.tau*((1-d13C.co2.phi^2)/(1-d13C.co2.phi^(2*dt[i-1])))
  d13C.co2.eps[i] ~ dnorm(d13C.co2.eps[i-1]*(d13C.co2.phi^dt[i-1]), d13C.co2.pc[i])T(-2, 2)
  d13C.co2[i] <- d13C.co2[i-1] + d13C.co2.eps[i]
  
  # [PO4] (umol kg^-1)
  po4.pc[i] <- po4.tau*((1-po4.phi^2)/(1-po4.phi^(2*dt[i-1])))
  po4.eps[i] ~ dnorm(po4.eps[i-1]*(po4.phi^dt[i-1]), po4.pc[i])T(-0.1, 0.1)
  po4[i] <- po4[i-1] * (1 + po4.eps[i])
  
  # Mean cell radius (m)
  rm.pc[i] <- rm.tau*((1-rm.phi^2)/(1-rm.phi^(2*dt[i-1])))
  rm.eps[i] ~ dnorm(rm.eps[i-1]*(rm.phi^dt[i-1]), rm.pc[i])T(-1*10^-6, 1*10^-6)
  rm[i] <- rm[i-1] + rm.eps[i]
}
############################################################################################   
}
