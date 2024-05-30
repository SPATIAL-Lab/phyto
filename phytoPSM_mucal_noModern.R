
model {

# Likelihood function: this block evaluates the data against the modeled values
############################################################################################ 
  
# Time series calibration data
for (i in 1:length(d13Cmarker.data)){
  d13Cmarker.data[i] ~ dnorm(d13Cmarker[i], d13Cmarker.p[i])
  d13Cmarker.p[i] = 1/d13Cmarker.data.sd[i]^2 # Gaussian precision for d13Cmarker measurements from sd
}

for (i in 1:length(d13Cpf.data)){
  d13Cpf.data[i] ~ dnorm(d13Cpf[i], d13Cpf.p[i])
  d13Cpf.p[i] = 1/d13Cpf.data.sd[i]^2 # Gaussian precision for d13Cpf measurements from sd
}

for (i in 1:length(len.lith.data)){
  len.lith.data[i] ~ dnorm(len.lith[i], len.lith.p[i])
  len.lith.p[i] = 1/(len.lith.data.sd[i])^2  # Gaussian precision for length of coccolith measurements from sd
}

for (i in 1:length(Uk.data)){
  Uk.data[i] ~ dnorm(Uk[i], Uk.p[i])
  Uk.p[i] = 1/(Uk.data.sd[i])^2  # Gaussian precision for Uk'37 measurements from sd
}

for (i in 1:length(iceco2.data)){
  iceco2.data[i] ~ dnorm(pco2[i], 1/(6)^2)
}

  

# Constants
############################################################################################
# eps.f = max fractionation by RuBisCO at infinite CO2 (generally 25 to 28‰)
eps.f ~ dnorm(27, 1/0.3^2)

# epsilon of diffusive transport of CO2(aq) in water
eps.d ~ dnorm(0.7, 1/0.1^2)

# epsilon of biomass/biomarker (C isotope fractionation b/w algal biomass and biomarkers, or bob [b over b])
eps.bob ~ dnorm(4.80, 1/0.2^2)

# cell wall permeability to CO2(aq) in m/s taken from Zhang et al. (2020)
P.c ~ dnorm(5.09*10^-5, 1/(0.16*10^-5)^2)T(0,)

# Uk'37 temperature calibration (Conte et al., 2006; sediment - AnnO linear model) with 1.1C = reported se of estimation
Uk.sl <- 29.876
Uk.int <- 1.334
Uk.cal.se <- 1.1

# gas constant in J / K*mol
R.gc <- 8.3143

# pH value for calculating rk - held constant here; varying this has almost zero effect on the model
pH.const <- 8
hyd.const <- 10^(-pH.const)


# mui coefficient priors
############################################################################################ 

# Coefficient for mu(i) - PO4 multi linear regression 
coeff.po4 ~ dnorm(po4.co.lr, 1/(po4.se.lr)^2)
# Coefficient for radius - PO4 multi linear regression 
coeff.rm ~ dnorm(r.co.lr, 1/(r.se.lr)^2)
# Coefficient for y intercept for mu(i) multi linear regression 
mui.y.int ~ dnorm(y.int.lr, 1/(y.int.se.lr)^2)
############################################################################################ 

# Proxy System Model - GIG calibration data
############################################################################################
for (i in 1:length(d13Cmarker.data)){
  # Calculate Uk'37 from temperature
  Uk[i] <- (tempC[i] + Uk.int)/Uk.sl
  temp[i] <- tempC[i] + 273.15

  # Pull K0 and Ksw from driver look up tables
  t.index[i] <- min(max(round((tempC[i]-tempC.lb)/t.inc)+1, 0), 261)
  s.index[i] <- min(max(round((sal[i]-sal.lb)/s.inc)+1, 0), 181)
  K0[i] <- K0a[t.index[i], s.index[i]]
  Ksw.noP[i] <- Ksw_sta[t.index[i], s.index[i]]

  # Calculate aqueous CO2 from atmospheric CO2 with Henry's law (carb chem)
  fco2[i] <- pco2[i]*0.9968
  co2[i] <- fco2[i]*K0[i]*10^-3 # mol/m^3 (uM)

  # Calculate length of the coccolith and diameter of coccosphere from mean radius (rm)
  len.lith[i] <- lith.m*(rm[i]*10^6) + lith.b  # in um

  # Calculate cell carbon content
  cell.vol[i] <- 4/3*3.141593*((rm[i])^3) # in m^3
  cell.vol.um[i] <- 4/3*3.141593*((rm[i]*10^6)^3) # in um^3
  gam.c[i] <- 1.46*10^-14*cell.vol.um[i]

  # Calculate DT (temperature-sensitive diffusivity of C02(aq)in seawater) from SST using eqn. 8 of Rau et al. (1996)
  DTi[i] <- 5.019*(10^-6)*exp(-(19510/(R.gc*temp[i])))
  DT[i] <- DTi[i]*(0.9508 - 7.389*(10^-4)*tempC[i])

  # Calculate rK (reacto-diffusive length) from SSS and SST using eqns. 6 and 7 of Rau et al. (1996); rate coeffs. follow Zhang et al. (2020)
  k.one[i] <- 8718*(exp(-(62800/(R.gc*temp[i]))))
  k.two[i] <- (680.5-4.72*sal[i])*10^8*(exp(-(69400/(R.gc*temp[i]))))
  hydrox[i] <- Ksw.noP[i] / hyd.const
  k.pr[i] <- k.one[i]*hydrox[i] + k.two[i]
  rk[i] <- sqrt(DT[i]/k.pr[i])

  # Calculate foram d13C from env parm d13C.co2 (rearranged eqns 19 and 20)
  eps.aog[i] <- (-373/temp[i]) + 0.19  # fractionation b/w CO2 (aq) and CO2 (g)
  d13C.co2g[i] <- ((d13C.co2[i] + 1000) / (eps.aog[i]/1000 +1)) - 1000 # corrected; this appears to be incorrect when rearranged from supplement eqn 19b
  d13Cpf[i] <- ((11.98 - 0.12*tempC[i])/1000 + 1)*(d13C.co2g[i] + 1000) - 1000 # corrected; this appears to be incorrect when rearranged from supplement eqn 19b


  # Calculate instantaneous growth rate (mu,i) from [PO4] and rmean
  mui[i] <- coeff.po4*po4[i] + coeff.rm*rm[i] + mui.y.int
  # Calculate Qs (the co2 flux into the cell per unit surface area of the cell membrane)
  Qr[i] <- mui[i]/log(2) * gam.c[i]
  Qs[i] <- Qr[i] / (4*3.141593*(rm[i]^2))
  #	Calculate b in uM from Qs, rmean, rK, DT, Pc, eps.f and eps.d using eqn. 15 of Rau et al (1996)
  b[i] <- ((eps.f - eps.d) * Qs[i] * ((rm[i]/((1 + rm[i]/rk[i])*DT[i])) + 1/P.c)) * 10^3
  # Calculate eps.p from [CO2](aq), eps.f and b
  eps.p[i] <- eps.f - (b[i]/(co2[i]*10^3))
  # Calculate d13C.biomass from d13C.co2 and epsilon.p (eqn 2)
  d13C.biomass[i] <- ((d13C.co2[i] + 1000) / ((eps.p[i] / 1000) + 1)) - 1000
  # Calculate d13C.marker from d13C.biomass (eqn 22)
  d13Cmarker[i] <- ((d13C.biomass[i] + 1000) / ((eps.bob/1000)+1)) - 1000

}
############################################################################################ 


# Environmental prior distributions 
############################################################################################   
for (i in 1:length(d13Cmarker.data)){
# Temperature (degrees C)
tempC[i] ~ dnorm(tempC.m, tempC.p)
# Salintiy (ppt)
sal[i] ~ dnorm(sal.m, sal.p)T(0,)
# pCO2 (uatm)
pco2[i] ~ dnorm(pco2.m, pco2.p)T(125,375)
# d13C of aqueous CO2 (per mille)
d13C.co2[i] ~ dnorm(d13C.co2.m, d13C.co2.p)
# Concentration of phosphate (PO4; umol/kg)
po4[i] ~ dnorm(po4.m[site.index[i]], po4.p)T(0,2)
# Mean cell radius (m)
rm[i] ~ dnorm(rm.m, rm.p)T(0,5e-6)
}
############################################################################################   
}


