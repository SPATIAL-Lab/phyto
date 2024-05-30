
model {

# Likelihood function
############################################################################################ 
# Modern calibration data   
for (i in 1:length(mui.cd.data)){
  mui.cd.data[i] ~ dnorm(mui.cd[i], 1/(1e-6)^2)
}

for(i in 1:length(radius.cd)){
  radius.cd[i] ~ dnorm(rm.cd[i], 1/(0.25e-6)^2)
}

for(i in 1:length(po4.cd.data)){
  po4.cd.data[i] ~ dnorm(po4.cd[i], 1 /(0.05)^2)
}
############################################################################################ 

# mui coefficient priors
############################################################################################ 
# Coefficient for mu(i) - PO4 multi linear regression 
coeff.po4 ~ dnorm(po4.co.lr, 1/(po4.se.lr)^2)
# Coefficient for radius - PO4 multi linear regression 
coeff.rm ~ dnorm(r.co.lr, 1/(r.se.lr)^2)
# Coefficient for y intercept for mu(i) multi linear regression 
mui.y.int ~ dnorm(y.int.lr, 1/(y.int.se.lr)^2)
############################################################################################ 

# mui PSM - modern calibration data
for (i in 1:length(radius.cd)){
  # Mean cell radius (m)
  rm.cd[i] ~ dnorm(rm.m, rm.p)T(0,5e-6)

  po4.cd[i] ~ dnorm(po4.m.cd, po4.p)T(0,2)

  # Calculate instantaneous growth rate (mu,i) from [PO4] and rmean
  mui.cd[i] <- coeff.po4*po4.cd[i] + coeff.rm*rm.cd[i] + mui.y.int
}
############################################################################################ 
}


