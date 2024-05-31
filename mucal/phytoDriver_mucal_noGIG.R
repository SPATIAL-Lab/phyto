

# Phytoplankton forward PSM driver for generating matrix of mui coefficients and y intercepts
# Incorporates only modern observational data to calibrate mui = f(radius, po4)
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


# Prior distributions for cell radius and po4
############################################################################################    
# Concentration of phosphate (PO4; umol/kg)
po4.m.cd = 1.5
po4.p = 1/0.5^2 # 1/0.1^2

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
                 "radius.cd" = cal.df$radius,
                 "po4.cd.data" = cal.df$po4,
                 "mui.cd.data" = cal.df$mui,
                 "po4.m.cd" = po4.m.cd,
                 "po4.p" = po4.p,
                 "rm.m" = rm.m,
                 "rm.p" = rm.p) 
############################################################################################


# Parameters to save as output 
############################################################################################
parms = c("po4.cd", "rm.cd", "mui.cd", "coeff.po4", "coeff.rm", "mui.y.int")
############################################################################################


# Run the inversion using jags 
############################################################################################
inv.out = jags.parallel(data = data.pass, model.file = "phytoPSM_mucal_noGIG.R", parameters.to.save = parms,
                          inits = NULL, n.chains = 3, n.iter = 5e5,
                          n.burnin = 2e5, n.thin = 100)

############################################################################################
# Save model parameters governing mui relationship w/ size and po4
############################################################################################
coeff.mat <- cbind(inv.out[["BUGSoutput"]][["sims.list"]][["coeff.po4"]],
                   inv.out[["BUGSoutput"]][["sims.list"]][["coeff.rm"]],
                   inv.out[["BUGSoutput"]][["sims.list"]][["mui.y.int"]])
write.csv(coeff.mat, file = "model_out/coeff_mui_noGIG.csv")



View(inv.out$BUGSoutput$summary)



