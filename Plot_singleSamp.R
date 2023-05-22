
# This script takes output from phytoDriver.singleSamp.R, extracts and plots the parameter output w/ 95% CI error bars
# Dustin T. Harper
# 5 May, 2023


###########################################################################################
# Load libraries
###########################################################################################
library(tidyverse)

###########################################################################################
# Extract data for plotting
###########################################################################################
rm1 <- inv.out.2ka$BUGSoutput$summary["rm",]
rm2 <- inv.out.40.8ka$BUGSoutput$summary["rm",]
rm3 <- inv.out.348ka$BUGSoutput$summary["rm",]
rm4 <- inv.out.403ka$BUGSoutput$summary["rm",]
rm.out <- data.frame(rbind(rm1, rm2, rm3, rm4))

po41 <- inv.out.2ka$BUGSoutput$summary["po4",]
po42 <- inv.out.40.8ka$BUGSoutput$summary["po4",]
po43 <- inv.out.348ka$BUGSoutput$summary["po4",]
po44 <- inv.out.403ka$BUGSoutput$summary["po4",]
po4.out <- data.frame(rbind(po41, po42, po43, po44))

d13C.co21 <- inv.out.2ka$BUGSoutput$summary["d13C.co2",]
d13C.co22 <- inv.out.40.8ka$BUGSoutput$summary["d13C.co2",]
d13C.co23 <- inv.out.348ka$BUGSoutput$summary["d13C.co2",]
d13C.co24 <- inv.out.403ka$BUGSoutput$summary["d13C.co2",]
d13C.co2.out <- data.frame(rbind(d13C.co21, d13C.co22, d13C.co23, d13C.co24))

sal1 <- inv.out.2ka$BUGSoutput$summary["sal",]
sal2 <- inv.out.40.8ka$BUGSoutput$summary["sal",]
sal3 <- inv.out.348ka$BUGSoutput$summary["sal",]
sal4 <- inv.out.403ka$BUGSoutput$summary["sal",]
sal.out <- data.frame(rbind(sal1, sal2, sal3, sal4))

tempC1 <- inv.out.2ka$BUGSoutput$summary["tempC",]
tempC2 <- inv.out.40.8ka$BUGSoutput$summary["tempC",]
tempC3 <- inv.out.348ka$BUGSoutput$summary["tempC",]
tempC4 <- inv.out.403ka$BUGSoutput$summary["tempC",]
tempC.out <- data.frame(rbind(tempC1, tempC2, tempC3, tempC4))

pco21 <- inv.out.2ka$BUGSoutput$summary["pco2",]
pco22 <- inv.out.40.8ka$BUGSoutput$summary["pco2",]
pco23 <- inv.out.348ka$BUGSoutput$summary["pco2",]
pco24 <- inv.out.403ka$BUGSoutput$summary["pco2",]
pco2.out <- data.frame(rbind(pco21, pco22, pco23, pco24))

###########################################################################################
# 95% CI plots for each time step 
###########################################################################################

ages.prox <- c(2, 40.8, 348.4, 403.1)
xmax <- 500
xmin <- -10

# Plot parms of interest

ggplot() + 
  geom_errorbar(data = rm.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray", width=0) +
  geom_point(data = rm.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(0,5e-6) +
  xlim(xmin,xmax) +
  labs(x= "age (ka)", y = "mean radius") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = d13C.co2.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray", width=0) +
  geom_point(data = d13C.co2.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(-10,-5) +
  xlim(xmin,xmax) +
  labs(x= "age (ka)", y = "d13Cco2") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = po4.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray", width=0) +
  geom_point(data = po4.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(0,3) +
  xlim(xmin,xmax) +
  labs(x= "age (ka)", y = "[PO4]") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = sal.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray", width=0) +
  geom_point(data = sal.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(32,38) +
  xlim(xmin,xmax) +
  labs(x= "age (ka)", y = "salinity (ppt)") +  
  theme_bw()

ggplot() + 
  geom_errorbar(data = tempC.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray", width=0) +
  geom_point(data = tempC.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(15,35) +
  xlim(xmin,xmax) +
  labs(x= "age (ka)", y = "temp (C)") +
  theme_bw()

ggplot() + 
  geom_errorbar(data = pco2.out, mapping = aes(x=ages.prox, y=mean, ymin=X2.5., ymax=X97.5.), color="gray", width=0) +
  geom_point(data = pco2.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(150,350) +
  xlim(xmin,xmax) +
  labs(x= "age (ka)", y = "pCO2 (atm)") +   
  theme_bw()




