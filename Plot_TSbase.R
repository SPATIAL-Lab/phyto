
# This script takes output from the phytoDriver, extracts and plots the output for  
# each time step w/ 95% CI error bars
# Dustin T. Harper
# 2 May, 2023


###########################################################################################
# Load libraries
###########################################################################################
library(tidyverse)


###########################################################################################
# Extract data for plotting
###########################################################################################

step.vector <- seq(1, n.steps, by=1)
parms.out <- inv.out$BUGSoutput$summary
parms.out <- data.frame(parms.out)
ages.prox <- rev(ages.prox)

tempC.v <- paste("tempC[", step.vector, "]", sep="")
tempC.out <- parms.out[c(tempC.v),]
tempC.out <- data.frame(ages.prox,tempC.out)

pco2.v <- paste("pco2[", step.vector, "]", sep="")
pco2.out <- parms.out[c(pco2.v),]
pco2.out <- data.frame(ages.prox,pco2.out)

sal.v <- paste("sal[", step.vector, "]", sep="")
sal.out <- parms.out[c(sal.v),]
sal.out <- data.frame(ages.prox,sal.out)

d13C.co2.v <- paste("d13C.co2[", step.vector, "]", sep="")
d13C.co2.out <- parms.out[c(d13C.co2.v),]
d13C.co2.out <- data.frame(ages.prox,d13C.co2.out)

po4.v <- paste("po4[", step.vector, "]", sep="")
po4.out <- parms.out[c(po4.v),]
po4.out <- data.frame(ages.prox,po4.out)

rm.v <- paste("rm[", step.vector, "]", sep="")
rm.out <- parms.out[c(rm.v),]
rm.out <- data.frame(ages.prox,rm.out)

iccc <- read.csv("data/icecoreco2comp.csv")

###########################################################################################
# 95% CI plots for each time step 
###########################################################################################


xmax <- 800
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
  geom_line(data=iccc, mapping = aes(x=age, y=iceco2), color = "red") +
  geom_point(data = pco2.out, mapping = aes(x=ages.prox, y=mean), color = "black") +
  ylim(150,350) +
  xlim(xmin,xmax) +
  labs(x= "age (ka)", y = "pCO2 (uatm)") +   
  theme_bw()


# View inversion summary 
View(inv.out$BUGSoutput$summary)




