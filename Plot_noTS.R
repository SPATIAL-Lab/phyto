
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
tempC.out <- apply(inv.out[["BUGSoutput"]][["sims.list"]][["tempC"]], 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
pco2.out <- apply(inv.out[["BUGSoutput"]][["sims.list"]][["pco2"]], 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
sal.out <- apply(inv.out[["BUGSoutput"]][["sims.list"]][["sal"]], 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
d13C.co2.out <- apply(inv.out[["BUGSoutput"]][["sims.list"]][["d13C.co2"]], 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
po4.out <- apply(inv.out[["BUGSoutput"]][["sims.list"]][["po4"]], 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
rm.out <- apply(inv.out[["BUGSoutput"]][["sims.list"]][["rm"]], 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

ages.prox <- prox.in$age

iccc <- read.csv("data/icecoreco2comp.csv")

###########################################################################################
# 95% CI plots for each time step 
###########################################################################################


xmax <- 800
xmin <- -10

# Plot parms of interest

ggplot() + 
  geom_errorbar(mapping = aes(x=ages.prox, y=rm.out[2,], ymin=rm.out[1,], ymax=rm.out[3,]), color="gray", width=0) +
  geom_point(mapping = aes(x=ages.prox, y=rm.out[2,]), color = "black") +
  ylim(0,5e-6) +
  xlim(xmax,xmin) +
  labs(x= "age (ka)", y = "mean radius") +
  theme_bw()

ggplot() + 
  geom_errorbar(mapping = aes(x=ages.prox, y=d13C.co2.out[2,], ymin=d13C.co2.out[1,], ymax=d13C.co2.out[3,]), color="gray", width=0) +
  geom_point(mapping = aes(x=ages.prox, y=d13C.co2.out[2,]), color = "black") +
  ylim(-10,-5) +
  xlim(xmax,xmin) +
  labs(x= "age (ka)", y = "d13Cco2") +
  theme_bw()

ggplot() + 
  geom_errorbar(mapping = aes(x=ages.prox, y=po4.out[2,], ymin=po4.out[1,], ymax=po4.out[3,]), color="gray", width=0) +
  geom_point(mapping = aes(x=ages.prox, y=po4.out[2,]), color = "black") +
  ylim(0,3) +
  xlim(xmax,xmin) +
  labs(x= "age (ka)", y = "[PO4]") +
  theme_bw()

ggplot() + 
  geom_errorbar(mapping = aes(x=ages.prox, y=sal.out[2,], ymin=sal.out[1,], ymax=sal.out[3,]), color="gray", width=0) +
  geom_point(mapping = aes(x=ages.prox, y=sal.out[2,]), color = "black") +
  ylim(32,38) +
  xlim(xmax,xmin) +
  labs(x= "age (ka)", y = "salinity (ppt)") +  
  theme_bw()

ggplot() + 
  geom_errorbar(mapping = aes(x=ages.prox, y=tempC.out[2,], ymin=tempC.out[1,], ymax=tempC.out[3,]), color="black", width=1) +
  geom_errorbar(mapping = aes(x=ages.prox, y=Uk.data*29.876, ymin=Uk.data-Uk.data.sd*29.876, 
                                              ymax=Uk.data+Uk.data.sd*29.876), color="red", linetype=2, width=0) +
  geom_point(mapping = aes(x=ages.prox, y=tempC.out[2,]), color = "black") +
  geom_point(mapping = aes(x=ages.prox, y=Uk.data*29.876), color = "red") +
  scale_y_continuous(sec.axis=sec_axis(~./29.876, name = "Uk'37"), limits=c(20,35))+
  theme_bw()+
  theme(axis.text.y.right = element_text(color = "red"), axis.line.y.right = element_line(color = "red"), axis.ticks.y.right = element_line(color = "red")) +
  xlim(xmax,xmin) +
  labs(x= "age (ka)", y = "temp (C)") 

ggplot() + 
  geom_errorbar(mapping = aes(x=ages.prox, y=pco2.out[2,], ymin=pco2.out[1,], ymax=pco2.out[3,]), color="gray", width=0) +
  geom_line(data=iccc, mapping = aes(x=age, y=iceco2), color = "red") +
  geom_point(mapping = aes(x=ages.prox, y=pco2.out[2,]), color = "black") +
  ylim(50,500) +
  xlim(xmax,xmin) +
  labs(x= "age (ka)", y = "pCO2 (uatm)") +   
  theme_bw()


# View inversion summary 
View(inv.out$BUGSoutput$summary)



