library(tidyverse)
library(viridis)

fig.font <- "Helvetica"
fontsize.axislabels <- 11
fontsize.scalelabels <- 9
po4.seq <- seq(0,3,0.1)
rm.seq <- seq(0,5,0.1)



# multi-linear regression using culture data
cal.df <- read.csv('data/caldata_culture.csv')
cal.df <- subset(cal.df, is.na(cal.df$radius) | is.na(cal.df$po4) | cal.df$radius <= 10 & cal.df$po4 <= 2)

igr.model = lm(formula = cal.df$mui ~ cal.df$po4 + cal.df$r, data = cal.df)
igr.model.sum <- summary(igr.model)

po4.co.lr <- igr.model.sum$coefficients[2,1]
po4.se.lr <- igr.model.sum$coefficients[2,2]
r.co.lr <- igr.model.sum$coefficients[3,1]
r.se.lr <- igr.model.sum$coefficients[3,2]
y.int.lr <- igr.model.sum$coefficients[1,1]
y.int.se.lr <- igr.model.sum$coefficients[1,2]

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



# Bayesian inversion using culture data 
coeff <- read.csv("model_out/coeff_mat_culture.csv")

mlr.fun <- function(x, y) {final_value = mean(coeff[,4]) + x*mean(coeff[,2]) + y*mean(coeff[,3])} 

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



# Bayesian inversion using ice core data 
coeff <- read.csv("model_out/coeff_mat_ice.csv")
mlr.fun <- function(x, y) {final_value = mean(coeff[,4]) + x*mean(coeff[,2]) + y*mean(coeff[,3])}  

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



# Bayesian inversion using ice + culture data 
coeff <- read.csv("model_out/coeff_mat_ice_culture.csv")
mlr.fun <- function(x, y) {final_value = mean(coeff[,4]) + x*mean(coeff[,2]) + y*mean(coeff[,3])} 

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


