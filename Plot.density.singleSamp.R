

# Density plots of output for each sample
############################################################################################
hist(inv.out.348ka$BUGSoutput$sims.list$pco2, probability = TRUE, xlim=c(80,420))
lines(density(inv.out.348ka$BUGSoutput$sims.list$pco2), col = "blue")
abline(v=mean(inv.out.348ka$BUGSoutput$sims.list$pco2), col="blue")

hist(inv.out.403ka$BUGSoutput$sims.list$pco2, probability = TRUE, xlim=c(80,420))
lines(density(inv.out.403ka$BUGSoutput$sims.list$pco2), col = "red")
abline(v=mean(inv.out.403ka$BUGSoutput$sims.list$pco2), col="red")

u.x <- seq(80, 420, length=100)
u.y <- dunif(u.x, min=pco2.l, max=pco2.u)
plot(u.x, u.y, type = "l", ylim=c(0, 0.03))
lines(density(inv.out.348ka$BUGSoutput$sims.list$pco2), col = "blue")
lines(density(inv.out.403ka$BUGSoutput$sims.list$pco2), col = "red")

hist(inv.out.348ka$BUGSoutput$sims.list$rm, probability = TRUE, xlim=c(0,5*10^-6))
lines(density(inv.out.348ka$BUGSoutput$sims.list$rm), col = "blue")
abline(v=mean(inv.out.348ka$BUGSoutput$sims.list$rm), col="blue")

hist(inv.out.403ka$BUGSoutput$sims.list$rm, probability = TRUE, xlim=c(0,5*10^-6))
lines(density(inv.out.403ka$BUGSoutput$sims.list$rm), col = "red")
abline(v=mean(inv.out.403ka$BUGSoutput$sims.list$rm), col="red")

hist(inv.out.348ka$BUGSoutput$sims.list$po4, probability = TRUE, xlim=c(0,4))
lines(density(inv.out.348ka$BUGSoutput$sims.list$po4), col = "blue")
abline(v=mean(inv.out.348ka$BUGSoutput$sims.list$po4), col="blue")

hist(inv.out.403ka$BUGSoutput$sims.list$po4, probability = TRUE, xlim=c(0,4))
lines(density(inv.out.403ka$BUGSoutput$sims.list$po4), col = "red")
abline(v=mean(inv.out.403ka$BUGSoutput$sims.list$po4), col="red")

hist(inv.out.348ka$BUGSoutput$sims.list$b, probability = TRUE, xlim=c(0,200))
lines(density(inv.out.348ka$BUGSoutput$sims.list$b), col = "blue")
abline(v=mean(inv.out.348ka$BUGSoutput$sims.list$b), col="blue")

hist(inv.out.403ka$BUGSoutput$sims.list$b, probability = TRUE, xlim=c(0,200))
lines(density(inv.out.403ka$BUGSoutput$sims.list$b), col = "red")
abline(v=mean(inv.out.403ka$BUGSoutput$sims.list$b), col="red")

hist(inv.out.348ka$BUGSoutput$sims.list$tempC, probability = TRUE, xlim=c(10,40))
lines(density(inv.out.348ka$BUGSoutput$sims.list$tempC), col = "blue")
abline(v=mean(inv.out.348ka$BUGSoutput$sims.list$tempC), col="blue")

hist(inv.out.403ka$BUGSoutput$sims.list$tempC, probability = TRUE, xlim=c(10,40))
lines(density(inv.out.403ka$BUGSoutput$sims.list$tempC), col = "red")
abline(v=mean(inv.out.403ka$BUGSoutput$sims.list$tempC), col="red")
############################################################################################

print(mean(inv.out.2ka$BUGSoutput$sims.list$pco2))
print(mean(inv.out.2ka$BUGSoutput$sims.list$tempC))
print(mean(inv.out.40.8ka$BUGSoutput$sims.list$pco2))
print(mean(inv.out.40.8ka$BUGSoutput$sims.list$tempC))
print(mean(inv.out.348ka$BUGSoutput$sims.list$pco2))
print(mean(inv.out.348ka$BUGSoutput$sims.list$tempC))
print(mean(inv.out.403ka$BUGSoutput$sims.list$pco2))
print(mean(inv.out.403ka$BUGSoutput$sims.list$tempC))
