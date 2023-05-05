# Distribution scratch work

# gammas
x <- seq(1e12,1e14, by=1e11)
y<-dgamma(x,1e13,1)
plot(y)

x <- seq(0,250, by=1)
y<-dgamma(x,100, 10)
plot(y)

x <- seq(0,250, by=1)
y<-dgamma(x,0.01,0.1)
plot(y)

x <- seq(0,250, by=1)
y<-dgamma(x,0.1,1)
plot(y)

#betas
x <- seq(0,1, by=0.05)
y<-dbeta(x, 5,2)
plot(x,y)

