library(rjags)
library(reshape)
library(ggplot2)
library(SESYNCBayes)
set.seed(5)
y=N2OEmission
y.n.sites = length(unique(y$group))
qplot(n.input, emission, data=y, color =  group)

data = list(
  y.emission = log(y$emission),
  y.n.input = log(y$n.input) - mean(log(y$n.input)) #center the data to speed convergence and aid in interpretation. Can recover 0 intercept if needed.
)

inits = list(
  list(
    alpha = 0,
    beta = .5,
    sigma = 50
  ),
  list(
    alpha = 1,
    beta = 1.5,
    sigma = 10
  )
)


####Pooled model
{
sink("Pooled")
cat("
model{
#priors
alpha ~ dnorm(0,.0001)
beta ~ dnorm(0,.0001)
sigma ~ dunif(0,100)
tau.reg <- 1/sigma^2
#likelihood
 for(i in 1:length(y.emission)){
    mu[i] <- alpha + beta * y.n.input[i]
    y.emission[i] ~ dnorm(mu[i], tau.reg)
    #simulated data for posterior predictive checks
    y.emission.sim[i] ~ dnorm(mu[i], tau.reg) 
    sq.error.data[i] <- (y.emission[i]-mu[i])^2
    sq.error.sim[i] <- (y.emission.sim[i] - mu[i])^2
 }
#Bayesian P values
sd.data <- sd(y.emission)
sd.sim <- sd(y.emission.sim)
p.sd <- step(sd.sim-sd.data)

mean.data <- mean(y.emission)
mean.sim  <- mean(y.emission.sim)
p.mean <- step(mean.sim - mean.data)

discrep.data <- sum(sq.error.data)
discrep.sim <- sum(sq.error.sim)
p.discrep <- step(discrep.sim - discrep.data)

min.data <- min(y.emission)
min.sim <- min(y.emission.sim)
p.min <-step(min.sim-min.data)
}
    
",fill=TRUE)
sink()
}

n.update=10000
n.iter=10000
n.update = 3000
jm.pooled = jags.model("Pooled", data=data, n.adapt = 3000, inits=inits, n.chains=length(inits))

update(jm.pooled, n.iter = n.update)
#zc.pooled = coda.samples(jm.pooled, variable.names = c("alpha", "beta", "sigma"), n.iter=n.iter)
zj.pooled = jags.samples(jm.pooled, variable.names = c("alpha", "beta", "sigma", "p.sd", "p.mean", "p.discrep","p.min", "y.emission.sim"), n.iter=n.iter)
zj.pooled$p.sd

zj.pooled$p.mean

zj.pooled$p.discrep

zj.pooled$p.min

hist(data$y.emission, breaks=20, freq=FALSE, main="Simulated and real data", xlab=expression(paste("log(", N0[2], " emission)")), cex.lab=1.2) #note that this is the log transformed data
lines(density(zj.pooled$y.emission.sim), col="red")
legend(-10,.2,"simulated", col="red", lty="solid")
