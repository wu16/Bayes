#Exercise 1 Writing a DAG
#There is no x because we are assuming it is measured without error.

#Exercise 2 
#use priors that are informative whenever possible.

# Exercise 3 Using for loops
for(i in 1:5){
  b[i] ~ dnorm(0, .000001)
}

#Exercise 4
sink("LogisticJAGS.R")
cat(" 
model{
  # priors
  K ~ dunif(0, 4000)
  r ~ dunif (0, 2)
  sigma ~ dunif(0, 100) 
  tau <- 1/sigma^2
  
  # likelihood
  for(i in 1:n){
    mu[i] <- r - r/K * x[i]
    y[i] ~ dnorm(mu[i], tau)
  }
} 
",fill = TRUE)
sink()


rm(list = ls())
library(rjags)
library(SESYNCBayes)
set.seed(1)

Logistic <- Logistic[order(Logistic$PopulationSize),]

inits = list(
  list(K = 1500, r = .2, sigma = 1),
  list(K = 1000, r = .15, sigma = 5),
  list(K = 900, r = .3, sigma = 10))

data = list(
  n = nrow(Logistic),
  x = as.double(Logistic$PopulationSize),
  y = as.double(Logistic$GrowthRate))

n.adapt = 5000
n.update = 10000
n.iter = 10000

jm = jags.model("LogisticJAGS.R", data = data, inits = inits, n.chains = length(inits), n.adapt = n.adapt)
update(jm, n.iter = n.update)
zm = coda.samples(jm, variable.names = c("K", "r", "sigma", "tau"), n.iter = n.iter, n.thin = 1)

#Exercise 5
initFunc <- function (){
return(list(
  K = runif(1, 10, 2000),
  r = runif(1, .1, 1.6),
  sigma = runif(1, 1, 80)))}

initFunc()
initFunc()

library(parallel)
detectCores()
cl <- makeCluster(3) # Here we use three cores
clusterExport(cl, c("data", "initFunc", "n.adapt", "n.update", "n.iter")) 

ptm <- proc.time()
out <- clusterEvalQ(cl, {
  library(rjags)
  set.seed(1)
  jm = jags.model("LogisticJAGS.R", data = data, inits = initFunc(), 
  n.chains = 1, n.adapt = n.adapt)
  update(jm, n.iter = n.update)
  zmCore = coda.samples(jm, variable.names = c("K", "r", "sigma", "tau"), 
  n.iter = n.iter, thin = 1)
  return(as.mcmc(zmCore))
})
ParallelTime <- proc.time() - ptm
ParallelTime
stopCluster(cl)
zmP <- mcmc.list(out)

#We rerun the model sequentially and use proc.time again for comparison
ptm <- proc.time()
jm = jags.model("LogisticJAGS.R", data = data, inits = inits, n.chains = length(inits), n.adapt = n.adapt)
update(jm, n.iter = n.update)
zm = coda.samples(jm, variable.names = c("K", "r", "sigma", "tau"), n.iter = n.iter, n.thin = 1)
SequentialTime <- proc.time() - ptm
SequentialTime

# Exercise 6: Summarizing coda objects
a <- signif(as.data.frame(summary(zm)$stat[, 1:2]), digits = 3)
b <- signif(as.data.frame(summary(zm)$quantile[, c(1, 3, 5)]), digits = 3)
LogisticParameters <- cbind(rownames(a), a, b)
rownames(LogisticParameters) <- c()
names(LogisticParameters) <- c("parameter", "mean", "standardDeviation", "lower95", "median", "upper95")
LogisticParameters

write.csv(LogisticParameters, file = "LogisticParameters.csv")

# Exercise 7: Understanding coda objects
n.adapt = 500
n.update = 500
n.iter = 20

jm.short = jags.model("LogisticJAGS.R", data = data, inits = inits, n.chains = length(inits), n.adapt = n.adapt)
update(jm.short, n.iter = n.update)
zm.short = coda.samples(jm.short, variable.names = c("K", "r", "sigma", "tau"), n.iter = n.iter, n.thin = 1)
zm.short[[2]][3,3]
zm.short[[1]][,2]
zm.short