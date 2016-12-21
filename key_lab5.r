library(rjags)
library(reshape)
library(ggplot2)
library(SESYNCBayes)
set.seed(5)
y=N2OEmission
w=SiteCarbon
w$mean=w$mean/100  #transform % to proportion
y.n.sites = length(unique(y$group))
head(y)
qplot(n.input, emission, data=y, color =  group)
qplot(n.input, emission, data=y, color =  fertilizer)

#Pooled model

# set up data and initial values
data = list(
  y.emission = log(y$emission),
  y.n.input = log(y$n.input)-mean(log(y$n.input)) #center the data to speed convergence and aid in interpretation.

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

# build model
{
sink("PooledJAGS.R")
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
 }

}
    
",fill=TRUE)
sink()
}

#compile the model

  n.adapt=3000
  n.update=5000
  n.iter= 5000
jm.pooled = jags.model(file="PooledJAGS.R", data=data, n.adapt = n.adapt, inits=inits, n.chains=length(inits))
update(jm.pooled, n.iter = n.update)
zc.pooled = coda.samples(jm.pooled, variable.names = c("alpha", "beta", "sigma"), n.iter=n.iter)
zj.pooled = jags.samples(jm.pooled, variable.names = c("alpha", "beta", "sigma"), n.iter=n.iter)
hist(zj.pooled$beta, xlab = expression(beta), breaks=50, main = "MCMC output")

summary(zc.pooled)
gelman.diag(zc.pooled)

# Intercepts for each site

data = list(
  y.emission = log(y$emission),
  y.n.input = log(y$n.input) - mean(log(y$n.input)), #center the data to speed convergence and aid in interpretation. Can recover 0 intercept if needed.
  y.group = y$group.index,  #use j to index groups
  y.n.sites = length(unique(y$group))
)

inits = list(
  list(
    alpha = rep(0,y.n.sites),
    beta = .5,
    sigma = 50,
    mu.alpha= 0,
    sigma.alpha = 10
  ),
  list(
    alpha = rep(1,y.n.sites),
    beta = 1.5,
    sigma = 10,
    mu.alpha= -2,
    sigma.alpha = 20
  )
)

{ #note this opening { and the closing } are needed by R markdown but not by R
####Hierarchical model, site level intercept, no site covariate
sink("Hier_1")
cat("
    model{
    ##hyperpriors
    mu.alpha ~ dnorm(0,.00001)
    sigma.alpha ~ dunif(0,200) #notated varsigma in model documentation
    tau.alpha <- 1/sigma.alpha^2

    ###priors
    for(j in 1:y.n.sites){
        alpha[j] ~ dnorm(mu.alpha,tau.alpha)
      }
    beta ~ dnorm(0,.0001)
    sigma ~ dunif(0,100)
    tau.reg <- 1/sigma^2

    ####
    #likelihood
    for(i in 1:length(y.emission)){
        mu[i] <- alpha[y.group[i]] + beta * y.n.input[i]
        y.emission[i] ~ dnorm(mu[i], tau.reg)
    }
    
    }
    
    ",fill=TRUE)
sink()
}

n.update=50000
n.iter=25000

jm.hier1 = jags.model("Hier_1", data=data, n.adapt = 3000, inits=inits, n.chains=length(inits))

update(jm.hier1, n.iter = n.update)
#You would wat to include alphas in check for convergnce but I eliminated them here to make output more compact.
zc.hier1 = coda.samples(jm.hier1, variable.names = c("sigma","beta", "mu.alpha", "sigma.alpha"), n.iter=n.iter)

summary(zc.hier1)

gelman.diag(zc.hier1)