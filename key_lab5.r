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

#We will also need a function to link the sequential indices used in JAGS to the groups (fertilizer and site) in the data. 
group_from_index = function(group, group.index, output ){
  #group is a vector of group names or numbers
  #group.index is vector of sequential indices to group names or numbers.  It is a sequence of integers 1 to length(group)
  #output is a matrix or dataframe of output with number of rows = length(group). Each row contains statistics, etc for each group.
  a = unique(as.vector(group)) 
  b = unique(group.index)
  group.key=as.data.frame(t(rbind(a,b))) #columns containing indices paired with group name or number
  names(group.key)= c(names(as.data.frame(group)), names(as.data.frame(group.index))) 
  link.column.name = names(group.key)[2] #name of column for merging output data with groups
  output2 = cbind(seq(1,nrow(output)),output) #give the output data sequential index the same as 
  colnames(output2)[1]=link.column.name
  group.data=as.data.frame(merge(group.key, output2, by = link.column.name )) #merge the output with the groups
  return(group.data)
}

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
    sigma ~ dunif(0,100)
    tau.reg <- 1/sigma^2
    ###priors
    for(j in 1:y.n.sites){
        alpha[j] ~ dnorm(mu.alpha,tau.alpha)
      }
    beta ~ dnorm(0,.0001)
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

# Intercepts vary with carbon level in site soils and slopes vary with fertilizer type

#######Hierarchical model, site level intercept predicted from carbon concentration covariate and slope varying with fertilizer type. 

w$mean=w$mean/100  #transform % to proportion
data = list(
  y.emission = log(y$emission),
  y.n.input = log(y$n.input)-mean(log(y$n.input)), #center the data to speed convergence and aid in interpretation
  y.group=  y$group.index,
  y.fert = y$fert.index,
  y.n.sites = length(unique(y$group)),
  y.n.fert = length(unique(y$fertilizer)),
  w = log(w$mean/(1-w$mean))   #logit of w$mean
)
y.n.sites = length(unique(y$group))
y.n.fert = length(unique(y$fertilizer))
inits = list(
  list(
    alpha = rep(0,y.n.sites),
    beta = rep(.5,y.n.fert),
    sigma = 50,
    sigma.alpha = 10,
    eta = .2,
    kappa = .5
  ),
  list(
    alpha = rep(-.2,y.n.sites),
    beta = rep(1.5, y.n.fert),
    sigma = 10,
    sigma.alpha = 20,
    eta = .2,
    kappa = 5
  )
)

{
sink("Hier_2")
cat("
    model{
    #priors for within site model######
    sigma ~ dunif(0,200)
    tau.reg <- 1/sigma^2
    
    #priors for intercept model#######
    kappa ~ dnorm(0,.00001)
    eta ~ dnorm(0, .000001)
    sigma.alpha ~ dunif(0,200)
    tau.alpha <- 1/sigma.alpha^2
    #hyper priors for slope model
    mu.beta ~ dnorm(0,.00001)
    sigma.beta ~ dunif(0,200)
    tau.beta <- 1/sigma.beta
  

    #likelihood for data, note that data are on log scale in data statement on R side
    for(i in 1:length(y.emission)){
      mu[i] <- alpha[y.group[i]] + beta[y.fert[i]] * y.n.input[i]
      y.emission[i] ~ dnorm(mu[i], tau.reg)
    }
    # carbon model for intercept
  for(j in 1:y.n.sites){
     #use normal because data are centered
      mu.alpha[j] <- kappa + eta *w[j]
      alpha[j] ~ dnorm(mu.alpha[j],tau.alpha)
  }
  #Allow slope to vary by fertilizer type
  for(k in 1:y.n.fert){
    beta[k] ~ dnorm(mu.beta, tau.beta)
  }
 } #end of model
    
    ",fill=TRUE)
sink()
  
}

n.update=50000
n.iter=25000

jm.hier2 = jags.model("Hier_2", data=data, n.adapt = 3000, inits=inits, n.chains=length(inits))

jm.hier2

update(jm.hier2, n.iter = n.update)
#You should run diagnostics on zc.hier coda object for intecepts and slopes but I eliminated them here to make output more compact
zc.hier2 = coda.samples(jm.hier2, variable.names = c("sigma","eta", "kappa", "mu.beta", "sigma.beta"), n.iter=n.iter)
zj.hier2 = jags.samples(jm.hier2, variable.names = c(variable.names = c("alpha", "beta", "sigma","eta", "kappa")), n.iter=n.iter)
gelman.diag(zc.hier2)
summary(zc.hier2)

slopes=t(summary(zj.hier2$beta,quantile, c(.025,.5,.957))$stat) #transpose, the t( ), is important to make next function work


#slopes as function of fetilizer type
group.data=group_from_index(group=y$fertilizer,group.index=y$fert.index,output=slopes)

#table with medians and credible intervals for slopes by fertilizer type
group.data

names(group.data)[3:5]=c("lower", "median", "upper")
library(ggplot2)
ggplot( group.data, aes(x = group, y = median)) +    geom_bar(position = position_dodge(), stat="identity", fill="red")  +   geom_errorbar(aes(ymin=lower, ymax=upper)) +   ggtitle("Medians of slopes by fertilizer type with 95% credible intervals") + # plot title 
  labs(x="Fertilizer", y=expression(beta)) +
  theme_bw() + # remove grey background (because Tufte said so)
  theme(panel.grid.major = element_blank()) # remove x and y major grid lines (because Tufte said so)

#Slope and intercepts vary by site
data = list(
  y.emission = log(y$emission),
  y.n.input = log(y$n.input)-mean(log(y$n.input)), #center the data to speed convergence and aid in interpretation-- there is no such thing as soil with 0 carbon
  y.group=  y$group.index,
  y.fert = y$fert.index,
  y.n.sites = length(unique(y$group)),
  y.n.fert = length(unique(y$fertilizer))
)
y.n.sites = length(unique(y$group))
B = matrix(nrow=y.n.sites, ncol=2)
B[,1]=0
B[,2]=1.5
inits = list(
  list(
    B=B,
    sigma = 50,
    mu.alpha = 0,
    mu.beta = 1.5,
    sigma.alpha = 10,
    sigma.beta = 10,
    rho=-.5
  ),
  list(
    B=B*.5,
    sigma = 20,
    mu.alpha = -.2,
    mu.beta = .8,
    sigma.alpha = 50,
    sigma.beta = 50,
    rho=.5
  )
)

{
sink("Hier_3")
cat("
    model{
    #priors for within site model######
    sigma ~ dunif(0,200)
    tau.reg <- 1/sigma^2
    
    #likelihood for data, note that data are on log scale in data statement on R side
    for(i in 1:length(y.emission)){
      mu[i] <- alpha[y.group[i]] + beta[y.group[i]] * y.n.input[i]
      y.emission[i] ~ dnorm(mu[i], tau.reg)
    }
    # Model for group intercept and slope:
    for(j in 1:y.n.sites){
        alpha[j] <- B[j,1]  #group level intercept
        beta[j]  <- B[j,2]  #group level slope
        B[j,1:2] ~ dmnorm(B.hat[j,1:2], Tau.B)  
        B.hat[j,1] <- mu.alpha  #required by JAGS syntax
        B.hat[j,2] <- mu.beta   #required by JAGS syntax
    }
    mu.alpha ~ dnorm(0,.0001)  #mean intercept
    mu.beta ~ dnorm(0, .0001)  #mean slope
    #Inverse of covariance matrix required by JAGS
    Tau.B[1:2,1:2] <- inverse(Sigma.B[1:2,1:2])
    #Elements of covariance matrix
    Sigma.B[1,1] <- sigma.alpha^2
    sigma.alpha ~ dunif(0,200)
    Sigma.B[2,2] <- sigma.beta^2
    sigma.beta ~ dunif(0,200)
    Sigma.B[1,2] <- rho*sigma.alpha*sigma.beta  # covariance is correlation coef. x product of variances
    Sigma.B[2,1] <- Sigma.B[1,2]
    rho ~ dunif(-1,1)
    } #end of model
    
    ",fill=TRUE)
sink()
}

n.update=50000
n.iter=10000

jm.hier3 = jags.model("Hier_3", data=data, n.adapt = 3000, inits=inits, n.chains=length(inits))

update(jm.hier3, n.iter = n.update)
#You should run diagnostics on zc.hier coda object for intecepts and slopes but I eliminated them here to make output more compact
zc.hier3 = coda.samples(jm.hier3, variable.names = c('mu.alpha', "mu.beta", "rho"), n.iter=n.iter)
zj.hier3 = jags.samples(jm.hier3, variable.names = c("alpha", "beta", 'mu.alpha', "mu.beta", "rho"), n.iter=n.iter)
gelman.diag(zc.hier3)

summary(zc.hier3)

#Make a vector to link sequential index (in ouput) to the group number.  Would work the same way if groups were character variables like names. t() is transpose. See group_from_index() function at top of file
slopes = t(summary(zj.hier3$beta, quantile, c(.025,.5,.975))$stat) #transpose is t() is important to make next function work
group.data=as.data.frame(group_from_index(group=y$group,group.index=y$group.index,output=slopes))
names(group.data)[3:5]=c("lower", "median", "upper")
library(ggplot2)
ggplot( group.data, aes(x = group, y = median)) +    geom_bar(position = position_dodge(), stat="identity", fill="red")  +   geom_errorbar(aes(ymin=lower, ymax=upper)) +   ggtitle("Medians of site-level slopes with 95% credible intervals") + # plot title 
  labs(x="Site", y=expression(beta)) +
  theme_bw() + # remove grey background (because Tufte said so)
  theme(panel.grid.major = element_blank()) # remove x and y major grid lines (because Tufte said so)

plot(density(zj.hier3$mu.beta), xlab = expression(mu[beta]), main="Posterior distribution of mean slope", cex.lab=1.25)



