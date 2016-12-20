#Scenario 1

aPrior <- 1
bPrior <- 1 
y <- 3
n <- 10
aPost <- aPrior + y
bPost <- bPrior + n - y

prior <- function(x) dbeta(x, aPrior, bPrior)
post <- function(x) dbeta(x, aPost, bPost)

plot(prior, 0, 1, ylim=c(0, 3), main = "Beta prior (dash) and  beta posterior (red)", xlab = expression(phi), ylab = "Density", lty = 2, lwd = 2)
plot(post, 0, 1, add = TRUE, col = "red", lwd = 2)

x <- seq(0,1,0.01)
plot(x, prior(x), ylim=c(0,3), type="l")
lines((x, post(x), col="red")

bounds <- c(0.025, 0.975)
bounds95 <- qbeta(bounds, aPost, bPost)
bounds95

1 - pbeta(.5, aPost, bPost)

# Scenario2

aPrior <- 4
bPrior <- 8 
y <- 38
n <- 100

aPost <- aPrior + y
bPost <- bPrior + n - y

prior <- function(x) dbeta(x, aPrior, bPrior)
post <- function(x) dbeta(x, aPost, bPost)

postFlat <- function(x) dbeta(x, 1 + y, 1 + n - y)

plot(prior, 0, 1, ylim=c(0, 10), main = "Beta prior (dash) and beta posterior (red, blue)", xlab = expression(phi), ylab = "Density", lty = 2, lwd = 2)
plot(post, 0, 1, add = TRUE, col = "red", lwd = 2)
plot(postFlat, 0, 1, add = TRUE, col = "blue", lwd = 2)

bounds <- c(0.025, 0.975)
bounds95 <- qbeta(bounds, aPost, bPost)
bounds95
1 - pbeta(.5, aPost, bPost)

#Scenario 3
mu <- 3
sigma <- 1.5

aPrior <- mu^2/sigma^2
bPrior <- mu/sigma^2

y <- c(15, 8, 6, 11, 4 )

aPost <- aPrior + sum(y)
bPost <- bPrior + length(y)

prior <- function(x) dgamma(x, aPrior, bPrior)
post <- function(x) dgamma(x, aPost, bPost)
plot(prior, 0, 20, main = "Gamma prior (dash) and gamma posterior (red)", xlab = expression(lambda), ylab = "Density", lty = 2, lwd = 2, ylim = c(0, 0.4))
plot(post, 0, 20, add = TRUE, col = "red", lwd = 2)

bounds <- c(0.025, 0.975)
bounds95 <- qgamma(bounds, aPost, bPost)
bounds95

#Scenario 4
setwd("F:/Wu/Bayes3/SESYNCBayes/")
install.packages("Packages/SESYNCBayes_0.1.0.tar.gz", repos = NULL, type = "source")

library(SESYNCBayes)
data(SolsticeTemp)

varKnown <- 15^2
n <- 50
varPrior <- 12^2
muPrior <- 30
y <- sum(SolsticeTemp$temp)
muPost <-((muPrior/varPrior) + y/varKnown)/((1/varPrior) + (n/varKnown))
varPost <- 1/((1/varPrior) + (n/varKnown))
prior <- function(x) dnorm(x, muPrior, varPrior^.5)
post <- function(x) dnorm(x, muPost, varPost^.5)
plot(prior, 0, 60, main = "Normal prior (dash) and normal posterior (red)", xlab = expression(mu),
ylab = "Density", lty = 2, ylim = c(0, 0.2), lwd = 2)
plot(post, 0, 60, add = TRUE, col = "red", lwd = 2)

bounds <- c(0.025, 0.975)
bounds95 <- qnorm(bounds, muPost, varPost^.5)
bounds95