# set random number generator
set.seed(3)

# simulate the data
y <- rpois(50, lambda = 6.4)

# histogram of the data
hist(y, freq = FALSE, breaks = 10, main = "Histogram of data", col = "gray")

# custom function for improved histogram of the data
library(plyr)

discrete_hist <- function(y) {
  z <- count(y)
  z$dens <- z$freq/sum(z$freq)
  plot(z$x, z$dens, type = "h", ylab = "Probability", xlab = "y", main = "Improved histogram of data", 
    frame = FALSE, xaxt = "n", lwd = 3, col = "blue")
  x <- seq(min(z$x), max(z$x), 1)
  axis(side = 1, at = x, labels = x)
}
discrete_hist(y)

# prior mean and standard deviation
mu.prior <- 10.2
sigma.prior <- 0.5

step <- .01
theta <- seq(0, 15, step)

prior <- function(theta, mu = mu.prior, sigma = sigma.prior) dgamma(theta, mu^2 / sigma^2, mu / sigma^2)

plot(theta, prior(theta), type = "l", ylab = expression(paste("[", theta, "]")), xlab = expression(theta),
  main = "Prior", xlim = c(5, 15))

sd(rgamma(100000, mu.prior^2 / sigma.prior^2, mu.prior / sigma.prior^2))
## [1] 0.5013619
mean(rgamma(100000, mu.prior^2 / sigma.prior^2, mu.prior / sigma.prior^2))
## [1] 10.20079

like <- function(theta, y){
  L = rep(0, length(theta))
    for(i in 1:length(theta)) L[i] = prod(dpois(y, theta[i], log = FALSE))
    return(L)
} 

plot(theta, like(theta, y = y), type = "l", xlim = c(5, 15), main = "Likelihood", xlab = expression(theta), 
  ylab = expression(paste("[y|", theta, "]")))

joint = function(theta) like(theta, y = y) * prior(theta)

plot(theta, joint(theta), type = "l",  main = "Joint", xlim = c(5, 15), xlab = expression(theta),
  ylab = expression(paste("[y|", theta, "] x [", theta, "]")))

(Py <- sum(like(theta, y = y) * prior(theta) * step))
## [1] 2.140143e-58

p.theta <- joint(theta) / Py

plot(theta, p.theta, typ = "l", xlim = c(5, 15), main = "Posterior", xlab = expression(theta), 
  ylab = expression(paste("[ ", theta, " | y]")))

par(mfrow = (c(2, 3)))
plot(theta, prior(theta), type = "l", ylab = expression(paste("[", theta, "]")), xlab = expression(theta),
  main = "Prior", xlim = c(5, 15))

hist(y, freq = FALSE, breaks = 10, main = "Histogram of data")
discrete_hist(y = y)

plot(theta, like(theta, y = y), type = "l", main = "Likelihood", xlim = c(5, 15), xlab = expression(theta),
  ylab = expression(paste("[y|", theta, "])")))
plot(theta, joint(theta), type = "l", main = "Likelihood", xlim = c(5, 15), xlab = expression(theta),
  ylab = expression(paste("[y | ", theta, "]) x [", theta, "]")))
plot(theta, p.theta, type = "l", xlim = c(5, 15), main = "Posterior", xlab = expression(theta),
  ylab = expression(paste("[ ", theta, " | y]")))
