## Correlated Data
library(MASS)
library(nimble, warn.conflicts = FALSE)

X <- seq(1,100,1)
n = 100
y <- mvrnorm(1, X*4, diag(n)*10)

# Simulated Data
tau <- 50.1
beta1 <- 2
beta2 <- 4
sigma1 <- 4
sigma2 <- 8

for(i in 1:n) {
  y[i] <- rnorm(1 , beta1 * X[i] * ( ( ( tau - i ) + abs(tau - i) ) / (2 * ( tau - i ) ) ) + beta2 *X[i] * 
                  ( ( ( i - tau ) + abs(i - tau) ) / (2 * ( i - tau ) ) ), 
                sigma1 * (((tau-i) + abs(tau-i)) / (2 * (tau-i)) + 
                  sigma2 * (((i-tau) + abs(i-tau)) / (2 * (i-tau))))) 
}
plot(X,y)
lines(X,y)


# Nimble Code
code <- nimbleCode({
  beta1 ~ dnorm(0, sd = 200)
  beta2 ~ dnorm(0, sd = 200) 
  tau ~ dunif(0,100)
  phiU ~ dunif(0,20)
  phiW ~ dunif(0,20)
  phiV ~ dunif(0,20)
  sigmaU[1:n,1:n] <- exp(-distmat[1:n,1:n]/phiU)
  sigmaV[1:n,1:n] <- exp(-distmat[1:n,1:n]/phiV)
  sigmaW[1:n,1:n] <- exp(-distmat[1:n,1:n]/phiW)
  U[1:n] ~ dmnorm(mu_zero[1:n],sigmaU[1:n,1:n])
  V[1:n] ~ dmnorm(mu_zero[1:n],sigmaV[1:n,1:n])
  W[1:n] ~ dmnorm(mu_zero[1:n],sigmaW[1:n,1:n])
  for(i in 1:n) {
    y[i] ~ dnorm(U[i] + (beta1 * x1[i] + V[i]) * ((tau-i) + abs(tau-i)) / (2 * (tau-i)) + 
                   (beta2 * x1[i] + W[i]) * ((i-tau) + abs(i-tau)) / (2 * (i-tau)), 
                 sd = 10) 
    
  }
})

# Setting Values
mu_zero <- rep(0, n)
distmat <- as.matrix(dist(1:n))
x1 <- X

# Create Model Object
constants <- list(n = n, x1 = x1, mu_zero = mu_zero)
data <- list(y = y, distmat = distmat)
inits <- list(beta1 = 1, beta2 = 2, tau = 10.1, phiU = 1, phiW = 1, phiV = 1, U = mu_zero, V = mu_zero, W = mu_zero)
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

# Compile Model, run MCMC
cmodel = compileNimble(model)
mcmc.out <- nimbleMCMC(code = code, constants = constants,
                       data = data, inits = inits,
                       nchains = 2, niter = 10000,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c('beta1', "beta2", "tau", "phiU", "phiW", "phiV", "U", "V", "W"))
# Model Summary
mcmc.out$summary



#IRRELEVANT NOW, JUST KEEPING IT JUST IN CASE
#mu[1:n] <- U[1:n] + beta1 * x1[1:n] * ((tau-i) + abs(tau-i)) / (2 * (tau-i)) + 
#beta2*x1[1:n] * ((i-tau) + abs(i-tau)) / (2 * (i-tau))
#covar_mat <- exp(-distmat[1:n,1:n]/phiW) * ((tau-i) + abs(tau-i)) / (2 * (tau-i)) + 
#exp(-distmat[1:n,1:n]/phiV) * ((i-tau) + abs(i-tau)) / (2 * (i-tau))
#y[1:n] ~ dmnorm(mu[1:n], cov = covar_mat) 
