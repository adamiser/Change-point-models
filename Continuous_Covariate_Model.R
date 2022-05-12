library(MASS)
library(nimble, warn.conflicts = FALSE)
library(tidyverse)

code <- nimbleCode({
  beta1death ~ dnorm(0, sd = 200)
  beta2death ~ dnorm(0, sd = 200) 
  beta1time ~ dnorm(0, sd = 200)
  beta2time ~ dnorm(0, sd = 200) 
  beta1vac ~ dnorm(0, sd = 200)
  beta2vac ~ dnorm(0, sd = 200) 
  
  tau ~ dunif(10, 90)
  sigma ~ dgamma(2, 0.5)
  
  phiU ~ dgamma(2, 2)
  phiW ~ dgamma(2,2)
  phiV ~ dgamma(2,2)
  cholU[1:n,1:n] <- chol(exp(-distmat[1:n,1:n]/phiU))
  cholW[1:n,1:n] <- chol(exp(-distmat[1:n,1:n]/phiW))
  cholV[1:n,1:n] <- chol(exp(-distmat[1:n,1:n]/phiV))
  U[1:n] ~ dmnorm(mu_zero[1:n], cholesky = cholU[1:n,1:n], prec_param = 0)
  W[1:n] ~ dmnorm(mu_zero[1:n], cholesky = cholW[1:n,1:n], prec_param = 0)
  V[1:n] ~ dmnorm(mu_zero[1:n], cholesky = cholV[1:n,1:n], prec_param = 0)
  
  for(i in 1:n){
    mu[i] <- U[i] + 
      (beta1time * t[i] + beta1death * deaths[i] + beta1vac * vac[i] + V[i]) * 
      ((tau-i) + abs(tau-i)) / (2 * (tau-i)) + 
    (beta2time * t[i] + beta2death * deaths[i] + beta2vac * vac[i] + W[i]) * 
      ((i-tau) + abs(i-tau)) / (2 * (i-tau))
    y[i] ~ dnorm(mu[i], sd = sigma)
    }
})

# Setting Values
albany <- albany[order(albany$time), ]
n <- dim(albany)[1]
y <- albany$new_cases_per_100k
t <- albany$time
deaths <- albany$prev_log_new_death
vac <- albany$Administered_Dose1_Pop_Pct
mu_zero <- rep(0, n)
distmat <- as.matrix(dist(1:n))

constants <- list(n = n,  mu_zero = mu_zero, distmat = distmat)
data <- list(deaths = deaths, t = t, vac = vac, y = y)
inits <- list(beta1time = 2, beta2time = 2, 
              beta1death = 2, beta2death = 2,
              beta1vac = 2, beta2vac = 2,
              tau = 10.1, sigma=1,
              phiU = 1, phiW = 1, phiV = 1, 
              U = mu_zero, W=mu_zero, V=mu_zero)
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

# Compile Model, run MCMC
cmodel = compileNimble(model)
mcmc.out <- nimbleMCMC(code = code, constants = constants,
                       data = data, inits = inits, nburnin=2000,
                       nchains = 5, niter = 200000,thin=20,
                       summary = TRUE,
                       monitors = c('beta1death', "beta2death", 
                                    "beta1time", "beta2time", 
                                    "beta1vac", "beta2vac",
                                    "tau", "phiU", "phiV", "phiW"),
                       samplesAsCodaMCMC = TRUE,
                       WAIC = TRUE)
# Model Summary
capture.output(mcmc.out$summary, file = "Model_output")

mcmc.chains <- mcmc(mcmc.out$samples)
mcmc.sims <- mcmc(as.matrix(mcmc.out$samples))
rhat = gelman.diag(mcmc.chains)
capture.output(rhat$psrf, file = "Model_gelman")
capture.output(autocorr(mcmc.sims), file = "Model_autocorr")
capture.output(mcmc.out$WAIC, file = "Model_WAIC")

