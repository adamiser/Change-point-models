library(MASS)
library(nimble, warn.conflicts = FALSE)
library(tidyverse)
library(lubridate)
library(coda)

new_data <- read.csv("joineddata.csv")
initial_data <- read.csv("US_COVID_weekly_NY_data.csv")
new_data$county <- gsub( " .*$", "", new_data$Recip_County)

simple_data <- new_data[c("Date", "time", "population", "Administered_Dose1_Pop_Pct",
                          "prev_log_new_death", "new_cases_per_100k", "category", "county")]

# Combining with data that is before vaccines
initial_data$week <- as.Date(initial_data$week)
initial_data$year <- year(initial_data$week)
initial_data$date <- initial_data$week
before_vaccine <- initial_data[initial_data$date < "2020-12-20",]
before_vaccine$Administered_Dose1_Pop_Pct = 0
before_vaccine$county <- gsub( ",.*$", "", before_vaccine$key)
colnames(before_vaccine)[26] = "Date"
before_vaccine <- before_vaccine[c("Date", "time", "population", "Administered_Dose1_Pop_Pct",
                          "prev_log_new_death", "new_cases_per_100k", "category", "county")]

vaccine_county_df <- rbind(simple_data, before_vaccine)

albany <- vaccine_county_df[vaccine_county_df$county =='Albany',]

code <- nimbleCode({
  betadeath ~ dnorm(0, sd = 200)
  betavac ~ dnorm(0, sd = 200)
  beta1time ~ dnorm(0, sd = 200)
  beta2time ~ dnorm(0, sd = 200) 
  
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
    mu[i] <- U[i] + betavac * vac[i] + betadeath * deaths[i] + 
      (beta1time * t[i] + V[i]) * 
      ((tau-i) + abs(tau-i)) / (2 * (tau-i)) + 
    (beta2time * t[i] + W[i]) * 
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
              betadeath = 2,betavac = 2,
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
                       monitors = c("betadeath", 
                                    "beta1time", "beta2time", 
                                    "betavac", 
                                    "tau", "phiU", "phiV", "phiW"),
                       samplesAsCodaMCMC = TRUE,
                       WAIC = TRUE)
# Model Summary
capture.output(mcmc.out$summary, file = "Model_output2")

mcmc.chains <- mcmc(mcmc.out$samples)
mcmc.sims <- mcmc(as.matrix(mcmc.out$samples))
rhat = gelman.diag(mcmc.chains)
capture.output(rhat$psrf, file = "Model_gelman2")
capture.output(autocorr(mcmc.sims), file = "Model_autocorr2")
capture.output(mcmc.out$WAIC, file = "Model_WAIC2")
capture.output(raftery.diag(mcmc.sims), file = "Model_raftery2")

