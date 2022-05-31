library(tidyverse)
library(nimble)
library(coda)
library(parallel)


dat <- read.csv("/Users/adamiser810/Desktop/Change-point-models/updatedjoineddata.csv")

counties <- unique(dat$Recip_County)
n.counties <- length(counties)

tausave <- matrix(0, ncol=6, nrow=n.counties)
beta1deathsave <- matrix(0, ncol=6, nrow = n.counties)
beta2deathsave <- matrix(0, ncol=6, nrow = n.counties)
beta1vaxsave <- matrix(0, ncol=6, nrow = n.counties)
beta2vaxsave <- matrix(0, ncol=6, nrow = n.counties)
phiUsave <- matrix(0, ncol=6, nrow = n.counties)
phiVsave <- matrix(0, ncol=6, nrow = n.counties)
phiWsave <- matrix(0, ncol=6, nrow = n.counties)
sigmasave <- matrix(0, ncol=6, nrow=n.counties)

waic <- matrix(0, ncol = 2, nrow = n.counties)
gelman <- vector('list', length = n.counties)


this_cluster <- makeCluster(4)

for (jj in 1:2) {
  
  thiscounty <- counties[jj]
  
  #Select just one county
  this <- dat %>%
    subset(Recip_County==thiscounty) %>%
    select(Date, category, log_new_death, new_vax_this_week, new_cases_per_100k.y)
  
run_MCMC_allcode <- function(seed, this) {
  library(nimble)
    
  

  ###############
  # Create data for Nimble Model
  # Note: Using rev() bc data appear to be in reverse time order so I'm putting them back in chronological order 
  z <- rev(this$category)
  x1 <- rev(this$log_new_death)
  x2 <- rev(this$new_vax_this_week)
  x2[1] <- 0 #no new vaccinations before the first date
  n <- length(z)
  y <- rev(this$new_cases_per_100k.y)

  ###############
  #Nimble Ordinal Model

  codeOrd <- nimbleCode({
    beta1death ~ dnorm(0, sd = 200)
    beta2death ~ dnorm(0, sd = 200) 
    beta1vax ~ dnorm(0, sd = 200)
    beta2vax ~ dnorm(0, sd = 200) 
    tau ~ dunif(5,n-5)
    sigma ~ dgamma(1, .1)
    alpha3 ~ dnorm(0,1)
    alpha4 ~ dnorm(0,1)
    gam[1] <- -9999999
    gam[2] <- 0
    gam[3] <- exp(alpha3)
    gam[4] <- exp(alpha4) + gam[3]
    gam[5] <- 9999999
  
    phiU ~ dunif(0,20)
    phiW ~ dunif(0,20)
    phiV ~ dunif(0,20)
    cholU[1:n,1:n] <- chol(exp(-distmat[1:n,1:n]/phiU))
    U[1:n] ~ dmnorm(mu_zero[1:n], cholesky = cholU[1:n,1:n], prec_param = 1)
    cholW[1:n,1:n] <- chol(exp(-distmat[1:n,1:n]/phiW))
    cholV[1:n,1:n] <- chol(exp(-distmat[1:n,1:n]/phiV))
    W[1:n] ~ dmnorm(mu_zero[1:n], cholesky=cholW[1:n,1:n], prec_param=1)
    V[1:n] ~ dmnorm(mu_zero[1:n], cholesky=cholV[1:n,1:n], prec_param=1)
  
    for(i in 1:n){
      mu[i] <- U[i] + (beta1death * x1[i] + beta1vax*x2[i] + V[i]) * ((tau-i) + abs(tau-i)) / (2 * (tau-i)) + 
      (beta2death * x1[i] +beta2vax*x2[i] + W[i])* ((i-tau) + abs(i-tau)) / (2 * (i-tau))
      z[i,] ~ dmulti(psi[i,1:4], size=1)
      for(j in 1:4){
        psi[i,j] <- iprobit((gam[j+1]-mu[i])/sigma) - iprobit((gam[j]-mu[i])/sigma)
      }
    } 
  })

  #Setting Values
  mu_zero <- rep(0, n)
  distmat <- as.matrix(dist(1:n))
  z2 <- matrix(0, nrow=n, ncol=4)
  for(i in 1:n){
    z2[i,z[i]] <- 1
  }

  constantsOrd <- list(n = n, x1 = x1,x2=x2, mu_zero = mu_zero)
  dataOrd <- list(z = z2, distmat=distmat)
  gaminit <- c(-2, 0, exp(0), exp(1) + exp(0),  5)
  y1 <- runif(length(z), min=gaminit[z], max=gaminit[z+1])
  xmat <- cbind(x1, x2)
  betainit <- c(solve(t(xmat)%*%xmat)%*%t(xmat)%*%y1)
  initsOrd <- list(beta1death = betainit[1], beta2death = betainit[1], 
              beta1vax=betainit[2], beta2vax=betainit[2], tau = n/2+.1, 
              sigma=1, alpha3 = 0, alpha4 = 1, phiU=1, phiW=1, phiV=1, U=mu_zero,
              W=mu_zero, V=mu_zero)
  modelOrd <- nimbleModel(codeOrd, constants = constantsOrd, data = dataOrd, inits = initsOrd)

  #Compile Model, run MCMC
  cmodelOrd = compileNimble(modelOrd)
  myMCMC = buildMCMC(cmodelOrd)
  ordMCMC = compileNimble(myMCMC)
  
  results <- runMCMC(ordMCMC, niter = 100000, nchains = 4,, nburnin = 5000, thin = 5,
                     summary = TRUE, samplesAsCodaMCMC = TRUE,
                     setSeed = seed)
  
  return(results)
  
}
}
  
  
  chain_output <- parLapply(cl = this_cluster,
                            X = 1:4,
                            fun = run_MCMC_allcode,
                            this = this)
  
  stopCluster(this_cluster)
  
  
  chain_output
  