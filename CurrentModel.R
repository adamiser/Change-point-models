library(tidyverse)
library(nimble)
library(coda)

#Load the data
dat <- read.csv("currentjoineddata.csv")

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
gamma3save <- matrix(0, ncol=3, nrow=n.counties)
gamma4save <- matrix(0, ncol=3, nrow=n.counties)

waic <- matrix(0, ncol = 2, nrow = n.counties)
gelman <- vector('list', length = n.counties)
raftery <- vector('list', length = n.counties)

for (jj in 1:n.counties) {
  thiscounty <- counties[jj]
  
  #Select just one county
  this <- dat %>%
    subset(Recip_County==thiscounty) %>%
    select(Date, category, prev_log_new_death, Administered_Dose1_Pop_Pct.x, log_new_death, 
           new_vax_this_week, new_cases_per_100k)

  ###############
  # Create data for Nimble Model
  # Note: Using rev() bc data appear to be in reverse time order so I'm putting them back in chronological order 
  z <- rev(this$category)
  x1 <- rev(this$prev_log_new_death)
  x2 <- rev(this$Administered_Dose1_Pop_Pct.x)
  x2[1] <- 0 #no new vaccinations before the first date
  n <- length(z)
  y <- rev(this$new_cases_per_100k)

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
    
    beta1deathstar <- beta1death / sigma
    beta2deathstar <- beta2death / sigma
    beta1vaxstar <- beta1vax / sigma
    beta2vaxstar <- beta2vax / sigma
  
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
  mcmc.outOrd <- nimbleMCMC(code = codeOrd, constants = constantsOrd,
                       data = dataOrd, inits = initsOrd, nburnin=5000,
                       nchains = 2, niter = 100000,thin=5,
                       summary = TRUE,
                       monitors = c('beta1death', "beta2death", "beta1vax", "beta2vax", 
                                    "alpha3","alpha4", "tau", "phiU", "phiV", "phiW", "sigma",
                                    "beta1deathstar", "beta2deathstar",
                                    "beta1vaxstar", "beta2vaxstar",
                                    "gam"),
                       samplesAsCodaMCMC = TRUE,
                       WAIC = TRUE)
  #Model Summary
  #mcmc.outOrd$summary

  tausave[jj,1:3] <- mcmc.outOrd$summary$all.chains[20,c(1,4,5)]
  beta1deathsave[jj, 1:3] <- mcmc.outOrd$summary$all.chains[4, c(1, 4, 5)]
  beta2deathsave[jj, 1:3] <- mcmc.outOrd$summary$all.chains[8, c(1, 4, 5)]
  beta1vaxsave[jj, 1:3] <- mcmc.outOrd$summary$all.chains[6, c(1, 4, 5)]
  beta2vaxsave[jj, 1:3] <- mcmc.outOrd$summary$all.chains[10, c(1, 4, 5)]
  phiUsave[jj, 1:3] <- mcmc.outOrd$summary$all.chains[16, c(1, 4, 5)]
  phiVsave[jj, 1:3] <- mcmc.outOrd$summary$all.chains[17, c(1, 4, 5)]
  phiWsave[jj, 1:3] <- mcmc.outOrd$summary$all.chains[18, c(1, 4, 5)]
  sigmasave[jj, 1:3] <- mcmc.outOrd$summary$all.chains[19, c(1, 4, 5)]
  gam3save[jj, 1:3] <- mcmc.outOrd$summary$all.chains[13, c(1, 2, 3)]
  gam4save[jj, 1:3] <- mcmc.outOrd$summary$all.chains[14, c(1, 2, 3)]
  
  waic[jj, 1] <- mcmc.outOrd$WAIC$WAIC
  # gelman[[jj]][[1]] <- gelman.diag(mcmc.outOrd$samples)
  raftery[[jj]][[1]] <- raftery.diag(mcmc.outOrd$samples)

  ###############################
  # Compare the change point when we estimate with continuous response

  codeCont <- nimbleCode({
    beta1death ~ dnorm(0, sd = 200)
    beta2death ~ dnorm(0, sd = 200) 
    beta1vax ~ dnorm(0, sd = 200)
    beta2vax ~ dnorm(0, sd = 200) 
    tau ~ dunif(5,n-5)
    sigma ~ dgamma(1, .1)
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
      y[i] ~ dnorm(mu[i], sigma)
    } 
  })

  # Setting Values
  mu_zero <- rep(0, n)
  distmat <- as.matrix(dist(1:n))

  constantsCont <- list(n = n, x1 = x1,x2=x2, mu_zero = mu_zero)
  dataCont <- list(y = y, distmat=distmat)
  xmat <- cbind(x1, x2)
  betainit <- c(solve(t(xmat)%*%xmat)%*%t(xmat)%*%y)
  initsCont <- list(beta1death = betainit[1], beta2death = betainit[1], 
                 beta1vax=betainit[2], beta2vax=betainit[2], tau = n/2+.1, 
                 sigma=1, phiU=1, phiW=1, phiV=1, U=mu_zero,
                 W=mu_zero, V=mu_zero)
  modelCont <- nimbleModel(codeCont, constants = constantsCont, data = dataCont, inits = initsCont)

  # Compile Model, run MCMC
  cmodelCont = compileNimble(modelCont)
  mcmc.outCont <- nimbleMCMC(code = codeCont, constants = constantsCont,
                          data = dataCont, inits = initsCont, nburnin=5000,
                          nchains = 2, niter = 100000,thin=5,
                          summary = TRUE,
                          monitors = c('beta1death', "beta2death", "beta1vax", "beta2vax", 
                                      "tau", "phiU", "phiV", "phiW", "sigma"),
                          samplesAsCodaMCMC = TRUE,
                          WAIC = TRUE)
  # Model Summary
  #mcmc.outCont$summary
  
  tausave[jj,4:6] <- mcmc.outCont$summary$all.chains[9,c(1,4,5)]
  beta1deathsave[jj, 4:6] <- mcmc.outCont$summary$all.chains[1, c(1, 4, 5)]
  beta2deathsave[jj, 4:6] <- mcmc.outCont$summary$all.chains[3, c(1, 4, 5)]
  beta1vaxsave[jj, 4:6] <- mcmc.outCont$summary$all.chains[2, c(1, 4, 5)]
  beta2vaxsave[jj, 4:6] <- mcmc.outCont$summary$all.chains[4, c(1, 4, 5)]
  phiUsave[jj, 4:6] <- mcmc.outCont$summary$all.chains[5, c(1, 4, 5)]
  phiVsave[jj, 4:6] <- mcmc.outCont$summary$all.chains[6, c(1, 4, 5)]
  phiWsave[jj, 4:6] <- mcmc.outCont$summary$all.chains[7, c(1, 4, 5)]
  sigmasave[jj, 4:6] <- mcmc.outCont$summary$all.chains[8, c(1, 4, 5)]
  
  waic[jj, 2] <- mcmc.outCont$WAIC$WAIC
  gelman[[jj]][[2]] <- gelman.diag(mcmc.outCont$samples)
  raftery[[jj]][[2]] <- raftery.diag(mcmc.outCont$samples)
  
  print(paste("Finished:", jj))
}


save(tausave, beta1deathsave, beta2deathsave, beta1vaxsave, beta2vaxsave,
     phiUsave, phiVsave, phiWsave, sigmasave, gam3save, gam4save, waic, gelman, raftery,
     counties, file="CurrentModel.Rdata")



        