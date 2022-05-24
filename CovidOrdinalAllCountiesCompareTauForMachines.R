library(tidyverse)
library(nimble)

#Load the data
dat <- read.csv("updatedjoineddata.csv")

counties <- unique(dat$Recip_County)
n.counties <- length(counties)

tausave <- matrix(0, ncol=6, nrow=n.counties)

for(jj in 1:n.counties){
  thiscounty <- counties[jj]
  
  #Select just one county
  this <- dat %>%
    subset(Recip_County==thiscounty) %>%
    select(Date, category, log_new_death, new_vax_this_week, new_cases_per_100k.y)

  ###############
  # Create data for Nimble Model
  # Note: Using rev() bc data appear to be in reverse time order so I'm putting them back in chronological order 
  z <- rev(this$category)
  x1 <- rev(this$log_new_death)
  x2 <- rev(this$new_vax_this_week)
  x2[1] <- 0 #no new vaccinations before the first date
  n <- length(z)
  y <- rev(log(this$new_cases_per_100k.y))

  ###############
  #Nimble Ordinal Model

  codeOrd <- nimbleCode({
    beta1death ~ dnorm(0, sd = 200)
    beta2death ~ dnorm(0, sd = 200) 
    beta1vax ~ dnorm(0, sd = 200)
    beta2vax ~ dnorm(0, sd = 200) 
    tau ~ dunif(5,n-5)
    sigma <- 1
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
  mcmc.outOrd <- nimbleMCMC(code = codeOrd, constants = constantsOrd,
                       data = dataOrd, inits = initsOrd, nburnin=2000,
                       nchains = 2, niter = 20000,thin=2,
                       summary = TRUE,
                       monitors = c('beta1death', "beta2death", "beta1vax", "beta2vax", 
                                    "alpha3","alpha4", "tau", "phiU", "phiV", "phiW"))
  #Model Summary
  #mcmc.outOrd$summary

  tausave[jj,1:3] <- mcmc.outOrd$summary$all.chains[10,c(1,4,5)]

  ###############################
  # Compare the change point when we estimate with continuous response

  codeCont <- nimbleCode({
    beta1death ~ dnorm(0, sd = 200)
    beta2death ~ dnorm(0, sd = 200) 
    beta1vax ~ dnorm(0, sd = 200)
    beta2vax ~ dnorm(0, sd = 200) 
    tau ~ dunif(5,n-5)
    sigma <- 1
  
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
                          data = dataCont, inits = initsCont, nburnin=2000,
                          nchains = 2, niter = 20000,thin=2,
                          summary = TRUE,
                          monitors = c('beta1death', "beta2death", "beta1vax", "beta2vax", 
                                      "tau", "phiU", "phiV", "phiW"))
  # Model Summary
  #mcmc.outCont$summary
  
  tausave[jj,4:6] <- mcmc.outCont$summary$all.chains[8,c(1,4,5)]
  print(paste("Finished:", jj))
}

save(tausave, counties, file="CovidOrdinalAllCountiesCompareTauForMachines.Rdata")


################
# I added this after I ran the above code (i.e., this wasn't part of the code when I ran it on the machine)
library(tidyverse)
setwd() #set directory where the files are saved
load('CovidOrdinalAllCountiesCompareTauForMachines.Rdata')
dat <- read.csv("updatedjoineddata.csv")

#Add latitude and longitude to the saved change point data set
cy <- dat %>% 
  group_by(Recip_County) %>%
  summarise_all(mean) %>%
  select(Recip_County, lat, long)
unique(counties)
tausave <- data.frame(tausave, counties)
names(tausave) <- c("O.mean", "O.low", "O.high", "C.mean", "C.low", "C.high", "Recip_County")
fuldat <- merge(tausave, cy)


#Quick plots to look at change points across NY
ggplot(aes(x=long, y=lat, colour=O.mean, cex=O.mean), data=fuldat) + geom_point() + labs(title="Change Point Mean - Ordinal Response")
ggplot(aes(x=long+.1, y=lat+.1, colour=C.mean, cex=C.mean), data=fuldat) + geom_point() + labs(title="Change Point Mean - Log New Cases Response")

#Wide plot to look at uncertainty in the change points 
nc <- length(counties)
for(i in 1:nc){
  if(i == 1){
    plot(c(i,i), tausave[i,2:3], col=rgb(.4, .2, .6, .5), lwd=4, xlim=c(1,nc), ylim=range(tausave), xaxt='n', ylab="Change Point", xlab="county", type='l')
    lines(c(i,i), tausave[i,5:6], col=rgb(.2, .3, .6, .7), lwd=4)
  }else{
    lines(c(i,i), tausave[i,2:3], col=rgb(.4, .2, .6, .5), lwd=4)
    lines(c(i,i), tausave[i,5:6], col=rgb(.2, .3, .6, .7), lwd=4)
  }
  text(i, 25, labels=counties[i], offset=0, col=rgb(.1, .1, .1, .5), srt=90)
}
abline(v=seq(.5, nc+.5), lty=2, col=rgb(.1, .1, .1, .1))

        