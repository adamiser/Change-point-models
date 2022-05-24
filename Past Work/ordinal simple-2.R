library(MASS)
library(nimble, warn.conflicts = FALSE)

n <- 100 #number of times

#Can loop through all of this if want to simulate values for multiple locations
x <- runif(n)-0.5 #covariate (centered at 0)
distmat <- as.matrix(dist(1:n)) #distance in time
x <- sort(x)

beta1 <- 2 #slope
beta2 <- 8

y <- mvrnorm(1, beta1*x[1:50], diag(50))
y2 <- mvrnorm(1, beta2*x[51:100], diag(50))
y <- c(y,y2)

z <- seq(1,100,1)
for(i in 1:100){
  if(y[i]< -1){
    z[i] <- 1
  }
  else if(y[i] < -0){
    z[i] <- 2
  }
  else if(y[i] < 1.1) {
    z[i] <- 3
  }
  else{
    z[i] <- 4
  }
}

#Look to see how it behaves
plot(x, z)
plot(x,y)
#plot(1:n, z)
#par(mfrow=c(2,1))

#--------------------------------------------------------
#Try just ordinal
#Change gammas around
#

code <- nimbleCode({
  beta1 ~ dnorm(0, sd = 200)
  beta2 ~ dnorm(0, sd = 200) 
  tau ~ dunif(0,100)
  sigma <- 1
  alpha3 ~ dnorm(0,1)
  alpha4 ~ dnorm(0,1)
  gam[1] <- -9999999
  gam[2] <- 0
  gam[3] <- exp(alpha3)
  gam[4] <- exp(alpha4) + gam[3]
  gam[5] <- 9999999
  
  for(i in 1:n){
    mu[i] <- beta1 * x1[i] * ((tau-i) + abs(tau-i)) / (2 * (tau-i)) + 
      beta2 * x1[i] * ((i-tau) + abs(i-tau)) / (2 * (i-tau))
    z[i,] ~ dmulti(psi[i,1:4], size=1)
    for(j in 1:4){
      psi[i,j] <- iprobit((gam[j+1]-mu[i])/sigma) - iprobit((gam[j]-mu[i])/sigma)
    }
  } 
})

# Setting Values
mu_zero <- rep(0, n)
distmat <- as.matrix(dist(1:n))
x1 <- x
z2 <- matrix(0, nrow=n, ncol=4)
for(i in 1:n){
  z2[i,z[i]] <- 1
}

constants <- list(n = n, x1 = x1, mu_zero = mu_zero)
data <- list(z = z2)
gaminit <- c(-2, 0, exp(0), exp(1) + exp(0),  5)
y1 <- runif(length(z), min=gaminit[z], max=gaminit[z+1])
betainit <- c(solve(t(x1)%*%x1)%*%t(x1)%*%y1)
inits <- list(beta1 = betainit, beta2 = betainit, tau = 50.1, sigma=2, alpha3 = 0, alpha4 = 1)
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

# Compile Model, run MCMC
cmodel = compileNimble(model)
mcmc.out <- nimbleMCMC(code = code, constants = constants,
                       data = data, inits = inits, nburnin=2000,
                       nchains = 2, niter = 20000,thin=2,
                       summary = TRUE,
                       monitors = c('beta1', "beta2", "alpha3","alpha4", "tau"))
# Model Summary
mcmc.out$summary

save(mcmc.out, file=paste("simulation", kk, ".Rdata", sep=''))