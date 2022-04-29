library(nimble, warn.conflicts = FALSE)
library(MASS)


#Simulate the Data
X <- seq(1,100,1)
n = 100
y <- mvrnorm(1, X*4, diag(n)*10)

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


#Nimble Code
code <- nimbleCode({
  beta1 ~ dnorm(0, sd = 200)
  beta2 ~ dnorm(0, sd = 200) 
  sigma1 ~ dgamma(0.1,0.1)
  sigma2 ~ dgamma(0.1,0.1)
  tau ~ dunif(0,100)
  for(i in 1:n) {
    y[i] ~ dnorm(beta1 * x1[i] * ((tau-i) + abs(tau-i)) / (2 * (tau-i)) + 
                   beta2*x1[i] * ((i-tau) + abs(i-tau)) / (2 * (i-tau)), 
                 sd = sigma1 * ((tau-i) + abs(tau-i)) / (2 * (tau-i)) + 
                   sigma2 * ((i-tau) + abs(i-tau)) / (2 * (i-tau))) 
    
  }
})


#Create Model
x1 <- X
constants <- list(n = n, x1 = x1, t = t)
data <- list(y = y)
inits <- list(beta1 = 1, beta2 = 2, tau = 10.1, sigma1 = 1, sigma2 = 1)
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

#Make sure it compiles and run MCMC
cmodel = compileNimble(model)
mcmc.out <- nimbleMCMC(code = code, constants = constants,
                       data = data, inits = inits,
                       nchains = 2, niter = 10000,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c('beta1', "sigma1",'sigma2', "beta2", "tau"))
mcmc.out$summary
