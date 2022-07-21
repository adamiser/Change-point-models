library(tidyverse)

dat <- read.csv("currentjoineddata.csv")
counties <- unique(dat$Recip_County)
n.counties <- length(counties)
populations <- unique(dat$population)
thiscounty <- counties[1]
n.weeks <- dim(this <- dat %>%
    subset(Recip_County==thiscounty))[1]

int <- 30

logdeath <- rep(0, n.weeks-int)
vax <- rep(0, n.weeks-int)
cat1 <- rep(0, n.weeks-int)
cat2 <- rep(0, n.weeks-int)
cat3 <- rep(0, n.weeks-int)
cat4 <- rep(0, n.weeks-int)

for (i in 1:(n.weeks-int)) {
  this <- dat %>%
    subset(Recip_County==thiscounty)
  lower <- i
  upper <- i + 30
  window <- this[lower:upper,]
  
  cat1[i] <- sum(window$category == 1) / int
  cat2[i] <- sum(window$category == 2)  / int
  cat3[i] <- sum(window$category == 3)  / int
  cat4[i] <- sum(window$category == 4)  / int
  
  logdeath[i] <- median(window$prev_log_new_death)
  vax[i] <- median(window$total_first_dose)
}


# Logdeath vs prob category
par(mfrow = c(2, 2))
plot(logdeath, cat1)
plot(logdeath, cat2)
plot(logdeath, cat3)
plot(logdeath, cat4)


# Vax vs prob category
plot(vax, cat1)
plot(vax, cat2)
plot(vax, cat3)
plot(vax, cat4)


# Time vs prob category
plot(1:92, cat1)
plot(1:92, cat2)
plot(1:92, cat3)
plot(1:92, cat4)

