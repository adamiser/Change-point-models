library(tidyverse)

dat <- read.csv("currentdata.csv")

counties <- unique(dat$county)
n.counties <- length(counties)

dat <- dat[order(dat$week), ]
dat <- dat[order(dat$county), ]


dat$Recip_County <- paste(dat$county, "County")
dat$Date <- dat$week
dat$Administered_Dose1_Pop_Pct.x <- dat$first_dose_prevalence


dat$new_vax_this_week <- rep(0, dim(dat)[1])
dat$new_vax_this_week[1] <- 0

thiscounty <- dat$county[1]
temp <- dat[dat$county == thiscounty, ]
size <- dim(temp)[1] - 1
for (i in 1:size) {
  temp$new_vax_this_week[i+1] <- dat$total_first_dose[i+1] - dat$total_first_dose[i]
  data <- temp
}

for (j in 2:n.counties) {
  thiscounty <- dat$county[j]
  temp <- dat[dat$county == thiscounty, ]
  size <- dim(temp)[1] - 1
  for (i in 1:size) {
  temp$new_vax_this_week[i+1] <- temp$total_first_dose[i+1] - temp$total_first_dose[i]
  }
  data <- rbind(data, temp)
}


write.csv(data, "currentjoineddata.csv")
