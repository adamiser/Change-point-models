### Understanding output from Ordinal Bayesian Model

### Step 1: Load the Rdata file "CovidOrdinalAllParameters.Rdata"
# This is found on the Github
# The easiest way that I have found to load the Rdata into your R studio is
# to open R studio, place the Rdata file on your Desktop and double click it;
# this will load the data into your Global Environment

### Step 2: Understanding the data
# Once you load the Rdata file, you'll have a bunch of variables in your Environment.
# Here are the main variables for you to worry about:
# counties - this is a vector of the 62 NY counties in the same order as the analysis
# tausave - these are the change points. These are the most important things to focus on for now


### Step 3: Explore the data
# Once you have the Rdata loaded, you can extract the valuable stuff, like this:

# ord_change_points is all of the posterior mean change points from our ordinal models for each county
ord_change_points <- tausave[,1]
### This is the most important information

# ord_change_point_intervals are the 95% credible intervals around each change point. These are cool
# and can likely be looked at later
ord_change_point_intervals <- tausave[,2:3]


# You can then explore relationships between change point and factors such as population
# Here's a very rough example:
library(tidyverse)
dat <- read.csv("updatedjoineddata.csv")
counties <- unique(dat$Recip_County)

populations <- rep(0, 62)
for (i in 1:length(counties)) {
  populations[i] <- dat[dat$Recip_County == counties[i],]$population[1]
}

plot(ord_change_points, populations)

