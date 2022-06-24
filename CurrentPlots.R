# Change Point Exploration

library(tidyverse)

dat <- read.csv("currentjoineddata.csv")
counties <- unique(dat$Recip_County)
n.counties <- length(counties)
populations <- unique(dat$population)


df <- data.frame("County" = counties, "Population" = populations,
                 "ChangePoint_Ord" = tausave[,1],
                 "ChangePoint_Cont" = tausave[,4])





# Relationship between population and change point
plot(df$Population, df$ChangePoint_Ord,
     xlab = "Population",
     ylab = "Change Point (Ordinal Response)",
     main = "Population vs. Change Point")




# Difference in change points (ordinal vs continuous)



# Change points on Google map of New York







