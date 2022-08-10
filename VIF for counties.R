install.packages("tidyverse")
install.packages("readxl")
install.packages("car")
library(car)
library(tidyverse)
library(readxl)
data <- read.csv("updatedjoineddata.csv")

view(data)
data$prev_log_death

data1<-data%>%
  select(Administered_Dose1_Recip.x,Administered_Dose1_Pop_Pct.x,population,death,time, 
         prev_week_total,new_cases,new_cases_per_100k.x,new_vax_this_week,prev_log_death)


colnames(data)
colnames(data1)

data1.lm<-lm(prev_log_death~., data=data1)
summary(data1.lm)
vif(data1.lm)

#everything with a vif over 5/10 is something we should look at. So here, we are seeing which variables have a high variance inflation factor, 
#or rather, which variables have similar effects on prev_log_death to other variables. Seeing that they are almomst all over 5 and many are over 10, 
#we need to take a closer look at how we're going to analyze the data. If it were linear regression, we would have to find the best model that 
#would cut out a lot of values.

#~~~~~~~~~~~~~~~~~~~~~~~~
#do VIF but for different counties, with a loop. 
#Refer to current model.R
#do a unique of county names, use that vector...

#look at Adam's Window R code

