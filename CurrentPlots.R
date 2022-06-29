# Change Point Exploration

library(tidyverse)

dat <- read.csv("currentjoineddata.csv")
counties <- unique(dat$Recip_County)
n.counties <- length(counties)
populations <- unique(dat$population)

load('CurrentModel.Rdata')

fipscodes <- read.table(header=F, sep=" ", text="
    36001        Albany County
    36003        Allegany County
    36005        Bronx County
    36007        Broome County
    36009        Cattaraugus County
    36011        Cayuga County
    36013        Chautauqua County
    36015        Chemung County
    36017        Chenango County
    36019        Clinton County
    36021        Columbia County
    36023        Cortland County
    36025        Delaware County
    36027        Dutchess County
    36029        Erie County
    36031        Essex County
    36033        Franklin County
    36035        Fulton County
    36037        Genesee County
    36039        Greene County
    36041        Hamilton County
    36043        Herkimer County
    36045        Jefferson County
    36047        Kings County
    36049        Lewis County
    36051        Livingston County
    36053        Madison County
    36055        Monroe County
    36057        Montgomery County
    36059        Nassau County
    36061        NewYork County
    36063        Niagara County
    36065        Oneida County
    36067        Onondaga County
    36069        Ontario County
    36071        Orange County
    36073        Orleans County
    36075        Oswego County
    36077        Otsego County
    36079        Putnam County
    36081        Queens County
    36083        Rensselaer County
    36085        Richmond County
    36087        Rockland County
    36089        StLawrence County
    36091        Saratoga County
    36093        Schenectady County
    36095        Schoharie County
    36097        Schuyler County
    36099        Seneca County
    36101        Steuben County
    36103        Suffolk County
    36105        Sullivan County
    36107        Tioga County
    36109        Tompkins County
    36111        Ulster County
    36113        Warren County
    36115        Washington County
    36117        Wayne County
    36119        Westchester County
    36121        Wyoming County
    36123        Yates County
")
fipscodes <- fipscodes[,c(5,13)]
names(fipscodes)=c("fips", "county")
fipscodes[45:50,] <- fipscodes[c(46:50,45),]
#double check to make sure ordering matches: 
#cbind(fipscodes, counties)

df <- data.frame("County"=counties, "Population"=populations, 
                 "ChangePoint_Ord"=tausave[,1], 
                 "ChangePoint_Cont"=tausave[,4],
                 "ChangePointRange_Ord"=tausave[,3]-tausave[,2],
                 "ChangePointRange_Cont"=tausave[,6]-tausave[,5],
                 "Beta1death_Ord"=beta1deathsave[,1],
                 "Beta2death_Ord"=beta2deathsave[,1],
                 "Beta1vax_Ord"=beta1vaxsave[,1],
                 "Beta2vax_Ord"=beta2vaxsave[,1],
                 "Beta1death_Cont"=beta1deathsave[,4],
                 "Beta2death_Cont"=beta2deathsave[,4],
                 "Beta1vax_Cont"=beta1vaxsave[,4],
                 "Beta2vax_Cont"=beta2vaxsave[,4],
                 "Beta1deathRange_Ord"=beta1deathsave[,3]-beta1deathsave[,2],
                 "Beta1deathRange_Cont"=beta1deathsave[,6]-beta1deathsave[,5],
                 "Beta1vaxRange_Ord"=beta1vaxsave[,3]-beta1vaxsave[,5],
                 "Beta1vaxRange_Cont"=beta1vaxsave[,6]-beta1vaxsave[,5],
                 "Beta2deathRange_Ord"=beta2deathsave[,3]-beta2deathsave[,2],
                 "Beta2deathRange_Cont"=beta2deathsave[,6]-beta2deathsave[,5],
                 "Beta2vaxRange_Ord"=beta2vaxsave[,3]-beta2vaxsave[,5],
                 "Beta2vaxRange_Cont"=beta2vaxsave[,6]-beta2vaxsave[,5],
                 "fips"=fipscodes$fips)




# Relationship between population and change point
plot(df$ChangePoint_Ord, log(df$Population), 
     ylab = "log Population",
     xlab = "Change Point (Ordinal Response)",
     main = "Population vs. Change Point")




# Change points on Google map of New York

library(ggplot2)
library(usmap)
library(gridExtra)


### To consider the question: where are change points across the state?  
### How are the estimated change points different for the ordinal response model and the continuous response model?

# Map of change points
plot_usmap(data=df, values="ChangePoint_Ord", include="New York") + scale_fill_continuous(name="Change Point") + 
  theme(legend.position='right') + 
  ggtitle("Mean Change Points Mapped Across NY")




#map of the range of the 95% CI bounds
range_ord_nyplot <- plot_usmap(data=df, values="ChangePointRange_Ord", include="New York") + scale_fill_continuous(name="Change Point") + 
  theme(legend.position='right') + 
  ggtitle("95% CI Range of Change Point - Ordinal Response")

#Maps of change points and 95% CI bounds for continuous response model
range_cont_nyplot <- plot_usmap(data=df, values="ChangePointRange_Cont", include="New York") + scale_fill_continuous(name="Change Point") + 
  theme(legend.position='right') + 
  ggtitle("95% CI Range of Change Point - Continuous Response")

grid.arrange(range_ord_nyplot, range_cont_nyplot, ncol = 2)

### OBSERVATIONS ###
# It seems that the ordinal response actually predicted the change point with less 
# variability than the continuous model. This may just show that the ordinal model is 
# better than the continuous model? 





### To consider the question: How are beta's behaving before and after change point?

# 2 maps: posterior mean of the beta's for death before the change and beta's for death after change
beta1death_ord_nyplot <- plot_usmap(data=df, values="Beta1death_Ord", include="New York") + 
  scale_fill_continuous(name="Change Point",
                        limits = c(0, 100)) + 
  theme(legend.position='right') + 
  ggtitle("Beta1death - Ordinal Response")

beta2death_ord_nyplot <- plot_usmap(data=df, values="Beta2death_Ord", include="New York") + 
  scale_fill_continuous(name="Change Point",
                        limits = c(0, 100)) + 
  theme(legend.position='right') + 
  ggtitle("Beta2death - Ordinal Response")

grid.arrange(beta1death_ord_nyplot, beta2death_ord_nyplot, ncol = 2)

### OBSERVATIONS ###
# The range of coefficient values across the counties is quite large. Most coefficients are small,
# but others are huge. This makes the plot of NY look rather boring.
# I've adjusted the limits to try and make it more interesting (I don't know that it's
# the best approach). It seems that the coefficients tend to be larger after the change point.



# 2 maps: beta's for vaccination instead of death
beta1vax_ord_nyplot <- plot_usmap(data=df, values="Beta1vax_Ord", include="New York") + 
  scale_fill_continuous(name="Change Point",
                        limits = c(0, 100)) + 
  theme(legend.position='right') + 
  ggtitle("Beta1vax - Ordinal Response")

beta2vax_ord_nyplot <- plot_usmap(data=df, values="Beta2vax_Ord", include="New York") + 
  scale_fill_continuous(name="Change Point",
                        limits = c(0, 1000)) + 
  theme(legend.position='right') + 
  ggtitle("Beta2vax - Ordinal Response")

grid.arrange(beta1vax_ord_nyplot, beta2vax_ord_nyplot, ncol = 2)

### OBSERVATIONS ###
# These look better than the death coefficient. It's clear that the coefficient dramatically 
# increases after the change point - this makes sense because vaccinations likely cause
# the change point, or at least the change point occurs around when vaccinations
# began to reach the public.


# conisder the ranges of credible intervals
remove_outlier <- df[-c(2,33,48,52,60),]
beta1deathRange_ord_nyplot <- plot_usmap(data=remove_outlier, values="Beta1deathRange_Ord", include="New York") + scale_fill_continuous(name="Beta1death") + 
  theme(legend.position='right') + 
  ggtitle("95% CI Range of Beta1death - Ordinal Response")
beta2deathRange_ord_nyplot <- plot_usmap(data=remove_outlier, values="Beta2deathRange_Ord", include="New York") + scale_fill_continuous(name="Beta2death") + 
  theme(legend.position='right') + 
  ggtitle("95% CI Range of Beta2death - Ordinal Response")

grid.arrange(beta1death_nyplot, beta2death_nyplot, ncol = 2)

### OBSERVATIONS ###
# I removed a few outliers that seemed to be disrupting the plots.





# Map the ranges for the beta's for the continuous response model (I expect the range to be a lot smaller)





# Map of population
df$logpop <- log(populations)
plot_usmap(data=df, values="logpop", include="New York") + scale_fill_continuous(name="log pop") + 
  theme(legend.position='right') + 
  ggtitle("Log Population Mapped Across NY")








### IGNORE FOR NOW
# Plots of weekly ordinal data with change points marked? 
ii <- 1 #for(ii in 1:n.counties){
thiscounty <- subset(dat, Recip_County==counties[ii])
plot(thiscounty$category, type='l', xlab="Week", ylab="Category", ylim=c(1, 4), main=counties[ii])
abline(v=tausave[ii,1], col='darkgray')
abline(v=tausave[ii,2:3], lty=2, col='darkgray') #I didn't look to see what these were but assumed they were bounds on the CI?
#}










### TEST ###


# 2 maps: posterior mean of the beta's for death before the change and beta's for death after change
df$logbeta1death_ord <- (df$Beta1death_Ord)^2
df$logbeta2death_ord <- log(df$Beta2death_Ord)^2
beta1death_ord_nyplot <- plot_usmap(data=df, values="logbeta1death_Ord", include="New York") + 
  scale_fill_continuous(name="Change Point") + 
  theme(legend.position='right') + 
  ggtitle("Beta1death - Ordinal Response")

beta2death_ord_nyplot <- plot_usmap(data=df, values="logbeta2death_Ord", include="New York") + 
  scale_fill_continuous(name="Change Point") + 
  theme(legend.position='right') + 
  ggtitle("Beta2death - Ordinal Response")

grid.arrange(beta1death_ord_nyplot, beta2death_ord_nyplot, ncol = 2)

