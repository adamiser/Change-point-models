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
                 "fips"=fipscodes$fips)




# Relationship between population and change point
plot(df$ChangePoint_Ord, log(df$Population), 
     ylab = "log Population",
     xlab = "Change Point (Ordinal Response)",
     main = "Population vs. Change Point")




# Difference in change points (ordinal vs continuous)



# Change points on Google map of New York

library(ggplot2)
library(usmap)

# Map of change points
plot_usmap(data=df, values="ChangePoint_Ord", include="New York") + scale_fill_continuous(name="Change Point") + 
  theme(legend.position='right') + 
  ggtitle("Mean Change Points Mapped Across NY")

# Map of population
df$logpop <- log(population)
plot_usmap(data=df, values="logpop", include="New York") + scale_fill_continuous(name="log pop") + 
  theme(legend.position='right') + 
  ggtitle("Log Population Mapped Across NY")




# Plots of weekly ordinal data with change points marked? 
ii <- 1 #for(ii in 1:n.counties){
thiscounty <- subset(dat, Recip_County==counties[ii])
plot(thiscounty$category, type='l', xlab="Week", ylab="Category", ylim=c(1, 4), main=counties[ii])
abline(v=tausave[ii,1], col='darkgray')
abline(v=tausave[ii,2:3], lty=2, col='darkgray') #I didn't look to see what these were but assumed they were bounds on the CI?
#}



