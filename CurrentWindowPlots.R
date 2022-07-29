library(tidyverse)
library(gridExtra)

dat <- read.csv("currentjoineddata.csv")
counties <- unique(dat$Recip_County)
n.counties <- length(counties)
populations <- unique(dat$population)
thiscounty <- counties[31]
n.weeks <- dim(this <- dat %>%
    subset(Recip_County==thiscounty))[1]

int <- 30

logdeath <- rep(0, n.weeks-int)
vax <- rep(0, n.weeks-int)
death <- rep(0, n.weeks - int)
cat1 <- rep(0, n.weeks-int)
cat2 <- rep(0, n.weeks-int)
cat3 <- rep(0, n.weeks-int)
cat4 <- rep(0, n.weeks-int)

data <- list()
logdeath_list <- list()
cat1_list <- list()
cat2_list <- list()
cat3_list <- list()
cat4_list <- list()
iter <- 1

for (kk in c(2, 15, 31)) {
  thiscounty <- counties[kk]
  
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
    
    logdeath[i] <- mean(window$prev_log_new_death)
    vax[i] <- median(window$Administered_Dose1_Pop_Pct.x)
    death[i] <- median(window$death)
    
    
  }
  
  df <- data.frame("Prob" = c(cat1, cat2, cat3, cat4), 
                   "Cat" = c(rep(1, n.weeks-int), rep(2, n.weeks-int), rep(3, n.weeks-int), rep(4, n.weeks-int)),
                   "Vax" = c(vax, vax, vax, vax),
                   "Death" = c(death, death, death, death),
                   "Logdeath" = c(logdeath, logdeath, logdeath, logdeath))
  
  logdeath_list[[iter]] <- logdeath
  cat1_list[[iter]] <- cat1
  cat2_list[[iter]] <- cat2
  cat3_list[[iter]] <- cat3
  cat4_list[[iter]] <- cat4
  
  data[[iter]] <- df
  iter <- iter + 1
  
}


df <- data.frame("Prob" = c(cat1, cat2, cat3, cat4), 
                 "Cat" = c(rep(1, n.weeks-int), rep(2, n.weeks-int), rep(3, n.weeks-int), rep(4, n.weeks-int)),
                 "Vax" = c(vax, vax, vax, vax),
                 "Death" = c(death, death, death, death),
                 "Logdeath" = c(logdeath, logdeath, logdeath, logdeath))

title_vax <- paste(thiscounty, "- Percent Vaccinated vs Ordinal Probability")

county <- gsub(" ", "", thiscounty)

point_vax <- ggplot(data = df, mapping = aes(x = Vax, y = Prob, 
                                group = factor(Cat), color = factor(Cat))) + 
  geom_point() + 
  ggtitle(title_vax) + 
  labs(color = "Category", group = "Category")

line_vax <- ggplot(data = df, mapping = aes(x = Vax, y = Prob, 
                                group = factor(Cat), color = factor(Cat))) + 
  geom_line(size = 2) + 
  ggtitle(title_vax) + 
  labs(color = "Category", group = "Category") + 
  xlab("First Dose Vaccination - Percent of Population") + 
  ylab("Probability")



line_vax

# grid.arrange(point_vax, line_vax, ncol = 2)

# grid.arrange(point_death, line_death, ncol = 2)

# grid.arrange(point_logdeath, line_logdeath, ncol = 2)

### Vax GGplot - 3 plots (New York City, Buffalo, Rural)
### Line plot - make lines thicker, make it look nicer
# Manhattan - County 31
# Buffalo - County 15
# Allegany - County 2












categories <- 1:4
pred_category <- list(rep(0, (n.weeks-int)),
                      rep(0, (n.weeks-int)),
                      rep(0, (n.weeks-int)))
for (kk in 1:3) {
  for (i in 1:(n.weeks-int)) {
    pred_category[[kk]][i] <- mean(sample(categories, size = 1000, 
                                       replace = TRUE, 
                                       prob = c(cat1_list[[kk]][i], cat2_list[[kk]][i], 
                                                cat3_list[[kk]][i], cat4_list[[kk]][i])))
  }
}

### Add 3 counties to the plot using different symbols
### Good for the "prev_new_log_death" variable

df <- data.frame("County" = c(rep("Allegany", 92), rep("Erie", 92), rep("NY", 92)),
                 "Logdeath" = c(logdeath_list[[1]], logdeath_list[[2]],
                                 logdeath_list[[3]]),
                 "Pred_Category" = c(pred_category[[1]], pred_category[[2]],
                                     pred_category[[3]]))

ggplot(data = df, mapping = aes(x = Logdeath, y = Pred_Category,
                                group = factor(County), color = factor(County))) +
  geom_point() + 
  labs(color = "County", group = "County") + 
  xlab("Previous Week Log(Death)") + 
  ylab("Expected Category in Interval") + 
  ggtitle("Previous Week Log(Death) vs Expected Category for 3 NY Counties")

















### Old Code ###




title_death <- paste(thiscounty, "- Death vs Ordinal Probability")

title_logdeath <- paste(thiscounty, "- Log Prev Death vs Ordinal Probability")

point_death <- ggplot(data = df, mapping = aes(x = Death, y = Prob, 
                                group = factor(Cat), color = factor(Cat))) + 
  geom_point() + 
  ggtitle(title_death) + 
  labs(color = "Category", group = "Category")

line_death <- ggplot(data = df, mapping = aes(x = Death, y = Prob, 
                                group = factor(Cat), color = factor(Cat))) + 
  geom_line() + 
  ggtitle(title_death) + 
  labs(color = "Category", group = "Category")

point_logdeath <- ggplot(data = df, mapping = aes(x = Logdeath, y = Prob, 
                                group = factor(Cat), color = factor(Cat))) + 
  geom_point() + 
  ggtitle(title_logdeath) + 
  labs(color = "Category", group = "Category")
  
line_logdeath <- ggplot(data = df, mapping = aes(x = Logdeath, y = Prob, 
                                group = factor(Cat), color = factor(Cat))) + 
  geom_line() + 
  ggtitle(title_logdeath) + 
  labs(color = "Category", group = "Category")

