library(tidyverse)
library(gridExtra)

dat <- read.csv("currentjoineddata.csv")
counties <- unique(dat$Recip_County)
n.counties <- length(counties)
populations <- unique(dat$population)
thiscounty <- counties[11]
n.weeks <- dim(this <- dat %>%
    subset(Recip_County==thiscounty))[1]

int <- 20

logdeath <- rep(0, n.weeks-int)
vax <- rep(0, n.weeks-int)
death <- rep(0, n.weeks - int)
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
  
  logdeath[i] <- mean(window$prev_log_new_death)
  vax[i] <- median(window$total_first_dose)
  death[i] <- median(window$death)
}


df <- data.frame("Prob" = c(cat1, cat2, cat3, cat4), 
                 "Cat" = c(rep(1, n.weeks-int), rep(2, n.weeks-int), rep(3, n.weeks-int), rep(4, n.weeks-int)),
                 "Vax" = c(vax, vax, vax, vax),
                 "Death" = c(death, death, death, death),
                 "Logdeath" = c(logdeath, logdeath, logdeath, logdeath))

title_vax <- paste(thiscounty, "- Vax vs Ordinal Probability")

title_death <- paste(thiscounty, "- Death vs Ordinal Probability")

title_logdeath <- paste(thiscounty, "- Log Prev Death vs Ordinal Probability")

county <- gsub(" ", "", thiscounty)

point_vax <- ggplot(data = df, mapping = aes(x = Vax, y = Prob, 
                                group = factor(Cat), color = factor(Cat))) + 
  geom_point() + 
  ggtitle(title_vax) + 
  labs(color = "Category", group = "Category")

line_vax <- ggplot(data = df, mapping = aes(x = Vax, y = Prob, 
                                group = factor(Cat), color = factor(Cat))) + 
  geom_line() + 
  ggtitle(title_vax) + 
  labs(color = "Category", group = "Category")

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

grid.arrange(point_vax, line_vax, ncol = 2)

grid.arrange(point_death, line_death, ncol = 2)

grid.arrange(point_logdeath, line_logdeath, ncol = 2)












categories <- 1:4
pred_category <- rep(0, n.weeks-int)
for (i in 1:(n.weeks-int)) {
  pred_category[i] <- mean(sample(categories, size = 1000, replace = TRUE, prob = c(cat1[i], cat2[i], cat3[i], cat4[i])))
}

par(mfrow = c(2, 2))

plot(1:(n.weeks-int), pred_category,
     xlab = "Interval",
     ylab = "Average Predicted Category",
     main = "Average Predicted Category - Sampled Data")

plot(vax, pred_category,
     xlab = "Vax",
     ylab = "Average Predicted Category",
     main = "Average Predicted Category - Sampled Data")

plot(death, pred_category,
     xlab = "Death",
     ylab = "Average Predicted Category",
     main = "Average Predicted Category - Sampled Data")

plot(logdeath, pred_category,
     xlab = "Prev_new_log_death",
     ylab = "Average Predicted Category",
     main = "Average Predicted Category - Sampled Data")



single_category <- rep(0, (n.weeks-int))
for (i in 1:(n.weeks-int)) {
  single_category[i] <- mean(sample(categories, size = 1, replace = TRUE, prob = c(cat1[i], cat2[i], cat3[i], cat4[i])))
}
plot(1:(n.weeks-int), single_category)



