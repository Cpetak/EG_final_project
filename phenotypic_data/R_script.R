library(tidyr)
library(ggplot2)
library(tidyverse)
library(car)

setwd("~/Desktop/EG_final_project/EG_final_project/phenotypic_data")


#analysing data for edge vs core in common garden
#height data

initial_ht <- read.csv(file="initialHt.csv", header=TRUE)
final_ht <- read.csv(file="final_ht_2019.csv", header=TRUE)

#let's compare all individuals from edge with all individuals from core
#preparing data

data_init <- initial_ht[,c("plant_ID", "Family", "Region","initial_ht" )]
data_init <-  data_init[-grep('M', data_init$Region), ]

data_fin <- final_ht[,c("plant_ID", "Family", "Region","Height_Fall2019" )]
data_fin <-  data_fin[-grep('M', data_fin$Region), ]
data_fin <- data_fin[,c("plant_ID","Height_Fall2019" )]

total <- merge(data_init,data_fin,by="plant_ID")
colnames(total) <- c("plant_ID", "Family", "Region", "1", "2")
total$Change <- total$`2`- total$`1`
ave <- total %>% 
  group_by(Family, Region) %>% 
    summarise(average = mean(Change))

#total <- pivot_longer(total, c("1", "2"), names_to = "Time")
#total$Time <- as.numeric(total$Time)


ave$Region <- as.character(ave$Region)
colnames(ave) <- c("Family","Region", "Change")


#two-way ANOVA

model = lm(Change ~Region,data=ave)
anova(model) 

hist(residuals(model),
     col="darkgray")

attach(ave)

bartlett.test(Change ~ interaction(Region)) # for heterogenous varience
shapiro.test(resid(aov(Change ~ Region))) #for normality

boxplot(Change ~ Region,
        data = ave,
        xlab = "Genotype x Sex",
        ylab = "MPI Activity")


