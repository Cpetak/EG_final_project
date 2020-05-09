library(tidyr)
library(ggplot2)
library(tidyverse)
library(car)

setwd("~/Desktop/EG_final_project/EG_final_project/phenotypic_data")


#analysing data for edge vs core in common garden
#height data
#-------------------------------------------------


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


#ANOVA

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

detach(ave)


#-----------------------------------
#budding
bud <- read.csv(file="BudSet2019.csv", header=TRUE)
data_bud <- bud[,c("plant_ID", "Family", "Region","BudSet" )]
data_bud <-  data_bud[-grep('M', data_bud$Region), ]
data_bud <- na.omit(data_bud)

ave_bud <- data_bud %>% 
  group_by(Family, Region) %>% 
  summarise(average = mean(BudSet))

ave_bud$Region <- as.character(ave_bud$Region)
colnames(ave_bud) <- c("Family","Region", "aBud")

#ANOVA

model = lm(aBud ~Region,data=ave_bud)
anova(model) 

hist(residuals(model),
     col="darkgray")

attach(ave_bud)

bartlett.test(aBud ~ interaction(Region)) # for heterogenous varience
shapiro.test(resid(aov(aBud ~ Region))) #for normality
# ----> not normally distributed residuals, ANOVA won't do

kruskal.test(aBud ~ Region, data = ave_bud) 


boxplot(aBud ~ Region,
        data = ave_bud,
        xlab = "Genotype x Sex",
        ylab = "MPI Activity")

detach(ave_bud)

#-----------------------------------------------
#analysing data for CW vs HD
#Harvest
harvest <- read.csv(file="harvest_data.csv", header=TRUE)

get_stats <- function(cols, treatment){
  result = "messed_up"
  harvest_green <- harvest[,c("Family", "Group","Treatment", cols )]
  harvest_green_10 <-  harvest_green[grep(treatment, harvest_green$Treatment), ]
  colnames(harvest_green_10) <- c("Family","Group", "Treatment", "Temp")
  ave_harv <- harvest_green_10 %>% 
    group_by(Family, Group) %>% 
    summarise(average = mean(Temp))
  
  attach(ave_harv)
  
  b=bartlett.test(average ~ interaction(Group))[3] # for heterogenous varience
  s = shapiro.test(resid(aov(average ~ Group)))[2] # for normal distribution
  
  if (b > 0.05 & s > 0.05){
    #print("PASSED_parametric")
    p = summary(aov(average~Group, data = ave_harv))[[1]][["Pr(>F)"]][1]
    if (p < 0.05){
      #print("Significant")
      p = as.character(p)
      result <- c("significant", p)
    } else {
      p = as.character(p)
      result <- c("not", p)
    }
    #print(p)
  } else {
      k = kruskal.test(average ~ Group, data = ave_harv) 
      #print(k[3])
      if (k[3] < 0.05){
        p = as.character(k[3])
        result <- c("significant, k", p)
      } else {
        p = as.character(k[3])
        result <- c("not", p)
      }
  }
  
  boxplot(average ~ Group,
          data = ave_harv,
          xlab = "Source",
          ylab = cols)
  
  detach(ave_harv)
  
  return(ave_harv)
}


mlist = c("control_0", "control_5", "control_10", "drought/heat_0", "drought/heat_5", "drought/heat_10", "heat_0", "heat_5", "heat_10")
allday_list = c("control","drought/heat", "heat")


for (el in allday_list){
  a <- get_stats("stem_wt", el)
  print(el)
  print(a)
}

#PCA data
pdata <- harvest[,c("Family", "Group","Treatment", "green_wt", "stem_wt", "needle_wt", "total_dry_wt", "root2shoot_ratio", "root_wt", "live_crown_ht")]

tpdata <- pdata %>% separate(Treatment, c("Treatment", "Day"), "_")

pdata.pr <- prcomp(tpdata[c(5:11)], center = TRUE, scale = TRUE)
summary(pdata.pr)
plot(pdata.pr$x[,1],pdata.pr$x[,2], xlab="PC1 (77.4%)", ylab = "PC2 (16.2%)", main = "PC1 / PC2 - plot")

library("factoextra")
p <- fviz_pca_ind(pdata.pr, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = tpdata$Day, 
             #col.ind = "black", 
             addEllipses = TRUE,
             label = "var",
             col.ind = tpdata$Day,
             mean.point = FALSE,
             ellipse.level = 0.95,
             repel = TRUE,
             legend.title = "Day") +
  theme(plot.title = element_text(hjust = 0.5))

p
p + scale_color_manual(values=c("blue", "red"))

#some more boxplots
a=get_stats("needle_wt", "drought/heat_0")
b=get_stats("needle_wt", "drought/heat_5")
c=get_stats("needle_wt", "drought/heat_10")
d=get_stats("needle_wt", "heat_0")
e=get_stats("needle_wt", "heat_5")
f=get_stats("needle_wt", "heat_10")
g=get_stats("needle_wt", "control_0")
h=get_stats("needle_wt", "control_5")
i=get_stats("needle_wt", "control_10")

attach(mtcars)
par(mfrow=c(3,3))
boxplot(average~Group, a, xlab="Drought/heat 0")
boxplot(average~Group, b, xlab="Drought/heat 5")
boxplot(average~Group, c, xlab="Drought/heat 10")
boxplot(average~Group, d, xlab="Heat 0")
boxplot(average~Group, e, xlab="Heat 5")
boxplot(average~Group, f, xlab="Heat 10")
boxplot(average~Group, g, xlab="Control 0")
boxplot(average~Group, h, xlab="Control 5")
boxplot(average~Group, i, xlab="Control 10")
