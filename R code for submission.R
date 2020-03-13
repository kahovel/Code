### ANALYSES FOR "JOINT EFFECTS OF PATCH EDGES AND HABITAT DEGRADATION ON FAUNAL PREDATION RISK 
### IN A WIDESPREAD MARINE FOUNDATION SPECIES" 
# Script author: Kevin Hovel 
# Contact: khovel@sdsu.edu
# Last updated: 13 March 2020

#clear R
rm(list = ls())

# Load libraries
library(effsize)
library(MuMIn)
library(bbmle)
library(dplyr)
library(ggplot2)
library(lme4)
library(nlme)
library(devtools)
library(reshape2)
library(car)
library(ggthemes) #allows custom plot themes
library(emmeans)
library(phia)
library(multcompView)
library(agricolae)
library(broom)
library(tidyverse)
library(Rmisc)
library(ggrepel) #allows option of having data labels in graphs not overlap one another
library(ARPobservation) #allows log response ratio  to be calculated

geco <- read.csv("D:/Research 2015 ZEN/ECO/R analysis/R data files/Final data for ms/predators_mod.csv", header=TRUE) 
alloddsRR <- read.csv("D:/Research 2015 ZEN/ECO/R analysis/R data files/Final data for ms/edge_effect_size_RR.csv", header=TRUE)
reductRR <- read.csv("D:/Research 2015 ZEN/ECO/R analysis/R data files/Final data for ms/reduct_effect_size_RR.csv", header=TRUE)

# Load data on outcome of tethering experiments, with associated environmental variables
geco <- read.csv("XXXX", header=TRUE) 

#Omit missing data
geco <- na.omit(geco)
alloddsRR <- na.omit(alloddsRR)
reductRR <- na.omit(reductRR)

#ensure factors are recognized and create log transformed versions of independent variables
geco$Shoot.reduction <- factor(geco$Shoot.reduction) #make sure shoot reduction is a factor
geco$Location <- factor(geco$Location) #make sure location is a factor
geco$Crust.biomass <- as.numeric(as.character(geco$Crust.biomass))
geco$Epiphyte.biomass <- as.numeric(as.character(geco$Epiphyte.biomass))
geco$Total.biomass <- as.numeric(as.character(geco$Total.biomass))
geco$Other.biomass <- as.numeric(as.character(geco$Other.biomass))
geco$Shtdens <- as.numeric(as.character(geco$Shtdens))
geco$Calc.shoot.dens <- as.numeric(as.character(geco$Calc.shoot.dens))

geco$log.Starting.shoot.dens <- log10(geco$Starting.shoot.dens + 1)
geco$log.Crust.biomass <- log10(geco$Crust.biomass + 1e-3)
geco$log.Epiphyte.biomass <-  log10(geco$Epiphyte.biomass + 1e-4)
geco$log.Calc.shoot.dens <- log10(geco$Calc.shoot.dens + 1)
geco$log.Shtdens <- log10(geco$Shtdens + 1)

#Subset the data so that predictions can be made separately for the different levels of habitat degradation
geco.ambient <- subset(geco, Shoot.reduction == "Amb") #use to select ambient plots (0), 50% plots (1), or 80% plots (2)
geco.fifty <- subset(geco, Shoot.reduction == "Fifty") #subsets the data to yield only plots with 50% removal
geco.eighty <- subset(geco, Shoot.reduction == "Eighty") #subsets the data to yield only plots with 80% removal

geco.ambient <- na.omit(geco.ambient)
geco.eighty <- na.omit(geco.eighty)
geco.eighty <- na.omit(geco.eighty)

##########PART 1: ANALYSES AT INDIVIDUAL SITES########

#Run logistic regressions for each site that calculate odds of mortality based on edge proximity, degradation, 
#shoot density, crustacean biomass, etc. Below, compare models using AIC.
m1 <- glm(Status.bi ~  Location*Shoot.reduction + log.Calc.shoot.dens + log.Crust.biomass + log.Epiphyte.biomass, 
          data = subset(geco, Site == "XX"), #replace "XX" with appropriate site code (see supplementary material for a list)
          family = "binomial")
odds <- exp(cbind(OR = coef(m1), confint(m1, level = 0.95)))
print(odds)
summary(m1)

m2 <- glm(Status.bi ~  Location*Shoot.reduction + log.Calc.shoot.dens + log.Crust.biomass, 
          data = subset(geco, Site == "XX"), 
          family = "binomial")
odds <- exp(cbind(OR = coef(m2), confint(m2, level = 0.95)))
print(odds)
summary(m2)

m3 <- glm(Status.bi ~  Location*Shoot.reduction + log.Calc.shoot.dens, 
          data = subset(geco, Site == "XX"), 
          family = "binomial")
odds <- exp(cbind(OR = coef(m3), confint(m3, level = 0.95)))
print(odds)
summary(m3)

m4 <- glm(Status.bi ~  Location*Shoot.reduction + log.Crust.biomass, 
          data = subset(geco, Site == "XX"), 
          family = "binomial")
odds <- exp(cbind(OR = coef(m4), confint(m4, level = 0.95)))
print(odds)
summary(m4)

m5 <- glm(Status.bi ~  Location*Shoot.reduction + Epiphyte.biomass, 
          data = subset(geco, Site == "XX"), 
          family = "binomial")
odds <- exp(cbind(OR = coef(m5), confint(m5, level = 0.95)))
print(odds)
summary(m5)

m6 <- glm(Status.bi ~  Location*Shoot.reduction, 
          data = subset(geco, Site == "XX"), 
          family = "binomial")
odds <- exp(cbind(OR = coef(m6), confint(m6, level = 0.95)))
print(odds)
summary(m6)

m7 <- glm(Status.bi ~  Location + Shoot.reduction, 
          data = subset(geco, Site == "XX"), 
          family = "binomial")
odds <- exp(cbind(OR = coef(m7), confint(m7, level = 0.95)))
print(odds)
summary(m7)

m8 <- glm(Status.bi ~  Location, 
          data = subset(geco, Site == "XX"), 
          family = "binomial")
odds <- exp(cbind(OR = coef(m8), confint(m8, level = 0.95)))
print(odds)
summary(m8)

m9 <- glm(Status.bi ~  Shoot.reduction, 
          data = subset(geco, Site == "XX"), 
          family = "binomial")
odds <- exp(cbind(OR = coef(m9), confint(m9, level = 0.95)))
print(odds)
summary(m9)

m10 <- glm(Status.bi ~  1, 
          data = subset(geco, Site == "XX"), 
          family = "binomial")
summary(m10)

#Rank by AIC
ModComp <- AICctab(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, base = T, weights = T, 
                   nobs = length(geco))
ModComp

#####ANALYSES INCORPORATING ALL SITES##########

#First we need effect sizes for edge vs. interior and for ambient vs. 50% habitat degradation, and 
# ambient vs. 80% habitat degradation. These are stored in the data set "alloddsRR" which is used below.

with(subset(geco.ambient, Site == "XX"), #replace "XX" with appropriate site code as listed in supplementary material
     logRespRatio(observations = Crust.biomass, phase = Location, base_level = "Interior", bias_correct=TRUE))
     #do the same with Epiphyte.biomass and Calc.shoot.dens

with(subset(geco.fifty, Site == "XX"), #replace "XX" with appropriate site code as listed in supplementary material
     logRespRatio(observations = Crust.biomass, phase = Location, base_level = "Interior", bias_correct=TRUE))
     #do the same with Epiphyte.biomass and Calc.shoot.dens

with(subset(geco.eighty, Site == "XX"), #replace "XX" with appropriate site code as listed in supplementary material
     logRespRatio(observations = Crust.biomass, phase = Location, base_level = "Interior", bias_correct=TRUE))
     #do the same with Epiphyte.biomass and Calc.shoot.dens

#Generalized linear mixed model testing how edge effect size for environmental variables affects predation risk
#edge effect sizes are expressed at response ratios for crustacean biomass, shoot density, and epibiont biomass
#predation risk is expressed as the odds ratio generated for each site by code above

M1e <- lme(log.Odds ~ Shoot.reduction*RespRaCrust + RespRaShoots + RespRaEpi, random = ~ 1 | Site, data=alloddsRR)
anova(M1e)
summary(M1e)
AIC(M1e)

M2e <- lme(log.Odds ~ Shoot.reduction + RespRaCrust + RespRaShoots + RespRaEpi, random = ~ 1 | Site, data=alloddsRR)
anova(M2e)
summary(M2e)
AIC(M2e)

M3e <- lme(log.Odds ~ Shoot.reduction*RespRaEpi + RespRaCrust + RespRaShoots, random = ~ 1 | Site, data=alloddsRR)
anova(M3e)
summary(M3e)
AIC(M3e)

M4e <- lme(log.Odds ~ Shoot.reduction*RespRaShoots + RespRaCrust + RespRaEpi, random = ~ 1 | Site, data=alloddsRR)
anova(M4e)
summary(M4e)
AIC(M4e)

M5e <- lme(log.Odds ~ Shoot.reduction*RespRaCrust + RespRaEpi, random = ~ 1 | Site, data=alloddsRR)
anova(M5e)
summary(M5e)
AIC(M5e)

M6e <- lme(log.Odds ~ Shoot.reduction*RespRaCrust + RespRaShoots, random = ~ 1 | Site, data=alloddsRR)
anova(M6e)
summary(M6e)
AIC(M6e)

M7e <- lme(log.Odds ~ Shoot.reduction*RespRaCrust, random = ~ 1 | Site, data=alloddsRR)
anova(M7e)
summary(M7e)
AIC(M7e)

M8e <- lme(log.Odds ~ 1, random = ~ 1 | Site, data=alloddsRR)
anova(M8e)
summary(M8e)
AIC(M8e)

#Rank by AIC
ModComp2 <- AICctab(m1e, m2e, m3e, m4e, m5e, m6e, m7e, m8e, base = T, weights = T, 
                   nobs = length(geco))
ModComp2

#Generalized linear mixed model testing how habitat degradation effect size for environmental variables affects predation risk
#habitat degradation effect sizes are expressed at response ratios for crustacean biomass and epibiont biomass
#ambient shoot density also is used in the model
#predation risk is expressed as the odds ratio generated for each site by code above

#First set of 6 models are for 50% habitat degradation

M1h50 <- lme(OR.fifty ~ RespRaCrustFif + RespRaEpiFif + Mn.shoots, random = ~ 1 | Site, data=reductRR) 
anova(M1h50)
summary(M1h50)
AIC(M1h50)

M2h50 <- lme(OR.fifty ~ RespRaCrustFif + Mn.shoots, random = ~ 1 | Site, data=reductRR) 
anova(M2h50)
summary(M2h50)
AIC(M2h50)

M3h50 <- lme(OR.fifty ~ RespRaEpiFif + Mn.shoots, random = ~ 1 | Site, data=reductRR) 
anova(M3h50)
summary(M3h50)
AIC(M3h50)

M4h50 <- lme(OR.fifty ~ RespRaCrustFif + RespRaEpiFif, random = ~ 1 | Site, data=reductRR) 
anova(M4h50)
summary(M4h50)
AIC(M4h50)

M5h50 <- lme(OR.fifty ~ Mn.shoots, random = ~ 1 | Site, data=reductRR) 
anova(M5h50)
summary(M5h50)
AIC(M5h50)

M6h50 <- lme(OR.fifty ~ 1, random = ~ 1 | Site, data=reductRR) 
anova(M6h50)
summary(M6h50)
AIC(M6h50)

#Next set of models are for 80% habitat degradation

M1h80 <- lme(OR.eighty ~ RespRaCrustEig + RespRaEpiEig + Mn.shoots, random = ~ 1 | Site, data=reductRR) 
anova(M1h80)
summary(M1h80)
AIC(M1h80)

M2h80 <- lme(OR.eighty ~ RespRaCrustEig + Mn.shoots, random = ~ 1 | Site, data=reductRR) 
anova(M2h80)
summary(M2h80)
AIC(M2h80)

M3h80 <- lme(OR.eighty ~ RespRaEpiEig + Mn.shoots, random = ~ 1 | Site, data=reductRR) 
anova(M3h80)
summary(M3h80)
AIC(M3h80)

M4h80 <- lme(OR.eighty ~ RespRaCrustEig + RespRaEpiEig, random = ~ 1 | Site, data=reductRR) 
anova(M4h80)
summary(M4h80)
AIC(M4h80)

M5h80 <- lme(OR.eighty ~ Mn.shoots, random = ~ 1 | Site, data=reductRR) 
anova(M5h80)
summary(M5h80)
AIC(M5h80)

M6h80 <- lme(OR.eighty ~ 1, random = ~ 1 | Site, data=reductRR) 
anova(M6h80)
summary(M6h80)
AIC(M6h80)


####PLOTS FOR PREDATION RISK EDGE VS. INTERIOR EFFECT SIZES####

#Make scatterplots for predation risk vs. edge effect sizes

#First create new variables that allow points to be color coded by basin or shoot reduction
alloddsRR$Basin1 <- factor(alloddsRR$Basin, levels = c("A","P"))
# create a character vector of colornames
colr <- as.character(unique(alloddsRR$Basin1))

#Shoot reduction color coding.
alloddsRR$SR1 <- factor(alloddsRR$Shoot.reduction, levels = c("0% loss","50% loss", "80% loss"))
# create a character vector of colornames
colr <- as.character(unique(alloddsRR$SR1))

#Make the scatterplot, with options
#Use first line of this graphing code for separate plots by shoot reduction; use the 2nd line for all data together
oddsplot <- ggplot(alloddsRR, aes(x = RespRaShoots, y = log.Odds, label=Site, legend.plot=FALSE)) +  
  geom_point(colour="black", size=12) +
  #Use the next two lines to color by basin; comment this out to color by shoot reduction
  #  geom_point(aes(shape="circle", color=Basin), size=9) +
  #  scale_color_manual(values = c("#0000CC", "#009900")) +
  #Use the next two lines to color by shoot reduction, and comment out the two lines above
  geom_point(aes(shape="circle", color=Shoot.reduction), size=9) +
  scale_color_manual(values = c("#519cb2", "#f5f503", "#286F44")) +
  geom_smooth(method = 'lm', col="black", linetype="solid", size=2, se=TRUE) +
  theme(panel.grid = element_blank()) +
  #geom_smooth(span=1, colour="cornflowerblue") +
  ylab("Edge effect strength (log odds ratio)\n") +
  xlab("\nEdge effect size: xxxxxxxxxx") + #put in correct axis label
  theme_bw() +
  theme(text=element_text(size = 30)) +
  theme(legend.position="none") +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text.x = element_blank()) +
  # legend.position = c(0.18, 0.79) +
  theme(axis.text.x = element_text(color="Black", size=30),
        axis.text.y = element_text(color="Black", 
                                   size=30)) +
  theme( axis.line = element_line(colour = "Black", 
                                  size = 1.0, linetype = "solid")) +
  theme( axis.line.x.top = element_line(colour = "Black", size = 1.4, linetype = "solid")) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(axis.ticks.x = element_line(size = 1.5)) +
  theme(axis.ticks.y = element_line(size = 1.5))
oddsplot

#To make sure that the correct number of decimal places are placed on axis tick labels
oddsplotdecimal <- oddsplot + scale_x_continuous(labels = scales::number_format(accuracy = 0.1))
oddsplotdecimal

#to label individual data points
oddsplotlabels <- oddsplotdecimal + geom_text_repel(show.legend = FALSE, hjust=0.5, vjust=2, size=9)
oddsplotlabels

#Use to adjust plot dimensions; change numbers in parentheses
oddsplotlabelsexpand <- oddsplotlabels + expand_limits(x = c(-2.2, 1.4))
oddsplotlabelsexpand

#To view separately by shoot reduction
oddsplotsep <- oddsplot + facet_grid(~Shoot.reduction3 ~ ., scales="free")
oddsplotsep

#To get both labels and separate plots
oddsplotlabsep <- oddsplotlabels + facet_grid(~Shoot.reduction3 ~ .)  #, scales="free"
oddsplotlabsep


#Make a boxplot showing odds ratios for each level of habitat degradation
alloddsRR$Shoot.reduction3 <- factor(alloddsRR$Shoot.red.per., levels = c("0% loss", "50% loss", "80% loss"))
oddsplot2 <- ggplot(data = alloddsRR, aes(x = Shoot.red.per., y = log.Odds,fill=Shoot.reduction3)) + 
  geom_boxplot(lwd=1.5, outlier.size = 2)+ geom_point(size=3, colour="Black", position=position_jitterdodge(jitter.width=0)) + #pch=21
  scale_fill_manual(values = c("#519cb2", "#f5f503", "#286F44")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("\nHabitat degradation treatment") +
  ylab("Edge effect strength (log odds ratio)\n") +
  theme_bw() +theme(axis.line = element_line(colour = "black", size = 1.4, linetype = "solid")) +
  theme(text=element_text(size = 30)) +
  theme(legend.position="none") +
  scale_x_discrete(limits=c("0% loss","50% loss","80% loss")) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text.x = element_blank()) +
  # legend.position = c(0.18, 0.79) +
  theme(axis.text.x = element_text(color="Black", size=30),
        axis.text.y = element_text(color="Black", size=30)) +
  theme(axis.ticks.length=unit(0.4,"cm")) +
  theme(axis.ticks.x = element_line(size = 1.5)) +
  theme(axis.ticks.y = element_line(size = 1.5))
oddsplot2


