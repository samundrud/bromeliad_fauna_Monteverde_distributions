### Analyses for Amundrud and Srivastava - Global Ecology and Biogeorgaphy
## do species interactions explain species elevational distributions?
## Sarah Amundrud


library(lavaan)
library(semPlot)
library(dplyr)
library(tidyr)
library(car)




## get survey data of insect abundances inside bromeliads
survey_melt<-read.csv(file.choose(),header=T) ## bromeliad_insect_abd.csv
names(survey_melt)
head(survey_melt)
tail(survey_melt)
str(survey_melt)




## rearrange data to get one column for each species

survey_melt <- survey_melt
names(survey_melt)

survey <- survey_melt %>%
  spread(taxa,abundance)
head(survey)



#### What is the role of species interactions in affecting species distributions? #####

## SEM for bromeliad food web (focal species) along elevation gradient in MV
## Model indirect effects (mediation analysis)
# response (Y) = detritivore (Culex_erethyzonfer / Scirtidae / Polypedilum)
# predictor (X) = site_elevation
# mediator (M) = Tipulidae, Damselfly (Atlantic only)

names(survey)
# calculate relative abundances of focal species
survey$insects <- survey$Culex_erethyzonfer + survey$Polypedilum + survey$Scirtidae + survey$Tipulidae +
  survey$M_modesta + survey$Dytisctid
survey$Culex_erethyzonfer_rel <- survey$Culex_erethyzonfer / survey$insects
survey$Polypedilum_rel <- survey$Polypedilum / survey$insects
survey$Scirtidae_rel <- survey$Scirtidae / survey$insects
survey$Tipulidae_rel <- survey$Tipulidae / survey$insects
survey$M_modesta_rel <- survey$M_modesta / survey$insects
survey$Dytisctid_rel <- survey$Dytisctid / survey$insects
survey$TopPred_rel <- (survey$M_modesta + survey$Dytisctid) / survey$insects







###### SEM on Pacific side #######
# subset Pacific side only
Pacific <- survey[which(survey$side == "Pacific"),]
Pacific$side <- factor(Pacific$side)
Pacific$site_ID <- factor(Pacific$site_ID)

# How many bromeliads do we have at each side (sample size)
table(Pacific$site_ID)



#################################------------------------------------
## it is important that variances between the different model components are roughly equal (within 10 to 100 fold)
## to check, I make up some random model that includes all components

names(Pacific)
for (i in c(5,17:20)){
  print(var(Pacific[i]))
}


## variances are very unequal!!
################------------------------------------------






### transform to make variances equal
## it is important that variances between the different model components are roughly equal (within 10 to 100 fold or so)
## recode so variances match (need to play around a bit until good transformations are found)
Pacific$site_elevation <- Pacific$site_elevation/150
Pacific$Culex_erethyzonfer_rel <- Pacific$Culex_erethyzonfer_rel*20
Pacific$Scirtidae_rel <- Pacific$Scirtidae_rel*5
Pacific$Polypedilum_rel <- Pacific$Polypedilum_rel*5
Pacific$Tipulidae_rel <- Pacific$Tipulidae_rel*10

## check variances again
for (i in c(5,17:20)){
  print(var(Pacific[i]))
}
## now the variances are roughly equal :)




###### SEMs for the three focal prey species

## Culex
mod_Culex <- ' # direct effect
Culex_erethyzonfer_rel ~ c * site_elevation
# mediator
Tipulidae_rel ~ a * site_elevation
Culex_erethyzonfer_rel ~ b* Tipulidae_rel
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'

fit <- sem(mod_Culex, data = Pacific)
varTable(fit)
summary(fit)
semPaths(fit)



## Scirtidae
mod_Scirtid <- ' # direct effect
Scirtidae_rel ~ c * site_elevation
# mediator
Tipulidae_rel ~ a * site_elevation
Scirtidae_rel ~ b* Tipulidae_rel
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'

fit <- sem(mod_Scirtid, data = Pacific)
varTable(fit)
summary(fit)
semPaths(fit)




## Polypedilum

mod_Poly <- ' # direct effect
Polypedilum_rel ~ c * site_elevation
# mediator
Tipulidae_rel ~ a * site_elevation
Polypedilum_rel ~ b* Tipulidae_rel
# indirect effect (a*b)
ab := a*b
# total effect
total := c + (a*b)
'

fit <- sem(mod_Poly, data = Pacific)
varTable(fit)
summary(fit)
semPaths(fit)

#######





###### SEM on Atlantic side #######
# subset Atlantic side only
Atlantic <- survey[which(survey$side == "Atlantic"),]
Atlantic$side <- factor(Atlantic$side)
Atlantic$site_ID <- factor(Atlantic$site_ID)

# How many bromeliads do we have at each side (sample size)
table(Atlantic$site_ID)

# check variances
names(Atlantic)
for (i in c(5,17:23)){
  print(var(Atlantic[i]))
}
# very unequal!

### transform to make variances equal
## recode so variances match (need to play around a bit until good transformations are found)
Atlantic2 <- Atlantic
Atlantic2$site_elevation <- Atlantic2$site_elevation/150
Atlantic2$Culex_erethyzonfer_rel <- Atlantic2$Culex_erethyzonfer_rel*15
Atlantic2$Scirtidae_rel <- Atlantic2$Scirtidae_rel*5
Atlantic2$Polypedilum_rel <- Atlantic2$Polypedilum_rel*5
Atlantic2$Tipulidae_rel <- Atlantic2$Tipulidae_rel*15
Atlantic2$M_modesta_rel <- Atlantic2$M_modesta_rel*100
Atlantic2$Dytisctid_rel <- Atlantic2$Dytisctid_rel*80
Atlantic2$TopPred_rel <- Atlantic2$TopPred_rel*80


# check variances again
for (i in c(5,17:23)){
  print(var(Atlantic2[i]))
}
# now they look good!




###### SEMs for the three focal prey species - use Damselfly as predator (no Dytisctid) ######

## Culex
mod_Culex <- ' # direct effect
Culex_erethyzonfer_rel ~ f * site_elevation
# mediators
M_modesta_rel ~ a * site_elevation
Culex_erethyzonfer_rel ~ b * M_modesta_rel
Tipulidae_rel ~ c * M_modesta_rel
Culex_erethyzonfer_rel ~ d * Tipulidae_rel
Tipulidae_rel ~ e * site_elevation
# indirect effects of elevation (a*b)
ab := a*b # of Damsel
ed := e*d # of Tipulid
acd := a*c*d # of Damsels and Tipulids

# total effects
total := f + (a*b) + (e*d) + (a*c*d)
'

fit <- sem(mod_Culex, data = Atlantic2)
varTable(fit)
summary(fit)
semPaths(fit)



## Scirtidae
mod_Scirtid <- ' # direct effect
Scirtidae_rel ~ f * site_elevation
# mediators
M_modesta_rel ~ a * site_elevation
Scirtidae_rel ~ b * M_modesta_rel
Tipulidae_rel ~ c * M_modesta_rel
Scirtidae_rel ~ d * Tipulidae_rel
Tipulidae_rel ~ e * site_elevation
# indirect effects of elevation (a*b)
ab := a*b # of Damsel
ed := e*d # of Tipulid
acd := a*c*d # of Damsels and Tipulids

# total effects
total := f + (a*b) + (e*d) + (a*c*d)
'

fit <- sem(mod_Scirtid, data = Atlantic2)
varTable(fit)
summary(fit)
semPaths(fit)




## Polypedilum

mod_Poly <- ' # direct effect
Polypedilum_rel ~ f * site_elevation
# mediators
M_modesta_rel ~ a * site_elevation
Polypedilum_rel ~ b * M_modesta_rel
Tipulidae_rel ~ c * M_modesta_rel
Polypedilum_rel ~ d * Tipulidae_rel
Tipulidae_rel ~ e * site_elevation
# indirect effects of elevation (a*b)
ab := a*b # of Damsel
ed := e*d # of Tipulid
acd := a*c*d # of Damsels and Tipulids

# total effects
total := f + (a*b) + (e*d) + (a*c*d)
'

fit <- sem(mod_Poly, data = Atlantic2)
varTable(fit)
summary(fit)
semPaths(fit)


#######











###### SEMs for the three focal prey species - use obligate predator (M. modesta and Dytisctids combined) #######

## Culex
mod_Culex <- ' # direct effect
Culex_erethyzonfer_rel ~ f * site_elevation
# mediators
TopPred_rel ~ a * site_elevation
Culex_erethyzonfer_rel ~ b * TopPred_rel
Tipulidae_rel ~ c * TopPred_rel
Culex_erethyzonfer_rel ~ d * Tipulidae_rel
Tipulidae_rel ~ e * site_elevation
# indirect effects of elevation (a*b)
ab := a*b # of Top predator
ed := e*d # of Tipulid
acd := a*c*d # of Damsels and Tipulids

# total effects
total := f + (a*b) + (e*d) + (a*c*d)
'

fit <- sem(mod_Culex, data = Atlantic2)
varTable(fit)
summary(fit)
semPaths(fit)



## Scirtidae
mod_Scirtid <- ' # direct effect
Scirtidae_rel ~ f * site_elevation
# mediators
TopPred_rel ~ a * site_elevation
Scirtidae_rel ~ b * TopPred_rel
Tipulidae_rel ~ c * TopPred_rel
Scirtidae_rel ~ d * Tipulidae_rel
Tipulidae_rel ~ e * site_elevation
# indirect effects of elevation (a*b)
ab := a*b # of top predator
ed := e*d # of Tipulid
acd := a*c*d # of top predator and Tipulids

# total effects
total := f + (a*b) + (e*d) + (a*c*d)
'

fit <- sem(mod_Scirtid, data = Atlantic2)
varTable(fit)
summary(fit)
semPaths(fit)




## Polypedilum

mod_Poly <- ' # direct effect
Polypedilum_rel ~ f * site_elevation
# mediators
TopPred_rel ~ a * site_elevation
Polypedilum_rel ~ b * TopPred_rel
Tipulidae_rel ~ c * TopPred_rel
Polypedilum_rel ~ d * Tipulidae_rel
Tipulidae_rel ~ e * site_elevation
# indirect effects of elevation (a*b)
ab := a*b # of top predator
ed := e*d # of Tipulid
acd := a*c*d # of top predator and Tipulids

# total effects
total := f + (a*b) + (e*d) + (a*c*d)
'

fit <- sem(mod_Poly, data = Atlantic2)
varTable(fit)
summary(fit)
semPaths(fit)

#######

