## Analysis of heat tolerances (CTmax) of bromeliad insect species
## Sarah Amundrud


# load libraries
library(car)


## import data for CTmax (CTmax_individuals.csv)
dataCT <-read.csv(file.choose(),header=T)
head(dataCT)
tail(dataCT)
str(dataCT)
levels(dataCT$type)
levels(dataCT$species)



#### does average CTmax depend on mountain side?

# Tipulid
model <- lm(dataCT$CTmax[which(dataCT$type == "Tipulid")] ~ dataCT$side[which(dataCT$type == "Tipulid")])
summary(model)
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
Anova(model)
plot(dataCT$CTmax[which(dataCT$type == "Tipulid")] ~ dataCT$side[which(dataCT$type == "Tipulid")])

# Scirtid
model <- lm(dataCT$CTmax[which(dataCT$type == "Scirtid")] ~ dataCT$side[which(dataCT$type == "Scirtid")])
summary(model)
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
Anova(model)
plot(dataCT$CTmax[which(dataCT$type == "Scirtid")] ~ dataCT$side[which(dataCT$type == "Scirtid")])

#Poly
model <- lm(dataCT$CTmax[which(dataCT$type == "Polypedilum")] ~ dataCT$side[which(dataCT$type == "Polypedilum")])
summary(model)
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
Anova(model)
plot(dataCT$CTmax[which(dataCT$type == "Polypedilum")] ~ dataCT$side[which(dataCT$type == "Polypedilum")])

# Culex
model <- lm(dataCT$CTmax[which(dataCT$type == "Culex")] ~ dataCT$side[which(dataCT$type == "Culex")])
summary(model)
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
Anova(model)
plot(dataCT$CTmax[which(dataCT$type == "Culex")] ~ dataCT$side[which(dataCT$type == "Culex")])


## mountain side does not affect CTmax of any focal species
## note that Damselflies were only present on Atlantic slope, so did not do analyisis for Damselflies



## calculate mean CTmax for each species to get species-specific heat tolerances

Meancomb<-tapply(dataCT$CTmax, dataCT$type, mean); Meancomb
SEcomb <- tapply(dataCT$CTmax, dataCT$type, sd)/
  sqrt(tapply(dataCT$CTmax, dataCT$type, length)); SEcomb

## add means and SEs to data frame
CTmaxMeans <- data.frame(species = c("Culex", "Damsel", "Poly", "Scirtid", "Tipulid"))
CTmaxMeans$CTmax_mean <- as.numeric(Meancomb)
CTmaxMeans$CTmax_SE <- as.numeric(SEcomb)
CTmaxMeans

## extract the species specific heat tolerances to be used in further analysis
#write.csv(CTmaxMeans, file = "CTmax_species.csv")

