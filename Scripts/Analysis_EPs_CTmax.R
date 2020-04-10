### Analyses for Amundrud and Srivastava - Global Ecology and Biogeorgaphy
## 1. how does insect abundance change along elevation gradients on Pacific and Atlantic slopes
## 2. does heat tolerance (CTmax) predict elevation preference (EP)
## Sarah Amundrud

# load packages
library(plyr)
library(dispmod)
library(dplyr)
library(deming)
library(car)


#### function for error bars for bar graph
y.error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
} 

x.error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(y) != length(x) | length(x) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x+upper,y, x-lower, y, angle=90, code=3, length=length, ...)
} 


## get survey data of insect abundances inside bromeliads
survey_melt<-read.csv(file.choose(),header=T) ## bromeliad_insect_abdp.csv
names(survey_melt)
head(survey_melt)
tail(survey_melt)
str(survey_melt)

# get species-specific CTmax data
CT<-read.csv(file.choose(),header=T) ## CTmax_species.csv
CT









## 1.1. How does insect abundance change along elevation gradients on Pacific and Atlantic slopes? ######
## Poisson models on insect abundances (abd ~ elevation * side + bromeliad maxVol)


### each focal species seperately
par(mfrow=c(3,2))

# pick focal species
species <- survey_melt[(survey_melt$taxa == "Culex_erethyzonfer"),]
species <- survey_melt[(survey_melt$taxa == "Polypedilum"),]
species <- survey_melt[(survey_melt$taxa == "Scirtidae"),]
species <- survey_melt[(survey_melt$taxa == "Tipulidae"),]


mod_abd <- glm(species$abundance ~ species$site_elevation * species$side + species$Max_vol, poisson)
summary(mod_abd) # overdispersion
disp<-glm.poisson.disp(mod_abd)
mod_abd<-update(disp,.~.-z:all)
par(mfrow=c(2,2))
plot(mod_abd)
summary(mod_abd)
Anova(mod_abd)




## for damselflies, don't include mountain side (only present on Atlantic slope)
species <- survey_melt[(survey_melt$taxa == "M_modesta" & survey_melt$side == "Atlantic"),]
mod_abd <- glm(species$abundance ~ species$site_elevation + species$Max_vol, poisson)
summary(mod_abd) # overdispersion
disp<-glm.poisson.disp(mod_abd)
mod_abd<-update(disp,.~.-z:all)
par(mfrow=c(2,2))
plot(mod_abd)
summary(mod_abd)
Anova(mod_abd)

#######



## 1.2. Extract the slopes of the Poisson models as indices of elevation preference #####

EPs <- function(DF){
  
  DF$brom_ID <- factor(DF$brom_ID)
  DF$side <- factor(DF$side)
  DF$taxa <- factor(DF$taxa)
  
  # calculate observed slope
  # Poisson model
  abd_Atl_obs <- DF$abundance
  elev_Atl_obs <- DF$site_elevation
  maxvol <- DF$Max_vol
  
  mod_abd_Atl_obs <- glm(abd_Atl_obs ~ elev_Atl_obs + maxvol, poisson)
  summary(mod_abd_Atl_obs)
  disp<-glm.poisson.disp(mod_abd_Atl_obs)
  mod_abd_Atl_obs<-update(disp,.~.-z:all)
  # extract slope
  slope_obs <- summary(mod_abd_Atl_obs)$coefficients[2,1] 
  slope_obs_SE <- summary(mod_abd_Atl_obs)$coefficients[2,2]
  
  # n species
  n_species <- sum(DF$abundance)
  
  # put it all in one dataframe
  data.frame(n_species = n_species, slope_obs = slope_obs, slope_obs_SE = slope_obs_SE)
  
}



## EPs for each side
MEOs <- ddply(survey_melt,.(side, taxa), EPs); MEOs

## EPs for both sides combined
MEOs_comb <- ddply(survey_melt,.(taxa), EPs); MEOs_comb

MEOs_comb <- cbind(data.frame(side = (rep("combined", 6))), MEOs_comb)

# combine
MEOs_slopes <- rbind(MEOs, MEOs_comb)


## define which slopes we use for which species for each side
MEOs_slopes$Atlantic_use <- rep(0)
MEOs_slopes$Pacific_use <- rep(0)

MEOs_slopes$Atlantic_use[which(MEOs_slopes$side == "Atlantic")] <- c(1,0,1,0,0,1)
MEOs_slopes$Atlantic_use[which(MEOs_slopes$side == "combined")] <- c(0,0,0,1,1,0)
MEOs_slopes$Pacific_use[which(MEOs_slopes$side == "Pacific")] <- c(1,0,0,0,0,1)
MEOs_slopes$Pacific_use[which(MEOs_slopes$side == "combined")] <- c(0,0,0,1,1,0)
MEOs_slopes

######





## 1.3. Calculate the mean elevation of occurrence (MEO) of each species to check for robustness of analysis (Appendix) ######

# use this to test commands inside function

meos <- function(df){
  # calculate the MEO (based on abundances of species)
  MEO_abd <- sum(df$abundance * df$site_elevation) / sum(df$abundance)
  # SD
  MEO_abd_SD <- sd(rep(df$site_elevation, df$abundance))
  
  data.frame(MEO_abd = MEO_abd, MEO_abd_SD = MEO_abd_SD)
}

MEOs <- ddply(survey_melt,.(side, taxa), meos); MEOs



## remove Dytisctids (no meaningful estimate of EP possible as they are likely linked to lagoon)

MEOs <- MEOs[which(MEOs$taxa != "Dytisctid"),]; MEOs
MEOs$taxa <- factor(MEOs$taxa)
levels(MEOs$taxa) <- c("Culex", "Damsel", "Poly", "Scirtid", "Tipulid")


## plot elevation preferences for Atlantic side
MEOsAt <- MEOs[which(MEOs$side == "Atlantic"),]; MEOsAt

par(mfrow=c(1,2))

plot(MEOsAt$MEO_abd ~ MEOsAt$taxa, xlab = "", ylab = "MEO (m)", border = "white",
     ylim = c(min(MEOsAt$MEO_abd) - max(MEOsAt$MEO_abd_SD), 
              max(MEOsAt$MEO_abd) + max(MEOsAt$MEO_abd_SD)), main = "Atlantic side", cex.axis = 0.9, las=2)
points(as.numeric(MEOsAt$MEO_abd) ~ MEOsAt$taxa, pch=16, cex=1, col = "black")

# add SD
arrows(as.numeric(MEOsAt$taxa), MEOsAt$MEO_abd + MEOsAt$MEO_abd_SD, 
       as.numeric(MEOsAt$taxa), MEOsAt$MEO_abd - MEOsAt$MEO_abd_SD, 
       length = 0.1, angle = 90, code = 3, col="black", lty = 2)


## plot elevation preferences for Pacific side
MEOsPc <- MEOs[which(MEOs$side == "Pacific"),][c(1,3,4,5),]; MEOsPc
MEOsPc$taxa <- factor(MEOsPc$taxa)

plot(MEOsPc$MEO_abd ~ MEOsPc$taxa, xlab = "", ylab = "MEO (m)", border = "white",
     ylim = c(min(MEOsPc$MEO_abd) - max(MEOsPc$MEO_abd_SD), 
              max(MEOsPc$MEO_abd) + max(MEOsPc$MEO_abd_SD)), main = "Pacific side", cex.axis = 0.9, las=2)
points(as.numeric(MEOsPc$MEO_abd) ~ MEOsPc$taxa, pch=16, cex=1, col = "black")

# add SD
arrows(as.numeric(MEOsPc$taxa), MEOsPc$MEO_abd + MEOsPc$MEO_abd_SD, 
       as.numeric(MEOsPc$taxa), MEOsPc$MEO_abd - MEOsPc$MEO_abd_SD, 
       length = 0.1, angle = 90, code = 3, col="black", lty = 2)

######









###---------------- 2. does heat tolerance (CTmax) predict elevation preference (EP) --------------------------

## 2.1. use EP as slope from Poisson model #######
# Atlantic side

data_Atl <- cbind(CT, MEOs_slopes[which(MEOs_slopes$Atlantic_use == 1),][c(1,2,4,5,3),]); data_Atl

## Atlantic side
EP_Atl <- data_Atl$slope_obs
EP_Atl_SE <- data_Atl$slope_obs_SE
CT_Atl <- data_Atl$CTmax_mean
CT_Atl_SE <- data_Atl$CTmax_SE
xlabel <- "Mean CTmax"



# Deming regression with estimated standard errors in both x and y
Atl_dem <- deming(EP_Atl ~ CT_Atl, xstd=CT_Atl_SE, ystd=EP_Atl_SE)
print(Atl_dem)

# calculate p-value from 95% CI
SE <- Atl_dem$ci[4] - Atl_dem$ci[2] /(2*1.96) # calculate SE from 95% CI
Z <- Atl_dem$coefficients[2]/SE; Z # calculate Z from XE
p <- exp(0.717*Z - 0.416*Z^2); p # calculate p from Z
2*pnorm(-abs(Z))

par(mfrow=c(1,1))
par(mar=c(4,5,2,1)+0.1)
plot(EP_Atl ~ CT_Atl, xlim = c(38, 46), ylim = c(-0.005, 0.006),
     main = "Atlantic side", ylab = "Elevation preferrence", xlab = xlabel)
x.error.bar(CT_Atl, EP_Atl, CT_Atl_SE, length = 0.05) 
y.error.bar(CT_Atl, EP_Atl, EP_Atl_SE, length = 0.05) 
text(EP_Atl ~ CT_Atl, labels = data_Atl$species, 
     adj = c(-0.2, -0.6), col = "black", cex = 0.6)
#abline(mod_Atlantic, col = "blue", lwd=2)
abline(Atl_dem, col = "red", lwd=2)
#legend("bottomleft", legend = c("Deming with measured errors","Weighted regression"),  lty = c(1, 1:2), col = c("red", "blue"), cex=0.8, lwd=2)



## Pacific side
data_Pac <- cbind(CT[c(1,3,4,5),], MEOs_slopes[which(MEOs_slopes$Pacific_use == 1),][c(1,3,4,2),]); data_Pac

EP_Pac <- data_Pac$slope_obs
EP_Pac_SE <- data_Pac$slope_obs_SE
CT_Pac <- data_Pac$CTmax_mean
CT_Pac_SE <- data_Pac$CTmax_SE
xlabel <- "Mean CTmax"



## Deming regression with estimated standard errors in both x and y
Pac_dem <- deming(EP_Pac ~ CT_Pac, xstd=CT_Pac_SE, ystd=EP_Pac_SE)
print(Pac_dem)

# calculate p-value from 95% CI
SE <- Pac_dem$ci[4] - Pac_dem$ci[2] /(2*1.96) # calculate SE from 95% CI
Z <- Pac_dem$coefficients[2]/SE; Z # calculate Z from SE
p <- exp(0.717*Z - 0.416*Z^2); p # calculate p from Z
2*pnorm(-abs(Z))

par(mfrow=c(1,1))
par(mar=c(4,5,2,1)+0.1)
plot(EP_Pac ~ CT_Pac, xlim = c(38, 43), ylim = c(-0.004, 0.004),
     main = "Pacific side", ylab = "Elevation preferrence", xlab = xlabel)
x.error.bar(CT_Pac, EP_Pac, CT_Pac_SE, length = 0.05) 
y.error.bar(CT_Pac, EP_Pac, EP_Pac_SE, length = 0.05) 
text(EP_Pac ~ CT_Pac, labels = data_Pac$species, 
     adj = c(-0.2, -0.6), col = "red", cex = 0.6)
#abline(Pac_dem, col = "red", lwd=2)

###############




## 2. use MEO instead of EP (Appendix) ######


# Atlantic side
data_Atl <- cbind(CT, MEOsAt); data_Atl

## Atlantic side
#EP_Atl <- data_Atl$Z_abd; title <- "EP (Z-score)"
EP_Atl <- data_Atl$MEO_abd; title <- "MEO (m)"
CT_Atl <- data_Atl$CTmax_mean
CT_Atl_SE <- data_Atl$CTmax_SE
xlabel <- "Mean CTmax"


mod_Atlantic <- lm(EP_Atl ~ CT_Atl, weights = 1/(CT_Atl_SE^2))
par(mfrow=c(2,2))
plot(mod_Atlantic)
summary(mod_Atlantic)
Anova(mod_Atlantic)



par(mfrow=c(1,1))
par(mar=c(4,5,2,1)+0.1)
plot(EP_Atl ~ CT_Atl, xlim = c(37, 46),
     main = "Atlantic side", ylab = title, xlab = xlabel)
x.error.bar(CT_Atl, EP_Atl, CT_Atl_SE, length = 0.05) 
text(EP_Atl ~ CT_Atl, labels = data_Atl$taxa, 
     adj = c(-0.2, -0.6), col = "black", cex = 0.6)
abline(mod_Atlantic, col = "blue", lwd=2)


# Pacific side
data_Pac <- cbind(CT[c(1,3,4,5),], MEOsPc); data_Pac
data_Pac$taxa <- factor(data_Pac$taxa)



## Pacific side
#EP_Atl <- data_Pac$Z_abd; title <- "EP (Z-score)"
EP_Atl <- data_Pac$MEO_abd; title <- "MEO (m)"
CT_Atl <- data_Pac$CTmax_mean
CT_Atl_SE <- data_Pac$CTmax_SE
xlabel <- "Mean CTmax"


mod_Atlantic <- lm(EP_Atl ~ CT_Atl, weights = 1/(CT_Atl_SE^2))
par(mfrow=c(2,2))
plot(mod_Atlantic)
summary(mod_Atlantic)
Anova(mod_Atlantic)


par(mfrow=c(1,1))
par(mar=c(4,5,2,1)+0.1)
plot(EP_Atl ~ CT_Atl, xlim = c(38, 43),
     main = "Pacific side", ylab = title, xlab = xlabel)
x.error.bar(CT_Atl, EP_Atl, CT_Atl_SE, length = 0.05) 
text(EP_Atl ~ CT_Atl, labels = data_Pac$taxa, 
     adj = c(-0.2, -0.6), col = "black", cex = 0.6)
#abline(mod_Atlantic, col = "blue", lwd=2)

########


