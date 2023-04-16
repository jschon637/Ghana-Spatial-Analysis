#change working directory
setwd("C:/Users/jscho/Documents/Ghana Spatial Work")

library(spdep)
library(maptools)
library(foreign)
library(rgdal)

library(CARBayes)
library(miceadds)
library(tidyverse)

##################################################################
### Maps, GWR, and hotspot analysis were all done in ArcMap
##################################################################

### Load Data Files

#####Districts 2008
dist2008<- readShapeSpatial("Population Weighted/Districts_2008/Volta Variable/20170226_Districts")

#Neighborhood File
dist2008.nb<- poly2nb(dist2008)

#### Convert neighbor list to spatial weights
ghana.dist2008.weights<- nb2listw(dist2008.nb, style="W", zero.policy=T)

################# Translating variables

###Districts 2008
#Count_ is Total Solar
#pov_p_2008 is Percent in poverty
#gini_2008 is Gini Index
#ferat_2008 is Female Ratio
#p_share is NDC 2008
#p_shvol is NDC vote share volatility
#p_turn is Turnout 2008
#p_tuvol is Turnout volatility
#p_ethfr is Ethnic Fractionalization
#Count_3 is World Bank projects
#Density_RD is Road density
#Pop_Densit is population density
#Count_4 is health facilities
#literacy is literacy
#grid_densi is grid density
#grid_perCa is grid per capita
#Count_5 is Lake Volta dummy

#In ArcMap, I used Intersect between districts and electric grid.
#Then, I used Dissolve for grid based on districts. I projected to 1984 WCS Mercator projection system.
#Then, I used "Calculate Geometry" to calculate line lengths.
#Read in the resultant data below.
newer.grid.length<- read.csv("Ghana Electricity Lines/Grid_line_lengths_2.csv")

dist2008.newer<- merge(dist2008,newer.grid.length,by.x = "DIST_2008",
                       by.y= "DIST_2008",all.x=T)

dist2008.newer@data$length[is.na(dist2008.newer@data$length)]<- 0

dist2008.newer@data$gpc<- dist2008.newer@data$length/dist2008.newer@data$pop_2008

################# SLDV models (Table 2)
### Model local, spill-over, and total effects of change in IV
# local effects are the "Direct" effects
# spill-over is the "Indirect"
# Total is the sum of the "Direct" and "Indirect" effects
#############
# Interpret the Direct Effects like OLS coefficients.

dist2008$turninv<- (1-dist2008$p_turn)

#For the robustness check with newer electric grid data
#dist2008.newer$turninv<- (1-dist2008.newer$p_turn)

#Model 1
dist2008.ols.controls<- lm.cluster(Count_~ pov_p_2008 + gini_2008
                                   + ferat_2008 + p_share
                                   + p_shvol + turninv
                                   + p_tuvol + p_ethfr
                                   + Count_3 + Density_RD
                                   + Pop_Densit + Count_4
                                   + literacy + grid_perCa,
                                   data=dist2008,cluster = "REGION")

#Model 2
dist2008.ols.controls.2<- lm.cluster(Count_~ pov_p_2008 + gini_2008
                                     + ferat_2008 + p_share:turninv
                                     + p_share + turninv
                                     + p_shvol 
                                     + p_tuvol + p_ethfr
                                     + Count_3 + Density_RD
                                     + Pop_Densit + Count_4
                                     + literacy + grid_perCa,
                                     data=dist2008,cluster = "REGION")


#Model 3
ghana.dist2008.sldv.3<- lagsarlm(Count_~ pov_p_2008 + gini_2008
                                 + ferat_2008 + p_share
                                 + p_shvol + turninv
                                 + p_tuvol + p_ethfr
                                 + Count_3 + Density_RD
                                 + Pop_Densit + Count_4
                                 + literacy + grid_perCa,
                                 data=dist2008, ghana.dist2008.weights)

#Robustness check with newer grid data
#ghana.dist2008.sldv.altgrid<- lagsarlm(Count_~ pov_p_2008 + gini_2008
#                                       + ferat_2008 + p_share
#                                       + p_shvol + turninv
#                                       + p_tuvol + p_ethfr
#                                       + Count_3 + Density_RD
#                                       + Pop_Densit + Count_4
#                                       + literacy + gpc,
#                                       data=dist2008.newer, ghana.dist2008.weights)


dist2008.impacts.3<- impacts(ghana.dist2008.sldv.3, 
                             listw=ghana.dist2008.weights,
                             R=1000,zstats=TRUE,useHESS = T)

#dist2008.impacts.altgrid<- impacts(ghana.dist2008.sldv.altgrid, 
#                                   listw=ghana.dist2008.weights,
#                                   R=1000,zstats=TRUE,useHESS = T)

summary(dist2008.impacts.3, zstats=TRUE)
#summary(dist2008.impacts.altgrid, zstats=TRUE)






################# Summary statistics (Table A1)
## Min, Max, Mean, N
## Standard deviations were calculated in Excel

#Districts 2008
summary(dist2008)

##### Diagnostic for model selection (Table A2)
# Spatial error vs. Spatially Lagged Dependent Variable

dist2008.ols.2<- lm(Count_~ pov_p_2008 + gini_2008
                    + ferat_2008 + p_share
                    + p_shvol + p_turn
                    + p_tuvol + p_ethfr
                    + Count_3 + Density_RD
                    + Pop_Densit + Count_4
                    + literacy + grid_perCa,
                    data=dist2008)

dist2008.diag.2<- lm.LMtests(dist2008.ols.2, ghana.dist2008.weights, test="all")

#This information goes into Table A2 in the Appendix
summary(dist2008.diag.2)
# SLDV wins


#######################################################
### Not reported

#ghana.dist2008.sldv.1.Durbin<- spatialreg::lagsarlm(Count_~ p_share
#                                      + p_shvol + turninv
#                                      + p_tuvol + p_ethfr
#                                      + Count_3 + Density_RD
#                                      + Pop_Densit + grid_perCa,
#                                      data=dist2008, Durbin=T,
#                                      ghana.dist2008.weights,tol.solve = 1e-15 )

### LM test for Durbin model does not reject the null hypothesis
### This means that spatially lagged independent variables are unnecessary
### Therefore, the spatial lag model is the best fit.
### For more on this, see:
### J. Paul Elhorst (2010) Applied Spatial Econometrics: Raising the Bar, Spatial
### Economic Analysis, 5:1, 9-28, DOI: 10.1080/17421770903541772
#########################################################

##############################################################
### Robustness checks for SLDV in the Appendix (Table A3)
##############################################################

# Robustness check with whether turnout increased from 2004 to 2008
dist2008$turn_inc[dist2008$dist_tur_1>dist2008$dist_turno]<- 1
dist2008$turn_inc[dist2008$dist_tur_1<dist2008$dist_turno]<- 0

#Model 1
ghana.dist2008.sldv.pop2<- lagsarlm(Count_~ pov_p_2008 + gini_2008
                                    + ferat_2008 + p_share
                                    + p_shvol + p_turn
                                    + p_tuvol + p_ethfr
                                    + Count_3 + Density_RD
                                    + Count_4
                                    + literacy + grid_perCa,
                                    data=dist2008, ghana.dist2008.weights)

#Model 2
ghana.dist2008.sldv.volta.2<- lagsarlm(Count_~ Count_5 + pov_p_2008 + gini_2008
                                       + ferat_2008 + p_share
                                       + p_shvol + p_turn
                                       + p_tuvol + p_ethfr
                                       + Count_3 + Density_RD
                                       + Count_4
                                       + literacy + grid_perCa,
                                       data=dist2008, ghana.dist2008.weights)

#Model 3
ghana.dist2008.sldv.volta<- lagsarlm(Count_~ Count_5 + pov_p_2008 + gini_2008
                                     + ferat_2008 + p_share
                                     + p_shvol + p_turn
                                     + p_tuvol + p_ethfr
                                     + Count_3 + Density_RD
                                     + Pop_Densit + Count_4
                                     + literacy + grid_perCa,
                                     data=dist2008, ghana.dist2008.weights)

#Model 4
ghana.dist2008.sldv.inc2<- lagsarlm(Count_~ pov_p_2008 + gini_2008
                                    + ferat_2008 + p_share
                                    + p_shvol + p_turn
                                    + p_tuvol + p_ethfr
                                    + Count_3 + Density_RD
                                    + Pop_Densit + Count_4
                                    + literacy + grid_perCa
                                    + turn_inc,
                                    data=dist2008, ghana.dist2008.weights)

#Model 5
#Robustness check with population
ghana.dist2008.sldv.pop.roaddensity<- lagsarlm(Count_~ pov_p_2008 + gini_2008
                                               + ferat_2008 + p_share
                                               + p_shvol + p_turn
                                               + p_tuvol + p_ethfr
                                               + Count_3 + Density_RD
                                               + pop_2008 + Count_4
                                               + literacy + grid_perCa, tol.solve = 1e-19,
                                               data=dist2008, ghana.dist2008.weights)


dist2008.impacts.r1<- impacts(ghana.dist2008.sldv.volta, listw=ghana.dist2008.weights,
                              R=1000,zstats=TRUE,useHESS = TRUE)
dist2008.impacts.r2<- impacts(ghana.dist2008.sldv.volta.2, listw=ghana.dist2008.weights,
                              R=1000,zstats=TRUE,useHESS = TRUE)
dist2008.impacts.r3<- impacts(ghana.dist2008.sldv.inc2, listw=ghana.dist2008.weights,
                              R=1000,zstats=TRUE,useHESS = TRUE)
dist2008.impacts.r4<- impacts(ghana.dist2008.sldv.pop2, listw=ghana.dist2008.weights,
                              R=1000,zstats=TRUE,useHESS = TRUE)
dist2008.impacts.pop.roaddensity<- impacts(ghana.dist2008.sldv.pop.roaddensity,
                                           listw=ghana.dist2008.weights,
                                           R=1000,zstats=TRUE,useHESS = F)

summary(dist2008.impacts.r1, zstats=TRUE)
summary(dist2008.impacts.r2, zstats=TRUE)
summary(dist2008.impacts.r3, zstats=TRUE)
summary(dist2008.impacts.r4, zstats=TRUE)
summary(dist2008.impacts.pop.roaddensity, zstats=TRUE)


##############################################################
### Corresponding results table for OLS in the Appendix (Table A4)
##############################################################

dist2008.ols.nocontrols<- lm.cluster(Count_~ p_share
                                     + p_shvol + turninv
                                     + p_tuvol + p_ethfr
                                     + Count_3 + Density_RD
                                     + Pop_Densit + grid_perCa,
                                     data=dist2008,cluster = "REGION")

dist2008.ols.controls<- lm.cluster(Count_~ pov_p_2008 + gini_2008
                                   + ferat_2008 + p_share
                                   + p_shvol + turninv
                                   + p_tuvol + p_ethfr
                                   + Count_3 + Density_RD
                                   + Pop_Densit + Count_4
                                   + literacy + grid_perCa,
                                   data=dist2008,cluster = "REGION")

dist2008.ols.nocontrols.2<- lm.cluster(Count_~ p_share:turninv
                                       + p_share + turninv
                                       + p_shvol
                                       + p_tuvol + p_ethfr
                                       + Count_3 + Density_RD
                                       + Pop_Densit + grid_perCa,
                                       data=dist2008,cluster = "REGION")

dist2008.ols.controls.2<- lm.cluster(Count_~ pov_p_2008 + gini_2008
                                     + ferat_2008 + p_share:turninv
                                     + p_share + turninv
                                     + p_shvol 
                                     + p_tuvol + p_ethfr
                                     + Count_3 + Density_RD
                                     + Pop_Densit + Count_4
                                     + literacy + grid_perCa,
                                     data=dist2008,cluster = "REGION")
#Model 1
summary(dist2008.ols.nocontrols)
#Model 2
summary(dist2008.ols.nocontrols.2)
#Model 3
summary(dist2008.ols.controls)
#Model 4
summary(dist2008.ols.controls.2)


#############################################################
### Spatial Count Model (Table A5)
#############################################################

CARBayes.model1.dist<- MVS.CARleroux(Count_~ dist_NDC_3 
                                     + dist_volat + dist_tur_1
                                     + dist_vol_1 + dist_ethfr
                                     + Density_RD + Count_3
                                     + Pop_Densit + pov_p_2008
                                     + gini_2008 + ferat_2008
                                     + Count_4 + literacy
                                     + grid_perCa, data=dist2008, trials = NULL,
                                     W=W.dist,
                                     family = "poisson", burnin = 50000,
                                     n.sample = 150000,
                                     verbose = TRUE)

print(CARBayes.model1.dist)



#####################################################################
### Bivariate regressions, not reported but code shared for interested readers
#####################################################################

bivariate.1.sldv<- lagsarlm(Count_~ literacy,data=dist2008, ghana.dist2008.weights)
bivariate.1.ols<- lm.cluster(Count_~ literacy,data=dist2008,cluster = "REGION")
#negative and significant in both
summary(bivariate.1.ols)
summary(bivariate.1.sldv)

bivariate.2.sldv<- lagsarlm(Count_~ pov_p_2008,data=dist2008, ghana.dist2008.weights)
bivariate.2.ols<- lm.cluster(Count_~ pov_p_2008,data=dist2008,cluster = "REGION")
#positive and significant in OLS; not significant in SLDV
summary(bivariate.2.ols)
summary(bivariate.2.sldv)

bivariate.3.sldv<- lagsarlm(Count_~ gini_2008,data=dist2008, ghana.dist2008.weights)
bivariate.3.ols<- lm.cluster(Count_~ gini_2008,data=dist2008,cluster = "REGION")
#not significant in either
summary(bivariate.3.ols)
summary(bivariate.3.sldv)

bivariate.4.sldv<- lagsarlm(Count_~ ferat_2008,data=dist2008, ghana.dist2008.weights)
bivariate.4.ols<- lm.cluster(Count_~ ferat_2008,data=dist2008,cluster = "REGION")
#negative and significant in both
summary(bivariate.4.ols)
summary(bivariate.4.sldv)

bivariate.5.sldv<- lagsarlm(Count_~ p_share,data=dist2008, ghana.dist2008.weights)
bivariate.5.ols<- lm.cluster(Count_~ p_share,data=dist2008,cluster = "REGION")
#positive and significant in OLS; not significant in SLDV
summary(bivariate.5.ols)
summary(bivariate.5.sldv)

bivariate.6.sldv<- lagsarlm(Count_~ turninv,data=dist2008, ghana.dist2008.weights)
bivariate.6.ols<- lm.cluster(Count_~ turninv,data=dist2008,cluster = "REGION")
#not significant in either
summary(bivariate.6.ols)
summary(bivariate.6.sldv)

bivariate.7.sldv<- lagsarlm(Count_~ p_shvol,data=dist2008, ghana.dist2008.weights)
bivariate.7.ols<- lm.cluster(Count_~ p_shvol,data=dist2008,cluster = "REGION")
#not significant in either
summary(bivariate.7.ols)
summary(bivariate.7.sldv)

bivariate.8.sldv<- lagsarlm(Count_~ p_tuvol,data=dist2008, ghana.dist2008.weights)
bivariate.8.ols<- lm.cluster(Count_~ p_tuvol,data=dist2008,cluster = "REGION")
#positive and significant in both
summary(bivariate.8.ols)
summary(bivariate.8.sldv)

bivariate.9.sldv<- lagsarlm(Count_~ p_ethfr,data=dist2008, ghana.dist2008.weights)
bivariate.9.ols<- lm.cluster(Count_~ p_ethfr,data=dist2008,cluster = "REGION")
#positive and significant in both
summary(bivariate.9.ols)
summary(bivariate.9.sldv)

bivariate.10.sldv<- lagsarlm(Count_~ Count_3,data=dist2008, ghana.dist2008.weights)
bivariate.10.ols<- lm.cluster(Count_~ Count_3,data=dist2008,cluster = "REGION")
#negative and significant only in OLS
summary(bivariate.10.ols)
summary(bivariate.10.sldv)

bivariate.11.sldv<- lagsarlm(Count_~ Density_RD,data=dist2008, ghana.dist2008.weights)
bivariate.11.ols<- lm.cluster(Count_~ Density_RD,data=dist2008,cluster = "REGION")
#negative and significant in both
summary(bivariate.11.ols)
summary(bivariate.11.sldv)

bivariate.12.sldv<- lagsarlm(Count_~ Pop_Densit,data=dist2008, ghana.dist2008.weights)
bivariate.12.ols<- lm.cluster(Count_~ Pop_Densit,data=dist2008,cluster = "REGION")
#negative and significant only in OLS 
summary(bivariate.12.ols)
summary(bivariate.12.sldv)

bivariate.13.sldv<- lagsarlm(Count_~ Count_4,data=dist2008, ghana.dist2008.weights)
bivariate.13.ols<- lm.cluster(Count_~ Count_4,data=dist2008,cluster = "REGION")
#not significant in either
summary(bivariate.13.ols)
summary(bivariate.13.sldv)

bivariate.14.sldv<- lagsarlm(Count_~ grid_perCa,data=dist2008, ghana.dist2008.weights)
bivariate.14.ols<- lm.cluster(Count_~ grid_perCa,data=dist2008,cluster = "REGION")
#not significant in either
summary(bivariate.14.ols)
summary(bivariate.14.sldv)






