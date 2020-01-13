######################################################################################
##########################Beta diversity vs Predictor variables##################
######################################################################################
# load package
library("vegan")
library("dplyr")
library("ggplot2")
library("purrr")
library("tidyr")
library("data.table")
library("gvlma") #offers a way to check the important assumptions on a given linear model

#################################################
# Data loading and preparation
#################################################

# Set working directory 
setwd("C:/Users/DELL/Dropbox (Personal)/anamalais.data/data analysis")


####metadata (predictors) reorganizing, 
#keeping only numerical variable
meta_data.num <-meta_data[,-c(1:5,7:10)]
meta_data.num$logFragSize <- log(meta_data$FragSize)
meta_data.num <-na.omit(meta_data.num) #Removing data with NA 

#Community data
#Adjusting rows of parasite communities data to match with metada
comm.dat.adj<-comm.dat[c(2:4,6:10,12,14,20,21),]

#Geographic distance data 
geo.dat <- read.csv("ana.dist2.csv",row.names=1, na.strings=c(""," ","-9", "-99","NA"))


###############
#Mantel test###
###############

comm.dist <- vegdist(comm.dat.adj, method = "bray")
disturb.dist <- vegdist(meta_data.num, method = "bray")
geo.dist<-vegdist(geo.dat, method="euclidean")

plot(geo.dist,comm.dist,pch=16,cex=0.7,col="black",bty="l",xlab="Spatial distance",ylab="Parasite composition dissimilarity")

#Distance decay through Mantel test
#https://www.rapidtables.com/convert/number/degrees-minutes-seconds-to-degrees.html

mantel(comm.dist, geo.dist, method="pearson", permutations=999) #ns r: 0.2825, p: 0.057

### A partial mantel test to test if parasite communities have more different 
#species compositions as the disturbance dissimilarity in the sites get 
#increasingly different. We will do that while removing any possible spatial 
#autocorrelation. So, the third distance matrix in the analysis will be the 
#spatial distance between sites.
mantel.partial(comm.dist, disturb.dist, geo.dist, method="pearson", permutations=999)

# Based on our results, there is no significant relationship between parasite 
# composition and disturbance composition dissimilarities. Mantel statistic 
#r: -0.01506 p: 0.524 



data(varespec)

## Bray-Curtis distances between samples
dis <- vegdist(varespec)

## First 16 sites grazed, remaining 8 sites ungrazed
groups <- factor(c(rep(1,16), rep(2,8)), labels = c("grazed","ungrazed"))

## Calculate multivariate dispersions
mod <- betadisper(dis, groups)
mod

## Perform test
anova(mod)

## Permutation test for F
permutest(mod, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))
plot(mod.HSD)

## Plot the groups and distances to centroids on the
## first two PCoA axes
plot(mod)
