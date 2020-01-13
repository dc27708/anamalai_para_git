######################################################################################
############## Host-Parasite distribution and networks ################
######################################################################################
# load package
library("vegan")
library("dplyr")
library("ggplot2")
library("igraph")

#################################################
# Data loading and preparation
#################################################

# Set working directory 
setwd("C:/Users/DELL/Dropbox (Personal)/anamalais.data/data analysis")

# Read data into R
p_data <- read.csv("data_parasites.csv",row.names=1, na.strings=c(""," ","-9", "-99","NA")) #Loading data
p_data <- p_data[,-11] #Removing grey junglefowl data

#decostand(x,"pa") counts presence/absence
p_data1 <- (decostand(p_data,"pa"))
p_data1$total<-rowSums(p_data1)
#Removing parasites with no host "Echinococcus_sp" , "Neorickettsia_sp" and "helminthosa" 
p_data2<- p_data1[-which(p_data1$total==0),]
rownames(p_data2)<-sub("Sum of ", "", rownames(p_data2))

write.csv(p_data2, "parasite_distribution")

m <- as.matrix(p_data2[,-27])
dimnames(m) <- list(rownames(p_data2), names(p_data2[,-27]))

# 2: Use iGraph to convert it to a network

n=graph.incidence(m)

# define color and shape mappings.
col <- c("steelblue", "orange")
shape <- c("circle", "square")

plot(n, layout=layout.bipartite,
     vertex.color = col[as.numeric(V(n)$type)+1],
     vertex.shape = shape[as.numeric(V(n)$type)+1],
     asm=0.2
)

# Get simple edge list
edgelist<-as_edgelist(n)
e.dat<-as.data.frame(edgelist)
names(e.dat)<-c("parasites", "hosts")

write.csv(e.dat, "parasite_distr")

# Project
pr <- bipartite.projection(n)
plot(pr[[1]])
plot(pr[[2]])

