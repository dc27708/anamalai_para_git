######################################################################################
##########################Parasite and host Richness models - Anamalai Data##################
######################################################################################
# load package
library("vegan")
library("dplyr")
library("tidyr")
library("tidyverse") 
library("MBESS") #To calculate ci around effect size (adj R2)
#http://yatani.jp/teaching/doku.php?id=hcistats:linearregression
library(MASS) # v7.3-50
library(car) # v3.0-2
library(lme4) # v1.1-18-1
library(mlmRev) # v1.0-6
library(ggplot2) # v3.0.0
library("glmmADMB")
library(MuMIn)
library(jtools)
library(summarytools)

#################################################
# Data loading and preparation
#################################################

# Set working directory 
setwd("C:/Users/DELL/Dropbox (Personal)/anamalais.data/data analysis")

# Read data into R
ana_data <- read.csv("data_v1.csv", na.strings=c(""," ","-9", "-99","NA")) #Loading data
ana_data <- ana_data[!ana_data$HostSp1 == "grey.junglefowl", ] #Removing grey junglefowl data
ana_data <- ana_data[!ana_data$HostSp1 == "cattle", ] #Removing cattle data
ana_data <- ana_data[!ana_data$HostSp1 == "domestic.goat", ] #Removing grey junglefowl data

ana_data<- ana_data[,-c(2,5,9)] #removing RangeName and VegType

ana_data <-na.omit(ana_data) #Removing data with NA #dim(ana_data) [1] 3615,48
ana_data[,8:48]<-ifelse(ana_data[,8:48]>0,1,0)
ana_data$HostSp1<-factor(ana_data$HostSp1)
ana_data<-ana_data[complete.cases(ana_data),]

#Creating new valiable for host and fragment combination
for(i in 1:nrow(ana_data)){
  
  ana_data$host_frag[i]<-paste(ana_data$HostSp1[i], ana_data$FragName[i], sep = "-")  
}
ana_data$host_frag<-as.factor(ana_data$host_frag)

#Creating parasite abundance matrix from raw data (sample level)
abun.mat <- ana_data[,8:48]  #dim(abun.mat) [1] 3615,41

# Read data into R
#abun.host <- read.csv("host1.csv",row.names=1, na.strings=c(""," ","-9", "-99","NA")) #Loading host data
ana_host<-ana_data[,-c(8:48)]
ana_host$value <-rowSums(ana_data[,c(8:48)])
ana_host$AbrFragSamID<-factor(ana_host$AbrFragSamID)

#Fragmentwise host abundance
abun.host <- spread(ana_host,HostSp1, value)
abun.host[is.na(abun.host)] <- 0
abun.host[,8:30]<-ifelse(abun.host[,8:30]>0,1,0)

#Read Metadata in to R
#Metadata that include predictor variables 
meta_data <- read.csv("meta_dat1.csv",row.names=1, na.strings=c(""," ","-9", "-99","NA"))
meta_data<-meta_data[,c("FragSize", "AvgIsoDist", "FragType", "VegType", "cattle", "plantation", "settlement")]
meta_data$frag<-rownames(meta_data)

#####################################################################################
############### Extrapolated Richness Estimators (nonparametric)#####################
#####################################################################################

#Sampling effort(sample size) may influence observed richness values 
#so extrapolated richness were estimated (using specpool) 
#to account for rare parasites in each fragment

####Bootstrap estimate of parasite richness for each fragment
par.pool<-with(ana_data, specpool(abun.mat, FragName)) #hostwise Parasite richness for each fragment

boot.p<-par.pool[,c("boot", "boot.se", "n")]
boot.p<-na.omit(boot.p)
colnames(boot.p)<- c("boot.p.rich", "boot.p.se", "p.n")
boot.p$lower.ci<-boot.p$boot.p.rich - qt(1 - (0.05 / 2), boot.p$p.n - 1) * boot.p$boot.p.se
boot.p$upper.ci<-boot.p$boot.p.rich + qt(1 - (0.05 / 2), boot.p$p.n - 1) * boot.p$boot.p.se
boot.p$frag<-rownames(boot.p)


#merging boot.p and meta_data
boot.p<-merge(boot.p, meta_data, by = "frag", all = FALSE)

write.csv(boot.p, file="boot.p.csv")###To calculate ci for parasite richness manually in excel
####New file saved as boot.p.ci

####Bootstrap estimate of parasite richness for each host species in each fragment

fspool<-with(ana_data, specpool(abun.mat, host_frag)) #hostwise Parasite richness for each fragment

boot.fspool<-fspool[,(7:8)]
boot.fspool$host_frag<-factor(rownames(boot.fspool))
boot.fspool<-separate(boot.fspool,"host_frag",sep = "-", c("host","frag"))
colnames(boot.fspool)<- c("par.rich", "par.rich.se","host", "frag")

####Bootstrap estimate of host richness for each fragment
hspool<-with(ana_data, specpool(abun.host, FragName)) #Host richness in each fragment


##Creating fragment column 
hspool$frag<-rownames(hspool)
###Subsetting boot data
boot.hspool<-hspool[, c("boot","boot.se","frag")]
colnames(boot.hspool)<-c("host.rich","host.rich.se","frag")

####Writing host richness estimates to file
write.csv(hspool, file = "Estimates of host species richness.csv")


#merging boot.hspool and boot.fspool to creat tot_data
tot_data<-merge(boot.fspool, boot.hspool, by = "frag", all = FALSE)

#merging tot_data and meta_data
total_data<-merge(tot_data, meta_data, by = "frag", all = FALSE)

total_data <-na.omit(total_data)

########################Mixed Models for Parasite Richness#################
###########################################################################
#http://rstudio-pubs-static.s3.amazonaws.com/263877_d811720e434d47fb8430b8f0bb7f7da4.html
#Visualization https://strengejacke.wordpress.com/2014/11/18/visualizing-generalized-linear-mixed-effects-models-part-2-rstats-lme4/
#http://atyre2.github.io/2017/06/16/rebutting_cade.html


hist(total_data$par.rich, freq=F) #Looks like Poisson distribution

dat<-total_data[, c("par.rich","host.rich","FragSize","AvgIsoDist", "frag", "host", "cattle", "plantation", "settlement")]
dat$host<-as.factor(dat$host)
dat$frag<-as.factor(dat$frag)

#Subsetting body size variable
body<-ana_data[,c("HostSp1", "SizeInG")]
body<-na.omit(body)
colnames(body)<-c("host", "mass")
body<-body[!duplicated(body$host),]
body$in.kg<-body$mass/1000

####adding body size variable to dat
dat1<-left_join(dat,body)

###Linear mixed model selection
#https://uoftcoders.github.io/rcourse/lec11-model-selection.html
#https://sites.google.com/site/rforfishandwildlifegrads/home/mumin_usage_examples
#https://www.andrew.cmu.edu/user/achoulde/94842/lectures/lecture10/lecture10-94842.html

###Dredging models
lmer_full<- lme4::lmer(par.rich ~ host.rich+in.kg+settlement+cattle+plantation+FragSize+
                   AvgIsoDist+(1|host)+(1|frag), data =dat1, REML=FALSE)
options(na.action = "na.fail") # Require for dredge to run
lmer_full_dredge <- dredge(lmer_full, beta = F, evaluate = T, rank = AICc)
options(na.action = "na.omit") # set back to default
head(lmer_full_dredge, 10)
nrow(lmer_full_dredge)

plot(lmer_full_dredge)
importance(lmer_full_dredge)

top_model <- get.models(lmer_full_dredge, subset = 1)[[1]]
top_model
# Summarize top model
summary(top_model)

#Model averaging
#http://atyre2.github.io/2017/06/16/rebutting_cade.html
top_model_avg<-model.avg(lmer_full_dredge, subset = delta <= 2)
summary(model.avg(lmer_full_dredge, subset = delta <= 2))
confint(top_model_avg)

#Diagnostics
shapiro.test(resid(top_model_avg))
plot(fitted(top_model),resid(top_model,type="pearson"))# this will create the plot
abline(0,0, col="red")
qqnorm(resid(top_model)) 
qqline(resid(top_model), col= "red") # add a perfect fit line
#The residuals are normally distributed.

######Output from jtools
#https://cran.r-project.org/web/packages/jtools/vignettes/summ.html
#https://strengejacke.wordpress.com/2014/10/26/visualizing-generalized-linear-mixed-effects-models-with-ggplot-rstats-lme4/

###Alternative model selection
top_model.d<-lme4::lmer(par.rich ~ cattle + settlement+ plantation+ host.rich+in.kg+(1|host) + (1|frag),REML = FALSE, dat1)
top_model.f<-lme4::lmer(par.rich ~ FragSize+AvgIsoDist+host.rich+in.kg+(1|host) + (1|frag),REML = FALSE, dat1)

summary(top_model.d)
summary(top_model.f)

out.put.p<-model.sel(top_model.d,top_model.f, rank = AIC)
export_summs(top_model.d,top_model.f,error_format = "[{conf.low}, {conf.high}]", digits = 3)
plot_summs(top_model.d,top_model.f,
           coefs = c("Livestock presence" = "cattle", "Human settlement presence" = "settlement",
                     "Plantation presence" = "plantation", 
                     "Host richness" = "host.rich", "Host bodymass (kg)" = "in.kg",
                     "Fragment size" = "FragSize", "Isolation index" = "AvgIsoDist"),
           model.names = c("Human disturbance model", "Habitat fragmentation model"),
           scale = TRUE, inner_ci_level = .9, legend.title = "Parasite richness models")


# coerce the object out.put into a data frame elements 6-10 in out.put have what we want
sel.table<-as.data.frame(out.put.p)[9:13]
sel.table
# a little clean-up, lets round things a bit
sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)
# that's better 
# how about a little renaming columns to fit proper conventions
# number of parameters (df) should be K
names(sel.table)[1] = "K"

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)

# replace Model name with formulas little tricky so be careful
for(i in 1:nrow(sel.table)) sel.table$Model[i]<- as.character(formula(paste(sel.table$Model[i])))[2]
#little reordering of columns
sel.table<-sel.table[,c(6,1,2,3,4,5)]
# let's see what is in there
sel.table 
# write to a file, here a comma separated values format

write.csv(sel.table,"par.rich model selection table.csv", row.names = T)

####Effect size
r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

r2.corr.mer(top_model1) #0.4817457
r2.corr.mer(top_model2) #0.4819008
r2.corr.mer(top_model3) #0.4837024
r2.corr.mer(top_model4) #0.4838355

r2.corr.mer(top_model.d) #0.477413
r2.corr.mer(top_model.b) #0.4974764
r2.corr.mer(top_model.h) #0.4913387
r2.corr.mer(top_model.f) #0.4960647

#Plot
ggplot(dat1, aes(x=FragSize, y=par.rich,color=plantation)) +
  geom_point(size=2) +
  labs(x = "FragSize",y = "Bootstrap estimate of parasite taxon richness")+
  geom_smooth(method=lm,se=T)+
  theme_classic2()



#######################################Linear Models for Host Richness#################
#######################################################################################

#merging boot.hspool and meta_data to creat tot_data
dath<-merge(boot.hspool, meta_data, by = "frag", all = FALSE)
dath<-dath[,c("host.rich", "FragSize", "AvgIsoDist", "cattle", "plantation", "settlement")]
dath<- na.omit(dath)

hist(dath$host.rich, freq=F) #Looks like normally distributed

lm_full<- lm(host.rich~FragSize+AvgIsoDist+settlement+cattle+plantation, data=dath)
options(na.action = "na.fail") # Require for dredge to run
lm_full_dredge <- dredge(lm_full, beta = F, evaluate = T, rank = AICc)
options(na.action = "na.omit") # set back to default

head(lm_full_dredge)
nrow(lm_full_dredge)

top_model_lm <- get.models(lm_full_dredge, subset = 1)[[1]]
top_model_lm
# Summarize top model
summary(top_model_lm)
#Model averaging
top_model_avg_lm<-model.avg(lm_full_dredge, subset = delta <= 2)
summary(top_model_avg_lm)
confint(top_model_avg_lm)

#Diagnostics
shapiro.test(resid(top_model_lm)) # Test for the assumption of normally distributed residuals
qqnorm(resid(top_model_lm)) # plot of the normal theoretical quantiles against the exponential data quantiles
qqline(resid(top_model_lm)) # the residuals should fall along this line if they are normally distributed
#The residuals are close to normally distributed.

###Top two models, not one
top.models.lm<-get.models(lm_full_dredge, subset = delta <= 2)
top.models.lm

top_model.lm1<-lm(host.rich ~ AvgIsoDist + cattle + plantation, dath)
top_model.lm2<-lm(host.rich ~ cattle + FragSize + plantation, dath)

export_summs(top_model.lm1,top_model.lm2,error_format = "[{conf.low}, {conf.high}]", digits = 3)
plot_summs(top_model.lm1,top_model.lm2, scale = TRUE)

###Alternative model selection strategy
lm_f<- lm(host.rich~FragSize+AvgIsoDist, data=dath)
lm_d<- lm(host.rich~settlement+cattle+plantation, data=dath)

out.put.h<-model.sel(lm_f,lm_d, rank = AIC)
export_summs(lm_f,lm_d,error_format = "[{conf.low}, {conf.high}]", digits = 3)
plot_summs(lm_f,lm_d,
           coefs = c("Livestock presence" = "cattle", "Human settlement presence" = "settlement",
                     "Plantation presence" = "plantation","Fragment size" = "FragSize", "Isolation index" = "AvgIsoDist"),
           model.names = c("Human disturbance model", "Habitat fragmentation model"),
           scale = TRUE, inner_ci_level = .9, legend.title = "Host richness models")

# coerce the object out.put into a data frame elements 6-10 in out.put have what we want
sel.table.h<-as.data.frame(out.put.h)[7:11]
sel.table.h
# a little clean-up, lets round things a bit
sel.table.h[,2:3]<- round(sel.table.h[,2:3],2)
sel.table.h[,4:5]<- round(sel.table.h[,4:5],3)
# that's better 
# how about a little renaming columns to fit proper conventions
# number of parameters (df) should be K
names(sel.table.h)[1] = "K"

## lets be sure to put the model names in a column
sel.table.h$Model<-rownames(sel.table.h)

# replace Model name with formulas little tricky so be careful
for(i in 1:nrow(sel.table.h)) sel.table.h$Model[i]<- as.character(formula(paste(sel.table.h$Model[i])))[2]
#little reordering of columns
sel.table.h<-sel.table.h[,c(6,1,2,3,4,5)]
# let's see what is in there
sel.table.h 
# write to a file, here a comma separated values format

write.csv(sel.table.h,"host.rich model selection table.csv", row.names = T)




ggplot(dath, aes(x = cattle, y = host.rich)) +
  geom_boxplot(fill = "grey80", colour = "blue") +
  scale_x_discrete() + xlab("Settlement") +
  ylab("Bootstrap estimate of host species richness")+ggtitle("host.rich~cattle")+
  theme_classic()



###########Missing parasites in undisturbed fragments#########################
##############################################################################

# Read data into R
mis.dat <- read.csv("missing.par.csv", na.strings=c(""," ","-9", "-99","NA")) #Loading data

#Converting to 0 and 1
mis.dat[,5:45]<-ifelse(mis.dat[,5:45]>0,1,0)
mis.dat<-na.omit(mis.dat) #Removing NA



#Subsetting plantation data
mis.plant<-mis.dat[,c(2,5:45)]
mis.plant<-mis.plant[, colSums(mis.plant != 0) > 0]#Removing columns with all 0.
##Names of parasites only in disturbed (plantation)
mis.plant.abs<-filter(mis.plant, plantation =="absent")
names(mis.plant.abs)[sapply(mis.plant.abs, setequal, 0)]

##Names of parasites only in undisturbed (plantation)
mis.plant.pres<-filter(mis.plant, plantation =="present")
names(mis.plant.pres)[sapply(mis.plant.pres, setequal, 0)]

#Subsetting livestock data
mis.cat<-mis.dat[,c(3,5:45)]
mis.cat<-mis.cat[, colSums(mis.cat != 0) > 0]#Removing columns with all 0.
##Names of parasites only in disturbed (livestock)
mis.cat.abs<-filter(mis.cat, livestock =="absent")
names(mis.cat.abs)[sapply(mis.cat.abs, setequal, 0)]
##Names of parasites only in undisturbed (livestock)
mis.cat.pres<-filter(mis.cat, livestock =="present")
names(mis.cat.pres)[sapply(mis.cat.pres, setequal, 0)] #paragonimus_sp

#Subsetting settlement data
mis.set<-mis.dat[,c(4,5:45)]
mis.set<-mis.set[, colSums(mis.set != 0) > 0]#Removing columns with all 0.
##Names of parasites only in disturbed (settlement)
mis.set.abs<-filter(mis.set, settlement =="absent")
names(mis.set.abs)[sapply(mis.set.abs, setequal, 0)]
##Names of parasites only in undisturbed (settlement)
mis.set.pres<-filter(mis.set, settlement =="present")
names(mis.set.pres)[sapply(mis.set.pres, setequal, 0)] #Sarcocystis_sp


#####Host species present in disturbed fragments
host.dat<-dat[c("host", "cattle")]
host.dat.pres<-filter(host.dat,cattle=="present") #indian.porcupine, hanuman.langur
unique(host.dat.pres$host)
host.dat.abs<-filter(host.dat,cattle=="absent")#spotted.deer, otter, sambar.deer
unique(host.dat.abs$host)

host.set<-dat[c("host", "settlement")]
host.set.pres<-filter(host.set,settlement=="present") #
unique(host.set.pres$host)
host.set.abs<-filter(host.set,settlement=="absent") #
unique(host.set.abs$host)





################Coinfection in primates
pri.host<-ana_data[,-c(8:48)]
pri.host$value <-rowSums(ana_data[,c(8:48)])
pri.host$AbrFragSamID<-factor(pri.host$AbrFragSamID)

pri.host %>%
  filter(HostSp1 == c("bonnet.macaque","liontailed.macaque","nilgiri.langur")) %>%
  ggplot(aes(x = value, fill = HostSp1))+
  geom_histogram(binwidth=1)+
  labs(x = "No. of parasite taxa/sample", y = "Frequency")+
  facet_wrap(HostSp1 ~ .)+
  geom_vline(aes(xintercept=mean(value)),color="blue", linetype="dashed", size=1)+
  theme_nice()+
  theme(legend.position = "none")


  















#######Scatter Plotting host and parasite richness for different factors
###############################################################

#With plantation
p1 <- ggplot(total_data, aes(x=host.rich, y=par.rich, color=plantation,palette = "jco")) +
  geom_point(size=2) +
  labs(x = "Bootstrap estimate of host species richness",y = "Bootstrap estimate of parasite taxon richness", title = "3a")+
  scale_color_manual(values=c("green", "brown"))+
  guides(color=guide_legend(title="Plantation"))+
  geom_smooth(method=lm,   
              se=T)
#With livestock
p3 <- ggplot(total_data, aes(x=host.rich, y=par.rich, color=cattle)) +
  geom_point(size=2) +
  labs(x = "Bootstrap estimate of host species richness",y = "Bootstrap estimate of parasite taxon richness", title = "3b")+
  scale_color_manual(values=c("green", "brown"))+
  guides(color=guide_legend(title="Livestock"))+
  geom_smooth(method=lm,   # Add linear regression line, ci shaded region
              se=T)
#With settlements
p2 <- ggplot(total_data, aes(x=host.rich, y=par.rich, color=settlement)) +
  geom_point(size=2) +
  labs(x = "Bootstrap estimate of host species richness",y = "Bootstrap estimate of parasite taxon richness", title = "3c")+
  scale_color_manual(values=c("green", "brown"))+
  guides(color=guide_legend(title="Settlement"))+
  geom_smooth(method=lm,   # Add linear regression line, no ci shaded region
              se=T)  



#With Fragment Size
p4 <- ggplot(total_data, aes(x=FragSize, y=par.rich)) +
  geom_point(size=2) +
  labs(x = "Fragment size (ha)",y = "Bootstrap estimate of parasite taxon richness")+
  geom_smooth(method=lm,   # Add linear regression line, no ci shaded region
              se=T)

ggarrange(p1, p2, p3, p4+ rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

#With Isolation distance
p5 <- ggplot(total_data, aes(x=AvgIsoDist, y=par.rich)) +
  geom_point(size=2) +
  labs(x = "Isolation index",y = "Bootstrap estimate of parasite taxon richness")+
  geom_smooth(method=lm,   # Add linear regression line, no ci shaded region
              se=T)

####Boxplotting Categorical variables for both host and parasite richness
#Finally, using boxplots, checking relationships 
#between response variables and categorical predictor variables
#One way ANOVA to test significance
#http://www.sthda.com/english/wiki/one-way-anova-test-in-r
ggplot(total_data, aes(x = plantation, y = host.rich)) +
  geom_boxplot(fill = "grey80", colour = "blue") +
  scale_x_discrete() + xlab("Plantation") +
  ylab("Bootstrap estimate of host species richness")+ggtitle("a")

ggplot(total_data, aes(x = cattle, y = host.rich)) +
  geom_boxplot(fill = "grey80", colour = "blue") +
  scale_x_discrete() + xlab("Livestock") +
  ylab("Bootstrap estimate of host species richness")+ggtitle("b")

ggplot(total_data, aes(x = settlement, y = host.rich)) +
  geom_boxplot(fill = "grey80", colour = "blue") +
  scale_x_discrete() + xlab("Settlement") +
  ylab("Bootstrap estimate of host species richness")+ggtitle("c")

#########################Parasite Diversity Plott#################################
##################################################################################

#reading frag mean and ci para.rich--to shade in geom_ribbon 
boot.p.ci <- read.csv("boot.p.ci.csv", na.strings=c(""," ","-9", "-99","NA"))


####plotting mean and CI par.rich for plantation
boot.p%>%
  ggplot(aes(x=frag, y=boot.p.rich,color=plantation)) + 
  labs(x="Forest fragments", y = "Bootstrap estimate of parasite taxon richness", title = "S1a plantation")+
  geom_point()+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.4) +
  geom_hline(yintercept = c(boot.p.ci$frag.mean[3], boot.p.ci$frag.mean[4]))+
  geom_text(aes(2,10,label = "mean (Absent)"), color = "black")+
  geom_text(aes(1,32,label = "mean (Present)"), color = "black")+
  coord_flip()+
  theme_bw()

####plotting mean and CI par.rich for Livestock
boot.p%>%
  ggplot(aes(x=frag, y=boot.p.rich, color=cattle)) + 
  labs(x="Forest fragments", y = "Bootstrap estimate of parasite taxon richness", title = "S1b Livestock")+
  geom_point()+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.4) +
  geom_hline(yintercept = c(boot.p.ci$frag.mean[1], boot.p.ci$frag.mean[2]))+
  geom_text(aes(1,12,label = "mean (Absent)"), color = "black")+
  geom_text(aes(1,32,label = "mean (Present)"), color = "black")+
  coord_flip()+
  theme_bw()

#geom_ribbon(aes(ymin=boot.p.ci$frag.lower.ci[1], ymax=boot.p.ci$frag.upper.ci[1]) ,fill="green", alpha=0.2)+ #absent
#geom_ribbon(aes(ymin=boot.p.ci$frag.lower.ci[2], ymax=boot.p.ci$frag.upper.ci[2]),fill="brown", alpha=0.2)+ #present


####plotting mean and CI par.rich for settlement
boot.p%>%
  ggplot(aes(x=frag, y=boot.p.rich, color=settlement)) + 
  labs(x="Forest fragments", y = "Bootstrap estimate of parasite taxon richness", title = "S1c Settlement")+
  geom_point()+
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.4) +
  geom_hline(yintercept = c(boot.p.ci$frag.mean[5], boot.p.ci$frag.mean[6]))+
  geom_text(aes(1,10,label = "Overall mean (Absent)"), color = "black")+
  geom_text(aes(1,32,label = "Overall mean (Present)"), color = "black")+
  coord_flip()+
  theme_bw()




























############################Fitting Multiple Linear Regression Model ###################
###############################################################################
##https://tutorials.iq.harvard.edu/R/Rstatistics/Rstatistics.html

#######Parasite richness model
#Variable selection
#http://r-statistics.co/Model-Selection-in-R.html
#http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/154-stepwise-regression-essentials-in-r/

# Fit the full model 
full.model.para <- lm(par.rich~host.rich*FragSize*plantation*cattle*settlement+FragSize, total_data)
summary(full.model.para)

confint(full.model.para, level = 0.95)

#Plotting of confidence intervals of coefficients
y = coef(full.model.para)
x = seq_along(y)
ci = confint(full.model.para)
xlim = range(x) + c(-0.5,0.2)
ylim = range(ci)
ylim = ylim + 0.1*c(-1,+1)*diff(ylim) # extend it a little
ylab = bquote(hat(beta))
xlab = "coefficient"
par(mar=c(4.5,12,1,1), las=1)
plot(y, pch=16, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", bty="n")
axis(1, at=x, labels=names(y), tick=FALSE)
abline(h=0, lty=3)
arrows(x,ci[,1],x,ci[,2], code=3, angle=90, length=0.05)
title(main = "Parasite.Richness.model")
### horizontal layout:
plot(x=y, y=x, pch=16, xlim=ylim, ylim=xlim, xlab=ylab, ylab="", yaxt="n", bty="n")
axis(2, at=x, labels=names(y), tick=FALSE)
abline(v=0, lty=3)
arrows(ci[,1],x,ci[,2],x, code=3, angle=90, length=0.05)
title(main = "Parasite.Richness.model")


#The unstandardized coefficients of the model are themselves 
#measures of effect size (and arguably better than standardized 
#ones). They tell you how big a change in Y is associated with 
#a change in X of 1 unit. For dummy variables this becomes a 
#difference in means between two groups (or the two values of 
#the dummy variable).
#Standardized effect size can be misleading even when it 
#appears to put effects on a common scale. For example, study 
#A has a standardized effect of 0.5 and study B 0.4. Can one 
#with confidence state that the effect in A is bigger than in B? No, not in general. For example if the SD of B was 2 and the SD of A was 1.2 then the mean difference in A is 0.60 and the mean difference in B is 0.8. (This matters because the SD in each study is influenced by factors that have no bearing on how big the effect is: measurement error, range restriction and so on).
#Standardized effect size is more useful for power calculations 
#but generally the effect size used for power calculations 
#includes all sources of error in the study, whereas the effect 
#size relevant to practical significance requires the sources of
#error to be constant for all comparisons (i.e., not based on 
#the characteristics of a study but on the characteristics of 
#the population of interest).

# Stepwise regression model
step.model.para <- step(full.model.para, direction = "forward", trace = FALSE)
summary(step.model.para)

#Multicollinearity and Statistical Significance
all_vifs <- car::vif(step.model.para)
print(all_vifs)

#Removing conlinear variables (VIF>4)
signif_all <- names(all_vifs)
# Remove vars with VIF> 4 and re-build model until none of VIFs don't exceed 4.
while(any(all_vifs > 4)){
  var_with_max_vif <- names(which(all_vifs == max(all_vifs)))  # get the var with max vif
  signif_all <- signif_all[!(signif_all) %in% var_with_max_vif]  # remove
  myForm <- as.formula(paste("boot ~ ", paste (signif_all, collapse=" + "), sep=""))  # new formula
  step.model.para <- lm(myForm, data=total_data)  # re-build model with new formula
  all_vifs <- car::vif(step.model.para)
}
#The best model had Adjusted R-squared:  0.2441
##Presence of plantation had significant effect (P=0.012)

#One way ANOVA to test significance
#http://www.sthda.com/english/wiki/one-way-anova-test-in-r


##################################Host Richness Model#########
##############################################################
####This is fine as well.
# Fit the full model 
full.model.host <- lm(boot.host~FragSize+plantation+cattle+settlement+AvgIsoDist, total_data)
summary(full.model.host)

confint(full.model.host, level = 0.95)

#Plotting of confidence intervals of coefficients
y = coef(full.model.host)
x = seq_along(y)
ci = confint(full.model.host)
xlim = range(x) + c(-0.5,0.2)
ylim = range(ci)
ylim = ylim + 0.1*c(-1,+1)*diff(ylim) # extend it a little
ylab = bquote(hat(beta))
xlab = "coefficient"
par(mar=c(4.5,12,1,1), las=1)
plot(y, pch=16, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", bty="n")
axis(1, at=x, labels=names(y), tick=FALSE)
abline(h=0, lty=3)
arrows(x,ci[,1],x,ci[,2], code=3, angle=90, length=0.05)
title(main = "Host.Richness.model")
### horizontal layout:
plot(x=y, y=x, pch=16, xlim=ylim, ylim=xlim, xlab=ylab, ylab="", yaxt="n", bty="n")
axis(2, at=x, labels=names(y), tick=FALSE)
abline(v=0, lty=3)
arrows(ci[,1],x,ci[,2],x, code=3, angle=90, length=0.05)
title(main = "Host.Richness.model")

# Stepwise regression model
step.model.host <- step(full.model.host, direction = "forward", trace = FALSE)
summary(step.model.host)

#Multicollinearity and Statistical Significance
all_vifs <- car::vif(step.model.host)
print(all_vifs)

#Removing conlinear variables (VIF>4)
signif_all <- names(all_vifs)
# Remove vars with VIF> 4 and re-build model until none of VIFs don't exceed 4.
while(any(all_vifs > 4)){
  var_with_max_vif <- names(which(all_vifs == max(all_vifs)))  # get the var with max vif
  signif_all <- signif_all[!(signif_all) %in% var_with_max_vif]  # remove
  myForm <- as.formula(paste("boot ~ ", paste (signif_all, collapse=" + "), sep=""))  # new formula
  step.model.host <- lm(myForm, data=total_data)  # re-build model with new formula
  all_vifs <- car::vif(step.model.host)
}
summary(step.model.host)




































###########################################old disturbance data##########
########################################################################
###Parasite richness
#Partial model with old disturbance data
total_data.boot1<- total_data[,c("boot","MeanNumTrees", "MeanHghtTrees",
                                 "MeanPcentCanCov" , "MeanBasalAreaPerPlot", 
                                 "MeanPcentShrubCov",    "MeanNumTreesCut", 
                                 "MeanNumLianas")]
total_data.boot1<-na.omit(total_data.boot1)
# Fit the full model 
full.model.boot1 <- lm(boot ~., data = total_data.boot1)
confint(full.model.boot1, level = 0.95)
# Stepwise regression model
step.model.boot1 <- step(full.model.boot1, direction = "forward", trace = FALSE)
summary(step.model.boot1)
#The best model had Adjusted R-squared:  0.5069
##None of the predictor variables had significant effect


#Plotting of confidence intervals of coefficients
y = coef(full.model.boot1)
x = seq_along(y)
ci = confint(full.model.boot1)
xlim = range(x) + c(-0.5,0.2)
ylim = range(ci)
ylim = ylim + 0.1*c(-1,+1)*diff(ylim) # extend it a little
ylab = bquote(hat(beta))
xlab = "coefficientl"
par(mar=c(4.5,12,1,1), las=1)
plot(y, pch=16, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", bty="n")
axis(1, at=x, labels=names(y), tick=FALSE)
abline(h=0, lty=3)
arrows(x,ci[,1],x,ci[,2], code=3, angle=90, length=0.05)
title(main = "Parasite.VegDisturb.model")
### horizontal layout:
plot(x=y, y=x, pch=16, xlim=ylim, ylim=xlim, xlab=ylab, ylab="", yaxt="n", bty="n")
axis(2, at=x, labels=names(y), tick=FALSE)
abline(v=0, lty=3)
arrows(ci[,1],x,ci[,2],x, code=3, angle=90, length=0.05)
title(main = "Parasite.VegDisturb.model")

###Host richness
#Partial model with old disturbance data
total_data.host1<- total_data[,c("logHostRich","MeanNumTrees", "MeanHghtTrees",
                                 "MeanPcentCanCov" , "MeanBasalAreaPerPlot", 
                                 "MeanPcentShrubCov",    "MeanNumTreesCut", 
                                 "MeanNumLianas")]
total_data.host1<-na.omit(total_data.host1)
# Fit the full model 
full.model.host1 <- lm(logHostRich ~., data = total_data.host1)
confint(full.model.host1, level = 0.95)
# Stepwise regression model
step.model.host1 <- step(full.model.host1, direction = "forward", trace = FALSE)
summary(step.model.host1)

#Plotting of confidence intervals of coefficients
y = coef(full.model.host1)
x = seq_along(y)
ci = confint(full.model.host1)
xlim = range(x) + c(-0.5,0.2)
ylim = range(ci)
ylim = ylim + 0.1*c(-1,+1)*diff(ylim) # extend it a little
ylab = bquote(hat(beta))
xlab = "coefficientl"
par(mar=c(4.5,12,1,1), las=1)
plot(y, pch=16, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", bty="n")
axis(1, at=x, labels=names(y), tick=FALSE)
abline(h=0, lty=3)
arrows(x,ci[,1],x,ci[,2], code=3, angle=90, length=0.05)
title(main = "Host.VegDisturb.model")
### horizontal layout:
plot(x=y, y=x, pch=16, xlim=ylim, ylim=xlim, xlab=ylab, ylab="", yaxt="n", bty="n")
axis(2, at=x, labels=names(y), tick=FALSE)
abline(v=0, lty=3)
arrows(ci[,1],x,ci[,2],x, code=3, angle=90, length=0.05)
title(main = "Host.VegDisturb.model")


##################END#####################
##########################################







































































###################################################################################
####################################################################################
#Extra codes

#Checking scatterplots of response to predictor variables
#to detect linear relationships graphically between continuous variables
#Smoothing applied.  method = "auto" the smoothing method is chosen based on the size 
#of the largest group (across all panels).
#https://drsimonj.svbtle.com/plot-some-variables-against-many-others

#How to interprete smooth scatterplots
#https://stats.stackexchange.com/questions/301428/how-do-i-interpret-this-scatter-plot

total_data.bootVegdisturb <- total_data[,c(8, 21:28)]
total_data.logFragSize <- total_data[,c(3,5,7,8,29)]
total_data.AvgIsoDist <- total_data[,c(3,5,7,8,12)]

#scatter plot: boot_v_veg disturbance variables
total_data.bootVegdisturb %>%
  gather(-boot, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = boot)) +
  geom_point(colour = "#D55E00", size= 2) +
  stat_smooth() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

#scatter plot: logFragSize_parasite rich estimates
total_data.logFragSize %>%
  gather(-logFragSize, key = "var", value = "value") %>%
  ggplot(aes(x = logFragSize, y = value)) +
  geom_point(colour = "#D55E00", size= 2) +
  stat_smooth() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

#scatter plot: AvgIsoDist_parasite rich estimates
total_data.AvgIsoDist %>%
  gather(-AvgIsoDist, key = "var", value = "value") %>%
  ggplot(aes(x = AvgIsoDist, y = value)) +
  geom_point(colour="#D55E00",size = 2) +
  stat_smooth() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

##Fitting linear models

lm.MeanNumTrees <- lm(boot ~ MeanNumTrees, total_data)
summary(lm.MeanNumTrees)#ns p=0.066
ci.R2(R2=0.3, N=12, K=1, conf.level=0.95) #95% CI: [0, 0.7] for effect size 0.2165

#Diagnostics
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(lm.MeanNumTrees, las = 1)      # Residuals, Fitted, ...

lm.MeanHghtTrees <- lm(boot ~ MeanHghtTrees, total_data)
summary(lm.MeanHghtTrees)#ns p=0.91
#Diagnostics
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(lm.MeanHghtTrees, las = 1)      # Residuals, Fitted, ...

lm.MeanPcentCanCov <- lm(boot ~ MeanPcentCanCov, total_data)
summary(lm.MeanPcentCanCov)#ns p=0.44
#Diagnostics
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(lm.MeanPcentCanCov, las = 1)      # Residuals, Fitted, ...

lm.MeanPcentShrubCov <- lm(boot ~ MeanPcentShrubCov, total_data)
summary(lm.MeanPcentShrubCov)#ns p=0.044 Adjusted R-squared (effect size):  0.28 
ci.R2(R2=0.28, N=12, K=1, conf.level=0.95) #95% CI: [0, 0.68] for effect size 0.28

#Diagnostics
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(lm.MeanPcentShrubCov, las = 1)      # Residuals, Fitted, ...
gvlma(lm.MeanPcentShrubCov)

lm.logFragSize <- lm(boot ~ logFragSize, total_data)
summary(lm.logFragSize)#ns p=0.38
ci.R2(R2=0.0456, N=19, K=1) #95% CI: [0, 0.35] for effect size 0.028
#Diagnostics
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(lm.logFragSize, las = 1)      # Residuals, Fitted, ...

####Categorical variables
#Finally, using boxplots, checking relationships 
#between response variables and categorical predictor variables

#One way ANOVA to test significance
#http://www.sthda.com/english/wiki/one-way-anova-test-in-r

par(mar=c(3,4,3,2), mfrow=c(2,2))
colors = terrain.colors(6)[5:1]
plot(chao~road, total_data,pch = 16, cex = 1, 
     col = "gray", ylab = "Chao")
plot(jack1~road, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "jack-knife1")
plot(jack2~road, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "jack-knife2")
plot(boot~road, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "bootstrap")
title( "Road", outer = T, line=-1)#Overall title
#One-way ANOVA
road.chao_anova <- aov(chao~road, data = total_data)
road.jack1_anova <- aov(jack1~road, data = total_data)
road.jack2_anova <- aov(jack2~road, data = total_data)
road.boot_anova <- aov(boot~road, data = total_data)
# Summary of the analysis
summary(road.chao_anova) #not significant p=0.73
summary(road.jack1_anova) #not significant p=0.81
summary(road.jack2_anova) #not significant p=0.91
summary(road.boot_anova) #not significant p=0.68

par(mar=c(3,4,3,2), mfrow=c(2,2))
colors = terrain.colors(6)[5:1]
plot(chao~plantation, total_data,pch = 16, cex = 1, 
     col = "gray", ylab = "Chao")
plot(jack1~plantation, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "jack-knife1")
plot(jack2~plantation, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "jack-knife2")
plot(boot~plantation, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "bootstrap")
title( "Plantation", outer = T, line=-1)#Overall title
#One-way ANOVA
plantation.chao_anova <- aov(chao~plantation, data = total_data)
plantation.jack1_anova <- aov(jack1~plantation, data = total_data)
plantation.jack2_anova <- aov(jack2~plantation, data = total_data)
plantation.boot_anova <- aov(boot~plantation, data = total_data)
# Summary of the analysis
summary(plantation.chao_anova) #not significant p=0.16
summary(plantation.jack1_anova) #not significant p=0.057
summary(plantation.jack2_anova) #not significant p=0.13
summary(plantation.boot_anova) #significant p=0.036
eta_sq(plantation.boot_anova,partial=TRUE, ci.lvl = 0.9) #partial Eta-squared Effect size=0.215
#95% CI: [0, 0.5] for effect size 0.23
#90% CI: [0.009, 0.46] for effect size 0.23

par(mar=c(3,4,3,2), mfrow=c(2,2))
colors = terrain.colors(6)[5:1]
plot(chao~cattle, total_data,pch = 16, cex = 1, 
     col = "gray", ylab = "Chao")
plot(jack1~cattle, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "jack-knife1")
plot(jack2~cattle, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "jack-knife2")
plot(boot~cattle, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "bootstrap")
title( "Cattle", outer = T, line=-1)#Overall title
#One-way ANOVA
cattle.chao_anova <- aov(chao~cattle, data = total_data)
cattle.jack1_anova <- aov(jack1~cattle, data = total_data)
cattle.jack2_anova <- aov(jack2~cattle, data = total_data)
cattle.boot_anova <- aov(boot~cattle, data = total_data)
# Summary of the analysis
summary(cattle.chao_anova) #not significant p=0.763
summary(cattle.jack1_anova) #not significant p=0.862
summary(cattle.jack2_anova) #not significant p=0.979
summary(cattle.boot_anova) #not significant p=0.726

par(mar=c(3,4,3,2), mfrow=c(2,2))
colors = terrain.colors(6)[5:1]
plot(chao~settlement, total_data,pch = 16, cex = 1, 
     col = "gray", ylab = "Chao")
plot(jack1~settlement, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "jack-knife1")
plot(jack2~settlement, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "jack-knife2")
plot(boot~settlement, total_data, pch = 16, cex = 1, 
     col = "gray",ylab = "bootstrap")
title( "Human settlement", outer = T, line=-1)#Overall title
#One-way ANOVA
settlement.chao_anova <- aov(chao~settlement, data = total_data)
settlement.jack1_anova <- aov(jack1~settlement, data = total_data)
settlement.jack2_anova <- aov(jack2~settlement, data = total_data)
settlement.boot_anova <- aov(boot~settlement, data = total_data)
# Summary of the analysis
summary(settlement.chao_anova) #not significant p=0.548
summary(settlement.jack1_anova) #not significant p=0.378
summary(settlement.jack2_anova) #not significant p=0.654
summary(settlement.boot_anova) #not significant p=0.219


######Based on HostRange################################################
#######################################################################
#Sampling effort(sample size) may influence observed richness values 
#so extrapolated richness were estimated (using specpool) 
#to account for rare parasites in each fragment
fspool_HRange<-with(ana_data, specpool(abun.mat, HostRange))

HRange<-fspool_HRange[,c(7,9)]
HRange<-HRange$HostRange
hr<-row.names(fspool_HRange)
HRange$HostRange<-hr

lm.HRange<-lm(boot~HostRange, HRange)
summary(lm.HRange)


ana_data.HR <- select(ana_data, FragName, HostRange, HostDiet)
head(ana_data.HR)

large<-filter(ana_data.HR, HostRange=="large")
med<-filter(ana_data.HR, HostRange=="med")
small<-filter(ana_data.HR, HostRange=="small")

levels(large$FragName)
levels(med$FragName)
levels(small$FragName)



##############################################################################################################################

##############################################################################################################################
# Calculating Richness with equal sample sizes (=Smallest sample size)
##############################################################################################################################

#Creating new data frame 16 random samples from each frag size pool
rr<-ana_data %>% group_by(FragName) %>% sample_n(size = 16, replace = F)

#Creating new abundance matrix with equal sample sizes
rr.mat <- rr[,5:45]

fspool.rr<-with(rr, specpool(rr.mat, FragName))
fspool.rr <-na.omit(fspool.rr) #Removing data with NA

total_data$rr.Species<-fspool.rr$Species
total_data$rr.chao<-fspool.rr$chao
total_data$rr.jack1<-fspool.rr$jack1
total_data$rr.jack2 <- fspool.rr$jack2
total_data$rr.boot <- fspool.rr$boot

cols.rrvars <- c("rr.species", "rr.chao", "rr.jack1", "rr.jack2", "rr.boot")

total_data %>%
  select(!!cols.rrvars) %>%               # select only useful, numerical variables
  gather() %>%                           # Convert to key-value pairs
  ggplot(aes(value)) +                   # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_density()+                         # as density
  theme_minimal()
#####Distribution of all response variables are close to normality




######Interaction between host richness and disturbance
######################################################
lmm.0 <- lmer(log(par.rich) ~ host.rich  + (1 | host), data = total_data,
              REML = FALSE)
summary(lmm.0)  #AIC 519.2
Anova(lmm.0)

lmm.1 <- lmer(log(par.rich) ~ host.rich*plantation  +  (1 | host), data = total_data,
              REML = FALSE)
summary(lmm.1) #514.3
Anova(lmm.1)

lmm.2 <- lmer(log(par.rich) ~ host.rich*plantation*cattle+ (1 | host), data = total_data,
              REML = FALSE)
summary(lmm.2) #516.6
Anova(lmm.2)

lmm.3 <- lmer(log(par.rich) ~ host.rich*plantation*cattle*settlement  + (1 | host), data = total_data,
              REML = FALSE)
summary(lmm.3) #520.5
Anova(lmm.3)

lmm.4 <- lmer(log(par.rich) ~ host.rich*plantation*cattle+FragSize +(1 |host), data = total_data,
              REML = FALSE)
summary(lmm.4) #513.4
Anova(lmm.4)
confint(lmm.4, oldNames=FALSE)


ggplot(total_data,aes(x=host.rich, y=par.rich, group = interaction(plantation, cattle))) + 
  geom_smooth(method="lm", alpha = 0.1, aes(lty=plantation, color=cattle))


lmm.5 <- lmer(log(par.rich) ~ host.rich*plantation*cattle+FragSize+AvgIsoDist  + (1 | host), data = total_data,
              REML = FALSE)
summary(lmm.5) #514.7
Anova(lmm.5)

lmm.6 <- lmer(log(par.rich) ~ host.rich*cattle  + (1 | host), data = total_data,
              REML = FALSE)
summary(lmm.6) #520.6
Anova(lmm.6)

lmm.7 <- lmer(par.rich ~ host.rich*settlement  + (1 | host), data = total_data,
              REML = FALSE)
summary(lmm.7) #1264
Anova(lmm.7)
confint(lmm.7, oldNames=FALSE)

lmm.4.1 <- lmer(par.rich ~ (host.rich+plantation+cattle+settlement)^2+FragSize+(1 | host), data = total_data,
                REML = FALSE)
summary(lmm.4.1) #1353





##Linear Mixed Effects Model for interaction between hostdi and disturbance ####
###############################################################################

###To compare richness between fished and unfished islands
#we ran a mixed-effects linear model with disturbance status (plantation, livestock, settlement) 
#as fixed effects and fragment and host species as random effects,where replicates were 
#bootstrap-estimated parasite taxon diversity for each host-island combination.
#We included the covariates FragSize, Isolation index to ensure that these factors were not 
#driving patterns in parasite diversity. Covariates were removed from the model through 
#backwards elimination if they were not significant at alpha=0.05. 

mod1<- lmer(par.rich~cattle*settlement+FragSize+AvgIsoDist+
              (1|host)+(1|frag), data=total_data, REML=FALSE)
summary(mod1) #AIC 1351.2
Anova(mod1)

mod2<- lmer(par.rich~cattle*settlement+FragSize+(1|host)+(1|frag), data=total_data, REML=FALSE)
summary(mod2) #AIC 1349.2
Anova(mod2)
confint(mod2, oldNames=FALSE)



##############Host diversity interaction##########################
##################################################################

full.model.para.plant <- lm(boot.para~plantation*boot.host, total_data)

summary(full.model.para)
confint(full.model.para)
#Plotting of confidence intervals of coefficients
y = coef(full.model.para.plant)
x = seq_along(y)
ci = confint(full.model.para.plant)
xlim = range(x) + c(-0.5,0.2)
ylim = range(ci)
ylim = ylim + 0.1*c(-1,+1)*diff(ylim) # extend it a little
ylab = bquote(hat(beta))
xlab = "coefficient"
par(mar=c(4.5,12,1,1), las=1)
plot(y, pch=16, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", bty="n")
axis(1, at=x, labels=names(y), tick=FALSE)
abline(h=0, lty=3)
arrows(x,ci[,1],x,ci[,2], code=3, angle=90, length=0.05)
title(main = "Para.Richness.plantation")
### horizontal layout:
plot(x=y, y=x, pch=16, xlim=ylim, ylim=xlim, xlab=ylab, ylab="", yaxt="n", bty="n")
axis(2, at=x, labels=names(y), tick=FALSE)
abline(v=0, lty=3)
arrows(ci[,1],x,ci[,2],x, code=3, angle=90, length=0.05)
title(main = "Para.Richness.plantation")

full.model.para.settlement <- lm(boot.para~settlement*boot.host, total_data)

summary(full.model.para.settlement)
#Plotting of confidence intervals of coefficients
y = coef(full.model.para.settlement)
x = seq_along(y)
ci = confint(full.model.para.settlement)
xlim = range(x) + c(-0.5,0.2)
ylim = range(ci)
ylim = ylim + 0.1*c(-1,+1)*diff(ylim) # extend it a little
ylab = bquote(hat(beta))
xlab = "coefficient"
par(mar=c(4.5,12,1,1), las=1)
plot(y, pch=16, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", bty="n")
axis(1, at=x, labels=names(y), tick=FALSE)
abline(h=0, lty=3)
arrows(x,ci[,1],x,ci[,2], code=3, angle=90, length=0.05)
title(main = "Para.Richness.settlement")
### horizontal layout:
plot(x=y, y=x, pch=16, xlim=ylim, ylim=xlim, xlab=ylab, ylab="", yaxt="n", bty="n")
axis(2, at=x, labels=names(y), tick=FALSE)
abline(v=0, lty=3)
arrows(ci[,1],x,ci[,2],x, code=3, angle=90, length=0.05)
title(main = "Para.Richness.settlement")


full.model.para.cattle <- lm(boot.para~cattle*boot.host, total_data)
summary(full.model.para.cattle)
#Plotting of confidence intervals of coefficients
y = coef(full.model.para.cattle)
x = seq_along(y)
ci = confint(full.model.para.cattle)
xlim = range(x) + c(-0.5,0.2)
ylim = range(ci)
ylim = ylim + 0.1*c(-1,+1)*diff(ylim) # extend it a little
ylab = bquote(hat(beta))
xlab = "coefficient"
par(mar=c(4.5,12,1,1), las=1)
plot(y, pch=16, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n", bty="n")
axis(1, at=x, labels=names(y), tick=FALSE)
abline(h=0, lty=3)
arrows(x,ci[,1],x,ci[,2], code=3, angle=90, length=0.05)
title(main = "Para.Richness.cattle")
### horizontal layout:
plot(x=y, y=x, pch=16, xlim=ylim, ylim=xlim, xlab=ylab, ylab="", yaxt="n", bty="n")
axis(2, at=x, labels=names(y), tick=FALSE)
abline(v=0, lty=3)
arrows(ci[,1],x,ci[,2],x, code=3, angle=90, length=0.05)
title(main = "Para.Richness.cattle")


#################################################################
#################################################################
#exploratory data analysys
#https://www.r-bloggers.com/exploratory-data-analysis-in-r-introduction/
basic_eda <- function(total_data)
{
  glimpse(total_data)
  df_status(total_data)
  freq(total_data) 
  profiling_num(total_data)
  plot_num(total_data)
  describe(total_data)
}

basic_eda(total_data)


#Testing normal dist
qqp(dat$par.rich, "norm")
qqp(dat$host.rich, "norm")

# lnorm means lognormal
qqp(dat$par.rich, "lnorm")
qqp(dat$host.rich, "lnorm")

# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. One can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.
nbinom <- fitdistr(dat$par.rich, "Negative Binomial") ###nbinom fitted
qqp(dat$par.rich, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

nbinom <- fitdistr(dat$host.rich, "Negative Binomial") ###nbinom fitted
qqp(dat$par.rich, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])


poisson <- fitdistr(dat$par.rich, "Poisson")
qqp(dat$host.rich, "pois", lambda = poisson$estimate)#####Poisson selected

gamma <- fitdistr(dat$host.rich, "gamma")
qqp(dat$host.rich, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

summary(dat$par.rich)


###Linear models
mod0<-lm(log(par.rich)~1, data=dat)
add1(mod0,scope = dat, test ="F")

m1<-lm(log(par.rich)~host.rich, data=dat)
summary(m1)
add1(m1,scope = dat, test ="F")

m2<-lm(log(par.rich)~host.rich+plantation, data=dat)
summary(m2)
add1(m2,scope = dat, test ="F")

#####To test difference between fragments with different tratments
######Since treatment groups are colinear, specificlly settlement presence, 
#each treatment will be analysed separately
###########ANOVA for host richness
##################################
vif(lm(host.rich~plantation+cattle+settlement, data = total_data)) #Testing multicolinearity

host1 <- lm(host.rich~plantation, data = total_data)
# Summary of the analysis
summary(host1) 
anova(host1)
confint(host1, conf.level = 0.95)

host2 <- lm(host.rich~cattle, data = total_data)
# Summary of the analysis
summary(host2) 
anova(host2)
confint(host2, conf.level = 0.95)

host3 <- lm(host.rich~settlement, data = total_data)
# Summary of the analysis
summary(host3) 
anova(host3)
confint(host3, conf.level = 0.95)

####One-way ANOVA for parasite richness
para_aov <- aov(par.rich~plantation*cattle*settlement, data = total_data)
# Summary of the analysis
summary(para_aov) #p<0.05
TukeyHSD(para_aov, conf.level = 0.95)


FragType.para_aov <- aov(par.rich~FragType, data = total_data)
summary(FragType.para_aov) #p>0.05
TukeyHSD(FragType.para_aov, conf.level = 0.95)


##Effect size_aov
#plot(eta_sq(plant.host_aov,partial=TRUE, ci.lvl = 0.95)) #partial Eta-squared Effect size=0.215
#95% CI: [0, 0.5] for effect size 0.23
#90% CI: [0.009, 0.46] for effect size 0.23