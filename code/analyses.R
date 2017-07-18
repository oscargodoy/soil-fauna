#Analyze soil fauna and its determinants Alejandro Jimenez TFM
#May 22th 2017

#loading libraries
#in case it is not loaded write the following install.packages('vegan'), change name accordingly to the library you want to download
library(vegan)
library(corrgram) # for correlograms

# NMDS----
# Load the data
d <- read.csv("data/original_database_r.csv", header=T, sep=",")

#obtain species (sp) separated from environmental variables (env)
sp <- d[,27:42]


#Following Santonja 2017 JECol square-root transformed and then standardized using a Wisconsin double standardization
sp <- sqrt(sp)
sp <-wisconsin(sp)

# put all units of litter in biomass for the three main species to have same units in (g) as total litter biomass and wood
env <- d[,c(7:9,11:26)]
env[,11:15] <- env[,10]* (env[,11:15]/100)

#Log the biomass data because differ an order of magnitude. 
env[,c(10:12, 15:19)] <- log(env[,c(10:12, 15:19)])
env[7,14] <- log(env[7,14])

#For analyses eliminate which plot it is, also OLEA because there is not litter in Marrufo forest. 
env <- env[,c(-1, -13)]

#Before performing analyses check whether some variables are high strongly correlated. 

corrgram(env, order=TRUE, lower.panel=panel.conf, upper.panel=panel.pts, text.panel=panel.txt,
         diag.panel=panel.density, pch=16, lty=1, main="Environmental correlations")


#Explore the data in a muldimensional way to know the main axes of variation
#code obtained from http://menugget.blogspot.com.es/2011/06/clarke-and-ainsworths-bioenv-and-bvstep.html

res<-bioenv(sp, env=env)
res

#Eliminate wood and non-id litter
env <- env[,c(-13, -17)]
res2<-bioenv(sp, env=env)
res2

source('code/bio.env.R')
res2.2<-bio.env(sp, env, fix.dist.method="bray", var.dist.method="euclidean", scale.fix=FALSE, scale.var=TRUE) 
res2.2

#lets plot these graphs 
source('code/bv.step.R')

#This analysis serve to know which gropus are best exaplining patterns
res.bv.step.biobio <- bv.step(sp, sp,  
                              fix.dist.method="bray", var.dist.method="bray", 
                              scale.fix=FALSE, scale.var=FALSE,  
                              max.rho=0.95, min.delta.rho=0.001, 
                              random.selection=TRUE, 
                              prop.selected.var=0.3, 
                              num.restarts=50, 
                              output.best=10, 
                              var.always.include=NULL)
res.bv.step.biobio

#Do the same analysis including those groups important
res.bv.step.biobio2  <- bv.step(sp, sp,  
                                fix.dist.method="bray", var.dist.method="bray", 
                                scale.fix=FALSE, scale.var=FALSE,  
                                max.rho=0.95, min.delta.rho=0.001, 
                                random.selection=TRUE, 
                                prop.selected.var=0.2, 
                                num.restarts=50, 
                                output.best=10, 
                                var.always.include=c(9,12,14,15)) 
res.bv.step.biobio2 


### Check how these groups relate to the each of the axes of the NDMS
MDS_res=metaMDS(sp, distance = "bray", k = 3, trymax = 50) 
bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio2$order.by.best$var.incl[1], ","))) 
bio.fit <- envfit(MDS_res, sp[,bio.keep], perm = 999) 
bio.fit 



#Do the same with the environmental variables, 
env.fit <- envfit(MDS_res, env[,c(1,5,15)], perm = 999) 
env.fit 

#Plot the graphs
pdf(file="results/NMDS_soil_fauna_AppendixS1.pdf", width=10, height=6) 
par(mfcol=c(1,2), omi=c(0.1, 0.1, 0.1, 0.1), mar=c(4, 4, 1, 1), ps=8) 
#first plot related to soil variables
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1, 1)) 
plot(env.fit, col="gray10", cex=1.2, font=4) 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]), max(MDS_res$points[,2]), paste("Stress =",round(MDS_res$stress, 2)), pos=4, font=3, col="gray30") 

#second plot related to soil fauna groups
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1, 1)) 
plot(bio.fit, col="gray10", cex=1.2, font=4) 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]), max(MDS_res$points[,2]), paste("Stress =",round(MDS_res$stress, 2)), pos=4, font=3, col="gray30") 
dev.off()





####2
# Load the data
d <- read.csv("data/original_database_r.csv", header=T, sep=",")

# put all units of litter in biomass for the three main species to have same units in (g) as total litter biomass and wood
env <- d[,c(7:9,11:26)]
env[,11:15] <- env[,10]* (env[,11:15]/100)

#Log the biomass data because differ an order of magnitude. 
env[,c(10:12, 15:19)] <- log(env[,c(10:12, 15:19)])
env[7,14] <- log(env[7,14])

#For analyses eliminate which plot it is, also OLEA because there is not litter in Marrufo forest. 
env <- env[,c(-1, -13)]

#Eliminate wood and non-id litter and others
env <- env[,c(-12, -13, -17)]

#light plot of samples 21 and 22 have much more light that the others. Therefore que arcsin sqrt transform them to avoid overdispersion in variance
env[,1]<-asin(sqrt(env[,1]/100))
##CHECK HOW TO TRANSFORM THE DATA TO AVOID OVERDISPERSION!!!


#obtain species (sp) separated from environmental variables (env)
sp <- d[,27:42]
# focus on those groups that are related to litter decomposition and nutrient cycling processes
sp2<- sp[,c(1:10, 14:16)]

sp2 <- sqrt(sp2)
sp2 <-wisconsin(sp2)

###Do the analyses again grouping insects by groups directly linked to litter. !Eliminate those groups that do not play a direct role. 

source('code/bio.env.R')
res2.3<-bio.env(sp2, env, fix.dist.method="bray", var.dist.method="euclidean", scale.fix=FALSE, scale.var=TRUE) 
res2.3
#Second analyses says that litter_decomposition 3 and litter_biomass are significant but they are highly correlated (eliminate litter_mass)

#This analysis serve to know which gropus are best exaplining patterns
source('code/bv.step.R')
res.bv.step.biobio <- bv.step(sp2, sp2,  
                              fix.dist.method="bray", var.dist.method="bray", 
                              scale.fix=FALSE, scale.var=FALSE,  
                              max.rho=0.95, min.delta.rho=0.001, 
                              random.selection=TRUE, 
                              prop.selected.var=0.3, 
                              num.restarts=50, 
                              output.best=10, 
                              var.always.include=NULL)
res.bv.step.biobio

#Do the same analysis including those groups important

res.bv.step.biobio2  <- bv.step(sp2, sp2,  
                                fix.dist.method="bray", var.dist.method="bray", 
                                scale.fix=FALSE, scale.var=FALSE,  
                                max.rho=0.95, min.delta.rho=0.001, 
                                random.selection=TRUE, 
                                prop.selected.var=0.2, 
                                num.restarts=50, 
                                output.best=10, 
                                var.always.include=c(1,9,10)) 
res.bv.step.biobio2 


### Check how these groups relate to the each of the axes of the NDMS
MDS_res=metaMDS(sp2, distance = "bray", k = 2, trymax = 50) 
bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio2$order.by.best$var.incl[1], ","))) 
bio.fit <- envfit(MDS_res, sp2[,bio.keep], perm = 999) 
bio.fit 

#Do the same with the environmental variables, 
env.fit <- envfit(MDS_res, env[,1], perm = 999) 
env.fit 

#Plot the graphs
pdf(file="results/NMDS_soil_fauna_specific_groups_Fig1.pdf", width=10, height=6) 
par(mfcol=c(1,2), omi=c(0.1, 0.1, 0.1, 0.1), mar=c(4, 4, 1, 1), ps=8) 
#first plot related to soil variables
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1, 1)) 
plot(env.fit, col="gray10", cex=1.2, font=4) 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]), max(MDS_res$points[,2]), paste("Stress =",round(MDS_res$stress, 2)), pos=2, font=3, col="gray30") 
text(x=0.05, y=-0.8, labels = "light_plot", col="gray10", cex=1.2, font=4)

#second plot related to soil fauna groups
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1, 1)) 
plot(bio.fit, col="gray10", cex=1.2, font=4) 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]), max(MDS_res$points[,2]), paste("Stress =",round(MDS_res$stress, 2)), pos=2, font=3, col="gray30") 
dev.off()




####3 decomposer and predators

# Load the data
d <- read.csv("data/original_database_r.csv", header=T, sep=",")

# put all units of litter in biomass for the three main species to have same units in (g) as total litter biomass and wood
env <- d[,c(7:9,11:26)]
env[,11:15] <- env[,10]* (env[,11:15]/100)

#Log the biomass data because differ an order of magnitude. 
env[,c(10:12, 15:19)] <- log(env[,c(10:12, 15:19)])
env[7,14] <- log(env[7,14])

#For analyses eliminate which plot it is, also OLEA because there is not litter in Marrufo forest. 
env <- env[,c(-1, -13)]

#Eliminate wood and non-id litter and others
env <- env[,c(-12, -13, -17)]

#light plot of samples 21 and 22 have much more light that the others. Therefore que arcsin sqrt transform them to avoid overdispersion in variance
env[,1]<-asin(sqrt(env[,1]/100))

env <- env[-5, ] # to remove an outlier


#obtain species (sp) separated from environmental variables (env)
sp <- d[,27:42]
# focus on those groups that are related to litter decomposition and nutrient cycling processes
sp2<- sp[,c(1:10, 14:16)]

#merge database between decomposers and predators
decomp <- sp2[, c(1,6:13)]
predators <-sp2[, c(2:5)]

decomp <- sqrt(decomp)
decomp <-wisconsin(decomp)

decomp <- decomp[-5, ] # to remove an outlier

predators <- sqrt(predators)
predators <-wisconsin(predators)

predators <- predators[-5, ] # to remove an outlier

#decomp_sum <- as.vector(apply(decomp,1,FUN=sum))
#predators_sum <- as.vector(apply(predators,1,FUN=sum))

#light plot of samples 21 and 22 have much more light that the others. Therefore que arcsin sqrt transform them to avoid overdispersion in variance


###Do the analyses again grouping insects by groups directly linked to litter. !Eliminate those groups that do not play a direct role. 

source('code/bio.env.R')
res2.4<-bio.env(decomp, env, fix.dist.method="bray", var.dist.method="euclidean", scale.fix=FALSE, scale.var=TRUE) 
res2.4
#Second analyses says that litter_decomposition 3 and litter_biomass are significant but they are highly correlated (eliminate litter_mass)

#This analysis serve to know which gropus are best exaplining patterns
source('code/bv.step.R')
res.bv.step.biobio <- bv.step(decomp, decomp,  
                              fix.dist.method="bray", var.dist.method="bray", 
                              scale.fix=FALSE, scale.var=FALSE,  
                              max.rho=0.95, min.delta.rho=0.001, 
                              random.selection=TRUE, 
                              prop.selected.var=0.3, 
                              num.restarts=50, 
                              output.best=10, 
                              var.always.include=NULL)
res.bv.step.biobio

#Do the same analysis including those groups important

res.bv.step.biobio2  <- bv.step(decomp, decomp,  
                                fix.dist.method="bray", var.dist.method="bray", 
                                scale.fix=FALSE, scale.var=FALSE,  
                                max.rho=0.95, min.delta.rho=0.001, 
                                random.selection=TRUE, 
                                prop.selected.var=0.2, 
                                num.restarts=50, 
                                output.best=10, 
                                var.always.include=c(5,8)) 
res.bv.step.biobio2 


### Check how these groups relate to the each of the axes of the NDMS
MDS_res=metaMDS(decomp, distance = "bray", k = 2, trymax = 50) #optimal number will be three
bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio2$order.by.best$var.incl[1], ","))) 
bio.fit <- envfit(MDS_res, decomp[,bio.keep], perm = 999) 
bio.fit 



#Do the same with the environmental variables, 
env.fit <- envfit(MDS_res, env[,c(1,5,8)], perm = 999) 
env.fit 

#Plot the graphs
pdf(file="results/NMDS_soil_fauna_decomposers_Fig2.pdf", width=10, height=6) 
par(mfcol=c(1,2), omi=c(0.1, 0.1, 0.1, 0.1), mar=c(4, 4, 1, 1), ps=8) 
#first plot related to soil variables
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1, 1), ylim = c(-1,0.5)) 
plot(env.fit, col="gray10", cex=1.2, font=4) 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]), max(MDS_res$points[,2]), paste("Stress =",round(MDS_res$stress, 2)), pos=2, font=3, col="gray30") 

#second plot related to soil fauna groups
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1, 1), ylim = c(-1,0.5)) 
plot(bio.fit, col="gray10", cex=1.2, font=4) 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]), max(MDS_res$points[,2]), paste("Stress =",round(MDS_res$stress, 2)), pos=2, font=3, col="gray30") 
dev.off()


###This is for predators, we are going to include abundance of decomposers as a explanatory variable.

decomposers <- sp2[, c(1,6:13)]
sum_decomp <- apply(decomposers, 1, FUN="sum")

##remove outlier #5 
sum_decomp <-sum_decomp[-5]

env2<-cbind(env, sum_decomp)

###Do the analyses again grouping insects by groups directly linked to litter. !Eliminate those groups that do not play a direct role. 

source('code/bio.env.R')
res2.5<-bio.env(predators, env2, fix.dist.method="bray", var.dist.method="euclidean", scale.fix=FALSE, scale.var=TRUE) 
res2.5
#Second analyses says that litter_decomposition 3 and litter_biomass are significant but they are highly correlated (eliminate litter_mass)

#This analysis serve to know which gropus are best exaplining patterns
source('code/bv.step.R')
res.bv.step.biobio <- bv.step(predators, predators,  
                              fix.dist.method="bray", var.dist.method="bray", 
                              scale.fix=FALSE, scale.var=FALSE,  
                              max.rho=0.95, min.delta.rho=0.001, 
                              random.selection=TRUE, 
                              prop.selected.var=0.3, 
                              num.restarts=50, 
                              output.best=10, 
                              var.always.include=NULL)
res.bv.step.biobio

#Do the same analysis including those groups important

res.bv.step.biobio2  <- bv.step(predators, predators,  
                                fix.dist.method="bray", var.dist.method="bray", 
                                scale.fix=FALSE, scale.var=FALSE,  
                                max.rho=0.95, min.delta.rho=0.001, 
                                random.selection=TRUE, 
                                prop.selected.var=0.2, 
                                num.restarts=50, 
                                output.best=10, 
                                var.always.include=c(1)) 
res.bv.step.biobio2 


### Check how these groups relate to the each of the axes of the NDMS
MDS_res=metaMDS(predators, distance = "bray", k = 2, trymax = 50) 
bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio2$order.by.best$var.incl[1], ","))) 
bio.fit <- envfit(MDS_res, decomp[,bio.keep], perm = 999) 
bio.fit 



#Do the same with the environmental variables, 
env.fit <- envfit(MDS_res, env2[,c(1, 6, 11)], perm = 999)  
env.fit 

#Plot the graphs
pdf(file="results/NMDS_soil_fauna_predators_Fig3.pdf", width=10, height=6) 
par(mfcol=c(1,2), omi=c(0.1, 0.1, 0.1, 0.1), mar=c(4, 4, 1, 1), ps=8) 
#first plot related to soil variables
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1.7, 1.7)) 
plot(env.fit, col="gray10", cex=1.2, font=4) 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]), max(MDS_res$points[,2]), paste("Stress =",round(MDS_res$stress, 2)), pos=1, font=3, col="gray30") 

#second plot related to soil fauna groups
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1.7, 1.7)) 
plot(bio.fit, col="gray10", cex=1.2, font=4) 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]), max(MDS_res$points[,2]), paste("Stress =",round(MDS_res$stress, 2)), pos=1, font=3, col="gray30") 
text(x = 0.2, y = 0.50, labels = "mesostigmata", col="gray10", cex=1.2, font=4)
dev.off()


#SEM analyses---- (Structural equation modelling)
#to see direct and indirect effect of soil fauna richness and abundances. 
library(lavaan)
library(qgraph)











##These are glm analyses to check better the broad scale varaition in soil fauna abundance


##Other analyses glm to know the ffect of climatic treatment of abundances of the different groups. 
# Load the data
d <- read.csv("data/original_database_r.csv", header=T, sep=",")
# put all units of litter in biomass for the three main species to have same units in (g) as total litter biomass and wood
env <- d[,c(5,8,9,12:17,25)]
env$litter_mass <- log(env$litter_mass)

#light plot of samples 21 and 22 have much more light that the others. Therefore que arcsin sqrt transform them to avoid overdispersion in variance
env$light_plot<-asin(sqrt(env$light_plot/100))

env_glm <- env[-5,] #remove this site with no predators

#obtain species (sp) separated from environmental variables (env)
sp <- d[,27:42]
# focus on those groups that are related to litter decomposition and nutrient cycling processes
sp2<- sp[,c(1:10, 14:16)]

#merge database between decomposers and predators
decomp <- sp2[, c(1,6:13)]
predators <-sp2[, c(2:5)]

decomp <- decomp[-5, ] # to remove an outlier
predators <- predators[-5, ] # to remove an outlier

decomp_glm <- as.vector(apply(decomp,1,FUN=sum))
predators_glm <- as.vector(apply(predators,1,FUN=sum))

ratio_decomp_predators <- decomp_glm / predators_glm

#First know which type of data distribution we have 
par(mfrow=c(1,2))
hist(decomp_glm, main="Decomposers", xlab= "N. individuals")
hist(predators_glm, main="Predators", xlab= "N. individuals")

# model selection by AIC using glm with poisson distribution. Akaike information criterion (AIC)
#First only main predictors 

corrgram(env_glm, order=TRUE, lower.panel=panel.conf, upper.panel=panel.pts, text.panel=panel.txt,
         diag.panel=panel.density, pch=16, lty=1, main="Environmental correlations")

#litter mass high correlated with decomposition level3 eliminate the latter
#tree size high correlated with litter depth, eliminate tree size. 
#tree defoliation and light plot high correlated, eliminate tree defoliation

env_glm <- env_glm[,c(-4,-10)]

#merge databases
decomp <- cbind(decomp_glm, env_glm)
predators <- cbind(predators_glm, env_glm)
decomp_predators <- cbind(ratio_decomp_predators, env_glm)


#NEED TO DO MODEL SELECTION AS A ROUTINE TO SEARCH FOR THE BEST MODEL USING AIC. 

fit1 <- glm(decomp_glm ~ ., data=decomp, family=poisson((link = "log")))
summary(fit1)

fit2 <- glm(predators_glm ~ ., data=predators, family=poisson((link = "log")))
summary(fit2)

fit3 <- glm(ratio_decomp_predators ~ ., data=decomp_predators, family=poisson((link = "log"))) ##BUSCAR MEJOR FUNCION DE UNION. 
summary(fit3)


##PATH ANALYSIS. 



