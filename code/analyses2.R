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
sp <- d[,27:45]


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
sp <- d[,27:45]
# focus on those groups that are related to litter decomposition and nutrient cycling processes
sp2<- sp[,c(1:10, 14:16, 18, 19)]

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
                                var.always.include=c(3,7,9)) 
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
env <- d[,c(8:10,12:26)]
env[,11:15] <- env[,10]* (env[,10:14]/100)

#Log the biomass data because differ an order of magnitude. 
env[,c(9:12, 15:18)] <- log(env[,c(9:12, 15:18)])
env[7,14] <- log(env[7,14]) # do the log as this is the only value in this column different from zero.

#For analyses eliminate which plot it is, also OLEA because there is not litter in Marrufo forest. 
env <- env[,-12]

#Eliminate wood and non-id litter and others
env <- env[,c(-12, -13, -17)]

#light plot of samples 21 and 22 have much more light that the others. Therefore que arcsin sqrt transform them to avoid overdispersion in variance
env[,1]<-asin(sqrt(env[,1]/100))

env <- env[-5, ] # to remove an outlier because the humidity level is to high


#obtain species (sp) separated from environmental variables (env)
sp <- d[,27:45]
# focus on those groups that are related to litter decomposition and nutrient cycling processes
sp2<- sp[,c(1:10, 14:16, 18, 19)]

#merge database between decomposers and predators
decomp <- sp2[, c(1,6:15)]
predators <-sp2[, c(2:5)]

decomp <- sqrt(decomp)
decomp <-wisconsin(decomp)

decomp <- decomp[-5, ] # to remove an outlier #DOUBLE CHECK WHETHER THIS IS TRUE. 


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
                                var.always.include=c(3,5,8)) # these are the groups that were included in the frist random run of bv.step analyses. 
res.bv.step.biobio2 


### Check how these groups relate to the each of the axes of the NDMS
MDS_res=metaMDS(decomp, distance = "bray", k = 2, trymax = 50) #optimal number will be three
bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio2$order.by.best$var.incl[1], ","))) 
bio.fit <- envfit(MDS_res, decomp[,bio.keep], perm = 999) 
bio.fit 



#Do the same with the environmental variables, 
env.fit <- envfit(MDS_res, env[,c(1,5,14)], perm = 999) 
env.fit 

#Plot the graphs
pdf(file="results/NMDS_soil_fauna_decomposers_Fig2.pdf", width=10, height=6) 
par(mfcol=c(1,2), omi=c(0.1, 0.1, 0.1, 0.1), mar=c(4, 4, 1, 1), ps=8) 
#first plot related to soil variables
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1.1, 1.1), ylim = c(-0.6,0.6)) 
plot(env.fit, col="gray10", cex=1.2, font=4,  labels = c("Light availability", "Tree defoliation", "High decomposed litter")) 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]) , max(MDS_res$points[,2]-0.05), paste("Stress =",round(MDS_res$stress, 2)), pos=4, font=3, col="gray30") 

#second plot related to soil fauna groups
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1.1, 1.1), ylim = c(-0.6,0.6)) 
plot(bio.fit, col="gray10", cex=1.2, font=4, labels = c("Oribatida", "Entomobryomorpha", "Poduromorpha",
                                                        "Myriapoda", "Enchytraeidae")) 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]), max(MDS_res$points[,2]-0.05), paste("Stress =",round(MDS_res$stress, 2)), pos=4, font=3, col="gray30") 
dev.off()


###This is for predators, we are going to include abundance of decomposers as a explanatory variable.

decomposers <- sp2[, c(1,6:15)]
sum_decomp <- apply(decomposers, 1, FUN="sum")

##remove outlier #5 
sum_decomp <-sum_decomp[-5]

env2<-cbind(env, sum_decomp)
#env2<- env2[-1, ] # to remove an outlier
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
                                var.always.include=c(2)) 
res.bv.step.biobio2 


### Check how these groups relate to the each of the axes of the NDMS
MDS_res=metaMDS(predators, distance = "bray", k = 2, trymax = 50) 
bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio2$order.by.best$var.incl[1], ","))) 
bio.fit <- envfit(MDS_res, predators[,bio.keep], perm = 999) 
bio.fit 



#Do the same with the environmental variables, 
env.fit <- envfit(MDS_res, env2[,c(5,7,9,13)], perm = 999)  
env.fit 

#Plot the graphs
pdf(file="results/NMDS_soil_fauna_predators_Fig3.pdf", width=10, height=6) 
par(mfcol=c(1,2), omi=c(0.1, 0.1, 0.1, 0.1), mar=c(4, 4, 1, 1), ps=8) 
#first plot related to soil variables
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1.7, 1.7)) 
plot(env.fit, col="gray10", cex=1.2, font=4, labels = c("Tree defoliation", "Litter depth", 
                                                        "Litter mass", "Medium decomposed litter")) 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]), max(MDS_res$points[,2]), paste("Stress =",round(MDS_res$stress, 2)), pos=2, font=3, col="gray30") 

#second plot related to soil fauna groups
plot(MDS_res$points, t="n",xlab="NMDS1", ylab="NMDS2", xlim=c(-1.7, 1.7)) 
plot(bio.fit, col="gray10", cex=1.2, font=4, labels="Prostigmata") 
text(MDS_res$points, as.character(1:length(MDS_res$points[,1])), cex=1.2, col="red") 
text(min(MDS_res$points[,1]), max(MDS_res$points[,2]), paste("Stress =",round(MDS_res$stress, 2)), pos=2, font=3, col="gray30") 
dev.off()


#SEM (Structural equation modelling) analyses---- 

# Load the data
d <- read.csv("data/original_database_r.csv", header=T, sep=",")

# put all units of litter in biomass for the three main species to have same units in (g) as total litter biomass and wood
env <- d[,c(8:10,12:26)]
env[,11:15] <- env[,10]* (env[,10:14]/100)

#Log the biomass data because differ an order of magnitude. 
env[,c(9:12, 15:18)] <- log(env[,c(9:12, 15:18)])
env[7,14] <- log(env[7,14]) # do the log as this is the only value in this column different from zero.

#For analyses eliminate which plot it is, also OLEA because there is not litter in Marrufo forest. 
env <- env[,-12]

#Eliminate wood and non-id litter and others
env <- env[,c(-12, -13, -17)]

#light plot of samples 21 and 22 have much more light that the others. Therefore que arcsin sqrt transform them to avoid overdispersion in variance
env[,1]<-asin(sqrt(env[,1]/100))

env <- env[-5, ] # to remove an outlier because the humidity level is to high


#obtain species (sp) separated from environmental variables (env)
sp <- d[,27:45]
# focus on those groups that are related to litter decomposition and nutrient cycling processes
sp2<- sp[,c(1:10, 14:16, 18, 19)]

sp2 <- sp2[-5, ]

#merge database between decomposers and predators
decomp <- sp2[, c(1,6:15)]
predators <-sp2[, c(2:5)]


#to see direct and indirect effect of soil fauna richness and abundances. 
library(lavaan)
library(qgraph)

#load env and sp2
#calculate richness it does not longer applies
#richness <- diversity(sp2, index = "shannon", MARGIN = 1, base = exp(1))
#richness_decomp <- diversity(decomp, index = "shannon", MARGIN = 1, base = exp(1))
#richness_predators <- diversity(predators, index = "shannon", MARGIN = 1, base = exp(1))

#calculate abundances
abundance_decomp <- apply(decomp, 1, FUN = "sum")
abundance_predators <- apply(predators, 1, FUN = "sum")
abundance_oribatida <- decomp[,1]
abundance_entomobryomorpha <- decomp[,3]
abundance_poduromorpha <- decomp[,5]
abundance_miryapoda <- decomp[,8]
abundance_mesostigmata <-predators[, 1]
abundance_prostigmata <-predators[, 2]
abundance_predator_mites <- apply(predators[,1:2], 1, FUN = "sum")


semd <-cbind(env, abundance_decomp, 
             abundance_predators, abundance_oribatida, abundance_entomobryomorpha, abundance_poduromorpha, abundance_miryapoda, abundance_mesostigmata,
             abundance_prostigmata, abundance_predator_mites)
             

##model fitting says pathogen variance is too high, log transform for solving 
semd$pathogen[17]<-1 # this is to do not have -Inf when transforming
semd$pathogen <- log(semd$pathogen)
semd$light_plot <- 100*semd$light_plot
semd$tree_size <- log(semd$tree_size)
semd$abundance_decomp <- log(semd$abundance_decomp)
semd$abundance_predators <- log(semd$abundance_predators)

# Model that includes the overall effect of the determinants of abundances among decomposers and predators. 

model1 <- '
# regressions

abundance_decomp ~ tree_size
abundance_decomp ~ light_plot
abundance_decomp ~~ pathogen

abundance_predators ~ light_plot
abundance_predators ~ tree_size
abundance_predators ~~ pathogen

# correlation
abundance_predators ~~ abundance_decomp
light_plot ~~ tree_size
pathogen ~~ tree_size


#latent variables

#intercept

'

fit <- sem(model1, data=semd, fixed.x=FALSE)
summary(fit, standardized=TRUE)
#this function is to obtain the fit form the model see http://lavaan.ugent.be/tutorial/tutorial.pdf
fitMeasures(fit, c("cfi","rmsea","srmr", "pvalue"))

corrgram(semd[,1:16], order=TRUE, lower.panel=panel.conf, upper.panel=panel.pts, text.panel=panel.txt,
         diag.panel=panel.density, pch=16, lty=1, main="Environmental correlations")

##represent the results. 
library(qgraph)
library(htmlwidgets)
library(semPlot)

#Standarized plot diagram http://sachaepskamp.com/semPlot/examples for examples
semPaths(fit, what = "std", whatLabels = 'std', style= "ram", layout= "spring", residuals= FALSE, 
         sizeInt =500, edge.label.cex = 1.5, node.width=2, node.height=2, label.cex= 3, curvePivot = TRUE)








#To avoid problems when logging
semd$abundance_oribatida <- log(semd$abundance_oribatida)
semd$abundance_mesostigmata[5] <- 1
abundance_entomobryomorpha <- decomp[,3]
abundance_poduromorpha <- decomp[,5]
abundance_miryapoda <- decomp[,8]
abundance_prostigmata <-predators[, 2]

#model 1 describes the best SEM model for decomposers
model1 <- '
# regressions
abundance_decomp ~ litter_moisture
abundance_decomp ~ light_plot
litter_moisture ~ light_plot
litter_moisture ~ litter_depth
litter_depth ~ tree_size

# residual correlations
#litter_moisture ~~ decomposition_level_3
'

fit <- sem(model1, data=semd, fixed.x=FALSE)
summary(fit, standardized=TRUE)

#model 2 describes the best SEM model for oribatida
model2 <- '
# regressions
abundance_oribatida ~ litter_moisture
abundance_oribatida ~ light_plot
litter_moisture ~ light_plot
litter_moisture ~ litter_depth
litter_depth ~ tree_size

# residual correlations
#litter_moisture ~~ decomposition_level_3
'

fit <- sem(model2, data=semd, fixed.x=FALSE)
summary(fit, standardized=TRUE)

#model 3 describes the best SEM model for poduromorpha
model3 <- '
# regressions
abundance_poduromorpha ~ litter_moisture
abundance_poduromorpha ~ light_plot
litter_moisture ~ light_plot
litter_moisture ~ litter_depth
litter_depth ~ tree_size

# residual correlations
#litter_moisture ~~ decomposition_level_3
'

fit <- sem(model3, data=semd, fixed.x=FALSE)
summary(fit, standardized=TRUE)

#model 4 describes the best SEM model for miryapoda
model4 <- '
# regressions
abundance_miryapoda ~ litter_moisture
abundance_miryapoda ~ light_plot
litter_moisture ~ light_plot
litter_moisture ~ litter_depth
litter_depth ~ tree_size

# residual correlations
#litter_moisture ~~ decomposition_level_3
'

fit <- sem(model4, data=semd, fixed.x=FALSE)
summary(fit, standardized=TRUE)

#model 5 describes the best SEM model for miryapoda
model5 <- '
# regressions
abundance_predators ~ pathogen
abundance_predators~ abundance_decomp
abundance_predators ~ light_plot
abundance_decomp ~ light_plot

#correlation

pathogen~~ tree_size
'

fit <- sem(model5, data=semd, fixed.x=FALSE)
summary(fit, standardized=TRUE)

####HASTA AQUI!!!

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



