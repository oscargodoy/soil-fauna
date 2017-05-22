#Analyze soil fauna and its determinants Alejandro Jimenez TFM
#May 22th 2017

#loading libraries
#in case it is not loaded write the following install.packages('vegan'), change name accordingly to the library you want to download
library(vegan)
library(corrgram) # for correlograms

# Load the data
d <- read.csv("data/original_database_r.csv", header=T, sep=",")

#obtain species (sp) separated from environmental variables (env)
sp <- d[,27:42]

#Following Santonja 2017 JECol square-root transformed and then standardized using a Wisconsin double standardization
sp <- sqrt(sp)
sp <-wisconsin(sp)

# put all units in biomass for the three main species 
env <- d[,c(7:9,11:26)]
env[,11:15] <- env[,10]* (env[,11:15]/100)

#Log the biomass data because differ an order of magnitude. 
env[,c(10:12, 15:19)] <- log(env[,c(10:12, 15:19)])
env[7,14] <- log(env[7,14])

#For analyses eliminate which plot it is, also OLEA because there is not litter
env <- env[,c(-1, -13)]

#IMPORTANTE!!!!!! Eliminate now the variables that we do not have data. 
env <- env[,c(-6, -8)]

##IMPORTANTE PONER VALORES QUE FALTAN DE HUMEDAD

env[c(4,11,18),1]<-0.33 #THIS values are not TRUEEEEEEE

#Before performing analyses check whether some variables are high strongly correlated. 

corrgram(env, order=TRUE, lower.panel=panel.conf, upper.panel=panel.pts, text.panel=panel.txt,
         diag.panel=panel.density, pch=16, lty=1, main="Environmental correlations")


#Explore the data in a muldimensional way to know the main axes of variation
#code obtained from http://menugget.blogspot.com.es/2011/06/clarke-and-ainsworths-bioenv-and-bvstep.html

res<-bioenv(sp, env=env)
res

#First results says that highly decomposed material is a strong predictor. Let's eliminate one of the variables
env <- env[,-11]
res2<-bioenv(sp, env=env)
res2

source('bio.env.R')
res2.2<-bio.env(sp, env, fix.dist.method="bray", var.dist.method="euclidean", scale.fix=FALSE, scale.var=TRUE) 
res2.2
#Second analyses says that litter_decomposition 3 and litter_biomass are significant but they are highly correlated (eliminate litter_mass)
env <- env[,-7]
res3<-bioenv(sp, env=env)
res3
res3.2<-bio.env(sp, env, fix.dist.method="bray", var.dist.method="euclidean", scale.fix=FALSE, scale.var=TRUE) 
res3.2

#lets plot these graphs 
source('bv.step.R')

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
                                var.always.include=c(9,11,12,15)) 
res.bv.step.biobio2 


### Check how these groups relate to the each of the axes of the NDMS
MDS_res=metaMDS(sp, distance = "bray", k = 3, trymax = 50) 
bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio2$order.by.best$var.incl[1], ","))) 
bio.fit <- envfit(MDS_res, sp[,bio.keep], perm = 999) 
bio.fit 



#Do the same with the environmental variables, 
env.fit <- envfit(MDS_res, env[,c(1,2,5,12)], perm = 999) 
env.fit 

#Plot the graphs
pdf(file="NMDS_soil_fauna.pdf", width=6, height=3) 
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

###Do the analyses again grouping insect by trophic groups!!!!Eliminate those groups that do not play a direct role. 

