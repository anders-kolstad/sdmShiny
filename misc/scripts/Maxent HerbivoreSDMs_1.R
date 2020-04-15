#Maxent herbivore SDM models
rm(list=ls())
require(raster)
require(dismo)
require(rasterVis)
require(ENMeval)

predvars<-stack("~/SDMs/InputData/predvars1k")
soilph<-raster('~/SDMs/InputData/geonode_phihox_m_sl2_250m.tif') 
soilphno<-mask(crop(projectRaster(soilph,predvars[[1]]),predvars[[1]]),predvars[[1]])/10


#Sample randomly 5000 points across Norway
#s1<-sampleRandom(bioclimNor,5000)
herbivoredensSp<-subset(herbivoredens,1:9) #Remove the summed data and beitelag data that has NAs
#s2<-sampleRandom(herbivoredensSp,5000)

#pcaClim<-prcomp(s1,scale=T)
#dput(pcaClim,"S:\\SDMs\\PCAs\\ClimatePca.csv")
pcaClim<-dget("S:\\SDMs\\PCAs\\ClimatePca.csv")
#pcaHerb<-prcomp(s2,scale=T)
#dput(pcaHerb,"S:\\SDMs\\PCAs\\HerbivorePca.csv")
pcaHerb<-dget("S:\\SDMs\\PCAs\\HerbivorePca.csv")

x11(12,8)
par(mfrow=c(2,3))
mar=c(4,4,2,2))
barplot((unclass(summary(pcaClim))$importance)[2,],xlab="Axis",ylab="Proportion of variance",ylim=c(0,0.65),las=1,main='Climate')
biplot(pcaClim,xlabs=rep(".", 5000),cex=c(1,0.7),las=1)
biplot(pcaClim,xlabs=rep(".", 5000),cex=c(1,0.7),choices=c(2,3),las=1)

barplot((unclass(summary(pcaHerb))$importance)[2,],xlab="Axis",ylab="Proportion of variance",ylim=c(0,0.65),las=1,main='Herbivores')
biplot(pcaHerb,xlabs=rep(".", 5000),cex=c(1,0.7))
biplot(pcaHerb,xlabs=rep(".", 5000),cex=c(1,0.7),choices=c(2,3))

#Climate data can be reduced to 2-3 variables (could use axes, but will use MST, MAP and PrecipSeasonality for interpretation ease)
#Herbivore data cannot be reduced in dimensions this way. Will sum herbivore metabolic biomass by main habitat

#s3<-sampleRandom(predvars,5000)
#pcaAll<-prcomp(s3,scale=T)
#dput(pcaAll,"S:\\SDMs\\PCAs\\AllvarsPca.csv")
pcaAll<-dget("S:\\SDMs\\PCAs\\AllvarsPca.csv")
barplot((unclass(summary(pcaAll))$importance)[2,],xlab="Axis",ylab="Proportion of variance",ylim=c(0,0.65),las=1,main='All variables')
biplot(pcaAll,xlabs=rep(".", nrow(s3)),cex=c(1,0.7))
biplot(pcaAll,xlabs=rep(".", nrow(s3)),cex=c(1,0.7),choices=c(2,3))

#Selected variables
as.data.frame(names(predvars))
selectvars1<-subset(predvars,c(1,14,25,27,30,35))
#pairs(selectvars1) #High correlation between MST and elevation so drop elevation (r=-0.74)
selectvars<-subset(predvars,c(14,25,27,30,35))
selectvars$AR50TYPE100m_ProjectRaster1[selectvars$AR50TYPE100m_ProjectRaster1>81]<-NA #Change sea and not mapped to NA
selectvars<-stack(selectvars,soilphno)
#Selected varaibles without herbivores
selectvarsnoherb<-subset(selectvars,2:5)

selectvarsP<-selectvars
names(selectvarsP)<-c("Tundra_herbivores","MST","MAP","Precip.season.", "Land_Cover","Soil pH")
selectvarsP$Land_Cover<-as.factor(selectvarsP$Land_Cover)
ratlc<- levels(selectvarsP$Land_Cover)[[1]]
ratlc[["Land_Cover"]] <- list("Built up","Agricultural","Forest","Natural vegetation","Mires","Glaciers/Ice/Snow","Freshwater")
levels(selectvarsP$Land_Cover)<-ratlc
levelplot(selectvarsP$Land_Cover)

p1<-levelplot(selectvarsP$MST/10,margin=F,colorkey=list(height=0.5),main=(expression(bold(paste(MST ~(~degree~C))))))
p2<-levelplot(selectvarsP$MAP,margin=F,colorkey=list(height=0.5),main="MAP (mm)")
p3<-levelplot(selectvarsP$Precip.season.,margin=F,colorkey=list(height=0.5),main="Precipitation seasonality")
p4<-levelplot(selectvarsP$Tundra_herbivores+1,zscaleLog=T,margin=F,colorkey=list(height=0.5),main=(expression(bold(paste("Tundra herbivores (kg km"^-2*')')))))
p5<-levelplot(selectvarsP$Land_Cover,margin=F,colorkey=list(height=0.5),main="Land cover")
names(p1$legend)<-'inside'
p1$legend$inside$x <- 0.6
p1$legend$inside$y <- 0.3
names(p2$legend)<-'inside'
p2$legend$inside$x <- 0.6
p2$legend$inside$y <- 0.3
names(p3$legend)<-'inside'
p3$legend$inside$x <- 0.6
p3$legend$inside$y <- 0.3
names(p4$legend)<-'inside'
p4$legend$inside$x <- 0.6
p4$legend$inside$y <- 0.3
names(p5$legend)<-'inside'
p5$legend$inside$x <- 1.1
p5$legend$inside$y <- 0.3
print(p1,split=c(1,1,3,2))
print(p2,split=c(2,1,3,2),newpage=F)
print(p3,split=c(3,1,3,2),newpage=F)
print(p4,split=c(1,2,3,2),newpage=F)
print(p5,split=c(2,2,3,2),newpage=F)

#Species data
#GBIF.org (9th February 2016) GBIF Occurrence Download http://doi.org/10.15468/dl.aisn6d
#Vascular plant species on the 2010 Norwegian red list, classificed as CR, EN, VU or NT. 
#Main habitat in alpine
#Habitat impact as the reason for red listing.
#alpspdat1<-read.table("S:\\SDMs\\SpeciesData\\RedListSpQC.csv",header=T,sep=",")
alpspdat1<-read.table("~/SDMs/SpeciesData/RedListSpQC.csv",header=T,sep=",")

with(alpspdat1,tapply(species,species,length))
alpspdat<-alpspdat1[!is.na(alpspdat1$Xutm) & !is.na(alpspdat1$Yutm),]
with(alpspdat,tapply(species,species,length))

alpspdat_new<-alpspdat[alpspdat$year>=1995,]

Botlan<-alpspdat_new[alpspdat_new$species=='Botrychium lanceolatum',]
Comten<-alpspdat_new[alpspdat_new$species=='Comastoma tenellum',]
Gencam<-alpspdat_new[alpspdat_new$species=='Gentianella campestris',]
Kobsim<-alpspdat_new[alpspdat_new$species=='Kobresia simpliciuscula',]
Prisca<-alpspdat_new[alpspdat_new$species=='Primula scandinavica',]
Psealb<-alpspdat_new[alpspdat_new$species=='Pseudorchis albida',]
Pulver<-alpspdat_new[alpspdat_new$species=='Pulsatilla vernalis',]

allselectedsp<-droplevels(rbind(Botlan,Comten,Gencam,Kobsim,Prisca,Psealb,Pulver))
#write.csv(allselectedsp[,c(13,2,3)],"S:\\SDMs\\SDM_maxent\\MIAT\\AllSelectedSp.csv",row.names=F)
with(allselectedsp,tapply(species,species,length))

#Make a bias file
biasfile<-raster('~/SpeciesObservations/BiasFiles/VascularPlantBias_allrecords.asc')
bg<-randomPoints(biasfile,10000)
#bg_BC<-randomPoints(biasfile,10000,prob=T) #Weighted selection of background by biasfile
#write.table(bg_BC,"~/SpeciesObservations/SDMs/SpeciesData/BackgroundBiasCorrected.txt")
bg_BC<-read.table("~/SDMs/SpeciesData/BackgroundBiasCorrected.txt",header=T)
plot(biasfile)
points(bg,pch='.')
points(bg_BC,pch='.',col=2)

#Tuning
#Priscaeval_bg<-ENMevaluate(Prisca[,c(2:3)],selectvars,a=bg,method="randomkfold",kfolds=5,categoricals="AR50TYPE100m_ProjectRaster1",RMvalues=c(0.5,1,1.5,2,3,4,6,8,10,16))
#Priscaeval_bgBC<-ENMevaluate(Prisca[,c(2:3)],selectvars,a=bg_BC,method="randomkfold",kfolds=5,categoricals="AR50TYPE100m_ProjectRaster1",RMvalues=c(0.5,1,1.5,2,3,4,6,8,10,16))
#Priscaeval_bg_noherb<-ENMevaluate(Prisca[,c(2:3)],selectvarsnoherb,a=bg,method="randomkfold",kfolds=5,categoricals="AR50TYPE100m_ProjectRaster1",RMvalues=c(0.5,1,1.5,2,3,4,6,8,10,16))
#Priscaeval_bgBC_noherb<-ENMevaluate(Prisca[,c(2:3)],selectvarsnoherb,a=bg_BC,method="randomkfold",kfolds=5,categoricals="AR50TYPE100m_ProjectRaster1",RMvalues=c(0.5,1,1.5,2,3,4,6,8,10,16))

#Evaluation with background bias corrected data
#################These can be loaded below if required
#Botlaneval<-ENMevaluate(Botlan[,c(2:3)],selectvars,a=bg,method="randomkfold",categoricals="AR50TYPE100m_ProjectRaster1",kfolds=5)
#Comteneval<-ENMevaluate(Comten[,c(2:3)],selectvars,a=bg,method="randomkfold",categoricals="AR50TYPE100m_ProjectRaster1",kfolds=5)
#Gencameval<-ENMevaluate(Gencam[,c(2:3)],selectvars,a=bg,method="randomkfold",categoricals="AR50TYPE100m_ProjectRaster1",kfolds=5)
#Kobsimeval<-ENMevaluate(Kobsim[,c(2:3)],selectvars,a=bg,method="randomkfold",categoricals="AR50TYPE100m_ProjectRaster1",kfolds=5)
#Priscaeval<-ENMevaluate(Prisca[,c(2:3)],selectvars,a=bg,method="randomkfold",categoricals="AR50TYPE100m_ProjectRaster1",kfolds=5)
#Psealbeval<-ENMevaluate(Psealb[,c(2:3)],selectvars,a=bg,method="randomkfold",categoricals="AR50TYPE100m_ProjectRaster1",kfolds=5)
#Pulvereval<-ENMevaluate(Pulver[,c(2:3)],selectvars,a=bg,method="randomkfold",categoricals="AR50TYPE100m_ProjectRaster1",kfolds=5)
#save.image('~/SDMs/R1_ENMevaluations.RData')

#load('~/SDMs/ENMevaluations.RData')
load('~/SDMs/R1_ENMevaluations.RData')
#tiff('~/SDMs/R1_Tuning.tif',width=7,height=10,units='cm',res=600)
par(mfrow=c(4,2))
eval.plot(Botlaneval@results)
title(main="Botrychium lanceolatum",line=0.2,font.main=4)
eval.plot(Comteneval@results)
title(main='Comastoma tenellum',line=0.2,font.main=4)
eval.plot(Gencameval@results)
title(main='Gentianella campestris',line=0.2,font.main=4)
eval.plot(Kobsimeval@results)
title(main='Kobresia simpliciuscula',line=0.2,font.main=4)
eval.plot(Priscaeval@results)
title(main='Primula scandinavica',line=0.2,font.main=4)
eval.plot(Psealbeval@results)
title(main='Pseudorchis albida',line=0.2,font.main=4)
eval.plot(Pulvereval@results)
title(main='Pulsatilla vernalis',line=0.2,font.main=4)
#dev.off()


#Tuned parameters - Replicated
argsBotlan=c("-P","-J","replicates=5","betamultiplier=1.5")
argsComten=c("-P","-J","replicates=5","betamultiplier=2.0")
argsGencam=c("-P","-J","replicates=5","betamultiplier=2.5")
argsKobsim=c("-P","-J","replicates=5","betamultiplier=1.5")
argsPrisca=c("-P","-J","replicates=5","betamultiplier=0.5") 
argsPsealb=c("-P","-J","replicates=5","betamultiplier=1.0")
argsPulver=c("-P","-J","replicates=5","betamultiplier=2.0")

#Tuned parameters - single
argsBotlan_1=c("-P","-J","replicates=1","betamultiplier=1.5")
argsComten_1=c("-P","-J","replicates=1","betamultiplier=2.0")
argsGencam_1=c("-P","-J","replicates=1","betamultiplier=2.5")
argsKobsim_1=c("-P","-J","replicates=1","betamultiplier=1.5")
argsPrisca_1=c("-P","-J","replicates=1","betamultiplier=0.5") 
argsPsealb_1=c("-P","-J","replicates=1","betamultiplier=1.0")
argsPulver_1=c("-P","-J","replicates=1","betamultiplier=2.0")


#Model and prediction plots using parameters selected above (unreplicated)
meBotlan<-maxent(selectvars,Botlan[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=argsBotlan_1)
meGencam<-maxent(selectvars,Gencam[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=argsGencam_1)
meComten<-maxent(selectvars,Comten[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=argsComten_1)
meKobsim<-maxent(selectvars,Kobsim[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=argsKobsim_1)
mePrisca<-maxent(selectvars,Prisca[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=argsPrisca_1)
mePsealb<-maxent(selectvars,Psealb[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=argsPsealb_1)
mePulver<-maxent(selectvars,Pulver[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=argsPulver_1)

pBotlan<-predict(meBotlan,selectvars)
pGencam<-predict(meGencam,selectvars)
pComten<-predict(meComten,selectvars)
pKobsim<-predict(meKobsim,selectvars)
pPrisca<-predict(mePrisca,selectvars)
pPsealb<-predict(mePsealb,selectvars)
pPulver<-predict(mePulver,selectvars)

#Species occurances spatial dataframes
Botlanspdf<-SpatialPointsDataFrame(Botlan[,2:3],Botlan,proj4string=CRS("+proj=utm +zone=32"))
Gencamspdf<-SpatialPointsDataFrame(Gencam[,2:3],Gencam,proj4string=CRS("+proj=utm +zone=32"))
Comtenspdf<-SpatialPointsDataFrame(Comten[,2:3],Comten,proj4string=CRS("+proj=utm +zone=32"))
Kobsimspdf<-SpatialPointsDataFrame(Kobsim[,2:3],Kobsim,proj4string=CRS("+proj=utm +zone=32"))
Priscaspdf<-SpatialPointsDataFrame(Prisca[,2:3],Prisca,proj4string=CRS("+proj=utm +zone=32"))
Psealbspdf<-SpatialPointsDataFrame(Psealb[,2:3],Psealb,proj4string=CRS("+proj=utm +zone=32"))
Pulverspdf<-SpatialPointsDataFrame(Pulver[,2:3],Pulver,proj4string=CRS("+proj=utm +zone=32"))
#List these
spocc<-list(Botlanspdf,Comtenspdf,Gencamspdf,Kobsimspdf,Priscaspdf,Psealbspdf,Pulverspdf)

#Plot predictions on the same axis by first stacking them up
sp<-stack(pBotlan,pComten,pGencam,pKobsim,pPrisca,pPsealb,pPulver)
names(sp)<-c('Botrychium lanceolatum','Comastoma tenellum','Gentianella campestris','Kobresia simpliciuscula'
             ,'Primula scandinavica','Pseudorchis albida','Pulsatilla vernalis')
levelplot(sp,main="Model Predictions") +
  layer(sp.points(spocc[[panel.number()]],
                  pch='.', cex=0.1, col="black"))

#Replicated maxent runs
meBotlan_R<-maxent(selectvars,Botlan[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=argsBotlan,path="~/SDMs/HerbivoreSDM_final/Botlan")
meGencam_R<-maxent(selectvars,Gencam[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=c(argsGencam),path="~/SDMs/HerbivoreSDM_final/Gencam")
meComten_R<-maxent(selectvars,Comten[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=c(argsComten),path="~/SDMs/HerbivoreSDM_final/Comten")
meKobsim_R<-maxent(selectvars,Kobsim[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=argsKobsim,path="~/SDMs/HerbivoreSDM_final/Kobsim")
mePrisca_R<-maxent(selectvars,Prisca[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=c(argsPrisca),path="~/SDMs/HerbivoreSDM_final/Prisca")
mePsealb_R<-maxent(selectvars,Psealb[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=c(argsPsealb),path="~/SDMs/HerbivoreSDM_final/Psealb")
mePulver_R<-maxent(selectvars,Pulver[,c(2:3)],a=bg_BC,factors="AR50TYPE100m_ProjectRaster1",args=c(argsPulver),path="~/SDMs/HerbivoreSDM_final/Pulver")


#MaxEnt result files
merBotlan<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\maxentResults.csv",sep=",",header=T,comment.char="")
merGencam<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\maxentResults.csv",sep=",",header=T,comment.char="")
merComten<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\maxentResults.csv",sep=",",header=T,comment.char="")
merKobsim<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\maxentResults.csv",sep=",",header=T,comment.char="")
merPrisca<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\maxentResults.csv",sep=",",header=T,comment.char="")
merPsealb<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\maxentResults.csv",sep=",",header=T,comment.char="")
merPulver<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\maxentResults.csv",sep=",",header=T,comment.char="")

#Percentage contribution figure
merallsp<-rbind(merBotlan,merComten,merGencam,merKobsim,merPrisca,merPsealb,merPulver)
merallsp$speciesname<-rep(c('Botrychium lanceolatum','Comastoma tenellum','Gentianella campestris','Kobresia simpliciuscula'
                            ,'Primula scandinavica','Pseudorchis albida','Pulsatilla vernalis'),each=6)
merallspreps<-merallsp[merallsp$Species!='species (average)',c(12:(12+nlayers(selectvars)-1),120)]
meanvarcont<-aggregate(merallspreps[,1:nlayers(selectvars)],list(merallspreps$speciesname),mean)
sem<-function(x)sd(x,na.rm=T)/sqrt(length(!is.na(x)))
meanvarsem<-aggregate(merallspreps[,1:nlayers(selectvars)],list(merallspreps$speciesname),sd)
m1<-as.matrix(t(meanvarcont[,2:6]))
s1<-as.matrix(t(meanvarsem[,2:6]))
colnames(m1)<-levels(as.factor(merallsp$speciesname))
rownames(m1)<-c("Land cover","MAP", "MST","Precipitation seasonality","Herbivore density")
par(mar=c(10,5,3,1))
b1<-barplot((m1),beside=T,legend.text=rownames(m1),las=2,args.legend=list(x='topl',ncol=2,cex=0.8),ylab='Percentage contribution',ylim=c(0,80))
arrows(b1,m1+s1,b1,m1-s1,length=0.05,code=3,angle=90)
aucmean<-with(merallsp[merallsp$Species!='species (average)',],tapply(Test.AUC,speciesname,mean))
mtext(side=3,at=b1[3,],format(aucmean,trim=T,digits=2))

#Perm importance
merallsp<-rbind(merBotlan,merComten,merGencam,merKobsim,merPrisca,merPsealb,merPulver)
merallsp$speciesname<-rep(c('Botrychium lanceolatum','Comastoma tenellum','Gentianella campestris','Kobresia simpliciuscula'
                            ,'Primula scandinavica','Pseudorchis albida','Pulsatilla vernalis'),each=6)
merallspreps<-merallsp[merallsp$Species!='species (average)',c(17:(17+nlayers(selectvars)-1),120)]
meanvarcont<-aggregate(merallspreps[,1:nlayers(selectvars)],list(merallspreps$speciesname),mean)
sem<-function(x)sd(x,na.rm=T)/sqrt(length(!is.na(x)))
meanvarsem<-aggregate(merallspreps[,1:nlayers(selectvars)],list(merallspreps$speciesname),sd)
m1<-as.matrix(t(meanvarcont[,2:6]))
s1<-as.matrix(t(meanvarsem[,2:6]))
colnames(m1)<-levels(as.factor(merallsp$speciesname))
rownames(m1)<-c("Land cover","MAP", "MST","Precipitation seasonality","Herbivore density")
par(mar=c(10,5,3,1))
b1<-barplot((m1),beside=T,legend.text=rownames(m1),las=2,args.legend=list(x='topl',ncol=2,cex=0.8),ylab='Permutation importance',ylim=c(0,80))
arrows(b1,m1+s1,b1,m1-s1,length=0.05,code=3,angle=90)
aucmean<-with(merallsp[merallsp$Species!='species (average)',],tapply(Test.AUC,speciesname,mean))
mtext(side=3,at=b1[3,],format(aucmean,trim=T,digits=2))


#################################################################################################################
#Make a limiting factors map
namedf<-cbind(shortname=c('Botlan','Comten','Gencam','Kobsim','Prisca','Psealb','Pulver')
              ,fullname=c('Botrychium lanceolatum','Comastoma tenellum','Gentianella campestris','Kobresia simpliciuscula'  ,'Primula scandinavica','Pseudorchis albida','Pulsatilla vernalis'))
                                                                                                

#Write predictor variables to file into folder where maxent will be run through command prompt 
#(only rerun if predictor variables changed)
#writeRaster(selectvars,bylayer=T,file.path('C:\\test\\botlan\\current', names(selectvars)),,format="ascii",overwrite=T)
#writeRaster(selectvars,bylayer=T,file.path('C:\\test\\Gencam\\current', names(selectvars)),,format="ascii",overwrite=T)
#writeRaster(selectvars,bylayer=T,file.path('C:\\test\\Comten\\current', names(selectvars)),,format="ascii",overwrite=T)
#writeRaster(selectvars,bylayer=T,file.path('C:\\test\\Kobsim\\current', names(selectvars)),,format="ascii",overwrite=T)
#writeRaster(selectvars,bylayer=T,file.path('C:\\test\\Prisca\\current', names(selectvars)),,format="ascii",overwrite=T)
#writeRaster(selectvars,bylayer=T,file.path('C:\\test\\Psealb\\current', names(selectvars)),,format="ascii",overwrite=T)
#writeRaster(selectvars,bylayer=T,file.path('C:\\test\\Pulver\\current', names(selectvars)),,format="ascii",overwrite=T)

#Save an unreplicated maxent run to the same folder
meBotlan<-maxent(selectvars,Botlan[,c(2:3)],factors="AR50TYPE100m_ProjectRaster1",args=c(argsBotlan,"randomtestpoints=5"),path="C:\\test\\botlan\\model")
meGencam<-maxent(selectvars,Gencam[,c(2:3)],factors="AR50TYPE100m_ProjectRaster1",args=c(argsGencam,"randomtestpoints=5"),path="C:\\test\\Gencam\\model")
meComten<-maxent(selectvars,Comten[,c(2:3)],factors="AR50TYPE100m_ProjectRaster1",args=c(argsComten,"randomtestpoints=5"),path="C:\\test\\Comten\\model")
meKobsim<-maxent(selectvars,Kobsim[,c(2:3)],factors="AR50TYPE100m_ProjectRaster1",args=c(argsKobsim,"randomtestpoints=5"),path="C:\\test\\Kobsim\\model")
mePrisca<-maxent(selectvars,Prisca[,c(2:3)],factors="AR50TYPE100m_ProjectRaster1",args=c(argsPrisca,"randomtestpoints=5"),path="C:\\test\\Prisca\\model")
mePsealb<-maxent(selectvars,Psealb[,c(2:3)],factors="AR50TYPE100m_ProjectRaster1",args=c(argsPsealb,"randomtestpoints=5"),path="C:\\test\\Psealb\\model")
mePulver<-maxent(selectvars,Pulver[,c(2:3)],factors="AR50TYPE100m_ProjectRaster1",args=c(argsPulver,"randomtestpoints=5"),path="C:\\test\\Pulver\\model")

#Now go to the command prompt
#Paste the following lines in
cd c:\test\botlan
java -cp maxent.jar density.tools.LimitingFactor model\species.lambdas current lf.asc
cd c:\test\Gencam
java -cp maxent.jar density.tools.LimitingFactor model\species.lambdas current lf.asc
cd c:\test\comten
java -cp maxent.jar density.tools.LimitingFactor model\species.lambdas current lf.asc
cd c:\test\kobsim
java -cp maxent.jar density.tools.LimitingFactor model\species.lambdas current lf.asc
cd c:\test\prisca
java -cp maxent.jar density.tools.LimitingFactor model\species.lambdas current lf.asc
cd c:\test\psealb
java -cp maxent.jar density.tools.LimitingFactor model\species.lambdas current lf.asc
cd c:\test\pulver
java -cp maxent.jar density.tools.LimitingFactor model\species.lambdas current lf.asc

#Read in the limiting factor rasters as factors
botlan_lf<-as.factor(raster("C:\\test\\botlan\\lf.asc"))
Gencam_lf<-as.factor(raster("C:\\test\\Gencam\\lf.asc"))
Comten_lf<-as.factor(raster("C:\\test\\Comten\\lf.asc"))
Kobsim_lf<-as.factor(raster("C:\\test\\Kobsim\\lf.asc"))
Prisca_lf<-as.factor(raster("C:\\test\\Prisca\\lf.asc"))
Psealb_lf<-as.factor(raster("C:\\test\\Psealb\\lf.asc"))
Pulver_lf<-as.factor(raster("C:\\test\\Pulver\\lf.asc"))

#Make a rat table for all the limiting factor maps and relevel the factors
ratlf<- levels(botlan_lf)[[1]]
ratlf[["Limiting_factor"]] <- c("Land Cover","MAP","MST","Precip season.","Tundra herbivores")
levels(botlan_lf)<-ratlf
levels(Gencam_lf)<-ratlf
levels(Comten_lf)<-ratlf
levels(Kobsim_lf)<-ratlf
levels(Prisca_lf)<-ratlf
levels(Psealb_lf)<-ratlf
levels(Pulver_lf)<-ratlf

norway<-getData('GADM', country='Norway', level=0)
norwayP<-spTransform(norway,"+proj=utm +zone=32")
myColors <- (brewer.pal(5, "YlOrRd"))
myKey <- list(text=list(lab=ratlf$Limiting_factor), rectangles=list(col = myColors),space='inside')
levelplot(botlan_lf,col.regions=myColors,colorkey=F,key=myKey) #+
  layer(sp.polygons(norwayP,lwd=0.1))

#Stack together
stacklf<-stack(botlan_lf,Comten_lf,Gencam_lf,Kobsim_lf,Prisca_lf,Psealb_lf,Pulver_lf)
#Note that the raster table output places the predictor variables in alphabetical order (0:n(layers))
lpp<-levelplot(stacklf,names=namedf[,2],col.regions=myColors,colorkey=F,key=myKey,main='Limiting Factors')
#lpp<-levelplot(stacklf,names=namedf[,2],colorkey=F)
  names(lpp$legend)<-'inside'
lpp$legend$inside$x <- 0.75
lpp$legend$inside$y <- 0.3
lpp+
  layer(sp.polygons(norwayP,,alpha=0.5))
  #levelplot(stacklf,names=namedf[,2],col.regions=myColors,colorkey=F,key=myKey)+
  #layer(sp.polygons(norwayP,,alpha=0.5))

lfdf<-matrix(summary(as.factor(getValues(botlan_lf)))[1:5],ncol=5)
lfdf<-rbind(lfdf,summary(as.factor(getValues(Comten_lf)))[1:5])
lfdf<-rbind(lfdf,summary(as.factor(getValues(Gencam_lf)))[1:5])
lfdf<-rbind(lfdf,summary(as.factor(getValues(Kobsim_lf)))[1:5])
lfdf<-rbind(lfdf,summary(as.factor(getValues(Prisca_lf)))[1:5])
lfdf<-rbind(lfdf,summary(as.factor(getValues(Psealb_lf)))[1:5])
lfdf<-rbind(lfdf,summary(as.factor(getValues(Pulver_lf)))[1:5])
colnames(lfdf)<-c("Land Cover","MAP","MST","Precip season.","Tundra herbivores")
rownames(lfdf)<-namedf[,2]
lfdf<-as.data.frame(lfdf)
lfdf$total<-rowSums(lfdf)
lfdfRel<-lfdf[,1:5]/lfdf[,6]
round(lfdf,2)
round(lfdfRel,2)
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#Response curves
#Note there are a lot of lines!

#Order variables by average contribution to models
sort(rowMeans(m1))

x11(11,8,pointsize=6)
par(mfcol=c(5,7))

namedf<-cbind(shortname=c('Botlan','Gencam','Comten','Kobsim','Prisca','Psealb','Pulver'),fullname=c('Botrychium lanceolatum','Gentianella campestris','Comastoma tenellum','Kobresia simpliciuscula'
                                                                               ,'Primula scandinavica','Pseudorchis albida','Pulsatilla vernalis'))
#Botlan
#MAP
Botlan_map0<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_0_AnnualPrecip.dat",sep=",",header=T)
Botlan_map1<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_1_AnnualPrecip.dat",sep=",",header=T)
Botlan_map2<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_2_AnnualPrecip.dat",sep=",",header=T)
Botlan_map3<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_3_AnnualPrecip.dat",sep=",",header=T)
Botlan_map4<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_4_AnnualPrecip.dat",sep=",",header=T)
Botlan_map_meanresp<-apply(cbind(Botlan_map0$y,Botlan_map1$y,Botlan_map2$y,Botlan_map3$y,Botlan_map4$y),1,mean)
Botlan_map_sdresp<-apply(cbind(Botlan_map0$y,Botlan_map1$y,Botlan_map2$y,Botlan_map3$y,Botlan_map4$y),1,sd)
Botlan_map_respdf<-data.frame(cbind(Botlan_map0$x,Botlan_map_meanresp,Botlan_map_sdresp))
names(Botlan_map_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Botlan_map_respdf[Botlan_map_respdf$x >= 0,],plot(x,Botlan_map_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Botlan'),2]
                                         ,xlab="Mean annual precipitation (mm)",ylab="Logistic output (probability of pressence)"))
with(Botlan_map_respdf[Botlan_map_respdf$x >= 0,],lines(x,Botlan_map_meanresp+Botlan_map_sdresp,type="l",lty=2))
with(Botlan_map_respdf[Botlan_map_respdf$x >= 0,],lines(x,Botlan_map_meanresp-Botlan_map_sdresp,type="l",lty=2))

#MST
Botlan_mst0<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_0_MeanTempWarmQuart.dat",sep=",",header=T)
Botlan_mst1<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_1_MeanTempWarmQuart.dat",sep=",",header=T)
Botlan_mst2<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_2_MeanTempWarmQuart.dat",sep=",",header=T)
Botlan_mst3<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_3_MeanTempWarmQuart.dat",sep=",",header=T)
Botlan_mst4<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_4_MeanTempWarmQuart.dat",sep=",",header=T)
Botlan_mst_meanresp<-apply(cbind(Botlan_mst0$y,Botlan_mst1$y,Botlan_mst2$y,Botlan_mst3$y,Botlan_mst4$y),1,mean)
Botlan_mst_sdresp<-apply(cbind(Botlan_mst0$y,Botlan_mst1$y,Botlan_mst2$y,Botlan_mst3$y,Botlan_mst4$y),1,sd)
Botlan_mst_respdf<-data.frame(cbind(Botlan_mst0$x,Botlan_mst_meanresp,Botlan_mst_sdresp))
names(Botlan_mst_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Botlan_mst_respdf[Botlan_mst_respdf$x >= 0,],plot(x/10,Botlan_mst_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Botlan'),2]
                                         ,xlab=expression(Temperature~(degree~C)) ,ylab="Logistic output (probability of pressence)"))
with(Botlan_mst_respdf[Botlan_mst_respdf$x >= 0,],lines(x/10,Botlan_mst_meanresp+Botlan_mst_sdresp,type="l",lty=2))
with(Botlan_mst_respdf[Botlan_mst_respdf$x >= 0,],lines(x/10,Botlan_mst_meanresp-Botlan_mst_sdresp,type="l",lty=2))

#Tundra herbivores(continuous)
Botlan_tundraherb0<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_0_TundraHerbivores.dat",sep=",",header=T)
Botlan_tundraherb1<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_1_TundraHerbivores.dat",sep=",",header=T)
Botlan_tundraherb2<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_2_TundraHerbivores.dat",sep=",",header=T)
Botlan_tundraherb3<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_3_TundraHerbivores.dat",sep=",",header=T)
Botlan_tundraherb4<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_4_TundraHerbivores.dat",sep=",",header=T)
Botlan_tundraherb_meanresp<-apply(cbind(Botlan_tundraherb0$y,Botlan_tundraherb1$y,Botlan_tundraherb2$y,Botlan_tundraherb3$y,Botlan_tundraherb4$y),1,mean)
Botlan_tundraherb_sdresp<-apply(cbind(Botlan_tundraherb0$y,Botlan_tundraherb1$y,Botlan_tundraherb2$y,Botlan_tundraherb3$y,Botlan_tundraherb4$y),1,sd)
Botlan_tundraherb_respdf<-data.frame(cbind(Botlan_tundraherb0$x,Botlan_tundraherb_meanresp,Botlan_tundraherb_sdresp))
names(Botlan_tundraherb_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Botlan_tundraherb_respdf[Botlan_tundraherb_respdf$x >= 0,],plot(x,Botlan_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Botlan'),2]
                                                       ,xlab=expression("Herbivore density (kg km"^-2*")"),ylab="Logistic output (probability of pressence)"))
with(Botlan_tundraherb_respdf[Botlan_tundraherb_respdf$x >= 0,],lines(x,Botlan_tundraherb_meanresp+Botlan_tundraherb_sdresp,type="l",lty=2))
with(Botlan_tundraherb_respdf[Botlan_tundraherb_respdf$x >= 0,],lines(x,Botlan_tundraherb_meanresp-Botlan_tundraherb_sdresp,type="l",lty=2))

#Precipitation seasonality
Botlan_precipseas0<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_0_PrecipSeasonality.dat",sep=",",header=T)
Botlan_precipseas1<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_1_PrecipSeasonality.dat",sep=",",header=T)
Botlan_precipseas2<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_2_PrecipSeasonality.dat",sep=",",header=T)
Botlan_precipseas3<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_3_PrecipSeasonality.dat",sep=",",header=T)
Botlan_precipseas4<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_4_PrecipSeasonality.dat",sep=",",header=T)
Botlan_precipseas_meanresp<-apply(cbind(Botlan_precipseas0$y,Botlan_precipseas1$y,Botlan_precipseas2$y,Botlan_precipseas3$y,Botlan_precipseas4$y),1,mean)
Botlan_precipseas_sdresp<-apply(cbind(Botlan_precipseas0$y,Botlan_precipseas1$y,Botlan_precipseas2$y,Botlan_precipseas3$y,Botlan_precipseas4$y),1,sd)
Botlan_precipseas_respdf<-data.frame(cbind(Botlan_precipseas0$x,Botlan_precipseas_meanresp,Botlan_precipseas_sdresp))
names(Botlan_precipseas_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Botlan_precipseas_respdf[Botlan_precipseas_respdf$x >= 0,],plot(x,Botlan_precipseas_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Botlan'),2]
                                                                     ,xlab="Precipitation Seasonality (CV)",ylab="Logistic output (probability of pressence)"))
with(Botlan_precipseas_respdf[Botlan_precipseas_respdf$x >= 0,],lines(x,Botlan_precipseas_meanresp+Botlan_precipseas_sdresp,type="l",lty=2))
with(Botlan_precipseas_respdf[Botlan_precipseas_respdf$x >= 0,],lines(x,Botlan_precipseas_meanresp-Botlan_precipseas_sdresp,type="l",lty=2))

#Land cover (factor)
Botlan_landcover0<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_0_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Botlan_landcover1<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_1_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Botlan_landcover2<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_2_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Botlan_landcover3<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_3_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Botlan_landcover4<-read.table("S:\\SDMs\\SDM_maxent\\Botlan\\plots\\species_4_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Botlan_landcover_meanresp<-apply(cbind(Botlan_landcover0$y,Botlan_landcover1$y,Botlan_landcover2$y,Botlan_landcover3$y,Botlan_landcover4$y),1,mean)
Botlan_landcover_sdresp<-apply(cbind(Botlan_landcover0$y,Botlan_landcover1$y,Botlan_landcover2$y,Botlan_landcover3$y,Botlan_landcover4$y),1,sd)
Botlan_landcover_respdf<-data.frame(cbind(Botlan_landcover0$x,Botlan_landcover_meanresp,Botlan_landcover_sdresp))
names(Botlan_landcover_respdf)[1]<-"x"
par(mar=c(8,5,3,1))
b1<-with(Botlan_landcover_respdf,barplot(Botlan_landcover_meanresp,ylim=c(0,1),main=namedf[which(namedf[,1]=='Botlan'),2],las=2
                                  ,names.arg=c("Built up","Agricultural","Forest","Natural vegetation \n (not forest)","Wetland","Ice/Snow","Freshwater")
                                  ,xlab="",ylab="Logistic output (probability of pressence)"))
with(Botlan_landcover_respdf,arrows(b1,Botlan_landcover_meanresp+Botlan_landcover_sdresp,b1,Botlan_landcover_meanresp-Botlan_landcover_sdresp
                             ,angle=90,code=3,length=0.05))

#Comten
#MAP
Comten_map0<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_0_AnnualPrecip.dat",sep=",",header=T)
Comten_map1<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_1_AnnualPrecip.dat",sep=",",header=T)
Comten_map2<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_2_AnnualPrecip.dat",sep=",",header=T)
Comten_map3<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_3_AnnualPrecip.dat",sep=",",header=T)
Comten_map4<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_4_AnnualPrecip.dat",sep=",",header=T)
Comten_map_meanresp<-apply(cbind(Comten_map0$y,Comten_map1$y,Comten_map2$y,Comten_map3$y,Comten_map4$y),1,mean)
Comten_map_sdresp<-apply(cbind(Comten_map0$y,Comten_map1$y,Comten_map2$y,Comten_map3$y,Comten_map4$y),1,sd)
Comten_map_respdf<-data.frame(cbind(Comten_map0$x,Comten_map_meanresp,Comten_map_sdresp))
names(Comten_map_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Comten_map_respdf[Comten_map_respdf$x >= 0,],plot(x,Comten_map_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Comten'),2]
                                                       ,xlab="Mean annual precipitation (mm)",ylab="Logistic output (probability of pressence)"))
with(Comten_map_respdf[Comten_map_respdf$x >= 0,],lines(x,Comten_map_meanresp+Comten_map_sdresp,type="l",lty=2))
with(Comten_map_respdf[Comten_map_respdf$x >= 0,],lines(x,Comten_map_meanresp-Comten_map_sdresp,type="l",lty=2))

#MST
Comten_mst0<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_0_MeanTempWarmQuart.dat",sep=",",header=T)
Comten_mst1<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_1_MeanTempWarmQuart.dat",sep=",",header=T)
Comten_mst2<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_2_MeanTempWarmQuart.dat",sep=",",header=T)
Comten_mst3<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_3_MeanTempWarmQuart.dat",sep=",",header=T)
Comten_mst4<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_4_MeanTempWarmQuart.dat",sep=",",header=T)
Comten_mst_meanresp<-apply(cbind(Comten_mst0$y,Comten_mst1$y,Comten_mst2$y,Comten_mst3$y,Comten_mst4$y),1,mean)
Comten_mst_sdresp<-apply(cbind(Comten_mst0$y,Comten_mst1$y,Comten_mst2$y,Comten_mst3$y,Comten_mst4$y),1,sd)
Comten_mst_respdf<-data.frame(cbind(Comten_mst0$x,Comten_mst_meanresp,Comten_mst_sdresp))
names(Comten_mst_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Comten_mst_respdf[Comten_mst_respdf$x >= 0,],plot(x/10,Comten_mst_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Comten'),2]
                                                       ,xlab=expression(Temperature~(degree~C)) ,ylab="Logistic output (probability of pressence)"))
with(Comten_mst_respdf[Comten_mst_respdf$x >= 0,],lines(x/10,Comten_mst_meanresp+Comten_mst_sdresp,type="l",lty=2))
with(Comten_mst_respdf[Comten_mst_respdf$x >= 0,],lines(x/10,Comten_mst_meanresp-Comten_mst_sdresp,type="l",lty=2))

#Tundra herbivores(continuous)
Comten_tundraherb0<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_0_TundraHerbivores.dat",sep=",",header=T)
Comten_tundraherb1<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_1_TundraHerbivores.dat",sep=",",header=T)
Comten_tundraherb2<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_2_TundraHerbivores.dat",sep=",",header=T)
Comten_tundraherb3<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_3_TundraHerbivores.dat",sep=",",header=T)
Comten_tundraherb4<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_4_TundraHerbivores.dat",sep=",",header=T)
Comten_tundraherb_meanresp<-apply(cbind(Comten_tundraherb0$y,Comten_tundraherb1$y,Comten_tundraherb2$y,Comten_tundraherb3$y,Comten_tundraherb4$y),1,mean)
Comten_tundraherb_sdresp<-apply(cbind(Comten_tundraherb0$y,Comten_tundraherb1$y,Comten_tundraherb2$y,Comten_tundraherb3$y,Comten_tundraherb4$y),1,sd)
Comten_tundraherb_respdf<-data.frame(cbind(Comten_tundraherb0$x,Comten_tundraherb_meanresp,Comten_tundraherb_sdresp))
names(Comten_tundraherb_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Comten_tundraherb_respdf[Comten_tundraherb_respdf$x >= 0,],plot(x,Comten_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Comten'),2]
                                                                     ,xlab=expression("Herbivore density (kg km"^-2*")"),ylab="Logistic output (probability of pressence)"))
with(Comten_tundraherb_respdf[Comten_tundraherb_respdf$x >= 0,],lines(x,Comten_tundraherb_meanresp+Comten_tundraherb_sdresp,type="l",lty=2))
with(Comten_tundraherb_respdf[Comten_tundraherb_respdf$x >= 0,],lines(x,Comten_tundraherb_meanresp-Comten_tundraherb_sdresp,type="l",lty=2))

#Precipitation seasonality
Comten_precipseas0<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_0_PrecipSeasonality.dat",sep=",",header=T)
Comten_precipseas1<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_1_PrecipSeasonality.dat",sep=",",header=T)
Comten_precipseas2<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_2_PrecipSeasonality.dat",sep=",",header=T)
Comten_precipseas3<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_3_PrecipSeasonality.dat",sep=",",header=T)
Comten_precipseas4<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_4_PrecipSeasonality.dat",sep=",",header=T)
Comten_precipseas_meanresp<-apply(cbind(Comten_precipseas0$y,Comten_precipseas1$y,Comten_precipseas2$y,Comten_precipseas3$y,Comten_precipseas4$y),1,mean)
Comten_precipseas_sdresp<-apply(cbind(Comten_precipseas0$y,Comten_precipseas1$y,Comten_precipseas2$y,Comten_precipseas3$y,Comten_precipseas4$y),1,sd)
Comten_precipseas_respdf<-data.frame(cbind(Comten_precipseas0$x,Comten_precipseas_meanresp,Comten_precipseas_sdresp))
names(Comten_precipseas_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Comten_precipseas_respdf[Comten_precipseas_respdf$x >= 0,],plot(x,Comten_precipseas_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Comten'),2]
                                                                     ,xlab="Precipitation Seasonality (CV)",ylab="Logistic output (probability of pressence)"))
with(Comten_precipseas_respdf[Comten_precipseas_respdf$x >= 0,],lines(x,Comten_precipseas_meanresp+Comten_precipseas_sdresp,type="l",lty=2))
with(Comten_precipseas_respdf[Comten_precipseas_respdf$x >= 0,],lines(x,Comten_precipseas_meanresp-Comten_precipseas_sdresp,type="l",lty=2))

#Land cover (factor)
Comten_landcover0<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_0_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Comten_landcover1<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_1_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Comten_landcover2<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_2_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Comten_landcover3<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_3_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Comten_landcover4<-read.table("S:\\SDMs\\SDM_maxent\\Comten\\plots\\species_4_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Comten_landcover_meanresp<-apply(cbind(Comten_landcover0$y,Comten_landcover1$y,Comten_landcover2$y,Comten_landcover3$y,Comten_landcover4$y),1,mean)
Comten_landcover_sdresp<-apply(cbind(Comten_landcover0$y,Comten_landcover1$y,Comten_landcover2$y,Comten_landcover3$y,Comten_landcover4$y),1,sd)
Comten_landcover_respdf<-data.frame(cbind(Comten_landcover0$x,Comten_landcover_meanresp,Comten_landcover_sdresp))
names(Comten_landcover_respdf)[1]<-"x"
par(mar=c(8,5,3,1))
b1<-with(Comten_landcover_respdf,barplot(Comten_landcover_meanresp,ylim=c(0,1),main=namedf[which(namedf[,1]=='Comten'),2],las=2
                                         ,names.arg=c("Built up","Agricultural","Forest","Natural vegetation \n (not forest)","Wetland","Ice/Snow","Freshwater")
                                         ,xlab="",ylab="Logistic output (probability of pressence)"))
with(Comten_landcover_respdf,arrows(b1,Comten_landcover_meanresp+Comten_landcover_sdresp,b1,Comten_landcover_meanresp-Comten_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))

#Gencam
#MAP
Gencam_map0<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_0_AnnualPrecip.dat",sep=",",header=T)
Gencam_map1<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_1_AnnualPrecip.dat",sep=",",header=T)
Gencam_map2<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_2_AnnualPrecip.dat",sep=",",header=T)
Gencam_map3<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_3_AnnualPrecip.dat",sep=",",header=T)
Gencam_map4<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_4_AnnualPrecip.dat",sep=",",header=T)
Gencam_map_meanresp<-apply(cbind(Gencam_map0$y,Gencam_map1$y,Gencam_map2$y,Gencam_map3$y,Gencam_map4$y),1,mean)
Gencam_map_sdresp<-apply(cbind(Gencam_map0$y,Gencam_map1$y,Gencam_map2$y,Gencam_map3$y,Gencam_map4$y),1,sd)
Gencam_map_respdf<-data.frame(cbind(Gencam_map0$x,Gencam_map_meanresp,Gencam_map_sdresp))
names(Gencam_map_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Gencam_map_respdf[Gencam_map_respdf$x >= 0,],plot(x,Gencam_map_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Gencam'),2]
                                                       ,xlab="Mean annual precipitation (mm)",ylab="Logistic output (probability of pressence)"))
with(Gencam_map_respdf[Gencam_map_respdf$x >= 0,],lines(x,Gencam_map_meanresp+Gencam_map_sdresp,type="l",lty=2))
with(Gencam_map_respdf[Gencam_map_respdf$x >= 0,],lines(x,Gencam_map_meanresp-Gencam_map_sdresp,type="l",lty=2))

#MST
Gencam_mst0<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_0_MeanTempWarmQuart.dat",sep=",",header=T)
Gencam_mst1<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_1_MeanTempWarmQuart.dat",sep=",",header=T)
Gencam_mst2<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_2_MeanTempWarmQuart.dat",sep=",",header=T)
Gencam_mst3<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_3_MeanTempWarmQuart.dat",sep=",",header=T)
Gencam_mst4<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_4_MeanTempWarmQuart.dat",sep=",",header=T)
Gencam_mst_meanresp<-apply(cbind(Gencam_mst0$y,Gencam_mst1$y,Gencam_mst2$y,Gencam_mst3$y,Gencam_mst4$y),1,mean)
Gencam_mst_sdresp<-apply(cbind(Gencam_mst0$y,Gencam_mst1$y,Gencam_mst2$y,Gencam_mst3$y,Gencam_mst4$y),1,sd)
Gencam_mst_respdf<-data.frame(cbind(Gencam_mst0$x,Gencam_mst_meanresp,Gencam_mst_sdresp))
names(Gencam_mst_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Gencam_mst_respdf[Gencam_mst_respdf$x >= 0,],plot(x/10,Gencam_mst_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Gencam'),2]
                                                       ,xlab=expression(Temperature~(degree~C)) ,ylab="Logistic output (probability of pressence)"))
with(Gencam_mst_respdf[Gencam_mst_respdf$x >= 0,],lines(x/10,Gencam_mst_meanresp+Gencam_mst_sdresp,type="l",lty=2))
with(Gencam_mst_respdf[Gencam_mst_respdf$x >= 0,],lines(x/10,Gencam_mst_meanresp-Gencam_mst_sdresp,type="l",lty=2))

#Tundra herbivores(continuous)
Gencam_tundraherb0<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_0_TundraHerbivores.dat",sep=",",header=T)
Gencam_tundraherb1<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_1_TundraHerbivores.dat",sep=",",header=T)
Gencam_tundraherb2<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_2_TundraHerbivores.dat",sep=",",header=T)
Gencam_tundraherb3<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_3_TundraHerbivores.dat",sep=",",header=T)
Gencam_tundraherb4<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_4_TundraHerbivores.dat",sep=",",header=T)
Gencam_tundraherb_meanresp<-apply(cbind(Gencam_tundraherb0$y,Gencam_tundraherb1$y,Gencam_tundraherb2$y,Gencam_tundraherb3$y,Gencam_tundraherb4$y),1,mean)
Gencam_tundraherb_sdresp<-apply(cbind(Gencam_tundraherb0$y,Gencam_tundraherb1$y,Gencam_tundraherb2$y,Gencam_tundraherb3$y,Gencam_tundraherb4$y),1,sd)
Gencam_tundraherb_respdf<-data.frame(cbind(Gencam_tundraherb0$x,Gencam_tundraherb_meanresp,Gencam_tundraherb_sdresp))
names(Gencam_tundraherb_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Gencam_tundraherb_respdf[Gencam_tundraherb_respdf$x >= 0,],plot(x,Gencam_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Gencam'),2]
                                                                     ,xlab=expression("Herbivore density (kg km"^-2*")"),ylab="Logistic output (probability of pressence)"))
with(Gencam_tundraherb_respdf[Gencam_tundraherb_respdf$x >= 0,],lines(x,Gencam_tundraherb_meanresp+Gencam_tundraherb_sdresp,type="l",lty=2))
with(Gencam_tundraherb_respdf[Gencam_tundraherb_respdf$x >= 0,],lines(x,Gencam_tundraherb_meanresp-Gencam_tundraherb_sdresp,type="l",lty=2))

#Precipitation seasonality
Gencam_precipseas0<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_0_PrecipSeasonality.dat",sep=",",header=T)
Gencam_precipseas1<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_1_PrecipSeasonality.dat",sep=",",header=T)
Gencam_precipseas2<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_2_PrecipSeasonality.dat",sep=",",header=T)
Gencam_precipseas3<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_3_PrecipSeasonality.dat",sep=",",header=T)
Gencam_precipseas4<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_4_PrecipSeasonality.dat",sep=",",header=T)
Gencam_precipseas_meanresp<-apply(cbind(Gencam_precipseas0$y,Gencam_precipseas1$y,Gencam_precipseas2$y,Gencam_precipseas3$y,Gencam_precipseas4$y),1,mean)
Gencam_precipseas_sdresp<-apply(cbind(Gencam_precipseas0$y,Gencam_precipseas1$y,Gencam_precipseas2$y,Gencam_precipseas3$y,Gencam_precipseas4$y),1,sd)
Gencam_precipseas_respdf<-data.frame(cbind(Gencam_precipseas0$x,Gencam_precipseas_meanresp,Gencam_precipseas_sdresp))
names(Gencam_precipseas_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Gencam_precipseas_respdf[Gencam_precipseas_respdf$x >= 0,],plot(x,Gencam_precipseas_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Gencam'),2]
                                                                     ,xlab="Precipitation Seasonality (CV)",ylab="Logistic output (probability of pressence)"))
with(Gencam_precipseas_respdf[Gencam_precipseas_respdf$x >= 0,],lines(x,Gencam_precipseas_meanresp+Gencam_precipseas_sdresp,type="l",lty=2))
with(Gencam_precipseas_respdf[Gencam_precipseas_respdf$x >= 0,],lines(x,Gencam_precipseas_meanresp-Gencam_precipseas_sdresp,type="l",lty=2))

#Land cover (factor)
Gencam_landcover0<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_0_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Gencam_landcover1<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_1_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Gencam_landcover2<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_2_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Gencam_landcover3<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_3_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Gencam_landcover4<-read.table("S:\\SDMs\\SDM_maxent\\Gencam\\plots\\species_4_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Gencam_landcover_meanresp<-apply(cbind(Gencam_landcover0$y,Gencam_landcover1$y,Gencam_landcover2$y,Gencam_landcover3$y,Gencam_landcover4$y),1,mean)
Gencam_landcover_sdresp<-apply(cbind(Gencam_landcover0$y,Gencam_landcover1$y,Gencam_landcover2$y,Gencam_landcover3$y,Gencam_landcover4$y),1,sd)
Gencam_landcover_respdf<-data.frame(cbind(Gencam_landcover0$x,Gencam_landcover_meanresp,Gencam_landcover_sdresp))
names(Gencam_landcover_respdf)[1]<-"x"
par(mar=c(8,5,3,1))
b1<-with(Gencam_landcover_respdf,barplot(Gencam_landcover_meanresp,ylim=c(0,1),main=namedf[which(namedf[,1]=='Gencam'),2],las=2
                                         ,names.arg=c("Built up","Agricultural","Forest","Natural vegetation \n (not forest)","Wetland","Ice/Snow","Freshwater")
                                         ,xlab="",ylab="Logistic output (probability of pressence)"))
with(Gencam_landcover_respdf,arrows(b1,Gencam_landcover_meanresp+Gencam_landcover_sdresp,b1,Gencam_landcover_meanresp-Gencam_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))

#Kobsim
#MAP
Kobsim_map0<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_0_AnnualPrecip.dat",sep=",",header=T)
Kobsim_map1<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_1_AnnualPrecip.dat",sep=",",header=T)
Kobsim_map2<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_2_AnnualPrecip.dat",sep=",",header=T)
Kobsim_map3<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_3_AnnualPrecip.dat",sep=",",header=T)
Kobsim_map4<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_4_AnnualPrecip.dat",sep=",",header=T)
Kobsim_map_meanresp<-apply(cbind(Kobsim_map0$y,Kobsim_map1$y,Kobsim_map2$y,Kobsim_map3$y,Kobsim_map4$y),1,mean)
Kobsim_map_sdresp<-apply(cbind(Kobsim_map0$y,Kobsim_map1$y,Kobsim_map2$y,Kobsim_map3$y,Kobsim_map4$y),1,sd)
Kobsim_map_respdf<-data.frame(cbind(Kobsim_map0$x,Kobsim_map_meanresp,Kobsim_map_sdresp))
names(Kobsim_map_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Kobsim_map_respdf[Kobsim_map_respdf$x >= 0,],plot(x,Kobsim_map_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Kobsim'),2]
                                                       ,xlab="Mean annual precipitation (mm)",ylab="Logistic output (probability of pressence)"))
with(Kobsim_map_respdf[Kobsim_map_respdf$x >= 0,],lines(x,Kobsim_map_meanresp+Kobsim_map_sdresp,type="l",lty=2))
with(Kobsim_map_respdf[Kobsim_map_respdf$x >= 0,],lines(x,Kobsim_map_meanresp-Kobsim_map_sdresp,type="l",lty=2))

#MST
Kobsim_mst0<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_0_MeanTempWarmQuart.dat",sep=",",header=T)
Kobsim_mst1<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_1_MeanTempWarmQuart.dat",sep=",",header=T)
Kobsim_mst2<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_2_MeanTempWarmQuart.dat",sep=",",header=T)
Kobsim_mst3<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_3_MeanTempWarmQuart.dat",sep=",",header=T)
Kobsim_mst4<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_4_MeanTempWarmQuart.dat",sep=",",header=T)
Kobsim_mst_meanresp<-apply(cbind(Kobsim_mst0$y,Kobsim_mst1$y,Kobsim_mst2$y,Kobsim_mst3$y,Kobsim_mst4$y),1,mean)
Kobsim_mst_sdresp<-apply(cbind(Kobsim_mst0$y,Kobsim_mst1$y,Kobsim_mst2$y,Kobsim_mst3$y,Kobsim_mst4$y),1,sd)
Kobsim_mst_respdf<-data.frame(cbind(Kobsim_mst0$x,Kobsim_mst_meanresp,Kobsim_mst_sdresp))
names(Kobsim_mst_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Kobsim_mst_respdf[Kobsim_mst_respdf$x >= 0,],plot(x/10,Kobsim_mst_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Kobsim'),2]
                                                       ,xlab=expression(Temperature~(degree~C)) ,ylab="Logistic output (probability of pressence)"))
with(Kobsim_mst_respdf[Kobsim_mst_respdf$x >= 0,],lines(x/10,Kobsim_mst_meanresp+Kobsim_mst_sdresp,type="l",lty=2))
with(Kobsim_mst_respdf[Kobsim_mst_respdf$x >= 0,],lines(x/10,Kobsim_mst_meanresp-Kobsim_mst_sdresp,type="l",lty=2))

#Tundra herbivores(continuous)
Kobsim_tundraherb0<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_0_TundraHerbivores.dat",sep=",",header=T)
Kobsim_tundraherb1<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_1_TundraHerbivores.dat",sep=",",header=T)
Kobsim_tundraherb2<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_2_TundraHerbivores.dat",sep=",",header=T)
Kobsim_tundraherb3<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_3_TundraHerbivores.dat",sep=",",header=T)
Kobsim_tundraherb4<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_4_TundraHerbivores.dat",sep=",",header=T)
Kobsim_tundraherb_meanresp<-apply(cbind(Kobsim_tundraherb0$y,Kobsim_tundraherb1$y,Kobsim_tundraherb2$y,Kobsim_tundraherb3$y,Kobsim_tundraherb4$y),1,mean)
Kobsim_tundraherb_sdresp<-apply(cbind(Kobsim_tundraherb0$y,Kobsim_tundraherb1$y,Kobsim_tundraherb2$y,Kobsim_tundraherb3$y,Kobsim_tundraherb4$y),1,sd)
Kobsim_tundraherb_respdf<-data.frame(cbind(Kobsim_tundraherb0$x,Kobsim_tundraherb_meanresp,Kobsim_tundraherb_sdresp))
names(Kobsim_tundraherb_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Kobsim_tundraherb_respdf[Kobsim_tundraherb_respdf$x >= 0,],plot(x,Kobsim_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Kobsim'),2]
                                                                     ,xlab=expression("Herbivore density (kg km"^-2*")"),ylab="Logistic output (probability of pressence)"))
with(Kobsim_tundraherb_respdf[Kobsim_tundraherb_respdf$x >= 0,],lines(x,Kobsim_tundraherb_meanresp+Kobsim_tundraherb_sdresp,type="l",lty=2))
with(Kobsim_tundraherb_respdf[Kobsim_tundraherb_respdf$x >= 0,],lines(x,Kobsim_tundraherb_meanresp-Kobsim_tundraherb_sdresp,type="l",lty=2))

#Precipitation seasonality
Kobsim_precipseas0<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_0_PrecipSeasonality.dat",sep=",",header=T)
Kobsim_precipseas1<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_1_PrecipSeasonality.dat",sep=",",header=T)
Kobsim_precipseas2<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_2_PrecipSeasonality.dat",sep=",",header=T)
Kobsim_precipseas3<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_3_PrecipSeasonality.dat",sep=",",header=T)
Kobsim_precipseas4<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_4_PrecipSeasonality.dat",sep=",",header=T)
Kobsim_precipseas_meanresp<-apply(cbind(Kobsim_precipseas0$y,Kobsim_precipseas1$y,Kobsim_precipseas2$y,Kobsim_precipseas3$y,Kobsim_precipseas4$y),1,mean)
Kobsim_precipseas_sdresp<-apply(cbind(Kobsim_precipseas0$y,Kobsim_precipseas1$y,Kobsim_precipseas2$y,Kobsim_precipseas3$y,Kobsim_precipseas4$y),1,sd)
Kobsim_precipseas_respdf<-data.frame(cbind(Kobsim_precipseas0$x,Kobsim_precipseas_meanresp,Kobsim_precipseas_sdresp))
names(Kobsim_precipseas_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Kobsim_precipseas_respdf[Kobsim_precipseas_respdf$x >= 0,],plot(x,Kobsim_precipseas_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Kobsim'),2]
                                                                     ,xlab="Precipitation Seasonality (CV)",ylab="Logistic output (probability of pressence)"))
with(Kobsim_precipseas_respdf[Kobsim_precipseas_respdf$x >= 0,],lines(x,Kobsim_precipseas_meanresp+Kobsim_precipseas_sdresp,type="l",lty=2))
with(Kobsim_precipseas_respdf[Kobsim_precipseas_respdf$x >= 0,],lines(x,Kobsim_precipseas_meanresp-Kobsim_precipseas_sdresp,type="l",lty=2))

#Land cover (factor)
Kobsim_landcover0<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_0_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Kobsim_landcover1<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_1_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Kobsim_landcover2<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_2_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Kobsim_landcover3<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_3_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Kobsim_landcover4<-read.table("S:\\SDMs\\SDM_maxent\\Kobsim\\plots\\species_4_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Kobsim_landcover_meanresp<-apply(cbind(Kobsim_landcover0$y,Kobsim_landcover1$y,Kobsim_landcover2$y,Kobsim_landcover3$y,Kobsim_landcover4$y),1,mean)
Kobsim_landcover_sdresp<-apply(cbind(Kobsim_landcover0$y,Kobsim_landcover1$y,Kobsim_landcover2$y,Kobsim_landcover3$y,Kobsim_landcover4$y),1,sd)
Kobsim_landcover_respdf<-data.frame(cbind(Kobsim_landcover0$x,Kobsim_landcover_meanresp,Kobsim_landcover_sdresp))
names(Kobsim_landcover_respdf)[1]<-"x"
par(mar=c(8,5,3,1))
b1<-with(Kobsim_landcover_respdf,barplot(Kobsim_landcover_meanresp,ylim=c(0,1),main=namedf[which(namedf[,1]=='Kobsim'),2],las=2
                                         ,names.arg=c("Built up","Agricultural","Forest","Natural vegetation \n (not forest)","Wetland","Ice/Snow","Freshwater")
                                         ,xlab="",ylab="Logistic output (probability of pressence)"))
with(Kobsim_landcover_respdf,arrows(b1,Kobsim_landcover_meanresp+Kobsim_landcover_sdresp,b1,Kobsim_landcover_meanresp-Kobsim_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))

#Prisca
#MAP
Prisca_map0<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_0_AnnualPrecip.dat",sep=",",header=T)
Prisca_map1<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_1_AnnualPrecip.dat",sep=",",header=T)
Prisca_map2<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_2_AnnualPrecip.dat",sep=",",header=T)
Prisca_map3<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_3_AnnualPrecip.dat",sep=",",header=T)
Prisca_map4<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_4_AnnualPrecip.dat",sep=",",header=T)
Prisca_map_meanresp<-apply(cbind(Prisca_map0$y,Prisca_map1$y,Prisca_map2$y,Prisca_map3$y,Prisca_map4$y),1,mean)
Prisca_map_sdresp<-apply(cbind(Prisca_map0$y,Prisca_map1$y,Prisca_map2$y,Prisca_map3$y,Prisca_map4$y),1,sd)
Prisca_map_respdf<-data.frame(cbind(Prisca_map0$x,Prisca_map_meanresp,Prisca_map_sdresp))
names(Prisca_map_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Prisca_map_respdf[Prisca_map_respdf$x >= 0,],plot(x,Prisca_map_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Prisca'),2]
                                                       ,xlab="Mean annual precipitation (mm)",ylab="Logistic output (probability of pressence)"))
with(Prisca_map_respdf[Prisca_map_respdf$x >= 0,],lines(x,Prisca_map_meanresp+Prisca_map_sdresp,type="l",lty=2))
with(Prisca_map_respdf[Prisca_map_respdf$x >= 0,],lines(x,Prisca_map_meanresp-Prisca_map_sdresp,type="l",lty=2))

#MST
Prisca_mst0<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_0_MeanTempWarmQuart.dat",sep=",",header=T)
Prisca_mst1<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_1_MeanTempWarmQuart.dat",sep=",",header=T)
Prisca_mst2<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_2_MeanTempWarmQuart.dat",sep=",",header=T)
Prisca_mst3<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_3_MeanTempWarmQuart.dat",sep=",",header=T)
Prisca_mst4<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_4_MeanTempWarmQuart.dat",sep=",",header=T)
Prisca_mst_meanresp<-apply(cbind(Prisca_mst0$y,Prisca_mst1$y,Prisca_mst2$y,Prisca_mst3$y,Prisca_mst4$y),1,mean)
Prisca_mst_sdresp<-apply(cbind(Prisca_mst0$y,Prisca_mst1$y,Prisca_mst2$y,Prisca_mst3$y,Prisca_mst4$y),1,sd)
Prisca_mst_respdf<-data.frame(cbind(Prisca_mst0$x,Prisca_mst_meanresp,Prisca_mst_sdresp))
names(Prisca_mst_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Prisca_mst_respdf[Prisca_mst_respdf$x >= 0,],plot(x/10,Prisca_mst_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Prisca'),2]
                                                       ,xlab=expression(Temperature~(degree~C)) ,ylab="Logistic output (probability of pressence)"))
with(Prisca_mst_respdf[Prisca_mst_respdf$x >= 0,],lines(x/10,Prisca_mst_meanresp+Prisca_mst_sdresp,type="l",lty=2))
with(Prisca_mst_respdf[Prisca_mst_respdf$x >= 0,],lines(x/10,Prisca_mst_meanresp-Prisca_mst_sdresp,type="l",lty=2))

#Tundra herbivores(continuous)
Prisca_tundraherb0<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_0_TundraHerbivores.dat",sep=",",header=T)
Prisca_tundraherb1<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_1_TundraHerbivores.dat",sep=",",header=T)
Prisca_tundraherb2<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_2_TundraHerbivores.dat",sep=",",header=T)
Prisca_tundraherb3<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_3_TundraHerbivores.dat",sep=",",header=T)
Prisca_tundraherb4<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_4_TundraHerbivores.dat",sep=",",header=T)
Prisca_tundraherb_meanresp<-apply(cbind(Prisca_tundraherb0$y,Prisca_tundraherb1$y,Prisca_tundraherb2$y,Prisca_tundraherb3$y,Prisca_tundraherb4$y),1,mean)
Prisca_tundraherb_sdresp<-apply(cbind(Prisca_tundraherb0$y,Prisca_tundraherb1$y,Prisca_tundraherb2$y,Prisca_tundraherb3$y,Prisca_tundraherb4$y),1,sd)
Prisca_tundraherb_respdf<-data.frame(cbind(Prisca_tundraherb0$x,Prisca_tundraherb_meanresp,Prisca_tundraherb_sdresp))
names(Prisca_tundraherb_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Prisca_tundraherb_respdf[Prisca_tundraherb_respdf$x >= 0,],plot(x,Prisca_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Prisca'),2]
                                                                     ,xlab=expression("Herbivore density (kg km"^-2*")"),ylab="Logistic output (probability of pressence)"))
with(Prisca_tundraherb_respdf[Prisca_tundraherb_respdf$x >= 0,],lines(x,Prisca_tundraherb_meanresp+Prisca_tundraherb_sdresp,type="l",lty=2))
with(Prisca_tundraherb_respdf[Prisca_tundraherb_respdf$x >= 0,],lines(x,Prisca_tundraherb_meanresp-Prisca_tundraherb_sdresp,type="l",lty=2))

#Precipitation seasonality
Prisca_precipseas0<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_0_PrecipSeasonality.dat",sep=",",header=T)
Prisca_precipseas1<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_1_PrecipSeasonality.dat",sep=",",header=T)
Prisca_precipseas2<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_2_PrecipSeasonality.dat",sep=",",header=T)
Prisca_precipseas3<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_3_PrecipSeasonality.dat",sep=",",header=T)
Prisca_precipseas4<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_4_PrecipSeasonality.dat",sep=",",header=T)
Prisca_precipseas_meanresp<-apply(cbind(Prisca_precipseas0$y,Prisca_precipseas1$y,Prisca_precipseas2$y,Prisca_precipseas3$y,Prisca_precipseas4$y),1,mean)
Prisca_precipseas_sdresp<-apply(cbind(Prisca_precipseas0$y,Prisca_precipseas1$y,Prisca_precipseas2$y,Prisca_precipseas3$y,Prisca_precipseas4$y),1,sd)
Prisca_precipseas_respdf<-data.frame(cbind(Prisca_precipseas0$x,Prisca_precipseas_meanresp,Prisca_precipseas_sdresp))
names(Prisca_precipseas_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Prisca_precipseas_respdf[Prisca_precipseas_respdf$x >= 0,],plot(x,Prisca_precipseas_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Prisca'),2]
                                                                     ,xlab="Precipitation Seasonality (CV)",ylab="Logistic output (probability of pressence)"))
with(Prisca_precipseas_respdf[Prisca_precipseas_respdf$x >= 0,],lines(x,Prisca_precipseas_meanresp+Prisca_precipseas_sdresp,type="l",lty=2))
with(Prisca_precipseas_respdf[Prisca_precipseas_respdf$x >= 0,],lines(x,Prisca_precipseas_meanresp-Prisca_precipseas_sdresp,type="l",lty=2))

#Land cover (factor)
Prisca_landcover0<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_0_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Prisca_landcover1<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_1_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Prisca_landcover2<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_2_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Prisca_landcover3<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_3_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Prisca_landcover4<-read.table("S:\\SDMs\\SDM_maxent\\Prisca\\plots\\species_4_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Prisca_landcover_meanresp<-apply(cbind(Prisca_landcover0$y,Prisca_landcover1$y,Prisca_landcover2$y,Prisca_landcover3$y,Prisca_landcover4$y),1,mean)
Prisca_landcover_sdresp<-apply(cbind(Prisca_landcover0$y,Prisca_landcover1$y,Prisca_landcover2$y,Prisca_landcover3$y,Prisca_landcover4$y),1,sd)
Prisca_landcover_respdf<-data.frame(cbind(Prisca_landcover0$x,Prisca_landcover_meanresp,Prisca_landcover_sdresp))
names(Prisca_landcover_respdf)[1]<-"x"
par(mar=c(8,5,3,1))
b1<-with(Prisca_landcover_respdf,barplot(Prisca_landcover_meanresp,ylim=c(0,1),main=namedf[which(namedf[,1]=='Prisca'),2],las=2
                                         ,names.arg=c("Built up","Agricultural","Forest","Natural vegetation \n (not forest)","Wetland","Ice/Snow","Freshwater")
                                         ,xlab="",ylab="Logistic output (probability of pressence)"))
with(Prisca_landcover_respdf,arrows(b1,Prisca_landcover_meanresp+Prisca_landcover_sdresp,b1,Prisca_landcover_meanresp-Prisca_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))

#Psealb
#MAP
Psealb_map0<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_0_AnnualPrecip.dat",sep=",",header=T)
Psealb_map1<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_1_AnnualPrecip.dat",sep=",",header=T)
Psealb_map2<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_2_AnnualPrecip.dat",sep=",",header=T)
Psealb_map3<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_3_AnnualPrecip.dat",sep=",",header=T)
Psealb_map4<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_4_AnnualPrecip.dat",sep=",",header=T)
Psealb_map_meanresp<-apply(cbind(Psealb_map0$y,Psealb_map1$y,Psealb_map2$y,Psealb_map3$y,Psealb_map4$y),1,mean)
Psealb_map_sdresp<-apply(cbind(Psealb_map0$y,Psealb_map1$y,Psealb_map2$y,Psealb_map3$y,Psealb_map4$y),1,sd)
Psealb_map_respdf<-data.frame(cbind(Psealb_map0$x,Psealb_map_meanresp,Psealb_map_sdresp))
names(Psealb_map_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Psealb_map_respdf[Psealb_map_respdf$x >= 0,],plot(x,Psealb_map_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Psealb'),2]
                                                       ,xlab="Mean annual precipitation (mm)",ylab="Logistic output (probability of pressence)"))
with(Psealb_map_respdf[Psealb_map_respdf$x >= 0,],lines(x,Psealb_map_meanresp+Psealb_map_sdresp,type="l",lty=2))
with(Psealb_map_respdf[Psealb_map_respdf$x >= 0,],lines(x,Psealb_map_meanresp-Psealb_map_sdresp,type="l",lty=2))

#MST
Psealb_mst0<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_0_MeanTempWarmQuart.dat",sep=",",header=T)
Psealb_mst1<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_1_MeanTempWarmQuart.dat",sep=",",header=T)
Psealb_mst2<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_2_MeanTempWarmQuart.dat",sep=",",header=T)
Psealb_mst3<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_3_MeanTempWarmQuart.dat",sep=",",header=T)
Psealb_mst4<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_4_MeanTempWarmQuart.dat",sep=",",header=T)
Psealb_mst_meanresp<-apply(cbind(Psealb_mst0$y,Psealb_mst1$y,Psealb_mst2$y,Psealb_mst3$y,Psealb_mst4$y),1,mean)
Psealb_mst_sdresp<-apply(cbind(Psealb_mst0$y,Psealb_mst1$y,Psealb_mst2$y,Psealb_mst3$y,Psealb_mst4$y),1,sd)
Psealb_mst_respdf<-data.frame(cbind(Psealb_mst0$x,Psealb_mst_meanresp,Psealb_mst_sdresp))
names(Psealb_mst_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Psealb_mst_respdf[Psealb_mst_respdf$x >= 0,],plot(x/10,Psealb_mst_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Psealb'),2]
                                                       ,xlab=expression(Temperature~(degree~C)) ,ylab="Logistic output (probability of pressence)"))
with(Psealb_mst_respdf[Psealb_mst_respdf$x >= 0,],lines(x/10,Psealb_mst_meanresp+Psealb_mst_sdresp,type="l",lty=2))
with(Psealb_mst_respdf[Psealb_mst_respdf$x >= 0,],lines(x/10,Psealb_mst_meanresp-Psealb_mst_sdresp,type="l",lty=2))

#Tundra herbivores(continuous)
Psealb_tundraherb0<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_0_TundraHerbivores.dat",sep=",",header=T)
Psealb_tundraherb1<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_1_TundraHerbivores.dat",sep=",",header=T)
Psealb_tundraherb2<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_2_TundraHerbivores.dat",sep=",",header=T)
Psealb_tundraherb3<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_3_TundraHerbivores.dat",sep=",",header=T)
Psealb_tundraherb4<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_4_TundraHerbivores.dat",sep=",",header=T)
Psealb_tundraherb_meanresp<-apply(cbind(Psealb_tundraherb0$y,Psealb_tundraherb1$y,Psealb_tundraherb2$y,Psealb_tundraherb3$y,Psealb_tundraherb4$y),1,mean)
Psealb_tundraherb_sdresp<-apply(cbind(Psealb_tundraherb0$y,Psealb_tundraherb1$y,Psealb_tundraherb2$y,Psealb_tundraherb3$y,Psealb_tundraherb4$y),1,sd)
Psealb_tundraherb_respdf<-data.frame(cbind(Psealb_tundraherb0$x,Psealb_tundraherb_meanresp,Psealb_tundraherb_sdresp))
names(Psealb_tundraherb_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Psealb_tundraherb_respdf[Psealb_tundraherb_respdf$x >= 0,],plot(x,Psealb_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Psealb'),2]
                                                                     ,xlab=expression("Herbivore density (kg km"^-2*")"),ylab="Logistic output (probability of pressence)"))
with(Psealb_tundraherb_respdf[Psealb_tundraherb_respdf$x >= 0,],lines(x,Psealb_tundraherb_meanresp+Psealb_tundraherb_sdresp,type="l",lty=2))
with(Psealb_tundraherb_respdf[Psealb_tundraherb_respdf$x >= 0,],lines(x,Psealb_tundraherb_meanresp-Psealb_tundraherb_sdresp,type="l",lty=2))

#Precipitation seasonality
Psealb_precipseas0<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_0_PrecipSeasonality.dat",sep=",",header=T)
Psealb_precipseas1<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_1_PrecipSeasonality.dat",sep=",",header=T)
Psealb_precipseas2<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_2_PrecipSeasonality.dat",sep=",",header=T)
Psealb_precipseas3<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_3_PrecipSeasonality.dat",sep=",",header=T)
Psealb_precipseas4<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_4_PrecipSeasonality.dat",sep=",",header=T)
Psealb_precipseas_meanresp<-apply(cbind(Psealb_precipseas0$y,Psealb_precipseas1$y,Psealb_precipseas2$y,Psealb_precipseas3$y,Psealb_precipseas4$y),1,mean)
Psealb_precipseas_sdresp<-apply(cbind(Psealb_precipseas0$y,Psealb_precipseas1$y,Psealb_precipseas2$y,Psealb_precipseas3$y,Psealb_precipseas4$y),1,sd)
Psealb_precipseas_respdf<-data.frame(cbind(Psealb_precipseas0$x,Psealb_precipseas_meanresp,Psealb_precipseas_sdresp))
names(Psealb_precipseas_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Psealb_precipseas_respdf[Psealb_precipseas_respdf$x >= 0,],plot(x,Psealb_precipseas_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Psealb'),2]
                                                                     ,xlab="Precipitation Seasonality (CV)",ylab="Logistic output (probability of pressence)"))
with(Psealb_precipseas_respdf[Psealb_precipseas_respdf$x >= 0,],lines(x,Psealb_precipseas_meanresp+Psealb_precipseas_sdresp,type="l",lty=2))
with(Psealb_precipseas_respdf[Psealb_precipseas_respdf$x >= 0,],lines(x,Psealb_precipseas_meanresp-Psealb_precipseas_sdresp,type="l",lty=2))

#Land cover (factor)
Psealb_landcover0<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_0_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Psealb_landcover1<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_1_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Psealb_landcover2<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_2_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Psealb_landcover3<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_3_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Psealb_landcover4<-read.table("S:\\SDMs\\SDM_maxent\\Psealb\\plots\\species_4_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Psealb_landcover_meanresp<-apply(cbind(Psealb_landcover0$y,Psealb_landcover1$y,Psealb_landcover2$y,Psealb_landcover3$y,Psealb_landcover4$y),1,mean)
Psealb_landcover_sdresp<-apply(cbind(Psealb_landcover0$y,Psealb_landcover1$y,Psealb_landcover2$y,Psealb_landcover3$y,Psealb_landcover4$y),1,sd)
Psealb_landcover_respdf<-data.frame(cbind(Psealb_landcover0$x,Psealb_landcover_meanresp,Psealb_landcover_sdresp))
names(Psealb_landcover_respdf)[1]<-"x"
par(mar=c(8,5,3,1))
b1<-with(Psealb_landcover_respdf,barplot(Psealb_landcover_meanresp,ylim=c(0,1),main=namedf[which(namedf[,1]=='Psealb'),2],las=2
                                         ,names.arg=c("Built up","Agricultural","Forest","Natural vegetation \n (not forest)","Wetland","Ice/Snow","Freshwater")
                                         ,xlab="",ylab="Logistic output (probability of pressence)"))
with(Psealb_landcover_respdf,arrows(b1,Psealb_landcover_meanresp+Psealb_landcover_sdresp,b1,Psealb_landcover_meanresp-Psealb_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))


#Pulver
#MAP
Pulver_map0<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_0_AnnualPrecip.dat",sep=",",header=T)
Pulver_map1<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_1_AnnualPrecip.dat",sep=",",header=T)
Pulver_map2<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_2_AnnualPrecip.dat",sep=",",header=T)
Pulver_map3<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_3_AnnualPrecip.dat",sep=",",header=T)
Pulver_map4<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_4_AnnualPrecip.dat",sep=",",header=T)
Pulver_map_meanresp<-apply(cbind(Pulver_map0$y,Pulver_map1$y,Pulver_map2$y,Pulver_map3$y,Pulver_map4$y),1,mean)
Pulver_map_sdresp<-apply(cbind(Pulver_map0$y,Pulver_map1$y,Pulver_map2$y,Pulver_map3$y,Pulver_map4$y),1,sd)
Pulver_map_respdf<-data.frame(cbind(Pulver_map0$x,Pulver_map_meanresp,Pulver_map_sdresp))
names(Pulver_map_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Pulver_map_respdf[Pulver_map_respdf$x >= 0,],plot(x,Pulver_map_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Pulver'),2]
                                                       ,xlab="Mean annual precipitation (mm)",ylab="Logistic output (probability of pressence)"))
with(Pulver_map_respdf[Pulver_map_respdf$x >= 0,],lines(x,Pulver_map_meanresp+Pulver_map_sdresp,type="l",lty=2))
with(Pulver_map_respdf[Pulver_map_respdf$x >= 0,],lines(x,Pulver_map_meanresp-Pulver_map_sdresp,type="l",lty=2))

#MST
Pulver_mst0<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_0_MeanTempWarmQuart.dat",sep=",",header=T)
Pulver_mst1<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_1_MeanTempWarmQuart.dat",sep=",",header=T)
Pulver_mst2<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_2_MeanTempWarmQuart.dat",sep=",",header=T)
Pulver_mst3<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_3_MeanTempWarmQuart.dat",sep=",",header=T)
Pulver_mst4<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_4_MeanTempWarmQuart.dat",sep=",",header=T)
Pulver_mst_meanresp<-apply(cbind(Pulver_mst0$y,Pulver_mst1$y,Pulver_mst2$y,Pulver_mst3$y,Pulver_mst4$y),1,mean)
Pulver_mst_sdresp<-apply(cbind(Pulver_mst0$y,Pulver_mst1$y,Pulver_mst2$y,Pulver_mst3$y,Pulver_mst4$y),1,sd)
Pulver_mst_respdf<-data.frame(cbind(Pulver_mst0$x,Pulver_mst_meanresp,Pulver_mst_sdresp))
names(Pulver_mst_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Pulver_mst_respdf[Pulver_mst_respdf$x >= 0,],plot(x/10,Pulver_mst_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Pulver'),2]
                                                       ,xlab=expression(Temperature~(degree~C)) ,ylab="Logistic output (probability of pressence)"))
with(Pulver_mst_respdf[Pulver_mst_respdf$x >= 0,],lines(x/10,Pulver_mst_meanresp+Pulver_mst_sdresp,type="l",lty=2))
with(Pulver_mst_respdf[Pulver_mst_respdf$x >= 0,],lines(x/10,Pulver_mst_meanresp-Pulver_mst_sdresp,type="l",lty=2))

#Tundra herbivores(continuous)
Pulver_tundraherb0<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_0_TundraHerbivores.dat",sep=",",header=T)
Pulver_tundraherb1<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_1_TundraHerbivores.dat",sep=",",header=T)
Pulver_tundraherb2<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_2_TundraHerbivores.dat",sep=",",header=T)
Pulver_tundraherb3<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_3_TundraHerbivores.dat",sep=",",header=T)
Pulver_tundraherb4<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_4_TundraHerbivores.dat",sep=",",header=T)
Pulver_tundraherb_meanresp<-apply(cbind(Pulver_tundraherb0$y,Pulver_tundraherb1$y,Pulver_tundraherb2$y,Pulver_tundraherb3$y,Pulver_tundraherb4$y),1,mean)
Pulver_tundraherb_sdresp<-apply(cbind(Pulver_tundraherb0$y,Pulver_tundraherb1$y,Pulver_tundraherb2$y,Pulver_tundraherb3$y,Pulver_tundraherb4$y),1,sd)
Pulver_tundraherb_respdf<-data.frame(cbind(Pulver_tundraherb0$x,Pulver_tundraherb_meanresp,Pulver_tundraherb_sdresp))
names(Pulver_tundraherb_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Pulver_tundraherb_respdf[Pulver_tundraherb_respdf$x >= 0,],plot(x,Pulver_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Pulver'),2]
                                                                     ,xlab=expression("Herbivore density (kg km"^-2*")"),ylab="Logistic output (probability of pressence)"))
with(Pulver_tundraherb_respdf[Pulver_tundraherb_respdf$x >= 0,],lines(x,Pulver_tundraherb_meanresp+Pulver_tundraherb_sdresp,type="l",lty=2))
with(Pulver_tundraherb_respdf[Pulver_tundraherb_respdf$x >= 0,],lines(x,Pulver_tundraherb_meanresp-Pulver_tundraherb_sdresp,type="l",lty=2))

#Precipitation seasonality
Pulver_precipseas0<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_0_PrecipSeasonality.dat",sep=",",header=T)
Pulver_precipseas1<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_1_PrecipSeasonality.dat",sep=",",header=T)
Pulver_precipseas2<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_2_PrecipSeasonality.dat",sep=",",header=T)
Pulver_precipseas3<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_3_PrecipSeasonality.dat",sep=",",header=T)
Pulver_precipseas4<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_4_PrecipSeasonality.dat",sep=",",header=T)
Pulver_precipseas_meanresp<-apply(cbind(Pulver_precipseas0$y,Pulver_precipseas1$y,Pulver_precipseas2$y,Pulver_precipseas3$y,Pulver_precipseas4$y),1,mean)
Pulver_precipseas_sdresp<-apply(cbind(Pulver_precipseas0$y,Pulver_precipseas1$y,Pulver_precipseas2$y,Pulver_precipseas3$y,Pulver_precipseas4$y),1,sd)
Pulver_precipseas_respdf<-data.frame(cbind(Pulver_precipseas0$x,Pulver_precipseas_meanresp,Pulver_precipseas_sdresp))
names(Pulver_precipseas_respdf)[1]<-"x"
par(mar=c(5,5,3,1))
with(Pulver_precipseas_respdf[Pulver_precipseas_respdf$x >= 0,],plot(x,Pulver_precipseas_meanresp,type="l",ylim=c(0,1),las=1,main=namedf[which(namedf[,1]=='Pulver'),2]
                                                                     ,xlab="Precipitation Seasonality (CV)",ylab="Logistic output (probability of pressence)"))
with(Pulver_precipseas_respdf[Pulver_precipseas_respdf$x >= 0,],lines(x,Pulver_precipseas_meanresp+Pulver_precipseas_sdresp,type="l",lty=2))
with(Pulver_precipseas_respdf[Pulver_precipseas_respdf$x >= 0,],lines(x,Pulver_precipseas_meanresp-Pulver_precipseas_sdresp,type="l",lty=2))

#Land cover (factor)
Pulver_landcover0<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_0_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Pulver_landcover1<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_1_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Pulver_landcover2<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_2_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Pulver_landcover3<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_3_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Pulver_landcover4<-read.table("S:\\SDMs\\SDM_maxent\\Pulver\\plots\\species_4_AR50TYPE100m_ProjectRaster1.dat",sep=",",header=T)
Pulver_landcover_meanresp<-apply(cbind(Pulver_landcover0$y,Pulver_landcover1$y,Pulver_landcover2$y,Pulver_landcover3$y,Pulver_landcover4$y),1,mean)
Pulver_landcover_sdresp<-apply(cbind(Pulver_landcover0$y,Pulver_landcover1$y,Pulver_landcover2$y,Pulver_landcover3$y,Pulver_landcover4$y),1,sd)
Pulver_landcover_respdf<-data.frame(cbind(Pulver_landcover0$x,Pulver_landcover_meanresp,Pulver_landcover_sdresp))
names(Pulver_landcover_respdf)[1]<-"x"
par(mar=c(8,5,3,1))
b1<-with(Pulver_landcover_respdf,barplot(Pulver_landcover_meanresp,ylim=c(0,1),main=namedf[which(namedf[,1]=='Pulver'),2],las=2
                                         ,names.arg=c("Built up","Agricultural","Forest","Natural vegetation \n (not forest)","Wetland","Ice/Snow","Freshwater")
                                         ,xlab="",ylab="Logistic output (probability of pressence)"))
with(Pulver_landcover_respdf,arrows(b1,Pulver_landcover_meanresp+Pulver_landcover_sdresp,b1,Pulver_landcover_meanresp-Pulver_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))



######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
#Alternative plot

x11(8,11,pointsize=6)
par(mfrow=c(7,5))
par(mar=c(0,0,0,0),oma=c(6,6,1,1))
par(tcl = -0.25)
par(xpd = NA)
#Botlan
with(Botlan_map_respdf[Botlan_map_respdf$x >= 0,],plot(x,Botlan_map_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
axis(2)
title(ylab="Logistic output (probability of pressence)")
box(col = "grey60")
with(Botlan_map_respdf[Botlan_map_respdf$x >= 0,],lines(x,Botlan_map_meanresp+Botlan_map_sdresp,type="l",lty=2))
with(Botlan_map_respdf[Botlan_map_respdf$x >= 0,],lines(x,Botlan_map_meanresp-Botlan_map_sdresp,type="l",lty=2))
text(1200,1,namedf[which(namedf[,1]=='Botlan'),2],cex=1.2)

with(Botlan_mst_respdf[Botlan_mst_respdf$x >= 0,],plot(x/10,Botlan_mst_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
with(Botlan_mst_respdf[Botlan_mst_respdf$x >= 0,],lines(x/10,Botlan_mst_meanresp+Botlan_mst_sdresp,type="l",lty=2))
with(Botlan_mst_respdf[Botlan_mst_respdf$x >= 0,],lines(x/10,Botlan_mst_meanresp-Botlan_mst_sdresp,type="l",lty=2))
box(col = "grey60")

with(Botlan_tundraherb_respdf[Botlan_tundraherb_respdf$x >= 0,],plot(x,Botlan_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Botlan_tundraherb_respdf[Botlan_tundraherb_respdf$x >= 0,],lines(x,Botlan_tundraherb_meanresp+Botlan_tundraherb_sdresp,type="l",lty=2))
with(Botlan_tundraherb_respdf[Botlan_tundraherb_respdf$x >= 0,],lines(x,Botlan_tundraherb_meanresp-Botlan_tundraherb_sdresp,type="l",lty=2))

with(Botlan_precipseas_respdf[Botlan_precipseas_respdf$x >= 0,],plot(x,Botlan_precipseas_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Botlan_precipseas_respdf[Botlan_precipseas_respdf$x >= 0,],lines(x,Botlan_precipseas_meanresp+Botlan_precipseas_sdresp,type="l",lty=2))
with(Botlan_precipseas_respdf[Botlan_precipseas_respdf$x >= 0,],lines(x,Botlan_precipseas_meanresp-Botlan_precipseas_sdresp,type="l",lty=2))

b1<-with(Botlan_landcover_respdf,barplot(Botlan_landcover_meanresp,ylim=c(0,1),las=2,axes=F,xlab="",ylab=""))
text(b1,0.2,c("Built up","Agricultural","Forest","Natural vegetation","Wetland","Ice/Snow","Freshwater"),srt=90)
with(Botlan_landcover_respdf,arrows(b1,Botlan_landcover_meanresp+Botlan_landcover_sdresp,b1,Botlan_landcover_meanresp-Botlan_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))
box(col = "grey60")
#Comten
with(Comten_map_respdf[Comten_map_respdf$x >= 0,],plot(x,Comten_map_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
axis(2)
title(ylab="Logistic output (probability of pressence)")
box(col = "grey60")
with(Comten_map_respdf[Comten_map_respdf$x >= 0,],lines(x,Comten_map_meanresp+Comten_map_sdresp,type="l",lty=2))
with(Comten_map_respdf[Comten_map_respdf$x >= 0,],lines(x,Comten_map_meanresp-Comten_map_sdresp,type="l",lty=2))
text(1200,1,namedf[which(namedf[,1]=='Comten'),2],cex=1.2)

with(Comten_mst_respdf[Comten_mst_respdf$x >= 0,],plot(x/10,Comten_mst_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
with(Comten_mst_respdf[Comten_mst_respdf$x >= 0,],lines(x/10,Comten_mst_meanresp+Comten_mst_sdresp,type="l",lty=2))
with(Comten_mst_respdf[Comten_mst_respdf$x >= 0,],lines(x/10,Comten_mst_meanresp-Comten_mst_sdresp,type="l",lty=2))
box(col = "grey60")

with(Comten_tundraherb_respdf[Comten_tundraherb_respdf$x >= 0,],plot(x,Comten_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Comten_tundraherb_respdf[Comten_tundraherb_respdf$x >= 0,],lines(x,Comten_tundraherb_meanresp+Comten_tundraherb_sdresp,type="l",lty=2))
with(Comten_tundraherb_respdf[Comten_tundraherb_respdf$x >= 0,],lines(x,Comten_tundraherb_meanresp-Comten_tundraherb_sdresp,type="l",lty=2))

with(Comten_precipseas_respdf[Comten_precipseas_respdf$x >= 0,],plot(x,Comten_precipseas_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Comten_precipseas_respdf[Comten_precipseas_respdf$x >= 0,],lines(x,Comten_precipseas_meanresp+Comten_precipseas_sdresp,type="l",lty=2))
with(Comten_precipseas_respdf[Comten_precipseas_respdf$x >= 0,],lines(x,Comten_precipseas_meanresp-Comten_precipseas_sdresp,type="l",lty=2))

b1<-with(Comten_landcover_respdf,barplot(Comten_landcover_meanresp,ylim=c(0,1),las=2,axes=F,xlab="",ylab=""))
text(b1,0.2,c("Built up","Agricultural","Forest","Natural vegetation","Wetland","Ice/Snow","Freshwater"),srt=90)
with(Comten_landcover_respdf,arrows(b1,Comten_landcover_meanresp+Comten_landcover_sdresp,b1,Comten_landcover_meanresp-Comten_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))
box(col = "grey60")
#Gencam
with(Gencam_map_respdf[Gencam_map_respdf$x >= 0,],plot(x,Gencam_map_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
axis(2)
title(ylab="Logistic output (probability of pressence)")
box(col = "grey60")
with(Gencam_map_respdf[Gencam_map_respdf$x >= 0,],lines(x,Gencam_map_meanresp+Gencam_map_sdresp,type="l",lty=2))
with(Gencam_map_respdf[Gencam_map_respdf$x >= 0,],lines(x,Gencam_map_meanresp-Gencam_map_sdresp,type="l",lty=2))
text(1200,1,namedf[which(namedf[,1]=='Gencam'),2],cex=1.2)

with(Gencam_mst_respdf[Gencam_mst_respdf$x >= 0,],plot(x/10,Gencam_mst_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
with(Gencam_mst_respdf[Gencam_mst_respdf$x >= 0,],lines(x/10,Gencam_mst_meanresp+Gencam_mst_sdresp,type="l",lty=2))
with(Gencam_mst_respdf[Gencam_mst_respdf$x >= 0,],lines(x/10,Gencam_mst_meanresp-Gencam_mst_sdresp,type="l",lty=2))
box(col = "grey60")

with(Gencam_tundraherb_respdf[Gencam_tundraherb_respdf$x >= 0,],plot(x,Gencam_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Gencam_tundraherb_respdf[Gencam_tundraherb_respdf$x >= 0,],lines(x,Gencam_tundraherb_meanresp+Gencam_tundraherb_sdresp,type="l",lty=2))
with(Gencam_tundraherb_respdf[Gencam_tundraherb_respdf$x >= 0,],lines(x,Gencam_tundraherb_meanresp-Gencam_tundraherb_sdresp,type="l",lty=2))

with(Gencam_precipseas_respdf[Gencam_precipseas_respdf$x >= 0,],plot(x,Gencam_precipseas_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Gencam_precipseas_respdf[Gencam_precipseas_respdf$x >= 0,],lines(x,Gencam_precipseas_meanresp+Gencam_precipseas_sdresp,type="l",lty=2))
with(Gencam_precipseas_respdf[Gencam_precipseas_respdf$x >= 0,],lines(x,Gencam_precipseas_meanresp-Gencam_precipseas_sdresp,type="l",lty=2))

b1<-with(Gencam_landcover_respdf,barplot(Gencam_landcover_meanresp,ylim=c(0,1),las=2,axes=F,xlab="",ylab=""))
text(b1,0.2,c("Built up","Agricultural","Forest","Natural vegetation","Wetland","Ice/Snow","Freshwater"),srt=90)
with(Gencam_landcover_respdf,arrows(b1,Gencam_landcover_meanresp+Gencam_landcover_sdresp,b1,Gencam_landcover_meanresp-Gencam_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))
box(col = "grey60")
#Kobsim
with(Kobsim_map_respdf[Kobsim_map_respdf$x >= 0,],plot(x,Kobsim_map_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
axis(2)
title(ylab="Logistic output (probability of pressence)")
box(col = "grey60")
with(Kobsim_map_respdf[Kobsim_map_respdf$x >= 0,],lines(x,Kobsim_map_meanresp+Kobsim_map_sdresp,type="l",lty=2))
with(Kobsim_map_respdf[Kobsim_map_respdf$x >= 0,],lines(x,Kobsim_map_meanresp-Kobsim_map_sdresp,type="l",lty=2))
text(1200,1,namedf[which(namedf[,1]=='Kobsim'),2],cex=1.2)

with(Kobsim_mst_respdf[Kobsim_mst_respdf$x >= 0,],plot(x/10,Kobsim_mst_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
with(Kobsim_mst_respdf[Kobsim_mst_respdf$x >= 0,],lines(x/10,Kobsim_mst_meanresp+Kobsim_mst_sdresp,type="l",lty=2))
with(Kobsim_mst_respdf[Kobsim_mst_respdf$x >= 0,],lines(x/10,Kobsim_mst_meanresp-Kobsim_mst_sdresp,type="l",lty=2))
box(col = "grey60")

with(Kobsim_tundraherb_respdf[Kobsim_tundraherb_respdf$x >= 0,],plot(x,Kobsim_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Kobsim_tundraherb_respdf[Kobsim_tundraherb_respdf$x >= 0,],lines(x,Kobsim_tundraherb_meanresp+Kobsim_tundraherb_sdresp,type="l",lty=2))
with(Kobsim_tundraherb_respdf[Kobsim_tundraherb_respdf$x >= 0,],lines(x,Kobsim_tundraherb_meanresp-Kobsim_tundraherb_sdresp,type="l",lty=2))

with(Kobsim_precipseas_respdf[Kobsim_precipseas_respdf$x >= 0,],plot(x,Kobsim_precipseas_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Kobsim_precipseas_respdf[Kobsim_precipseas_respdf$x >= 0,],lines(x,Kobsim_precipseas_meanresp+Kobsim_precipseas_sdresp,type="l",lty=2))
with(Kobsim_precipseas_respdf[Kobsim_precipseas_respdf$x >= 0,],lines(x,Kobsim_precipseas_meanresp-Kobsim_precipseas_sdresp,type="l",lty=2))

b1<-with(Kobsim_landcover_respdf,barplot(Kobsim_landcover_meanresp,ylim=c(0,1),las=2,axes=F,xlab="",ylab=""))
text(b1,0.2,c("Built up","Agricultural","Forest","Natural vegetation","Wetland","Ice/Snow","Freshwater"),srt=90)
with(Kobsim_landcover_respdf,arrows(b1,Kobsim_landcover_meanresp+Kobsim_landcover_sdresp,b1,Kobsim_landcover_meanresp-Kobsim_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))
box(col = "grey60")
#Prisca
with(Prisca_map_respdf[Prisca_map_respdf$x >= 0,],plot(x,Prisca_map_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
axis(2)
title(ylab="Logistic output (probability of pressence)")
box(col = "grey60")
with(Prisca_map_respdf[Prisca_map_respdf$x >= 0,],lines(x,Prisca_map_meanresp+Prisca_map_sdresp,type="l",lty=2))
with(Prisca_map_respdf[Prisca_map_respdf$x >= 0,],lines(x,Prisca_map_meanresp-Prisca_map_sdresp,type="l",lty=2))
text(1200,1,namedf[which(namedf[,1]=='Prisca'),2],cex=1.2)

with(Prisca_mst_respdf[Prisca_mst_respdf$x >= 0,],plot(x/10,Prisca_mst_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
with(Prisca_mst_respdf[Prisca_mst_respdf$x >= 0,],lines(x/10,Prisca_mst_meanresp+Prisca_mst_sdresp,type="l",lty=2))
with(Prisca_mst_respdf[Prisca_mst_respdf$x >= 0,],lines(x/10,Prisca_mst_meanresp-Prisca_mst_sdresp,type="l",lty=2))
box(col = "grey60")

with(Prisca_tundraherb_respdf[Prisca_tundraherb_respdf$x >= 0,],plot(x,Prisca_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Prisca_tundraherb_respdf[Prisca_tundraherb_respdf$x >= 0,],lines(x,Prisca_tundraherb_meanresp+Prisca_tundraherb_sdresp,type="l",lty=2))
with(Prisca_tundraherb_respdf[Prisca_tundraherb_respdf$x >= 0,],lines(x,Prisca_tundraherb_meanresp-Prisca_tundraherb_sdresp,type="l",lty=2))

with(Prisca_precipseas_respdf[Prisca_precipseas_respdf$x >= 0,],plot(x,Prisca_precipseas_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Prisca_precipseas_respdf[Prisca_precipseas_respdf$x >= 0,],lines(x,Prisca_precipseas_meanresp+Prisca_precipseas_sdresp,type="l",lty=2))
with(Prisca_precipseas_respdf[Prisca_precipseas_respdf$x >= 0,],lines(x,Prisca_precipseas_meanresp-Prisca_precipseas_sdresp,type="l",lty=2))

b1<-with(Prisca_landcover_respdf,barplot(Prisca_landcover_meanresp,ylim=c(0,1),las=2,axes=F,xlab="",ylab=""))
text(b1,0.2,c("Built up","Agricultural","Forest","Natural vegetation","Wetland","Ice/Snow","Freshwater"),srt=90)
with(Prisca_landcover_respdf,arrows(b1,Prisca_landcover_meanresp+Prisca_landcover_sdresp,b1,Prisca_landcover_meanresp-Prisca_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))
box(col = "grey60")
#Psealb
with(Psealb_map_respdf[Psealb_map_respdf$x >= 0,],plot(x,Psealb_map_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
axis(2)
title(ylab="Logistic output (probability of pressence)")
box(col = "grey60")
with(Psealb_map_respdf[Psealb_map_respdf$x >= 0,],lines(x,Psealb_map_meanresp+Psealb_map_sdresp,type="l",lty=2))
with(Psealb_map_respdf[Psealb_map_respdf$x >= 0,],lines(x,Psealb_map_meanresp-Psealb_map_sdresp,type="l",lty=2))
text(1200,1,namedf[which(namedf[,1]=='Psealb'),2],cex=1.2)

with(Psealb_mst_respdf[Psealb_mst_respdf$x >= 0,],plot(x/10,Psealb_mst_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
with(Psealb_mst_respdf[Psealb_mst_respdf$x >= 0,],lines(x/10,Psealb_mst_meanresp+Psealb_mst_sdresp,type="l",lty=2))
with(Psealb_mst_respdf[Psealb_mst_respdf$x >= 0,],lines(x/10,Psealb_mst_meanresp-Psealb_mst_sdresp,type="l",lty=2))
box(col = "grey60")

with(Psealb_tundraherb_respdf[Psealb_tundraherb_respdf$x >= 0,],plot(x,Psealb_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Psealb_tundraherb_respdf[Psealb_tundraherb_respdf$x >= 0,],lines(x,Psealb_tundraherb_meanresp+Psealb_tundraherb_sdresp,type="l",lty=2))
with(Psealb_tundraherb_respdf[Psealb_tundraherb_respdf$x >= 0,],lines(x,Psealb_tundraherb_meanresp-Psealb_tundraherb_sdresp,type="l",lty=2))

with(Psealb_precipseas_respdf[Psealb_precipseas_respdf$x >= 0,],plot(x,Psealb_precipseas_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Psealb_precipseas_respdf[Psealb_precipseas_respdf$x >= 0,],lines(x,Psealb_precipseas_meanresp+Psealb_precipseas_sdresp,type="l",lty=2))
with(Psealb_precipseas_respdf[Psealb_precipseas_respdf$x >= 0,],lines(x,Psealb_precipseas_meanresp-Psealb_precipseas_sdresp,type="l",lty=2))

b1<-with(Psealb_landcover_respdf,barplot(Psealb_landcover_meanresp,ylim=c(0,1),las=2,axes=F,xlab="",ylab=""))
text(b1,0.2,c("Built up","Agricultural","Forest","Natural vegetation","Wetland","Ice/Snow","Freshwater"),srt=90)
with(Psealb_landcover_respdf,arrows(b1,Psealb_landcover_meanresp+Psealb_landcover_sdresp,b1,Psealb_landcover_meanresp-Psealb_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))
box(col = "grey60")
#Pulver
with(Pulver_map_respdf[Pulver_map_respdf$x >= 0,],plot(x,Pulver_map_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
axis(2)
title(ylab="Logistic output (probability of pressence)")
axis(1)
title(xlab="Annual Precipitation (mm)")
box(col = "grey60")
with(Pulver_map_respdf[Pulver_map_respdf$x >= 0,],lines(x,Pulver_map_meanresp+Pulver_map_sdresp,type="l",lty=2))
with(Pulver_map_respdf[Pulver_map_respdf$x >= 0,],lines(x,Pulver_map_meanresp-Pulver_map_sdresp,type="l",lty=2))
text(1200,1,namedf[which(namedf[,1]=='Pulver'),2],cex=1.2)

with(Pulver_mst_respdf[Pulver_mst_respdf$x >= 0,],plot(x/10,Pulver_mst_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
with(Pulver_mst_respdf[Pulver_mst_respdf$x >= 0,],lines(x/10,Pulver_mst_meanresp+Pulver_mst_sdresp,type="l",lty=2))
with(Pulver_mst_respdf[Pulver_mst_respdf$x >= 0,],lines(x/10,Pulver_mst_meanresp-Pulver_mst_sdresp,type="l",lty=2))
axis(1)
title(xlab=expression("Mean summer temperature"~(degree~C)))
box(col = "grey60")

with(Pulver_tundraherb_respdf[Pulver_tundraherb_respdf$x >= 0,],plot(x,Pulver_tundraherb_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
axis(1)
title(xlab=expression("Herbivore density (kg km"^-2*")"))
with(Pulver_tundraherb_respdf[Pulver_tundraherb_respdf$x >= 0,],lines(x,Pulver_tundraherb_meanresp+Pulver_tundraherb_sdresp,type="l",lty=2))
with(Pulver_tundraherb_respdf[Pulver_tundraherb_respdf$x >= 0,],lines(x,Pulver_tundraherb_meanresp-Pulver_tundraherb_sdresp,type="l",lty=2))

with(Pulver_precipseas_respdf[Pulver_precipseas_respdf$x >= 0,],plot(x,Pulver_precipseas_meanresp,type="l",ylim=c(0,1),las=1,axes=F,xlab="",ylab=""))
box(col = "grey60")
with(Pulver_precipseas_respdf[Pulver_precipseas_respdf$x >= 0,],lines(x,Pulver_precipseas_meanresp+Pulver_precipseas_sdresp,type="l",lty=2))
with(Pulver_precipseas_respdf[Pulver_precipseas_respdf$x >= 0,],lines(x,Pulver_precipseas_meanresp-Pulver_precipseas_sdresp,type="l",lty=2))
axis(1)
title(xlab="Precipitation seasonality (CV)")

b1<-with(Pulver_landcover_respdf,barplot(Pulver_landcover_meanresp,ylim=c(0,1),las=2,axes=F,xlab="",ylab=""))
text(b1,0.2,c("Built up","Agricultural","Forest","Natural vegetation","Wetland","Ice/Snow","Freshwater"),srt=90)
with(Pulver_landcover_respdf,arrows(b1,Pulver_landcover_meanresp+Pulver_landcover_sdresp,b1,Pulver_landcover_meanresp-Pulver_landcover_sdresp
                                    ,angle=90,code=3,length=0.05))
title(xlab="Land cover")
box(col = "grey60")
