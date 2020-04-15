#Forest RedList species
require(rgbif)#Spatial data
require(raster)#Spatial data
require(rgdal)#Spatial data
require(sp)#Spatial data
require(data.table)#Reading in big tables
require(rasterVis)#Visualising spatial data
require(dismo)#Distribution modelling
require(ENMeval)#Tuning maxent
require(usdm)#VIFs

#Norway polygon
norway<-raster::getData('GADM', country='NOR', level=0)
norwayP<-sp::spTransform(norway,"+proj=utm +zone=32")
plot(norwayP)
plot(norway)


# # #Redlist -----------------------------------------------------------------
# # #Make species lists with those linked to grazing changes
# # #Data from Artsdatabanken Rødlist 2015. All redlisted spp, påvirkning factor: changed habitat use
# # #Only species (not subspecies)
# # #Filter by 'beite' 
# vasc<-read.table('Rodlista2015_vascplant.csv'
#                   ,header=T,sep=";",fileEncoding="UTF-16LE")
#  
# herbvasc<-droplevels(vasc[grep('*beit*',vasc$Påvirkningsfaktorer,ignore.case=T),])
# dim(herbvasc)
# herbvasc$Vitenskapelig.navn
# 
# moss<-read.table('Rodlista2015_mosses.csv'
#                  ,header=T,sep=";",fileEncoding="UTF-16LE")
# 
# herbmoss<-droplevels(moss[grep('*beit*',moss$Påvirkningsfaktorer,ignore.case=T),])
# dim(herbmoss)
# herbmoss$Vitenskapelig.navn
# 
# lichen<-read.table('Rodlista2015_lichens.csv'
#                  ,header=T,sep=";",fileEncoding="UTF-16LE")
# 
# herblichen<-droplevels(lichen[grep('*beit*',lichen$Påvirkningsfaktorer,ignore.case=T),])
# dim(herblichen)
# herblichen$Vitenskapelig.navn
# 
# 
# redlistsp_all<-rbind(vasc,moss,lichen)
# redlistsp_all$PlantGroup<-c(rep('Vascular',times=nrow(vasc)),rep('Bryphyte',times=nrow(moss)),rep('Lichen',times=nrow(lichen)))
# 
# View(redlistsp_all)
# 
# 
# # GBIF import -------------------------------------------------------------
# 
# #Keyed download with doi linkk
# keysL<-sapply(as.character(herblichen$Vitenskapelig.navn),function(x) name_backbone(x,rank='species')$speciesKey)
# paste(keysL,collapse=',')#Copy and paste to below
# odLichen<-occ_download('taxonKey = 2609133,2607387,3408876,2603661,2602561,2609349,3422592,3433876,2601243,2608140,5260765,7083298,2605872,7425899,5260770,9198278,5258394,8704544,3397573,3429282,3390552,2609409,2609180,5260761,2607725,8335421,5260565,2599762,3429909,3389602'
#                        ,'country = NO','hasCoordinate = TRUE',
#                        user='jamesspeed',pwd='*****',email='*****')
# occ_download_meta(odLichen)
# gbif_citation(occ_download_meta(odLichen))# GBIF Occurrence Download https://doi.org/10.15468/dl.pl144i Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2018-11-26"
# 
# keysM<-sapply(as.character(herbmoss$Vitenskapelig.navn),function(x) name_backbone(x,rank='species')$speciesKey)
# paste(keysM,collapse=',')
# odMoss<-occ_download('taxonKey = 5280585,4276928,5283214,7800507,2689413'
#                        ,'country = NO','hasCoordinate = TRUE',
#                        user='jamesspeed',pwd='*****',email='*****')
# occ_download_meta(odMoss)
# gbif_citation(occ_download_meta(odMoss))# "GBIF Occurrence Download https://doi.org/10.15468/dl.ydfbtt Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2018-11-26"
# 
# keysV<-sapply(as.character(herbvasc$Vitenskapelig.navn),function(x) name_backbone(x,rank='species')$speciesKey)
# paste(keysV,collapse=',')
# odVasc<-occ_download('taxonKey = 9139754,8558004,2704845,3025813,7931051,2792588,9020552,5405976,2849252,8231647,3001509,8917443,2820517,3012376,2975152,8277403,3033129,5361866,5284517,8869754,7270427,5410857,8915737,2926086,5409958,2787993,3111049,9177060,2927078,3012509,2926055,7660935,2914396,2975380,2925944'
#                      ,'country = NO','hasCoordinate = TRUE',
#                      user='jamesspeed',pwd='*****',email='*****')
# occ_download_meta(odVasc)
# gbif_citation(occ_download_meta(odVasc))# "GBIF Occurrence Download https://doi.org/10.15468/dl.7eqijw Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2018-11-26"
# 
# odLichdat<-occ_download_get(odLichen,overwrite=T)
# odMossdat<-occ_download_get(odMoss,overwrite=T)
# odVascdat<-occ_download_get(odVasc,overwrite=T)
# 
# odLichendata<-occ_download_import(odLichdat)
# odMossdata<-occ_download_import(odMossdat)
# odVascdata<-occ_download_import(odVascdat)
# 
# summary(as.factor(odVascdata$species))
# summary(as.factor(odMossdata$species))
# summary(as.factor(odLichendata$species))
# 
# #Convert to SPDF
# mossspdf<-SpatialPointsDataFrame(cbind(odMossdata$decimalLongitude,odMossdata$decimalLatitude),as.data.frame(odMossdata),proj4string = crs(norway))
# lichenspdf<-SpatialPointsDataFrame(cbind(odLichendata$decimalLongitude,odLichendata$decimalLatitude),as.data.frame(odLichendata),proj4string = crs(norway))
# vascspdf<-SpatialPointsDataFrame(cbind(odVascdata$decimalLongitude,odVascdata$decimalLatitude),as.data.frame(odVascdata),proj4string = crs(norway))
# 
# moss_utm<-spTransform(mossspdf,crs(norwayP))
# lichen_utm<-spTransform(lichenspdf,crs(norwayP))
# vasc_utm<-spTransform(vascspdf,crs(norwayP))
# 
# plot(norwayP)
# points(vasc_utm,pch=16,cex=0.1)
# 
# 
# # Clipping to Norway outline bufferd --------------------------------------
# 
# #Clip to remove points outside of Norway
# #Polygon of 1km buffer around Norway
# library(rgeos)
# norway1km<-gBuffer(norwayP, width=1000)
# #plot(norway1km)
# #plot(norwayP,add=T)
# 
# #Clip
# lichen_clip<-lichen_utm[which(!is.na(over(lichen_utm,norway1km))),]
# moss_clip<-moss_utm[which(!is.na(over(moss_utm,norway1km))),]
# vasc_clip<-vasc_utm[which(!is.na(over(vasc_utm,norway1km))),]
# plot(norwayP)
# points(lichen_utm,pch=16,col=2)
# points(lichen_clip,pch=16,col=3)
# plot(norwayP)
# points(moss_utm,pch=16,col=2)
# points(moss_clip,pch=16,col=3)
# plot(norwayP)
# points(vasc_utm,pch=16,col=2)
# points(vasc_clip,pch=16,col=3)
# 
# AllForestRedList<-rbind(lichen_clip,moss_clip,vasc_clip)
# AllForestRedList$PlantGroup<-c(rep('Lichen',times=nrow(lichen_clip)),rep('Bryophyte',times=nrow(moss_clip)),rep('Vascular',times=nrow(vasc_clip)))
# plot(norwayP)
# points(AllForestRedList[AllForestRedList$PlantGroup=='Lichen',])
# 
# #Merge with red list details
# match(AllForestRedList$species,redlistsp_all$Vitenskapelig.navn)
# ForestRedList_adb<-merge(AllForestRedList,redlistsp_all,by.x='species',by.y='Vitenskapelig.navn',all.x=T)
# 
# 
# # Final spp data ----------------------------------------------------------
# 
# write.table(lichen_clip,'Lichens_forest_redlisted_herbivory.csv')
# write.table(moss_clip,'Moss_forest_redlisted_herbivory.csv')
# write.table(vasc_clip,'Vascular_forest_redlisted_herbivory.csv')
# 
# fwrite(ForestRedList_adb@data,file='RedListedForestSpeciesNorwayBeite.csv')


# Environmental data ------------------------------------------------------


#Elevation
Norelev<-getData('alt',country='NOR')
plot(Norelev)

# #Climate
# Norbioclim<-getData('worldclim',var='bio',res=0.5,lon=5,lat=60)
# Norbioclim1<-getData('worldclim',var='bio',res=0.5,lon=5,lat=70)
# Norbioclim2<-getData('worldclim',var='bio',res=0.5,lon=40,lat=70)
# plot(Norbioclim[[1]])
# plot(Norbioclim1[[1]])
# mergclim<-merge(Norbioclim,Norbioclim1)
# mergclim1<-merge(mergclim,Norbioclim2)
# cropclim<-crop(mergclim1,Norelev)
# Norclimdat<-mask(cropclim,Norelev)
# plot(Norclimdat[[1]])
# 
# NorClimElev<-stack(Norclimdat,Norelev)
# names(NorClimElev)<-c(names(Norbioclim),names(Norelev))
# writeRaster(NorClimElev,'NorClimElev')
# 
# 
# NorClimElev<-stack('NorClimElev')
# NorClimElev
#  
# 
# #Land cover data
# landcover<-stack('landcoverstack')
# 
# 
# #Climate and elevation in same stack
# NorClimElev_utm<-projectRaster(NorClimElev,landcover[[1]],method='bilinear')
# 
# #Soil data
# soilph<-raster('geonode_phihox_m_sl2_250m.tif')
# soilph_utm<-projectRaster(soilph,landcover[[1]],method='bilinear')
# Norsoilph<-mask(crop(soilph_utm,landcover[[1]]),landcover[[1]])
# plot(Norsoilph)
# 
# #Herbivore data
# herbdat<-readOGR('KommuneMetabolicBiomass')
# herbdat1<-spTransform(herbdat,crs(landcover))
# moose1949<-rasterize(herbdat1,field='moose.1949',landcover[[1]])
# moose1959<-rasterize(herbdat1,field='moose.1959',landcover[[1]])
# moose1969<-rasterize(herbdat1,field='moose.1969',landcover[[1]])
# moose1979<-rasterize(herbdat1,field='moose.1979',landcover[[1]])
# moose1989<-rasterize(herbdat1,field='moose.1989',landcover[[1]])
# moose1999<-rasterize(herbdat1,field='moose.1999',landcover[[1]])
# moose2009<-rasterize(herbdat1,field='moose.2009',landcover[[1]])
# moose2015<-rasterize(herbdat1,field='moose.2015',landcover[[1]])
# moosestack<-stack(moose1949,moose1959,moose1969,moose1979,moose1989,moose1999,moose2009,moose2015)
# names(moosestack)<-c('moose1949','moose1959','moose1969','moose1979','moose1989','moose1999','moose2009','moose2015')
# writeRaster(moosestack,'Moose_metbio',overwrite=T)
# 
# red_deer1949<-rasterize(herbdat1,field='red_deer.1949',landcover[[1]])
# red_deer1959<-rasterize(herbdat1,field='red_deer.1959',landcover[[1]])
# red_deer1969<-rasterize(herbdat1,field='red_deer.1969',landcover[[1]])
# red_deer1979<-rasterize(herbdat1,field='red_deer.1979',landcover[[1]])
# red_deer1989<-rasterize(herbdat1,field='red_deer.1989',landcover[[1]])
# red_deer1999<-rasterize(herbdat1,field='red_deer.1999',landcover[[1]])
# red_deer2009<-rasterize(herbdat1,field='red_deer.2009',landcover[[1]])
# red_deer2015<-rasterize(herbdat1,field='red_deer.2015',landcover[[1]])
# red_deerstack<-stack(red_deer1949,red_deer1959,red_deer1969,red_deer1979,red_deer1989,red_deer1999,red_deer2009,red_deer2015)
# names(red_deerstack)<-c('red_deer1949','red_deer1959','red_deer1969','red_deer1979','red_deer1989','red_deer1999','red_deer2009','red_deer2015')
# writeRaster(red_deerstack,'red_deer_metbio',overwrite=T)
# 
# roe_deer1949<-rasterize(herbdat1,field='roe_deer.1949',landcover[[1]])
# roe_deer1959<-rasterize(herbdat1,field='roe_deer.1959',landcover[[1]])
# roe_deer1969<-rasterize(herbdat1,field='roe_deer.1969',landcover[[1]])
# roe_deer1979<-rasterize(herbdat1,field='roe_deer.1979',landcover[[1]])
# roe_deer1989<-rasterize(herbdat1,field='roe_deer.1989',landcover[[1]])
# roe_deer1999<-rasterize(herbdat1,field='roe_deer.1999',landcover[[1]])
# roe_deer2009<-rasterize(herbdat1,field='roe_deer.2009',landcover[[1]])
# roe_deer2015<-rasterize(herbdat1,field='roe_deer.2015',landcover[[1]])
# roe_deerstack<-stack(roe_deer1949,roe_deer1959,roe_deer1969,roe_deer1979,roe_deer1989,roe_deer1999,roe_deer2009,roe_deer2015)
# names(roe_deerstack)<-c('roe_deer1949','roe_deer1959','roe_deer1969','roe_deer1979','roe_deer1989','roe_deer1999','roe_deer2009','roe_deer2015')
# writeRaster(roe_deerstack,'roe_deer_metbio',overwrite=T)
# 
# #Stack together
# PredVars<-stack(NorClimElev_utm,landcover,Norsoilph,moosestack,red_deerstack,roe_deerstack)
# writeRaster(PredVars,'PredictorVariables',overwrite=T)


# Setup -------------------------------------------------------------------

PredVars<-stack('PredictorVariables')#Available here https://ntnu.box.com/s/wcmr0dgoyz2yu6ielw6er1pm7h0gaisa 
names(PredVars)[20:25]<-c('Elevation','Land_Cover','Forest_Type','Forest_Productivity','Vegetation_Type','SoilpH')


#Convert to factor variables
#Ratify land cover
PredVars$Land_Cover<-ratify(PredVars$Land_Cover)
ratlc<- levels(PredVars$Land_Cover)[[1]]
ratlc[["Land_Cover"]] <- c("Built-up","Agricultural","Forest","Open-natural vegetation","Mires","Glaciers/Ice/Snow","Freshwater","Sea","NA")
levels(PredVars$Land_Cover)<-ratlc
levelplot(PredVars$Land_Cover)

#Ratify forest productivty
PredVars$Forest_Productivity[PredVars$Forest_Productivity>18]<-NA #Class 99 ikke registrert. Gjelder skogområder som ligger utenfor AR5 kartleggingsområder
PredVars$Forest_Productivity<-ratify(PredVars$Forest_Productivity)
ratlcp<-levels(PredVars$Forest_Productivity)[[1]]
ratlcp[['Forest_Productivity']]<-c('Unproductive','Low','Medium','High')
levels(PredVars$Forest_Productivity)<-ratlcp
levelplot(PredVars$Forest_Productivity)

#Ratify forest type
PredVars$Forest_Type[PredVars$Forest_Type>33]<-NA #	99 Ikke registrert. Gjelder skogområder som ligger utenfor AR5 kartleggingsområder
PredVars$Forest_Type<-ratify(PredVars$Forest_Type)
ratlct<-levels(PredVars$Forest_Type)[[1]]
ratlct
ratlct[['ForestType']]<-c('Coniferous','Deciduous','Mixed')
levels(PredVars$Forest_Type)<-ratlct
levelplot(PredVars$Forest_Type)



#Correlation plots
pairs(PredVars)


#Reimport Species data
redlistforest<-fread('RedListedForestSpeciesNorwayBeite.csv')

#Convert to spdf
rlforsp<-SpatialPointsDataFrame(cbind(redlistforest$decimalLongitude,redlistforest$decimalLatitude),redlistforest,proj4string = crs(norway))
rlforsp_utm<-spTransform(rlforsp,crs(norwayP))

plot(norwayP)
points(rlforsp_utm[rlforsp_utm$PlantGroup.x=='Vascular',],col='green',pch=16,cex=0.2)
points(rlforsp_utm[rlforsp_utm$PlantGroup.x=='Lichen',],col='brown',pch=16,cex=0.2)
points(rlforsp_utm[rlforsp_utm$PlantGroup.x=='Bryophyte',],col='blue',pch=16,cex=0.2)

#Select points only with good geographic precision (coordinateuncertainty <1415m (sqrt(1000^2+1000^2)))
rlfor_use<-rlforsp_utm[!is.na(rlforsp_utm$coordinateUncertaintyInMeters) & rlforsp_utm$coordinateUncertaintyInMeters<=1415,]

#Mask out points outside of Norway coverage (worldclim data)
rlfor_useNor<-rlfor_use[!is.na(extract(PredVars[[10]],rlfor_use)),]

#Points per species
spNorwayGoodPrec<-with(rlfor_useNor@data,tapply(species,species,length))

#Select points Only in forest land-cover
extforest<-extract(PredVars$Land_Cover,rlfor_useNor)
extforest[is.na(extforest)]<-0
forestonly<-rlfor_use[extforest==30,]

#Points per species
spforocc<-with(forestonly@data,tapply(species,species,length))
spforocc
hist(spforocc)

#Only in forest with productivity
extforestprod<-extract(PredVars$Forest_Productivity,forestonly)
extforestprod[is.na(extforestprod)]<-0
forestprodonly<-forestonly[extforestprod>0,]
#Only in forest with forest type too
extforestprodtype<-extract(PredVars$Forest_Type,forestprodonly)
extforestprodtype[is.na(extforestprodtype)]<-0
forestprodonlytype<-forestprodonly[extforestprodtype>0,]

#Points per species - prod
spforprodocc<-with(forestprodonly@data,tapply(species,species,length))
spforprodocc
hist(spforprodocc)

#Points per species - prod and type
spforprodtypeocc<-with(forestprodonlytype@data,tapply(species,species,length))
spforprodtypeocc
hist(spforprodtypeocc)

#Counts of all poitns per species in different type
reccts<-as.data.frame(t(Reduce(function(...) merge(..., all=TRUE), list(t(as.data.frame(spNorwayGoodPrec)),
                                                t(as.data.frame(spforocc)),
                                                t(as.data.frame(spforprodocc)),
                                                t(as.data.frame(spforprodtypeocc))))))
                      
colnames(reccts)<-c('Records with forest type and productivity','Records with forest productivity','Records in forest','Records with good precisicion within worldclim')
reccts$Group<-rlfor_use$PlantGroup.x[match(rownames(reccts),rlfor_use$species)]
write.csv(reccts[,5:1],'RecordCountsperSpecies.csv')


#Points per higher taxa
tapply(forestonly$PlantGroup.x,forestonly$PlantGroup.x,length)

#Plot all species
col<-colorRampPalette('white')
#Make a raster to plot
noralt<-getData('alt',country='NOR')
nornull<-crop(projectRaster(noralt,crs=crs(norwayP)),norwayP)
values(nornull)<-0
levelplot(nornull,margin=F,colorkey=F
          ,col.regions=col,main='Red listed forest species', xlab=NULL, ylab=NULL, scales=list(draw=FALSE))+
  layer(sp.polygons(norwayP))+
  layer(sp.points(forestonly[forestonly$PlantGroup.x=='Vascular',],pch=16,cex=0.2,col='green'))+
  layer(sp.points(forestonly[forestonly$PlantGroup.x=='Lichen',],pch=16,cex=0.2,col='blue'))+
  layer(sp.points(forestonly[forestonly$PlantGroup.x=='Bryophyte',],pch=16,cex=0.5,col='tan4'))


# Plot species ------------------------------------------------------------
#Lichens
lichen_use<-forestonly[forestonly$PlantGroup.x=='Lichen',]
p<-list()
for(i in 1:length(levels(as.factor(lichen_use$species)))){
  p[[i]]<-levelplot(nornull,margin=F,colorkey=F,col.regions=col, xlab=NULL, ylab=NULL, scales=list(draw=FALSE)
                    ,main=list(label=levels(as.factor(lichen_use$species))[i],cex=0.7,fontface=3))+
    layer(sp.polygons(norwayP))+
    layer(sp.points(lichen_use[lichen_use$species==levels(as.factor(lichen_use$species))[i],],pch=16,cex=0.5))
}

#Make locations to plot (28 lichen species)
rowi<-c(rep(1:5,times=6))
coli<-c(rep(1:6,each=5))
tiff(width=210,height=297,units='mm',res=300,'Lichen_distributions.tif')
lattice.options(  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
                  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0)))
for (i in 1:length(p)){
  print(p[[i]],split=c(rowi[i],coli[i],5,6),more=T)}
dev.off()

#Mosses
moss_use<-forestonly[forestonly$PlantGroup.x=='Bryophyte',]
m<-list()
for(i in 1:length(levels(as.factor(moss_use$species)))){
  m[[i]]<-levelplot(nornull,margin=F,colorkey=F,col.regions=col, xlab=NULL, ylab=NULL, scales=list(draw=FALSE)
                    ,main=list(label=levels(as.factor(moss_use$species))[i],cex=0.7,fontface=3))+
    layer(sp.polygons(norwayP))+
    layer(sp.points(moss_use[moss_use$species==levels(as.factor(moss_use$species))[i],],pch=16,cex=0.5))
}

#Make locations to plot (5 moss species)
rowi<-c(rep(1:5,times=1))
coli<-c(rep(1,each=5))
tiff(width=210,height=60,units='mm',res=300,'Bryophyte_distributions.tif')
lattice.options(  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
                  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0)))
for (i in 1:length(m)){
  print(m[[i]],split=c(rowi[i],1,5,1),more=T)}
dev.off()


#Vascular
vasc_use<-forestonly[forestonly$PlantGroup.x=='Vascular',]
v<-list()
for(i in 1:length(levels(as.factor(vasc_use$species)))){
  v[[i]]<-levelplot(nornull,margin=F,colorkey=F,col.regions=col, xlab=NULL, ylab=NULL, scales=list(draw=FALSE)
                    ,main=list(label=levels(as.factor(vasc_use$species))[i],cex=0.7,fontface=3))+
    layer(sp.polygons(norwayP))+
    layer(sp.points(vasc_use[vasc_use$species==levels(as.factor(vasc_use$species))[i],],pch=16,cex=0.5))
}

#Make locations to plot (31 vasc species)
rowi<-c(rep(1:5,times=7))
coli<-c(rep(1:7,each=5))
tiff(width=210,height=297,units='mm',res=150,'Vascular distributions.tif')
lattice.options(  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
                  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0)))
for (i in 1:length(v)){
  print(v[[i]],split=c(rowi[i],coli[i],5,7),more=T)}
dev.off()


# Distribution modelling --------------------------------------------------
#Currently only points with forest productivity and type data 


#Background data
background<-read.table('BackgroundBiasCorrected.txt',header=T)#Vascular plants bias file used for Speed & Austrheim et al. 2017
backsp<-SpatialPoints(background,proj4string = crs(norwayP))
plot(norwayP)
points(backsp,cex=0.1,pch=16)

#Background in forest
backforestext<-extract(PredVars$Land_Cover,backsp)
backforestext[is.na(backforestext)]<-0
backforest<-backsp[backforestext==30,]
  
#MaxEnt Tuning ####
#Always allow linear features

 args_sp<-vector("list", length(levels(as.factor(forestprodonlytype$species))))
 names(args_sp)<-levels(as.factor(forestprodonlytype$species))
 nrecs<-summary(as.factor(forestprodonlytype$species))

 # for(i in 1:length(levels(as.factor(forestprodonlytype$species)))){
#   print(i)
#   print(levels(as.factor(forestprodonlytype$species))[[i]])
#   print(nrecs[i])  
#   if(nrecs[i]>=15){
#   tuneparameters<-ENMevaluate(occ=forestprodonlytype@coords[forestprodonlytype$species==levels(as.factor(forestprodonlytype$species))[[i]],],
#                             env=PredVars[[c(10,12,15,33,41,49,22:23)]],#Many NA values for forest type and productivity
#                             RMvalues = c(0.5,1,1.5,2,2.5,3,3.5,4,6,8), 
#                             fc = c("L", "LQ","LQH", "LQHP", "LQHPT"),
#                             categoricals=c("Forest_Type","Forest_Productivity"),
#                             method="block",
#                             bg.coords=backforest)
# tuneparameters@results[which.min(tuneparameters@results$AICc),]
# 
# b<-tuneparameters@results$rm[which.min(tuneparameters@results$AICc)]
# lin <-grepl('L',tuneparameters@results$features[which.min(tuneparameters@results$AICc)])
# quad<-grepl('Q',tuneparameters@results$features[which.min(tuneparameters@results$AICc)])
# prod<-grepl('P',tuneparameters@results$features[which.min(tuneparameters@results$AICc)])
# hing<-grepl('H',tuneparameters@results$features[which.min(tuneparameters@results$AICc)]) 
# thres<-grepl('T',tuneparameters@results$features[which.min(tuneparameters@results$AICc)])
# 
# args_sp[[i]]<-as.character(c(paste0('betamultiplier=',b),
#                       paste0('linear=',lin),
#                       paste0('quadratic=',quad),
#                       paste0('product=',prod),
#                       paste0('hinge=',hing),
#                       paste0('threshold=',thres),
#                       "-P",
#                       "-J",
#                       'replicates=5'
# ))}
# }

#Write tuning as list
#saveRDS(args_sp,file='Tuning')
args_sp<-readRDS('Tuning')

#MaxEnt modelling ####
me_list<-list()
#for(i in 1:4){
for(i in 1:length(levels(as.factor(forestprodonlytype$species)))){
  print(i)
  print(levels(as.factor(forestprodonlytype$species))[[i]])
  print(args_sp[[i]])
  #nrecs<-nrow(forestprodonlytype@coords[forestprodonlytype$species==levels(as.factor(forestprodonlytype$species))[[i]],])
  print(nrecs[i])  
  if(nrecs[i]>=25){
  me_list[[i]]<-maxent(p=forestprodonlytype@coords[forestprodonlytype$species==levels(as.factor(forestprodonlytype$species))[[i]],],
            x=PredVars[[c(10,12,15,33,41,49,22,23)]],#Many NA values for forest type and productivity
            factors=c("Forest_Type","Forest_Productivity"),
            a=backforest,
            args=args_sp[[i]],
            path=paste0('MaxEnt/',levels(as.factor(forestprodonlytype$species))[[i]]))
}}
names(me_list)<-levels(as.factor(forestprodonlytype$species))[nrecs>=25]


#Predictions ####
predictions<-list()
#for(i in 1:4){
for(i in 1:length(levels(as.factor(forestprodonlytype$species)))){
  #nrecs<-nrow(forestprodonlytype@coords[forestprodonlytype$species==levels(as.factor(forestprodonlytype$species))[[i]],])
  print(nrecs[i])  
  if(nrecs[i]>=25){ predictions[[i]]<-predict(me_list[[i]],PredVars)
} }
predictions_sp<-predictions[nrecs>=25]
names(predictions_sp)<-levels(as.factor(forestprodonlytype$species))[nrecs>25]

#Averaging predictions across multiple runs
meanpreds<-lapply(predictions_sp,function(x)calc(x,mean))

allspeciespredictions<-stack(meanpreds)
writeRaster(allspeciespredictions,'ModelPredictions/.tif',format='GTiff',bylayer=T,suffix=names(allspeciespredictions))
writeRaster(allspeciespredictions,'ModelPredictions/AllSpeciesPredictions.tif',format='GTiff')




# Using sdm package ------------------------------------------------------
library(sdm)
#https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.01881

sem<-function(x)sd(x,na.rm=T)/length(!is.na(x))

#InstallAll()#One time to install all dependent packages
#Testing

ulmgla<-forestprodonlytype[forestprodonlytype@data$species=='Ulmus glabra',]#@coords
names(ulmgla)[1]<-'ulmgla'
ulmgla<-cbind(ulmgla[,1],ulmgla@coords)
ulmgla$ulmgla<-'ulmgla'#Problems with spaces in species names...

#ulmgla<-data.frame(cbind(ulmgla=rep('ulmgla',times=nrow(ulmgla)),ulmgla))
sdmdatasetug<-sdmData(ulmgla~roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16+f(Forest_Productivity),train=ulmgla,
                    predictors=PredVars,bg=list(n=1000,method='gRandom',remove=TRUE))
sdmdataset
plot(sdmdataset)

mod1<-sdm(ulmgla~roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16+Forest_Productivity,data=sdmdatasetug,
          methods=c('glm','gam','gbm','cart','fda','rf'),
          replication=c('cv'),cv.folds=5)

#mod1<-sdm(ulmgla~roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16,data=sdmdataset,
#          methods=c('gbm','tree','mda','fda'),replication=c('cv','boot'),cv.folds=5,n=10)

roc(mod1)
rcurve(mod1)
getVarImp(mod1,1)# 1 model at a time

e1<-ensemble(mod1,pv1,setting=list(method='weighted',stat='AUC'))

# Making good dataframe ---------------------------------------------------

listdf<-list()
for (i in 1: length(levels(as.factor(forestprodonlytype$species)))){
#  for (i in 1:3){
  a<-forestprodonlytype[forestprodonlytype$species==levels(as.factor(forestprodonlytype$species))[i],]
  b<-cbind(a[,1],a@coords)
 # names(b)[1]<-levels(as.factor(forestprodonlytype$species))[i]
  names(b)[1]<-'species'
  listdf[[i]]<-b
  }

  
sdmdataset<-sdmData(ulmgla~roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16,train=ulmgla,
                      predictors=PredVars,bg=list(n=1000,method='gRandom',remove=TRUE))
  
s1<-sdmData(Ajuga_reptans~roe_deer2015,train = listdf[[1]],predictors = PredVars,
            bg=list(n=1000,method='gRandom',remove=TRUE))

slist<-list()
for(i in 1:3){
  #s[[i]]<-sdmData(paste(levels(as.factor(forestprodonlytype$species))[i])~
  s[[1]]<-sdmData(species~
                    roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16,
    train=listdf[[1]],
    predictors=PredVars,bg=list(n=1000,method='gRandom',remove=TRUE))
}


AllSpp <- do.call("rbind", listdf)
#Remove space from species name to avoid errors
AllSpp$species<-sub(" ","_",AllSpp$species)

#Make the sdm dataset with all species and relevent environmental variables (specifiy factors)
#1000 random points in the forest area
bg<-sampleRandom(PredVars$Forest_Productivity,1000,sp=T)

sdmdataset<-sdmData(species~roe_deer2015+red_deer2015+moose2015
                    +bio10_16+bio12_16+bio15_16+SoilpH
                    +f(Forest_Type)+f(Forest_Productivity)
                    ,train=AllSpp,predictors=PredVars,bg=list(sample(backforest,1000),remove=T))
sdmdataset


# All speceis model -------------------------------------------------------
#Model with all species concurrently 
#All species with >=20 records in forest

# sdm_Allspp<-sdm( Ajuga_reptans               +Anastrophyllum_donnianum   
#                 +Arnica_montana              +Asperugo_procumbens          
#                 +Campanula_barbata           +Campanula_cervicaria        +Cetrelia_olivetorum        
#                 +Cinna_latifolia             +Collema_curtisporum        
#                 +Collema_occultatum          +Crepis_praemorsa            +Cypripedium_calceolus      
#                 +Dactylorhiza_sambucina      +Epipogium_aphyllum          +Galium_sterneri            
#                 +Gentianella_campestris      +Gyalecta_flotowii           +Gyalecta_truncigena        
#                 +Gyalecta_ulmi               +Hackelia_deflexa            +Herbertus_stramineus        +Heterodermia_speciosa      
#                 +Lithospermum_officinale     +Malus_sylvestris            +Menegazzia_subsimilis       +Menegazzia_terebrata       
#                 +Opegrapha_vermicellifera    +Ophrys_insectifera          +Pectenia_cyanoloma         
#                 +Phaeophyscia_kairamoi       +Physconia_detersa           +Pseudorchis_albida          +Ramalina_dilacerata        
#                 +Ramalina_sinensis           +Ramboldia_subcinnabarina    +Rinodina_disjuncta          
#                 +Schismatomma_graphidioides  +Scorzonera_humilis          +Sorbus_lancifolia           +Sorbus_subpinnata           +Staurolemma_omphalarioides  +Taxus_baccata              
#                 +Thalictrum_minus            +Thalictrum_simplex          +Thelotrema_macrosporum                  
#                 +Ulmus_glabra                +Vicia_cassubica                     
#           
#                 ~roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16+Forest_Type+Forest_Productivity+SoilpH,
#                 data=sdmdataset,
#           methods=c('glm','gam','rf','gbm','mda','fda','brt'),
#           replication=c('cv'),cv.folds=5)
#saveRDS(sdm_Allspp,'SDM package/SDMAllSpecies')
#sdm_Allspp@run.info

#Extract model evaluations
#modeval<-cbind(sdm_Allspp@run.info,getEvaluation(sdm_Allspp))
modeval<-merge(sdm_Allspp@run.info,getEvaluation(sdm_Allspp),by='modelID')
modeval
write.csv(modeval,'SDM package/AllModels_Evaluation.csv')

aucmean<-with(modeval,tapply(AUC,list(method,species),mean))
aucsd<-with(modeval,tapply(AUC,list(method,species),sd))
barplot(aucmean,beside=T)

auc_spmethod<-matrix(paste(round(aucmean,3),round(aucsd,3),sep=' +/- '),dim(aucmean),dimnames = dimnames(aucmean))
write.csv(auc_spmethod,'SDM package/AUC_SpeciesMethod.csv')

aucmean
aucgrandmean<-apply(aucmean,2,mean,na.rm=T)
aucgrandsd<-apply(aucmean,2,sd,na.rm=T)

auc_sp<-cbind(round(aucgrandmean,3),round(aucgrandsd,3))
auc_sp
write.csv(auc_sp,'SDM package/AUC_Species.csv')

b1<-barplot(aucgrandmean)
arrows(b1,aucgrandmean+aucgrandsd,b1,aucgrandmean-aucgrandsd,length=0.05,code=3,angle=90)

#Extract variable importances
varimplist<-list()

#Null Df for models where variable importance not extracted
df1<-getVarImp(sdm_Allspp,id=i)@varImportance
df1$corTest<-NA
df1$AUCtest<-NA

for (i in 1:max(sdm_Allspp@run.info$modelID)){
  ifelse(sdm_Allspp@run.info$success[i]==TRUE,
         {
           ifelse(!is.null(getVarImp(sdm_Allspp,id=i)),
                  varimplist[[i]]<-getVarImp(sdm_Allspp,id=i)@varImportance,
                  varimplist[[i]]<-df1)
           varimplist[[i]]$species<-sdm_Allspp@run.info$species[i]
           varimplist[[i]]$method<-sdm_Allspp@run.info$method[i]
           varimplist[[i]]$repid<-sdm_Allspp@run.info$replicationID[i]}
         ,print(paste('Model failiure run ',i)))
}

AllVarImp<-do.call('rbind',varimplist)
AllVarImp
write.csv(AllVarImp,'SDM package/AllModelsVariableImportance.csv')

#Plot
varimpmean<-with(AllVarImp,tapply(corTest,list(variables,species),mean,na.rm=T))
varimpsem<-with(AllVarImp,tapply(corTest,list(variables,species),sem,na.rm=T))
par(mar=c(5,12,1,1))
b1<-barplot(varimpmean,beside=T,horiz=T,las=1,legend.text=T)
arrows(varimpmean+varimpsem,b1,varimpmean-varimpsem,b1,code=3,angle=90,length=0.05)

varimpdata<-t(matrix(paste((round(varimpmean,3)),(round(varimpsem,3)),sep= ' +/- '),dim((varimpmean))
                     ,dimnames=dimnames(varimpmean)))
write.csv(varimpdata,'SDM package/VariableImportance_species.csv')

#Response curves
#GetResponseCruve function gives object that can be plotted rather than just plots
#For some reason, does not work with rf, mda, fda
responsecurvelist<-list()
for (i in 1:length(levels(as.factor(sdm_Allspp@run.info$species)))){
  responsecurvelist[[i]]<-getResponseCurve(sdm_Allspp,id=sdm_Allspp@run.info$modelID[sdm_Allspp@run.info$species==levels(as.factor(sdm_Allspp@run.info$species))[i]
                                                                                     &sdm_Allspp@run.info$method%in% c('glm','gam','brt')]
                                           ,mean=T,main=levels(as.factor(sdm_Allspp@run.info$species))[i])
}

responsecurvelist[[1]]

saveRDS(responsecurvelist,'SDM package/ResponseCurvesobj')

plot(responsecurvelist[[2]])


#Ensemble models for each species
ensemblelist<-list()

#Need to remove attibute tables from factor variables to allow ensemble to work
pv1<-subset(PredVars,names(PredVars)[names(PredVars)%in%rownames(varimpmean)])
pv1$Forest_Productivity<-setValues(raster(pv1$Forest_Productivity),pv1$Forest_Productivity[])
pv1$Forest_Type<-setValues(raster(pv1$Forest_Type),pv1$Forest_Type[])

for (i in 1:length(levels(as.factor(sdm_Allspp@run.info$species)))){
  #for(i in 1:2){
  ensemblelist[[i]]<-ensemble(sdm_Allspp,newdata=pv1,filename=paste0('SDM package/EnsemblePredictions/',levels(as.factor(sdm_Allspp@run.info$species))[i]),
                              setting=list(method='weighted',stat='AUC'
                                           ,id=sdm_Allspp@run.info$modelID[sdm_Allspp@run.info$species==levels(as.factor(sdm_Allspp@run.info$species))[i]]))
}

plot(ensemblelist[[1]],main=levels(as.factor(sdm_Allspp@run.info$species))[1])
points(AllSpp[AllSpp$species==levels(as.factor(sdm_Allspp@run.info$species))[1],])
plot(ensemblelist[[2]],main=levels(as.factor(sdm_Allspp@run.info$species))[2])
points(AllSpp[AllSpp$species==levels(as.factor(sdm_Allspp@run.info$species))[2],])


#Niches 
niche(PredVars,ensemblelist[[1]],n=c('bio16_16','moose2015'))

#Forest spp over 20 recs

sdm_Allspp_for20<-sdm(
Ajuga_reptans
+Allium_scorodoprasum
+Anastrophyllum_donnianum
+Arnica_montana
+Asperugo_procumbens
+Campanula_barbata
+Campanula_cervicaria
+Cetrelia_olivetorum
+Cinna_latifolia
#+Cladonia_callosa
#+Cladonia_krogiana
+Collema_curtisporum
+Collema_occultatum
+Cotoneaster_laxiflorus
+Crepis_praemorsa
+Cypripedium_calceolus
+Dactylorhiza_sambucina
+Epipogium_aphyllum
#Fissidens_exilis
+Galium_sterneri
+Gentianella_campestris
#Gyalecta_derivata
+Gyalecta_flotowii
+Gyalecta_truncigena
+Gyalecta_ulmi
+Hackelia_deflexa
+Herbertus_stramineus
+Heterodermia_speciosa
+Lithospermum_officinale
+Malus_sylvestris
+Menegazzia_subsimilis
+Menegazzia_terebrata
+Opegrapha_vermicellifera
+Ophrys_insectifera
+Pectenia_cyanoloma
+Phaeophyscia_kairamoi
+Physconia_detersa
+Pseudorchis_albida
+Ramalina_dilacerata
+Ramalina_sinensis
+Ramboldia_subcinnabarina
+Rinodina_disjuncta
+Schismatomma_graphidioides
+Scorzonera_humilis
+Sorbus_lancifolia
+Sorbus_subpinnata
+Staurolemma_omphalarioides
+Taxus_baccata
+Thalictrum_minus
+Thalictrum_simplex
+Thelotrema_macrosporum
+Ulmus_glabra
+Vicia_cassubica
+Vicia_orobus
~roe_deer2015+red_deer2015+moose2015+bio10_16+bio12_16+Forest_Type+Forest_Productivity+SoilpH,
data=sdmdataset,
methods=c('glm','gam','rf','gbm','mda','fda','brt'),
replication=c('cv'),cv.folds=5)


saveRDS(sdm_Allspp_for20,'SDM package/SDMAllSpecies_for20')
sdm_Allspp_for20@run.info

# Figures -----------------------------------------------------------------

#Set the order to plot spp (moss, lichen, vascular)
spgroups<-read.csv('SppGroups.csv',sep=';')
spgroups<-droplevels(spgroups[!spgroups$Species%in%c('Vicia_orobus','Allium_scorodoprasum','Cotoneaster_laxiflorus'),])
levels(spgroups$Species)[(levels(spgroups$Species)%in%'Collema_occultatum')]<-"Rostania_occultata"
levels(spgroups$Species)[(levels(spgroups$Species)%in%'Schismatomma_graphidioides')]<-"Schismatomma_pericleum"
levels(spgroups$Species)[(levels(spgroups$Species)%in%'Hackelia_deflexa')]<-"Lappula deflexa"
spgroups<-spgroups[order(spgroups$Group,spgroups$Species),]
#Extract model evaluations
#modeval<-cbind(sdm_Allspp@run.info,getEvaluation(sdm_Allspp))
modeval_for20<-merge(sdm_Allspp_for20@run.info,getEvaluation(sdm_Allspp_for20),by='modelID')
levels(modeval_for20$species)[(levels(modeval_for20$species)%in%'Collema_occultatum')]<-"Rostania_occultata"
levels(modeval_for20$species)[(levels(modeval_for20$species)%in%'Schismatomma_graphidioides')]<-"Schismatomma_pericleum"
levels(modeval_for20$species)[(levels(modeval_for20$species)%in%'Hackelia_deflexa')]<-"Lappula deflexa"

write.csv(modeval_for20,'SDM package/AllModels_Evaluation_for20.csv')

aucmean_for20<-with(modeval_for20,tapply(AUC,list(method,species),mean))
aucsd_for20<-with(modeval_for20,tapply(AUC,list(method,species),sd))
barplot(aucmean_for20,beside=T)

auc_spmethod_for20<-data.frame(t(matrix(paste(round(aucmean_for20,3),round(aucsd_for20,3),sep=' +/- '),dim(aucmean_for20),dimnames = dimnames(aucmean_for20))))
auc_spmethod_for20$Species<-rownames(auc_spmethod_for20)
auc_spmethod_for20<-merge(auc_spmethod_for20,spgroups,by.x='Species',by.y='Species')
write.csv(auc_spmethod_for20,'SDM package/AUC_SpeciesMethod_for20.csv')

aucmean_for20<-aucmean_for20[,colnames(aucmean_for20)%in%spgroups$Species]
aucgrandmean_for20<-apply(aucmean_for20,2,mean,na.rm=T)
aucgrandsd_for20<-apply(aucmean_for20,2,sd,na.rm=T)

auc_sp_for20<-cbind(round(aucgrandmean_for20,3),round(aucgrandsd_for20,3))
auc_sp_for20
write.csv(auc_sp_for20,'SDM package/AUC_Species_for20.csv')

m1<-merge(spgroups,data.frame(species=names(aucgrandmean_for20),AUC=aucgrandmean_for20),by.x='Species',by.y='species')
with(m1,tapply(aucgrandmean_for20,Group,mean))
with(m1,tapply(aucgrandmean_for20,Group,sd))

names(aucgrandmean_for20)<-sub("_", " ", names(aucgrandmean_for20))
sporder<-match(sub("_", " ",vf3$spp),names(aucgrandmean_for20))
sporder


tiff('Figures/AUC.tif',width=1200,height=800,units='px',pointsize = 20)
par(mar=c(11,5,3,1))
par(xpd=T)
b1<-barplot(aucgrandmean_for20[sporder],ylim=c(0,1.1),
            las=2,cex.names=0.9,font=3,yaxt='n',
            col=c(grey(0.2),grey(0.5),grey(0.8))[spgroups$Group])
axis(2,las=1)
title(ylab='AUC')
arrows(b1,aucgrandmean_for20[sporder]+aucgrandsd_for20[sporder],b1,aucgrandmean_for20[sporder]-aucgrandsd_for20[sporder],length=0.05,code=3,angle=90)
#abline(v=(b1[3]-b1[2])/2+b1[2],lty=2,lwd=2)
#abline(v=(b1[23]-b1[22])/2+b1[22],lty=2,lwd=2)
legend(45,1.25,fill=c(grey(0.2),grey(0.5),grey(0.8)),c('Bryophyte','Lichen','Vascular'),cex=0.9)
#text(b1[1.5],1.1,'B')
dev.off()

summary(aucgrandmean_for20)

#Extract variable importances
varimplist_for20<-list()

#Null Df for models where variable importance not extracted
df1_for20<-getVarImp(sdm_Allspp_for20,id=i)@varImportance
df1_for20$corTest<-NA
df1_for20$AUCtest<-NA

for (i in 1:max(sdm_Allspp_for20@run.info$modelID)){
  ifelse(sdm_Allspp_for20@run.info$success[i]==TRUE,
         {
           ifelse(!is.null(getVarImp(sdm_Allspp_for20,id=i)),
                  varimplist_for20[[i]]<-getVarImp(sdm_Allspp_for20,id=i)@varImportance,
                  varimplist_for20[[i]]<-df1)
           varimplist_for20[[i]]$species<-sdm_Allspp_for20@run.info$species[i]
           varimplist_for20[[i]]$method<-sdm_Allspp_for20@run.info$method[i]
           varimplist_for20[[i]]$repid<-sdm_Allspp_for20@run.info$replicationID[i]}
         ,print(paste('Model failiure run ',i)))
}

AllVarImp_for20<-do.call('rbind',varimplist_for20)
AllVarImp_for20

AllVarImp_for20group<-merge(AllVarImp_for20)
write.csv(AllVarImp_for20,'SDM package/AllModelsVariableImportance_for20.csv')

#Plot
varimpmean_for20<-with(AllVarImp_for20,tapply(corTest,list(variables,species),mean,na.rm=T))
varimpsem_for20<-with(AllVarImp_for20,tapply(corTest,list(variables,species),sem))
par(mar=c(5,12,1,1))
b1<-barplot(varimpmean_for20,beside=T,horiz=T,las=1,legend.text=T)
arrows(varimpmean_for20[sporder]+varimpsem_for20[sporder],b1,varimpmean_for20[sporder]-varimpsem_for20[sporder],b1,code=3,angle=90,length=0.05)

varimpdata_for20<-t(matrix(paste((round(varimpmean_for20,3)),(round(varimpsem_for20,3)),sep= ' +/- '),dim((varimpmean_for20))
                     ,dimnames=dimnames(varimpmean_for20)))


write.csv(varimpdata_for20,'SDM package/VariableImportance_species_for20.csv')
vim1<-as.data.frame(varimpdata_for20)
vim1$species<-rownames(varimpdata_for20)
vimpgroup<-merge(vim1,spgroups,by.x='species',by.y='Species')
write.csv(vimpgroup,'SDM package/VariableImportance_species_for20_Group.csv')

#Fix spp names
levels(AllVarImp_for20$species)%in% "Collema_occultatum" <-"Rostania_occultata"
levels(AllVarImp_for20$species)%in% "Schismatomma_graphidioides"<-"Schismatomma_pericleum"
levels(AllVarImp_for20$species)%in% "Hackelia_deflexa"<-"Lappula deflexa"


#Ordering by spp groups
vf1<-as.data.frame(t(varimpmean_for20))
vf1$spp<-rownames(vf1)
vf2<-merge(vf1,spgroups,by.x='spp',by.y='Species')
vf2$spp[vf2$spp=='Collema_occultatum']<-"Rostania_occultata"
vf2$spp[vf2$spp=='Schismatomma_graphidioides']<-"Schismatomma_pericleum"
vf2$spp[vf2$spp=='Hackelia_deflexa']<-"Lappula deflexa"
vf3<-vf2[order(vf2$Group,vf2$spp),]
vfse1<-as.data.frame(t(varimpsem_for20))
vfse1$spp<-rownames(vfse1)
vfse2<-merge(vfse1,spgroups,by.x='spp',by.y='Species')
vfse2$spp[vfse2$spp=='Collema_occultatum']<-"Rostania_occultata"
vfse2$spp[vfse2$spp=='Schismatomma_graphidioides']<-"Schismatomma_pericleum"
vfse2$spp[vfse2$spp=='Hackelia_deflexa']<-"Lappula deflexa"
vfse3<-vfse2[order(vfse2$Group,vfse2$spp),]

tiff('Figures/VarImp.tif',width=1200,height=800,units='px',pointsize = 20)
par(mar=c(12,5,1,1))
par(xpd=F)
b1<-barplot(t(as.matrix(vf3[,6:8])),beside=T,cex.names=0.8,horiz=F,las=1,legend.text = T,
            names.arg=sub("_"," ",vf3$spp),las=2,font=3,yaxt='n',ylim=c(0,0.67),
            args.legend = list(legend=c('Moose', 'Red deer', 'Roe deer')))
axis(2,las=1)
title(ylab='Variable importance')
arrows(b1,t(as.matrix(vf3[,6:8]))+t(as.matrix(vfse3[,6:8])),b1,t(as.matrix(vf3[,6:8]))-t(as.matrix(vfse3[,6:8])),code=3,length=0.05,angle=90)
abline(v=(b1[2,3]-b1[2,2])/2+b1[2,2],lty=2,lwd=2)
abline(v=(b1[2,23]-b1[2,22])/2+b1[2,22],lty=2,lwd=2)
text(1,0.6,'B')
text(50,0.6,'Lichens')
text(130,0.6,'Vascular')
#abline(h=1/6)
dev.off()

#Try as a two part figure
tiff('Figures/VarImp2.tif',width=1200,height=1600,units='px',pointsize = 20)
par(mfrow=c(2,1))
par(mar=c(12,5,2,1))
par(xpd=F)
b1<-barplot(t(as.matrix(vf3[1:22,6:8])),beside=T,cex.names=0.8,horiz=F,las=1,legend.text = T,
            names.arg=sub("_"," ",vf3$spp[1:22]),las=2,font=3,yaxt='n',ylim=c(0,0.67),
            args.legend = list(legend=c('Moose', 'Red deer', 'Roe deer')))
axis(2,las=1)
title(ylab='Variable importance')
title(main='Bryophytes and Lichens')
arrows(b1,t(as.matrix(vf3[1:22,6:8]))+t(as.matrix(vfse3[1:22,6:8])),b1,t(as.matrix(vf3[1:22,6:8]))-t(as.matrix(vfse3[1:22,6:8])),code=3,length=0.05,angle=90)
abline(v=(b1[2,3]-b1[2,2])/2+b1[2,2],lty=2,lwd=2)
#text(1,0.6,'B')
#text(50,0.6,'Lichens')

b1<-barplot(t(as.matrix(vf3[23:47,6:8])),beside=T,cex.names=0.8,horiz=F,las=1,legend=F,
            names.arg=sub("_"," ",vf3$spp[23:47]),las=2,font=3,yaxt='n',ylim=c(0,0.67))
            #args.legend = list(legend=c('Moose', 'Red deer', 'Roe deer')))
axis(2,las=1)
title(ylab='Variable importance')
title(main='Vascular plants')
arrows(b1,t(as.matrix(vf3[23:47,6:8]))+t(as.matrix(vfse3[23:47,6:8])),b1,t(as.matrix(vf3[23:47,6:8]))-t(as.matrix(vfse3[23:47,6:8])),code=3,length=0.05,angle=90)
dev.off()

#Correlation figure
tiff('Figures/VarImp_Corrs.tif',width=800,height=1200,units='px',pointsize = 20)
par(mfrow=c(2,1))
par(mar=c(10,5,1,1))
#Average varimp across spp
vf4<-vf3
names(vf4)<-c('Species','MST','MAP','Forest productivity','Forest type','Moose biomass','Red deer biomass','Roe deer biomass', 'Soil pH','Group')
par(cex=0.8)
boxplot(vf4[,2:9],las=2,cex=0.8,ylab='Variable importance')
points(1:8,apply(vf4[,2:9],2,mean),pch=3)
mtext('a',side=3,adj=-0.14)

#Correlation between variable importances
#plot(vf3$moose2015,vf3$red_deer2015)
#points(vf3$moose2015,vf3$roe_deer2015,pch=16)
corrs <- cor(vf4[,2:9])
corrplot(corrs, type = "upper",diag=F, tl.col = "black", tl.srt = 45)
mtext('b',side=3,adj=0)
dev.off()


#Response curves
#GetResponseCruve function gives object that can be plotted rather than just plots
#For some reason, does not work with rf, mda, fda
responsecurvelist_for20<-list()
for (i in 1:length(levels(as.factor(sdm_Allspp_for20@run.info$species)))){
  print(i)
  responsecurvelist_for20[[i]]<-getResponseCurve(sdm_Allspp_for20,id=sdm_Allspp_for20@run.info$modelID[sdm_Allspp_for20@run.info$species==levels(as.factor(sdm_Allspp_for20@run.info$species))[i]
                                                                                     &sdm_Allspp_for20@run.info$method%in% c('glm','gam','brt')]
                                           ,mean=T,main=levels(as.factor(sdm_Allspp_for20@run.info$species))[i])
}

responsecurvelist_for20[[1]]

saveRDS(responsecurvelist_for20,'SDM package/ResponseCurvesobj_for20')

plot(responsecurvelist_for20[[2]],main=levels(as.factor(sdm_Allspp_for20@run.info$species))[2])

#Plotting response curves against moose density
#par(mfrow=c(5,10))
for(i in sporder){
  print(i)
  plot(responsecurvelist_for20[[i]]@response$moose2015[,1],apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,mean)
     ,type='l',main=levels(modeval_for20$species)[i],xlab='Moose density',ylab='Response',las=1)
  lines(responsecurvelist_for20[[i]]@response$moose2015[,1],apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,mean)
        +apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,sem),lty=2)
  lines(responsecurvelist_for20[[i]]@response$moose2015[,1],apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,mean)
      -apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,sem),lty=2)
}

levels(modeval_for20$species)[levels(modeval_for20$species)=="Schismatomma_graphidioides"]<-"Schismatomma_pericleum"
#Plot response curves for species-herbivore combinations where VarImp herbivore >0.2
tiff('Figures/RespCurces.tif',width=1600,height=1600,units='px',pointsize = 20)
tiff('Figures/RespCurces2.tif',width=1600,height=1600,units='px',pointsize = 20)
par(mfcol=c(7,3))
par(mar = c(0, 4, 3, 0), oma = c(6, 3, 0.5, 0.5))
par(tcl = -0.25)
#Moose
for(i in which(varimpmean_for20[5,]>0.25)[c(2:6,1)]){
  print(i)
  print(colnames(varimpmean_for20)[i])
  plot(responsecurvelist_for20[[i]]@response$moose2015[,1],apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,mean)
       ,type='l',main=paste0(sub("_"," ",levels(modeval_for20$species)[i]),"\n (",round(varimpmean_for20[5,i],3),")"),xlab='',ylab='',las=1,cex.main=1.5,xaxt='n')
  lines(responsecurvelist_for20[[i]]@response$moose2015[,1],apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,mean)
        +apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,sem),lty=2)
  lines(responsecurvelist_for20[[i]]@response$moose2015[,1],apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,mean)
        -apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,sem),lty=2)
  }
axis(1,cex=1.5)
mtext(side=2,'Response',cex=1.5,outer=T)
#Dummy plot to fill 'gap'
plot(responsecurvelist_for20[[i]]@response$moose2015[,1],apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,mean)
     ,main="",xlab='',ylab='',las=1,cex.main=0.9,xaxt='n',type='n',axes=F)
mtext(side=1,expression('Moose' ~(kg~km^{-2})),outer=T,at=0.15,padj=T,line=2,cex=1.5)

#Red deer
for(i in which(varimpmean_for20[6,]>0.25)){
  print(i)
plot(responsecurvelist_for20[[i]]@response$red_deer2015[,1],apply(responsecurvelist_for20[[i]]@response$red_deer2015[,2:ncol(responsecurvelist_for20[[i]]@response$red_deer2015)],1,mean)
       ,type='l',main=paste0(sub("_"," ",levels(modeval_for20$species)[i]),"\n (",round(varimpmean_for20[6,i],3),")"),ylab='',las=1,cex.main=1.5,xaxt='n')
  lines(responsecurvelist_for20[[i]]@response$red_deer2015[,1],apply(responsecurvelist_for20[[i]]@response$red_deer2015[,2:ncol(responsecurvelist_for20[[i]]@response$red_deer2015)],1,mean)
        +apply(responsecurvelist_for20[[i]]@response$red_deer2015[,2:ncol(responsecurvelist_for20[[i]]@response$red_deer2015)],1,sem),lty=2)
  lines(responsecurvelist_for20[[i]]@response$red_deer2015[,1],apply(responsecurvelist_for20[[i]]@response$red_deer2015[,2:ncol(responsecurvelist_for20[[i]]@response$red_deer2015)],1,mean)
        -apply(responsecurvelist_for20[[i]]@response$red_deer2015[,2:ncol(responsecurvelist_for20[[i]]@response$red_deer2015)],1,sem),lty=2)
}
axis(1)
mtext(side=1,expression('Red deer' ~(kg~km^{-2})),outer=T,at=0.52,padj=1,line=2,cex=1.5)

#Roe deer
#par(mfrow=c(3,3))
for(i in which(varimpmean_for20[7,]>0.25)[c(4,3,5:6)]){
  print(i)
  plot(responsecurvelist_for20[[i]]@response$roe_deer2015[,1],apply(responsecurvelist_for20[[i]]@response$roe_deer2015[,2:ncol(responsecurvelist_for20[[i]]@response$roe_deer2015)],1,mean)
       ,type='l',main=paste0(sub("_"," ",levels(modeval_for20$species)[i]),"\n (",round(varimpmean_for20[7,i],3),")"),xlab='',ylab='',las=1,cex.main=1.5,xaxt='n')
  lines(responsecurvelist_for20[[i]]@response$roe_deer2015[,1],apply(responsecurvelist_for20[[i]]@response$roe_deer2015[,2:ncol(responsecurvelist_for20[[i]]@response$roe_deer2015)],1,mean)
        +apply(responsecurvelist_for20[[i]]@response$roe_deer2015[,2:ncol(responsecurvelist_for20[[i]]@response$roe_deer2015)],1,sem),lty=2)
  lines(responsecurvelist_for20[[i]]@response$roe_deer2015[,1],apply(responsecurvelist_for20[[i]]@response$roe_deer2015[,2:ncol(responsecurvelist_for20[[i]]@response$roe_deer2015)],1,mean)
        -apply(responsecurvelist_for20[[i]]@response$roe_deer2015[,2:ncol(responsecurvelist_for20[[i]]@response$roe_deer2015)],1,sem),lty=2)
}
axis(1)
#Dummy plot to fill 'gap'
plot(responsecurvelist_for20[[i]]@response$moose2015[,1],apply(responsecurvelist_for20[[i]]@response$moose2015[,2:ncol(responsecurvelist_for20[[i]]@response$moose2015)],1,mean)
     ,main="",xlab='',ylab='',las=1,cex.main=0.9,xaxt='n',type='n',axes=F)
mtext(side=1,expression('Roe deer' ~(kg~km^{-2})),outer=T,at=0.85,padj=1,line=2,cex=1.5)
dev.off()



#Need to remove attibute tables from factor variables to allow ensemble to work
pv1<-subset(PredVars,names(PredVars)[names(PredVars)%in%rownames(varimpmean)])
pv1$SoilpH<-mask(pv1$SoilpH,pv1$bio10_16)
pv1$Forest_Productivity<-setValues(raster(pv1$Forest_Productivity),pv1$Forest_Productivity[])
pv1$Forest_Type<-setValues(raster(pv1$Forest_Type),pv1$Forest_Type[])

#Plot env vars
pvplot1<-levelplot(pv1$bio10_16/10,margin=F,main="Mean summer temperature",scales=list(draw=F),
                   colorkey=list(title=expression(~degree~C)),par.settings='YlOrRdTheme')+
  layer(sp.polygons(norwayP),col=grey(0.5))
pvplot2<-levelplot(pv1$bio12_16,margin=F,main="Annual precipitation",scales=list(draw=F),
                   colorkey=list(title='mm'),par.settings='RdBuTheme')+
  layer(sp.polygons(norwayP),col=grey(0.5))
pvplot3<-levelplot(pv1$SoilpH/10,margin=F,main="Soil pH",scales=list(draw=F),
                   colorkey=list(title='pH'),par.settings='RdBuTheme')+
  layer(sp.polygons(norwayP),col=grey(0.5))
pvplot4<-levelplot(PredVars$Forest_Type,margin=F,main="Forest Type",scales=list(draw=F))+
  layer(sp.polygons(norwayP),col=grey(0.5))
pvplot5<-levelplot(PredVars$Forest_Productivity,margin=F,main="Forest Productivity",scales=list(draw=F))+
  layer(sp.polygons(norwayP),col=grey(0.5))
pvplot6<-levelplot(pv1$moose2015,margin=F,main="Moose metabolic biomass",scales=list(draw=F),
                   colorkey=list(title=expression(~ kg~km^{-2}),space='right'),par.settings='YlOrRdTheme')+
  layer(sp.polygons(norwayP),col=grey(0.5))
pvplot7<-levelplot(pv1$red_deer2015,margin=F,main="Red deer density",scales=list(draw=F),
                   colorkey=list(title=expression(~ kg~km^{-2}),space='right'),par.settings='YlOrRdTheme')+
  layer(sp.polygons(norwayP),col=grey(0.5))
pvplot8<-levelplot(pv1$roe_deer2015,margin=F,main="Roe deer density",scales=list(draw=F),
                   colorkey=list(title=expression(~kg~km^{-2}),space='right'),par.settings='YlOrRdTheme')+
  layer(sp.polygons(norwayP),col=grey(0.5))

pvplot1
pvplot2
pvplot3
pvplot4
pvplot5
pvplot6
pvplot7
pvplot8

tiff('Figures/PredVars.tif',width=800,height=1200,units='px',pointsize = 20)
grid.arrange(pvplot1,pvplot6,pvplot2,pvplot7,pvplot3,pvplot8,pvplot4,pvplot5,ncol=2)
dev.off()

ensemblelist_for20<-list()
for (i in 1:length(levels(as.factor(sdm_Allspp_for20@run.info$species)))){
#  for(i in 1:2){
  ensemblelist_for20[[i]]<-ensemble(sdm_Allspp_for20,newdata=pv1,filename=paste0('SDM package/EnsemblePredictions_for20/',levels(as.factor(sdm_Allspp_for20@run.info$species))[i]),
                              setting=list(method='weighted',stat='AUC'
                                           ,id=sdm_Allspp_for20@run.info$modelID[sdm_Allspp_for20@run.info$species==levels(as.factor(sdm_Allspp_for20@run.info$species))[i]]))
}

plot(ensemblelist_for20[[1]],main=levels(as.factor(sdm_Allspp_for20@run.info$species))[1])
points(AllSpp[AllSpp$species==levels(as.factor(sdm_Allspp_for20@run.info$species))[1],])
plot(ensemblelist_for20[[2]],main=levels(as.factor(sdm_Allspp_for20@run.info$species))[2])
points(AllSpp[AllSpp$species==levels(as.factor(sdm_Allspp_for20@run.info$species))[2],])

#pv1949 Replace 2015 with 1949 data
pv49<-pv1
pv49$moose2015<-PredVars$moose1949
pv49$red_deer2015 <-PredVars$red_deer1949
pv49$roe_deer2015<-PredVars$roe_deer1949


#Predictions from sdm
#Average within method
#Multiple spp

#Selected spp
selspp<-names(c(which(varimpmean_for20[5,]>0.25),which(varimpmean_for20[6,]>0.25),which(varimpmean_for20[7,]>0.25)[3:6]))
selsppun<-unique(selspp)

sppreds2015<-list()
sppreds1949<-list()
for(i in 1:length(selsppun)){
  print(i)
  sppreds2015[i]<-predict(sdm_Allspp_for20,pv1,filename=paste0('ModelPredictions/ModelPreds2015/sppreds2015_',selsppun[i]),
                   species=selsppun[i],mean=T)
  sppreds1949[i]<-predict(sdm_Allspp_for20,pv49,filename=paste0('ModelPredictions/ModelPreds1949/sppreds1949_',selsppun[i]),
                   species=selsppun[i],mean=T)
}


#Averages across methods
sppredictions2015<-lapply(sppreds2015,function(x)calc(x,mean))
sppredictions1949<-lapply(sppreds1949,function(x)calc(x,mean))
sppredictionsstack1949<-stack(sppredictions1949)
sppredictionsstack1949<-sppredictionsstack1949[[c(1:10,13:16)]]
names(sppredictionsstack1949)<-selsppun
sppredictionsstack2015<-stack(sppredictions2015)
sppredictionsstack2015<-sppredictionsstack2015[[c(1:10,13:16)]]
names(sppredictionsstack2015)<-selsppun
#Write
writeRaster(sppredictionsstack1949,filename=paste0('ModelPredictions/SelSppPredAverages/1949/',names(sppredictionsstack1949)),format='GTiff', bylayer=TRUE)
writeRaster(sppredictionsstack2015,filename=paste0('ModelPredictions/SelSppPredAverages/2015/',names(sppredictionsstack2015)),format='GTiff', bylayer=TRUE)

#Change in range
#Plot on diverging colour scale
diverge0 <- function(p, ramp) {
  # p: a trellis object resulting from rasterVis::levelplot
  # ramp: the name of an RColorBrewer palette (as character), a character 
  #       vector of colour names to interpolate, or a colorRampPalette.
  require(RColorBrewer)
  require(rasterVis)
  if(length(ramp)==1 && is.character(ramp) && ramp %in% 
     row.names(brewer.pal.info)) {
    ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
    ramp <- colorRampPalette(ramp)
  } else if(!is.function(ramp)) 
    stop('ramp should be either the name of a RColorBrewer palette, ', 
         'a vector of colours to be interpolated, or a colorRampPalette.')
  rng <- range(p$legend[[1]]$args$key$at)
  s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
  i <- findInterval(rng[which.min(abs(rng))], s)
  zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
  p$legend[[1]]$args$key$at <- s[zlim]
  p$par.settings$regions$col <- ramp(1000)[zlim[-length(zlim)]]
  p
}

predchanges<-sppredictionsstack1949
for(i in 1:nlayers(sppredictionsstack1949)){
  predchanges[[i]]<-sppredictions2015[[i]]-sppredictions1949[[i]]}
p1<-levelplot(predchanges,margin=F,scales=list(draw=F),names.attr=selsppun)
BuRd<-colorRampPalette(BuRdTheme()$regions$col)
diverge0(p1,BuRd)


#Try with alternative breaks
pc100<-predchanges*100
rng <- range(cellStats(pc100, range))
lim <- ceiling(log(abs(rng), 2))
b <- sort(c(0, unique(unlist(mapply(function(x, y) y*2^(0:x), lim, sign(rng))))))
b[1] <- rng[1]
b[length(b)] <- rng[2]

p.strip <- list(cex=0.5, lines=1,font=3)
p <- levelplot(pc100, par.settings=BuRdTheme(), at=b, 
                 colorkey=list(height=0.8, labels=list(at=b[c(1:4,8,12:16)], labels=round(b[c(1:4,8,12:16)], 0))),
               names.attr=sub("_"," ",selsppun),scales=list(draw=F),par.strip.text=p.strip)
p

#Remove spp with too low n
pc100<-pc100[[c(1:10,13:16)]]

sppredictionsstack1949<-sppredictionsstack1949[[c(1:10,13:16)]]
#Species data
plotord<-c(7,2,8,3,9,12,4,5,10,6,1,11,13,14)
listsppts<-list()
selsppunOrder<-selsppun[plotord]
for (i in plotord){
  print(i)
  listsppts[i]<-forestprodonlytype[forestprodonlytype$species%in%sub("_"," ",selsppunOrder[i]),]}

selsppun[selsppun=="Schismatomma_graphidioides"]<-"Schismatomma_pericleum"
tiff('Figures/RangeChange.tif',width=297,height=210,units='mm',pointsize = 20,res=300)
p.strip <- list(cex=0.75, lines=1,font=3)
p <- levelplot(pc100[[plotord]], par.settings=BuRdTheme(), at=b, 
               colorkey=list(height=0.6,title='Change in predictions\n (%)\n \n ', title.gpar = list(cex = 1), labels=list(at=b[c(1:4,8,12:16)], labels=round(b[c(1:4,8,12:16)], 0),cex=1)),
               names.attr=sub("_"," ",selsppun[plotord]),scales=list(draw=F),par.strip.text=p.strip)+
  layer(sp.polygons(norwayP,lwd=0.1))+
  layer(sp.points(listsppts[[panel.number()]],pch=1,cex=0.5,col=1,alpha=0.8))
p
dev.off()

tiff('Figures/Range2015.tif',width=297,height=210,units='mm',pointsize = 20,res=300)
p.strip <- list(cex=0.75, lines=1,font=3)
p <- levelplot(sppredictionsstack2015[[plotord]]*100,par.settings=YlOrRdTheme(), 
               colorkey=list(height=0.6,title='Prediction 2015\n (%)\n \n', title.gpar = list(cex = 1), cex=1),
               names.attr=sub("_"," ",selsppun[plotord]),scales=list(draw=F),par.strip.text=p.strip)+
  layer(sp.polygons(norwayP,lwd=0.1))+
  layer(sp.points(listsppts[[panel.number()]],pch=1,cex=0.5,col=1,alpha=0.8))
p
dev.off()

tiff('Figures/Range1949.tif',width=297,height=210,units='mm',pointsize = 20,res=300)
p.strip <- list(cex=0.75, lines=1,font=3)
p <- levelplot(sppredictionsstack1949[[plotord]]*100, par.settings=YlOrRdTheme(), 
               colorkey=list(height=0.6,title='Prediction 1949\n (%)\n \n', title.gpar = list(cex = 1), cex=1),
               names.attr=sub("_"," ",selsppun[plotord]),scales=list(draw=F),par.strip.text=p.strip)+
  layer(sp.polygons(norwayP,lwd=0.1))+
  layer(sp.points(listsppts[[panel.number()]],pch=1,cex=0.5,col=1,alpha=0.8))
p
dev.off()


pvP<-pv1
names(pvP)<-c('MST','MAP','Forest type','Forest Prod.','Soil pH','Moose biomass','Red deer biomass','Roe deer biomass')
tiff('Figures/PredPairs.tif',width=1600,height=1000,units='px',pointsize = 20)
pairs(pvP)
dev.off()

#Plot herbivore differences

pdiffX<-levelplot(stack(pv1$moose2015-pv49$moose2015,pv1$red_deer2015-pv49$red_deer2015,pv1$roe_deer2015-pv49$roe_deer2015),names.attr=c('Moose','Red deer','Roe deer'),
                  main='Change in metabolic biomass',colorkey=list(title=expression(~kg~km^{2})),scales=list(draw=F))+
  layer(sp.polygons(norwayP))


pvdiffM<-levelplot(pv1$moose2015-pv49$moose2015,margin=F,main="Moose",scales=list(draw=F),
                   colorkey=list(title=expression(~kg~km^{-2}),space='right'))+
    layer(sp.polygons(norwayP))
#pvdiffM
pvdiffRed<-levelplot(pv1$red_deer2015-pv49$red_deer2015,margin=F,main="Red deer",scales=list(draw=F),
                   colorkey=list(title=expression(~kg~km^{-2}),space='right'))+
  layer(sp.polygons(norwayP))
pvdiffRoe<-levelplot(pv1$roe_deer2015-pv49$roe_deer2015,margin=F,main="Roe deer",scales=list(draw=F),
                   colorkey=list(title=expression(~kg~km^{-2}),space='right'))+
  layer(sp.polygons(norwayP))

tiff('Figures/CervidDiff.tif',width=1200,height=500,units='px',pointsize=20)
#diverge0(pdiffX,BuRd)
grid.arrange(diverge0(pvdiffM,BuRd),diverge0(pvdiffRed,BuRd),diverge0(pvdiffRoe,BuRd),ncol=3)
dev.off()


#Range changes as KMZ
pc100[[plotord]]
rangechangesp<-stack(pc100[[plotord]])
names(rangechangesp)<-sub("_"," ",selsppun[plotord])
rangechangespLL<-projectRaster(rangechangesp,crs=crs(norway))
p <- levelplot(pc100, par.settings=BuRdTheme(), at=b, 
               colorkey=list(height=0.8, labels=list(at=b[c(1:4,8,12:16)], labels=round(b[c(1:4,8,12:16)], 0))),
               names.attr=sub("_"," ",selsppun),scales=list(draw=F),par.strip.text=p.strip)

KML(rangechangespLL,'HabitatSuitabilityChange',col=rev(brewer.pal(11,'RdBu')),overwrite=T)

KML(rangechangespLL,'HabitatSuitabilityChange',col=p$par.settings$regions$col,overwrite=T)

KML(projectRaster(sppredictionsstack2015[[c(1:10,13:16)]][[plotord]]*100,crs=crs(norway)),'Predictions2015',col=brewer.pal(9, 'YlOrRd'),overwrite=T)
KML(projectRaster(sppredictionsstack1949[[c(1:10,13:16)]][[plotord]]*100,crs=crs(norway)),'Predictions1949',col=brewer.pal(9, 'YlOrRd'),overwrite=T)


kml(rangechangespLL[[1]],file='HabSuitChange')


#VIF
vif(cbind(sdmdataset@features[c(2:6,8)],cbind(as.numeric(sdmdataset@features$Forest_Type),as.numeric(sdmdataset@features$Forest_Productivity))))
vif(as.numeric(sdmdataset@features))
