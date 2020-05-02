# program

# same as the readme file

# reindeerSheep <- raster::stack('data/large/reindeerSheep.grd')
# PredVars <- raster::stack("data/large/PredictorVariables.grd")
# names(PredVars)[20:25]<-c('Elevation','Land_Cover','Forest_Type','Forest_Productivity','Vegetation_Type','SoilpH')
# PredVars$SoilpH <- PredVars$SoilpH/10
# PredVars <- PredVars[[c(21:23, 25, 31, 39, 47)]]
# names(PredVars)
# NorClimElev <- raster::stack('data/NorClimElev.grd')
# NorClimElev$temp <- NorClimElev$temp/10
# newproj <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# cats               <- raster::stack(PredVars$Land_Cover, PredVars$Forest_Type, PredVars$Forest_Productivity) 
# cats               <- raster::projectRaster(cats, crs = newproj, method='ngb')
# num                <- raster::stack(PredVars$SoilpH,
#                                     PredVars$moose1999,
#                                     PredVars$red_deer1999,
#                                     PredVars$roe_deer1999)
# num                <- raster::projectRaster(num, crs = newproj, method='bilinear')
# reindeerSheep2     <- raster::projectRaster(reindeerSheep, num[[1]], method='bilinear')  
# NorClimElev2       <- raster::projectRaster(NorClimElev, num[[1]], method='bilinear')
# IV                 <- raster::stack(cats, num, reindeerSheep2, NorClimElev2)
# IV$SoilpH               <- raster::mask(IV$SoilpH,              IV$roe_deer1999)
# IV$Land_Cover           <- raster::mask(IV$Land_Cover,          IV$roe_deer1999)
# IV$Forest_Type          <- raster::mask(IV$Forest_Type,         IV$roe_deer1999)
# IV$Forest_Productivity  <- raster::mask(IV$Forest_Productivity, IV$roe_deer1999)
# writeRaster(IV, 'data/IV', overwrite=TRUE) # 76 MB
# rm(cats, IV, NorClimElev, Norclimdat, Norbioclim2, Norbioclim, Norbioclim1, NorClimElev2, num, PredVars, reindeerSheep, reindeerSheep2, newproj)



# Run top section ####
myIVs              <- raster:: stack('data/IV.grd')
names(myIVs)

myIVs$Forest_Productivity[myIVs$Forest_Productivity>18]<-NA
myIVs$Forest_Productivity <-   raster::ratify(myIVs$Forest_Productivity)
ratlcp <- raster::levels(myIVs$Forest_Productivity)[[1]]
ratlcp[['Forest_Productivity']] <- 
  c('Unproductive',
    'Low',
    'Medium',
    'High')
levels(myIVs$Forest_Productivity) <- ratlcp

myIVs$Land_Cover[myIVs$Land_Cover == 10 |
                   myIVs$Land_Cover == 70 |
                   myIVs$Land_Cover == 81 |
                   myIVs$Land_Cover == 82] <- 98
myIVs$Land_Cover <- raster::ratify(myIVs$Land_Cover)
ratlc               <- raster::levels(myIVs$Land_Cover)[[1]]
ratlc$Land_Cover <- c(
  "Agricultural",
  "Forest",
  "Open-natural vegetation",
  "Mires",
  "others",
  "NA")
levels(myIVs$Land_Cover) <- ratlc

myIVs$Forest_Type[myIVs$Forest_Type>33]<-NA
myIVs$Forest_Type<-raster::ratify(myIVs$Forest_Type)
ratlct<-raster::levels(myIVs$Forest_Type)[[1]]
ratlct[['ForestType']] <-
  c('Coniferous','Deciduous','Mixed')
levels(myIVs$Forest_Type) <- ratlct

source("./R/spList.R")
mySpecies <- sl()


# End top section ####





# TEST DOWNLOAD   ####
mySpecies2 <- mySpecies[1:10]

# Let's do a test loop without downloading anything, just seeing how many records there are.

nOccurences_df <- data.frame(species = mySpecies2, 
                  nOccurences = as.numeric(NA))

for(i in 1:length(mySpecies2)){
  myName  <- mySpecies2[i]
  myName2 <- stringr::str_split(myName, " ")[[1]]
  nOccurences_df$nOccurences[i] <- dismo::gbif(myName2[1], myName2[2], download = F) 
}

nOccurences_df



# SELECT SPECIES ####
mySpecies3 <- mySpecies[mySpecies == c("Primula scandinavica", "Kobresia simpliciuscula")]
# write.csv(mySpecies3, 'data/species_list.csv')




# DONT RUN - GET OCCURENCES  ####

#for(i in 1:length(mySpecies3)){
#  myName  <- mySpecies3[i]
#  myName2 <- stringr::str_split(myName, " ")[[1]]
#  
#  assign(
#    sub(' ', '_', mySpecies3[i]), 
#    dismo::gbif(myName2[1], myName2[2], 
#                download = T,
#                geo = T, 
#                sp = F) 
#  )
#}
#
#
#
##Now we can turn the dataframes into spatialPointsDataFrames, define the CRS, and plot the points. The dataset comes as #lonlat.
#
#for(i in 1:length(mySpecies3)){
#  
#  d <- get(
#    sub(' ', '_', mySpecies3[i]))
#  
#  sp::coordinates(d) <- ~lon + lat
#  sp::proj4string(d) <- sp::proj4string(raster::raster())
#  
#  assign(
#    sub(' ', '_', mySpecies3[i]),  d)
#}
#mapview::mapview(Kobresia_simpliciuscula, 
#                 map.types = c("Esri.WorldShadedRelief",
#                               "Esri.WorldImagery"),
#                 cex = 5, lwd = 0,
#                 alpha.regions = 0.5,
#                 col.regions = "blue")
#
#
#
#mapview::mapview(Primula_scandinavica, 
#                 map.types = c("Esri.WorldShadedRelief",
#                               "Esri.WorldImagery"),
#                 cex = 5, lwd = 0,
#                 alpha.regions = 0.5,
#                 col.regions = "blue")
#
## First, notice that Kobresia is called Carex in GBIF, but Kobresia in ADB. This is not a problem and they are #recognised as synonyms. The Kobresia is a widespread species, whereas the Primula is endemic to Norway and Sweden. We #only need the points that fall on Norway. First we need something to clip against, so we'll get an outline of Norway. 
#
## outline <- norway()
##saveRDS(outline, "data/large/outline_Norway.RData") # 1.8MB
#outline <- readRDS("data/large/outline_Norway.RData")
#raster::plot(outline)
#
#
## Now to clip away occurences outside this polygon (can take a few minutes)
#
#
#for(i in 1:length(mySpecies3)){
#
#d <- get(
#sub(' ', '_', mySpecies3[i]))
#
#d <- raster::crop(d, outline)
#
#assign(
#sub(' ', '_', mySpecies3[i]),  d)
#}
#
#
#
## Lets see it it worked.
#
#
#mapview::mapview(Kobresia_simpliciuscula, 
#                 map.types = c("Esri.WorldShadedRelief",
#                               "Esri.WorldImagery"),
#                 cex = 5, lwd = 0,
#                 alpha.regions = 0.5,
#                 col.regions = "blue")
#
#
#
#mapview::mapview(Primula_scandinavica, 
#                 map.types = c("Esri.WorldShadedRelief",
#                               "Esri.WorldImagery"),
#                 cex = 5, lwd = 0,
#                 alpha.regions = 0.5,
#                 col.regions = "blue")
#
##Looks like it. Now we just need to get this over to UTM32 to match the IV data, and save it on file.
#
#for(i in 1:length(mySpecies3)){
#  
#  d <- get(
#    sub(' ', '_', mySpecies3[i]))
#  
#  d <- sp::spTransform(d, myIVs[[1]]@crs)
#  
#  assign(
#    sub(' ', '_', mySpecies3[i]),  d)
#}
#
#oDat <- get(sub(' ', '_', mySpecies3[1]))
#for(i in 2:length(mySpecies3)){
#  oDat <- rbind(oDat, get(sub(' ', '_', mySpecies3[i])))
#}
#saveRDS(oDat, 'data/large/oDat.RData')
#rm(oDat)




# LOAD OCCURENCES ####
oDat <- readRDS('data/large/oDat.RData')
raster::plot(myIVs$Forest_Productivity)
raster::plot(oDat,add=T)


# CHECK SAMPLE SIZES

t <- myIVs$elev
df <- data.frame(species = NA,
                points = as.numeric(NA),
                unique = as.numeric(NA))
for(i in 1:length(unique(oDat$species))){
  s       <- unique(oDat$species)[i]
  df[i,1] <- paste(s)
  t1      <- oDat[oDat$species == s,]
  df[i,2] <- length(t1)
  u       <- raster::rasterize(t1, t, 1, fun = "count")
  df[i,3] <- length(u[u>0])
}
df


# SDM DATA OBJECTS

library(sdm)

### SDM-data

for(i in 1:length(mySpecies3)){
  rm(s, s2, d, dat)
  s    <- unique(oDat$species)[i]
  s2   <- paste0(s, "_d")
  d    <- oDat[oDat$species == s,]
  dat  <- sdm::sdmData(species~temp+prec+SoilpH+TundraHerbivores,
                       train = d,
                       predictors = myIVs,
                       bg = list(n=1000, method = "gRandom"))
  
  assign(s2, dat)
  sdm::write.sdm(dat, paste0("models/sdmData/", s2), overwrite=TRUE)
}



Primula_scandinavica_d
Carex_simpliciuscula_d
