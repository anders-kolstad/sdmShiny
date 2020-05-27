# TOPP ####

# Fitting SDM models  
# Anders L. Kolstad may 2020


# Here we attempt to fit SDM models for multiple red-listed Norwegian species.

# This code chunk will get us up to speed.
library(sdm)
source("./R/spList.R")
myIVs      <- raster::stack('data/IV.grd')
# removing categorical layers - they dont add much and cause trouble with NAs
myIVs <- raster::stack(myIVs$SoilpH, myIVs$moose1999, 
                       myIVs$red_deer1999, myIVs$roe_deer1999,
                        myIVs$TundraHerbivores, myIVs$temp, myIVs$prec, myIVs$elev)
names(myIVs)
#oDat       <- readRDS('data/large/oDat.RData')   # 2 species only as in the Rdm documentation file
# oDat       <- readRDS('data/allOccurences.RData') # All species
dim(oDat)
oDat <- oDat[oDat$year>1989,]
dim(oDat)

#mySpecies  <- c("Primula scandinavica", "Kobresia simpliciuscula")

mySpecies <- unique(oDat$species) # 132
mySpecies2  <- sl(df=T)
myAlpine <- mySpecies2$mySpList[mySpecies2$type == "alpine"]
myForest <- mySpecies2$mySpList[mySpecies2$type == "forest"]
myAlpine <- sub(' ', '_', myAlpine)
myForest <- sub(' ', '_', myForest)

comb <- c(as.character(myAlpine), as.character(myForest))
(comb <- comb[which(duplicated(comb))])
# shold these be modelled as forets plants or as alpine plants? Perhaps a combination, adding sheep and reindeer as IV


myAlpine <- mySpecies[mySpecies %in% myAlpine]
myAlpine <- myAlpine[!myAlpine %in% comb]

myForest <- mySpecies[mySpecies %in% myForest]
myForest <- myForest[!myForest %in% comb]


# ratify
#myIVs$Forest_Type[myIVs$Forest_Type>33]<-NA
#myIVs$Forest_Type<-raster::ratify(myIVs$Forest_Type)
#ratlct<-raster::levels(myIVs$Forest_Type)[[1]]
#ratlct[['ForestType']] <-
#  c('Coniferous','Deciduous','Mixed')
#levels(myIVs$Forest_Type) <- ratlct
#cols3 <- colorRampPalette(c("darkgreen", "orange", "blue" ))
#rasterVis::levelplot(myIVs$Forest_Type, main= "Forest type", col.regions = cols3)


# Background data ####
# This uses a bias file.
#bg <- readRDS('data/background.RData')
#bg <- bg[!is.na(bg$decimalLatitude),]
#bg <- bg[!is.na(bg$decimalLongitude),]
#bg <- select(bg, decimalLatitude, decimalLongitude)
#sp::coordinates(bg) <- ~decimalLongitude + decimalLatitude
#sp::proj4string(bg) <- sp::proj4string(raster::raster())
#bg <- spTransform(bg,         crs(myIVs))
##source("R/norway.R")
##nor <- norway(lonlat = TRUE)
##plot(nor)
##points(bg)
#bg2 <- as.data.frame(bg)
#rasValue <- raster::extract(myIVs, bg2)
#combinePointValue <- cbind(bg2,rasValue)
#head(combinePointValue)
#t <- combinePointValue[complete.cases(combinePointValue),]
#names(t)[1] <- 'lon'
#names(t)[2] <- 'lat'
#head(t)


# background data v2
library(dismo)
bg <- dismo::randomPoints(myIVs$elev, 5000)
bg <- as.data.frame(bg)
colnames(bg) <- c("lon", "lat")
sp::coordinates(bg) <- ~lon + lat
sp::proj4string(bg) <- sp::proj4string(myIVs)
#mapview::mapview(bg, 
#                 map.types = c("Esri.WorldShadedRelief",
#                               "Esri.WorldImagery"),
#                 cex = 5, lwd = 0,
#                 alpha.regions = 0.5,
#                 col.regions = "blue")
rasValue <- raster::extract(myIVs, bg)
bg2 <- as.data.frame(bg)
bg3 <- cbind(bg2,rasValue)
head(bg3)
#bg3x <- bg3[complete.cases(bg3),] # alsost the same
bg4 <- bg3[sample(1:nrow(bg3), 1000),]
rm(bg, bg2, bg3, rasValue)




## SDM-data object ####

for(i in 1:length(myAlpine)){
  s    <- myAlpine[i]
  s2   <- paste0(s, "_d")
  d    <- oDat[oDat$species == s,]
  
  dat  <- sdm::sdmData(species~temp+prec+SoilpH+TundraHerbivores,
                       train = d,
                       predictors = myIVs,
                       bg = bg4)
  #assign(s2, dat)
  sdm::write.sdm(dat, paste0("models/sdmData/", s2), overwrite=TRUE)
  #rm(s, s2, d, dat)
}



for(i in 1:length(myForest)){
  s    <- myForest[i]
  s2   <- paste0(s, "_d")
  d    <- oDat[oDat$species == s,]
  if(length(d)>20){
    dat  <- sdm::sdmData(species~
                           temp+
                           prec+
                           SoilpH+
                           #Forest_Type+
                           moose1999+
                           red_deer1999+
                           roe_deer1999,
                         train = d,
                         predictors = myIVs,
                         bg = bg4)
  
  #assign(s2, dat)
  sdm::write.sdm(dat, paste0("models/sdmData/", s2), overwrite=TRUE)
  }
}



for(i in 1:length(comb)){
  s    <- comb[i]
  s2   <- paste0(s, "_d")
  d    <- oDat[oDat$species == s,]
  if(length(d)>20){
    dat  <- sdm::sdmData(species~
                           temp+
                           prec+
                           SoilpH+
                           #f(Forest_Type)+
                           moose1999+
                           red_deer1999+
                           roe_deer1999+
                           TundraHerbivores,
                         train = d,
                         predictors = myIVs,
                         bg = bg4)
    
    #assign(s2, dat)
    sdm::write.sdm(dat, paste0("models/sdmData/", s2), overwrite=TRUE)
  }
  #rm(s, s2, d, dat)
}
#The 'number of records' is not equal to n+1000, even if remove=F is added. 
#Not sure what is happening.



# get a list of the species we have sdmData objects for
myS <- list.files("models/sdmData/", pattern = ".sdd")
myS2 <- as.list(NA)
for(i in 1:length(myS)){
  myS2[i] <- 
    paste(
      stringr::str_split(myS[i], "_")[[1]][1],
      stringr::str_split(myS[i], "_")[[1]][2],
      collapse = "_",
      sep="_")
}
mySpecies <- unique(as.character(myS2))
length(mySpecies) # 85


# get a list of the species we have sdmModel objects for
myS <- list.files("shiny/sdmModels/", pattern = ".sdm")
myS2 <- as.list(NA)
for(i in 1:length(myS)){
  myS2[i] <- 
    paste(
      stringr::str_split(myS[i], "_")[[1]][1],
      stringr::str_split(myS[i], "_")[[1]][2],
      collapse = "_",
      sep="_")
}
mySpeciesM <- unique(as.character(myS2))
length(mySpeciesM) # 


# GAM ####
# There is  compromise between model accuracy and model object file size (need to be low to allow many species in the shinyapps bundle), I have ended up using a single method (gam, one of the most successful in trils) with 5 replicates (because gams sometimes fail, and varImp variesa lot between runs). They should  differ more when the number of observations is very low because then you could get 'unlucky' with the partitioning. 
mySpecies2 <- mySpecies[!mySpecies %in% mySpeciesM]
length(mySpecies2) # 85
mySpecies <- mySpecies2
mySpecies <- mySpecies[1:10]
# Aphanes_australis is extremly rare
for(i in 1:length(mySpecies)){
  
  s       <- mySpecies[i]
  file1   <- paste0("models/sdmData/", s, "_d.sdd")
  obj     <- paste0(s, "_3gams")
  file2   <- paste0("shiny/sdmModels/", obj)
  d      <- sdm::read.sdm(file1)
  
  mod <- sdm::sdm(.~.,
                  data = d, 
                  methods = c('gam'),
                  replication = c('boot'), n=3)     
  if(any(mod@run.info$success)){
  sdm::write.sdm(mod, file2, overwrite=TRUE)
  }
}

# These models are very simple and overfitted, but it's not such a big deal for the purpuse they are intended for.
# However, we can manually remove those that are way off by simply deleting the model object from the folder shiny/sdmModels. 
# I'm going to do that for Stellaria hebecalyx.


# END #### -----------------------------------------------------------------------------
# The stuff below is just for future reference


# pick the best one out of the three
for(i in 12:length(mySpecies)){
  s <- mySpecies[i]
  mod <- sdm::read.sdm(paste0("models/sdmModels/", s, "_maxent3.sdm"))
  if(any(mod@run.info$success)){
    s3 <- getEvaluation(mod)
    file <- paste0("shiny/sdmModels/", s, "_maxent")
    best <- s3$modelID[which.max(s3$AUC)]
    mod2 <- mod[[best]]
    
    # save the model object for easy lading in the shiny app.
    sdm::write.sdm(mod2, file, overwrite=TRUE)
  }
}


## 5 x 3 models ####
#We don't have any independent test data so I'll use bootstrapping to partition test data, and I'll do that 5 times. I shouldn't do much less because there are very few data points for some of the species. Records from inside the same 1km cells will be counted as duplicates and removed. I will use 3 methods as well, resulting in 3 x 5 = 15 models per species. 

mySpecies2 <- mySpecies
mySpecies <- mySpecies[1:30]

for(i in 1:length(mySpecies)){
  
  s       <- mySpecies[i]
  file1   <- paste0("models/sdmData/", s, "_d.sdd")
  obj     <- paste0(s, "_m")
  file2   <- paste0("models/sdmModels/", obj)
  d      <- sdm::read.sdm(file1)
  
  mod <- sdm::sdm(.~.,
                  data = d, 
                  methods = c('glm', 'gam', 'mars'),
                  replication = c('boot'), n=5)     
  
  sdm::write.sdm(mod, file2, overwrite=TRUE)
  
  #rm(s, file1, obj, file2, d, mod)
}


## Best candidate model ####

#To make predictions from the 5 x 3 models we need to run the ensamble function again, this time with altered IVs in the newdata argument. The function is probably too slow to run on the fly. This I need to test on the io server later. We could choose one method, eg maxent, and just use that for all species, and with only one replication. Then we could use the predict function in raster which will be quicker. Or we can look at these 15 models we have generated and choose the best one from there. That would be a bit safer.



for(i in 12:length(mySpecies)){
  s <- mySpecies[i]
  s2 <- paste0(s, "_m")
  mod <- sdm::read.sdm(paste0("models/sdmModels/", s, "_m.sdm"))
  if(any(mod@run.info$success)){
  s3 <- getEvaluation(mod)
  file <- paste0("shiny/sdmModels/", s, "_bcm")
  best <- s3$modelID[which.max(s3$AUC)]
  mod2 <- mod[[best]]
  
  # save the model object for easy lading in the shiny app.
  sdm::write.sdm(mod2, file, overwrite=TRUE)
  
  mod3 <- raster::predict(mod2,
                          object = myIVs, 
                          filename = paste0("shiny/predictions/", s, "_bcm.img"), 
                          overwrite=TRUE)
  
  #assign(paste0(s, "_best"), 
  #       mod3)
  
  #rm(s, s2, s3, file, best, mod, mod2, mod3)
  }
}



# Plotting - just change the ending of the file name and rerun line with a few different species
#p <- stack("shiny/predictions/Carex_bicolor_bcm.img")
#p <- p > -Inf
#plot(p)
#rm(p)






# Let's bring the models back in to the environment. They're 3-20 20 MB each on file, or 17-55 when unzipped in the environment.
#for(i in 1:length(mySpecies)){
#  s      <- unique(oDat$species)[i]
#  s2     <- paste0(s, "_m")
#  file   <- paste0("models/sdmModels/", s2, ".sdm")
#  
#  assign(s2, sdm::read.sdm(file))
#}




### Ensemble ####
#Let's put the 5 x 3 models together and make a map of the current habitat suitability, using ensamble. 

for(i in 1:length(mySpecies)){
  s <- unique(oDat$species)[i]
  d <- sdm::read.sdm(paste0("models/sdmModels/", s, "_m.sdm"))
  #
  
  # using this because overwrite=TRUE dont work
  fn <- paste0("models/predictions/", s, "_ens.img")
  if(file.exists(fn)){file.remove(fn)}
  
  mod <- sdm::ensemble(d,
                       newdata = myIVs, 
                       filename = fn, 
                       overwrite=TRUE,   
                       setting = list(method='weighted', stat = 'AUC')
  )     
  
  #assign(paste0(s, "_ens"), 
  #      mod)
  
  rm(s, d, fn, mod)
}


# It a good idea to plot some  of them here.
#p <- stack("models/predictions/Primula_scandinavica_ens.img")
#plot(p)
#rm(p)






## Replicated single method  ####
#This approach chooses one method, maxent, and does five replicated partitionings, and then the raster::predict function to make maps. I'll keep these models to to a speed test later in the actuall Shiny app on the io server.


for(i in 1:length(mySpecies)){
  
  s       <- unique(oDat$species)[i]
  file1   <- paste0("models/sdmData/", s, "_d.sdd")
  
  obj    <- paste0(s, "_5maxent")
  file2   <- paste0("models/sdmModels/", obj)
  d       <- sdm::read.sdm(file1)
  
  mod <- sdm::sdm(.~.,
                  data = d, 
                  methods = 'maxent',   
                  replication = c('boot'), n=5)     
  
  sdm::write.sdm(mod, file2, overwrite=TRUE)
  #assign(obj, mod)
  
  mod2 <- raster::predict(mod,
                          object = myIVs,
                          mean = T,
                          filename = paste0("models/predictions/", obj, ".img"), 
                          overwrite=TRUE)
  #assign(obj2, mod2)

}
 

# Plotting - just change the ending of the file name and rerun line with a few different species
#p <- stack("models/predictions/Carex_simpliciuscula_5maxent.img")
#plot(p)


# Option for presence-absence maps instead:
# ev     <- sdm::getEvaluation(mod, stat = "threshold")
#  th     <- mean(ev$threshold)
#  s      <- paste0(unique(oDat$species)[i], "_5maxent_pa")
#  mod2[] <- ifelse(mod2[] >= th, 1, 0)
#  assign(s, mod2)
  



#Some runs gives white spots on this map. I think I have traced this problem to sdmData function - the 'number of records' ois not the sum of the pseudo-absences plus unique spatial points in the occurence data. Usually I get a fairly close number when I run the for-loop with an empty environment (eg 800-1060 records for Primula). But if I run it again straight after I get maybe around 200. It's easier to check the coverage when there is just one colour.






