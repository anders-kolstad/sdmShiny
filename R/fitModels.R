# TOPP ####

# Fitting SDM models  
# Anders L. Kolstad may 2020


# Here we attempt to fit SDM models for multiple red-listed Norwegian species.

# This code chunk will get us up to speed.
library(sdm)
source("./R/spList.R")
myIVs      <- raster::stack('data/IV.grd')
#oDat       <- readRDS('data/large/oDat.RData')   # 2 species only as in the Rdm documentation file
# oDat       <- readRDS('data/allOccurences.RData') # All species



#mySpecies  <- c("Primula scandinavica", "Kobresia simpliciuscula")
#mySpecies  <- sl()
#mySpecies <- 
  unique(oDat$species)
# UPDATE NEEDED: remove records prior to 1990

## SDM-data object ####

for(i in 1:length(mySpecies)){
  s    <- unique(oDat$species)[i]
  s2   <- paste0(s, "_d")
  d    <- oDat[oDat$species == s,]
  
  dat  <- sdm::sdmData(species~temp+prec+SoilpH+TundraHerbivores,
                       train = d,
                       predictors = myIVs,
                       bg = list(n=length(d)+1000, method = "gRandom"))
  
  assign(s2, dat)
  sdm::write.sdm(dat, paste0("models/sdmData/", s2), overwrite=TRUE)
  rm(s, s2, d, dat)
}


#This takes the names over to the gbif approved names (eg Kobresia to Carex). 
#The 'number of records' is not equal to n+1000, even if remove=F is added. 
#Not sure what is happening.


## 5 x 3 models
#We don't have any independent test data so I'll use bootstrapping to partition test data, and I'll do that 5 times. I shouldn't do much less because there are very few data points for some of the species. I think records from inside the same 1km cells will be counted as duplicates and removed. I will use 3 methods as well, resulting in 3 x 5 = 15 models per species. 


for(i in 1:length(mySpecies)){
  
  s       <- unique(oDat$species)[i]
  file1   <- paste0("models/sdmData/", s, "_d.sdd")
  obj     <- paste0(s, "_m")
  file2   <- paste0("models/sdmModels/", obj)
  d      <- sdm::read.sdm(file1)
  
  mod <- sdm::sdm(.~.,
                  data = d, 
                  methods = c('glm', 'gam', 'maxent'),
                  replication = c('boot'), n=5)     
  
  sdm::write.sdm(mod, file2, overwrite=TRUE)
  
  rm(s, file1, obj, file2, d, mod)
}



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



## Best candidate model ####

#To make predictions from the 5 x 3 models we need to run the ensamble function again, this time with altered IVs in the newdata argument. The function is probably too slow to run on the fly. This I need to test on the io server later. We could choose one method, eg maxent, and just use that for all species, and with only one replication. Then we could use the predict function in raster which will be quicker. Or we can look at these 15 models we have generated and choose the best one from there. That would be a bit safer.



for(i in 1:length(mySpecies)){
  s <- unique(oDat$species)[i]
  s2 <- paste0(s, "_m")
  mod <- sdm::read.sdm(paste0("models/sdmModels/", s, "_m.sdm"))
  s3 <- getEvaluation(mod)
  file <- paste0("models/sdmModels/", s, "_bcm")
  best <- s3$modelID[which.max(s3$AUC)]
  mod2 <- mod[[best]]
  
  # save the model object for easy lading in the shiny app.
  sdm::write.sdm(mod2, file, overwrite=TRUE)
  
  
  
  mod3 <- raster::predict(mod2,
                          object = myIVs, 
                          filename = paste0("models/predictions/", s, "_bcm.img"), 
                          overwrite=TRUE)
  
  #assign(paste0(s, "_best"), 
  #       mod3)
  
  rm(s, s2, s3, file, best, mod, mod2, mod3)
}



# Plotting - just change the ending of the file name and rerun line with a few different species
#p <- stack("models/predictions/Primula_scandinavica_bcm.img")
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






