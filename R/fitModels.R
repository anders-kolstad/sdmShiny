# Fitting SDM models


# This script attempt to fit SDM models for multiple red-listed Norwegian species.

# This code chunk will get us up to speed.

library(sdm)
source("./R/spList.R")
myIVs      <- raster::stack('data/IV.grd')
mySpecies  <- sl()
oDat2       <- readRDS('data/large/oDat.RData')   # 2 species only as in the Rdm documentation file
oDat       <- readRDS('data/large/allOccurences.RData') # All species


## SDM-data object

for(i in 1:length(mySpecies)){
  s    <- unique(oDat$species)[i]
  s2   <- paste0(s, "_d")
  d    <- oDat[oDat$species == s,]
  dat  <- sdm::sdmData(species~temp+prec+SoilpH+TundraHerbivores,
                       train = d,
                       predictors = myIVs,
                       bg = list(n=1000, method = "gRandom"))    # change bd
  
  assign(s2, dat)
  sdm::write.sdm(dat, paste0("models/sdmData/", s2), overwrite=TRUE)
  rm(s, s2, d, dat)
}

```
#This takes the names over to the gbif approved names (eg Kobresia to Carex). 




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
for(i in 1:length(mySpecies3)){
  s      <- unique(oDat$species)[i]
  s2     <- paste0(s, "_m")
  file   <- paste0("temp/models/sdmModels/", s2, ".sdm")
  
  assign(s2, sdm::read.sdm(file))
}




### Ensemble
#Let's put the 5 x 3 models together and make a map of the current habitat suitability, using ensamble. 

for(i in 1:length(mySpecies)){
  s <- unique(oDat$species)[i]
  d <- get(paste0(s, "_m"))
  
  # using this because overwrite=TRUE dont work
  fn <- paste0("models/predictions/", s, "_ens.img")
  if(file.exists(fn)){file.remove(fn)}
  
  mod <- sdm::ensemble(d,
                       newdata = myIVs, 
                       filename = fn, overwrite=TRUE,   
                       setting = list(method='weighted', stat = 'AUC')
  )     
  assign(paste0(s, "_ens"), 
         mod)
  
  rm(s, d, fn, mod)
}

# Lets plot them.





## Best candidate model

#To make predictions from the 5 x 3 models we need to run the ensamble function again, this time with altered IVs in the newdata argument. The function is probably too slow to run on the fly. This I need to test on the io server later. We could choose one method, eg maxent, and just use that for all species, and with only one replication. Then we could use the predict function in raster which will be quicker. Or we can look at these 15 models we have generated and choose the best one from there. That would be a bit safer.



```{r}
for(i in 1:length(mySpecies)){
  s <- unique(oDat$species)[i]
  s2 <- paste0(s, "_m")
  s3 <- get(paste0(s, "_ev"))
  file <- paste0("models/sdmModels/", s, "_bcm")
  best <- s3$modelID[s3$AUC == max(s3$AUC)]
  mod <- get(s2)
  mod2 <- mod[[best]]
  
  # save the model object for easy lading in the shiny app.
  sdm::write.sdm(mod2, file, overwrite=TRUE)
  
  
  
  mod3 <- raster::predict(mod2,
                          object = myIVs, 
                          filename = paste0("models/predictions/", s, "_best.img"), 
                          overwrite=TRUE)
  
  #assign(paste0(s, "_best"), 
  #       mod3)
  
  rm(s, s2, s3, file, best, mod, mod2, mod3)
}

## Replicated single method
#This approach chooses one method, maxent, and does five replicated partitionings, and then the raster::predict function to make maps. I'll keep these models to to a speed test later in the actuall Shiny app on the io server.


for(i in 1:length(mySpecies)){
  
  s       <- unique(oDat$species)[i]
  file1   <- paste0("models/sdmData/", s, "_d.sdd")
  obj     <- paste0(s, "_ms")
  obj2    <- paste0(s, "_5maxent")
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
                          filename = paste0("models/predictions/", obj2, ".img"), 
                          overwrite=TRUE)
  assign(obj2, mod2)
  
  ev     <- sdm::getEvaluation(mod, stat = "threshold")
  th     <- mean(ev$threshold)
  s      <- paste0(unique(oDat$species)[i], "_5maxent_pa")
  mod2[] <- ifelse(mod2[] >= th, 1, 0)
  assign(s, mod2)
  
}


#Some runs gives white spots on this map. I think I have traced this problem to sdmData function - the 'number of records' ois not the sum of the pseudo-absences plus unique spatial points in the occurence data. Usually I get a fairly close number when I run the for-loop with an empty environment (eg 800-1060 records for Primula). But if I run it again straight after I get maybe around 200. It's easier to check the coverage when there is just one colour.



## Response curves
#Let make response curve plots for each species and save those as well. These are the mean responses for the 5 x 3 models.

for(i in 1:length(mySpecies3)){
  s      <- unique(oDat$species)[i]
  s2     <- paste0(s, "_m")
  s3     <- paste0("models/rcurves/", s, "_rcurves.png")
  d      <- get(s2)
  d2     <- d
  
  p <- sdm::rcurve(d2, ylab="Habitategnethet\nHabitat suitability", 
                   xlab = "",
                   main = "")
  # not sure how to rename the variables:
  #labs <- c("Temperature", "Precipitation", 
  #           "Soil pH", "Sheep and reindeer")
  # p2 <- p + facet_grid(labeller = labeller(variable = labs))
  
  png(filename = s3,
      width = 480, height = 380, units = "px")
  print(p)
  dev.off()
  
}


## Variable importance
#I would also like to get a plot of the variable importance. This script uses the combined/averaged impirtances from the 5 x 3 models.
# create empty list
varimp <- list()

# list of IVs (different for alpine and forest species). 
# Make sure to put them in the correct order. 
# Find the order with > getVarImp(yourModel, id = 1)@varImportance
IV <- c("temp", "prec", "SoilpH", "TundraHerbivores")

# Create an empty varimp table
df1 <- data.frame(variables = IV,
                  corTest = as.numeric(NA),
                  AUCtest = as.numeric(NA))

for(i in 1:length(mySpecies3)){
  s      <- unique(oDat$species)[i]
  s2     <- paste0(s, "_m")
  d      <- get(s2)
  tab    <- d@run.info
  
  for(t in 1:max(tab$modelID)){
    r <- length(varimp)+1
    ifelse(tab$success[t]==TRUE,
           varimp[[r]]          <-sdm::getVarImp(d,id=t)@varImportance,
           varimp[[r]]          <- df1)
    
    varimp[[r]]$species  <-tab$species[t]
    varimp[[r]]$method   <-tab$method[t]
    varimp[[r]]$repid    <-tab$replicationID[t]      
    
    #if(tab$success[t]==FALSE) return(print(paste('Model failiure run ',t)))
  }
}

varimp<-do.call('rbind',varimp)
rm(s, s2, d, tab, df1)
source("R/se.R")
varimpmean <- aggregate(data = varimp,
                        corTest ~ species + variables,
                        FUN = function(x) c(mean = mean(x, na.rm=T), se = se(x)))
varimpmean <- do.call(data.frame, varimpmean)
head(varimpmean)


#Changing the variable names for axis tick labels
varimpmean$variables <- plyr::revalue(varimpmean$variables, c(
  prec         = "Årlig nedbørsmengde\n
  Annual precipitation",
  SoilpH       = "pH i jorden\nSoil pH",
  temp         = "Gjennomsnittemperatur i varmeste kvartal\n
  Mean temperature in warmenst quarter",
  TundraHerbivores = "Tetthet av sau og reinsdyr\n
  Sheep and reindeer densities"))


#We can make one plot per species and show it in the app. That way the user will know what slider should induce the biggest effect.
library(ggplot2)

for(i in 1:length(unique(varimpmean$species))){
  
  d <- varimpmean[varimpmean$species==unique(varimpmean$species)[i],]
  s <- paste0("temp/models/varimp/",
              unique(varimpmean$species)[i],
              ".png")
  
  p <- ggplot2::ggplot(data = d)+
    geom_bar(aes(y = corTest.mean, 
                 x = variables, 
                 fill = variables), 
             stat = "identity", 
             colour = "black")+
    coord_flip()+
    ylab("Variabelviktigheten\nVariable importance")+
    xlab("")+
    theme_minimal()+
    theme(axis.text.y = element_text(size = 10))+
    theme(legend.position="none")+
    geom_errorbar(aes(x = variables, ymin=corTest.mean-corTest.se, ymax=corTest.mean+corTest.se), width=.2)
  
  png(filename = s,
      width = 480, height = 380, units = "px")
  print(p)
  dev.off()
}
