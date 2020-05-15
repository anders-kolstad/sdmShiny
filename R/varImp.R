## Variable importance



# Get a list of species for each IV combo
source("./R/spList.R")
mySpecies2  <- sl(df=T)
myAlpine <- mySpecies2$mySpList[mySpecies2$type == "alpine"]
myForest <- mySpecies2$mySpList[mySpecies2$type == "forest"]
myAlpine <- sub(' ', '_', myAlpine)
myForest <- sub(' ', '_', myForest)

comb <- c(as.character(myAlpine), as.character(myForest))
(comb <- comb[which(duplicated(comb))])
myAlpine <- myAlpine[!myAlpine %in% comb]
myForest <- myForest[!myForest %in% comb]

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
mySpecies <- unique(as.character(myS2))
#length(mySpecies) # 85

comb2 <- mySpecies[mySpecies %in% comb]
myAlpine2 <- mySpecies[mySpecies %in% myAlpine]
myForest2 <- mySpecies[mySpecies %in% myForest]

#mySpecies  <- c("Primula scandinavica", "Kobresia simpliciuscula")
#mySpecies  <- sl()
#oDat       <- readRDS('data/large/oDat.RData')   # 2 species only as in the Rdm documentation file
# oDat2       <- readRDS('data/large/allOccurences.RData') # All species
# oDat2 <- oDat



# list of IVs (different for alpine and forest species). 
# Make sure to put them in the correct order. 
# Find the order with > getVarImp(yourModel, id = 1)@varImportance
IV_alpine <- c("temp", "prec", "SoilpH", "TundraHerbivores")
IV_forest <- c("temp", "prec", "SoilpH", "moose1999", "red_deer1999", "roe_deer1999")

#x      <- read.sdm("shiny/sdmModels/Allium_scorodoprasum_bcm.sdm")
#getVarImp(x)@varImportance


# Create an empty varimp table
df1 <- data.frame(variables = IV_forest,
                  corTest = as.numeric(NA),
                  AUCtest = as.numeric(NA))


# create empty list
varimp <- list()

for(i in 1:length(myForest2)){
  s      <- unique(myForest2)[i]
  s2     <- paste0(s, "_bcm")
  d      <- read.sdm(paste0("shiny/sdmModels/", s2, ".sdm"))
  tab    <- d@run.info
  varimp[[i]]          <-sdm::getVarImp(d)@varImportance 
    varimp[[i]]$species  <-tab$species
    varimp[[i]]$method   <-tab$method
    #varimp[[i]]$repid    <-tab$replicationID[i]      
    
    #if(tab$success[t]==FALSE) return(print(paste('Model failiure run ',t)))
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
  moose1999 = "Tetthet av elg\n
  Moose densities",
  red_deer1999 = "Tetthet av hjort\nDeer densities",
  roe_deer1999 = "Tetthet av rådyr\nRoe deer densities"))
##Changing the variable names for axis tick labels
#varimpmean$variables <- plyr::revalue(varimpmean$variables, c(
#  prec         = "Årlig nedbørsmengde\n
#  Annual precipitation",
#  SoilpH       = "pH i jorden\nSoil pH",
#  temp         = "Gjennomsnittemperatur i varmeste kvartal\n
#  Mean temperature in warmenst quarter",
#  TundraHerbivores = "Tetthet av sau og reinsdyr\n
#  Sheep and reindeer densities"))


#We can make one plot per species and show it in the app. That way the user will know what slider should induce the biggest effect.
library(ggplot2)

for(i in 1:length(unique(varimpmean$species))){
  
  d <- varimpmean[varimpmean$species==unique(varimpmean$species)[i],]
  s <- paste0("shiny/varimp/",         
              unique(varimpmean$species)[i],
              ".png")
  #s <- paste0("models/varimp/",         
  #            unique(varimpmean$species)[i],
  #            ".png")
  
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

#***************************************END **********************************⋅####
# This script works for getting var imp from the 5*3 models
# Create an empty varimp table
df1 <- data.frame(variables = IV,
                  corTest = as.numeric(NA),
                  AUCtest = as.numeric(NA))




for(i in 1:length(mySpecies)){
  s      <- unique(oDat$species)[i]
  s2     <- paste0(s, "_m")
  d      <- read.sdm(paste0("models/sdmModels/", s2, ".sdm"))
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
  #s <- paste0("models/varimp/",         
  #            unique(varimpmean$species)[i],
  #            ".png")
  
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
