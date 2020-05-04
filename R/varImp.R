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