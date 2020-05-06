# Response curves

#Let make response curve plots for each species and save those as well. These are the mean responses for the 5 x 3 models.
mySpecies  <- c("Primula scandinavica", "Kobresia simpliciuscula")
#mySpecies  <- sl()
oDat       <- readRDS('data/large/oDat.RData')   # 2 species only as in the Rdm documentation file
# oDat2       <- readRDS('data/large/allOccurences.RData') # All species
# oDat2 <- oDat



for(i in 1:length(mySpecies)){
  s      <- unique(oDat$species)[i]
  s2     <- paste0(s, "_m")
  s3     <- paste0("models/rcurves/", s, "_rcurves.png")
  d      <- read.sdm(paste0("models/sdmModels/", s2, ".sdm"))
  
  
  p <- sdm::rcurve(d, ylab="Habitategnethet\nHabitat suitability", 
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
