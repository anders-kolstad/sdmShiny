

library(shiny)
library(raster)
#library(dismo)
library(sdm)
library(rasterVis)



# Model objects
# 5*3 model = _m
# Replicated single method = _ms
# Best candidate model = _bcm

# Species
# Primula_scandinavica
# Carex_simpliciuscula

# IVs
# Temperature
# Precipitation
# Soil pH
# Sheep and reindeer

# import species list
s <- read.csv('data/species_list.csv')
s <- as.vector(s$x)

# Bring all models into the environment
files <- list.files("models/sdmModels/", pattern = ".sdm")
mlist <- list()
for(i in 1:length(files)){
  assign(files[i], read.sdm(paste0("models/sdmModels/", files[i])))
  }


# Get IVs
IV <- raster::stack("data/IV.grd")


shinyServer(
  function(input, output) {
    
    #modify IVs with selectable parameters
    output$map <- renderPlot({
      IV2      <- IV
      IV2$TundraHerbivores <-IV2$TundraHerbivores * (input$herbivory/100 + 1)
      IV2$temp             <-IV2$temp             + input$temperature         # won't work because temp is x10 !!
      IV2$prec             <-IV2$prec             + input$precipitation
      
      #names(IV2)<-names(IV)
      
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      
      if(input$`modelling approach` == "Ensemble"){
        m <- get(paste0(ss2, "_m.sdm"))
        current <- sdm::ensemble(m, newdata = IV, 
                                 filename = paste0("models/predictions/app/", s, "_ens.img"), 
                                 overwrite=TRUE,   
                                 setting = list(method='weighted', stat = 'AUC'))
        pred    <- sdm::ensemble(m, newdata = IV2, 
                                 filename = paste0("models/predictions/app/", s, "_ens2.img"), 
                                 overwrite=TRUE,   
                                 setting = list(method='weighted', stat = 'AUC'))
      } else{
      
      ending <- ifelse(input$`modelling approach` == "Replicated maxent (n=5)", "_ms", "_bcm")
      m <- get(paste0(ss2, ending, ".sdm"))
      
      current  <-  predict(m, IV)
      pred     <-  predict(m, IV2)
      }
      
      predp<-rasterViz::levelplot(raster::stack(p1,p2),
                                  margin=F,
                                  main="Habitat suitability \n Primula scandinavica",
                                  names.attr=c("Current","Changed"))
      return(print(predp))
      
    })})

