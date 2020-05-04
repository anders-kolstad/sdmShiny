

library(shiny)
library(raster)
library(sdm)
library(rasterVis)
library(gridExtra)


# * * * *  INFO * * * * #

# Model objects:
# 5*3 model = _m
# Replicated single method = _ms
# Best candidate model = _bcm

# Species:
# Primula_scandinavica
# Carex_simpliciuscula

# IVs:
# Temperature
# Precipitation
# Soil pH
# Sheep and reindeer


# * * * *  Start * * * * #

# import species list
#s <- read.csv('data/species_list.csv')
#s <- as.vector(s$x)   # ADB names

# Bring all models into the environment
files <- list.files("models/sdmModels/", pattern = ".sdm")
for(i in 1:length(files)){
  assign(files[i], read.sdm(paste0("models/sdmModels/", files[i])))
  }


# Get IVs
IV <- raster::stack("models/IVapp.grd")



shinyServer(
  function(input, output) {
    
    #modify IVs with selectable parameters
    output$map <- renderPlot({
      
      IV2      <- IV
      IV2$TundraHerbivores <-IV2$TundraHerbivores * (input$herbivory/100 + 1)
      IV2$temp             <-IV2$temp             + input$temperature         
      IV2$prec             <-IV2$prec             + input$precipitation
      
      
      
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      
      if(input$`modelling approach` == "Ensemble"){
        m <- read.sdm(paste0("models/sdmModels/", ss2, "_m.sdm"))
        
        current <- raster::stack(paste0("models/predictions/", ss2, "_ens.img"))
        
        if(input$herbivory == 0 & input$temperature == 0 & input$precipitation == 0){
          pred <- current} else{
          
          fn <- paste0("models/predictions/app/", ss2, "_ens2.img")
          if(file.exists(fn)){file.remove(fn)}
          pred    <- sdm::ensemble(m, newdata = IV2, 
                                 filename = paste0("models/predictions/app/", ss2, "_ens2.img"), 
                                 overwrite=TRUE,   
                                 setting = list(method='weighted', stat = 'AUC'))
          }
      } 
      
      
      if(input$`modelling approach` == "Replicated maxent (n=5)"){
        
        m       <- read.sdm(paste0("models/sdmModels/", ss2, "_5maxent.sdm"))
        current <- stack(paste0("models/predictions/", ss2, "_5maxent.img"))
        if(input$herbivory == 0 & input$temperature == 0 & input$precipitation == 0){
          pred <- current
        } else{
          pred    <- raster::predict(object = m, newdata=IV2, mean=TRUE,
                                     filename = "models/predictions/app/predHSM.img",
                                     overwrite=TRUE)
        }
      }
      
      
      if(input$`modelling approach` == "Best candidate model"){
        
        m       <- read.sdm(paste0("models/sdmModels/", ss2, "_bcm.sdm"))
        current <- stack(paste0("models/predictions/", ss2, "_bcm.img"))
        if(input$herbivory == 0 & input$temperature == 0 & input$precipitation == 0){
          pred <- current
        } else{
        
        pred    <- raster::predict(object = m, newdata=IV2, mean=TRUE,
                                   filename = "models/predictions/app/predHSM.img",
                                   overwrite=TRUE)
        }
      }
      
      cols <- colorRampPalette(c("beige", "darkgreen" ))
      cols2 <- colorRampPalette(c("red", "white", "darkgreen" ))
      
      
      diff <- pred-current
      
     
      
      p1 <- rasterVis::levelplot(current,
                           margin=F,
                           main="Estimated current\nhabitat suitability",
                           col.regions = cols,
                           at=seq(0, 1,len=19),
                           scales=list(draw=FALSE))
      
      p2 <- rasterVis::levelplot(diff,
                           margin=F,
                           main="Predicted change in\nhabitat suitability",
                           col.regions = cols2,
                           at=seq(-1, 1,len=19),
                           scales=list(draw=FALSE))
      predp <- grid.arrange(p1, p2, ncol = 2, top = paste(ss))
      
      
      return(print(predp))
      
    })})

