# Load packages ####

library(shiny)
library(sdm)
library(rasterVis)
library(gridExtra)
library(rgdal)
library(raster)
library(grDevices)
# Possible error:
# Make sure that all necessary packages are installed for the same user set under run_as in your Shiny Server configuration file.
# https://support.rstudio.com/hc/en-us/articles/231249288-Why-does-my-app-work-locally-but-not-on-my-Shiny-Server-




# * * * *  INFO * * * * ####

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


# * * * *  Start * * * * #####

# import species list
#s <- read.csv('data/species_list.csv')
#s <- as.vector(s$x)   # ADB names

# Get all model object
# saves time to have this outside the server function
#fn <- list.files("models/sdmModels/", pattern = ".sdm")
#for(i in 1:length(fn)){
#  d <- read.sdm(paste0("models/sdmModels/", fn[i]))
#  assign(fn[i], d)
#}


shinyServer(
  function(input, output) {
    
    
    # Get IVs ####
    IV <- raster::stack("./models/IVapp.grd")
    cat(stderr(), "load raster")
    
    
    #modify IVs ####
    #with selectable parameters
    output$map <- renderPlot({
      
      IV2                  <- IV
      IV2$TundraHerbivores <-IV2$TundraHerbivores * (input$herbivory/100 + 1)
      IV2$temp             <-IV2$temp             + input$temperature         
      IV2$prec             <-IV2$prec             + input$precipitation
      
      cat(stderr(), "change IV")
      
      # Get species name  ####
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      cat(stderr(), "alter input species names")
      
      # ensemble    ####
      if(input$`modelling approach` == "Ensemble"){
        m <- sdm::read.sdm(paste0("./models/sdmModels/", ss2, "_m.sdm"))
        cat(stderr(), "read ensemble model")
        current <- raster::stack(paste0("./models/predictions/", ss2, "_ens.img"))
        cat(stderr(), "stack current")
        
        if(input$herbivory == 0 & input$temperature == 0 & input$precipitation == 0){
          pred <- current} else{
            
            fn <- paste0("./models/predictions/app/", ss2, "_ens2.img")
            if(file.exists(fn)){file.remove(fn)}
            fn2 <- paste0("./models/predictions/app/", ss2, "_ens2.img.aux.xml")
            if(file.exists(fn2)){file.remove(fn2)}
            
            pred    <- sdm::ensemble(m, newdata = IV2, 
                                     filename = paste0("./models/predictions/app/", ss2, "_ens2.img"), 
                                     overwrite=TRUE,   
                                     setting = list(method='weighted', stat = 'AUC'))
          }
      } 
      
      
      
      # Maxent ####
      if(input$`modelling approach` == "Replicated maxent (n=5)"){
        
        m       <- sdm::read.sdm(paste0("./models/sdmModels/", ss2, "_5maxent.sdm"))
        current <- stack(paste0("./models/predictions/", ss2, "_5maxent.img"))
        if(input$herbivory == 0 & input$temperature == 0 & input$precipitation == 0){
          pred <- current
        } else{
          pred    <- raster::predict(m, IV2, mean=TRUE,
                                     filename = "./models/predictions/app/predHSM.img",
                                     overwrite=TRUE)
        }
      }
      
      
      
      # Best candidate ####
      if(input$`modelling approach` == "Best candidate model"){
        
        m       <- read.sdm(paste0("models/sdmModels/", ss2, "_bcm.sdm"))
        current <- stack(paste0("models/predictions/", ss2, "_bcm.img"))
        if(input$herbivory == 0 & input$temperature == 0 & input$precipitation == 0){
          pred <- current
        } else{
          
          pred    <- predict(m, IV2, mean=TRUE,
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
      
      
      
      
    }) # render plot
    
    output$rcurves <- renderImage({
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      outfile1 <- paste0("models/rcurves/", ss2, "_rcurves.png")
      
      list(src = outfile1,
           contentType = 'image/png',
           width = 400,
           height = 300,
           alt = "Responce curves")
      
      
    }, deleteFile = FALSE) #renderImage
    
    
    
    output$varimp <- renderImage({
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      outfile2 <- paste0("models/varimp/", ss2, ".png")
      # here we can generate on the fly perhaps.
      
      
      list(src = outfile2,
           contentType = 'image/png',
           width = 400,
           height = 300,
           alt = "Variable importance")
      
      
    }, deleteFile = FALSE) #renderImage
    
    output$varimptext <- renderText({
      ss <- input$species
      paste("This figure shows which variables are most important, according to this mode, in determening the distribution of", ss)
    })
    
    })  # app

