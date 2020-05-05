library(shiny)

require(sdm)
require(rasterVis)
require(gridExtra)
require(rgdal)
require(raster)

# Possible error:
# Make sure that all necessary packages are installed for the same user set under run_as in your Shiny Server configuration file.
# https://support.rstudio.com/hc/en-us/articles/231249288-Why-does-my-app-work-locally-but-not-on-my-Shiny-Server-


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




shinyServer(
  function(input, output) {
    
    
    # Bring all models into the environment
    files <- list.files("./models/sdmModels/", pattern = ".sdm")
    for(i in 1:length(files)){
      assign(files[i], read.sdm(paste0("./models/sdmModels/", files[i])))
    }
    # Possible error - assigning string to variable. Use get() for each time?
    # https://stackoverflow.com/questions/37190516/shiny-dowloadhandler-error-replacement-has-length-zero
    
    # Get IVs
    IV <- raster::stack("./models/IVapp.grd")
    
    
    #modify IVs with selectable parameters
    output$map <- renderPlot({
      
      IV2                  <- IV
      IV2$TundraHerbivores <-IV2$TundraHerbivores * (input$herbivory/100 + 1)
      IV2$temp             <-IV2$temp             + input$temperature         
      IV2$prec             <-IV2$prec             + input$precipitation
      
      
      
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      
      if(input$`modelling approach` == "Ensemble"){
        m <- read.sdm(paste0("./models/sdmModels/", ss2, "_m.sdm"))
        
        current <- raster::stack(paste0("./models/predictions/", ss2, "_ens.img"))
        
        if(input$herbivory == 0 & input$temperature == 0 & input$precipitation == 0){
          pred <- current} else{
          
          fn <- paste0("./models/predictions/app/", ss2, "_ens2.img")
          if(file.exists(fn)){file.remove(fn)}
          fn <- paste0("./models/predictions/app/", ss2, "_ens2.img.aux.xml")
          if(file.exists(fn)){file.remove(fn)}
          
          pred    <- sdm::ensemble(m, newdata = IV2, 
                                 filename = paste0("./models/predictions/app/", ss2, "_ens2.img"), 
                                 overwrite=TRUE,   
                                 setting = list(method='weighted', stat = 'AUC'))
          }
      } 
      
      
      if(input$`modelling approach` == "Replicated maxent (n=5)"){
        
        m       <- read.sdm(paste0("./models/sdmModels/", ss2, "_5maxent.sdm"))
        current <- stack(paste0("./models/predictions/", ss2, "_5maxent.img"))
        if(input$herbivory == 0 & input$temperature == 0 & input$precipitation == 0){
          pred <- current
        } else{
          pred    <- predict(m, IV2, mean=TRUE,
                                     filename = "./models/predictions/app/predHSM.img",
                                     overwrite=TRUE)
        }
      }
      
      
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
      
    })})

