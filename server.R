library(shiny)
library(raster)
library(sdm)
library(rasterVis)
library(gridExtra)
library(rgdal)
library(grDevices)
library(dismo)



IV <- raster::stack("./models/IVapp.grd")
Kobresia_new.sdm <- read.sdm("models/sdmModels/Kobresia_new.sdm")
p <- raster::predict(Kobresia_new.sdm, IV, 
                     mean=TRUE, 
                     filename  = "models/predictions/app/Kobresia.img",
                     overwrite=TRUE)

shinyServer(
  function(input, output) {
    output$map <- renderPlot({
   
    plot(p)

         }) #renderPlot
    
    
   
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
          }) # app


