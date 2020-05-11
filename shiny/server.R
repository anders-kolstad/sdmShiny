library(shiny)
library(sdm)
library(rasterVis)
library(gridExtra)
library(rgdal)


#library(raster)
#library(grDevices)
#library(dismo)


# Independenet variables
IV <- raster::stack("./models/IVapp.grd")

# Define colour palette
cols <- colorRampPalette(c("beige", "darkgreen" ))
cols2 <- colorRampPalette(c("red", "white", "darkgreen" ))



shinyServer(
  function(input, output) {
    
    # PLOTS ####
    output$map <- renderPlot({
   
      #modify IVs
      #with selectable parameters
        IV2                  <- IV
        IV2$TundraHerbivores <-IV2$TundraHerbivores * (input$herbivory/100 + 1)
        IV2$temp             <-IV2$temp             + input$temperature         
        IV2$prec             <-IV2$prec             * (input$precipitation/100+1)

        
        withProgress(message = 
            "Working on it ... please wait" , value = 0, {
                
      # Get species name  
        ss <- input$species
        ss2 <- sub(' ', '_', ss)

      # Predictions
        # Get model object (best (out of 15) candidate model)
        m       <- read.sdm(paste0("models/sdmModels/", ss2, "_bcm.sdm"))
        
      # -get img with predicted current habitat suitability
        current <- stack(paste0("models/predictions/", ss2, "_bcm.img"))
        
      # -make new predictions based o user input
        if(input$herbivory == 0 & input$temperature == 0 & input$precipitation == 0){
          pred <- current
        } else{
          pred    <- predict(m, IV2, mean=TRUE,
                             filename = "models/predictions/app/predicted.img",
                             overwrite=TRUE)
        }
        diff <- pred-current
        
      # Plot  
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
        
        }) # progress
        
        return(print(predp))

         }) #renderPlot
    
    
    
    # Response curves ####
   
    output$rcurves <- renderImage({
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      outfile1 <- paste0("models/rcurves/", ss2, "_rcurves.png")
      
      list(src = outfile1,
           contentType = 'image/png',
           width = 800,
           height = 600,
           alt = "Responce curves")
      
      
    }, deleteFile = FALSE) #renderImage
    
  
    # Variable importance ####
    
    output$varimp <- renderImage({
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      outfile2 <- paste0("models/varimp/", ss2, ".png")
      # here we can generate on the fly perhaps.
      
      
      list(src = outfile2,
           contentType = 'image/png',
           width = 800,
           height = 600,
           alt = "Variable importance")
      
      
    }, deleteFile = FALSE) #renderImage
    
    
    output$varimptext <- renderText({
      ss <- input$species
      paste("This figure shows which variables are most important, according to this mode, in determening the distribution of", ss)
    })  # renderText
    
    
    # Pitures ####
    
    output$pic <- renderImage({
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      fn <- list.files("www/", pattern = ss2)
      outfile3 <- paste0("www/", fn)
      
      list(src = outfile3,
           contentType = 'image/png',
           width = 400,
           height = 300,
           alt = fn)
      
      
    }, deleteFile = FALSE) #renderImage
    
    output$credits <- renderText({
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      fn <- list.files("www/", pattern = ss2)
      creds <- stringr::str_split(fn, "-")[[1]]
      creds <- stringr::str_split(creds[2], "\\.")[[1]]
      creds2 <- stringr::str_split(creds[1], "\\|")[[1]]
      
      print(paste("Photo credits:", paste(creds2, collapse = " ")))
      
      
    })  # renderText
    
    
    # Help ####
    
    output$help <- renderText({
      
      paste("
Norsk\n
Begynn med å velge deg en art i nedtrekksmenyen til venstre. Disse plantene er, eller har nylig vært, på den norske rødlista. Det vil si de er truede arter. Det venstre kartet viser habitat-egnetheten til denne arten rundt om i Norge. Dette er områder hvor modellene våre tilsier at det finnes egnede vekstforhold. Nå kan du prøve å endre på noen av miljøparametrene ved å flytte på skyvelinjalene. Kartet til høyre vil etter en stund oppdatere seg og vise hvordan dette vil endre habitat-egnetheten for denne arten. Under fanen 'Variable importance' kan du se hvilke av de aktuelle variablene som er viktigst i å bestemme utbredelsen til arten du har valgt, og under 'Response curves' kan du se hvordan det faktiske foholdet er mellom utbredelsen til artene og hver av miljøvariablene. Lykke til.

English\n
Start by selecting a species in the dropdown meny in the left. These are (currently or previously) red-listed plant species in Norway. The left map indicates the habitat suitability for this species across the country. These are areas where or models suggest that suitable habitat exists. Now, you can try make some changes to some og the key varaiblesthat determine habitat suitability for plants by adjusting the sliders on the left. The map on the right will soon update (it might take a few seconds) and show how the habitat suitability is predicted to change in 'your scenario'. Under the pane 'Varaible importance' you can see which of the included variables are most important in determening the distribution of the species you have chosen, and under 'Response curvs' you can even see the actuall relationship between the distribution of the species and each of the environmental variables. Good luck!\n\n
         

            ")
    })  # renderText
    
    
          }) # app


