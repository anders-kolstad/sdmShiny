# Run setwd to test locally, but don't include in server version
#setwd("shiny/")

#Packages ####

library(shiny)
library(sdm)
library(sp)
library(raster)
library(dismo)
library(class)
library(mgcv)
library(nlme)
library(Matrix)
library(earth)
library(Formula)
library(plotrix)
library(rJava)
library(Rcpp)
library(rasterVis)
library(gridExtra)
library(rgdal)
library(shinyjs)
library(lattice)
library(latticeExtra)
library(shiny.i18n)
library(ggplot2)




# IVs ####
IV <- raster::stack("IVapp.grd")

# Colours  ####
cols <- colorRampPalette(c("beige", "darkgreen" ))
cols2 <- colorRampPalette(c("red", "white", "darkgreen" ))



shinyServer(
  function(input, output, session) {
    
    # ConditionalPanel Outputs ####
    
    output$cond1 <- renderText({
        con <- data.frame(c1 = as.character(NA), c2 = as.character(NA))
        ss <- input$species
        ss2 <- sub(' ', '_', ss)
      # Get model object to know what IVs to include as sliders
        m       <- read.sdm(paste0("sdmModels/", ss2, "_3gams.sdm"))
      if('TundraHerbivores' %in% names(m@data@features)){print("yes")}
    })
    output$cond2 <- renderText({
      con <- data.frame(c1 = as.character(NA), c2 = as.character(NA))
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      # Get model object to know what IVs to include as sliders
      m       <- read.sdm(paste0("sdmModels/", ss2, "_3gams.sdm"))
      if('moose1999' %in% names(m@data@features)){print("yes")}
    })
    
    
    
    # translated ui ####
    isolate(
    output$contr <- renderUI({
      i18n$set_translation_language(input$lan)
      
     
      box(
        shinyjs::useShinyjs(),
        id = "side-panel",
        title = i18n$t("Controls"), 
        status = "primary", 
        solidHeader = TRUE,
        width = NULL,
          
          sliderInput("temperature", 
                      label = i18n$t("Mean summer temperature (\u00B0C)"),
                      min = -5, max = 5, value = c(0)),
          sliderInput("precipitation", 
                      label = i18n$t("Annual precipitation (%):"),
                      min = -50, max = 50, value = c(0)),
          conditionalPanel(
            condition = "output.cond1 == 'yes'",
          sliderInput("herbivory", 
                      label = i18n$t("Sheep and reindeer density (%):"),
                      min = -50, max = 50, value = c(0))),
          conditionalPanel(
            condition = "output.cond2 == 'yes'",
            sliderInput("moose", 
                        label = i18n$t("Moose density (%):"),
                        min = -50, max = 50, value = c(0)),
            sliderInput("deer", 
                        label = i18n$t("Red deer density (%):"),
                        min = -50, max = 50, value = c(0))),
            sliderInput("roedeer", 
                        label = i18n$t("Roe deer density (%):"),
                        min = -50, max = 50, value = c(0)),
          actionButton("reset_input", i18n$t("Reset inputs"))
      ) 
    } )
    )
    
    
# Top banner ####    
    output$top <- renderUI({
      i18n$set_translation_language(input$lan)
      box(
        width = NULL, background = "yellow",
        i18n$t("Make changes to the climate and herbivore density variables to the right and see how that affects the habitat suitability of your selected plant species. You can look at the 'variable importance' at the bottom right of the page to see which variables are having the biggest effect for this species"))})

    
# Bottom banners ####    
    output$credits <- renderText({
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      fn <- list.files("www/", pattern = ss2)
      if(length(fn)>0){
        creds <- stringr::str_split(fn, "-")[[1]]
        creds <- stringr::str_split(creds[2], "\\.")[[1]]
        creds2 <- stringr::str_split(creds[1], "\\|")[[1]]
        print(paste(creds2, collapse = " "))
      } else{
        print("There is no picture for this species. Enjoy this view over Muddus National Park instead")
      }
    })  # renderText
    
    output$credUI <- renderUI({
      i18n$set_translation_language(input$lan)
      box(
        width = 2, height = 200, background = "aqua",
        i18n$t("Photo credits:"),
        textOutput("credits"))})
   
    output$usUI <- renderUI({
      i18n$set_translation_language(input$lan)
      box(
        width = 4, height = 200, background = "aqua",
        i18n$t("This Shyny app was made by Anders Lorentzen Kolstad and James D. M. Speed, both at the NTNU Univerity Museum, Trondheim, Norway"))})
    
    output$dissclaimer <- renderUI({
      i18n$set_translation_language(input$lan)
      box(
        width = 4, height = 200, background = "aqua",
        i18n$t("The habitat suitability models presented here have had to be made quite simple and are therefore not the most accurate. See the following  published papers for more precise predictions:"))
                     })
    
    output$refs <- renderUI({
    box(i18n$t("References"),width = 2, height = 200, background = "aqua",
        tags$a(href="https://www.sciencedirect.com/science/article/abs/pii/S0006320716309168", "Paper 1"))})
    
    
   
    
    
    
    # Response curves ####
     output$rcurves <- renderPlot({
       
       
       ss2     <- sub(' ', '_', input$species)
       # Get model object 
       m       <- read.sdm(paste0("sdmModels/", ss2, "_3gams.sdm"))
       p       <- sdm::rcurve(m, ylab="Habitategnethet\nHabitat suitability", 
                        xlab = "",
                        main = "")
       # not sure how to rename the variables:
       #labs <- c("Temperature", "Precipitation", 
       #           "Soil pH", "Sheep and reindeer")
       # p2 <- p + facet_grid(labeller = labeller(variable = labs))
      
       print(p)
       
     })
#    output$rcurves <- renderImage({
#      ss <- input$species
#      ss2 <- sub(' ', '_', ss)
#      outfile1 <- paste0("models/rcurves/", ss2, "_rcurves.png")
#      
#      list(src = outfile1,
#           contentType = 'image/png',
#           width = 800,
#           height = 600,
#           alt = "Responce curves")
#      
#      
#    }, deleteFile = FALSE) #renderImage
    
  
    # Variable importance ####
    
    output$varimp <- renderPlot({
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      df1 <- data.frame(variables = names(IV),
                        corTest = as.numeric(NA),
                        AUCtest = as.numeric(NA))
      varimp <- list()
      
        d      <- read.sdm(paste0("sdmModels/", ss2, "_3gams.sdm"))
        tab    <- d@run.info
        
        for(t in 1:max(tab$modelID)){
          r <- length(varimp)+1
          ifelse(tab$success[t]==TRUE,
                 varimp[[r]]          <-sdm::getVarImp(d,id=t)@varImportance,
                 varimp[[r]]          <- df1)
          
          varimp[[r]]$species  <-tab$species[t]
          #varimp[[r]]$method   <-tab$method[t]
          #varimp[[r]]$repid    <-tab$replicationID[t]      
          #if(tab$success[t]==FALSE) return(print(paste('Model failiure run ',t)))
        }
      
      varimp<-do.call('rbind',varimp)
      #rm(s, s2, d, tab, df1)
      #source("R/se.R")
      varimpmean <- aggregate(data = varimp,
                              corTest ~ species + variables,
                              FUN = mean, na.rm=T)
      varimpmean <- do.call(data.frame, varimpmean)
     
      #Changing the variable names for axis tick labels
 #     varimpmean$variables <- plyr::revalue(varimpmean$variables, c(
 #       prec         = "Årlig nedbørsmengde\n
 #       Annual precipitation",
 #       SoilpH       = "pH i jorden\nSoil pH",
 #       temp         = "Gjennomsnittemperatur i varmeste kvartal\n
 #       Mean temperature in warmenst quarter",
 #       TundraHerbivores = "Tetthet av sau og reinsdyr\n
 #       Sheep and reindeer densities"))
        
        p <- ggplot2::ggplot(data = varimpmean)+
          geom_bar(aes(y = corTest, 
                       x = variables, 
                       fill = variables), 
                   stat = "identity", 
                   colour = "black")+
          coord_flip()+
          ylab("Variabelviktigheten\nVariable importance")+
          xlab("")+
          theme_minimal()+
          theme(axis.text.y = element_text(size = 10))+
          theme(legend.position="none")
          
        
        print(p)
       
      })

    
    
    output$varimptext <- renderText({
      ss <- input$species
      paste("This figure shows which variables are most important, according to this model, in determening the distribution of", ss)
    })  # renderText
    
    
    # Pictures ####
    
    output$pic <- renderImage({
      myWidth  <- session$clientData$output_pic_width
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      fn <- list.files("www/", pattern = ss2)
      
      ifelse(length(fn)>0, outfile3 <- paste0("www/", fn), outfile3 <- "www/example.png")
      
      list(src = outfile3,
           contentType = 'image/png',
           width = myWidth,
           height = 'auto',
           alt = fn)
      
      
    }, deleteFile = FALSE) #renderImage
    
   
    
    
    
    # AUC ####
    output$AUC <- renderText({
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      m       <- read.sdm(paste0("sdmModels/", ss2, "_3gams.sdm"))
      
      print(round(mean(getEvaluation(m)[,2], na.rm=T), digits = 2))
      })
    
    # N = ####
    output$n <- renderText({
      ss <- input$species
      ss2 <- sub(' ', '_', ss)
      m       <- read.sdm(paste0("sdmModels/", ss2, "_3gams.sdm"))
      t <- m@data@species[[1]]
      length(t@presence)
    })
    
    # PLOTS ####
    output$map <- renderPlot({
      
      #modify IVs
      #with selectable parameters
      IV2                  <- IV
      IV2$TundraHerbivores <- IV2$TundraHerbivores * (input$herbivory/100 + 1)
      IV2$temp             <- IV2$temp             + input$temperature         
      IV2$prec             <- IV2$prec             * (input$precipitation/100+1)
      
      
      withProgress(message = 
                     "Working on it ... please wait" , value = 0, {
                       
                       # Get species name  
                       ss <- input$species
                       ss2 <- sub(' ', '_', ss)
                       
                       # Predictions
                       # Get model object 
                       m       <- read.sdm(paste0("sdmModels/", ss2, "_3gams.sdm"))
                       
                       # -get img with predicted current habitat suitability
                       # current <- stack(paste0("predictions/", ss2, "_bcm.img"))
                       
                       current <- predict(m, IV,
                                          filename = "predictions/current.img",
                                          overwrite=TRUE,
                                          mean=TRUE)
                       
                       # -make new predictions based on user input
                       
                       # OBS add moose, deer etc....
                       if(input$herbivory == 0 & input$temperature == 0 & input$precipitation == 0){
                         pred <- current
                       } else{
                         pred    <- predict(m, IV2,
                                            filename = "predictions/pred.img",
                                            overwrite=TRUE,
                                            mean=TRUE)
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
    # Help ####
    
    output$help <- renderText({
      
      paste("
Norsk\n
Begynn med å velge deg en art i nedtrekksmenyen til venstre. Disse plantene er, eller har nylig vært, på den norske rødlista. Det vil si de er truede arter. Det venstre kartet viser habitat-egnetheten til denne arten rundt om i Norge. Dette er områder hvor modellene våre tilsier at det finnes egnede vekstforhold. Nå kan du prøve å endre på noen av miljøparametrene ved å flytte på skyvelinjalene. Kartet til høyre vil etter en stund oppdatere seg og vise hvordan dette vil endre habitat-egnetheten for denne arten. Under fanen 'Variable importance' kan du se hvilke av de aktuelle variablene som er viktigst i å bestemme utbredelsen til arten du har valgt, og under 'Response curves' kan du se hvordan det faktiske foholdet er mellom utbredelsen til artene og hver av miljøvariablene. Lykke til.

English\n
Start by selecting a species in the dropdown meny in the left. These are (currently or previously) red-listed plant species in Norway. The left map indicates the habitat suitability for this species across the country. These are areas where or models suggest that suitable habitat exists. Now, you can try make some changes to some og the key varaiblesthat determine habitat suitability for plants by adjusting the sliders on the left. The map on the right will soon update (it might take a few seconds) and show how the habitat suitability is predicted to change in 'your scenario'. Under the pane 'Varaible importance' you can see which of the included variables are most important in determening the distribution of the species you have chosen, and under 'Response curvs' you can even see the actuall relationship between the distribution of the species and each of the environmental variables. Good luck!\n\n
         

            ")
    })  # renderText
    
    observeEvent(input$reset_input, {
      shinyjs::reset("side-panel")
    })
    outputOptions(output, "cond2", suspendWhenHidden = FALSE)
    outputOptions(output, "cond1", suspendWhenHidden = FALSE)
    
          }) # app


