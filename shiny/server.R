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
    
    # prepare species name ####
  theName <- reactive({
      myName <- sub(' ', '_', input$species)
      return(myName)
    })
    
    # Get model ####
  theMod <- reactive({
      myName <- theName()
      m       <- read.sdm(paste0("sdmModels/", myName, "_3gams.sdm"))
      return(m)
    })

# New IVs based on user input
  
  IV2 <- reactive({
    #modify IVs ####
    #with selectable parameters
    IV2                  <- IV
    IV2$TundraHerbivores <- IV2$TundraHerbivores * (input$herbivory/100 + 1)
    IV2$moose1999        <- IV2$moose1999        * (input$moose/100 + 1)
    IV2$red_deer1999     <- IV2$red_deer1999     * (input$deer/100 + 1)
    IV2$roe_deer1999     <- IV2$roe_deer1999     * (input$roedeer/100 + 1)
    IV2$SoilpH           <- IV2$SoilpH           +  input$ph
    IV2$temp             <- IV2$temp             + input$temperature         
    IV2$prec             <- IV2$prec             * (input$precipitation/100+1)
    return(IV2)
  })
  
#ConditionalPanel Outputs ####
    
    output$cond1 <- renderText({
      m <- theMod()
              if('TundraHerbivores' %in% names(m@data@features)){print("yes")}
    })
    
    
    output$cond2 <- renderText({
      m       <- theMod()
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
                      label = i18n$t("Mean summer temperature (\u00B0C):"),
                      min = -5, max = 5, value = c(0)),
          sliderInput("precipitation", 
                      label = i18n$t("Annual precipitation (%):"),
                      min = -50, max = 50, value = c(0)),
        sliderInput("ph", 
                    label = i18n$t("Soil pH (in pH units):"),
                    min = -0.6, max = 0.6, value = c(0)),
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
                        min = -50, max = 50, value = c(0)),
            sliderInput("roedeer", 
                        label = i18n$t("Roe deer density (%):"),
                        min = -50, max = 50, value = c(0))),
            
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
      i18n$set_translation_language(input$lan)
      fn <- list.files("www/", pattern = theName())
      if(length(fn)>0){
        creds <- stringr::str_split(fn, "-")[[1]]
        creds <- stringr::str_split(creds[2], "\\.")[[1]]
        creds2 <- stringr::str_split(creds[1], "\\|")[[1]]
        print(paste(creds2, collapse = " "))
      } else{
        print(i18n$t("There is no picture for this species. Enjoy this view over Muddus National Park instead"))
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
       
       m <- theMod()
       
       p       <- sdm::rcurve(m, ylab="Habitategnethet\nHabitat suitability", 
                        xlab = "",
                        main = "")
       # not sure how to rename the variables:
       #labs <- c("Temperature", "Precipitation", 
       #           "Soil pH", "Sheep and reindeer")
       # p2 <- p + facet_grid(labeller = labeller(variable = labs))
      
       print(p)
       
     })

    
    
    
    
    # Variable importance ####
    
    output$varimp <- renderPlot({
      i18n$set_translation_language(input$lan)
      df1 <- data.frame(variables = names(IV),
                        corTest = as.numeric(NA),
                        AUCtest = as.numeric(NA))
      varimp <- list()
      
        d      <- read.sdm(paste0("sdmModels/", theName(), "_3gams.sdm"))
        tab    <- d@run.info
        
        for(t in 1:max(tab$modelID)){
          r <- length(varimp)+1
          ifelse(tab$success[t]==TRUE,
                 varimp[[r]]          <-sdm::getVarImp(d,id=t)@varImportance,
                 varimp[[r]]          <- df1)
          
          varimp[[r]]$species  <-tab$species[t]
          
        }
      
      varimp<-do.call('rbind',varimp)
    
      varimpmean <- aggregate(data = varimp,
                              corTest ~ species + variables,
                              FUN = mean, na.rm=T)
      varimpmean <- do.call(data.frame, varimpmean)
     
      #Changing the variable names for axis tick labels
      varimpmean$variables <- plyr::revalue(varimpmean$variables, c(
        prec         = i18n$t("Annual precipitation"),
        SoilpH       = i18n$t("Soil pH"),
        temp         = i18n$t("Mean summer temperature"),
        TundraHerbivores = i18n$t("Sheep and reindeer density"),
        moose1999    = i18n$t("Moose density"),
        red_deer1999 = i18n$t("Red deer density"),
        roe_deer1999 = i18n$t("Roe deer density")
        ))
        
        p <- ggplot2::ggplot(data = varimpmean)+
          geom_bar(aes(y = corTest, 
                       x = variables, 
                       fill = variables), 
                   stat = "identity", 
                   colour = "black")+
          coord_flip()+
          ylab(i18n$t("Variable importance"))+
          xlab("")+
          theme_minimal()+
          theme(axis.text.y = element_text(size = 10))+
          theme(legend.position="none")
          
        
        print(p)
       
      })

    
    
    output$varimptext <- renderText({
      i18n$set_translation_language(input$lan)
      # ss <- input$species
      print(i18n$t("This figure shows which variables are most important, according to this model, in determening the distribution of the selected species"))  # DON'T work for some reason
    })  # renderText
    
    
    # Pictures ####
    
    output$pic <- renderImage({
      myWidth  <- session$clientData$output_pic_width
      height <- session$clientData$output_plot2_height
      
      
      fn <- list.files("www/", pattern = theName())
      
      ifelse(length(fn)>0, outfile3 <- paste0("www/", fn), outfile3 <- "www/example.png")
       
      list(src = outfile3,
           contentType = 'image/png',
           #width =   myWidth, 
           height = 300,
           alt = fn)
 
    }, deleteFile = FALSE) #renderImage
    
   
    
    
    
    # AUC ####
    output$AUC <- renderText({
      
      m       <- read.sdm(paste0("sdmModels/", theName(), "_3gams.sdm"))
      
      print(round(mean(getEvaluation(m)[,2], na.rm=T), digits = 2))
      })
    
    # N = ####
    output$n <- renderText({
      m       <- theMod()
      t <- m@data@species[[1]]
      length(t@presence)
    })
    
    # PLOTS ####
    
    current <- reactive({
      
      m   <- theMod()
      withProgress(message =   "Mapping current habitat suitability ... please wait" , value = 0, {
      predict(m, IV,
                           filename = "predictions/current.img",
                           overwrite=TRUE,
                           mean=TRUE)
      })        
    })
    
    

    output$map <- renderPlot({
      
         # -make new predictions based on user input
         IV2                  <- IV2()
         m <- theMod()
         
         if(
           identical(IV2,IV)
           ){
           pred <- current()
         } else{
           withProgress(message = 
    "Mapping predicted changes in habitat suitability ... please wait" , value = 0, {
                        pred    <- predict(m, IV2,
                                           filename = "predictions/pred.img",
                                           overwrite=TRUE,
                                           mean=TRUE)
                          })}
                      
           diff <- pred-current()
           p1 <- rasterVis::levelplot(current(),
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
            
           predp <- grid.arrange(p1, p2, ncol = 2, top = paste(input$species))   
            
           return(print(predp))
           
    }) #renderPlot    
  
                       
                    

 
    
    observeEvent(input$reset_input, {
      shinyjs::reset("side-panel")
    })
    outputOptions(output, "cond2", suspendWhenHidden = FALSE)
    outputOptions(output, "cond1", suspendWhenHidden = FALSE)

    
          }) # app


