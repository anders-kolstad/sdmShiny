# Run setwd to test locally, but don't include in server version
#setwd("shiny/")

#Packages ####

library(shiny)
library(sdm)
library(sp)
library(raster)
library(dismo)
#library(class)
#library(mgcv)
library(nlme)
library(Matrix)
library(earth)
library(Formula)
library(plotrix)
#library(rJava)
#library(Rcpp)
library(rasterVis)
library(gridExtra)
library(rgdal)
library(shinyjs)
library(lattice)
library(latticeExtra)
library(shiny.i18n)
library(ggplot2)
library(readr)
library(mapview)
library(leaflet)
library(plyr)
library(leafsync)
library(shinyWidgets)

# Get species list ####
namelist <- readRDS('namelist.RData') # se R/vernacular.R
namelist$sci        <- as.character(namelist$sci)       # scientific names
namelist$vern_nor   <- as.character(namelist$vern_nor)  # norwegian names
namelist$vern_eng   <- as.character(namelist$vern_eng)  # english names with scientific names as placeholders
namelist$vern_eng2  <- as.character(namelist$vern_eng2) # english names


# IVs ####
IV <- raster::stack("IVapp.grd")


# Occurences
occ <- readRDS("allOccurences.RData")

# Colours  ####
cols <- colorRampPalette(c("beige", "darkgreen" ))
cols2 <- colorRampPalette(c("red", "yellow", "white", "white", "green", "darkgreen" ))
pal <- colorNumeric(c("#FFFFFF", "#157300"), c(0:1), 
                    na.color = "transparent", reverse = F)
pal_rev <- colorNumeric(c("#FFFFFF", "#157300"), c(0:1), 
                        na.color = "transparent", reverse = T)

pal2 <- colorNumeric(c("#FF0000", "#FFFFFF", "#157300"), c(-1:1), 
                     na.color = "transparent", reverse = F)
pal2_rev <- colorNumeric(c("#FF0000", "#FFFFFF", "#157300"), c(-1:1), 
                         na.color = "transparent", reverse = T)




shinyServer(
  function(input, output, session) {
    
    # translated ui ####
    
    output$spLanguage <- renderUI({
    i18n$set_translation_language(input$lan)
      norVal <- c("Vitenskapelige navn", "Norske navn", "Engelske navn")
      engVal <- c("Scientific names", "Norwegian names", "English names")
      vals <- switch(input$lan,
                     "en" = engVal,
                     "no" = norVal)
    radioGroupButtons("spLang",
                 #i18n$t("Show species by:"),
                 choiceNames = vals,
                 choiceValues = c("sci", "vern_nor", "vern_eng"),
                 direction = "vertical")
    })
    
    
    # select species
    output$spNames <- renderUI({
      i18n$set_translation_language(input$lan)
      
      whatNames <- switch(input$spLang,
                          "sci"      = namelist$sci,
                          "vern_nor" = namelist$vern_nor,
                          "vern_eng" = namelist$vern_eng,
                          namelist$sci)
      pickerInput("species",
                   i18n$t("Pick a species"),
                  choices = whatNames )
                  #choiceNames = whatNames,
                   #choiceValues = namelist$sci) # choiceNames and choiceValues are not options for the pickerInput
    })
    
    
    # prepare species name ####
  theName <- reactive({
      temp <- which(namelist == input$species, arr.ind = T)
      row <- temp[1,1]
      myName <- sub(' ', '_', namelist$sci[row])
      return(myName)
    })
    
    # Get model ####
  theMod <- reactive({
      myName <- theName()
      m       <- sdm::read.sdm(paste0("sdmModels/", myName, "_3glms.sdm"))
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
    #IV2$SoilpH           <- IV2$SoilpH           +  input$ph
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
    
    
    

    
    isolate(
    output$contr <- renderUI({
      i18n$set_translation_language(input$lan)
      
     
      box(
        shinyjs::useShinyjs(),
        id = "side-panel",
        title = i18n$t("Control Panel"), 
        status = "primary", 
        solidHeader = TRUE,
        width = NULL,
          
          sliderInput("temperature", 
                      label = i18n$t("Mean summer temperature (\u00B0C):"),
                      min = -5, max = 5, value = c(0)),
          sliderInput("precipitation", 
                      label = i18n$t("Annual precipitation (%):"),
                      min = -50, max = 50, value = c(0)),
        #sliderInput("ph", 
        #            label = i18n$t("Soil pH (in pH units):"),
        #            min = -0.6, max = 0.6, value = c(0)),
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
    
    
    output$myTitle <- renderText({
      i18n$set_translation_language(input$lan)
      i18n$t("Interactive distribution modelling")
      })
    
# Top banner ####    
    output$top <- renderUI({
      i18n$set_translation_language(input$lan)
      box(
        width = NULL, background = "yellow",
        i18n$t("Make changes to the climate and herbivore density variables in the Control Panel and see how that affects the habitat suitability of your selected plant species. You can look at the 'variable importance' below to see which variables are having the biggest effect for this species. Start by playing around, or take 'the Challenge'."))})

    output$selectedSp <- renderText({
      temp <- which(namelist == input$species, arr.ind = T)
      row <- temp[1,1]
      #theSci <- input$species
      paste(
        #theSci, 
        #namelist$vern_nor[match(theSci, namelist$sci)],
        #namelist$vern_eng2[match(theSci, namelist$sci)],
        
        namelist$sci[row],
        namelist$vern_nor[row],
        namelist$vern_eng[row],
        sep = " | "
      )
    })
   
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
        width = NULL,  #background = "aqua",
        i18n$t("Photo credits:"),
        textOutput("credits"))})
   
    output$usUI <- renderUI({
      i18n$set_translation_language(input$lan)
      box(
        width = NULL, #background = "aqua",
        i18n$t("This Shiny app was made by Anders Lorentzen Kolstad and James D. M. Speed, both at the NTNU Univerity Museum, Trondheim, Norway. To see where the data comes from, please visit the"), 
        a("project homepage.", href = "https://anders-kolstad.github.io/sdmShiny/"),
        i18n$t("The habitat suitability models presented here have had to be made quite simple and are therefore not the most accurate. See the following  published papers for more precise predictions:"),
        a("Paper 1", href = "https://www.sciencedirect.com/science/article/abs/pii/S0006320716309168"),
        i18n$t("This work was financed by the Research Council of Norway, project ID 262064.")
        )
       })
    
  
    
    # Response curves ####
     output$rcurves <- renderPlot({
       i18n$set_translation_language(input$lan)
       
       m <- theMod()
       
       p       <- sdm::rcurve(m, ylab=i18n$t("Habitat suitability"), 
                        xlab = "",
                        main = "",
                        size = 100) # defult size   
       
       p$data$variable <- plyr::revalue(p$data$variable, c(
         prec         = i18n$t("Annual\nprecipitation"),
         SoilpH       = i18n$t("Soil pH"),
         temp         = i18n$t("Mean summer\ntemperature"),
         TundraHerbivores = i18n$t("Sheep and\nreindeer density"),
         moose1999    = i18n$t("Moose\ndensity"),
         red_deer1999 = i18n$t("Red deer\ndensity"),
         roe_deer1999 = i18n$t("Roe deer\ndensity")
       ))
       # not sure how to rename the variables:
       #labs <- c("Temperature", "Precipitation", 
       #           "Soil pH", "Sheep and reindeer")
       # p2 <- p + facet_grid(labeller = labeller(variable = labs))
      
       print(p)
       
     })
    
    output$rcurcetext <- renderText({
      i18n$set_translation_language(input$lan)
      
      print(i18n$t("This figure shows how the habitat suitability for your selected plant (vertical axis) changes with each of the explanatory variables (horizontal axis). The units for herbivore density is metabolic biomass (kg) per squre km. This represent the feeding requirements for the species and not the number of individuals."))  # DON'T work for some reason
    })  # renderText

    
    
    
    
    # Variable importance ####
    
    output$varimp <- renderPlot({
      i18n$set_translation_language(input$lan)
      df1 <- data.frame(variables = names(IV),
                        corTest = as.numeric(NA),
                        AUCtest = as.numeric(NA))
      varimp <- list()
      
        d      <- theMod()
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
    
    output$obstext <- renderText({
      i18n$set_translation_language(input$lan)
      print(i18n$t("This is a map of the occurence data used to fit the model. For plants they are herbarium specimens."))  
    })
    
    output$rcurvetext <- renderText({
      i18n$set_translation_language(input$lan)
      print(i18n$t("This figure how the habitat suitability (on the vertical axis) changes along gradients of the explanatory variables (horizontal axes)."))  
    })  # renderText
    
    
    
    # Occurences ####
    output$occurenceMap <- renderLeaflet({
      dat <- occ[occ$species == theName(), ]
      theMap <- mapview::mapview(dat,
                                 layer.name = theName(),
                                 map.types = c("Esri.WorldShadedRelief",
                                               "Esri.WorldImagery"),
                                 cex = 5, lwd = 0,
                                 alpha.regions = 0.5,
                                 col.regions = "blue")
      theMap@map
      
    })
    
    
    # Combine responce curces, var imp and occurences ####
    output$explore <- renderUI({
      i18n$set_translation_language(input$lan)
       tabBox(width = NULL, id = 'tabset1', selected = NULL, #height = 600,
                             
                             tabPanel(i18n$t("Variable importance"),
                                      plotOutput("varimp"),
                                      textOutput("varimptext")),
              
                             tabPanel(i18n$t("Observations"),
                                      leafletOutput('occurenceMap'),
                                      textOutput("obstext")),
                                      
                             tabPanel(i18n$t("Response curves"),
                                      imageOutput("rcurves"),
                                      textOutput("rcurvetext")))
    })
    # Pictures ####
    
    output$pic <- renderImage({
      myWidth  <- session$clientData$output_pic_width
      #height <- session$clientData$output_plot2_height
      
      
      fn <- list.files("www/", pattern = theName())
      
      ifelse(length(fn)>0, outfile3 <- paste0("www/", fn), outfile3 <- "www/example.png")
       
      list(src = outfile3,
           contentType = 'image/png',
           width =   myWidth, 
           height = 'auto',
           alt = fn)
 
    }, deleteFile = FALSE) #renderImage
    
   
    
    
    
    # AUC ####
    output$AUC <- renderText({
      
      m       <- theMod()
      
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
      i18n$set_translation_language(input$lan)
      m   <- theMod()
      withProgress(message =   i18n$t("Mapping current habitat suitability ... please wait") , value = 0, {
      predict(m, IV,
                           filename = "predictions/current.img",
                           overwrite=TRUE,
                           mean=TRUE)
      })        
    })
    
    

   
    diff <- reactive({   
     i18n$set_translation_language(input$lan)
      # Outline of Norway for plotting
      outline <- readRDS("outlineNorway.RData")
      
         # -make new predictions based on user input
         IV2                  <- IV2()
         m <- theMod()
         
         if(
           identical(IV2,IV)
           ){
           pred <- current()
         } else{
           withProgress(message = 
    i18n$t("Mapping predicted changes in habitat suitability ... please wait") , value = 0, {
                        pred    <- predict(m, IV2,
                                           filename = "predictions/pred.img",
                                           overwrite=TRUE,
                                           mean=TRUE)
                          })}
                      
           diff <- pred-current()
           return(diff)
           
    }) #diff    
  

    
# here I want to make use of leafletProxy to keep the map zzom from updating when one changes the input of the control panel:
    #https://stackoverflow.com/questions/28393310/how-to-prevent-leaflet-map-from-resetting-zoom-in-shiny-app
    
#  map2 <- 
#      leaflet() %>% 
#      addProviderTiles(providers$Stamen.TonerLite) 
#    
#   observe({
#     diff <- diff()
#     
#     leafletProxy('map2')%>%
#       clearShapes() %>%
#       addRasterImage(diff, layerId= "layer1") 
#     #%>%
#     #  addLegend(pal = pal2_rev, values = c(-1:1),
#     #                            title = i18n$t("Relative change"),
#     #                            labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
#    })

   m <- reactive(
     leaflet() %>% 
     addProviderTiles(providers$Stamen.TonerLite) %>%
     addRasterImage(current(), opacity = 0.8, color = pal) %>%
     addLegend(pal = pal_rev, values = c(0:1),
               title = i18n$t("Habitat suitability"),
               labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
   )
    
   
   m2 <- reactive(
     leaflet() %>% 
             addProviderTiles(providers$Stamen.TonerLite) %>%
             addRasterImage(diff(), opacity = 0.9, color = pal2) %>%
             addLegend(pal = pal2_rev, values = c(-1:1),
                       title = i18n$t("Relative change"),
                       labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
   )
    output$mapOut <- renderUI({
      
      leafsync::sync(m(), m2() )
    })
    
  
    
#    output$mapOut <- renderUI({
#      
#      m2 <- leaflet() %>% 
#        addProviderTiles(providers$Stamen.TonerLite) %>%
#        addRasterImage(diff(), opacity = 0.9, color = pal2) %>%
#        addLegend(pal = pal2_rev, values = c(-1:1),
#                  title = i18n$t("Relative change"),
#                  labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
#      
#      leafsync::sync(m, m2) 
#    })
    
#  output$currentMap <- renderLeaflet({
#    m <- leaflet() %>% 
#      addProviderTiles(providers$Stamen.TonerLite) %>%
#      addRasterImage(current(), opacity = 0.8, color = pal) %>%
#      addLegend(pal = pal_rev, values = c(0:1),
#                title = i18n$t("Habitat suitability"),
#                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
#    #m2 <- leaflet() %>% 
#    #  addProviderTiles(providers$Stamen.TonerLite) %>%
#    #  addRasterImage(diff(), opacity = 0.9, color = pal2) %>%
#    #  addLegend(pal = pal2_rev, values = c(-1:1),
#    #            title = i18n$t("Relative change"),
#    #            labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
#    
#    m
#    #latticeView(m, m2)    # don't work
#    #leafsync::sync(m, m2) # not supported yet
#  })
#  
#  output$predMap <- renderLeaflet({
#    m2 <- leaflet() %>% 
#      addProviderTiles(providers$Stamen.TonerLite) %>%
#      addRasterImage(diff(), opacity = 0.9, color = pal2) %>%
#      addLegend(pal = pal2_rev, values = c(-1:1),
#                title = i18n$t("Relative change"),
#                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
#    m2
#  })
  
  
    
    # Plot IVs ####
    output$temp <- renderPlot({
      raster::plot(IV$temp, main = "Mean temperature of warmest quarter")
    })
    output$SoilpH <- renderPlot({
      raster::plot(IV$SoilpH, main = "Soil pH at 5cm depth")
    })
    output$moose1999 <- renderPlot({
      raster::plot(IV$moose1999, main = "Moose densities (metabolic biomass) in year 1999")
    })
    output$red_deer1999 <- renderPlot({
      raster::plot(IV$red_deer1999, main = "Red deer densities (metabolic biomass) in year 1999")
    })
    output$roe_deer1999 <- renderPlot({
      raster::plot(IV$roe_deer1999, main = "Roe deer densities (metabolic biomass) in year 1999")
    })
    output$TundraHerbivores <- renderPlot({
      raster::plot(IV$TundraHerbivores, main = "Combined sheep and reindeer densities\n(metabolic biomass) in year 1999")
    })
    output$prec <- renderPlot({
      raster::plot(IV$prec, main = "Annual precipitation (mm)")
    })
    
    output$expvar <- renderUI({
      i18n$set_translation_language(input$lan)
      
      box(title = i18n$t('Explanatory variables'), 
                      width = NULL, 
                      collapsible = T, collapsed = T,
                     
                    tabBox(width = NULL, id = 'tabset2', selected = NULL,
                         
                         tabPanel(i18n$t('Temperature'),
                                  plotOutput('temp')),
                         tabPanel(i18n$t('Precipitation'),
                                  plotOutput('prec')),
                         tabPanel(i18n$t('Soil pH'),
                                  plotOutput('SoilpH')),
                         tabPanel(i18n$t('Moose'),
                                  plotOutput('moose1999')),
                         tabPanel(i18n$t('Red deer'),
                                  plotOutput('red_deer1999')),
                         tabPanel(i18n$t('Roe deer'),
                                  plotOutput('roe_deer1999')),
                         tabPanel(i18n$t('Sheep and reindeer'),
                                  plotOutput('TundraHerbivores'))
                         ))
    })

  #The Challenge ####
    output$challenge <- renderUI({
      i18n$set_translation_language(input$lan)
        box(title = i18n$t("Take 'The Challange'"),
          width = NULL, collapsible = T, collapsed = F,
                i18n$t("Global temperatures are increasing and a warming of 2 (\u00B0C) is probable within the near future. Try increasing the temperature by 2 (\u00B0C) and see what happens to the predicted habitat suitability of you selected species. Then you challenge is: can you counteract some of these changes by modifying the herbivore densities?"))
    })
    
    
    
  
    
    observeEvent(input$reset_input, {
      shinyjs::reset("side-panel")
    })
    outputOptions(output, "cond2", suspendWhenHidden = FALSE)
    outputOptions(output, "cond1", suspendWhenHidden = FALSE)

    
    
          }) # app


