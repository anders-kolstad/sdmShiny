library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinydashboard)
library(dashboardthemes)
library(leaflet)
library(readr)
#
library(sdm)
library(sp)
library(raster)
library(dismo)
library(rgdal)
library(lattice)
library(latticeExtra)
library(plyr)


# Species list ------------------------------------------------------------

namelist <- readRDS('namelist.RData') # se R/vernacular.R
# Picking species that have some relevance to Trøndela
sel  <- c("Ulmus glabra",
          #"Androsace septentrionalis",
          #"Arnica montana",
          "Gentianella campestris",
          "Primula scandinavica",
          #"Pulsatilla vernalis",
          "Campanula cervicaria",
          #"Comastoma tenellum",
          "Cypripedium calceolus",
          "Pseudorchis albida",
          "Schoenus ferrugineus"
)

namelist <- namelist[namelist$sci %in% sel,]
namelist$sci        <- as.character(namelist$sci)       # scientific names
namelist$vern_nor   <- as.character(namelist$vern_nor)  # norwegian names
namelist$vern_eng   <- as.character(namelist$vern_eng)  # english names with scientific names as placeholders
namelist$vern_eng2  <- as.character(namelist$vern_eng2) # english names


# Data -----------------------------------------------------------------------


# IVs ####
IV <- raster::stack("IVapp.grd")


# Occurences
occ <- readRDS("allOccurences.RData")
occ <- sp::spTransform(occ, "+proj=longlat +datum=WGS84 +no_defs")



# Colours -----------------------------------------------------------------

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



# Species details ---------------------------------------------------------
Ulmus_glabra_i <- "Alm er en av de såkalte edelløvtrærne vår. Det vil si at den liker det litt varmt. Den er vurdert som 'sårbar' i norsk natur på grunn av sterk tilbakegang i senere år. Dette skyldes både en sykdom (almesyke) og beiteskader av blant annet hjortedyr. Utenfor Trondheim kan man se en flott almeskog i Lauglolia, nord for Gaulosen, men alm finnes også lengre nord til Nordland."
Gentianella_campestris_i <- "Bakkesøte er en engart og den krever også noe kalk i bakken. Den er vurdert som nær truet i norsk natur på grunn av nedgang i den tradisjonelle slåtte- og beitebruken slik at habitatet dens blir mer sjeldent. Den vokser for eksempel på Lianjordene i Trondheim."
Primula_scandinavica_i <- "Fjellnøkleblom er endemisk til Skandinavia. Det vil si den finnes bare her. Og her vokser den på litt utsatte steder, som i fjellenger og ved kysten, helst på litt kalkrike plasser. Fjellnøkleblom er ikke på rødlisten (!), men den var det forrige gang det ble laget en rødliste, i 2010. Grunnen til at den nå er vurdert som livskraftig er at reduksjonen som først skjedde som en følge av opphørt beite og seterdrift nå ser ut til å ha stabilisert seg."
Campanula_cervicaria_i <- "Stavklokke går omtrendt så langt nord som til Trondheimsfjorden. Den trives i varme skogslier som gjerne beites eller holdes åpne på andre måter. Arter er nær truet på grunn av at den er i tilbakegang."
Cypripedium_calceolus_i <- "Marisko er en fantastisk fin orkidee som vokser i kalkrike skoger og i rasmark. Arter er i pågående tilbakegang og er vurdert som nær truet i norsk natur. Marisko har en noe østlig utbredelse, men går fra sør til nord i landet."
Pseudorchis_albida_i <- "Hvitkurle er en orkide som er knyttet til slåttemarker, beitemarker og beitet skog. Den trues blant annet av gjengroing og er vurdert som nær truet i norsk natur. Hvitkurle er fortsatt tallrik på Lianjordene  i Trondheim."
Schoenus_ferrugineus_i <- "Brunskjene vokser på ekstremrike myrer og er vurdert som truet i norsk natur fordi disse myrene blir mer og mer sjeldne, blant annet på grunn av grøfting og planting av skog. I Trondheim kan arten ses blant annet på Rørmyra ved Granåsen."

infoknapp <- "Denne appen er utviklet av Anders L. Kolstad og James Speed ved NTNU Vitenskapsmuseet i år 2020. planteR lar deg endre på klimavariabler og beitetrykket, 
og en artsutbredelsesmodell vil kjøre i sanntid for å beregne hvordan dette vil slå ut for den arten du har valgt. Grønne områder i kartet betyr at endringene i klima og/eller beitetrykket som du har gjort virker positivt og øker sansyneligheten for at arten du har valgt vil kunne vokse i disse områdene. Røde felter betyr det motsatte. 
Merk at modellen er en forenkling og er assosiert med stor usikkerhet, men den gir oss en pekepin om hva som kan komme til å skje i fremtiden."
# UI ----------------------------------------------------------------------


ui <- dashboardPage(title = "planteR",
              
              # HEADER ####
              dashboardHeader(disable = T),
              
              # SIDEBAR ####    
              dashboardSidebar(disable = T
              ), # end sidebar
              
              
              # BODY ####
              dashboardBody(
                # Hide error messages
                tags$style(type="text/css",
                           ".shiny-output-error { visibility: hidden; }",
                           ".shiny-output-error:before { visibility: hidden; }"),
                
                
                
                
                
                tags$style(type = "text/css", ".map-container {position:absolute; top:0; bottom:0; right:0; left:0;}"),
                tags$div(
                  class = "map-container",
                  leafletOutput("mapOut", width = "100%", height = "100%")
                ),
                
                tags$style(type = "text/css", "#dropdown {margin-top: 180px; margin-left: 80px;}"),
                dropdownButton(
                  inputId = "dropdown",
                  circle = T, status = "success",
                  icon = icon("leaf"), width = "200px",
                  tooltip = tooltipOptions(title = "Velg en plante"),
                  radioGroupButtons(
                    inputId = "species",
                    choices = c(namelist$vern_nor),selected = namelist$vern_nor[2],
                    justified = T,
                    direction = 'vertical'
                  )
                ),
                
                tags$style(type = "text/css", "#dropdown2 {margin-top: 20px; margin-left: 80px;}"),
                dropdownButton(
                  inputId = "dropdown2",
                  circle = T, status = "warning",
                  icon = icon("info"), width = "300px",
                  tooltip = tooltipOptions(title = "Se info"),
                 uiOutput("details")
                  
                ),
                
                tags$style(type = "text/css", "#dropdown3 {margin-top: 20px; margin-left: 80px;}"),
                dropdownButton(
                  inputId = "dropdown3",
                  circle = T, status = "info",
                  icon = icon("question"), width = "300px",
                  tooltip = tooltipOptions(title = "Se info"),
                  h4(infoknapp)
                  
                ),
                
                absolutePanel(
                  id = "controls", 
                  class = "panel panel-default", 
                  fixed = F,
                  draggable = TRUE, 
                  top = 50, left = 'auto', right = 50, bottom = "auto",
                  width = 300, height = "auto",
                  uiOutput("contr")
                ),
                
                absolutePanel(
                  id = "pic", 
                  class = "panel panel-default", 
                  fixed = F,
                  draggable = TRUE, 
                  top = 'auto', left = 50, right = 'auto', bottom = 50,
                  width = 400, height = 400,
                  imageOutput("pic")  ),
                
                absolutePanel(id = "name",class = "panel panel-title",top  = 80, left  = 80, uiOutput('selectedSp'),draggable = F)
              
                
                
                
              )# Body
) # Page


# END UI ---------------------------------------------------------------------


# Server ------------------------------------------------------------------


server <- function(input, output, session) {
    
    # translated ui ####
    
    
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
      IV2$moose1999        <- IV2$moose1999        * (input$herbivory/100 + 1)
      IV2$red_deer1999     <- IV2$red_deer1999     * (input$herbivory/100 + 1)
      IV2$roe_deer1999     <- IV2$roe_deer1999     * (input$herbivory/100 + 1)
      IV2$temp             <- IV2$temp             + input$temperature         
      IV2$prec             <- IV2$prec             * (input$precipitation/100+1)
      return(IV2)
    })
 
    
    
    
    
    output$contr <- renderUI({
      
      box(
        shinyjs::useShinyjs(), # to enable reset function
        id = "side-panel",
        status = "primary", 
        solidHeader = TRUE,
        width = NULL,
        
        sliderInput("temperature", 
                    label = "Temperatur (\u00B0C):",
                    min = -5, max = 5, value = c(0)),
        sliderInput("precipitation", 
                    label = "Nedbør (%):",
                    min = -50, max = 50, value = c(0)),
        sliderInput("herbivory", 
                    label = "Ville beitedyr (%):",
                    min = -100, max = 100, value = c(0)),
        
        
        actionButton("reset_input", "Nullstill"),
        br(),br(),
        shiny::HTML("<p>I kartet blir <span style='color: red'>røde områder</span> dårlige, og 
                    <span style='color: green'> grønne områder </span> bedre egnet som voksesteder. 
                    <span style='color: blue'> Blå prikker </span> er hvor planten vokser i dag</p>")
      ) 
    } )
    
    
    
    
    output$selectedSp <- renderUI({
      
      temp <- which(namelist == input$species, arr.ind = T)
      row <- temp[1,1]
      HTML(
       paste0(
        h1(namelist$vern_nor[row]), 
        " (",
        em(namelist$sci[row]), 
        ")")
       )
    })
    
    
    
    
    # PLOTS ####
    
    current <- reactive({
      raster::stack(paste0("predictions/", theName(), ".img"))
    })
    
    
    #
    
    diff <- reactive({   
      
      # -make new predictions based on user input
      IV2                  <- IV2()
      m <- theMod()
      
      if(
        identical(IV2,IV)
      ){
        pred <- current()
      } else{
        withProgress(message = 
                       "Jobber med saken ... " , value = 0, {
                         pred    <- predict(m, IV2,
                                            filename = "predictions/pred.img",
                                            overwrite=TRUE,
                                            mean=TRUE)
                       })}
      
      diff <- pred-current()
      return(diff)
      
    }) #diff    
    
    
    
    dat <- reactive({
      
      occ[occ$species == theName(), ]
      
    })
    
    output$mapOut <- renderLeaflet({
      
      leaflet() %>% 
        addProviderTiles(providers$Stamen.TonerLite) %>%    # boring basemap, but doesn't mask the raster layer as much
        addRasterImage(diff(), opacity = 1, color = pal2) %>%
        addCircleMarkers(data = dat(), radius =1, fill=F)
      
    })
    
    
    output$details <- renderUI({
      det <- get(paste0(theName(), "_i"))
      h4(det)
    })
    
    
    output$pic <- renderImage({
      myWidth  <- session$clientData$output_pic_width
      outfile3 <- paste0("www/", theName(), ".JPG")
      list(src = outfile3,
           contentType = 'image/jpeg',
           width =   myWidth, 
           height = 'auto')
    }, deleteFile = FALSE) #renderImage
    
    observeEvent(input$reset_input, {
      shinyjs::reset("side-panel")
    })
    
 
    
    
  } 

shinyApp(ui = ui, server = server)
