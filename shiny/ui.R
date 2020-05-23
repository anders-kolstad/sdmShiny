# Run setwd to test locally, but don't include in server version
#setwd("shiny/")

library(shiny)
library(shinyjs)
library(shinydashboard)
library(dashboardthemes)
library(leaflet)


            

dashboardPage(title = "sdmShiny",
  
# HEADER ####
  dashboardHeader(title = textOutput('myTitle'),
                  
                  titleWidth = 400,
                  
                  tags$li(div(img(src = "logoNFRw.png",
                                  height = "30px"),
                              style = "padding-top:10px; padding-bottom:10px;"),
                          class = "dropdown"),
                  tags$li(div(img(src = "NTNUlong.png",
                                  height = "30px"),
                              style = "padding-top:10px; padding-bottom:10px; padding-right:10px; padding-left:10px"),
                          class = "dropdown"),
                  dropdownMenu(type = "messages",
                  messageItem(
                    from = "AUC",
                    message = textOutput("AUC"),
                    icon("bullseye")),
                  messageItem(
                    from = "Method",
                    message = "GAM",
                    icon("bullseye")),
                  messageItem(
                    from = "# unique samples",
                    message = textOutput("n"),
                    icon("bullseye"))
                  )
  ),
              
# SIDEBAR ####    
  dashboardSidebar( 
    
    radioButtons("lan", "",
                 c("English" = "en",
                   "Norsk" = "no")),
    
    uiOutput('spLanguage'),
#    radioButtons("spLang",
#                 "List species by:",
#                 c("Scientific names" = "sci", "Norwegian names" = "com")
#                      ),
    uiOutput('spNames')
    ), # end sidebar
  
  
  # BODY ####
  dashboardBody(
    shinyDashboardThemes(
      theme = "onenote"
    ),
    
    # Aligning plots to  center of boxes
    tags$head(tags$style(HTML(".grViz { width:100%!important;}"))),
    
    # top row (orange box)
    fluidRow(width=12, uiOutput("top")),
    
    
    fluidRow(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"),
      
      # First column ####
      column(width = 8,
        box(width = NULL, 
          plotOutput("map")),
        
        uiOutput('challenge'),
        
        
        tabBox(width = NULL, id = 'tabset1', selected = "Records", #height = 600,
               tabPanel("Records",
                        leafletOutput('occurenceMap')),
               tabPanel("Variable importance",
                        plotOutput("varimp"),
                        textOutput("varimptext")),
               tabPanel("Response curves",
                        imageOutput("rcurves"))),
        box(title = 'Explanatory variables', 
            width = NULL, 
            collapsible = T, collapsed = T,
            #status = "primary", 
            #solidHeader = TRUE,
          tabBox(width = NULL, id = 'tabset2', selected = "Temperature",
               
               tabPanel('Temperature',
                        plotOutput('temp')),
               tabPanel('Precipitation',
                        plotOutput('prec')),
               tabPanel('soil pH',
                        plotOutput('SoilpH')),
               tabPanel('Moose',
                        plotOutput('moose1999')),
               tabPanel('Red deer',
                        plotOutput('red_deer1999')),
               tabPanel('Roe deer',
                        plotOutput('roe_deer1999')),
               tabPanel('Sheep and reindeer',
                        plotOutput('TundraHerbivores'))
               )),
               
        uiOutput("usUI")
        #,
        #uiOutput("credUI")
                ),  # end of column
      column( width = 4,
        box(width = NULL,
         uiOutput("contr")),
        #box(width = NULL, align = "center", height = 800,
            imageOutput("pic"),
            uiOutput("credUI")
             ))   # end of row
      
    

  ) # Body
) # Page

  