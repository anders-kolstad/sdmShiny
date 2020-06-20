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
                    message = "mean of 3 GLM models",
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

    uiOutput('spNames')
    ), # end sidebar
  
  
  # BODY ####
  dashboardBody(
    shinyDashboardThemes(
      theme = "onenote"
    ),
    
    # Aligning plots to  center of boxes
    tags$head(tags$style(HTML(".grViz { width:100%!important;}"))),
    
    # Center text
    tags$head(tags$style(HTML("
                                #selectedSp {
                                  text-align: center;
                                }"))),
    
    # top row (orange box)
    fluidRow( uiOutput("top")),
    
                 
    
    fluidRow(
      # Hide error messages
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"),
      
      # First column ####
      column(width = 8,
             
            textOutput("selectedSp"),
            
            box(width = NULL,
             uiOutput('mapOut')),
            
            box(width=NULL,
                 uiOutput('challenge')),
        
                  
        
            uiOutput('explore'),


          uiOutput('expvar'),

               
        uiOutput("usUI")
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

  