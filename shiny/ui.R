# Run setwd to test locally, but don't include in server version
#setwd("shiny/")

library(shiny)
library(shinyjs)
library(shinydashboard)
library(dashboardthemes)
library(leaflet)


            

dashboardPage(
  
# HEADER ####
  dashboardHeader(title = textOutput('myTitle'),
                  
                  titleWidth = 450,
                  
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
    
    tags$head(tags$style(HTML(".grViz { width:100%!important;}"))),
    fluidRow(width=12, uiOutput("top")),
    
    fluidRow(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"),
      column(width = 8,
        box(width = NULL, 
          plotOutput("map")),
        tabBox(width = NULL, id = 'tabset1', selected = "Records", #height = 600,
               tabPanel("Records",
                        leafletOutput('occurenceMap')),
               tabPanel("Variable importance",
                        plotOutput("varimp"),
                        textOutput("varimptext")),
               tabPanel("Response curves",
                        imageOutput("rcurves"))),
               
        uiOutput("usUI"),
        uiOutput("credUI")
                ),  # end of column
      column( width = 4,
        box(width = NULL,
         uiOutput("contr")),
        #box(width = NULL, align = "center", height = 800,
            imageOutput("pic")
             ))   # end of row
      
    

  ) # Body
) # Page

  