# Run setwd to test locally, but don't include in server version
#setwd("shiny/")

library(shiny)
library(shinyjs)
library(shinydashboard)


# translations

# Get species list ####
myS <- list.files("sdmModels/", pattern = ".sdm")
myS2 <- as.list(NA)
for(i in 1:length(myS)){
  myS2[i] <- 
    paste(
      stringr::str_split(myS[i], "_")[[1]][1],
      stringr::str_split(myS[i], "_")[[1]][2],
      collapse = " ")
}
myS3 <- unique(as.character(myS2))


dashboardPage(
  
# HEADER ####
  dashboardHeader(title = i18n$t("Interactive distribution modelling"),
                  
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
    radioButtons("lan", "Choose language:",
                              c("English" = "en",
                                "Norsk" = "no")),
    radioButtons("species",
                 "Pick a species",
                 choiceNames = myS3,
                 choiceValues = myS3
                 )),
  
  
  # BODY ####
  dashboardBody(
    tags$head(tags$style(HTML(".grViz { width:100%!important;}"))),
    fluidRow(width=12, uiOutput("top")),
    
    fluidRow(
      box(width = 8, height = 650,
          plotOutput("map")),
          
      box(width = 4, height = 650,
      uiOutput("contr"))
    ),
      
    
    
    fluidRow(
      box(width = 4, align = "center", height = 600,
        imageOutput("pic")
              ),
      tabBox(width = 8, id = 'tabset1', selected = "Variable importance", height = 600,
             tabPanel("Variable importance",
                      plotOutput("varimp"),
                      textOutput("varimptext")),
             tabPanel("Response curves",
                      imageOutput("rcurves"))
            
      )
    ),
    
# Photo credits ####
  fluidRow(
             uiOutput("credUI"),
             uiOutput("usUI"),
             uiOutput("dissclaimer"),
             uiOutput("refs"))


  ) # Body
) # Page

  