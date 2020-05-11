# SDM-Shiny

library(shiny)
  

# Get species list
myS <- list.files("models/sdmModels/", pattern = ".sdm")
myS2 <- as.list(NA)
for(i in 1:length(myS)){
  myS2[i] <- 
    paste(
      stringr::str_split(myS[i], "_")[[1]][1],
      stringr::str_split(myS[i], "_")[[1]][2],
      collapse = " ")
}
myS3 <- unique(as.character(myS2))


shinyUI(fluidPage(
  titlePanel("Interactive distribution modelling"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Explore effects of changing climate and herbivore densities on distribution of rare plant species."),
      
      
      # Possible improvement: https://shiny.rstudio.com/articles/selectize.html
      selectInput("species", 
                  label = "Choose a species to display",
                  choices = myS3, 
                  selected = myS3[1]),
      
      #selectInput("modelling approach", 
      #            label = "Choose modelling approach",
      #            choices = c("Ensemble","Replicated maxent (n=5)", "Best candidate model"),
      #            selected = "Ensemble"),
      
      sliderInput("temperature", 
                  label = "Change in mean summer temperature (\u00B0C)",
                  min = -5, max = 5, value = c(0)),
      
      sliderInput("herbivory", 
                  label = "Change in sheep and reindeer density (%):",
                  min = -50, max = 50, value = c(0)),
      
      sliderInput("precipitation", 
                  label = "Change in annual precipitation (%):",
                  min = -50, max = 50, value = c(0))
    ) ,
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plots", 
                 plotOutput("map"),
                 imageOutput("pic"),
                 textOutput("credits")), 
        tabPanel("Variable importance", 
                 textOutput("varimptext"),
                   imageOutput("varimp")), 
        tabPanel("Response curves", imageOutput("rcurves")),
        tabPanel("Help", 
                 textOutput("help"),
                 tags$style(type="text/css", "#help {white-space: pre-wrap;}")
        )
      ))
    
  )))