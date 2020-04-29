# SDM-Shiny

library(shiny)

shinyUI(fluidPage(
  titlePanel("Interactive distribution modelling"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Explore effects of changing climate and herbivore densities on distribution of rare plant species."),
      
      selectInput("species", 
                  label = "Choose a species to display",
                  choices = c("Primula scandinavica","Carex simpliciuscula"),
                  selected = "Primula scandinavica"),
      
      sliderInput("temperature", 
                  label = "Change in summer temperature degC:",
                  min = -5, max = 5, value = c(0)),
      
      sliderInput("herbivory", 
                  label = "Change in sheep and reindeer density %:",
                  min = -50, max = 50, value = c(0))
      
      sliderInput("precipitation", 
                  label = "Change in annual precipitation mm/year:",
                  min = -1000, max = 1000, value = c(0))
    ) ,
    
    mainPanel(plotOutput("map")
              #,imageOutput('image')
    )
    
  )))