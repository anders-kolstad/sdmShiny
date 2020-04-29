

library(shiny)
library(raster)
library(dismo)
library(sdm)



# Model objects
# 5*3 model = _m
# Replicated single method = _ms
# Best candidate model = _bcm

# Species
# Primula_scandinavica
# Carex_simpliciuscula

# IVs
# Temperature
# Precipitation
# Soil pH
# Sheep and reindeer

# import species list
s <- read.csv('data/species_list.csv')
s <- as.vector(s$x)

# Bring all models into the environment
files <- list.files("models/sdmModels/", pattern = ".sdm")



Priscamaxent <- readRDS("data/prisca_maxent.rds")


shinyServer(
  function(input, output) {
    #modify env vars by selectable parameters
    output$map <- renderPlot({
      selectvars_use<-selectvarsshiny
      selectvars_use[[1]]<-selectvarsshiny[[1]]+(selectvarsshiny[[1]]*(input$herbivory/100))
      selectvars_use[[2]]<-selectvarsshiny[[2]]+input$temperature
      names(selectvars_use)<-names(selectvarsshiny)
      
      if(input$species=='Primula scandinavica'){
        p1<-predict(Priscamaxent,selectvarsshiny)
        p2<-predict(Priscamaxent,selectvars_use)
        diff<-p2-p1
        diffp <- diverge0(levelplot(diff,margin=F),ramp='RdBu')
        predp<-levelplot(stack(p1,p2),margin=F,main="Habitat suitability \n Primula scandinavica",names.attr=c("Current","Changed"))
        return(print(predp  )#, split=c(1,1,1,2),more=T)
               # print(diffp,split=c(1,2,1,2))
        )}
      
      else if(input$species=='Botrychium lanceolatum'){
        p1<-predict(Botlanmaxent,selectvarsshiny)
        p2<-predict(Botlanmaxent,selectvars_use) 
        return(  levelplot(stack(p1,p2),margin=F,main="Habitat suitability \n Botrychium lanceolatum",names.attr=c("Current","Changed"))
        )}
    })
  }
)


# predvars<-stack("S:\\SDMs\\InputData\\predvars1k")
# as.data.frame(names(predvars))
# selectvarsShiny<-subset(predvars,c(14,25))
# selectvarsShiny<-aggregate(selectvarsShiny,10)
# writeRaster(selectvarsShiny,file='data/selectvars',overwrite=T) 
# alpspdat1<-read.table("S:\\SDMs\\SpeciesData\\RedListSpQC.csv",header=T,sep=",")
# alpspdat<-alpspdat1[!is.na(alpspdat1$Xutm) & !is.na(alpspdat1$Yutm),]
# Prisca<-alpspdat[alpspdat$species=='Primula scandinavica',]
# mePrisca<-maxent(selectvarsShiny,Prisca[,c(2:3)])
# saveRDS(Prisca,file="data/prisca_obs.rds")
# saveRDS(mePrisca,file='data/prisca_maxent.rds')
#Botlan<-alpspdat[alpspdat$species=='Botrychium lanceolatum',]
#meBotlan<-maxent(selectvarsShiny,Botlan[,c(2:3)])
#saveRDS(Botlan,file="data/botlan_obs.rds")
#saveRDS(meBotlan,file='data/botlan_maxent.rds')

#prisca<-readRDS('ShinyApps/testSDM/data/prisca_obs.rds')
#selecvatrs1 <- readRDS("ShinyApps/testSDM/data/selectvars.rds")
#Priscamaxent <- readRDS("ShinyApps/testSDM/data/prisca_maxent.rds")
prisca<-readRDS('data/prisca_obs.rds')
Priscamaxent <- readRDS("data/prisca_maxent.rds")
botlan<-readRDS('data/botlan_obs.rds')
Botlanmaxent <- readRDS("data/botlan_maxent.rds")

selectvarsshiny<-stack('data/selectvars')

#diverge0 <- function(p, ramp) {
#  # p: a trellis object resulting from rasterVis::levelplot
#  # ramp: the name of an RColorBrewer palette (as character), a character 
#  #       vector of colour names to interpolate, or a colorRampPalette.
#  if(length(ramp)==1 && is.character(ramp) && ramp %in% 
#     row.names(brewer.pal.info)) {
#    ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
#  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
#    ramp <- colorRampPalette(ramp)
#  } else if(!is.function(ramp)) 
#    stop('ramp should be either the name of a RColorBrewer palette, ', 
#         'a vector of colours to be interpolated, or a colorRampPalette.')
#  rng <- range(p$legend[[1]]$args$key$at)
#  s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
#  i <- findInterval(rng[which.min(abs(rng))], s)
#  zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
#  p$legend[[1]]$args$key$at <- s[zlim]
#  p$par.settings$regions$col <- ramp(1000)[zlim[-length(zlim)]]
#  p
#}

#source("helpers.R")
#library(rasterVis)
#require(RColorBrewer)
#library(maps)
#library(mapproj)

