#### TOP #####

## Download occurence data from GBIF

# This script is for downloading occurence data from GBIF for multiple species at ones. 
# It is similar to the Rmd example that used only 2 species: https://anders-kolstad.github.io/sdmShiny/occurences
# One should only need to run this script once and it may take a while to run.

# Anders L. Kolstad May 2020

##  Get species list ####
#This function  produces a list of species that we will later use to harvest occurence data from gbif. 


source("./R/spList.R")
mySpecies <- sl()
mySpeciesdf <- sl(df = T)

head(mySpecies)
head(mySpeciesdf)


## Occurence data
#To get occurence data I will use the gbif function in the dismo package. 
#It can only handle one species at the time, so I will need to make a for-loop.

## Extent ####
# We dont want to download all records, just those for Norway.
source("R/norway.R")
nor <- norway(lonlat = TRUE)
ext <- raster::extent(nor)


### Test run ####
#Let's do a test loop without downloading anything, just seeing how many records there are.
#This are all the records, not only from Norway.
nOccurences_df <- data.frame(species = mySpecies, 
                             nOccurences = as.numeric(NA))

for(i in 1:length(mySpecies)){
myName  <- mySpecies[i]
myName2 <- stringr::str_split(myName, " ")[[1]]
nOccurences_df$nOccurences[i] <- dismo::gbif(myName2[1], myName2[2], download = F) 
}

View(nOccurences_df)
#plot(nOccurences_df$nOccurences)
summary(nOccurences_df$nOccurences)

#Let's do the same, but with an extent argument
nOccurences_df2 <- data.frame(species = mySpecies, 
                             nOccurences = as.numeric(NA))

for(i in 1:length(mySpecies[1:80])){
  myName  <- mySpecies[i]
  myName2 <- stringr::str_split(myName, " ")[[1]]
  nOccurences_df2$nOccurences[i] <- dismo::gbif(myName2[1], myName2[2], download = F, ext = ext) 
}

View(nOccurences_df2)
plot(nOccurences_df2$nOccurences)
summary(nOccurences_df2$nOccurences)

sum(nOccurences_df$nOccurences)
sum(nOccurences_df2$nOccurences)
## Download  ####
#For real this time:

# BIG JOB ALERT # ///////////////////////////////////////////

for(i in 1:length(mySpecies[1:80])){
  myName  <- mySpecies[i]
  myName2 <- stringr::str_split(myName, " ")[[1]]
  
  assign(
    sub(' ', '_', mySpecies[i]), 
    dismo::gbif(myName2[1], myName2[2], 
                download = T,
                geo = T, 
                sp = F,
                ext = ext) 
  )
}
# 8:42 - 
# BIG JOB FINISHED # ///////////////////////////////////////////

# basisOfRecord == "PRESERVED_SPECIMEN"
# !is.na(year)
# keep year
# country == "Norway" or just clip top be sure

#Two new dataframes are put in the environment. They have a lot of columns to start with, so lets get rid of som to make the objects smaller. I only need the species names and the coordinates (perhaps some more, but I can add those later). 



qc <- data.frame(Species = mySpecies[1:80],
                 lon_NA           =     as.numeric(NA),
                 lat_NA           =     as.numeric(NA),
                 lon_zero         =     as.numeric(NA),
                 lat_zero         =     as.numeric(NA),
                 year_NA          =     as.numeric(NA),
                 unvalidated      =     as.numeric(NA),
                 original_length  =     as.numeric(NA),
                 new_length       =     as.numeric(NA),
                 deleted          =     as.numeric(NA))

for(i in 1:80){
  
 
  d <- get(    sub(' ', '_', mySpecies[i]))
  
  if(!is.null(d)){
  d <- d[,c("species","lat","lon", "year", "basisOfRecord", "occurrenceID")]
  
  # remove spaces in names (it clogs up the sdm function)
  d$species <- sub(' ', '_', d$species)
  
  # number of records:
  n <- nrow(d)
  
  # remove NA's
  w1 <- d$occurrenceID[which(is.na(d$lon))]
  w2 <- d$occurrenceID[which(is.na(d$lat))]
  
  # remove those with coordinates equal to zero
  w3 <- d$occurrenceID[which(d$lon == 0)]
  w4 <- d$occurrenceID[which(d$lat == 0)]
  
  # remove those with no year
  w5 <- d$occurrenceID[which(is.na(d$year))]
  
  # remove 'HUMAN OBSERVATIONS'
  w6 <- d$occurrenceID[which(d$basisOfRecord == "HUMAN_OBSERVATION")]
  
  w <- c(w1, w2, w3, w4, w5, w6)

  if(length(w) != 0) {d <- d[!d$occurrenceID %in% w,]}
  
  
  # remaining records
  n2 <- nrow(d)
  
  # deleted
  n3 <- n-n2
  
  
  
  assign(
    sub(' ', '_', mySpecies[i]),  d)
  
  qc[i,2] <- length(w1)
  qc[i,3] <- length(w2)
  qc[i,4] <- length(w3)
  qc[i,5] <- length(w4)
  qc[i,6] <- length(w5)
  qc[i,7] <- length(w6)
  qc[i,8] <- n
  qc[i,9] <- n2
  qc[i,10] <- n3
  
  } else{
    name <- as.name(sub(' ', '_', mySpecies[i]))
    rm(name)}
}


#A dataframe called qc tells us what has happened.
View(qc)

# OBS. fitModels should not use mySpecies as before, since species get dropped is there are no records.

#Now we can turn the dataframes into spatialPointsDataFrames, define the CRS, and plot the points. The dataset comes as lonlat.


for(i in 1:80){
  
  name <- sub(' ', '_', mySpecies[i])
  if(exists(name)){
    
  d <- get( name   )
  
  if(is.data.frame(d)){
   if(nrow(d)>30){
  sp::coordinates(d) <- ~lon + lat
  sp::proj4string(d) <- sp::proj4string(raster::raster())
  
  assign(
    sub(' ', '_', mySpecies[i]),  d)
  } else{
    name <- sub(' ', '_', mySpecies[i])
    rm(list = name)
    } # rm if <30 records
  } # is.data.frame
  } # if exists
}



# If the above loop didn't work..:
#for(i in 1:80){
#  
#  name <- sub(' ', '_', mySpecies[i])
#  
#  if(exists(name)){
#    d <- get(  name)}
#  
#  if(class(d) == "SpatialPointsDataFrame"){
#  
#  
#  if(length(d) < 30){
#        rm(list = name)}
#  } # if spatial points
#  }


# CLIP
# First we need something to clip against, so we'll get an outline of Norway. 
outline <- readRDS("data/large/outline_Norway.RData")
raster::plot(outline)

# Now to clip away occurences outside this polygon

# BIG JOB ALERT # ///////////////////////////////////////////


for(i in 1:80){

  name <- sub(' ', '_', mySpecies[i])
  
  if(exists(name)){
    
    d <- get( name   )
    d <- raster::crop(d, outline)

   assign(name,  d)
  } # if exists
 }

# >30min

# Let's see it it worked.
#raster::plot(nor)
#raster::plot(Ajuga_reptans, add=T)
#raster::plot(Ulmus_glabra, add=T)

#mapview::mapview(Ajuga_reptans, 
#                 map.types = c("Esri.WorldShadedRelief",
#                               "Esri.WorldImagery"),
#                 cex = 5, lwd = 0,
#                 alpha.regions = 0.5,
#                 col.regions = "blue")


# Now we just need to get this over to UTM32 to match the IV data, and save it on file.
myIVs      <- raster::stack('data/IV.grd')

for(i in 1:80){
  name <- sub(' ', '_', mySpecies[i])
  if(exists(name)){
  d <- get(name)
  
  d <- sp::spTransform(d, myIVs[[1]]@crs)
  
  assign(name,  d)
}}


oDat <- get(sub(' ', '_', mySpecies[1]))
for(i in 2:80){
    name <- sub(' ', '_', mySpecies[i])
  if(exists(name)){
    d    <- get(name)
    oDat <- rbind(oDat, d)
}}


#nor2 <- norway(lonlat = FALSE)
#par(mfrow=c(1,2))
#raster::plot(nor2)
#raster::plot(oDat, add=T)
#raster::plot(nor2)
#raster::plot(oDat[oDat$species == "Ulmus_glabra",], add=T)

saveRDS(oDat, 'data/allOccurences1-80.RData')

#testImp <- readRDS('data/large/allOccurences1-80.RData')



## Check sample sizes
#Let's see how mny point there are for each species, and how many of these that fall on the same 1x1 km grid cells.

t <- myIVs$elev
df <- data.frame(species = NA,
                points = as.numeric(NA),
                unique = as.numeric(NA))
for(i in 1:length(unique(oDat$species))){
  s       <- unique(oDat$species)[i]
  df[i,1] <- paste(s)
  t1      <- oDat[oDat$species == s,]
  df[i,2] <- length(t1)
  u       <- raster::rasterize(t1, t, 1, fun = "count")
  df[i,3] <- length(u[u>0])
}
df

