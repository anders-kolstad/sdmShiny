#### TOP #####

## Download occurence data from GBIF

# This script is for downloading occurence data from GBIF for multiple species at ones. 
# It is similar to the Rmd example that used only 2 species: https://anders-kolstad.github.io/sdmShiny/occurences
# One should only need to run this script ones and it may take a while to run.

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
plot(nOccurences_df$nOccurences)
summary(nOccurences_df$nOccurences)


## Download  ####
#For real this time:

# BIG JOB ALERT # ///////////////////////////////////////////

for(i in 1:length(mySpecies)){
  myName  <- mySpecies[i]
  myName2 <- stringr::str_split(myName, " ")[[1]]
  
  assign(
    sub(' ', '_', mySpecies[i]), 
    dismo::gbif(myName2[1], myName2[2], 
                download = T,
                geo = T, 
                sp = F) 
  )
}

# BIG JOB FINISHED # ///////////////////////////////////////////

#Two new dataframes are put in the environment. They have a lot of columns to start with, so lets get rid of som to make the objects smaller. I only need the species names and the coordinates (perhaps some more, but I can add those later). 



qc <- data.frame(Species = mySpecies,
                 lon_is_NA =                        as.numeric(NA),
                 lat_NA_when_lon_not          =     as.numeric(NA),
                 lon_is_zero =                      as.numeric(NA),
                 lat_zero_when_lon_not = as.numeric(NA))

for(i in 1:length(mySpecies)){
  
  
  d <- get(
    sub(' ', '_', mySpecies[i]))
  d <- d[,c("species","lat","lon")]
  
  # remove spaces in names (it clogs up the sdm function)
  d$species <- sub(' ', '_', d$species)
  
  # remove NA's
  w1 <- which(is.na(d$lon))
  if(length(w1) != 0) d <- d[-w1,]
  w2 <- which(is.na(d$lat))
  if(length(w2) != 0) d <- d[-w2,]
  
  # remove those with coordinates equal to zero
  w3 <- which(d$lon == 0)
  if(length(w3) != 0) d <- d[-w3,]
  w4 <- which(d$lat == 0)
  if(length(w4) != 0) d <- d[-w4,]
  
  assign(
    sub(' ', '_', mySpecies[i]),  d)
  
  qc[i,2] <- length(w1)
  qc[i,3] <- length(w2)
  qc[i,4] <- length(w3)
  qc[i,5] <- length(w4)
  
}


#A dataframe called qc tells us what has happened.
qc

#Now we can turn the dataframes into spatialPointsDataFrames, define the CRS, and plot the points. The dataset comes as lonlat.

for(i in 1:length(mySpecies)){
  
  d <- get(
    sub(' ', '_', mySpecies[i]))
  
  sp::coordinates(d) <- ~lon + lat
  sp::proj4string(d) <- sp::proj4string(raster::raster())
  
  assign(
    sub(' ', '_', mySpecies[i]),  d)
}


# CLIP
# First we need something to clip against, so we'll get an outline of Norway. 
outline <- readRDS("data/large/outline_Norway.RData")
raster::plot(outline)

# Now to clip away occurences outside this polygon

# BIG JOB ALERT # ///////////////////////////////////////////


for(i in 1:length(mySpecies)){

d <- get(
sub(' ', '_', mySpecies[i]))

d <- raster::crop(d, outline)

assign(
sub(' ', '_', mySpecies[i]),  d)
}


# Let's see it it worked.


# Now we just need to get this over to UTM32 to match the IV data, and save it on file.
myIVs      <- raster::stack('data/IV.grd')

for(i in 1:length(mySpecies)){
  
  d <- get(
    sub(' ', '_', mySpecies[i]))
  
  d <- sp::spTransform(d, myIVs[[1]]@crs)
  
  assign(
    sub(' ', '_', mySpecies[i]),  d)
}

oDat <- get(sub(' ', '_', mySpecies[1]))
for(i in 2:length(mySpecies)){
  oDat <- rbind(oDat, get(sub(' ', '_', mySpecies[i])))
}

saveRDS(oDat, 'data/large/allOccurences.RData')




oDat <- readRDS('data/large/oDat.RData')



# Check that they are inside Norway
raster::plot(myIVs$Forest_Productivity)
raster::plot(oDat,add=T)



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

