
-   [sdmShiny](#sdmshiny)
    -   [Installation](#installation)
    -   [Example application](#example-application)
-   [Documentation](#documentation)
    -   [Environmental data](#environmental-data)
        -   [Overview](#overview)
        -   [Existing raster stacks](#existing-raster-stacks)
            -   [The 'alpine paper'](#the-alpine-paper)
            -   [The 'forest' paper](#the-forest-paper)
        -   [Worldclim data and DTM](#worldclim-data-and-dtm)
    -   [Change projection and stack](#change-projection-and-stack)
        -   [Modify IVs](#modify-ivs)
            -   [Clip](#clip)
            -   [Ratify categorical layers](#ratify-categorical-layers)
    -   [Write IV file](#write-iv-file)
        -   [Ratify v2](#ratify-v2)
            -   [Forest productivity](#forest-productivity)
            -   [Land cover](#land-cover)
            -   [Forest type](#forest-type)
    -   [Get species list](#get-species-list)
    -   [Occurence data](#occurence-data)
    -   [SDM](#sdm)
        -   [SDM-data](#sdm-data)
        -   [5 x 3 models](#x-3-models)
            -   [Ensemble](#ensemble)
            -   [Presence-absence map](#presence-absence-map)
        -   [Best candidate model](#best-candidate-model)
        -   [Replicated single method](#replicated-single-method)
        -   [Response curves](#response-curves)
        -   [Variable importance](#variable-importance)

<!-- README.md is generated from README.Rmd. Please edit that file -->
sdmShiny
========

<!-- badges: start -->
<!-- badges: end -->
Last update:

``` r
Sys.time()
#> [1] "2020-04-29 14:50:14 CEST"
```

This project is for disseminating the species distribution modeling work done in James Speed's group at the NTNU University Museum. We will use web-based Shiny apps to present distribution maps of several species and allow these to change with the predictions of the SDM as the user tweaks the parameters for climate and herbivory. The Shiny app will look something like this: ![The Shiny app will look something like this](figures/app.png)

Installation
------------

Developers should clone the repo and work from there.

Example application
-------------------

Follow this link <https://anderskolstad.shinyapps.io/demoSDM/>

Documentation
=============

This section explains the workflow that ended up the the shiny app. Large files around 100mb or greater, or unessential raster files etc., are in the data/large/ folder whihc is not pushed (it's in the .gitignore file) and therefore only exists locally with Anders. The same for RData files intill I find a way to load them without getting magic number errors (readRDS don't work with knitr).

Environmental data
------------------

### Overview

The alpine red-listed species paper used the following IVs, in the approximate order of importance:

-   bio10 (temperature)
-   bio12 (precipitation)
-   soil ph
-   Tundraherbivores - the combined metabolic biomass of sheep aand reindeer tundraherbivores
-   bio15 (precipitation seasonality)
-   AR50 categorical land-use classes (how many levels I don't know)

The 'forest paper' used these:

-   bio10
-   bio12
-   Forest type (three levels, from AR50)
-   soil pH
-   Cervid densities...

### Existing raster stacks

#### The 'alpine paper'

From the alpine red-listed species paper, I got this sent over from James:

``` r
IV <- raster::stack("data/large/selectvars.grd")
names(IV)
#> [1] "TundraHerbivores"  "MeanTempWarmQuart"
```

The first layer is the combined metabolic biomas of reindeer (wild and semi-domesticated) and sheep from the year 1999. I will use this later, but not the worldclim variable (see below). The resolution is 10km, which is larger than for the following layers.

``` r
#reindeerSheep <- IV[[1]]
#writeRaster(reindeerSheep, 'data/large/reindeerSheep')
reindeerSheep <- raster::stack('data/large/reindeerSheep.grd')
rm(IV)
raster::plot(reindeerSheep)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

#### The 'forest' paper

From another of James' projects there is a file already collated with environmental data. It was downloaded from NTNU box: <https://ntnu.app.box.com/s/wcmr0dgoyz2yu6ielw6er1pm7h0gaisa/file/393633279036>

``` r
PredVars <- raster::stack("data/large/PredictorVariables.grd")
names(PredVars)[20:25]<-c('Elevation','Land_Cover','Forest_Type','Forest_Productivity','Vegetation_Type','SoilpH')

# geonode is soil pH (soilgrids.org). The units are ph * 10.
PredVars$SoilpH <- PredVars$SoilpH/10

# Subset and keep only the layers we'll need.
PredVars <- PredVars[[c(21:23, 25, 31, 39, 47)]]
names(PredVars)
#> [1] "Land_Cover"          "Forest_Type"         "Forest_Productivity"
#> [4] "SoilpH"              "moose1999"           "red_deer1999"       
#> [7] "roe_deer1999"
```

Info

-   ar50 maps are land use classes, including forest productivity (skogbon) and dominating forest tree species (treslag)
-   bio1 to bio19 are worldclim variables. Only bio 10 and12 are used in the publication (bio15 is used in the 'alpine paper'). I'm remaking these variables belowe. so don't need to keep these.
-   resolution is 1km

### Worldclim data and DTM

Worldclim was updated jan 2020, so I can get the bioclim variables again. This is a dataset of interpolated climate variables for the whole world at high resolution (0.5 arc minutes). It is build on data from 1970 onwards and so is not representing any one year. I need to download it as three tiles before merging these together. I will save each tile, but only the two variable bio10 and bio12, Mean Temperature of Warmest Quarter and Annual Precipitation, respectively.

``` r
# first tile
#Norbioclim<-getData('worldclim',var='bio',res=0.5,lon=5,lat=60) # approx 3 min
#Norbioclim <- Norbioclim[[c(10,12)]]
#writeRaster(Norbioclim,'data/large/Norbioclim')
Norbioclim <- raster::stack("data/large/Norbioclim.grd")
raster::plot(Norbioclim)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

``` r
# second tile
#Norbioclim1<-getData('worldclim',var='bio',res=0.5,lon=5,lat=70)
#Norbioclim1 <- Norbioclim1[[c(10,12)]]
#writeRaster(Norbioclim1,'data/large/Norbioclim1')
Norbioclim1 <- raster::stack("data/large/Norbioclim1.grd")
raster::plot(Norbioclim1)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

``` r
#third tile
#Norbioclim2<-getData('worldclim',var='bio',res=0.5,lon=40,lat=70)
#Norbioclim2 <- Norbioclim2[[c(10,12)]]
#writeRaster(Norbioclim2,'data/large/Norbioclim2')
Norbioclim2 <- raster::stack("data/large/Norbioclim2.grd")
raster::plot(Norbioclim2)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

Then I merge these together.

``` r
mergclim<-raster::merge(Norbioclim,Norbioclim1)
mergclim1<-raster::merge(mergclim,Norbioclim2)
raster::plot(mergclim1)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

Now I get a DTM for Norway to be used as an IV, but also to crop the wordclim data.

``` r
#Norelev<-getData('alt',country='NOR', res = 0.5) # 0.86 km2
#names(Norelev) # "NOR_msk_alt"
#writeRaster(Norelev, "data/large/Norelev") # 20mb
Norelev <- raster::stack("data/large/Norelev.grd")
raster::plot(Norelev)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

Then I crop the worldclim data

``` r
cropclim<-raster::crop(mergclim1,Norelev)
raster::plot(cropclim)
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

That took care of the extent. Now I want to put all cells that are outside the DTM as NA also in the climate layers

``` r
Norclimdat<-raster::mask(cropclim,Norelev)
raster::plot(Norclimdat)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

I can put these two together.

``` r
# NorClimElev<-stack(Norclimdat,Norelev)
# names(NorClimElev)<-c("temp", "prec", "elev")
# writeRaster(NorClimElev,'data/NorClimElev') # only 60 mb so shold be able to go on GitHub
NorClimElev <- raster::stack('data/NorClimElev.grd')
NorClimElev
#> class      : RasterStack 
#> dimensions : 1608, 3192, 5132736, 3  (nrow, ncol, ncell, nlayers)
#> resolution : 0.008333333, 0.008333333  (x, y)
#> extent     : 4.6, 31.2, 57.9, 71.3  (xmin, xmax, ymin, ymax)
#> crs        : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 
#> names      : temp, prec, elev 
#> min values :   11,  381,  -27 
#> max values :  164, 2958, 2292
```

Temperature seems to be x10. However, the algebra stuff crashes in knitr for some reason

``` r
#NorClimElev$temp <- NorClimElev$temp/10
#plot(NorClimElev$temp)
```

Change projection and stack
---------------------------

I'm going to use UTM32 projection because the maps look better (more familiar) than with latlon. For the worldclim data and the reindeerSheep data there will also be a resampling. I'm also going to delete some layers I don't need and save it all as IV (overwriting previous name).

``` r
newproj <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# I need to use nearest neighbour for the categorical layers.
cats               <- raster::stack(PredVars$Land_Cover, PredVars$Forest_Type, PredVars$Forest_Productivity) 
cats               <- raster::projectRaster(cats, crs = newproj, method='ngb')

num                <- raster::stack(PredVars$SoilpH,
                                    PredVars$moose1999,
                                    PredVars$red_deer1999,
                                    PredVars$roe_deer1999)
num                <- raster::projectRaster(num, crs = newproj, method='bilinear')

# use 'num' as template to ensure same extent when resampling
reindeerSheep2     <- raster::projectRaster(reindeerSheep, num[[1]], method='bilinear')  
NorClimElev2       <- raster::projectRaster(NorClimElev, num[[1]], method='bilinear')

# Combine
IV                 <- raster::stack(cats, num, reindeerSheep2, NorClimElev2)
```

### Modify IVs

Let's first look at some of them.

``` r
par(mfrow = c(2,2))
raster::plot(IV$Land_Cover, main = "Land cover")
raster::plot(IV$SoilpH, main = "Soil pH")
raster::plot(IV$Forest_Type, main = "Forest type")
raster::plot(IV$Forest_Productivity, main = "Forest productivity")
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%" /> These all need trimming and all exept pH need to be ratified.

``` r
par(mfrow = c(2,2))
raster::plot(IV$moose1999, main = "Moose")
raster::plot(IV$red_deer1999, main = "Red deer")
raster::plot(IV$roe_deer1999, main = "Roe deer")
```

<img src="man/figures/README-unnamed-chunk-17-1.png" width="100%" /> Cervid data seems fine.

#### Clip

We can use any of the cervid data to mask the other layers.

``` r
IV$SoilpH               <- raster::mask(IV$SoilpH,              IV$roe_deer1999)
IV$Land_Cover           <- raster::mask(IV$Land_Cover,          IV$roe_deer1999)
IV$Forest_Type          <- raster::mask(IV$Forest_Type,         IV$roe_deer1999)
IV$Forest_Productivity  <- raster::mask(IV$Forest_Productivity, IV$roe_deer1999)
```

These should now be trimmed to the outline of Norway.

``` r
par(mfrow = c(2,2))
raster::plot(IV$Land_Cover, main = "Land cover")
raster::plot(IV$SoilpH, main = "Soil pH")
raster::plot(IV$Forest_Type, main = "Forest type")
raster::plot(IV$Forest_Productivity, main = "Forest productivity")
```

<img src="man/figures/README-unnamed-chunk-19-1.png" width="100%" /> Reset par.

``` r
par(mfrow = c(1,1))
```

#### Ratify categorical layers

This should naturally happen at this stage, but the RAT, raster attributes table, doesn't carry with the writeRaster function, so I need to do it after reading the file from disc.

Write IV file
-------------

``` r
# writeRaster(IV, 'data/IV', overwrite=TRUE) # 76 MB
myIVs              <- raster:: stack('data/IV.grd')
rm(cats, IV, NorClimElev, Norclimdat, Norbioclim2, Norbioclim, Norbioclim1, NorClimElev2, num, PredVars, reindeerSheep, reindeerSheep2)
names(myIVs)
#>  [1] "Land_Cover"          "Forest_Type"         "Forest_Productivity"
#>  [4] "SoilpH"              "moose1999"           "red_deer1999"       
#>  [7] "roe_deer1999"        "TundraHerbivores"    "temp"               
#> [10] "prec"                "elev"
```

### Ratify v2

So now I can ratify.

#### Forest productivity

``` r
raster::levels(raster::ratify(myIVs$Forest_Productivity))
#> [[1]]
#>   ID
#> 1 11
#> 2 12
#> 3 13
#> 4 18
#> 5 98
#> 6 99
```

We don't want levels 98 and 99. Class 99 = 'ikke registrert'.

``` r
myIVs$Forest_Productivity[myIVs$Forest_Productivity>18]<-NA
myIVs$Forest_Productivity <-   raster::ratify(myIVs$Forest_Productivity)

ratlcp <- raster::levels(myIVs$Forest_Productivity)[[1]]
ratlcp[['Forest_Productivity']] <- 
  c('Unproductive',
    'Low',
    'Medium',
    'High')

levels(myIVs$Forest_Productivity) <- ratlcp
cols2 <- colorRampPalette(c("lightgreen", "darkgreen" ))(4)
rasterVis::levelplot(myIVs$Forest_Productivity, main = "Forest productivity", col.regions = cols2)
```

<img src="man/figures/README-unnamed-chunk-23-1.png" width="100%" /> Ares of NA where there's no forest.

#### Land cover

For this one, let's also look at some statistics

``` r
Land_Cover_stats <- raster::ratify(myIVs$Land_Cover, count = T)
ratlc               <- raster::levels(Land_Cover_stats)[[1]]
ratlc$Land_Cover <- c("Built-up",
                                    "Agricultural",
                                    "Forest",
                                    "Open-natural vegetation",
                                    "Mires",
                                    "Glaciers/Ice/Snow",
                                    "Freshwater",
                                    "Sea",
                                    "NA")
ratlc
#>   ID  COUNT              Land_Cover
#> 1 10   1435                Built-up
#> 2 20   9482            Agricultural
#> 3 30 138689                  Forest
#> 4 50 144319 Open-natural vegetation
#> 5 60  10278                   Mires
#> 6 70   2443       Glaciers/Ice/Snow
#> 7 81  11452              Freshwater
#> 8 82  12467                     Sea
#> 9 99   3237                      NA
```

Norway is mostly open natural vegetation and forest. There's very little buildt up areas and glaciers. For the sake of modeling plant distributions we can group all the obviously unsuitable areas like water, glaciers, and also 'Built-up' I think (red-liste plants are not found on parking lots that often). This should simplify the models considerably. I'll keep the NA as they are, but put the others in a clas == 98 wich I'll call 'other'.

``` r
myIVs$Land_Cover[myIVs$Land_Cover == 10 |
                   myIVs$Land_Cover == 70 |
                   myIVs$Land_Cover == 81 |
                   myIVs$Land_Cover == 82] <- 98
myIVs$Land_Cover <- raster::ratify(myIVs$Land_Cover)
ratlc               <- raster::levels(myIVs$Land_Cover)[[1]]
ratlc$Land_Cover <- c(
                                    "Agricultural",
                                    "Forest",
                                    "Open-natural vegetation",
                                    "Mires",
                                    
                                    
                                    "others",
                                    "NA")
ratlc
#>   ID              Land_Cover
#> 1 20            Agricultural
#> 2 30                  Forest
#> 3 50 Open-natural vegetation
#> 4 60                   Mires
#> 5 98                  others
#> 6 99                      NA
```

``` r

levels(myIVs$Land_Cover) <- ratlc
cols <- colorRampPalette(c(
                           "yellow", 
                           "darkgreen",
                           "tan",
                           "lightblue",
                           "black",
                           "white"
                           ))
rasterVis::levelplot(myIVs$Land_Cover, main = "Land cover", col.regions = cols)
```

<img src="man/figures/README-unnamed-chunk-26-1.png" width="100%" />

#### Forest type

``` r
raster::levels(raster::ratify(myIVs$Forest_Type))
#> [[1]]
#>   ID
#> 1 31
#> 2 32
#> 3 33
#> 4 39
#> 5 98
#> 6 99
```

Deleting class 98 and 99 as above, but also 39 although I', not sure what that is...

``` r
myIVs$Forest_Type[myIVs$Forest_Type>33]<-NA
myIVs$Forest_Type<-raster::ratify(myIVs$Forest_Type)
ratlct<-raster::levels(myIVs$Forest_Type)[[1]]
ratlct[['ForestType']] <-
  c('Coniferous','Deciduous','Mixed')
levels(myIVs$Forest_Type) <- ratlct
cols3 <- colorRampPalette(c("darkgreen", "orange", "blue" ))
rasterVis::levelplot(myIVs$Forest_Type, main= "Forest type", col.regions = cols3)
```

<img src="man/figures/README-unnamed-chunk-28-1.png" width="100%" />

Get species list
----------------

This function produces a list of species that we will later use to harvest occurence data from gbif. More on this later.

``` r
source("./R/spList.R")
mySpecies <- sl()
head(mySpecies)
#> [1] "Botrychium lanceolatum"  "Comastoma tenellum"     
#> [3] "Gentianella campestris"  "Kobresia simpliciuscula"
#> [5] "Primula scandinavica"    "Pseudorchis albida"

# alternatively
# mySpecies(df = TRUE)
```

Occurence data
--------------

To get occurence data I will use the gbif function in the dismo package. It can only handle one species at the time, so I will need to make a for-loop.

``` r
head(mySpecies)
#> [1] "Botrychium lanceolatum"  "Comastoma tenellum"     
#> [3] "Gentianella campestris"  "Kobresia simpliciuscula"
#> [5] "Primula scandinavica"    "Pseudorchis albida"
```

This is my species list with correct spelling (direct from ADB).

To test the functions I will use a shorter list of 10 species.

``` r
mySpecies2 <- mySpecies[1:10]
```

Let's do a test loop without downloading anything, just seeing how many records there are.

``` r
nOccurences_df <- data.frame(species = mySpecies2, 
                  nOccurences = as.numeric(NA))

for(i in 1:length(mySpecies2)){
  myName  <- mySpecies2[i]
  myName2 <- stringr::str_split(myName, " ")[[1]]
  nOccurences_df$nOccurences[i] <- dismo::gbif(myName2[1], myName2[2], download = F) 
}
#> Loading required namespace: jsonlite

nOccurences_df
#>                    species nOccurences
#> 1   Botrychium lanceolatum        5053
#> 2       Comastoma tenellum        6321
#> 3   Gentianella campestris       53391
#> 4  Kobresia simpliciuscula        3392
#> 5     Primula scandinavica         321
#> 6       Pseudorchis albida       20958
#> 7      Pulsatilla vernalis       23679
#> 8    Buglossoides arvensis       42200
#> 9       Anisantha sterilis      113335
#> 10       Sorbus lancifolia         101
```

This shows some of the bias in the occurence data. For example that B arvensis, a super rare plant found almost only on Hovedøya, an island outside of Oslo, has 37k records, whereas P. scandinavica, a relatively common plant, has 321. However, I'm pretty sure occurences within the same 1x1 km cell will only count as one (duplicates removed by the sdm function.)

For the next part I will use two species with a quite low number of records to reduce processing time and test potentias problems due to low sample sizes.

``` r
mySpecies3 <- mySpecies[mySpecies == c("Primula scandinavica", "Kobresia simpliciuscula")]
# write.csv(mySpecies3, 'data/species_list.csv')
```

For fun. lets see what these plants look like.

``` r
list.files("./figures/plants")
```

![Alt text](figures/plants/Kobresia_simpliciuscula_Andrey_Zharkikh_CCBY2.jpg) Picture: *Kobresia simpliciuscula* (Andrey Zharkikh CC-BY 2.0)

![Alt text](figures/plants/Primula_scandinavica_Anders_Kolstad_CCBY4.JPG) Picture: *Primula scandinavica* (Anders Kolstad CC-BY 4.0)

(Note: The picture sizes are 250p and 400p, respectively)

For real this time:

``` r
for(i in 1:length(mySpecies3)){
  myName  <- mySpecies3[i]
  myName2 <- stringr::str_split(myName, " ")[[1]]
  
  assign(
    sub(' ', '_', mySpecies3[i]), 
         dismo::gbif(myName2[1], myName2[2], 
                                        download = T,
                                        geo = T, 
                                        sp = F) 
  )
}
#> 3392 records found
#> 0-300-600-900-1200-1500-1800-2100-2400-2700-3000-3300-3392 records downloaded
#> 321 records found
#> 0-300-321 records downloaded
```

Two new dataframes are put in the environment. They have a lot of columns to start with, so lets get rid of som to make the objects smaller. I only need the species names and the coordinates (perhaps some more, but I can add those later).

``` r
qc <- data.frame(Species = mySpecies3,
                 lon_is_NA =                        as.numeric(NA),
                 lat_NA_when_lon_not          =     as.numeric(NA),
                 lon_is_zero =                      as.numeric(NA),
                 lat_zero_when_lon_not = as.numeric(NA))

for(i in 1:length(mySpecies3)){
  
  
  d <- get(
           sub(' ', '_', mySpecies3[i]))
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
      sub(' ', '_', mySpecies3[i]),  d)
  
  qc[i,2] <- length(w1)
  qc[i,3] <- length(w2)
  qc[i,4] <- length(w3)
  qc[i,5] <- length(w4)
  
}
```

A dataframe called qc tells us what has happened.

``` r
qc
#>                   Species lon_is_NA lat_NA_when_lon_not lon_is_zero
#> 1 Kobresia simpliciuscula       786                   0           0
#> 2    Primula scandinavica       172                   0           0
#>   lat_zero_when_lon_not
#> 1                     0
#> 2                     0
```

We have cut 787 rows from the Kobresia and 172 from Primula, due to missing coordinates.

Now we can turn the dataframes into spatialPointsDataFrames, define the CRS, and plot the points. The dataset comes as lonlat.

``` r
for(i in 1:length(mySpecies3)){
  
  d <- get(
           sub(' ', '_', mySpecies3[i]))
  
  sp::coordinates(d) <- ~lon + lat
  sp::proj4string(d) <- sp::proj4string(raster::raster())
  
  assign(
      sub(' ', '_', mySpecies3[i]),  d)
}
```

``` r
mapview::mapview(Kobresia_simpliciuscula, 
                 map.types = c("Esri.WorldShadedRelief",
                               "Esri.WorldImagery"),
                 cex = 5, lwd = 0,
                 alpha.regions = 0.5,
                 col.regions = "blue")
```

<img src="man/figures/README-unnamed-chunk-39-1.png" width="100%" />

``` r
mapview::mapview(Primula_scandinavica, 
                 map.types = c("Esri.WorldShadedRelief",
                               "Esri.WorldImagery"),
                 cex = 5, lwd = 0,
                 alpha.regions = 0.5,
                 col.regions = "blue")
```

<img src="man/figures/README-unnamed-chunk-40-1.png" width="100%" /> First, notice that Kobresia is called Carex in GBIF, but Kobresia in ADB. This is not a problem and they are recognised as synonyms. The Kobresia is a widespread species, whereas the Primula is endemic to Norway and Sweden. We only need the points that fall on Norway. First we need something to clip against, so we'll get an outline of Norway.

``` r
# outline <- norway()
#saveRDS(outline, "data/large/outline_Norway.RData") # 1.8MB
outline <- readRDS("data/large/outline_Norway.RData")
raster::plot(outline)
```

<img src="man/figures/README-unnamed-chunk-41-1.png" width="100%" />

Now to clip away occurences outside this polygon (can take a few minutes)

``` r

for(i in 1:length(mySpecies3)){
  
  d <- get(
           sub(' ', '_', mySpecies3[i]))
  
  d <- raster::crop(d, outline)
  
  assign(
      sub(' ', '_', mySpecies3[i]),  d)
}
```

Let's see it it worked.

``` r
mapview::mapview(Kobresia_simpliciuscula, 
                 map.types = c("Esri.WorldShadedRelief",
                               "Esri.WorldImagery"),
                 cex = 5, lwd = 0,
                 alpha.regions = 0.5,
                 col.regions = "blue")
```

<img src="man/figures/README-unnamed-chunk-43-1.png" width="100%" />

``` r
mapview::mapview(Primula_scandinavica, 
                 map.types = c("Esri.WorldShadedRelief",
                               "Esri.WorldImagery"),
                 cex = 5, lwd = 0,
                 alpha.regions = 0.5,
                 col.regions = "blue")
```

<img src="man/figures/README-unnamed-chunk-44-1.png" width="100%" /> Looks like it. Now we just need to get this over to UTM32 to match the IV data, and save it on file.

``` r
for(i in 1:length(mySpecies3)){
  
  d <- get(
           sub(' ', '_', mySpecies3[i]))
  
  d <- sp::spTransform(d, myIVs[[1]]@crs)
  
  assign(
      sub(' ', '_', mySpecies3[i]),  d)
}

oDat <- get(sub(' ', '_', mySpecies3[1]))
for(i in 2:length(mySpecies3)){
  oDat <- rbind(oDat, get(sub(' ', '_', mySpecies3[i])))
  }
saveRDS(oDat, 'data/large/oDat.RData')
rm(oDat)
```

``` r
oDat <- readRDS('data/large/oDat.RData')
```

``` r
raster::plot(myIVs$Forest_Productivity)
raster::plot(oDat,add=T)
```

<img src="man/figures/README-unnamed-chunk-47-1.png" width="100%" />

SDM
---

Now we have all we need to make a model. I will save all model objects and use them for making predictions live in the application. The two test species are alpine red listed plants, so I'll use the IV from the alpine paper. Actually, I'll drop bio15 as it was so unimportant. I'll create 1000 random pseudo absences across the goegraphical area (Norway I suppose). To avoid struggling too much with subsetting these strange sdm objects, I will run one set of models for each species, although I know that in sdm::sdm you can specify multiple species.

``` r
library(sdm)
#> Loading required package: sp
#> sdm 1.0-82 (2020-02-03)
```

### SDM-data

``` r
for(i in 1:length(mySpecies3)){
  s    <- unique(oDat$species)[i]
  s2   <- paste0(s, "_d")
  d    <- oDat[oDat$species == s,]
  dat  <- sdm::sdmData(species~temp+prec+SoilpH+TundraHerbivores,
                   train = d,
                   predictors = myIVs,
                   bg = list(n=1000, method = "gRandom"))
  
  # assign(s2, dat)
  sdm::write.sdm(dat, paste0("models/sdmData/", s2), overwrite=TRUE)
}
```

This takes the names over to the gbif approved names (Kobresia to Carex). Perhaps it's for the best. (Hmm, no, that'll make it more difficult to find the model object in the folders when there are tens of models. )

This code could bring the sdmData back in into the environment, but I don't need it now:

``` r
for(i in 1:length(mySpecies3)){
  s    <- unique(oDat$species)[i]
  s2   <- paste0(s, "_d")
  s3   <- paste0("models/sdmData/", s, "_d.sdd")
  assign(s2, sdm::read.sdm(s3))
}
```

### 5 x 3 models

We don't have any independent test data so I'll use bootstrapping to partition test data, and I'll do that 5 times. I shouldn't do much less because there are very few data points for some of the species, like P scandinavica. I think records from inside the same 1km cells will be counted as duplicates and removed. I will use 3 methods as well, resulting in 3 x 5 = 15 models per species. It takes about 1 min to run this.

``` r
for(i in 1:length(mySpecies3)){
  
  s       <- unique(oDat$species)[i]
  file1   <- paste0("models/sdmData/", s, "_d.sdd")
  obj     <- paste0(s, "_m")
  file2   <- paste0("models/sdmModels/", obj)
  d      <- sdm::read.sdm(file1)
  
  mod <- sdm::sdm(.~.,
              data = d, 
              methods = c('glm', 'gam', 'maxent'),   
              replication = c('boot'), n=5)     
  
  sdm::write.sdm(mod, file2, overwrite=TRUE)
  
}
```

Depending on methods, these are some common warning: Warning 1: *prediction from a rank-deficient fit may be misleading* Warning 2: *The response has five or fewer unique values. Are you sure you want to do regression?* I think this was mainly due to random forest method on small sample sizes.

*Technical note: I seems to not be possible to print characters inside the sdm function formula, and therefore we neede seperate sdmData files for each species. I.e. the following fails.*

``` r
for(i in 1:length(mySpecies3)){
  s <- unique(oDat$species)[i]
  d <- get(paste0(s, "_d"))
  mod <- sdm::sdm(paste(s)~.,     # also tried noquote, eval, and  print(quote = F)...
              data = d, 
              methods = c('glm', 'gam', 'rf'),   
              replication = c('boot'), n=3)
    assign(paste0(s, "_mxxx"), 
         mod)
}
```

Let's bring the models back in to the environment. They're 3-20 20 MB each on file, or 17-55 when unzipped in the environment.

``` r
for(i in 1:length(mySpecies3)){
  s      <- unique(oDat$species)[i]
  s2     <- paste0(s, "_m")
  file   <- paste0("models/sdmModels/", s2, ".sdm")
  
  assign(s2, sdm::read.sdm(file))
}
```

``` r
Carex_simpliciuscula_m
#> Loading required package: dismo
#> Loading required package: raster
#> Loading required package: gbm
#> Loaded gbm 2.1.5
#> Loading required package: tree
#> Registered S3 method overwritten by 'tree':
#>   method     from
#>   print.tree cli
#> Loading required package: mda
#> Loading required package: class
#> Loaded mda 0.4-10
#> Loading required package: mgcv
#> Loading required package: nlme
#> 
#> Attaching package: 'nlme'
#> The following object is masked from 'package:raster':
#> 
#>     getData
#> This is mgcv 1.8-28. For overview type 'help("mgcv-package")'.
#> Loading required package: glmnet
#> Loading required package: Matrix
#> Loaded glmnet 3.0-2
#> Loading required package: earth
#> Loading required package: Formula
#> Loading required package: plotmo
#> Loading required package: plotrix
#> Loading required package: TeachingDemos
#> Loading required package: rJava
#> Loading required package: RSNNS
#> Loading required package: Rcpp
#> Loading required package: randomForest
#> randomForest 4.6-14
#> Type rfNews() to see new features/changes/bug fixes.
#> Loading required package: rpart
#> Loading required package: kernlab
#> 
#> Attaching package: 'kernlab'
#> The following objects are masked from 'package:raster':
#> 
#>     buffer, rotated
#> class                                 : sdmModels 
#> ======================================================== 
#> number of species                     :  1 
#> number of modelling methods           :  3 
#> names of modelling methods            :  glm, gam, maxent 
#> replicate.methods (data partitioning) :  bootstrap 
#> number of replicates (each method)    :  5 
#> toral number of replicates per model  :  5 (per species) 
#> ------------------------------------------
#> model run success percentage (per species)  :
#> ------------------------------------------
#> method          Carex_simpliciuscula     
#> ------------------------------ 
#> glm        :            100   %
#> gam        :            100   %
#> maxent     :            100   %
#> 
#> ###################################################################
#> model Mean performance (per species), using test dataset (generated using partitioning):
#> -------------------------------------------------------------------------------
#> 
#>  ## species   :  Carex_simpliciuscula 
#> =========================
#> 
#> methods    :     AUC     |     COR     |     TSS     |     Deviance 
#> -------------------------------------------------------------------------
#> glm        :     0.75    |     0.52    |     0.52    |     1.08     
#> gam        :     0.88    |     0.7     |     0.66    |     0.87     
#> maxent     :     0.9     |     0.71    |     0.69    |     0.9
```

``` r
Primula_scandinavica_m
#> class                                 : sdmModels 
#> ======================================================== 
#> number of species                     :  1 
#> number of modelling methods           :  3 
#> names of modelling methods            :  glm, gam, maxent 
#> replicate.methods (data partitioning) :  bootstrap 
#> number of replicates (each method)    :  5 
#> toral number of replicates per model  :  5 (per species) 
#> ------------------------------------------
#> model run success percentage (per species)  :
#> ------------------------------------------
#> method          Primula_scandinavica     
#> ------------------------------ 
#> glm        :            100   %
#> gam        :            100   %
#> maxent     :            100   %
#> 
#> ###################################################################
#> model Mean performance (per species), using test dataset (generated using partitioning):
#> -------------------------------------------------------------------------------
#> 
#>  ## species   :  Primula_scandinavica 
#> =========================
#> 
#> methods    :     AUC     |     COR     |     TSS     |     Deviance 
#> -------------------------------------------------------------------------
#> glm        :     0.77    |     0.3     |     0.49    |     0.52     
#> gam        :     0.87    |     0.56    |     0.66    |     0.81     
#> maxent     :     0.87    |     0.48    |     0.64    |     0.7
```

GAMs fail when land class is included. RF is not perhaps suitable for such low sample sizes, and when I change from rf to maxent as the third model I removed large white areas in the ensamble map.

``` r
(ev1 <- sdm::getEvaluation(Primula_scandinavica_m, stat = c("AUC", "threshold")))
#>    modelID   AUC  threshold
#> 1        1 0.791 0.17428732
#> 2        2 0.756 0.08328167
#> 3        3 0.714 0.14566953
#> 4        4 0.761 0.24131305
#> 5        5 0.804 0.09338707
#> 6        6 0.889 0.34990089
#> 7        7 0.907 0.01904153
#> 8        8 0.920 0.02030771
#> 9        9 0.865 0.03371827
#> 10      10 0.768 0.23451925
#> 11      11 0.899 0.19887266
#> 12      12 0.863 0.50493407
#> 13      13 0.858 0.36808664
#> 14      14 0.836 0.25700560
#> 15      15 0.904 0.66059893
```

Note the large variation in the threshold value. Density plots show many models have problems discriminating. The maxent models 11-15 are perhaps the most stable. Note also that you can have a 'good model' with high AUC and at the same time a threshold close to 1. Then something is wrong. A very variable threshold value makes it risky to make presense absense maps as I do below.

``` r
(ev2 <- sdm::getEvaluation(Carex_simpliciuscula_m, stat = c("AUC", "threshold")))
#>    modelID   AUC threshold
#> 1        1 0.732 0.4640570
#> 2        2 0.777 0.5327341
#> 3        3 0.744 0.5760212
#> 4        4 0.747 0.5755588
#> 5        5 0.739 0.5536726
#> 6        6 0.861 0.5699635
#> 7        7 0.882 0.6817512
#> 8        8 0.899 0.7605437
#> 9        9 0.873 0.6225270
#> 10      10 0.900 0.6673734
#> 11      11 0.880 0.5872504
#> 12      12 0.893 0.5812584
#> 13      13 0.921 0.6306559
#> 14      14 0.897 0.5963424
#> 15      15 0.898 0.5239558
```

The Carex models are much more stable. Lets compare.

``` r
df <- rbind(ev1, ev2)
df <- data.frame(species   = rep(c("Primula", "Carex"), each = 15),
                 threshold = df$threshold,
                 AUC       = df$AUC )
par(mfrow = c(1,2))
boxplot(threshold~species, data = df)
boxplot(AUC~species, data = df)
```

<img src="man/figures/README-unnamed-chunk-58-1.png" width="100%" /> I would not make a presence absence map for the Primula with threshold values being so fluctuating. Either I make probability maps for all, or I drop some species out. The AUCs is more stable, but some replicated partitioning seems warranted. The trade off here is with processing time for this thing to be able to run live predictions.

#### Ensemble

Lets put the 5 x 3 models together and make a map of the current habitat suitability, using ensamble. Runtime approx. 1 min.

``` r
for(i in 1:length(mySpecies3)){
  s <- unique(oDat$species)[i]
  d <- get(paste0(s, "_m"))
  
  mod <- sdm::ensemble(d,
              newdata = myIVs, 
              filename = paste0("models/predictions/", s, "_ens.img"), overwrite=TRUE,   
              setting = list(method='weighted', stat = 'AUC')
              )     
  assign(paste0(s, "_ens"), 
         mod)
}
```

Lets plot them.

``` r
Primula_scandinavica_ens <- raster::raster("models/predictions/Primula_scandinavica_ens.img")
raster::plot(Primula_scandinavica_ens, main = "Habitat suitability for Primula scandinavica")
raster::plot(oDat[oDat$species == "Primula_scandinavica",], add=T, cex = 0.8, pch = 1)
```

<img src="man/figures/README-unnamed-chunk-60-1.png" width="100%" /> The realised nishe is considerably smaller then the desirable or fundamental nishe.

``` r
Carex_simpliciuscula_ens <- raster::raster("models/predictions/Carex_simpliciuscula_ens.img")
raster::plot(Carex_simpliciuscula_ens, main = "Habitat suitability for Carex simpliciuscula")
raster::plot(oDat[oDat$species == "Carex_simpliciuscula",], add=T, cex = 0.2)
```

<img src="man/figures/README-unnamed-chunk-61-1.png" width="100%" />

For most people I think presence-absence is more understandable than probability of occurence (but see discussion about threshold values above). Let's use mean threshold value for defining p-a and plot a discrete map instead. Also, I like to think of this as habitat suitability rather than distribution per se, and we can name the categories 'suitable' or 'unsuitable' rather then 'present' and 'absent'. Then the next map of 'collonisation' and 'extinction' could be classes as 'new habitat' and 'lost habitat'.

#### Presence-absence map

``` r
for(i in 1:length(mySpecies3)){
  s <- as.name(paste0(unique(oDat$species)[i], "_m"))
  ev <- paste0(unique(oDat$species)[i], "_ev")
  assign(ev, sdm::getEvaluation(eval(s), 
                   stat = c('AUC', 'threshold', 'TSS'), opt = 2)
  )
}
```

``` r
for(i in 1:length(mySpecies3)){
  ens    <- get(paste0(unique(oDat$species)[i], "_ens"))
  ev     <- get(paste0(unique(oDat$species)[i], "_ev"))
  s      <- paste0(unique(oDat$species)[i], "_pa")
  
  
  ens[]   <- ifelse(ens[] >= mean(ev$threshold), 1, 0)
  assign(s, ens)
}
```

``` r
Primula_scandinavica_pa2 <- raster::ratify(Primula_scandinavica_pa)
lev <- raster::levels(Primula_scandinavica_pa2)[[1]] 
lev$pa <- c("Unsuitable habitat", "Suitable habitat")
levels(Primula_scandinavica_pa2) <- lev
source("R/norway.R")
nor <- norway(lonlat = FALSE)
library(latticeExtra)
#> Loading required package: lattice
rasterVis::levelplot(Primula_scandinavica_pa2, 
                   margin=F,
                   main="Primula scandinavica",
                   scales=list(draw=F),
                   col.regions= c("grey", "blue"))+
  layer(sp::sp.polygons(nor),col=grey(0.5))+
  layer(sp::sp.points(oDat[oDat$species == "Primula_scandinavica",], pch = 1, cex = 0.8, col = "black", lwd = 0.5, alpha = 0.7))
```

<img src="man/figures/README-unnamed-chunk-64-1.png" width="100%" /> This map has changes drastically between runs!!, indicating that 5\*3 is not enough for it to stabilise.

``` r
Carex_simpliciuscula_pa2 <- raster::ratify(Carex_simpliciuscula_pa)
lev <- raster::levels(Carex_simpliciuscula_pa2)[[1]] 
lev$pa <- c("Unsuitable habitat", "Suitable habitat")
levels(Carex_simpliciuscula_pa2) <- lev

rasterVis::levelplot(Carex_simpliciuscula_pa2, 
                   margin=F,
                   main="Carex simpliciuscula",
                   scales=list(draw=F),
                   col.regions= c("grey", "blue"))+
  layer(sp::sp.polygons(nor),col=grey(0.5))+
layer(sp::sp.points(oDat[oDat$species == "Carex_simpliciuscula",], pch = 1, cex = 0.8, col = "black", lwd = 0.5, alpha = 0.7))
```

<img src="man/figures/README-unnamed-chunk-65-1.png" width="100%" /> Nice map. Easy to read and understand to laymen I think.

### Best candidate model

To make predictions from the 5 x 3 models we need to run the ensamble function again, this time with altered IVs in the newdata argument. The function is probably too slow to run on the fly. This I need to test on the io server later. We could choose one method, eg maxent, and just use that for all species, and with only one replication. Then we could use the predict function in raster which will be quicker. Or we can look at these 15 models we have generated and choose the best one from there. That would be a bit safer.

``` r
Primula_scandinavica_ev
#>    modelID   AUC  threshold   TSS
#> 1        1 0.791 0.17428732 0.509
#> 2        2 0.756 0.08328167 0.485
#> 3        3 0.714 0.14566953 0.408
#> 4        4 0.761 0.24131305 0.482
#> 5        5 0.804 0.09338707 0.584
#> 6        6 0.889 0.34990089 0.627
#> 7        7 0.907 0.01904153 0.663
#> 8        8 0.920 0.02030771 0.708
#> 9        9 0.865 0.03371827 0.686
#> 10      10 0.768 0.23451925 0.604
#> 11      11 0.899 0.19887266 0.662
#> 12      12 0.863 0.50493407 0.657
#> 13      13 0.858 0.36808664 0.623
#> 14      14 0.836 0.25700560 0.521
#> 15      15 0.904 0.66059893 0.727
```

Based on AUC, model nr ...

``` r
which.max(Primula_scandinavica_ev$AUC)
#> [1] 8
```

is the best for the Primula. But the threshold is super low.

``` r
Carex_simpliciuscula_ev
#>    modelID   AUC threshold   TSS
#> 1        1 0.732 0.4640570 0.519
#> 2        2 0.777 0.5327341 0.551
#> 3        3 0.744 0.5760212 0.485
#> 4        4 0.747 0.5755588 0.537
#> 5        5 0.739 0.5536726 0.509
#> 6        6 0.861 0.5699635 0.615
#> 7        7 0.882 0.6817512 0.670
#> 8        8 0.899 0.7605437 0.710
#> 9        9 0.873 0.6225270 0.630
#> 10      10 0.900 0.6673734 0.699
#> 11      11 0.880 0.5872504 0.658
#> 12      12 0.893 0.5812584 0.696
#> 13      13 0.921 0.6306559 0.754
#> 14      14 0.897 0.5963424 0.672
#> 15      15 0.898 0.5239558 0.681
```

Based on AUC, model nr ...

    #> [1] 13

is the best for the Carex.

Let's compare the ensemble with the best candidate model.

``` r
for(i in 1:length(mySpecies3)){
  s <- unique(oDat$species)[i]
  s2 <- paste0(s, "_m")
  s3 <- get(paste0(s, "_ev"))
  file <- paste0("models/sdmModels/", s, "_bcm")
  best <- s3$modelID[s3$AUC == max(s3$AUC)]
  mod <- get(s2)
  mod2 <- mod[[best]]
  sdm::write.sdm(mod2, file, overwrite=TRUE)

  
  
  mod3 <- raster::predict(mod2,
              object = myIVs, 
              filename = paste0("models/predictions/", s, "_best.img"), 
              overwrite=TRUE)
                
  assign(paste0(s, "_best"), 
         mod3)
}
```

``` r
par(mfrow=c(2,2))
raster::plot(Primula_scandinavica_best, main = "Best candidate model\nfor Primula")   
raster::plot(Primula_scandinavica_ens, main = "Ensemble model\nfor Primula")
raster::plot(Carex_simpliciuscula_best, main = "Best candidate model\nfor Carex")   
raster::plot(Carex_simpliciuscula_ens, main = "Ensemble model\nfor Carex")
```

<img src="man/figures/README-unnamed-chunk-71-1.png" width="100%" /> The models change considerably between runs, sometimes ending up with incomplete coverage (white areas on the map). The strange thing is that when that happens, it's the same for both species.

``` r
tt  <- Carex_simpliciuscula_best > -Inf
tt2 <- Primula_scandinavica_best > -Inf
par(mfrow = c(1,2))
raster::plot(tt, main = "Carex")
raster::plot(tt2, main = "Primula")
```

<img src="man/figures/README-unnamed-chunk-72-1.png" width="100%" />

The following script makes these comparisons with presence absence map, but thats not so important.

``` r
for(i in 1:length(mySpecies3)){
  p    <- get(paste0(unique(oDat$species)[i], "_best"))
  ev     <- get(paste0(unique(oDat$species)[i], "_ev"))
  th     <- ev$threshold[ev$AUC == max(s3$AUC)]
  s      <- paste0(unique(oDat$species)[i], "_pab")
  
  
  p[]   <- ifelse(p[] >= mean(ev$threshold), 1, 0)
  assign(s, p)
}


# Find divergent areas for the Primula

ps <- Primula_scandinavica_pa2 - Primula_scandinavica_pab
ps <- raster::ratify(ps)
ratlcp <- raster::levels(ps)[[1]]
ratlcp[['errors']] <- 
  c('new',
    'same',
    'lost')
levels(ps) <- ratlcp


# Find divergent areas for the Carex

cs <- Carex_simpliciuscula_pa2 - Carex_simpliciuscula_pab
cs <- raster::ratify(cs)
ratlcp <- raster::levels(cs)[[1]]
ratlcp[['errors']] <- 
  c('new',
    'same',
    'lost')
levels(cs) <- ratlcp


cols <- colorRampPalette(c("lightgreen", "white", "red" ))

# save as pictures to save knitr time
png("figures/comparisons/Primula.png", height = 380, width = 380, units = "px")
rasterVis::levelplot(ps, col.regions = cols, main = "Primula")+
  layer(sp::sp.polygons(nor),col=grey(0.5))
dev.off()
png("figures/comparisons/Carex.png", height = 380, width = 380, units = "px")
rasterVis::levelplot(cs, col.regions = cols, main = "Carex")+
  layer(sp::sp.polygons(nor),col=grey(0.5))
dev.off()
```

![Alt txt](figures/comparisons/Primula.png) ![Alt txt](figures/comparisons/Carex.png) In conclusion I think we ideally should to replicate the partitioning because of the low sample sizes, but if that takes too long to process it not impossible to do the 'single best' approach, ie to find the best out of 15 candidate models and use that. This is better than choosing on method, ie maxent, and one of a few replicates, because our test showed that on some runs actually the GAMs were better. Then should use the probabilities and not the presence-absence maps, because of the large variation in the threshold value between model runs. If the 'single best models differs from the ensemble that may not be a problem as the change in predictions when altering the IVs (the input values) should be qualitatively similar, which is good enough for our purpose. Also, the ensemble prediction includes a lot of crap models so it's not neccessarily better.

### Replicated single method

This approach chooses one method, maxent, and does five replicated partitionings, and then the raster::predict function to make maps. I'll keep these models to to a speed test later in the actuall Shiny app on the io server.

``` r
for(i in 1:length(mySpecies3)){
  
  s       <- unique(oDat$species)[i]
  file1   <- paste0("models/sdmData/", s, "_d.sdd")
  obj     <- paste0(s, "_ms")
  obj2     <- paste0(s, "_5maxent")
  file2   <- paste0("models/sdmModels/", obj)
  d      <- sdm::read.sdm(file1)
  
  mod <- sdm::sdm(.~.,
              data = d, 
              methods = 'maxent',   
              replication = c('boot'), n=5)     
  
  sdm::write.sdm(mod, file2, overwrite=TRUE)
  assign(obj, mod)
  
  mod2 <- raster::predict(mod,
              object = myIVs,
              mean = T,
              filename = paste0("models/predictions/", obj2, ".img"), 
              overwrite=TRUE)
  assign(obj2, mod2)

  ev     <- sdm::getEvaluation(mod, stat = "threshold")
  th     <- mean(ev$threshold)
  s      <- paste0(unique(oDat$species)[i], "_5maxent_pa")
  mod2[] <- ifelse(mod2[] >= th, 1, 0)
  assign(s, mod2)

}
#> Loading required package: parallel
```

Check coverage. Some runs gives white spots on this map.

``` r
tt <- Primula_scandinavica_5maxent > -Inf
raster::plot(tt)
```

<img src="man/figures/README-unnamed-chunk-75-1.png" width="100%" />

``` r
raster::plot(Primula_scandinavica_5maxent, main = "Habitat suitability\nPrimula, 5x maxent")
```

<img src="man/figures/README-unnamed-chunk-76-1.png" width="100%" />

``` r
getEvaluation(Primula_scandinavica_ms)
#>   modelID   AUC   COR  Deviance   TSS
#> 1       1 0.973 0.560 0.6288563 0.885
#> 2       2 0.889 0.513 0.5720289 0.729
#> 3       3 0.908 0.464 0.7771705 0.717
#> 4       4 0.938 0.509 0.6519923 0.773
#> 5       5 0.923 0.489 0.5709049 0.707
```

Compare 'best candidate' with 'five best maxent'

``` r
diffpa <- Primula_scandinavica_best - Primula_scandinavica_5maxent_pa
raster::plot(diffpa)
```

<img src="man/figures/README-unnamed-chunk-78-1.png" width="100%" /> Obs, plot changes between updates.

### Response curves

Let make response curve plots for each species and save those as well

``` r


for(i in 1:length(mySpecies3)){
  s      <- unique(oDat$species)[i]
  s2     <- paste0(s, "_m")
  s3     <- paste0("models/rcurves/", s, "_rcurves.png")
  d      <- get(s2)
  d2     <- d
 
p <- sdm::rcurve(d2, ylab="Habitategnethet\nHabitat suitability", 
                 xlab = "",
                 main = "")
# not sure how to rename the variables:
#labs <- c("Temperature", "Precipitation", 
#           "Soil pH", "Sheep and reindeer")
# p2 <- p + facet_grid(labeller = labeller(variable = labs))

png(filename = s3,
    width = 480, height = 380, units = "px")
  print(p)
  dev.off()

}
#> Loading required package: ggplot2
#> 
#> Attaching package: 'ggplot2'
#> The following object is masked from 'package:latticeExtra':
#> 
#>     layer
#> The following object is masked from 'package:kernlab':
#> 
#>     alpha
#> The following object is masked from 'package:randomForest':
#> 
#>     margin
#> The id argument is not specified; The modelIDs of 15 successfully fitted models are assigned to id...!
#> The id argument is not specified; The modelIDs of 15 successfully fitted models are assigned to id...!
```

![Alt text](models/rcurves/Primula_scandinavica_rcurves.png) ![Alt text](models/rcurves/Carex_simpliciuscula_rcurves.png)

### Variable importance

I would also like to get a plot of the variable importance. This script uses the combined/averaged impirtances from the 5\*3 models.

``` r
# create empty list
varimp <- list()

# list of IVs (different for alpine and forest species). 
# Make sure to put them in the correct order. 
# Find the order with > getVarImp(yourModel, id = 1)@varImportance
IV <- c("temp", "prec", "SoilpH", "TundraHerbivores")

# Create an empty varimp table
df1 <- data.frame(variables = IV,
                  corTest = as.numeric(NA),
                  AUCtest = as.numeric(NA))
                  
for(i in 1:length(mySpecies3)){
  s      <- unique(oDat$species)[i]
  s2     <- paste0(s, "_m")
  d      <- get(s2)
  tab    <- d@run.info

for(t in 1:max(tab$modelID)){
  r <- length(varimp)+1
  ifelse(tab$success[t]==TRUE,
           varimp[[r]]          <-sdm::getVarImp(d,id=t)@varImportance,
           varimp[[r]]          <- df1)
         
           varimp[[r]]$species  <-tab$species[t]
           varimp[[r]]$method   <-tab$method[t]
           varimp[[r]]$repid    <-tab$replicationID[t]      
         
   #if(tab$success[t]==FALSE) return(print(paste('Model failiure run ',t)))
}
}

varimp<-do.call('rbind',varimp)
rm(s, s2, d, tab, df1)
source("R/se.R")
varimpmean <- aggregate(data = varimp,
                        corTest ~ species + variables,
                        FUN = function(x) c(mean = mean(x, na.rm=T), se = se(x)))
varimpmean <- do.call(data.frame, varimpmean)
head(varimpmean)
#>                species variables corTest.mean  corTest.se
#> 1 Carex_simpliciuscula      prec    0.3956200 0.015210185
#> 2 Primula_scandinavica      prec    0.5914267 0.046246521
#> 3 Carex_simpliciuscula    SoilpH    0.0354000 0.005247094
#> 4 Primula_scandinavica    SoilpH    0.2610533 0.038254301
#> 5 Carex_simpliciuscula      temp    0.3300467 0.020798081
#> 6 Primula_scandinavica      temp    0.1030133 0.023270557
```

Changing the variable names for axis tick labels

``` r
varimpmean$variables <- plyr::revalue(varimpmean$variables, c(
              prec         = "Årlig nedbørsmengde\n
                              Annual precipitation",
              SoilpH       = "pH i jorden\nSoil pH",
              temp         = "Gjennomsnittemperatur i varmeste kvartal\n
                              Mean temperature in warmenst quarter",
              TundraHerbivores = "Tetthet av sau og reinsdyr\n
                                  Sheep and reindeer densities"))
```

We can make one plot per species and show it in the app. That way the user will know what slider should induce the biggest effect.

``` r
library(ggplot2)

for(i in 1:length(unique(varimpmean$species))){
  
  d <- varimpmean[varimpmean$species==unique(varimpmean$species)[i],]
  s <- paste0("./models/varimp/",
              unique(varimpmean$species)[i],
              ".png")
  
  p <- ggplot2::ggplot(data = d)+
  geom_bar(aes(y = corTest.mean, 
               x = variables, 
               fill = variables), 
           stat = "identity", 
           colour = "black")+
  coord_flip()+
  ylab("Variabelviktigheten\nVariable importance")+
  xlab("")+
  theme_minimal()+
  theme(axis.text.y = element_text(size = 10))+
  theme(legend.position="none")+
  geom_errorbar(aes(x = variables, ymin=corTest.mean-corTest.se, ymax=corTest.mean+corTest.se), width=.2)

  png(filename = s,
    width = 480, height = 380, units = "px")
  print(p)
  dev.off()
}
```

![Alt text](models/varimp/Primula_scandinavica.png) ![Alt text](models/varimp/Carex_simpliciuscula.png)
