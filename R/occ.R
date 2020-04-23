#' Get occurence data from gbif
#' 
#' Take a species list and harvest gbif data
#' @examples
#' \dontrun{
#'        x <- y
#'        # put example here
#'     }

#' @import dismo
#' @import sp
#' @import mapview




# Forest plant occurences from GBIF

#file <- "/home/anders/Documents/R/Git projects/sdmShiny/raw data/0008843-181108115102211.zip"
#out <- finch::dwca_read # ERROR finch wont download

#keysL  <-sapply(as.character(herblichen$Vitenskapelig.navn),function(x) name_backbone(x,rank='species')$speciesKey)

#paste(keysL, collapse = ',')#Copy and paste to below

#odLichen <- occ_download('taxonKey = 2609133,2607387,3408876,2603661,2602561,2609349,3422592,3433876,2601243,2608140,5260765,7083298,2605872,7425899,5260770,9198278,5258394,8704544,3397573,3429282,3390552,2609409,2609180,5260761,2607725,8335421,5260565,2599762,3429909,3389602'
#                        ,'country = NO','hasCoordinate = TRUE',
#                      user='jamesspeed',pwd='*****',email='*****')
#occ_download_meta(odLichen)
#gbif_citation(occ_download_meta(odLichen))# GBIF Occurrence Download https://doi.org/10.15468/dl.pl144i Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2018-11-26"





# test download
#dismo::gbif("Botrychium", "lanceolatum", download = F)
#dismo::gbif("Primula", "scandinavica", download = F)

#ps <- dismo::gbif("Primula", "scandinavica", download = T, geo = T, sp = F)

# this will require a for loop

#class(ps)

# backup
#PS <- ps

# remove NA's
#w <- which(is.na(ps$lon))
#if(length(w) != 0) ps <- ps[-w,]
#w <- which(is.na(ps$lat))
#if(length(w) != 0) ps <- ps[-w,]
#w <- which(ps$lon == 0)
#if(length(w) != 0) ps <- ps[-w,]
#w <- which(ps$lat == 0)
#if(length(w) != 0) ps <- ps[-w,]


#ps$species <- 'Primula scandinavica'
#ps <- ps[,c("lon", "lat","species")]
#head(ps)
#sp::coordinates(ps) <- ~lon + lat
#sp::proj4string(ps) <- sp::proj4string(raster::raster())
##library(mapview)
#mapview::mapview(ps, 
#                 map.types = c("Esri.WorldShadedRelief",
#                               "Esri.WorldImagery"),
#                 cex = 5, lwd = 0,
#                 alpha.regions = 0.5,
#                 col.regions = "blue")#

#library(margrittr) # forward pipe
#leaflet::leaflet(data = bl) %>% 
#  addTiles(group = "OSM", 
#  options = providerTileOptions(minZoom = 2, maxZoom = 100)) # %>% 
#  addCircleMarkers(lat = ~lat, lng = ~lon,
#                   color = "blue",
#                   radius = 0.05)


