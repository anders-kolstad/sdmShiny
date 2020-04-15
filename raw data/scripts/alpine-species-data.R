#' Get species data
#' 
#' This scripts collects and collates spatial data on red listed alpine planst in Norway.
#' 
#' 
#' 


#library(dismo) # (>= 1.1.4)


# Species list
sl <- c("Botrychium lanceolatum",
  'Comastoma tenellum',
  'Gentianella campestris',
  'Kobresia simpliciuscula',
  'Primula scandinavica',
  'Pseudorchis albida',
  'Pulsatilla vernalis')

# Any synonyms? 
# Gbif should recognize synonyms through the 'Catalogue of Life. 
# For example, P. vernalis has 5 synonyms and three infraspecies, 
# but all data appear when searching for P. vernalis.

# test download
BL <- dismo::gbif("Botrychium", "lanceolatum", download = T, geo = T, sp = F)
class(BL)

# backup
bl <- BL

# remove NA's
w <- which(is.na(BL$lon))
if(length(w) != 0) BL <- BL[-w,]
BL <- BL[-w,]
w <- which(is.na(BL$lat))
if(length(w) != 0) BL <- BL[-w,]


BL$species <- "Botrychium lanceolatum"
BL <- BL[,c("lon", "lat","species")]
head(BL)
sp::coordinates(BL) <- ~lon + lat
library(mapview)
mapview::mapview(BL, map.types = c("Esri.WorldShadedRelief", "OpenStreetMap.DE"))
