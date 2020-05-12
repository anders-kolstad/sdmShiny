#background_data


# Get 5000 random physical plant records

rgbif::occ_count(country = "NO", 
                 taxonKey=7707728,              
                 basisOfRecord = "PRESERVED_SPECIMEN", 
                 year='1989,2020'
                 )

mybg <- rgbif::occ_search(country = "NO", 
                          phylumKey=7707728,
                          basisOfRecord = "PRESERVED_SPECIMEN", 
                          year='1989,2020', 
                          return='data',
                          hasGeospatialIssue=FALSE,
                          hasCoordinate=TRUE,
                          limit=99990) # max 100.000

saveRDS(mybg, 'data/background.RData')

plot(mybg$decimalLongitude, mybg$decimalLatitude)
summary(mybg$decimalLatitude      )

hist(mybg$year)
sp::coordinates(mybg) <- ~decimalLongitude + decimalLatitude
sp::proj4string(test) <- sp::proj4string(raster::raster())
mapview::mapview(test, 
                 map.types = c("Esri.WorldShadedRelief",
                               "Esri.WorldImagery"),
                 cex = 5, lwd = 0,
                 alpha.regions = 0.5,
                 col.regions = "blue")
