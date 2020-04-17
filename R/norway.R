#' Create an outline of Norway
#' 
#' \code{norway} returs the outline of Norway
#' 
#' Here comes the details...
#' 
#' @param cs Logical. Should the funtion return a raster in the lonlat coordinate system (defult) or "UTM32"
#' @return Returns a Large SpatialPolygonsDataFrame
#' @import raster
#' @import sp
#' @examples 
#' \dontrun{
#' norway()
#' }
#' @export


norway <- function(cs = "lonlat"){
  norway_lonlat <- raster::getData('GADM', country='NOR', level=0)
  norway_UTM32   <-sp::spTransform(norway_lonlat,"+proj=utm +zone=32")
  ifelse(cs == "lonlat", return(norway_lonlat), return(norway_UTM32))
  }