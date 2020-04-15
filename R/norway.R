#' Create an outline of Norway
#' 
#' \code{norway} returs the outline of Norway
#' 
#' Here comes the details...
#' 
#' @param 
#' @return Returns a Large SpatialPolygonsDataFrame
#' @import raster
#' @import sp
#' @examples 
#' \dontrun{
#' norway()
#' }
#' @export


norway <- function(){
  norway_longlat <<- raster::getData('GADM', country='NOR', level=0)
  norway_UTM32   <<-sp::spTransform(norway_longlat,"+proj=utm +zone=32")
}