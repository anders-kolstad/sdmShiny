#' Standard error of the mean
#' 
#' Calculates SE
#' 
#' @examples 
#' \dontrun{
#' # craete a vector with sd = 1 and n = 150. SE should be 1/sqrt(100) = 0.08
#' a <- rnorm(150,5,1)
#' se() 
#' 
#' # it automatically removed NA's
#' a2 <- a
#' a2[3] <- NA
#' se(a2)
#' }
#' @export

se <- function(x) sd(x, na.rm = T)/sqrt(length(x))