#' Get species lists
#' 
#' This scripts produces a list of the species that we will quiery gbif for occurence data.
#' \code{sl()} produces a species list and a dataframe with the same list pluss some more
#' @param df Logical. Whether to return a dataframe which also includes the species group in a second column. Default is df = FALSE which returns acharacter string of names (type: vector).
#' 
#' @examples
#' \dontrun{
#'        myList <- sl()
#'     }



#'@import readxl 
#'@export


sl <- function(df = FALSE){
# for the alpine plants, it a short list, from Speed and Austrheim 2017 Biological Conservation:
alp_sl <- c("Botrychium lanceolatum",
        'Comastoma tenellum',
        'Gentianella campestris',
        'Kobresia simpliciuscula',   # NT in 2010, but LC in 2015
        'Primula scandinavica',
        'Pseudorchis albida',
        'Pulsatilla vernalis')
# All these, I think, are on the red list due to reduced grazing. None are here because of overgrazing -watch out for negative parameter estimates for cervid densities!


# What about synonyms? 
# Gbif should recognize synonyms through the 'Catalogue of Life. 
# For example, P. vernalis has 5 synonyms and three infraspecies, 
# but all data appear when searching for P. vernalis.


# Forest plants and lichens species (Speed et al submitted 2020 FEM)
# This is not actually neccessary as I can also download the DwC archive from the doi cited in th paper. However, that makes it harder to retrace our steps in the future. 

# Note: ABD also offer shapefiles that are updated daily. But we will stick with the more generic solution with gbif which isn't limited to redlisted species.

# The redlist is downloaded from here:
# https://artsdatabanken.no/Rodliste2015/sok?vurderings%u00e5r=2015

# First, the (v)ascular (f)orest plants (karplanter) in the (r)ed (l)ist:
vf_rl <- readxl::read_excel("./raw data/vasc_rl.xlsx")
vf_sl<-droplevels(vf_rl[grep('*beit*',vf_rl$Påvirkningsfaktorer,ignore.case=T),])
vf_sl <- vf_sl[,'Vitenskapelig navn']

#View(vf_rl) 
# Still quite a long list

mySpList <- c(alp_sl, vf_sl$`Vitenskapelig navn`)
mySpListDF <- data.frame(mySpList, "type" = c(rep("alpine", length(alp_sl)), rep("forest", nrow(vf_sl))))

ifelse(df == FALSE, return(mySpList), return(mySpListDF))

  
}
