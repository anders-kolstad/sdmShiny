
<!-- README.md is generated from README.Rmd. Please edit that file -->
sdmShiny
========

<!-- badges: start -->
<!-- badges: end -->
This project is for disseminating the species distribution modeling work done in James Speed's group at the NTNU University Museum. We will use web-based Shiny apps to present distribution maps of several species and allow these to change with the predictions of the SDM as the user tweaks the parameters for climate and herbivory. ![The Shiny app will look something like this](figures/app.png)

Installation
------------

Install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
#devtools::install_github("anders-kolstad/sdmShiny")
```

Example
-------

Here is an example application: <https://anderskolstad.shinyapps.io/demoSDM/>

Documentation
=============

This section explains the workflow that ended up the the shiny app. The r project is arranged as an r package. All r code such as functions are in the 'data' folder.

1) Get species list
-------------------

This is the first function. It produces a list of species that we will later use to harvest occurence data from gbif. ...

``` r
mySpecies <- sl()
head(mySpecies)
# alternatively
# mySpecies(df = TRUE)
```

Get environmental data (explanatory variables)
==============================================

A .gri file with already collated environmental data was downloaded from NTNU box: <https://ntnu.app.box.com/s/wcmr0dgoyz2yu6ielw6er1pm7h0gaisa/file/393633279036>

``` r
env <- raster::stack("data/large/PredictorVariables.grd")
names(env)
#>  [1] "bio1_16"                   "bio2_16"                  
#>  [3] "bio3_16"                   "bio4_16"                  
#>  [5] "bio5_16"                   "bio6_16"                  
#>  [7] "bio7_16"                   "bio8_16"                  
#>  [9] "bio9_16"                   "bio10_16"                 
#> [11] "bio11_16"                  "bio12_16"                 
#> [13] "bio13_16"                  "bio14_16"                 
#> [15] "bio15_16"                  "bio16_16"                 
#> [17] "bio17_16"                  "bio18_16"                 
#> [19] "bio19_16"                  "NOR_msk_alt"              
#> [21] "ar50_type"                 "ar50_treslag"             
#> [23] "ar50_skogbon"              "ar50_veget"               
#> [25] "geonode_phihox_m_sl2_250m" "moose1949"                
#> [27] "moose1959"                 "moose1969"                
#> [29] "moose1979"                 "moose1989"                
#> [31] "moose1999"                 "moose2009"                
#> [33] "moose2015"                 "red_deer1949"             
#> [35] "red_deer1959"              "red_deer1969"             
#> [37] "red_deer1979"              "red_deer1989"             
#> [39] "red_deer1999"              "red_deer2009"             
#> [41] "red_deer2015"              "roe_deer1949"             
#> [43] "roe_deer1959"              "roe_deer1969"             
#> [45] "roe_deer1979"              "roe_deer1989"             
#> [47] "roe_deer1999"              "roe_deer2009"             
#> [49] "roe_deer2015"
```
