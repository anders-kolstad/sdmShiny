
-   [sdmShiny](#sdmshiny)
    -   [Project homepage](#project-homepage)
    -   [Installation](#installation)
    -   [About the repo](#about-the-repo)

<!-- README.md is generated from README.Rmd. Please edit that file -->
sdmShiny
========

<!-- badges: start -->
<!-- badges: end -->
Last update:

``` r
Sys.time()
#> [1] "2020-05-20 12:49:56 CEST"
```

Follow this link to go to the Shiny application: <https://anderskolstad.shinyapps.io/sdmShiny/>

Project homepage
----------------

All the documentation for this project is found on [this webpage](https://anders-kolstad.github.io/sdmShiny/).

Installation
------------

Developers should clone the repo and work from there.

About the repo
--------------

Large files around 100mb or greater, or unessential raster files etc., are in the data/large/ folder whihc is not pushed (it's in the .gitignore file) and therefore only exists locally with Anders. The same for RData files intill I find a way to load them without getting magic number errors (readRDS don't work with knitr). R markdown files in the root folder produces html files that goes to the docs folder - these make up the github pages website (except the README file which makes up this page).
