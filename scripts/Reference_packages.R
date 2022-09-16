################################################################################
## R script to export R packages citation used in the Hyperspectral Soil Mapping     
## Version: 1, Sep 2022                                                         
## Author: Yuri Andrei Gelsleichter
## License: CC-BY-NC-SA                                                         
################################################################################

#' ### Auto setwd
library(rstudioapi); current_path <- getActiveDocumentContext()$path; setwd(dirname(current_path)) ### https://eranraviv.com/r-tips-and-tricks-working-directory/
getwd()
gc(); rm(list=ls())

# For import, export dataset:  
library(readr)
library(fst)
library(rstudioapi)

# Data management:  
library(dplyr)
library(magrittr)
library(reshape)

# Graphs: 
library(ggplot2)
library(RStoolbox)
library(ggspatial)
library(RColorBrewer)
library(hrbrthemes)
library(viridis)
library(gridExtra)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# Tune Random Forest model: 
library(tidymodels)
library(themis)
library(caret)

# Manage raster data:  
library(rgdal)
library(raster)
library(terra)

library(RSAGA)

# Model and goodness of fit measures: 
library(randomForest)
library(ithir)

pk <- c(
  "readr",
  "fst",
  "rstudioapi",
  "dplyr",
  "magrittr",
  "reshape",
  "ggplot2",
  "RStoolbox",
  "ggspatial",
  "RColorBrewer",
  "hrbrthemes",
  "viridis",
  "gridExtra",
  "DiagrammeR",
  "DiagrammeRsvg",
  "rsvg",
  "tidymodels",
  "themis",
  "caret",
  "rgdal",
  "raster",
  "terra",
  "RSAGA",
  "randomForest",
  "ithir")

# knitr::write_bib(c(.packages()), "packages.bib")
knitr::write_bib(c(pk), "manuscript_packages.bib")


# For import, export dataset: readr, fst, rstudioapi; data management: dplyr, magrittr, reshape;  
# graphs: ggplot2, RStoolbox, ggspatial, RColorBrewer, hrbrthemes, viridis, gridExtra, DiagrammeR, DiagrammeRsvg, rsvg; 
# tune Random Forest model: tidymodels, themis, caret; manage raster, grid data: raster, rgdal, terra;
# derivations of terrain attributes from a digital elevation model: RSAGA; 
# model and goodness of fit measures: randomForest and ithir.


vpk <- paste0(pk, " \\cite{R-", pk, "} ")
# gsub("\\", "", vpk, fixed= TRUE)
dput(vpk)


