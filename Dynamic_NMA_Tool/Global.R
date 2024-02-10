packages = c("shiny", "ggplot2", "xtable", "forestplot",
              "grDevices", "igraph", "grid", "gridBase", 
              "gridExtra", "DT", "formattable", "htmlwidgets",
              "dplyr", "kableExtra", "scales", "ggthemes", 
             "shinyjs", "RColorBrewer", "markdown", "plotly")


testin = function(package){
  if(!(package %in% installed.packages())){
    install.packages(package)
  }  
}


#sapply(packages, testin)
lapply(packages, require, character.only = TRUE)
library(gridExtra)
library(plotly)
library(igraph)




source("NMA Functions Baseline Risk.R")
