rm(list=ls())

# The idea of this script is to generate a reasonable distribution model for each
#   species of Coenonympha. This will allow us to get an idea of their geographic
#   extent, and visually represent their range/habitat overlaps.

#------------------------#
#### Loading Packages ####
#------------------------#

#install the pacman package manager if necessary
if(!requireNamespace("pacman", quietly = TRUE)){
  install.packages("pacman")
}


pacman::p_load(here, # to provide document paths relative to the R project
               dplyr,
               spThin,
               rgeos,
               sp,
               ggplot2,
               maptools,
               raster,
               ecospat,
               nicheROVER,
               RColorBrewer,
               PerformanceAnalytics
)

