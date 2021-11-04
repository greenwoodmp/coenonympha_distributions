rm(list=ls())

setwd('C:/Users/Matt/Documents/Coenonympha_Data/Results')


# from https://rspatial.org/raster/sdm/4_sdm_envdata.html

library(rgdal)
library(dismo)
library(raster)
library(sf)
library(dplyr)
library(spThin)
library(tidyverse)


#path to the folder containing the bioclim variables
pathb <- file.path('C:/Users/Matt/Documents/Coenonympha_Data/BioClim/2.5m/bioclim')
biofiles <- list.files(pathb, pattern='bil$', full.names=TRUE)
biofiles

#path to the folder containing the elevation data
pathe <- file.path('C:/Users/Matt/Documents/Coenonympha_Data/BioClim/2.5m/elev')
elefiles <- list.files(pathe, pattern='tif$', full.names=TRUE)
elefiles


predictorsb <- stack(biofiles)
predictorse <- stack(elefiles)

# Check alignment of predictor maps - keep in mind occurence data matched wrld simpl perfectly
#library(maptools)
#data(wrld_simpl)
#plot(predictorsb,1)
#plot(wrld_simpl, add=TRUE)
#points(occur[ , c("decimalLongitude", "decimalLatitude")], pch = ".")
#reading in a WWF terrestrial biomes classification file
# following https://stackoverflow.com/questions/35096133/converting-shapefile-to-raster
# and https://gis.stackexchange.com/questions/376280/converting-a-component-of-a-shapefile-into-a-geotiff-file-in-r

#library(fasterize)
#library(dplyr)
#library(sf)
#library(tidyverse)

# Make a BIOM tif file first
#biomefiles <- st_read('C:/Users/Matt/Documents/Coenonympha_Data/BioClim/official_teow/official/wwf_terr_ecos.shp', stringsAsFactors = F)
# the raster cell size is 2.5m for the other data - therefore resolution is 2.5/60
#raster(resolution = 2.5/60, crs = "+proj=longlat +datum=WGS84 +no_defs", vals = 0) %>% 
#  fasterize(biomefiles, ., field = "BIOME") %>% 
#  writeRaster(., "C:/Users/Matt/Documents/Coenonympha_Data/BioClim/official_teow/official/biomes.tif")

pathbiome <- file.path('C:/Users/Matt/Documents/Coenonympha_Data/BioClim/official_teow/official')
biomefiles <- list.files(pathbiome, pattern='tif$', full.names=TRUE)
biomefiles

predictorsbiome <- raster(biomefiles)

#remove the 98 (lakes) and 99 (rock/ice) assignments
predictorsbiome[predictorsbiome==98] <- NA
predictorsbiome[predictorsbiome==99] <- NA

# Read in the occurrence data 
occur <- read.csv('final occurrence data.csv', sep=';', dec=',') #need to let R know that commas are the decimals in the original dataset - will convert to periods
occur <- occur[,-1] #a dummy column of integers was added for some reason so I remove it

occur_points <- occur[c(5,4)]

extracted_valuesb <- extract(predictorsb, occur_points)
#the extracted bios correspond to annual temp (1), annual precipitation (12) and seasonality in temp (4) and precipitation (15)
extracted_valuese <- extract(predictorse, occur_points)
extracted_valuesbiome <- extract(predictorsbiome, occur_points)

#looking at points on rasters
samp <- sample(nrow(occur), 50000)
plot(predictorsbiome)
points(occur_points[samp,], pch='.')

final_dataset <- data.frame(occur$species, occur_points, extracted_valuesb, extracted_valuese, extracted_valuesbiome)
colnames(final_dataset) <- c('Species','decimalLongitude','decimalLatitude',
                             'MeanTemp','MeanTempWarmQuart','MeanTempColdQuart',
                             'MeanPrecip','PrecipWet','PrecipDry','PrecipSeasonality',
                             'PrecipWetQuart','PrecipDryQuart','PrecipWarmQuart','PrecipColdQuart',
                             'MeanDiurnal','Isothermality','TempSeasonality','MaxTempWarm','MaxTempCold',
                             'TempAnnualRange','MeanTempWetQuart','MeanTempDryQuart',
                             'Elevation','Biome')



# Filter out locations with outlier environmental variables/values
# Maybe filter for the other variables now?
final_dataset <- final_dataset %>%
  drop_na() %>% #there are some NA assignments from the rasters to remove
  group_by(Species) %>%
  filter(between(MeanTemp, mean(MeanTemp, na.rm=TRUE) - (2.5 * sd(MeanTemp, na.rm=TRUE)), 
                 mean(MeanTemp, na.rm=TRUE) + (2.5 * sd(MeanTemp, na.rm=TRUE)))) %>%
  filter(between(TempSeasonality, mean(TempSeasonality, na.rm=TRUE) - (2.5 * sd(TempSeasonality, na.rm=TRUE)), 
                 mean(TempSeasonality, na.rm=TRUE) + (2.5 * sd(TempSeasonality, na.rm=TRUE)))) %>%
  filter(between(MeanPrecip, mean(MeanPrecip, na.rm=TRUE) - (2.5 * sd(MeanPrecip, na.rm=TRUE)), 
                 mean(MeanPrecip, na.rm=TRUE) + (2.5 * sd(MeanPrecip, na.rm=TRUE)))) %>%
  filter(between(PrecipSeasonality, mean(PrecipSeasonality, na.rm=TRUE) - (2.5 * sd(PrecipSeasonality, na.rm=TRUE)), 
                 mean(PrecipSeasonality, na.rm=TRUE) + (2.5 * sd(PrecipSeasonality, na.rm=TRUE)))) %>%
  filter(between(Elevation, mean(Elevation, na.rm=TRUE) - (2.5 * sd(Elevation, na.rm=TRUE)), 
                 mean(Elevation, na.rm=TRUE) + (2.5 * sd(Elevation, na.rm=TRUE))))


write.csv(final_dataset, 'Final dataset.csv')

