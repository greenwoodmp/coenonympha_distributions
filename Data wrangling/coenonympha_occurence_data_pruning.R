rm(list=ls())

setwd('C:/Users/Matt/Documents/Coenonympha_Data/Smaller Occurence Data')


# from http://rstudio-pubs-static.s3.amazonaws.com/415459_61fae98ff3984241bd5f7f31da425c98.html#download-occurrence-data-from-gbif

#######################################
## DOWNLOAD AND CLEAN DATA FROM GBIF ##
#######################################

library(rgbif)
library(scrubr)
library(maps)
library(data.table)
library(dplyr)
library(rgeos)
library(HH)
library(tidyverse)
library(scales)
library(spThin)

# When doing this properly, do this for each Coenonympha species and combine final result
data <- as.data.frame(occ_data(scientificName = "Coenonympha pamphilus", limit = 1000)[[2]])

## Removing missing data rows, and restricting the dataset to the Palearctic
data <- data[!is.na(data$decimalLongitude) & !is.na(data$decimalLatitude)] #remove entries with NA coordinates
data <- data[!which(data$species=='')] #remove entries without species-level ID
data <- data[data$decimalLatitude > 0 & data$decimalLongitude > -30] #rough European co-ordinates

## Removing duplicate rows

recs.dups <- duplicated(data %>% dplyr::select(decimalLongitude, decimalLatitude))
occurrences <- occurrences[!recs.dups, ]

## Applying spatial thinning


# map the occurrence data:
#map("world", xlim = range(data$decimalLongitude), ylim = range(data$decimalLatitude))  # if the map doesn't appear right at first, run this command again
#points(data[ , c("decimalLongitude", "decimalLatitude")], pch = ".")
# you may notice (especially if you zoom in, e.g. by specifying a smaller range of coordinates under 'xlim' and 'ylim' above) that many points are too regularly spaced to be exact locations of species sightings; rather, such points are likely to be centroids of (relatively large) grid cells on which particular surveys was based, so remember to adjust the spatial resolution of your analysis accordingly!
# also, these data are likely to contain species absences and location errors, so jump to "CLEAN THE DATASET" section below - this is VERY IMPORTANT!!!


# CLEAN THE DATASET! ----
# mind that data often contain errors, so careful inspection and cleaning are necessary! 
# here we'll first remove records of absence or zero-abundance (if any):
#names(myspecies_coords)
#sort(unique(data$individualCount))  # notice if some points correspond to zero abundance
#sort(unique(data$occurrenceStatus))  # check for different indications of "absent", which could be in different languages! and remember that R is case-sensitive
absence_rows <- which(data$individualCount == 0 | data$occurrenceStatus %in% c("absent", "Absent", "ABSENT", "ausente", "Ausente", "AUSENTE"))
length(absence_rows)
if (length(absence_rows) > 0) {
  data <- data[-absence_rows, ]
}

# let's do some further data cleaning with functions of the 'scrubr' package (but note this cleaning is not exhaustive!)
#nrow(data)
data <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(data))))
#nrow(data)

# map the cleaned occurrence data:
#map("world", xlim = range(data$decimalLongitude), ylim = range(data$decimalLatitude))  # if the map doesn't appear right at first, run this command again
#points(data[ , c("decimalLongitude", "decimalLatitude")], col = as.factor(data$species), pch = ".")
# possible erroneous points e.g. on the Equator (lat and lon = 0) should have disappeared now
# also eliminate presences with reported coordinate uncertainty (location error, spatial resolution) larger than 5 km (5000 m):
data <- coord_uncertain(data, coorduncertainityLimit = 5000)
#nrow(data)
# but note that this will only get rid of records where coordinate uncertainty is adequately reported, which may not always be the case! Careful mapping and visual inspection is necessary
# map the cleaned occurrence records with a different colour on top of the raw ones:
#points(data[ , c("decimalLongitude", "decimalLatitude")], pch = 20, cex = 0.5, col = "turquoise")

final_occurrence_data <- data

library(dplyr)

final_occurrence_data <- final_occurrence_data %>% 
  arrange(species)

# Remove outliers? Unsure if this is correct
final_occurrence_data <- final_occurrence_data %>%
  group_by(species) %>%
  filter(between(decimalLatitude, mean(decimalLatitude, na.rm=TRUE) - (2.5 * sd(decimalLatitude, na.rm=TRUE)), 
                 mean(decimalLatitude, na.rm=TRUE) + (2.5 * sd(decimalLatitude, na.rm=TRUE)))) %>%
  filter(between(decimalLongitude, mean(decimalLongitude, na.rm=TRUE) - (2.5 * sd(decimalLongitude, na.rm=TRUE)), 
                 mean(decimalLongitude, na.rm=TRUE) + (2.5 * sd(decimalLongitude, na.rm=TRUE))))

# Removing dense geographic data which may be an indicator of bias - using a avg distance of 5 km between samples - arbitrary 
# see https://rmacroecology.netlify.app/2019/01/21/niche-overlap-update-to-silva-et-al-2014-supplementary-matterial/

write.csv2(final_occurrence_data, "Final occurrence data.csv")

