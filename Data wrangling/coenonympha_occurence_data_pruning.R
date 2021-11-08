rm(list=ls())

# this cleaning process largely follows http://rstudio-pubs-static.s3.amazonaws.com/415459_61fae98ff3984241bd5f7f31da425c98.html#download-occurrence-data-from-gbif

######################
## LOADING PACKAGES ##
######################

#install the pacman package manager if necessary
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")

# this was necessary for the scrubr package to install as it is incompatible with 4.0 and above
#devtools::install_version("scrubr", version = "0.4.0", repos = "http://cran.us.r-project.org")

pacman::p_load(here, # to provide document paths relative to the R project
               data.table, # for the fread function
               rgbif,
               maps,
               scrubr,
               spThin)

#######################################
## DOWNLOAD AND CLEAN DATA FROM GBIF ##
#######################################

# GBIF data (search term "Coenonympha") downloaded at https://doi.org/10.15468/dl.xprfux; here, the unzipped Darwin Core Archive file is placed straight into the Raw data folder and used

c_data <- fread(here("./Data wrangling/Raw data/0285527-200613084148143/occurrence.txt"), #look at fread() as it is faster than read.csv()
                header=TRUE,
                sep="\t",
                select=c('species',
                         'individualCount',
                         'occurrenceStatus',
                         'decimalLatitude',
                         'decimalLongitude'
                         ,"coordinateUncertaintyInMeters",
                         "institutionCode")) #change retained columns if necessary

c_data <- c_data %>% 
  dplyr::filter(!is.na(decimalLongitude), # only retain rows with provided longitude coordinates
         !is.na(decimalLatitude), # only retain rows with provided latitude coordinates
         species!='', # remove rows with no species-level taxonomic assignments
         !occurrenceStatus %in% c("absent", "Absent", "ABSENT", "ausente", "Ausente", "AUSENTE"), # remove any entries with absent notes
         coordinateUncertaintyInMeters <= (5000) # set maximum coordinate uncertainty for an arbitrary but reasonable grid value like 5km
         #decimalLatitude > 0, decimalLongitude > -30 # can be used to constrain coordinates to, e.g., the Palearctic
         )

# Clean the data further with functions of the 'scrubr' package
c_data <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(c_data))))

# Thin the data points to reduce geographic sampling bias. Here, the previously selected aribtrary uncertainty theshold (5km) is used
c_data.thin <- thin(c_data, 
                        verbose = FALSE, 
                        lat.col = "decimalLatitude",
                        long.col = "decimalLongitude",
                        spec.col = "species",
                        thin.par = 5,
                        reps = 100, 
                        write.files = FALSE,
                        write.log.file = FALSE,
                        locs.thinned.list.return = TRUE)

final_occurrence_data <- final_occurrence_data %>%
  group_by(species) %>%
  filter(between(decimalLatitude, mean(decimalLatitude, na.rm=TRUE) - (2.5 * sd(decimalLatitude, na.rm=TRUE)), 
                 mean(decimalLatitude, na.rm=TRUE) + (2.5 * sd(decimalLatitude, na.rm=TRUE)))) %>%
  filter(between(decimalLongitude, mean(decimalLongitude, na.rm=TRUE) - (2.5 * sd(decimalLongitude, na.rm=TRUE)), 
                 mean(decimalLongitude, na.rm=TRUE) + (2.5 * sd(decimalLongitude, na.rm=TRUE))))

# Removing dense geographic data which may be an indicator of bias - using a avg distance of 5 km between samples - arbitrary 
# see https://rmacroecology.netlify.app/2019/01/21/niche-overlap-update-to-silva-et-al-2014-supplementary-matterial/


##############################
## QUICK DATA VISUALIZATION ##
##############################

map("world", xlim = range(c_data$decimalLongitude), ylim = range(c_data$decimalLatitude))  # if the map doesn't appear right at first, run this command again
points(c_data[ , c("decimalLongitude", "decimalLatitude")], pch = ".")


#########################
## EXPORT CLEANED DATA ##
#########################

c_data <- c_data %>% 
  dplyr::arrange(species)

write.csv2(c_data, here("./Data wrangling/Processed data/coenonympha_occurrence_data.csv"))
