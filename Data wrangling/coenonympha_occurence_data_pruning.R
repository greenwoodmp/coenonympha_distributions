rm(list=ls())

# this data cleaning process largely follows http://rstudio-pubs-static.s3.amazonaws.com/415459_61fae98ff3984241bd5f7f31da425c98.html#download-occurrence-data-from-gbif


#------------------------#
#### Loading Packages ####
#------------------------#

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


#------------------------------------------------#
#### Downloading and Preparing Data from GBIF ####
#------------------------------------------------#

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
         coordinateUncertaintyInMeters <= 5000 # set maximum coordinate uncertainty for an arbitrary but reasonable grid value like 5km
         #decimalLatitude > 0, decimalLongitude > -30 # can be used to constrain coordinates to, e.g., the Palearctic
         )

# Clean the data further with functions of the 'scrubr' package
c_data <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(c_data))))

# In a later script, one occurrence point per raster cell will be selected to avoid geographic sampling bias


#---------------------------------#
#### Simple Data Visualization ####
#---------------------------------#

map("world", xlim = range(c_data$decimalLongitude), ylim = range(c_data$decimalLatitude))  # if the map doesn't appear right at first, run this command again
points(c_data[ , c("decimalLongitude", "decimalLatitude")], pch = ".")


#------------------------------#
#### Exporting Cleaned Data ####
#------------------------------#

write.csv2(c_data, 
           here("./Data wrangling/Processed data/coenonympha_occurrence_data.csv"),
           row.names=FALSE
           )