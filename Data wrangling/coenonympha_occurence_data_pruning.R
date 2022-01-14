rm(list=ls())

# this data cleaning process largely follows http://rstudio-pubs-static.s3.amazonaws.com/415459_61fae98ff3984241bd5f7f31da425c98.html#download-occurrence-data-from-gbif

# As data from collaborators has amassed, this document has been edited to 
#   collate data from multiple sources. All additional data are manipulated to
#   conform to GBIF naming conventions before assembling a large occurrence
#   dataset. Final data filtering (e.g. removal by coordinate uncertainty) is
#   performed on the final dataset.

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
               spThin,
               dplyr)


#------------------------------------------------#
#### Downloading and Preparing Data from GBIF ####
#------------------------------------------------#

# GBIF data (search term "Coenonympha") downloaded at https://doi.org/10.15468/dl.xprfux; here, the unzipped Darwin Core Archive file is placed straight into the Raw data folder and used

gbif_data <- fread(here("./Data wrangling/Raw data/0285527-200613084148143/occurrence.txt"), #look at fread() as it is faster than read.csv()
                header=TRUE,
                sep="\t",
                select=c('species',
                         #'individualCount',
                         'occurrenceStatus',
                         'decimalLatitude',
                         'decimalLongitude'
                         ,"coordinateUncertaintyInMeters",
                         "institutionCode")) #change retained columns if necessary

gbif_data <- gbif_data %>% 
  dplyr::filter(!is.na(decimalLongitude), # only retain rows with provided longitude coordinates
         !is.na(decimalLatitude), # only retain rows with provided latitude coordinates
         species!='', # remove rows with no species-level taxonomic assignments
         !occurrenceStatus %in% c("absent", "Absent", "ABSENT", "ausente", "Ausente", "AUSENTE"), # remove any entries with absent notes
         ) %>%
  dplyr::select(species, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters)

# It should be noted that Coenonympha darwiniana is not present in the gbif dataset 
#   directly as it is synonomised with Coenonympha macromma
# As these species are highly localized to Swiss Alps (C. darwiniana) and French
#   Alps (C. macromma), a latitude cut-off will be used to recover C. darwiniana 
#   occurrences from the total list of C. macromma observations. 
gbif_data <- gbif_data %>% 
  mutate(species = case_when(species=="Coenonympha macromma" & decimalLatitude >= 45.5 ~ "Coenonympha darwiniana",
                             TRUE ~ species))

# Check the darwiana change on a map
gbif_test <- gbif_data[gbif_data$species %in% c("Coenonympha macromma","Coenonympha darwiniana")]
map("world", xlim = range(gbif_test$decimalLongitude), ylim = range(gbif_test$decimalLatitude))  # if the map doesn't appear right at first, run this command again
points(gbif_test[ , c("decimalLongitude", "decimalLatitude")], pch = ".", cex=3, col=as.factor(gbif_test$species))
abline(h=45.5) # plot the latitude cut-off line

rm(gbif_test)


#-------------------------------------------#
#### Preparing Data from CSCF Info Fauna ####
#-------------------------------------------#

cscf_data <- read.csv(here("./Data wrangling/Raw data/Collaborator occurrence data/cropped 21_0078 LDesprés Uni Grenoble Coenonympha 20210217.csv"))

cscf_data <- cscf_data %>%
  select('Taxon','LAT','LON','Rayon') %>%
  rename('species'='Taxon', # change variable names to match GBIF format
         'decimalLatitude'='LAT',
         'decimalLongitude'='LON',
         'coordinateUncertaintyInMeters'='Rayon')

# In this file, the coordinate uncertainty is expressed as a range. We'll keep 
#   only the maximum uncertainty
cscf_data$coordinateUncertaintyInMeters <- as.numeric(gsub(".*- ","",cscf_data$coordinateUncertaintyInMeters))

# Remove the gardetta species entries with unclear taxonomy
cscf_data <- subset(cscf_data, species!="Coenonympha gardetta aggr.")         
           

#--------------------------------------#
#### Preparing Data from PCORRADINI ####
#--------------------------------------#

# Find out who this collaborator is

pcor_data <- read.csv(here("./Data wrangling/Raw data/Collaborator occurrence data/Data_Coll_Coenonympha_200821_PCORRADINI.txt"), sep='\t') #look at fread() as it is faster than read.csv()
                   
pcor_data <- pcor_data %>%
  select('Lien_Taxons..c_Binome_latin','Lien_Stations..Coord_Lat','Lien_Stations..Coord_Long') %>%
  rename('species'='Lien_Taxons..c_Binome_latin', # change variable names to match GBIF format
         'decimalLatitude'='Lien_Stations..Coord_Lat',
         'decimalLongitude'='Lien_Stations..Coord_Long')

# First, ensure that the location data uses a period as the decimal marker, and then convert the data to a numeric format
pcor_data$decimalLatitude <- as.numeric(gsub(",",".", pcor_data$decimalLatitude))
pcor_data$decimalLongitude <- as.numeric(gsub(",",".", pcor_data$decimalLongitude))


#-------------------------------------------------#
#### Collating Data with Uncertainty Filtering ####
#-------------------------------------------------#

c_data <- rbind(gbif_data, cscf_data)
  
  
c_data <- c_data %>% 
  dplyr::filter(!is.na(decimalLongitude), # only retain rows with provided longitude coordinates
                !is.na(decimalLatitude), # only retain rows with provided latitude coordinates
                species!='', # remove rows with no species-level taxonomic assignment
                coordinateUncertaintyInMeters < 5000 # set maximum coordinate uncertainty for an arbitrary but reasonable grid value like 5km
                #decimalLatitude > 0, decimalLongitude > -30 # can be used to constrain coordinates to, e.g., the Palearctic
  )

# Count the number of occurrences left per species
data_count <- c_data %>% 
  count(species)

# Clean the data further with functions of the 'scrubr' package
#c_data <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(c_data))))

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