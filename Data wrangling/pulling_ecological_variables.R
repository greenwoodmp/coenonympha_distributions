rm(list=ls())

# This data preparation process largely follows:
#   https://rspatial.org/raster/sdm/4_sdm_envdata.html

#------------------------#
#### Loading Packages ####
#------------------------#

#install the pacman package manager if necessary
if(!requireNamespace("pacman", quietly = TRUE)){
  install.packages("pacman")
  }


pacman::p_load(here, # to provide document paths relative to the R project
               rgdal, # for the fread function
               dismo, # for the gridSample function
               raster,
               sf,
               dplyr,
               spThin,
               tidyverse,
               sf, # for reading in shape files - st_read()
               fasterize,
               ggplot2 # for graph plotting
               )

#---------------------------------#
#### Reading in Worldclim Data ####
#---------------------------------#

climate <- getData('worldclim',
                   var='bio',
                   path=here("./Data wrangling/Raw data/Bioclim/"),
                   res=2.5
                   )

altitude <- getData('worldclim',
                    var='alt',
                    path=here("./Data wrangling/Raw data/Bioclim/"),
                    res=2.5
                    )


predictorsb <- stack(biofiles)
predictorse <- stack(elefiles)


#------------------------------#
#### Reading in Biome Data  ####
#------------------------------#

# Unfortunately, adding a raster for categorical biome data is not as simple as
#   as adding the worldcilm rasters. The biome file I am using can be retrieved
#   from https://ecoregions.appspot.com by downloading the Ecoregions2017.zip.
#   The downloaded file should be unzipped and placed into the
#   "/Data wrangling/Raw data/" folder

# The unzipped file we need is a shape file (.shp), which is converted to a 
#   raster following guidelines stipulated in the tutorial below:
#   https://gis.stackexchange.com/questions/376280/converting-a-component-of-a-shapefile-into-a-geotiff-file-in-r

# The below check is used to avoid doing this conversion every time the script is run
#   Essentially, we look for the name of the raster to be produced (Ecoregions2017.tif)
#   and, if it already exists, skip the raster production steps

if(file.exists(here("./Data wrangling/Raw data/Ecoregions2017/Ecoregions2017.tif")) == FALSE){
  
  # Read in the shape file
  biomefile <- st_read(here("./Data wrangling/Raw data/Ecoregions2017/Ecoregions2017.shp"), stringsAsFactors = F)
  
  # Convert the shape file to a raster
  # Note that the raster cell size is 2.5 arc minutes for the other raster data.
  #   Therefore, the required resolution (in degrees) is (2.5 arc mins)/(60 arc mins)
  raster(resolution = 2.5/60, 
         crs = "+proj=longlat +datum=WGS84 +no_defs", 
         vals = 0) %>%
    fasterize(biomefile, ., field = "BIOME_NUM") %>%
    writeRaster(., here("./Data wrangling/Raw data/Ecoregions2017/Ecoregions2017.tif"))
  
  rm(biomefile) # clean the biome shape file from memory
  
}else{
  message("The converted raster file for the biome shape file already exists")
}

# Load in the prepared biome raster data
biome <- raster(here("./Data wrangling/Raw data/Ecoregions2017/Ecoregions2017.tif"))


#----------------------------------#
#### Checking Raster Alignments ####
#----------------------------------#

# Check that the resolution of each map matches
res(climate)
res(altitude)
res(biome)

# Check that the extent of each map starts and ends in compatible locations
extent(climate)
extent(altitude)
extent(biome) 
# note that the ymin of biome is lower than the other datasets - the extents 
#   should otherwise match perfectly

# Alignments can also be checked by graphing the rasters over eachother with ggplot2
# Just a small section of the rasters are plotted to inspect cell overlap, which
#   is more challenging to observe when the entire raster is plotted
plot(climate[[1]],
     xlim=range(-5:10), # The x and y lims used are arbitrary
     ylim=range(42:51),
     alpha=0.3
     )
plot(biome, 
     alpha=0.3,
     xlim=range(-5:10),
     ylim=range(42:51),
     add=TRUE # ensures that this raster is placed above the first one
     )

# For the most part, the alignment is okay. The biome raster appears to generally
#   fall a little short of the edges of the climate and altitude raster data
#   (e.g., at xlim=range(-55:-50); ylim=range(48:51) or xlim=range(-5:10); ylim=range(48:51) )
#   but the populated cells always overlap
# It may be useful to look into cropping the original biome shapefile by the
#   worldclim rasters for added accuracy (see, as a start, https://stackoverflow.com/questions/23073669/clipping-raster-using-shapefile-in-r-but-keeping-the-geometry-of-the-shapefile)


#--------------------------------------------------------------#
#### Assigning Ecological Raster Data to Occurrence Records ####
#--------------------------------------------------------------#

# Read in the occurrence data 
# Need to indicate that commas are the decimals in this dataset - these will be converted to periods
occur <- read.csv(here("./Data wrangling/Processed data/coenonympha_occurrence_data.csv"),
                  sep=';',
                  dec=','
                  ) 

# Pull only the latitude and longitude data from the file
occur_points <- occur[c("decimalLongitude", "decimalLatitude")]

extracted_climate <- raster::extract(climate, occur_points, cellnumbers=TRUE)
extracted_altitude <- raster::extract(altitude, occur_points, cellnumbers=TRUE)
extracted_biome <- raster::extract(biome, occur_points, cellnumbers=TRUE)

# Test to see if the raster cells IDed in each case of extraction are the same
# If they are identical, the datasets can be simply appended to eachother
identical(extracted_climate[,1], extracted_altitude[,1])
identical(extracted_biome[,1], extracted_climate[,1])
identical(extracted_biome[,1], extracted_altitude[,1])

# Generate a dataframe containing species names and all the climate/ecology rasters
#   mapped to the occurrence data.
# The raster cell calls are retained for a point filtering step coming up next
dataset <- data.frame(occur$species, 
                      occur_points, 
                      extracted_climate, 
                      extracted_altitude[,-1], # The raster cell ID is removed as it is already captured in extracted_climate
                      extracted_biome[,-1] # The raster cell ID is removed as it is already captured in extracted_climate
                      )


#-------------------------------------------------#
#### Subsampling Occurrence Data Using Rasters ####
#-------------------------------------------------#

# There is almost certainly sampling bias in the occurrence dataset, with
#   repeated ID of the same species in some regions (e.g., in close proximity 
#   to urban regions) and sparse sampling in others (e.g. uninhabited locales)
# To reduce the impact of this "geographic sampling bias" on the occurrence data,
#   the occurrence points are resampled to ensure that, for each taxon, presence
#   within an environmental raster cell is only counted once.
# To do this, the raster cell number calls will be used - for each raster cell
#   per species, only a single occurrence point will be retained, preventing
#   the overrepresentation of a species within a close sampling region

subsampled_dataset <- dataset %>%
  group_by(occur.species, cells) %>%
  drop_na() %>% # Remove any rows for which any raster has an NA before subsampling
  slice(1) # Slice will select only 1 occurrence entry per raster location per species

# Confirm raster subsampling on plot
plot(altitude[[1]],
     xlim=range(-10:10), # The x and y lims used are arbitrary
     ylim=range(30:45),
     alpha=0.3
     )
points(x=as.numeric(unlist(subsampled_dataset["decimalLongitude"])), # the numeric values must be unlisted to plot
       y=as.numeric(unlist(subsampled_dataset["decimalLatitude"])), # the numeric values must be unlisted to plot
       pch="."
       )

# Remove the raster ID column from the dataset
subsampled_dataset <- subset(subsampled_dataset, select=-c(cells))

# Look at the remaining occurrence point counts per species
counts <- table(subsampled_dataset$occur.species)


#------------------------------------------------------------------#
#### Preparing the Final Dataset for Spatial and Niche Analysis ####
#------------------------------------------------------------------#

colnames(subsampled_dataset) <- c('Species','decimalLongitude','decimalLatitude',
                             'MeanTemp','MeanTempWarmQuart','MeanTempColdQuart',
                             'MeanPrecip','PrecipWet','PrecipDry','PrecipSeasonality',
                             'PrecipWetQuart','PrecipDryQuart','PrecipWarmQuart','PrecipColdQuart',
                             'MeanDiurnal','Isothermality','TempSeasonality','MaxTempWarm','MaxTempCold',
                             'TempAnnualRange','MeanTempWetQuart','MeanTempDryQuart',
                             'Altitude','Biome')

data.table::fwrite(subsampled_dataset,
          here("./Data wrangling/Processed data/coenonympha_spatial_data.csv"),
          sep=",",
          dec="."
          )

