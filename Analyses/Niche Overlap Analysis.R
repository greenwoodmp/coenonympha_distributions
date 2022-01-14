rm(list=ls())

# This data preparation process largely follows:
#   https://rmacroecology.netlify.app/2019/01/21/niche-overlap-update-to-silva-et-al-2014-supplementary-matterial/


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


#----------------------------------#
#### Filtering the Species List ####
#----------------------------------#

data <- read.csv(here("./Data wrangling/Processed data/coenonympha_spatial_data.csv"))


# List of the species we have genomic data for
species.list <- c("Coenonympha amaryllis",
                  "Coenonympha dorus",
                  "Coenonympha leander",
                  "Coenonympha rhodopensis",
                  "Coenonympha arcania",
                  "Coenonympha hero",
                  "Coenonympha macromma",
                  "Coenonympha oedippus",
                  "Coenonympha saadi",
                  "Coenonympha thyrsis",
                  "Coenonympha arcanioides",
                  "Coenonympha gardetta",
                  #"Coenonympha myops",
                  #"Coenonympha orientalis",
                  #"Coenonympha semenovi",
                  "Coenonympha tullia",
                  "Coenonympha corinna",
                  "Coenonympha glycerion",
                  #"Coenonympha nolckeni",
                  "Coenonympha pamphilus",
                  #"Coenonympha sunbecca",
                  "Coenonympha vaucheri")

# Filter the list for only species that we have genomic data for
data <- data %>% #need C phryne, L. myops, C. elbana, iphoides - where is our orientalis
  filter(Species %in% species.list)


#-----------------------------------#
#### Preparing Data for Analyses ####
#-----------------------------------#

# Define the number of groups to be tested. Here, this will equal the number of
#   species we wish to assess
n.groups <- length(unique(data$Species))

# To perform niche analysis, we need a background environmental spaces for each
#   species
# To get this, a minimum convex polygon (MCP) for each species will be generated
buffer.size <- 1

# A function for defining the MCP, pulled from https://github.com/ndimhypervol/wallace
mcp <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
}

# Getting the ecological variables; these only need to download once
# The 19 Worldclim bioclimatic variables
climate <- getData('worldclim',
                   var='bio',
                   path=here("./Data wrangling/Raw data/Bioclim/"),
                   res=2.5
                   )

# Altitudinal worldclim data
altitude <- getData('worldclim',
                    var='alt',
                    path=here("./Data wrangling/Raw data/Bioclim/"),
                    res=2.5
                    )

variables <- stack(climate, altitude)


#-------------------------#
#### Group Assignments ####
#-------------------------#

# Now we generate MCPs for the data
# Keep in mind, the occurrence records are already thinned by the worldclim rasters
#   in the pulling_ecological_variables.R script, so they are good to go

# Union of the world map
lps <- getSpPPolygonsLabptSlots(wrld_simpl)
IDFourBins <- cut(lps[,1], range(lps[,1]), include.lowest=TRUE)
world <- unionSpatialPolygons(wrld_simpl, IDFourBins)


# Empty objects
g.assign <- numeric(length(data[[1]]))
xy.mcp <- list()
back.env <- list()
spec.env <- list()
row.sp <- list()

# Plot map
data(wrld_simpl)
plot(wrld_simpl)

# Loop
for(i in 1:n.groups) {
  
  # Define groups
  cut1 <- data[[1]][, 1] >= group.long[i]
  cut2 <- data[[1]][, 1] < group.long[i + 1]
  g.limit <- cut1 & cut2
  
  # Save row numbers per species
  row.sp[[i]] <- which(g.limit)
  g.assign[g.limit] <- g.names[i]
  
  # Background polygon
  mcp.occ <- mcp(as.matrix(data[[1]][g.limit, ]))
  xy.mcp.i <- gBuffer(mcp.occ, width = buffer.size)
  proj4string(xy.mcp.i) <- proj4string(world)
  xy.mcp[[i]] <- gIntersection(xy.mcp.i, world, byid=TRUE, drop_lower_td=TRUE)
  # Background environment
  back.env[[i]] <- na.exclude(do.call(rbind.data.frame, extract(variables, xy.mcp[[i]])))
  # Species environment
  spec.env[[i]] <- na.exclude(extract(variables, data[[1]][g.limit, ]))
  
  # Plot
  points(data[[1]][g.limit, ], col = g.colors[i],
         pch = 20, cex = 0.7)
  plot(xy.mcp[[i]], add = TRUE, border = g.colors[i], lwd = 2)
  
}


#-----------------------------------------#
#### Trying a Different ecospat Method ####
#-----------------------------------------#
rm(list=ls())
# see https://rsh249.github.io/spatial_bioinformatics/niche_overlap.html

pacman::p_load(here, # to provide document paths relative to the R project
               spocc,
               raster,
               viridis,
               ENMTools,
               dplyr,
               tibble,
               data.table,
               ecospat
)

data <- read.csv(here("./Data wrangling/Processed data/coenonympha_spatial_data.csv"))

# Get the extent for the species occurrence data, which has already been thinned
#  by the environmental rasters cells in pulling_ecological_variables.py
ext=extent(c(min(data$decimalLongitude)-5, max(data$decimalLongitude)+5, min(data$decimalLatitude)-5, max(data$decimalLatitude)+5))

# Get environmental rasters
# Worldclim bioclimatic variables
climate <- getData('worldclim',
                   var='bio',
                   path=here("./Data wrangling/Raw data/Bioclim/"),
                   res=2.5
)

# Altitudinal worldclim data
altitude <- getData('worldclim',
                    var='alt',
                    path=here("./Data wrangling/Raw data/Bioclim/"),
                    res=2.5
)

variables <- stack(climate, altitude)

Env = crop(variables, ext)


# Writing a function for pairwise species comarisons

# Note that ncores = detectCores() -1 could be a useful parameter to set to speed up the niche equivalence and similarity tests
niche_compare <- function(species1, species2) {
  
  #background by radius
  sp1 <- data[data$Species == species1, c("decimalLongitude","decimalLatitude")]
  sp2 <- data[data$Species == species2, c("decimalLongitude","decimalLatitude")]
  
  # Some datasets are extremely large, and will be subsampled for computational efficiency
  # Setting a cut-off size
  n <- 1000
  set.seed(12345)
  sp1 <- sp1 %>%
    sample_n(if(n() < n) n() else n)
  
  set.seed(12345)
  sp2 <- sp2 %>%
    sample_n(if(n() < n) n() else n)
  
  # Simulate background environmental points
  bgsp1 <- background.points.buffer(sp1, 
                                    radius = 200000, 
                                    n = 10*nrow(sp1), 
                                    mask = Env[[1]])
  bgsp2 <- background.points.buffer(sp2, 
                                    radius = 200000, 
                                    n = 10*nrow(sp2), 
                                    mask = Env[[1]])
  
  # Get environmental data
  extractsp1 <- na.omit(cbind(sp1[,1:2], extract(Env, sp1[,1:2]), rep(1, nrow(sp1))))
  extractsp2 = na.omit(cbind(sp2[,1:2], extract(Env, sp2[,1:2]), rep(1, nrow(sp2))))
  
  colnames(extractsp1)[ncol(extractsp1)] = 'occ'
  colnames(extractsp2)[ncol(extractsp2)] = 'occ'
  
  extbgsp1 = na.omit(cbind(bgsp1, extract(Env, bgsp1), rep(0, nrow(bgsp1))))
  extbgsp2 = na.omit(cbind(bgsp2, extract(Env, bgsp2), rep(0, nrow(bgsp2))))
  
  colnames(extbgsp1)[ncol(extbgsp1)] = 'occ'
  colnames(extbgsp2)[ncol(extbgsp2)] = 'occ'
  
  # merge occ and bg data 
  datsp1 = rbind(extractsp1, extbgsp1)
  datsp2 = rbind(extractsp2, extbgsp2)
  
  
  pca.env <- dudi.pca(
    rbind(datsp1, datsp2)[,3:22],
    scannf=FALSE,
    nf=2
  )
  
  #Variable contribution
  #ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
  
  scores.globclim<-pca.env$li # PCA scores for the whole study area (all points)
  
  scores.sp1 <- suprow(pca.env,
                         extractsp1[which(extractsp1[,23]==1),3:22])$li # PCA scores for the species 1 distribution
  scores.sp2 <- suprow(pca.env,
                         extractsp2[which(extractsp2[,23]==1),3:22])$li # PCA scores for the species 1 distribution
  
  scores.climsp1 <- suprow(pca.env,datsp1[,3:22])$li # PCA scores for the whole native study area
  scores.climsp2 <- suprow(pca.env,datsp2[,3:22])$li # PCA scores for the whole native study area
  
  
  grid.climsp1 <- ecospat.grid.clim.dyn(
    glob = scores.globclim,
    glob1 = scores.climsp1,
    sp = scores.sp1,
    R = 100,
    th.sp = 0
  )
  grid.climsp2 <- ecospat.grid.clim.dyn(
    glob = scores.globclim,
    glob1 = scores.climsp2,
    sp = scores.sp2,
    R = 100,
    th.sp = 0
  )
  
  #D.overlap <- ecospat.niche.overlap (grid.climama, grid.climarc, cor=T)$D 
  
  eq.test <- ecospat.niche.equivalency.test(grid.climsp1, grid.climsp2,
                                            ncores = 5, # change this based on core availability
                                            rep=1000, alternative = "lower") ##rep = 1000 recommended for operational runs
  
  sim.test <- ecospat.niche.similarity.test(grid.climsp1, grid.climsp2,
                                            ncores = 5, # change this based on core availability
                                            rep=1000, alternative = "greater",
                                            rand.type=1) 
  output <- data.frame(eq.test$obs$D, eq.test$p.D, sim.test$p.D)
  colnames(output) <- c("Observed_D","Equivalence_p", "Similarity_p")
  

  return(output)

}

# List of the species we have genomic data for
species.list <- c("Coenonympha amaryllis",
                  "Coenonympha dorus",
                  "Coenonympha leander",
                  "Coenonympha rhodopensis",
                  "Coenonympha arcania",
                  "Coenonympha hero",
                  "Coenonympha macromma",
                  "Coenonympha oedippus",
                  "Coenonympha saadi",
                  "Coenonympha thyrsis",
                  "Coenonympha arcanioides",
                  "Coenonympha gardetta",
                  #"Coenonympha myops",
                  #"Coenonympha orientalis",
                  #"Coenonympha semenovi",
                  "Coenonympha tullia",
                  "Coenonympha corinna",
                  "Coenonympha glycerion",
                  #"Coenonympha nolckeni",
                  "Coenonympha pamphilus",
                  #"Coenonympha sunbecca",
                  "Coenonympha vaucheri")

# Compare all species pairwise using the niche_compare function

# First, make a dataframe containing all pairwise combinations of species names
pairwise <- expand.grid(species.list,species.list)

# Add empty columns that will be filled with niche test data
pairwise.niche <- pairwise %>%
  add_column(Observed_D = NA, Equivalence_p = NA, Similarity_p = NA)

# Run the function iteratively through the dataset, filling columns as it goes
# This will take a LONG time to run - the niche tests are slow, and this is iterating
#   through hundreds of them
# Try and write the output to a file to avoid having to perform the calculations
#   again
for(i in 1:length(pairwise.niche$Var1)){
  holder <- niche_compare(pairwise.niche$Var1[i],pairwise.niche$Var2[i])
  pairwise.niche$Observed_D[i] <- holder$Observed_D
  pairwise.niche$Equivalence_p[i] <- holder$Equivalence_p
  pairwise.niche$Similarity_p[i] <- holder$Similarity_p
  
  message(paste0("Loop ",i," completed; ",round((i/length(pairwise.niche$Var1))*100, 2),"% done."))
}

# Write the niche comparison data to an output file

data.table::fwrite(pairwise.niche,
                   here("./Data wrangling/Processed data/coenonympha_niche_comparisons.csv"),
                   sep=",",
                   dec="."
)















#-------------------------------------#
#### Trying nicheRover Temporarily ####
#-------------------------------------#

# random draws from posterior distribution with default prior
nsamples <- 500

# Define the data to be used

# nicheROVER needs us to use tapply, which can't handle correlation in data sets
cor <- chart.Correlation(data[, 4:23])
sel.dat <- c(4,7,10,16,17,18,19,21,22,23)
data.par <- tapply(1:nrow(data), data$Species,
                   function(ii) niw.post(nsamples = nsamples, X = data[ii,c(5,12)]))

# display p(mu | X) and p(Sigma | X) for each butterfly species
clrs <- colorRampPalette(brewer.pal(11, "Paired"))(length(unique(data$Species))) # colors for each species
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(data.par, col = clrs)
legend(x = "topright", legend = names(data.par), fill = clrs)

# 2-d projections of 10 niche regions
nsamples <- 10
data.par <- tapply(1:nrow(data), data$Species,
                   function(ii) niw.post(nsamples = nsamples, X = data[ii,4:9]))

# format data for plotting function
data.data <- tapply(1:nrow(data), data$Species, function(ii) X = data[ii,4:9])

niche.plot(niche.par = data.par, niche.data = data.data, pfrac = .05,
           col = clrs, xlab = expression("Isotope Ratio (\u2030)"))

# niche overlap plots for 95% niche region sizes

# overlap calculation.  use nsamples = nprob = 1e4 for higher accuracy.
nsamples <- 1e3
over.stat <- overlap(data.par, nreps = nsamples, nprob = nsamples, alpha = .95)

# overlap plot
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise",
             equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")
