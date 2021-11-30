rm(list=ls())

# This data preparation process largely follows:
#   https://github.com/KarasiewiczStephane/WitOMI


#------------------------#
#### Loading Packages ####
#------------------------#

#install the pacman package manager if necessary
if(!requireNamespace("pacman", quietly = TRUE)){
  install.packages("pacman")
}


pacman::p_load(here, # to provide document paths relative to the R project
               dplyr,
               subniche,
               ecospat,
               nicheROVER,
               ggplot2,
               factoextra,
               NbClust,
               mclust,
               ggforce, # for the covex hull drawing; need to run devtools::install_github("joelgombin/concaveman") for this package to actually work
               concaveman
               )

data <- read.csv(here("./Data wrangling/Processed data/coenonympha_spatial_data.csv"))


#----------------------------------#
#### Filtering the Species List ####
#----------------------------------#

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

# Reduce dataframe size for testing purposes
#n <- 100 
#set.seed(123)
#data <- data %>% 
#  group_by(Species) %>%
  #filter(n() >=100) %>% #remove species with less than 100 entries
#  sample_n(if(n() < 100) n() else n) # see https://stackoverflow.com/questions/55676368/custom-grouped-dplyr-function-sample-n/55676429#55676429
                                        # Here, I subsample if a group has more than 100 entries

#data <- filter(data, !grepl('Coenonympha sunbecca', Species)) #remove it because it's skewing things
data$Species <- factor(data$Species) #setting factor levels as dropped species were still in factor list
## Make a dataframe for species presence by site ########################

# It's important to realise every coordinate is a unique site
# We therefore need to think differently about assigning presence/absence values
# First, we can make the skeleton of the new dataframe
# To do this, we have columns for each species
# In the original dataset (data), each unique geographic location is associated with a species
# We'll exploit this propery to get a p/a matrix

# Replicate the species names and place them under the column names for each different species
data2 <- data.frame(replicate(length(unique(data$Species)), data[,1]))
colnames(data2) <- unique(data$Species)
data2[] <- lapply(data2, as.character) #get rid of the pesky introduced factor-ship

# Now, every time a species in the rows matches the designated species in the column, count it as a presence (i.e, 1)
for(i in unique(data$Species)){
  data2[,i][data2[,i]==i] <- 1
} 

# Now, without using a time-costly loop, just change any non-1 values to zero as they will not match the designated species in the column
data2[data2 != 1] <- 0

## Perform OMI ##########################################################

#env <- data[,-c(1:3,24)] #This is an environmental dataset with the actual coordinates and species names trimmed out - also removed seasonality as it skewed a lot
#the above has all the environmental variables - below we select only a few
env <- data[,-c(1:3,24)] # remove non_environmental variables and the biome variable
spe <- data.frame(lapply(data2, as.numeric)) # This is the adjusted 'site-specific' dataset we created above

# Note, the way that data2 was made should leave the sites matched up between data and data2
omi <- list(spe=spe, env=env) # an omi list containing these dataframes


dudi1 <- dudi.pca(omi$env, 
                  scale = TRUE, 
                  scan = FALSE, 
                  nf = length(colnames(env))) # retain as many axes as there are environmental variables
#scatter(dudi1)

nic <- niche(dudi1, 
             omi$spe, 
             scannf=FALSE,
             nf = length(colnames(env))) # retain as many axes as there are environmental variables

# Get the relative contribution of each axis to the total variation
var.exp <- (nic$eig/sum(nic$eig))*100

# Retrieve a dataset containing the location of each habitat point 
# niche.points <- data.frame(niche.x=nic$ls[,1], niche.y =  nic$ls[,2])
niche.points <- data.frame(nic$ls)

# Retrieve a dataset containing the location of each species point, which is
#   essentially the central coordinate for all habitat points for a species across all axes
# species.points <- data.frame(spe.x=nic$li[,1], spe.y =  nic$li[,2]) 
species.points <- data.frame(nic$li)

# Make an additional dataframe to add "spider points" for each species, i.e., drawing
#   lines from the PCA point of each species to the species centroid to visualise
#   niche dispersion
# First, make a dummy dataframe where the pca coordinates are assigned to their 
#   representative species
dummy <- data.frame(niche.points, Species=as.character(data$Species))
dummy$Species <- gsub(" ", ".", dummy$Species) # Replace the spaces with a period to match the species.points dataset row names

# Make the species.points rownames a data column
species.points.dummy <- data.frame(species.points, Species=rownames(species.points))
colnames(species.points.dummy) <- gsub("Axis", "CentreAxis", colnames(species.points.dummy)) # Rename the columns so they are distinguishable from the pca data axes titles

# Merge the two datasets by the Species column so that each species occurrence 
#   point is matched to the species centroid value
spider <- merge(dummy, species.points.dummy, by="Species")

# Clean the dummy sets from memory
rm(dummy, species.points.dummy)


# Plot the geographic data with the species niche centroid positions
                                                            # bring in the biome variable to colour position points
ggplot(niche.points, aes(x=CS1, y=CS2, colour=as.factor(data$Species)))+
  theme(
    panel.background=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_rect(colour='black', fill=NA, size=0.5),
    #legend.title=element_blank(),
    legend.position='none',
    legend.key=element_rect(fill='white'),
    axis.title=element_text(size=14),
    axis.text=element_text(size=12),
    legend.title=element_text(size=12)
  )+
  geom_point(alpha=0.3)+
  stat_ellipse(aes(x=CS1, y=CS2, group=as.factor(data$Species)), level=0.95, col="black")+
  #geom_segment(aes(x=spider$CS1,
  #                 xend=spider$CentreAxis1,
  #                 y=spider$CS2,
  #                 yend=spider$CentreAxis2)) # Drawing the spider lines - linking species points to centroids 
  #geom_mark_hull(expand=0.01,aes(fill=data$Species))+
  geom_point(data=species.points, aes(x=Axis1,y=Axis2), col='black', size=0.1)+
  geom_text(data=species.points, aes(x=Axis1,y=Axis2, label=colnames(data2)), col="black")+
  xlab(paste("PC1 (", round(var.exp[1], digits = 3), " % explained)", sep=""))+
  ylab(paste("PC2 (", round(var.exp[2], digits = 3), " % explained)", sep=""))


ggplot(niche.points, aes(x=CS1, y=CS3, colour=as.factor(data$Species)))+
  theme(legend.position = 'none',
        panel.background = element_blank())+
  geom_point(alpha=0.3)+
  stat_ellipse(aes(x=CS1, y=CS3, group=as.factor(data$Species)), level=0.95, col="black")+
  #geom_segment(aes(x=spider$CS1,
  #                 xend=spider$CentreAxis1,
  #                 y=spider$CS2,
  #                 yend=spider$CentreAxis2)) # Drawing the spider lines - linking species points to centroids 
  #geom_mark_hull(expand=0.01,aes(fill=data$Species))+
  geom_point(data=species.points, aes(x=Axis1,y=Axis3), col='black', size=0.1)+
  geom_text(data=species.points, aes(x=Axis1,y=Axis3, label=colnames(data2)), col="black")+
  xlab(paste("PC1 (", round(var.exp[1], digits = 3), " % explained)", sep=""))+
  ylab(paste("PC3 (", round(var.exp[3], digits = 3), " % explained)", sep=""))


ggplot(niche.points, aes(x=CS2, y=CS3, colour=as.factor(data$Species)))+
  theme(legend.position = 'none',
        panel.background = element_blank())+
  geom_point(alpha=0.3)+
  stat_ellipse(aes(x=CS2, y=CS3, group=as.factor(data$Species)), level=0.95, col="black")+
  #geom_segment(aes(x=spider$CS1,
  #                 xend=spider$CentreAxis1,
  #                 y=spider$CS2,
  #                 yend=spider$CentreAxis2)) # Drawing the spider lines - linking species points to centroids 
  #geom_mark_hull(expand=0.01,aes(fill=data$Species))+
  geom_point(data=species.points, aes(x=Axis2,y=Axis3), col='black', size=0.1)+
  geom_text(data=species.points, aes(x=Axis2,y=Axis3, label=colnames(data2)), col="black")+
  xlab(paste("PC2 (", round(var.exp[2], digits = 3), " % explained)", sep=""))+
  ylab(paste("PC3 (", round(var.exp[3], digits = 3), " % explained)", sep=""))
 

#---------------------------------------#
#### Interspecific Niche Comparisons ####
#---------------------------------------#

# This largely follows the tutorial below:
#   https://plantarum.ca/2021/07/29/ecospat/

# We want to be able to statistically test the hypothesis that different species
#   are found in different environmental niches.
# While this provides a simple overview of niche diversification, it also serves
#   as a starting block for assessing factors like niche conservatism within the 
#   genus, which is something interesting to consider in an ecological speciation
#   framework

# First, make a grid object for all species being compared

#We need a partitioned dataframe of occurrence and environmental variables for 
#   each species
global.scores <- dudi1$li[,1:2]

# C. amaryllis
ama <- filter(data, Species=="Coenonympha amaryllis")
amaLS <- ama[,c("decimalLongitude","decimalLatitude")]
amaEnv <- ama[,!(names(ama) %in% c("Species","decimalLongitude","decimalLatitude", "Biome"))]

amaLS.scores <- suprow(dudi1, data.frame(ama)[, colnames(env)])$li[,1:2]
amaEnv.scores <- suprow(dudi1, amaEnv)$li[,1:2]

ama.grid <- ecospat.grid.clim.dyn(global.scores,
                                  amaEnv.scores,
                                  amaLS.scores)

rm(ama, amaLS, amaEnv, amaLS.scores, amaEnv.scores) 

# C. pamphilus
pam <- filter(data, Species=="Coenonympha pamphilus")
pamLS <- pam[,c("decimalLongitude","decimalLatitude")]
pamEnv <- pam[,!(names(pam) %in% c("Species","decimalLongitude","decimalLatitude", "Biome"))]

pamLS.scores <- suprow(dudi1, data.frame(pam)[, colnames(env)])$li[,1:2]
pamEnv.scores <- suprow(dudi1, pamEnv)$li[,1:2]
  
pam.grid <- ecospat.grid.clim.dyn(global.scores,
                                  pamEnv.scores,
                                  pamLS.scores)

rm(pam, pamLS, pamEnv, pamLS.scores, pamEnv.scores)

ecospat.plot.niche.dyn(ama.grid, pam.grid,
                       quant = 0.05)

# Test niche equivalency
eqv <- ecospat.niche.equivalency.test(ama.grid, pam.grid,
                               rep=100, alternative = "greater")

sim <- ecospat.niche.similarity.test(ama.grid, pam.grid,
                                     rep=1000, alternative = "lower",
                                     rand.type=2)
ecospat.plot.overlap.test(eqv, "D", "Equivalency")
ecospat.plot.overlap.test(sim, "D", "Similarity")


  #--------------------------------#
  #### Environmental Clustering ####
  #--------------------------------#
  
  # Looking at how species ecological boundaries could be delineated using an 
  #   a priori clustering approach
  
  # Ensure reproducibility of clustering procedure
  set.seed(123)
  
  # First scale the env variables so that they are all comparable
  scaled_env <- scale(env)
  
  # Compute and plot wss for k = 2 to k = 15. This allows for the selection of
  #   an optimal number of groups (k) for points to be clustered into
  k.max <- 15
  wss <- sapply(1:k.max, 
                function(k){kmeans(scaled_env, k, nstart=50,iter.max = 15)$tot.withinss})
  
  plot(1:k.max, wss,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  
  # Try a different approach to k selection
  d_clust <- Mclust(as.matrix(scaled_env), G=1:15, 
                    modelNames = mclust.options("emModelNames"))
  d_clust$BIC
  plot(d_clust)
  
  # Together the above suggests there should be 3 or 4 clusters
  # Cluster the points 
  clusters <- kmeans(scaled_env, 4) # When pulling the cluster assignments (clusters$cluster), the data will be parallel with the scaled_env data and niche.points
  
  # Go back and cluster the ggplot by these assignments, e.g. colour=as.factor(clusters$cluster)
  

# Figuring out how to plot niche ellipses - checked source code for niche: https://rdrr.io/cran/ade4/src/R/niche.R
bs <- ecospat.plot.contrib(contrib=dudi1$co, eigen=dudi1$eig)

# From nicheROVER https://cran.r-project.org/web/packages/nicheROVER/vignettes/ecol-vignette.html
# generate parameter draws from the 'default' posteriors of each fish
library(nicheROVER)

nsamples <- 1000
system.time({
  data.par <- tapply(1:nrow(data), data$Species, function(ii) niw.post(nsamples = nsamples, 
                                                                       X = data[ii, c(4,7,10,17,23)]))
})

# set colours later
clrs <- c(rep('black',13))

# mu plots
par(mar=c(1,1,1,1))
niche.par.plot(data.par, col=clrs, plot.mu=TRUE)

#2-d projecttions of niches
data.data <- tapply(1:nrow(data), data$Species, function(ii) X = data[ii, c(4,7,10,17,23)])

niche.plot(data.par, data.data, pfrac=0.05)
x<-s.distri(nic$ls, eval.parent(as.list(nic$call)[[3]]), 
         cstar = 0, axesell = FALSE, cellipse = 1, sub = "Niches", csub = 2,
         cpoint=0)
# Note the above ellipse is NOT a confidence ellipse; it can be adjusted by changing cellipse.

## Check the ecospat package following: http://www.ecography.org/sites/ecography.org/files/appendix/ecog-02671.pdf
library(ecospat)

# we can use the previously generated dudi.pca
contri