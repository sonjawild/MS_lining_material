
# ## LINING MATERIAL ------------------------------------------------------

# 1) Load libraries -------------------------------------------------------

library(asnipe)
library(sp)
library(geosphere)
library(Imap)
library(data.table)
library(ggplot2)
library(gridExtra)


# 2) Foraging associations ------------------------------------------------


# 2.1. Read network data --------------------------------------------------

load("C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/Foraging network/gmm.spring.RData")


# 2.2. Calculate association network --------------------------------------

foraging_network <-
  get_network(
    association_data = gmm.spring$gbi,
    data_format = "GBI",
    association_index = "SRI"
  )

dim(foraging_network)


# 3) Distance matrices ----------------------------------------------------


# 3.1. Read GPS data ------------------------------------------------------

GPS_data <- read.csv("C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/GPS_data/coordinates_boxes+dispensers_ALL.csv")



# 3.2. Load functions to calculate distances in m -------------------------

##  GeoDistanceInMetresMatrix() function that generates the distance matrix
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}


# 3.3. calculate distances between nest boxes and dispensers --------------

head(GPS_data)

coordinates <- data.frame(name = GPS_data$Box,
                          lat  = GPS_data$Lat,
                          lon  = GPS_data$Long)

distance_matrix<-round(GeoDistanceInMetresMatrix(coordinates))
head(distance_matrix)

# double check the coordinates make sens
plot(coordinates$lon, coordinates$lat)


# 3.4. Extract closest dispenser ------------------------------------------

closest_dispenser <- distance_matrix[,c("D1","D2","D3","D4","D5")] 

closest_disp <- data.frame(cbind(colnames(closest_dispenser)[apply(closest_dispenser, 1, which.min)]),
                    apply(closest_dispenser, 1, min))
colnames(closest_disp) <- c("Closest dispenser", "Distance to dispenser")
head(closest_disp)

# 3.5. Extract which boxes are within 200m of dispensers ------------------
neighbour_matrix <- distance_matrix
# m contains the boxes within 200m of the dispensers
m <- neighbour_matrix[(length(rownames(neighbour_matrix))-4):length(rownames(neighbour_matrix)),1:(length(rownames(neighbour_matrix))-5)]<=200

boxes.to.include <- unique(colnames(m)[col(m)[which(m)]])
length(boxes.to.include)
# 175 boxes are within a 200m radius of the dispensers

# 3.6. Calculate the neighbour matrix -------------------------------------

# we set all values above 50m to 0
# we use the inverted distance

neighbour_matrix <- distance_matrix

neighbour_matrix[neighbour_matrix>=50] <- 0

neighbour_matrix <- 1/neighbour_matrix

# set all infinity values to 0
neighbour_matrix[neighbour_matrix==Inf] <- 0
hist(neighbour_matrix[neighbour_matrix!=0])


# 4) Extract data from dispensers -----------------------------------------

test.tags  <- c("0110178F23", 
                "0700EDB728",
                "0110178F7B",
                "01101790B3",
                "0110179023",
                "0700EDFA90",
                "0700EDBDD2",
                "0110178F75",
                "01101792A2")

# read in/prep dispenser data
read.dispenser.data <- function(path, location) {
  list<-
    list.files(
      path = path,
      all.files = FALSE,
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = FALSE,
      include.dirs = FALSE,
      no.. = FALSE
    )
  
  
  dispenser.data <-
    as.data.frame(rbindlist(
      sapply(
        list,
        read.delim,
        header = F,
        sep = "" ,
        simplify = FALSE
      ),
      use.names = TRUE,
      idcol = FALSE
    ))
  dispenser.data$location <- location
  colnames(dispenser.data) <- c("date.time", "Antenna", "PIT", "Location")
  dispenser.data <- unique(dispenser.data)
  # remove test tags
  dispenser.data <- subset(dispenser.data, !(dispenser.data$PIT %in% test.tags))
  return(dispenser.data)
  
}

# Read in Mill data
dispenser.data.1 <- read.dispenser.data(path="C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/Lining dispenser data/D1", location="D1")
dispenser.data.2 <- read.dispenser.data(path="C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/Lining dispenser data/D2", location="D2")
dispenser.data.3 <- read.dispenser.data(path="C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/Lining dispenser data/D3", location="D3")
dispenser.data.4 <- read.dispenser.data(path="C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/Lining dispenser data/D4", location="D4")
dispenser.data.5 <- read.dispenser.data(path="C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/Lining dispenser data/D5", location="D5")

# 4.2. Extract the number of visits ---------------------------------------

# make a function that extracts separate visits
extract.visits <- function(file, time.window){
  # make a new column called 'visit'
  file$visit <- NA
  # assign a 1 to the very first data point
  file[1,"visit"] <-1
  # write a loop that goes through every line
  for(i in 1:(length(file[,1])-1)){
    # it first extracts the time difference between each two consecutive lines in seconds
    if(difftime(
      as.POSIXct(as.character(file[i + 1, 1]), format = "%y%m%d%H%M%S"),
      as.POSIXct(as.character(file[i, 1]), format = "%y%m%d%H%M%S"),
      units = "secs")<=time.window & file[i+1, 3]==file[i, 3]){
      # if the time difference is no more than 1 second and the ID of the PIT tag is the same
      # it considers it as the same visit
      file[i+1,"visit"] <- file[i, "visit"]
    } else { # otherwise it considers it a next visit
      file[i+1,"visit"] <- file[i, "visit"]+1
    }
  }
  
  # now we can subset the file so that we only keep one row per 'visit'
  file.visits <-   file[match(unique(file$visit), file$visit),]
  return(file.visits)
}

# We can run this function on each file separately
dispenser.data.1 <- extract.visits(file=dispenser.data.1, time.window=15)
dispenser.data.2 <- extract.visits(file=dispenser.data.2, time.window=15)
dispenser.data.3 <- extract.visits(file=dispenser.data.3, time.window=15)
dispenser.data.4 <- extract.visits(file=dispenser.data.4, time.window=15)
dispenser.data.5 <- extract.visits(file=dispenser.data.5, time.window=15)


#Frequency of visits to determine which NA birds are male/female
d1.visits <- as.data.frame(table(dispenser.data.1$PIT))
d2.visits <- as.data.frame(table(dispenser.data.2$PIT)) 
d3.visits <- as.data.frame(table(dispenser.data.3$PIT)) 
d4.visits <- as.data.frame(table(dispenser.data.4$PIT)) 
d5.visits <- as.data.frame(table(dispenser.data.5$PIT)) 


d1.visits$Location <- "D1"
d2.visits$Location <- "D2"
d3.visits$Location <- "D3"
d4.visits$Location <- "D4"
d5.visits$Location <- "D5"

dispenser.visits.all <- rbind(d1.visits, d2.visits, d3.visits, d4.visits, d5.visits)
# write.table(dispenser.visits.all, file="dispenser.visits.txt")

dispenser.visits.all <- read.delim("C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/dispenser.visits.txt")


# Create species and list of females
species_sex <- read.delim("C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/species_sex.txt")

females <- subset(species_sex$Pit, species_sex$Sex.confirmed=="F" | species_sex$Sex.suspected=="F")
# add those inferred manually from dispenser visits
females <- unique(c(females, dispenser.visits.all$PIT[dispenser.visits.all$Sex=="F"]))
# note that this list includes both confirmed and suspected females
females <- females[!females==""]
# 4.2. Extract the order of visit to the dispensers --------------------------

dispenser.data.comb <- rbind(dispenser.data.1, 
                             dispenser.data.2,
                             dispenser.data.3,
                             dispenser.data.4,
                             dispenser.data.5)

# How many different birds have visited
length(unique(dispenser.data.comb$PIT))
# 39

# How many of those are females?

dispenser.data.comb.sub <- subset(dispenser.data.comb, dispenser.data.comb$PIT %in% females)
length(unique(dispenser.data.comb.sub$PIT))
# 25 (note that this also includes suspected females)


# 5) Prepare ILVs ---------------------------------------------------------

breeders <- read.delim("C:/Users/swild/Desktop/Konstanz/Breeding Season 2021/Breeders 2021.txt")
nestgrid <- read.delim("C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/nestgrid.txt", sep="\t")


ILVs <- as.data.frame(matrix(NA, nrow=length(GPS_data[,1])-5, ncol=7))
ILVs[,1] <- GPS_data$Box[1:(length(GPS_data$Box)-5)] # remove the Dispensers (last 5)

colnames(ILVs) <- c("Box", "PIT_f", "Species", "Age", "Lay.Date", "Hatch.Date", "First.visit")


# 5.1. Add the PIT tags of the females and the species ------------------------------------

for(i in ILVs$Box){
PIT.fem <- subset(breeders$Tag, breeders$Box==i & breeders$Who=="Female")
PIT.m <- subset(breeders$Tag, breeders$Box==i & breeders$Who=="Male")
species <- unique(subset(nestgrid$Species, nestgrid$Box_No==i))
if(length(species)==0){
  species <- NA
}
if(length(PIT.fem)==0){
  PIT.fem <- NA
}
if(length(PIT.m)==0){
  PIT.m <- NA
}
# in cases where it is not clear who is the female, check the sex that was manually indicated on the dispenser.visit.all data frame
if(PIT.fem %in% females & PIT.m %in% females){
 # ILVs[which(ILVs$Box==i), "check"] <- "check"
  num.disp.visits.f <- sum(subset(dispenser.visits.all$Freq, dispenser.visits.all$PIT==PIT.fem))
  num.disp.visits.m <- sum(subset(dispenser.visits.all$Freq, dispenser.visits.all$PIT==PIT.m))

 if((num.disp.visits.m-num.disp.visits.f)>=3){ # if the supposed male has visited the dispenser 5 more times than the supposed female, we assume the sex is the other way round
   PIT.temp <- PIT.fem
   PIT.fem <- PIT.m
   PIT.m <- PIT.temp
 }
  
  }

ILVs[which(ILVs$Box==i), "PIT_f"] <- PIT.fem
#ILVs[which(ILVs$Box==i), "PIT_m"] <- PIT.m
ILVs[which(ILVs$Box==i), "Species"] <- species
}

dispenser.visits.all
length(dispenser.visits.all[dispenser.visits.all$Sex=="F",1])
# 24 tagged females visited the dispenser

# 5.2. Add lay date/incubation date/hatch date ----------------------------

for( i in ILVs$Box){
  lay.date <- subset(nestgrid$Lay.date, nestgrid$Box_No==i)
  if(length(lay.date)==0){lay.date <- NA}
  hatch.date <- subset(nestgrid$Observed.hatch.date, nestgrid$Box_No==i)
  if(length(hatch.date)==0){hatch.date <- NA}
  ILVs[which(ILVs$Box==i), "Lay.Date"] <- lay.date 
  ILVs[which(ILVs$Box==i), "Hatch.Date"] <- hatch.date
}


# 5.3. Closest dispenser --------------------------------------------------

ILVs <- cbind.data.frame(ILVs, closest_disp[1:(length(closest_disp[,1])-5),])
colnames(ILVs) <- c("Box", "PIT_f", "Species", "Age", "Lay.Date", "Hatch.Date", "First.visit", "Closest.dispenser", "Distance.to.dispenser")
rownames(ILVs) <- NULL

# extract the maximum distance travelled between nest box and dispenser
max(subset(ILVs, ILVs$PIT_f %in% dispenser.visits.all$PIT)[,"Distance.to.dispenser"])
#178
# extract the average distance travelled
mean(subset(ILVs, ILVs$PIT_f %in% dispenser.visits.all$PIT)[,"Distance.to.dispenser"])
# 88.3

# 5.4. The maximum distance travelled between nest box and dispenser was 178m -----------------------------
# subset to within 200 m (allow for inaccuracy of GPS +-10m of nestboxes and dispensers)
# to tits only
# to boxes with PIT tagged females
# and with a maximum lay date of 40 (which corresponds to the 11th of May when the experiment ended)
ILVs.sub <- subset(
  ILVs,
  ILVs$Distance.to.dispenser < 200 &
    !(is.na(ILVs$Species)) &
    !(is.na(ILVs$PIT_f)) & ILVs$Species != "NUTHA" & ILVs$Lay.Date <= 40
)


# 5.5. Add the females who only visited the dispenser ---------------------

# some females have visited the dispenser, but did not breed in any of the nest boxes
# add those to the ILV file
add.PIT <- subset(dispenser.visits.all, dispenser.visits.all$Sex=="F" & !(dispenser.visits.all$PIT%in% ILVs.sub$PIT_f))
ILVs.add <- ILVs.sub[FALSE,]
ILVs.add[12,] <- NA
ILVs.add$PIT_f <- add.PIT$PIT
ILVs.add$Box <- "no_box"
ILVs.add$Distance.to.dispenser <- 88.3 # this corresponds to the average distance between nest boxes with coloured wool and their closest dispenser

# add species
for(i in ILVs.add$PIT){
  species <- unique(subset(species_sex$Species, species_sex$Pit==i)) 
  ILVs.add[which(ILVs.add$PIT_f==i),"Species"] <- species
  
  closest.disp.i <- subset(dispenser.visits.all$Location, dispenser.visits.all$PIT==i) 
  ILVs.add[which(ILVs.add$PIT_f==i),"Closest.dispenser"] <- closest.disp.i
}


ILVs.combined <- rbind.data.frame(ILVs.sub, ILVs.add)

# 5.6. Add age ------------------------------------------------------------

age <- read.delim("C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/age.data.spring.2021.txt", sep=",", row.names = 1)

for(i in ILVs.combined$PIT_f[!(is.na(ILVs.combined$PIT_f)) & ILVs.combined$PIT_f!="" & ILVs.combined$PIT_f!="xxxxxx"]){
  ILVs.combined[which(ILVs.combined$PIT_f==i), "Age"] <- unique(subset(age$Age.in.2021, age$Pit==i)) 
}



# 5.7. Add info on whether birds visited dispensers the previous year --------



birds.with.prev.exp <- c("011017811F",
"0700EDA1F4",
"0700EDEFDA",
"01101793C7",
"0110178BF4",
"01101798BE",
"0110177F58",
"01101772C2",
"01101788AB",
"0110179303",
"01101775E5",
"011017691B",
"0700EE42F4",
"01101791DC",
"01101760A6",
"0110179961",
"0110178AA6",
"0700EDAC58",
"0110178C37",
"01101793B2",
"0110177ACE",
"0110178427",
"01101789A9",
"0110178B1D",
"0110175B4D",
"011017A3C1",
"01101783E9")

ILVs.combined$Prev.exp <- NA

for(i in ILVs.combined$PIT_f){
  if(i %in% birds.with.prev.exp){
    ILVs.combined[which(ILVs.combined$PIT_f==i),"Prev.exp"] <- "yes"
    
  } else {
    ILVs.combined[which(ILVs.combined$PIT_f==i),"Prev.exp"] <- "no"
  }
}


# 5.8. Add first visit to dispenser -------------------------------------

for(i in ILVs.combined$PIT_f){
  suppressWarnings(
    first.v <- min(subset(dispenser.data.comb$date.time, dispenser.data.comb$PIT==i))  
  )
if(is.infinite(first.v)){
  first.v <- NA
}  
  ILVs.combined[which(ILVs.combined$PIT_f==i), "First.visits"] <- first.v

}


# 5.9. Add whether in 200m radius of each dispenser -----------------------
ILVs.combined$D1 <- NA
ILVs.combined$D2 <- NA
ILVs.combined$D3 <- NA
ILVs.combined$D4 <- NA
ILVs.combined$D5 <- NA

for(i in ILVs.combined$Box[ILVs.combined$Box!="no_box"]){
  
  for(j in c("D1", "D2", "D3", "D4", "D5")){
    dist <- distance_matrix[i, j]
    if(dist<=200){
      ILVs.combined[which(ILVs.combined$Box==i), j] <- 1
    } else {
      ILVs.combined[which(ILVs.combined$Box==i), j] <- 0
    }
  }
  # also check whether the females have visited other dispensers and add them to those areas
  PIT.i <- subset(ILVs.combined$PIT_f, ILVs.combined$Box==i)
  dips.visited <- unique(subset(dispenser.data.comb$Location, dispenser.data.comb$PIT==PIT.i))
  for(l in dips.visited){
    ILVs.combined[which(ILVs.combined$PIT_f==PIT.i), l] <- 1
  }
    
}

# for those not breeding in boxes
# check which dispensers they have visited and add them to those radii

for(i in ILVs.combined$PIT_f[ILVs.combined$Box=="no_box"]){
  dips.visited <- subset(dispenser.data.comb$Location, dispenser.data.comb$PIT==i)
  for(j in dips.visited){
    ILVs.combined[which(ILVs.combined$PIT_f==i), j] <- 1
  }

}


# 5.9. Add information on the use of wool ------------------------------

ILVs.combined$Wool.found

ILVs.combined$First.col

ILVs.combined$Second.col

ILVs.combined$First.provided

ILVs.combined$Demo


# 6. NBDA -----------------------------------------------------------------


# 6.2. Social information use to find lining material ---------------------
# here we include all females that are part of our ILVs_combined



# 6.2.1. prepare matrices -------------------------------------------------


IDs.to.include.in.NBDA <- intersect(ILVs.combined$PIT_f, rownames(foraging_network))
# create a second where we only include the females breeding in the boxes
# IDs.to.include.in.NBDA <- intersect(ILVs.combined$PIT_f[ILVs.combined$Box!="no_box"], rownames(foraging_network))
length(IDs.to.include.in.NBDA)
# 46
foraging.network.NBDA <- foraging_network[rownames(foraging_network) %in% IDs.to.include.in.NBDA, colnames(foraging_network) %in% IDs.to.include.in.NBDA]
foraging.network.NBDA <- foraging.network.NBDA[order(rownames(foraging.network.NBDA)), order(colnames(foraging.network.NBDA))]
dim(foraging.network.NBDA)
hist(foraging.network.NBDA)


boxes.to.include.in.NBDA <- unique(subset(ILVs.combined$Box, ILVs.combined$PIT_f %in% IDs.to.include.in.NBDA))
boxes.to.include.in.NBDA <- boxes.to.include.in.NBDA[boxes.to.include.in.NBDA!="no_box"]
new.names.all <- NULL
neighbour_matrix.NBDA <- neighbour_matrix[rownames(neighbour_matrix) %in% boxes.to.include.in.NBDA, colnames(neighbour_matrix) %in% boxes.to.include.in.NBDA]
for(i in rownames(neighbour_matrix.NBDA)){
  new.name <- subset(ILVs.combined$PIT_f, ILVs.combined$Box==i)
  new.names.all[which(rownames(neighbour_matrix.NBDA)==i)] <- new.name
}
rownames(neighbour_matrix.NBDA) <- new.names.all
colnames(neighbour_matrix.NBDA) <- new.names.all
dim(neighbour_matrix.NBDA)
# this matrix only contains 36 individuals (excludes the ones that weren't breeding in nest boxes)
# we add those individuals and assign the average distance

neighbour_matrix.NBDA.new <- matrix(mean(neighbour_matrix.NBDA), ncol=length(rownames(foraging.network.NBDA)), nrow=length(rownames(foraging.network.NBDA)), dimnames = list(rownames(foraging.network.NBDA),rownames(foraging.network.NBDA)))
# now fill in the real distance values for those available
cols <- colnames(neighbour_matrix.NBDA.new)[colnames(neighbour_matrix.NBDA.new) %in% colnames(neighbour_matrix.NBDA)]
rows <- rownames(neighbour_matrix.NBDA.new)[rownames(neighbour_matrix.NBDA.new) %in% rownames(neighbour_matrix.NBDA)]

neighbour_matrix.NBDA.new[rows, cols] <- neighbour_matrix.NBDA[rows,cols]
dim(neighbour_matrix.NBDA.new)
hist(neighbour_matrix.NBDA.new)
  
# to ensure they are on a similar scale, we multiply the values in the neighbour matrix *10
neighbour_matrix.NBDA.new <- neighbour_matrix.NBDA.new*10
hist(neighbour_matrix.NBDA.new)

# mantel(neighbour_matrix.NBDA.new, foraging.network.NBDA, permutations = 999)
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = neighbour_matrix.NBDA.new, ydis = foraging.network.NBDA,      permutations = 999) 
# 
# Mantel statistic r: 0.02521 
#       Significance: 0.204 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0452 0.0660 0.0825 0.0947 
# Permutation: free
# Number of permutations: 999


# 6.2.1. Prepare networks into array ------------------------------------------------------------------
assMatrix <- array(data = c(foraging.network.NBDA, neighbour_matrix.NBDA.new), dim=c(nrow(foraging.network.NBDA), ncol(foraging.network.NBDA), 2))
dim(assMatrix)

# 6.2.2. Prepare ILVs ------------------------------------------------------------------

prepare.NBDA.data <- function(dispenser.data){
  dispenser.data <- dispenser.data[order(dispenser.data$date.time),] # ensure it is sorted according to date
  location <- unique(dispenser.data$Location)
  
  
  
  
  # subset to those in the correct dispenser area and those who were seen at the network feeders
  ILVs.sub.disp <- subset(ILVs.combined, ILVs.combined[,location]==1 & ILVs.combined$PIT_f %in% IDs.to.include.in.NBDA)
  
  unique(dispenser.data$PIT) %in% ILVs.sub.disp$PIT_f
  
  IDs.included <- ILVs.sub.disp$PIT_f
  
  dispenser.data <- subset(dispenser.data, dispenser.data$PIT %in% IDs.included)
  
  # order data ascending
  ILVs.sub.disp <- ILVs.sub.disp[order(ILVs.sub.disp$PIT_f),] # order ascending
  
  # subset the two networks to those IDs
  forage.net <- foraging.network.NBDA[rownames(foraging.network.NBDA) %in% IDs.included, colnames(foraging.network.NBDA) %in% IDs.included]
  neighbour.net <- neighbour_matrix.NBDA.new[rownames(neighbour_matrix.NBDA.new) %in% IDs.included, colnames(neighbour_matrix.NBDA.new) %in% IDs.included]
  
  assMatrix.nbda <- array(data=c(forage.net, neighbour.net), dim=c(nrow(forage.net), ncol(forage.net), 2))
  
  # create objects in the global environment for each ILV
  species.nbda <- as.matrix(ILVs.sub.disp$Species) 
  species.nbda[species.nbda!="GRETI"] <- -0.5 
  species.nbda[species.nbda=="GRETI"] <- 0.5 
  species.nbda <- as.matrix(as.numeric(species.nbda))
  
  age.nbda <- as.matrix(ILVs.sub.disp$Age) 
  age.nbda[age.nbda=="first.year"] <- -0.5 
  age.nbda[age.nbda=="adult"] <- 0.5 
  age.nbda <- as.matrix(as.numeric(age.nbda))
  
  distance.nbda <- as.matrix(as.numeric(scale(ILVs.sub.disp$Distance.to.dispenser))) 
  
  prev.exp.nbda <- ILVs.sub.disp$Prev.exp
  prev.exp.nbda[prev.exp.nbda=="yes"] <- 0.5
  prev.exp.nbda[prev.exp.nbda=="no"] <- -0.5
  prev.exp.nbda <- as.matrix(as.numeric(prev.exp.nbda))
  
  assign(paste("species", location, sep="_"), species.nbda, envir = .GlobalEnv)
  assign(paste("age", location, sep="_"), age.nbda, envir = .GlobalEnv)
  assign(paste("distance", location, sep="_"), distance.nbda, envir = .GlobalEnv)
  assign(paste("prev.exp", location, sep="_"), prev.exp.nbda, envir = .GlobalEnv)
  
  ILVs <- paste(c("species", "age", "distance", "prev.exp"), location, sep="_")
#  ILVs <- paste(c("distance"), location, sep="_")
 # ILVs <- paste(c("species", "distance"), location, sep="_")
  assign(paste("ILVs", location, sep="_"), ILVs, envir = .GlobalEnv)
  
  # extract the order of finding
  t.first <- dispenser.data[match(unique(dispenser.data$PIT), dispenser.data$PIT),]
  order <- NULL
  time <- NULL
  num.visits <- NULL
  for( i in t.first$PIT){
    order[which(t.first$PIT==i)] <- which(rownames(forage.net)==i)
    time[which(t.first$PIT==i)] <- as.POSIXct(as.character(subset(t.first$date.time, t.first$PIT==i)), format="%y%m%d%H%M%S", origin="1970-01-01")-as.POSIXct("21032612300000", format="%y%m%d%H%M%S") # difference in days
    num.visits[which(t.first$PIT==i)] <- as.numeric(subset(dispenser.visits.all$Freq, dispenser.visits.all$PIT==i & dispenser.visits.all$Location==location))
    # 26.03.21 12:30 CEST
  }
  object <- NULL
  object$forage.net <- forage.net
  object$neighbour.net <- neighbour.net
  object$ILVs.full <- ILVs.sub.disp
  object$assMatrix <- assMatrix.nbda
  
  object$OAc <- order
  object$TAc <- time
  object$transmission.weights <- rep(0, length(rownames(forage.net)))
  object$transmission.weights[order] <- num.visits
  
  return(object)
}

d1.NBDA <- prepare.NBDA.data(dispenser.data = dispenser.data.1)
d2.NBDA <- prepare.NBDA.data(dispenser.data = dispenser.data.2)
d3.NBDA <- prepare.NBDA.data(dispenser.data = dispenser.data.3)
d4.NBDA <- prepare.NBDA.data(dispenser.data = dispenser.data.4)
d5.NBDA <- prepare.NBDA.data(dispenser.data = dispenser.data.5)


# 6.2.4. Prepare NBDA Data Objects ----------------------------------------
#install.packages("devtools")
#install.packages("rtools40")
#install.packages("ade4")
library(devtools)
# install_github("whoppitt/NBDA")
library(NBDA)


nbdaData_D1 <- nbdaData(label="D1",                        # specify an informative label
                            assMatrix=d1.NBDA$assMatrix,           # our array with the matrices
                            asoc_ilv=get(paste("ILVs", "D1", sep="_")),            # we specify that ILVs can influence asocial learning, if no ILV: "ILVabsent"
                            int_ilv=get(paste("ILVs", "D1", sep="_")),             # we specify that our ILVs can influence social learning, if no ILV: "ILVabsent" 
                            multi_ilv="ILVabsent",        # we specify that our ILVs can influence asocial and social learning to the same extent, if no ILV: "ILVabsent" 
                            orderAcq=d1.NBDA$OAc,          # vector with the order of acquisition 
                            timeAcq=d1.NBDA$TAc,           # numerical vector giving the time at which each individual acquired the target behaviour, given in the order matching orderAcqv
                            endTime=41                    # numeric giving the time at which the diffusion ended. (11th of May = day 40)
                            )

nbdaData_D2 <- nbdaData(label="D2",                        
                        assMatrix=d2.NBDA$assMatrix,          
                        asoc_ilv=get(paste("ILVs", "D2", sep="_")),            
                        int_ilv=get(paste("ILVs", "D2", sep="_")),            
                        multi_ilv="ILVabsent",        
                        orderAcq=d2.NBDA$OAc,          
                        timeAcq=d2.NBDA$TAc,           
                        endTime=41
)


nbdaData_D3 <- nbdaData(label="D3",                        
                        assMatrix=d3.NBDA$assMatrix,          
                        asoc_ilv=get(paste("ILVs", "D3", sep="_")),            
                        int_ilv=get(paste("ILVs", "D3", sep="_")),            
                        multi_ilv="ILVabsent",        
                        orderAcq=d3.NBDA$OAc,          
                        timeAcq=d3.NBDA$TAc,           
                        endTime=41
)



# skip 4

nbdaData_D5 <- nbdaData(label="D5",                        
                        assMatrix=d5.NBDA$assMatrix,          
                        asoc_ilv=get(paste("ILVs", "D5", sep="_")),            
                        int_ilv=get(paste("ILVs", "D5", sep="_")),            
                        multi_ilv="ILVabsent",        
                        orderAcq=d5.NBDA$OAc,          
                        timeAcq=d5.NBDA$TAc,           
                        endTime=41
)

# create function for making constraintsVectorMatrix

## 03.08.2020: this version of the code can create the matrix for unconstrained and/or multiplicative models
# it takes all the info straight from the NBDA data object
# I have not tried multiple networks yet - will implement that in the next version


# 6.2.5. Create constraints Vector MAtrix ---------------------------------


create.constraints.Vect.Matrix <- function(NBDA_data_object, num_networks, num_ILVs){
  suppressWarnings(
    if(NBDA_data_object@asoc_ilv=="ILVabsent"){
      num.ILV.asoc <- 0
    } else {num.ILV.asoc <- length(NBDA_data_object@asoc_ilv)})
  
  suppressWarnings(
    if(NBDA_data_object@int_ilv=="ILVabsent"){
      num.ILV.int<- 0
    } else {num.ILV.int<- length(NBDA_data_object@int_ilv)})
  
  suppressWarnings(
    if(NBDA_data_object@multi_ilv=="ILVabsent"){
      num.ILV.multi <- 0
    } else {num.ILV.multi <- length(NBDA_data_object@multi_ilv)})
  
  vector <- seq(1:(num_networks+num.ILV.asoc+num.ILV.int+num.ILV.multi))
  
  count <- 0 # create an object 'count', which starts on 0
  
  constraintsVect <- matrix(nrow = 10000000, ncol=(num_networks+num.ILV.asoc+num.ILV.int+num.ILV.multi)) # create a matrix to save the combination of parameters in
  constraintsVect[1,] <- vector # the first row gets filled with a sequence from 1:8 (all parameters will be estimated, none are set to 0)
  
  for (i in 1:(length(vector)-1)){ # a loop for each number of parameters to be estimated
    array <- combn(vector, i, FUN = NULL, simplify = TRUE) # for each number of paramters to be estiamted (e.g. 2) create all possible combinations of numbers between 1:12 (e.g. 2&8, 1&5 etc)
    
    for (j in 1:length(array[1,])){ # for each of those combinations
      vector2 <- seq(1:((num_networks+(num.ILV.asoc+num.ILV.int+num.ILV.multi))-i)) # create a second vector with 11-i free spaces
      position <- array[,j] # for each created combination
      count <- count+1 # add +1 to the count
      
      for (k in position){ # at each possible position
        vector2 <- append(vector2, 0, after=k-1) # add a 0 (e.g. 1 0 2 3 ...; 1 2 0 3 4 5 ...; 1 2 3 0 4 5 ....)
      }
      constraintsVect[count+1,] <- vector2 # and save the resulting order in a matrix
    }
  }
  
  
  constraintsVect <- na.omit(constraintsVect) # remove all NAs from the matrix
  
  # extract which columns are networks
  col.networks <- c(1:num_networks)
  
  col.names <- NULL
  
  if(num.ILV.asoc!=0){
    col.names <- rep("asoc", num.ILV.asoc)
  }
  
  if(num.ILV.int!=0){
    col.names <- c(col.names, rep("int", num.ILV.int))
  }
  
  if(num.ILV.multi!=0){
    col.names <- c(col.names, rep("multi", num.ILV.multi))
  }
  
  colnames(constraintsVect) <- c(rep("network", num_networks), col.names)
  
  constraintsVect <- as.matrix(as.data.frame(constraintsVect))
  
  # extract the models containing any social network
  
  social.models <- rep(NA, length(constraintsVect[,1]))
  
  for (k in 1:length(constraintsVect[,1])){
    sum <- sum(constraintsVect[k,1:num_networks])
    if(sum!=0){
      social.models[k] <- k
    }
  }
  social.models <- as.vector(na.omit(social.models))
  
  social.models.matrix <- constraintsVect[social.models,]
  
  # if multiplicative models are fit, we need to adjust the matrix
  # if the multiplicative slots are filled, it automatically fits the parameter for asoc and social (just constrained to be the same)
  # meaning that we can remove it from the asoc and int slot
  
  if(num.ILV.multi!=0){
    social.models.retain <- rep(NA, length(social.model.matrix[,1]))
    multi.models <- rep(NA, length(social.models.matrix[,1]))
    for (k in 1:length(social.models.matrix[,1])){
      sum <- sum(social.models.matrix[k,which(colnames(social.models.matrix)=="multi")])
      sum2 <- sum(social.models.matrix[k, c(which(colnames(social.models.matrix)=="asoc"),which(colnames(social.models.matrix)=="int"))])
      if(sum!=0 & sum2==0){ # if multi models are fit and int and asoc are set to 0
        multi.models[k] <- k # then retain the model
      } else if (sum==0){
        social.models.retain[k] <- k
      }
    }
    
    multi.models <- as.vector(na.omit(multi.models))
    social.models.retain <- as.vector(na.omit(social.models.retain))
    
    models.to.retain <- c(multi.models, social.models.retain)
    
    # these models are retained
    retain.matrix.soc <- social.models.matrix[models.to.retain,]
    
    social.models.matrix <- retain.matrix.soc
  }
  
  # extract the models containing no social network
  
  asocial.models <- rep(NA, length(constraintsVect[,1]))
  
  for (k in 1:length(constraintsVect[,1])){
    sum <- sum(constraintsVect[k,1:num_networks])
    if(sum==0){
      asocial.models[k] <- k
    }
  }
  asocial.models <- as.vector(na.omit(asocial.models))
  
  asocial.models.matrix <- constraintsVect[asocial.models,]
  
  cols.asoc <- which(colnames(constraintsVect)=="asoc")
  
  asocial.retain <- rep(NA, length(asocial.models))
  for (k in 1:length(asocial.models)){
    sum <- sum(asocial.models.matrix[k,which(colnames(constraintsVect)!="asoc")])
    if(sum==0){
      asocial.retain[k] <- k
    }
  }
  
  
  asocial.retain <- as.vector(na.omit(asocial.retain))
  
  asocial.models.to.retain <- asocial.models.matrix[asocial.retain, ]
  asocial.models.to.retain.matrix <- as.matrix(asocial.models.to.retain)
  constraintsVectMatrix <- rbind(social.models.matrix,asocial.models.to.retain)
  
  # add the Null model (without social learning, and no ILVs)
  constraintsVectMatrix <- rbind(constraintsVectMatrix, rep(0, length(constraintsVectMatrix[1,])))
  
  row.names(constraintsVectMatrix) <- NULL
  return(constraintsVectMatrix)
}


constraintsVectMatrix <- create.constraints.Vect.Matrix(NBDA_data_object = nbdaData_D2, num_networks = 2, num_ILVs = 4)

###########

nbdaData_D1_noILV <- nbdaData(label="D1",                        # specify an informative label
                        assMatrix=d1.NBDA$assMatrix,           # our array with the matrices
                        asoc_ilv="ILVabsent",            # we specify that ILVs can influence asocial learning, if no ILV: "ILVabsent"
                        int_ilv="ILVabsent",             # we specify that our ILVs can influence social learning, if no ILV: "ILVabsent" 
                        multi_ilv="ILVabsent",        # we specify that our ILVs can influence asocial and social learning to the same extent, if no ILV: "ILVabsent" 
                        orderAcq=d1.NBDA$OAc,          # vector with the order of acquisition 
                        timeAcq=d1.NBDA$TAc,           # numerical vector giving the time at which each individual acquired the target behaviour, given in the order matching orderAcqv
                        endTime=41                    # numeric giving the time at which the diffusion ended. (11th of May = day 40)
)

nbdaData_D2_noILV <- nbdaData(label="D2",                        
                        assMatrix=d2.NBDA$assMatrix,          
                        asoc_ilv="ILVabsent",            
                        int_ilv="ILVabsent",            
                        multi_ilv="ILVabsent",        
                        orderAcq=d2.NBDA$OAc,          
                        timeAcq=d2.NBDA$TAc,           
                        endTime=41
)


nbdaData_D3_noILV <- nbdaData(label="D3",                        
                        assMatrix=d3.NBDA$assMatrix,          
                        asoc_ilv="ILVabsent",            
                        int_ilv="ILVabsent",            
                        multi_ilv="ILVabsent",        
                        orderAcq=d3.NBDA$OAc,          
                        timeAcq=d3.NBDA$TAc,           
                        endTime=41
)



# skip 4

nbdaData_D5_noILV <- nbdaData(label="D5",                        
                        assMatrix=d5.NBDA$assMatrix,          
                        asoc_ilv="ILVabsent",            
                        int_ilv="ILVabsent",            
                        multi_ilv="ILVabsent",        
                        orderAcq=d5.NBDA$OAc,          
                        timeAcq=d5.NBDA$TAc,           
                        endTime=41
)


TADA.finding.noILV <-
  tadaFit(
    nbdadata = list(
      nbdaData_D1_noILV,
      nbdaData_D2_noILV,
      nbdaData_D3_noILV,
      nbdaData_D5_noILV),
    type="social")

TADA.finding.noILV@optimisation

start.vals <- TADA.finding.noILV@outputPar


# 6.2.6. Run TADA ---------------------------------------------------------
# get start values



TADA.finding.all.start.val <-
  tadaAICtable(
    nbdadata = list(
      nbdaData_D1,
      nbdaData_D2,
      nbdaData_D3,
      nbdaData_D5),
    constraintsVectMatrix = constraintsVectMatrix, 
    writeProgressFile = F,
    startValue = c(start.vals, rep(0, 8))
    ) # start values correspond to: baseline, s1, s2, 3xasoc_ILV, 3xsoc_ILV

# view the results with 'print'
print(TADA.finding.full@printTable) # with all 4 ILVs
TADA.finding.full@printTable <- subset(TADA.finding.full@printTable, TADA.finding.full@printTable$convergence==T)


print(TADA.finding.all.start.val@printTable) # with 4 ILVs and start values
subset(TADA.finding.all.start.val@printTable, TADA.finding.all.start.val@printTable$convergence==T)

subset(TADA.finding@printTable, TADA.finding@printTable$convergence==T)

remove.unfitted <- function(object){
  # Create a new object with a printTable that excludes unfitted models
  newobject<-object
  newobject@printTable<-object@printTable[!is.nan(object@printTable$aicc)&!is.na(object@printTable$aicc),]
  
  # remove those that did not converge
  newobject@printTable <- subset(newobject@printTable, newobject@printTable$convergence==T)
  
  
  # recalculate model support only including models that were fitted
  newobject@aicc<-object@aicc[!is.nan(object@aicc)&!is.na(object@aicc)]
  newobject@MLEs<-object@MLEs[!is.nan(object@aicc)&!is.na(object@aicc),]
  newobject@MLEilv<-object@MLEilv[!is.nan(object@aicc)&!is.na(object@aicc),]
  newobject@MLEint<-object@MLEint[!is.nan(object@aicc)&!is.na(object@aicc),]
  
  newobject@printTable<-newobject@printTable[order(newobject@printTable$aicc),]
  newobject@printTable$deltaAICc<-newobject@printTable$aicc-newobject@printTable$aicc[1]
  
  newobject@printTable$RelSupport<- exp(-0.5*newobject@printTable$deltaAICc)
  newobject@printTable$AkaikeWeight<-newobject@printTable$RelSupport/sum(newobject@printTable$RelSupport)
  
  
  newobject@deltaAIC<-newobject@aicc-min(newobject@aicc)
  
  # calculate model support and akaike weights for each model
  newobject@RelSupport<- exp(-0.5*newobject@deltaAIC)
  newobject@AkaikeWeight<-newobject@RelSupport/sum(newobject@RelSupport)
  
  
  # extract the number of unfitted models that were removed. 
  dim(object@printTable)[1]-dim(newobject@printTable)[1]
  dim(object@printTable)[1]
  ## In our data set, 902 models could not be fitted put of 17'984
  # they likely have too many parameters for the data set 
  

  
  return(newobject)
  
}


TADA.full.start.val.new <- remove.unfitted(object=TADA.finding.all.start.val)

TADA.full.start.val.new@printTable
subset(TADA.full.start.val.new@printTable, TADA.full.start.val.new@printTable$convergence==T)


# network support
networksSupport <- networksSupport(TADA.full.start.val.new)
round(networksSupport, 2)

# variable support
variable_support <- variableSupport(TADA.full.start.val.new, includeAsocial = T)
round(variable_support,3)

# model averaged estimates
MLE_med  <- modelAverageEstimates(TADA.full.start.val.new , averageType = "median")
round(MLE_med,2)


# 6.3. Extract efect sizes ------------------------------------------------

bestModelData1 <- constrainedNBDAdata(nbdadata=nbdaData_D1,constraintsVect =constraintsVectMatrix[112,])
bestModelData2 <- constrainedNBDAdata(nbdadata=nbdaData_D2,constraintsVect =constraintsVectMatrix[112,])
bestModelData3 <- constrainedNBDAdata(nbdadata=nbdaData_D3,constraintsVect =constraintsVectMatrix[112,])
#bestModelData4 <- constrainedNBDAdata(nbdadata=nbdaData_D4,constraintsVect =constraintsVectMatrix[28,])
bestModelData5 <- constrainedNBDAdata(nbdadata=nbdaData_D5,constraintsVect =constraintsVectMatrix[112,])


model.best.social <-
  tadaFit(
    list(
      bestModelData1,
      bestModelData2,
      bestModelData3,
    #  bestModelData4,
      bestModelData5
    )
  )
model.best.social@outputPar
model.best.social@varNames
# [1] "Scale (1/rate):"         "1 Social transmission 1" "2 Social: ILV_dist" 
# [1]     925.3983540   7.6147362  -0.6511193

model.best.social@optimisation

# extract the % of events occurring through social learning
prop.solve.social.byevent <-
  oadaPropSolveByST.byevent(
    nbdadata = list(
      bestModelData1,
      bestModelData2,
      bestModelData3,
  #   bestModelData4,
      bestModelData5
    ),
    model = model.best.social
  )
prop.solve.social.byevent 

# this extracts the overall percentage that have learned socially
prop.solve.social <-
  oadaPropSolveByST(
    nbdadata = list(
      bestModelData1,
      bestModelData2,
      bestModelData3,
 #     bestModelData4,
      bestModelData5
    ),
    model = model.best.social
  )
prop.solve.social # P=0.24

# extract profile likelihood. which=1 extracts the first parameter (in this case s for the vertical network)
# (in this case the s parameter for vertical social learning)
plotProfLik(which=1,model=model.best.social,range=c(0,200), resolution=10) # start with large range
plotProfLik(which=1,model=model.best.social,range=c(0,30), resolution=5) # zoom in on lower range
CIs <- profLikCI(which=1,model=model.best.social,lowerRange=c(5,10), upperRange = c(20,30)) # extract confidence intervals
CIs
# Lower CI  Upper CI 
# 5.000078 24.787620 
#which = refers to the number of models (see model@varNames)

######### lower and upper bound in %

#To get the estimates for the lower bound we should really find the corresponding value of the other parameters to plug in when s1 is constrained to this value
bestModelDataS1LowerBound.D1 <- constrainedNBDAdata(
  nbdadata =
    nbdaDataWOOL_D1,
  constraintsVect = constraintsVectMatrix[58, ],
  offset = c(CIs[1] , rep(0, 6))
)

bestModelDataS1LowerBound.D2 <- constrainedNBDAdata(
  nbdadata =
    nbdaDataWOOL_D2,
  constraintsVect = constraintsVectMatrix[58, ],
  offset = c(CIs[1] , rep(0, 6))
)

bestModelDataS1LowerBound.D3 <- constrainedNBDAdata(
  nbdadata =
    nbdaDataWOOL_D3,
  constraintsVect = constraintsVectMatrix[58, ],
  offset = c(CIs[1] , rep(0, 6))
)

bestModelDataS1LowerBound.D4 <- constrainedNBDAdata(
  nbdadata =
    nbdaDataWOOL_D4,
  constraintsVect = constraintsVectMatrix[58, ],
  offset = c(CIs[1] , rep(0, 6))
)

bestModelDataS1LowerBound.D5 <- constrainedNBDAdata(
  nbdadata =
    nbdaDataWOOL_D5,
  constraintsVect = constraintsVectMatrix[58, ],
  offset = c(CIs[1] , rep(0, 6))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then the value of s at the lower bound is added to s as an offset
bestModelS1LowerBound <-
  tadaFit(
    list(
      bestModelDataS1LowerBound.D1,
      bestModelDataS1LowerBound.D2,
      bestModelDataS1LowerBound.D3,
      bestModelDataS1LowerBound.D4,
      bestModelDataS1LowerBound.D5
    ) ,
    type = "asocial"
  )
bestModelS1LowerBound@outputPar
# [1] 822.5877818  -0.7606382
#Now plug into the prop solve function
prop.solve.social.lower.dial <-
  oadaPropSolveByST(
    model = bestModelS1LowerBound.dial,
    nbdadata = list(
      bestModelDataS1LowerBound.D1,
      bestModelDataS1LowerBound.D2,
      bestModelDataS1LowerBound.D3,
      bestModelDataS1LowerBound.D4,
      bestModelDataS1LowerBound.D5
    )
  )
prop.solve.social.lower.dial
# lower bound for % of birds having learned the dial task through social learning is 44.8%





# 7) Wool choice ----------------------------------------------------------








