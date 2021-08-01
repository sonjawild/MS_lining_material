
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




# 5.3. Subset -----------------------------
# to tits only
# to boxes within 200m of each dispenser
# to boxes with PIT tagged females
# and with a maximum lay date of 40 (which corresponds to the 11th of May when the experiment ended)
# or no lay date at all (for those who have just built the nest but did not lay eggs)
ILVs.sub <- subset(
  ILVs,
    !(is.na(ILVs$Species)) &
    !(is.na(ILVs$PIT_f)) & ILVs$Species != "NUTHA" 
)

# subset it to those with lay date before or on the 11th of May, plus those who never laid eggs (NAs)
 ILVs.sub <- rbind(subset(ILVs.sub, ILVs.sub$Lay.Date<=41), subset(ILVs.sub, is.na(ILVs.sub$Lay.Date)))

# 5.4. Add the females who only visited the dispenser ---------------------

# some females have visited the dispenser, but did not breed in any of the nest boxes
# add those to the ILV file
add.PIT <- subset(dispenser.visits.all, dispenser.visits.all$Sex=="F" & !(dispenser.visits.all$PIT%in% ILVs.sub$PIT_f))
ILVs.add <- ILVs.sub[FALSE,]
ILVs.add[11,] <- NA
ILVs.add$PIT_f <- add.PIT$PIT
ILVs.add$Box <- "no_box"

# add species
for(i in ILVs.add$PIT){
  species <- unique(subset(species_sex$Species, species_sex$Pit==i)) 
  ILVs.add[which(ILVs.add$PIT_f==i),"Species"] <- species

}


ILVs.combined <- rbind.data.frame(ILVs.sub, ILVs.add)

# add species manually for two birds
ILVs.combined[ILVs.combined$PIT_f%in% c("0700EDAD6F", "0700EDEFDA"),"Species"] <- "GRETI"


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


ILVs.combined$D1.visited <- 0
ILVs.combined$D2.visited <- 0
ILVs.combined$D3.visited <- 0
ILVs.combined$D4.visited <- 0
ILVs.combined$D5.visited <- 0


for(i in ILVs.combined$Box[ILVs.combined$Box!="no_box"]){
  
  for(j in c("D1", "D2", "D3", "D4", "D5")){
    dist <- distance_matrix[i, j]
    ILVs.combined[which(ILVs.combined$Box==i), j] <- dist
    
  }

    
}

# for those not breeding in boxes
# check which dispensers they have visited and add them to those radii

for(i in ILVs.combined$PIT_f){
  disp.visited <- unique(subset(dispenser.data.comb$Location, dispenser.data.comb$PIT==i))
  
  if(length(disp.visited)!=0){
    cols <- paste(disp.visited, ".visited", sep="")  
    ILVs.combined[which(ILVs.combined$PIT_f==i), cols] <- 1}

}


# extract the maximum distance travelled between nest box and dispenser
distances.all <- NULL
for(i in c("D1", "D2", "D3", "D4", "D5")){
  distances.i <- subset(ILVs.combined,  ILVs.combined[,paste(i, ".visited", sep="")]==1 & !is.na(ILVs.combined[,i]))[,i]
distances.all <- c(distances.all, distances.i)
  }
hist(distances.all)
max(distances.all)

# subset ILVs to those that are within a maximum of 200m of one of the dispensers
ILVs.combined$closest.dispenser <- NA
colnames(closest_disp) <- c("closest_disp", "distance_to_closest_disp")

for(i in ILVs.combined$Box[ILVs.combined$Box!="no_box"]){
  dist <- subset(closest_disp$distance_to_closest_disp, rownames(closest_disp)==i)
  ILVs.combined[which(ILVs.combined$Box==i), "closest.dispenser"] <- dist
  
}
head(ILVs.combined)
tail(ILVs.combined)
ILVs.combined <- subset(ILVs.combined, ILVs.combined$closest.dispenser<=200 | is.na(ILVs.combined$closest.dispenser))

# 6. NBDA - social information to use to find lining material -----------------------------------------------------------------

#install.packages("devtools")
#install.packages("rtools40")
#install.packages("ade4")
library(devtools)
# install_github("whoppitt/NBDA")
library(NBDA)

# 6.1. prepare matrices -------------------------------------------------

# extract females that were breeding in our ILV data frame and were part of the foraging network
IDs.to.include.in.NBDA <- intersect(ILVs.combined$PIT_f, rownames(foraging_network))
length(IDs.to.include.in.NBDA)
# 46
# subset the foraging network to these females
foraging.network.NBDA <- foraging_network[rownames(foraging_network) %in% IDs.to.include.in.NBDA, colnames(foraging_network) %in% IDs.to.include.in.NBDA]
foraging.network.NBDA <- foraging.network.NBDA[order(rownames(foraging.network.NBDA)), order(colnames(foraging.network.NBDA))]
dim(foraging.network.NBDA)
hist(foraging.network.NBDA)

# subset the neighbour netowrk to these boxes
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
# this matrix only contains 43 individuals (excludes the ones that weren't breeding in nest boxes)
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

# 6.2. Check for correlation ----------------------------------------------

library(vegan)
 mantel(neighbour_matrix.NBDA.new, foraging.network.NBDA, permutations = 9999)
#  Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = neighbour_matrix.NBDA.new, ydis = foraging.network.NBDA,      permutations = 9999) 
# 
# Mantel statistic r: 0.05238 
#       Significance: 0.0651 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0404 0.0595 0.0758 0.0980 
# Permutation: free
# Number of permutations: 9999


# 6.3. Prepare ILVs ------------------------------------------------------------------

prepare.NBDA.data <- function(dispenser.data, include.all, ILVs.include){
  dispenser.data <- dispenser.data[order(dispenser.data$date.time),] # ensure it is sorted according to date
  location <- unique(dispenser.data$Location)
  # if(include.all == TRUE){ # here, we include all individuals (including those that were only seen at the dispensers)
  #   # subset to those in the correct dispenser area and those who were seen at the network feeders
  #   ILVs.sub.disp <- subset(ILVs.combined, ILVs.combined[,location]<=200 | ILVs.combined[, paste(location, ".visited", sep="")] )
  #   ILVs.sub.disp <- subset(ILVs.sub.disp, ILVs.sub.disp$PIT_f %in% IDs.to.include.in.NBDA)  
  # } else { # or we can only subset it to females breeding in known nest boxes
  #   ILVs.sub.disp <- subset(ILVs.combined, ILVs.combined[,location]<=200 | ILVs.combined[, paste(location, ".visited", sep="")] )
  #   ILVs.sub.disp <- subset(ILVs.sub.disp, ILVs.sub.disp$PIT_f %in% IDs.to.include.in.NBDA & ILVs.sub.disp$Box!="no_box")  
  # }
  # 
  # 
  # IDs.included <- ILVs.sub.disp$PIT_f
  
  if(include.all==TRUE){
    IDs.included <- IDs.to.include.in.NBDA

    
  } else if(include.all ==FALSE){
    IDs.included <- subset(ILVs.combined$PIT_f, ILVs.combined$Box!="no_box" & ILVs.combined$PIT_f %in% IDs.to.include.in.NBDA)


  }
  # we remove boxes D04, R06, G33 (females breeding in two boxes - we retaiend the ones closer to the dispenser)
  ILVs.sub.disp <- subset(ILVs.combined, ILVs.combined$PIT_f %in% IDs.included & !(ILVs.combined$Box %in% c("D04", "R06", "G33")))

  dispenser.data <- subset(dispenser.data, dispenser.data$PIT %in% IDs.to.include.in.NBDA)
  
  # order data ascending
  ILVs.sub.disp <- ILVs.sub.disp[order(ILVs.sub.disp$PIT_f),] # order ascending
  
  # subset the two networks to those IDs
  forage.net <- foraging.network.NBDA[rownames(foraging.network.NBDA) %in% IDs.included, colnames(foraging.network.NBDA) %in% IDs.included]
  neighbour.net <- neighbour_matrix.NBDA.new[rownames(neighbour_matrix.NBDA.new) %in% IDs.included, colnames(neighbour_matrix.NBDA.new) %in% IDs.included]
 # if(network=="forage"){
    assMatrix.nbda <- array(data=c(forage.net, neighbour.net), dim=c(nrow(forage.net), ncol(forage.net), 2))
    
#  } else if(network=="neighbour"){
 #   assMatrix.nbda <- array(data=c(neighbour.net), dim=c(nrow(forage.net), ncol(forage.net), 1))
    
#  }

  
  # create objects in the global environment for each ILV
  species.nbda <- as.matrix(ILVs.sub.disp$Species) 
  species.nbda[species.nbda!="GRETI"] <- -0.5 
  species.nbda[species.nbda=="GRETI"] <- 0.5 
  species.nbda <- as.matrix(as.numeric(species.nbda))
  
  age.nbda <- as.matrix(ILVs.sub.disp$Age) 
  age.nbda[age.nbda=="first.year"] <- -0.5 
  age.nbda[age.nbda=="adult"] <- 0.5 
  age.nbda <- as.matrix(as.numeric(age.nbda))
  
  distance.nbda <- as.matrix(as.numeric(ILVs.sub.disp[,location]))
  distance.nbda[is.na(distance.nbda)] <- 94.1 # for those not breeding in boxes, add the average distance
  distance.nbda <- as.matrix(as.numeric(scale(distance.nbda)))
  
  prev.exp.nbda <- ILVs.sub.disp$Prev.exp
  prev.exp.nbda[prev.exp.nbda=="yes"] <- 0.5
  prev.exp.nbda[prev.exp.nbda=="no"] <- -0.5
  prev.exp.nbda <- as.matrix(as.numeric(prev.exp.nbda))
  
  assign(paste("species", location, sep="_"), species.nbda, envir = .GlobalEnv)
  assign(paste("age", location, sep="_"), age.nbda, envir = .GlobalEnv)
  assign(paste("distance", location, sep="_"), distance.nbda, envir = .GlobalEnv)
  assign(paste("prev.exp", location, sep="_"), prev.exp.nbda, envir = .GlobalEnv)
  
 # ILVs <- paste(c("prev.exp"), location, sep="_")
#  ILVs <- paste(c("species", "age", "distance"), location, sep="_")
#  ILVs <- paste(c("species", "distance"), location, sep="_")
#  ILVs <- paste(c("species", "distance", "prev.exp"), location, sep="_")
#  ILVs <- paste(c("distance"), location, sep="_")
 # ILVs <- paste(c("species", "distance"), location, sep="_")
  ILVs <- paste(ILVs.include, location, sep="_")
#  ILVs <- paste(c("age"), location, sep="_")
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
  

  
  object2 <- nbdaData(label=location,                        
           assMatrix=object$assMatrix,          
           asoc_ilv=get(paste("ILVs", location, sep="_")),            
           int_ilv=get(paste("ILVs", location, sep="_")),            
           multi_ilv="ILVabsent",        
           orderAcq=object$OAc,          
           timeAcq=object$TAc,           
           endTime=41
  )
  
  return(object2)
}


# 6.4. Including all females ----------------------------------------------

 
# we first prepare NBDA data objects including all females (those breeding in boxes + the ones who have visited the dispenser) 
nbdaData_D1.all <- prepare.NBDA.data(dispenser.data = dispenser.data.1, include.all = TRUE, ILVs.include = c("age", "species", "distance", "prev.exp"))
nbdaData_D2.all <- prepare.NBDA.data(dispenser.data = dispenser.data.2, include.all = TRUE, ILVs.include = c("age", "species", "distance", "prev.exp"))
nbdaData_D3.all <- prepare.NBDA.data(dispenser.data = dispenser.data.3, include.all = TRUE, ILVs.include = c("age", "species", "distance", "prev.exp"))
# nbdaData_D4.all <- prepare.NBDA.data(dispenser.data = dispenser.data.4, include.all = TRUE, ILVs.include = c("age", "species", "distance", "prev.exp"))
nbdaData_D5.all <- prepare.NBDA.data(dispenser.data = dispenser.data.5, include.all = TRUE, ILVs.include = c("age", "species", "distance", "prev.exp"))

# number of learners:
# D1: 5
# D2: 7
# D3: 5
# D4: only one learner - function will not work
# D5: 4


# 6.5. Create constraints Vector MAtrix ---------------------------------


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

# it does not matter which NBDA data object we choose - the matrix is the same for all
constraintsVectMatrix <- create.constraints.Vect.Matrix(NBDA_data_object = nbdaData_D2.all, num_networks = 2, num_ILVs = 4)


# 6.6. Run TADA on all females  --------------------------------------------------------------------

TADA.finding.all <-
  tadaAICtable(
    nbdadata = list(
      nbdaData_D1.all,
      nbdaData_D2.all,
      nbdaData_D3.all,
      nbdaData_D5.all),
    constraintsVectMatrix = constraintsVectMatrix, 
    writeProgressFile = F
)
    
print(TADA.finding.all@printTable)
networksSupport(TADA.finding.all)
variableSupport(TADA.finding.all)

# we can see from the AIC table and the variable support that the models struggle to converge and the output is pretty much 'nonsense'
# it appears that certain combinations of ILVs cause troubles.
# We run TADA with each ILV alone to determine to then only include the ones that appear important for model fitting

# we prepare a new constraintsVectMatrix with only one ILV
constraintsVectMatrix <- create.constraints.Vect.Matrix(NBDA_data_object = nbdaData_D2.age, num_networks = 2, num_ILVs = 1)

# run TADA with just 'age'
# prepare NBDA data objects
nbdaData_D1.age <- prepare.NBDA.data(dispenser.data = dispenser.data.1, include.all = TRUE, ILVs.include = c("age"))
nbdaData_D2.age <- prepare.NBDA.data(dispenser.data = dispenser.data.2, include.all = TRUE, ILVs.include = c("age"))
nbdaData_D3.age <- prepare.NBDA.data(dispenser.data = dispenser.data.3, include.all = TRUE, ILVs.include = c("age"))
# nbdaData_D4.age <- prepare.NBDA.data(dispenser.data = dispenser.data.4, include.all = TRUE, ILVs.include = c("age"))
nbdaData_D5.age <- prepare.NBDA.data(dispenser.data = dispenser.data.5, include.all = TRUE, ILVs.include = c("age"))

# run TADA
TADA.finding.age <-
  tadaAICtable(
    nbdadata = list(
      nbdaData_D1.age,
      nbdaData_D2.age,
      nbdaData_D3.age,
      nbdaData_D5.age),
    constraintsVectMatrix = constraintsVectMatrix, 
    writeProgressFile = F
  )

print(TADA.finding.age@printTable)
variableSupport(TADA.finding.age)
# s1        s2 ASOC:age_D1 SOCIAL:age_D1
# support 0.9999818 0.1898923   0.1936806     0.2721888
# we have evidence for an influence of age on social or asocial learning (all support <0.5)


# run TADA with just 'species'
# prepare NBDA data objects
nbdaData_D1.species <- prepare.NBDA.data(dispenser.data = dispenser.data.1, include.all = TRUE, ILVs.include = c("species"))
nbdaData_D2.species <- prepare.NBDA.data(dispenser.data = dispenser.data.2, include.all = TRUE, ILVs.include = c("species"))
nbdaData_D3.species <- prepare.NBDA.data(dispenser.data = dispenser.data.3, include.all = TRUE, ILVs.include = c("species"))
# nbdaData_D4.species <- prepare.NBDA.data(dispenser.data = dispenser.data.4, include.all = TRUE, ILVs.include = c("species"))
nbdaData_D5.species <- prepare.NBDA.data(dispenser.data = dispenser.data.5, include.all = TRUE, ILVs.include = c("species"))

# run TADA
TADA.finding.species <-
  tadaAICtable(
    nbdadata = list(
      nbdaData_D1.species,
      nbdaData_D2.species,
      nbdaData_D3.species,
      nbdaData_D5.species),
    constraintsVectMatrix = constraintsVectMatrix, 
    writeProgressFile = F
  )

print(TADA.finding.species@printTable)
variableSupport(TADA.finding.species)
# s1        s2 ASOC:species_D1 SOCIAL:species_D1
# support 0.9999799 0.1917265       0.2046157         0.1923134
# we also have no evidence for species influencing social or asocial learning rate


# run TADA with just 'distance'
# prepare NBDA data objects
nbdaData_D1.distance <- prepare.NBDA.data(dispenser.data = dispenser.data.1, include.all = TRUE, ILVs.include = c("distance"))
nbdaData_D2.distance <- prepare.NBDA.data(dispenser.data = dispenser.data.2, include.all = TRUE, ILVs.include = c("distance"))
nbdaData_D3.distance <- prepare.NBDA.data(dispenser.data = dispenser.data.3, include.all = TRUE, ILVs.include = c("distance"))
# nbdaData_D4.distance <- prepare.NBDA.data(dispenser.data = dispenser.data.4, include.all = TRUE, ILVs.include = c("distance"))
nbdaData_D5.distance <- prepare.NBDA.data(dispenser.data = dispenser.data.5, include.all = TRUE, ILVs.include = c("distance"))

# run TADA
TADA.finding.distance <-
  tadaAICtable(
    nbdadata = list(
      nbdaData_D1.distance,
      nbdaData_D2.distance,
      nbdaData_D3.distance,
      nbdaData_D5.distance),
    constraintsVectMatrix = constraintsVectMatrix, 
    writeProgressFile = F
  )

print(TADA.finding.distance@printTable)
variableSupport(TADA.finding.distance)
# s1        s2 ASOC:distance_D1 SOCIAL:distance_D1
# support 0.9958509 0.1579205        0.9979706          0.6569978
# we support for distance influencing both social and asocial learning rate (both AIC > 0.5)

# run TADA with just 'prev.exp'
# prepare NBDA data objects
nbdaData_D1.prev.exp <- prepare.NBDA.data(dispenser.data = dispenser.data.1, include.all = TRUE, ILVs.include = c("prev.exp"))
nbdaData_D2.prev.exp <- prepare.NBDA.data(dispenser.data = dispenser.data.2, include.all = TRUE, ILVs.include = c("prev.exp"))
nbdaData_D3.prev.exp <- prepare.NBDA.data(dispenser.data = dispenser.data.3, include.all = TRUE, ILVs.include = c("prev.exp"))
# nbdaData_D4.prev.exp <- prepare.NBDA.data(dispenser.data = dispenser.data.4, include.all = TRUE, ILVs.include = c("prev.exp"))
nbdaData_D5.prev.exp <- prepare.NBDA.data(dispenser.data = dispenser.data.5, include.all = TRUE, ILVs.include = c("prev.exp"))

# run TADA
TADA.finding.prev.exp <-
  tadaAICtable(
    nbdadata = list(
      nbdaData_D1.prev.exp,
      nbdaData_D2.prev.exp,
      nbdaData_D3.prev.exp,
      nbdaData_D5.prev.exp),
    constraintsVectMatrix = constraintsVectMatrix, 
    writeProgressFile = F
  )

print(TADA.finding.prev.exp@printTable)
variableSupport(TADA.finding.prev.exp)
# s1        s2 ASOC:prev.exp_D1 SOCIAL:prev.exp_D1
# support 0.9999871 0.1848491         0.244346          0.4118598
# no evidence for previous experience influencing social or asocial learning rate


# run TADA with age, species and distance

nbdaData_D1.3.ILVs <- prepare.NBDA.data(dispenser.data = dispenser.data.1, include.all = TRUE, ILVs.include = c("age", "species", "distance"))
nbdaData_D2.3.ILVs <- prepare.NBDA.data(dispenser.data = dispenser.data.2, include.all = TRUE, ILVs.include = c("age", "species", "distance"))
nbdaData_D3.3.ILVs <- prepare.NBDA.data(dispenser.data = dispenser.data.3, include.all = TRUE, ILVs.include = c("age", "species", "distance"))
# nbdaData_D4.3.ILVs <- prepare.NBDA.data(dispenser.data = dispenser.data.4, include.all = TRUE, ILVs.include = c("age", "species", "distance"))
nbdaData_D5.3.ILVs <- prepare.NBDA.data(dispenser.data = dispenser.data.5, include.all = TRUE, ILVs.include = c("age", "species", "distance"))

constraintsVectMatrix <- create.constraints.Vect.Matrix(NBDA_data_object = nbdaData_D2.3.ILVs, num_networks = 2, num_ILVs = 3)

# run TADA
TADA.finding.3.ILVs <-
  tadaAICtable(
    nbdadata = list(
      nbdaData_D1.3.ILVs,
      nbdaData_D2.3.ILVs,
      nbdaData_D3.3.ILVs,
      nbdaData_D5.3.ILVs),
    constraintsVectMatrix = constraintsVectMatrix, 
    writeProgressFile = F
  )

print(TADA.finding.3.ILVs@printTable)
variableSupport(TADA.finding.3.ILVs)
# again, troubles converting

# we run it with just distance and species
nbdaData_D1.2.ILVs <- prepare.NBDA.data(dispenser.data = dispenser.data.1, include.all = TRUE, ILVs.include = c("species", "distance"))
nbdaData_D2.2.ILVs <- prepare.NBDA.data(dispenser.data = dispenser.data.2, include.all = TRUE, ILVs.include = c("species", "distance"))
nbdaData_D3.2.ILVs <- prepare.NBDA.data(dispenser.data = dispenser.data.3, include.all = TRUE, ILVs.include = c("species", "distance"))
# nbdaData_D4.2.ILVs <- prepare.NBDA.data(dispenser.data = dispenser.data.4, include.all = TRUE, ILVs.include = c("species", "distance"))
nbdaData_D5.2.ILVs <- prepare.NBDA.data(dispenser.data = dispenser.data.5, include.all = TRUE, ILVs.include = c("species", "distance"))

constraintsVectMatrix <- create.constraints.Vect.Matrix(NBDA_data_object = nbdaData_D2.2.ILVs, num_networks = 2, num_ILVs = 2)

# run TADA
TADA.finding.2.ILVs <-
  tadaAICtable(
    nbdadata = list(
      nbdaData_D1.2.ILVs,
      nbdaData_D2.2.ILVs,
      nbdaData_D3.2.ILVs,
      nbdaData_D5.2.ILVs),
    constraintsVectMatrix = constraintsVectMatrix, 
    writeProgressFile = F
  )

print(TADA.finding.2.ILVs@printTable)
networksSupport(TADA.finding.2.ILVs)
# support numberOfModels
# 0:0 0.005163674              4
# 0:1 0.001749869             16
# 1:0 0.848641864             16
# 1:2 0.144444593             16

# most evidential support for transmission along the foraging network (0.85), 
# followed by transmission through both the neighbour and the foraging network (0.15)
# asocial models and transmission through the neighbour network alone get very little support (<0.005)

variableSupport(TADA.finding.2.ILVs)
#             s1        s2       ASOC:species_D1   ASOC:distance_D1  SOCIAL:species_D1 SOCIAL:distance_D1
# support 0.9930865 0.1461945       0.2166048        0.9983375         0.1473031          0.7179743

# support for distance influencing both asocial (0.99) and social learning rate (0.63)
# no evidence for species influencing social or asocial learning rate (both <0.5)

# model averaged estimates
modelAverageEstimates(TADA.finding.2.ILVs , averageType = "median")
#    s1                 s2          ASOCIALspecies_D1  ASOCIALdistance_D1   SOCIALspecies_D1  SOCIALdistance_D1 
# 26.110377           0.000000           0.000000          -2.258964           0.000000          -1.472468 


# 6.7. Extract effect sizes ------------------------------------------------

# the best model is model 29 (top model in AIC table)
constraintsVectMatrix[29,]
# network network    asoc    asoc     int     int 
# 1       0       0       2       0       3 

# it contains the foraging network and distance influencing both asocial and social learning rate

# we create constrained NBDA Data Objects for that specific model
bestModelData1 <- constrainedNBDAdata(nbdadata=nbdaData_D1.2.ILVs,constraintsVect =constraintsVectMatrix[29,])
bestModelData2 <- constrainedNBDAdata(nbdadata=nbdaData_D2.2.ILVs,constraintsVect =constraintsVectMatrix[29,])
bestModelData3 <- constrainedNBDAdata(nbdadata=nbdaData_D3.2.ILVs,constraintsVect =constraintsVectMatrix[29,])
bestModelData5 <- constrainedNBDAdata(nbdadata=nbdaData_D5.2.ILVs,constraintsVect =constraintsVectMatrix[29,])


# and run TADA on the best model
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
cbind.data.frame(model.best.social@varNames, model.best.social@outputPar)

# model.best.social@varNames model.best.social@outputPar
# 1            Scale (1/rate):                 2960.292745
# 2    1 Social transmission 1                   26.110377
# 3     2 Asocial: distance_D1                   -2.258964
# 4      3 Social: distance_D1                   -1.472468


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
prop.solve.social # P=0.64696

# this means that 64.7% of birds have found the dispensers through social learning 
# the remaining 35.3% have done so through asocial learning

# extract profile likelihood. which=1 extracts the first parameter 
# (in this case s for the foraging network)
plotProfLik(which=1,model=model.best.social,range=c(0,800), resolution=10) 
# we zoom into the lower and upper range (where the profile likelihood crosses the dotted line)
plotProfLik(which=1,model=model.best.social,range=c(0,10), resolution=10) 
plotProfLik(which=1,model=model.best.social,range=c(600,650), resolution=10) 
# we check where the profile likelihood crosses the dotted line to get the
# range for the lower and upper interval - set the ranges accordingly
CIs <- profLikCI(which=1,model=model.best.social, lowerRange = c(0,10), upperRange = c(600, 650)) # extract confidence intervals
CIs
# Lower CI   Upper CI 
# 1.925461 622.794178

######### lower and upper bound in %

#To get the estimates for the lower bound we should really find the corresponding value of the other parameters to plug in when s1 is constrained to this value
bestModelDataS1LowerBound.D1 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D1.2.ILVs,
  constraintsVect = constraintsVectMatrix[29, ],
  offset = c(CIs[1] , rep(0, 5))
)

bestModelDataS1LowerBound.D2 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D2.2.ILVs,
  constraintsVect = constraintsVectMatrix[29, ],
  offset = c(CIs[1] , rep(0, 5))
)

bestModelDataS1LowerBound.D3 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D3.2.ILVs,
  constraintsVect = constraintsVectMatrix[29, ],
  offset = c(CIs[1] , rep(0, 5))
)

bestModelDataS1LowerBound.D5 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D5.2.ILVs,
  constraintsVect = constraintsVectMatrix[29, ],
  offset = c(CIs[1] , rep(0, 5))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then the value of s at the lower bound is added to s as an offset
bestModelS1LowerBound <-
  tadaFit(
    list(
      bestModelDataS1LowerBound.D1,
      bestModelDataS1LowerBound.D2,
      bestModelDataS1LowerBound.D3,
#      bestModelDataS1LowerBound.D4,
      bestModelDataS1LowerBound.D5
    ) ,
    type = "asocial"
  )
bestModelS1LowerBound@outputPar
# [1] 1158.293029   -1.602609   -2.754137
#Now plug into the prop solve function
prop.solve.social.lower <-
  oadaPropSolveByST(
    model = bestModelS1LowerBound,
    nbdadata = list(
      bestModelDataS1LowerBound.D1,
      bestModelDataS1LowerBound.D2,
      bestModelDataS1LowerBound.D3,
      bestModelDataS1LowerBound.D4,
      bestModelDataS1LowerBound.D5
    )
  )
prop.solve.social.lower
# P(S offset) 
# 0.50618
# lower bound for % of birds having learned the dial task through social learning is 50.6%

#To get the estimates for the upper bound we should really find the corresponding value of the other parameters to plug in when s1 is constrained to this value
bestModelDataS1upperBound.D1 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D1.2.ILVs,
  constraintsVect = constraintsVectMatrix[29, ],
  offset = c(CIs[2] , rep(0, 5))
)

bestModelDataS1upperBound.D2 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D2.2.ILVs,
  constraintsVect = constraintsVectMatrix[29, ],
  offset = c(CIs[2] , rep(0, 5))
)

bestModelDataS1upperBound.D3 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D3.2.ILVs,
  constraintsVect = constraintsVectMatrix[29, ],
  offset = c(CIs[2] , rep(0, 5))
)

bestModelDataS1upperBound.D5 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D5.2.ILVs,
  constraintsVect = constraintsVectMatrix[29, ],
  offset = c(CIs[2] , rep(0, 5))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then the value of s at the upper bound is added to s as an offset
bestModelS1upperBound <-
  tadaFit(
    list(
      bestModelDataS1upperBound.D1,
      bestModelDataS1upperBound.D2,
      bestModelDataS1upperBound.D3,
      #      bestModelDataS1upperBound.D4,
      bestModelDataS1upperBound.D5
    ) ,
    type = "asocial"
  )
bestModelS1upperBound@outputPar
# [1] 42733.558644    -4.364275    -1.127624
#Now plug into the prop solve function
prop.solve.social.upper <-
  oadaPropSolveByST(
    model = bestModelS1upperBound,
    nbdadata = list(
      bestModelDataS1upperBound.D1,
      bestModelDataS1upperBound.D2,
      bestModelDataS1upperBound.D3,
      bestModelDataS1upperBound.D4,
      bestModelDataS1upperBound.D5
    )
  )
prop.solve.social.upper
# P(S offset) 
# 0.74311
# upper bound for % of birds having learned the dial task through social learning is 74.3%



# 7) Wool choice ----------------------------------------------------------


wool.choice <- read.delim("C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/Wool_boxes.txt", sep="\t")

start.colours <- rbind.data.frame(cbind("D1", "Pi"),
cbind("D2", "Pu"),
cbind("D3", "O"),
cbind("D4", "B"),
cbind("D5", "Pi"))

colnames(start.colours) <- c("dispenser", "start.colour")


# add the closest dispenser to all ILVs

# extract which dispenser was closest 
for(i in wool.choice$Box){
  closest <- subset(closest_disp$`Closest dispenser`, rownames(closest_disp)==i)
  wool.choice[which(wool.choice$Box==i), "dispenser"] <- closest
}


# extract those that accessed the dispenser before the second colour was made available 
demos <- c("B07",
           "C08",
           "D03",
           "H08",
           "H24",
           "S16",
           "T01",
           "V02")


wool.choice.learners <- subset(wool.choice, !(wool.choice$Box%in% demos))



for(i in wool.choice.learners$dispenser){
  col <- subset(start.colours$start.colour, start.colours$dispenser==i)
  wool.choice.learners[which(wool.choice.learners$dispenser==i), "start.col"] <- col
}


# Fishers test

fisher <- fisher.test(wool.choice.learners$first_color, wool.choice.learners$start.col, alternative = "greater")
fisher

# Fisher's Exact Test for Count Data
# 
# data:  wool.choice.learners$first_color and wool.choice.learners$start.col
# p-value = 0.02498
# alternative hypothesis: greater


# 8) Visualization --------------------------------------------------------

# create a network with polgygons around dispenser areas

library(igraph)

# for each PIT tag, extract whether it was seen at a network feeder
col.plot <- NULL

for(i in rownames(foraging.network.NBDA)){
  box <- subset(ILVs.combined$Box, ILVs.combined$PIT_f==i)
  if(length(box)>1){
    box <- subset(box, box %in% wool.choice$Box)
  } 
  if(length(box)==0){
    col <- "none"
    box <- "no"
  }
    if(box=="no_box"){
    wool.used <- "unknown"
  } else {
    wool.used <- subset(wool.choice$first_color, wool.choice$Box==box)  
    
  }
  if(length(wool.used)==0){
    wool.used <- "none"
  }
    col.plot[which(rownames(foraging.network.NBDA)==i)] <- wool.used
}

col.plot[col.plot=="unknown"] <- "#575656"
col.plot[col.plot=="none"] <- "#faf7f7"
col.plot[col.plot=="Pi"] <- "#F056B3"
col.plot[col.plot=="O"] <- "#F49527"
col.plot[col.plot=="Pu"] <- "#873EE0"
col.plot[col.plot=="B"] <- "#33CCFF"


g.net <- graph_from_adjacency_matrix(foraging.network.NBDA, mode = "undirected",
                                     weighted = TRUE, diag = FALSE)

E(g.net)$width <- E(g.net)$weight

V(g.net)$colour <- col.plot
l <- layout_with_kk(g.net)
l <- layout_with_fr(g.net )
l <- layout_with_gem(g.net)
# for each individual, extract which dispenser was closest
# for those not breeding in boxes, we assign the dispenser they visited
all.D <- NULL
suppressWarnings(
  for(i in rownames(foraging.network.NBDA)){
    box <- subset(ILVs.combined$Box, ILVs.combined$PIT_f==i)
    if(box=="no_box"){
      sub <- subset(ILVs.combined, ILVs.combined$PIT_f==i)
      D <- c("D1", "D2", "D3", "D4", "D5")[which(sub[,c("D1.visited", "D2.visited", "D3.visited", "D4.visited", "D5.visited"),]==1)]
    } else {
      D <- unique(subset(closest_disp$`Closest dispenser`, rownames(closest_disp) %in% box))
    }
    D <- unique(D)
    all.D[which(rownames(foraging.network.NBDA)==i)] <- D  
  }
  
)

length(all.D)
set.seed(12)

 
list(c(14, 24, 27, 28, 36, 41, 43), 
     c(1,  2,  6,  7, 10, 15, 16, 17, 18, 21, 22, 23, 25, 26, 33, 37, 39, 42, 46),
     c(3,  9, 11, 20, 29, 30, 32, 34, 38, 40, 45),
     c(4, 12, 13, 35, 4),
     c( 5,  8, 19, 31))

# set transparent polgygons
col.adj <- grDevices::adjustcolor(c("#a08f00","#62dab9", "#b738bd", "#7a3f63", "#c31910"), alpha=0.15)










png( "network.png", units="in", width=12, height=4, res=400)


igraph::plot.igraph( g.net,
      vertex.size = 6,
      edge.curved = 0.2,
      edge.color =  "#8c8989",
      vertex.color = V(g.net)$colour,
      vertex.label = NA,
      vertex.frame.colour = "black",
      edge.width = E(g.net)$width*5,
      frame = FALSE,
      layout=l,
      asp = 1,
  #    rescale = TRUE,
      mark.groups = list(c(14, 24, 27, 28, 36, 41, 43), 
                         c(1,  2,  6,  7, 10, 15, 16, 17, 18, 21, 22, 23, 25, 26, 33, 37, 39, 42, 46),
                         c(3,  9, 11, 20, 29, 30, 32, 34, 38, 40, 45),
                         c(4, 12, 13, 35, 4),
                         c( 5,  8, 19, 31)),
      mark.border =NA,
    mark.col=col.adj
      
      
)







#abd2c1
#e6b8b3
#98d4e4
#d4d7b2
#c4bedf







dev.off()











