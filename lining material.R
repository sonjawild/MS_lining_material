
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
for(i in ILVs.combined$Box[ILVs.combined$Box!="no_box"]){
  dist <- subset(closest_disp$`Distance to dispenser`, rownames(closest_disp)==i)
  ILVs.combined[which(ILVs.combined$Box==i), "closest.dispenser"] <- dist
  
}
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
#       Significance: 0.0709 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0425 0.0617 0.0778 0.0974 
# Permutation: free
# Number of permutations: 9999

# ILVs.combined <- subset(ILVs.combined, ILVs.combined$Box!="no_box")


# 6.3. Prepare ILVs ------------------------------------------------------------------

prepare.NBDA.data <- function(dispenser.data, include.all){
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
    ILVs.sub.disp <- ILVs.combined
    
  } else if(include.all ==FALSE){
    IDs.included <- subset(ILVs.combined$PIT_f, ILVs.combined$Box!="no_box" & ILVs.combined$PIT_f %in% IDs.to.include.in.NBDA)
    ILVs.sub.disp <- subset(ILVs.combined, ILVs.combined$PIT_f %in% IDs.included)

  }
  

  
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
#  assign(paste("prev.exp", location, sep="_"), prev.exp.nbda, envir = .GlobalEnv)
  
#  ILVs <- paste(c("prev.exp"), location, sep="_")
  ILVs <- paste(c("species", "age", "distance"), location, sep="_")
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


# 6.4. Including all females ----------------------------------------------

 
# we first prepare NBDA data objects including all females (those breeding in boxes + the ones who have visited the dispenser) 
d1.NBDA.all <- prepare.NBDA.data(dispenser.data = dispenser.data.1, include.all = TRUE)
d2.NBDA.all <- prepare.NBDA.data(dispenser.data = dispenser.data.2, include.all = TRUE)
d3.NBDA.all <- prepare.NBDA.data(dispenser.data = dispenser.data.3, include.all = TRUE)
d4.NBDA.all <- prepare.NBDA.data(dispenser.data = dispenser.data.4, include.all = TRUE)
d5.NBDA.all <- prepare.NBDA.data(dispenser.data = dispenser.data.5, include.all = TRUE)

# number of learners:
# D1: 5
# D2: 7
# D3: 5
# D4: 1
# D5: 4

nbdaData_D1.all <- nbdaData(label="D1",                        # specify an informative label
                            assMatrix=d1.NBDA.all$assMatrix,           # our array with the matrices
                            asoc_ilv=get(paste("ILVs", "D1", sep="_")),            # we specify that ILVs can influence asocial learning, if no ILV: "ILVabsent"
                            int_ilv=get(paste("ILVs", "D1", sep="_")),             # we specify that our ILVs can influence social learning, if no ILV: "ILVabsent" 
                            multi_ilv="ILVabsent",        # we specify that our ILVs can influence asocial and social learning to the same extent, if no ILV: "ILVabsent" 
                            orderAcq=d1.NBDA.all$OAc,          # vector with the order of acquisition 
                            timeAcq=d1.NBDA.all$TAc,           # numerical vector giving the time at which each individual acquired the target behaviour, given in the order matching orderAcqv
                            endTime=41                    # numeric giving the time at which the diffusion ended. (11th of May = day 40)
                            )

nbdaData_D2.all <- nbdaData(label="D2",                        
                        assMatrix=d2.NBDA.all$assMatrix,          
                        asoc_ilv=get(paste("ILVs", "D2", sep="_")),            
                        int_ilv=get(paste("ILVs", "D2", sep="_")),            
                        multi_ilv="ILVabsent",        
                        orderAcq=d2.NBDA.all$OAc,          
                        timeAcq=d2.NBDA.all$TAc,           
                        endTime=41
)


nbdaData_D3.all <- nbdaData(label="D3",                        
                        assMatrix=d3.NBDA.all$assMatrix,          
                        asoc_ilv=get(paste("ILVs", "D3", sep="_")),            
                        int_ilv=get(paste("ILVs", "D3", sep="_")),            
                        multi_ilv="ILVabsent",        
                        orderAcq=d3.NBDA.all$OAc,          
                        timeAcq=d3.NBDA.all$TAc,           
                        endTime=41
)



# skip 4

nbdaData_D5.all <- nbdaData(label="D5",                        
                        assMatrix=d5.NBDA.all$assMatrix,          
                        asoc_ilv=get(paste("ILVs", "D5", sep="_")),            
                        int_ilv=get(paste("ILVs", "D5", sep="_")),            
                        multi_ilv="ILVabsent",        
                        orderAcq=d5.NBDA.all$OAc,          
                        timeAcq=d5.NBDA.all$TAc,           
                        endTime=41
)

# create function for making constraintsVectorMatrix



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
constraintsVectMatrix <- create.constraints.Vect.Matrix(NBDA_data_object = nbdaData_D2.all, num_networks = 2, num_ILVs = 3)


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
    

# view the results with 'print'
print(TADA.finding.all@printTable) # with all females, 3 ILVs

# network support
networksSupport <- networksSupport(TADA.finding.all)
round(networksSupport, 2)

# variable support
variable_support <- variableSupport(TADA.finding.all, includeAsocial = T)
round(variable_support,3)

# model averaged estimates
MLE_med  <- modelAverageEstimates(TADA.finding.all , averageType = "median")
round(MLE_med,2)


# 6.7. Extract effect sizes ------------------------------------------------

bestModelData1 <- constrainedNBDAdata(nbdadata=nbdaData_D1.all,constraintsVect =constraintsVectMatrix[160,])
bestModelData2 <- constrainedNBDAdata(nbdadata=nbdaData_D2.all,constraintsVect =constraintsVectMatrix[160,])
bestModelData3 <- constrainedNBDAdata(nbdadata=nbdaData_D3.all,constraintsVect =constraintsVectMatrix[160,])
bestModelData5 <- constrainedNBDAdata(nbdadata=nbdaData_D5.all,constraintsVect =constraintsVectMatrix[160,])


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
# [1] "Scale (1/rate):"         "1 Social transmission 1" "2 Asocial: distance_D1"  "3 Social: distance_D1" 
# [1] 823.7992383   5.7684467   0.7607148  -1.8509040

# best model is showing 'relative convergence'
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
prop.solve.social # P=0.55

# extract profile likelihood. which=1 extracts the first parameter (in this case s for the vertical network)
# (in this case the s parameter for vertical social learning)
plotProfLik(which=1,model=model.best.social,range=c(0,100), resolution=10) # start with large range
CIs <- profLikCI(which=1,model=model.best.social, lowerRange = c(0,5), upperRange = c(40, 60)) # extract confidence intervals
CIs
#   Lower CI   Upper CI 
# 0.5396941 45.9442201  

######### lower and upper bound in %

#To get the estimates for the lower bound we should really find the corresponding value of the other parameters to plug in when s1 is constrained to this value
bestModelDataS1LowerBound.D1 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D1.all,
  constraintsVect = constraintsVectMatrix[160, ],
  offset = c(CIs[1] , rep(0, 7))
)

bestModelDataS1LowerBound.D2 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D2.all,
  constraintsVect = constraintsVectMatrix[160, ],
  offset = c(CIs[1] , rep(0, 7))
)

bestModelDataS1LowerBound.D3 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D3.all,
  constraintsVect = constraintsVectMatrix[160, ],
  offset = c(CIs[1] , rep(0, 7))
)

bestModelDataS1LowerBound.D5 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D5.all,
  constraintsVect = constraintsVectMatrix[160, ],
  offset = c(CIs[1] , rep(0, 7))
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
# [1] 591.956573   0.604114  -3.265292
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
# lower bound for % of birds having learned the dial task through social learning is 48.1%

#To get the estimates for the upper bound we should really find the corresponding value of the other parameters to plug in when s1 is constrained to this value
bestModelDataS1upperBound.D1 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D1.all,
  constraintsVect = constraintsVectMatrix[160, ],
  offset = c(CIs[2] , rep(0, 7))
)

bestModelDataS1upperBound.D2 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D2.all,
  constraintsVect = constraintsVectMatrix[160, ],
  offset = c(CIs[2] , rep(0, 7))
)

bestModelDataS1upperBound.D3 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D3.all,
  constraintsVect = constraintsVectMatrix[160, ],
  offset = c(CIs[2] , rep(0, 7))
)

bestModelDataS1upperBound.D5 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D5.all,
  constraintsVect = constraintsVectMatrix[160, ],
  offset = c(CIs[2] , rep(0, 7))
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
# [1] 1787.6480094    1.0372553   -0.7438817
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
# upper bound for % of birds having learned the dial task through social learning is 69.8%

# 6.8. Run TADA on females in nest boxes ----------------------------------

# we first prepare NBDA data objects including all females (those breeding in boxes + the ones who have visited the dispenser) 
d1.NBDA.sub <- prepare.NBDA.data(dispenser.data = dispenser.data.1, include.all = FALSE)
d2.NBDA.sub <- prepare.NBDA.data(dispenser.data = dispenser.data.2, include.all = FALSE)
d3.NBDA.sub <- prepare.NBDA.data(dispenser.data = dispenser.data.3, include.all = FALSE)
d4.NBDA.sub <- prepare.NBDA.data(dispenser.data = dispenser.data.4, include.all = FALSE)
d5.NBDA.sub <- prepare.NBDA.data(dispenser.data = dispenser.data.5, include.all = FALSE)

# number of learners:
# D1: 1
# D2: 3
# D3: 3
# D4: 0
# D5: 4

# nbdaData_D1.sub <- nbdaData(label="D1",                        # specify an informative label
#                             assMatrix=d1.NBDA.sub$assMatrix,           # our array with the matrices
#                             asoc_ilv=get(paste("ILVs", "D1", sep="_")),            # we specify that ILVs can influence asocial learning, if no ILV: "ILVabsent"
#                             int_ilv=get(paste("ILVs", "D1", sep="_")),             # we specify that our ILVs can influence social learning, if no ILV: "ILVabsent" 
#                             multi_ilv="ILVabsent",        # we specify that our ILVs can influence asocial and social learning to the same extent, if no ILV: "ILVabsent" 
#                             orderAcq=d1.NBDA.sub$OAc,          # vector with the order of acquisition 
#                             timeAcq=d1.NBDA.sub$TAc,           # numerical vector giving the time at which each individual acquired the target behaviour, given in the order matching orderAcqv
#                             endTime=41                    # numeric giving the time at which the diffusion ended. (11th of May = day 40)
# )

nbdaData_D2.sub <- nbdaData(label="D2",                        
                            assMatrix=d2.NBDA.sub$assMatrix,          
                            asoc_ilv=get(paste("ILVs", "D2", sep="_")),            
                            int_ilv=get(paste("ILVs", "D2", sep="_")),            
                            multi_ilv="ILVabsent",        
                            orderAcq=d2.NBDA.sub$OAc,          
                            timeAcq=d2.NBDA.sub$TAc,           
                            endTime=41
)


nbdaData_D3.sub <- nbdaData(label="D3",                        
                            assMatrix=d3.NBDA.sub$assMatrix,          
                            asoc_ilv=get(paste("ILVs", "D3", sep="_")),            
                            int_ilv=get(paste("ILVs", "D3", sep="_")),            
                            multi_ilv="ILVabsent",        
                            orderAcq=d3.NBDA.sub$OAc,          
                            timeAcq=d3.NBDA.sub$TAc,           
                            endTime=41
)



# skip 4

nbdaData_D5.sub <- nbdaData(label="D5",                        
                            assMatrix=d5.NBDA.sub$assMatrix,          
                            asoc_ilv=get(paste("ILVs", "D5", sep="_")),            
                            int_ilv=get(paste("ILVs", "D5", sep="_")),            
                            multi_ilv="ILVabsent",        
                            orderAcq=d5.NBDA.sub$OAc,          
                            timeAcq=d5.NBDA.sub$TAc,           
                            endTime=41
)






TADA.finding.sub <-
  tadaAICtable(
    nbdadata = list(
    #  nbdaData_D1.sub,
      nbdaData_D2.sub,
      nbdaData_D3.sub,
      nbdaData_D5.sub),
    constraintsVectMatrix = constraintsVectMatrix, 
    writeProgressFile = F
  )


# view the results with 'print'
print(TADA.finding.sub@printTable) # with sub females, 3 ILVs

# network support
networksSupport <- networksSupport(TADA.finding.sub)
round(networksSupport, 2)

# variable support
variable_support <- variableSupport(TADA.finding.sub, includeAsocial = T)
round(variable_support,3)

# model averaged estimates
MLE_med  <- modelAverageEstimates(TADA.finding.sub , averageType = "median")
round(MLE_med,2)


# 6.7. Extract effect sizes ------------------------------------------------

bestModelData1 <- constrainedNBDAdata(nbdadata=nbdaData_D1.all,constraintsVect =constraintsVectMatrix[192,])
bestModelData2 <- constrainedNBDAdata(nbdadata=nbdaData_D2.all,constraintsVect =constraintsVectMatrix[192,])
bestModelData3 <- constrainedNBDAdata(nbdadata=nbdaData_D3.all,constraintsVect =constraintsVectMatrix[192,])
bestModelData5 <- constrainedNBDAdata(nbdadata=nbdaData_D5.all,constraintsVect =constraintsVectMatrix[192,])


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
# [1] "Scale (1/rate):"         "1 Social transmission 1"
# [1]   125.065973   2.124156

# best model is showing 'relative convergence'
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
prop.solve.social # P=0.34

# extract profile likelihood. which=1 extracts the first parameter (in this case s for the vertical network)
# (in this case the s parameter for vertical social learning)
plotProfLik(which=1,model=model.best.social,range=c(0,20), resolution=10) # start with large range
CIs <- profLikCI(which=1,model=model.best.social, upperRange = c(7,11)) # extract confidence intervals
CIs
# Lower CI Upper CI 
# 0.00000  9.00025  

######### lower and upper bound in %

#To get the estimates for the lower bound we should really find the corresponding value of the other parameters to plug in when s1 is constrained to this value
bestModelDataS1LowerBound.D1 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D1.all,
  constraintsVect = constraintsVectMatrix[192, ],
  offset = c(CIs[1] , rep(0, 7))
)

bestModelDataS1LowerBound.D2 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D2.all,
  constraintsVect = constraintsVectMatrix[192, ],
  offset = c(CIs[1] , rep(0, 7))
)

bestModelDataS1LowerBound.D3 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D3.all,
  constraintsVect = constraintsVectMatrix[192, ],
  offset = c(CIs[1] , rep(0, 7))
)

bestModelDataS1LowerBound.D5 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D5.all,
  constraintsVect = constraintsVectMatrix[192, ],
  offset = c(CIs[1] , rep(0, 7))
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
# [1] 90.27884
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
# lower bound for % of birds having learned the dial task through social learning is 0%

#To get the estimates for the upper bound we should really find the corresponding value of the other parameters to plug in when s1 is constrained to this value
bestModelDataS1upperBound.D1 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D1.all,
  constraintsVect = constraintsVectMatrix[192, ],
  offset = c(CIs[2] , rep(0, 7))
)

bestModelDataS1upperBound.D2 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D2.all,
  constraintsVect = constraintsVectMatrix[192, ],
  offset = c(CIs[2] , rep(0, 7))
)

bestModelDataS1upperBound.D3 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D3.all,
  constraintsVect = constraintsVectMatrix[192, ],
  offset = c(CIs[2] , rep(0, 7))
)

bestModelDataS1upperBound.D5 <- constrainedNBDAdata(
  nbdadata =
    nbdaData_D5.all,
  constraintsVect = constraintsVectMatrix[192, ],
  offset = c(CIs[2] , rep(0, 7))
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
# [1] 237.6753
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
# upper bound for % of birds having learned the dial task through social learning is 59.4%







# 7) Wool choice ----------------------------------------------------------


wool.choice <- read.delim("C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/Wool_boxes.txt", sep="\t")

start.colours <- rbind.data.frame(cbind("D1", "Pi"),
cbind("D2", "Pu"),
cbind("D3", "O"),
cbind("D4", "B"),
cbind("D5", "Pi"))

colnames(start.colours) <- c("dispenser", "start.colour")

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

# extract which dispenser was closest 
for(i in wool.choice.learners$Box){
  closest.disp <- subset(ILVs$Closest.dispenser, ILVs$Box==i)
  wool.choice.learners[which(wool.choice.learners$Box==i), "dispenser"] <- closest.disp 
}

for(i in wool.choice.learners$dispenser){
  col <- subset(start.colours$start.colour, start.colours$dispenser==i)
  wool.choice.learners[which(wool.choice.learners$dispenser==i), "start.col"] <- col
}


# Fishers test
fisher <- fisher.test(table(learners$provided_color, learners$first_col), alternative = "greater")


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
which.disp.plot <- NULL

for(i in rownames(foraging.network.NBDA)){
  disp <- unique(subset(ILVs.combined$Closest.dispenser, ILVs.combined$PIT_f==i))
  which.disp.plot[which(rownames(foraging.network.NBDA)==i)] <- disp
}

which.disp.plot[which.disp.plot=="D1"] <- "#7A7978"
which.disp.plot[which.disp.plot=="D2"] <- "#87CBAC"
which.disp.plot[which.disp.plot=="D3"] <- "#90FFDC"
which.disp.plot[which.disp.plot=="D4"] <- "#8DE4FF"
which.disp.plot[which.disp.plot=="D5"] <- "#8AC4FF"

length(rownames(foraging.network.NBDA))

g.net <- graph_from_adjacency_matrix(foraging.network.NBDA, mode = "undirected",
                                     weighted = TRUE, diag = FALSE)

E(g.net)$width <- E(g.net)$weight

V(g.net)$colour <- which.disp.plot
l <- layout_with_fr(g.net )
l <- layout_in_circle(g.net )

plot( g.net,
      vertex.size = 5,
      edge.curved = 0.2,
      edge.color =  "darkgrey",
      vertex.color = V(g.net)$coSlour,
      vertex.label = NA,
      vertex.frame.colour = "black",
      edge.width = E(g.net)$width*10,
      frame = TRUE,
      layout=l,
      asp = 1,
      rescale = TRUE
)



plot.network <- function(net, demos, learners.list, experiment, n, solves, diffusion, position, simple.list){
  # create a colour vector that distinguishes between learners.list, naive birds and demonstrators
  cat.all <- NULL
  for (i in rownames(net)){
    if(i %in% demos){
      cat <- "#F0DC6A"
    } else if(i %in% simple.list){
      cat <- "#33CCCC"  # blue for simple learners
    } else if (i %in% learners.list$RING){
      cat <- "#660099"  # purple for complex learners
    } 
    else {cat <- "#FFFFFF"} #white for non-learners
    cat.all[which(rownames(net)==i)] <- cat    
  }
  
  # create graph object
  g.net <- graph_from_adjacency_matrix(net, mode = "undirected",
                                       weighted = TRUE, diag = FALSE)
  V(g.net)$colour <- cat.all
  
  E(g.net)$width <- E(g.net)$weight
  E(g.net)$edge.color <- "gray50"
  # remove edge weights below 0.005
  g.net <- delete_edges(g.net, E(g.net)[weight<0.05])
  # layout fruchtermann reingold
  l <- layout_with_fr(g.net )
  title <-  paste(experiment, paste("(N=", n, ")", sep =
                                      ""), sep = " ")
  
  
  plot( g.net,
        vertex.size = 6,
        edge.curved = 0.2,
        edge.color =  E(g.net)$edge.color,
        vertex.color = V(g.net)$colour,
        vertex.label = NA,
        vertex.frame.colour = "black",
        edge.width = E(g.net)$width*10,
        frame = TRUE,
        layout=l,
        #   main= title,
        adj = c(0,-1),
        margin=c(0.02,0.0,0.0,0.0),
        asp = 1,
        rescale = TRUE
  )
  
  title(title, line=0.5, adj=0.0, cex.main=1.5, font =1, family = "sans")
  if(position!="none"){
    legend(x=-1.15, y=-0.7, c("naive","simple", "demonstrator", "complex"), pch=21,
           
           col="#777777", pt.bg=c("#FFFFFF", "#33CCCC", "#F0DC6A", "#660099"), pt.cex=1.5, cex=1.2, bty="n", ncol=1, 
           y.intersp=0.8)
  }
  
  
}
set.seed(10)

png( "networks.png", units="in", width=12, height=4, res=400)
par(mfrow=c(1,3), mai = c(0.0, 0.0, 0.2, 0.07))


p1 <-
  plot.network(
    net = dial.net,
    demos = demos.dial,
    learners.list = learners.list.dial,
    n = BB_NBDA_DATA$num.birds,
    experiment = "A) Dial diffusion",
    solves = "dial",
    diffusion = "dial_diffusion",
    position = "left",
    simple.list = learners.list.dial$RING
  )


p2 <-
  plot.network(
    net = complex1st.net,
    demos = NA,
    learners.list = learners.list.complex1st,
    n = BB_complex_progressive$num.birds,
    experiment = "B) Complex 1st generation",
    solves = "complex",
    diffusion = "complex_prog",
    position = "none", 
    simple.list = learners.complex1st.simple
  )

p3 <-
  plot.network(
    net = complex2nd.net,
    demos = demos.complex2nd,
    learners.list = learners.list.complex2nd,
    n = BB_complex_next_gen$num.birds,
    experiment = "C) Complex 2nd generation",
    solves = "complex",
    diffusion = "nextgen",
    position = "none",
    simple.list = learners.complex2nd.simple
  )


dev.off()











