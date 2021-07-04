
# ## LINING MATERIAL ------------------------------------------------------

# 1) Load libraries -------------------------------------------------------

library(asnipe)
library(sp)
library(geosphere)
library(Imap)



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


# 3.5. Calculate the neighbour matrix -------------------------------------

# we set all values above 50m to 0
# we use the inverted squared distance

neighbour_matrix <- distance_matrix

neighbour_matrix[neighbour_matrix>=50] <- 0

neighbour_matrix <- 1/neighbour_matrix

# set all infinity values to 0
neighbour_matrix[neighbour_matrix==Inf] <- 0
hist(neighbour_matrix[neighbour_matrix!=0])


# 4) Prepare ILVs ---------------------------------------------------------

breeders <- read.delim("C:/Users/swild/Desktop/Konstanz/Breeding Season 2021/Breeders 2021.txt")

ILVs <- as.data.frame(matrix(NA, nrow=length(GPS_data[,1])-5, ncol=5))
ILVs[,1] <- GPS_data$Box[1:(length(GPS_data$Box)-5)] # remove the Dispensers (last 5)

colnames(ILVs) <- c("Box", "PIT_fem", "Species", "Age", "NA")

# 4.1. Add the PIT tags of the females and the species ------------------------------------


for(i in ILVs$Box){
PIT.fem <- subset(breeders$Tag, breeders$Box==i & breeders$Who=="Female")
species <- subset(breeders$Species, breeders$Box==i & breeders$Who=="Female")
if(length(PIT.fem)==0){
  PIT.fem <- NA
  species <- NA
}

ILVs[which(ILVs$Box==i), "PIT_fem"] <- PIT.fem
ILVs[which(ILVs$Box==i), "Species"] <- species
}

# remove the boxes with no breeders
ILVs <- subset(ILVs, !(is.na(ILVs$PIT_fem)))

# 4.2. Add age ------------------------------------------------------------

age <- read.delim("C:/Users/swild/Desktop/Konstanz/Lining material - social learning/Raw data/age.data.spring.2021.txt", sep=",", row.names = 1)

for(i in ILVs$PIT_fem){
    ILVs[which(ILVs$PIT_fem==i), "Age"] <- unique(subset(age$Age.in.2021, age$Pit==i)) 
}


