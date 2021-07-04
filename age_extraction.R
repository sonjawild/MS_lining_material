## extracts hatch date - but still has an error for the 2020 chicks


age.data <- read.delim("C:/Users/swild/Desktop/Konstanz/Restricted puzzle box/Age extraction/age data.txt", sep="\t", header=TRUE)


age.cat <- sort(as.vector(unique(age.data$Age)))
age.data$Date <- as.Date(age.data$Date, '%d/%m/%Y')

age.data$PIT.year <- as.numeric(format(age.data$Date,'%Y'))

age.data$hatch.year <- NA

age.data <- subset(age.data, subset=!is.na(age.data$Age))


for (i in 1:length(age.data[, 1])) {
  if (age.data[i, "Age"] == 1) {
    age.data[i, "hatch.year"] <- age.data[i, "PIT.year"]
  } else if (age.data[i, "Age"] == 2) {
    age.data[i, "hatch.year"] <- age.data[i, "PIT.year"]
  } else if (age.data[i, "Age"] == 3) {
    age.data[i, "hatch.year"] <- age.data[i, "PIT.year"]
  } else if (age.data[i, "Age"] == 4) {
    age.data[i, "hatch.year"] <- age.data[i, "PIT.year"] - 1
  } else if (age.data[i, "Age"] == 5) {
    age.data[i, "hatch.year"] <- age.data[i, "PIT.year"] - 1
  } else if (age.data[i, "Age"] == 6) {
    age.data[i, "hatch.year"] <- age.data[i, "PIT.year"] - 2
  } else if (age.data[i, "Age"] == 7) {
    age.data[i, "hatch.year"] <- age.data[i, "PIT.year"] - 2
  } else if (age.data[i, "Age"] == 8) {
    age.data[i, "hatch.year"] <- age.data[i, "PIT.year"] - 3
  } else if (age.data[i, "Age"] == 9) {
    age.data[i, "hatch.year"] <- age.data[i, "PIT.year"] - 3
  } 
}

age.data$Age.in.2021 <- NA

for (i in 1:length(age.data[,1])){
  if(age.data[i,"hatch.year"]==2020){
    age.data[i,"Age.in.2021"] <- "first.year"
  } else {
    age.data[i,"Age.in.2021"] <- "adult"
  }
}


write.csv(age.data, file="age.data.spring.2021.txt")
