# Cleans the environment
rm(list = ls())
dev.off()
cat("\014")

# Import libraries
require("ggplot2")


# # Import the .csv file with all the stations of Meteo France
# 
# MF_Stations <- read.csv("MeteoFrance/listeStations_Metro-OM_PackRadome.csv", sep = ";")
# 
# # Coordinates SOLA (approximately: 42°28'56.7"N 3°07'58.4"E)
# # Converted SOLA coordinates Latitude: 42.4824167, Longitude: 3.132888888888889
# 
# SOLA_lat <- 42.4824167
# SOLA_lon <- 3.132888888888889
# MF_Stations[which(abs(MF_Stations$Latitude - SOLA_lat) == min(abs(MF_Stations$Latitude - SOLA_lat))),]
# MF_Stations[which(abs(MF_Stations$Longitude - SOLA_lon) == min(abs(MF_Stations$Longitude - SOLA_lon))),]
# 
# # Calculate the distance of every station to SOLA. I make the assumption that a degree of Longitude is constant whatever the longitude in France.
# MF_Stations[["DistToSola"]] <- sqrt(((abs(MF_Stations$Latitude - SOLA_lat))**2)+((abs(MF_Stations$Latitude - SOLA_lat))**2))
# MF_Stations[which(MF_Stations[["DistToSola"]] == min(MF_Stations[["DistToSola"]])),]

# Cap Bear station files (one sample per day)

CapBear2013 <- read.csv("MeteoFrance/CapBearAPIClimatologie/CapBear2013.csv", sep = ";", dec = ",")
CapBear2014 <- read.csv("MeteoFrance/CapBearAPIClimatologie/CapBear2014.csv", sep = ";", dec = ",")
CapBear2015 <- read.csv("MeteoFrance/CapBearAPIClimatologie/CapBear2015.csv", sep = ";", dec = ",")
CapBear2016 <- read.csv("MeteoFrance/CapBearAPIClimatologie/CapBear2016.csv", sep = ";", dec = ",")
CapBear2017 <- read.csv("MeteoFrance/CapBearAPIClimatologie/CapBear2017.csv", sep = ";", dec = ",")
InfoVariables <- read.csv("MeteoFrance/MeteoFranceAPIClimatologieVariablesInformations/api_clim_table_parametres_quotidiens_20240103_354.csv", sep = ";")
# RR: HAUTEUR DE PRECIPITATIONS QUOTIDIENNE en MILLIMETRES ET 1/10
# DRR: DUREE DES PRECIPITATIONS QUOTIDIENNES en MINUTES
# FF2M: MOYENNE DES VITESSES DU VENT A 2 METRES QUOTIDIENNE en M/S ET 1/10
# FFM: MOYENNE DES VITESSES DU VENT A 10M QUOTIDIENNE en M/S ET 1/10

# Check column names are the same and in the same order

identical(colnames(CapBear2013), colnames(CapBear2014))
identical(colnames(CapBear2013), colnames(CapBear2015))
identical(colnames(CapBear2013), colnames(CapBear2016))
identical(colnames(CapBear2013), colnames(CapBear2017))

identical(colnames(CapBear2014), colnames(CapBear2015))
identical(colnames(CapBear2014), colnames(CapBear2016))
identical(colnames(CapBear2014), colnames(CapBear2017))

identical(colnames(CapBear2015), colnames(CapBear2016))
identical(colnames(CapBear2015), colnames(CapBear2017))

identical(colnames(CapBear2016), colnames(CapBear2017))

# rbind the dataframes # 2013, 2014, 2015,  2016, 2017: OK

CapBear2013_17 <- rbind(CapBear2013, CapBear2014, CapBear2015, CapBear2016, CapBear2017)

# Add date in the date format

CapBear2013_17$Date_YYYY <- substr(CapBear2013_17$DATE, 1, 4)
CapBear2013_17$Date_mm <- substr(CapBear2013_17$DATE, 5, 6)
CapBear2013_17$Date_dd <- substr(CapBear2013_17$DATE, 7, 8)
CapBear2013_17$Date_YYYYmmdd <- paste(CapBear2013_17$Date_YYYY, CapBear2013_17$Date_mm, CapBear2013_17$Date_dd, sep = "_")
CapBear2013_17$Date_num <- as.Date(CapBear2013_17$Date_YYYYmmdd, "%Y_%m_%d")

ggplot(CapBear2013_17, aes(x = Date_num, y = RR)) +
  geom_point() +
  geom_line()

ggplot(CapBear2013_17, aes(x = Date_num, y = FFM)) +
  geom_point() +
  geom_line()

# Export the rbind dataframe to file

write.csv(x = CapBear2013_17, file = "MeteoFrance/CapBearAPIClimatologie/CapBearConcatenated20132017.csv", row.names = FALSE)
