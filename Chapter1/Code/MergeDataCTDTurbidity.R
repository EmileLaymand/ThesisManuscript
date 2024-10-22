# Cleans the environment
rm(list = ls())
dev.off()
cat("\014")

# Import libraries

require("ggplot2")
require("phyloseq")
library(microViz)
library(stringr)
library(tidyverse)
library(gridExtra)
library(microbiomeSeq)
library(hues)
library(microbiome)
library(gt)

# library(raster)
# library(adehabitatHR)
# 
# r <- raster("../CTD_CBE19_turbidite_2013_2017/2013/somlit_20131009s.asc")
# r <- raster("../CTD_CBE19_turbidite_2013_2017/2013/somlit_20131009s.asc")
# 
# read.csv("../CTD_CBE19_turbidite_2013_2017/Test/somlit_20131009s.asc", sep = "\t", header = TRUE)

# Import libraries



# Import files

# CTDFiles <- list.files(pattern="CTDallTSVOneFolder/\\.tsv$")

PathCTD <- "CTDallTSVOneFolder/"
CTDFiles <- list.files(path = PathCTD, pattern = "*.tsv")
ListCTDDF <- lapply(CTDFiles, function(x) read.csv(paste(PathCTD, x, sep = ""),  sep = "\t"))
# lapply(ListCTDDF, colnames)
# MergedCTDDF <- ListCTDDF %>% reduce(full_join)
MergedCTDDF <- bind_rows(ListCTDDF) 

# Check if the number of all observations (i.e., dates and each depth for each date) is the same between MergedCTDDF and ListCTDDF
sum(sapply(ListCTDDF, nrow))
nrow(MergedCTDDF) # Seems that it works

# Check number of columns in MergedCTDDF
ncol(MergedCTDDF)
#View(MergedCTDDF[is.na(MergedCTDDF[["Sal00.1"]]) == FALSE,]) # OK, some columns look duplicated (e.g., Sal00.1) but it is because they are actually duplicated in the original files.

# Add a turbidity column with the values of all columns with turbidity (e.g., Turbidite, turbidite, etc.)

MergedCTDDF$Turbidity <- NA

ReceivingCol <- "Turbidity"
Col2Add <- "Turbidite"
MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]

ReceivingCol <- "Turbidity"
Col2Add <- "turbidite"
MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]

ReceivingCol <- "Turbidity"
Col2Add <- "SeaTurbMtr"
MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]


# Add a TheoreticalDepth column (with the values of columns: Profondeur.1, then PrdM, then profondeur, then DepSM.1, then Profondeur)
# PrSM is suspicious, as there are values below 27 (e.g., 87) which are not supposed to be present, as the depth below SOLA is only 27 m.
# There will be some NA, that I will remove afterwards

MergedCTDDF$TheoreticalDepth <- NA

ReceivingCol <- "TheoreticalDepth"
Col2Add <- "Profondeur.1"
MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]

ReceivingCol <- "TheoreticalDepth"
Col2Add <- "PrdM"
MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]

ReceivingCol <- "TheoreticalDepth"
Col2Add <- "profondeur"
MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]

ReceivingCol <- "TheoreticalDepth"
Col2Add <- "DepSM.1"
MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]

ReceivingCol <- "TheoreticalDepth"
Col2Add <- "Profondeur"
MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]

# ReceivingCol <- "TheoreticalDepth"
# Col2Add <- "PrSM"
# MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]

# Add a Date_num column
#----------------------

# Fill the DDMMMYYYY column

MergedCTDDF$MMMnum <- NA
MergedCTDDF[which(MergedCTDDF[["MMM"]] == "Oct"),"MMMnum"] <- 10
MergedCTDDF[which(MergedCTDDF[["MMM"]] == "Nov"),"MMMnum"] <- 11
MergedCTDDF[which(MergedCTDDF[["MMM"]] == "Dec"),"MMMnum"] <- 12

MergedCTDDF$DDMMMYYYY <- paste(MergedCTDDF$DD, MergedCTDDF$MMMnum, MergedCTDDF$YYYY, sep="_")
MergedCTDDF[MergedCTDDF[["DDMMMYYYY"]] == "NA_NA_NA","DDMMMYYYY"] <- NA
MergedCTDDF$DDMMMYYYY_date <- as.Date(MergedCTDDF$DDMMMYYYY, "%d_%m_%Y")

#as.Date(EukaAll_Metadata$Date_Euk, "%d_%m_%Y")
# MergedCTDDF$Date_num <- NA

# Fill the mm.dd.yyyy_date column

MergedCTDDF$mm.dd.yyyy_date <- as.Date(MergedCTDDF$mm.dd.yyyy, "%m/%d/%Y")

# 0-Replace all spaces by dashes 

MergedCTDDF$Date <- gsub(" ", "-", MergedCTDDF$Date)
MergedCTDDF$Date <- gsub("/", "-", MergedCTDDF$Date)

# 1-Fill the DateCol_num column

MergedCTDDF$Date_d <- sapply(strsplit(MergedCTDDF$Date, split = "-"), function(x) x[1])
MergedCTDDF$Date_m <- sapply(strsplit(MergedCTDDF$Date, split = "-"), function(x) x[2])
MergedCTDDF$Date_y <- sapply(strsplit(MergedCTDDF$Date, split = "-"), function(x) x[3])

# 2-Replace the months in names by numbers

MergedCTDDF$Date_m_num <- NA
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "Jan"),"Date_m_num"] <- 01
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "Feb"),"Date_m_num"] <- 02
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "Mar"),"Date_m_num"] <- 03
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "Apr"),"Date_m_num"] <- 04
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "May"),"Date_m_num"] <- 05
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "Jun"),"Date_m_num"] <- 06
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "Jul"),"Date_m_num"] <- 07
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "Aug"),"Date_m_num"] <- 08
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "Sep"),"Date_m_num"] <- 09
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "Oct"),"Date_m_num"] <- 10
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "Nov"),"Date_m_num"] <- 11
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "Dec"),"Date_m_num"] <- 12

# 2bis-Sometimes, there is a number for the month so keep the number

MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "01"),"Date_m_num"] <- 01
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "02"),"Date_m_num"] <- 02
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "03"),"Date_m_num"] <- 03
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "04"),"Date_m_num"] <- 04
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "05"),"Date_m_num"] <- 05
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "06"),"Date_m_num"] <- 06
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "07"),"Date_m_num"] <- 07
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "08"),"Date_m_num"] <- 08
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "09"),"Date_m_num"] <- 09
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "10"),"Date_m_num"] <- 10
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "11"),"Date_m_num"] <- 11
MergedCTDDF[which(MergedCTDDF[["Date_m"]] == "12"),"Date_m_num"] <- 12

# 3-Replace the format 14 by 2014

MergedCTDDF$Date_y_long <- NA
MergedCTDDF[which(MergedCTDDF[["Date_y"]] == 13),"Date_y_long"] <- 2013
MergedCTDDF[which(MergedCTDDF[["Date_y"]] == 14),"Date_y_long"] <- 2014
MergedCTDDF[which(MergedCTDDF[["Date_y"]] == 15),"Date_y_long"] <- 2015
MergedCTDDF[which(MergedCTDDF[["Date_y"]] == 16),"Date_y_long"] <- 2016
MergedCTDDF[which(MergedCTDDF[["Date_y"]] == 17),"Date_y_long"] <- 2017

MergedCTDDF[which(MergedCTDDF[["Date_y"]] == 2013),"Date_y_long"] <- 2013
MergedCTDDF[which(MergedCTDDF[["Date_y"]] == 2014),"Date_y_long"] <- 2014
MergedCTDDF[which(MergedCTDDF[["Date_y"]] == 2015),"Date_y_long"] <- 2015
MergedCTDDF[which(MergedCTDDF[["Date_y"]] == 2016),"Date_y_long"] <- 2016
MergedCTDDF[which(MergedCTDDF[["Date_y"]] == 2017),"Date_y_long"] <- 2017

# 4-Merge the columns

MergedCTDDF$DateColMerged <- paste(MergedCTDDF$Date_y_long, MergedCTDDF$Date_m_num, MergedCTDDF$Date_d, sep = "-")
MergedCTDDF[which(MergedCTDDF[["DateColMerged"]] == "NA-NA-NA"), "DateColMerged"] <- NA

# 5-Convert to Date format

MergedCTDDF$DateCol_num <- as.Date(MergedCTDDF$DateColMerged, "%Y-%m-%d")

# 6-Check it worked
MergedCTDDF[, c("Date", "DateCol_num")]

# Put all the date_num columns together --> Does not work

MergedCTDDF$Date_num <- NA
MergedCTDDF$Date_num <- as.Date(MergedCTDDF$Date_num)


ReceivingCol <- "Date_num"
Col2Add <- "DDMMMYYYY_date"
MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]

ReceivingCol <- "Date_num"
Col2Add <- "mm.dd.yyyy_date"
MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]

ReceivingCol <- "Date_num"
Col2Add <- "DateCol_num"
MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), ReceivingCol] <- MergedCTDDF[is.na(MergedCTDDF[[ReceivingCol]]) & (is.na(MergedCTDDF[[Col2Add]]) == FALSE), Col2Add]

# Keep only the Date_num, Turbidity and TheoreticalDepth columns for export

ShortDF <- MergedCTDDF[, c("Date_num", "TheoreticalDepth", "Turbidity")]

# Quick checks 

# length(unique(ShortDF[!complete.cases(ShortDF),"Date_num"]))
# ShortDF[which(!complete.cases(ShortDF)),"Date_num"]
which(is.na(MergedCTDDF[["Date_num"]])) # If integer(0): means that there is no missing value in the Date_num column

# Count number of complete cases

length(unique(ShortDF[["Date_num"]])) # Two elements seem to miss
ShortDF[["Date_num"]][match(unique(ShortDF[["Date_num"]]), ShortDF[["Date_num"]])]

# Find which are the two missing samples

sapply(ListCTDDF, function(x) dim(x)[1])
length(sapply(ListCTDDF, function(x) dim(x)[1])) # length 204

sort(sapply(ListCTDDF, function(x) dim(x)[1]))
sort(table(ShortDF[["Date_num"]]))

sort(table(ShortDF[["Date_num"]]))[!(sort(table(ShortDF[["Date_num"]])) %in% sort(sapply(ListCTDDF, function(x) dim(x)[1])))]

# Look for the indexes of the two missing samples

MergedCTDDF[MergedCTDDF[["Date_num"]] == "2017-02-15",]
which(MergedCTDDF[["Date_num"]] == "2017-02-15") # Problem solved: two files were actually duplicated

# See what is the best depth (i.e., deep enough and present in enough samples to be relevant)

ShortDFCompleteCases <- ShortDF[complete.cases(ShortDF),]
table(ShortDFCompleteCases[["TheoreticalDepth"]]) # Best is maybe 3m and 20m

# Format the dataframe to have one column with turbidity at 3m, and another with turbidity at 20m

ShortDFCompleteCases_3m <- ShortDFCompleteCases[ShortDFCompleteCases[["TheoreticalDepth"]] == 3,]
colnames(ShortDFCompleteCases_3m)[colnames(ShortDFCompleteCases_3m) == "Turbidity"] <- "Turbidity_3m"
ShortDFCompleteCases_3m$TheoreticalDepth <- NULL

ShortDFCompleteCases_20m <- ShortDFCompleteCases[ShortDFCompleteCases[["TheoreticalDepth"]] == 20,]
colnames(ShortDFCompleteCases_20m)[colnames(ShortDFCompleteCases_20m) == "Turbidity"] <- "Turbidity_20m"
ShortDFCompleteCases_20m$TheoreticalDepth <- NULL

ShortDFCompleteCases_3m_20m <- merge(x = ShortDFCompleteCases_3m, y = ShortDFCompleteCases_20m, by = "Date_num", all = FALSE)
plot(ShortDFCompleteCases_3m_20m$Turbidity_3m, ShortDFCompleteCases_3m_20m$Turbidity_20m)
cor(ShortDFCompleteCases_3m_20m$Turbidity_3m, ShortDFCompleteCases_3m_20m$Turbidity_20m) # Correlation is very good, 0.9667534.

# Plot turbidity against time

ggplot(ShortDFCompleteCases_3m_20m, aes(x = Date_num, y = Turbidity_3m)) +
  geom_line() +
  geom_point()

ggplot(ShortDFCompleteCases_3m_20m, aes(x = Date_num, y = Turbidity_20m)) +
  geom_line() +
  geom_point()

# Export to csv format # /!\ Data quality must be checked out (script and with Paul)

write.csv(ShortDFCompleteCases_3m_20m, file = "CTDFromPaulMerged/Turbidity3m20mMerged.csv", row.names = FALSE)
