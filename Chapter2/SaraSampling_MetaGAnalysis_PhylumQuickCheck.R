#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  GENERAL CLEANING AND SETTING  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Cleans the environment
rm(list = ls())
dev.off()
cat("\014")

## Set the working directory
#setwd("/home/emilelaymand/Documents/Science/Master2_MES/StageM2")

library(plyr)
require("ggplot2")
require("phyloseq")
library(microViz)
library(stringr)
library(tidyverse)
library(gridExtra)
library(microbiomeSeq)
library(hues)
library(microbiome)
library(plotly)
library(ggfortify)
library(FactoMineR)
library("factoextra")
library("corrplot")
library(missMDA)
library(reshape2)
library(multcomp)
library(FSA)
library(plotrix)
library(ggbreak)
library(ggforce)
library(rstatix)
library(vegan)
library(lubridate)

source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/Generics.R")
source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/Pathnames.R")
source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/phyloseqExtended.R")

#==========================================#
#            Import all datasets           #
#==========================================#

# Load the metadata
#------------------

# Import the correspondance between SampleIDs and the AGTUIDs from Fasteris

DF_IDCorrespondance <- read.csv("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/EnvironmentalParameters/Correspondance_SampleName_AGTUIndex.csv")
DF_IDCorrespondance$SampleShort <- substr(DF_IDCorrespondance$SampleName, 3, 7)
# DF_IDCorrespondance$Site <- substr(DF_IDCorrespondance$SampleName, 3, 5)

# Import metadata

DF_Metadata <- read.csv("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/EnvironmentalParameters/EnvironmentalParameters.csv")
DF_Metadata$SampleShort <- DF_Metadata$sampleID
DF_Metadata$sampleID <- NULL

# Merge both dataframes without removing the duplicates (i.e., samples that have been sequenced twice, that do have the same environmental characteristics but not the same IDs)

DF_MergedMetadata <- join(DF_IDCorrespondance, DF_Metadata, type="full") # /!\ Hard check (really check)

# Load the Overview dataframes, and for each remove the percentage signs
#-----------------------------------------------------------------------

# Confidence 0.0 

DF_MetaG_Overview_Conf_0.0 <- read.csv("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaG/SortiesKraken2/Pavian_Confidence_0.0/Analysis/Overview/SaraSampling_MetaG_Kraken2_Confidence_0.0_Overview.csv")

DF_MetaG_Overview_Conf_0.0_woPerc <- DF_MetaG_Overview_Conf_0.0

for (colID in c("Classified.reads", "Chordate.reads", "Artificial.reads", "Unclassified.reads", "Microbial.reads", "Bacterial.reads", "Viral.reads", "Fungal.reads", "Protozoan.reads")) {
  DF_MetaG_Overview_Conf_0.0_woPerc[[colID]] <- as.numeric(substr(DF_MetaG_Overview_Conf_0.0_woPerc[[colID]], 1, nchar(DF_MetaG_Overview_Conf_0.0_woPerc[[colID]])-1))
}

DF_MetaG_Overview_Conf_0.0_woPerc$Confidence <- "0.0"
DF_MetaG_Overview_Conf_0.0_woPerc$AGTUIndex <- substr(DF_MetaG_Overview_Conf_0.0_woPerc$Name, 22, nchar(DF_MetaG_Overview_Conf_0.0_woPerc$Name))

# Confidence 0.1

DF_MetaG_Overview_Conf_0.1 <- read.csv("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaG/SortiesKraken2/Pavian_Confidence_0.1/Analysis/Overview/SaraSampling_MetaG_Kraken2_Confidence_0.1_Overview.csv")

DF_MetaG_Overview_Conf_0.1_woPerc <- DF_MetaG_Overview_Conf_0.1

for (colID in c("Classified.reads", "Chordate.reads", "Artificial.reads", "Unclassified.reads", "Microbial.reads", "Bacterial.reads", "Viral.reads", "Fungal.reads", "Protozoan.reads")) {
  DF_MetaG_Overview_Conf_0.1_woPerc[[colID]] <- as.numeric(substr(DF_MetaG_Overview_Conf_0.1_woPerc[[colID]], 1, nchar(DF_MetaG_Overview_Conf_0.1_woPerc[[colID]])-1))
}

DF_MetaG_Overview_Conf_0.1_woPerc$Confidence <- "0.1"
DF_MetaG_Overview_Conf_0.1_woPerc$AGTUIndex <- substr(DF_MetaG_Overview_Conf_0.1_woPerc$Name, 22, nchar(DF_MetaG_Overview_Conf_0.1_woPerc$Name))

# Confidence 0.6

DF_MetaG_Overview_Conf_0.6 <- read.csv("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaG/SortiesKraken2/Pavian_Confidence_0.6/Analysis/Overview/SaraSampling_MetaG_Kraken2_Confidence_0.6_Overview.csv")

DF_MetaG_Overview_Conf_0.6_woPerc <- DF_MetaG_Overview_Conf_0.6

for (colID in c("Classified.reads", "Chordate.reads", "Artificial.reads", "Unclassified.reads", "Microbial.reads", "Bacterial.reads", "Viral.reads", "Fungal.reads", "Protozoan.reads")) {
  DF_MetaG_Overview_Conf_0.6_woPerc[[colID]] <- as.numeric(substr(DF_MetaG_Overview_Conf_0.6_woPerc[[colID]], 1, nchar(DF_MetaG_Overview_Conf_0.6_woPerc[[colID]])-1))
}

DF_MetaG_Overview_Conf_0.6_woPerc$Confidence <- "0.6"
DF_MetaG_Overview_Conf_0.6_woPerc$AGTUIndex <- substr(DF_MetaG_Overview_Conf_0.6_woPerc$Name, 22, nchar(DF_MetaG_Overview_Conf_0.6_woPerc$Name))

# rbind the dataframes with the three confidences, then merge the dataframes with the metadata
#---------------------------------------------------------------------------------------------

# Check colnames are in the same order

identical(colnames(DF_MetaG_Overview_Conf_0.0_woPerc), colnames(DF_MetaG_Overview_Conf_0.1_woPerc)) # Must be TRUE
identical(colnames(DF_MetaG_Overview_Conf_0.1_woPerc), colnames(DF_MetaG_Overview_Conf_0.6_woPerc)) # Must be TRUE
identical(colnames(DF_MetaG_Overview_Conf_0.0_woPerc), colnames(DF_MetaG_Overview_Conf_0.6_woPerc)) # Must be TRUE

# rbind the Pavian reports

DF_MetaG_Overview_rbind <- rbind(DF_MetaG_Overview_Conf_0.0_woPerc, DF_MetaG_Overview_Conf_0.1_woPerc, DF_MetaG_Overview_Conf_0.6_woPerc)

# Merge with the metadata

DF_MetaG_Overview_rbind

DF_MergedMetadataPavian <- join(DF_MetaG_Overview_rbind, DF_MergedMetadata, type="full") # /!\ Hard check (really check)

# Add a column to specify if the sample is from France or Chile

DF_MergedMetadataPavian$Country <- substr(DF_MergedMetadataPavian$SampleShort, 1, 1)

# Make the siteID column complete

DF_MergedMetadataPavian$siteID <- substr(DF_MergedMetadataPavian$SampleShort, 1, 3)

#======================================================================#
# Import the files with the number of reads attributed to each species # # Copy-pasted from SaraSampling_MetaGMetaTComparison.R on 2024/05/13. Check everything
#======================================================================#

# Import the files
#-----------------

# Conf 0.0

DF_PavianSpecies_Conf_0.0 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaG/SortiesKraken2/Pavian_Confidence_0.0/Analysis/RawPavianReports/SaraSampling_MetaG_Kraken2_Confidence_0.0_Fungi_Phylum.tsv", sep = "\t", header = TRUE)

colnames(DF_PavianSpecies_Conf_0.0)[5:114] <- gsub("\\.", "-", substr(colnames(DF_PavianSpecies_Conf_0.0)[5:114], 23, nchar(colnames(DF_PavianSpecies_Conf_0.0)[5:114])-11))

# Conf 0.1

DF_PavianSpecies_Conf_0.1 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaG/SortiesKraken2/Pavian_Confidence_0.1/Analysis/RawPavianReports/SaraSampling_MetaG_Kraken2_Confidence_0.1_Fungi_Phylum.tsv", sep = "\t", header = TRUE)

colnames(DF_PavianSpecies_Conf_0.1)[5:114] <- gsub("\\.", "-", substr(colnames(DF_PavianSpecies_Conf_0.1)[5:114], 23, nchar(colnames(DF_PavianSpecies_Conf_0.1)[5:114])-11))

# Conf 0.6

DF_PavianSpecies_Conf_0.6 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaG/SortiesKraken2/Pavian_Confidence_0.6/Analysis/RawPavianReports/SaraSampling_MetaG_Kraken2_Confidence_0.6_Fungi_Phylum.tsv", sep = "\t", header = TRUE)

colnames(DF_PavianSpecies_Conf_0.6)[5:114] <- gsub("\\.", "-", substr(colnames(DF_PavianSpecies_Conf_0.6)[5:114], 23, nchar(colnames(DF_PavianSpecies_Conf_0.6)[5:114])-11))

# Normalize read counts by the number of reads in the sample
#-----------------------------------------------------------

# Conf 0.0

DF_PavianSpecies_Conf_0.0_Norm <- DF_PavianSpecies_Conf_0.0

for (i in 5:114) {
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.0_Norm)[i]
  DF_PavianSpecies_Conf_0.0_Norm[[NameOfCol]] <- DF_PavianSpecies_Conf_0.0_Norm[[NameOfCol]]/DF_MetaG_Overview_Conf_0.0_woPerc[DF_MetaG_Overview_Conf_0.0_woPerc[["AGTUIndex"]] == NameOfCol,"Number.of.raw.reads"]
} # To be checked

# Conf 0.1 

DF_PavianSpecies_Conf_0.1_Norm <- DF_PavianSpecies_Conf_0.1

for (i in 5:114) {
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.1_Norm)[i]
  DF_PavianSpecies_Conf_0.1_Norm[[NameOfCol]] <- DF_PavianSpecies_Conf_0.1_Norm[[NameOfCol]]/DF_MetaG_Overview_Conf_0.1_woPerc[DF_MetaG_Overview_Conf_0.1_woPerc[["AGTUIndex"]] == NameOfCol,"Number.of.raw.reads"]
} # To be checked

# Conf 0.6 

DF_PavianSpecies_Conf_0.6_Norm <- DF_PavianSpecies_Conf_0.6

for (i in 5:114) {
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.6_Norm)[i]
  DF_PavianSpecies_Conf_0.6_Norm[[NameOfCol]] <- DF_PavianSpecies_Conf_0.6_Norm[[NameOfCol]]/DF_MetaG_Overview_Conf_0.6_woPerc[DF_MetaG_Overview_Conf_0.6_woPerc[["AGTUIndex"]] == NameOfCol,"Number.of.raw.reads"]
} # To be checked

# Rename rows according to the species name

row.names(DF_PavianSpecies_Conf_0.0_Norm) <- DF_PavianSpecies_Conf_0.0_Norm$name
row.names(DF_PavianSpecies_Conf_0.1_Norm) <- DF_PavianSpecies_Conf_0.1_Norm$name
row.names(DF_PavianSpecies_Conf_0.6_Norm) <- DF_PavianSpecies_Conf_0.6_Norm$name

# Replace colnames by sample names # Added on 2024/05/13

# Conf 0.0

NewColnames <- c()
for (Colname_i in colnames(DF_PavianSpecies_Conf_0.0_Norm)[5:114]){
  NewColnames <- c(NewColnames, DF_IDCorrespondance[DF_IDCorrespondance[["AGTUIndex"]] == Colname_i, "SampleName"])
}

# Check there is no duplicated value

identical(length(unique(NewColnames)), length(NewColnames))

# Replace colnames # To be checked

colnames(DF_PavianSpecies_Conf_0.0_Norm)[5:114] <- NewColnames

# Sort colnames to have all samples from the same site close in the dataframe

DF_PavianSpecies_Conf_0.0_Norm_ColSort <- DF_PavianSpecies_Conf_0.0_Norm[, c("name", "taxRank", "taxID", "Max", sort(colnames(DF_PavianSpecies_Conf_0.0_Norm)[5:114]), "lineage")]

DF_PavianSpecies_Conf_0.0_Norm <- DF_PavianSpecies_Conf_0.0_Norm_ColSort

# Conf 0.1

NewColnames <- c()
for (Colname_i in colnames(DF_PavianSpecies_Conf_0.1_Norm)[5:114]){
  NewColnames <- c(NewColnames, DF_IDCorrespondance[DF_IDCorrespondance[["AGTUIndex"]] == Colname_i, "SampleName"])
}

# Check there is no duplicated value

identical(length(unique(NewColnames)), length(NewColnames))

# Replace colnames # To be checked

colnames(DF_PavianSpecies_Conf_0.1_Norm)[5:114] <- NewColnames

# Sort colnames to have all samples from the same site close in the dataframe

DF_PavianSpecies_Conf_0.1_Norm_ColSort <- DF_PavianSpecies_Conf_0.1_Norm[, c("name", "taxRank", "taxID", "Max", sort(colnames(DF_PavianSpecies_Conf_0.1_Norm)[5:114]), "lineage")]

DF_PavianSpecies_Conf_0.1_Norm <- DF_PavianSpecies_Conf_0.1_Norm_ColSort

# Conf 0.6

NewColnames <- c()
for (Colname_i in colnames(DF_PavianSpecies_Conf_0.6_Norm)[5:114]){
  NewColnames <- c(NewColnames, DF_IDCorrespondance[DF_IDCorrespondance[["AGTUIndex"]] == Colname_i, "SampleName"])
}

# Check there is no duplicated value

identical(length(unique(NewColnames)), length(NewColnames))

# Replace colnames # To be checked

colnames(DF_PavianSpecies_Conf_0.6_Norm)[5:114] <- NewColnames

# Sort colnames to have all samples from the same site close in the dataframe

DF_PavianSpecies_Conf_0.6_Norm_ColSort <- DF_PavianSpecies_Conf_0.6_Norm[, c("name", "taxRank", "taxID", "Max", sort(colnames(DF_PavianSpecies_Conf_0.6_Norm)[5:114]), "lineage")]

DF_PavianSpecies_Conf_0.6_Norm <- DF_PavianSpecies_Conf_0.6_Norm_ColSort

# Look at the most abundant species in each sample 
#-------------------------------------------------

# Conf 0.0

ListSorted_Conf_0.0 <- list() # Check this out thoroughly
for (i in 5:114){
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.0_Norm)[i]
  DF_Sorted <- DF_PavianSpecies_Conf_0.0_Norm[order(DF_PavianSpecies_Conf_0.0_Norm[[NameOfCol]], decreasing = TRUE), c("name", NameOfCol)] # Check the ordering of names is good
  colnames(DF_Sorted) <- c(paste("name_", NameOfCol, sep = ""), NameOfCol)
  ListSorted_Conf_0.0[[NameOfCol]] <- DF_Sorted
}

SortedSpecies_Conf_0.0 <- bind_cols(ListSorted_Conf_0.0)
row.names(SortedSpecies_Conf_0.0) <- NULL

write.csv(SortedSpecies_Conf_0.0, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/MetaG/Phylum/RelativeAbundance/RelativeAbundancePhylum_Conf_0.0.csv", row.names = FALSE)

# Conf 0.1

ListSorted_Conf_0.1 <- list() # Check this out thoroughly
for (i in 5:114){
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.1_Norm)[i]
  DF_Sorted <- DF_PavianSpecies_Conf_0.1_Norm[order(DF_PavianSpecies_Conf_0.1_Norm[[NameOfCol]], decreasing = TRUE), c("name", NameOfCol)] # Check the ordering of names is good
  colnames(DF_Sorted) <- c(paste("name_", NameOfCol, sep = ""), NameOfCol)
  ListSorted_Conf_0.1[[NameOfCol]] <- DF_Sorted
}

SortedSpecies_Conf_0.1 <- bind_cols(ListSorted_Conf_0.1)
row.names(SortedSpecies_Conf_0.1) <- NULL

write.csv(SortedSpecies_Conf_0.1, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/MetaG/Phylum/RelativeAbundance/RelativeAbundancePhylum_Conf_0.1.csv", row.names = FALSE)

# Conf 0.6

ListSorted_Conf_0.6 <- list() # Check this out thoroughly
for (i in 5:114){
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.6_Norm)[i]
  DF_Sorted <- DF_PavianSpecies_Conf_0.6_Norm[order(DF_PavianSpecies_Conf_0.6_Norm[[NameOfCol]], decreasing = TRUE), c("name", NameOfCol)] # Check the ordering of names is good
  colnames(DF_Sorted) <- c(paste("name_", NameOfCol, sep = ""), NameOfCol)
  ListSorted_Conf_0.6[[NameOfCol]] <- DF_Sorted
}

SortedSpecies_Conf_0.6 <- bind_cols(ListSorted_Conf_0.6)
row.names(SortedSpecies_Conf_0.6) <- NULL

write.csv(SortedSpecies_Conf_0.6, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/MetaG/Phylum/RelativeAbundance/RelativeAbundancePhylum_Conf_0.6.csv", row.names = FALSE)

#===========================================================================# 
# End of copy-paste from SaraSampling_MetaGAnalysis_Species.R on 2024/05/16 #
#===========================================================================# 

# Make a boxplot of the relative abundance of each phylum 

# Conf 0.0
#---------

# Select only the good columns 
DF_PavianSpecies_Conf_0.0_Norm_ForBox <- DF_PavianSpecies_Conf_0.0_Norm
DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol <- DF_PavianSpecies_Conf_0.0_Norm_ForBox[, !(colnames(DF_PavianSpecies_Conf_0.0_Norm_ForBox) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]

# Select only samples from France

ColToKeep_0.0 <- colnames(DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol)[which(substr(colnames(DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol), 3, 5) %in% c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF"))]

DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc <- DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol[, colnames(DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol) %in% ColToKeep_0.0]
DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc <- DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc %>% replace(is.na(.), 0)
DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc_t <- data.frame(t(DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc))
DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc_t$Confidence <- "0.0"

# Conf 0.1
#---------

# Select only the good columns 
DF_PavianSpecies_Conf_0.1_Norm_ForBox <- DF_PavianSpecies_Conf_0.1_Norm
DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol <- DF_PavianSpecies_Conf_0.1_Norm_ForBox[, !(colnames(DF_PavianSpecies_Conf_0.1_Norm_ForBox) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]

# Select only samples from France

ColToKeep_0.1 <- colnames(DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol)[which(substr(colnames(DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol), 3, 5) %in% c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF"))]

DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol_To_calc <- DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol[, colnames(DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol) %in% ColToKeep_0.1]
DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol_To_calc <- DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol_To_calc %>% replace(is.na(.), 0)
DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol_To_calc_t <- data.frame(t(DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol_To_calc))
DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol_To_calc_t$Confidence <- "0.1"

# Conf 0.6
#---------

# Select only the good columns 
DF_PavianSpecies_Conf_0.6_Norm_ForBox <- DF_PavianSpecies_Conf_0.6_Norm
DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol <- DF_PavianSpecies_Conf_0.6_Norm_ForBox[, !(colnames(DF_PavianSpecies_Conf_0.6_Norm_ForBox) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]

# Select only samples from France

ColToKeep_0.6 <- colnames(DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol)[which(substr(colnames(DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol), 3, 5) %in% c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF"))]

DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol_To_calc <- DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol[, colnames(DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol) %in% ColToKeep_0.6]
DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol_To_calc <- DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol_To_calc %>% replace(is.na(.), 0)
DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol_To_calc_t <- data.frame(t(DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol_To_calc))
DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol_To_calc_t$Confidence <- "0.6"

# rbind the dataframes, convert to long, and make the plot

# Add missing columns full of 0 (I replace fornow with NAs)

# 0.1
MissingIn0.1 <- colnames(DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc_t)[which(!(colnames(DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc_t) %in% colnames(DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol_To_calc_t)))]

for (i in MissingIn0.1){
  # DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol_To_calc_t[[i]] <- 0
  DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol_To_calc_t[[i]] <- NA
}

# 0.6
MissingIn0.6 <- colnames(DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc_t)[which(!(colnames(DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc_t) %in% colnames(DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol_To_calc_t)))]

for (i in MissingIn0.6){
  # DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol_To_calc_t[[i]] <- 0
  DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol_To_calc_t[[i]] <- NA
}

# Sort the columns in each dataframe

DF_ForBoxPhylum_0.0 <- DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc_t[,sort(colnames(DF_PavianSpecies_Conf_0.0_Norm_ForBox_GoodCol_To_calc_t))]
DF_ForBoxPhylum_0.1 <- DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol_To_calc_t[,sort(colnames(DF_PavianSpecies_Conf_0.1_Norm_ForBox_GoodCol_To_calc_t))]
DF_ForBoxPhylum_0.6 <- DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol_To_calc_t[,sort(colnames(DF_PavianSpecies_Conf_0.6_Norm_ForBox_GoodCol_To_calc_t))]

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

# 0.0

DF_ForBoxPhylum_0.0 <- DF_ForBoxPhylum_0.0[!(row.names(DF_ForBoxPhylum_0.0) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

DF_ForBoxPhylum_0.0[["SiteID"]] <- substr(row.names(DF_ForBoxPhylum_0.0), 3, 5)

DF_ForBoxPhylum_0.0 <- DF_ForBoxPhylum_0.0[DF_ForBoxPhylum_0.0[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

DF_ForBoxPhylum_0.0$SiteID <- NULL

# 0.1

DF_ForBoxPhylum_0.1 <- DF_ForBoxPhylum_0.1[!(row.names(DF_ForBoxPhylum_0.1) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

DF_ForBoxPhylum_0.1[["SiteID"]] <- substr(row.names(DF_ForBoxPhylum_0.1), 3, 5)

DF_ForBoxPhylum_0.1 <- DF_ForBoxPhylum_0.1[DF_ForBoxPhylum_0.1[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

DF_ForBoxPhylum_0.1$SiteID <- NULL

# 0.6 

DF_ForBoxPhylum_0.6 <- DF_ForBoxPhylum_0.6[!(row.names(DF_ForBoxPhylum_0.6) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),]
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

DF_ForBoxPhylum_0.6[["SiteID"]] <- substr(row.names(DF_ForBoxPhylum_0.6), 3, 5)

DF_ForBoxPhylum_0.6 <- DF_ForBoxPhylum_0.6[DF_ForBoxPhylum_0.6[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

DF_ForBoxPhylum_0.6$SiteID <- NULL

# End of added on 2024-06-03

# Check the columns are in the same order in the three dataframes

identical(colnames(DF_ForBoxPhylum_0.0), colnames(DF_ForBoxPhylum_0.1)) # Must be TRUE
identical(colnames(DF_ForBoxPhylum_0.0), colnames(DF_ForBoxPhylum_0.6)) # Must be TRUE

# rbind

DF_ForBoxPhylum_rbind <- rbind(DF_ForBoxPhylum_0.0, DF_ForBoxPhylum_0.1, DF_ForBoxPhylum_0.6)

# Added on 2024-06-17

# Replace NAs by 0

DF_ForBoxPhylum_rbind <- DF_ForBoxPhylum_rbind %>% replace(is.na(.), 0)

# End of added on 2024-06-17

DF_ForBoxPhylum_rbind_long <- DF_ForBoxPhylum_rbind %>%
  pivot_longer(
    cols = colnames(DF_ForBoxPhylum_0.0)[which(colnames(DF_ForBoxPhylum_0.0) != "Confidence")],
    names_to = "Phylum",
    values_to = "RelAbundance"
  )

# Make the plot

ggplot(data = DF_ForBoxPhylum_rbind_long, aes(x = Phylum, y = RelAbundance, color = Confidence)) +
  #stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  #ylab(label = "Percent of Chytridiomycota amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )

ggplot(data = DF_ForBoxPhylum_rbind_long, aes(x = Phylum, y = RelAbundance)) +
  facet_grid(rows = vars(Confidence), scales = "free") +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  #scale_y_continuous(trans='log10') +
  #ylab(label = "Percent of Chytridiomycota amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

# Convert all 0 to a very small value (1e-12)
# Note: minimum value except 0 is 4.281023e-09 # Edit 2024-06-17: with the selected 7 sites, the min is 6.109345e-09
# min(DF_ForBoxPhylum_rbind_long$RelAbundance[(DF_ForBoxPhylum_rbind_long$RelAbundance != 0) & (is.na(DF_ForBoxPhylum_rbind_long$RelAbundance) == FALSE)])
DF_ForBoxPhylum_rbind_long_zero_replaced <- DF_ForBoxPhylum_rbind_long
DF_ForBoxPhylum_rbind_long_zero_replaced[DF_ForBoxPhylum_rbind_long_zero_replaced == 0] <- 0.0000000001

# Added on 2024-06-17. Convert all values to percentage by multiplying all values by 100

DF_ForBoxPhylum_rbind_long_zero_replaced$RelAbundance <- DF_ForBoxPhylum_rbind_long_zero_replaced$RelAbundance*100

# End of added on 2024-06-17. Convert all values to percentage by multiplying all values by 100

p1 <- ggplot(data = DF_ForBoxPhylum_rbind_long_zero_replaced, aes(x = Phylum, y = RelAbundance)) +
  facet_grid(rows = vars(Confidence), scales = "free") +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  scale_y_continuous(trans='log10') +
  annotation_logticks(sides = "l") +
  ylab(label = "Percent of phylum amongst all reads") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/PercentPhylum_Log10_0replaced_Perc_AllConf.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/PercentPhylum_Log10_0replaced_Perc_AllConf.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

# Exact same, but with only Confidence 0.1 on the plot

DF_ForBoxPhylum_rbind_long_zero_replaced_0.1 <- DF_ForBoxPhylum_rbind_long_zero_replaced[DF_ForBoxPhylum_rbind_long_zero_replaced[["Confidence"]] == "0.1", ]

# Remove the phyla that were added to the plot because they were present at confidence 0.0 or that were present in other sites

DF_ForBoxPhylum_rbind_long_zero_replaced_0.1 <- DF_ForBoxPhylum_rbind_long_zero_replaced_0.1[!(DF_ForBoxPhylum_rbind_long_zero_replaced_0.1[["Phylum"]] %in% c("Olpidiomycota", "Sanchytriomycota")),]

# Convert to factor

DF_ForBoxPhylum_rbind_long_zero_replaced_0.1$Phylum <- factor(DF_ForBoxPhylum_rbind_long_zero_replaced_0.1$Phylum, levels = c("Ascomycota", "Basidiomycota", "Mucoromycota", "Chytridiomycota", "Microsporidia", "Zoopagomycota", "Cryptomycota", "Blastocladiomycota"))

# Make the plot

p1 <- ggplot(data = DF_ForBoxPhylum_rbind_long_zero_replaced_0.1, aes(x = Phylum, y = RelAbundance)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  scale_y_continuous(trans='log10') +
  annotation_logticks(sides = "l") +
  ylab(label = "Percent of phylum amongst all reads") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/PercentPhylum_Log10_0replaced_Perc_0.1.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/PercentPhylum_Log10_0replaced_Perc_0.1.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

# Same plot but with values in frequency and not in percentage

DF_ForBoxPhylum_rbind_long_zero_replaced_0.1_Fraction <- DF_ForBoxPhylum_rbind_long_zero_replaced_0.1
DF_ForBoxPhylum_rbind_long_zero_replaced_0.1_Fraction$RelAbundance <- DF_ForBoxPhylum_rbind_long_zero_replaced_0.1_Fraction$RelAbundance / 100

p1 <- ggplot(data = DF_ForBoxPhylum_rbind_long_zero_replaced_0.1_Fraction, aes(x = Phylum, y = RelAbundance)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1, outlier.size = 6) +
  scale_y_continuous(trans='log10') +
  annotation_logticks(sides = "l") +
  ylab(label = "Relative abundance of phylum amongst all reads") +
  geom_hline(yintercept = 0.000000001, color = "black", linetype = "dashed") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/FractionPhylum_Log10_0replaced_Perc_0.1.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/FractionPhylum_Log10_0replaced_Perc_0.1.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

# Subset per Confidence level to make small calculations

DF_ForBoxPhylum_rbind_0.1_QuickCalculations <- DF_ForBoxPhylum_rbind[DF_ForBoxPhylum_rbind[["Confidence"]] == "0.1",]
DF_ForBoxPhylum_rbind_0.1_QuickCalculations$Confidence <- NULL
sapply(DF_ForBoxPhylum_rbind_0.1_QuickCalculations, function(x) mean(x)*100)
sort(sapply(DF_ForBoxPhylum_rbind_0.1_QuickCalculations, function(x) mean(x)*100), decreasing = TRUE)

# Make the same plots but per site

DF_ForBoxPhylum_0.1_ForPlotPerSite <- DF_ForBoxPhylum_0.1
DF_ForBoxPhylum_0.1_ForPlotPerSite$Olpidiomycota <- NULL 
DF_ForBoxPhylum_0.1_ForPlotPerSite$Sanchytriomycota <- NULL
DF_ForBoxPhylum_0.1_ForPlotPerSite <- DF_ForBoxPhylum_0.1_ForPlotPerSite %>% replace(is.na(.), 0)
DF_ForBoxPhylum_0.1_ForPlotPerSite$SiteID <- substr(row.names(DF_ForBoxPhylum_0.1_ForPlotPerSite), 3, 5)
DF_ForBoxPhylum_0.1_ForPlotPerSite$Confidence <- NULL

DF_ForBoxPhylum_0.1_ForPlotPerSite_long <- DF_ForBoxPhylum_0.1_ForPlotPerSite %>%
  pivot_longer(
    cols = colnames(DF_ForBoxPhylum_0.1_ForPlotPerSite)[which(colnames(DF_ForBoxPhylum_0.1_ForPlotPerSite) != "SiteID")],
    names_to = "Phylum",
    values_to = "RelAbundance"
  )

DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced <- DF_ForBoxPhylum_0.1_ForPlotPerSite_long
DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced[DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced == 0] <- 0.0000000001

DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced$SiteID <- factor(DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced$SiteID, levels = c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"))
DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced$Phylum <- factor(DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced$Phylum, levels = c("Ascomycota", "Basidiomycota", "Mucoromycota", "Chytridiomycota", "Microsporidia", "Zoopagomycota", "Cryptomycota", "Blastocladiomycota"))

dev.new()
p1 <- ggplot(data = DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced, aes(x = Phylum, y = RelAbundance, fill = SiteID)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1, outlier.size = 6) +
  scale_y_continuous(trans='log10') +
  annotation_logticks(sides = "l") +
  #scale_fill_manual(values = c("deepskyblue", "gold", "gold2", "darkolivegreen1", "darkolivegreen3", "darkolivegreen4", "forestgreen")) +
  scale_fill_manual(values = c("darkorchid2", "gold", "goldenrod4", "deepskyblue", "dodgerblue4", "darkolivegreen1", "forestgreen")) +
  ylab(label = "Relative abundance amongst all reads") +
  geom_hline(yintercept = 0.000000001, color = "black", linetype = "dashed") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/PercentPhylum_Log10_0replaced_Perc_0.1_PerSite.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/PercentPhylum_Log10_0replaced_Perc_0.1_PerSite.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

# Same but with the fraction of reads instead of the percentage

# DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced_FracPlot <- DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced
# DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced_FracPlot$RelAbundance <- DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced_FracPlot$RelAbundance / 100
# 
# dev.new()
# p1 <- ggplot(data = DF_ForBoxPhylum_0.1_ForPlotPerSite_long_zero_replaced_FracPlot, aes(x = Phylum, y = RelAbundance, fill = SiteID)) +
#   geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
#   scale_y_continuous(trans='log10') +
#   annotation_logticks(sides = "l") +
#   #scale_fill_manual(values = c("deepskyblue", "gold", "gold2", "darkolivegreen1", "darkolivegreen3", "darkolivegreen4", "forestgreen")) +
#   scale_fill_manual(values = c("darkorchid2", "gold", "goldenrod4", "deepskyblue", "dodgerblue4", "darkolivegreen1", "forestgreen")) +
#   ylab(label = "Percent of phylum amongst all reads") +
#   geom_hline(yintercept = 0.000000001, color = "black", linetype = "dashed") +
#   theme_bw() +
#   theme(
#     axis.text = element_text(size = 20),
#     axis.title = element_text(size = 20),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#     legend.text = element_text(size = 20),
#     legend.title = element_text(size = 20)
#   )
# 
# print(p1)
# 
# ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/PercentPhylum_Log10_0replaced_Perc_0.1_PerSite_Fraction.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
# ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/PercentPhylum_Log10_0replaced_Perc_0.1_PerSite_Fraction.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
# dev.off()
# 
