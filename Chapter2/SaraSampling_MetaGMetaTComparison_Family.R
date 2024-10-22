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

source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/Generics.R")
source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/Pathnames.R")
source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/phyloseqExtended.R")

#==========================================#
#            Import all datasets           #
#==========================================#

# Load the metadata
#------------------

# Import the correspondance between SampleIDs and the AGTUIDs from Fasteris # /!\ I changed the Correspondance file, check everything is right

DF_IDCorrespondance <- read.csv("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/CorrespondenceMetaTMetaG/CorrespondanceRealNameFasterisIndex.csv")
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

# Load the Overview dataframes, and for each remove the percentage signs
#-----------------------------------------------------------------------

# Confidence 0.0 

DF_MetaG_Overview_Conf_0.0 <- read.csv("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SortiesKraken2/Confidence_0.0/Analysis/Overview/SaraSampling_metaT_Kraken2_Confidence_0.0_Overview.csv")

DF_MetaG_Overview_Conf_0.0_woPerc <- DF_MetaG_Overview_Conf_0.0

for (colID in c("Classified.reads", "Chordate.reads", "Artificial.reads", "Unclassified.reads", "Microbial.reads", "Bacterial.reads", "Viral.reads", "Fungal.reads", "Protozoan.reads")) {
  DF_MetaG_Overview_Conf_0.0_woPerc[[colID]] <- as.numeric(substr(DF_MetaG_Overview_Conf_0.0_woPerc[[colID]], 1, nchar(DF_MetaG_Overview_Conf_0.0_woPerc[[colID]])-1))
}

DF_MetaG_Overview_Conf_0.0_woPerc$Confidence <- "0.0"
DF_MetaG_Overview_Conf_0.0_woPerc$AGTUIndex <- substr(DF_MetaG_Overview_Conf_0.0_woPerc$Name, 22, nchar(DF_MetaG_Overview_Conf_0.0_woPerc$Name))

# Confidence 0.1

DF_MetaG_Overview_Conf_0.1 <- read.csv("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SortiesKraken2/Confidence_0.1/Analysis/Overview/SaraSampling_metaT_Kraken2_Confidence_0.1_Overview.csv")

DF_MetaG_Overview_Conf_0.1_woPerc <- DF_MetaG_Overview_Conf_0.1

for (colID in c("Classified.reads", "Chordate.reads", "Artificial.reads", "Unclassified.reads", "Microbial.reads", "Bacterial.reads", "Viral.reads", "Fungal.reads", "Protozoan.reads")) {
  DF_MetaG_Overview_Conf_0.1_woPerc[[colID]] <- as.numeric(substr(DF_MetaG_Overview_Conf_0.1_woPerc[[colID]], 1, nchar(DF_MetaG_Overview_Conf_0.1_woPerc[[colID]])-1))
}

DF_MetaG_Overview_Conf_0.1_woPerc$Confidence <- "0.1"
DF_MetaG_Overview_Conf_0.1_woPerc$AGTUIndex <- substr(DF_MetaG_Overview_Conf_0.1_woPerc$Name, 22, nchar(DF_MetaG_Overview_Conf_0.1_woPerc$Name))

# Confidence 0.6

DF_MetaG_Overview_Conf_0.6 <- read.csv("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SortiesKraken2/Confidence_0.6/Analysis/Overview/SaraSampling_metaT_Kraken2_Confidence_0.6_Overview.csv")

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

# Add a column with MetaG or MetaT

DF_MergedMetadataPavian$ColForMetaTG <- substr(DF_MergedMetadataPavian$AGTUIndex, 1, 4)
DF_MergedMetadataPavian$MetaType <- NA
DF_MergedMetadataPavian[which(DF_MergedMetadataPavian[["ColForMetaTG"]] == "AGTU"),"MetaType"] <- "MetaG"
DF_MergedMetadataPavian[which(DF_MergedMetadataPavian[["ColForMetaTG"]] == "AXWN"),"MetaType"] <- "MetaT"

# Remove all samples that do not have metaT (and that are here due to earlier full join)

DF_MergedMetadataPavian_GoodSamples <- DF_MergedMetadataPavian[is.na(DF_MergedMetadataPavian$Name) == FALSE,]

#======================================================================#
# Plot percent of reads assigned according to environmental parameters #
#======================================================================#

ggplot(data=DF_MergedMetadataPavian_GoodSamples, aes(x = SampleShort, y = Fungal.reads)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point(aes(color = MetaType)) +
  #scale_colour_manual(values = rainbow(12)) +
  theme_bw()

ggplot(data=DF_MergedMetadataPavian_GoodSamples, aes(x = SampleShort, y = Classified.reads)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point(aes(color = MetaType)) +
  #scale_colour_manual(values = rainbow(12)) +
  theme_bw()

#======================================================================#
# Import the files with the number of reads attributed to each species #
#======================================================================#

# Import the files
#-----------------

# Conf 0.0

DF_PavianSpecies_Conf_0.0 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SortiesKraken2/Confidence_0.0/Analysis/RawPavianReports/SaraSampling_metaT_Kraken2_Confidence_0.0_Fungi_Family.tsv", sep = "\t", header = TRUE)

colnames(DF_PavianSpecies_Conf_0.0)[5:14] <- gsub("\\.", "-", substr(colnames(DF_PavianSpecies_Conf_0.0)[5:14], 23, nchar(colnames(DF_PavianSpecies_Conf_0.0)[5:14])-11))

# Conf 0.1

DF_PavianSpecies_Conf_0.1 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SortiesKraken2/Confidence_0.1/Analysis/RawPavianReports/SaraSampling_metaT_Kraken2_Confidence_0.1_Fungi_Family.tsv", sep = "\t", header = TRUE)

colnames(DF_PavianSpecies_Conf_0.1)[5:14] <- gsub("\\.", "-", substr(colnames(DF_PavianSpecies_Conf_0.1)[5:14], 23, nchar(colnames(DF_PavianSpecies_Conf_0.1)[5:14])-11))

# Conf 0.6

DF_PavianSpecies_Conf_0.6 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SortiesKraken2/Confidence_0.6/Analysis/RawPavianReports/SaraSampling_metaT_Kraken2_Confidence_0.6_Fungi_Family.tsv", sep = "\t", header = TRUE)

colnames(DF_PavianSpecies_Conf_0.6)[5:14] <- gsub("\\.", "-", substr(colnames(DF_PavianSpecies_Conf_0.6)[5:14], 23, nchar(colnames(DF_PavianSpecies_Conf_0.6)[5:14])-11))

# Normalize read counts by the number of reads in the sample
#-----------------------------------------------------------

# Conf 0.0

DF_PavianSpecies_Conf_0.0_Norm <- DF_PavianSpecies_Conf_0.0

for (i in 5:14) {
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.0_Norm)[i]
  DF_PavianSpecies_Conf_0.0_Norm[[NameOfCol]] <- DF_PavianSpecies_Conf_0.0_Norm[[NameOfCol]]/DF_MetaG_Overview_Conf_0.0_woPerc[DF_MetaG_Overview_Conf_0.0_woPerc[["AGTUIndex"]] == NameOfCol,"Number.of.raw.reads"]
} # To be checked

# Conf 0.1 

DF_PavianSpecies_Conf_0.1_Norm <- DF_PavianSpecies_Conf_0.1

for (i in 5:14) {
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.1_Norm)[i]
  DF_PavianSpecies_Conf_0.1_Norm[[NameOfCol]] <- DF_PavianSpecies_Conf_0.1_Norm[[NameOfCol]]/DF_MetaG_Overview_Conf_0.1_woPerc[DF_MetaG_Overview_Conf_0.1_woPerc[["AGTUIndex"]] == NameOfCol,"Number.of.raw.reads"]
} # To be checked

# Conf 0.6 

DF_PavianSpecies_Conf_0.6_Norm <- DF_PavianSpecies_Conf_0.6

for (i in 5:14) {
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.6_Norm)[i]
  DF_PavianSpecies_Conf_0.6_Norm[[NameOfCol]] <- DF_PavianSpecies_Conf_0.6_Norm[[NameOfCol]]/DF_MetaG_Overview_Conf_0.6_woPerc[DF_MetaG_Overview_Conf_0.6_woPerc[["AGTUIndex"]] == NameOfCol,"Number.of.raw.reads"]
} # To be checked

# Rename rows according to the species name

row.names(DF_PavianSpecies_Conf_0.0_Norm) <- DF_PavianSpecies_Conf_0.0_Norm$name
row.names(DF_PavianSpecies_Conf_0.1_Norm) <- DF_PavianSpecies_Conf_0.1_Norm$name
row.names(DF_PavianSpecies_Conf_0.6_Norm) <- DF_PavianSpecies_Conf_0.6_Norm$name


# Look at the most abundant species in each sample 
#-------------------------------------------------

# Conf 0.0

ListSorted_Conf_0.0 <- list() # Check this out thoroughly
for (i in 5:14){
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.0_Norm)[i]
  DF_Sorted <- DF_PavianSpecies_Conf_0.0_Norm[order(DF_PavianSpecies_Conf_0.0_Norm[[NameOfCol]], decreasing = TRUE), c("name", NameOfCol)] # Check the ordering of names is good
  colnames(DF_Sorted) <- c(paste("name_", NameOfCol, sep = ""), NameOfCol)
  ListSorted_Conf_0.0[[NameOfCol]] <- DF_Sorted
}

SortedSpecies_Conf_0.0 <- bind_cols(ListSorted_Conf_0.0)
row.names(SortedSpecies_Conf_0.0) <- NULL

# Conf 0.1

ListSorted_Conf_0.1 <- list() # Check this out thoroughly
for (i in 5:14){
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.1_Norm)[i]
  DF_Sorted <- DF_PavianSpecies_Conf_0.1_Norm[order(DF_PavianSpecies_Conf_0.1_Norm[[NameOfCol]], decreasing = TRUE), c("name", NameOfCol)] # Check the ordering of names is good
  colnames(DF_Sorted) <- c(paste("name_", NameOfCol, sep = ""), NameOfCol)
  ListSorted_Conf_0.1[[NameOfCol]] <- DF_Sorted
}

SortedSpecies_Conf_0.1 <- bind_cols(ListSorted_Conf_0.1)
row.names(SortedSpecies_Conf_0.1) <- NULL

# Conf 0.6

ListSorted_Conf_0.6 <- list() # Check this out thoroughly
for (i in 5:14){
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.6_Norm)[i]
  DF_Sorted <- DF_PavianSpecies_Conf_0.6_Norm[order(DF_PavianSpecies_Conf_0.6_Norm[[NameOfCol]], decreasing = TRUE), c("name", NameOfCol)] # Check the ordering of names is good
  colnames(DF_Sorted) <- c(paste("name_", NameOfCol, sep = ""), NameOfCol)
  ListSorted_Conf_0.6[[NameOfCol]] <- DF_Sorted
}

SortedSpecies_Conf_0.6 <- bind_cols(ListSorted_Conf_0.6)
row.names(SortedSpecies_Conf_0.6) <- NULL

#==============================================================================#
#         Check differential relative abundance between metaT and metaG        #
#==============================================================================#

# Divide the relative abundance of each species from the metaT sample by the relative abundance of the metaG sample

# Reminder for the samples correspondance
# SampleNameMetaG_1	MetaG_1	SampleNameMetaG_2	MetaG_2	SampleNameMetaT_1	MetaT_1	SampleNameMetaT_2	MetaT_2
# G_FSS09.1_G1	    AGTU-34	     NA	            NA	     SS09_U1-1	     AXWN-1	    SS09_U1-2	     AXWN-2
# G_FCT10.1_G1	    AGTU-61	   G_FCT10.1_G2	 AGTU-263	   CT10_U2-2	     AXWN-3	    CT10_U2-3	     AXWN-4
# G_FHB10.1_G1	    AGTU-67	     NA	            NA	     HB10_U1-1	     AXWN-5	    HB10_U2-1	     AXWN-6

# Conf 0.0 # Check the calculations do not return any error, especially with zeros
#---------

DF_PavianSpecies_Conf_0.0_Norm_Diff <- DF_PavianSpecies_Conf_0.0_Norm
DF_PavianSpecies_Conf_0.0_Norm_Diff[["SS09DiffExpr_1"]] <- DF_PavianSpecies_Conf_0.0_Norm_Diff[["AXWN-1"]]/DF_PavianSpecies_Conf_0.0_Norm_Diff[["AGTU-34"]]
DF_PavianSpecies_Conf_0.0_Norm_Diff[["SS09DiffExpr_2"]] <- DF_PavianSpecies_Conf_0.0_Norm_Diff[["AXWN-2"]]/DF_PavianSpecies_Conf_0.0_Norm_Diff[["AGTU-34"]]

DF_PavianSpecies_Conf_0.0_Norm_Diff[["CT10DiffExpr_1"]] <- DF_PavianSpecies_Conf_0.0_Norm_Diff[["AXWN-3"]]/DF_PavianSpecies_Conf_0.0_Norm_Diff[["AGTU-263"]]
DF_PavianSpecies_Conf_0.0_Norm_Diff[["CT10DiffExpr_2"]] <- DF_PavianSpecies_Conf_0.0_Norm_Diff[["AXWN-4"]]/DF_PavianSpecies_Conf_0.0_Norm_Diff[["AGTU-263"]]

DF_PavianSpecies_Conf_0.0_Norm_Diff[["HB10DiffExpr_1"]] <- DF_PavianSpecies_Conf_0.0_Norm_Diff[["AXWN-5"]]/DF_PavianSpecies_Conf_0.0_Norm_Diff[["AGTU-67"]]
DF_PavianSpecies_Conf_0.0_Norm_Diff[["HB10DiffExpr_2"]] <- DF_PavianSpecies_Conf_0.0_Norm_Diff[["AXWN-6"]]/DF_PavianSpecies_Conf_0.0_Norm_Diff[["AGTU-67"]]

# Add a column with the sum of the differential expressions to sort on them 

# Replace all NAs by 0

DF_PavianSpecies_Conf_0.0_Norm_Diff <- DF_PavianSpecies_Conf_0.0_Norm_Diff %>% replace(is.na(.), 0)

# Make the column to sort on

DF_PavianSpecies_Conf_0.0_Norm_Diff$ColToSort <- DF_PavianSpecies_Conf_0.0_Norm_Diff[["SS09DiffExpr_1"]] + DF_PavianSpecies_Conf_0.0_Norm_Diff[["SS09DiffExpr_2"]] + DF_PavianSpecies_Conf_0.0_Norm_Diff[["CT10DiffExpr_1"]] + DF_PavianSpecies_Conf_0.0_Norm_Diff[["CT10DiffExpr_2"]] + DF_PavianSpecies_Conf_0.0_Norm_Diff[["HB10DiffExpr_1"]] + DF_PavianSpecies_Conf_0.0_Norm_Diff[["HB10DiffExpr_2"]]

# Order the dataframe according to the new column

DF_PavianSpecies_Conf_0.0_Norm_Diff_Sorted <- DF_PavianSpecies_Conf_0.0_Norm_Diff[order(DF_PavianSpecies_Conf_0.0_Norm_Diff[["ColToSort"]], decreasing = TRUE),]

# Conf 0.1
#---------

DF_PavianSpecies_Conf_0.1_Norm_Diff <- DF_PavianSpecies_Conf_0.1_Norm
DF_PavianSpecies_Conf_0.1_Norm_Diff[["SS09DiffExpr_1"]] <- DF_PavianSpecies_Conf_0.1_Norm_Diff[["AXWN-1"]]/DF_PavianSpecies_Conf_0.1_Norm_Diff[["AGTU-34"]]
DF_PavianSpecies_Conf_0.1_Norm_Diff[["SS09DiffExpr_2"]] <- DF_PavianSpecies_Conf_0.1_Norm_Diff[["AXWN-2"]]/DF_PavianSpecies_Conf_0.1_Norm_Diff[["AGTU-34"]]

DF_PavianSpecies_Conf_0.1_Norm_Diff[["CT10DiffExpr_1"]] <- DF_PavianSpecies_Conf_0.1_Norm_Diff[["AXWN-3"]]/DF_PavianSpecies_Conf_0.1_Norm_Diff[["AGTU-263"]]
DF_PavianSpecies_Conf_0.1_Norm_Diff[["CT10DiffExpr_2"]] <- DF_PavianSpecies_Conf_0.1_Norm_Diff[["AXWN-4"]]/DF_PavianSpecies_Conf_0.1_Norm_Diff[["AGTU-263"]]

DF_PavianSpecies_Conf_0.1_Norm_Diff[["HB10DiffExpr_1"]] <- DF_PavianSpecies_Conf_0.1_Norm_Diff[["AXWN-5"]]/DF_PavianSpecies_Conf_0.1_Norm_Diff[["AGTU-67"]]
DF_PavianSpecies_Conf_0.1_Norm_Diff[["HB10DiffExpr_2"]] <- DF_PavianSpecies_Conf_0.1_Norm_Diff[["AXWN-6"]]/DF_PavianSpecies_Conf_0.1_Norm_Diff[["AGTU-67"]]

# Add a column with the sum of the differential expressions to sort on them 

# Replace all NAs by 0

DF_PavianSpecies_Conf_0.1_Norm_Diff <- DF_PavianSpecies_Conf_0.1_Norm_Diff %>% replace(is.na(.), 0)

# Make the column to sort on

DF_PavianSpecies_Conf_0.1_Norm_Diff$ColToSort <- DF_PavianSpecies_Conf_0.1_Norm_Diff[["SS09DiffExpr_1"]] + DF_PavianSpecies_Conf_0.1_Norm_Diff[["SS09DiffExpr_2"]] + DF_PavianSpecies_Conf_0.1_Norm_Diff[["CT10DiffExpr_1"]] + DF_PavianSpecies_Conf_0.1_Norm_Diff[["CT10DiffExpr_2"]] + DF_PavianSpecies_Conf_0.1_Norm_Diff[["HB10DiffExpr_1"]] + DF_PavianSpecies_Conf_0.1_Norm_Diff[["HB10DiffExpr_2"]]

# Order the dataframe according to the new column

DF_PavianSpecies_Conf_0.1_Norm_Diff_Sorted <- DF_PavianSpecies_Conf_0.1_Norm_Diff[order(DF_PavianSpecies_Conf_0.1_Norm_Diff[["ColToSort"]], decreasing = TRUE),]

# Conf 0.6
#---------

DF_PavianSpecies_Conf_0.6_Norm_Diff <- DF_PavianSpecies_Conf_0.6_Norm
DF_PavianSpecies_Conf_0.6_Norm_Diff[["SS09DiffExpr_1"]] <- DF_PavianSpecies_Conf_0.6_Norm_Diff[["AXWN-1"]]/DF_PavianSpecies_Conf_0.6_Norm_Diff[["AGTU-34"]]
DF_PavianSpecies_Conf_0.6_Norm_Diff[["SS09DiffExpr_2"]] <- DF_PavianSpecies_Conf_0.6_Norm_Diff[["AXWN-2"]]/DF_PavianSpecies_Conf_0.6_Norm_Diff[["AGTU-34"]]

DF_PavianSpecies_Conf_0.6_Norm_Diff[["CT10DiffExpr_1"]] <- DF_PavianSpecies_Conf_0.6_Norm_Diff[["AXWN-3"]]/DF_PavianSpecies_Conf_0.6_Norm_Diff[["AGTU-263"]]
DF_PavianSpecies_Conf_0.6_Norm_Diff[["CT10DiffExpr_2"]] <- DF_PavianSpecies_Conf_0.6_Norm_Diff[["AXWN-4"]]/DF_PavianSpecies_Conf_0.6_Norm_Diff[["AGTU-263"]]

DF_PavianSpecies_Conf_0.6_Norm_Diff[["HB10DiffExpr_1"]] <- DF_PavianSpecies_Conf_0.6_Norm_Diff[["AXWN-5"]]/DF_PavianSpecies_Conf_0.6_Norm_Diff[["AGTU-67"]]
DF_PavianSpecies_Conf_0.6_Norm_Diff[["HB10DiffExpr_2"]] <- DF_PavianSpecies_Conf_0.6_Norm_Diff[["AXWN-6"]]/DF_PavianSpecies_Conf_0.6_Norm_Diff[["AGTU-67"]]

# Add a column with the sum of the differential expressions to sort on them 

# Replace all NAs by 0

DF_PavianSpecies_Conf_0.6_Norm_Diff <- DF_PavianSpecies_Conf_0.6_Norm_Diff %>% replace(is.na(.), 0)

# Make the column to sort on

DF_PavianSpecies_Conf_0.6_Norm_Diff$ColToSort <- DF_PavianSpecies_Conf_0.6_Norm_Diff[["SS09DiffExpr_1"]] + DF_PavianSpecies_Conf_0.6_Norm_Diff[["SS09DiffExpr_2"]] + DF_PavianSpecies_Conf_0.6_Norm_Diff[["CT10DiffExpr_1"]] + DF_PavianSpecies_Conf_0.6_Norm_Diff[["CT10DiffExpr_2"]] + DF_PavianSpecies_Conf_0.6_Norm_Diff[["HB10DiffExpr_1"]] + DF_PavianSpecies_Conf_0.6_Norm_Diff[["HB10DiffExpr_2"]]

# Order the dataframe according to the new column

DF_PavianSpecies_Conf_0.6_Norm_Diff_Sorted <- DF_PavianSpecies_Conf_0.6_Norm_Diff[order(DF_PavianSpecies_Conf_0.6_Norm_Diff[["ColToSort"]], decreasing = TRUE),]

