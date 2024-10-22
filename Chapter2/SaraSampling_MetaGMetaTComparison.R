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

p1 <- ggplot(data=DF_MergedMetadataPavian_GoodSamples, aes(x = SampleShort, y = Fungal.reads)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point(aes(color = MetaType)) +
  #scale_colour_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/ComparisonMetaGMetaT/FigMetaTMetaG_FungalReads.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/ComparisonMetaGMetaT/FigMetaTMetaG_FungalReads.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_GoodSamples, aes(x = SampleShort, y = Classified.reads)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point(aes(color = MetaType)) +
  #scale_colour_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/ComparisonMetaGMetaT/FigMetaTMetaG_ClassifiedReads.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/ComparisonMetaGMetaT/FigMetaTMetaG_ClassifiedReads.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

#======================================================================#
# Import the files with the number of reads attributed to each species #
#======================================================================#

# Import the files
#-----------------

# Conf 0.0

DF_PavianSpecies_Conf_0.0 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SortiesKraken2/Confidence_0.0/Analysis/RawPavianReports/SaraSampling_metaT_Kraken2_Confidence_0.0_Fungi_Species.tsv", sep = "\t", header = TRUE)

colnames(DF_PavianSpecies_Conf_0.0)[5:14] <- gsub("\\.", "-", substr(colnames(DF_PavianSpecies_Conf_0.0)[5:14], 23, nchar(colnames(DF_PavianSpecies_Conf_0.0)[5:14])-11))

# Conf 0.1

DF_PavianSpecies_Conf_0.1 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SortiesKraken2/Confidence_0.1/Analysis/RawPavianReports/SaraSampling_metaT_Kraken2_Confidence_0.1_Fungi_Species.tsv", sep = "\t", header = TRUE)

colnames(DF_PavianSpecies_Conf_0.1)[5:14] <- gsub("\\.", "-", substr(colnames(DF_PavianSpecies_Conf_0.1)[5:14], 23, nchar(colnames(DF_PavianSpecies_Conf_0.1)[5:14])-11))

# Conf 0.6

DF_PavianSpecies_Conf_0.6 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SortiesKraken2/Confidence_0.6/Analysis/RawPavianReports/SaraSampling_metaT_Kraken2_Confidence_0.6_Fungi_Species.tsv", sep = "\t", header = TRUE)

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

write.csv(SortedSpecies_Conf_0.0, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/ComparisonMetaGMetaT/Species/RelativeAbundance/RelativeAbundanceSpecies_Conf_0.0.csv", row.names = FALSE)

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

write.csv(SortedSpecies_Conf_0.1, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/ComparisonMetaGMetaT/Species/RelativeAbundance/RelativeAbundanceSpecies_Conf_0.1.csv", row.names = FALSE)

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

write.csv(SortedSpecies_Conf_0.6, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/ComparisonMetaGMetaT/Species/RelativeAbundance/RelativeAbundanceSpecies_Conf_0.6.csv", row.names = FALSE)

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

# Added on 2024-06-27 to calculate the number of species most abundant in metaT compared to metaG

DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount <- DF_PavianSpecies_Conf_0.1_Norm_Diff

SupMetaTMetaG_0.1_SS09_1 <- DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-1"]] > DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AGTU-34"]]), ] # Check this
SupMetaTMetaG_0.1_SS09_2 <- DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-2"]] > DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AGTU-34"]]), ] # Check this

SupMetaTMetaG_0.1_SS09_Names <- row.names(SupMetaTMetaG_0.1_SS09_1)[which(row.names(SupMetaTMetaG_0.1_SS09_1) %in% row.names(SupMetaTMetaG_0.1_SS09_2))]
length(SupMetaTMetaG_0.1_SS09_Names)

RowNames_0.1_SS09_MetaT1 <- row.names(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount)[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-1"]] > 0)]
RowNames_0.1_SS09_MetaT2 <- row.names(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount)[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-2"]] > 0)]
RowNames_0.1_SS09_MetaG <- row.names(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount)[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AGTU-34"]] > 0)]

length(unique(c(RowNames_0.1_SS09_MetaT1, RowNames_0.1_SS09_MetaT2, RowNames_0.1_SS09_MetaG))) # Check

SupMetaTMetaG_0.1_CT10_1 <- DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-3"]] > DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AGTU-263"]]), ] # Check this
SupMetaTMetaG_0.1_CT10_2 <- DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-4"]] > DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AGTU-263"]]), ] # Check this

SupMetaTMetaG_0.1_CT10_Names <- row.names(SupMetaTMetaG_0.1_CT10_1)[which(row.names(SupMetaTMetaG_0.1_CT10_1) %in% row.names(SupMetaTMetaG_0.1_CT10_2))]
length(SupMetaTMetaG_0.1_CT10_Names)

RowNames_0.1_CT10_MetaT1 <- row.names(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount)[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-3"]] > 0)]
RowNames_0.1_CT10_MetaT2 <- row.names(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount)[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-4"]] > 0)]
RowNames_0.1_CT10_MetaG <- row.names(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount)[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AGTU-263"]] > 0)]

length(unique(c(RowNames_0.1_CT10_MetaT1, RowNames_0.1_CT10_MetaT2, RowNames_0.1_CT10_MetaG))) # Check

SupMetaTMetaG_0.1_HB10_1 <- DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-5"]] > DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AGTU-67"]]), ] # Check this
SupMetaTMetaG_0.1_HB10_2 <- DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-6"]] > DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AGTU-67"]]), ] # Check this

SupMetaTMetaG_0.1_HB10_1_Names <- row.names(SupMetaTMetaG_0.1_HB10_1)[which(row.names(SupMetaTMetaG_0.1_HB10_1) %in% row.names(SupMetaTMetaG_0.1_HB10_2))]
length(SupMetaTMetaG_0.1_HB10_1_Names)

RowNames_0.1_HB10_MetaT1 <- row.names(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount)[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-5"]] > 0)]
RowNames_0.1_HB10_MetaT2 <- row.names(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount)[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AXWN-6"]] > 0)]
RowNames_0.1_HB10_MetaG <- row.names(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount)[which(DF_PavianSpecies_Conf_0.1_Norm_Diff_ForCount[["AGTU-67"]] > 0)]

length(unique(c(RowNames_0.1_HB10_MetaT1, RowNames_0.1_HB10_MetaT2, RowNames_0.1_HB10_MetaG))) # Check

# End of added on 2024-06-27 to calculate the number of species most abundant in metaT compared to metaG

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

#-------------------------------------------------------------#
# List most differentially expressed species in for all sites #
#-------------------------------------------------------------#

# Conf 0.0

ListSorted_Conf_0.0 <- list() # Check this out thoroughly
for (i in 16:21){
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.0_Norm_Diff)[i]
  DF_Sorted <- DF_PavianSpecies_Conf_0.0_Norm_Diff[order(DF_PavianSpecies_Conf_0.0_Norm_Diff[[NameOfCol]], decreasing = TRUE), c("name", NameOfCol)] # Check the ordering of names is good
  colnames(DF_Sorted) <- c(paste("name_", NameOfCol, sep = ""), NameOfCol)
  ListSorted_Conf_0.0[[NameOfCol]] <- DF_Sorted
}

SortedSpecies_Conf_0.0_DiffRelAb <- bind_cols(ListSorted_Conf_0.0)
row.names(SortedSpecies_Conf_0.0_DiffRelAb) <- NULL

write.csv(SortedSpecies_Conf_0.0_DiffRelAb, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/ComparisonMetaGMetaT/Species/DifferentialExpression/MostDifferentialyAbundantSpecies_Conf_0.0.csv", row.names = FALSE)

# Conf 0.1

ListSorted_Conf_0.1 <- list() # Check this out thoroughly
for (i in 16:21){
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.1_Norm_Diff)[i]
  DF_Sorted <- DF_PavianSpecies_Conf_0.1_Norm_Diff[order(DF_PavianSpecies_Conf_0.1_Norm_Diff[[NameOfCol]], decreasing = TRUE), c("name", NameOfCol)] # Check the ordering of names is good
  colnames(DF_Sorted) <- c(paste("name_", NameOfCol, sep = ""), NameOfCol)
  ListSorted_Conf_0.1[[NameOfCol]] <- DF_Sorted
}

SortedSpecies_Conf_0.1_DiffRelAb <- bind_cols(ListSorted_Conf_0.1)
row.names(SortedSpecies_Conf_0.1_DiffRelAb) <- NULL

write.csv(SortedSpecies_Conf_0.1_DiffRelAb, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/ComparisonMetaGMetaT/Species/DifferentialExpression/MostDifferentialyAbundantSpecies_Conf_0.1.csv", row.names = FALSE)

# Conf 0.6

ListSorted_Conf_0.6 <- list() # Check this out thoroughly
for (i in 16:21){
  NameOfCol <- colnames(DF_PavianSpecies_Conf_0.6_Norm_Diff)[i]
  DF_Sorted <- DF_PavianSpecies_Conf_0.6_Norm_Diff[order(DF_PavianSpecies_Conf_0.6_Norm_Diff[[NameOfCol]], decreasing = TRUE), c("name", NameOfCol)] # Check the ordering of names is good
  colnames(DF_Sorted) <- c(paste("name_", NameOfCol, sep = ""), NameOfCol)
  ListSorted_Conf_0.6[[NameOfCol]] <- DF_Sorted
}

SortedSpecies_Conf_0.6_DiffRelAb <- bind_cols(ListSorted_Conf_0.6)
row.names(SortedSpecies_Conf_0.6_DiffRelAb) <- NULL

write.csv(SortedSpecies_Conf_0.6_DiffRelAb, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/ComparisonMetaGMetaT/Species/DifferentialExpression/MostDifferentialyAbundantSpecies_Conf_0.6.csv", row.names = FALSE)

#======================================================================================================================================================#
# Make a plot of the species that are present in both the metaG and the metaT, list the species that are present only in the metaG, only in the metaT, #
# and 
#======================================================================================================================================================#

# Confidence 0.1 - Version with complete cases

# Remove rows that have NA

DF_PavianSpecies_Conf_0.1_Norm_FSS09 <- DF_PavianSpecies_Conf_0.1_Norm[, c("AGTU-34", "AXWN-1")]
DF_PavianSpecies_Conf_0.1_Norm_FSS09 <- DF_PavianSpecies_Conf_0.1_Norm_FSS09[complete.cases(DF_PavianSpecies_Conf_0.1_Norm_FSS09),]

# SS09 A

ggplot(data = DF_PavianSpecies_Conf_0.1_Norm_FSS09, aes(x = `AGTU-34`, y = `AXWN-1`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black", linewidth=0.5) +
  geom_vline(xintercept = 0.00001) +
  geom_hline(yintercept = 0.00001) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  xlab(label = "FSS09 MetaG") +
  ylab(label = "FSS09 MetaT") +
  theme_bw()

ggplot(data = DF_PavianSpecies_Conf_0.1_Norm_FSS09, aes(x = `AGTU-34`, y = `AXWN-1`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black", linewidth=0.5) +
  geom_vline(xintercept = 0.00001) +
  geom_hline(yintercept = 0.00001) +
  xlab(label = "FSS09 MetaG") +
  ylab(label = "FSS09 MetaT") +
  theme_bw()

# FCT10

DF_PavianSpecies_Conf_0.1_Norm_FCT10 <- DF_PavianSpecies_Conf_0.1_Norm[, c("AGTU-263", "AXWN-3")]
DF_PavianSpecies_Conf_0.1_Norm_FCT10 <- DF_PavianSpecies_Conf_0.1_Norm_FCT10[complete.cases(DF_PavianSpecies_Conf_0.1_Norm_FCT10),]

ggplot(data = DF_PavianSpecies_Conf_0.1_Norm_FCT10, aes(x = `AGTU-263`, y = `AXWN-3`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black", linewidth=0.5) +
  geom_vline(xintercept = 0.00001) +
  geom_hline(yintercept = 0.00001) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  xlab(label = "FCT10 MetaG") +
  ylab(label = "FCT10 MetaT") +
  theme_bw()

ggplot(data = DF_PavianSpecies_Conf_0.1_Norm_FCT10, aes(x = `AGTU-263`, y = `AXWN-3`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black", linewidth=0.5) +
  geom_vline(xintercept = 0.00001) +
  geom_hline(yintercept = 0.00001) +
  xlab(label = "FCT10 MetaG") +
  ylab(label = "FCT10 MetaT") +
  theme_bw()

# HB10

DF_PavianSpecies_Conf_0.1_Norm_FHB10 <- DF_PavianSpecies_Conf_0.1_Norm[, c("AGTU-67", "AXWN-5")]
DF_PavianSpecies_Conf_0.1_Norm_FHB10 <- DF_PavianSpecies_Conf_0.1_Norm_FHB10[complete.cases(DF_PavianSpecies_Conf_0.1_Norm_FHB10),]

ggplot(data = DF_PavianSpecies_Conf_0.1_Norm_FHB10, aes(x = `AGTU-67`, y = `AXWN-5`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black", linewidth=0.5) +
  geom_vline(xintercept = 0.00001) +
  geom_hline(yintercept = 0.00001) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  xlab(label = "FHB10 MetaG") +
  ylab(label = "FHB10 MetaT") +
  theme_bw()

ggplot(data = DF_PavianSpecies_Conf_0.1_Norm_FHB10, aes(x = `AGTU-67`, y = `AXWN-5`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black", linewidth=0.5) +
  geom_vline(xintercept = 0.00001) +
  geom_hline(yintercept = 0.00001) +
  xlab(label = "FHB10 MetaG") +
  ylab(label = "FHB10 MetaT") +
  theme_bw()

# Confidence 0.1 - Version with the two replicates on the same plot

# SS09

DF_PavianSpecies_Conf_0.1_Norm_FSS09 <- DF_PavianSpecies_Conf_0.1_Norm[, c("AGTU-34", "AXWN-1", "AXWN-2")]
DF_PavianSpecies_Conf_0.1_Norm_FSS09_long <- DF_PavianSpecies_Conf_0.1_Norm_FSS09 %>%
  pivot_longer(
    c("AXWN-1", "AXWN-2"),
    names_to = c("metaT"),
    values_to = "Count"
  )

ggplot(data = DF_PavianSpecies_Conf_0.1_Norm_FSS09_long, aes(x = `AGTU-34`, y = Count, color=metaT)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black", linewidth=0.5) +
  geom_vline(xintercept = 0.00001) +
  geom_hline(yintercept = 0.00001) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  xlab(label = "FSS09 MetaG") +
  ylab(label = "FSS09 MetaT") +
  theme_bw()

ggplot(data = DF_PavianSpecies_Conf_0.1_Norm_FSS09_long, aes(x = `AGTU-34`, y = Count, color=metaT)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black", linewidth=0.5) +
  geom_vline(xintercept = 0.00001) +
  geom_hline(yintercept = 0.00001) +
  xlab(label = "FSS09 MetaG") +
  ylab(label = "FSS09 MetaT") +
  theme_bw()

# With three panels
#------------------

# 0.0

# Select good columns 
DF_PavianSpecies_Conf_0.0_Norm_ThreePan <- DF_PavianSpecies_Conf_0.0_Norm
DF_PavianSpecies_Conf_0.0_Norm_ThreePan <- DF_PavianSpecies_Conf_0.0_Norm_ThreePan[, c("name", "AGTU-263", "AGTU-34", "AGTU-67", "AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6")]

# Find the minimum, to set all NA below the lowest value
min(sapply(DF_PavianSpecies_Conf_0.0_Norm_ThreePan[, c("AGTU-263", "AGTU-34", "AGTU-67", "AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6")], function(x) min(x, na.rm = TRUE))) # The min is 2.140426e-08

# Check there is no NA in names, then set all NA to 10-9 (one ordre of magnitude below the lowest value)
which(is.na(DF_PavianSpecies_Conf_0.0_Norm_ThreePan$name)) # Should be integer(0)

DF_PavianSpecies_Conf_0.0_Norm_ThreePan <- DF_PavianSpecies_Conf_0.0_Norm_ThreePan %>% replace(is.na(.), 0.00000001)

# Pivot to long, then add columns to make groups

DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long <- DF_PavianSpecies_Conf_0.0_Norm_ThreePan %>%
  pivot_longer(
    c("AGTU-263", "AGTU-34", "AGTU-67", "AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6"),
    names_to = "Sample",
    values_to = "RelAb"
  )

# Add columns to mention if the sample is metaG or metaT and from which site it is

DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long$Site <- NA
DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-34", "AXWN-1", "AXWN-2"), "Site"] <- "FSS09"
DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-263", "AXWN-3", "AXWN-4"), "Site"] <- "FCT10"
DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-67", "AXWN-5", "AXWN-6"), "Site"] <- "FHB10"

DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long$TypeSam <- NA
DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-34", "AGTU-263", "AGTU-67"), "TypeSam"] <- "MetaG"
DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[["Sample"]] %in% c("AXWN-1", "AXWN-3", "AXWN-5"), "TypeSam"] <- "MetaT_Rep1"
DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long[["Sample"]] %in% c("AXWN-2", "AXWN-4", "AXWN-6"), "TypeSam"] <- "MetaT_Rep2"

DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long$Sample <- NULL
DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long_Wide <- DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long %>%
  pivot_wider(names_from = TypeSam, values_from = RelAb)

DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long_Wide_Long <- DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long_Wide %>%
  pivot_longer(
    c("MetaT_Rep1", "MetaT_Rep2"),
    names_to = "MetaT_Rep",
    values_to = "MetaT"
  )

p1 <- ggplot(DF_PavianSpecies_Conf_0.0_Norm_ThreePan_Long_Wide_Long, aes(x = MetaG, y = MetaT, color = MetaT_Rep)) +
  facet_grid(cols = vars(Site)) +
  geom_point(shape = 1) +
  geom_abline(intercept = 0, slope = 1, color="black", linewidth=0.5) +
  geom_vline(xintercept = 0.00001, color = "chartreuse3") +
  geom_hline(yintercept = 0.00001, color = "chartreuse3") +
  geom_vline(xintercept = 0.00000002, linetype="dotted") +
  geom_hline(yintercept = 0.00000002, linetype="dotted") +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  annotation_logticks(sides = "bl") +
  xlab(label = "Relative abundance in metagenomics") +
  ylab(label = "Relative abundance in metatranscriptomics") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text.x = element_text(size = 15)
  )

plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/MetaT_MetaG_Species_RelAb0.0.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/MetaT_MetaG_Species_RelAb0.0.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

# 0.1

# Select good columns 
DF_PavianSpecies_Conf_0.1_Norm_ThreePan <- DF_PavianSpecies_Conf_0.1_Norm
DF_PavianSpecies_Conf_0.1_Norm_ThreePan <- DF_PavianSpecies_Conf_0.1_Norm_ThreePan[, c("name", "AGTU-263", "AGTU-34", "AGTU-67", "AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6")]

# Find the minimum, to set all NA below the lowest value
min(sapply(DF_PavianSpecies_Conf_0.1_Norm_ThreePan[, c("AGTU-263", "AGTU-34", "AGTU-67", "AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6")], function(x) min(x, na.rm = TRUE))) # The min is 2.140426e-08

# Check there is no NA in names, then set all NA to 10-9 (one ordre of magnitude below the lowest value)
which(is.na(DF_PavianSpecies_Conf_0.1_Norm_ThreePan$name)) # Should be integer(0)

DF_PavianSpecies_Conf_0.1_Norm_ThreePan <- DF_PavianSpecies_Conf_0.1_Norm_ThreePan %>% replace(is.na(.), 0.00000001)

# Pivot to long, then add columns to make groups

DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long <- DF_PavianSpecies_Conf_0.1_Norm_ThreePan %>%
  pivot_longer(
    c("AGTU-263", "AGTU-34", "AGTU-67", "AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6"),
    names_to = "Sample",
    values_to = "RelAb"
  )

# Add columns to mention if the sample is metaG or metaT and from which site it is

DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long$Site <- NA
DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-34", "AXWN-1", "AXWN-2"), "Site"] <- "FSS09"
DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-263", "AXWN-3", "AXWN-4"), "Site"] <- "FCT10"
DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-67", "AXWN-5", "AXWN-6"), "Site"] <- "FHB10"

DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long$TypeSam <- NA
DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-34", "AGTU-263", "AGTU-67"), "TypeSam"] <- "MetaG"
DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[["Sample"]] %in% c("AXWN-1", "AXWN-3", "AXWN-5"), "TypeSam"] <- "MetaT_Rep1"
DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long[["Sample"]] %in% c("AXWN-2", "AXWN-4", "AXWN-6"), "TypeSam"] <- "MetaT_Rep2"

DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long$Sample <- NULL
DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long_Wide <- DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long %>%
  pivot_wider(names_from = TypeSam, values_from = RelAb)

DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long_Wide_Long <- DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long_Wide %>%
  pivot_longer(
    c("MetaT_Rep1", "MetaT_Rep2"),
    names_to = "MetaT_Rep",
    values_to = "MetaT"
  )

p1 <- ggplot(DF_PavianSpecies_Conf_0.1_Norm_ThreePan_Long_Wide_Long, aes(x = MetaG, y = MetaT, color = MetaT_Rep)) +
  facet_grid(cols = vars(Site)) +
  geom_point(shape = 1) +
  geom_abline(intercept = 0, slope = 1, color="black", linewidth=0.5) +
  geom_vline(xintercept = 0.00001, color = "chartreuse3") +
  geom_hline(yintercept = 0.00001, color = "chartreuse3") +
  geom_vline(xintercept = 0.00000002, linetype="dotted") +
  geom_hline(yintercept = 0.00000002, linetype="dotted") +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  annotation_logticks(sides = "bl") +
  xlab(label = "Relative abundance in metagenomics") +
  ylab(label = "Relative abundance in metatranscriptomics") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text.x = element_text(size = 15)
  )

plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/MetaT_MetaG_Species_RelAb0.1.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/MetaT_MetaG_Species_RelAb0.1.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

# 0.6

# Select good columns 
DF_PavianSpecies_Conf_0.6_Norm_ThreePan <- DF_PavianSpecies_Conf_0.6_Norm
DF_PavianSpecies_Conf_0.6_Norm_ThreePan <- DF_PavianSpecies_Conf_0.6_Norm_ThreePan[, c("name", "AGTU-263", "AGTU-34", "AGTU-67", "AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6")]

# Find the minimum, to set all NA below the lowest value
min(sapply(DF_PavianSpecies_Conf_0.6_Norm_ThreePan[, c("AGTU-263", "AGTU-34", "AGTU-67", "AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6")], function(x) min(x, na.rm = TRUE))) # The min is 2.140426e-08

# Check there is no NA in names, then set all NA to 10-9 (one ordre of magnitude below the lowest value)
which(is.na(DF_PavianSpecies_Conf_0.6_Norm_ThreePan$name)) # Should be integer(0)

DF_PavianSpecies_Conf_0.6_Norm_ThreePan <- DF_PavianSpecies_Conf_0.6_Norm_ThreePan %>% replace(is.na(.), 0.00000001)

# Pivot to long, then add columns to make groups

DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long <- DF_PavianSpecies_Conf_0.6_Norm_ThreePan %>%
  pivot_longer(
    c("AGTU-263", "AGTU-34", "AGTU-67", "AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6"),
    names_to = "Sample",
    values_to = "RelAb"
  )

# Add columns to mention if the sample is metaG or metaT and from which site it is

DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long$Site <- NA
DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-34", "AXWN-1", "AXWN-2"), "Site"] <- "FSS09"
DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-263", "AXWN-3", "AXWN-4"), "Site"] <- "FCT10"
DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-67", "AXWN-5", "AXWN-6"), "Site"] <- "FHB10"

DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long$TypeSam <- NA
DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[["Sample"]] %in% c("AGTU-34", "AGTU-263", "AGTU-67"), "TypeSam"] <- "MetaG"
DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[["Sample"]] %in% c("AXWN-1", "AXWN-3", "AXWN-5"), "TypeSam"] <- "MetaT_Rep1"
DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long[["Sample"]] %in% c("AXWN-2", "AXWN-4", "AXWN-6"), "TypeSam"] <- "MetaT_Rep2"

DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long$Sample <- NULL
DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long_Wide <- DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long %>%
  pivot_wider(names_from = TypeSam, values_from = RelAb)

DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long_Wide_Long <- DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long_Wide %>%
  pivot_longer(
    c("MetaT_Rep1", "MetaT_Rep2"),
    names_to = "MetaT_Rep",
    values_to = "MetaT"
  )

p1 <- ggplot(DF_PavianSpecies_Conf_0.6_Norm_ThreePan_Long_Wide_Long, aes(x = MetaG, y = MetaT, color = MetaT_Rep)) +
  facet_grid(cols = vars(Site)) +
  geom_point(shape = 1) +
  geom_abline(intercept = 0, slope = 1, color="black", linewidth=0.5) +
  geom_vline(xintercept = 0.00001, color = "chartreuse3") +
  geom_hline(yintercept = 0.00001, color = "chartreuse3") +
  geom_vline(xintercept = 0.00000002, linetype="dotted") +
  geom_hline(yintercept = 0.00000002, linetype="dotted") +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  annotation_logticks(sides = "bl") +
  xlab(label = "Relative abundance in metagenomics") +
  ylab(label = "Relative abundance in metatranscriptomics") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text.x = element_text(size = 15)
  )

plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/MetaT_MetaG_Species_RelAb0.6.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/MetaT_MetaG_Species_RelAb0.6.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

#----------------------------------------------------------------------------------------------#
# Check how many species are above the thresholds confidence 0.6 and relative abundance > 10-5 #
#----------------------------------------------------------------------------------------------#

ThreshSpec <- 0.00001

# MetaT
#------

# Conf 0.0

DF_PavianSpecies_Conf_0.0_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.0_Norm
DF_PavianSpecies_Conf_0.0_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.0_Norm_AboveThresh[, c("AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6")]

ListSpecAboveThresh_0.0 <- list()
for (i in (1:(ncol(DF_PavianSpecies_Conf_0.0_Norm_AboveThresh)))) {
  TheColName <- colnames(DF_PavianSpecies_Conf_0.0_Norm_AboveThresh)[i]
  SpecCounts <- DF_PavianSpecies_Conf_0.0_Norm_AboveThresh[[TheColName]][which(DF_PavianSpecies_Conf_0.0_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  SpecNames <- row.names(DF_PavianSpecies_Conf_0.0_Norm_AboveThresh)[which(DF_PavianSpecies_Conf_0.0_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  names(SpecCounts) <- SpecNames
  ListSpecAboveThresh_0.0[[TheColName]] <- SpecCounts
} # Check thoroughly

ListSpecAboveThresh_0.0

unlist(sapply(ListSpecAboveThresh_0.0, function(x) names(x)))
unique(unlist(sapply(ListSpecAboveThresh_0.0, function(x) names(x))))
sort(unique(unlist(sapply(ListSpecAboveThresh_0.0, function(x) names(x)))))

DF_ListSpecAboveThresh_0.0 <- sort(unique(unlist(sapply(ListSpecAboveThresh_0.0, function(x) names(x)))))
DF_ListSpecAboveThresh_0.0 <- data.frame(ListSpecies = DF_ListSpecAboveThresh_0.0)
DF_ListSpecAboveThresh_0.0$MetaT_0.0 <- "Yes"

# Conf 0.1

DF_PavianSpecies_Conf_0.1_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.1_Norm
DF_PavianSpecies_Conf_0.1_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.1_Norm_AboveThresh[, c("AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6")]

ListSpecAboveThresh_0.1 <- list()
for (i in (1:(ncol(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh)))) {
  TheColName <- colnames(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh)[i]
  SpecCounts <- DF_PavianSpecies_Conf_0.1_Norm_AboveThresh[[TheColName]][which(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  SpecNames <- row.names(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh)[which(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  names(SpecCounts) <- SpecNames
  ListSpecAboveThresh_0.1[[TheColName]] <- SpecCounts
} # Check thoroughly

ListSpecAboveThresh_0.1

# The list is : Malassezia restricta, Mortierella globalpina, Alternaria alternata, Malassezia globosa, Mortierella globalpina, Aureobasidium pullulans, uncultured fungus (not more precision, only that it is a Fungi)

unlist(sapply(ListSpecAboveThresh_0.1, function(x) names(x)))
unique(unlist(sapply(ListSpecAboveThresh_0.1, function(x) names(x))))
sort(unique(unlist(sapply(ListSpecAboveThresh_0.1, function(x) names(x)))))

DF_ListSpecAboveThresh_0.1 <- sort(unique(unlist(sapply(ListSpecAboveThresh_0.1, function(x) names(x)))))
DF_ListSpecAboveThresh_0.1 <- data.frame(ListSpecies = DF_ListSpecAboveThresh_0.1)
DF_ListSpecAboveThresh_0.1$MetaT_0.1 <- "Yes"

# Conf 0.6

DF_PavianSpecies_Conf_0.6_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.6_Norm
DF_PavianSpecies_Conf_0.6_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.6_Norm_AboveThresh[, c("AXWN-1", "AXWN-2", "AXWN-3", "AXWN-4", "AXWN-5", "AXWN-6")]

ListSpecAboveThresh_0.6 <- list()
for (i in (1:(ncol(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh)))) {
  TheColName <- colnames(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh)[i]
  SpecCounts <- DF_PavianSpecies_Conf_0.6_Norm_AboveThresh[[TheColName]][which(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  SpecNames <- row.names(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh)[which(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  names(SpecCounts) <- SpecNames
  ListSpecAboveThresh_0.6[[TheColName]] <- SpecCounts
} # Check thoroughly

ListSpecAboveThresh_0.6
# The only Fungi which passes the threshold is a Mucoromycota (Mortierella globalpina)

unlist(sapply(ListSpecAboveThresh_0.6, function(x) names(x)))
unique(unlist(sapply(ListSpecAboveThresh_0.6, function(x) names(x))))
sort(unique(unlist(sapply(ListSpecAboveThresh_0.6, function(x) names(x)))))

DF_ListSpecAboveThresh_0.6 <- sort(unique(unlist(sapply(ListSpecAboveThresh_0.6, function(x) names(x)))))
DF_ListSpecAboveThresh_0.6 <- data.frame(ListSpecies = DF_ListSpecAboveThresh_0.6)
DF_ListSpecAboveThresh_0.6$MetaT_0.6 <- "Yes"

# Merge all three dataframes 

DF_ListSpecAboveThresh_AllConf_MetaT <- list(DF_ListSpecAboveThresh_0.0, DF_ListSpecAboveThresh_0.1, DF_ListSpecAboveThresh_0.6) %>% reduce(full_join) # To be checked

DF_ListSpecAboveThresh_AllConf_MetaG <- read.csv("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SpeciesAboveThreshold/ListSpecAboveThresh_AllConf_MetaG.csv")

DF_ListSpecAboveThresh_AllConf_MetaT_MetaG <- list(DF_ListSpecAboveThresh_AllConf_MetaG, DF_ListSpecAboveThresh_AllConf_MetaT) %>% reduce(full_join) # To be checked

# Replace "Yes" by 1, and NA by 0

DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final <- DF_ListSpecAboveThresh_AllConf_MetaT_MetaG
DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final[DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final == "Yes"] <- 1
DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final <- DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final %>% replace(is.na(.), 0)

cols.num <- c("MetaG_0.0", "MetaG_0.1", "MetaG_0.6", "MetaT_0.0", "MetaT_0.1", "MetaT_0.6")
DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final[,cols.num] <- sapply(DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final[,cols.num], as.numeric)

#DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final$RowSumAll <- rowSums(DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final[,cols.num])
DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final_Sorted <- DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final[order(DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final$MetaT_0.6, DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final$MetaT_0.1, DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final$MetaG_0.6, DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final$MetaG_0.1, DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final$MetaT_0.0, DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final$MetaG_0.0, decreasing = TRUE),]

write.csv(DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final_Sorted, "/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SpeciesAboveThreshold/ListSpecAboveThresh_AllConf_MetaT_MetaG.csv", row.names = FALSE)

# Export the same table, but without the lines containing the 0.0 thresholds

DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106 <- DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final[, c("ListSpecies", "MetaG_0.1", "MetaG_0.6", "MetaT_0.1", "MetaT_0.6")]

DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106_Sorted <- DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106[order(DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106$MetaT_0.6, DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106$MetaT_0.1, DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106$MetaG_0.6, DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106$MetaG_0.1, decreasing = TRUE),]

NonEmptyRowIndices <- which(rowSums(DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106_Sorted[,!(colnames(DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106_Sorted) %in% c("ListSpecies"))]) != 0)

DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106_Sorted_Nonempty <- DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106_Sorted[NonEmptyRowIndices,]

write.csv(DF_ListSpecAboveThresh_AllConf_MetaT_MetaG_Final0106_Sorted_Nonempty, "/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SpeciesAboveThreshold/ListSpecAboveThresh_AllConf_MetaT_MetaG_0.1_0.6_only.csv", row.names = FALSE)

# MetaG
#------

# Conf 0.1

DF_PavianSpecies_Conf_0.1_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.1_Norm
DF_PavianSpecies_Conf_0.1_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.1_Norm_AboveThresh[, c("AGTU-34", "AGTU-263", "AGTU-67")]

ListSpecAboveThresh_0.1 <- list()
for (i in (1:(ncol(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh)))) {
  TheColName <- colnames(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh)[i]
  SpecCounts <- DF_PavianSpecies_Conf_0.1_Norm_AboveThresh[[TheColName]][which(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  SpecNames <- row.names(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh)[which(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  names(SpecCounts) <- SpecNames
  ListSpecAboveThresh_0.1[[TheColName]] <- SpecCounts
} # Check thoroughly

ListSpecAboveThresh_0.1 # Only Epichloe typhina in AGTU-34

# Conf 0.6

DF_PavianSpecies_Conf_0.6_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.6_Norm
DF_PavianSpecies_Conf_0.6_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.6_Norm_AboveThresh[, c("AGTU-34", "AGTU-263", "AGTU-67")]

ListSpecAboveThresh_0.6 <- list()
for (i in (1:(ncol(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh)))) {
  TheColName <- colnames(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh)[i]
  SpecCounts <- DF_PavianSpecies_Conf_0.6_Norm_AboveThresh[[TheColName]][which(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  SpecNames <- row.names(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh)[which(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  names(SpecCounts) <- SpecNames
  ListSpecAboveThresh_0.6[[TheColName]] <- SpecCounts
} # Check thoroughly

ListSpecAboveThresh_0.6 # Nothing

