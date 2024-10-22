# 2024/06/03: I  have removed five samples as they were resequenced by Fasteris as they were considered of a too low quality.
# Check this does not disturb the code (especially when using index numbers in dataframes).
# 2024/10/01: I corrected the Kraken2 reports for 4 samples (AGTU-3, AGTU-6, AGTU-9 and AGTU-25). I replaced the corresponding Kraken2 reports and Pavian outputs  
# directly in the source folder.

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
library(pheatmap)
library(paletteer)
library(pairwiseAdonis) # Added on 2024-09-25
library(indicspecies) # Added on 2024-09-25

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

DF_PavianSpecies_Conf_0.0 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaG/SortiesKraken2/Pavian_Confidence_0.0/Analysis/RawPavianReports/SaraSampling_MetaG_Kraken2_Confidence_0.0_Fungi_Species.tsv", sep = "\t", header = TRUE)

colnames(DF_PavianSpecies_Conf_0.0)[5:114] <- gsub("\\.", "-", substr(colnames(DF_PavianSpecies_Conf_0.0)[5:114], 23, nchar(colnames(DF_PavianSpecies_Conf_0.0)[5:114])-11))

# Conf 0.1

DF_PavianSpecies_Conf_0.1 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaG/SortiesKraken2/Pavian_Confidence_0.1/Analysis/RawPavianReports/SaraSampling_MetaG_Kraken2_Confidence_0.1_Fungi_Species.tsv", sep = "\t", header = TRUE)

colnames(DF_PavianSpecies_Conf_0.1)[5:114] <- gsub("\\.", "-", substr(colnames(DF_PavianSpecies_Conf_0.1)[5:114], 23, nchar(colnames(DF_PavianSpecies_Conf_0.1)[5:114])-11))

# Conf 0.6

DF_PavianSpecies_Conf_0.6 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaG/SortiesKraken2/Pavian_Confidence_0.6/Analysis/RawPavianReports/SaraSampling_MetaG_Kraken2_Confidence_0.6_Fungi_Species.tsv", sep = "\t", header = TRUE)

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

write.csv(SortedSpecies_Conf_0.0, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/MetaG/Species/RelativeAbundance/RelativeAbundanceSpecies_Conf_0.0.csv", row.names = FALSE)

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

write.csv(SortedSpecies_Conf_0.1, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/MetaG/Species/RelativeAbundance/RelativeAbundanceSpecies_Conf_0.1.csv", row.names = FALSE)

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

write.csv(SortedSpecies_Conf_0.6, file = "/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/Tables/MetaG/Species/RelativeAbundance/RelativeAbundanceSpecies_Conf_0.6.csv", row.names = FALSE)

#==========================================================#
#                Make alpha diversity plot                 # # Check I did not make mistakes when copy-pasting the code for 0.1 and 0.6
#==========================================================#

# Conf 0.0

DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha <- DF_PavianSpecies_Conf_0.0_Norm
DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha <- DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha %>% replace(is.na(.), 0)

DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol <- DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha[,!(colnames(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]

ShannonInd_Conf_0.0 <- sapply(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol, function(x) vegan::diversity(x, index = "shannon"))
ShannonInd_Conf_0.0_v <- unlist(ShannonInd_Conf_0.0)
ShannonConf0.0_DF <- data.frame("Shannon"=ShannonInd_Conf_0.0_v)
ShannonConf0.0_DF[["SiteID"]] <- substr(row.names(ShannonConf0.0_DF), 3, 5)
ShannonConf0.0_DF$SiteID <- factor(ShannonConf0.0_DF$SiteID, levels = c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF", "FCM", "FLM", "FLA", "CDI", "CPL", "CHI", "CBD", "CDE", "CHD", "CYH", "CCO", "CXH", "CBS", "CMH", "CHS", "CRA"))

ShannonConf0.0_DF$Country <- substr(ShannonConf0.0_DF$SiteID, 1, 1)
ShannonConf0.0_DF_F <- ShannonConf0.0_DF[ShannonConf0.0_DF[["Country"]] == "F",]

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

ShannonConf0.0_DF_F <- ShannonConf0.0_DF_F[!(row.names(ShannonConf0.0_DF_F) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # Added G_FSS01.1_G1 on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

ShannonConf0.0_DF_F <- ShannonConf0.0_DF_F[ShannonConf0.0_DF_F[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

# End of added on 2024-06-03

p1 <- ggplot(data = ShannonConf0.0_DF_F, aes(x=SiteID, y=Shannon)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Shannon index (H')") +
  #ylim(c(0, 3.5)) +
  #scale_fill_manual(values = c("chocolate1", "darkolivegreen3", "gold", "deepskyblue")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position="none"
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Shannon_Conf_0.0.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Shannon_Conf_0.0.png", plot=p1, device =  "png", dpi = 300, width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Shannon_Conf_0.0.jpeg", plot=p1, device =  "jpeg", dpi = 300, width=21, height=20)
dev.off()

SimpsonInd_Conf_0.0 <- sapply(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol, function(x) vegan::diversity(x, index = "simpson"))
SimpsonInd_Conf_0.0_v <- unlist(SimpsonInd_Conf_0.0)
SimpsonConf0.0_DF <- data.frame("Simpson"=SimpsonInd_Conf_0.0_v)
SimpsonConf0.0_DF[["SiteID"]] <- substr(row.names(SimpsonConf0.0_DF), 3, 5)
SimpsonConf0.0_DF$SiteID <- factor(SimpsonConf0.0_DF$SiteID, levels = c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF", "FCM", "FLM", "FLA", "CDI", "CPL", "CHI", "CBD", "CDE", "CHD", "CYH", "CCO", "CXH", "CBS", "CMH", "CHS", "CRA"))

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

SimpsonConf0.0_DF <- SimpsonConf0.0_DF[!(row.names(SimpsonConf0.0_DF) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 addedon 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

SimpsonConf0.0_DF <- SimpsonConf0.0_DF[SimpsonConf0.0_DF[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

# End of added on 2024-06-03

p1 <- ggplot(data = SimpsonConf0.0_DF, aes(x=SiteID, y=Simpson)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Gini-Simpson index (1-D)") +
  #ylim(c(0, 3.5)) +
  #scale_fill_manual(values = c("chocolate1", "darkolivegreen3", "gold", "deepskyblue")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position="none"
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Simpson_Conf_0.0.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Simpson_Conf_0.0.png", plot=p1, device =  "png", dpi = 300, width=21, height=20)
dev.off()

#=========================
# Statistical tests

# Shannon 
kruskal_test(data = ShannonConf0.0_DF_F, formula = Shannon ~ SiteID) # Significant
print(dunn_test(data = ShannonConf0.0_DF_F, formula = Shannon ~ SiteID, p.adjust.method = "bonferroni"), n = 30) # Significant bewteen FGR and FCT

# Simpson 
kruskal_test(data = SimpsonConf0.0_DF, formula = Simpson ~ SiteID) # Significant
print(dunn_test(data = SimpsonConf0.0_DF, formula = Simpson ~ SiteID, p.adjust.method = "bonferroni"), n = 30)  # Significant between FMS and FCT, and FSS and FCT

#=========================

# Conf 0.1

DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha <- DF_PavianSpecies_Conf_0.1_Norm
DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha %>% replace(is.na(.), 0)

DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha[,!(colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]

ShannonInd_Conf_0.1 <- sapply(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol, function(x) vegan::diversity(x, index = "shannon"))
ShannonInd_Conf_0.1_v <- unlist(ShannonInd_Conf_0.1)
ShannonConf0.1_DF <- data.frame("Shannon"=ShannonInd_Conf_0.1_v)
ShannonConf0.1_DF[["SiteID"]] <- substr(row.names(ShannonConf0.1_DF), 3, 5)
ShannonConf0.1_DF$SiteID <- factor(ShannonConf0.1_DF$SiteID, levels = c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF", "FCM", "FLM", "FLA", "CDI", "CPL", "CHI", "CBD", "CDE", "CHD", "CYH", "CCO", "CXH", "CBS", "CMH", "CHS", "CRA"))

ShannonConf0.1_DF$Country <- substr(ShannonConf0.1_DF$SiteID, 1, 1)
ShannonConf0.1_DF_F <- ShannonConf0.1_DF[ShannonConf0.1_DF[["Country"]] == "F",]

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

ShannonConf0.1_DF_F <- ShannonConf0.1_DF_F[!(row.names(ShannonConf0.1_DF_F) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

ShannonConf0.1_DF_F <- ShannonConf0.1_DF_F[ShannonConf0.1_DF_F[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

# End of added on 2024-06-03

p1 <- ggplot(data = ShannonConf0.1_DF_F, aes(x=SiteID, y=Shannon)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1, outlier.size = 6) +
  ylab(label = "Shannon index (H')") +
  xlab(label = "Site") +
  #ylim(c(0, 3.5)) +
  #scale_fill_manual(values = c("chocolate1", "darkolivegreen3", "gold", "deepskyblue")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 40),
    axis.title = element_text(size = 40),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.position="none"
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Shannon_Conf_0.1.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Shannon_Conf_0.1.png", plot=p1, device =  "png", dpi = 300, width=21, height=20)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/FigAlphaDivFranceSpecies_Shannon_Conf_0.1.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/FigAlphaDivFranceSpecies_Shannon_Conf_0.1.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

SimpsonInd_Conf_0.1 <- sapply(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol, function(x) vegan::diversity(x, index = "simpson"))
SimpsonInd_Conf_0.1_v <- unlist(SimpsonInd_Conf_0.1)
SimpsonConf0.1_DF <- data.frame("Simpson"=SimpsonInd_Conf_0.1_v)
SimpsonConf0.1_DF[["SiteID"]] <- substr(row.names(SimpsonConf0.1_DF), 3, 5)
SimpsonConf0.1_DF$SiteID <- factor(SimpsonConf0.1_DF$SiteID, levels = c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF", "FCM", "FLM", "FLA", "CDI", "CPL", "CHI", "CBD", "CDE", "CHD", "CYH", "CCO", "CXH", "CBS", "CMH", "CHS", "CRA"))

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

SimpsonConf0.1_DF <- SimpsonConf0.1_DF[!(row.names(SimpsonConf0.1_DF) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

SimpsonConf0.1_DF <- SimpsonConf0.1_DF[SimpsonConf0.1_DF[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

# End of added on 2024-06-03

p1 <- ggplot(data = SimpsonConf0.1_DF, aes(x=SiteID, y=Simpson)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1, outlier.size = 6) +
  ylab(label = "Gini-Simpson index (1-D)") +
  xlab(label = "Site") +
  #ylim(c(0, 3.5)) +
  #scale_fill_manual(values = c("chocolate1", "darkolivegreen3", "gold", "deepskyblue")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 40),
    axis.title = element_text(size = 40),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.position="none"
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Simpson_Conf_0.1.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Simpson_Conf_0.1.png", plot=p1, device = "png", dpi = 300, width=21, height=20)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/FigAlphaDivFranceSpecies_Simpson_Conf_0.1.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/FigAlphaDivFranceSpecies_Simpson_Conf_0.1.png", plot=p1, device = "png", dpi = 300, width=21, height=20)

dev.off()

#=========================
# Statistical tests

# Shannon 
kruskal_test(data = ShannonConf0.1_DF_F, formula = Shannon ~ SiteID) # Not significant
print(dunn_test(data = ShannonConf0.1_DF_F, formula = Shannon ~ SiteID, p.adjust.method = "bonferroni"), n = 30) # Nothing significant with adjusted p-value

# Simpson 
kruskal_test(data = SimpsonConf0.1_DF, formula = Simpson ~ SiteID) # Significant
print(dunn_test(data = SimpsonConf0.1_DF, formula = Simpson ~ SiteID, p.adjust.method = "bonferroni"), n = 30) # Significant between FSS and FHB

#=========================

# Conf 0.6

DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha <- DF_PavianSpecies_Conf_0.6_Norm
DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha <- DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha %>% replace(is.na(.), 0)

DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol <- DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha[,!(colnames(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]

ShannonInd_Conf_0.6 <- sapply(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol, function(x) vegan::diversity(x, index = "shannon"))
ShannonInd_Conf_0.6_v <- unlist(ShannonInd_Conf_0.6)
ShannonConf0.6_DF <- data.frame("Shannon"=ShannonInd_Conf_0.6_v)
ShannonConf0.6_DF[["SiteID"]] <- substr(row.names(ShannonConf0.6_DF), 3, 5)
ShannonConf0.6_DF$SiteID <- factor(ShannonConf0.6_DF$SiteID, levels = c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF", "FCM", "FLM", "FLA", "CDI", "CPL", "CHI", "CBD", "CDE", "CHD", "CYH", "CCO", "CXH", "CBS", "CMH", "CHS", "CRA"))

ShannonConf0.6_DF$Country <- substr(ShannonConf0.6_DF$SiteID, 1, 1)
ShannonConf0.6_DF_F <- ShannonConf0.6_DF[ShannonConf0.6_DF[["Country"]] == "F",]

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

ShannonConf0.6_DF_F <- ShannonConf0.6_DF_F[!(row.names(ShannonConf0.6_DF_F) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

ShannonConf0.6_DF_F <- ShannonConf0.6_DF_F[ShannonConf0.6_DF_F[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

# End of added on 2024-06-03

p1 <- ggplot(data = ShannonConf0.6_DF_F, aes(x=SiteID, y=Shannon)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Shannon index (H')") +
  #ylim(c(0, 3.5)) +
  #scale_fill_manual(values = c("chocolate1", "darkolivegreen3", "gold", "deepskyblue")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position="none"
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Shannon_Conf_0.6.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Shannon_Conf_0.6.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

SimpsonInd_Conf_0.6 <- sapply(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol, function(x) vegan::diversity(x, index = "simpson"))
SimpsonInd_Conf_0.6_v <- unlist(SimpsonInd_Conf_0.6)
SimpsonConf0.6_DF <- data.frame("Simpson"=SimpsonInd_Conf_0.6_v)
SimpsonConf0.6_DF[["SiteID"]] <- substr(row.names(SimpsonConf0.6_DF), 3, 5)
SimpsonConf0.6_DF$SiteID <- factor(SimpsonConf0.6_DF$SiteID, levels = c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF", "FCM", "FLM", "FLA", "CDI", "CPL", "CHI", "CBD", "CDE", "CHD", "CYH", "CCO", "CXH", "CBS", "CMH", "CHS", "CRA"))

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

SimpsonConf0.6_DF <- SimpsonConf0.6_DF[!(row.names(SimpsonConf0.6_DF) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

SimpsonConf0.6_DF <- SimpsonConf0.6_DF[SimpsonConf0.6_DF[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

# End of added on 2024-06-03

p1 <- ggplot(data = SimpsonConf0.6_DF, aes(x=SiteID, y=Simpson)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Gini-Simpson index (1-D)") +
  #ylim(c(0, 3.5)) +
  #scale_fill_manual(values = c("chocolate1", "darkolivegreen3", "gold", "deepskyblue")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position="none"
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Simpson_Conf_0.6.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/AlphaDiversity/Species/FigAlphaDivFranceSpecies_Simpson_Conf_0.6.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

#=========================
# Statistical tests

# Shannon 
kruskal_test(data = ShannonConf0.6_DF_F, formula = Shannon ~ SiteID) # Significant
print(dunn_test(data = ShannonConf0.6_DF_F, formula = Shannon ~ SiteID, p.adjust.method = "bonferroni"), n = 30) # Significant between FHB and FLC, and FHB and FLP

# Simpson 
kruskal_test(data = SimpsonConf0.6_DF, formula = Simpson ~ SiteID) # Not significant
print(dunn_test(data = SimpsonConf0.6_DF, formula = Simpson ~ SiteID, p.adjust.method = "bonferroni"), n = 30) # Not significant

#=========================

#================#
# Beta-diversity #
#================#

CustomPalette <- c("black", "#DF536B", "#61D04F", "#2297E6", "#28E2E5", "#CD0BBC", "#F5C710", "gray62", "#8DD3C7",  "#FB8072", "#BEBADA")

# Conf 0.0
#---------

# Take only France samples

DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t <- as.data.frame(t(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol))
DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t$Country <- substr(row.names(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t), 3, 3)
DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F <- DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t[DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t[["Country"]] == "F",]
DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F$Country <- NULL

# Added on 2024-05-31 to remove the site we decided not to keep

DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F$Site <- substr(row.names(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F), 3, 5)
DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F <- DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F[DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F[["Site"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]
DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F$Site <- NULL

# End of added on 2024-05-31

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F <- DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F[!(row.names(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08

# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

# End of added on 2024-06-03

NMDS_Species_0.0 <- metaMDS(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_GoodCol_t_F, distance = "bray", k = 2, autotransform = TRUE, trymax=100)
data.scores_0.0 <- as.data.frame(scores(NMDS_Species_0.0)[["sites"]])
data.scores_0.0$SiteID <- substr(row.names(data.scores_0.0), 3, 5)

data.scores_0.0$SiteID <- factor(data.scores_0.0$SiteID, levels = c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF", "FCM", "FLM", "FLA"))

hull_data_0.0 <- 
  data.scores_0.0 %>%
  drop_na() %>%
  group_by(SiteID) %>% 
  slice(chull(NMDS1, NMDS2))

ggplot(data=data.scores_0.0, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point() +
  scale_colour_manual(values = CustomPalette) +
  geom_polygon(data = hull_data_0.0,
               aes(colour = SiteID),
               alpha = 0.0,
               show.legend = FALSE) +
  #stat_ellipse() +
  theme_bw()


ggplot(data=data.scores_0.0, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point() +
  scale_colour_manual(values = CustomPalette) +
  stat_ellipse() +
  theme_bw()

# Tester si les points se regroupent par moment dans l'annee

data.scores_0.0$SampleName <- row.names(data.scores_0.0)
data.scores_0.0_metadata <- left_join(data.scores_0.0, DF_MergedMetadata) # /!\ Check 
data.scores_0.0_metadata$Date_num <- as.Date(data.scores_0.0_metadata$date, "%Y-%m-%d")
data.scores_0.0_metadata$JulianDay <- yday(data.scores_0.0_metadata$Date_num) # Check

PlotTest <- ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = JulianDay)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(10)) +
  theme_bw()

# ggplot() +
#   geom_text(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = JulianDay, label=SiteID)) +
#   scale_color_gradientn(colours = rainbow(10)) +
#   geom_polygon(data = hull_data_0.0,
#                aes(x = NMDS1, y = NMDS2), colour = "black",
#                #aes(colour = SiteID),
#                alpha = 0.0,
#                show.legend = FALSE) +
#   theme_bw()

plot(PlotTest)

# PlotTest <- ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, fill = JulianDay, color = SiteID)) +
#   geom_point(shape = 21, stroke = 0.01) +
#   scale_fill_gradientn(colours = rainbow(10)) +
#   scale_colour_manual(values = CustomPalette) +
#   theme_bw()
# 
# ggplotly(PlotTest)

# Check grouping by salinity 

ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = salinity)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = `T`)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = BGE)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = BP)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = RES)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = chla.mean)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = DOC.mean)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = protein.like.DOM)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = humic.like.DOM)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

# Idem, but shape as a function of the season

data.scores_0.0_metadata$Month <- month(data.scores_0.0_metadata$Date_num)

# Make the seasons column

data.scores_0.0_metadata$Season <- NA
for (i in 1:dim(data.scores_0.0_metadata)[1]) {
  if (data.scores_0.0_metadata[i,"Month"] %in% c(1, 2, 3)) {
    data.scores_0.0_metadata[i, "Season"] <- "Winter"
  } else if (data.scores_0.0_metadata[i,"Month"] %in% c(4, 5, 6)) {
    data.scores_0.0_metadata[i, "Season"] <- "Spring"
  } else if (data.scores_0.0_metadata[i,"Month"] %in% c(7, 8, 9)) {
    data.scores_0.0_metadata[i, "Season"] <- "Summer"
  } else if (data.scores_0.0_metadata[i,"Month"] %in% c(10, 11, 12)) {
    data.scores_0.0_metadata[i, "Season"] <- "Autumn"
  } else {}
}

ggplot(data=data.scores_0.0_metadata, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point(aes(shape = Season)) +
  scale_colour_manual(values = CustomPalette) +
  geom_polygon(data = hull_data_0.0,
               aes(colour = SiteID),
               alpha = 0.0,
               show.legend = FALSE) +
  #stat_ellipse() +
  theme_bw()

#
# A completer en faisant n NMDS, une pour chaque site
#

# Conf 0.1
#---------

# Take only France samples

DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t <- as.data.frame(t(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol))
DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t$Country <- substr(row.names(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t), 3, 3)
DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t[DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t[["Country"]] == "F",]
DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F$Country <- NULL

# Added on 2024-05-31 to remove the site we decided not to keep

DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F$Site <- substr(row.names(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F), 3, 5)
DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F[DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F[["Site"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]
DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F$Site <- NULL

# End of added on 2024-05-31

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F[!(row.names(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08

# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

# End of added on 2024-06-03

NMDS_Species_0.1 <- metaMDS(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F, distance = "bray", k = 2, autotransform = TRUE, trymax=100)
data.scores_0.1 <- as.data.frame(scores(NMDS_Species_0.1)[["sites"]])
data.scores_0.1$SiteID <- substr(row.names(data.scores_0.1), 3, 5)

data.scores_0.1$SiteID <- factor(data.scores_0.1$SiteID, levels = c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF", "FCM", "FLM", "FLA"))

hull_data_0.1 <- 
  data.scores_0.1 %>%
  drop_na() %>%
  group_by(SiteID) %>% 
  slice(chull(NMDS1, NMDS2))

ggplot(data=data.scores_0.1, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point() +
  scale_colour_manual(values = CustomPalette) +
  geom_polygon(data = hull_data_0.1,
               aes(colour = SiteID),
               alpha = 0.0,
               show.legend = FALSE) +
  #stat_ellipse() +
  theme_bw()


ggplot(data=data.scores_0.1, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point() +
  scale_colour_manual(values = CustomPalette) +
  stat_ellipse() +
  theme_bw()

# Tester si les points se regroupent par moment dans l'annee

data.scores_0.1$SampleName <- row.names(data.scores_0.1)
data.scores_0.1_metadata <- left_join(data.scores_0.1, DF_MergedMetadata) # /!\ Check 
data.scores_0.1_metadata$Date_num <- as.Date(data.scores_0.1_metadata$date, "%Y-%m-%d")
data.scores_0.1_metadata$JulianDay <- yday(data.scores_0.1_metadata$Date_num) # Check

PlotTest <- ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = JulianDay)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(10)) +
  theme_bw()

# ggplot() +
#   geom_text(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = JulianDay, label=SiteID)) +
#   scale_color_gradientn(colours = rainbow(10)) +
#   geom_polygon(data = hull_data_0.1,
#                aes(x = NMDS1, y = NMDS2), colour = "black",
#                #aes(colour = SiteID),
#                alpha = 0.0,
#                show.legend = FALSE) +
#   theme_bw()

plot(PlotTest)

# PlotTest <- ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, fill = JulianDay, color = SiteID)) +
#   geom_point(shape = 21, stroke = 0.01) +
#   scale_fill_gradientn(colours = rainbow(10)) +
#   scale_colour_manual(values = CustomPalette) +
#   theme_bw()
# 
# ggplotly(PlotTest)

# Check grouping by salinity 

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = salinity)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = `T`)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = BGE)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = BP)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = RES)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = chla.mean)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = DOC.mean)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = protein.like.DOM)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = humic.like.DOM)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

# Idem, but shape as a function of the season

data.scores_0.1_metadata$Month <- month(data.scores_0.1_metadata$Date_num)

# Make the seasons column

data.scores_0.1_metadata$Season <- NA
for (i in 1:dim(data.scores_0.1_metadata)[1]) {
  if (data.scores_0.1_metadata[i,"Month"] %in% c(1, 2, 3)) {
    data.scores_0.1_metadata[i, "Season"] <- "Winter"
  } else if (data.scores_0.1_metadata[i,"Month"] %in% c(4, 5, 6)) {
    data.scores_0.1_metadata[i, "Season"] <- "Spring"
  } else if (data.scores_0.1_metadata[i,"Month"] %in% c(7, 8, 9)) {
    data.scores_0.1_metadata[i, "Season"] <- "Summer"
  } else if (data.scores_0.1_metadata[i,"Month"] %in% c(10, 11, 12)) {
    data.scores_0.1_metadata[i, "Season"] <- "Autumn"
  } else {}
}

p1 <- ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point(aes(shape = Season)) +
  #scale_color_manual(values = c("deepskyblue", "gold", "gold2", "darkolivegreen1", "darkolivegreen3", "darkolivegreen4", "forestgreen")) + # Added on 2024-06-28
  scale_color_manual(values = c("darkorchid2", "gold", "goldenrod4", "deepskyblue", "dodgerblue4", "darkolivegreen1", "forestgreen")) +
  #scale_colour_manual(values = CustomPalette) +
  geom_polygon(data = hull_data_0.1,
               aes(colour = SiteID),
               alpha = 0.0,
               show.legend = FALSE) +
  #stat_ellipse() +
  theme_bw()

plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/NMDS_ColorSite_0.1.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/NMDS_ColorSite_0.1.png", plot=p1, device = "png", dpi = 300, width=10, height=10)
dev.off()

#
# A completer en faisant n NMDS, une pour chaque site
#

#===========================================#
#            Added on 2024-06-28            #
#===========================================#

# Make the same NMDS, but normalize Fungi to the percent of classified Fungal species

DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForNormSpec <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F
DF_0.1_ForNormSpecNew_20240628 <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForNormSpec/rowSums(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForNormSpec) # To be checked

NMDS_Species_0.1 <- metaMDS(DF_0.1_ForNormSpecNew_20240628, distance = "bray", k = 2, autotransform = TRUE, trymax=100)
data.scores_0.1 <- as.data.frame(scores(NMDS_Species_0.1)[["sites"]])
data.scores_0.1$SiteID <- substr(row.names(data.scores_0.1), 3, 5)

data.scores_0.1$SiteID <- factor(data.scores_0.1$SiteID, levels = c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF", "FCM", "FLM", "FLA"))

hull_data_0.1 <- 
  data.scores_0.1 %>%
  drop_na() %>%
  group_by(SiteID) %>% 
  slice(chull(NMDS1, NMDS2))

ggplot(data=data.scores_0.1, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point() +
  scale_colour_manual(values = CustomPalette) +
  geom_polygon(data = hull_data_0.1,
               aes(colour = SiteID),
               alpha = 0.0,
               show.legend = FALSE) +
  #stat_ellipse() +
  theme_bw()


ggplot(data=data.scores_0.1, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point() +
  scale_colour_manual(values = CustomPalette) +
  stat_ellipse() +
  theme_bw()

# Tester si les points se regroupent par moment dans l'annee

data.scores_0.1$SampleName <- row.names(data.scores_0.1)
data.scores_0.1_metadata <- left_join(data.scores_0.1, DF_MergedMetadata) # /!\ Check 
data.scores_0.1_metadata$Date_num <- as.Date(data.scores_0.1_metadata$date, "%Y-%m-%d")
data.scores_0.1_metadata$JulianDay <- yday(data.scores_0.1_metadata$Date_num) # Check

PlotTest <- ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = JulianDay)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(10)) +
  theme_bw()

# ggplot() +
#   geom_text(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = JulianDay, label=SiteID)) +
#   scale_color_gradientn(colours = rainbow(10)) +
#   geom_polygon(data = hull_data_0.1,
#                aes(x = NMDS1, y = NMDS2), colour = "black",
#                #aes(colour = SiteID),
#                alpha = 0.0,
#                show.legend = FALSE) +
#   theme_bw()

plot(PlotTest)

# PlotTest <- ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, fill = JulianDay, color = SiteID)) +
#   geom_point(shape = 21, stroke = 0.01) +
#   scale_fill_gradientn(colours = rainbow(10)) +
#   scale_colour_manual(values = CustomPalette) +
#   theme_bw()
# 
# ggplotly(PlotTest)

# Check grouping by salinity 

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = salinity)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = `T`)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = BGE)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = BP)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = RES)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = chla.mean)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = DOC.mean)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = protein.like.DOM)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = humic.like.DOM)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

# Idem, but shape as a function of the season

data.scores_0.1_metadata$Month <- month(data.scores_0.1_metadata$Date_num)

# Make the seasons column

data.scores_0.1_metadata$Season <- NA
for (i in 1:dim(data.scores_0.1_metadata)[1]) {
  if (data.scores_0.1_metadata[i,"Month"] %in% c(1, 2, 3)) {
    data.scores_0.1_metadata[i, "Season"] <- "Winter"
  } else if (data.scores_0.1_metadata[i,"Month"] %in% c(4, 5, 6)) {
    data.scores_0.1_metadata[i, "Season"] <- "Spring"
  } else if (data.scores_0.1_metadata[i,"Month"] %in% c(7, 8, 9)) {
    data.scores_0.1_metadata[i, "Season"] <- "Summer"
  } else if (data.scores_0.1_metadata[i,"Month"] %in% c(10, 11, 12)) {
    data.scores_0.1_metadata[i, "Season"] <- "Autumn"
  } else {}
}

p1 <- ggplot(data=data.scores_0.1_metadata, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point(aes(shape = Season)) +
  #scale_color_manual(values = c("deepskyblue", "gold", "gold2", "darkolivegreen1", "darkolivegreen3", "darkolivegreen4", "forestgreen")) + # Added on 2024-06-28
  scale_color_manual(values = c("darkorchid2", "gold", "goldenrod4", "deepskyblue", "dodgerblue4", "darkolivegreen1", "forestgreen")) +
  #scale_colour_manual(values = CustomPalette) +
  geom_polygon(data = hull_data_0.1,
               aes(colour = SiteID),
               alpha = 0.0,
               show.legend = FALSE) +
  #stat_ellipse() +
  theme_bw()

plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/NMDS_ColorSite_0.1_NormAssignedFungalSpecies.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/NMDS_ColorSite_0.1_NormAssignedFungalSpecies.png", plot=p1, device = "png", dpi = 300, width=10, height=10)
dev.off()


#===========================================#
#         End of added on 2024-06-28        #
#===========================================#

# Conf 0.6
#---------

# Take only France samples

DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t <- as.data.frame(t(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol))
DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t$Country <- substr(row.names(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t), 3, 3)
DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F <- DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t[DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t[["Country"]] == "F",]
DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F$Country <- NULL

# Added on 2024-05-31 to remove the site we decided not to keep

DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F$Site <- substr(row.names(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F), 3, 5)
DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F <- DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F[DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F[["Site"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]
DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F$Site <- NULL

# End of added on 2024-05-31

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F <- DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F[!(row.names(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08

# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

# End of added on 2024-06-03

NMDS_Species_0.6 <- metaMDS(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_GoodCol_t_F, distance = "bray", k = 2, autotransform = TRUE, trymax=100)
data.scores_0.6 <- as.data.frame(scores(NMDS_Species_0.6)[["sites"]])
data.scores_0.6$SiteID <- substr(row.names(data.scores_0.6), 3, 5)

data.scores_0.6$SiteID <- factor(data.scores_0.6$SiteID, levels = c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF", "FCM", "FLM", "FLA"))

hull_data_0.6 <- 
  data.scores_0.6 %>%
  drop_na() %>%
  group_by(SiteID) %>% 
  slice(chull(NMDS1, NMDS2))

ggplot(data=data.scores_0.6, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point() +
  scale_colour_manual(values = CustomPalette) +
  geom_polygon(data = hull_data_0.6,
               aes(colour = SiteID),
               alpha = 0.0,
               show.legend = FALSE) +
  #stat_ellipse() +
  theme_bw()


ggplot(data=data.scores_0.6, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point() +
  scale_colour_manual(values = CustomPalette) +
  stat_ellipse() +
  theme_bw()

# Tester si les points se regroupent par moment dans l'annee

data.scores_0.6$SampleName <- row.names(data.scores_0.6)
data.scores_0.6_metadata <- left_join(data.scores_0.6, DF_MergedMetadata) # /!\ Check 
data.scores_0.6_metadata$Date_num <- as.Date(data.scores_0.6_metadata$date, "%Y-%m-%d")
data.scores_0.6_metadata$JulianDay <- yday(data.scores_0.6_metadata$Date_num) # Check

PlotTest <- ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = JulianDay)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(10)) +
  theme_bw()

# ggplot() +
#   geom_text(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = JulianDay, label=SiteID)) +
#   scale_color_gradientn(colours = rainbow(10)) +
#   geom_polygon(data = hull_data_0.6,
#                aes(x = NMDS1, y = NMDS2), colour = "black",
#                #aes(colour = SiteID),
#                alpha = 0.0,
#                show.legend = FALSE) +
#   theme_bw()

plot(PlotTest)

# PlotTest <- ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, fill = JulianDay, color = SiteID)) +
#   geom_point(shape = 21, stroke = 0.01) +
#   scale_fill_gradientn(colours = rainbow(10)) +
#   scale_colour_manual(values = CustomPalette) +
#   theme_bw()
# 
# ggplotly(PlotTest)

# Check grouping by salinity 

ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = salinity)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = `T`)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = BGE)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = BP)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = RES)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = chla.mean)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = DOC.mean)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = protein.like.DOM)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = humic.like.DOM)) +
  geom_point() +
  #scale_color_gradientn() +
  theme_bw()

# Idem, but shape as a function of the season

data.scores_0.6_metadata$Month <- month(data.scores_0.6_metadata$Date_num)

# Make the seasons column

data.scores_0.6_metadata$Season <- NA
for (i in 1:dim(data.scores_0.6_metadata)[1]) {
  if (data.scores_0.6_metadata[i,"Month"] %in% c(1, 2, 3)) {
    data.scores_0.6_metadata[i, "Season"] <- "Winter"
  } else if (data.scores_0.6_metadata[i,"Month"] %in% c(4, 5, 6)) {
    data.scores_0.6_metadata[i, "Season"] <- "Spring"
  } else if (data.scores_0.6_metadata[i,"Month"] %in% c(7, 8, 9)) {
    data.scores_0.6_metadata[i, "Season"] <- "Summer"
  } else if (data.scores_0.6_metadata[i,"Month"] %in% c(10, 11, 12)) {
    data.scores_0.6_metadata[i, "Season"] <- "Autumn"
  } else {}
}

ggplot(data=data.scores_0.6_metadata, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
  geom_point(aes(shape = Season)) +
  scale_colour_manual(values = CustomPalette) +
  geom_polygon(data = hull_data_0.6,
               aes(colour = SiteID),
               alpha = 0.0,
               show.legend = FALSE) +
  #stat_ellipse() +
  theme_bw()

#
# A completer en faisant n NMDS, une pour chaque site
#

#=================================================================#
# Line plots with the 10 most abundant species overtime per site #
#=================================================================#

# Number of species to plot

Nspecies <- 10

# Conf 0.0
#---------

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_SamRm <- DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha[,!(colnames(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1"))] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

# End of added on 2024-06-03

# Make an empty list to store the plots

ListPlots <- list()

# Get the list of all sites

ListSites <- unique(DF_MergedMetadata$siteID)[which(is.na(unique(DF_MergedMetadata$siteID)) == FALSE)]

# Remove FLA and FLM for now as they do not have metagenomic data

ListSites <- c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT") # ListSites[which(!(ListSites %in% c("FLA", "FLM")))]

# Loop on all sites

for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the 10 most abundant species names
  
  SubDF <- DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  NamesSpecies <- names(sort(rowSums(SubDF), decreasing = TRUE))[1:Nspecies]

  # Extract the right species in the SubDF, transpose the dataframe and merge with the appropriate metadata
  
  SubDF_Species <- SubDF[row.names(SubDF) %in% NamesSpecies, , drop = FALSE] # Check that drop does not disturb
  SubDF_Species_t <- data.frame(t(SubDF_Species))
  SubDF_Species_t$SampleName <- row.names(SubDF_Species_t)
  
  SubDF_Species_t_metadata <- left_join(SubDF_Species_t, DF_MergedMetadata)
  
  # Convert the date to type date
  
  SubDF_Species_t_metadata$Date_num <- as.Date(SubDF_Species_t_metadata$date, "%Y-%m-%d")
  
  # Keep only the columns for the plot, then convert to long format 
  
  SubDF_Species_t_metadata_ForPlot <- SubDF_Species_t_metadata[, colnames(SubDF_Species_t_metadata) %in% c("Date_num", "SampleName", gsub(" ", ".", NamesSpecies))]
  
  SubDF_Species_t_metadata_ForPlot_Long <- SubDF_Species_t_metadata_ForPlot %>%
    pivot_longer(
      cols = gsub(" ", ".", NamesSpecies),
      names_to = "Species",
      values_to = "RelAbundance",
    )
  
  # Make the plot
  
  print(ListSites[i])
  ListPlots[[ListSites[i]]] <- ggplot(data = SubDF_Species_t_metadata_ForPlot_Long, aes(x = Date_num, y = RelAbundance, color = Species)) +
    geom_point() +
    geom_line() +
    scale_colour_manual(values = CustomPalette) +
    scale_x_date(date_breaks="1 month", date_labels="%b\n%Y", limits = c(min = as.Date("2021-12-02", "%Y-%m-%d"), max = as.Date("2022-12-06", "%Y-%m-%d"))) +
    labs(x = "Date", y = "Relative abundance", title = ListSites[i]) +
    theme_bw()
  
}

p1 <- gridExtra::grid.arrange(grobs = ListPlots)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_Conf_0.0.pdf", plot=p1, device = cairo_pdf(), width=49, height=49)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_Conf_0.0.png", plot=p1, device = "png", dpi = 300, width=49, height=49)
dev.off()

#
# Add the plots for confidence 0.1 and 0.6
#



#====================================================#
# Barplots with the n most abundant species per site #
#====================================================#

# Number of species to plot

Nspecies <- 100

# Conf 0.0
#---------

# Make an empty list to store the plots

ListPlots <- list()

# Get the list of all sites

ListSites <- unique(DF_MergedMetadata$siteID)[which(is.na(unique(DF_MergedMetadata$siteID)) == FALSE)]

# Remove FLA and FLM for now as they do not have metagenomic data

ListSites <- ListSites[which(!(ListSites %in% c("FLA", "FLM", "FBB", "FCF", "FCM", "FMD")))] # Before: ListSites <- ListSites[which(!(ListSites %in% c("FLA", "FLM")))]

# Loop on all sites

for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  MostAbSpecies <- sort(rowSums(SubDF), decreasing = TRUE)[1:Nspecies]

  # Plot the most abundant species in the form of a barplot
  
  MostAbSpecies_DF <- data.frame("Species" = names(MostAbSpecies), "RelAb" = MostAbSpecies)
  MostAbSpecies_DF$Species <- factor(MostAbSpecies_DF$Species, levels = names(MostAbSpecies))

  ListPlots[[ListSites[i]]] <- ggplot(data = MostAbSpecies_DF, aes(x = Species, y = RelAb)) +
    geom_bar(stat="identity") +
    labs(x = "Species", y = "Summed relative abundance", title = ListSites[i]) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    theme_bw()

}
  
p1 <- gridExtra::grid.arrange(grobs = ListPlots)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_SummedAllSamples_Conf_0.0.pdf", plot=p1, device = cairo_pdf(), width=49, height=49)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_SummedAllSamples_Conf_0.0.png", plot=p1, device = "png", dpi = 300, width=49, height=49)
dev.off()
  
# Take the first 10 species per site

Nspecies <- 10

ListMostAbundSpecies <- list()
for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  MostAbSpecies <- sort(rowSums(SubDF), decreasing = TRUE)[1:Nspecies]
  ListMostAbundSpecies[[ListSites[i]]] <- MostAbSpecies
}

# List all the names of the species that appear within the 10 most abundant species in at least 1 site
ListNmostAbAllSites <- unique(unlist(lapply(ListMostAbundSpecies, names)))

# For each site extract the average relative abundance of these species

ListMostAbundSpecies <- list()
for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.0_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  #MostAbSpecies <- sort(rowSums(SubDF), decreasing = TRUE)[1:Nspecies]
  MostAbSpecies <- SubDF[ListNmostAbAllSites,]
  NumberSamples <- ncol(MostAbSpecies)
  MostAbSpeciesSums <- rowSums(MostAbSpecies)/NumberSamples
  ListMostAbundSpecies[[ListSites[i]]] <- MostAbSpeciesSums
}

# Convert the list to dataframe, then convert to long format

DF_RelAbMostAbSpecForPlot <- as.data.frame(ListMostAbundSpecies)
DF_RelAbMostAbSpecForPlot_t <- data.frame(t(DF_RelAbMostAbSpecForPlot))
NamesFactor <- names(sort(colSums(DF_RelAbMostAbSpecForPlot_t), decreasing = TRUE))
DF_RelAbMostAbSpecForPlot_t$Site <- row.names(DF_RelAbMostAbSpecForPlot_t)
DF_RelAbMostAbSpecForPlot_t_longer <- pivot_longer(DF_RelAbMostAbSpecForPlot_t, cols = !Site, names_to = "Species", values_to = "RelAbundance")
DF_RelAbMostAbSpecForPlot_t_longer$Species <- factor(DF_RelAbMostAbSpecForPlot_t_longer$Species, levels = NamesFactor)
DF_RelAbMostAbSpecForPlot_t_longer$Site <- factor(DF_RelAbMostAbSpecForPlot_t_longer$Site, levels = c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"))

# Convert relative abundance to % of all reads

DF_RelAbMostAbSpecForPlot_t_longer$RelAbundance <- DF_RelAbMostAbSpecForPlot_t_longer$RelAbundance * 100

# Make the plot

p1 <- ggplot(data = DF_RelAbMostAbSpecForPlot_t_longer, aes(x = Species, y = RelAbundance, fill = Site)) +
  geom_bar(position='dodge', stat='identity') +
  #scale_fill_manual(values = c("deepskyblue", "gold", "gold2", "darkolivegreen1", "darkolivegreen3", "darkolivegreen4", "forestgreen")) +
  scale_fill_manual(values = c("darkorchid2", "gold", "goldenrod4", "deepskyblue", "dodgerblue4", "darkolivegreen1", "forestgreen")) +
  labs(x = "Species", y = "Mean relative abundance (% of all reads)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))

plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/Fig10MostAbTaxaPerSite_Species_Conf_0.0.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/Fig10MostAbTaxaPerSite_Species_Conf_0.0.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

# Confidence 0.1
#===============

# Number of species to plot

Nspecies <- 10

# Conf 0.1
#---------

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha[,!(colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1"))] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

# End of added on 2024-06-03

# Make an empty list to store the plots

ListPlots <- list()

# Get the list of all sites

ListSites <- unique(DF_MergedMetadata$siteID)[which(is.na(unique(DF_MergedMetadata$siteID)) == FALSE)]

# Remove FLA and FLM for now as they do not have metagenomic data

ListSites <- c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT") # ListSites[which(!(ListSites %in% c("FLA", "FLM")))]

# Loop on all sites

for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the 10 most abundant species names
  
  SubDF <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  NamesSpecies <- names(sort(rowSums(SubDF), decreasing = TRUE))[1:Nspecies]
  
  # Extract the right species in the SubDF, transpose the dataframe and merge with the appropriate metadata
  
  SubDF_Species <- SubDF[row.names(SubDF) %in% NamesSpecies, , drop = FALSE] # Check that drop does not disturb
  SubDF_Species_t <- data.frame(t(SubDF_Species))
  SubDF_Species_t$SampleName <- row.names(SubDF_Species_t)
  
  SubDF_Species_t_metadata <- left_join(SubDF_Species_t, DF_MergedMetadata)
  
  # Convert the date to type date
  
  SubDF_Species_t_metadata$Date_num <- as.Date(SubDF_Species_t_metadata$date, "%Y-%m-%d")
  
  # Keep only the columns for the plot, then convert to long format 
  
  SubDF_Species_t_metadata_ForPlot <- SubDF_Species_t_metadata[, colnames(SubDF_Species_t_metadata) %in% c("Date_num", "SampleName", gsub(" ", ".", NamesSpecies))]
  
  SubDF_Species_t_metadata_ForPlot_Long <- SubDF_Species_t_metadata_ForPlot %>%
    pivot_longer(
      cols = gsub(" ", ".", NamesSpecies),
      names_to = "Species",
      values_to = "RelAbundance",
    )
  
  # Make the plot
  
  print(ListSites[i])
  ListPlots[[ListSites[i]]] <- ggplot(data = SubDF_Species_t_metadata_ForPlot_Long, aes(x = Date_num, y = RelAbundance, color = Species)) +
    geom_point() +
    geom_line() +
    scale_colour_manual(values = CustomPalette) +
    scale_x_date(date_breaks="1 month", date_labels="%b\n%Y", limits = c(min = as.Date("2021-12-02", "%Y-%m-%d"), max = as.Date("2022-12-06", "%Y-%m-%d"))) +
    labs(x = "Date", y = "Relative abundance", title = ListSites[i]) +
    theme_bw()
  
}

p1 <- gridExtra::grid.arrange(grobs = ListPlots)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_Conf_0.1.pdf", plot=p1, device = cairo_pdf(), width=49, height=49)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_Conf_0.1.png", plot=p1, device = "png", dpi = 300, width=49, height=49)
dev.off()

#
# Add the plots for confidence 0.1 and 0.6
#



#====================================================#
# Barplots with the n most abundant species per site #
#====================================================#

# Number of species to plot

Nspecies <- 100

# Conf 0.1
#---------

# Make an empty list to store the plots

ListPlots <- list()

# Get the list of all sites

ListSites <- unique(DF_MergedMetadata$siteID)[which(is.na(unique(DF_MergedMetadata$siteID)) == FALSE)]

# Remove FLA and FLM for now as they do not have metagenomic data

ListSites <- ListSites[which(!(ListSites %in% c("FLA", "FLM", "FBB", "FCF", "FCM", "FMD")))] # Before: ListSites <- ListSites[which(!(ListSites %in% c("FLA", "FLM")))]

# Loop on all sites

for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  MostAbSpecies <- sort(rowSums(SubDF), decreasing = TRUE)[1:Nspecies]
  
  # Plot the most abundant species in the form of a barplot
  
  MostAbSpecies_DF <- data.frame("Species" = names(MostAbSpecies), "RelAb" = MostAbSpecies)
  MostAbSpecies_DF$Species <- factor(MostAbSpecies_DF$Species, levels = names(MostAbSpecies))
  
  ListPlots[[ListSites[i]]] <- ggplot(data = MostAbSpecies_DF, aes(x = Species, y = RelAb)) +
    geom_bar(stat="identity") +
    labs(x = "Species", y = "Summed relative abundance", title = ListSites[i]) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme_bw()
  
}

p1 <- gridExtra::grid.arrange(grobs = ListPlots)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_SummedAllSamples_Conf_0.1.pdf", plot=p1, device = cairo_pdf(), width=49, height=49)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_SummedAllSamples_Conf_0.1.png", plot=p1, device = "png", dpi = 300, width=49, height=49)
dev.off()

# Take the first 10 species per site

Nspecies <- 10

ListMostAbundSpecies <- list()
for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  MostAbSpecies <- sort(rowSums(SubDF), decreasing = TRUE)[1:Nspecies]
  ListMostAbundSpecies[[ListSites[i]]] <- MostAbSpecies
}

# List all the names of the species that appear within the 10 most abundant species in at least 1 site
ListNmostAbAllSites <- unique(unlist(lapply(ListMostAbundSpecies, names)))

# For each site extract the average relative abundance of these species

ListMostAbundSpecies <- list()
for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  #MostAbSpecies <- sort(rowSums(SubDF), decreasing = TRUE)[1:Nspecies]
  MostAbSpecies <- SubDF[ListNmostAbAllSites,]
  NumberSamples <- ncol(MostAbSpecies)
  MostAbSpeciesSums <- rowSums(MostAbSpecies)/NumberSamples
  ListMostAbundSpecies[[ListSites[i]]] <- MostAbSpeciesSums
}

# Convert the list to dataframe, then convert to long format

DF_RelAbMostAbSpecForPlot <- as.data.frame(ListMostAbundSpecies)
DF_RelAbMostAbSpecForPlot_t <- data.frame(t(DF_RelAbMostAbSpecForPlot))
NamesFactor <- names(sort(colSums(DF_RelAbMostAbSpecForPlot_t), decreasing = TRUE))
DF_RelAbMostAbSpecForPlot_t$Site <- row.names(DF_RelAbMostAbSpecForPlot_t)
DF_RelAbMostAbSpecForPlot_t_longer <- pivot_longer(DF_RelAbMostAbSpecForPlot_t, cols = !Site, names_to = "Species", values_to = "RelAbundance")
DF_RelAbMostAbSpecForPlot_t_longer$Species <- factor(DF_RelAbMostAbSpecForPlot_t_longer$Species, levels = NamesFactor)
DF_RelAbMostAbSpecForPlot_t_longer$Site <- factor(DF_RelAbMostAbSpecForPlot_t_longer$Site, levels = c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"))

# Convert relative abundance to % of all reads

DF_RelAbMostAbSpecForPlot_t_longer$RelAbundance <- DF_RelAbMostAbSpecForPlot_t_longer$RelAbundance * 100

# Make the plot

p1 <- ggplot(data = DF_RelAbMostAbSpecForPlot_t_longer, aes(x = Species, y = RelAbundance, fill = Site)) +
  geom_bar(position='dodge', stat='identity') +
  #scale_fill_manual(values = c("deepskyblue", "gold", "gold2", "darkolivegreen1", "darkolivegreen3", "darkolivegreen4", "forestgreen")) +
  scale_fill_manual(values = c("darkorchid2", "gold", "goldenrod4", "deepskyblue", "dodgerblue4", "darkolivegreen1", "forestgreen")) +
  labs(x = "Species", y = "Mean relative abundance (% of all reads)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/Fig10MostAbTaxaPerSite_Species_Conf_0.1.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/Fig10MostAbTaxaPerSite_Species_Conf_0.1.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

# Same, but in fraction of all reads instead of in percents

DF_RelAbMostAbSpecForPlot_t_longer$RelAbundance_frac <- DF_RelAbMostAbSpecForPlot_t_longer$RelAbundance / 100

# Make the plot

p1 <- ggplot(data = DF_RelAbMostAbSpecForPlot_t_longer, aes(x = Species, y = RelAbundance_frac, fill = Site)) +
  geom_bar(position='dodge', stat='identity') +
  #scale_fill_manual(values = c("deepskyblue", "gold", "gold2", "darkolivegreen1", "darkolivegreen3", "darkolivegreen4", "forestgreen")) +
  scale_fill_manual(values = c("darkorchid2", "gold", "goldenrod4", "deepskyblue", "dodgerblue4", "darkolivegreen1", "forestgreen")) +
  labs(x = "Species", y = "Mean relative abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/Fig10MostAbTaxaPerSite_Species_Conf_0.1_fraction.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/Fig10MostAbTaxaPerSite_Species_Conf_0.1_fraction.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

# Added on 2024-06-26

#=========================================================================================#
# Make a heatmap using only the 36 species that are amongst the top 10 of at least 1 site #
#=========================================================================================#

# List species to select

ListOfSpeciesHeatmap_0.1 <- ListNmostAbAllSites

# Make dataframe 

DF_ForHeatmap_0.1 <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm
DF_ForHeatmap_0.1_t <- as.data.frame(t(DF_ForHeatmap_0.1))

# Select good sites and species

DF_ForHeatmap_0.1_t_GoodCol <- DF_ForHeatmap_0.1_t[!(row.names(DF_ForHeatmap_0.1_t) %in% c("name", "taxRank", "taxID", "Max", "lineage")), colnames(DF_ForHeatmap_0.1_t) %in% ListOfSpeciesHeatmap_0.1]

DF_ForHeatmap_0.1_t_GoodCol$Site <- substr(row.names(DF_ForHeatmap_0.1_t_GoodCol), 3, 5)

DF_ForHeatmap_0.1_t_GoodCol_GoodSite <- DF_ForHeatmap_0.1_t_GoodCol[DF_ForHeatmap_0.1_t_GoodCol[["Site"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]
DF_ForHeatmap_0.1_t_GoodCol_GoodSite$Site <- NULL

# Convert all columns to numeric

DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num <- DF_ForHeatmap_0.1_t_GoodCol_GoodSite %>% mutate_if(is.character,as.numeric) # Check
pheatmap(DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num, cutree_rows = 10, cluster_cols = TRUE, color = paletteer_c("grDevices::Purples 3", 30, direction = -1))

DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num_sqrt <- sqrt(DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num)
pheatmap(DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num_sqrt, cutree_rows = 10, cluster_cols = TRUE, color = paletteer_c("grDevices::Purples 3", 30, direction = -1))

DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num_sd <- data.frame(lapply(DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num, function(x) x/sd(x)))
row.names(DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num_sd) <- row.names(DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num) # To be absolutely checked
pheatmap(DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num_sd, cutree_rows = 10, cluster_cols = TRUE, color = paletteer_c("grDevices::Purples 3", 30, direction = -1))

DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num_NormRow <- DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num/rowSums(DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num) # Check this thoroughly
row.names(DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num_NormRow) <- substr(row.names(DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num_NormRow), 3, 7)
p1 <- pheatmap(DF_ForHeatmap_0.1_t_GoodCol_GoodSite_Num_NormRow, cutree_rows = 10, cluster_cols = TRUE, color = paletteer_c("grDevices::Purples 3", 30, direction = -1), fontsize = 20)
p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/HeatmapNorm36spec_Species_Conf_0.1.pdf", plot=p1, device = cairo_pdf(), width=20, height=30)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/HeatmapNorm36spec_Species_Conf_0.1.png", plot=p1, device = "png", dpi = 300, width=20, height=20)
dev.off()

#=====================================================================================================# 
# Calculate the 10 most abundant species per sample instead of per site, according to Louis's comment #
#=====================================================================================================# 

Nspecies <- 10

ListMostAbundSpecies_PerSample <- list()

for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  
  for (j in 1:length(colnames(SubDF))) {
    SampleVect <- SubDF[[colnames(SubDF)[j]]] # To be checked
    names(SampleVect) <- row.names(SubDF) # To be checked
    SampleVect_sorted <- sort(SampleVect, decreasing = TRUE)[1:Nspecies]
    ListMostAbundSpecies_PerSample[[colnames(SubDF)[j]]] <- SampleVect_sorted
  }
}

DF_names_MostAbundSpecies_PerSample <- lapply(ListMostAbundSpecies_PerSample, names) # Check
Vect_names_MostAbundSpecies_PerSample <- do.call(c, DF_names_MostAbundSpecies_PerSample) # Check
unique(Vect_names_MostAbundSpecies_PerSample)

#=============================================================================================# 
# Calculate the 10 most abundant species per sample per site but using median instead of mean #
#=============================================================================================# 

# Take the first 10 species per site

Nspecies <- 10

ListMostAbundSpecies <- list()
for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  SubDF_t <- data.frame(t(SubDF))
  MostAbSpecies <- sort(sapply(SubDF_t, median), decreasing = TRUE)[1:Nspecies] # Check
  ListMostAbundSpecies[[ListSites[i]]] <- MostAbSpecies
}

# List all the names of the species that appear within the 10 most abundant species in at least 1 site
ListNmostAbAllSites <- unique(unlist(lapply(ListMostAbundSpecies, names)))

# For each site extract the average relative abundance of these species

ListMostAbundSpecies <- list()
for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  SubDF_t <- data.frame(t(SubDF))
  MostAbSpecies <- SubDF_t[,ListNmostAbAllSites]
  
  MostAbSpeciesSums <- sapply(MostAbSpecies, median) # Check
  ListMostAbundSpecies[[ListSites[i]]] <- MostAbSpeciesSums
}

# Convert the list to dataframe, then convert to long format

DF_RelAbMostAbSpecForPlot <- as.data.frame(ListMostAbundSpecies)
DF_RelAbMostAbSpecForPlot_t <- data.frame(t(DF_RelAbMostAbSpecForPlot))
NamesFactor <- names(sort(colSums(DF_RelAbMostAbSpecForPlot_t), decreasing = TRUE))
DF_RelAbMostAbSpecForPlot_t$Site <- row.names(DF_RelAbMostAbSpecForPlot_t)
DF_RelAbMostAbSpecForPlot_t_longer <- pivot_longer(DF_RelAbMostAbSpecForPlot_t, cols = !Site, names_to = "Species", values_to = "RelAbundance")
DF_RelAbMostAbSpecForPlot_t_longer$Species <- factor(DF_RelAbMostAbSpecForPlot_t_longer$Species, levels = NamesFactor)
DF_RelAbMostAbSpecForPlot_t_longer$Site <- factor(DF_RelAbMostAbSpecForPlot_t_longer$Site, levels = c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"))

# Convert relative abundance to % of all reads

DF_RelAbMostAbSpecForPlot_t_longer$RelAbundance <- DF_RelAbMostAbSpecForPlot_t_longer$RelAbundance * 100

# Make the plot

p1 <- ggplot(data = DF_RelAbMostAbSpecForPlot_t_longer, aes(x = Species, y = RelAbundance, fill = Site)) +
  geom_bar(position='dodge', stat='identity') +
  #scale_fill_manual(values = c("deepskyblue", "gold", "gold2", "darkolivegreen1", "darkolivegreen3", "darkolivegreen4", "forestgreen")) +
  scale_fill_manual(values = c("darkorchid2", "gold", "goldenrod4", "deepskyblue", "dodgerblue4", "darkolivegreen1", "forestgreen")) +
  labs(x = "Species", y = "Median relative abundance (% of all reads)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))

plot(p1)

# End of added on 2024-06-26

# Confidence 0.6
#===============

# Number of species to plot

Nspecies <- 10

# Conf 0.6
#---------

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_SamRm <- DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha[,!(colnames(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1"))] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

# End of added on 2024-06-03

# Make an empty list to store the plots

ListPlots <- list()

# Get the list of all sites

ListSites <- unique(DF_MergedMetadata$siteID)[which(is.na(unique(DF_MergedMetadata$siteID)) == FALSE)]

# Remove FLA and FLM for now as they do not have metagenomic data

ListSites <- c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT") # ListSites[which(!(ListSites %in% c("FLA", "FLM")))]

# Loop on all sites

for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the 10 most abundant species names
  
  SubDF <- DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  NamesSpecies <- names(sort(rowSums(SubDF), decreasing = TRUE))[1:Nspecies]
  
  # Extract the right species in the SubDF, transpose the dataframe and merge with the appropriate metadata
  
  SubDF_Species <- SubDF[row.names(SubDF) %in% NamesSpecies, , drop = FALSE] # Check that drop does not disturb
  SubDF_Species_t <- data.frame(t(SubDF_Species))
  SubDF_Species_t$SampleName <- row.names(SubDF_Species_t)
  
  SubDF_Species_t_metadata <- left_join(SubDF_Species_t, DF_MergedMetadata)
  
  # Convert the date to type date
  
  SubDF_Species_t_metadata$Date_num <- as.Date(SubDF_Species_t_metadata$date, "%Y-%m-%d")
  
  # Keep only the columns for the plot, then convert to long format 
  
  # Added as some species names had [] or () in their names # Check this does not disturb the rest
  
  NamesSpecies <- gsub("\\[", "X.", NamesSpecies)
  NamesSpecies <- gsub("\\]", ".", NamesSpecies)
  NamesSpecies <- gsub("\\(", ".", NamesSpecies)
  NamesSpecies <- gsub("\\)", ".", NamesSpecies)
  
  # End of Added as some species names had [] or () in their names
  
  SubDF_Species_t_metadata_ForPlot <- SubDF_Species_t_metadata[, colnames(SubDF_Species_t_metadata) %in% c("Date_num", "SampleName", gsub(" ", ".", NamesSpecies))]
  
  SubDF_Species_t_metadata_ForPlot_Long <- SubDF_Species_t_metadata_ForPlot %>%
    pivot_longer(
      cols = gsub(" ", ".", NamesSpecies),
      names_to = "Species",
      values_to = "RelAbundance",
    )
  
  # Make the plot
  
  print(ListSites[i])
  ListPlots[[ListSites[i]]] <- ggplot(data = SubDF_Species_t_metadata_ForPlot_Long, aes(x = Date_num, y = RelAbundance, color = Species)) +
    geom_point() +
    geom_line() +
    scale_colour_manual(values = CustomPalette) +
    scale_x_date(date_breaks="1 month", date_labels="%b\n%Y", limits = c(min = as.Date("2021-12-02", "%Y-%m-%d"), max = as.Date("2022-12-06", "%Y-%m-%d"))) +
    labs(x = "Date", y = "Relative abundance", title = ListSites[i]) +
    theme_bw()
  
}

p1 <- gridExtra::grid.arrange(grobs = ListPlots)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_Conf_0.6.pdf", plot=p1, device = cairo_pdf(), width=49, height=49)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_Conf_0.6.png", plot=p1, device = "png", dpi = 300, width=49, height=49)
dev.off()

#
# Add the plots for confidence 0.1 and 0.6
#



#====================================================#
# Barplots with the n most abundant species per site #
#====================================================#

# Number of species to plot

Nspecies <- 100

# Conf 0.6
#---------

# Make an empty list to store the plots

ListPlots <- list()

# Get the list of all sites

ListSites <- unique(DF_MergedMetadata$siteID)[which(is.na(unique(DF_MergedMetadata$siteID)) == FALSE)]

# Remove FLA and FLM for now as they do not have metagenomic data

ListSites <- ListSites[which(!(ListSites %in% c("FLA", "FLM", "FBB", "FCF", "FCM", "FMD")))] # Before: ListSites <- ListSites[which(!(ListSites %in% c("FLA", "FLM")))]

# Loop on all sites

for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  MostAbSpecies <- sort(rowSums(SubDF), decreasing = TRUE)[1:Nspecies]
  
  # Plot the most abundant species in the form of a barplot
  
  MostAbSpecies_DF <- data.frame("Species" = names(MostAbSpecies), "RelAb" = MostAbSpecies)
  MostAbSpecies_DF$Species <- factor(MostAbSpecies_DF$Species, levels = names(MostAbSpecies))
  
  ListPlots[[ListSites[i]]] <- ggplot(data = MostAbSpecies_DF, aes(x = Species, y = RelAb)) +
    geom_bar(stat="identity") +
    labs(x = "Species", y = "Summed relative abundance", title = ListSites[i]) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  theme_bw()
  
}

p1 <- gridExtra::grid.arrange(grobs = ListPlots)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_SummedAllSamples_Conf_0.6.pdf", plot=p1, device = cairo_pdf(), width=49, height=49)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/MostAbundantTaxa/Species/FigMostAbTaxaAmongstAllReads_Species_SummedAllSamples_Conf_0.6.png", plot=p1, device = "png", dpi = 300, width=49, height=49)
dev.off()

# Take the first 10 species per site

Nspecies <- 10

ListMostAbundSpecies <- list()
for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  MostAbSpecies <- sort(rowSums(SubDF), decreasing = TRUE)[1:Nspecies]
  ListMostAbundSpecies[[ListSites[i]]] <- MostAbSpecies
}

# List all the names of the species that appear within the 10 most abundant species in at least 1 site
ListNmostAbAllSites <- unique(unlist(lapply(ListMostAbundSpecies, names)))

# For each site extract the average relative abundance of these species

ListMostAbundSpecies <- list()
for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe, sum per species, and extract the n most abundant species 
  
  SubDF <- DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.6_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  #MostAbSpecies <- sort(rowSums(SubDF), decreasing = TRUE)[1:Nspecies]
  MostAbSpecies <- SubDF[ListNmostAbAllSites,]
  NumberSamples <- ncol(MostAbSpecies)
  MostAbSpeciesSums <- rowSums(MostAbSpecies)/NumberSamples
  ListMostAbundSpecies[[ListSites[i]]] <- MostAbSpeciesSums
}

# Convert the list to dataframe, then convert to long format

DF_RelAbMostAbSpecForPlot <- as.data.frame(ListMostAbundSpecies)
DF_RelAbMostAbSpecForPlot_t <- data.frame(t(DF_RelAbMostAbSpecForPlot))
NamesFactor <- names(sort(colSums(DF_RelAbMostAbSpecForPlot_t), decreasing = TRUE))
DF_RelAbMostAbSpecForPlot_t$Site <- row.names(DF_RelAbMostAbSpecForPlot_t)
DF_RelAbMostAbSpecForPlot_t_longer <- pivot_longer(DF_RelAbMostAbSpecForPlot_t, cols = !Site, names_to = "Species", values_to = "RelAbundance")
DF_RelAbMostAbSpecForPlot_t_longer$Species <- factor(DF_RelAbMostAbSpecForPlot_t_longer$Species, levels = NamesFactor)
DF_RelAbMostAbSpecForPlot_t_longer$Site <- factor(DF_RelAbMostAbSpecForPlot_t_longer$Site, levels = c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"))

# Convert relative abundance to % of all reads

DF_RelAbMostAbSpecForPlot_t_longer$RelAbundance <- DF_RelAbMostAbSpecForPlot_t_longer$RelAbundance * 100

# Make the plot

p1 <- ggplot(data = DF_RelAbMostAbSpecForPlot_t_longer, aes(x = Species, y = RelAbundance, fill = Site)) +
  geom_bar(position='dodge', stat='identity') +
  #scale_fill_manual(values = c("deepskyblue", "gold", "gold2", "darkolivegreen1", "darkolivegreen3", "darkolivegreen4", "forestgreen")) +
  scale_fill_manual(values = c("darkorchid2", "gold", "goldenrod4", "deepskyblue", "dodgerblue4", "darkolivegreen1", "forestgreen")) +
  labs(x = "Species", y = "Mean relative abundance (% of all reads)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))

plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/Fig10MostAbTaxaPerSite_Species_Conf_0.6.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/Fig10MostAbTaxaPerSite_Species_Conf_0.6.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

#============================================================================================================#
# Calculate how many species are required to gather at least xx % of the total of fungal reads for each site #
#============================================================================================================#

# Conf 0.1
#---------

# Create a list to hold for each site the number of secies required to gather up 50 and 90 % of the fungal reads of a site

Listn50n90 <- list()

# Get the list of all sites

ListSites <- unique(DF_MergedMetadata$siteID)[which(is.na(unique(DF_MergedMetadata$siteID)) == FALSE)]

# Remove FLA and FLM for now as they do not have metagenomic data

ListSites <- c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT") # ListSites[which(!(ListSites %in% c("FLA", "FLM")))]

# Loop on all sites

for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe
  
  SubDF <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  
  # Calculate total relative abundance of Fungi amongst all samples from the good site
  
  SumEachSample <- colSums(SubDF)
  Names_SumEachSample <- names(SumEachSample)
  SumPercSite <- sum(DF_MergedMetadataPavian[(DF_MergedMetadataPavian[["SampleName"]] %in% Names_SumEachSample) & (DF_MergedMetadataPavian[["Confidence"]] == "0.1"), "Fungal.reads"]/100) # Check this
  
  # Calculate the sum for each species
  
  SortedSpecSums <- sort(rowSums(SubDF), decreasing = TRUE)
  SortedSpecSumsPerc <- 100*(SortedSpecSums/SumPercSite)
  
  CumSumSpec <- cumsum(SortedSpecSumsPerc)
  n50 <- length(CumSumSpec[which(CumSumSpec < 50)]) +1 # +1 because we want the mimimun of species to represent 50 % of fungal reads
  n90 <- length(CumSumSpec[which(CumSumSpec < 90)]) +1 # +1 because we want the mimimun of species to represent 90 % of fungal reads
  
  VectOut <- c(n50, n90)
  names(VectOut) <- c("n50", "n90")
  
  Listn50n90[[ListSites[i]]] <- VectOut
}

bind_rows(Listn50n90)

# Same, but sample per sample to make barplots
#---------------------------------------------

DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm_ColNum <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm[,!(colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]

ListCumSums_0.1 <- list()
for (i in 1:length(colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm_ColNum))) {
  NameSample <- colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm_ColNum)[i]
  RelAbFungi <- (DF_MergedMetadataPavian[(DF_MergedMetadataPavian[["SampleName"]] %in% NameSample) & (DF_MergedMetadataPavian[["Confidence"]] == "0.1"), "Fungal.reads"])/100
  SortedSpeciesSample <- sort(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm_ColNum[[NameSample]], decreasing = TRUE)
  SortedSpeciesSample_Div <- SortedSpeciesSample/RelAbFungi
  SortedSpeciesSample_Div_CumSum <- cumsum(SortedSpeciesSample_Div)
  ListCumSums_0.1[[NameSample]] <- SortedSpeciesSample_Div_CumSum
}

# Convert to dataframe
DF_ListCumSums_0.1 <- as.data.frame(ListCumSums_0.1)

# Add a Rank column, then convert to long
DF_ListCumSums_0.1$Rank <- row.names(DF_ListCumSums_0.1)

DF_ListCumSums_0.1_Long <- DF_ListCumSums_0.1 %>%
  pivot_longer(
    cols = !Rank,
    names_to = "Sample",
    values_to = "CumRelAbundance",
  )

# Select good sites
DF_ListCumSums_0.1_Long$Site <- substr(DF_ListCumSums_0.1_Long$Sample, 3, 5)
DF_ListCumSums_0.1_Long_GoodCol <- DF_ListCumSums_0.1_Long[DF_ListCumSums_0.1_Long[["Site"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT") ,]

# Calculate mean by site and by rank

DF_ListCumSums_0.1_Long_GoodCol_mean <- aggregate(DF_ListCumSums_0.1_Long_GoodCol$CumRelAbundance, list(DF_ListCumSums_0.1_Long_GoodCol$Site, DF_ListCumSums_0.1_Long_GoodCol$Rank), FUN=mean) # Hard check
colnames(DF_ListCumSums_0.1_Long_GoodCol_mean) <- c("Site", "Rank", "Mean_CumRelAb")

DF_ListCumSums_0.1_Long_GoodCol_sd <- aggregate(DF_ListCumSums_0.1_Long_GoodCol$CumRelAbundance, list(DF_ListCumSums_0.1_Long_GoodCol$Site, DF_ListCumSums_0.1_Long_GoodCol$Rank), FUN=sd) # Hard check
colnames(DF_ListCumSums_0.1_Long_GoodCol_sd) <- c("Site", "Rank", "Sd_CumRelAb")

DF_ListCumSums_0.1_Long_GoodCol_median <- aggregate(DF_ListCumSums_0.1_Long_GoodCol$CumRelAbundance, list(DF_ListCumSums_0.1_Long_GoodCol$Site, DF_ListCumSums_0.1_Long_GoodCol$Rank), FUN=median) # Hard check
colnames(DF_ListCumSums_0.1_Long_GoodCol_median) <- c("Site", "Rank", "Median_CumRelAb")

DF_ListCumSums_0.1_Long_GoodCol_mean_sd <- list(DF_ListCumSums_0.1_Long_GoodCol_mean, DF_ListCumSums_0.1_Long_GoodCol_sd, DF_ListCumSums_0.1_Long_GoodCol_median) %>% reduce(full_join) # To be checked
DF_ListCumSums_0.1_Long_GoodCol_mean_sd$Rank <- as.numeric(DF_ListCumSums_0.1_Long_GoodCol_mean_sd$Rank) # Check

ggplot(DF_ListCumSums_0.1_Long_GoodCol_mean_sd, aes(x=Rank, y = Mean_CumRelAb, group = Site)) +
  geom_line(aes(color = Site)) +
  #geom_point(aes(color = Site)) +
  geom_pointrange(aes(ymin=Mean_CumRelAb-Sd_CumRelAb, ymax=Mean_CumRelAb+Sd_CumRelAb, color = Site)) +
  xlim(1, 60) +
  theme_bw()

ggplot(DF_ListCumSums_0.1_Long_GoodCol_mean_sd, aes(x=Rank, y = Median_CumRelAb, group = Site)) +
  geom_line(aes(color = Site)) +
  geom_point(aes(color = Site)) +
  #geom_pointrange(aes(ymin=Mean_CumRelAb-Sd_CumRelAb, ymax=Median_CumRelAb+Sd_CumRelAb, color = Site)) +
  xlim(1, 60) +
  theme_bw()

ggplot(DF_ListCumSums_0.1_Long_GoodCol_mean_sd, aes(x=Rank, y = Mean_CumRelAb, group = Site)) +
  geom_line(aes(color = Site)) +
  geom_point(aes(color = Site)) +
  #geom_pointrange(aes(ymin=Mean_CumRelAb-Sd_CumRelAb, ymax=Median_CumRelAb+Sd_CumRelAb, color = Site)) +
  xlim(1, 60) +
  theme_bw()

ggplot(DF_ListCumSums_0.1_Long_GoodCol_mean_sd, aes(x=Rank, y = Mean_CumRelAb, group = Site)) +
  geom_line(aes(color = Site)) +
  geom_point(aes(color = Site)) +
  geom_line(aes(x=Rank, y=Median_CumRelAb+Sd_CumRelAb, color = Site), linetype = "dashed") +
  geom_line(aes(x=Rank, y=Median_CumRelAb-Sd_CumRelAb, color = Site), linetype = "dashed") +
  #geom_pointrange(aes(ymin=Mean_CumRelAb-Sd_CumRelAb, ymax=Median_CumRelAb+Sd_CumRelAb, color = Site)) +
  xlim(1, 60) +
  theme_bw()

DF_ListCumSums_0.1_Long_GoodCol_summary <- aggregate(DF_ListCumSums_0.1_Long_GoodCol$CumRelAbundance, list(DF_ListCumSums_0.1_Long_GoodCol$Site, DF_ListCumSums_0.1_Long_GoodCol$Rank), FUN=summary) # Hard check
DF_ListCumSums_0.1_Long_GoodCol_summary <- DF_ListCumSums_0.1_Long_GoodCol_summary %>%
  mutate_at("x", as.data.frame) |> 
  unnest(x) # Check this thoroughly
colnames(DF_ListCumSums_0.1_Long_GoodCol_summary) <- c("Site", "Rank", "Min", "FirstQuartile", "Median", "Mean", "ThirdQuartile", "Max")
DF_ListCumSums_0.1_Long_GoodCol_summary$Rank <- as.numeric(DF_ListCumSums_0.1_Long_GoodCol_summary$Rank) # Check
DF_ListCumSums_0.1_Long_GoodCol_summary$Site <- factor(DF_ListCumSums_0.1_Long_GoodCol_summary$Site, levels = c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"))

p1 <- ggplot(DF_ListCumSums_0.1_Long_GoodCol_summary, aes(x=Rank, y = Median, group = Site)) +
  geom_line(aes(color = Site)) +
  geom_point(aes(color = Site)) +
  geom_line(aes(x=Rank, y=FirstQuartile, color = Site), linetype = "dotted") +
  geom_line(aes(x=Rank, y=ThirdQuartile, color = Site), linetype = "dashed") +
  #scale_color_manual(values = c("deepskyblue", "gold", "gold2", "darkolivegreen1", "darkolivegreen3", "darkolivegreen4", "forestgreen")) +
  scale_color_manual(values = c("darkorchid2", "gold", "goldenrod4", "deepskyblue", "dodgerblue4", "darkolivegreen1", "forestgreen")) +
  geom_hline(yintercept = 0.5, color = "black") +
  ylab(label = "Fraction of the fungal reads") +
  xlab(label = "Number of fungal species (ordered by decreasing relative abundance)") +
  #geom_pointrange(aes(ymin=Mean_CumRelAb-Sd_CumRelAb, ymax=Median_CumRelAb+Sd_CumRelAb, color = Site)) +
  xlim(1, 65) +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)
  )

plot(p1)
#ggplotly(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/NumberCumulSpecies_0.1.pdf", plot=p1, device = cairo_pdf(), width=12, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/NumberCumulSpecies_0.1.png", plot=p1, device = "png", dpi = 300, width=12, height=10)
dev.off()

# Calculate the number of reads that were affected at the species level
#======================================================================

# Conf 0.1
#---------

# Create a list to hold for each site the number of secies required to gather up 50 and 90 % of the fungal reads of a site

ListDFPercFungi <- list()

# Get the list of all sites

ListSites <- unique(DF_MergedMetadata$siteID)[which(is.na(unique(DF_MergedMetadata$siteID)) == FALSE)]

# Remove FLA and FLM for now as they do not have metagenomic data

ListSites <- c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT") # ListSites[which(!(ListSites %in% c("FLA", "FLM")))]

# Loop on all sites

for (i in 1:length(ListSites)) {
  
  # Get the columns that are only from given site
  
  ListColnames <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] == ListSites[i],"SampleName"]
  ListColnames <- ListColnames[which(is.na(ListColnames) == FALSE)]
  
  # Extract the columns from the right site in the dataframe
  
  SubDF <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm[,colnames(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_SamRm) %in% ListColnames, drop = FALSE] # Check that drop does not disturb
  
  # Calculate total relative abundance of Fungi amongst all samples from the good site
  
  SumEachSample <- colSums(SubDF)
  Names_SumEachSample <- names(SumEachSample)
  
  DF_SumFungiSamples <- data.frame(FungalReadsSpecies = SumEachSample, SampleName = names(SumEachSample))
  colSums(SubDF)
  DF_NumberAttributedAllFungi <- DF_MergedMetadataPavian[(DF_MergedMetadataPavian[["SampleName"]] %in% Names_SumEachSample) & (DF_MergedMetadataPavian[["Confidence"]] == "0.1"), c("Fungal.reads", "SampleName")]
  
  DF_CompReadsFungalSpecAll <- merge(DF_SumFungiSamples, DF_NumberAttributedAllFungi, by = "SampleName")
  DF_CompReadsFungalSpecAll$FungalReadsAll <- DF_CompReadsFungalSpecAll$Fungal.reads/100
  DF_CompReadsFungalSpecAll$Fungal.reads <- NULL
  DF_CompReadsFungalSpecAll$PercReadsSpecAmongstFung <- (DF_CompReadsFungalSpecAll$FungalReadsSpecies/DF_CompReadsFungalSpecAll$FungalReadsAll)*100
  
  ListDFPercFungi[[ListSites[i]]] <- DF_CompReadsFungalSpecAll
}

ListDFPercFungi_AllSites <- bind_rows(ListDFPercFungi, .id = "column_label")
mean(ListDFPercFungi_AllSites$PercReadsSpecAmongstFung)
min(ListDFPercFungi_AllSites$PercReadsSpecAmongstFung)
max(ListDFPercFungi_AllSites$PercReadsSpecAmongstFung)
sd(ListDFPercFungi_AllSites$PercReadsSpecAmongstFung)

#======================================================================================================================#
# Find the number of species that are above a certain number of reads in at least one sample for each confidence level #
# (to be more accurate, the best option is to have a curve that for n thresholds with n high on x axis displays the    #
# number of species equal or above this threshold)                                                                     #
#======================================================================================================================#

# Number of values to evaluate 

Xseq <- seq(from = 1, to = 10000, by=1) # seq(from = 0, to = 100000, by=10)

# Conf 0.0
#---------

DF_PavianSpecies_Conf_0.0_ForNumbSpecies <- DF_PavianSpecies_Conf_0.0
row.names(DF_PavianSpecies_Conf_0.0_ForNumbSpecies) <- DF_PavianSpecies_Conf_0.0_ForNumbSpecies$name
DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols <- DF_PavianSpecies_Conf_0.0_ForNumbSpecies[, !(colnames(DF_PavianSpecies_Conf_0.0_ForNumbSpecies) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]
DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t <- data.frame(t(DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols))
DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t$AGTUIndex <- row.names(DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t)
DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged <- left_join(DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t, DF_MergedMetadata)
DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged$Country <- substr(DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged$SampleShort, 1, 1)
DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged_F <- DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged[DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged[["Country"]] == "F",]
row.names(DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged_F) <- DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged_F$SampleName # Check names are unique
DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged_F_GoodCols <- DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged_F[, !(colnames(DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged_F) %in% c(colnames(DF_MergedMetadata), "Country"))]

DF_ForSpeciesNumber_0.0 <- DF_PavianSpecies_Conf_0.0_ForNumbSpecies_GoodCols_t_merged_F_GoodCols %>% replace(is.na(.), 0)

MaxPerSpecies_0.0 <- sapply(DF_ForSpeciesNumber_0.0, max)

NumbSpecies_0.0 <- integer(0) 

for (i in 1:length(Xseq)) {
  NumbSpecies_0.0[i] <- length(which(MaxPerSpecies_0.0 >= Xseq[i]))
}

DF_NumbSpec_0.0 <- data.frame("Threshold"=Xseq, "NumberSpecies"=NumbSpecies_0.0)
DF_NumbSpec_0.0$Confidence <- "0.0"

# Conf 0.1
#---------

DF_PavianSpecies_Conf_0.1_ForNumbSpecies <- DF_PavianSpecies_Conf_0.1
row.names(DF_PavianSpecies_Conf_0.1_ForNumbSpecies) <- DF_PavianSpecies_Conf_0.1_ForNumbSpecies$name
DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols <- DF_PavianSpecies_Conf_0.1_ForNumbSpecies[, !(colnames(DF_PavianSpecies_Conf_0.1_ForNumbSpecies) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]
DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t <- data.frame(t(DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols))
DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t$AGTUIndex <- row.names(DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t)
DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged <- left_join(DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t, DF_MergedMetadata)
DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged$Country <- substr(DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged$SampleShort, 1, 1)
DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged_F <- DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged[DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged[["Country"]] == "F",]
row.names(DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged_F) <- DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged_F$SampleName # Check names are unique
DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged_F_GoodCols <- DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged_F[, !(colnames(DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged_F) %in% c(colnames(DF_MergedMetadata), "Country"))]

DF_ForSpeciesNumber_0.1 <- DF_PavianSpecies_Conf_0.1_ForNumbSpecies_GoodCols_t_merged_F_GoodCols %>% replace(is.na(.), 0)

MaxPerSpecies_0.1 <- sapply(DF_ForSpeciesNumber_0.1, max)

NumbSpecies_0.1 <- integer(0) 

for (i in 1:length(Xseq)) {
  NumbSpecies_0.1[i] <- length(which(MaxPerSpecies_0.1 >= Xseq[i]))
}

DF_NumbSpec_0.1 <- data.frame("Threshold"=Xseq, "NumberSpecies"=NumbSpecies_0.1)
DF_NumbSpec_0.1$Confidence <- "0.1"


# Conf 0.6
#---------

DF_PavianSpecies_Conf_0.6_ForNumbSpecies <- DF_PavianSpecies_Conf_0.6
row.names(DF_PavianSpecies_Conf_0.6_ForNumbSpecies) <- DF_PavianSpecies_Conf_0.6_ForNumbSpecies$name
DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols <- DF_PavianSpecies_Conf_0.6_ForNumbSpecies[, !(colnames(DF_PavianSpecies_Conf_0.6_ForNumbSpecies) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]
DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t <- data.frame(t(DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols))
DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t$AGTUIndex <- row.names(DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t)
DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged <- left_join(DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t, DF_MergedMetadata)
DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged$Country <- substr(DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged$SampleShort, 1, 1)
DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged_F <- DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged[DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged[["Country"]] == "F",]
row.names(DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged_F) <- DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged_F$SampleName # Check names are unique
DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged_F_GoodCols <- DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged_F[, !(colnames(DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged_F) %in% c(colnames(DF_MergedMetadata), "Country"))]

DF_ForSpeciesNumber_0.6 <- DF_PavianSpecies_Conf_0.6_ForNumbSpecies_GoodCols_t_merged_F_GoodCols %>% replace(is.na(.), 0)

MaxPerSpecies_0.6 <- sapply(DF_ForSpeciesNumber_0.6, max)

NumbSpecies_0.6 <- integer(0) 

for (i in 1:length(Xseq)) {
  NumbSpecies_0.6[i] <- length(which(MaxPerSpecies_0.6 >= Xseq[i]))
}

DF_NumbSpec_0.6 <- data.frame("Threshold"=Xseq, "NumberSpecies"=NumbSpecies_0.6)
DF_NumbSpec_0.6$Confidence <- "0.6"

# rbind the 3 dataframes, then plot

DF_NumbSpec_rbind <- rbind(DF_NumbSpec_0.0, DF_NumbSpec_0.1, DF_NumbSpec_0.6)

ggplot(data = DF_NumbSpec_rbind, aes(x = Threshold, y = NumberSpecies, color = Confidence)) +
  geom_line()+
  geom_vline(xintercept = 100) +
  geom_vline(xintercept = 1000) +
  annotate("text", x=1500, y=7500, label= sprintf("At 100 reads:\n-Conf 0.0: %d\n-Conf 0.1: %d\n-Conf 0.6: %d", DF_NumbSpec_0.0[DF_NumbSpec_0.0[["Threshold"]] == 100, "NumberSpecies"], DF_NumbSpec_0.1[DF_NumbSpec_0.1[["Threshold"]] == 100, "NumberSpecies"], DF_NumbSpec_0.6[DF_NumbSpec_0.6[["Threshold"]] == 100, "NumberSpecies"]), size = 6, hjust = 0) +
  annotate("text", x=2000, y=7500, label= sprintf("At 1000 reads:\n-Conf 0.0: %d\n-Conf 0.1: %d\n-Conf 0.6: %d", DF_NumbSpec_0.0[DF_NumbSpec_0.0[["Threshold"]] == 1000, "NumberSpecies"], DF_NumbSpec_0.1[DF_NumbSpec_0.1[["Threshold"]] == 1000, "NumberSpecies"], DF_NumbSpec_0.6[DF_NumbSpec_0.6[["Threshold"]] == 1000, "NumberSpecies"]), size = 6, hjust = 0) +
  annotate("text", x=2500, y=7500, label= sprintf("At 4000 reads:\n-Conf 0.0: %d\n-Conf 0.1: %d\n-Conf 0.6: %d", DF_NumbSpec_0.0[DF_NumbSpec_0.0[["Threshold"]] == 4000, "NumberSpecies"], DF_NumbSpec_0.1[DF_NumbSpec_0.1[["Threshold"]] == 4000, "NumberSpecies"], DF_NumbSpec_0.6[DF_NumbSpec_0.6[["Threshold"]] == 4000, "NumberSpecies"]), size = 6, hjust = 0) +
  xlim(0, 4000) +
  theme_bw()

######################################################
#                        PCA                         #
######################################################

DF_MergedMetadata_GoodSites <- DF_MergedMetadata[DF_MergedMetadata[["siteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]
DF_MergedMetadata_GoodSites <- DF_MergedMetadata_GoodSites[!(DF_MergedMetadata_GoodSites[["AGTUIndex"]] %in% c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59", "AGTU-48")),] # AGTU-48 added on 2024-07-08

DF_MergedMetadata_GoodSites_missMDA <- DF_MergedMetadata_GoodSites[,colnames(DF_MergedMetadata_GoodSites) %in% c("T", "salinity", "chla.mean", "DOC.mean", "protein.like.DOM", "humic.like.DOM")]

# Extrapolate missing values

nb <- estim_ncpPCA(DF_MergedMetadata_GoodSites_missMDA,ncp.max=5) # Output is 1
res.comp_DF_MergedMetadata_GoodSites_missMDA = imputePCA(DF_MergedMetadata_GoodSites_missMDA,ncp=1)
Extrapolated_DF_metaMDA <- data.frame(res.comp_DF_MergedMetadata_GoodSites_missMDA$completeObs)

# Center and scale the columns, then make the PCA

Extrapolated_DF_metaMDA_CentRed <- as.data.frame(lapply(Extrapolated_DF_metaMDA, function (x) (x-mean(x))/sd(x)))

PCA_3umNice <- PCA(Extrapolated_DF_metaMDA_CentRed, scale.unit = FALSE, ncp = 10, graph = TRUE) #, quanti.sup = 11:22)

# Check the metadata dataframe and the extrapolated dataframes are in the same order (using DOC.mean that had no missing values)

identical(Extrapolated_DF_metaMDA$DOC.mean, DF_MergedMetadata_GoodSites$DOC.mean) # Must be TRUE

# Make a more complete plot with sites as colors and seasons as shapes

DF_MergedMetadata_GoodSites$Date_num <- as.Date(DF_MergedMetadata_GoodSites$date, "%Y-%m-%d")
DF_MergedMetadata_GoodSites$Month <- month(DF_MergedMetadata_GoodSites$Date_num)

DF_MergedMetadata_GoodSites$Season <- NA
for (i in 1:dim(DF_MergedMetadata_GoodSites)[1]) {
  if (DF_MergedMetadata_GoodSites[i,"Month"] %in% c(1, 2, 3)) {
    DF_MergedMetadata_GoodSites[i, "Season"] <- "Winter"
  } else if (DF_MergedMetadata_GoodSites[i,"Month"] %in% c(4, 5, 6)) {
    DF_MergedMetadata_GoodSites[i, "Season"] <- "Spring"
  } else if (DF_MergedMetadata_GoodSites[i,"Month"] %in% c(7, 8, 9)) {
    DF_MergedMetadata_GoodSites[i, "Season"] <- "Summer"
  } else if (DF_MergedMetadata_GoodSites[i,"Month"] %in% c(10, 11, 12)) {
    DF_MergedMetadata_GoodSites[i, "Season"] <- "Autumn"
  } else {}
}

PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", #col.var = "contrib",
                                                  legend.title = "Sites", 
                                                  #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"),
                                                  title = "") +
  geom_point(aes(x = x, y = y), colour = "white", fill = "white", size = 2) +
  geom_point(aes(shape = factor(DF_MergedMetadata_GoodSites$Season, levels=c("Spring", "Summer", "Autumn", "Winter")), colour = factor(DF_MergedMetadata_GoodSites$siteID, levels = c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT")))) + 
  stat_ellipse(aes (x = x, y = y, colour = DF_MergedMetadata_GoodSites$siteID)) +
  #scale_color_manual(values = c("deepskyblue", "gold", "gold2", "darkolivegreen1", "darkolivegreen3", "darkolivegreen4", "forestgreen")) + # Added on 2024-06-28
  scale_color_manual(values = c("darkorchid2", "gold", "goldenrod4", "deepskyblue", "dodgerblue4", "darkolivegreen1", "forestgreen")) +
  #scale_colour_continuous(aes(colour = FullExtrapDF_withCyto$PercFungi))) +
  #scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  guides(shape = guide_legend(title = "Season"))#,
         #colour = guide_colorbar(title = "Site"))

plot(PCA_Figure3_CytoFlux_PercFungi)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/PCA_Figure3_FinalPCA.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=11, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/PCA_Figure3_FinalPCA.png", plot=PCA_Figure3_CytoFlux_PercFungi, device = "png", dpi = 300, width=11, height=10)
dev.off()

# hull_data_PCA_0.0 <- 
#   PCA_3umNice %>%
#   drop_na() %>%
#   group_by(SiteID) %>% 
#   slice(chull(NMDS1, NMDS2))
# 
# ggplot(data=data.scores_0.0, aes(x = NMDS1, y = NMDS2, color = SiteID)) +
#   geom_point() +
#   scale_colour_manual(values = CustomPalette) +
#   geom_polygon(data = hull_data_0.0,
#                aes(colour = SiteID),
#                alpha = 0.0,
#                show.legend = FALSE) +
#   #stat_ellipse() +
#   theme_bw()








# /!\
#############################################################################################################
# /!\ Section of PCA just below not working. Beware, something might have been deleted somewhere before /!\ # 
#############################################################################################################
# /!\

# Merge with the big dataframe

DF_3um_ForMergingCyto <- DF_SampleDataForPCA_rmCol_sortedDate_3um

DF_3um_withcyto <- merge(DF_3um_ForMergingCyto, CytometryFileTrunc, by = "Date_num", all.x = TRUE) # Check this merging

# Extrapolate the values with missMDA

DF_3um_withcyto_missMDA <- DF_3um_withcyto[,colnames(DF_3um_withcyto) %in% c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4",  "BaillauryHeigthIn_m", "CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC")]
nb <- estim_ncpPCA(DF_3um_withcyto_missMDA,ncp.max=5) # Output is 1
res.comp_DF_3um_withcyto_missMDA = imputePCA(DF_3um_withcyto_missMDA,ncp=1)
Extrapolated_DF_withCyto <- data.frame(res.comp_DF_3um_withcyto_missMDA$completeObs)




# With percentage of Fungi as colors # /!\ Check the values are in the right order
PCA_3umNice <- PCA(DF_3um_CentRed_woRR7daysDerBaill_10ASVs, scale.unit = FALSE, ncp = 10, graph = TRUE, quanti.sup = 11:22)
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto$PercFungi, #col.var = "contrib",
                                                  legend.title = "Percent Fungi", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"),
                                                  title = "") +
  geom_point(aes(x = x, y = y), colour = "white", fill = "white", size = 2) +
  geom_point(aes(shape = factor(FullExtrapDF_withCyto_ForSeason$Season, levels=c("Spring", "Summer", "Autumn", "Winter")), colour = FullExtrapDF_withCyto$PercFungi)) + 
  #scale_colour_continuous(aes(colour = FullExtrapDF_withCyto$PercFungi))) +
  #scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  guides(shape = guide_legend(title = "Season"),
         colour = guide_colorbar(title = "Percent of Fungi\nAmongst eukaryotes"))

# Modify the PCA to include the season as shape
# PCA_Figure3_CytoFlux_PercFungi$data$Season <- factor(FullExtrapDF_withCyto_ForSeason$Season)
# PCA_Figure3_CytoFlux_PercFungi$layers[[1]]$data$Season <- factor(FullExtrapDF_withCyto_ForSeason$Season)
# PCA_Figure3_CytoFlux_PercFungi$layers[[1]]$mapping <- aes(x, y, colour = Col., shape = Season)
#PCA_Figure3_CytoFlux_PercFungi$layers[[1]]$aes_params$size <- 3
#PCA_Figure3_CytoFlux_PercFungi <- PCA_Figure3_CytoFlux_PercFungi + labs(shape = "Season")
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/PCA_Figure3_FinalPCA.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/PCA_Figure3_FinalPCA.png", plot=PCA_Figure3_CytoFlux_PercFungi, device = "png", dpi = 300, width=15, height=10)
dev.off()


###################################
# List of species above threshold #
###################################

ThreshSpec <- 0.00001

# 0.0

DF_PavianSpecies_Conf_0.0_Norm_Thresh <- DF_PavianSpecies_Conf_0.0_Norm
DF_PavianSpecies_Conf_0.0_Norm_Thresh <- DF_PavianSpecies_Conf_0.0_Norm_Thresh[,!(colnames(DF_PavianSpecies_Conf_0.0_Norm_Thresh) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]
DF_PavianSpecies_Conf_0.0_Norm_Thresh <- as.data.frame(t(DF_PavianSpecies_Conf_0.0_Norm_Thresh))

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

DF_PavianSpecies_Conf_0.0_Norm_Thresh <- DF_PavianSpecies_Conf_0.0_Norm_Thresh[!(row.names(DF_PavianSpecies_Conf_0.0_Norm_Thresh) %in% c("G_CCO01.1_G1", "G_FCT10.6_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

DF_PavianSpecies_Conf_0.0_Norm_Thresh$SiteID <- substr(row.names(DF_PavianSpecies_Conf_0.0_Norm_Thresh), 3, 5)

DF_PavianSpecies_Conf_0.0_Norm_Thresh <- DF_PavianSpecies_Conf_0.0_Norm_Thresh[DF_PavianSpecies_Conf_0.0_Norm_Thresh[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

DF_PavianSpecies_Conf_0.0_Norm_Thresh$SiteID <- NULL

DF_PavianSpecies_Conf_0.0_Norm_Thresh <- as.data.frame(t(DF_PavianSpecies_Conf_0.0_Norm_Thresh))

# End of added on 2024-06-03

DF_PavianSpecies_Conf_0.0_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.0_Norm_Thresh

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
DF_ListSpecAboveThresh_0.0$MetaG_0.0 <- "Yes"

# 0.1

DF_PavianSpecies_Conf_0.1_Norm_Thresh <- DF_PavianSpecies_Conf_0.1_Norm
DF_PavianSpecies_Conf_0.1_Norm_Thresh <- DF_PavianSpecies_Conf_0.1_Norm_Thresh[,!(colnames(DF_PavianSpecies_Conf_0.1_Norm_Thresh) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]
DF_PavianSpecies_Conf_0.1_Norm_Thresh <- as.data.frame(t(DF_PavianSpecies_Conf_0.1_Norm_Thresh))

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

DF_PavianSpecies_Conf_0.1_Norm_Thresh <- DF_PavianSpecies_Conf_0.1_Norm_Thresh[!(row.names(DF_PavianSpecies_Conf_0.1_Norm_Thresh) %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

DF_PavianSpecies_Conf_0.1_Norm_Thresh$SiteID <- substr(row.names(DF_PavianSpecies_Conf_0.1_Norm_Thresh), 3, 5)

DF_PavianSpecies_Conf_0.1_Norm_Thresh <- DF_PavianSpecies_Conf_0.1_Norm_Thresh[DF_PavianSpecies_Conf_0.1_Norm_Thresh[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

DF_PavianSpecies_Conf_0.1_Norm_Thresh$SiteID <- NULL

DF_PavianSpecies_Conf_0.1_Norm_Thresh <- as.data.frame(t(DF_PavianSpecies_Conf_0.1_Norm_Thresh))

# End of added on 2024-06-03

DF_PavianSpecies_Conf_0.1_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.1_Norm_Thresh

ListSpecAboveThresh_0.1 <- list()
for (i in (1:(ncol(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh)))) {
  TheColName <- colnames(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh)[i]
  SpecCounts <- DF_PavianSpecies_Conf_0.1_Norm_AboveThresh[[TheColName]][which(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  SpecNames <- row.names(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh)[which(DF_PavianSpecies_Conf_0.1_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  names(SpecCounts) <- SpecNames
  ListSpecAboveThresh_0.1[[TheColName]] <- SpecCounts
} # Check thoroughly

ListSpecAboveThresh_0.1

unlist(sapply(ListSpecAboveThresh_0.1, function(x) names(x)))
unique(unlist(sapply(ListSpecAboveThresh_0.1, function(x) names(x))))
sort(unique(unlist(sapply(ListSpecAboveThresh_0.1, function(x) names(x)))))

DF_ListSpecAboveThresh_0.1 <- sort(unique(unlist(sapply(ListSpecAboveThresh_0.1, function(x) names(x)))))
DF_ListSpecAboveThresh_0.1 <- data.frame(ListSpecies = DF_ListSpecAboveThresh_0.1)
DF_ListSpecAboveThresh_0.1$MetaG_0.1 <- "Yes"

# 0.6

DF_PavianSpecies_Conf_0.6_Norm_Thresh <- DF_PavianSpecies_Conf_0.6_Norm
DF_PavianSpecies_Conf_0.6_Norm_Thresh <- DF_PavianSpecies_Conf_0.6_Norm_Thresh[,!(colnames(DF_PavianSpecies_Conf_0.6_Norm_Thresh) %in% c("name", "taxRank", "taxID", "Max", "lineage"))]
DF_PavianSpecies_Conf_0.6_Norm_Thresh <- as.data.frame(t(DF_PavianSpecies_Conf_0.6_Norm_Thresh))

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

DF_PavianSpecies_Conf_0.6_Norm_Thresh <- DF_PavianSpecies_Conf_0.6_Norm_Thresh[!(row.names(DF_PavianSpecies_Conf_0.6_Norm_Thresh) %in% c("G_CCO01.1_G1", "G_FCT10.6_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

DF_PavianSpecies_Conf_0.6_Norm_Thresh$SiteID <- substr(row.names(DF_PavianSpecies_Conf_0.6_Norm_Thresh), 3, 5)

DF_PavianSpecies_Conf_0.6_Norm_Thresh <- DF_PavianSpecies_Conf_0.6_Norm_Thresh[DF_PavianSpecies_Conf_0.6_Norm_Thresh[["SiteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

DF_PavianSpecies_Conf_0.6_Norm_Thresh$SiteID <- NULL

DF_PavianSpecies_Conf_0.6_Norm_Thresh <- as.data.frame(t(DF_PavianSpecies_Conf_0.6_Norm_Thresh))

# End of added on 2024-06-03

DF_PavianSpecies_Conf_0.6_Norm_AboveThresh <- DF_PavianSpecies_Conf_0.6_Norm_Thresh

ListSpecAboveThresh_0.6 <- list()
for (i in (1:(ncol(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh)))) {
  TheColName <- colnames(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh)[i]
  SpecCounts <- DF_PavianSpecies_Conf_0.6_Norm_AboveThresh[[TheColName]][which(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  SpecNames <- row.names(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh)[which(DF_PavianSpecies_Conf_0.6_Norm_AboveThresh[[TheColName]]>ThreshSpec)]
  names(SpecCounts) <- SpecNames
  ListSpecAboveThresh_0.6[[TheColName]] <- SpecCounts
} # Check thoroughly

ListSpecAboveThresh_0.6

unlist(sapply(ListSpecAboveThresh_0.6, function(x) names(x)))
unique(unlist(sapply(ListSpecAboveThresh_0.6, function(x) names(x))))
sort(unique(unlist(sapply(ListSpecAboveThresh_0.6, function(x) names(x)))))

DF_ListSpecAboveThresh_0.6 <- sort(unique(unlist(sapply(ListSpecAboveThresh_0.6, function(x) names(x)))))
DF_ListSpecAboveThresh_0.6 <- data.frame(ListSpecies = DF_ListSpecAboveThresh_0.6)
DF_ListSpecAboveThresh_0.6$MetaG_0.6 <- "Yes"

# Merge the three dataframes and export it

DF_ListSpecAboveThresh_AllConf <- list(DF_ListSpecAboveThresh_0.0, DF_ListSpecAboveThresh_0.1, DF_ListSpecAboveThresh_0.6) %>% reduce(full_join) # To be checked
write.csv(DF_ListSpecAboveThresh_AllConf, "/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaT/SpeciesAboveThreshold/ListSpecAboveThresh_AllConf_MetaG.csv", row.names = FALSE)

# Check relative abundance of the species that passed the threshold in metatranscriptomics

DF_PavianSpecies_Conf_0.1_Norm_Thresh_SpeciesMetaT <- DF_PavianSpecies_Conf_0.1_Norm_Thresh[row.names(DF_PavianSpecies_Conf_0.1_Norm_Thresh) %in% c("Mortierella globalpina", "Alternaria alternata", "Aureobasidium pullulans", "Malassezia restricta", "Malassezia globosa"),]

DF_PavianSpecies_Conf_0.1_Norm_Thresh_SpeciesMetaT_t <- as.data.frame(t(DF_PavianSpecies_Conf_0.1_Norm_Thresh_SpeciesMetaT))
sapply(DF_PavianSpecies_Conf_0.1_Norm_Thresh_SpeciesMetaT_t, function(x) max(x, na.rm = TRUE))

#====================================================================#
# Statistics to check if the communities are different between sites #
#====================================================================#

# 0.1 
#----

# Normalized to all reads
#------------------------

# Prepare the dataframe

DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F
DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat_WithGroup <- DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat
DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat_WithGroup$Group <- substr(row.names(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat_WithGroup), 3, 5)

ANOSIM_0.1_species_PercFungiAllReads <- anosim(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat, grouping = DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat_WithGroup$Group, permutations = 999, distance = "bray", strata = NULL)

Adonis_0.1_species_PercFungiAllReads <- adonis2(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat ~ DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat_WithGroup$Group, permutations = 999, by = "terms", distance = "bray", strata = NULL)

pairwise.adonis(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat, DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat_WithGroup$Group, perm = 999, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "bonferroni")

Multipatt_0.1_species_PercFungiAllReads <- multipatt(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat, cluster = DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat_WithGroup$Group, func = "IndVal.g", duleg = FALSE, restcomb = NULL, min.order = 1, max.order = NULL, control = how(), permutations = NULL, print.perm = FALSE)

# Check for beta dispersion (which may be the cause of the significance of the difference between groups)

DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat_DIST <- vegdist(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat, method = "bray")

Betadisper_0.1_species_PercFungiAllReads <- betadisper(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat_DIST, group = DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F_ForStat_WithGroup$Group) # /!\ I am not sure the groups are in the right order (they should be)

permutest(Betadisper_0.1_species_PercFungiAllReads) # Dispersion is significant, but not sure this is the good test
permutest(Betadisper_0.1_species_PercFungiAllReads, pairwise = TRUE) # Dispersion is not significat for most group, not sure this is the good test

# Normalized to the number of reads assigned to Fungi
#----------------------------------------------------

# Prepare the dataframe

DF_0.1_ForNormSpecNew_20240628_ForStat <- DF_0.1_ForNormSpecNew_20240628
DF_0.1_ForNormSpecNew_20240628_ForStat_WithGroup <- DF_0.1_ForNormSpecNew_20240628_ForStat
DF_0.1_ForNormSpecNew_20240628_ForStat_WithGroup$Group <- substr(row.names(DF_0.1_ForNormSpecNew_20240628_ForStat_WithGroup), 3, 5)

ANOSIM_0.1_species_PercFungiAllFungi <- anosim(DF_0.1_ForNormSpecNew_20240628_ForStat, grouping = DF_0.1_ForNormSpecNew_20240628_ForStat_WithGroup$Group, permutations = 999, distance = "bray", strata = NULL)

Adonis_0.1_species_PercFungiAllFungi <- adonis2(DF_0.1_ForNormSpecNew_20240628_ForStat ~ DF_0.1_ForNormSpecNew_20240628_ForStat_WithGroup$Group, permutations = 999, by = "terms", distance = "bray", strata = NULL)

pairwise.adonis(DF_0.1_ForNormSpecNew_20240628_ForStat, DF_0.1_ForNormSpecNew_20240628_ForStat_WithGroup$Group, perm = 999, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "bonferroni")

# Check for beta dispersion (which may be the cause of the significance of the difference between groups)

DF_0.1_ForNormSpecNew_20240628_ForStat_DIST <- vegdist(DF_0.1_ForNormSpecNew_20240628_ForStat, method = "bray")

Betadisper_0.1_species_PercFungiAllFungi <- betadisper(DF_0.1_ForNormSpecNew_20240628_ForStat_DIST, group = DF_0.1_ForNormSpecNew_20240628_ForStat_WithGroup$Group) # /!\ I am not sure the groups are in the right order (they should be)

permutest(Betadisper_0.1_species_PercFungiAllFungi) # Dispersion is not significant, but not sure this is the good test
permutest(Betadisper_0.1_species_PercFungiAllFungi, pairwise = TRUE) # Dispersion is not significant for most group, not sure this is the good test


#======================================================================================#
# Calculate the number of species, and the number of species above a certain threshold #
#======================================================================================#

# 0.1
#----

# All species including those that have a relative abundance of 0 (because they were in samples that were removed)

length(colSums(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F)) # 898

# All species that are present in at lease one sample of the good samples

length(colSums(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F)[which(colSums(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F) > 0)]) # 795

# All species that are above a certain threshold

length(colSums(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F)[which(sapply(DF_PavianSpecies_Conf_0.1_Norm_ForDivAlpha_GoodCol_t_F, max) > 0.000001)])

