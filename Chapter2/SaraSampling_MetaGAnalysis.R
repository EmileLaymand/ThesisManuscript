# 2024/06/03: I  have removed five samples as they were resequenced by Fasteris as they were considered of a too low quality.
# Check this does not disturb the code (especially when using index numbers in dataframes).

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
# Plot percent of reads assigned according to environmental parameters #
#======================================================================#

# # Remove samples from Chile
# 
# DF_MergedMetadataPavian_F <- DF_MergedMetadataPavian[DF_MergedMetadataPavian[["Country"]] == "F",]

# # Remove rows with NA /!\ This will need to be modified when the data arrives
# 
# DF_MergedMetadataPavian_F_NArm <- DF_MergedMetadataPavian_F[is.na(DF_MergedMetadataPavian_F[["Classified.reads"]]) == FALSE,]

# Order the factors in the site column so that the sites appear in the right order on the plot

DF_MergedMetadataPavian$siteID <- factor(DF_MergedMetadataPavian$siteID, levels = c("FMD", "FMS", "FSS", "FBB", "FHB", "FLC", "FGR", "FLP", "FCT", "FCF", "FCM", "FLM", "FLA", "CDI", "CPL", "CHI", "CBD", "CDE", "CHD", "CYH", "CCO", "CXH", "CBS", "CMH", "CHS", "CRA"))

# Remove rows with NA /!\ This will need to be modified when the data arrives

DF_MergedMetadataPavian_NArm <- DF_MergedMetadataPavian[is.na(DF_MergedMetadataPavian[["Classified.reads"]]) == FALSE,]

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm, aes(x = sal.sd, y = Classified.reads)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point() +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Classified.reads_sal.sd_FranceAndChile.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Classified.reads_sal.sd_FranceAndChile.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm, aes(x = sal.sd, y = Fungal.reads)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point() +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_sal.sd_FranceAndChile.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_sal.sd_FranceAndChile.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm, aes(x = siteID, y = Fungal.reads)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_boxplot() +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_siteID_FranceAndChile.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_siteID_FranceAndChile.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

# Keep only samples from France

DF_MergedMetadataPavian_NArm_F <- DF_MergedMetadataPavian_NArm[DF_MergedMetadataPavian_NArm[["Country"]] == "F",]
DF_MergedMetadataPavian_NArm_F$Month <- month(DF_MergedMetadataPavian_NArm_F$date)
DF_MergedMetadataPavian_NArm_F$MonthChar <- as.character(DF_MergedMetadataPavian_NArm_F$Month)
DF_MergedMetadataPavian_NArm_F$MonthChar <- factor(DF_MergedMetadataPavian_NArm_F$MonthChar, levels = c("1", "2", "3", "4", "5" ,"6", "7", "8", "9", "10", "11", "12"))

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm_F, aes(x = siteID, y = Fungal.reads)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point(aes(color = MonthChar), position=position_dodge(width=0.3)) +
  scale_colour_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_siteID_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_siteID_France.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm_F, aes(x = salinity, y = Fungal.reads)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point(aes(color = MonthChar)) +
  scale_colour_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_salinity_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_salinity_France.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm_F, aes(x = humic.like.DOM, y = Fungal.reads)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point(aes(color = MonthChar)) +
  scale_colour_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_humic.like.DOM_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_humic.like.DOM_France.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm_F, aes(x = protein.like.DOM, y = Fungal.reads)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point(aes(color = MonthChar)) +
  scale_colour_manual(values = rainbow(12)) +
  theme_bw()

p1 

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_protein.like.DOM_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/CLassifiedReadsEnvParameters/Fig_Fungal.reads_protein.like.DOM_France.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

# Check the number of missing values (divide by 3 as all rows are duplicated three times as there are 3 confidence levels)

table(DF_MergedMetadataPavian_NArm_F[is.na(DF_MergedMetadataPavian_NArm_F$DOC.mean),"siteID"])/3
table(DF_MergedMetadataPavian_NArm_F[is.na(DF_MergedMetadataPavian_NArm_F$chla.mean),"siteID"])/3


# Get the salinity standard deviation per site, and sort it

DF_MergedMetadataPavian_NArm_F_Sal <- DF_MergedMetadataPavian_NArm_F[is.na(DF_MergedMetadataPavian_NArm_F$salinity) == FALSE,]
DF_salinity <- aggregate(DF_MergedMetadataPavian_NArm_F_Sal$salinity, list(DF_MergedMetadataPavian_NArm_F_Sal$siteID), FUN=sd) 
Salinity_vect <- DF_salinity$x
names(Salinity_vect) <- DF_salinity$Group.1
sort(Salinity_vect, decreasing = TRUE)

# Plot environmental parameters per site

DF_MergedMetadataPavian_NArm_F_Conf0.0 <- DF_MergedMetadataPavian_NArm_F[DF_MergedMetadataPavian_NArm_F$Confidence == "0.0",]

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm_F_Conf0.0, aes(x = siteID, y = chla.mean)) +
  geom_point(aes(fill = MonthChar), shape = 21, colour = "black", position=position_dodge(width=0.3)) +
  scale_fill_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_Chla_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_Chla_France.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm_F_Conf0.0, aes(x = siteID, y = humic.like.DOM)) +
  geom_point(aes(fill = MonthChar), shape = 21, colour = "black", position=position_dodge(width=0.3)) +
  scale_fill_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_humic.like.DOM_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_humic.like.DOM_France.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm_F_Conf0.0, aes(x = siteID, y = protein.like.DOM)) +
  geom_point(aes(fill = MonthChar), shape = 21, colour = "black", position=position_dodge(width=0.3)) +
  scale_fill_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_protein.like.DOM_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_protein.like.DOM_France.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm_F_Conf0.0, aes(x = siteID, y = BP)) +
  geom_point(aes(fill = MonthChar), shape = 21, colour = "black", position=position_dodge(width=0.3)) +
  scale_fill_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_BP_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_BP_France.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm_F_Conf0.0, aes(x = siteID, y = RES)) +
  geom_point(aes(fill = MonthChar), shape = 21, colour = "black", position=position_dodge(width=0.3)) +
  scale_fill_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_RES_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_RES_France.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm_F_Conf0.0, aes(x = siteID, y = BGE)) +
  geom_point(aes(fill = MonthChar), shape = 21, colour = "black", position=position_dodge(width=0.3)) +
  scale_fill_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_BGE_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/FigEnvVarPerSite_BGE_France.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

p1 <- ggplot(data=DF_MergedMetadataPavian_NArm_F_Conf0.0, aes(x = humic.like.DOM, y = chla.mean)) +
  geom_point(aes(fill = MonthChar), shape = 21, colour = "black") +
  scale_fill_manual(values = rainbow(12)) +
  theme_bw()

p1

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/Fig_humic.like.DOM_chla.mean_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/EnvironmentalVariablesPerSite/Fig_humic.like.DOM_chla.mean_France.png", plot=p1, device = "png", dpi = 300, width=21, height=10)
dev.off()

#=====================================================================================#
# Make plots of classified reads and fungal reads as a function of the date per group #
#=====================================================================================#

cc <- c(palette("Set3"), palette("Dark2"), palette("Accent"), palette("Set1"))
#CustomPalette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "darkorange1", "black")
CustomPalette <- c("black", "#DF536B", "#61D04F", "#2297E6", "#28E2E5", "#CD0BBC", "#F5C710", "gray62", "#8DD3C7",  "#FB8072", "#BEBADA")

DF_MergedMetadataPavian_NArm_F$Date_num <- as.Date(DF_MergedMetadataPavian_NArm_F$date, "%Y-%m-%d")
p1 <- ggplot(data = DF_MergedMetadataPavian_NArm_F, aes(x = Date_num, y = Classified.reads, color = siteID)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = CustomPalette) +
  scale_x_date(date_breaks="1 month", date_labels="%b\n%Y") +
  labs(x = "Date", y = "Percent of classified reads") +
  theme_bw()

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/PercentClassifiedReads/FigClassifiedReads_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/PercentClassifiedReads/FigClassifiedReads_France.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

p1 <- ggplot(data = DF_MergedMetadataPavian_NArm_F, aes(x = Date_num, y = Fungal.reads, color = siteID)) +
  facet_grid(rows = vars(Confidence), scales = "free_y") +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = CustomPalette) +
  scale_x_date(date_breaks="1 month", date_labels="%b\n%Y") +
  labs(x = "Date", y = "Percent of reads classified as Fungi") +
  theme_bw()

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/PercentClassifiedReads/FigFungalReads_France.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/RawFigures/PercentClassifiedReads/FigFungalReads_France.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

# Boxplot of the relative abundance of Fungi per site
#----------------------------------------------------

DF_MergedMetadataPavian_ForBoxplot <- DF_MergedMetadataPavian

# Select only the right sites

RightSites <- c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT")
DF_MergedMetadataPavian_ForBoxplotRightSites <- DF_MergedMetadataPavian_ForBoxplot[DF_MergedMetadataPavian_ForBoxplot[["siteID"]] %in% RightSites,]
DF_MergedMetadataPavian_ForBoxplotRightSites$siteID <- factor(DF_MergedMetadataPavian_ForBoxplotRightSites$siteID, levels = RightSites)

# Remove two sites with no metagenomes

DF_MergedMetadataPavian_ForBoxplotRightSites <- DF_MergedMetadataPavian_ForBoxplotRightSites[which(is.na(DF_MergedMetadataPavian_ForBoxplotRightSites$Name) == FALSE),]

# Added on 2024-06-03 to remove the badly sequenced metagenomes when sequencing was reperformed by Fasteris

DF_MergedMetadataPavian_ForBoxplotRightSites <- DF_MergedMetadataPavian_ForBoxplotRightSites[!(DF_MergedMetadataPavian_ForBoxplotRightSites[["SampleName"]] %in% c("G_CCO01.1_G1", "G_FCT10.1_G1", "G_FGR08.1_G1", "G_FGR09.1_G1", "G_FMS04.1_G1", "G_FSS01.1_G1")),] # G_FSS01.1_G1 added on 2024-07-08
# c("AGTU-31", "AGTU-61", "AGTU-12", "AGTU-6", "AGTU-59")

DF_MergedMetadataPavian_ForBoxplotRightSites <- DF_MergedMetadataPavian_ForBoxplotRightSites[DF_MergedMetadataPavian_ForBoxplotRightSites[["siteID"]] %in% c("FMS", "FSS", "FHB", "FLC", "FGR", "FLP", "FCT"),]

# End of added on 2024-06-03

# Add color according to site

DF_MergedMetadataPavian_ForBoxplotRightSites$Coastality <- NA
DF_MergedMetadataPavian_ForBoxplotRightSites[DF_MergedMetadataPavian_ForBoxplotRightSites$siteID == "FMS","Coastality"]  <- "Offshore"
DF_MergedMetadataPavian_ForBoxplotRightSites[DF_MergedMetadataPavian_ForBoxplotRightSites$siteID %in% c("FHB", "FSS"),"Coastality"]  <- "Coastal"
DF_MergedMetadataPavian_ForBoxplotRightSites[DF_MergedMetadataPavian_ForBoxplotRightSites$siteID %in% c("FLC", "FGR", "FCT", "FLP"),"Coastality"]  <- "Lagoon"
DF_MergedMetadataPavian_ForBoxplotRightSites$Coastality <- factor(DF_MergedMetadataPavian_ForBoxplotRightSites$Coastality, levels = c("Offshore", "Coastal", "Lagoon"))

# Make the plot

p1 <- ggplot(data = DF_MergedMetadataPavian_ForBoxplotRightSites, aes(x = siteID, y = Classified.reads)) +
  facet_grid(cols = vars(Confidence)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(aes(fill = Coastality), outlier.colour = "black", outlier.shape = 1, outlier.size = 6) +
  #scale_fill_manual(values = c("deepskyblue", "darkgoldenrod1", "darkolivegreen1")) +
  scale_fill_manual(values = c("darkorchid2", "darkgoldenrod1", "darkolivegreen1")) +
  ylab(label = "Percent of reads classified by Kraken2") +
  xlab(label = "Site") +
  theme_bw() + 
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    strip.text = element_text(size = 20)
  )

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/AllClassifiedReads_7sites_AllConf.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/AllClassifiedReads_7sites_AllConf.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

p1 <- ggplot(data = DF_MergedMetadataPavian_ForBoxplotRightSites, aes(x = siteID, y = Fungal.reads)) +
    facet_grid(cols = vars(Confidence)) +
    stat_boxplot(geom ='errorbar', width = 0.5) +
    geom_boxplot(aes(fill = Coastality)) +
    #scale_fill_manual(values = c("deepskyblue", "darkgoldenrod1", "darkolivegreen1")) +
    scale_fill_manual(values = c("darkorchid2", "darkgoldenrod1", "darkolivegreen1")) +
    ylab(label = "Percent of reads classified as Fungi by Kraken2") +
    xlab(label = "Site") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20)
    )
  
print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/FungiClassifiedReads_7sites_AllConf_Column.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/FungiClassifiedReads_7sites_AllConf_Column.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

p1 <- ggplot(data = DF_MergedMetadataPavian_ForBoxplotRightSites, aes(x = siteID, y = Fungal.reads)) +
    facet_grid(rows = vars(Confidence), scales = "free") +
    stat_boxplot(geom ='errorbar', width = 0.5) +
    geom_boxplot(aes(fill = Coastality), outlier.colour = "black", outlier.shape = 1, outlier.size = 6) +
    #scale_fill_manual(values = c("deepskyblue", "darkgoldenrod1", "darkolivegreen1")) +
    scale_fill_manual(values = c("darkorchid2", "darkgoldenrod1", "darkolivegreen1")) +
    ylab(label = "Percent of reads classified as Fungi by Kraken2") +
    xlab(label = "Site") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      strip.text = element_text(size = 20)
    )
  
print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/FungiClassifiedReads_7sites_AllConf_Rows.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/SupplementaryFigures/FungiClassifiedReads_7sites_AllConf_Rows.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

# Select only samples with Confidence 0.1

DF_MergedMetadataPavian_ForBoxplotRightSites_0.1 <- DF_MergedMetadataPavian_ForBoxplotRightSites[DF_MergedMetadataPavian_ForBoxplotRightSites[["Confidence"]] == 0.1,]

p1 <- ggplot(data = DF_MergedMetadataPavian_ForBoxplotRightSites_0.1, aes(x = siteID, y = Classified.reads)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(aes(fill = Coastality), outlier.colour = "black", outlier.shape = 1, outlier.size = 6) +
  #scale_fill_manual(values = c("deepskyblue", "darkgoldenrod1", "darkolivegreen1")) +
  scale_fill_manual(values = c("darkorchid2", "darkgoldenrod1", "darkolivegreen1")) +
  ylab(label = "Percent of reads classified by Kraken2") +
  xlab(label = "Site") +
  theme_bw() + 
  theme(
    axis.text = element_text(size = 40),
    axis.title = element_text(size = 40),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.position="none"
  )

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/AllClassifiedReads_7sites_0.1.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/AllClassifiedReads_7sites_0.1.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

p1 <- ggplot(data = DF_MergedMetadataPavian_ForBoxplotRightSites_0.1, aes(x = siteID, y = Fungal.reads)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(aes(fill = Coastality), outlier.colour = "black", outlier.shape = 1, outlier.size = 6) +
  #scale_fill_manual(values = c("deepskyblue", "darkgoldenrod1", "darkolivegreen1")) +
  scale_fill_manual(values = c("darkorchid2", "darkgoldenrod1", "darkolivegreen1")) +
  ylab(label = "Percent of reads classified as Fungi by Kraken2") +
  xlab(label = "Site") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 40),
    axis.title = element_text(size = 40),
    legend.text = element_text(size = 40),
    legend.title = element_text(size = 40),
    legend.position="none"
  )

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/FungiClassifiedReads_7sites_0.1.pdf", plot=p1, device = cairo_pdf(), width=21, height=20)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre2/MainFigures/FungiClassifiedReads_7sites_0.1.png", plot=p1, device = "png", dpi = 300, width=21, height=20)
dev.off()

# Take only confidence 0.0 and 0.6 for short calculations

DF_MergedMetadataPavian_ForBoxplotRightSites_0.0 <- DF_MergedMetadataPavian_ForBoxplotRightSites[DF_MergedMetadataPavian_ForBoxplotRightSites[["Confidence"]] == "0.0",]
DF_MergedMetadataPavian_ForBoxplotRightSites_0.6 <- DF_MergedMetadataPavian_ForBoxplotRightSites[DF_MergedMetadataPavian_ForBoxplotRightSites[["Confidence"]] == 0.6,]

#======================================================================================================#
# Added on 2024-09-25. Calculate the proportion of each domain (Bacteria, Eukaryota, Archaea, Viruses) #
#======================================================================================================#

DF_MetaG_Domains_Conf_0.1 <- read.table("/home/emilelaymand/Documents/Science/These/Sara_Sampling_2022/metaG/SortiesKraken2/Pavian_Confidence_0.1/Analysis/RawPavianReports/SaraSampling_MetaG_Kraken2_Confidence_0.1_Domains.tsv", sep = "\t", header = TRUE)

# Rename colomns to keep only the AGTU index

colnames(DF_MetaG_Domains_Conf_0.1)[5:114] <- gsub("\\.", "-", substr(colnames(DF_MetaG_Domains_Conf_0.1)[5:114], 23, nchar(colnames(DF_MetaG_Domains_Conf_0.1)[5:114])-11))

# Remove useless columns

DF_MetaG_Domains_Conf_0.1$lineage <- NULL
DF_MetaG_Domains_Conf_0.1$taxRank <- NULL
DF_MetaG_Domains_Conf_0.1$taxID <- NULL
DF_MetaG_Domains_Conf_0.1$Max <- NULL

# Rename the rows with the domains

row.names(DF_MetaG_Domains_Conf_0.1) <- DF_MetaG_Domains_Conf_0.1$name
DF_MetaG_Domains_Conf_0.1$name <- NULL

# Transpose the dataframe

DF_MetaG_Domains_Conf_0.1_t <- data.frame(t(DF_MetaG_Domains_Conf_0.1))

# Add a column that is the AGTU index of the sample

DF_MetaG_Domains_Conf_0.1_t$AGTUIndex <- row.names(DF_MetaG_Domains_Conf_0.1_t)

# Merge both dataframes, keeping only the good samples

DF_MetaG_Domains_Conf_0.1_Merged <- join(DF_MergedMetadataPavian_ForBoxplotRightSites_0.1, DF_MetaG_Domains_Conf_0.1_t, type="left") # /!\ Hard check (really check)

# Calculate the percent of all reads assigned to each domain

DF_MetaG_Domains_Conf_0.1_Merged$BacteriaPerc <- (DF_MetaG_Domains_Conf_0.1_Merged$Bacteria/DF_MetaG_Domains_Conf_0.1_Merged$Number.of.raw.reads)*100
DF_MetaG_Domains_Conf_0.1_Merged$EukaryotaPerc <- (DF_MetaG_Domains_Conf_0.1_Merged$Eukaryota/DF_MetaG_Domains_Conf_0.1_Merged$Number.of.raw.reads)*100
DF_MetaG_Domains_Conf_0.1_Merged$ArchaeaPerc <- (DF_MetaG_Domains_Conf_0.1_Merged$Archaea/DF_MetaG_Domains_Conf_0.1_Merged$Number.of.raw.reads)*100
DF_MetaG_Domains_Conf_0.1_Merged$VirusesPerc <- (DF_MetaG_Domains_Conf_0.1_Merged$Viruses/DF_MetaG_Domains_Conf_0.1_Merged$Number.of.raw.reads)*100
