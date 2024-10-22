#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  GENERAL CLEANING AND SETTING  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Cleans the environment
rm(list = ls())
dev.off()
cat("\014")

## Set the working directory
#setwd("/home/emilelaymand/Documents/Science/Master2_MES/StageM2")

require("ggplot2")
require("phyloseq")
library(microViz)
library(stringr)
library(tidyverse)
library(gridExtra)
library(microbiomeSeq)
library(hues)
library(microbiome)

source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/Generics.R")
source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/Pathnames.R")
source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/phyloseqExtended.R")

#=========================#
#  Convenient functions   #
#=========================#

plot_bar_2 <-  function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

#================================================#
# 1st section: Eukaryotes 0.2 micrometers import #
#================================================#

#--------------------------------------------------------
# Step 1 : Import the data and build the phyloseq object
#--------------------------------------------------------

# Import the file containing sequences, OTU abundances per sample and taxonomy for each OTU

Euka02_OTU_taxo_Seq <- read.csv("SOLA/wetransfer_sola_2021-10-18_1242/2013-2017_OTU_DADA2_Euk_02_modified.csv")

# Create a dummy OTU name for each sequence

Euka02_OTU_taxo_Seq$ASV <- paste("ASV", c(1:length(row.names(Euka02_OTU_taxo_Seq))), sep = "")

# Separate the sequences, the taxonomy, and the abundances in 3 different dataframes. Keep column "ASV" in all dataframes so that you have correspondences between dataframes

## Sequences

Euka02_Seq <- Euka02_OTU_taxo_Seq[ ,c("ASV", "X")]

## Taxonomy

Euka02_taxo <- Euka02_OTU_taxo_Seq[, c("ASV", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

## OTU abundances (all columns but previous)

Euka02_OTU <- Euka02_OTU_taxo_Seq[, !(colnames(Euka02_OTU_taxo_Seq) %in% c( "X", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"))]

# Make the ASV column of all these dataframes the row names (necessary to build the Phyloseq object)

row.names(Euka02_Seq) <- Euka02_Seq$ASV

row.names(Euka02_taxo) <- Euka02_taxo$ASV
Euka02_taxo$ASV <- NULL
Euka02_taxo_mat <- as.matrix(Euka02_taxo)

row.names(Euka02_OTU) <- Euka02_OTU$ASV
Euka02_OTU$ASV <- NULL

# Import the file containing metadata (sampling date, chlorophyll concentration, depth, GPS coordinates, ...)

Euka02_Metadata <- read.csv("SOLA/wetransfer_sola_2021-10-18_1242/Equivalence_name_date_3um_modified.csv")

# Crate the right sample names, then make the sample names the row names (necessary to build the phyloseq object)
Euka02_Metadata$Name_Euk02 <- paste(str_sub(Euka02_Metadata$Name_Pro, 1, -4), "Euk", sep = "")
Euka02_Metadata <- Euka02_Metadata[Euka02_Metadata[["Name_Euk02"]] %in% colnames(Euka02_OTU),]

row.names(Euka02_Metadata) <- Euka02_Metadata$Name_Euk02

# Merge the data into a phyloseq object

OTU_02 <- otu_table(Euka02_OTU, taxa_are_rows = TRUE)

TAX_02 <- tax_table(Euka02_taxo_mat)

META_02 <- sample_data(Euka02_Metadata)

Phyloseq_Euk02 <- phyloseq(OTU_02, TAX_02, META_02)

#==============================================#
# 2nd section: Eukaryotes 3 micrometers import #
#==============================================#

#--------------------------------------------------------
# Step 1 : Import the data and build the phyloseq object
#--------------------------------------------------------

# Import the file containing sequences, OTU abundances per sample and taxonomy for each OTU

Euka3_OTU_taxo_Seq <- read.csv("SOLA/wetransfer_sola_2021-10-18_1242/2013-2017_OTU_DADA2_Euk_3um.csv")

# Create a dummy OTU name for each sequence

Euka3_OTU_taxo_Seq$ASV <- paste("ASV", c(1:length(row.names(Euka3_OTU_taxo_Seq))), sep = "")

# Separate the sequences, the taxonomy, and the abundances in 3 different dataframes. Keep column "ASV" in all dataframes so that you have correspondences between dataframes

## Sequences

Euka3_Seq <- Euka3_OTU_taxo_Seq[ ,c("ASV", "X")]

## Taxonomy

Euka3_taxo <- Euka3_OTU_taxo_Seq[, c("ASV", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

## OTU abundances (all columns but previous)

Euka3_OTU <- Euka3_OTU_taxo_Seq[, !(colnames(Euka3_OTU_taxo_Seq) %in% c( "X", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"))]

# Make the ASV column of all these dataframes the row names (necessary to build the Phyloseq object)

row.names(Euka3_Seq) <- Euka3_Seq$ASV

row.names(Euka3_taxo) <- Euka3_taxo$ASV
Euka3_taxo$ASV <- NULL
Euka3_taxo_mat <- as.matrix(Euka3_taxo)

row.names(Euka3_OTU) <- Euka3_OTU$ASV
Euka3_OTU$ASV <- NULL

# Import the file containing metadata (sampling date, chlorophyll concentration, depth, GPS coordinates, ...)

Euka3_Metadata <- read.csv("SOLA/wetransfer_sola_2021-10-18_1242/Equivalence_name_date_3um_modified.csv")

# Make the sample names the row names (necessary to build the phyloseq object)

row.names(Euka3_Metadata) <- Euka3_Metadata$Name_Euk

# Merge the data into a phyloseq object

OTU_3 <- otu_table(Euka3_OTU, taxa_are_rows = TRUE)

TAX_3 <- tax_table(Euka3_taxo_mat)

META_3 <- sample_data(Euka3_Metadata)

Phyloseq_Euk3 <- phyloseq(OTU_3, TAX_3, META_3)

########################################################################################
# CORRELATION BETWEEN THE NUMBER OF READS IN A SAMPLE AND THE FUNGI RELATIVE ABUNDANCE #
########################################################################################

# Number of reads per sample
#---------------------------

sort(sample_sums(Phyloseq_Euk02)) ## SMALL SIZE FRACTION
ReadCount02Sorted <- sample_sums(Phyloseq_Euk02)[sort(names(sample_sums(Phyloseq_Euk02)))]

sort(sample_sums(Phyloseq_Euk3)) ## BIG SIZE FRACTION
ReadCount3Sorted <- sample_sums(Phyloseq_Euk3)[sort(names(sample_sums(Phyloseq_Euk3)))]


# Relative abundance of Fungi
#----------------------------

# Extract Fungi relative abundance and sort it
Phyloseq_Euk02_perc <- taxa_percentize(Phyloseq_Euk02, TaxLevel = "Division")
ASVnameFungi02 <- row.names(tax_table(Phyloseq_Euk02_perc)[tax_table(Phyloseq_Euk02_perc)[,"Division"] %in% c("Fungi"),])
RelAbFungi02 <- otu_table(Phyloseq_Euk02_perc)[ASVnameFungi02,]
RelAbFungi02V <- as.vector(RelAbFungi02)
names(RelAbFungi02V) <- colnames(RelAbFungi02)
RelAbFungi02Sorted <- RelAbFungi02V[sort(names(RelAbFungi02V))]

Phyloseq_Euk3_perc <- taxa_percentize(Phyloseq_Euk3, TaxLevel = "Division")
ASVnameFungi3 <- row.names(tax_table(Phyloseq_Euk3_perc)[tax_table(Phyloseq_Euk3_perc)[,"Division"] %in% c("Fungi"),])
RelAbFungi3 <- otu_table(Phyloseq_Euk3_perc)[ASVnameFungi3,]
RelAbFungi3V <- as.vector(RelAbFungi3)
names(RelAbFungi3V) <- colnames(RelAbFungi3)
RelAbFungi3Sorted <- RelAbFungi3V[sort(names(RelAbFungi3V))]

# Check if the sample order is the same between read counts and fungal relative abundance
identical(names(ReadCount02Sorted), names(RelAbFungi02Sorted)) ## SMALL SIZE FRACTION

identical(names(ReadCount3Sorted), names(RelAbFungi3Sorted))

# Plot read counts again fungal relative abundance
plot(RelAbFungi02Sorted ~ ReadCount02Sorted)

plot(RelAbFungi3Sorted ~ ReadCount3Sorted)

##### CONCLUSION #####
# Some samples with low read counts (especially the one sample with the lowest read count) have a high fungal relative abundance, but overall the number of 
# reads and the relative abundance of Fungi is not correlated.

# Check if Ascomycota dominate the low counts samples 
#----------------------------------------------------

cc <- c(palette("Set3"), palette("Dark2"), palette("Accent"), palette("Set1"))

FungiPhyloseq_Euk02 <- subset_taxa(Phyloseq_Euk02, Division == "Fungi")
FungiPhyloseq_Euk02 <- prune_taxa(taxa_sums(FungiPhyloseq_Euk02) > 0, FungiPhyloseq_Euk02)
FungiPhyloseq_Euk02_perc <- taxa_percentize(FungiPhyloseq_Euk02, TaxLevel = "Class")
plot_bar(FungiPhyloseq_Euk02_perc, x="Name_Euk", fill="Class") + scale_fill_manual(values = cc) # Have an idea of present Kingdoms
# -> 50 % Basidio, 2 % Asco, 48 % Chytridio

FungiPhyloseq_Euk3 <- subset_taxa(Phyloseq_Euk3, Division == "Fungi")
FungiPhyloseq_Euk3 <- prune_taxa(taxa_sums(FungiPhyloseq_Euk3) > 0, FungiPhyloseq_Euk3)
FungiPhyloseq_Euk3_perc <- taxa_percentize(FungiPhyloseq_Euk3, TaxLevel = "Class")
plot_bar(FungiPhyloseq_Euk3_perc, x="Name_Euk", fill="Class") + scale_fill_manual(values = cc) # Have an idea of present Kingdoms
# -> 3/4 Basidio, 1/4 Asco