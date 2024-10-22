# /!\ I added a prune sample to remove sample PF256Euk from data.
# I must check this does not disturb the code.

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
library(plotly)
library(ggfortify)
library(FactoMineR)
library("factoextra")
library("corrplot")
library(missMDA)
library(reshape2)
library(multcomp)
library(FSA)
library(vegan)
library(rstatix)

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

#======================================================#
# 3rd section: Eukaryotes 0.2 and 3 micrometers import #
#======================================================#

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The idea with this section is to merge data from 0.2 and 3 micrometers data so that the ASV names are the same
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Import the file containing sequences, OTU abundances per sample and taxonomy for each OTU

Euka02_OTU_taxo_Seq_2merge <- read.csv("/home/emilelaymand/Documents/Science/These/SOLA_Galand/SOLA/wetransfer_sola_2021-10-18_1242/2013-2017_OTU_DADA2_Euk_02_modified.csv")

Euka3_OTU_taxo_Seq_2merge <- read.csv("/home/emilelaymand/Documents/Science/These/SOLA_Galand/SOLA/wetransfer_sola_2021-10-18_1242/2013-2017_OTU_DADA2_Euk_3um.csv")

# Remove the taxonomy columns to avoid duplicates of taxa that do not have exactly the same taxonomy

Euka02_OTU_taxo_Seq_2merge_wotax <- Euka02_OTU_taxo_Seq_2merge[, !(colnames(Euka02_OTU_taxo_Seq_2merge) %in% c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"))] 

Euka3_OTU_taxo_Seq_2merge_wotax <- Euka3_OTU_taxo_Seq_2merge[, !(colnames(Euka3_OTU_taxo_Seq_2merge) %in% c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"))] 

# Merge both dataframes

EukaAll_OTU_Seq <- list(Euka02_OTU_taxo_Seq_2merge_wotax, Euka3_OTU_taxo_Seq_2merge_wotax) %>% reduce(full_join) # Hard check

# ASVs that are in one filter size and not the other will create NAs in the latter. Replace all NAs by 0

EukaAll_OTU_Seq <- EukaAll_OTU_Seq %>% replace(is.na(.), 0)

# That that merging was good (no modification in taxonomy, not twice the same sequence with two different taxonomies)

length(unique(Euka02_OTU_taxo_Seq_2merge[["X"]]))

length(Euka02_OTU_taxo_Seq_2merge[["X"]])

length(unique(Euka3_OTU_taxo_Seq_2merge[["X"]]))
length(Euka3_OTU_taxo_Seq_2merge[["X"]])
sum((unique(Euka02_OTU_taxo_Seq_2merge[["X"]]) %in% unique(Euka3_OTU_taxo_Seq_2merge[["X"]])), na.rm = TRUE)
length(EukaAll_OTU_Seq[["X"]])
length(unique(EukaAll_OTU_Seq[["X"]]))

# For the taxa that would be duplicated keep only the taxonomy of the Eukaryotes 3 um (this is arbitrary)

## Take all sequences from 3 um

Euka3_taxo_seq_2merge <- Euka3_OTU_taxo_Seq_2merge[, c("X", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

## Take sequences that are in 2 um but not in 3 um

Euka02_taxo_seq_2merge <- Euka02_OTU_taxo_Seq_2merge[!(Euka02_OTU_taxo_Seq_2merge[["X"]] %in% Euka3_taxo_seq_2merge[["X"]]), c("X", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

## rbind these two dataframes to have the final taxonomy

EukaAll_taxo_Seq <- rbind(Euka3_taxo_seq_2merge, Euka02_taxo_seq_2merge)

## Check dimensions

length(unique(EukaAll_taxo_Seq[["X"]]))

# Merge the synthetic taxonomy and the OTU dataframes to make sure everything is in the right order

EukaAll_OTU_taxo_Seq <- list(EukaAll_OTU_Seq, EukaAll_taxo_Seq) %>% reduce(full_join) # Hard check

# Check dimensions

dim(EukaAll_OTU_taxo_Seq)

# Check for some sequences that the taxonomy is good

## 1st check

EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["X"]] == "AGCTCCAATAGCGTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCGGTTGAGAACGGCCGGTCCGCCGTTTGGTGTGCACTGGCTGGTTTCAACTTCCTGTAGAGGACGCGCTCTGGCTTCATCGCTGGACGCGGAGTCTACGTGGTTACTTTGAAAAAATTAGAGTGTTCAAAGCGGGCTTACGCTTGAATATTTCAGCATGGAATAACACTATAGGACTCCTGTCCTATTTCGTTGGTCTCGGGACGGGAGTAATGATTAAGAGGAACAGTTGGGGGCATTCGTATTTCATTGTCAGAGGTGAAATTCTTGGATTTATGAAAGACGAACTTCTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka02_OTU_taxo_Seq_2merge[Euka02_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCGGTTGAGAACGGCCGGTCCGCCGTTTGGTGTGCACTGGCTGGTTTCAACTTCCTGTAGAGGACGCGCTCTGGCTTCATCGCTGGACGCGGAGTCTACGTGGTTACTTTGAAAAAATTAGAGTGTTCAAAGCGGGCTTACGCTTGAATATTTCAGCATGGAATAACACTATAGGACTCCTGTCCTATTTCGTTGGTCTCGGGACGGGAGTAATGATTAAGAGGAACAGTTGGGGGCATTCGTATTTCATTGTCAGAGGTGAAATTCTTGGATTTATGAAAGACGAACTTCTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka3_OTU_taxo_Seq_2merge[Euka3_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCGGTTGAGAACGGCCGGTCCGCCGTTTGGTGTGCACTGGCTGGTTTCAACTTCCTGTAGAGGACGCGCTCTGGCTTCATCGCTGGACGCGGAGTCTACGTGGTTACTTTGAAAAAATTAGAGTGTTCAAAGCGGGCTTACGCTTGAATATTTCAGCATGGAATAACACTATAGGACTCCTGTCCTATTTCGTTGGTCTCGGGACGGGAGTAATGATTAAGAGGAACAGTTGGGGGCATTCGTATTTCATTGTCAGAGGTGAAATTCTTGGATTTATGAAAGACGAACTTCTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

## 2nd check

EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATTTCTGTTGAGGACGGCCGGTCCGCATTATGTGCGTGTATCTGGTTCGGCCTTGGCATCCTCCAGGGGAACGTTCCTGCGCTTCGCTGCGTGGGACGGTATTCTGGACTTTTACTTTGAGGAAATTAGAGTGTTCACGGCAGGCAGTCGCCTTGAATATATTAGCATGGAATAATAATATAGGACCTTGGTTCTATTTTGTTGGTTTCTAGAACTTAGGTAATGATTAATAGGGACAATTGGGGGCATTTGTATTAACACGTCAGAGGTGAAATTCTTGGATTGTGTTACGACAAACTACTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka02_OTU_taxo_Seq_2merge[Euka02_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATTTCTGTTGAGGACGGCCGGTCCGCATTATGTGCGTGTATCTGGTTCGGCCTTGGCATCCTCCAGGGGAACGTTCCTGCGCTTCGCTGCGTGGGACGGTATTCTGGACTTTTACTTTGAGGAAATTAGAGTGTTCACGGCAGGCAGTCGCCTTGAATATATTAGCATGGAATAATAATATAGGACCTTGGTTCTATTTTGTTGGTTTCTAGAACTTAGGTAATGATTAATAGGGACAATTGGGGGCATTTGTATTAACACGTCAGAGGTGAAATTCTTGGATTGTGTTACGACAAACTACTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka3_OTU_taxo_Seq_2merge[Euka3_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATTTCTGTTGAGGACGGCCGGTCCGCATTATGTGCGTGTATCTGGTTCGGCCTTGGCATCCTCCAGGGGAACGTTCCTGCGCTTCGCTGCGTGGGACGGTATTCTGGACTTTTACTTTGAGGAAATTAGAGTGTTCACGGCAGGCAGTCGCCTTGAATATATTAGCATGGAATAATAATATAGGACCTTGGTTCTATTTTGTTGGTTTCTAGAACTTAGGTAATGATTAATAGGGACAATTGGGGGCATTTGTATTAACACGTCAGAGGTGAAATTCTTGGATTGTGTTACGACAAACTACTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

# 3rd check in 02 but not in 3

EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGAGTTCTGCCTGGTGACGCCCGTCCGCCCAAGTGGTGTGTACAGGGTGTACATCTGGCCCTTTCAAGGGGAACGTATCTGCACTTAACTGTGCGGTGCGAGATCCTTGACTTTTACTTTGAGGAAATCAGAGTGTTCCAAGCAGGCTCTCGTCGTGCATGTTTCAGCATGGAATAATAGCATTGGACCTCGTCTCTCAGCTGTTGGTTGCAAGAAGCGAGGTAATGATGAAGAGGGATAGTTGGGGGCATTCGTATTTAACTGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGGACTACTGCGAAAGCATCTGCCATGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka02_OTU_taxo_Seq_2merge[Euka02_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGAGTTCTGCCTGGTGACGCCCGTCCGCCCAAGTGGTGTGTACAGGGTGTACATCTGGCCCTTTCAAGGGGAACGTATCTGCACTTAACTGTGCGGTGCGAGATCCTTGACTTTTACTTTGAGGAAATCAGAGTGTTCCAAGCAGGCTCTCGTCGTGCATGTTTCAGCATGGAATAATAGCATTGGACCTCGTCTCTCAGCTGTTGGTTGCAAGAAGCGAGGTAATGATGAAGAGGGATAGTTGGGGGCATTCGTATTTAACTGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGGACTACTGCGAAAGCATCTGCCATGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka3_OTU_taxo_Seq_2merge[Euka3_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGAGTTCTGCCTGGTGACGCCCGTCCGCCCAAGTGGTGTGTACAGGGTGTACATCTGGCCCTTTCAAGGGGAACGTATCTGCACTTAACTGTGCGGTGCGAGATCCTTGACTTTTACTTTGAGGAAATCAGAGTGTTCCAAGCAGGCTCTCGTCGTGCATGTTTCAGCATGGAATAATAGCATTGGACCTCGTCTCTCAGCTGTTGGTTGCAAGAAGCGAGGTAATGATGAAGAGGGATAGTTGGGGGCATTCGTATTTAACTGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGGACTACTGCGAAAGCATCTGCCATGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

# 4th check in both but different taxonomy

Test_tax <- list(Euka02_OTU_taxo_Seq_2merge, Euka3_OTU_taxo_Seq_2merge) %>% reduce(full_join) # Hard check 
Test_tax[duplicated(Test_tax[["X"]]),]

EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["X"]] == "AGCTCCAATAGTGTATATTAAAGTTGCTGCGGTTAAAATGCTCGTAGTTTAATTTCTGCTGAGGGTAATCGGTCCACCCACTGGGTGAGTAATTGTTTGCCCTTTGCATTTTGTGAACAATGCATCTGCACTTGATTGTGCTGTGTGTCCTGTTCATGCAATTTACTTTGAGGAAATTAGAGTGTTTCGAGCATGCATATGCACCGAGCACATTAGCATGGAATAATTTACACTGATCGTCGCTGTATTTGTTGGTTTTTAGGGCTTCGATAATGACTGATAGGGATAATTGGTGGTATTCGCATTTAATAGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGAACTGCTGCGAAAGCATTTGCCAACTATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka02_OTU_taxo_Seq_2merge[Euka02_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGTGTATATTAAAGTTGCTGCGGTTAAAATGCTCGTAGTTTAATTTCTGCTGAGGGTAATCGGTCCACCCACTGGGTGAGTAATTGTTTGCCCTTTGCATTTTGTGAACAATGCATCTGCACTTGATTGTGCTGTGTGTCCTGTTCATGCAATTTACTTTGAGGAAATTAGAGTGTTTCGAGCATGCATATGCACCGAGCACATTAGCATGGAATAATTTACACTGATCGTCGCTGTATTTGTTGGTTTTTAGGGCTTCGATAATGACTGATAGGGATAATTGGTGGTATTCGCATTTAATAGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGAACTGCTGCGAAAGCATTTGCCAACTATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka3_OTU_taxo_Seq_2merge[Euka3_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGTGTATATTAAAGTTGCTGCGGTTAAAATGCTCGTAGTTTAATTTCTGCTGAGGGTAATCGGTCCACCCACTGGGTGAGTAATTGTTTGCCCTTTGCATTTTGTGAACAATGCATCTGCACTTGATTGTGCTGTGTGTCCTGTTCATGCAATTTACTTTGAGGAAATTAGAGTGTTTCGAGCATGCATATGCACCGAGCACATTAGCATGGAATAATTTACACTGATCGTCGCTGTATTTGTTGGTTTTTAGGGCTTCGATAATGACTGATAGGGATAATTGGTGGTATTCGCATTTAATAGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGAACTGCTGCGAAAGCATTTGCCAACTATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

# Separate again the OTU, sequences and taxonomy into 3 different dataframes

# Create a dummy OTU name for each sequence

EukaAll_OTU_taxo_Seq$ASV <- paste("ASV", c(1:length(row.names(EukaAll_OTU_taxo_Seq))), sep = "")

# Separate the sequences, the taxonomy, and the abundances in 3 different dataframes. Keep column "ASV" in all dataframes so that you have correspondences between dataframes

## Sequences

EukaAll_Seq <- EukaAll_OTU_taxo_Seq[ ,c("ASV", "X")]

## Taxonomy

EukaAll_taxo <- EukaAll_OTU_taxo_Seq[, c("ASV", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

## OTU abundances (all columns but previous)

EukaAll_OTU <- EukaAll_OTU_taxo_Seq[, !(colnames(EukaAll_OTU_taxo_Seq) %in% c( "X", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"))]

# Make the ASV column of all these dataframes the row names (necessary to build the Phyloseq object)

row.names(EukaAll_Seq) <- EukaAll_Seq$ASV

row.names(EukaAll_taxo) <- EukaAll_taxo$ASV
EukaAll_taxo$ASV <- NULL
EukaAll_taxo_mat <- as.matrix(EukaAll_taxo)

row.names(EukaAll_OTU) <- EukaAll_OTU$ASV
EukaAll_OTU$ASV <- NULL

# Import the file containing metadata (sampling date, chlorophyll concentration, depth, GPS coordinates, ...)

EukaAll_Metadata <- read.csv("/home/emilelaymand/Documents/Science/These/SOLA_Galand/SOLA/wetransfer_sola_2021-10-18_1242/Equivalence_name_date_3um_concatenated.csv")

# Crate the right sample names, then make the sample names the row names (necessary to build the phyloseq object)

EukaAll_Metadata[str_sub(EukaAll_Metadata[["Name_Euk"]], -3, -1) == "Pro", "Name_Euk"] <- paste(str_sub(EukaAll_Metadata[str_sub(EukaAll_Metadata[["Name_Euk"]], -3, -1) == "Pro", "Name_Euk"], 1, -4), "Euk", sep = "")

row.names(EukaAll_Metadata) <- EukaAll_Metadata$Name_Euk

# Add one more column with a date that is actually considered as a date

EukaAll_Metadata$Date_num <- as.Date(EukaAll_Metadata$Date_Euk, "%d_%m_%Y")

# Import a second metadata dataframe, that contains details on nutrients, chl a, etc.

EukaAll_Metadata_nut <- read.csv("/home/emilelaymand/Documents/Science/These/SOLA_Galand/SOLA/wetransfer_sola_2021-10-18_1242/climato_2007-2017_modified.csv")

EukaAll_Metadata_nut$Date_num <- as.Date(EukaAll_Metadata_nut$DATE, "%d/%m/%Y")

# Merge both dataframes, then keep only rows that are present in the first dataframe 

EukaAll_Metadata_merged <- list(EukaAll_Metadata, EukaAll_Metadata_nut) %>% reduce(full_join) # Hard check 

EukaAll_Metadata_merged_red <- EukaAll_Metadata_merged[EukaAll_Metadata_merged[["Date_num"]] %in% unique(EukaAll_Metadata[["Date_num"]]), ]

# Check dimensions

dim(EukaAll_Metadata)

dim(EukaAll_Metadata_merged_red)

# Replace EukaAll_Metadata with the new, more complete EukaAll_Metadata_merged_red

EukaAll_Metadata <- EukaAll_Metadata_merged_red

row.names(EukaAll_Metadata) <- EukaAll_Metadata$Name_Euk

# Merge the data into a phyloseq object

OTU_All <- otu_table(EukaAll_OTU, taxa_are_rows = TRUE)

TAX_All <- tax_table(EukaAll_taxo_mat)

META_All <- sample_data(EukaAll_Metadata)

Phyloseq_EukAll <- phyloseq(OTU_All, TAX_All, META_All) 

# Check for ASV 500 (to compare to the raw .csv files)

Seq_ASV500 <- EukaAll_Seq[EukaAll_Seq[["ASV"]] == "ASV500", "X"]
tax_table(Phyloseq_EukAll)["ASV500",]
otu_table(Phyloseq_EukAll)["ASV500",]

# Check 2: present in both but different taxonomy

Seq_ASV346 <- EukaAll_Seq[EukaAll_Seq[["ASV"]] == "ASV346", "X"]
tax_table(Phyloseq_EukAll)["ASV346",]
otu_table(Phyloseq_EukAll)["ASV346",]

# Check 3: present in one dataset, but not the other

Seq_ASV1420 <- EukaAll_Seq[EukaAll_Seq[["ASV"]] == "ASV1420", "X"]
tax_table(Phyloseq_EukAll)["ASV1420",]
otu_table(Phyloseq_EukAll)["ASV1420",]

#===================================================================================================#
# End of the copy from the original MetaB_Sola.R file. Below is the new stuff from MetaB_SOLA_NMDS. #
#===================================================================================================#

sample_data(Phyloseq_EukAll)[["Month"]] <- substr(sample_data(Phyloseq_EukAll)[["Date_num"]], start = 6, stop = 7)
sample_data(Phyloseq_EukAll)[["Year"]] <- substr(sample_data(Phyloseq_EukAll)[["Date_num"]], start = 1, stop = 4)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# I remove the sample PF256Euk as it has a very low number of reads and a strange composition # # Added on 2024_02_23
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

Phyloseq_EukAll_woPF256Euk <- subset_samples(Phyloseq_EukAll, Name_Euk != "PF256Euk")
Phyloseq_EukAll <- Phyloseq_EukAll_woPF256Euk

# Define "outstanding samples" for easier spotting on the plots (outstanding = the relative abundance of FUngi is high)

OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)

# Here is the modification to select only Fungi in the dataset
OTUPhyloseq_NormDomain <- subset_taxa(OTUPhyloseq_NormDomain, Division == "Fungi")
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain) # Some samples seem to contain 0 Fungi, and it causes issues with converting to percent as it divides by 0
# End of modification

OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "OTU")
# Phyloseq_ForOutstanding <- tax_glom(OTUPhyloseq_NormDomain, "Division", NArm = FALSE)
# ASVnameFungi <- row.names(tax_table(Phyloseq_ForOutstanding)[tax_table(Phyloseq_ForOutstanding)[,"Division"] %in% c("Fungi"),])
# 
# Above1 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 1])
# Above5 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 5])
# Above10 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 10])
# Above15 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 15])
# Above20 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 20])
# Above30 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 30])
# Above40 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 40])
# Bet1And5 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 1 & otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] <= 5])
# Bet5And10 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 5  & otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] <= 10])
# Bet10And15 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 10  & otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] <= 15])
# Bet15And20 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 15  & otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] <= 20])
# Bet20And30 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 20  & otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] <= 30])
# Above30FungiRange <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 30])
# 
# sample_data(Phyloseq_EukAll)[["Above1"]] <- FALSE
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above1, "Above1"] <- TRUE
# 
# sample_data(Phyloseq_EukAll)[["Above5"]] <- FALSE
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above5, "Above5"] <- TRUE
# 
# sample_data(Phyloseq_EukAll)[["Above10"]] <- FALSE
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above10, "Above10"] <- TRUE
# 
# sample_data(Phyloseq_EukAll)[["Above15"]] <- FALSE
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above15, "Above15"] <- TRUE
# 
# sample_data(Phyloseq_EukAll)[["Above20"]] <- FALSE
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above20, "Above20"] <- TRUE
# 
# sample_data(Phyloseq_EukAll)[["Above30"]] <- FALSE
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above30, "Above30"] <- TRUE
# 
# sample_data(Phyloseq_EukAll)[["Above40"]] <- FALSE
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above40, "Above40"] <- TRUE
# 
# sample_data(Phyloseq_EukAll)[["FungiRange"]] <- "Below1"
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet1And5, "FungiRange"] <- "Bet1And5"
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet5And10, "FungiRange"] <- "Bet5And10"
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet10And15, "FungiRange"] <- "Bet10And15"
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet15And20, "FungiRange"] <- "Bet15And20"
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet20And30, "FungiRange"] <- "Bet20And30"
# sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above30FungiRange, "FungiRange"] <- "Above30"

# # Prepare to add the percent of Fungi in each sample to the metadata
# 
# DF_PercFungi <- data.frame(t(data.frame(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,])))
# colnames(DF_PercFungi) <- "PercFungi"
# DF_PercFungi[["Name_Euk"]] <- row.names(DF_PercFungi)

# Prepare to add the percent of each phylum to the metadata 

Phyloseq_ForOutstandingClass <- tax_glom(OTUPhyloseq_NormDomain, "Class", NArm = FALSE)
ASVnameFungiClass <- row.names(tax_table(Phyloseq_ForOutstandingClass)[tax_table(Phyloseq_ForOutstandingClass)[,"Division"] %in% c("Fungi"),])
DF_PercFungiClass <- data.frame(otu_table(Phyloseq_ForOutstandingClass)[ASVnameFungiClass,])

DF_PercFungiClassTax <- data.frame(tax_table(Phyloseq_ForOutstandingClass)[tax_table(Phyloseq_ForOutstandingClass)[,"Division"] %in% c("Fungi"),])
DF_PercFungiClassTax[is.na(DF_PercFungiClassTax[["Class"]]), "Class"] <- "Fungi_NA"
DF_PercFungiClass[["ASVcol"]] <- row.names(DF_PercFungiClass)
DF_PercFungiClassTax[["ASVcol"]] <- row.names(DF_PercFungiClassTax)
DF_PercFungiClassMerged <- merge(DF_PercFungiClass, DF_PercFungiClassTax, by = "ASVcol")
row.names(DF_PercFungiClassMerged) <- DF_PercFungiClassMerged[["Class"]]
DF_PercFungiClassMerged[, c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species", "ASVcol")] <- NULL
DF_PercFungiClassMerged_t <- data.frame(t(DF_PercFungiClassMerged))
DF_PercFungiClassMerged_t[["Name_Euk"]] <- row.names(DF_PercFungiClassMerged_t)

# Check if the samples are in the same order in all dataframes prior to merging
#identical(sample_data(Phyloseq_EukAll)[["Name_Euk"]], DF_PercFungi[["Name_Euk"]]) # It works
identical(sample_data(OTUPhyloseq_NormDomain)[["Name_Euk"]], DF_PercFungiClassMerged_t[["Name_Euk"]]) # It works

# Add columns to the sample_data dataframe

#sample_data(Phyloseq_EukAll)[["PercFungi"]] <- DF_PercFungi[["PercFungi"]]

sample_data(OTUPhyloseq_NormDomain)[["PercMucoromycota"]] <- DF_PercFungiClassMerged_t[["Mucoromycota"]]
sample_data(OTUPhyloseq_NormDomain)[["PercBlastocladiomycota"]] <- DF_PercFungiClassMerged_t[["Blastocladiomycota"]]
sample_data(OTUPhyloseq_NormDomain)[["PercCryptomycota"]] <- DF_PercFungiClassMerged_t[["Cryptomycota"]]
sample_data(OTUPhyloseq_NormDomain)[["PercAscomycota"]] <- DF_PercFungiClassMerged_t[["Ascomycota"]]
sample_data(OTUPhyloseq_NormDomain)[["PercChytridiomycota"]] <- DF_PercFungiClassMerged_t[["Chytridiomycota"]]
sample_data(OTUPhyloseq_NormDomain)[["PercBasidiomycota"]] <- DF_PercFungiClassMerged_t[["Basidiomycota"]]
sample_data(OTUPhyloseq_NormDomain)[["PercFungi_NA"]] <- DF_PercFungiClassMerged_t[["Fungi_NA"]]
sample_data(OTUPhyloseq_NormDomain)[["PercEntomophthoromycota"]] <- DF_PercFungiClassMerged_t[["Entomophthoromycota"]]
sample_data(OTUPhyloseq_NormDomain)[["PercFungi_X"]] <- DF_PercFungiClassMerged_t[["Fungi_X"]]

#Check if the sum of all classes equals the percent of Fungi

sample_data(OTUPhyloseq_NormDomain)[["SumPercClasses"]] <- sample_data(OTUPhyloseq_NormDomain)[["PercMucoromycota"]] + 
  sample_data(OTUPhyloseq_NormDomain)[["PercBlastocladiomycota"]] +
  sample_data(OTUPhyloseq_NormDomain)[["PercCryptomycota"]] +
  sample_data(OTUPhyloseq_NormDomain)[["PercAscomycota"]] +
  sample_data(OTUPhyloseq_NormDomain)[["PercChytridiomycota"]] +
  sample_data(OTUPhyloseq_NormDomain)[["PercBasidiomycota"]] +
  sample_data(OTUPhyloseq_NormDomain)[["PercFungi_NA"]] +
  sample_data(OTUPhyloseq_NormDomain)[["PercEntomophthoromycota"]] +
  sample_data(OTUPhyloseq_NormDomain)[["PercFungi_X"]] # It works

# #### Added in this script, not present in MetaB_SOLA_PCA.R ####
# 
# # Prepare to add the relative abundance of the 40 most abundant ASVs
# 
# Phyloseq_ForOutstandingASVs <- OTUPhyloseq_NormDomain
# Phyloseq_ForOutstandingASVs <- subset_taxa(Phyloseq_ForOutstandingASVs, Division == "Fungi")
# 
# topN <- 40
# most_abundant_taxa <- sort(taxa_sums(Phyloseq_ForOutstandingASVs), TRUE)[1:topN]
# Phyloseq_ForOutstandingASVs39ASV <- prune_taxa(names(most_abundant_taxa), Phyloseq_ForOutstandingASVs)
# 
# Phyloseq_ForOutstandingASVs39ASV_t <- data.frame(t(otu_table(Phyloseq_ForOutstandingASVs39ASV)))
# 
# # Check the order of columns 
# 
# identical(sample_data(Phyloseq_EukAll)[["Name_Euk"]], row.names(Phyloseq_ForOutstandingASVs39ASV_t)) # It works
# 
# sample_data(Phyloseq_EukAll) <- cbind(sample_data(Phyloseq_EukAll), Phyloseq_ForOutstandingASVs39ASV_t) # /!\ To be checked
# 
# #### End of added in this script, not present in MetaB_SOLA_PCA.R ####
# 
# #=============================================================================================#
# # Import the metadata from Meteo France, and the Baillaury level, and the flux cytometry data #
# #=============================================================================================#
# 
# #=============================================================================================
# # From MetaB_SOLA_PlusTurbidityMeteoFranceBaillaury3__RR7daysOffset.R
# 
# #--------------------------------------------------------------------------------------------------------------------------
# # 6- Load the Turbidity and MeteoFrance data, and merge them with the dataframe
# #--------------------------------------------------------------------------------------------------------------------------
# 
# # Import MeteoFrance data
# # RR: HAUTEUR DE PRECIPITATIONS QUOTIDIENNE en MILLIMETRES ET 1/10
# # DRR: DUREE DES PRECIPITATIONS QUOTIDIENNES en MINUTES
# # FF2M: MOYENNE DES VITESSES DU VENT A 2 METRES QUOTIDIENNE en M/S ET 1/10
# # FFM: MOYENNE DES VITESSES DU VENT A 10M QUOTIDIENNE en M/S ET 1/10
# 
# MeteoFranceMetadata <- read.csv("MeteoFrance/CapBearAPIClimatologie/CapBearConcatenated20132017.csv", sep = ",")
# MeteoFranceMetadata <- MeteoFranceMetadata[, c("Date_num", "RR", "FFM")]
# MeteoFranceMetadata$Date_num2 <- as.Date(MeteoFranceMetadata$Date_num, "%Y-%m-%d")
# MeteoFranceMetadata$Date_num <- MeteoFranceMetadata$Date_num2
# MeteoFranceMetadata$Date_num2 <- NULL
# 
# # Add a sum of the n previous days ## Modified from MetaB_SOLA_PlusTurbidityMeteoFranceBaillaury2.R to include the 7 days before instead
# 
# ndays <- 7
# SumndaysV <- rep(NA, nrow(MeteoFranceMetadata))
# for (i in (ndays:nrow(MeteoFranceMetadata))){
#   SumndaysV[i] <- sum(MeteoFranceMetadata[(i-ndays):(i-1),"RR"]) # Hence, we have the value for day i corresponding to day-7 to day-1. It may be more relevant as the samples are taken in the morning, hence the rain of day i is not so relevant.
# }
# MeteoFranceMetadata$RR7days <- SumndaysV
# 
# # Import turbidity data
# 
# TurbidityMetadata <- read.csv("CTDFromPaulMerged/Turbidity3m20mMerged.csv", sep = ",")
# TurbidityMetadata$Date_num2 <- as.Date(TurbidityMetadata$Date_num, "%Y-%m-%d")
# TurbidityMetadata$Date_num <- TurbidityMetadata$Date_num2
# TurbidityMetadata$Date_num2 <- NULL
# 
# # Import Baillaury data # 2024-02-20. I found data for the heigth of Baillaury on HydroPortail, let's use them.
# #--------------------------------------------------------------------------------------------------------------
# 
# BaillauryHeigth <- read.csv("/home/emilelaymand/Documents/Science/These/SOLA_Galand/SOLA_project_Emile_Laymand/OtherMetadata/BaillauryHeigth/Y010522001_HIXnJ(n=1_non-glissant).csv")
# 
# # Add date
# 
# BaillauryHeigth$Date <- substr(BaillauryHeigth$Date..Europe.Paris., 1, 10)
# BaillauryHeigth$Date_num <- as.Date(BaillauryHeigth$Date, "%Y-%m-%d")
# 
# # Check if there are some dubious data
# 
# unique(BaillauryHeigth$Statut)
# unique(BaillauryHeigth$Qualification)
# unique(BaillauryHeigth$Méthode)
# unique(BaillauryHeigth$Continuité)
# # Looks like there is no dubious data
# 
# # Remove columns that are useless for the rest of the analysis
# 
# BaillauryHeigth$Statut <- NULL
# BaillauryHeigth$Qualification <- NULL
# BaillauryHeigth$Méthode <- NULL
# BaillauryHeigth$Continuité <- NULL
# BaillauryHeigth$Date.de.la.mesure.du.min.max..Europe.Paris. <- NULL
# BaillauryHeigth$Date..Europe.Paris. <- NULL
# BaillauryHeigth$Date <- NULL
# colnames(BaillauryHeigth)[colnames(BaillauryHeigth) == "Valeur..en.m."] <- "BaillauryHeigthIn_m"
# 
# # Quick look at the data
# #ggplot(BaillauryHeigth, aes(x=Date_num, y=BaillauryHeigthIn_m)) +
# #  geom_point(size=0.7) +
# #  geom_line()
# 
# # Calculate derivative of BaillauryHeigthIn_m
# 
# # Check values are in the right order
# identical(BaillauryHeigth$Date_num, sort(BaillauryHeigth$Date_num)) # TRUE means we can go on
# 
# DerivHeigth <- rep(NA, nrow(BaillauryHeigth))
# for (i in (2:nrow(BaillauryHeigth))){
#   DerivHeigth[i] <- (BaillauryHeigth$BaillauryHeigthIn_m[i]-BaillauryHeigth$BaillauryHeigthIn_m[i-1])/(as.numeric(BaillauryHeigth$Date_num[i]-BaillauryHeigth$Date_num[i-1]))
# }
# BaillauryHeigth$DerivHeigth_m_j <- DerivHeigth
# 
# # Add 10*BaillauryHeigthIn_m and 10*DerivHeigth_m_j so that it is correctly scaled for plots
# 
# BaillauryHeigth$BaillauryHeigthIn_m_x10 <- 10*BaillauryHeigth$BaillauryHeigthIn_m
# BaillauryHeigth$DerivHeigth_m_j_x10 <- 10*BaillauryHeigth$DerivHeigth_m_j
# 
# # End of section imported from MetaB_SOLA_PlusTurbidityMeteoFranceBaillaury3__RR7daysOffset.R
# #====================================================================================================
# 
# # Add the Meteo France and the Baillaury data to the sample_data of Phyloseq_EukAll
# 
# # Take only the dates that are present in the metabarcoding dataset
# 
# DatesToKeep <- sample_data(Phyloseq_EukAll)[["Date_num"]]
# 
# MeteoFranceMetadata_OnlyDates <- MeteoFranceMetadata[MeteoFranceMetadata[["Date_num"]] %in% DatesToKeep,]
# dim(MeteoFranceMetadata_OnlyDates)
# 
# TurbidityMetadata_OnlyDates <- TurbidityMetadata[TurbidityMetadata[["Date_num"]] %in% DatesToKeep,]
# dim(TurbidityMetadata_OnlyDates)
# 
# BaillauryHeigth_OnlyDates <- BaillauryHeigth[BaillauryHeigth[["Date_num"]] %in% DatesToKeep,]
# dim(BaillauryHeigth_OnlyDates)
# 
# # Append the dataframes to phyloseq_EukAll by adding columns full of NAs, then filling the columns with the values from the corresponding dates
# 
# sample_data(Phyloseq_EukAll)[["RR"]] <- NA
# sample_data(Phyloseq_EukAll)[["FFM"]] <- NA
# sample_data(Phyloseq_EukAll)[["RR7days"]] <- NA
# sample_data(Phyloseq_EukAll)[["Turbidity_3m"]] <- NA
# sample_data(Phyloseq_EukAll)[["Turbidity_20m"]] <- NA
# sample_data(Phyloseq_EukAll)[["BaillauryHeigthIn_m"]] <- NA
# sample_data(Phyloseq_EukAll)[["DerivHeigth_m_j"]] <- NA
# sample_data(Phyloseq_EukAll)[["BaillauryHeigthIn_m_x10"]] <- NA
# sample_data(Phyloseq_EukAll)[["DerivHeigth_m_j_x10"]] <- NA
# 
# # For Meteo France and the Baillaury # /!\ Hard Check --> Checked: looks like it works
# 
# for (i in (1:nrow(sample_data(Phyloseq_EukAll)))){
#   sample_data(Phyloseq_EukAll)[i,"RR"] <- MeteoFranceMetadata_OnlyDates[which(MeteoFranceMetadata_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "RR"]
#   
#   sample_data(Phyloseq_EukAll)[i,"FFM"] <- MeteoFranceMetadata_OnlyDates[which(MeteoFranceMetadata_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "FFM"]
#   
#   sample_data(Phyloseq_EukAll)[i,"RR7days"] <- MeteoFranceMetadata_OnlyDates[which(MeteoFranceMetadata_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "RR7days"]
#   
#   sample_data(Phyloseq_EukAll)[i,"BaillauryHeigthIn_m"] <- BaillauryHeigth_OnlyDates[which(BaillauryHeigth_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "BaillauryHeigthIn_m"]
#   
#   sample_data(Phyloseq_EukAll)[i,"DerivHeigth_m_j"] <- BaillauryHeigth_OnlyDates[which(BaillauryHeigth_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "DerivHeigth_m_j"]
#   
#   sample_data(Phyloseq_EukAll)[i,"BaillauryHeigthIn_m_x10"] <- BaillauryHeigth_OnlyDates[which(BaillauryHeigth_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "BaillauryHeigthIn_m_x10"]
#   
#   sample_data(Phyloseq_EukAll)[i,"DerivHeigth_m_j_x10"] <- BaillauryHeigth_OnlyDates[which(BaillauryHeigth_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "DerivHeigth_m_j_x10"]
# }
# 
# # For Turbidity # /!\ Hard Check --> Checked: looks like it works
# 
# for (i in which(sample_data(Phyloseq_EukAll)[["Date_num"]] %in% TurbidityMetadata[["Date_num"]])){
#   sample_data(Phyloseq_EukAll)[i,"Turbidity_3m"] <- TurbidityMetadata_OnlyDates[which(TurbidityMetadata_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "Turbidity_3m"]
#   
#   sample_data(Phyloseq_EukAll)[i,"Turbidity_20m"] <- TurbidityMetadata_OnlyDates[which(TurbidityMetadata_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "Turbidity_20m"]
#}

#===========================================#
#  Compute percent of each phylum per month #
#===========================================#

DF_ForBoxplot <- data.frame(sample_data(OTUPhyloseq_NormDomain))

# 3 um

DF_ForBoxplot_3um <- DF_ForBoxplot[DF_ForBoxplot[["Filter_Euk"]] == "3µM",]

# Ascomycota 
boxplot(DF_ForBoxplot_3um$PercAscomycota ~ DF_ForBoxplot_3um$Month)

p1 <- ggplot(data = DF_ForBoxplot_3um, aes(x=Month, y=PercAscomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Ascomycota amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercAscomycota_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Basidiomycota 
boxplot(DF_ForBoxplot_3um$PercBasidiomycota ~ DF_ForBoxplot_3um$Month)

p1 <- ggplot(data = DF_ForBoxplot_3um, aes(x=Month, y=PercBasidiomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Basidiomycota amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercBasidiomycota_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Chytridiomycota 
boxplot(DF_ForBoxplot_3um$PercChytridiomycota ~ DF_ForBoxplot_3um$Month)

p1 <- ggplot(data = DF_ForBoxplot_3um, aes(x=Month, y=PercChytridiomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Chytridiomycota amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercChytridiomycota_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Cryptomycota 
boxplot(DF_ForBoxplot_3um$PercCryptomycota ~ DF_ForBoxplot_3um$Month)

p1 <- ggplot(data = DF_ForBoxplot_3um, aes(x=Month, y=PercCryptomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Cryptomycota amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercCryptomycota_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Unknown phylum 

DF_ForBoxplot_3um$PercFungiUnknownPhyla <- DF_ForBoxplot_3um$PercFungi_NA + DF_ForBoxplot_3um$PercFungi_X

p1 <- ggplot(data = DF_ForBoxplot_3um, aes(x=Month, y=PercFungiUnknownPhyla)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi from unkwnown phylum amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercFungiFromUnkwnonPhylum_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Other phyla 

DF_ForBoxplot_3um$PercOtherPhyla <- DF_ForBoxplot_3um$PercBlastocladiomycota + DF_ForBoxplot_3um$PercEntomophthoromycota + DF_ForBoxplot_3um$PercMucoromycota

p1 <- ggplot(data = DF_ForBoxplot_3um, aes(x=Month, y=PercOtherPhyla)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi from other phyla amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercFungiFromOtherPhyla_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# 02 um

DF_ForBoxplot_02um <- DF_ForBoxplot[DF_ForBoxplot[["Filter_Euk"]] == "Sterivex",]

# Ascomycota 
boxplot(DF_ForBoxplot_02um$PercAscomycota ~ DF_ForBoxplot_02um$Month)

p1 <- ggplot(data = DF_ForBoxplot_02um, aes(x=Month, y=PercAscomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Ascomycota amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercAscomycota_02um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Basidiomycota 
boxplot(DF_ForBoxplot_02um$PercBasidiomycota ~ DF_ForBoxplot_02um$Month)

p1 <- ggplot(data = DF_ForBoxplot_02um, aes(x=Month, y=PercBasidiomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Basidiomycota amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercBasidiomycota_02um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Chytridiomycota 
boxplot(DF_ForBoxplot_02um$PercChytridiomycota ~ DF_ForBoxplot_02um$Month)

p1 <- ggplot(data = DF_ForBoxplot_02um, aes(x=Month, y=PercChytridiomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Chytridiomycota amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercChytridiomycota_02um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Cryptomycota 
boxplot(DF_ForBoxplot_02um$PercCryptomycota ~ DF_ForBoxplot_02um$Month)

p1 <- ggplot(data = DF_ForBoxplot_02um, aes(x=Month, y=PercCryptomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Cryptomycota amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercCryptomycota_02um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Unknown phylum 

DF_ForBoxplot_02um$PercFungiUnknownPhyla <- DF_ForBoxplot_02um$PercFungi_NA + DF_ForBoxplot_02um$PercFungi_X

p1 <- ggplot(data = DF_ForBoxplot_02um, aes(x=Month, y=PercFungiUnknownPhyla)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi from unkwnown phylum amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercFungiFromUnknownPhylum_02um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Other phyla 

DF_ForBoxplot_02um$PercOtherPhyla <- DF_ForBoxplot_02um$PercBlastocladiomycota + DF_ForBoxplot_02um$PercEntomophthoromycota + DF_ForBoxplot_02um$PercMucoromycota

p1 <- ggplot(data = DF_ForBoxplot_02um, aes(x=Month, y=PercOtherPhyla)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi from other phyla amongst Fungi") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2_NormalisedByFungi/Fig2PercFungiFromOtherPhyla_02um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

#======================================#
# Compute alpha diversity within Fungi #
#======================================#

DF_ForAlphaDiv <- data.frame(otu_table(OTUPhyloseq_NormDomain))
ShannonInd <- sapply(DF_ForAlphaDiv, function(x) diversity(x, index = "shannon"))
SimpsonInd <- sapply(DF_ForAlphaDiv, function(x) diversity(x, index = "simpson"))

ShannonInd_DF <- data.frame(ShannonInd)
ShannonInd_DF[["Name_Euk"]] <- row.names(ShannonInd_DF)

SimpsonInd_DF <- data.frame(SimpsonInd)
SimpsonInd_DF[["Name_Euk"]] <- row.names(SimpsonInd_DF)

# Check columns are in the same order

identical(sample_data(OTUPhyloseq_NormDomain)[["Name_Euk"]], ShannonInd_DF[["Name_Euk"]]) # It works
identical(sample_data(OTUPhyloseq_NormDomain)[["Name_Euk"]], SimpsonInd_DF[["Name_Euk"]]) # It works

# Add the new columns

sample_data(OTUPhyloseq_NormDomain)[["ShannonInd"]] <- ShannonInd_DF$ShannonInd
sample_data(OTUPhyloseq_NormDomain)[["SimpsonInd"]] <- SimpsonInd_DF$SimpsonInd

# Extract dataframe and make plots

DF_ForAlphaDivPlot <- data.frame(sample_data(OTUPhyloseq_NormDomain))

# Make the seasons column

DF_ForAlphaDivPlot$Season <- NA
for (i in 1:dim(DF_ForAlphaDivPlot)[1]) {
  if (DF_ForAlphaDivPlot[i,"Month"] %in% c("01", "02", "03")) {
    DF_ForAlphaDivPlot[i, "Season"] <- "Winter"
  } else if (DF_ForAlphaDivPlot[i,"Month"] %in% c("04", "05", "06")) {
    DF_ForAlphaDivPlot[i, "Season"] <- "Spring"
  } else if (DF_ForAlphaDivPlot[i,"Month"] %in% c("07", "08", "09")) {
    DF_ForAlphaDivPlot[i, "Season"] <- "Summer"
  } else if (DF_ForAlphaDivPlot[i,"Month"] %in% c("10", "11", "12")) {
    DF_ForAlphaDivPlot[i, "Season"] <- "Autumn"
  } else {}
}

# 3 um

DF_ForAlphaDivPlot_3um <- DF_ForAlphaDivPlot[DF_ForAlphaDivPlot[["Filter_Euk"]] == "3µM",]

boxplot(DF_ForAlphaDivPlot_3um$ShannonInd ~ DF_ForAlphaDivPlot_3um$Month)
boxplot(DF_ForAlphaDivPlot_3um$SimpsonInd ~ DF_ForAlphaDivPlot_3um$Month)

p1 <- ggplot(data = DF_ForAlphaDivPlot_3um, aes(x=Month, y=ShannonInd, fill = Season)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Shannon index (H')") +
  ylim(c(0, 3.5)) +
  scale_fill_manual(values = c("chocolate1", "darkolivegreen3", "gold", "deepskyblue")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position="none"
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_AlphaDiversite/FigAlphaDiv_Shannon_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

p1 <- ggplot(data = DF_ForAlphaDivPlot_3um, aes(x=Month, y=SimpsonInd, fill = Season)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Gini-Simpson index (1-D)") +
  ylim(c(0, 1)) +
  scale_fill_manual(values = c("chocolate1", "darkolivegreen3", "gold", "deepskyblue")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position="none"
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_AlphaDiversite/FigAlphaDiv_GiniSimpson_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# 02 um

DF_ForAlphaDivPlot_02um <- DF_ForAlphaDivPlot[DF_ForAlphaDivPlot[["Filter_Euk"]] == "Sterivex",]

p1 <- ggplot(data = DF_ForAlphaDivPlot_02um, aes(x=Month, y=ShannonInd, fill = Season)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Shannon index (H')") +
  ylim(c(0, 3.5)) +
  scale_fill_manual(values = c("chocolate1", "darkolivegreen3", "gold", "deepskyblue")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position="none"
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_AlphaDiversite/FigAlphaDiv_Shannon_02um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

p1 <- ggplot(data = DF_ForAlphaDivPlot_02um, aes(x=Month, y=SimpsonInd, fill = Season)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Gini-Simpson index (1-D)") +
  ylim(c(0, 1)) +
  scale_fill_manual(values = c("chocolate1", "darkolivegreen3", "gold", "deepskyblue")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position="none"
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_AlphaDiversite/FigAlphaDiv_GiniSimpson_02um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

#==================#
# Stats per season #
#==================#

# Make the statistics

# Check normality
# 3 um

DF_ForAlphaDivPlot_3um %>% group_by(Season) %>% shapiro_test(SimpsonInd) # Not normal distribution
DF_ForAlphaDivPlot_3um %>% group_by(Season) %>% shapiro_test(ShannonInd) # Normal distribution

# 02 um

DF_ForAlphaDivPlot_02um %>% group_by(Season) %>% shapiro_test(SimpsonInd) # Not normal distribution (except for Spring)
DF_ForAlphaDivPlot_02um %>% group_by(Season) %>% shapiro_test(ShannonInd) # Normal distribution


#Kruskal-Wallis test 

# 3um
kruskal_test(data = DF_ForAlphaDivPlot_3um, formula = SimpsonInd ~ Season) # Significant
kruskal_test(data = DF_ForAlphaDivPlot_3um, formula = ShannonInd ~ Season) # Significant

# 02 um
kruskal_test(data = DF_ForAlphaDivPlot_02um, formula = SimpsonInd ~ Season) # Not significant, I think we do not have enough points for the test to be significant
kruskal_test(data = DF_ForAlphaDivPlot_02um, formula = ShannonInd ~ Season) # Not significant, I think we do not have enough points for the test to be significant

# Post-hoc test --> I use a Dunn test

# 3um
# Simpson
dunn_test(data = DF_ForAlphaDivPlot_3um, formula = SimpsonInd ~ Season, p.adjust.method = "bonferroni") # Autumn Summer
dunn_test(data = DF_ForAlphaDivPlot_3um, formula = SimpsonInd ~ Season, p.adjust.method = "BH") # Autumn Summer

# Shannon
dunn_test(data = DF_ForAlphaDivPlot_3um, formula = ShannonInd ~ Season, p.adjust.method = "bonferroni") # Summer Winter
dunn_test(data = DF_ForAlphaDivPlot_3um, formula = ShannonInd ~ Season, p.adjust.method = "BH") # Autumn Summer + Summer Winter

# At the month scale
#-------------------

#Kruskal-Wallis test 

# 3um
kruskal_test(data = DF_ForAlphaDivPlot_3um, formula = SimpsonInd ~ Month) # Significant
kruskal_test(data = DF_ForAlphaDivPlot_3um, formula = ShannonInd ~ Month) # Significant

# 02 um
kruskal_test(data = DF_ForAlphaDivPlot_02um, formula = SimpsonInd ~ Month) # Not significant, I think we do not have enough points for the test to be significant
kruskal_test(data = DF_ForAlphaDivPlot_02um, formula = ShannonInd ~ Month) # Not significant, I think we do not have enough points for the test to be significant

# Post-hoc test --> I use a Dunn test

# 3um
# Simpson
View(dunn_test(data = DF_ForAlphaDivPlot_3um, formula = SimpsonInd ~ Month, p.adjust.method = "bonferroni")) # Nothing significant
View(dunn_test(data = DF_ForAlphaDivPlot_3um, formula = SimpsonInd ~ Month, p.adjust.method = "BH")) # Nothing significant

# Shannon
View(dunn_test(data = DF_ForAlphaDivPlot_3um, formula = ShannonInd ~ Month, p.adjust.method = "bonferroni")) # 
View(dunn_test(data = DF_ForAlphaDivPlot_3um, formula = ShannonInd ~ Month, p.adjust.method = "BH")) # 

# Between size fractions
#-----------------------

kruskal_test(data = DF_ForAlphaDivPlot, formula = SimpsonInd ~ Filter_Euk) # 
kruskal_test(data = DF_ForAlphaDivPlot, formula = ShannonInd ~ Filter_Euk) # 


# # Use permanova to test the differences between groups --> Only one response variable, should not be used
# 
# adonis2(data = DF_ForAlphaDivPlot_3um, formula = DF_ForAlphaDivPlot_3um$ShannonInd ~ DF_ForAlphaDivPlot_3um$Season)
# adonis2(data = DF_ForAlphaDivPlot_3um, formula = DF_ForAlphaDivPlot_3um$SimpsonInd ~ DF_ForAlphaDivPlot_3um$Season)
# adonis2(data = DF_ForAlphaDivPlot_3um, formula = DF_ForAlphaDivPlot_3um$Season ~ DF_ForAlphaDivPlot_3um$)
