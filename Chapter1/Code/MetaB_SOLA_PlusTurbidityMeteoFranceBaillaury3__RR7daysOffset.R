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

####################################################
#                    NEW SECTION                   #
####################################################

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# I remove the sample PF256Euk as it has a very low number of reads and a strange composition #
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

Phyloseq_EukAll_woPF256Euk <- subset_samples(Phyloseq_EukAll, Name_Euk != "PF256Euk")
Phyloseq_EukAll <- Phyloseq_EukAll_woPF256Euk

#--------------------------------------------------------------------------------------------------------------------------
# 1- Put data in percent of domain, then keep only Fungi
#--------------------------------------------------------------------------------------------------------------------------

# Convert to percent at the OTU level

Phyloseq_EukAll_perc_OTU <- taxa_percentize(Phyloseq_EukAll, TaxLevel = "OTU")

# Extract Fungi from that global dataset

FungiOTUPercPhyloseq_EukAll <- subset_taxa(Phyloseq_EukAll_perc_OTU, Division == "Fungi")
FungiOTUPercPhyloseq_EukAll <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll) > 0, FungiOTUPercPhyloseq_EukAll)
FungiOTUPercPhyloseq_EukAll <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll) > 0, FungiOTUPercPhyloseq_EukAll)

#--------------------------------------------------------------------------------------------------------------------------
# 2- Take the 20 most abundant ASVs, the total Fungi, extract them and merge them with the metadata in a dataframe
#--------------------------------------------------------------------------------------------------------------------------

# Only keep the 20 most abundant Fungi

topN <- 41 # Before, it was 20 # 2024/05/22: switched from 40 to 41
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll20ASV <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll)

#--------------------------------------------------------------------------------------------------------------------------
# 3- psmelt the phyloseq object with only the 20 most abundant taxa
#--------------------------------------------------------------------------------------------------------------------------

DFmelt20ASV <- psmelt(FungiOTUPercPhyloseq_EukAll20ASV)

#--------------------------------------------------------------------------------------------------------------------------
# 4- Agglomerate all Fungi at the Kingdom level
#--------------------------------------------------------------------------------------------------------------------------

PhyloseqFungiDomainGlom <- tax_glom(FungiOTUPercPhyloseq_EukAll, "Division", NArm = FALSE)

#--------------------------------------------------------------------------------------------------------------------------
# 5- psmelt the phyoseq object with Fungi agglomerated at the Kingdom level
#--------------------------------------------------------------------------------------------------------------------------

DFmeltFungiDomainGlom <- psmelt(PhyloseqFungiDomainGlom)
DFmeltFungiDomainGlom$OTU <- "AllFungi"

#--------------------------------------------------------------------------------------------------------------------------
# 6- Merge the two psmelted objets
#--------------------------------------------------------------------------------------------------------------------------

MergedASVDF_AllFungi20ASV <- list(DFmelt20ASV, DFmeltFungiDomainGlom) %>% reduce(full_join) # Hard check

#--------------------------------------------------------------------------------------------------------------------------
# 6- Load the Turbidity and MeteoFrance data, and merge them with the dataframe
#--------------------------------------------------------------------------------------------------------------------------

# Import MeteoFrance data
# RR: HAUTEUR DE PRECIPITATIONS QUOTIDIENNE en MILLIMETRES ET 1/10
# DRR: DUREE DES PRECIPITATIONS QUOTIDIENNES en MINUTES
# FF2M: MOYENNE DES VITESSES DU VENT A 2 METRES QUOTIDIENNE en M/S ET 1/10
# FFM: MOYENNE DES VITESSES DU VENT A 10M QUOTIDIENNE en M/S ET 1/10

MeteoFranceMetadata <- read.csv("MeteoFrance/CapBearAPIClimatologie/CapBearConcatenated20132017.csv", sep = ",")
MeteoFranceMetadata <- MeteoFranceMetadata[, c("Date_num", "RR", "FFM")]
MeteoFranceMetadata$Date_num2 <- as.Date(MeteoFranceMetadata$Date_num, "%Y-%m-%d")
MeteoFranceMetadata$Date_num <- MeteoFranceMetadata$Date_num2
MeteoFranceMetadata$Date_num2 <- NULL

# Add a sum of the n previous days ## Modified from MetaB_SOLA_PlusTurbidityMeteoFranceBaillaury2.R to include the 7 days before instead

ndays <- 7
SumndaysV <- rep(NA, nrow(MeteoFranceMetadata))
for (i in (ndays:nrow(MeteoFranceMetadata))){
  SumndaysV[i] <- sum(MeteoFranceMetadata[(i-ndays):(i-1),"RR"]) # Hence, we have the value for day i corresponding to day-7 to day-1. It may be more relevant as the samples are taken in the morning, hence the rain of day i is not so relevant.
}
MeteoFranceMetadata$RR7days <- SumndaysV

# Import turbidity data

TurbidityMetadata <- read.csv("CTDFromPaulMerged/Turbidity3m20mMerged.csv", sep = ",")
TurbidityMetadata$Date_num2 <- as.Date(TurbidityMetadata$Date_num, "%Y-%m-%d")
TurbidityMetadata$Date_num <- TurbidityMetadata$Date_num2
TurbidityMetadata$Date_num2 <- NULL

# Import Baillaury data # 2024-02-20. I found data for the heigth of Baillaury on HydroPortail, let's use them.
#--------------------------------------------------------------------------------------------------------------

BaillauryHeigth <- read.csv("/home/emilelaymand/Documents/Science/These/SOLA_Galand/SOLA_project_Emile_Laymand/OtherMetadata/BaillauryHeigth/Y010522001_HIXnJ(n=1_non-glissant).csv")

# Add date

BaillauryHeigth$Date <- substr(BaillauryHeigth$Date..Europe.Paris., 1, 10)
BaillauryHeigth$Date_num <- as.Date(BaillauryHeigth$Date, "%Y-%m-%d")

# Check if there are some dubious data

unique(BaillauryHeigth$Statut)
unique(BaillauryHeigth$Qualification)
unique(BaillauryHeigth$Méthode)
unique(BaillauryHeigth$Continuité)
# Looks like there is no dubious data

# Remove columns that are useless for the rest of the analysis

BaillauryHeigth$Statut <- NULL
BaillauryHeigth$Qualification <- NULL
BaillauryHeigth$Méthode <- NULL
BaillauryHeigth$Continuité <- NULL
BaillauryHeigth$Date.de.la.mesure.du.min.max..Europe.Paris. <- NULL
BaillauryHeigth$Date..Europe.Paris. <- NULL
BaillauryHeigth$Date <- NULL
colnames(BaillauryHeigth)[colnames(BaillauryHeigth) == "Valeur..en.m."] <- "BaillauryHeigthIn_m"

# Quick look at the data
ggplot(BaillauryHeigth, aes(x=Date_num, y=BaillauryHeigthIn_m)) +
  geom_point(size=0.7) +
  geom_line()

# Calculate derivative of BaillauryHeigthIn_m

# Check values are in the right order
identical(BaillauryHeigth$Date_num, sort(BaillauryHeigth$Date_num)) # TRUE means we can go on

DerivHeigth <- rep(NA, nrow(BaillauryHeigth))
for (i in (2:nrow(BaillauryHeigth))){
  DerivHeigth[i] <- (BaillauryHeigth$BaillauryHeigthIn_m[i]-BaillauryHeigth$BaillauryHeigthIn_m[i-1])/(as.numeric(BaillauryHeigth$Date_num[i]-BaillauryHeigth$Date_num[i-1]))
}
BaillauryHeigth$DerivHeigth_m_j <- DerivHeigth

# Add 10*BaillauryHeigthIn_m and 10*DerivHeigth_m_j so that it is correctly scaled for plots

BaillauryHeigth$BaillauryHeigthIn_m_x10 <- 10*BaillauryHeigth$BaillauryHeigthIn_m
BaillauryHeigth$DerivHeigth_m_j_x10 <- 10*BaillauryHeigth$DerivHeigth_m_j

# Merge these new dataframes with the original metadataDF
#--------------------------------------------------------

#MergedASVDF_AllFungi20ASV_MeteoTurb <- list(MergedASVDF_AllFungi20ASV, MeteoFranceMetadata, TurbidityMetadata) %>% reduce(full_join) # Hard check # Original line 

MergedASVDF_AllFungi20ASV_MeteoTurb <- list(MergedASVDF_AllFungi20ASV, MeteoFranceMetadata, TurbidityMetadata, BaillauryHeigth) %>% reduce(full_join) # New line 2024-02-20

#--------------------------------------------------------------------------------------------------------------------------
# 7- Make the plot with 3um, 0.2um, metadata as in MetaB_SOLA.R
#--------------------------------------------------------------------------------------------------------------------------

# Set the dataframe to long format, and add a column that will have three levels (0.2 um, 3 um, metadata) that will be used in facet_grid

Col2Long <- c("Abundance", "T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA", "RR", "RR7days", "FFM", "Turbidity_3m", "Turbidity_20m", "BaillauryHeigthIn_m", "DerivHeigth_m_j", "BaillauryHeigthIn_m_x10", "DerivHeigth_m_j_x10")
MergedASVDF_AllFungi20ASV_MeteoTurb_Long <- MergedASVDF_AllFungi20ASV_MeteoTurb %>%
  pivot_longer(all_of(Col2Long), names_to = "Variable", values_to = "Count") # Check

MergedASVDF_AllFungi20ASV_MeteoTurb_Long$Facet_ID <- NA
MergedASVDF_AllFungi20ASV_MeteoTurb_Long[MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Variable"]] == "Abundance", "Facet_ID"] <- MergedASVDF_AllFungi20ASV_MeteoTurb_Long[MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Variable"]] == "Abundance", "Filter_Euk"]
MergedASVDF_AllFungi20ASV_MeteoTurb_Long[which(is.na(MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Facet_ID"]])), "Facet_ID"] <- "Metadata" # Check these 3 last lines

MergedASVDF_AllFungi20ASV_MeteoTurb_Long$Group_ID <- NA
MergedASVDF_AllFungi20ASV_MeteoTurb_Long[MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Variable"]] == "Abundance", "Group_ID"] <- "Fungi"
MergedASVDF_AllFungi20ASV_MeteoTurb_Long[which(is.na(MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Group_ID"]])), "Group_ID"] <- "Metadata" # Check these 3 last lines

# Change the "Variable" column: replace the term "Abundance" by the ASV name # /!\ Check this is the right method!

MergedASVDF_AllFungi20ASV_MeteoTurb_Long[MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Variable"]] == "Abundance","Variable"] <- MergedASVDF_AllFungi20ASV_MeteoTurb_Long[MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Variable"]] == "Abundance","OTU"]

# Remove all lines for which "Count" is na (as we do not need them)

MergedASVDF_AllFungi20ASV_MeteoTurb_Long_NaRemoved <- MergedASVDF_AllFungi20ASV_MeteoTurb_Long[is.na(MergedASVDF_AllFungi20ASV_MeteoTurb_Long$Count) == FALSE,] # Check


# Change the name of the dataframe to avoid modifying all the plot in command (this is the original name in MetaB_Sola.R)

Fungi_otu_metadata_perc_OTU_Long <- MergedASVDF_AllFungi20ASV_MeteoTurb_Long_NaRemoved

# Make the plot

cc <- c(palette("Set3"), palette("Dark2"), palette("Accent"), palette("Set1"))

ggplot(Fungi_otu_metadata_perc_OTU_Long, aes(x = Date_num, y = Count, colour = Variable)) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

# Use a different legend for all facet_grid

xs <- split(Fungi_otu_metadata_perc_OTU_Long, f = Fungi_otu_metadata_perc_OTU_Long$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same graph, but remove S, T, pH, O, "CHLA", "NH4", "NO2", "NO3", "PO4", "SIO4" from the metadata graph

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "S", "pH", "O", "CHLA", "NH4", "NO2", "NO3", "PO4", "SIO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  #geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# # Plot the same graph, but with only rain in the metadata graph

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "S", "pH", "O", "CHLA", "NH4", "NO2", "NO3", "PO4", "SIO4", "FFM", "Turbidity_3m", "Turbidity_20m", "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  #geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same graph, but with only wind in the metadata graph

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "S", "pH", "O", "CHLA", "NH4", "NO2", "NO3", "PO4", "SIO4", "RR", "Turbidity_3m", "Turbidity_20m", "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  #geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same graph, but with only Turbidity_3m in the metadata graph

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "S", "pH", "O", "CHLA", "NH4", "NO2", "NO3", "PO4", "SIO4", "RR", "FFM", "Turbidity_20m", "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  #geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same graph, but with only RR7days in the metadata graph

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "S", "pH", "O", "CHLA", "NH4", "NO2", "NO3", "PO4", "SIO4", "RR", "FFM", "Turbidity_20m", "Turbidity_3m")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  #geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "pH", "O", "CHLA", "NH4", "NO2", "NO3", "PO4", "SIO4", "RR", "FFM", "Turbidity_20m", "Turbidity_3m")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  #geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

#-------------------------------------------------------------------
# Plots ASV per ASV
#-------------------------------------------------------------------

#========================
# ASV 224
#========================

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("ASV224", "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("ASV224", "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("ASV224", "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("ASV224", "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with the Height of Baillaury

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("ASV224", "BaillauryHeigthIn_m")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain and river

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("ASV224", "RR", "BaillauryHeigthIn_m_x10")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

#======================================================# 
#         No Metadata -- Added on 2024/04/18           #
#======================================================# 

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("ASV224", "RR", "BaillauryHeigthIn_m_x10")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )


print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/01_ASV224.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV275
#=========================================

ASVToTest <- "ASV275"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
  #annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/02_ASV275.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV409
#=========================================

ASVToTest <- "ASV409"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/03_ASV409.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV1613
#=========================================

ASVToTest <- "ASV1613"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/04_ASV1613.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV1016
#=========================================

ASVToTest <- "ASV1016"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/05_ASV1016.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV166
#=========================================

ASVToTest <- "ASV166"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/06_ASV166.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()


#=========================================
# ASV847
#=========================================

ASVToTest <- "ASV847"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/07_ASV847.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV594
#=========================================

ASVToTest <- "ASV594"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem with the heigth of Baillaury

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "BaillauryHeigthIn_m")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain and river

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR", "BaillauryHeigthIn_m_x10")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/08_ASV594.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV367
#=========================================

ASVToTest <- "ASV367"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/09_ASV367.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV2336
#=========================================

ASVToTest <- "ASV2336"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/10_ASV2336.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV553
#=========================================

ASVToTest <- "ASV553"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/11_ASV553.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV2777
#=========================================

ASVToTest <- "ASV2777"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/12_ASV2777.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

# #=========================================
# # ASV9123
# #=========================================
# 
# ASVToTest <- "ASV9123"
# 
# # Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)
# 
# Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]
# 
# xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
# p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
#   geom_point(size = 0.7) +
#   geom_line() +
#   scale_colour_manual(values = cc) +
#   facet_grid(rows = vars(Facet_ID), scales = "free") +
#   scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
#   theme_bw()
# 
# p2 <- p1 %+% xs[["Metadata"]]
# 
# gA <- ggplotGrob(p1)
# gB <- ggplotGrob(p2)
# grid::grid.newpage()
# grid::grid.draw(rbind(gA, gB))
# 
# # Idem, but with rain
# 
# Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]
# 
# xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
# p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
#   geom_point(size=0.7) +
#   geom_line() +
#   scale_colour_manual(values = cc) +
#   facet_grid(rows = vars(Facet_ID), scales = "free") +
#   scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
#   theme_bw()
# 
# p2 <- p1 %+% xs[["Metadata"]]
# 
# gA <- ggplotGrob(p1)
# gB <- ggplotGrob(p2)
# grid::grid.newpage()
# grid::grid.draw(rbind(gA, gB))
# 
# # Idem, but with NO3
# 
# Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]
# 
# xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
# p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
#   geom_point(size=0.7) +
#   geom_line() +
#   scale_colour_manual(values = cc) +
#   facet_grid(rows = vars(Facet_ID), scales = "free") +
#   scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
#   theme_bw()
# 
# p2 <- p1 %+% xs[["Metadata"]]
# 
# gA <- ggplotGrob(p1)
# gB <- ggplotGrob(p2)
# grid::grid.newpage()
# grid::grid.draw(rbind(gA, gB))
# 
# # Idem, but with PO4
# 
# Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]
# 
# xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
# p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
#   geom_point(size=0.7) +
#   geom_line() +
#   scale_colour_manual(values = cc) +
#   facet_grid(rows = vars(Facet_ID), scales = "free") +
#   scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
#   theme_bw()
# 
# p2 <- p1 %+% xs[["Metadata"]]
# 
# gA <- ggplotGrob(p1)
# gB <- ggplotGrob(p2)
# grid::grid.newpage()
# grid::grid.draw(rbind(gA, gB))
# 
# # Without metadata
# 
# xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
# 
# xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
# xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"
# 
# p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
#   geom_point(size=0.7) +
#   geom_line() +
#   scale_colour_manual(values = "black") +
#   facet_grid(rows = vars(Facet_ID), scales = "free") +
#   scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
#                date_minor_breaks = "1 month", date_labels="%Y") +
#   #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
#   #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
#   labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
#   theme_bw() +
#   theme(legend.position = "none", 
#         axis.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         strip.text.y = element_text(size = 20)
#   )
# #annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")
# 
# print(p1)
# 
# ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/13_ASV9123.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
# dev.off()

#=========================================
# ASV774
#=========================================

ASVToTest <- "ASV774"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/14_ASV774.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV5211
#=========================================

ASVToTest <- "ASV5211"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/15_ASV5211.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV843
#=========================================

ASVToTest <- "ASV843"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/16_ASV843.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV324
#=========================================

ASVToTest <- "ASV324"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/17_ASV324.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV1724
#=========================================

ASVToTest <- "ASV1724"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/18_ASV1724.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV6458
#=========================================

ASVToTest <- "ASV6458"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/19_ASV6458.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV420
#=========================================

ASVToTest <- "ASV420"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/20_ASV420.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV1991
#=========================================

ASVToTest <- "ASV1991"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/21_ASV1991.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV1633
#=========================================

ASVToTest <- "ASV1633"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/22_ASV1633.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV8230
#=========================================

ASVToTest <- "ASV8230"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/23_ASV8230.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV460
#=========================================

ASVToTest <- "ASV460"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/24_ASV460.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV505
#=========================================

ASVToTest <- "ASV505"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/25_ASV505.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV5782
#=========================================

ASVToTest <- "ASV5782"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/26_ASV5782.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

# #=========================================
# # ASV9513
# #=========================================
# 
# ASVToTest <- "ASV9513"
# 
# # Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)
# 
# Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]
# 
# xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
# p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
#   geom_point(size = 0.7) +
#   geom_line() +
#   scale_colour_manual(values = cc) +
#   facet_grid(rows = vars(Facet_ID), scales = "free") +
#   scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
#   theme_bw()
# 
# p2 <- p1 %+% xs[["Metadata"]]
# 
# gA <- ggplotGrob(p1)
# gB <- ggplotGrob(p2)
# grid::grid.newpage()
# grid::grid.draw(rbind(gA, gB))
# 
# # Idem, but with rain
# 
# Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]
# 
# xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
# p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
#   geom_point(size=0.7) +
#   geom_line() +
#   scale_colour_manual(values = cc) +
#   facet_grid(rows = vars(Facet_ID), scales = "free") +
#   scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
#   theme_bw()
# 
# p2 <- p1 %+% xs[["Metadata"]]
# 
# gA <- ggplotGrob(p1)
# gB <- ggplotGrob(p2)
# grid::grid.newpage()
# grid::grid.draw(rbind(gA, gB))
# 
# # Idem, but with NO3
# 
# Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]
# 
# xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
# p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
#   geom_point(size=0.7) +
#   geom_line() +
#   scale_colour_manual(values = cc) +
#   facet_grid(rows = vars(Facet_ID), scales = "free") +
#   scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
#   theme_bw()
# 
# p2 <- p1 %+% xs[["Metadata"]]
# 
# gA <- ggplotGrob(p1)
# gB <- ggplotGrob(p2)
# grid::grid.newpage()
# grid::grid.draw(rbind(gA, gB))
# 
# # Idem, but with PO4
# 
# Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]
# 
# xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
# p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
#   geom_point(size=0.7) +
#   geom_line() +
#   scale_colour_manual(values = cc) +
#   facet_grid(rows = vars(Facet_ID), scales = "free") +
#   scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
#   theme_bw()
# 
# p2 <- p1 %+% xs[["Metadata"]]
# 
# gA <- ggplotGrob(p1)
# gB <- ggplotGrob(p2)
# grid::grid.newpage()
# grid::grid.draw(rbind(gA, gB))
# 
# # Without metadata
# 
# xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
# 
# xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
# xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"
# 
# p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
#   geom_point(size=0.7) +
#   geom_line() +
#   scale_colour_manual(values = "black") +
#   facet_grid(rows = vars(Facet_ID), scales = "free") +
#   scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
#                date_minor_breaks = "1 month", date_labels="%Y") +
#   #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
#   #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
#   labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
#   theme_bw() +
#   theme(legend.position = "none", 
#         axis.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         strip.text.y = element_text(size = 20)
#   )
# #annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")
# 
# print(p1)
# 
# ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/27_ASV9513.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
# dev.off()

#=========================================
# ASV1163
#=========================================

ASVToTest <- "ASV1163"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/28_ASV1163.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV8944
#=========================================

ASVToTest <- "ASV8944"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/29_ASV8944.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV2737
#=========================================

ASVToTest <- "ASV2737"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/30_ASV2737.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV866
#=========================================

ASVToTest <- "ASV866"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/31_ASV866.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV3073
#=========================================

ASVToTest <- "ASV3073"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/32_ASV3073.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV1421
#=========================================

ASVToTest <- "ASV1421"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/33_ASV1421.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV8958
#=========================================

ASVToTest <- "ASV8958"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/34_ASV8958.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV2651
#=========================================

ASVToTest <- "ASV2651"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/35_ASV2651.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV8967
#=========================================

ASVToTest <- "ASV8967"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/36_ASV8967.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV8969
#=========================================

ASVToTest <- "ASV8969"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain and river

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR", "BaillauryHeigthIn_m_x10")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/37_ASV8969.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV6352
#=========================================

ASVToTest <- "ASV6352"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/38_ASV6352.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV3729
#=========================================

ASVToTest <- "ASV3729"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/39_ASV3729.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV971
#=========================================

ASVToTest <- "ASV971"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/40_ASV971.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV9009
#=========================================

ASVToTest <- "ASV9009"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/41_ASV9009.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV5590
#=========================================

ASVToTest <- "ASV5590"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/42_ASV5590.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

#=========================================
# ASV1367
#=========================================

ASVToTest <- "ASV1367"

# Plot the same graph, but with only RR7days and salinity in the metadata graph (to see if there is a correlation between these variables)

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with rain

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "RR7days")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with NO3

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "NO3")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Idem, but with PO4

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c(ASVToTest, "PO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Without metadata

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)

xs$Fungi[xs$Fungi[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Fungi[xs$Fungi[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size=0.7) +
  geom_line() +
  scale_colour_manual(values = "black") +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks = "1 year", 
               date_minor_breaks = "1 month", date_labels="%Y") +
  #geom_label(aes(hjust = 1, vjust = 1, label=paste(ASVToTest))) +
  #scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01")), date_breaks="1 month", date_labels="%Y") +
  labs(x = "Date", y = "Relative abundance (% of eukaryotes)") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.y = element_text(size = 20)
  )
#annotate("label", x=as.Date("2013-01-01"),y=Inf,hjust=0.5,vjust=0.7, label="test")

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/FigureSup_Profiles40FirstASVs/43_ASV1367.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()

# Correlations

# 

ListDFASVs <- split(MergedASVDF_AllFungi20ASV_MeteoTurb, f = MergedASVDF_AllFungi20ASV_MeteoTurb$OTU)

# sapply(ListDFASVs, function(x) cor(x[["Abundance"]], x[["RR"]])**2)
# sapply(ListDFASVs, function(x) cor(x[["Abundance"]], x[["RR"]], method = "spearman")**2)

VariablesToTest <- c("T", "S", "pH", "O", "CHLA", "NH4", "NO2", "NO3", "PO4", "SIO4", "RR", "FFM", "Turbidity_20m", "Turbidity_3m")

CorAllMetadata <- sapply(VariablesToTest, function(y) sapply(ListDFASVs, function(x) cor(x[["Abundance"]], x[[y]], use = "pairwise.complete.obs")**2))
#CorAllMetadata <- sapply(VariablesToTest, function(y) sapply(ListDFASVs, function(x) cor(x[["Abundance"]], x[[y]], use = "complete.obs")**2))

which(CorAllMetadata > 0.1)
CorAllMetadata > 0.1

# Correlation between rain and salinity --> No correlation between rain and salinity

cor(MergedASVDF_AllFungi20ASV_MeteoTurb$FFM, MergedASVDF_AllFungi20ASV_MeteoTurb$S, use = "pairwise.complete.obs") # 0.01514427
cor(MergedASVDF_AllFungi20ASV_MeteoTurb$RR, MergedASVDF_AllFungi20ASV_MeteoTurb$S, use = "pairwise.complete.obs") # 0.05092254
cor(MergedASVDF_AllFungi20ASV_MeteoTurb$RR7days, MergedASVDF_AllFungi20ASV_MeteoTurb$S, use = "pairwise.complete.obs") # 0.03629584
cor(MergedASVDF_AllFungi20ASV_MeteoTurb$RR7days, MergedASVDF_AllFungi20ASV_MeteoTurb$S, use = "pairwise.complete.obs", method = "spearman") # -0.1126421
cor(MergedASVDF_AllFungi20ASV_MeteoTurb$RR, MergedASVDF_AllFungi20ASV_MeteoTurb$S, use = "pairwise.complete.obs", method = "spearman") # 0.05804371

# Correlation between salinity and temperature --> There is a positive correlation

cor(MergedASVDF_AllFungi20ASV_MeteoTurb$T, MergedASVDF_AllFungi20ASV_MeteoTurb$S, use = "pairwise.complete.obs", method = "spearman") # 0.0729181
cor(MergedASVDF_AllFungi20ASV_MeteoTurb$T, MergedASVDF_AllFungi20ASV_MeteoTurb$S, use = "pairwise.complete.obs", method = "pearson") # 0.5351823 









#########################################################################################################
#          Investigate links between some groups' decrease (namely, dinoflagellates) and Fungi          #
#########################################################################################################

#--------------------------------------------------------------------------------------------------------------------------
# 1- Add a phyloseq object with ASVs glommed to the Domain level
#--------------------------------------------------------------------------------------------------------------------------

Phyloseq_EukAll_perc_OTU_glom_domain <- tax_glom(Phyloseq_EukAll_perc_OTU, "Division", NArm = FALSE)
MeltedDomainDF <- psmelt(Phyloseq_EukAll_perc_OTU_glom_domain)
MeltedDomainDF$OTU <- MeltedDomainDF$Division # Check

MergedASVDF_AllDomainGlom20ASVFungi <- list(DFmelt20ASV, MeltedDomainDF) %>% reduce(full_join) # Very Hard check

# MergedASVDF_AllDomainGlom20ASVFungi_MeteoTurb <- list(MergedASVDF_AllDomainGlom20ASVFungi, MeteoFranceMetadata, TurbidityMetadata) %>% reduce(full_join) # Very Hard check
MergedASVDF_AllDomainGlom20ASVFungi_MeteoTurb <- list(MergedASVDF_AllDomainGlom20ASVFungi, MeteoFranceMetadata, TurbidityMetadata, BaillauryHeigth) %>% reduce(full_join) # Very Hard check # Added on 2024-02-20. Must be checked

#--------------------------------------------------------------------------------------------------------------------------
# 2- Change name for compatibility with code upper in the script
#--------------------------------------------------------------------------------------------------------------------------

MergedASVDF_AllFungi20ASV_MeteoTurb <- MergedASVDF_AllDomainGlom20ASVFungi_MeteoTurb

#--------------------------------------------------------------------------------------------------------------------------
# 7- Make the plot with 3um, 0.2um, metadata as in MetaB_SOLA.R
#--------------------------------------------------------------------------------------------------------------------------

# Set the dataframe to long format, and add a column that will have three levels (0.2 um, 3 um, metadata) that will be used in facet_grid

Col2Long <- c("Abundance", "T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA", "RR", "RR7days", "FFM", "Turbidity_3m", "Turbidity_20m", "BaillauryHeigthIn_m", "DerivHeigth_m_j", "BaillauryHeigthIn_m_x10", "DerivHeigth_m_j_x10")
MergedASVDF_AllFungi20ASV_MeteoTurb_Long <- MergedASVDF_AllFungi20ASV_MeteoTurb %>%
  pivot_longer(all_of(Col2Long), names_to = "Variable", values_to = "Count") # Check

MergedASVDF_AllFungi20ASV_MeteoTurb_Long$Facet_ID <- NA
MergedASVDF_AllFungi20ASV_MeteoTurb_Long[MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Variable"]] == "Abundance", "Facet_ID"] <- MergedASVDF_AllFungi20ASV_MeteoTurb_Long[MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Variable"]] == "Abundance", "Filter_Euk"]
MergedASVDF_AllFungi20ASV_MeteoTurb_Long[which(is.na(MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Facet_ID"]])), "Facet_ID"] <- "Metadata" # Check these 3 last lines

MergedASVDF_AllFungi20ASV_MeteoTurb_Long$Group_ID <- NA
MergedASVDF_AllFungi20ASV_MeteoTurb_Long[MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Variable"]] == "Abundance", "Group_ID"] <- "Fungi"
MergedASVDF_AllFungi20ASV_MeteoTurb_Long[which(is.na(MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Group_ID"]])), "Group_ID"] <- "Metadata" # Check these 3 last lines

# Change the "Variable" column: replace the term "Abundance" by the ASV name # /!\ Check this is the right method!

MergedASVDF_AllFungi20ASV_MeteoTurb_Long[MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Variable"]] == "Abundance","Variable"] <- MergedASVDF_AllFungi20ASV_MeteoTurb_Long[MergedASVDF_AllFungi20ASV_MeteoTurb_Long[["Variable"]] == "Abundance","OTU"]

# Remove all lines for which "Count" is na (as we do not need them)

MergedASVDF_AllFungi20ASV_MeteoTurb_Long_NaRemoved <- MergedASVDF_AllFungi20ASV_MeteoTurb_Long[is.na(MergedASVDF_AllFungi20ASV_MeteoTurb_Long$Count) == FALSE,] # Check


# Change the name of the dataframe to avoid modifying all the plot in command (this is the original name in MetaB_Sola.R)

Fungi_otu_metadata_perc_OTU_Long <- MergedASVDF_AllFungi20ASV_MeteoTurb_Long_NaRemoved

# Make the plot

cc <- c(palette("Set3"), palette("Dark2"), palette("Accent"), palette("Set1"))

ggplot(Fungi_otu_metadata_perc_OTU_Long, aes(x = Date_num, y = Count, colour = Variable)) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

# Use a different legend for all facet_grid

xs <- split(Fungi_otu_metadata_perc_OTU_Long, f = Fungi_otu_metadata_perc_OTU_Long$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same graph, but remove S, T, pH, O, "CHLA", "NH4", "NO2", "NO3", "PO4", "SIO4" from the metadata graph

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("Dinoflagellata", "Metazoa", "Fungi", "ASV224", "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  #geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

#-----------------------------
# ASV594
#-----------------------------

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("Dinoflagellata", "Metazoa", "Fungi", "ASV594", "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  #geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))


#================================================================================
# Quick test: plot salinity, rain and Baillaury heigth to test their correlation 
#================================================================================

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("Dinoflagellata", "Metazoa", "Fungi", "ASV594", "RR7days", "BaillauryHeigthIn_m", "S", "DerivHeigth_m_j")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  #geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  scale_x_date(limits = as.Date(c("2013-01-01", "2017-06-01"))) +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

MergedASVDF_AllFungi20ASV_MeteoTurb_ForPlot <- MergedASVDF_AllFungi20ASV_MeteoTurb
MergedASVDF_AllFungi20ASV_MeteoTurb_ForPlot$BaillauryHeigthIn_m10 <- 10*MergedASVDF_AllFungi20ASV_MeteoTurb_ForPlot$BaillauryHeigthIn_m
MergedASVDF_AllFungi20ASV_MeteoTurb_ForPlot$DerivHeigth_m_j_50 <- 10*MergedASVDF_AllFungi20ASV_MeteoTurb_ForPlot$DerivHeigth_m_j
  
ggplot(MergedASVDF_AllFungi20ASV_MeteoTurb_ForPlot, aes(x=Date_num)) +
  geom_line(aes(y=BaillauryHeigthIn_m10, color="Baillaury Heigth 10")) +
  geom_line(aes(y=S, color = "Salinity")) +
  geom_point(aes(y=S, color = "Salinity")) +
  geom_line(aes(y=RR7days, color = "RR7days")) +
  geom_line(aes(y=DerivHeigth_m_j_50, color = "DerivHeigth_m_j_50"))

ggplot(MergedASVDF_AllFungi20ASV_MeteoTurb_ForPlot, aes(x=Date_num)) +
  geom_line(aes(y=BaillauryHeigthIn_m10, color="Baillaury Heigth 10")) +
  geom_line(aes(y=S, color = "Salinity")) +
  geom_point(aes(y=S, color = "Salinity")) +
  geom_line(aes(y=RR, color = "RR")) +
  geom_line(aes(y=DerivHeigth_m_j_50, color = "DerivHeigth_m_j_50"))
  
  
ggplot(MergedASVDF_AllFungi20ASV_MeteoTurb_ForPlot, aes(x=Date_num)) +
  geom_line(aes(y=BaillauryHeigthIn_m10, color="Baillaury Heigth")) +
  geom_line(aes(y=S, color = "Salinity")) +
  geom_point(aes(y=S, color = "Salinity")) +
  geom_line(aes(y=RR7days, color = "RR"))

