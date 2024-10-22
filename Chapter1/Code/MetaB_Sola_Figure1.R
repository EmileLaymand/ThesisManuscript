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

######################################
# Section imported from metaB_Sola.R # 
######################################

cc <- c(palette("Set3"), palette("Dark2"), palette("Accent"), palette("Set1"))

#====
# Division level 
#====

OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)

OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "Division")

# Extract OTU table and sample_data from the phyloseq object to plot them as geom_point

OTUPhyloseq_NormDomain_DF <- psmelt(OTUPhyloseq_NormDomain)

# Set the dataframe to long format, and add a column that will have three levels (0.2 um, 3 um, metadata) that will be used in facet_grid

Col2Long <- c("Abundance", "T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")
OTUPhyloseq_NormDomain_DF_Long <- OTUPhyloseq_NormDomain_DF %>%
  pivot_longer(all_of(Col2Long), names_to = "Variable", values_to = "Count")

OTUPhyloseq_NormDomain_DF_Long$Facet_ID <- NA
OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Facet_ID"] <- OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Filter_Euk"]
OTUPhyloseq_NormDomain_DF_Long[which(is.na(OTUPhyloseq_NormDomain_DF_Long[["Facet_ID"]])), "Facet_ID"] <- "Metadata"

OTUPhyloseq_NormDomain_DF_Long$Group_ID <- NA
OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Group_ID"] <- "Species"
OTUPhyloseq_NormDomain_DF_Long[which(is.na(OTUPhyloseq_NormDomain_DF_Long[["Group_ID"]])), "Group_ID"] <- "Metadata"

# Rename "Abundance" in the "Variable" column to be the Domain name

OTUPhyloseq_NormDomain_DF_Long_Div <- OTUPhyloseq_NormDomain_DF_Long
OTUPhyloseq_NormDomain_DF_Long_Div[OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] == "Abundance", "Variable"] <- OTUPhyloseq_NormDomain_DF_Long_Div[OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] == "Abundance", "Division"]


xs <- split(OTUPhyloseq_NormDomain_DF_Long_Div, f = OTUPhyloseq_NormDomain_DF_Long_Div$Group_ID)
p1 <- ggplot(xs[["Species"]] ,aes(x = Date_num, y = Count, colour = Variable)) +
  geom_point(size = 0.5) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Same but remove pH as a metadata as one point is very far from mean

OTUPhyloseq_NormDomain_DF_Long_Div_butpH <- OTUPhyloseq_NormDomain_DF_Long_Div[!(OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] %in% c("pH")),]

xs <- split(OTUPhyloseq_NormDomain_DF_Long_Div_butpH, f = OTUPhyloseq_NormDomain_DF_Long_Div_butpH$Group_ID)
p1 <- ggplot(xs[["Species"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.5) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Same but remove pH as a metadata as one point is very far from mean and keep only Fungi --> Nice

OTUPhyloseq_NormDomain_DF_Long_Div_butpH_Fung <- OTUPhyloseq_NormDomain_DF_Long_Div[!(OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] %in% c("pH", "Dinoflagellata", "Metazoa", "Stramenopiles_X", "Chlorophyta", "Ochrophyta", "Ciliophora", "Cercozoa", "Cryptophyta", "Radiolaria", "Mesomycetozoa", "Haptophyta", "Telonemia", "Choanoflagellida", "Foraminifera", "Rhodophyta", "Opisthokonta_X", "Streptophyta", "Picozoa", "Katablepharidophyta", "Apicomplexa", "Centroheliozoa", "Discoba", "Apusomonadidae", "Conosa", "Metamonada", "Hilomonadea", "Perkinsea", "Alveolata_X", "Lobosa", "Breviatea", NA)),]

xs <- split(OTUPhyloseq_NormDomain_DF_Long_Div_butpH_Fung, f = OTUPhyloseq_NormDomain_DF_Long_Div_butpH_Fung$Group_ID)
p1 <- ggplot(xs[["Species"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.5) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))


# Same but only NO2 as metadata

OTUPhyloseq_NormDomain_DF_Long_Div_NO2 <- OTUPhyloseq_NormDomain_DF_Long_Div[!(OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] %in% c("T", "S", "pH", "O", "NH4", "NO3", "PO4", "CHLA", "SIO4")),]

xs <- split(OTUPhyloseq_NormDomain_DF_Long_Div_NO2, f = OTUPhyloseq_NormDomain_DF_Long_Div_NO2$Group_ID)
p1 <- ggplot(xs[["Species"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.5) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

#====
# Fungi Class (actually phylum) level 
#====

OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)

OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "OTU")

# Z-score the sample_data columns (T, Sal, PO4, etc.)

sample_data(OTUPhyloseq_NormDomain)[, c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")] <- lapply(sample_data(OTUPhyloseq_NormDomain)[, c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")], function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))

# Select only Fungi

OTUPhyloseq_NormDomain <- subset_taxa(OTUPhyloseq_NormDomain, Division == "Fungi")
OTUPhyloseq_NormDomain <- tax_glom(OTUPhyloseq_NormDomain, "Class", NArm = FALSE)

# Extract OTU table and sample_data from the phyloseq object to plot them as geom_point

OTUPhyloseq_NormDomain_DF <- psmelt(OTUPhyloseq_NormDomain)

# Set the dataframe to long format, and add a column that will have three levels (0.2 um, 3 um, metadata) that will be used in facet_grid

Col2Long <- c("Abundance", "T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")
OTUPhyloseq_NormDomain_DF_Long <- OTUPhyloseq_NormDomain_DF %>%
  pivot_longer(all_of(Col2Long), names_to = "Variable", values_to = "Count")

OTUPhyloseq_NormDomain_DF_Long$Facet_ID <- NA
OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Facet_ID"] <- OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Filter_Euk"]
OTUPhyloseq_NormDomain_DF_Long[which(is.na(OTUPhyloseq_NormDomain_DF_Long[["Facet_ID"]])), "Facet_ID"] <- "Metadata"

OTUPhyloseq_NormDomain_DF_Long$Group_ID <- NA
OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Group_ID"] <- "Species"
OTUPhyloseq_NormDomain_DF_Long[which(is.na(OTUPhyloseq_NormDomain_DF_Long[["Group_ID"]])), "Group_ID"] <- "Metadata"

# Rename "Abundance" in the "Variable" column to be the Class name

OTUPhyloseq_NormDomain_DF_Long_Div <- OTUPhyloseq_NormDomain_DF_Long
OTUPhyloseq_NormDomain_DF_Long_Div[OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] == "Abundance", "Variable"] <- OTUPhyloseq_NormDomain_DF_Long_Div[OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] == "Abundance", "Class"]

# Remove pH as a metadata as one point is very far from mean --> Nice

OTUPhyloseq_NormDomain_DF_Long_Div_butpH <- OTUPhyloseq_NormDomain_DF_Long_Div[!(OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] %in% c("pH")),]

xs <- split(OTUPhyloseq_NormDomain_DF_Long_Div_butpH, f = OTUPhyloseq_NormDomain_DF_Long_Div_butpH$Group_ID)
p1 <- ggplot(xs[["Species"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.5) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

##########################################################################################################
#       Added on 2024/04/16 - Merge low abundance phyla together, and Fungi_X and Fungi NA together      #
##########################################################################################################

IndexNA <- which(is.na(tax_table(OTUPhyloseq_NormDomain)[,"Class"]))
IndexFungiX <- which(tax_table(OTUPhyloseq_NormDomain)[,"Class"] == "Fungi_X")

tax_table(OTUPhyloseq_NormDomain)[IndexNA,"Class"] <- "Unknown phylum"
tax_table(OTUPhyloseq_NormDomain)[IndexFungiX,"Class"] <- "Unknown phylum"

IndexEntomophthoromycota <- which(tax_table(OTUPhyloseq_NormDomain)[,"Class"] == "Entomophthoromycota")
IndexMucoromycota <- which(tax_table(OTUPhyloseq_NormDomain)[,"Class"] == "Mucoromycota")
IndexBlastocladiomycota <- which(tax_table(OTUPhyloseq_NormDomain)[,"Class"] == "Blastocladiomycota")

tax_table(OTUPhyloseq_NormDomain)[IndexEntomophthoromycota,"Class"] <- "Other phylum"
tax_table(OTUPhyloseq_NormDomain)[IndexMucoromycota,"Class"] <- "Other phylum"
tax_table(OTUPhyloseq_NormDomain)[IndexBlastocladiomycota,"Class"] <- "Other phylum"


OTUPhyloseq_NormDomain <- tax_glom(OTUPhyloseq_NormDomain, "Class", NArm = FALSE)

# Extract OTU table and sample_data from the phyloseq object to plot them as geom_point

OTUPhyloseq_NormDomain_DF <- psmelt(OTUPhyloseq_NormDomain)

# Set the dataframe to long format, and add a column that will have three levels (0.2 um, 3 um, metadata) that will be used in facet_grid

Col2Long <- c("Abundance", "T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")
OTUPhyloseq_NormDomain_DF_Long <- OTUPhyloseq_NormDomain_DF %>%
  pivot_longer(all_of(Col2Long), names_to = "Variable", values_to = "Count")

OTUPhyloseq_NormDomain_DF_Long$Facet_ID <- NA
OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Facet_ID"] <- OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Filter_Euk"]
OTUPhyloseq_NormDomain_DF_Long[which(is.na(OTUPhyloseq_NormDomain_DF_Long[["Facet_ID"]])), "Facet_ID"] <- "Metadata"

OTUPhyloseq_NormDomain_DF_Long$Group_ID <- NA
OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Group_ID"] <- "Species"
OTUPhyloseq_NormDomain_DF_Long[which(is.na(OTUPhyloseq_NormDomain_DF_Long[["Group_ID"]])), "Group_ID"] <- "Metadata"

# Rename "Abundance" in the "Variable" column to be the Class name

OTUPhyloseq_NormDomain_DF_Long_Div <- OTUPhyloseq_NormDomain_DF_Long
OTUPhyloseq_NormDomain_DF_Long_Div[OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] == "Abundance", "Variable"] <- OTUPhyloseq_NormDomain_DF_Long_Div[OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] == "Abundance", "Class"]

# Remove pH as a metadata as one point is very far from mean --> Nice

OTUPhyloseq_NormDomain_DF_Long_Div_butpH <- OTUPhyloseq_NormDomain_DF_Long_Div[!(OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] %in% c("pH")),]

xs <- split(OTUPhyloseq_NormDomain_DF_Long_Div_butpH, f = OTUPhyloseq_NormDomain_DF_Long_Div_butpH$Group_ID)
p1 <- ggplot(xs[["Species"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.5) +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

CustomPalette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "darkorange1", "black")

##########################################################################################################
#                                       End of Added on 2024/04/16                                       #
##########################################################################################################


# Only Fungi sequences, not metadata

xs <- split(OTUPhyloseq_NormDomain_DF_Long_Div_butpH, f = OTUPhyloseq_NormDomain_DF_Long_Div_butpH$Group_ID)

xs$Species[xs$Species[["Facet_ID"]] == "3µM","Facet_ID"] <- "> 3 µm"
xs$Species[xs$Species[["Facet_ID"]] == "Sterivex","Facet_ID"] <- "0.2 - 3 µm"

p1 <- ggplot(xs[["Species"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point(size = 0.5) +
  geom_line() +
  scale_colour_manual(values = CustomPalette) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  labs(x = "Date", y = "Relative abundance of fungal phyla (% of eukaryotes)") +
  scale_x_date(date_breaks="1 month", date_labels="%b\n%Y") +
  theme_bw()

print(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/LinePlot_Figure1.pdf", plot=p1, device = cairo_pdf(), width=20, height=10)
dev.off()
