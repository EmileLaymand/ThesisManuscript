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

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# I remove the sample PF256Euk as it has a very low number of reads and a strange composition # # Added on 2024_03_13
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

Phyloseq_EukAll_woPF256Euk <- subset_samples(Phyloseq_EukAll, Name_Euk != "PF256Euk")
Phyloseq_EukAll <- Phyloseq_EukAll_woPF256Euk

####################################################
#                    NEW SECTION                   #
####################################################

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

topN <- 20
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

# Add a sum of the n previous days

ndays <- 7
SumndaysV <- rep(NA, nrow(MeteoFranceMetadata))
for (i in (ndays:nrow(MeteoFranceMetadata))){
  SumndaysV[i] <- sum(MeteoFranceMetadata[(i-ndays+1):i,"RR"])
}
MeteoFranceMetadata$RR7days <- SumndaysV

# Import turbidity data

TurbidityMetadata <- read.csv("CTDFromPaulMerged/Turbidity3m20mMerged.csv", sep = ",")
TurbidityMetadata$Date_num2 <- as.Date(TurbidityMetadata$Date_num, "%Y-%m-%d")
TurbidityMetadata$Date_num <- TurbidityMetadata$Date_num2
TurbidityMetadata$Date_num2 <- NULL

# Import Baillaury data # Too bad for now, let's not use it

# TO BE ADDED

# Merge these new dataframes with the original metadataDF

MergedASVDF_AllFungi20ASV_MeteoTurb <- list(MergedASVDF_AllFungi20ASV, MeteoFranceMetadata, TurbidityMetadata) %>% reduce(full_join) # Hard check

########################################################################
# End of copy-paste from MetaB_SOLA_PlusTurbidityMeteoFranceBaillaury2.R
########################################################################

#----------------------------------------------------------------------------------#
#                        Add the cytometry values from SOMLIT                      #
#----------------------------------------------------------------------------------#

# Import the file

CytometryFile <- read.csv("Cytometry/Somlit_Extraction_piconano_2013_2017_modified.csv", sep = ";")

# Add a numeric date column 
CytometryFile$Date_num <- as.Date(CytometryFile$DATE, "%Y-%m-%d")

############################################################################
# End of copy-paste fromMetaB_Sola_with_Prokaryotes_BaillauryMeteoFrance.R #
############################################################################

#----------------------------------------------------------------------------------------------# 
# Extract specific groups at the Kingdom level, and correlate them with specific ASVs of Fungi # 
#----------------------------------------------------------------------------------------------# 

# 0- Select only 3 um samples
#----------------------------

Phyloseq_EukAll_perc_OTU_3umForPlot <- subset_samples(Phyloseq_EukAll_perc_OTU, Filter_Euk=="3µM")

# 1- Select Maxillopoda
#----------------------

# Glom the dataset to Kingdom level, keep only metazoa

MetazoaPhylo <- subset_taxa(Phyloseq_EukAll_perc_OTU_3umForPlot, Division == "Metazoa")

# Select only "Maxillopoda" level to get only Copepods and relatives, then glom at the genus level to see which are the most abundant groups (to see if Copepoda really dominate Maxillopoda)

MaxillopodaPhylo <- subset_taxa(MetazoaPhylo, Family == "Maxillopoda")
MaxillopodaPhylo_GenusGlom <- tax_glom(MaxillopodaPhylo, "Genus", NArm = FALSE)

sort(taxa_sums(MaxillopodaPhylo_GenusGlom), decreasing = TRUE)
tax_table(MaxillopodaPhylo_GenusGlom)[names(sort(taxa_sums(MaxillopodaPhylo_GenusGlom), decreasing = TRUE)),]

# Merge all Maxillopoda

MaxillopodaPhylo_FamilyGlom <- tax_glom(MaxillopodaPhylo, "Family", NArm = FALSE)

# 2- Select only the interesting ASVs of Fungi
#---------------------------------------------

FungiOTUPercPhyloseq_EukAll_3umForPlot <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")


InterestingFungi <- prune_taxa(c("ASV224", "ASV1613"), FungiOTUPercPhyloseq_EukAll_3umForPlot)

# 3- Merge the phyloseq objects
#------------------------------

MaxillopodaFungi <- merge_phyloseq(MaxillopodaPhylo_FamilyGlom, InterestingFungi) # /!\ Hard check this

MaxillopodaFungi_melt <- psmelt(MaxillopodaFungi)

ggplot(MaxillopodaFungi_melt, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

# 4- Plot all the Genera of Maxillopoda 1 by 1 with ggplotly
#------------------------------

MaxillopodaFungiGenus <- merge_phyloseq(MaxillopodaPhylo_GenusGlom, InterestingFungi) # /!\ Hard check this

MaxillopodaFungi_melt_Genus <- psmelt(MaxillopodaFungiGenus)

MaxillopodaFungi_melt_Genus_p <- ggplot(MaxillopodaFungi_melt_Genus, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(MaxillopodaFungi_melt_Genus_p)

# Identify the genera that were most notable genera

MostNotableGenera <- c("Oithona", "Centropages", "Acartia", "Paracalanus", "Clausocalanus", "Temora", "Calanus", "Maxillopoda_X", "Pleuromamma", "Monstrilla", "Undinula")

# Groups to keep on the plot
#ASV1     "Maxillopoda"     "Oithona"         NA                      
#ASV2     "Maxillopoda"     "Centropages"     NA                      
#ASV6     "Maxillopoda"     "Acartia"         NA                      
#ASV7     "Maxillopoda"     "Paracalanus"     NA                      
#ASV11    "Maxillopoda"     "Clausocalanus"   NA                      
#ASV21    "Maxillopoda"     "Temora"          NA                      
#ASV29    "Maxillopoda"     "Calanus"         NA                      
#ASV90    "Maxillopoda"     "Maxillopoda_X"   NA                      
#ASV431   "Maxillopoda"     "Pleuromamma"     NA                      
#ASV8935  "Maxillopoda"     "Monstrilla"      NA                      
#ASV8949  "Maxillopoda"     "Undinula"        NA     

# Plot only these genera together

MaxillopodaPhylo_GenusGlom_MostRelGenera <- subset_taxa(MaxillopodaPhylo_GenusGlom, Genus %in% MostNotableGenera)

MaxillopodaFungiGenus11 <- merge_phyloseq(MaxillopodaPhylo_GenusGlom_MostRelGenera, InterestingFungi) # /!\ Hard check this
MaxillopodaFungi_melt_Genus11 <- psmelt(MaxillopodaFungiGenus11)

MaxillopodaFungi_melt_Genus11_p <- ggplot(MaxillopodaFungi_melt_Genus11, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(MaxillopodaFungi_melt_Genus11_p)

# Look at the dynamics of Sordariomycetes (the very wide group containing Cordycipitacaea)

SordariomycetesPhylo <- subset_taxa(FungiOTUPercPhyloseq_EukAll_3umForPlot, Family == "Sordariomycetes")
SordariomycetesPhylo_glom <- tax_glom(SordariomycetesPhylo, "Family", NArm = FALSE)
taxa_names(SordariomycetesPhylo_glom) <- "Sordariomycetes" # /!\ Hard check this

SordarioMerge <- merge_phyloseq(SordariomycetesPhylo_glom, InterestingFungi) # /!\ Hard check this
SordarioMerge_melt <- psmelt(SordarioMerge)

SordarioMerge_melt_p <- ggplot(SordarioMerge_melt, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordarioMerge_melt_p)
# Conclusion: ASVs "ASV224" and "ASV1613" are a very good proxy of Sordariomycetes

# Correlation of Sordariomycetes with all families withis the dataset
#--------------------------------------------------------------------

PhyloAllFamilies <- tax_glom(Phyloseq_EukAll_perc_OTU_3umForPlot, "Family", NArm = FALSE)
PhyloAllFamiliesDFt <- data.frame(t(otu_table(PhyloAllFamilies)))
Sordariomycetes_Vect <- PhyloAllFamiliesDFt$ASV224

CorrToSordario <- sapply(PhyloAllFamiliesDFt, function(x) cor(Sordariomycetes_Vect, x))
head(sort(abs(CorrToSordario), decreasing = TRUE))
head(sort(CorrToSordario, decreasing = TRUE))

NamesMostCorr <- names(head(sort(abs(CorrToSordario), decreasing = TRUE), 100))
tax_table(PhyloAllFamilies)[NamesMostCorr,]

# Before
# ASV1463 Eukaryota Archaeplastida Rhodophyta Florideophyceae Ceramiales Wrangeliaceae
# ASV16905 Eukaryota Archaeplastida Rhodophyta Florideophyceae Colaconematales Colaconematales_X
# ASV1800 Eukaryota Stramenopiles Stramenopiles_X MAST MAST-9 MAST-9D
# ASV12345 Eukaryota Amoebozoa Lobosa Tubulinea Echinamoebida Vermamoebidae
# ASV2298 Eukaryota Opisthokonta Metazoa Craniata Craniata_X Teleostei

# After 
#ASV224   "Eukaryota" "Opisthokonta"   "Fungi"            "Ascomycota"            "Pezizomycotina"          "Sordariomycetes"                          NA    NA     
#ASV8230  "Eukaryota" "Opisthokonta"   "Fungi"            "Ascomycota"            "Pezizomycotina"          "Eurotiomycetes"                           NA    NA     
#ASV1463  "Eukaryota" "Archaeplastida" "Rhodophyta"       "Florideophyceae"       "Ceramiales"              "Wrangeliaceae"                            NA    NA     
#ASV16905 "Eukaryota" "Archaeplastida" "Rhodophyta"       "Florideophyceae"       "Colaconematales"         "Colaconematales_X"                        NA    NA     
#ASV1800  "Eukaryota" "Stramenopiles"  "Stramenopiles_X"  "MAST"                  "MAST-9"                  "MAST-9D"                                  NA    NA     
#ASV12345 "Eukaryota" "Amoebozoa"      "Lobosa"           "Tubulinea"             "Echinamoebida"           "Vermamoebidae"                            NA    NA     
#ASV2298  "Eukaryota" "Opisthokonta"   "Metazoa"          "Craniata"              "Craniata_X"              "Teleostei"                                NA    NA     
#ASV581   "Eukaryota" "Opisthokonta"   "Metazoa"          "Mollusca"              "Gastropoda"              "Heterobranchia"                           NA    NA     
#ASV2777  "Eukaryota" "Opisthokonta"   "Fungi"            "Ascomycota"            "Pezizomycotina"          "Dothideomycetes"

# Plot these most correlated groups with Sordariomycetes

SordariomycetesCorPhylo <- subset_taxa(PhyloAllFamilies, Family %in% c("Sordariomycetes", "Wrangeliaceae", "Colaconematales_X", "MAST-9D", "Vermamoebidae", "Teleostei", "Heterobranchia", "Dothideomycetes"))
SordariomycetesCorPhylo_melt <- psmelt(SordariomycetesCorPhylo)

SordariomycetesCorPhylo_melt_p <- ggplot(SordariomycetesCorPhylo_melt, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)
# Those groups are in reality not correlated, most of them are just barely present.

# With spearman correlation

CorrToSordarioSpear <- sapply(PhyloAllFamiliesDFt, function(x) cor(Sordariomycetes_Vect, x, method = "spearman"))
head(sort(abs(CorrToSordarioSpear), decreasing = TRUE))
head(sort(CorrToSordarioSpear, decreasing = TRUE))
NamesMostCorr <- names(head(sort(abs(CorrToSordarioSpear), decreasing = TRUE), 100))
tax_table(PhyloAllFamilies)[NamesMostCorr,]

# Plot these most correlated groups with Sordariomycetes

SordariomycetesCorPhylo <- subset_taxa(PhyloAllFamilies, Family %in% c("Sordariomycetes", "Eurotiomycetes", "Leotiomycetes", "Dothideomycetes", "Cryptomycotina_X", "Saccharomycetales", "Arachnida", "Dino-Group-II-Clade-23", "Vermamoebidae", "Dinophysiaceae", "Strombidiidae_M", "Vetigastropoda", "Vetigastropoda", "Cercozoa_XXX"))
SordariomycetesCorPhylo_melt <- psmelt(SordariomycetesCorPhylo)

SordariomycetesCorPhylo_melt_p <- ggplot(SordariomycetesCorPhylo_melt, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

#=============================================#
#          With the ASV level                 #
#=============================================#

DF_ASV_3umForCorAllAgAll <- data.frame(t(otu_table(Phyloseq_EukAll_perc_OTU_3umForPlot)))
CorAllAgAll <- cor(DF_ASV_3umForCorAllAgAll)
CorAllAgAll <- data.frame(CorAllAgAll)
class(CorAllAgAll)

# ASV594
Vect594 <- CorAllAgAll$ASV594
names(Vect594) <- row.names(CorAllAgAll)
head(sort(abs(Vect594), decreasing = TRUE), 100)
NamesAllAgAll594 <- names(head(sort(abs(Vect594), decreasing = TRUE), 100))
NamesAllAgAll594_20 <- names(head(sort(abs(Vect594), decreasing = TRUE), 20))
tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[NamesAllAgAll594,]

PhyloForCor <- Phyloseq_EukAll_perc_OTU_3umForPlot
tax_table(PhyloForCor) <- cbind(tax_table(PhyloForCor), ASV=row.names(tax_table(PhyloForCor)))

SordariomycetesCorPhylo <- subset_taxa(PhyloForCor, ASV %in% NamesAllAgAll594_20)
SordariomycetesCorPhylo_melt <- psmelt(SordariomycetesCorPhylo)

SordariomycetesCorPhylo_melt_p <- ggplot(SordariomycetesCorPhylo_melt, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

# ASV224

Vect594 <- CorAllAgAll$ASV224
names(Vect594) <- row.names(CorAllAgAll)
head(sort(abs(Vect594), decreasing = TRUE), 100)
NamesAllAgAll594 <- names(head(sort(abs(Vect594), decreasing = TRUE), 100))
NamesAllAgAll594_20 <- names(head(sort(abs(Vect594), decreasing = TRUE), 20))
tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[NamesAllAgAll594,]

PhyloForCor <- Phyloseq_EukAll_perc_OTU_3umForPlot
tax_table(PhyloForCor) <- cbind(tax_table(PhyloForCor), ASV=row.names(tax_table(PhyloForCor)))

SordariomycetesCorPhylo <- subset_taxa(PhyloForCor, ASV %in% NamesAllAgAll594_20)
SordariomycetesCorPhylo_melt <- psmelt(SordariomycetesCorPhylo)

SordariomycetesCorPhylo_melt_p <- ggplot(SordariomycetesCorPhylo_melt, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

# ASV1613

Vect594 <- CorAllAgAll$ASV1613
names(Vect594) <- row.names(CorAllAgAll)
head(sort(abs(Vect594), decreasing = TRUE), 100)
NamesAllAgAll594 <- names(head(sort(abs(Vect594), decreasing = TRUE), 100))
NamesAllAgAll594_20 <- names(head(sort(abs(Vect594), decreasing = TRUE), 20))
tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[NamesAllAgAll594,]

PhyloForCor <- Phyloseq_EukAll_perc_OTU_3umForPlot
tax_table(PhyloForCor) <- cbind(tax_table(PhyloForCor), ASV=row.names(tax_table(PhyloForCor)))

SordariomycetesCorPhylo <- subset_taxa(PhyloForCor, ASV %in% NamesAllAgAll594_20)
SordariomycetesCorPhylo_melt <- psmelt(SordariomycetesCorPhylo)

SordariomycetesCorPhylo_melt_p <- ggplot(SordariomycetesCorPhylo_melt, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

# Same, but spearman

DF_ASV_3umForCorAllAgAll <- data.frame(t(otu_table(Phyloseq_EukAll_perc_OTU_3umForPlot)))
CorAllAgAll <- cor(DF_ASV_3umForCorAllAgAll, method = "spearman")
CorAllAgAll <- data.frame(CorAllAgAll)
class(CorAllAgAll)

# ASV 594 
Vect594 <- CorAllAgAll$ASV594
names(Vect594) <- row.names(CorAllAgAll)
head(sort(abs(Vect594), decreasing = TRUE), 100)
NamesAllAgAll594 <- names(head(sort(abs(Vect594), decreasing = TRUE), 100))
NamesAllAgAll594_20 <- names(head(sort(abs(Vect594), decreasing = TRUE), 20))
tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[NamesAllAgAll594_20,]

PhyloForCor <- Phyloseq_EukAll_perc_OTU_3umForPlot
tax_table(PhyloForCor) <- cbind(tax_table(PhyloForCor), ASV=row.names(tax_table(PhyloForCor)))

SordariomycetesCorPhylo <- subset_taxa(PhyloForCor, ASV %in% NamesAllAgAll594_20)
SordariomycetesCorPhylo_melt <- psmelt(SordariomycetesCorPhylo)

SordariomycetesCorPhylo_melt_p <- ggplot(SordariomycetesCorPhylo_melt, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

# ASV 224
Vect594 <- CorAllAgAll$ASV224
names(Vect594) <- row.names(CorAllAgAll)
head(sort(abs(Vect594), decreasing = TRUE), 100)
NamesAllAgAll594 <- names(head(sort(abs(Vect594), decreasing = TRUE), 100))
NamesAllAgAll594_20 <- names(head(sort(abs(Vect594), decreasing = TRUE), 20))
tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[NamesAllAgAll594_20,]

PhyloForCor <- Phyloseq_EukAll_perc_OTU_3umForPlot
tax_table(PhyloForCor) <- cbind(tax_table(PhyloForCor), ASV=row.names(tax_table(PhyloForCor)))

SordariomycetesCorPhylo <- subset_taxa(PhyloForCor, ASV %in% NamesAllAgAll594_20)
SordariomycetesCorPhylo_melt <- psmelt(SordariomycetesCorPhylo)

SordariomycetesCorPhylo_melt_p <- ggplot(SordariomycetesCorPhylo_melt, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

# ASV 1613
Vect594 <- CorAllAgAll$ASV1613
names(Vect594) <- row.names(CorAllAgAll)
head(sort(abs(Vect594), decreasing = TRUE), 100)
NamesAllAgAll594 <- names(head(sort(abs(Vect594), decreasing = TRUE), 100))
NamesAllAgAll594_20 <- names(head(sort(abs(Vect594), decreasing = TRUE), 20))
tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[NamesAllAgAll594_20,]

PhyloForCor <- Phyloseq_EukAll_perc_OTU_3umForPlot
tax_table(PhyloForCor) <- cbind(tax_table(PhyloForCor), ASV=row.names(tax_table(PhyloForCor)))

SordariomycetesCorPhylo <- subset_taxa(PhyloForCor, ASV %in% NamesAllAgAll594_20)
SordariomycetesCorPhylo_melt <- psmelt(SordariomycetesCorPhylo)

SordariomycetesCorPhylo_melt_p <- ggplot(SordariomycetesCorPhylo_melt, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

#================================================================================# 
# Correlations between ASVs in the big size fraction and the small size fraction

# Merge both dataframes
Phyloseq_EukAll_perc_OTU_02umForPlot <- subset_samples(Phyloseq_EukAll_perc_OTU, Filter_Euk=="Sterivex")
Phyloseq_EukAll_perc_OTU_02umForPlot_DF <- data.frame(t(otu_table(Phyloseq_EukAll_perc_OTU_02umForPlot)))
Phyloseq_EukAll_perc_OTU_02umForPlot_DF[["Name_Euk"]] <- row.names(Phyloseq_EukAll_perc_OTU_02umForPlot_DF)
Phyloseq_EukAll_perc_OTU_02umForPlot_SamDat <- data.frame(sample_data(Phyloseq_EukAll_perc_OTU_02umForPlot))
Phyloseq_EukAll_perc_OTU_02umForPlot_SamDat <- Phyloseq_EukAll_perc_OTU_02umForPlot_SamDat[, c("Name_Euk","Date_num")]
DF_02_DateASV <- merge(Phyloseq_EukAll_perc_OTU_02umForPlot_DF, Phyloseq_EukAll_perc_OTU_02umForPlot_SamDat, by = "Name_Euk", all.x = TRUE)
DF_02_DateASV$Name_Euk <- NULL
colnames(DF_02_DateASV) <- paste(colnames(DF_02_DateASV), "_02", sep = "")
colnames(DF_02_DateASV)[colnames(DF_02_DateASV) == "Date_num_02"] <- "Date_num"

Phyloseq_EukAll_perc_OTU_3um_BIS_ForPlot <- subset_samples(Phyloseq_EukAll_perc_OTU, Filter_Euk=="3µM")
Phyloseq_EukAll_perc_OTU_3um_BIS_ForPlot_DF <- data.frame(t(otu_table(Phyloseq_EukAll_perc_OTU_3um_BIS_ForPlot)))
Phyloseq_EukAll_perc_OTU_3um_BIS_ForPlot_DF[["Name_Euk"]] <- row.names(Phyloseq_EukAll_perc_OTU_3um_BIS_ForPlot_DF)
Phyloseq_EukAll_perc_OTU_3um_BIS_ForPlot_SamDat <- data.frame(sample_data(Phyloseq_EukAll_perc_OTU_3um_BIS_ForPlot))
Phyloseq_EukAll_perc_OTU_3um_BIS_ForPlot_SamDat <- Phyloseq_EukAll_perc_OTU_3um_BIS_ForPlot_SamDat[, c("Name_Euk","Date_num")]
DF_3_DateASV <- merge(Phyloseq_EukAll_perc_OTU_3um_BIS_ForPlot_DF, Phyloseq_EukAll_perc_OTU_3um_BIS_ForPlot_SamDat, by = "Name_Euk", all.x = TRUE)
DF_3_DateASV$Name_Euk <- NULL
colnames(DF_3_DateASV) <- paste(colnames(DF_3_DateASV), "_3", sep = "")
colnames(DF_3_DateASV)[colnames(DF_3_DateASV) == "Date_num_3"] <- "Date_num"

MergedDF_02_3 <- merge(DF_02_DateASV,DF_3_DateASV, by = "Date_num", all = FALSE) # Check this out
row.names(MergedDF_02_3) <- MergedDF_02_3$Date_num
MergedDF_02_3$Date_num <- NULL

# Make some RAM space 

DF_ASV_3umForCorAllAgAll <- NULL
CorAllAgAll <- NULL

# Make correlations
#--------------------

# Spearman

#CorAllAgAll <- cor(MergedDF_02_3, method = "spearman") # Demands too much computer resources
# CorAllAgAll <- data.frame(CorAllAgAll)
#class(CorAllAgAll)

# ASV 594 

# Pearson 

Sordariomycetes_Vect <- MergedDF_02_3$ASV594_3

CorrToSordario <- sapply(MergedDF_02_3, function(x) cor(Sordariomycetes_Vect, x))
head(sort(abs(CorrToSordario), decreasing = TRUE), 100)
head(sort(CorrToSordario, decreasing = TRUE))

NamesAllAgAll594 <- names(head(sort(abs(CorrToSordario), decreasing = TRUE), 100))
NamesAllAgAll594_20 <- names(head(sort(abs(CorrToSordario), decreasing = TRUE), 20))
Names_trunc <- unique(sapply(str_split(NamesAllAgAll594_20, pattern = "_"), function(x) x[1]))
tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[Names_trunc,]

MergedDF_02_3_ForPlot <- merge(DF_02_DateASV,DF_3_DateASV, by = "Date_num", all = FALSE) # Check this out

NamesToKeep <- c(NamesAllAgAll594_20, "Date_num")
MergedDF_02_3_ForPlot_short <- MergedDF_02_3_ForPlot[,colnames(MergedDF_02_3_ForPlot) %in% NamesToKeep]
MergedDF_02_3_ForPlot_short_pivotlong <- pivot_longer(MergedDF_02_3_ForPlot_short, cols = all_of(NamesAllAgAll594_20), names_to = "OTU", values_to = "Abundance")

SordariomycetesCorPhylo_melt_p <- ggplot(MergedDF_02_3_ForPlot_short_pivotlong, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

# Spearman

Sordariomycetes_Vect <- MergedDF_02_3$ASV594_3

CorrToSordario <- sapply(MergedDF_02_3, function(x) cor(Sordariomycetes_Vect, x, method = "spearman"))
head(sort(abs(CorrToSordario), decreasing = TRUE), 100)
head(sort(CorrToSordario, decreasing = TRUE))

NamesAllAgAll594 <- names(head(sort(abs(CorrToSordario), decreasing = TRUE), 100))
#NamesAllAgAll594_20 <- names(head(sort(abs(CorrToSordario), decreasing = TRUE), 20))
NamesAllAgAll594_20 <- names(head(sort(CorrToSordario, decreasing = TRUE), 20))
Names_trunc <- unique(sapply(str_split(NamesAllAgAll594_20, pattern = "_"), function(x) x[1]))
tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[Names_trunc,]

MergedDF_02_3_ForPlot <- merge(DF_02_DateASV,DF_3_DateASV, by = "Date_num", all = FALSE) # Check this out

NamesToKeep <- c(NamesAllAgAll594_20, "Date_num")
MergedDF_02_3_ForPlot_short <- MergedDF_02_3_ForPlot[,colnames(MergedDF_02_3_ForPlot) %in% NamesToKeep]
MergedDF_02_3_ForPlot_short_pivotlong <- pivot_longer(MergedDF_02_3_ForPlot_short, cols = all_of(NamesAllAgAll594_20), names_to = "OTU", values_to = "Abundance")

SordariomycetesCorPhylo_melt_p <- ggplot(MergedDF_02_3_ForPlot_short_pivotlong, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

NamesToKeep <- c(NamesAllAgAll594_20, "ASV186_3", "Date_num") # I added this section as ASV186_02 correlated fairly well
NamesAllAgAll594_20 <- c(NamesAllAgAll594_20, "ASV186_3")
MergedDF_02_3_ForPlot_short <- MergedDF_02_3_ForPlot[,colnames(MergedDF_02_3_ForPlot) %in% NamesToKeep]
MergedDF_02_3_ForPlot_short_pivotlong <- pivot_longer(MergedDF_02_3_ForPlot_short, cols = all_of(NamesAllAgAll594_20), names_to = "OTU", values_to = "Abundance")

SordariomycetesCorPhylo_melt_p <- ggplot(MergedDF_02_3_ForPlot_short_pivotlong, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

# Look for Mamiellophyceae
ASVMamiell <- unique(row.names(tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[,"Class"] == "Mamiellophyceae",]))
ASVMamiell <- ASVMamiell[!(is.na(ASVMamiell))]
ASVMamiell <- c(ASVMamiell, "ASV224", "ASV1613", "ASV594")
ASVMamiell02 <- paste(ASVMamiell, "_02", sep = "")
ASVMamiell3 <- paste(ASVMamiell, "_3", sep = "")
ASVMamiell02_3 <- c(ASVMamiell02, ASVMamiell3)
ASVMamiellPlusDate <- c(ASVMamiell02_3, "Date_num")
MergedDF_02_3_ForPlot_short <- MergedDF_02_3_ForPlot[,colnames(MergedDF_02_3_ForPlot) %in% ASVMamiellPlusDate]
MergedDF_02_3_ForPlot_short_pivotlong <- pivot_longer(MergedDF_02_3_ForPlot_short, cols = all_of(ASVMamiell02_3), names_to = "OTU", values_to = "Abundance")

SordariomycetesCorPhylo_melt_p <- ggplot(MergedDF_02_3_ForPlot_short_pivotlong, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

# # Look for metazoa
# ASVMamiell <- unique(row.names(tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[,"Division"] == "Metazoa",]))
# ASVMamiell <- ASVMamiell[!(is.na(ASVMamiell))]
# ASVMamiell <- c(ASVMamiell, "ASV224", "ASV1613", "ASV594")
# ASVMamiell02 <- paste(ASVMamiell, "_02", sep = "")
# ASVMamiell3 <- paste(ASVMamiell, "_3", sep = "")
# ASVMamiell02_3 <- c(ASVMamiell02, ASVMamiell3)
# ASVMamiellPlusDate <- c(ASVMamiell02_3, "Date_num")
# MergedDF_02_3_ForPlot_short <- MergedDF_02_3_ForPlot[,colnames(MergedDF_02_3_ForPlot) %in% ASVMamiellPlusDate]
# MergedDF_02_3_ForPlot_short_pivotlong <- pivot_longer(MergedDF_02_3_ForPlot_short, cols = all_of(ASVMamiell02_3), names_to = "OTU", values_to = "Abundance")
# 
# SordariomycetesCorPhylo_melt_p <- ggplot(MergedDF_02_3_ForPlot_short_pivotlong, aes(x = Date_num, y = Abundance, group = OTU)) +
#   geom_point(aes(col = OTU)) +
#   geom_line(aes(col = OTU))
# 
# ggplotly(SordariomycetesCorPhylo_melt_p)

# ASV 553

ASVToInvest <- "ASV553_3"

# Pearson 

Sordariomycetes_Vect <- MergedDF_02_3[[ASVToInvest]]

CorrToSordario <- sapply(MergedDF_02_3, function(x) cor(Sordariomycetes_Vect, x))
head(sort(abs(CorrToSordario), decreasing = TRUE), 100)
head(sort(CorrToSordario, decreasing = TRUE))

NamesAllAgAll594 <- names(head(sort(abs(CorrToSordario), decreasing = TRUE), 100))
NamesAllAgAll594_20 <- names(head(sort(abs(CorrToSordario), decreasing = TRUE), 20))
Names_trunc <- unique(sapply(str_split(NamesAllAgAll594_20, pattern = "_"), function(x) x[1]))
tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[Names_trunc,]

MergedDF_02_3_ForPlot <- merge(DF_02_DateASV,DF_3_DateASV, by = "Date_num", all = FALSE) # Check this out

NamesToKeep <- c(NamesAllAgAll594_20, "Date_num")
MergedDF_02_3_ForPlot_short <- MergedDF_02_3_ForPlot[,colnames(MergedDF_02_3_ForPlot) %in% NamesToKeep]
MergedDF_02_3_ForPlot_short_pivotlong <- pivot_longer(MergedDF_02_3_ForPlot_short, cols = all_of(NamesAllAgAll594_20), names_to = "OTU", values_to = "Abundance")

SordariomycetesCorPhylo_melt_p <- ggplot(MergedDF_02_3_ForPlot_short_pivotlong, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

# Spearman

Sordariomycetes_Vect <- MergedDF_02_3[[ASVToInvest]]

CorrToSordario <- sapply(MergedDF_02_3, function(x) cor(Sordariomycetes_Vect, x, method = "spearman"))
head(sort(abs(CorrToSordario), decreasing = TRUE), 100)
head(sort(CorrToSordario, decreasing = TRUE))

NamesAllAgAll594 <- names(head(sort(abs(CorrToSordario), decreasing = TRUE), 100))
#NamesAllAgAll594_20 <- names(head(sort(abs(CorrToSordario), decreasing = TRUE), 20))
NamesAllAgAll594_20 <- names(head(sort(CorrToSordario, decreasing = TRUE), 20))
Names_trunc <- unique(sapply(str_split(NamesAllAgAll594_20, pattern = "_"), function(x) x[1]))
tax_table(Phyloseq_EukAll_perc_OTU_3umForPlot)[Names_trunc,]

MergedDF_02_3_ForPlot <- merge(DF_02_DateASV,DF_3_DateASV, by = "Date_num", all = FALSE) # Check this out

NamesToKeep <- c(NamesAllAgAll594_20, "Date_num")
MergedDF_02_3_ForPlot_short <- MergedDF_02_3_ForPlot[,colnames(MergedDF_02_3_ForPlot) %in% NamesToKeep]
MergedDF_02_3_ForPlot_short_pivotlong <- pivot_longer(MergedDF_02_3_ForPlot_short, cols = all_of(NamesAllAgAll594_20), names_to = "OTU", values_to = "Abundance")

SordariomycetesCorPhylo_melt_p <- ggplot(MergedDF_02_3_ForPlot_short_pivotlong, aes(x = Date_num, y = Abundance, group = OTU)) +
  geom_point(aes(col = OTU)) +
  geom_line(aes(col = OTU))

ggplotly(SordariomycetesCorPhylo_melt_p)

################################################################################################################
# I will now try to prove relative abundance is irrelevant, and that absolute abundance should be used instead #
# To this purpose, I will use pico-plankton measurements for flow cytometry, and try to correlate it to mamiellophyca

#Phyloseq_EukAll_perc_OTU_02umForMamiell <- subset_samples(Phyloseq_EukAll_perc_OTU, Filter_Euk=="Sterivex")
