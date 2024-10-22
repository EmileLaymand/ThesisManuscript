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
library(gt)

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

#=============================================================================#
# End of the copy from the original MetaB_Sola.R file. Below is the new stuff.
#=============================================================================#

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
# I remove the sample PF256Euk as it has a very low number of reads and a strange composition # # Added on 2024-02-21
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#

Phyloseq_EukAll_woPF256Euk <- subset_samples(Phyloseq_EukAll, Name_Euk != "PF256Euk")
Phyloseq_EukAll <- Phyloseq_EukAll_woPF256Euk

# Select only Fungi
#------------------

### ADDED ON 2023 06 12 TO turn the salinities and temperatures equal to 0 into NA ###
sample_data(Phyloseq_EukAll)[sample_data(Phyloseq_EukAll)[["Name_Euk"]] %in% c("S109Euk", "S115Euk", "S36Euk", "PF281Euk", "PF287Euk"),"S"] <- NA
sample_data(Phyloseq_EukAll)[sample_data(Phyloseq_EukAll)[["Name_Euk"]] %in% c("S109Euk", "S115Euk", "S36Euk", "PF281Euk", "PF287Euk"),"T"] <- NA
### END OF ADDED ON 2023 06 12 TO turn the salinities and temperatures equal to 0 into NA ###

# Convert to percent at the OTU level

Phyloseq_EukAll_perc_OTU <- taxa_percentize(Phyloseq_EukAll, TaxLevel = "OTU")

# Extract Fungi from that global dataset

FungiOTUPercPhyloseq_EukAll <- subset_taxa(Phyloseq_EukAll_perc_OTU, Division == "Fungi")
FungiOTUPercPhyloseq_EukAll <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll) > 0, FungiOTUPercPhyloseq_EukAll)
FungiOTUPercPhyloseq_EukAll <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll) > 0, FungiOTUPercPhyloseq_EukAll)

# Only keep the 20 most abundant Fungi
#-------------------------------------
#-------------------------------------

topN <- 20
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll20ASV <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll)

# Both size fractions together
#-----------------------------

SumPerASV <- taxa_sums(FungiOTUPercPhyloseq_EukAll20ASV)
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll))
SumPerASVNormalised <- (SumPerASV/SumAllFungi)*100
sum(SumPerASVNormalised) # These 20 ASVs gather 82 % of normalised sequences

# Only 3 um size fraction
#------------------------

FungiOTUPercPhyloseq_EukAll_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")
FungiOTUPercPhyloseq_EukAll_3um <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)
FungiOTUPercPhyloseq_EukAll_3um <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)


topN <- 20
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll_3um20ASV <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll_3um)

# Both size fractions together
#-----------------------------

SumPerASV <- taxa_sums(FungiOTUPercPhyloseq_EukAll_3um20ASV)
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll_3um))
SumPerASVNormalised <- (SumPerASV/SumAllFungi)*100
sum(SumPerASVNormalised) # These 20 ASVs gather ~85 % of normalised sequences

# Only keep the 30 most abundant Fungi
#-------------------------------------
#-------------------------------------

topN <- 30
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll20ASV <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll)

# Both size fractions together
#-----------------------------

SumPerASV <- taxa_sums(FungiOTUPercPhyloseq_EukAll20ASV)
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll))
SumPerASVNormalised <- (SumPerASV/SumAllFungi)*100
sum(SumPerASVNormalised) # These 20 ASVs gather 87 % of normalised sequences

# Only 3 um size fraction
#------------------------

FungiOTUPercPhyloseq_EukAll_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")
FungiOTUPercPhyloseq_EukAll_3um <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)
FungiOTUPercPhyloseq_EukAll_3um <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)


topN <- 30
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll_3um20ASV <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll_3um)

# Both size fractions together
#-----------------------------

SumPerASV <- taxa_sums(FungiOTUPercPhyloseq_EukAll_3um20ASV)
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll_3um))
SumPerASVNormalised <- (SumPerASV/SumAllFungi)*100
sum(SumPerASVNormalised) # These 20 ASVs gather ~89 % of normalised sequences

# Only keep the 50 most abundant Fungi
#-------------------------------------
#-------------------------------------

topN <- 50
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll20ASV <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll)

# Both size fractions together
#-----------------------------

SumPerASV <- taxa_sums(FungiOTUPercPhyloseq_EukAll20ASV)
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll))
SumPerASVNormalised <- (SumPerASV/SumAllFungi)*100
sum(SumPerASVNormalised) # These 20 ASVs gather 92 % of normalised sequences

# Only 3 um size fraction
#------------------------

FungiOTUPercPhyloseq_EukAll_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")
FungiOTUPercPhyloseq_EukAll_3um <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)
FungiOTUPercPhyloseq_EukAll_3um <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)


topN <- 50
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll_3um20ASV <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll_3um)

# Both size fractions together
#-----------------------------

SumPerASV <- taxa_sums(FungiOTUPercPhyloseq_EukAll_3um20ASV)
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll_3um))
SumPerASVNormalised <- (SumPerASV/SumAllFungi)*100
sum(SumPerASVNormalised) # These 20 ASVs gather 93 % of normalised sequences

# Sort all ASVs, and plot the cumulative sum of relative overall abundance
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

# All samples # At ASV ranked 39th, more than 90 % of relative fungal quantities are explained
sort(taxa_sums(FungiOTUPercPhyloseq_EukAll), TRUE)

CumulSum <- cumsum(sort(taxa_sums(FungiOTUPercPhyloseq_EukAll), TRUE))
SumTotal <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll))
CumulSumPercent <- (CumulSum/SumTotal)*100
plot(CumulSumPercent, ylim = c(0,100))

CumulSumPercentDF <- data.frame(CumulSumPercent)
CumulSumPercentDF$ASVs <- row.names(CumulSumPercentDF)
CumulSumPercentDF$Rank <- 1:nrow(CumulSumPercentDF)
ggplot(CumulSumPercentDF, aes(x = Rank, y = CumulSumPercent)) +
  geom_point() + theme(axis.text.x = element_text(angle=90, hjust = 1)) + 
  scale_x_discrete(labels = CumulSumPercentDF$ASVs) + xlab("ASVs")

# Plot the total abundance of Fungi in each sample, and the cumulative abundance of the 20 most abundant ASVs
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

# 3 um
#-----

# 10 ASVs # 

FungiOTUPercPhyloseq_EukAll_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")
FungiOTUPercPhyloseq_EukAll_3um <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)
FungiOTUPercPhyloseq_EukAll_3um <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)


topN <- 10
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll_3um20ASV <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll_3um)

# Check row orders in both dataframes are the same, so that values can be directly added to the total Fungi dataframe
identical(row.names(FungiOTUPercPhyloseq_EukAll_3um), row.names(FungiOTUPercPhyloseq_EukAll_3um20ASV)) # TRUE = we can proceed

sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["FungiTotalSum"]] <- sample_sums(FungiOTUPercPhyloseq_EukAll_3um)
sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["Fungi20ASVSum"]] <- sample_sums(FungiOTUPercPhyloseq_EukAll_3um20ASV)
sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["OtherASVsSum"]] <- sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["FungiTotalSum"]] - sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["Fungi20ASVSum"]]

DFFungiPlot3um <- sample_data(FungiOTUPercPhyloseq_EukAll_3um)
ggplot(DFFungiPlot3um, aes(x=Date_num)) +
  geom_line(aes(y = FungiTotalSum, color = "FungiTotalSum"))  +
  geom_line(aes(y = Fungi20ASVSum, color = "Fungi20ASVSum")) +
  theme_bw()

DFFungiPlot3um <- sample_data(FungiOTUPercPhyloseq_EukAll_3um)
ggplot(DFFungiPlot3um, aes(x=Date_num)) +
  geom_line(aes(y = FungiTotalSum, color = "FungiTotalSum"))  +
  geom_line(aes(y = OtherASVsSum, color = "OtherASVsSum")) +
  theme_bw()

# 20 ASVs #

FungiOTUPercPhyloseq_EukAll_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")
FungiOTUPercPhyloseq_EukAll_3um <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)
FungiOTUPercPhyloseq_EukAll_3um <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)


topN <- 20
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll_3um20ASV <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll_3um)

# Check row orders in both dataframes are the same, so that values can be directly added to the total Fungi dataframe
identical(row.names(FungiOTUPercPhyloseq_EukAll_3um), row.names(FungiOTUPercPhyloseq_EukAll_3um20ASV)) # TRUE = we can proceed

sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["FungiTotalSum"]] <- sample_sums(FungiOTUPercPhyloseq_EukAll_3um)
sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["Fungi20ASVSum"]] <- sample_sums(FungiOTUPercPhyloseq_EukAll_3um20ASV)
sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["OtherASVsSum"]] <- sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["FungiTotalSum"]] - sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["Fungi20ASVSum"]]


DFFungiPlot3um <- sample_data(FungiOTUPercPhyloseq_EukAll_3um)
ggplot(DFFungiPlot3um, aes(x=Date_num)) +
  geom_line(aes(y = FungiTotalSum, color = "FungiTotalSum"))  +
  geom_line(aes(y = Fungi20ASVSum, color = "Fungi20ASVSum")) +
  theme_bw()

DFFungiPlot3um <- sample_data(FungiOTUPercPhyloseq_EukAll_3um)
ggplot(DFFungiPlot3um, aes(x=Date_num)) +
  geom_line(aes(y = FungiTotalSum, color = "FungiTotalSum"))  +
  geom_line(aes(y = OtherASVsSum, color = "OtherASVsSum")) +
  theme_bw()

# 33 ASVs # --> The 33 most abundant ASVs represent 90 % of normalized sequences across samples

FungiOTUPercPhyloseq_EukAll_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")
FungiOTUPercPhyloseq_EukAll_3um <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)
FungiOTUPercPhyloseq_EukAll_3um <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)

topN <- 33
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll_3um20ASV <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll_3um)

# Check row orders in both dataframes are the same, so that values can be directly added to the total Fungi dataframe
identical(row.names(FungiOTUPercPhyloseq_EukAll_3um), row.names(FungiOTUPercPhyloseq_EukAll_3um20ASV)) # TRUE = we can proceed

sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["FungiTotalSum"]] <- sample_sums(FungiOTUPercPhyloseq_EukAll_3um)
sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["Fungi20ASVSum"]] <- sample_sums(FungiOTUPercPhyloseq_EukAll_3um20ASV)
sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["OtherASVsSum"]] <- sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["FungiTotalSum"]] - sample_data(FungiOTUPercPhyloseq_EukAll_3um)[["Fungi20ASVSum"]]

DFFungiPlot3um <- sample_data(FungiOTUPercPhyloseq_EukAll_3um)
ggplot(DFFungiPlot3um, aes(x=Date_num)) +
  geom_line(aes(y = FungiTotalSum, color = "FungiTotalSum"))  +
  geom_line(aes(y = Fungi20ASVSum, color = "Fungi20ASVSum")) +
  theme_bw()

DFFungiPlot3um <- sample_data(FungiOTUPercPhyloseq_EukAll_3um)
ggplot(DFFungiPlot3um, aes(x=Date_num)) +
  geom_line(aes(y = FungiTotalSum, color = "FungiTotalSum"))  +
  geom_line(aes(y = OtherASVsSum, color = "OtherASVsSum")) +
  theme_bw()

# Number Of ASVs to explain X percent of cumulative abundance -- 02 um
#---------------------------------------------------------------------

FungiOTUPercPhyloseq_EukAll_02um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="Sterivex")
FungiOTUPercPhyloseq_EukAll_02um <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll_02um) > 0, FungiOTUPercPhyloseq_EukAll_02um)
FungiOTUPercPhyloseq_EukAll_02um <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll_02um) > 0, FungiOTUPercPhyloseq_EukAll_02um)

CumulSum02 <- cumsum(sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_02um), decreasing = TRUE))
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll_02um))
CumulSum02Norm <- (CumulSum02/SumAllFungi)*100 ## ASVs explain 90 % of the cumulative relative abundance 

# Number Of ASVs to explain X percent of cumulative abundance -- 3 um
#---------------------------------------------------------------------

FungiOTUPercPhyloseq_EukAll_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")
FungiOTUPercPhyloseq_EukAll_3um <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)
FungiOTUPercPhyloseq_EukAll_3um <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll_3um) > 0, FungiOTUPercPhyloseq_EukAll_3um)

CumulSum3 <- cumsum(sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_3um), decreasing = TRUE))
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll_3um))
CumulSum3Norm <- (CumulSum3/SumAllFungi)*100 ## ASVs explain 90 % of the cumulative relative abundance 

# Number Of ASVs to explain X percent of cumulative abundance -- All samples
#---------------------------------------------------------------------------

CumulSumAll <- cumsum(sort(taxa_sums(FungiOTUPercPhyloseq_EukAll), decreasing = TRUE))
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll))
CumulSumAllNorm <- (CumulSumAll/SumAllFungi)*100 ## 41 ASVs needed according to this calculation

#############################################################################
# Calculate the cumulative relative abundance of fungal phyla in all samples #
#############################################################################

# All samples
#------------

FungiOTUPercPhyloseq_EukAll_Phyla <- tax_glom(FungiOTUPercPhyloseq_EukAll, "Class", NArm = FALSE)
# DF_PercFungiClassTax[is.na(DF_PercFungiClassTax[["Class"]]), "Class"] <- "Fungi_NA"
SumPerASV <- taxa_sums(FungiOTUPercPhyloseq_EukAll_Phyla)
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll))
SumPerASVNormalised <- (SumPerASV/SumAllFungi)*100
tax_table(FungiOTUPercPhyloseq_EukAll_Phyla)
sort(SumPerASVNormalised)

# 3 um
#------
  
FungiOTUPercPhyloseq_EukAll_Phyla_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll_Phyla, Filter_Euk=="3µM")
FungiOTUPercPhyloseq_EukAll_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")

SumPerASV <- taxa_sums(FungiOTUPercPhyloseq_EukAll_Phyla_3um)
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll_3um))
SumPerASVNormalised <- (SumPerASV/SumAllFungi)*100
tax_table(FungiOTUPercPhyloseq_EukAll_Phyla_3um)
sort(SumPerASVNormalised)

# 0.2 um 
#-------

FungiOTUPercPhyloseq_EukAll_Phyla_02um <- subset_samples(FungiOTUPercPhyloseq_EukAll_Phyla, Filter_Euk=="Sterivex")
FungiOTUPercPhyloseq_EukAll_02um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="Sterivex")

SumPerASV <- taxa_sums(FungiOTUPercPhyloseq_EukAll_Phyla_02um)
SumAllFungi <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll_02um))
SumPerASVNormalised <- (SumPerASV/SumAllFungi)*100
tax_table(FungiOTUPercPhyloseq_EukAll_Phyla_02um)
sort(SumPerASVNormalised)

############################################################################################
# Number of samples required to explain 90 % of the cumulative relative abundance of Fungi #
############################################################################################

# 3 µm
#-----

FungiOTUPercPhyloseq_EukAll_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")
sample_sums(FungiOTUPercPhyloseq_EukAll_3um)

CumulSum <- cumsum(sort(sample_sums(FungiOTUPercPhyloseq_EukAll_3um), TRUE))
SumTotal <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll_3um))
CumulSumPercent <- (CumulSum/SumTotal)*100

# 0.2 µm
#-------

FungiOTUPercPhyloseq_EukAll_02um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="Sterivex")
sample_sums(FungiOTUPercPhyloseq_EukAll_02um)

CumulSum <- cumsum(sort(sample_sums(FungiOTUPercPhyloseq_EukAll_02um), TRUE))
SumTotal <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll_02um))
CumulSumPercent <- (CumulSum/SumTotal)*100

#############################################################
# List how many ASVs dominate the samples where Fungi burst #
#############################################################

sort(sample_sums(FungiOTUPercPhyloseq_EukAll), decreasing = TRUE)
FungiOTUPercPhyloseq_EukAllAbove5Perc <- subset_samples(FungiOTUPercPhyloseq_EukAll, sample_sums(FungiOTUPercPhyloseq_EukAll) >= 5)
DFlapply <- data.frame(otu_table(FungiOTUPercPhyloseq_EukAllAbove5Perc))
CountASVsBurst <- sapply(X = DFlapply, function(x) length(which(x >= 1)))
table(CountASVsBurst)

################################################################################################################################
# List the ASVs that explain 90 % of the cumulative relative abundance, sort them in decreasing order, and find their taxonomy #
################################################################################################################################

CumulSum <- cumsum(sort(taxa_sums(FungiOTUPercPhyloseq_EukAll), TRUE))
SumTotal <- sum(sample_sums(FungiOTUPercPhyloseq_EukAll))
CumulSumPercent <- (CumulSum/SumTotal)*100
which(CumulSumPercent < 90.2) # 90.1 as the one sample that allows to go over 90 % is 90.01315 % and the next one 90.20946 %  # Modified to 90.2 due to removal of 1 sample
NamesASVs90PercCumRelAb <- names(which(CumulSumPercent < 90.2))
TableForExport <- data.frame(tax_table(FungiOTUPercPhyloseq_EukAll)[NamesASVs90PercCumRelAb,])
CumulRelAbDF <- as.data.frame((sort(taxa_sums(FungiOTUPercPhyloseq_EukAll), TRUE)/SumTotal)*100)
CumulRelAbDF <- head(CumulRelAbDF, 41) # Before 2024-02-21, it was 39 # Before 2024-05-01 it was 40
colnames(CumulRelAbDF) <- "CumulativeRelativeAbundance"
identical(row.names(CumulRelAbDF), row.names(TableForExport)) # Must be TRUE for merging
TableForExport$ASV <- row.names(TableForExport)
TableForExport$CumulativeRelativeAbundance <- CumulRelAbDF$CumulativeRelativeAbundance # CHECK ALL RESULTS
TableForExportGT <- gt(TableForExport)
gtsave(TableForExportGT, filename = "TableTaxoCumulRelAb.pdf", path = "Tables")
#write.csv(x = ) 

##################################################################################
# Get sequences of the most abundant ASVs to BLAST them and check their taxonomy #
##################################################################################

EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["ASV"]] == "ASV224","X"]
EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["ASV"]] == "ASV275",] # -> ASV 275: really seems to be Candida Parapsilopsis, or to a lesser extent, Candida metapsilosis (that is apparently part of the Candida parapsilosis species complex)
EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["ASV"]] == "ASV409",]

###########################################################
# List ASVs and sequences to be blasted to check taxonomy #
###########################################################

Sequences90percASVs <- EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["ASV"]] %in% TableForExport$ASV, c("ASV", "X")]
write.csv(Sequences90percASVs, "BlastFungi/RawList40ASVsBlastNew.csv", row.names = FALSE)

#####################################################################################
# Add the percent of abundance in the 3 µm and the percent of abundance in the 3 µm #
#####################################################################################

# 3 um

FungiOTUPercPhyloseq_EukAll_3umForTable <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")
CumRelAbVect_3um <- (sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_3umForTable), decreasing = TRUE)/sum(sample_sums(FungiOTUPercPhyloseq_EukAll_3umForTable)))*100

DF_CumRelAbVect_3um <- as.data.frame(CumRelAbVect_3um)
DF_CumRelAbVect_3um$ASV <- row.names(DF_CumRelAbVect_3um)

DF_CumRelAbVect_3um_only40 <- DF_CumRelAbVect_3um[row.names(TableForExport),]

# 0.2 um

FungiOTUPercPhyloseq_EukAll_02umForTable <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="Sterivex")
CumRelAbVect_02um <- (sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_02umForTable), decreasing = TRUE)/sum(sample_sums(FungiOTUPercPhyloseq_EukAll_02umForTable)))*100

DF_CumRelAbVect_02um <- as.data.frame(CumRelAbVect_02um)
DF_CumRelAbVect_02um$ASV <- row.names(DF_CumRelAbVect_02um)

DF_CumRelAbVect_02um_only40 <- DF_CumRelAbVect_02um[row.names(TableForExport),]

# Add these dataframes to the main DF for export

TableForExport_merged <- list(TableForExport, DF_CumRelAbVect_3um_only40, DF_CumRelAbVect_02um_only40) %>% reduce(full_join) # Hard check 

# Export table to csv, to be modified in LibreOffice (to add curation, Presence in Marine Fungi, Putative role, etc.)

write.csv(TableForExport_merged, "BlastFungi/RelativeAbundanceTableRaw.csv", row.names = FALSE)

