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

# Define "outstanding samples" for easier spotting on the plots (outstanding = the relative abundance of FUngi is high)

OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)
OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "OTU")
Phyloseq_ForOutstanding <- tax_glom(OTUPhyloseq_NormDomain, "Division", NArm = FALSE)
ASVnameFungi <- row.names(tax_table(Phyloseq_ForOutstanding)[tax_table(Phyloseq_ForOutstanding)[,"Division"] %in% c("Fungi"),])

Above1 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 1])
Above5 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 5])
Above10 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 10])
Above15 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 15])
Above20 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 20])
Above30 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 30])
Above40 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 40])
Bet1And5 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 1 & otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] <= 5])
Bet5And10 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 5  & otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] <= 10])
Bet10And15 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 10  & otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] <= 15])
Bet15And20 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 15  & otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] <= 20])
Bet20And30 <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 20  & otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] <= 30])
Above30FungiRange <- colnames(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi, otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,] > 30])

sample_data(Phyloseq_EukAll)[["Above1"]] <- FALSE
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above1, "Above1"] <- TRUE

sample_data(Phyloseq_EukAll)[["Above5"]] <- FALSE
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above5, "Above5"] <- TRUE

sample_data(Phyloseq_EukAll)[["Above10"]] <- FALSE
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above10, "Above10"] <- TRUE

sample_data(Phyloseq_EukAll)[["Above15"]] <- FALSE
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above15, "Above15"] <- TRUE

sample_data(Phyloseq_EukAll)[["Above20"]] <- FALSE
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above20, "Above20"] <- TRUE

sample_data(Phyloseq_EukAll)[["Above30"]] <- FALSE
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above30, "Above30"] <- TRUE

sample_data(Phyloseq_EukAll)[["Above40"]] <- FALSE
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above40, "Above40"] <- TRUE

sample_data(Phyloseq_EukAll)[["FungiRange"]] <- "Below1"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet1And5, "FungiRange"] <- "Bet1And5"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet5And10, "FungiRange"] <- "Bet5And10"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet10And15, "FungiRange"] <- "Bet10And15"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet15And20, "FungiRange"] <- "Bet15And20"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet20And30, "FungiRange"] <- "Bet20And30"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above30FungiRange, "FungiRange"] <- "Above30"

# Prepare to add the percent of Fungi in each sample to the metadata 

DF_PercFungi <- data.frame(t(data.frame(otu_table(Phyloseq_ForOutstanding)[ASVnameFungi,])))
colnames(DF_PercFungi) <- "PercFungi"
DF_PercFungi[["Name_Euk"]] <- row.names(DF_PercFungi)

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
identical(sample_data(Phyloseq_EukAll)[["Name_Euk"]], DF_PercFungi[["Name_Euk"]]) # It works
identical(sample_data(Phyloseq_EukAll)[["Name_Euk"]], DF_PercFungiClassMerged_t[["Name_Euk"]]) # It works

# Add columns to the sample_data dataframe

sample_data(Phyloseq_EukAll)[["PercFungi"]] <- DF_PercFungi[["PercFungi"]]

sample_data(Phyloseq_EukAll)[["PercMucoromycota"]] <- DF_PercFungiClassMerged_t[["Mucoromycota"]]
sample_data(Phyloseq_EukAll)[["PercBlastocladiomycota"]] <- DF_PercFungiClassMerged_t[["Blastocladiomycota"]]
sample_data(Phyloseq_EukAll)[["PercCryptomycota"]] <- DF_PercFungiClassMerged_t[["Cryptomycota"]]
sample_data(Phyloseq_EukAll)[["PercAscomycota"]] <- DF_PercFungiClassMerged_t[["Ascomycota"]]
sample_data(Phyloseq_EukAll)[["PercChytridiomycota"]] <- DF_PercFungiClassMerged_t[["Chytridiomycota"]]
sample_data(Phyloseq_EukAll)[["PercBasidiomycota"]] <- DF_PercFungiClassMerged_t[["Basidiomycota"]]
sample_data(Phyloseq_EukAll)[["PercFungi_NA"]] <- DF_PercFungiClassMerged_t[["Fungi_NA"]]
sample_data(Phyloseq_EukAll)[["PercEntomophthoromycota"]] <- DF_PercFungiClassMerged_t[["Entomophthoromycota"]]
sample_data(Phyloseq_EukAll)[["PercFungi_X"]] <- DF_PercFungiClassMerged_t[["Fungi_X"]]

#Check if the sum of all classes equals the percent of Fungi

sample_data(Phyloseq_EukAll)[["SumPercClasses"]] <- sample_data(Phyloseq_EukAll)[["PercMucoromycota"]] + 
  sample_data(Phyloseq_EukAll)[["PercBlastocladiomycota"]] +
  sample_data(Phyloseq_EukAll)[["PercCryptomycota"]] +
  sample_data(Phyloseq_EukAll)[["PercAscomycota"]] +
  sample_data(Phyloseq_EukAll)[["PercChytridiomycota"]] +
  sample_data(Phyloseq_EukAll)[["PercBasidiomycota"]] +
  sample_data(Phyloseq_EukAll)[["PercFungi_NA"]] +
  sample_data(Phyloseq_EukAll)[["PercEntomophthoromycota"]] +
  sample_data(Phyloseq_EukAll)[["PercFungi_X"]] # It works

# Add the anomaly of salinity as a variable. The anomaly of salinity is here defined as the median of salinity for each month over all the years
  
# Create a dataframe with the closest values of debit prior to the day of sampling
#---------------------------------------------------------------------------------

HydroBaillaury <- read.csv("OtherMetadata/BaillauryMaillolHydroportailAccessed20231205DateSplit.csv", sep = ",", dec = ".")

# Check if there are the exact dates
PhyloseqForBaillaury <- data.frame(sample_data(Phyloseq_EukAll))
PhyloseqForBaillaury[["Date_num"]]
HydroBaillauryOnlyDatesSampling <- HydroBaillaury[HydroBaillaury[["Date_TU"]] %in% PhyloseqForBaillaury[["Date_num"]],]
length(unique(HydroBaillauryOnlyDatesSampling[["Date_TU"]])) # 144
length(PhyloseqForBaillaury[["Date_num"]]) # 308


HydroBaillauryShlagEdit <- read.csv("OtherMetadata/BaillauryMaillolHydroportailAccessed20231205DateSplitForWorkLibreOfficeOnlySampleDatesOneDateDatesReplaced.csv", sep = ",", dec = ".")
colnames(HydroBaillauryShlagEdit)[colnames(HydroBaillauryShlagEdit) %in% c("Date_TU")] <- "Date_num"
PhyloseqShlagHydro <- merge(PhyloseqForBaillaury, HydroBaillauryShlagEdit, by = "Date_num")

# 3 um 
#-----

PhyloseqShlagHydro3um <- PhyloseqShlagHydro[PhyloseqShlagHydro[["Filter_Euk"]] == "3µM", ]
plot(PhyloseqShlagHydro3um$PercFungi ~ PhyloseqShlagHydro3um$Valeur_m3_s)
cor(y = PhyloseqShlagHydro3um$PercFungi, x =PhyloseqShlagHydro3um$Valeur_m3_s, method = "spearman")

PhyloseqShlagHydroNo0TS3um <- PhyloseqShlagHydro3um[!(PhyloseqShlagHydro3um[["Name_Euk"]] %in% c("PF281Euk", "PF287Euk")), ]
plot(PhyloseqShlagHydroNo0TS3um$S ~ PhyloseqShlagHydroNo0TS3um$Valeur_m3_s)
ggplot(PhyloseqShlagHydroNo0TS3um, aes(x= Valeur_m3_s, y = S)) +
  geom_point(aes(colour = Month))
cor(y = PhyloseqShlagHydro3um$PercFungi, x =PhyloseqShlagHydro3um$S, use = "complete.obs", method = "spearman") # Not sure if complete.obs is the right option

#Correlation salinity with the rest
plot(PhyloseqShlagHydroNo0TS3um$PercFungi ~ PhyloseqShlagHydroNo0TS3um$S)
plot(PhyloseqShlagHydroNo0TS3um$PercMucoromycota ~ PhyloseqShlagHydroNo0TS3um$S)
plot(PhyloseqShlagHydroNo0TS3um$PercChytridiomycota ~ PhyloseqShlagHydroNo0TS3um$S)
plot(PhyloseqShlagHydroNo0TS3um$PercCryptomycota ~ PhyloseqShlagHydroNo0TS3um$S)
plot(PhyloseqShlagHydroNo0TS3um$PercCryptomycota ~ PhyloseqShlagHydroNo0TS3um$CHLA)
plot(PhyloseqShlagHydroNo0TS3um$PercAscomycota ~ PhyloseqShlagHydroNo0TS3um$S)
plot(PhyloseqShlagHydroNo0TS3um$PercBasidiomycota ~ PhyloseqShlagHydroNo0TS3um$S)
plot(PhyloseqShlagHydroNo0TS3um$NO3 ~ PhyloseqShlagHydroNo0TS3um$S)
plot(PhyloseqShlagHydroNo0TS3um$NO2 ~ PhyloseqShlagHydroNo0TS3um$S)
plot(PhyloseqShlagHydroNo0TS3um$PO4 ~ PhyloseqShlagHydroNo0TS3um$S)
plot(PhyloseqShlagHydroNo0TS3um$SIO4 ~ PhyloseqShlagHydroNo0TS3um$S)


ggplot(PhyloseqShlagHydroNo0TS3um, aes(x= S, y = PercBasidiomycota)) +
  geom_point()

plot(PhyloseqShlagHydroNo0TS3um$PercFungi ~ PhyloseqShlagHydroNo0TS3um$Month)
plot(PhyloseqShlagHydroNo0TS3um$PercAscomycota ~ PhyloseqShlagHydroNo0TS3um$Month)
plot(PhyloseqShlagHydroNo0TS3um$PercAscomycota ~ PhyloseqShlagHydroNo0TS3um$Month)

boxplot(PhyloseqShlagHydroNo0TS3um$PercFungi ~ PhyloseqShlagHydroNo0TS3um$Month)
boxplot(PhyloseqShlagHydroNo0TS3um$PercAscomycota ~ PhyloseqShlagHydroNo0TS3um$Month)
boxplot(PhyloseqShlagHydroNo0TS3um$PercBasidiomycota ~ PhyloseqShlagHydroNo0TS3um$Month)
boxplot(PhyloseqShlagHydroNo0TS3um$PercChytridiomycota ~ PhyloseqShlagHydroNo0TS3um$Month)
boxplot(PhyloseqShlagHydroNo0TS3um$PercCryptomycota ~ PhyloseqShlagHydroNo0TS3um$Month)
boxplot(PhyloseqShlagHydroNo0TS3um$Valeur_m3_s ~ PhyloseqShlagHydroNo0TS3um$Month)

boxplot(PhyloseqShlagHydroNo0TS3um$S ~ PhyloseqShlagHydroNo0TS3um$Month)


#===================================================================================#
# End of the copy from the original MetaB_SOLA_NMDS.R file. Below is the new stuff. #
#===================================================================================#

#######
# PCA #
#######

MetadataAllSamples <- data.frame(sample_data(Phyloseq_EukAll))

# Keep only 3um samples
#----------------------

MetadataAll3um <- MetadataAllSamples[MetadataAllSamples[["Filter_Euk"]] == "3µM",]
MetadataAll3umOnlySomeVar <- MetadataAll3um[,c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA", "PercFungi", "PercMucoromycota", "PercBlastocladiomycota", "PercCryptomycota", "PercAscomycota", "PercChytridiomycota", "PercBasidiomycota", "PercFungi_NA", "PercEntomophthoromycota", "PercFungi_X")]

# Remove samples with at least one NA for a variable
MetadataAll3umOnlySomeVarNoNA <- MetadataAll3umOnlySomeVar[complete.cases(MetadataAll3umOnlySomeVar), ]
dim(MetadataAll3umOnlySomeVar)
dim(MetadataAll3umOnlySomeVarNoNA) # 16 samples lack at least one parameter, ~1/10 of the samples

# Remove samples that have weird values
plot(MetadataAll3umOnlySomeVarNoNA$T)
MetadataAll3umOnlySomeVarNoNA[MetadataAll3umOnlySomeVarNoNA[["T"]] < 5,] # For two samples, the Temperature and salinity is 0: they must be removed
MetadataAll3umOnlySomeVarNoNA <- MetadataAll3umOnlySomeVarNoNA[!(row.names(MetadataAll3umOnlySomeVarNoNA) %in% c("PF281Euk", "PF287Euk")), ]
dim(MetadataAll3umOnlySomeVarNoNA) 
plot(MetadataAll3umOnlySomeVarNoNA$T)

plot(MetadataAll3umOnlySomeVarNoNA$S)
plot(MetadataAll3umOnlySomeVarNoNA$O)
plot(MetadataAll3umOnlySomeVarNoNA$pH)
plot(MetadataAll3umOnlySomeVarNoNA$NH4)
plot(MetadataAll3umOnlySomeVarNoNA$NO3) # Kinda cyclic (and I think the samples are not exactly in the right order)
plot(MetadataAll3umOnlySomeVarNoNA$NO2) # Kinda cyclic(and I think the samples are not exactly in the right order)
plot(MetadataAll3umOnlySomeVarNoNA$PO4) # Kinda cyclic, but less clear (and I think the samples are not exactly in the right order)
plot(MetadataAll3umOnlySomeVarNoNA$SIO4) 
plot(MetadataAll3umOnlySomeVarNoNA$CHLA) # Kinda cyclic, but way less (and I think the samples are not exactly in the right order)

# PCAs
#-----

# MetadataAll3umOnlyGoodSamples <- MetadataAll3um[row.names(MetadataAll3umOnlySomeVarNoNA),]
# 
# PCA_3um <- prcomp(MetadataAll3umOnlySomeVarNoNA, center = TRUE, scale. = TRUE)
# PCA_3um$x <- -1*PCA_3um$x  # /!\ Check these rotations must be done
# PCA_3um$rotation <- -1*PCA_3um$rotation # /!\ Check these rotations must be done
# plot(PCA_3um)
# biplot(PCA_3um, scale = 0)
# PCA_3um$sdev^2 / sum(PCA_3um$sdev^2)
# p <- autoplot(PCA_3um, data = MetadataAll3umOnlyGoodSamples, colour = "FungiRange", loadings = TRUE, loadings.colour = "blue", 
#               loadings.label = TRUE, loadings.label.size = 3, label = TRUE, label.label = "Name_Euk")
# ggplotly(p) # Check the loadings are the right thing to plot

# PCA with FactoMineR

PCA_3umNice <- PCA(MetadataAll3umOnlySomeVarNoNA, scale.unit = TRUE, ncp = 5, graph = TRUE, quanti.sup = 11:20)
PCA_3umNice.eig.val <- get_eigenvalue(PCA_3umNice)
fviz_eig(PCA_3umNice, addlabels = TRUE, ylim = c(0, 50))
PCA_3umNice.var <- get_pca_var(PCA_3umNice)
fviz_pca_var(PCA_3umNice, col.var = "black") # Correlation circle
corrplot(PCA_3umNice.var$cos2, is.corr=FALSE)
fviz_pca_var(PCA_3umNice, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
corrplot(PCA_3umNice.var$contrib, is.corr=FALSE) # Contribution of each variables to each dimension
fviz_pca_var(PCA_3umNice, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
corrplot(PCA_3umNice$quanti.sup$cos2, is.corr=FALSE)

#note: there is a possibility to add ellipses to plots.
#note: there is the possibility to plot the samples on the graph

# Idem, but with all 10 dimensions displayed on the plots

PCA_3umNice <- PCA(MetadataAll3umOnlySomeVarNoNA, scale.unit = TRUE, ncp = 10, graph = TRUE, quanti.sup = 11:20)
PCA_3umNice.eig.val <- get_eigenvalue(PCA_3umNice)
fviz_eig(PCA_3umNice, addlabels = TRUE, ylim = c(0, 50))
PCA_3umNice.var <- get_pca_var(PCA_3umNice)
fviz_pca_var(PCA_3umNice, col.var = "black") # Correlation circle
corrplot(PCA_3umNice.var$cos2, is.corr=FALSE)
fviz_pca_var(PCA_3umNice, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
corrplot(PCA_3umNice.var$contrib, is.corr=FALSE) # Contribution of each variables to each dimension
fviz_pca_var(PCA_3umNice, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
corrplot(PCA_3umNice$quanti.sup$cos2, is.corr=FALSE)

# Test if salinity can reasonably be extrapolated (i.e., easily predicted by a sine function)
MetadataAll3umRmOutliers <- MetadataAll3um[!(row.names(MetadataAll3um) %in% c("PF281Euk", "PF287Euk")), ] # Remove the two outliers
ggplot(data = MetadataAll3umRmOutliers, aes(x=Date_num, y=S)) +
  geom_point(na.rm=TRUE) +
  geom_line()

# Test if the samples with most Fungi are overrepresented within the ANR samples # NOT FINISHED
Phyloseq_EukAll3um <- subset_samples(Phyloseq_EukAll, Filter_Euk=="3µM")

###############################################
# Statistics per phylum and per size fraction #
###############################################

MetadataAllSamplesStatistics <- MetadataAllSamples

# Rename columns with an _ to no _ because it is required as separator later on

colnames(MetadataAllSamplesStatistics)[colnames(MetadataAllSamplesStatistics) == 'PercFungi_NA'] <- 'PercFungiNA'
colnames(MetadataAllSamplesStatistics)[colnames(MetadataAllSamplesStatistics) == 'PercFungi_X'] <- 'PercFungiX'

# All samples 
#------------

Col2Long <- c("PercFungi", "PercMucoromycota", "PercBlastocladiomycota", "PercCryptomycota", "PercAscomycota", "PercChytridiomycota", "PercBasidiomycota", "PercFungiNA", "PercEntomophthoromycota", "PercFungiX")
MetadataAllSamplesStatistics_DF_Long <- MetadataAllSamplesStatistics %>%
  pivot_longer(all_of(Col2Long), names_to = "Variable", values_to = "Count")

# Synthetic plot, but unreadable (too much information)
ggplot(MetadataAllSamplesStatistics_DF_Long, aes(x=Variable,y=Count, fill=factor(Filter_Euk))) +
  geom_boxplot() + 
  labs(fill = "Filter_Euk") + 
  geom_point(position=position_jitterdodge(),alpha=0.3) +
  theme_bw(base_size = 16)

MetadataAllSamplesStatisticsSomeVars <- MetadataAllSamplesStatistics[, c("PercFungi", "PercMucoromycota", "PercBlastocladiomycota", "PercCryptomycota", "PercAscomycota", "PercChytridiomycota", "PercBasidiomycota", "PercFungiNA", "PercEntomophthoromycota", "PercFungiX")]


DFStatsAllSamples <- MetadataAllSamplesStatisticsSomeVars %>% summarise(across(where(is.numeric), .fns = 
                          list(min = min,
                               median = median,
                               mean = mean,
                               stdev = sd,
                               q25 = ~quantile(., 0.25),
                               q75 = ~quantile(., 0.75),
                               max = max))) %>%
  pivot_longer(everything(), names_sep='_', names_to=c('variable', '.value'))

# Only 3 um
#----------

MetadataAllSamplesStatistics <- MetadataAllSamplesStatistics[MetadataAllSamplesStatistics[["Filter_Euk"]] == "3µM", ]

Col2Long <- c("PercFungi", "PercMucoromycota", "PercBlastocladiomycota", "PercCryptomycota", "PercAscomycota", "PercChytridiomycota", "PercBasidiomycota", "PercFungiNA", "PercEntomophthoromycota", "PercFungiX")
MetadataAllSamplesStatistics_DF_Long <- MetadataAllSamplesStatistics %>%
  pivot_longer(all_of(Col2Long), names_to = "Variable", values_to = "Count")

# Synthetic plot, but unreadable (too much information)
ggplot(MetadataAllSamplesStatistics_DF_Long, aes(x=Variable,y=Count, fill=factor(Filter_Euk))) +
  geom_boxplot() + 
  labs(fill = "Filter_Euk") + 
  geom_point(position=position_jitterdodge(),alpha=0.3) +
  theme_bw(base_size = 16)

MetadataAllSamplesStatisticsSomeVars <- MetadataAllSamplesStatistics[, c("PercFungi", "PercMucoromycota", "PercBlastocladiomycota", "PercCryptomycota", "PercAscomycota", "PercChytridiomycota", "PercBasidiomycota", "PercFungiNA", "PercEntomophthoromycota", "PercFungiX")]


DFStatsAllSamples <- MetadataAllSamplesStatisticsSomeVars %>% summarise(across(where(is.numeric), .fns = 
                                                                                 list(min = min,
                                                                                      median = median,
                                                                                      mean = mean,
                                                                                      stdev = sd,
                                                                                      q25 = ~quantile(., 0.25),
                                                                                      q75 = ~quantile(., 0.75),
                                                                                      max = max))) %>%
  pivot_longer(everything(), names_sep='_', names_to=c('variable', '.value'))

# Only 2 um
#----------

# Need to reset the dataset, as it was modified to contain only 3 um samples

MetadataAllSamplesStatistics <- MetadataAllSamples

# Rename columns with an _ to no _ because it is required as separator later on

colnames(MetadataAllSamplesStatistics)[colnames(MetadataAllSamplesStatistics) == 'PercFungi_NA'] <- 'PercFungiNA'
colnames(MetadataAllSamplesStatistics)[colnames(MetadataAllSamplesStatistics) == 'PercFungi_X'] <- 'PercFungiX'

MetadataAllSamplesStatistics <- MetadataAllSamplesStatistics[MetadataAllSamplesStatistics[["Filter_Euk"]] == "Sterivex", ]

Col2Long <- c("PercFungi", "PercMucoromycota", "PercBlastocladiomycota", "PercCryptomycota", "PercAscomycota", "PercChytridiomycota", "PercBasidiomycota", "PercFungiNA", "PercEntomophthoromycota", "PercFungiX")
MetadataAllSamplesStatistics_DF_Long <- MetadataAllSamplesStatistics %>%
  pivot_longer(all_of(Col2Long), names_to = "Variable", values_to = "Count")

# Synthetic plot, but unreadable (too much information)
ggplot(MetadataAllSamplesStatistics_DF_Long, aes(x=Variable,y=Count, fill=factor(Filter_Euk))) +
  geom_boxplot() + 
  labs(fill = "Filter_Euk") + 
  geom_point(position=position_jitterdodge(),alpha=0.3) +
  theme_bw(base_size = 16)

MetadataAllSamplesStatisticsSomeVars <- MetadataAllSamplesStatistics[, c("PercFungi", "PercMucoromycota", "PercBlastocladiomycota", "PercCryptomycota", "PercAscomycota", "PercChytridiomycota", "PercBasidiomycota", "PercFungiNA", "PercEntomophthoromycota", "PercFungiX")]


DFStatsAllSamples <- MetadataAllSamplesStatisticsSomeVars %>% summarise(across(where(is.numeric), .fns = 
                                                                                 list(min = min,
                                                                                      median = median,
                                                                                      mean = mean,
                                                                                      stdev = sd,
                                                                                      q25 = ~quantile(., 0.25),
                                                                                      q75 = ~quantile(., 0.75),
                                                                                      max = max))) %>%
  pivot_longer(everything(), names_sep='_', names_to=c('variable', '.value'))
