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

ggplot(data = CytometryFile, aes(x=Date_num)) +
  geom_line(aes(y=SYNC, color = "Synechococcus")) +
  geom_point(aes(y=SYNC, color = "Synechococcus")) +
  geom_line(aes(y=PROC, color = "Prochlorococcus")) +
  geom_point(aes(y=PROC, color = "Prochlorococcus"))

# Check if the scale is right

ggplot(data = CytometryFile, aes(x=Date_num)) +
  geom_line(aes(y=PROC, color = "Prochlorococcus")) +
  geom_point(aes(y=PROC, color = "Prochlorococcus")) +
  geom_line(aes(y=SYNC, color = "Synechococcus")) +
  geom_point(aes(y=SYNC, color = "Synechococcus")) 

# Add vertical lines depending on the month

ggplot(data = CytometryFile, aes(x=Date_num)) +
  geom_line(aes(y=SYNC, color = "Synechococcus")) +
  geom_point(aes(y=SYNC, color = "Synechococcus")) +
  geom_line(aes(y=PROC, color = "Prochlorococcus")) +
  geom_point(aes(y=PROC, color = "Prochlorococcus")) +
  geom_vline(xintercept = as.numeric(as.Date("2013-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2013-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2014-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2014-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2015-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2015-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2016-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2016-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2017-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2017-08-31")), linetype=4)

# Add picoeukaryotes and nanoeukaryotes, and cryptophytes

ggplot(data = CytometryFile, aes(x=Date_num)) +
  geom_line(aes(y=SYNC, color = "Synechococcus")) +
  geom_point(aes(y=SYNC, color = "Synechococcus")) +
  geom_line(aes(y=PROC, color = "Prochlorococcus")) +
  geom_point(aes(y=PROC, color = "Prochlorococcus")) +
  geom_line(aes(y=PICOEC, color = "Pico eukaryotes")) +
  geom_point(aes(y=PICOEC, color = "Pico eukaryotes")) +
  geom_line(aes(y=NANOEC, color = "Nano eukaryotes")) +
  geom_point(aes(y=NANOEC, color = "Nano eukaryotes")) +
  geom_line(aes(y=CRYC, color = "Cryptophytes")) +
  geom_point(aes(y=CRYC, color = "Cryptophytes")) +
  scale_color_manual(values=c("black", "purple", "green", "red", "blue")) +
  geom_vline(xintercept = as.numeric(as.Date("2013-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2013-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2014-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2014-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2015-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2015-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2016-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2016-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2017-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2017-08-31")), linetype=4)

# Plot fluorescence of synechococcus and synechococcus abundance

CytometryFile$SYNC_scaled <- CytometryFile$SYNC / 10000

ggplot(data = CytometryFile, aes(x=Date_num)) +
  geom_line(aes(y=SYNC_scaled, color = "Synechococcus")) +
  geom_point(aes(y=SYNC_scaled, color = "Synechococcus")) +
  geom_line(aes(y=SYNFLR, color = "Fluo red chlorophylle")) +
  geom_point(aes(y=SYNFLR, color = "Fluo red chlorophylle")) +
  geom_line(aes(y=SYNFLO, color = "Fluo orange phycoerythrin")) +
  geom_point(aes(y=SYNFLO, color = "Fluo orange phycoerythrin")) +
  scale_color_manual(values=c("orange", "red", "green")) +
  geom_vline(xintercept = as.numeric(as.Date("2013-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2013-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2014-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2014-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2015-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2015-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2016-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2016-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2017-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2017-08-31")), linetype=4)

# Plot fluorescences of all groups

ggplot(data = CytometryFile, aes(x=Date_num)) +
  geom_line(aes(y=SYNFLR, color = "Synechococcus fluo red")) +
  geom_point(aes(y=SYNFLR, color = "Synechococcus fluo red")) +
  geom_line(aes(y=PROFLR, color = "Prochlorococcus fluo red")) +
  geom_point(aes(y=PROFLR, color = "Prochlorococcus fluo red")) +
  geom_line(aes(y=PICOEFLR, color = "Pico eukaryotes fluo red")) +
  geom_point(aes(y=PICOEFLR, color = "Pico eukaryotes fluo red")) +
  geom_line(aes(y=NANOEFLR, color = "Nano eukaryotes fluo red")) +
  geom_point(aes(y=NANOEFLR, color = "Nano eukaryotes fluo red")) +
  geom_line(aes(y=CRYFLR, color = "Cryptophytes fluo red")) +
  geom_point(aes(y=CRYFLR, color = "Cryptophytes fluo red")) +
  scale_color_manual(values=c("black", "purple", "green", "red", "blue")) +
  geom_vline(xintercept = as.numeric(as.Date("2013-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2013-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2014-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2014-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2015-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2015-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2016-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2016-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2017-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2017-08-31")), linetype=4) +
  ylim(0, 300)

# # Multiply the abundance of synechococcus by the fluorescence of its chlorophyl
# # Do the same for prochlorococcus, picoeukaryotes and nano eukaryotes
# 
# CytometryFile$SYNC_NormChloro <- CytometryFile$SYNC * CytometryFile$SYNFLR
# CytometryFile$PROC_NormChloro <- CytometryFile$PROC * CytometryFile$PROFLR
# CytometryFile$PICOEC_NormChloro <- CytometryFile$PICOEC * CytometryFile$PICOEFLR
# CytometryFile$NANOEC_NormChloro <- CytometryFile$NANOEFLR * CytometryFile$NANOEFLR
# CytometryFile$CRYC_NormChloro <- CytometryFile$CRYC * CytometryFile$CRYFLR
# 
# 
# ggplot(data = CytometryFile, aes(x=Date_num)) +
#   geom_line(aes(y=SYNC_NormChloro, color = "Synechococcus norm chlorophyl")) +
#   geom_point(aes(y=SYNC_NormChloro, color = "Synechococcus norm chlorophyl")) +
#   geom_line(aes(y=PROC_NormChloro, color = "Prochlorococcus norm chlorophyl")) +
#   geom_point(aes(y=PROC_NormChloro, color = "Prochlorococcus norm chlorophyl")) +
#   geom_line(aes(y=PICOEC_NormChloro, color = "Picoeukaryotes norm chlorophyl")) +
#   geom_point(aes(y=PICOEC_NormChloro, color = "Picoeukaryotes norm chlorophyl")) +
#   geom_line(aes(y=NANOEC_NormChloro, color = "Nanoeukaryotes norm chlorophyl")) +
#   geom_point(aes(y=NANOEC_NormChloro, color = "Nanoeukaryotes norm chlorophyl")) +
#   geom_line(aes(y=CRYC_NormChloro, color = "Cryptophytes norm chlorophyl")) +
#   geom_point(aes(y=CRYC_NormChloro, color = "Cryptophytes norm chlorophyl")) +
#   scale_color_manual(values=c("black", "purple", "green", "red", "blue")) +
#   geom_vline(xintercept = as.numeric(as.Date("2013-08-01")), linetype=4) +
#   geom_vline(xintercept = as.numeric(as.Date("2013-08-31")), linetype=4) +
#   geom_vline(xintercept = as.numeric(as.Date("2014-08-01")), linetype=4) +
#   geom_vline(xintercept = as.numeric(as.Date("2014-08-31")), linetype=4) +
#   geom_vline(xintercept = as.numeric(as.Date("2015-08-01")), linetype=4) +
#   geom_vline(xintercept = as.numeric(as.Date("2015-08-31")), linetype=4) +
#   geom_vline(xintercept = as.numeric(as.Date("2016-08-01")), linetype=4) +
#   geom_vline(xintercept = as.numeric(as.Date("2016-08-31")), linetype=4) +
#   geom_vline(xintercept = as.numeric(as.Date("2017-08-01")), linetype=4) +
#   geom_vline(xintercept = as.numeric(as.Date("2017-08-31")), linetype=4)

# Divide  the fluorescence of its chlorophyl by the abundance of synechococcus
# Do the same for prochlorococcus, picoeukaryotes and nano eukaryotes

CytometryFile$SYNC_NormChloro <- CytometryFile$SYNFLR / CytometryFile$SYNC 
CytometryFile$SYNC_NormChloro[which(is.infinite(CytometryFile$SYNC_NormChloro))] <- NA

CytometryFile$PROC_NormChloro <- CytometryFile$PROFLR / CytometryFile$PROC 
CytometryFile$PROC_NormChloro[which(is.infinite(CytometryFile$PROC_NormChloro))] <- NA

CytometryFile$PICOEC_NormChloro <- CytometryFile$PICOEFLR / CytometryFile$PICOEC 
CytometryFile$PICOEC_NormChloro[which(is.infinite(CytometryFile$PICOEC_NormChloro))] <- NA

CytometryFile$NANOEC_NormChloro <- CytometryFile$NANOEFLR / CytometryFile$NANOEFLR
CytometryFile$NANOEC_NormChloro[which(is.infinite(CytometryFile$NANOEC_NormChloro))] <- NA

CytometryFile$CRYC_NormChloro <- CytometryFile$CRYFLR / CytometryFile$CRYC
CytometryFile$CRYC_NormChloro[which(is.infinite(CytometryFile$CRYC_NormChloro))] <- NA

ggplot(data = CytometryFile, aes(x=Date_num)) +
  geom_line(aes(y=SYNC_NormChloro, color = "Synechococcus norm chlorophyl")) +
  geom_point(aes(y=SYNC_NormChloro, color = "Synechococcus norm chlorophyl")) +
  geom_line(aes(y=PROC_NormChloro, color = "Prochlorococcus norm chlorophyl")) +
  geom_point(aes(y=PROC_NormChloro, color = "Prochlorococcus norm chlorophyl")) +
  # geom_line(aes(y=PICOEC_NormChloro, color = "Picoeukaryotes norm chlorophyl")) +
  # geom_point(aes(y=PICOEC_NormChloro, color = "Picoeukaryotes norm chlorophyl")) +
  # geom_line(aes(y=NANOEC_NormChloro, color = "Nanoeukaryotes norm chlorophyl")) +
  # geom_point(aes(y=NANOEC_NormChloro, color = "Nanoeukaryotes norm chlorophyl")) +
  # geom_line(aes(y=CRYC_NormChloro, color = "Cryptophytes norm chlorophyl")) +
  # geom_point(aes(y=CRYC_NormChloro, color = "Cryptophytes norm chlorophyl")) +
  scale_color_manual(values=c("black", "purple", "green", "red", "blue")) +
  geom_vline(xintercept = as.numeric(as.Date("2013-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2013-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2014-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2014-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2015-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2015-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2016-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2016-08-31")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2017-08-01")), linetype=4) +
  geom_vline(xintercept = as.numeric(as.Date("2017-08-31")), linetype=4)

# Correlation between phycoerythrin and chlorophyll in Synechococcus, and cell concentration

cor(CytometryFile$SYNFLR, CytometryFile$SYNFLO)
cor(CytometryFile$SYNC_scaled, CytometryFile$SYNFLO)
cor(CytometryFile$SYNFLR, CytometryFile$SYNC_scaled)


#----------------------------------------------------------------------------------#
#          Add 0.2 um size fraction Prokaryotes from Lambert et al. 2021           #
#----------------------------------------------------------------------------------#

# Import prokaryotes dataset

Prok02_OTU_taxo_Seq <- read.csv("/home/emilelaymand/Documents/Science/These/SOLA_Galand/SOLA_project_Emile_Laymand/SOLA/wetransfer_sola_2021-10-18_1242/ProkaryotesLambertEtAl2021_2015_2017.csv", sep = "\t")

Prok02_OTU_taxo_Seq$ASV <- paste("Prok_ASV", c(1:length(row.names(Prok02_OTU_taxo_Seq))), sep = "")
row.names(Prok02_OTU_taxo_Seq) <- Prok02_OTU_taxo_Seq$ASV

# Separate the dataframe in three, according if it is the sequence, the taxonomy or the abundance

Prok02_OTU_tax <- Prok02_OTU_taxo_Seq[, c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")]

Prok02_OTU_sequences <- Prok02_OTU_taxo_Seq[, c("ASV", "Sequence")]

Prok02_OTU_abund <- Prok02_OTU_taxo_Seq[, !(colnames(Prok02_OTU_taxo_Seq) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Sequence"))]
