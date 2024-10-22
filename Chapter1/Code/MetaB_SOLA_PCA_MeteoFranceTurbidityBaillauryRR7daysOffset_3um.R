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
library(plotrix)
library(ggbreak)
library(ggforce)
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

#### Added in this script, not present in MetaB_SOLA_PCA.R ####

# Prepare to add the relative abundance of the 40 most abundant ASVs

Phyloseq_ForOutstandingASVs <- OTUPhyloseq_NormDomain
Phyloseq_ForOutstandingASVs <- subset_taxa(Phyloseq_ForOutstandingASVs, Division == "Fungi")

topN <- 40
most_abundant_taxa <- sort(taxa_sums(Phyloseq_ForOutstandingASVs), TRUE)[1:topN]
Phyloseq_ForOutstandingASVs39ASV <- prune_taxa(names(most_abundant_taxa), Phyloseq_ForOutstandingASVs)

Phyloseq_ForOutstandingASVs39ASV_t <- data.frame(t(otu_table(Phyloseq_ForOutstandingASVs39ASV)))

# Check the order of columns 

identical(sample_data(Phyloseq_EukAll)[["Name_Euk"]], row.names(Phyloseq_ForOutstandingASVs39ASV_t)) # It works

sample_data(Phyloseq_EukAll) <- cbind(sample_data(Phyloseq_EukAll), Phyloseq_ForOutstandingASVs39ASV_t) # /!\ To be checked

#### End of added in this script, not present in MetaB_SOLA_PCA.R ####

#=============================================================================================#
# Import the metadata from Meteo France, and the Baillaury level, and the flux cytometry data #
#=============================================================================================#

#=============================================================================================
# From MetaB_SOLA_PlusTurbidityMeteoFranceBaillaury3__RR7daysOffset.R

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
#ggplot(BaillauryHeigth, aes(x=Date_num, y=BaillauryHeigthIn_m)) +
#  geom_point(size=0.7) +
#  geom_line()

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

# End of section imported from MetaB_SOLA_PlusTurbidityMeteoFranceBaillaury3__RR7daysOffset.R
#====================================================================================================

# Add the Meteo France and the Baillaury data to the sample_data of Phyloseq_EukAll

# Take only the dates that are present in the metabarcoding dataset

DatesToKeep <- sample_data(Phyloseq_EukAll)[["Date_num"]]

MeteoFranceMetadata_OnlyDates <- MeteoFranceMetadata[MeteoFranceMetadata[["Date_num"]] %in% DatesToKeep,]
dim(MeteoFranceMetadata_OnlyDates)

TurbidityMetadata_OnlyDates <- TurbidityMetadata[TurbidityMetadata[["Date_num"]] %in% DatesToKeep,]
dim(TurbidityMetadata_OnlyDates)

BaillauryHeigth_OnlyDates <- BaillauryHeigth[BaillauryHeigth[["Date_num"]] %in% DatesToKeep,]
dim(BaillauryHeigth_OnlyDates)

# Append the dataframes to phyloseq_EukAll by adding columns full of NAs, then filling the columns with the values from the corresponding dates

sample_data(Phyloseq_EukAll)[["RR"]] <- NA
sample_data(Phyloseq_EukAll)[["FFM"]] <- NA
sample_data(Phyloseq_EukAll)[["RR7days"]] <- NA
sample_data(Phyloseq_EukAll)[["Turbidity_3m"]] <- NA
sample_data(Phyloseq_EukAll)[["Turbidity_20m"]] <- NA
sample_data(Phyloseq_EukAll)[["BaillauryHeigthIn_m"]] <- NA
sample_data(Phyloseq_EukAll)[["DerivHeigth_m_j"]] <- NA
sample_data(Phyloseq_EukAll)[["BaillauryHeigthIn_m_x10"]] <- NA
sample_data(Phyloseq_EukAll)[["DerivHeigth_m_j_x10"]] <- NA

# For Meteo France and the Baillaury # /!\ Hard Check --> Checked: looks like it works

for (i in (1:nrow(sample_data(Phyloseq_EukAll)))){
  sample_data(Phyloseq_EukAll)[i,"RR"] <- MeteoFranceMetadata_OnlyDates[which(MeteoFranceMetadata_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "RR"]
  
  sample_data(Phyloseq_EukAll)[i,"FFM"] <- MeteoFranceMetadata_OnlyDates[which(MeteoFranceMetadata_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "FFM"]
  
  sample_data(Phyloseq_EukAll)[i,"RR7days"] <- MeteoFranceMetadata_OnlyDates[which(MeteoFranceMetadata_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "RR7days"]
  
  sample_data(Phyloseq_EukAll)[i,"BaillauryHeigthIn_m"] <- BaillauryHeigth_OnlyDates[which(BaillauryHeigth_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "BaillauryHeigthIn_m"]
  
  sample_data(Phyloseq_EukAll)[i,"DerivHeigth_m_j"] <- BaillauryHeigth_OnlyDates[which(BaillauryHeigth_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "DerivHeigth_m_j"]
  
  sample_data(Phyloseq_EukAll)[i,"BaillauryHeigthIn_m_x10"] <- BaillauryHeigth_OnlyDates[which(BaillauryHeigth_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "BaillauryHeigthIn_m_x10"]
  
  sample_data(Phyloseq_EukAll)[i,"DerivHeigth_m_j_x10"] <- BaillauryHeigth_OnlyDates[which(BaillauryHeigth_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "DerivHeigth_m_j_x10"]
}

# For Turbidity # /!\ Hard Check --> Checked: looks like it works

for (i in which(sample_data(Phyloseq_EukAll)[["Date_num"]] %in% TurbidityMetadata[["Date_num"]])){
  sample_data(Phyloseq_EukAll)[i,"Turbidity_3m"] <- TurbidityMetadata_OnlyDates[which(TurbidityMetadata_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "Turbidity_3m"]
  
  sample_data(Phyloseq_EukAll)[i,"Turbidity_20m"] <- TurbidityMetadata_OnlyDates[which(TurbidityMetadata_OnlyDates[["Date_num"]] == as.vector(sample_data(Phyloseq_EukAll)[i,"Date_num"])) , "Turbidity_20m"]
}

#======================================================#
#                      Make the PCA                    #
#======================================================#

DF_SampleDataForPCA <- data.frame(sample_data(Phyloseq_EukAll))

# 1- Remove all columns that are not useful for PCA

DF_SampleDataForPCA_rmCol <- DF_SampleDataForPCA[,!(colnames(DF_SampleDataForPCA) %in% c("Name_Euk", "BP_Euk", "Date_Euk", "Step", "DATE", "ANR", "Month", "Year", "Above1", "Above5", "Above10", "Above15", "Above20", "Above30", "Above40", "FungiRange", "SumPercClasses", "Turbidity_3m", "Turbidity_20m", "BaillauryHeigthIn_m_x10", "DerivHeigth_m_j_x10"))]

# 2- Extrapolate the missing values by using number of days

# Order the dataframe by date # /!\ Hard check

DF_SampleDataForPCA_rmCol_sortedDate <- DF_SampleDataForPCA_rmCol[order(DF_SampleDataForPCA_rmCol$Date_num), ]

# Remove samples from 0.2 um size fraction

DF_SampleDataForPCA_rmCol_sortedDate_3um <- DF_SampleDataForPCA_rmCol_sortedDate[DF_SampleDataForPCA_rmCol_sortedDate[["Filter_Euk"]] == "3µM",]

# Add a column number of days before and number of days after 

DF_SampleDataForPCA_rmCol_sortedDate_3um$DaysBefore <- NA
for (i in (2:nrow(DF_SampleDataForPCA_rmCol_sortedDate_3um))){
  DF_SampleDataForPCA_rmCol_sortedDate_3um$DaysBefore[i] <- as.numeric(DF_SampleDataForPCA_rmCol_sortedDate_3um$Date_num[i]-DF_SampleDataForPCA_rmCol_sortedDate_3um$Date_num[i-1])
}

DF_SampleDataForPCA_rmCol_sortedDate_3um$DaysAfter <- NA
for (i in (1:(nrow(DF_SampleDataForPCA_rmCol_sortedDate_3um)-1))){
  DF_SampleDataForPCA_rmCol_sortedDate_3um$DaysAfter[i] <- as.numeric(DF_SampleDataForPCA_rmCol_sortedDate_3um$Date_num[i+1]-DF_SampleDataForPCA_rmCol_sortedDate_3um$Date_num[i])
}

# Extrapolate the values --> Values may be missing for T, S, O, pH, NH4, NO3, NO2, PO4, SIO4 and CHLA
# Extrapolated values is the mean of the closest value, ponderated by the temporal closeness to the closest actual points
# If two points are missing in a row, the value is extrapolated with the value two points away, and so on. 

# FOR T AND S: replace dubious 0's by NAs # Check there is no dubious 0's in other variables

DF_SampleDataForPCA_rmCol_sortedDate_3um[which(DF_SampleDataForPCA_rmCol_sortedDate_3um[["T"]] == 0), "T"] <- NA

DF_SampleDataForPCA_rmCol_sortedDate_3um[which(DF_SampleDataForPCA_rmCol_sortedDate_3um[["S"]] == 0), "S"] <- NA

# Make a function that calculates the extrapolated values

extrapolate_values <- function(DF_ToExtrap, V_ValToExtrap){
  # This function extrapolates the missing values in DF_ToExtrap for all variables stored in V_ValToExtrap
  
  DF_SampleDataForPCA_rmCol_sortedDate_3um <- DF_ToExtrap
  
  for (Variable in V_ValToExtrap){
    
    VarToExtrap <- Variable
    print(VarToExtrap)
    
    IndexToExtrapolate <- which(is.na(DF_SampleDataForPCA_rmCol_sortedDate_3um[[VarToExtrap]]))
    VectExtrapolatedValues <- numeric(0)
    IndexWhereAdd <- numeric(0)
    
    # Calculate all the extrapolated values
    for (i in IndexToExtrapolate){
      if ((i == 1) | (i == length(DF_SampleDataForPCA_rmCol_sortedDate_3um[[VarToExtrap]]))) { # Do nothing, as a point with no point farther must not be extrapolated
      } else {
        NonNABefore <- 1
        NonNAAfter <- 1
        while (is.na(DF_SampleDataForPCA_rmCol_sortedDate_3um[[VarToExtrap]][i-NonNABefore])){
          NonNABefore <- NonNABefore +1
        }
        while (is.na(DF_SampleDataForPCA_rmCol_sortedDate_3um[[VarToExtrap]][i+NonNAAfter])){
          NonNAAfter <- NonNAAfter +1
        }
        NdaysBefore <- as.numeric(DF_SampleDataForPCA_rmCol_sortedDate_3um$Date_num[i]-DF_SampleDataForPCA_rmCol_sortedDate_3um$Date_num[i-NonNABefore])
        NdaysAfter <- as.numeric(DF_SampleDataForPCA_rmCol_sortedDate_3um$Date_num[i+NonNAAfter]-DF_SampleDataForPCA_rmCol_sortedDate_3um$Date_num[i])
        
        ExtrapVal <- ((1-((NdaysBefore)/(NdaysBefore+NdaysAfter))) * DF_SampleDataForPCA_rmCol_sortedDate_3um[[VarToExtrap]][i-NonNABefore]) + ((1-((NdaysAfter)/(NdaysBefore+NdaysAfter))) * DF_SampleDataForPCA_rmCol_sortedDate_3um[[VarToExtrap]][i+NonNAAfter])
        
        VectExtrapolatedValues <- c(VectExtrapolatedValues, ExtrapVal)
        IndexWhereAdd <- c(IndexWhereAdd, i)
        # print("New Line")
        # print(paste("NonNABefore:", NonNABefore))
        # print(paste("NonNAAfter:", NonNAAfter))
        # print(paste("ExtrapVal:",ExtrapVal))
        # print(paste("IndexWhereAdd:", i))
        # print(paste("NumberBefore", DF_SampleDataForPCA_rmCol_sortedDate_3um[[VarToExtrap]][i-NonNABefore]))
        # print(paste("NumberAfter:", DF_SampleDataForPCA_rmCol_sortedDate_3um[[VarToExtrap]][i+NonNAAfter]))
      }
    }
    
    # Add the extrapolated values to the original vector # /!\ Hard check on the results /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
    for (i in (1:length(VectExtrapolatedValues))){
      DF_SampleDataForPCA_rmCol_sortedDate_3um[[VarToExtrap]][IndexWhereAdd[i]] <- VectExtrapolatedValues[i]
    }
  }
  
  return(DF_SampleDataForPCA_rmCol_sortedDate_3um)
}

V_VarToExtrap <- c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA", "FFM")
DF_ToExtrap <- DF_SampleDataForPCA_rmCol_sortedDate_3um

DF_SampleDataForPCA_rmCol_sortedDate_3um_Extrap <- extrapolate_values(DF_ToExtrap, V_VarToExtrap)

# 3- Remove some columns that I will not use for the computation of the PCA

DF_SampleDataForPCA_rmCol_sortedDate_3um_Extrap_rmCol <- DF_SampleDataForPCA_rmCol_sortedDate_3um_Extrap[, !(colnames (DF_SampleDataForPCA_rmCol_sortedDate_3um_Extrap) %in% c("Filter_Euk", "Date_num", "DaysBefore", "DaysAfter"))]

# 4- Sort columns to have first all the columns used to compute the PCA, and then the supplementary variables
DF_SampleDataForPCA_rmCol_Reorder <- DF_SampleDataForPCA_rmCol_sortedDate_3um_Extrap_rmCol[,c(1:9,61:65,10:60)]

# 5- Remove rows that still contain NAs

dim(DF_SampleDataForPCA_rmCol_Reorder)
DF_SampleDataForPCA_rmCol_Reorder_Complete <- DF_SampleDataForPCA_rmCol_Reorder[complete.cases(DF_SampleDataForPCA_rmCol_Reorder),]
dim(DF_SampleDataForPCA_rmCol_Reorder_Complete)

#======================================================#
#                  Make the PCA itself                 #
#======================================================#

# 1- Center and reduce all variables

DF_3um_CentRed <- data.frame(lapply(DF_SampleDataForPCA_rmCol_Reorder_Complete, function(x) (x-mean(x))/sd(x) ))

# Check it worked 

lapply(DF_3um_CentRed, mean) # The mean gets close to zero but is not exactly 0. Check if this is an issue
lapply(DF_3um_CentRed, sd)

# 2- Remove columns that have NaN

DF_3um_CentRed$ASV971 <- NULL
DF_3um_CentRed$PercEntomophthoromycota <- NULL

# 3- Make the PCA

PCA_3umNice <- PCA(DF_3um_CentRed, scale.unit = FALSE, ncp = 10, graph = TRUE, quanti.sup = 15:63)
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

# Add the plot of the samples

# Percent Fungi
fviz_pca_ind(PCA_3umNice, label="none", col.ind = DF_SampleDataForPCA_rmCol_Reorder_Complete$ASV224,
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# Idem but remove Chla from the plot (to improve the color scale of the Cos2 of the supplementary variables)
#-----------------------------------------------------------------------------------------------------------

DF_3um_CentRed_woChla <- DF_3um_CentRed
DF_3um_CentRed_woChla$CHLA <- NULL

PCA_3umNice <- PCA(DF_3um_CentRed_woChla, scale.unit = FALSE, ncp = 10, graph = TRUE, quanti.sup = 15:62)
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

# Add the plot of the samples

# Percent Fungi
fviz_pca_ind(PCA_3umNice, label="none", col.ind = DF_SampleDataForPCA_rmCol_Reorder_Complete$ASV224,
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# Idem but remove the derivate of Baillaury Heigth and RR7days
#-------------------------------------------------------------

DF_3um_CentRed_woRR7daysDerBaill <- DF_3um_CentRed
DF_3um_CentRed_woRR7daysDerBaill$DerivHeigth_m_j <- NULL
DF_3um_CentRed_woRR7daysDerBaill$RR7days <- NULL

PCA_3umNice <- PCA(DF_3um_CentRed_woRR7daysDerBaill, scale.unit = FALSE, ncp = 10, graph = TRUE, quanti.sup = 13:61)
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

# Add the plot of the samples

# Percent Fungi
fviz_pca_ind(PCA_3umNice, label="none", col.ind = DF_SampleDataForPCA_rmCol_Reorder_Complete$ASV224,
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# Idem, but remove the derivate of Baillaury Heigth and RR7days, and keep only the 10 first ASVs as supplementary variables

DF_3um_CentRed_woRR7daysDerBaill_10ASVs <- DF_3um_CentRed
DF_3um_CentRed_woRR7daysDerBaill_10ASVs$DerivHeigth_m_j <- NULL
DF_3um_CentRed_woRR7daysDerBaill_10ASVs$RR7days <- NULL

DF_3um_CentRed_woRR7daysDerBaill_10ASVs <- DF_3um_CentRed_woRR7daysDerBaill_10ASVs[,colnames(DF_3um_CentRed_woRR7daysDerBaill_10ASVs) %in% c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "RR", "FFM",  "BaillauryHeigthIn_m", "ASV224", "ASV275", "ASV1613", "ASV1016", "ASV409", "ASV166", "ASV594", "ASV847", "ASV367", "ASV553","CHLA")]

PCA_3umNice <- PCA(DF_3um_CentRed_woRR7daysDerBaill_10ASVs, scale.unit = FALSE, ncp = 10, graph = TRUE, quanti.sup = 13:22)
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

# Add the plot of the samples

# Percent Fungi
fviz_pca_ind(PCA_3umNice, label="none", col.ind = DF_SampleDataForPCA_rmCol_Reorder_Complete$ASV224,
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, label="none", col.ind = DF_SampleDataForPCA_rmCol_Reorder_Complete$ASV594,
             legend.title = "ASV594", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

#=====================================#
#   Test missMDA to extrapolate data  #
#=====================================#

# 1- Extrapolate the values according to the values projected on the first component as returned by estimncpPCA
#--------------------------------------------------------------------------------------------------------------

DF_3um_ToExtrap_missMDA <- DF_SampleDataForPCA_rmCol_sortedDate_3um[,colnames(DF_SampleDataForPCA_rmCol_sortedDate_3um) %in% c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "RR", "FFM",  "BaillauryHeigthIn_m", "CHLA")]
nb <- estim_ncpPCA(DF_3um_ToExtrap_missMDA,ncp.max=5) # Output is 1
res.comp_3um_ToExtrap = imputePCA(DF_3um_ToExtrap_missMDA,ncp=1)
Extrapolated_DF <- data.frame(res.comp_3um_ToExtrap$completeObs)

# 2-Replace the columns in the original dataset with the extrapolated columns
#----------------------------------------------------------------------------

# Check row names are in the same order
identical(row.names(DF_SampleDataForPCA_rmCol_sortedDate_3um), row.names(Extrapolated_DF))

# Replace the original columns with the complete ones 

DF_3um_NonExtrapCols <- DF_SampleDataForPCA_rmCol_sortedDate_3um[,!(colnames(DF_SampleDataForPCA_rmCol_sortedDate_3um) %in% c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "RR", "FFM",  "BaillauryHeigthIn_m", "CHLA"))]

FullExtrapDF <- cbind(Extrapolated_DF, DF_3um_NonExtrapCols) # /!\ Hard check

# 3- Center and reduce the variables
#-----------------------------------

# 0- Remove columns that will not be used for the PCA 

FullExtrapDF$Filter_Euk <- NULL
FullExtrapDF$Date_num <- NULL 
FullExtrapDF$DaysBefore <- NULL
FullExtrapDF$DaysAfter <- NULL

# 1- Center and reduce all variables

DF_3um_CentRed <- data.frame(lapply(FullExtrapDF, function(x) (x-mean(x))/sd(x) ))

# Check it worked 

lapply(DF_3um_CentRed, mean) # The mean gets close to zero but is not exactly 0. Check if this is an issue
lapply(DF_3um_CentRed, sd)

# 2- Remove columns that have NaN

DF_3um_CentRed$ASV971 <- NULL
DF_3um_CentRed$PercEntomophthoromycota <- NULL

# 4- Compute the PCA
#-------------------

FullExtrapDF_PlusLog <- FullExtrapDF
FullExtrapDF_PlusLog$LogPercFungi <- log(FullExtrapDF_PlusLog$PercFungi)

# Idem, but remove the derivate of Baillaury Heigth and RR7days, and keep only the 10 first ASVs as supplementary variables
#==========================================================================================================================

DF_3um_CentRed_woRR7daysDerBaill_10ASVs <- DF_3um_CentRed
DF_3um_CentRed_woRR7daysDerBaill_10ASVs$DerivHeigth_m_j <- NULL
DF_3um_CentRed_woRR7daysDerBaill_10ASVs$RR7days <- NULL

DF_3um_CentRed_woRR7daysDerBaill_10ASVs <- DF_3um_CentRed_woRR7daysDerBaill_10ASVs[,c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "RR", "FFM",  "BaillauryHeigthIn_m", "ASV224", "ASV275", "ASV1613", "ASV1016", "ASV409", "ASV166", "ASV594", "ASV847", "ASV367", "ASV553","CHLA")]

PCA_3umNice <- PCA(DF_3um_CentRed_woRR7daysDerBaill_10ASVs, scale.unit = FALSE, ncp = 10, graph = TRUE, quanti.sup = 13:23)
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

# Add the plot of the samples

# Percent Fungi
fviz_pca_ind(PCA_3umNice, label="none", col.ind = FullExtrapDF$ASV224,
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, label="none", col.ind = FullExtrapDF$ASV594,
             legend.title = "ASV594", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_PlusLog$PercFungi, #col.var = "contrib",
                legend.title = "Percent Fungi", 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                #gradient.cols = c("white", "blue", "red"),
                #gradient.cols = palette("rainbow"),
                repel = TRUE, geom.ind = c("point", "text"))# +
                #scale_color_continuous(breaks = c(0, 5, 50))

PCA_Figure3 <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_PlusLog$ASV594, #col.var = "contrib",
                               legend.title = "ASV594", 
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                               #gradient.cols = c("white", "blue", "red"),
                               #gradient.cols = palette("rainbow"),
                               repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))


ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/PCA_Figure3.pdf", plot=PCA_Figure3, device = cairo_pdf(), width=15, height=10)
dev.off()

# Same, but remove the derivate of Baillaury Heigth and RR7days, and keep all ASVs 
#=================================================================================

DF_3um_CentRed_woRR7daysDerBaill_10ASVs <- DF_3um_CentRed
DF_3um_CentRed_woRR7daysDerBaill_10ASVs$DerivHeigth_m_j <- NULL
DF_3um_CentRed_woRR7daysDerBaill_10ASVs$RR7days <- NULL

DF_3um_CentRed_woRR7daysDerBaill_10ASVs_Reordered <- DF_3um_CentRed_woRR7daysDerBaill_10ASVs[, c(1:9, 11:61, 10)]

PCA_3umNice <- PCA(DF_3um_CentRed_woRR7daysDerBaill_10ASVs_Reordered, scale.unit = FALSE, ncp = 10, graph = TRUE, quanti.sup = 13:61)
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
corrplot(PCA_3umNice$quanti.sup$cor, is.corr=FALSE)


# Add the plot of the samples

# Percent Fungi
fviz_pca_ind(PCA_3umNice, label="none", col.ind = FullExtrapDF$ASV224,
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, label="none", col.ind = FullExtrapDF$ASV594,
             legend.title = "ASV594", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))


# Same, but remove the derivate of Baillaury Heigth and RR7days and CHLA, and keep all ASVs
#=================================================================================

DF_3um_CentRed_woRR7daysDerBaill_10ASVs <- DF_3um_CentRed
DF_3um_CentRed_woRR7daysDerBaill_10ASVs$DerivHeigth_m_j <- NULL
DF_3um_CentRed_woRR7daysDerBaill_10ASVs$RR7days <- NULL

DF_3um_CentRed_woRR7daysDerBaill_10ASVs_Reordered <- DF_3um_CentRed_woRR7daysDerBaill_10ASVs[, c(1:9, 11:61)]

PCA_3umNice <- PCA(DF_3um_CentRed_woRR7daysDerBaill_10ASVs_Reordered, scale.unit = FALSE, ncp = 10, graph = TRUE, quanti.sup = 13:60)
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
corrplot(PCA_3umNice$quanti.sup$cor, is.corr=FALSE)


# Add the plot of the samples

# Percent Fungi
fviz_pca_ind(PCA_3umNice, label="none", col.ind = FullExtrapDF$ASV224,
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, label="none", col.ind = FullExtrapDF$ASV594,
             legend.title = "ASV594", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

#############################################################
#                  Add the flux cytometry                   #
#############################################################

# Import the file

CytometryFile <- read.csv("Cytometry/Somlit_Extraction_piconano_2013_2017_modified.csv", sep = ";")

# Add a numeric date column 
CytometryFile$Date_num <- as.Date(CytometryFile$DATE, "%Y-%m-%d")

# Keep only SYNC, PROC, PICOEC, NANOEC, CRYC as columns

CytometryFileTrunc <- CytometryFile[,c("Date_num", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC")]

# Merge with the big dataframe

DF_3um_ForMergingCyto <- DF_SampleDataForPCA_rmCol_sortedDate_3um

DF_3um_withcyto <- merge(DF_3um_ForMergingCyto, CytometryFileTrunc, by = "Date_num", all.x = TRUE) # Check this merging

# Extrapolate the values with missMDA

DF_3um_withcyto_missMDA <- DF_3um_withcyto[,colnames(DF_3um_withcyto) %in% c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "RR", "FFM",  "BaillauryHeigthIn_m", "CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC")]
nb <- estim_ncpPCA(DF_3um_withcyto_missMDA,ncp.max=5) # Output is 1
res.comp_DF_3um_withcyto_missMDA = imputePCA(DF_3um_withcyto_missMDA,ncp=1)
Extrapolated_DF_withCyto <- data.frame(res.comp_DF_3um_withcyto_missMDA$completeObs)

# 2-Replace the columns in the original dataset with the extrapolated columns
#----------------------------------------------------------------------------

# Check row names are in the same order
identical(row.names(DF_3um_withcyto), row.names(Extrapolated_DF_withCyto))

# Replace the original columns with the complete ones 

DF_3um_withcyto_NonExtrapCols <- DF_3um_withcyto[,!(colnames(DF_3um_withcyto) %in% c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "RR", "FFM",  "BaillauryHeigthIn_m", "CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC"))]

FullExtrapDF_withCyto <- cbind(Extrapolated_DF_withCyto, DF_3um_withcyto_NonExtrapCols) # /!\ Hard check




# 3- Center and reduce the variables
#-----------------------------------

# 0- Remove columns that will not be used for the PCA 

FullExtrapDF_withCyto$Filter_Euk <- NULL
FullExtrapDF_withCyto$Date_num <- NULL 
FullExtrapDF_withCyto$DaysBefore <- NULL
FullExtrapDF_withCyto$DaysAfter <- NULL

# 1- Center and reduce all variables

FullExtrapDF_withCyto_3um_CentRed <- data.frame(lapply(FullExtrapDF_withCyto, function(x) (x-mean(x))/sd(x) ))

# Check it worked 

lapply(FullExtrapDF_withCyto_3um_CentRed, mean) # The mean gets close to zero but is not exactly 0. Check if this is an issue
lapply(FullExtrapDF_withCyto_3um_CentRed, sd)

# 2- Remove columns that have NaN

FullExtrapDF_withCyto_3um_CentRed$ASV971 <- NULL
FullExtrapDF_withCyto_3um_CentRed$PercEntomophthoromycota <- NULL


#=================================================#
#                PCAs with Flow cytometry         #
#=================================================#

DF_3um_CentRed_woRR7daysDerBaill_10ASVs <- FullExtrapDF_withCyto_3um_CentRed
DF_3um_CentRed_woRR7daysDerBaill_10ASVs$DerivHeigth_m_j <- NULL
#DF_3um_CentRed_woRR7daysDerBaill_10ASVs$RR7days <- NULL

DF_3um_CentRed_woRR7daysDerBaill_10ASVs <- DF_3um_CentRed_woRR7daysDerBaill_10ASVs[,c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "RR", "FFM", "RR7days", "BaillauryHeigthIn_m", "ASV224", "ASV275", "ASV1613", "ASV1016", "ASV409", "ASV166", "ASV594", "ASV847", "ASV367", "ASV553","CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC")]

PCA_3umNice <- PCA(DF_3um_CentRed_woRR7daysDerBaill_10ASVs, scale.unit = FALSE, ncp = 10, graph = TRUE, quanti.sup = 14:29)
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

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_3um_PrincipalCos2.pdf",         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")
corrplot(PCA_3umNice.var$cos2, is.corr=FALSE)
dev.off()


pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_3um_PrincipalCor.pdf",         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")
corrplot(PCA_3umNice.var$cor, is.corr=FALSE)
dev.off()

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_3um_SupplementaryCos2.pdf",         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")
corrplot(PCA_3umNice$quanti.sup$cos2, is.corr=FALSE)
dev.off()

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_3um_SupplementaryCor.pdf",         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")
corrplot(PCA_3umNice$quanti.sup$cor, is.corr=FALSE)
dev.off()

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_3um_VarianceAxes.pdf",         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")
fviz_eig(PCA_3umNice, addlabels = TRUE, ylim = c(0, 40))
dev.off()


# Add the plot of the samples

# Percent Fungi
fviz_pca_ind(PCA_3umNice, label="none", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV224,
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, axes = c(2,3),label="none", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV224, # Axes 2,3
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, axes = c(3,4),label="none", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV224, # Axes 2,3
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, axes = c(4,5),label="none", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV224, # Axes 2,3
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, axes = c(5,6),label="none", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV224, # Axes 2,3
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, axes = c(6,7),label="none", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV224, # Axes 2,3
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, axes = c(7,8),label="none", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV224, # Axes 2,3
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, axes = c(8,9),label="none", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV224, # Axes 2,3
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, axes = c(9,10),label="none", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV224, # Axes 2,3
             legend.title = "ASV224", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

fviz_pca_ind(PCA_3umNice, label="none", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV594,
             legend.title = "ASV594", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

PCA_Figure3_CytoFlux <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$PercFungi, #col.var = "contrib",
                               legend.title = "Percentage Fungi\n(centered-scaled)", 
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                               #gradient.cols = c("white", "blue", "red"),
                               #gradient.cols = palette("rainbow"),
                               repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/PCA_Figure3_CytoFlux.pdf", plot=PCA_Figure3_CytoFlux, device = cairo_pdf(), width=15, height=10)
dev.off()

PCA_Figure3_CytoFlux_ASV224 <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV224, #col.var = "contrib",
                                        legend.title = "Percent ASV224 \n(centered-scaled)", 
                                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                        #gradient.cols = c("white", "blue", "red"),
                                        #gradient.cols = palette("rainbow"),
                                        repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV224.pdf", plot=PCA_Figure3_CytoFlux_ASV224, device = cairo_pdf(), width=15, height=10)
dev.off()

# Same plot but with percent of Fungi as the color, then all phyla
#-------------------------------------------------

PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$PercFungi, #col.var = "contrib",
                                               legend.title = "Percent Fungi \n(centered-scaled)", 
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               #gradient.cols = c("white", "blue", "red"),
                                               #gradient.cols = palette("rainbow"),
                                               repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorPercFungi.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# Same plot, but with the first 10 ASVs as colors
#------------------------------------------------

# ASV224 
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV224, #col.var = "contrib",
                                                  legend.title = "Percent ASV224 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV224.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV275 
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV275, #col.var = "contrib",
                                                  legend.title = "Percent ASV275 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV275.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV1613
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV1613, #col.var = "contrib",
                                                  legend.title = "Percent ASV1613 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV1613.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV1016
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV1016, #col.var = "contrib",
                                                  legend.title = "Percent ASV1016 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV1016.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV409
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV409, #col.var = "contrib",
                                                  legend.title = "Percent ASV409 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV409.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV166
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV166, #col.var = "contrib",
                                                  legend.title = "Percent ASV166 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV166.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV594
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV594, #col.var = "contrib",
                                                  legend.title = "Percent ASV594 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV594.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV847
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV847, #col.var = "contrib",
                                                  legend.title = "Percent ASV847 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV847.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV367
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV367, #col.var = "contrib",
                                                  legend.title = "Percent ASV367 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV367.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV553
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV553, #col.var = "contrib",
                                                  legend.title = "Percent ASV553 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV553.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV774
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV774, #col.var = "contrib",
                                                  legend.title = "Percent ASV774 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV774.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV843
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV843, #col.var = "contrib",
                                                  legend.title = "Percent ASV843 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV843.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# ASV420
PCA_Figure3_CytoFlux_PercFungi <- fviz_pca_biplot(PCA_3umNice, label="var", col.ind = FullExtrapDF_withCyto_3um_CentRed$ASV420, #col.var = "contrib",
                                                  legend.title = "Percent ASV420 \n(centered-scaled)", 
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  #gradient.cols = c("white", "blue", "red"),
                                                  #gradient.cols = palette("rainbow"),
                                                  repel = TRUE, geom.ind = c("point", "text"))# +
#scale_color_continuous(breaks = c(0, 5, 50))
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/PCA_Figure3_CytoFlux_ColorASV420.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

#==================================================================================================================#
# Added on 20240504, PCA with only the variables I will use for the final PCA                                      #
# Note: I just copied the code without modifying the variable names, so go back to the upper sections with caution #
#==================================================================================================================#

# Merge with the big dataframe

DF_3um_ForMergingCyto <- DF_SampleDataForPCA_rmCol_sortedDate_3um

DF_3um_withcyto <- merge(DF_3um_ForMergingCyto, CytometryFileTrunc, by = "Date_num", all.x = TRUE) # Check this merging

# Extrapolate the values with missMDA

DF_3um_withcyto_missMDA <- DF_3um_withcyto[,colnames(DF_3um_withcyto) %in% c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4",  "BaillauryHeigthIn_m", "CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC")]
nb <- estim_ncpPCA(DF_3um_withcyto_missMDA,ncp.max=5) # Output is 1
res.comp_DF_3um_withcyto_missMDA = imputePCA(DF_3um_withcyto_missMDA,ncp=1)
Extrapolated_DF_withCyto <- data.frame(res.comp_DF_3um_withcyto_missMDA$completeObs)

# 2-Replace the columns in the original dataset with the extrapolated columns
#----------------------------------------------------------------------------

# Check row names are in the same order
identical(row.names(DF_3um_withcyto), row.names(Extrapolated_DF_withCyto))

# Replace the original columns with the complete ones 

DF_3um_withcyto_NonExtrapCols <- DF_3um_withcyto[,!(colnames(DF_3um_withcyto) %in% c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4",  "BaillauryHeigthIn_m", "CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC"))]

FullExtrapDF_withCyto <- cbind(Extrapolated_DF_withCyto, DF_3um_withcyto_NonExtrapCols) # /!\ Hard check

# 3- Center and reduce the variables
#-----------------------------------

# 0- Remove columns that will not be used for the PCA 

FullExtrapDF_withCyto$Filter_Euk <- NULL
FullExtrapDF_withCyto$Date_num <- NULL 
FullExtrapDF_withCyto$DaysBefore <- NULL
FullExtrapDF_withCyto$DaysAfter <- NULL

# Add a column that groups Fungi_NA and Fungi_X

FullExtrapDF_withCyto$PercFungiUnknownPhylum <- FullExtrapDF_withCyto$PercFungi_NA + FullExtrapDF_withCyto$PercFungi_X

# 1- Center and reduce all variables

FullExtrapDF_withCyto_3um_CentRed <- data.frame(lapply(FullExtrapDF_withCyto, function(x) (x-mean(x))/sd(x) ))

# Check it worked 

lapply(FullExtrapDF_withCyto_3um_CentRed, mean) # The mean gets close to zero but is not exactly 0. Check if this is an issue
lapply(FullExtrapDF_withCyto_3um_CentRed, sd)

# 2- Remove columns that have NaN

FullExtrapDF_withCyto_3um_CentRed$ASV971 <- NULL
FullExtrapDF_withCyto_3um_CentRed$PercEntomophthoromycota <- NULL


#=================================================#
#                PCAs with Flow cytometry         #
#=================================================#

DF_3um_CentRed_woRR7daysDerBaill_10ASVs <- FullExtrapDF_withCyto_3um_CentRed
DF_3um_CentRed_woRR7daysDerBaill_10ASVs$DerivHeigth_m_j <- NULL
#DF_3um_CentRed_woRR7daysDerBaill_10ASVs$RR7days <- NULL

DF_3um_CentRed_woRR7daysDerBaill_10ASVs <- DF_3um_CentRed_woRR7daysDerBaill_10ASVs[,c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "BaillauryHeigthIn_m", "CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC", "PercFungi", "PercCryptomycota", "PercAscomycota", "PercChytridiomycota", "PercBasidiomycota", "PercFungiUnknownPhylum")]

#DF_3um_CentRed_woRR7daysDerBaill_10ASVs <- DF_3um_CentRed_woRR7daysDerBaill_10ASVs[,c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "BaillauryHeigthIn_m", "ASV224", "ASV275", "ASV1613", "ASV1016", "ASV409", "ASV166", "ASV594", "ASV847", "ASV367", "ASV553","CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC")]

# Recompute the original FullExtrapDF_withCyto to have the column with the date and be able to add the season

FullExtrapDF_withCyto_ForSeason <- cbind(Extrapolated_DF_withCyto, DF_3um_withcyto_NonExtrapCols)
FullExtrapDF_withCyto_ForSeason$Month <- month(FullExtrapDF_withCyto_ForSeason$Date_num)

# Make the seasons column

FullExtrapDF_withCyto_ForSeason$Season <- NA
for (i in 1:dim(FullExtrapDF_withCyto_ForSeason)[1]) {
  if (FullExtrapDF_withCyto_ForSeason[i,"Month"] %in% c(1, 2, 3)) {
    FullExtrapDF_withCyto_ForSeason[i, "Season"] <- "Winter"
  } else if (FullExtrapDF_withCyto_ForSeason[i,"Month"] %in% c(4, 5, 6)) {
    FullExtrapDF_withCyto_ForSeason[i, "Season"] <- "Spring"
  } else if (FullExtrapDF_withCyto_ForSeason[i,"Month"] %in% c(7, 8, 9)) {
    FullExtrapDF_withCyto_ForSeason[i, "Season"] <- "Summer"
  } else if (FullExtrapDF_withCyto_ForSeason[i,"Month"] %in% c(10, 11, 12)) {
    FullExtrapDF_withCyto_ForSeason[i, "Season"] <- "Autumn"
  } else {}
}

# Check the values are in the right order
SeasonBaillauryHeight <- (FullExtrapDF_withCyto_ForSeason$BaillauryHeigthIn_m - mean(FullExtrapDF_withCyto_ForSeason$BaillauryHeigthIn_m))/sd(FullExtrapDF_withCyto_ForSeason$BaillauryHeigthIn_m)
CytoBaillauryHeight <- (FullExtrapDF_withCyto$BaillauryHeigthIn_m - mean(FullExtrapDF_withCyto$BaillauryHeigthIn_m))/sd(FullExtrapDF_withCyto$BaillauryHeigthIn_m)
identical(SeasonBaillauryHeight, CytoBaillauryHeight) # Must be TRUE
identical(SeasonBaillauryHeight, DF_3um_CentRed_woRR7daysDerBaill_10ASVs$BaillauryHeigthIn_m) # Must be TRUE

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
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/PCA_Figure3_FinalPCA.pdf", plot=PCA_Figure3_CytoFlux_PercFungi, device = cairo_pdf(), width=15, height=10)
dev.off()

# Make a big heatmap to correlate all variables with all variables including the PCA dimensions (the heatmap production is at the end of the script)

PCA_3umNice
CoordIndivPCA <- PCA_3umNice$ind$coord
DF_MetadataWithDimesions <- cbind(FullExtrapDF_withCyto_ForSeason, CoordIndivPCA)

# Check the order of lines is good
DimensionsBaillauryNorm <- (DF_MetadataWithDimesions$BaillauryHeigthIn_m - mean(DF_MetadataWithDimesions$BaillauryHeigthIn_m))/sd(DF_MetadataWithDimesions$BaillauryHeigthIn_m)
identical(PCA_3umNice$call$X$BaillauryHeigthIn_m, DimensionsBaillauryNorm) # Must be TRUE

#########################################################################
#                         Correlations                                  #
#########################################################################

# With centered and scaled variables

AllAgainstAllCor <- cor(DF_3um_CentRed)
heatmap(AllAgainstAllCor)
AllAgainstAllCor[,"ASV224"]
sort(AllAgainstAllCor[,"ASV594"])
plot_ly(z = AllAgainstAllCor, type = "heatmap", )  

#===============================================
# WITH FLUX CYTOMETRY ADDED

AllAgainstAllCor_FluxCyto_sortedCols <- AllAgainstAllCor_FluxCyto[, c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "RR", "FFM", "RR7days", "BaillauryHeigthIn_m", "CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC", "PercFungi", "PercMucoromycota", "PercBlastocladiomycota", "PercCryptomycota", "PercAscomycota", "PercChytridiomycota", "PercBasidiomycota", "PercFungi_NA", "PercFungi_X" ,"ASV224", "ASV275", "ASV1613", "ASV1016", "ASV409", "ASV166", "ASV594", "ASV847", "ASV367", "ASV553")]

AllAgainstAllCor_FluxCyto <- cor(AllAgainstAllCor_FluxCyto_sortedCols)
heatmap(AllAgainstAllCor_FluxCyto, Rowv = NA, Colv = NA, scale = "none", revC = TRUE)

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(AllAgainstAllCor_FluxCyto)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

CorMatrixExtrapolatedData_3um <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   size = 10, hjust = 1))+
  coord_fixed()

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/CorMatrixExtrapolatedData_3um.pdf", plot=CorMatrixExtrapolatedData_3um, device = cairo_pdf(), width=15, height=10)
dev.off()

#################################################################
# Make the barplots of the relative abundance of Fungi by month #
#################################################################

DF_ForBarplot <- sample_data(Phyloseq_EukAll)

# 3 um
#----------

DF_ForBarplot_3um <- DF_ForBarplot[DF_ForBarplot[["Filter_Euk"]] == "3µM",]

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/BoxPlot_Figure4_3um_PercentFungiMonth.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
    #paper = "A4")
boxplot(DF_ForBarplot_3um$PercFungi ~ DF_ForBarplot_3um$Month, xlab = "Month", ylab = "Percent of Fungi amongst eukaryotes")
dev.off()

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_3um, aes(x=Month, y=PercFungi)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi amongst eukaryotes") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercFungi_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

p1 <- ggplot(data = DF_ForBarplot_3um, aes(x=Month, y=PercFungi)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi amongst eukaryotes") +
  facet_zoom(ylim = c(0, 14.8), zoom.size = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercFungi_3um_alt.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
dev.off()

# Ascomycota

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/BoxPlot_3um_PercentAscomycotaMonth.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_3um$PercAscomycota ~ DF_ForBarplot_3um$Month, xlab = "Month", ylab = "Percent of Ascomycota amongst eukaryotes")
dev.off()

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_3um, aes(x=Month, y=PercAscomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Ascomycota amongst eukaryotes") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercAscomycota_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Basidiomycota

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/BoxPlot_3um_PercentBasidiomycotaMonth.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_3um$PercBasidiomycota ~ DF_ForBarplot_3um$Month, xlab = "Month", ylab = "Percent of Basidiomycota amongst eukaryotes")
dev.off()

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_3um, aes(x=Month, y=PercBasidiomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Basidiomycota amongst eukaryotes") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercBasidiomycota_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Chytridiomycota --> Y axis must be cut

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/BoxPlot_3um_PercentChytridiomycotaMonth.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_3um$PercChytridiomycota ~ DF_ForBarplot_3um$Month, xlab = "Month", ylab = "Percent of Chytridiomycota amongst eukaryotes")
# gap.boxplot(DF_ForBarplot_3um$PercChytridiomycota ~ DF_ForBarplot_3um$Month, gap=list(top=c(2.7,4.8),bottom=c(NA,NA)), gap.axis="y", axes = FALSE)
# abline(h=seq(2.69,2.8,.001), col="white")
# axis.break(axis=2, breakpos=2.74,style="slash")
# #axis(2, labels=c(0:2.5,5:6), at=c(0:2.5,(5-(5-2.1)):(6-(5-2.1))))
# #axis(2, labels=c(0:2.5,5:6), at=c(0:2.5,(5-2.1),(6-2.1)))# at=c(0:2.5,(5-(5-2.1)):(6-(5-2.1))))
# axis(2, labels=c(0,1,2,2.5,5,6), at=c(0,1,2,2.5, 5-2.1, 6-2.1))
dev.off()

#dev.new()
# ggplot(data = DF_ForBarplot_3um, aes(x=Month, y=PercChytridiomycota)) +
# stat_boxplot(geom ='errorbar', width = 0.5) +
# geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
#   scale_y_break(c(2.3, 5.4), scales = "fixed", ticklabels = c(5.5,6)) +
#   ylim(c(0,6)) +
#   theme_bw()

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_3um, aes(x=Month, y=PercChytridiomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Chytridiomycota amongst eukaryotes") +
  facet_zoom(ylim = c(0, 1), zoom.size = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
    )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercChytridiomycota_3um.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
dev.off()

# Cryptomycota --> Y axis must be cut

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/BoxPlot_3um_PercentCryptomycotaMonth.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_3um$PercCryptomycota ~ DF_ForBarplot_3um$Month, xlab = "Month", ylab = "Percent of Cryptomycota amongst eukaryotes")
dev.off()

# ggplot(data = DF_ForBarplot_3um, aes(x=Month, y=PercCryptomycota)) +
#   stat_boxplot(geom ='errorbar', width = 0.5) +
#   geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
#   scale_y_break(c(2.1,7), scales = "fixed") + #, ticklabels = c(6,8,10)) +
#   ylim(c(0,10)) +
#   theme_bw()

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_3um, aes(x=Month, y=PercCryptomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Cryptomycota amongst eukaryotes") +
  facet_zoom(ylim = c(0, 1), zoom.size = 1) +
  theme_bw()+
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercCryptomycota_3um.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
dev.off()

# Unknown phyla

DF_ForBarplot_3um$PercFungiUnknownPhyla <- DF_ForBarplot_3um$PercFungi_NA + DF_ForBarplot_3um$PercFungi_X

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_3um, aes(x=Month, y=PercFungiUnknownPhyla)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi from unkwnown phylum amongst eukaryotes") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercUnknownPhylum_3um.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()

# Other phyla

DF_ForBarplot_3um$OtherPhyla <- DF_ForBarplot_3um$PercBlastocladiomycota + DF_ForBarplot_3um$PercEntomophthoromycota + DF_ForBarplot_3um$PercMucoromycota

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_3um, aes(x=Month, y=OtherPhyla)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi from other phyla amongst eukaryotes") +
  facet_zoom(ylim = c(0, 0.0000005), zoom.size = 1) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercOtherPhyla_3um.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
dev.off()

# ASVs

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/BoxPlot_3um_PercentASV224Month.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_3um$ASV224 ~ DF_ForBarplot_3um$Month, xlab = "Month", ylab = "Percent of ASV224 amongst eukaryotes")
dev.off()

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/BoxPlot_3um_PercentASV1016Month.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_3um$ASV1016 ~ DF_ForBarplot_3um$Month, xlab = "Month", ylab = "Percent of ASV1016 amongst eukaryotes")
dev.off()

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/BoxPlot_3um_PercentASV594Month.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_3um$ASV594 ~ DF_ForBarplot_3um$Month, xlab = "Month", ylab = "Percent of ASV594 amongst eukaryotes")
dev.off()



# 0.2 um
#----------

DF_ForBarplot_02um <- DF_ForBarplot[DF_ForBarplot[["Filter_Euk"]] == "Sterivex",]

# PercFungi --> Y axis must be cut

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/BoxPlot_Figure4_02um_PercentFungiMonth.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_02um$PercFungi ~ DF_ForBarplot_02um$Month, xlab = "Month", ylab = "Percent of Fungi amongst eukaryotes")
dev.off()

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_02um, aes(x=Month, y=PercFungi)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi amongst eukaryotes") +
  facet_zoom(ylim = c(0, 1.7), zoom.size = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercFungi_02um.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
dev.off()

p1 <- ggplot(data = DF_ForBarplot_02um, aes(x=Month, y=PercFungi)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi amongst eukaryotes") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercFungi_02um_alt.pdf", plot=p1, device = cairo_pdf(), width=10, height=10)
dev.off()


# Ascomycota 

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/BoxPlot_02um_PercentAscomycotaMonth.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_02um$PercAscomycota ~ DF_ForBarplot_02um$Month, xlab = "Month", ylab = "Percent of Ascomycota amongst eukaryotes")
dev.off()

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_02um, aes(x=Month, y=PercAscomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Ascomycota amongst eukaryotes") +
  facet_zoom(ylim = c(0, 0.75), zoom.size = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercAscomycota_02um.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
dev.off()

# Basidiomycota --> Y axis must be cut

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/BoxPlot_02um_PercentBasidiomycotaMonth.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_02um$PercBasidiomycota ~ DF_ForBarplot_02um$Month, xlab = "Month", ylab = "Percent of Basidiomycota amongst eukaryotes")
dev.off()

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_02um, aes(x=Month, y=PercBasidiomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Basidiomycota amongst eukaryotes") +
  facet_zoom(ylim = c(0, 0.15), zoom.size = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercBasidiomycota_02um.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
dev.off()

# Chytridiomycota

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/BoxPlot_02um_PercentChytridiomycotaMonth.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_02um$PercChytridiomycota ~ DF_ForBarplot_02um$Month, xlab = "Month", ylab = "Percent of Chytridiomycota amongst eukaryotes")
dev.off()

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_02um, aes(x=Month, y=PercChytridiomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Chytridiomycota amongst eukaryotes") +
  facet_zoom(ylim = c(0, 0.55), zoom.size = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercChytridiomycota_02um.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
dev.off()

# Cryptomycota --> Y axis must be cut

pdf("/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/BoxPlot_02um_PercentCryptomycotaMonth.pdf",         # File name
    width = 8, height = 8, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk") #,    # Color model (cmyk is required for most publications)
#paper = "A4")
boxplot(DF_ForBarplot_02um$PercCryptomycota ~ DF_ForBarplot_02um$Month, xlab = "Month", ylab = "Percent of Cryptomycota amongst eukaryotes")
dev.off()

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_02um, aes(x=Month, y=PercCryptomycota)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Cryptomycota amongst eukaryotes") +
  facet_zoom(ylim = c(0, 1.2), zoom.size = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercCryptomycota_02um.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
dev.off()

# Fungi from unknown phylum

DF_ForBarplot_02um$PercFungiUnknownPhylum <- DF_ForBarplot_02um$PercFungi_NA + DF_ForBarplot_02um$PercFungi_X

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_02um, aes(x=Month, y=PercFungiUnknownPhylum)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi from unkwnown phylum amongst eukaryotes") +
  facet_zoom(ylim = c(0, 0.4), zoom.size = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercUnknownPhylum_02um.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
dev.off()

# Fungi from other phyla 

DF_ForBarplot_02um$OtherPhyla <- DF_ForBarplot_02um$PercBlastocladiomycota + DF_ForBarplot_02um$PercEntomophthoromycota + DF_ForBarplot_02um$PercMucoromycota

#dev.new()
p1 <- ggplot(data = DF_ForBarplot_02um, aes(x=Month, y=OtherPhyla)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1) +
  ylab(label = "Percent of Fungi from other phyla amongst eukaryotes") +
  facet_zoom(ylim = c(0, 0.0000005), zoom.size = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20)
  )
plot(p1)

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Figure2/Fig2PercOtherPhyla_02um.pdf", plot=p1, device = cairo_pdf(), width=21, height=10)
dev.off()

################################################################################################
# Statistical tests to check if the percent of Fungi is significantly different between months #
################################################################################################

# /!\ AT LEAST homosedasticity and normality must be checked (Shapiro-Wilk test and )

# 3 um 
#-------

DF_ForBarplot_3um_RealDF <- data.frame(DF_ForBarplot_3um)
one.way_3um <- aov(PercFungi ~ Month, data = DF_ForBarplot_3um_RealDF)
summary(one.way_3um)
# tukey.one.way_3um <-TukeyHSD(one.way_3um)
# summary(tukey.one.way_3um)
post_test <- glht(one.way_3um,
                  linfct = mcp(Month = "Tukey") # Does not work
)

summary(post_test)

Tukey.one.way <- TukeyHSD(one.way_3um, conf.level=.95)
Tukey.one.way_DF <- data.frame(Tukey.one.way$Month)
Tukey.one.way_DF[order(Tukey.one.way_DF$p.adj), ]
pairwise.t.test(one.way_3um$model$PercFungi, one.way_3um$model$Month, p.adj='bonferroni')

# Test normality
#---------------

# Make a list of the values for each Month

ListValMon <- list()
for (Mon  in unique(DF_ForBarplot_3um_RealDF$Month)){
  V_temp <- DF_ForBarplot_3um_RealDF[DF_ForBarplot_3um_RealDF[["Month"]] == Mon,"PercFungi"]
  ListValMon[[Mon]] <- V_temp
} # Check results

lapply(ListValMon, shapiro.test) # The distribution is not normal, I should use a Kruskal-Wallis test instead

# Kruskal Wallis

kw_3um <- kruskal.test(PercFungi ~ Month, data = DF_ForBarplot_3um_RealDF)

# Kruskal-Wallis rank sum test
# 
# data:  PercFungi by Month
# Kruskal-Wallis chi-squared = 27.503, df = 11, p-value = 0.003856

# Pairwise comparisons

# pairwise.wilcox.test(DF_ForBarplot_3um_RealDF$PercFungi, DF_ForBarplot_3um_RealDF$Month,
#                      p.adjust.method = "BH") # Apparently not the good method

dunn_3um <- dunnTest(PercFungi ~ Month,
         data=DF_ForBarplot_3um_RealDF,
         method="bonferroni") # Nothing significant
View(dunn_3um$res[order(dunn_3um$res$P.adj), ])

dunn_3um_bh <- dunnTest(PercFungi ~ Month,
                     data=DF_ForBarplot_3um_RealDF,
                     method="bh")
dunn_3um_bh$res[order(dunn_3um_bh$res$P.adj), ]

# With functions from the rstatix package

# 3 um
kruskal_test(data = DF_ForBarplot_3um_RealDF, formula = PercFungi ~ Month) # Significant

# 02 um
DF_ForBarplot_02um_RealDF <- data.frame(DF_ForBarplot_02um)
kruskal_test(data = DF_ForBarplot_02um_RealDF, formula = PercFungi ~ Month) # Not significant

# Post-hoc test --> I use a Dunn test

# 3um
View(dunn_test(data = DF_ForBarplot_3um_RealDF, formula = PercFungi ~ Month, p.adjust.method = "bonferroni")) #
View(dunn_test(data = DF_ForBarplot_3um_RealDF, formula = PercFungi ~ Month, p.adjust.method = "BH")) # 

# Comparison between both size fractions

DF_ForBarplot_RealDF <- data.frame(DF_ForBarplot)
kruskal_test(data = DF_ForBarplot_RealDF, formula = PercFungi ~ Filter_Euk) # Not significant

#=======================================================================================#
#               Plot all environmental variables as facets in one big plot              # 
#=======================================================================================#

# Take the dataframe with everything but turbidity

DF_3um_withcyto_ForMetadataPlot <- DF_3um_withcyto

# Take the dataframe with turbidity, then select only 3 µm samples

DF_TUrbidity_ForMetadataPlot <- DF_SampleDataForPCA
DF_TUrbidity_ForMetadataPlot_3um <- DF_TUrbidity_ForMetadataPlot[DF_TUrbidity_ForMetadataPlot[["Filter_Euk"]] == "3µM",]
DF_TUrbidity_ForMetadataPlot_3um <- DF_TUrbidity_ForMetadataPlot_3um[,c("Date_num", "Turbidity_3m", "Turbidity_20m")]

# Merge the two dataframes

DF_Merged_ForMetadataPlot_3um <- merge(DF_3um_withcyto_ForMetadataPlot, DF_TUrbidity_ForMetadataPlot_3um, by = "Date_num") # Check

# Select variables to be plotted

DF_Merged_ForMetadataPlot_3um_SlctCol <- DF_Merged_ForMetadataPlot_3um[, c("Date_num", "T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "BaillauryHeigthIn_m", "RR", "FFM", "Turbidity_3m", "Turbidity_20m", "CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC")]

# Rename columns with their unit

colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "T")] <- " T (°C)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "S")] <- " S (PSU)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "O")] <- " O (mL/L)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "pH")] <- "pH"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "NH4")] <- "NH4 (µM)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "NO3")] <- "NO3 (µM)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "NO2")] <- "NO2 (µM)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "PO4")] <- "PO4 (µM)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "SIO4")] <- "SIO4 (µM)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "BaillauryHeigthIn_m")] <- "Baillaury height (m)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "RR")] <- "RR (mm)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "FFM")] <- "FFM (m/s)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "Turbidity_3m")] <- "Turb 3 m (NTU)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "Turbidity_20m")] <- "Turb 20 m (NTU)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "CHLA")] <- "CHLA (µg/L)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "SYNC")] <- "SYNC (cell/mL)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "PROC")] <- "PROC (cell/mL)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "PICOEC")] <- "PICOEC (cell/mL)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "NANOEC")] <- "NANOEC (cell/mL)"
colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[which(colnames(DF_Merged_ForMetadataPlot_3um_SlctCol) == "CRYC")] <- "CRYC (cell/mL)"

# Make the plot

DF_LongerMetadataForPlot <- DF_Merged_ForMetadataPlot_3um_SlctCol %>%
  pivot_longer(!Date_num, names_to = "EnvParam", values_to = "Values")
DF_LongerMetadataForPlot$EnvParamFactor <- factor(DF_LongerMetadataForPlot$EnvParam, levels=colnames(DF_Merged_ForMetadataPlot_3um_SlctCol)[2:length(DF_Merged_ForMetadataPlot_3um_SlctCol)])

p1 <- ggplot(data = DF_LongerMetadataForPlot, aes(x = Date_num, y = Values)) +
  geom_point()+
  geom_line() +
  facet_grid(rows = vars(EnvParamFactor), scales = "free") +
  theme_bw()

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/EnvironmentalParameters/FigEnvironmentalParameters.pdf", plot=p1, device = cairo_pdf(), width=21, height=50, limitsize = FALSE)
dev.off()

# Make correlations between all variables and all dimensions of the PCA

DF_MetadataWithDimesions_ColSelect <- DF_MetadataWithDimesions[,c("Date_num", "Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5", "Dim.6", "Dim.7", "Dim.8", "Dim.9", "Dim.10")]
DF_MetadataForHeatmap <- merge(DF_MetadataWithDimesions_ColSelect, DF_Merged_ForMetadataPlot_3um, by = "Date_num")
DF_MetadataForHeatmap$PercUnknownPhylum <- DF_MetadataForHeatmap$PercFungi_NA + DF_MetadataForHeatmap$PercFungi_X

# Add the 41st ASV

PhyloseqFor41stASV <- OTUPhyloseq_NormDomain
PhyloseqFor41stASV <- subset_taxa(PhyloseqFor41stASV, Division == "Fungi")

topN_41 <- 41
most_abundant_taxa_41 <- sort(taxa_sums(PhyloseqFor41stASV), TRUE)[1:topN_41]
PhyloseqFor41stASV_short <- prune_taxa(names(most_abundant_taxa_41), PhyloseqFor41stASV)

PhyloseqFor41stASV_short_t <- data.frame(t(otu_table(PhyloseqFor41stASV_short)))

# Check the order of columns 

identical(sample_data(Phyloseq_EukAll)[["Name_Euk"]], row.names(PhyloseqFor41stASV_short_t)) # It works

DF_For41stASV <- cbind(sample_data(Phyloseq_EukAll), PhyloseqFor41stASV_short_t) # /!\ To be checked

# Keep only 3 um size fraction

DF_For41stASV_3um <- DF_For41stASV[DF_For41stASV[["Filter_Euk"]] == "3µM",]

DF_For41stASV_3um_short <- DF_For41stASV[,c("Date_num", "ASV1367")]

# Merge both dataframes 

DF_MetadataForHeatmapPlus41 <- merge(DF_MetadataForHeatmap, DF_For41stASV_3um_short, by = "Date_num")
DF_MetadataForHeatmap <- DF_MetadataForHeatmapPlus41

# Remove columns I will not use
DF_MetadataForHeatmap$Filter_Euk <- NULL
DF_MetadataForHeatmap$Date_num <- NULL
DF_MetadataForHeatmap$DaysAfter <- NULL
DF_MetadataForHeatmap$DaysBefore <- NULL
DF_MetadataForHeatmap$RR7days <- NULL
DF_MetadataForHeatmap$PercFungi_NA <- NULL
DF_MetadataForHeatmap$PercFungi_X <- NULL
DF_MetadataForHeatmap$DerivHeigth_m_j <- NULL
DF_MetadataForHeatmap$ASV971 <- NULL # Every value is 0 in the > 3 um size fraction
DF_MetadataForHeatmap$PercEntomophthoromycota <- NULL # Every value is 0 in the > 3 um size fraction

# Rename columns for the final plot

colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Dim.1")] <- "PCA Dim 1"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Dim.2")] <- "PCA Dim 2"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Dim.3")] <- "PCA Dim 3"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Dim.4")] <- "PCA Dim 4"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Dim.5")] <- "PCA Dim 5"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Dim.6")] <- "PCA Dim 6"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Dim.7")] <- "PCA Dim 7"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Dim.8")] <- "PCA Dim 8"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Dim.9")] <- "PCA Dim 9"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Dim.10")] <- "PCA Dim 10"

colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "PercFungi")] <- "Fungi"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "PercAscomycota")] <- "Ascomycota"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "PercBasidiomycota")] <- "Basidiomycota"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "PercChytridiomycota")] <- "Chytridiomycota"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "PercCryptomycota")] <- "Cryptomycota"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "PercMucoromycota")] <- "Mucoromycota"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "PercBlastocladiomycota")] <- "Blastocladiomycota"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "BaillauryHeigthIn_m")] <- "Baillaury height"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Turbidity_3m")] <- "Turbidity 3 m"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "Turbidity_20m")] <- "Turbidity 20 m"
colnames(DF_MetadataForHeatmap)[which(colnames(DF_MetadataForHeatmap) == "PercUnknownPhylum")] <- "Fungi unknown phylum"

# Reorder columns as desired in the plot 

DF_MetadataForHeatmapOrder <- DF_MetadataForHeatmap[, c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "Baillaury height", "RR", "FFM", "Turbidity 3 m", "Turbidity 20 m", "CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC", "Fungi", "Ascomycota", "Basidiomycota", "Chytridiomycota", "Cryptomycota", "Fungi unknown phylum", "Mucoromycota", "Blastocladiomycota", "ASV166", "ASV224", "ASV275", "ASV324", "ASV367", "ASV409", "ASV420", "ASV460", "ASV505", "ASV553", "ASV594", "ASV774", "ASV843", "ASV847", "ASV866", "ASV1016", "ASV1163", "ASV1367", "ASV1421", "ASV1613", "ASV1633", "ASV1724", "ASV1991", "ASV2336", "ASV2651", "ASV2737", "ASV2777", "ASV3073", "ASV3729", "ASV5211", "ASV5590", "ASV5782", "ASV6352", "ASV6458", "ASV8230", "ASV8944", "ASV8958", "ASV8967", "ASV8969", "ASV9009", "PCA Dim 1", "PCA Dim 2", "PCA Dim 3", "PCA Dim 4", "PCA Dim 5", "PCA Dim 6", "PCA Dim 7", "PCA Dim 8", "PCA Dim 9", "PCA Dim 10")]

# Make the correlation and the plot

#CorrelationMatrixForHeatmap <- round(cor(DF_MetadataForHeatmapOrder, use = "pairwise.complete.obs", method = "pearson"), 2)
CorrelationMatrixForHeatmap <- cor_test(DF_MetadataForHeatmapOrder, use = "pairwise.complete.obs", method = "pearson")
CorrelationMatrixForHeatmap_rounded <- CorrelationMatrixForHeatmap
CorrelationMatrixForHeatmap_rounded$cor_round <- round(CorrelationMatrixForHeatmap_rounded$cor, 2)
CorrelationMatrixForHeatmap_rounded$p_round <- round(CorrelationMatrixForHeatmap_rounded$p, 3)
# CorrelationMatrixForHeatmap_rounded_short <- CorrelationMatrixForHeatmap_rounded[,c("var1", "var2", "cor_round", "p_round")]

# Put a new column where the value is the concatenation of cor and p if the p-value is below 0.05

CorrelationMatrixForHeatmap_rounded$corPvalue <- paste("r=",CorrelationMatrixForHeatmap_rounded$cor_round, "\nP=", CorrelationMatrixForHeatmap_rounded$p_round, sep = "")

CorrelationMatrixForHeatmap_rounded[which(CorrelationMatrixForHeatmap_rounded[["p"]] > 0.05), "corPvalue"] <- ""
CorrelationMatrixForHeatmap_rounded[which(is.na(CorrelationMatrixForHeatmap_rounded[["cor"]])), "corPvalue"] <- ""

CorrelationMatrixForHeatmap_rounded_DF <- data.frame(CorrelationMatrixForHeatmap_rounded)
CorrelationMatrixForHeatmap_rounded_DF$var1 <- factor(CorrelationMatrixForHeatmap_rounded_DF$var1, levels = c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "Baillaury height", "RR", "FFM", "Turbidity 3 m", "Turbidity 20 m", "CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC", "Fungi", "Ascomycota", "Basidiomycota", "Chytridiomycota", "Cryptomycota", "Fungi unknown phylum", "Mucoromycota", "Blastocladiomycota", "ASV166", "ASV224", "ASV275", "ASV324", "ASV367", "ASV409", "ASV420", "ASV460", "ASV505", "ASV553", "ASV594", "ASV774", "ASV843", "ASV847", "ASV866", "ASV1016", "ASV1163", "ASV1367", "ASV1421", "ASV1613", "ASV1633", "ASV1724", "ASV1991", "ASV2336", "ASV2651", "ASV2737", "ASV2777", "ASV3073", "ASV3729", "ASV5211", "ASV5590", "ASV5782", "ASV6352", "ASV6458", "ASV8230", "ASV8944", "ASV8958", "ASV8967", "ASV8969", "ASV9009", "PCA Dim 1", "PCA Dim 2", "PCA Dim 3", "PCA Dim 4", "PCA Dim 5", "PCA Dim 6", "PCA Dim 7", "PCA Dim 8", "PCA Dim 9", "PCA Dim 10"))
CorrelationMatrixForHeatmap_rounded_DF$var2 <- factor(CorrelationMatrixForHeatmap_rounded_DF$var2, levels = c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "Baillaury height", "RR", "FFM", "Turbidity 3 m", "Turbidity 20 m", "CHLA", "SYNC", "PROC", "PICOEC", "NANOEC", "CRYC", "Fungi", "Ascomycota", "Basidiomycota", "Chytridiomycota", "Cryptomycota", "Fungi unknown phylum", "Mucoromycota", "Blastocladiomycota", "ASV166", "ASV224", "ASV275", "ASV324", "ASV367", "ASV409", "ASV420", "ASV460", "ASV505", "ASV553", "ASV594", "ASV774", "ASV843", "ASV847", "ASV866", "ASV1016", "ASV1163", "ASV1367", "ASV1421", "ASV1613", "ASV1633", "ASV1724", "ASV1991", "ASV2336", "ASV2651", "ASV2737", "ASV2777", "ASV3073", "ASV3729", "ASV5211", "ASV5590", "ASV5782", "ASV6352", "ASV6458", "ASV8230", "ASV8944", "ASV8958", "ASV8967", "ASV8969", "ASV9009", "PCA Dim 1", "PCA Dim 2", "PCA Dim 3", "PCA Dim 4", "PCA Dim 5", "PCA Dim 6", "PCA Dim 7", "PCA Dim 8", "PCA Dim 9", "PCA Dim 10"))
#CorrelationMatrixForHeatmap_rounded_DF[which(is.na(CorrelationMatrixForHeatmap_rounded_DF[["cor"]])), "cor"] <- 0

# Make the plot /!\ Check the plot is good

#melted_corr_mat <- melt(CorrelationMatrixForHeatmap)
#head(melted_corr_mat)

#txtsize <- par('din')[2] / 2
pHeatmap <- ggplot(data = CorrelationMatrixForHeatmap_rounded_DF, aes(x=var1, y=var2, 
                                   fill=cor)) + 
  geom_tile(aes(height=1)) +
  #geom_raster() +
  geom_text(aes(var1, var2, label = corPvalue), # before it was var1, var2
            color = "black", size = 2) +
  #scale_fill_gradient2(limits = c(-1, 1)) +
  scale_fill_gradient2(midpoint = 0, low = "red", mid = "white",
                        high = "deepskyblue", limits = c(-1, 1)) +
  #coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) #, axis.text.y = element_text(vjust = 1))

pHeatmap

ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/EnvironmentalParameters/HeatmapEnvironmentalParameters.pdf", plot=pHeatmap, device = cairo_pdf(), width=30, height=21, limitsize = FALSE)
dev.off()

