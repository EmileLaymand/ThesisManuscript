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
library("pheatmap")
library(RColorBrewer)
library(viridis)
library(paletteer)

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

#################################################################################################
# End of section copy pasted from MetaB_SOLA_PCA_MeteoFranceTurbidityBaillauryRR7daysOffset_3um #
#################################################################################################

# Make the heatmap with only Fungi
#------------------------------------

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

topN <- 40 # Before, it was 20
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll20ASV <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll)

#--------------------------------------------------------------------------------------------------------------------------
# 3- Extract the otu table, rename the rows according to dates, then draw the heatmap
#--------------------------------------------------------------------------------------------------------------------------

FungiOTUPercPhyloseq_EukAll20ASV_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll20ASV, Filter_Euk=="3µM")

DF_ASV40Fungi_3um <- data.frame(t(data.frame(otu_table(FungiOTUPercPhyloseq_EukAll20ASV_3um))))
DF_ASV40Fungi_3um[["Name_Euk"]] <- row.names(DF_ASV40Fungi_3um)
DF_SampleData_3um <- data.frame(sample_data(FungiOTUPercPhyloseq_EukAll20ASV_3um))
DF_SampleData_3um <- DF_SampleData_3um[,c("Name_Euk", "Date_num")]

DF_ASV40Fungi_3um_dates <- merge(DF_ASV40Fungi_3um, DF_SampleData_3um, by = "Name_Euk", all.x = TRUE)
row.names(DF_ASV40Fungi_3um_dates) <- DF_ASV40Fungi_3um_dates$Date_num
DF_ASV40Fungi_3um_dates$Name_Euk <- NULL
DF_ASV40Fungi_3um_dates$Date_num <- NULL

DF_ASV40Fungi_3um_dates_t <- data.frame(t(DF_ASV40Fungi_3um_dates))

pheatmap(DF_ASV40Fungi_3um_dates_t, cutree_rows = 10, cluster_cols = FALSE, color = paletteer_c("grDevices::Purples 3", 30, direction = -1))
#pheatmap(DF_ASV40Fungi_3um_dates_t, cutree_rows = 10, cluster_cols = FALSE, color = paletteer_c("grDevices::Blue-Red 3", 30, direction = 1), breaks = c(0, 1, 25))

# Log10
DF_ASV40Fungi_3um_dates_t_NA <- DF_ASV40Fungi_3um_dates_t
DF_ASV40Fungi_3um_dates_t_NA[DF_ASV40Fungi_3um_dates_t_NA == 0] <- NA
DF_ASV40Fungi_3um_dates_t_NA_log10 <- log10(DF_ASV40Fungi_3um_dates_t_NA)

pheatmap(DF_ASV40Fungi_3um_dates_t_NA_log10, cutree_rows = 10, cluster_cols = FALSE, color = paletteer_c("grDevices::Purples 3", 30, direction = -1) )

#sqrt

DF_ASV40Fungi_3um_dates_t_sqrt <- sqrt(DF_ASV40Fungi_3um_dates_t)

pheatmap(DF_ASV40Fungi_3um_dates_t_sqrt, cutree_rows = 10, cluster_cols = FALSE, color = paletteer_c("grDevices::Purples 3", 30, direction = -1) )

# Each column divided by its own standard deviation 

DF_ASV40Fungi_3um_dates_t_sd <- data.frame(lapply(DF_ASV40Fungi_3um_dates_t, function(x) x/sd(x)))
row.names(DF_ASV40Fungi_3um_dates_t_sd) <- row.names(DF_ASV40Fungi_3um_dates_t) # TO BE ABSOLUTELY CHECKED

pheatmap(DF_ASV40Fungi_3um_dates_t_sd, cutree_rows = 10, cluster_cols = FALSE, color = paletteer_c("grDevices::Purples 3", 30, direction = -1) )
