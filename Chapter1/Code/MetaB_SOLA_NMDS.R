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
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above1, "Above40"] <- TRUE

sample_data(Phyloseq_EukAll)[["FungiRange"]] <- "Below1"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet1And5, "FungiRange"] <- "Bet1And5"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet5And10, "FungiRange"] <- "Bet5And10"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet10And15, "FungiRange"] <- "Bet10And15"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet15And20, "FungiRange"] <- "Bet15And20"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Bet20And30, "FungiRange"] <- "Bet20And30"
sample_data(Phyloseq_EukAll)[row.names(sample_data(Phyloseq_EukAll)) %in% Above30FungiRange, "FungiRange"] <- "Above30"

##########
## NMDS ##
##########

# Normalize to % of eukaryotes at the domain level
OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)
OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "Division")

# Total data (all kingdoms)
#--------------------------

# stat_plot_repel(Phyloseq_Euk3,
#                 StatMethod="NMDS",
#                 StatMethodDistance="bray",
#                 Normalization="relative",
#                 Colorparameter=NULL,
#                 Shapeparameter=NULL,
#                 OrderLegend=NULL,
#                 PointLabel="Date_Euk",
#                 Fillwith=NULL,
#                 Stat_ellipse=NULL,
#                 Title=NULL,
#                 PlotMethod="Repel",
#                 FilePath=NULL)

# Fungi only 
#-----------

# OTU LEVEL

# Check if the Fungi composition is different between 0.2 and 3 micrometers fractions 
#====================================================================================

OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)
OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "OTU")

FungiPhyloseq_EukAll <- subset_taxa(OTUPhyloseq_NormDomain, Division == "Fungi")
FungiPhyloseq_EukAll <- prune_taxa(taxa_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)
FungiPhyloseq_EukAll <- prune_samples(sample_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)

# stat_plot_repel(FungiPhyloseq_EukAll,
#                 StatMethod="NMDS",
#                 StatMethodDistance="bray",
#                 Normalization="Raw",
#                 Colorparameter=c("indianred4","darkorange","deepskyblue","#009300",
#                                  "#F5D11D", "#009999", "#990000","#FF9999", "#FFFF00", "black", "blue", "red"),
#                 Shapeparameter="Filter_Euk",
#                 OrderLegend=NULL,
#                 PointLabel="Name_Euk",
#                 Fillwith="Month",
#                 Stat_ellipse=NULL,
#                 Title=NULL,
#                 PlotMethod="Repel",
#                 FilePath=NULL)

# Only BIG size fraction
#-----------------------

FungiPhyloseq_3um <- subset_samples(FungiPhyloseq_EukAll, Filter_Euk=="3µM")
# stat_plot_repel(FungiPhyloseq_3um,
#                 StatMethod="NMDS",
#                 StatMethodDistance="bray",
#                 Normalization="Raw",
#                 Colorparameter=c("indianred4","darkorange","deepskyblue","#009300",
#                                  "#F5D11D", "#009999", "#990000","#FF9999", "#FFFF00", "black", "blue", "red"),
#                 Shapeparameter="Filter_Euk",
#                 OrderLegend=NULL,
#                 PointLabel="Name_Euk",
#                 Fillwith="Month",
#                 Stat_ellipse=NULL,
#                 Title=NULL,
#                 PlotMethod="Repel",
#                 FilePath=NULL)

FungiPhyloseq_3um_ord <- ordinate(FungiPhyloseq_3um, "NMDS", "bray")
p1 = plot_ordination(FungiPhyloseq_3um, FungiPhyloseq_3um_ord, type="samples", color="Month", shape="Year")
print(p1)
p1 = plot_ordination(FungiPhyloseq_3um, FungiPhyloseq_3um_ord, type="samples", color="FungiRange", shape="Year")
print(p1)
p1 = plot_ordination(FungiPhyloseq_3um, FungiPhyloseq_3um_ord, type="samples", color="Above20", shape="Year")
print(p1)

FungiPhyloseq_3um_ord <- ordinate(FungiPhyloseq_3um, "PCoA", "bray")
p1 = plot_ordination(FungiPhyloseq_3um, FungiPhyloseq_3um_ord, type="samples", color="Month", shape="Year")
print(p1)

#ord_explore(FungiPhyloseq_3um)

# Only SMALL size fraction
#-------------------------

OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)
OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "OTU")

FungiPhyloseq_EukAll <- subset_taxa(OTUPhyloseq_NormDomain, Division == "Fungi")
FungiPhyloseq_EukAll <- prune_taxa(taxa_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)
FungiPhyloseq_EukAll <- prune_samples(sample_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)
FungiPhyloseq_02um <- subset_samples(FungiPhyloseq_EukAll, Filter_Euk=="Sterivex")

FungiPhyloseq_02um_ord <- ordinate(FungiPhyloseq_02um, "NMDS", "bray")
p1 = plot_ordination(FungiPhyloseq_02um, FungiPhyloseq_02um_ord, type="samples", color="Month", shape="Year")
print(p1)

# CLASS (Phylum) LEVEL
#=====================

OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)
OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "OTU")

FungiPhyloseq_EukAll <- subset_taxa(OTUPhyloseq_NormDomain, Division == "Fungi")
FungiPhyloseq_EukAll <- prune_taxa(taxa_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)
FungiPhyloseq_EukAll <- prune_samples(sample_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)
FungiPhyloseq_EukAll <- tax_glom(FungiPhyloseq_EukAll, "Class", NArm = FALSE)

# Both size fractions

FungiPhyloseq_all_ord <- ordinate(FungiPhyloseq_EukAll, "NMDS", "bray")
p1 = plot_ordination(FungiPhyloseq_EukAll, FungiPhyloseq_all_ord, type="samples", color="Filter_Euk", shape="Year")
print(p1)

# Only 3um

OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)
OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "OTU")

FungiPhyloseq_EukAll <- subset_taxa(OTUPhyloseq_NormDomain, Division == "Fungi")
FungiPhyloseq_EukAll <- prune_taxa(taxa_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)
FungiPhyloseq_EukAll <- prune_samples(sample_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)
FungiPhyloseq_EukAll <- tax_glom(FungiPhyloseq_EukAll, "Class", NArm = FALSE)
FungiPhyloseq_3um <- subset_samples(FungiPhyloseq_EukAll, Filter_Euk=="3µM")

FungiPhyloseq_3um_ord <- ordinate(FungiPhyloseq_3um, "NMDS", "bray")
p1 = plot_ordination(FungiPhyloseq_3um, FungiPhyloseq_3um_ord, type="samples", color="Month", shape="Year")
print(p1)

# Only 0.2 um

OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)
OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "OTU")

FungiPhyloseq_EukAll <- subset_taxa(OTUPhyloseq_NormDomain, Division == "Fungi")
FungiPhyloseq_EukAll <- prune_taxa(taxa_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)
FungiPhyloseq_EukAll <- prune_samples(sample_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)
FungiPhyloseq_EukAll <- tax_glom(FungiPhyloseq_EukAll, "Class", NArm = FALSE)
FungiPhyloseq_02um <- subset_samples(FungiPhyloseq_EukAll, Filter_Euk=="Sterivex")

FungiPhyloseq_02um_ord <- ordinate(FungiPhyloseq_02um, "NMDS", "bray")
p1 = plot_ordination(FungiPhyloseq_02um, FungiPhyloseq_02um_ord, type="samples", color="Month", shape="Year")
print(p1)


#########################################################################################
#                              NMDS with everything but Fungi                           #
#########################################################################################

# The goal of this section is to make NMDS at the ASV and kingdom level with everything but Fungi, in order to study what is around Fungi

# Remove Fungi

Phyloseq_EukAll_woFungi <- subset_taxa(Phyloseq_EukAll, Division != "Fungi")

# Convert to percent 

Phyloseq_woFungi_Perc <- prune_taxa(taxa_sums(Phyloseq_EukAll_woFungi) > 0, Phyloseq_EukAll_woFungi)
Phyloseq_woFungi_Perc <- prune_samples(sample_sums(Phyloseq_woFungi_Perc) > 0, Phyloseq_woFungi_Perc)
Phyloseq_woFungi_Perc <- taxa_percentize(Phyloseq_woFungi_Perc, TaxLevel = "OTU")

# Keep only Fungi, and add to the sample data of Phyloseq_woFungi_Perc the percent of Fungi

OTUPhyloseq_ForPercFung <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_ForPercFung <- prune_samples(sample_sums(OTUPhyloseq_ForPercFung) > 0, OTUPhyloseq_ForPercFung)
OTUPhyloseq_ForPercFung <- taxa_percentize(OTUPhyloseq_ForPercFung, TaxLevel = "Division")
OTUPhyloseq_ForPercFung <- subset_taxa(OTUPhyloseq_ForPercFung, Division == "Fungi")
DFAbundFungi <- data.frame(t(otu_table(OTUPhyloseq_ForPercFung)["ASV224",]))
DFAbundFungi$Name_Euk <- row.names(DFAbundFungi)

identical(sample_data(Phyloseq_woFungi_Perc)[["Name_Euk"]], DFAbundFungi$Name_Euk) 
sample_data(Phyloseq_woFungi_Perc)[["PercFungi"]] <- DFAbundFungi$ASV224

# Check the addition of column is good. 20240118 Looks good

sample_data(Phyloseq_woFungi_Perc)[sample_data(Phyloseq_woFungi_Perc)[["Name_Euk"]] == "PF192Euk", "PercFungi"]
DFAbundFungi[DFAbundFungi$Name_Euk == "PF192Euk","ASV224"]

sample_data(Phyloseq_woFungi_Perc)[sample_data(Phyloseq_woFungi_Perc)[["Name_Euk"]] == "PF308Euk", "PercFungi"]
DFAbundFungi[DFAbundFungi$Name_Euk == "PF308Euk","ASV224"]

sample_data(Phyloseq_woFungi_Perc)[sample_data(Phyloseq_woFungi_Perc)[["Name_Euk"]] == "S08Euk", "PercFungi"]
DFAbundFungi[DFAbundFungi$Name_Euk == "S08Euk","ASV224"]

# Add the proportion of the 20 most abundant fungal ASVs to the sample_data

OTUPhyloseq_NormDomain_MostAbFungASVs <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain_MostAbFungASVs <- prune_samples(sample_sums(OTUPhyloseq_NormDomain_MostAbFungASVs) > 0, OTUPhyloseq_NormDomain_MostAbFungASVs)
OTUPhyloseq_NormDomain_MostAbFungASVs <- taxa_percentize(OTUPhyloseq_NormDomain_MostAbFungASVs, TaxLevel = "OTU")
OTUPhyloseq_NormDomain_MostAbFungASVs <- subset_taxa(OTUPhyloseq_NormDomain_MostAbFungASVs, Division == "Fungi")
MostAbFungASVs <- names(sort(taxa_sums(OTUPhyloseq_NormDomain_MostAbFungASVs), decreasing = TRUE)[1:20])
DF_MostAbFungASVs <- data.frame(t(otu_table(OTUPhyloseq_NormDomain_MostAbFungASVs)))
DF_MostAbFungASVs <- DF_MostAbFungASVs[,colnames(DF_MostAbFungASVs) %in% MostAbFungASVs]


# Make the NMDS #
#----------------

# 3 um

Phyloseq_woFungi_Perc_3um <- subset_samples(Phyloseq_woFungi_Perc, Filter_Euk=="3µM")

Phyloseq_woFungi_Perc_3um_ord <- ordinate(Phyloseq_woFungi_Perc_3um, "NMDS", "bray")
p1 = plot_ordination(Phyloseq_woFungi_Perc_3um, Phyloseq_woFungi_Perc_3um_ord, type="samples", color="PercFungi", shape="Year")
print(p1) # One sample goes really out, and seems to contain a lot of Fungi
ggplotly(p1)

Phyloseq_woFungi_Perc_3um <- subset_samples(Phyloseq_woFungi_Perc_3um, Date_Euk != "29_03_2016")

Phyloseq_woFungi_Perc_3um_ord <- ordinate(Phyloseq_woFungi_Perc_3um, "NMDS", "bray")
p1 = plot_ordination(Phyloseq_woFungi_Perc_3um, Phyloseq_woFungi_Perc_3um_ord, type="samples", color="PercFungi", shape="Year") +
  annotate("text", x=-1, y=0.6, label= sprintf("Stress = %s",as.character(round(Phyloseq_woFungi_Perc_3um_ord$stress, digits = 2))), size = 8) +
  scale_colour_gradient2() +
  geom_point(size=5) +
  scale_shape_manual(values=c(3, 7, 16, 15, 17))+
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
print(p1) # Samples with most Fungi seem to be in the center of all points
ggplotly(p1)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/NMDS_WithoutFungiASVLevel_3um_ColorPercFungi.pdf", plot=p1, device = cairo_pdf(), width=15, height=14)
dev.off()

p1 = plot_ordination(Phyloseq_woFungi_Perc_3um, Phyloseq_woFungi_Perc_3um_ord, type="samples", color="Month", shape="Year") +
  annotate("text", x=-1, y=0.6, label= sprintf("Stress = %s",as.character(round(Phyloseq_woFungi_Perc_3um_ord$stress, digits = 2))), size = 8) +
  geom_point(size=5) +
  scale_shape_manual(values=c(3, 7, 16, 15, 17))+
  theme_bw() +
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
print(p1) # It seems that samples from the same months really group together, that is really clear.
#ggplotly(p1)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/NMDS_WithoutFungiASVLevel_3um_ColorMonth.pdf", plot=p1, device = cairo_pdf(), width=15, height=14)
dev.off()

p1 = plot_ordination(Phyloseq_woFungi_Perc_3um, Phyloseq_woFungi_Perc_3um_ord, type="samples", color="Above5", shape="Year") +
  annotate("text", x=-1, y=0.6, label= sprintf("Stress = %s",as.character(round(Phyloseq_woFungi_Perc_3um_ord$stress, digits = 2))), size = 8) +
  geom_point(size=5) +
  scale_shape_manual(values=c(3, 7, 16, 15, 17))+
  theme_bw() +
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
print(p1) # It seems that samples from the same months really group together, that is really clear.
ggplotly(p1)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/NMDS_WithoutFungiASVLevel_3um_ColorAbove5PercFungi.pdf", plot=p1, device = cairo_pdf(), width=15, height=14)
dev.off()

# Just to test : which months gather the most samples with > 5 % of Fungi

Phyloseq_woFungi_Perc_3um_above5 <- subset_samples(Phyloseq_woFungi_Perc_3um, Above5 == TRUE)
TablePercFung <- table(sample_data(Phyloseq_woFungi_Perc_3um_above5)[["Month"]])
VectMissMonth <- c("05" = 0, "10" = 0, "11" = 0)
AllMonth <- c(TablePercFung, VectMissMonth)
AllMonthSortedNames <- AllMonth[sort(names(AllMonth))]
barplot(AllMonthSortedNames, xlab = " Month", ylab = "Number of samples with Fungi > 5 % of eukaryotes")
# Normalize by the number of samples in each month (otherwise, winter months are favorised over other months)
AllMonthSortedNamesAllSamples <- table(sample_data(Phyloseq_woFungi_Perc_3um)[["Month"]])
identical(names(AllMonthSortedNames), names(AllMonthSortedNamesAllSamples))
PercSamplesAbove5PercFungi <- (AllMonthSortedNames/AllMonthSortedNamesAllSamples)*100
barplot(PercSamplesAbove5PercFungi, xlab = " Month", ylab = "Percent of samples with Fungi > 5 % of eukaryotes")


# Idem but with the data glommed at the Kingdom level

Phyloseq_woFungi_Perc_3um_KingdomGlom <- tax_glom(Phyloseq_woFungi_Perc_3um, "Division", NArm = FALSE)

Phyloseq_woFungi_Perc_3um_KingdomGlom_ord <- ordinate(Phyloseq_woFungi_Perc_3um_KingdomGlom, "NMDS", "bray")
p1 = plot_ordination(Phyloseq_woFungi_Perc_3um_KingdomGlom, Phyloseq_woFungi_Perc_3um_KingdomGlom_ord, type="samples", color="PercFungi", shape="Year") +
  scale_colour_gradient2()
print(p1) # Samples with lots of Fungi do not really stand out from the NMDS 
ggplotly(p1)
p1 = plot_ordination(Phyloseq_woFungi_Perc_3um_KingdomGlom, Phyloseq_woFungi_Perc_3um_KingdomGlom_ord, type="samples", color="Month", shape="Year")
print(p1) # It seems that samples from the same months group together, less clear than at the ASV level
ggplotly(p1)

# 0.2 um

Phyloseq_woFungi_Perc_02um <- subset_samples(Phyloseq_woFungi_Perc, Filter_Euk=="Sterivex")

Phyloseq_woFungi_Perc_02um_ord <- ordinate(Phyloseq_woFungi_Perc_02um, "NMDS", "bray")
p1 = plot_ordination(Phyloseq_woFungi_Perc_02um, Phyloseq_woFungi_Perc_02um_ord, type="samples", color="PercFungi", shape="Year") +
  annotate("text", x=-1, y=1.8, label= sprintf("Stress = %s",as.character(round(Phyloseq_woFungi_Perc_02um_ord$stress, digits = 2))), size = 8) +
  scale_colour_gradient2()+
  geom_point(size=5) +
  scale_shape_manual(values=c(16, 15, 17))+
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
print(p1) # One sample goes really out, and seems to contain a lot of Fungi
ggplotly(p1)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/NMDS_WithoutFungiASVLevel_02um_ColorPercFungi.pdf", plot=p1, device = cairo_pdf(), width=15, height=14)
dev.off()

p1 = plot_ordination(Phyloseq_woFungi_Perc_02um, Phyloseq_woFungi_Perc_02um_ord, type="samples", color="Month", shape="Year")+
  annotate("text", x=-1, y=1.8, label= sprintf("Stress = %s",as.character(round(Phyloseq_woFungi_Perc_02um_ord$stress, digits = 2))), size = 8) +
  geom_point(size=5) +
  scale_shape_manual(values=c(16, 15, 17))+
  theme_bw() +
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
print(p1) # It seems that samples from the same months really group together, that is really clear.
ggplotly(p1)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/NMDS_WithoutFungiASVLevel_02um_ColorMonth.pdf", plot=p1, device = cairo_pdf(), width=15, height=14)
dev.off()

p1 = plot_ordination(Phyloseq_woFungi_Perc_02um, Phyloseq_woFungi_Perc_02um_ord, type="samples", color="Above5", shape="Year") +
  annotate("text", x=-1, y=1.8, label= sprintf("Stress = %s",as.character(round(Phyloseq_woFungi_Perc_02um_ord$stress, digits = 2))), size = 8) +
  geom_point(size=5) +
  scale_shape_manual(values=c(16, 15, 17))+
  theme_bw() +
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20))
print(p1) # It seems that samples from the same months really group together, that is really clear.
ggplotly(p1)
ggsave(file="/home/emilelaymand/Documents/Science/These/Manuscrit/PartiesManuscrit/Chapitre1/RawFigures/Supplementary/NMDS_WithoutFungiASVLevel_02um_ColorAbove5PercFungi.pdf", plot=p1, device = cairo_pdf(), width=15, height=14)
dev.off()

