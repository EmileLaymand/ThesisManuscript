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
library(lubridate)
library(lomb)

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

lsp_allspecies <- function(PhyObj, times = NULL, from = NULL, to = NULL, type = c("frequency", "period"),
                           ofac = 1, alpha = 0.01, normalize=c("standard","press")) {
  # This function calculates the lomb-scargale periodogram for all species in a phyloseq object
  
  # Melt the phyloseq object
  PhyObjMelt <- psmelt(PhyObj)
  
  # Calculate the number of days from the first day in the time series
  
  PhyObjMelt$DaysFromStart <- as.numeric(difftime(PhyObjMelt$Date_num, sort(PhyObjMelt$Date_num)[1], units = "days"))
  
  # Make the counters for the for loop
  counter <- 0
  MaxCounter <- length(unique(PhyObjMelt[["OTU"]]))
  
  # Make the empty output object 
  OutputDF <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(OutputDF) <-  c("ASV", "sig.level", "peak", "peak.at_1", "peak.at_2", "p.value", "normalize", "n", "type", "ofac", "alpha")
  
  # For every species name, calculate the lomb-scargale periodogram
  for (i in unique(PhyObjMelt[["OTU"]])) {
    
    # Calculate the lomb-scargale periodogram for species i
    PhyObjMelt_i <- PhyObjMelt[PhyObjMelt[["OTU"]] == i, c("DaysFromStart", "Abundance")]
    LSP_i <- lsp(PhyObjMelt_i, times = times, from = from, to = to, type = type,
                 ofac = ofac, alpha = alpha, normalize=normalize, plot = FALSE)
    
    # Put the output of lsp in the final object
    LSP_i_vector <- c(i, LSP_i$sig.level, LSP_i$peak, LSP_i$peak.at[1], LSP_i$peak.at[2], LSP_i$p.value, LSP_i$normalize, LSP_i$n, LSP_i$type, LSP_i$ofac, LSP_i$alpha)
    OutputDF <- rbind(OutputDF, LSP_i_vector)
    
    # Increase the counter and output progress
    counter <- counter + 1
    cat("\r", counter, " ASVs processed out of ", MaxCounter, sep = "")
  }
  
  # Change the names of the columns to be the right ones 
  colnames(OutputDF) <-  c("ASV", "sig.level", "peak", "peak.at_1", "peak.at_2", "p.value", "normalize", "n", "type", "ofac", "alpha")
  
  # Return the output object
  return(OutputDF)
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

# Select only 3um samples

FungiOTUPercPhyloseq_EukAll_3um <- subset_samples(FungiOTUPercPhyloseq_EukAll, Filter_Euk=="3µM")

# Merge the otu_table, the tax_table and the sample_data

FungiOTUPercPhyloseq_EukAll_3um_melt <- psmelt(FungiOTUPercPhyloseq_EukAll_3um)

# Take only ASV224

DFASV224 <- FungiOTUPercPhyloseq_EukAll_3um_melt[FungiOTUPercPhyloseq_EukAll_3um_melt[["OTU"]] == "ASV224", c("OTU", "Date_num", "Abundance")]
DFASV224$Days <- difftime(DFASV224$Date_num, sort(DFASV224$Date_num)[1], units = "days")
DFASV224$DaysNum <- as.numeric(DFASV224$Days)
DFASV224Lomb <- DFASV224[, c("DaysNum", "Abundance")]

# Make lsp

# LSPASV224 <- lsp(DFASV224Lomb, times = NULL, from = 2, to = 366, type = "frequency",
#     ofac = 1, alpha = 0.01, normalize="standard", plot = FALSE)
# summary.lsp(LSPASV224)

# Make lsp with sorted dataframe, in case this causes a miscalculation
DFASV224LombSorted <- DFASV224Lomb[order(DFASV224Lomb$DaysNum), ]
# LSPASV224Sorted <- lsp(DFASV224LombSorted, times = NULL, from = 2, to = 366, type = "frequency",
#                  ofac = 1, alpha = 0.01, normalize="standard", plot = FALSE)
LSPASV224Sorted <- lsp(DFASV224LombSorted, times = NULL, from = 2, to = 366, type = "period",
                       ofac = 1, alpha = 0.01, normalize="standard", plot = FALSE)
LSPASV224Sorted <- lsp(DFASV224LombSorted, times = NULL, from = 2, to = 366, type = "period",
                       ofac = 100, alpha = 0.01, normalize="standard", plot = FALSE)
LSPASV224Sorted <- lsp(DFASV224LombSorted, times = NULL, from = 2, to = 900, type = "period",
                       ofac = 1, alpha = 0.01, normalize="standard", plot = FALSE)
plot(LSPASV224Sorted)

# RandLSPASV224Sorted <- randlsp(repeats=1000, DFASV224LombSorted, from = 2, to = 366,
#                                type = "period", ofac = 10, alpha = 0.01, trace = TRUE)

# Test with the ibex and lynx datasets

data(ibex)
IbexLSP <- lsp(ibex[,2:3],ofac=5)
IbexLSP2 <- lsp(ibex$temp,times=ibex$hours,type='period',ofac=5)

data(lynx)
lynx.spec <- lsp(lynx,type='period',from=2,to=20,ofac=5)
summary(lynx.spec)

# Check if temperature is rythmic with a period of around 360 days to check lsp is working correctly

DFASVTemp <- FungiOTUPercPhyloseq_EukAll_3um_melt[FungiOTUPercPhyloseq_EukAll_3um_melt[["OTU"]] == "ASV224", c("OTU", "Date_num", "T")]
DFASVTemp$Days <- difftime(DFASVTemp$Date_num, sort(DFASVTemp$Date_num)[1], units = "days")
DFASVTemp$DaysNum <- as.numeric(DFASVTemp$Days)
DFASVTempLomb <- DFASVTemp[complete.cases(DFASVTemp), c("DaysNum", "T")]

LSPASVTemp <- lsp(DFASVTempLomb, times = NULL, from = 2, to = 900, type = "period",
                 ofac = 100, alpha = 0.01, normalize="standard", plot = FALSE)
plot(LSPASVTemp) # It works well!
summary(LSPASVTemp)

# Let's calculate the lomb-scargale periodogram to all species

LSPAllASVs <- lsp_allspecies(FungiOTUPercPhyloseq_EukAll_3um, times = NULL, from = 2, to = 900, type = "period",
                             ofac = 1, alpha = 0.01, normalize="standard")

LSPAllASVs$peak <- as.numeric(LSPAllASVs$peak)
LSPAllASVs[order(-LSPAllASVs$peak),]

LSPAllASVs_press <- lsp_allspecies(FungiOTUPercPhyloseq_EukAll_3um, times = NULL, from = 2, to = 900, type = "period",
                             ofac = 1, alpha = 0.01, normalize="press") # Does not work

LSPAllASVs$peak <- as.numeric(LSPAllASVs$peak)
LSPAllASVs[order(-LSPAllASVs$peak),]

# ASV8969 seems to be nearly periodic according to LSPAllASVs. 

DFASV8969 <- FungiOTUPercPhyloseq_EukAll_3um_melt[FungiOTUPercPhyloseq_EukAll_3um_melt[["OTU"]] == "ASV8969", c("OTU", "Date_num", "Abundance")]
DFASV8969$Days <- difftime(DFASV8969$Date_num, sort(DFASV8969$Date_num)[1], units = "days")
DFASV8969$DaysNum <- as.numeric(DFASV8969$Days)
DFASV8969Lomb <- DFASV8969[, c("DaysNum", "Abundance")]

# Make lsp

DFASV8969lsp <- lsp(DFASV8969Lomb, times = NULL, from = 2, to = 900, type = "period",
                 ofac = 1, alpha = 0.01, normalize="standard", plot = FALSE)
DFASV8969lsp <- lsp(DFASV8969Lomb, times = NULL, from = 2, to = 900, type = "period",
                 ofac = 1, alpha = 0.01, normalize="press", plot = FALSE)
plot(DFASV8969)
#abline(h = 10, col = "red")
summary.lsp(DFASV8969lsp)

# End of ASV8969 --> According to Max's standards, this should be rythmic

# Now, make the same calculation, but at the genus level, then family level, then order level, then class level, then Fungi all together

# Genus --> 1 genus significant

FungiOTUPercPhyloseq_EukAll_3um_genus <- tax_glom(FungiOTUPercPhyloseq_EukAll_3um, "Genus", NArm = FALSE)
LSPAllASVsGenus <- lsp_allspecies(FungiOTUPercPhyloseq_EukAll_3um_genus, times = NULL, from = 2, to = 900, type = "period",
                             ofac = 1, alpha = 0.01, normalize="standard")

LSPAllASVsGenus$peak <- as.numeric(LSPAllASVsGenus$peak)
LSPAllASVsGenus[order(-LSPAllASVsGenus$peak),]

# Family --> No family significant

FungiOTUPercPhyloseq_EukAll_3um_family <- tax_glom(FungiOTUPercPhyloseq_EukAll_3um, "Family", NArm = FALSE)
LSPAllASVsFamily <- lsp_allspecies(FungiOTUPercPhyloseq_EukAll_3um_family, times = NULL, from = 2, to = 900, type = "period",
                                  ofac = 1, alpha = 0.01, normalize="standard")

LSPAllASVsFamily$peak <- as.numeric(LSPAllASVsFamily$peak)
LSPAllASVsFamily[order(-LSPAllASVsFamily$peak),]

# Order --> No order significant

FungiOTUPercPhyloseq_EukAll_3um_order <- tax_glom(FungiOTUPercPhyloseq_EukAll_3um, "Order", NArm = FALSE)
LSPAllASVsOrder <- lsp_allspecies(FungiOTUPercPhyloseq_EukAll_3um_order, times = NULL, from = 2, to = 900, type = "period",
                                   ofac = 1, alpha = 0.01, normalize="standard")

LSPAllASVsOrder$peak <- as.numeric(LSPAllASVsOrder$peak)
LSPAllASVsOrder[order(-LSPAllASVsOrder$peak),]

# Class --> No class significant

FungiOTUPercPhyloseq_EukAll_3um_class <- tax_glom(FungiOTUPercPhyloseq_EukAll_3um, "Class", NArm = FALSE)
LSPAllASVsClass <- lsp_allspecies(FungiOTUPercPhyloseq_EukAll_3um_class, times = NULL, from = 2, to = 900, type = "period",
                                  ofac = 1, alpha = 0.01, normalize="standard")

LSPAllASVsClass$peak <- as.numeric(LSPAllASVsClass$peak)
LSPAllASVsClass[order(-LSPAllASVsClass$peak),]

# Fungi all together --> Fungi all together are not significant

FungiOTUPercPhyloseq_EukAll_3um_division <- tax_glom(FungiOTUPercPhyloseq_EukAll_3um, "Division", NArm = FALSE)
LSPAllASVsDivision <- lsp_allspecies(FungiOTUPercPhyloseq_EukAll_3um_division, times = NULL, from = 2, to = 900, type = "period",
                                  ofac = 1, alpha = 0.01, normalize="standard")

LSPAllASVsDivision$peak <- as.numeric(LSPAllASVsDivision$peak)
LSPAllASVsDivision[order(-LSPAllASVsDivision$peak),]
