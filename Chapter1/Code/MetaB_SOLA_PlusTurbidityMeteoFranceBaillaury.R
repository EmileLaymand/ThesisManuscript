##################################################################################
# Script for data analysis of the SOLA metabarcoding dataset                     #
# 1st section: Eukaryotes 0.2 micrometers (ADNr 18S ? What region)               #
# 2nd section: Eukaryotes 3 micrometers (ADNr 18S ? What region)                 #
# 03/02/2021                                                                     #
# Author: Emile Laymand                                                          #
# Affilation: - Atelier de Bioinformatique, Institut de Systematique, Evolution, # 
#               Biodiversite, Museum National d'Histoire Naturelle, Paris,       #
#               France                                                           #
#             - Laboratoire d'Oceanalogie Microbienne, Observatoire              #
#               Oceanologique de Banyuls sur Mer, Sorbonne Universite/CNRS,      #
#               Banyuls sur Mer, France                                          #
#             - Departement de Geosciences, Ecole Normale Superieure, Paris,     #
#               France                                                           #
##################################################################################

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

#==============#
#     Flags    #
#==============#

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

#---------------------------------
# Step 2 : First look at the data
#---------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Here, my goal is to see how many Fungi I have in the samples (as absolute and relative abundance)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Create a color palette so that adjacent colors are different enough to distinguish them

cc <- c(palette("Set3"), palette("Dark2"), palette("Accent"), palette("Set1"))

# Plot absolute abundances

#plot_bar(Phyloseq_Euk02, x="Name_Euk", fill="Division")  + scale_fill_manual(values = cc) # Have an idea of present Kingdoms
Phyloseq_Euk02_glommed_div <- tax_glom(Phyloseq_Euk02, taxrank = "Division")
plot_bar(Phyloseq_Euk02_glommed_div, x="Name_Euk02", fill ="Division") + scale_fill_manual(values = cc) # Same but without delineation of ASVs

# Plot relative abundances

Phyloseq_Euk02_perc <- taxa_percentize(Phyloseq_Euk02, TaxLevel = "Division")
plot_bar(Phyloseq_Euk02_perc, x="Name_Euk02", fill="Division") + scale_fill_manual(values = cc) # Have an idea of present Kingdoms

#--------------------------------------------------------------------------------------------------------------------------------

# Take only Fungi 

FungiPhyloseq_Euk02 <- subset_taxa(Phyloseq_Euk02, Division == "Fungi")
FungiPhyloseq_Euk02 <- prune_taxa(taxa_sums(FungiPhyloseq_Euk02) > 0, FungiPhyloseq_Euk02)

# Plot absolute abundances

plot_bar(FungiPhyloseq_Euk02, x="Name_Euk02", fill="Class")  + scale_fill_manual(values = cc) # Have an idea of present Kingdoms

# Plot relative abundances

FungiPhyloseq_Euk02_perc <- taxa_percentize(FungiPhyloseq_Euk02, TaxLevel = "Class")
FungiPhyloseq_Euk02_perc <- prune_samples(sample_sums(FungiPhyloseq_Euk02_perc) > 0, FungiPhyloseq_Euk02_perc)
plot_bar(FungiPhyloseq_Euk02_perc, x="Name_Euk02", fill="Class") + scale_fill_manual(values = cc) # Have an idea of present Kingdoms

# Explore data

## Remove empty samples (those not containing any fungal read)

FungiPhyloseq_Euk02_samp_pruned <- prune_samples(sample_sums(FungiPhyloseq_Euk02) > 0, FungiPhyloseq_Euk02)

#ord_explore(FungiPhyloseq_Euk02_samp_pruned)

# Make a Phyloseq with a modified tax table, with a dummy variable called "OTU"

FungiPhyloseq_Euk02_OTU <- FungiPhyloseq_Euk02 
tax_table(FungiPhyloseq_Euk02_OTU) <- cbind(tax_table(FungiPhyloseq_Euk02_OTU), ASV = row.names(tax_table(FungiPhyloseq_Euk02_OTU)))

## Plot absolute abundance
# Not working

plot_bar(FungiPhyloseq_Euk02_OTU, x="Name_Euk02", fill="ASV") + scale_fill_manual(values = cc) # Have an idea of present Kingdoms

#--------------------------------------------------------------------------------------------------------------------------------

# Take only Stramenopiles 

StramenopilesPhyloseq_Euk02 <- subset_taxa(Phyloseq_Euk02, Division == "Stramenopiles_X")
StramenopilesPhyloseq_Euk02 <- prune_taxa(taxa_sums(StramenopilesPhyloseq_Euk02) > 0, StramenopilesPhyloseq_Euk02)

# Plot absolute abundances

plot_bar(StramenopilesPhyloseq_Euk02, x="Name_Euk02", fill="Class")  + scale_fill_manual(values = cc) # Have an idea of present Kingdoms

# Plot relative abundances

StramenopilesPhyloseq_Euk02_perc <- taxa_percentize(StramenopilesPhyloseq_Euk02, TaxLevel = "Class")
StramenopilesPhyloseq_Euk02_perc <- prune_samples(sample_sums(StramenopilesPhyloseq_Euk02_perc) > 0, StramenopilesPhyloseq_Euk02_perc)
plot_bar(StramenopilesPhyloseq_Euk02_perc, x="Name_Euk02", fill="Class") + scale_fill_manual(values = cc) # Have an idea of present Kingdoms

# Explore data

## Remove empty samples (those not containing any fungal read)

StramenopilesPhyloseq_Euk02_samp_pruned <- prune_samples(sample_sums(StramenopilesPhyloseq_Euk02) > 0, StramenopilesPhyloseq_Euk02)

#ord_explore(StramenopilesPhyloseq_Euk02_samp_pruned)

# Make a Phyloseq with a modified tax table, with a dummy variable called "OTU"

StramenopilesPhyloseq_Euk02_OTU <- StramenopilesPhyloseq_Euk02 
tax_table(StramenopilesPhyloseq_Euk02_OTU) <- cbind(tax_table(StramenopilesPhyloseq_Euk02_OTU), ASV = row.names(tax_table(StramenopilesPhyloseq_Euk02_OTU)))




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

#---------------------------------
# Step 2 : First look at the data
#---------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Here, my goal is to see how many Fungi I have in the samples (as absolute and relative abundance)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Create a color palette so that adjacent colors are different enough to distinguish them

cc <- c(palette("Set3"), palette("Dark2"), palette("Accent"), palette("Set1"))

# Plot absolute abundances

#plot_bar(Phyloseq_Euk3, x="Name_Euk", fill="Division")  + scale_fill_manual(values = cc) # Have an idea of present Kingdoms
Phyloseq_Euk3_glommed_div <- tax_glom(Phyloseq_Euk3, taxrank = "Division")
plot_bar(Phyloseq_Euk3_glommed_div, x="Name_Euk", fill ="Division") + scale_fill_manual(values = cc) # Same but without delineation of ASVs

# Plot relative abundances

Phyloseq_Euk3_perc <- taxa_percentize(Phyloseq_Euk3, TaxLevel = "Division")
plot_bar(Phyloseq_Euk3_perc, x="Name_Euk", fill="Division") + scale_fill_manual(values = cc) # Have an idea of present Kingdoms

#--------------------------------------------------------------------------------------------------------------------------------

# Take only Fungi 

FungiPhyloseq_Euk3 <- subset_taxa(Phyloseq_Euk3, Division == "Fungi")
FungiPhyloseq_Euk3 <- prune_taxa(taxa_sums(FungiPhyloseq_Euk3) > 0, FungiPhyloseq_Euk3)

# Plot absolute abundances

plot_bar(FungiPhyloseq_Euk3, x="Name_Euk", fill="Class")  + scale_fill_manual(values = cc) # Have an idea of present Kingdoms

# Absolute number of fungal reads per sample

rowSums(otu_table(FungiPhyloseq_Euk3))

# Plot relative abundances

FungiPhyloseq_Euk3_perc <- taxa_percentize(FungiPhyloseq_Euk3, TaxLevel = "Class")
plot_bar(FungiPhyloseq_Euk3_perc, x="Name_Euk", fill="Class") + scale_fill_manual(values = cc) # Have an idea of present Kingdoms

# Plot relative abundances only for samples that have more than 100 fungal reads

FungiPhyloseq_Euk3_SortedRawCounts <- sort(colSums(otu_table(FungiPhyloseq_Euk3)), decreasing = TRUE)
FungiPhyloseq_Euk3_more1000reads <- names(FungiPhyloseq_Euk3_SortedRawCounts[FungiPhyloseq_Euk3_SortedRawCounts >= 1000])
FungiPhyloseq_Euk3_perc_more1000reads <- subset_samples(FungiPhyloseq_Euk3_perc, Name_Euk %in% FungiPhyloseq_Euk3_more1000reads)
plot_bar(FungiPhyloseq_Euk3_perc_more1000reads, x="Name_Euk", fill="Class") + scale_fill_manual(values = cc) # Have an idea of present Kingdoms

# Plot relative abundances

#ord_explore(FungiPhyloseq_Euk3)

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

EukaAll_OTU_Seq <- list(Euka02_OTU_taxo_Seq_2merge_wotax, Euka3_OTU_taxo_Seq_2merge_wotax) %>% reduce(full_join) # Hard check # 20240109 Seems OK

# ASVs that are in one filter size and not the other will create NAs in the latter. Replace all NAs by 0

EukaAll_OTU_Seq <- EukaAll_OTU_Seq %>% replace(is.na(.), 0)

# That that merging was good (no modification in taxonomy, not twice the same sequence with two different taxonomies)

length(unique(Euka02_OTU_taxo_Seq_2merge[["X"]]))

length(Euka02_OTU_taxo_Seq_2merge[["X"]])

length(unique(Euka3_OTU_taxo_Seq_2merge[["X"]]))
length(Euka3_OTU_taxo_Seq_2merge[["X"]])
sum((unique(Euka02_OTU_taxo_Seq_2merge[["X"]]) %in% unique(Euka3_OTU_taxo_Seq_2merge[["X"]])), na.rm = TRUE)
length(EukaAll_OTU_Seq[["X"]])
length(unique(EukaAll_OTU_Seq[["X"]])) # 20240109 Seems OK

# For the taxa that would be duplicated keep only the taxonomy of the Eukaryotes 3 um (this is arbitrary)

## Take all sequences from 3 um

Euka3_taxo_seq_2merge <- Euka3_OTU_taxo_Seq_2merge[, c("X", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

## Take sequences that are in 2 um but not in 3 um

Euka02_taxo_seq_2merge <- Euka02_OTU_taxo_Seq_2merge[!(Euka02_OTU_taxo_Seq_2merge[["X"]] %in% Euka3_taxo_seq_2merge[["X"]]), c("X", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

## rbind these two dataframes to have the final taxonomy

EukaAll_taxo_Seq <- rbind(Euka3_taxo_seq_2merge, Euka02_taxo_seq_2merge)

## Check dimensions

length(unique(EukaAll_taxo_Seq[["X"]]))
length(EukaAll_taxo_Seq[["X"]]) # 20240109 Seems OK

# Merge the synthetic taxonomy and the OTU dataframes to make sure everything is in the right order

EukaAll_OTU_taxo_Seq <- list(EukaAll_OTU_Seq, EukaAll_taxo_Seq) %>% reduce(full_join) # Hard check

# Check dimensions

dim(EukaAll_OTU_taxo_Seq) # 20240109 Seems OK

# Check for some sequences that the taxonomy is good

## 1st check # 20240109 Seems OK

EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["X"]] == "AGCTCCAATAGCGTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCGGTTGAGAACGGCCGGTCCGCCGTTTGGTGTGCACTGGCTGGTTTCAACTTCCTGTAGAGGACGCGCTCTGGCTTCATCGCTGGACGCGGAGTCTACGTGGTTACTTTGAAAAAATTAGAGTGTTCAAAGCGGGCTTACGCTTGAATATTTCAGCATGGAATAACACTATAGGACTCCTGTCCTATTTCGTTGGTCTCGGGACGGGAGTAATGATTAAGAGGAACAGTTGGGGGCATTCGTATTTCATTGTCAGAGGTGAAATTCTTGGATTTATGAAAGACGAACTTCTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka02_OTU_taxo_Seq_2merge[Euka02_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCGGTTGAGAACGGCCGGTCCGCCGTTTGGTGTGCACTGGCTGGTTTCAACTTCCTGTAGAGGACGCGCTCTGGCTTCATCGCTGGACGCGGAGTCTACGTGGTTACTTTGAAAAAATTAGAGTGTTCAAAGCGGGCTTACGCTTGAATATTTCAGCATGGAATAACACTATAGGACTCCTGTCCTATTTCGTTGGTCTCGGGACGGGAGTAATGATTAAGAGGAACAGTTGGGGGCATTCGTATTTCATTGTCAGAGGTGAAATTCTTGGATTTATGAAAGACGAACTTCTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka3_OTU_taxo_Seq_2merge[Euka3_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCGGTTGAGAACGGCCGGTCCGCCGTTTGGTGTGCACTGGCTGGTTTCAACTTCCTGTAGAGGACGCGCTCTGGCTTCATCGCTGGACGCGGAGTCTACGTGGTTACTTTGAAAAAATTAGAGTGTTCAAAGCGGGCTTACGCTTGAATATTTCAGCATGGAATAACACTATAGGACTCCTGTCCTATTTCGTTGGTCTCGGGACGGGAGTAATGATTAAGAGGAACAGTTGGGGGCATTCGTATTTCATTGTCAGAGGTGAAATTCTTGGATTTATGAAAGACGAACTTCTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

## 2nd check (# 20240109 Seems OK in 02 but not in 3)

EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATTTCTGTTGAGGACGGCCGGTCCGCATTATGTGCGTGTATCTGGTTCGGCCTTGGCATCCTCCAGGGGAACGTTCCTGCGCTTCGCTGCGTGGGACGGTATTCTGGACTTTTACTTTGAGGAAATTAGAGTGTTCACGGCAGGCAGTCGCCTTGAATATATTAGCATGGAATAATAATATAGGACCTTGGTTCTATTTTGTTGGTTTCTAGAACTTAGGTAATGATTAATAGGGACAATTGGGGGCATTTGTATTAACACGTCAGAGGTGAAATTCTTGGATTGTGTTACGACAAACTACTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka02_OTU_taxo_Seq_2merge[Euka02_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATTTCTGTTGAGGACGGCCGGTCCGCATTATGTGCGTGTATCTGGTTCGGCCTTGGCATCCTCCAGGGGAACGTTCCTGCGCTTCGCTGCGTGGGACGGTATTCTGGACTTTTACTTTGAGGAAATTAGAGTGTTCACGGCAGGCAGTCGCCTTGAATATATTAGCATGGAATAATAATATAGGACCTTGGTTCTATTTTGTTGGTTTCTAGAACTTAGGTAATGATTAATAGGGACAATTGGGGGCATTTGTATTAACACGTCAGAGGTGAAATTCTTGGATTGTGTTACGACAAACTACTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka3_OTU_taxo_Seq_2merge[Euka3_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATTTCTGTTGAGGACGGCCGGTCCGCATTATGTGCGTGTATCTGGTTCGGCCTTGGCATCCTCCAGGGGAACGTTCCTGCGCTTCGCTGCGTGGGACGGTATTCTGGACTTTTACTTTGAGGAAATTAGAGTGTTCACGGCAGGCAGTCGCCTTGAATATATTAGCATGGAATAATAATATAGGACCTTGGTTCTATTTTGTTGGTTTCTAGAACTTAGGTAATGATTAATAGGGACAATTGGGGGCATTTGTATTAACACGTCAGAGGTGAAATTCTTGGATTGTGTTACGACAAACTACTGCGAAAGCATTTGCCAAGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

# 3rd check in 02 but not in 3 # 20240109 Seems OK

EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGAGTTCTGCCTGGTGACGCCCGTCCGCCCAAGTGGTGTGTACAGGGTGTACATCTGGCCCTTTCAAGGGGAACGTATCTGCACTTAACTGTGCGGTGCGAGATCCTTGACTTTTACTTTGAGGAAATCAGAGTGTTCCAAGCAGGCTCTCGTCGTGCATGTTTCAGCATGGAATAATAGCATTGGACCTCGTCTCTCAGCTGTTGGTTGCAAGAAGCGAGGTAATGATGAAGAGGGATAGTTGGGGGCATTCGTATTTAACTGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGGACTACTGCGAAAGCATCTGCCATGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka02_OTU_taxo_Seq_2merge[Euka02_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGAGTTCTGCCTGGTGACGCCCGTCCGCCCAAGTGGTGTGTACAGGGTGTACATCTGGCCCTTTCAAGGGGAACGTATCTGCACTTAACTGTGCGGTGCGAGATCCTTGACTTTTACTTTGAGGAAATCAGAGTGTTCCAAGCAGGCTCTCGTCGTGCATGTTTCAGCATGGAATAATAGCATTGGACCTCGTCTCTCAGCTGTTGGTTGCAAGAAGCGAGGTAATGATGAAGAGGGATAGTTGGGGGCATTCGTATTTAACTGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGGACTACTGCGAAAGCATCTGCCATGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka3_OTU_taxo_Seq_2merge[Euka3_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGAGTTCTGCCTGGTGACGCCCGTCCGCCCAAGTGGTGTGTACAGGGTGTACATCTGGCCCTTTCAAGGGGAACGTATCTGCACTTAACTGTGCGGTGCGAGATCCTTGACTTTTACTTTGAGGAAATCAGAGTGTTCCAAGCAGGCTCTCGTCGTGCATGTTTCAGCATGGAATAATAGCATTGGACCTCGTCTCTCAGCTGTTGGTTGCAAGAAGCGAGGTAATGATGAAGAGGGATAGTTGGGGGCATTCGTATTTAACTGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGGACTACTGCGAAAGCATCTGCCATGGATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

# 4th check in both but different taxonomy # 20240109 Seems OK

Test_tax <- list(Euka02_OTU_taxo_Seq_2merge, Euka3_OTU_taxo_Seq_2merge) %>% reduce(full_join) # Hard check 
Test_tax[duplicated(Test_tax[["X"]]),]

EukaAll_OTU_taxo_Seq[EukaAll_OTU_taxo_Seq[["X"]] == "AGCTCCAATAGTGTATATTAAAGTTGCTGCGGTTAAAATGCTCGTAGTTTAATTTCTGCTGAGGGTAATCGGTCCACCCACTGGGTGAGTAATTGTTTGCCCTTTGCATTTTGTGAACAATGCATCTGCACTTGATTGTGCTGTGTGTCCTGTTCATGCAATTTACTTTGAGGAAATTAGAGTGTTTCGAGCATGCATATGCACCGAGCACATTAGCATGGAATAATTTACACTGATCGTCGCTGTATTTGTTGGTTTTTAGGGCTTCGATAATGACTGATAGGGATAATTGGTGGTATTCGCATTTAATAGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGAACTGCTGCGAAAGCATTTGCCAACTATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka02_OTU_taxo_Seq_2merge[Euka02_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGTGTATATTAAAGTTGCTGCGGTTAAAATGCTCGTAGTTTAATTTCTGCTGAGGGTAATCGGTCCACCCACTGGGTGAGTAATTGTTTGCCCTTTGCATTTTGTGAACAATGCATCTGCACTTGATTGTGCTGTGTGTCCTGTTCATGCAATTTACTTTGAGGAAATTAGAGTGTTTCGAGCATGCATATGCACCGAGCACATTAGCATGGAATAATTTACACTGATCGTCGCTGTATTTGTTGGTTTTTAGGGCTTCGATAATGACTGATAGGGATAATTGGTGGTATTCGCATTTAATAGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGAACTGCTGCGAAAGCATTTGCCAACTATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

Euka3_OTU_taxo_Seq_2merge[Euka3_OTU_taxo_Seq_2merge[["X"]] == "AGCTCCAATAGTGTATATTAAAGTTGCTGCGGTTAAAATGCTCGTAGTTTAATTTCTGCTGAGGGTAATCGGTCCACCCACTGGGTGAGTAATTGTTTGCCCTTTGCATTTTGTGAACAATGCATCTGCACTTGATTGTGCTGTGTGTCCTGTTCATGCAATTTACTTTGAGGAAATTAGAGTGTTTCGAGCATGCATATGCACCGAGCACATTAGCATGGAATAATTTACACTGATCGTCGCTGTATTTGTTGGTTTTTAGGGCTTCGATAATGACTGATAGGGATAATTGGTGGTATTCGCATTTAATAGTCAGAGGTGAAATTCTTGGATTTGTTAAAGACGAACTGCTGCGAAAGCATTTGCCAACTATGTTT", c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")]

# Separate again the OTU, sequences and taxonomy into 3 different dataframes

# Create a dummy OTU name for each sequence

EukaAll_OTU_taxo_Seq$ASV <- paste("ASV", c(1:length(row.names(EukaAll_OTU_taxo_Seq))), sep = "")
length(which(is.na(EukaAll_OTU_taxo_Seq$ASV) == TRUE)) # Must be 0 # 20240109 Seems OK
which(duplicated(EukaAll_OTU_taxo_Seq$ASV) == TRUE) # Must be integer(0) # 20240109 Seems OK

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
EukaAll_taxo_mat <- as.matrix(EukaAll_taxo) # 20240109 Seems OK

row.names(EukaAll_OTU) <- EukaAll_OTU$ASV
EukaAll_OTU$ASV <- NULL

# Import the file containing metadata (sampling date, chlorophyll concentration, depth, GPS coordinates, ...)

EukaAll_Metadata <- read.csv("/home/emilelaymand/Documents/Science/These/SOLA_Galand/SOLA/wetransfer_sola_2021-10-18_1242/Equivalence_name_date_3um_concatenated.csv")

# Crate the right sample names, then make the sample names the row names (necessary to build the phyloseq object)

EukaAll_Metadata[str_sub(EukaAll_Metadata[["Name_Euk"]], -3, -1) == "Pro", "Name_Euk"] <- paste(str_sub(EukaAll_Metadata[str_sub(EukaAll_Metadata[["Name_Euk"]], -3, -1) == "Pro", "Name_Euk"], 1, -4), "Euk", sep = "") # 20240109 Seems OK

row.names(EukaAll_Metadata) <- EukaAll_Metadata$Name_Euk

# Add one more column with a date that is actually considered as a date

EukaAll_Metadata$Date_num <- as.Date(EukaAll_Metadata$Date_Euk, "%d_%m_%Y")

# Import a second metadata dataframe, that contains details on nutrients, chl a, etc.

EukaAll_Metadata_nut <- read.csv("/home/emilelaymand/Documents/Science/These/SOLA_Galand/SOLA/wetransfer_sola_2021-10-18_1242/climato_2007-2017_modified.csv")

EukaAll_Metadata_nut$Date_num <- as.Date(EukaAll_Metadata_nut$DATE, "%d/%m/%Y")

# Merge both dataframes, then keep only rows that are present in the first dataframe 

EukaAll_Metadata_merged <- list(EukaAll_Metadata, EukaAll_Metadata_nut) %>% reduce(full_join) # Hard check # 20240109 Seems OK

EukaAll_Metadata_merged_red <- EukaAll_Metadata_merged[EukaAll_Metadata_merged[["Date_num"]] %in% unique(EukaAll_Metadata[["Date_num"]]), ]

# Check dimensions

dim(EukaAll_Metadata)

dim(EukaAll_Metadata_merged_red)

EukaAll_Metadata_merged_red[EukaAll_Metadata_merged_red[["Date_num"]] == "2015-01-05", ] # Metadata must be identical for the two rows # 20240109 Seems OK

identical(EukaAll_Metadata[duplicated(EukaAll_Metadata[["Date_num"]]),"Date_num"], EukaAll_Metadata_merged_red[duplicated(EukaAll_Metadata_merged_red[["Date_num"]]),"Date_num"]) # Must be identical # 20240109 Seems OK

# Replace EukaAll_Metadata with the new, more complete EukaAll_Metadata_merged_red

EukaAll_Metadata <- EukaAll_Metadata_merged_red

row.names(EukaAll_Metadata) <- EukaAll_Metadata$Name_Euk

# ########### NEW IN THIS FILE ###################
# 
# # Import the files for turbidity, MeteoFrance (Precipitations and wind) and the Baillaury river
# 
# # Import MeteoFrance data 
# # RR: HAUTEUR DE PRECIPITATIONS QUOTIDIENNE en MILLIMETRES ET 1/10
# # DRR: DUREE DES PRECIPITATIONS QUOTIDIENNES en MINUTES
# # FF2M: MOYENNE DES VITESSES DU VENT A 2 METRES QUOTIDIENNE en M/S ET 1/10
# # FFM: MOYENNE DES VITESSES DU VENT A 10M QUOTIDIENNE en M/S ET 1/10
# 
# MeteoFranceMetadata <- read.csv("MeteoFrance/CapBearAPIClimatologie/CapBearConcatenated20132017.csv", sep = ",")
# MeteoFranceMetadata <- MeteoFranceMetadata[, c("Date_num", "RR", "FFM")]
# 
# # Import turbidity data
# 
# TurbidityMetadata <- read.csv("CTDFromPaulMerged/Turbidity3m20mMerged.csv", sep = ",")
# 
# # Import Baillaury data # Too bad for now, let's not use it
# 
# # TO BE ADDED
# 
# # Merge these new dataframes with the original metadataDF
# 
# EukaAll_MetadataMeteo <- merge(EukaAll_Metadata, MeteoFranceMetadata, by = "Date_num", all = TRUE)
# EukaAll_MetadataMeteoTurb <- merge(EukaAll_MetadataMeteo, MeteoFranceMetadata, by = "Date_num", all = TRUE)




# Merge the data into a phyloseq object

OTU_All <- otu_table(EukaAll_OTU, taxa_are_rows = TRUE)

TAX_All <- tax_table(EukaAll_taxo_mat)

META_All <- sample_data(EukaAll_Metadata)

Phyloseq_EukAll <- phyloseq(OTU_All, TAX_All, META_All) 

# Check for ASV 500 (to compare to the raw .csv files) # 20240109 Seems OK

Seq_ASV500 <- EukaAll_Seq[EukaAll_Seq[["ASV"]] == "ASV500", "X"]
tax_table(Phyloseq_EukAll)["ASV500",]
otu_table(Phyloseq_EukAll)["ASV500",]
Euka3_OTU_taxo_Seq_2merge[Euka3_OTU_taxo_Seq_2merge[["X"]] == Seq_ASV500, ]
Euka02_OTU_taxo_Seq_2merge[Euka02_OTU_taxo_Seq_2merge[["X"]] == Seq_ASV500, ]

# Check 2: present in both but different taxonomy # 20240109 Seems OK

Seq_ASV346 <- EukaAll_Seq[EukaAll_Seq[["ASV"]] == "ASV346", "X"]
tax_table(Phyloseq_EukAll)["ASV346",]
otu_table(Phyloseq_EukAll)["ASV346",]
Euka3_OTU_taxo_Seq_2merge[Euka3_OTU_taxo_Seq_2merge[["X"]] == Seq_ASV346, ]
Euka02_OTU_taxo_Seq_2merge[Euka02_OTU_taxo_Seq_2merge[["X"]] == Seq_ASV346, ]

# Check 3: present in one dataset, but not the other # 20240109 Seems OK

Seq_ASV1420 <- EukaAll_Seq[EukaAll_Seq[["ASV"]] == "ASV1420", "X"]
tax_table(Phyloseq_EukAll)["ASV1420",]
otu_table(Phyloseq_EukAll)["ASV1420",]
Euka3_OTU_taxo_Seq_2merge[Euka3_OTU_taxo_Seq_2merge[["X"]] == Seq_ASV1420, ]
Euka02_OTU_taxo_Seq_2merge[Euka02_OTU_taxo_Seq_2merge[["X"]] == Seq_ASV1420, ]

#---------------------------------
# Step 2 : First look at the data
#---------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Here, my goal is to see how many Fungi I have in the samples (as absolute and relative abundance)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Create a color palette so that adjacent colors are different enough to distinguish them

cc <- c(palette("Set3"), palette("Dark2"), palette("Accent"), palette("Set1"))

# Plot absolute abundances

#plot_bar(Phyloseq_Euk3, x="Name_Euk", fill="Division")  + scale_fill_manual(values = cc) # Have an idea of present Kingdoms
Phyloseq_EukAll_glommed_div <- tax_glom(Phyloseq_EukAll, taxrank = "Division")
plot_bar(Phyloseq_EukAll_glommed_div, x="Name_Euk", fill ="Division") + 
  facet_grid(~Filter_Euk, scales = "free") +
  scale_fill_manual(values = cc) # Same but without delineation of ASVs

plot_bar(Phyloseq_EukAll_glommed_div, x="Date_num", fill ="Division") + 
  facet_grid(~Filter_Euk, scales = "free") +
  scale_fill_manual(values = cc) # Same but without delineation of ASVs and date as x label

# Plot relative abundances

Phyloseq_EukAll_perc <- taxa_percentize(Phyloseq_EukAll, TaxLevel = "Division")
plot_bar(Phyloseq_EukAll_perc, x="Name_Euk", fill="Division") + 
  facet_grid(~Filter_Euk, scales = "free") +
  scale_fill_manual(values = cc) # Have an idea of present Kingdoms

plot_bar(Phyloseq_EukAll_perc, x="Date_num", fill="Division") + 
  facet_grid(~Filter_Euk, scales = "free") +
  scale_fill_manual(values = cc) # Have an idea of present Kingdoms

plot_bar_2(Phyloseq_EukAll_perc, x="Date_num", fill="Division") + 
  facet_grid(~Filter_Euk, scales = "free") +
  scale_fill_manual(values = cc) # Have an idea of present Kingdoms)


# WARNING: stange behaviour with filling -> Is the value under or above?
ggplot(psmelt(Phyloseq_EukAll_perc), aes(x = Date_num, y = Abundance, fill = Division)) +
  facet_grid(~Filter_Euk, scales = "free") +
  geom_area(position = "stack") +
  geom_point() +
  scale_fill_manual(values = cc) # Have an idea of present Kingdoms)



#--------------------------------------------------------------------------------------------------------------------------------

# Take only Fungi 

FungiPhyloseq_EukAll <- subset_taxa(Phyloseq_EukAll, Division == "Fungi")
FungiPhyloseq_EukAll <- prune_taxa(taxa_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)
FungiPhyloseq_EukAll <- prune_samples(sample_sums(FungiPhyloseq_EukAll) > 0, FungiPhyloseq_EukAll)

# Plot absolute abundances

plot_bar(FungiPhyloseq_EukAll, x="Name_Euk", fill="Class")  + 
  facet_grid(~Filter_Euk, scales = "free") +
  scale_fill_manual(values = cc) # Have an idea of present Kingdoms

# Plot relative abundances

FungiPhyloseq_EukAll_perc <- taxa_percentize(FungiPhyloseq_EukAll, TaxLevel = "Class")
plot_bar(FungiPhyloseq_EukAll_perc, x="Name_Euk", fill="Class") + 
  facet_grid(~Filter_Euk, scales = "free") +
  scale_fill_manual(values = cc) # Have an idea of present Kingdoms

# Plot relative abundances

#ord_explore(FungiPhyloseq_EukAll)

# Check some ASVs

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV324",]



tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV224",]
# Kingdom     Supergroup     Division Class        Order            Family            Genus           Species                 
# ASV224 "Eukaryota" "Opisthokonta" "Fungi"  "Ascomycota" "Pezizomycotina" "Sordariomycetes" "Lecanicillium" "Lecanicillium_saksenae"
# Fugus used for its anti-nematodes and anti-insects properties

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV275",]
# Kingdom     Supergroup     Division Class        Order              Family              Genus     Species               
# ASV275 "Eukaryota" "Opisthokonta" "Fungi"  "Ascomycota" "Saccharomycotina" "Saccharomycetales" "Candida" "Candida_parapsilosis"
# Pathogenic yeast of skin for human and animals, including insects. Has been found in soils as well

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV409",]
# Kingdom     Supergroup     Division Class           Order                Family              Genus        Species               
# ASV409 "Eukaryota" "Opisthokonta" "Fungi"  "Basidiomycota" "Ustilaginomycotina" "Exobasidiomycetes" "Malassezia" "Malassezia_restricta"
# Pathogenic and skin yeast of human

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV1613",]
# Kingdom     Supergroup     Division Class        Order            Family            Genus Species
# ASV1613 "Eukaryota" "Opisthokonta" "Fungi"  "Ascomycota" "Pezizomycotina" "Sordariomycetes" NA    NA    

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV1016",]
# Kingdom     Supergroup     Division Class        Order            Family          Genus                 Species                  
# ASV1016 "Eukaryota" "Opisthokonta" "Fungi"  "Ascomycota" "Pezizomycotina" "Leotiomycetes" "Rhexocercosporidium" "Rhexocercosporidium_sp."
# Endophytic and plant parasite. Soil Fungus.

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV166",]
#Kingdom     Supergroup     Division Class          Order            Family             Genus               Species                
#ASV166 "Eukaryota" "Opisthokonta" "Fungi"  "Cryptomycota" "Cryptomycotina" "Cryptomycotina_X" "Cryptomycotina_XX" "Cryptomycotina_XX_sp."

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV1724",]
#Kingdom     Supergroup     Division Class        Order            Family            Genus        Species               
#ASV1724 "Eukaryota" "Opisthokonta" "Fungi"  "Ascomycota" "Pezizomycotina" "Dothideomycetes" "Ochroconis" "Ochroconis_mirabilis"

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV2336",]
#Kingdom     Supergroup     Division Class        Order            Family           Genus         Species                     
#ASV2336 "Eukaryota" "Opisthokonta" "Fungi"  "Ascomycota" "Pezizomycotina" "Eurotiomycetes" "Penicillium" "Penicillium_brevicompactum"

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV2777",]
#Kingdom     Supergroup     Division Class        Order            Family            Genus        Species                
#ASV2777 "Eukaryota" "Opisthokonta" "Fungi"  "Ascomycota" "Pezizomycotina" "Dothideomycetes" "Davidiella" "Davidiella_macrospora"

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV324",]
#Kingdom     Supergroup     Division Class             Order               Family             Genus            Species             
#ASV324 "Eukaryota" "Opisthokonta" "Fungi"  "Chytridiomycota" "Chytridiomycotina" "Chytridiomycetes" "Rhyzophidiales" "Rhyzophidiales_sp."

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV367",]
#Kingdom     Supergroup     Division Class             Order               Family             Genus              Species               
#ASV367 "Eukaryota" "Opisthokonta" "Fungi"  "Chytridiomycota" "Chytridiomycotina" "Chytridiomycetes" "Rhyzophidiales_X" "Rhyzophidiales_X_sp."

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV409",]
#Kingdom     Supergroup     Division Class           Order                Family              Genus        Species               
#ASV409 "Eukaryota" "Opisthokonta" "Fungi"  "Basidiomycota" "Ustilaginomycotina" "Exobasidiomycetes" "Malassezia" "Malassezia_restricta"

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV420",]
#Kingdom     Supergroup     Division Class             Order               Family             Genus                Species                 
#ASV420 "Eukaryota" "Opisthokonta" "Fungi"  "Chytridiomycota" "Chytridiomycotina" "Chytridiomycetes" "Chytridiomycetes_X" "Chytridiomycetes_X_sp."

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV5211",]
#Kingdom     Supergroup     Division Class        Order            Family            Genus          Species           
#ASV5211 "Eukaryota" "Opisthokonta" "Fungi"  "Ascomycota" "Pezizomycotina" "Dothideomycetes" "Cladosporium" "Cladosporium_sp."

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV553",]
#Kingdom     Supergroup     Division Class             Order               Family             Genus              Species               
#ASV553 "Eukaryota" "Opisthokonta" "Fungi"  "Chytridiomycota" "Chytridiomycotina" "Chytridiomycetes" "Rhyzophidiales_X" "Rhyzophidiales_X_sp."

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV594",]
#Kingdom     Supergroup     Division Class Order Family Genus Species
#ASV594 "Eukaryota" "Opisthokonta" "Fungi"  NA    NA    NA     NA    NA     

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV6458",]
#Kingdom     Supergroup     Division Class        Order              Family              Genus           Species                   
#ASV6458 "Eukaryota" "Opisthokonta" "Fungi"  "Ascomycota" "Saccharomycotina" "Saccharomycetales" "Saccharomyces" "Saccharomyces_cerevisiae"

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV774",]
#Kingdom     Supergroup     Division Class        Order            Family            Genus           Species            
#ASV774 "Eukaryota" "Opisthokonta" "Fungi"  "Ascomycota" "Pezizomycotina" "Dothideomycetes" "Aureobasidium" "Aureobasidium_sp."

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV843",]
#       Kingdom     Supergroup     Division Class Order Family Genus Species
#ASV843 "Eukaryota" "Opisthokonta" "Fungi"  NA    NA    NA     NA    NA     

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV847",]
#Kingdom     Supergroup     Division Class           Order                Family              Genus        Species             
#ASV847 "Eukaryota" "Opisthokonta" "Fungi"  "Basidiomycota" "Ustilaginomycotina" "Exobasidiomycetes" "Malassezia" "Malassezia_globosa"

tax_table(FungiPhyloseq_EukAll)[row.names(tax_table(FungiPhyloseq_EukAll)) == "ASV9123",]
#Kingdom     Supergroup     Division Class           Order                Family              Genus        Species         
#ASV9123 "Eukaryota" "Opisthokonta" "Fungi"  "Basidiomycota" "Ustilaginomycotina" "Exobasidiomycetes" "Malassezia" "Malassezia_sp."

#----------------------------------------------------------------------------------------
# Step 3: take a look with dot plot to co-occurrence between primary producers and Fungi
#----------------------------------------------------------------------------------------

# ### ADDED ON 2023 06 12 TO REMOVE THE SAMPLES WITH 0 in Temperature and Salinity ###
# Phyloseq_EukAll <- subset_samples(Phyloseq_EukAll, !(Name_Euk %in% c("S109Euk", "S115Euk", "S36Euk", "PF281Euk", "PF287Euk")))
# ### END OF ADDED ON 2023 06 12 TO REMOVE THE SAMPLES WITH 0 in Temperature and Salinity ###

### ADDED ON 2023 06 12 TO turn the salinities and temperatures equal to 0 into NA ###
sample_data(Phyloseq_EukAll)[sample_data(Phyloseq_EukAll)[["Name_Euk"]] %in% c("S109Euk", "S115Euk", "S36Euk", "PF281Euk", "PF287Euk"),"S"] <- NA
sample_data(Phyloseq_EukAll)[sample_data(Phyloseq_EukAll)[["Name_Euk"]] %in% c("S109Euk", "S115Euk", "S36Euk", "PF281Euk", "PF287Euk"),"T"] <- NA
### END OF ADDED ON 2023 06 12 TO turn the salinities and temperatures equal to 0 into NA ###

# Convert to percent at the OTU level

Phyloseq_EukAll_perc_OTU <- taxa_percentize(Phyloseq_EukAll, TaxLevel = "OTU")

# Extract diatoms from that global phyloseq

DiatomPhyloseq_EukAll <- subset_taxa(Phyloseq_EukAll_perc_OTU, Class == "Bacillariophyta")
DiatomPhyloseq_EukAll <- prune_taxa(taxa_sums(DiatomPhyloseq_EukAll) > 0, DiatomPhyloseq_EukAll)
DiatomPhyloseq_EukAll <- prune_samples(sample_sums(DiatomPhyloseq_EukAll) > 0, DiatomPhyloseq_EukAll)

# Extract Fungi from that global dataset

FungiOTUPercPhyloseq_EukAll <- subset_taxa(Phyloseq_EukAll_perc_OTU, Division == "Fungi")
FungiOTUPercPhyloseq_EukAll <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll) > 0, FungiOTUPercPhyloseq_EukAll)
FungiOTUPercPhyloseq_EukAll <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll) > 0, FungiOTUPercPhyloseq_EukAll)

# Only keep the x most abundant Fungi

topN <- 20
most_abundant_taxa <- sort(taxa_sums(FungiOTUPercPhyloseq_EukAll), TRUE)[1:topN]
FungiOTUPercPhyloseq_EukAll <- prune_taxa(names(most_abundant_taxa), FungiOTUPercPhyloseq_EukAll)

# Extract OTU table and sample_data from the phyloseq object to plot them as geom_point

Fungi_otu_perc_OTU <- otu_table(FungiOTUPercPhyloseq_EukAll)
Fungi_otu_perc_OTU <- data.frame(t(Fungi_otu_perc_OTU))
Fungi_otu_perc_OTU$Name_Euk <- row.names(Fungi_otu_perc_OTU)

Fungi_metadata_perc_OTU <- sample_data(FungiOTUPercPhyloseq_EukAll)

Fungi_otu_metadata_perc_OTU <- list(Fungi_otu_perc_OTU, Fungi_metadata_perc_OTU) %>% reduce(full_join) # Hard check 

dim(Fungi_otu_perc_OTU)
dim(Fungi_metadata_perc_OTU)
dim(Fungi_otu_metadata_perc_OTU)

# Set the dataframe to long format, and add a column that will have three levels (0.2 um, 3 um, metadata) that will be used in facet_grid

Col2Long <- c(colnames(Fungi_otu_perc_OTU)[colnames(Fungi_otu_perc_OTU) != "Name_Euk"], "T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")
Fungi_otu_metadata_perc_OTU_Long <- Fungi_otu_metadata_perc_OTU %>%
  pivot_longer(all_of(Col2Long), names_to = "Variable", values_to = "Count")

Fungi_otu_metadata_perc_OTU_Long$Facet_ID <- NA
Fungi_otu_metadata_perc_OTU_Long[Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% colnames(Fungi_otu_perc_OTU)[colnames(Fungi_otu_perc_OTU) != "Name_Euk"], "Facet_ID"] <- Fungi_otu_metadata_perc_OTU_Long[Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% colnames(Fungi_otu_perc_OTU)[colnames(Fungi_otu_perc_OTU) != "Name_Euk"], "Filter_Euk"]
Fungi_otu_metadata_perc_OTU_Long[which(is.na(Fungi_otu_metadata_perc_OTU_Long[["Facet_ID"]])), "Facet_ID"] <- "Metadata"

Fungi_otu_metadata_perc_OTU_Long$Group_ID <- NA
Fungi_otu_metadata_perc_OTU_Long[Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% colnames(Fungi_otu_perc_OTU)[colnames(Fungi_otu_perc_OTU) != "Name_Euk"], "Group_ID"] <- "Fungi"
Fungi_otu_metadata_perc_OTU_Long[which(is.na(Fungi_otu_metadata_perc_OTU_Long[["Group_ID"]])), "Group_ID"] <- "Metadata"

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
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same graph, but remove S, T, pH and O from the metadata graph

Fungi_otu_metadata_perc_OTU_Long_cut <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "S", "pH", "O")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_cut, f = Fungi_otu_metadata_perc_OTU_Long_cut$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same, but with only chl a as metadata

Fungi_otu_metadata_perc_OTU_Long_chla <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "S", "pH", "O", "NH4", "NO2", "NO3", "PO4", "SIO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_chla, f = Fungi_otu_metadata_perc_OTU_Long_chla$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same, but with only PO4 as metadata

Fungi_otu_metadata_perc_OTU_Long_PO4 <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "S", "pH", "O", "NH4", "NO2", "NO3", "CHLA", "SIO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_PO4, f = Fungi_otu_metadata_perc_OTU_Long_PO4$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same but with only NO3 as metadata

Fungi_otu_metadata_perc_OTU_Long_NO3 <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "S", "pH", "O", "NH4", "NO2", "PO4", "CHLA", "SIO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_NO3, f = Fungi_otu_metadata_perc_OTU_Long_NO3$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same but with only NO2 as metadata

Fungi_otu_metadata_perc_OTU_Long_NO2 <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "S", "pH", "O", "NH4", "NO3", "PO4", "CHLA", "SIO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_NO2, f = Fungi_otu_metadata_perc_OTU_Long_NO2$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

# Plot the same but with only Salinity (S) as metadata

Fungi_otu_metadata_perc_OTU_Long_NO2 <- Fungi_otu_metadata_perc_OTU_Long[!(Fungi_otu_metadata_perc_OTU_Long[["Variable"]] %in% c("T", "NO2", "pH", "O", "NH4", "NO3", "PO4", "CHLA", "SIO4")),]

xs <- split(Fungi_otu_metadata_perc_OTU_Long_NO2, f = Fungi_otu_metadata_perc_OTU_Long_NO2$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
  geom_point() +
  geom_line() +
  scale_colour_manual(values = cc) +
  facet_grid(rows = vars(Facet_ID), scales = "free") +
  theme_bw()

p2 <- p1 %+% xs[["Metadata"]]

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

#------------------------------------------------------------------------------------------------------
# Same graphs, but with only Fungi at different resolution levels (Phylum, family, genus, order, etc.)
#------------------------------------------------------------------------------------------------------

# Extract Fungi from that global dataset

FungiOTUPhyloseq_NormFung <- subset_taxa(Phyloseq_EukAll, Division == "Fungi")
FungiOTUPhyloseq_NormFung <- prune_taxa(taxa_sums(FungiOTUPhyloseq_NormFung) > 0, FungiOTUPhyloseq_NormFung)
FungiOTUPhyloseq_NormFung <- prune_samples(sample_sums(FungiOTUPhyloseq_NormFung) > 0, FungiOTUPhyloseq_NormFung)

#====
# Phylum level (it is actually indicated as Class for Fungi in this dataset)
#====

FungiOTUPhyloseq_NormFung_Clas <- taxa_percentize(FungiOTUPhyloseq_NormFung, TaxLevel = "Class")

# Extract OTU table and sample_data from the phyloseq object to plot them as geom_point

FungiOTUPhyloseq_NormFung_Clas_DF <- psmelt(FungiOTUPhyloseq_NormFung_Clas)

# Set the dataframe to long format, and add a column that will have three levels (0.2 um, 3 um, metadata) that will be used in facet_grid

Col2Long <- c("Abundance", "T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")
FungiOTUPhyloseq_NormFung_Clas_DF_Long <- FungiOTUPhyloseq_NormFung_Clas_DF %>%
  pivot_longer(all_of(Col2Long), names_to = "Variable", values_to = "Count")

FungiOTUPhyloseq_NormFung_Clas_DF_Long$Facet_ID <- NA
FungiOTUPhyloseq_NormFung_Clas_DF_Long[FungiOTUPhyloseq_NormFung_Clas_DF_Long[["Variable"]] == "Abundance", "Facet_ID"] <- FungiOTUPhyloseq_NormFung_Clas_DF_Long[FungiOTUPhyloseq_NormFung_Clas_DF_Long[["Variable"]] == "Abundance", "Filter_Euk"]
FungiOTUPhyloseq_NormFung_Clas_DF_Long[which(is.na(FungiOTUPhyloseq_NormFung_Clas_DF_Long[["Facet_ID"]])), "Facet_ID"] <- "Metadata"

FungiOTUPhyloseq_NormFung_Clas_DF_Long$Group_ID <- NA
FungiOTUPhyloseq_NormFung_Clas_DF_Long[FungiOTUPhyloseq_NormFung_Clas_DF_Long[["Variable"]] == "Abundance", "Group_ID"] <- "Fungi"
FungiOTUPhyloseq_NormFung_Clas_DF_Long[which(is.na(FungiOTUPhyloseq_NormFung_Clas_DF_Long[["Group_ID"]])), "Group_ID"] <- "Metadata"

# Rename "Abundance" in the "Variable" column to be the Class name

FungiOTUPhyloseq_NormFung_Clas_DF_Long[FungiOTUPhyloseq_NormFung_Clas_DF_Long[["Variable"]] == "Abundance", "Variable"] <- FungiOTUPhyloseq_NormFung_Clas_DF_Long[FungiOTUPhyloseq_NormFung_Clas_DF_Long[["Variable"]] == "Abundance", "Class"]


xs <- split(FungiOTUPhyloseq_NormFung_Clas_DF_Long, f = FungiOTUPhyloseq_NormFung_Clas_DF_Long$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
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

# Plot the same but with only NO2 as metadata

FungiOTUPhyloseq_NormFung_Clas_DF_Long_NO2 <- FungiOTUPhyloseq_NormFung_Clas_DF_Long[!(FungiOTUPhyloseq_NormFung_Clas_DF_Long[["Variable"]] %in% c("T", "S", "pH", "O", "NH4", "NO3", "PO4", "CHLA", "SIO4")),]

xs <- split(FungiOTUPhyloseq_NormFung_Clas_DF_Long_NO2, f = FungiOTUPhyloseq_NormFung_Clas_DF_Long_NO2$Group_ID)
p1 <- ggplot(xs[["Fungi"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
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

#------------------------------------------------------------------------------------------------------
# Same graphs, but the abundance is normalized at the domain level
#------------------------------------------------------------------------------------------------------

# # Extract Fungi from that global dataset
# 
# OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
# OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)
# 
# #====
# # Phylum level (it is actually indicated as Class for Fungi in this dataset)
# #====
# 
# OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "OTU")
# 
# # Extract OTU table and sample_data from the phyloseq object to plot them as geom_point
# 
# OTUPhyloseq_NormDomain_DF <- psmelt(OTUPhyloseq_NormDomain)
# 
# # Set the dataframe to long format, and add a column that will have three levels (0.2 um, 3 um, metadata) that will be used in facet_grid
# 
# Col2Long <- c("Abundance", "T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")
# OTUPhyloseq_NormDomain_DF_Long <- OTUPhyloseq_NormDomain_DF %>%
#   pivot_longer(all_of(Col2Long), names_to = "Variable", values_to = "Count")
# 
# OTUPhyloseq_NormDomain_DF_Long$Facet_ID <- NA
# OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Facet_ID"] <- OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Filter_Euk"]
# OTUPhyloseq_NormDomain_DF_Long[which(is.na(OTUPhyloseq_NormDomain_DF_Long[["Facet_ID"]])), "Facet_ID"] <- "Metadata"
# 
# OTUPhyloseq_NormDomain_DF_Long$Group_ID <- NA
# OTUPhyloseq_NormDomain_DF_Long[OTUPhyloseq_NormDomain_DF_Long[["Variable"]] == "Abundance", "Group_ID"] <- "Species"
# OTUPhyloseq_NormDomain_DF_Long[which(is.na(OTUPhyloseq_NormDomain_DF_Long[["Group_ID"]])), "Group_ID"] <- "Metadata"
# 
# # Rename "Abundance" in the "Variable" column to be the Domain name
# 
# OTUPhyloseq_NormDomain_DF_Long_Div <- OTUPhyloseq_NormDomain_DF_Long
# OTUPhyloseq_NormDomain_DF_Long_Div[OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] == "Abundance", "Variable"] <- OTUPhyloseq_NormDomain_DF_Long_Div[OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] == "Abundance", "Division"]
# 
# 
# xs <- split(OTUPhyloseq_NormDomain_DF_Long_Div, f = OTUPhyloseq_NormDomain_DF_Long_Div$Group_ID)
# p1 <- ggplot(xs[["Species"]] ,aes(x = Date_num, y = Count, colour = Variable)) + 
#   geom_point(size = 0.5) +
#   geom_line() +
#   scale_colour_manual(values = cc) +
#   facet_grid(rows = vars(Facet_ID), scales = "free") +
#   theme_bw()
# 
# p2 <- p1 %+% xs[["Metadata"]]
# 
# gA <- ggplotGrob(p1)
# gB <- ggplotGrob(p2)
# grid::grid.newpage()
# grid::grid.draw(rbind(gA, gB))


OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)

# Z-score the sample_data columns (T, Sal, PO4, etc.)

sample_data(OTUPhyloseq_NormDomain)[, c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")] <- lapply(sample_data(OTUPhyloseq_NormDomain)[, c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")], function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))

#====
# Division level 
#====

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


#====
# Fungi OTU (actually phylum) level (20 best)
#====

OTUPhyloseq_NormDomain <- prune_taxa(taxa_sums(Phyloseq_EukAll) > 0, Phyloseq_EukAll)
OTUPhyloseq_NormDomain <- prune_samples(sample_sums(OTUPhyloseq_NormDomain) > 0, OTUPhyloseq_NormDomain)

OTUPhyloseq_NormDomain <- taxa_percentize(OTUPhyloseq_NormDomain, TaxLevel = "OTU")

# Z-score the sample_data columns (T, Sal, PO4, etc.)

sample_data(OTUPhyloseq_NormDomain)[, c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")] <- lapply(sample_data(OTUPhyloseq_NormDomain)[, c("T", "S", "O", "pH", "NH4", "NO3", "NO2", "PO4", "SIO4", "CHLA")], function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))

# Select only Fungi

OTUPhyloseq_NormDomain <- subset_taxa(OTUPhyloseq_NormDomain, Division == "Fungi")

# Select most abundant fungal OTUs

topN <- 20
most_abundant_taxa <- sort(taxa_sums(OTUPhyloseq_NormDomain), TRUE)[1:topN]
OTUPhyloseq_NormDomain <- prune_taxa(names(most_abundant_taxa), OTUPhyloseq_NormDomain)

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
OTUPhyloseq_NormDomain_DF_Long_Div[OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] == "Abundance", "Variable"] <- OTUPhyloseq_NormDomain_DF_Long_Div[OTUPhyloseq_NormDomain_DF_Long_Div[["Variable"]] == "Abundance", "OTU"]

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

#--------------------------
# Investigate Fungi deeper
#--------------------------

FungiOTUPercPhyloseq_EukAll <- subset_taxa(Phyloseq_EukAll_perc_OTU, Division == "Fungi")
FungiOTUPercPhyloseq_EukAll <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll) > 0, FungiOTUPercPhyloseq_EukAll)
FungiOTUPercPhyloseq_EukAll <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll) > 0, FungiOTUPercPhyloseq_EukAll)

# Phylum level

FungiOTUPercPhyloseq_EukAll_percCla <- taxa_percentize(FungiOTUPercPhyloseq_EukAll, TaxLevel = "Class")

FungiOTUPercPhyloseq_EukAll_percCla <- prune_samples(sample_data(FungiOTUPercPhyloseq_EukAll_percCla)[["Filter_Euk"]] == "3µM", FungiOTUPercPhyloseq_EukAll_percCla)

topN_Fung <- 9
most_abundant_taxa_Fung <- tax_table(FungiOTUPercPhyloseq_EukAll_percCla)[names(sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_percCla), TRUE)[1:topN_Fung]), "Class"]
FungiOTUPercPhyloseq_EukAll_percCla <- prune_taxa(tax_table(FungiOTUPercPhyloseq_EukAll_percCla)[, "Class"] %in% most_abundant_taxa_Fung, FungiOTUPercPhyloseq_EukAll_percCla)

plot_bar(FungiOTUPercPhyloseq_EukAll_percCla, x="Name_Euk", fill="Class") + 
  scale_fill_manual(values =iwanthue(n = length(unique(tax_table(FungiOTUPercPhyloseq_EukAll_percCla)[,"Class"])))) # Have an idea of present Kingdoms

# Order level 

FungiOTUPercPhyloseq_EukAll_percCla <- taxa_percentize(FungiOTUPercPhyloseq_EukAll, TaxLevel = "Order")

FungiOTUPercPhyloseq_EukAll_percCla <- prune_samples(sample_data(FungiOTUPercPhyloseq_EukAll_percCla)[["Filter_Euk"]] == "3µM", FungiOTUPercPhyloseq_EukAll_percCla)

topN_Fung <- 9
most_abundant_taxa_Fung <- tax_table(FungiOTUPercPhyloseq_EukAll_percCla)[names(sort(taxa_sums(FungiOTUPercPhyloseq_EukAll_percCla), TRUE)[1:topN_Fung]), "Class"]
FungiOTUPercPhyloseq_EukAll_percCla <- prune_taxa(tax_table(FungiOTUPercPhyloseq_EukAll_percCla)[, "Genus"] %in% most_abundant_taxa_Fung, FungiOTUPercPhyloseq_EukAll_percCla)

plot_bar(FungiOTUPercPhyloseq_EukAll_percCla, x="Name_Euk", fill="Class") + 
  scale_fill_manual(values =iwanthue(n = length(unique(tax_table(FungiOTUPercPhyloseq_EukAll_percCla)[,"Class"])))) # Have an idea of present Kingdoms

#-------------------------
# Investigate dinophyceae
#-------------------------

DinoflagellataPhyloseq_EukAll <- subset_taxa(Phyloseq_EukAll, Division == "Dinoflagellata")
DinoflagellataPhyloseq_EukAll <- prune_taxa(taxa_sums(DinoflagellataPhyloseq_EukAll) > 0, DinoflagellataPhyloseq_EukAll)
DinoflagellataPhyloseq_EukAll <- prune_samples(sample_sums(DinoflagellataPhyloseq_EukAll) > 0, DinoflagellataPhyloseq_EukAll)

DinoflagellataPhyloseq_EukAll_percFam <- taxa_percentize(DinoflagellataPhyloseq_EukAll, TaxLevel = "Family")
#DinoflagellataPhyloseq_EukAll_percFam <- prune_samples(sample_data(DinoflagellataPhyloseq_EukAll_percFam)[["Filter_Euk"]] == "3µM", DinoflagellataPhyloseq_EukAll_percFam)

plot_bar(DinoflagellataPhyloseq_EukAll_percFam, x="Date_num", fill="Family") + 
  facet_grid(~Filter_Euk, scales = "free") +
  scale_fill_manual(values =iwanthue(n = length(unique(tax_table(DinoflagellataPhyloseq_EukAll_percFam)[,"Family"])))) # Have an idea of present Kingdoms

DinoflagellataPhyloseq_EukAll_percFam <- prune_samples(sample_data(DinoflagellataPhyloseq_EukAll_percFam)[["Filter_Euk"]] == "3µM", DinoflagellataPhyloseq_EukAll_percFam)

topN_Dino <- 10
most_abundant_taxa_Dino <- tax_table(DinoflagellataPhyloseq_EukAll_percFam)[names(sort(taxa_sums(DinoflagellataPhyloseq_EukAll_percFam), TRUE)[1:topN_Dino]), "Family"]
DinoflagellataPhyloseq_EukAll_percFam <- prune_taxa(tax_table(DinoflagellataPhyloseq_EukAll_percFam)[, "Family"] %in% most_abundant_taxa_Dino, DinoflagellataPhyloseq_EukAll_percFam)

plot_bar(DinoflagellataPhyloseq_EukAll_percFam, x="Name_Euk", fill="Family") + 
  scale_fill_manual(values =iwanthue(n = length(unique(tax_table(DinoflagellataPhyloseq_EukAll_percFam)[,"Family"])))) # Have an idea of present Kingdoms

# Family 

DinoflagellataPhyloseq_EukAll <- subset_taxa(Phyloseq_EukAll, Family == "Gymnodiniaceae")
DinoflagellataPhyloseq_EukAll <- prune_taxa(taxa_sums(DinoflagellataPhyloseq_EukAll) > 0, DinoflagellataPhyloseq_EukAll)
DinoflagellataPhyloseq_EukAll <- prune_samples(sample_sums(DinoflagellataPhyloseq_EukAll) > 0, DinoflagellataPhyloseq_EukAll)

DinoflagellataPhyloseq_EukAll_percFam <- taxa_percentize(DinoflagellataPhyloseq_EukAll, TaxLevel = "Genus")

DinoflagellataPhyloseq_EukAll_percFam <- prune_samples(sample_data(DinoflagellataPhyloseq_EukAll_percFam)[["Filter_Euk"]] == "3µM", DinoflagellataPhyloseq_EukAll_percFam)

topN_Dino <- 10
most_abundant_taxa_Dino <- tax_table(DinoflagellataPhyloseq_EukAll_percFam)[names(sort(taxa_sums(DinoflagellataPhyloseq_EukAll_percFam), TRUE)[1:topN_Dino]), "Genus"]
DinoflagellataPhyloseq_EukAll_percFam <- prune_taxa(tax_table(DinoflagellataPhyloseq_EukAll_percFam)[, "Genus"] %in% most_abundant_taxa_Dino, DinoflagellataPhyloseq_EukAll_percFam)

plot_bar(DinoflagellataPhyloseq_EukAll_percFam, x="Name_Euk", fill="Genus") + 
  scale_fill_manual(values =iwanthue(n = length(unique(tax_table(DinoflagellataPhyloseq_EukAll_percFam)[,"Genus"])))) # Have an idea of present Kingdoms

#=================================================#
# Correlation heatmap using microbiomeSeq package # /!\ microbiomeSeq issue: KMDA (required package formicrobiomeSeq) is not available for this version of R
#=================================================#

# Convert to percent at the OTU level

Phyloseq_EukAll_perc_OTU <- taxa_percentize(Phyloseq_EukAll, TaxLevel = "OTU")

# Extract Fungi from that global dataset

FungiOTUPercPhyloseq_EukAll <- subset_taxa(Phyloseq_EukAll_perc_OTU, Division == "Fungi")
FungiOTUPercPhyloseq_EukAll <- prune_taxa(taxa_sums(FungiOTUPercPhyloseq_EukAll) > 0, FungiOTUPercPhyloseq_EukAll)
FungiOTUPercPhyloseq_EukAll <- prune_samples(sample_sums(FungiOTUPercPhyloseq_EukAll) > 0, FungiOTUPercPhyloseq_EukAll)

# Make the heatmap (calls a lot of warnings because of stadard deviation being 0)
# it probably means that the ASV is only present in one sample, hence the correlation coefficient
# must be not defined or wrong, and the p-value artificially low enough to be significant

env.taxa.cor <- taxa.env.correlation(FungiOTUPercPhyloseq_EukAll, grouping_column="Filter_Euk", method="pearson", pvalue.threshold=0.05,
                                     padjust.method="BH", adjustment=5, num.taxa=28, select.variables=NULL)

p <- plot_taxa_env(env.taxa.cor)
print(p)

#================================================#
#             Ordination analysis                #
#================================================#

as.Date_origin <- function(x){
  as.Date(x, origin = '2013-05-21')
}

Phyloseq_EukAll_Hell <- transform_sample_counts(Phyloseq_EukAll, function(x) sqrt(x / sum(x)))
GP.ord.Phyloseq.EukAll.Hell <- ordinate(Phyloseq_EukAll_Hell, "NMDS", "bray")
plot_ordination(Phyloseq_EukAll_Hell, GP.ord.Phyloseq.EukAll.Hell, type="sample", color="Filter_Euk")

GP.ord.Phyloseq.EukAll.Hell <- ordinate(Phyloseq_EukAll_Hell, "PCoA", "bray")
plot_ordination(Phyloseq_EukAll_Hell, GP.ord.Phyloseq.EukAll.Hell, type="sample", color = "Date_num", shape="Filter_Euk") +
  scale_colour_gradientn(colours=c('red','green','blue'), labels=as.Date_origin)
