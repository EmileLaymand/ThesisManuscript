##################################################################################
# Script for data analysis of the SOLA metagenomics dataset                      #
# 1st section: Eukaryotes 0.2 micrometers (ADNr 18S ? What region)               #
# 2nd section: Eukaryotes 3 micrometers (ADNr 18S ? What region)                 #
# 14/02/2021                                                                     #
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
setwd("/home/emilelaymand/Documents/Science/Master2_MES/StageM2")

require("ggplot2")
require("phyloseq")
library(microViz)
library(stringr)
library(tidyverse)
library(gridExtra)

source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/Generics.R")
source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/Pathnames.R")
source("/home/emilelaymand/Documents/Science/Cadagno/S4PhyloseqDecontam/ClassFiles/phyloseqExtended.R")

#==============#
#     Flags    #
#==============#

#================================================#
# 1st section: Eukaryotes 0.2 micrometers import #
#================================================#

MetaG_DF <- read.csv("/home/emilelaymand/Documents/Science/These/SOLA_Galand/SOLA/wetransfer_sola_2021-10-18_1242/BBMOSOLA-common-dates-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.gene_name-kegg_annotation.tbl.hierarchy.txt", sep = "\t")
MetaG_DF <- separate(data = MetaG_DF, col = Gene.name..KEGG.annotation, into = c("Gene.name", "KEGG.annotation"), sep = ";")
