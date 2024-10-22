# The goal of this file is to analyse the output files of FlashWeave

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  GENERAL CLEANING AND SETTING  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Cleans the environment
rm(list = ls())
dev.off()
cat("\014")

# Load libraries

require("ggplot2")
require("phyloseq")
library(plyr)
library(igraph)
#library(pairwiseAdonis)
library(RVAideMemoire)
source("/home/emilelaymand/Documents/Science/Master2_MES/StageM2/Scripts/FunctionsScriptR/FunctionsGraphs.R")

# Define convenient functions

altconcat <- function(vector1,vector2) {
  ## This function alternates the elements of vector1 and vector2. Used to fit the syntax of the function E() ##
  mergedvector <- c()
  for (i in 1:length(vector1)) {
    mergedvector <- c(mergedvector, vector1[i],  vector2[i])
  }
  return(mergedvector)
}

findweightsadj <- function(Graph, Node){
  ## This functions takes a graph and a node within this graph, and returns the weights of the edges to the adjacent nodes ##
  # Avoid modifying objects
  GraphObj <- Graph
  MainNode <- Node
  # Find adjacent nodes
  AdjNodes <- names(adjacent_vertices(GraphObj, V(GraphObj))[[MainNode]])
  # Duplicate the name of the main node to fit the number of adjacent nodes
  MainNodeGoodLength <- rep(MainNode, length(AdjNodes))
  # Make a vector that alternate the main node name and the adjacent nodes names to fit the syntax of the function E()
  VectForListWeights <- altconcat(MainNodeGoodLength, AdjNodes)
  # Find the weights of the adjacent nodes
  VectorWeights <- E(GraphObj, P=VectForListWeights)$weight
  names(VectorWeights) <- AdjNodes
  # Return the weights to the adjacent nodes, with their names
  return(VectorWeights)
}

Graph_Free_Living_AllEukaryotesPerc3um_t <- read_graph("FlashWeaveOutput/AllEuk3umPerc.gml", format = "gml")

Graph_Free_Living_AllEukaryotesPerc3um_t_named <- index_as_vertices_names(Graph_Free_Living_AllEukaryotesPerc3um_t)
#degree(Graph_Free_Living_AllEukaryotesPerc3um_t_named)
#AdjMat <- get.adjacency(Graph_Free_Living_AllEukaryotesPerc3um_t_named)

Graph_Free_Living_AllEukaryotesPerc3um_t_named
is_weighted(Graph_Free_Living_AllEukaryotesPerc3um_t_named)
is_directed(Graph_Free_Living_AllEukaryotesPerc3um_t_named)
#E(Graph_Free_Living_AllEukaryotesPerc3um_t_named)$weight

# GraphObj <- Graph_Free_Living_AllEukaryotesPerc3um_t_named
# MainNode <- "ASV224"
# #AdjNodes <- adjacent_vertices(GraphObj, V(GraphObj))$ASV224
# AdjNodes <- names(adjacent_vertices(GraphObj, V(GraphObj))[[MainNode]])
# MainNodeGoodLength <- rep(MainNode, length(AdjNodes))
# VectForListWeights <- altconcat(MainNodeGoodLength, AdjNodes)
# VectorWeights <- E(GraphObj, P=VectForListWeights)$weight
# names(VectorWeights) <- AdjNodes
# #GraphObj_AdjWeight <- as_adj(GraphObj,attr="weight") #, edges = T,sparse = F)

###################################################################
# Import merged phyloseq object (see PrepareFilesForFlashWeave.R) #
###################################################################

Phyloseq_EukAll <- readRDS(file = "PrepareFilesForFlashWeaveOutput/PhyloseqData.rds")

######################################################
# Find the adjacent nodes to the most abundant Fungi #
######################################################

# ASV224

MainNode <- "ASV224"
AdjASVs <- findweightsadj(Graph_Free_Living_AllEukaryotesPerc3um_t_named, MainNode)
AdjASVs
tax_table(Phyloseq_EukAll)[MainNode,]
tax_table(Phyloseq_EukAll)[names(AdjASVs),]

# ASV275

MainNode <- "ASV275"
AdjASVs <- findweightsadj(Graph_Free_Living_AllEukaryotesPerc3um_t_named, MainNode)
AdjASVs
tax_table(Phyloseq_EukAll)[MainNode,]
tax_table(Phyloseq_EukAll)[names(AdjASVs),]

# ASV409 

MainNode <- "ASV409"
AdjASVs <- findweightsadj(Graph_Free_Living_AllEukaryotesPerc3um_t_named, MainNode)
AdjASVs
tax_table(Phyloseq_EukAll)[MainNode,]
tax_table(Phyloseq_EukAll)[names(AdjASVs),]
