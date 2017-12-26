#====================================================================
# To read graphData and visualize the correspondig graph
# By (Ameer) Amirhossein Afshinfard - a.afshinfard@gmail.com
# Bioinformatics Research Lab, Sharif Unversity of Technology.
#====================================================================
## read data
library(igraph)


dev.off()

sparseMatrix = read.table("/home/ameer/ExactSV/SV_out/graphData",sep = "\t",col.names=c("row","col","Weight"));
sparseMatrix
nodeLoci = read.table("/home/ameer/ExactSV/SV_out/graphData2",sep = "\t",col.names=c("X","Y"));
nodeLoci
nodeWeights = read.table("/home/ameer/ExactSV/SV_out/graphData3",sep = "\t",col.names=c("Weight"));
nodeWeights
graphDim = dim(nodeWeights)[1]

adjaMatrix = matrix(0,graphDim,graphDim)
adjaMatrix[as.matrix(sparseMatrix[,1:2])]=sparseMatrix[,3]
sum(as.logical(adjaMatrix))

gr <- graph_from_adjacency_matrix(adjaMatrix, mode = c("", "undirected",
                                                "max", "min", "upper", "lower", "plus"), weighted = TRUE, diag = TRUE,
                            add.colnames = NULL, add.rownames = NA)
gr <- graph_from_adjacency_matrix(adjaMatrix,weighted = TRUE,mode = "undirected",diag = TRUE)
E(gr)
plot(gr)

aM = as.data.frame(adjaMatrix)
plot(aM$index)
from_adjacency

#####################################################################################################
#####################################################################################################
#####################################################################################################

install.packages("GGally")
install.packages("network")
install.packages("sna")
install.packages("devtools")
devtools::install_github("briatte/ggnet")
library(ggnet)
library(network)
library(sna)
library(ggplot2)

net = network(adjaMatrix,directed = FALSE,ignore.eval = FALSE,names.eval = "weights")
network.vertex.names(net) = as.character(1:10000)
net %v% "x" = nodeLoci[,1]
net %v% "y" = nodeLoci[,2]
ggnet2(net, node.size = 5, mode = c("x", "y"),label = as.character(as.matrix(nodeWeights)), label.color = "grey",edge.label = "weights",edge.label.color = "darkred", edge.size = 1, edge.color = "grey",node.color = "black")


dev.off()
