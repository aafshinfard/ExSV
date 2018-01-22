#====================================================================
# To read graphData and visualize the correspondig graph
# By (Ameer) Amirhossein Afshinfard - a.afshinfard@gmail.com
# Bioinformatics Research Lab, Sharif Unversity of Technology.
# https://rstudio-pubs-static.s3.amazonaws.com/74248_3bd99f966ed94a91b36d39d8f21e3dc3.html
# http://users.cecs.anu.edu.au/~xlx/teaching/css2013/handout/igraph.pdf
# http://igraph.org/r/doc/tkplot.html
# https://stackoverflow.com/questions/5364264/how-to-control-the-igraph-plot-layout-with-fixed-positions
# https://github.com/igraph/rigraph/issues/63
# http://igraph.org/r/doc/subgraph.html
# http://www.shizukalab.com/toolkits/sna/plotting-networks-pt-2
# http://igraph.org/r/doc/decompose.html
# http://igraph.org/r/doc/components.html
# http://igraph.org/r/doc/layout_in_circle.html
# http://igraph.org/r/doc/edge_attr.html
# 
#====================================================================
## read data
library(igraph)
#<<<<<<< HEAD
#=======


#>>>>>>> 5e7bf2f71818386c65a6e73cd252a10aca92b1d5
dev.off()

sparseMatrix = read.table("/home/ameer/ExactSV/SV_out/graphData",sep = "\t",col.names=c("row","col","Weight"));
#sparseMatrix[,1:2]=sparseMatrix[,1:2]+1
nodeLoci = read.table("/home/ameer/ExactSV/SV_out/graphData2",sep = "\t",col.names=c("X","Y"));
nodeLoci
nodeWeights = read.table("/home/ameer/ExactSV/SV_out/graphData3",sep = "\t",col.names=c("Weight"));
nodeWeights
graphDim = dim(nodeWeights)[1]

adjaMatrix = matrix(0,graphDim,graphDim)
adjaMatrix[as.matrix(sparseMatrix[,1:2])]=sparseMatrix[,3]
sum(as.logical(adjaMatrix))
class(adjaMatrix)
dim(adjaMatrix)

myGraph = graph_from_adjacency_matrix(adjaMatrix, mode = c("max"), weighted = TRUE, diag = TRUE,
                            add.colnames = NULL, add.rownames = NA)
E(myGraph)
is_weighted(myGraph)
tkplot(myGraph)
plot(myGraph)
l <-layout.reingold.tilford(myGraph)
class(l )
class(nodeLoci)
l = as.matrix(nodeLoci)
l2 = l
l2[,1] = log(l[,1])
l2[32,1] = 14
plot(myGraph,
     layout = l)
tkplot(myGraph,
     layout = l2)
tkplot(myGraph,
       layout = l)
a= V(myGraph)
a[1] = 2
a = 1:56
a[1]=2
a[2]=1
a
l = layout_in_circle(myGraph, order = V(myGraph))
l = layout_in_circle(myGraph, order = a)
tkplot(myGraph,
       layout = l)
node
l
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

ggnet2(net, node.size = 5,label = as.character(as.matrix(nodeWeights)),
       label.color = "black",edge.label = "weights",edge.label.color = "darkred", edge.size = 1,
       edge.color = "black",node.color = "grey")


net %v% "x" = nodeLoci[,1]
net %v% "y" = as.numeric(as.matrix((runif(length(nodeLoci[,2]))/5)-0.1))

ggnet2(net, node.size = 5, mode = c("x", "y"),label = as.character(as.matrix(nodeWeights)),
       label.color = "black",edge.label = "weights",edge.label.color = "darkred", edge.size = 1,
       edge.color = "black",node.color = "grey")


net %v% "x" = nodeLoci[,1]
net %v% "y" = as.numeric(nodeLoci[,2]+as.matrix((runif(length(nodeLoci[,2]))/5)-0.1))

ggnet2(net, node.size = 5, mode = c("x", "y"),label = as.character(as.matrix(nodeWeights)),
       label.color = "black",edge.label = "weights",edge.label.color = "darkred", edge.size = 1,
       edge.color = "black",node.color = "grey")


dev.off()
