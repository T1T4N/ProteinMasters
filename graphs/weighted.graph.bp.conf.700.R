# Create graph with all metrics
#   Input: graphs_all file
#   Output: graph with edge attributes all metrics
library(igraph)
library(data.table)

g.sparse = as.data.table(read.table("data/graphs_all.txt", header = TRUE))
graph = graph.data.frame(g.sparse, directed = FALSE)

remove(g.sparse)
