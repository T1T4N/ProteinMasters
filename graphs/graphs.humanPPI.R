# Make graphs from HumanPPI700 set and HumanPPI900 set
#   Input: HumanPPI700, HumanPPI900
#   Output: humanPPI700.graph, humanPPI900.graph (in RAM)
library(data.table)
library(igraph)

create.interaction.graph = function(x, directed = FALSE) {
  graph = graph.data.frame(x, directed)
  graph = simplify(graph, edge.attr.comb = "first")
  graph
}

print("Reading databases...")
t1 = Sys.time()
humanPPI700 = as.data.table(read.table("data/HumanPPI700.txt", header = TRUE))
humanPPI900 = as.data.table(read.table("data/HumanPPI900.txt", header = TRUE))
setkey(humanPPI700, "protein1", "protein2")
setkey(humanPPI900, "protein1", "protein2")
t2 = Sys.time()
print(paste("Reading finish in ", t2 - t1))

print("Create graphs...")
t1 = Sys.time()
humanPPI700[, neighborhood := NULL]
humanPPI700[, fusion := NULL]
humanPPI700[, cooccurence := NULL]
humanPPI700[, coexpression := NULL]
humanPPI700[, experimental := NULL]
humanPPI700[, database := NULL]
humanPPI700[, textmining := NULL]
humanPPI700[, hi := NULL]
humanPPI700[, lit := NULL]
humanPPI700[, venkatesan := NULL]
humanPPI700[, yu := NULL]
humanPPI700[, combined_score := NULL]
humanPPI900[, neighborhood := NULL]
humanPPI900[, fusion := NULL]
humanPPI900[, cooccurence := NULL]
humanPPI900[, coexpression := NULL]
humanPPI900[, experimental := NULL]
humanPPI900[, database := NULL]
humanPPI900[, textmining := NULL]
humanPPI900[, hi := NULL]
humanPPI900[, lit := NULL]
humanPPI900[, venkatesan := NULL]
humanPPI900[, yu := NULL]
humanPPI900[, combined_score := NULL]
humanPPI700.graph = create.interaction.graph(humanPPI700)
humanPPI900.graph = create.interaction.graph(humanPPI900)
t2 = Sys.time()
print(paste("Graphs created in ", t2 - t1))

remove(humanPPI700, humanPPI900, t1, t2)
