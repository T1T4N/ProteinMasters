library(igraph)
library(data.table)

dt = as.data.table(read.table("data//graphs_all.txt", header = TRUE))
humanPPI700 = as.data.table(read.table("data/HumanPPI700.txt", header = TRUE))
setkey(humanPPI700, "protein1", "protein2")

dt[, neighborhood := humanPPI700[dt[, list(protein1, protein2)], neighborhood/1000]]
dt[, fusion := humanPPI700[dt[, list(protein1, protein2)], fusion/1000]]
dt[, cooccurence := humanPPI700[dt[, list(protein1, protein2)], cooccurence/1000]]
dt[, coexpression := humanPPI700[dt[, list(protein1, protein2)], coexpression/1000]]
dt[, experimental := humanPPI700[dt[, list(protein1, protein2)], experimental/1000]]
dt[, database := humanPPI700[dt[, list(protein1, protein2)], database/1000]]
dt[, textmining := humanPPI700[dt[, list(protein1, protein2)], textmining/1000]]
graph = graph.data.frame(dt, directed = FALSE)

remove(dt, humanPPI700)