library(data.table)
library(igraph)

reading_db = function(data, separator = "", head = FALSE, quote_char = "") {
  initial = read.table(data, sep = separator, header = head, nrow = 10)
  classes = sapply(initial, class)
  tmp = read.table(data, sep = separator, header = head, quote = quote_char, colClasses = classes, comment.char = "")
  tmp = as.data.table(tmp)
  tmp
}

create_interaction_graph = function(x, directed = FALSE) {
  graph = graph.data.frame(x, directed)
  graph = simplify(graph, edge.attr.comb = "first")
  graph
}

print("Reading databases...")
t1 = Sys.time()
humanPPI700 = reading_db("data/HumanPPI700.txt", head = TRUE, quote_char = "\"'")
humanPPI900 = reading_db("data/HumanPPI900.txt", head = TRUE, quote_char = "\"'")
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
humanPPI700_graph = create_interaction_graph(humanPPI700)
humanPPI900_graph = create_interaction_graph(humanPPI900)
t2 = Sys.time()
print(paste("Graphs created in ", t2 - t1))
