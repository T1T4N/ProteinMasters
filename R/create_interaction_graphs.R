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
  weight_vector = get.edge.attribute(graph, "combined_score")
  graph = set.edge.attribute(graph, "weight", value = weight_vector)
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
humanPPI700_graph = create_interaction_graph(humanPPI700)
humanPPI900_graph = create_interaction_graph(humanPPI900)
t2 = Sys.time()
print(paste("Graphs created in ", t2 - t1))
