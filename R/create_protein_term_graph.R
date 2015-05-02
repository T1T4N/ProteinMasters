source("create_interaction_graphs.R")

# Нашето множество е BP_conf за HumanPPI700
set = reading_db("data/HumanPPI700_GO_BP_conf.txt", head = TRUE, quote_char = "\"'")
setkey(set, GO, protein)
set[, c("V3", "V5", "V6", "V7", "conf", "count"):=NULL]

protein_term_graph = induced.subgraph(humanPPI700_graph, which(V(humanPPI700_graph)$name %in% set[, protein]))
protein_term_graph = decompose.graph(protein_term_graph)[[1]]
protein_term_graph = set.vertex.attribute(protein_term_graph, "type", value = "protein")
set = set[protein %in% V(protein_term_graph)$name]
protein_term_graph = add.vertices(protein_term_graph, nv = nlevels(set[, GO]), name = levels(set[, GO]), type = "term")


print("Create protein-term graph...")
t1 = Sys.time()
for(i in levels(set[, GO])){
  term = i;
  proteins = as.character(set[term, protein])
  protein_term_graph[term, proteins] = TRUE
}
t2 = Sys.time()
print(paste("Finish in ", t2 - t1))
