# Graph with two types of nodes: protein and term
# new edges: protein has-a term
#   Input: humanPPI700.graph, HumanPPI700_GO_BP_conf
#   Output: protein-term graph (RAM)
source("graphs/graphs.humanPPI.R")
remove(humanPPI900.graph)

set = as.data.table(read.table("data/HumanPPI700_GO_BP_conf.txt", header = TRUE))
setkey(set, GO, protein)
set[, c("V3", "V5", "V6", "V7", "conf", "count"):=NULL]

protein.term.graph = induced.subgraph(humanPPI700.graph, which(V(humanPPI700.graph)$name %in% set[, protein]))
protein.term.graph = decompose.graph(protein.term.graph)[[1]]
protein.term.graph = set.vertex.attribute(protein.term.graph, "type", value = "protein")
set = set[protein %in% V(protein.term.graph)$name]
protein.term.graph = add.vertices(protein.term.graph, nv = nlevels(set[, GO]), name = levels(set[, GO]), type = "term")

print("Create protein-term graph...")
t1 = Sys.time()
for(i in levels(set[, GO])){
  term = i;
  proteins = as.character(set[term, protein])
  protein.term.graph[term, proteins] = TRUE
}
t2 = Sys.time()
print(paste("Finish in ", t2 - t1))

remove(set, term, i, t1, t2, proteins)
