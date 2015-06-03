source("create.interaction.graphs.R")
remove(humanPPI900.graph)

# Нашето множество е BP_conf за HumanPPI700
set = as.data.table(read.table("data/HumanPPI700_GO_BP_conf.txt", header = TRUE))
setkey(set, GO, protein)
set[, c("V3", "V5", "V6", "V7", "conf", "count"):=NULL]

# Се извлекуват протеините кои ги има само во нашето множество
protein.term.graph = induced.subgraph(humanPPI700.graph, which(V(humanPPI700.graph)$name %in% set[, protein]))

# Се извлекува само најголемата сврзана компонента
protein.term.graph = decompose.graph(protein.term.graph)[[1]]

# Се означуваат овие јазли како протеини
protein.term.graph = set.vertex.attribute(protein.term.graph, "type", value = "protein")

# Се острануваат протеините кои ги нема во најголемата сврзана компонента
set = set[protein %in% V(protein.term.graph)$name]

# Се поставува нова ознака, анотација
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
