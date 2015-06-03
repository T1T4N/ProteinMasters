source("create.interaction.graphs.R")

humanPPI700.BP.conf = as.data.table(read.table("data/HumanPPI700_GO_BP_conf.txt", header = TRUE))
humanPPI900.BP.conf = as.data.table(read.table("data/HumanPPI900_GO_BP_conf.txt", header = TRUE))
setkey(humanPPI700.BP.conf, protein)
setkey(humanPPI900.BP.conf, protein)

bp.conf.700.graph = induced.subgraph(humanPPI700.graph, which(V(humanPPI700.graph)$name %in% humanPPI700.BP.conf[, protein]))
bp.conf.900.graph = induced.subgraph(humanPPI900.graph, which(V(humanPPI900.graph)$name %in% humanPPI900.BP.conf[, protein]))

t1 = Sys.time()
print("Setting GO terms for protein node...")
for(p in levels(humanPPI700.BP.conf[, protein])){
  bp.conf.700.graph = set.vertex.attribute(bp.conf.700.graph, name = "GO", index = which(p == V(bp.conf.700.graph)$name), 
                                           value = list(as.vector(humanPPI700.BP.conf[p, GO])))
}
for(p in levels(humanPPI900.BP.conf[, protein])){
  bp.conf.900.graph = set.vertex.attribute(bp.conf.900.graph, name = "GO", index = which(p == V(bp.conf.900.graph)$name), 
                                           value = list(as.vector(humanPPI900.BP.conf[p, GO])))
}
t2 = Sys.time()
print(t2 - t1)

remove(humanPPI700.BP.conf, humanPPI900.BP.conf, t1, t2, p)

