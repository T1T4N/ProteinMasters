source("create_interaction_graphs.R")

HumanPPI700_BP_conf = reading_db("data/HumanPPI700_GO_BP_conf.txt", head = TRUE, quote_char = "\"'")
HumanPPI900_BP_conf = reading_db("data/HumanPPI900_GO_BP_conf.txt", head = TRUE, quote_char = "\"'")
setkey(HumanPPI700_BP_conf, protein)
setkey(HumanPPI900_BP_conf, protein)

bp_conf_700_graph = induced.subgraph(humanPPI700_graph, which(V(humanPPI700_graph)$name %in% HumanPPI700_BP_conf[, protein]))
bp_conf_900_graph = induced.subgraph(humanPPI900_graph, which(V(humanPPI900_graph)$name %in% HumanPPI900_BP_conf[, protein]))

t1 = Sys.time()
print("Setting GO terms for protein node...")
for(p in levels(HumanPPI700_BP_conf[, protein])){
  bp_conf_700_graph = set.vertex.attribute(bp_conf_700_graph, name = "GO", index = which(p == V(bp_conf_700_graph)$name), 
                                           value = list(as.vector(HumanPPI700_BP_conf[p, GO])))
}
for(p in levels(HumanPPI900_BP_conf[, protein])){
  bp_conf_900_graph = set.vertex.attribute(bp_conf_900_graph, name = "GO", index = which(p == V(bp_conf_900_graph)$name), 
                                           value = list(as.vector(HumanPPI900_BP_conf[p, GO])))
}
t2 = Sys.time()
print(t2 - t1)

