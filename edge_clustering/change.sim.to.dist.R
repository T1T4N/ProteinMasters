# This script change transformed graphs' weights to be distance not similarity
# distance = 1 - similarity
#   Input: 9 transformed graphs
#   Output: 9 graphs in ncol format
library(data.table)
library(igraph)

jc = as.data.table(read.table("data/transformGraphs/jaccard_content.txt", header = TRUE))
jc_graph = graph.data.frame(jc, directed = TRUE)
jc_graph = set.edge.attribute(graph = jc_graph, name = "weight", 
                              value = 1- get.edge.attribute(graph = jc_graph, name = "weight"))
write.graph(jc_graph, "data/transformGraphs/jc.ncol", format = "ncol")

js = as.data.table(read.table("data/transformGraphs/jaccard_structure.txt", header = TRUE))
js_graph = graph.data.frame(js, directed = TRUE)
js_graph = set.edge.attribute(graph = js_graph, name = "weight", 
                              value = 1- get.edge.attribute(graph = js_graph, name = "weight"))
write.graph(js_graph, "data/transformGraphs/js.ncol", format = "ncol")

jh = as.data.table(read.table("data/transformGraphs/jaccard_hybrid.txt", header = TRUE))
jh_graph = graph.data.frame(jh, directed = TRUE)
jh_graph = set.edge.attribute(graph = jh_graph, name = "weight", 
                              value = 1- get.edge.attribute(graph = jh_graph, name = "weight"))
write.graph(jh_graph, "data/transformGraphs/jh.ncol", format = "ncol")

rc = as.data.table(read.table("data/transformGraphs/resnik_content.txt", header = TRUE))
rc_graph = graph.data.frame(rc, directed = TRUE)
rc_graph = set.edge.attribute(graph = rc_graph, name = "weight", 
                              value = 1- get.edge.attribute(graph = rc_graph, name = "weight"))
write.graph(rc_graph, "data/transformGraphs/rc.ncol", format = "ncol")

rs = as.data.table(read.table("data/transformGraphs/resnik_structure.txt", header = TRUE))
rs_graph = graph.data.frame(rs, directed = TRUE)
rs_graph = set.edge.attribute(graph = rs_graph, name = "weight", 
                              value = 1- get.edge.attribute(graph = rs_graph, name = "weight"))
write.graph(rs_graph, "data/transformGraphs/rs.ncol", format = "ncol")

rh = as.data.table(read.table("data/transformGraphs/resnik_hybrid.txt", header = TRUE))
rh_graph = graph.data.frame(rh, directed = TRUE)
rh_graph = set.edge.attribute(graph = rh_graph, name = "weight", 
                              value = 1- get.edge.attribute(graph = rh_graph, name = "weight"))
write.graph(rh_graph, "data/transformGraphs/rh.ncol", format = "ncol")

wc = as.data.table(read.table("data/transformGraphs/wang_content.txt", header = TRUE))
wc_graph = graph.data.frame(wc, directed = TRUE)
wc_graph = set.edge.attribute(graph = wc_graph, name = "weight", 
                              value = 1- get.edge.attribute(graph = wc_graph, name = "weight"))
write.graph(wc_graph, "data/transformGraphs/wc.ncol", format = "ncol")

ws = as.data.table(read.table("data/transformGraphs/wang_structure.txt", header = TRUE))
ws_graph = graph.data.frame(ws, directed = TRUE)
ws_graph = set.edge.attribute(graph = ws_graph, name = "weight", 
                              value = 1- get.edge.attribute(graph = ws_graph, name = "weight"))
write.graph(ws_graph, "data/transformGraphs/ws.ncol", format = "ncol")

wh = as.data.table(read.table("data/transformGraphs/wang_hybrid.txt", header = TRUE))
wh_graph = graph.data.frame(wh, directed = TRUE)
wh_graph = set.edge.attribute(graph = wh_graph, name = "weight", 
                              value = 1- get.edge.attribute(graph = wh_graph, name = "weight"))
write.graph(wh_graph, "data/transformGraphs/wh.ncol", format = "ncol")