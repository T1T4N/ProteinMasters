# This script will transform graph to graph with edges as vertices and
# edges will be calculated by Eq.5.38 inside paper
# Because of large number of vertices this script is copied on few instances
# and with variable current we split the calculation on R sessions
# (we not work with multicore processing because we need observation on our
# execution)
#   Input: original graph with all 9 metrics as edge attributes
#   Output: 9 transformed graphs with weights: metric calculations
source("graphs/weighted.graph.bp.conf.700.R")

jc = data.frame(edge1 = numeric(), edge2 = numeric(), weight = numeric(), stringsAsFactors = FALSE)
js = data.frame(edge1 = numeric(), edge2 = numeric(), weight = numeric(), stringsAsFactors = FALSE)
jh = data.frame(edge1 = numeric(), edge2 = numeric(), weight = numeric(), stringsAsFactors = FALSE)
rc = data.frame(edge1 = numeric(), edge2 = numeric(), weight = numeric(), stringsAsFactors = FALSE)
rs = data.frame(edge1 = numeric(), edge2 = numeric(), weight = numeric(), stringsAsFactors = FALSE)
rh = data.frame(edge1 = numeric(), edge2 = numeric(), weight = numeric(), stringsAsFactors = FALSE)
wc = data.frame(edge1 = numeric(), edge2 = numeric(), weight = numeric(), stringsAsFactors = FALSE)
ws = data.frame(edge1 = numeric(), edge2 = numeric(), weight = numeric(), stringsAsFactors = FALSE)
wh = data.frame(edge1 = numeric(), edge2 = numeric(), weight = numeric(), stringsAsFactors = FALSE)

current = 100000:ecount(graph)
t1 = Sys.time()
for(e in current) {
  t3 = Sys.time()
  edge = get.edge(graph, e)
  v1 = edge[1]
  v2 = edge[2]
  neig1 = get.edge.ids(graph, as.vector(t(get.edges(graph, E(graph)[from(v1)]))), directed = F)
  neig1 = neig1[-which(neig1 == e)]
  neig2 = get.edge.ids(graph, as.vector(t(get.edges(graph, E(graph)[from(v2)]))), directed = F)
  neig2 = neig2[-which(neig2 == e)]
  
  we.jc = get.edge.attribute(graph, "jacard_content", e)
  we.js = get.edge.attribute(graph, "jacard_structure", e)
  we.jh = get.edge.attribute(graph, "jacard_hybrid", e)
  we.rc = get.edge.attribute(graph, "resnik_content", e)
  we.rs = get.edge.attribute(graph, "resnik_structure", e)
  we.rh = get.edge.attribute(graph, "resnik_hybrid", e)
  we.wc = get.edge.attribute(graph, "wang_content", e)
  we.ws = get.edge.attribute(graph, "wang_structure", e)
  we.wh = get.edge.attribute(graph, "wang_hybrid", e)
  
  if(length(neig1) != 0) {
    w = get.edge.attribute(graph, "jacard_content", neig1) 
    c1.jc = we.jc/(sum(w)+we.jc-w)
    w = get.edge.attribute(graph, "jacard_structure", neig1)
    c1.js = we.js/(sum(w)+we.js-w)
    w = get.edge.attribute(graph, "jacard_hybrid", neig1)
    c1.jh = we.jh/(sum(w)+we.jh-w)
    w = get.edge.attribute(graph, "resnik_content", neig1)
    c1.rc = we.rc/(sum(w)+we.rc-w)
    w = get.edge.attribute(graph, "resnik_structure", neig1)
    c1.rs = we.rs/(sum(w)+we.rs-w)
    w = get.edge.attribute(graph, "resnik_hybrid", neig1)
    c1.rh = we.rh/(sum(w)+we.rh-w)
    w = get.edge.attribute(graph, "wang_content", neig1)
    c1.wc = we.wc/(sum(w)+we.wc-w)
    w = get.edge.attribute(graph, "wang_structure", neig1)
    c1.ws = we.ws/(sum(w)+we.ws-w)
    w = get.edge.attribute(graph, "wang_hybrid", neig1)
    c1.wh = we.wh/(sum(w)+we.wh-w)
    
    jc = rbind(jc, data.table(edge1 = e, edge2 = neig1, weight = c1.jc))
    js = rbind(js, data.table(edge1 = e, edge2 = neig1, weight = c1.js))
    jh = rbind(jh, data.table(edge1 = e, edge2 = neig1, weight = c1.jh))
    rc = rbind(rc, data.table(edge1 = e, edge2 = neig1, weight = c1.rc))
    rs = rbind(rs, data.table(edge1 = e, edge2 = neig1, weight = c1.rs))
    rh = rbind(rh, data.table(edge1 = e, edge2 = neig1, weight = c1.rh))
    wc = rbind(wc, data.table(edge1 = e, edge2 = neig1, weight = c1.wc))
    ws = rbind(ws, data.table(edge1 = e, edge2 = neig1, weight = c1.ws))
    wh = rbind(wh, data.table(edge1 = e, edge2 = neig1, weight = c1.wh))
  }
  
  if(length(neig2) != 0) {
    w = get.edge.attribute(graph, "jacard_content", neig2)
    c2.jc = we.jc/(sum(w)+we.jc-w)
    w = get.edge.attribute(graph, "jacard_structure", neig2)
    c2.js = we.js/(sum(w)+we.js-w)
    w = get.edge.attribute(graph, "jacard_hybrid", neig2)
    c2.jh = we.jh/(sum(w)+we.jh-w)
    w = get.edge.attribute(graph, "resnik_content", neig2)
    c2.rc = we.rc/(sum(w)+we.rc-w)
    w = get.edge.attribute(graph, "resnik_structure", neig2)
    c2.rs = we.rs/(sum(w)+we.rs-w)
    w = get.edge.attribute(graph, "resnik_hybrid", neig2)
    c2.rh = we.rh/(sum(w)+we.rh-w)
    w = get.edge.attribute(graph, "wang_content", neig2)
    c2.wc = we.wc/(sum(w)+we.wc-w)
    w = get.edge.attribute(graph, "wang_structure", neig2)
    c2.ws = we.ws/(sum(w)+we.ws-w)
    w = get.edge.attribute(graph, "wang_hybrid", neig2)
    c2.wh = we.wh/(sum(w)+we.wh-w)
    
    jc = rbind(jc, data.table(edge1 = e, edge2 = neig2, weight = c2.jc))
    js = rbind(js, data.table(edge1 = e, edge2 = neig2, weight = c2.js))
    jh = rbind(jh, data.table(edge1 = e, edge2 = neig2, weight = c2.jh))
    rc = rbind(rc, data.table(edge1 = e, edge2 = neig2, weight = c2.rc))
    rs = rbind(rs, data.table(edge1 = e, edge2 = neig2, weight = c2.rs))
    rh = rbind(rh, data.table(edge1 = e, edge2 = neig2, weight = c2.rh))
    wc = rbind(wc, data.table(edge1 = e, edge2 = neig2, weight = c2.wc))
    ws = rbind(ws, data.table(edge1 = e, edge2 = neig2, weight = c2.ws))
    wh = rbind(wh, data.table(edge1 = e, edge2 = neig2, weight = c2.wh))
  }
  t4 = Sys.time()
  print(paste(e, t4-t3, sep = ":")) 
}
t2 = Sys.time()
print(t2 - t1)

write.table(jc, paste0(paste("data/jaccard_content", current[1], current[length(current)], sep = "_"), ".txt"), row.names = FALSE, col.names = TRUE)
write.table(js, paste0(paste("data/jaccard_structure", current[1], current[length(current)], sep = "_"), ".txt"), row.names = FALSE, col.names = TRUE)
write.table(jh, paste0(paste("data/jaccard_hybrid", current[1], current[length(current)], sep = "_"), ".txt"), row.names = FALSE, col.names = TRUE)
write.table(rc, paste0(paste("data/resnik_content", current[1], current[length(current)], sep = "_"), "txt"), row.names = FALSE, col.names = TRUE)
write.table(rs, paste0(paste("data/resnik_structure", current[1], current[length(current)], sep = "_"), ".txt"), row.names = FALSE, col.names = TRUE)
write.table(rh, paste0(paste("data/resnik_hybrid", current[1], current[length(current)], sep = "_"), ".txt"), row.names = FALSE, col.names = TRUE)
write.table(wc, paste0(paste("data/wang_content", current[1], current[length(current)], sep = "_"), ".txt"), row.names = FALSE, col.names = TRUE)
write.table(ws, paste0(paste("data/wang_structure", current[1], current[length(current)], sep = "_"), ".txt"), row.names = FALSE, col.names = TRUE)
write.table(wh, paste0(paste("data/wang_hybrid", current[1], current[length(current)], sep = "_"), ".txt"), row.names = FALSE, col.names = TRUE)

write.table(as.data.table(get.edgelist(graph, names = FALSE)), file = "data/transformGraphs/edges.txt", 
            row.names = TRUE, col.names = FALSE)

