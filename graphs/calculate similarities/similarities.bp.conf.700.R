# Set content similarity, calculate and set structure and hybrid similarities
# (metrics: Jacard, Resnik, Wang) to bp.conf.700.graph
#   Input: bp.conf.700.graph, content similarity with all metrics in W_sim
#   Output: bp.conf.700.graph edgelist file with all metrics 
source("graphs/graphs.BP.conf.R")
library(doParallel)
registerDoParallel(cores = "3")

graph = decompose.graph(bp.conf.700.graph)[[1]]
edges = get.edges(graph, E(graph, directed = FALSE))

t1 = Sys.time()
W = as.data.table(read.table("data/W_sim", header = TRUE))
setkey(W, "protein1", "protein2")
t2 = Sys.time()
print(t2 - t1)

get.names = function(e) {
  c(V(graph)[edges[e, 2]]$name, V(graph)[edges[e, 1]]$name)
}

all.sim = function(e) {
  p1 = V(graph)[edges[e, 2]]$name
  p2 = V(graph)[edges[e, 1]]$name
  neighbors1 = neighbors(graph, p1)
  neighbors2 = neighbors(graph, p2)
  
  content = W[list(p1, p2), c(jacard, resnik, wang)]
  
  firstSum = c(0, 0, 0)
  for(i in neighbors1) {
    n = V(graph)[i]$name
    if(i > edges[e, 1]) {
      firstSum = firstSum + W[list(p2, n), c(jacard, resnik, wang)]
    }
    if(i < edges[e, 1]) {
      firstSum = firstSum + W[list(n, p2), c(jacard, resnik, wang)]
    }
  }
  firstSum = firstSum/length(neighbors1)
  secondSum = c(0, 0, 0)
  for(i in neighbors2) {
    n = V(graph)[i]$name
    if(i > edges[e, 2]) {
      secondSum = secondSum + W[list(p1, n), c(jacard, resnik, wang)]
    }
    if(i < edges[e, 2]) {
      secondSum = secondSum + W[list(n, p1), c(jacard, resnik, wang)]
    }
  }
  secondSum = secondSum/length(neighbors2)
  structure = (firstSum + secondSum)/2
  
  hybrid = (content + structure)/2
  
  sim = c(content, structure, hybrid)
}

t1 = Sys.time()
result = foreach(e = 1:nrow(edges), .combine = "rbind", .inorder = FALSE) %dopar% {
  s = all.sim(e)
  data.frame(protein1 = V(graph)[edges[e, 2]]$name, protein2 = V(graph)[edges[e, 1]]$name,
             jacard_content = s[1], resnik_content = s[2], wang_content = s[3],
             jacard_structure = s[4], resnik_structure = s[5], wang_structure = s[6],
             jacard_hybrid = s[7], resnik_hybrid = s[8], wang_hybrid = s[9])
}
t2 = Sys.time()
print(t2 - t1)
result = as.data.table(result)

write.table(result, file = "data/graphs_all.txt", col.names = TRUE, row.names = FALSE)

# normalize without observed protein in firstSum and secondSum
# write.table(result, "data/graphs.txt", col.names = TRUE, row.names = FALSE)

remove(edges, W, result, graph, t1, t2)


