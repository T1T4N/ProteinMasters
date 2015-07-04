# Support script for SRWapproach.R
# Write in txt file page ranks of random walk with restart for every node
#   Input: bp.conf.700 with all metrics
#   Output: txt files with vcount elements inside
source("graphs/weighted.graph.bp.conf.700.R")

# wang_hybrid or resnik_hybrid
metric = "wang_hybrid"
epsilon = 1e-12
restart.prob = 0.3
max.iter = 100
A = get.adjacency(graph, attr = metric)

for(v in V(graph)) {
  print(paste("Calculate page ranks of node", v))
  Q = (1-restart.prob) * (Matrix::Diagonal(x = (Matrix::rowSums(A))^(-1)) %*% A)
  Q[, v] = Q[, v] + restart.prob  
  p = rep(1/vcount(graph), vcount(graph))
  t = 0
  repeat{
    t = t + 1
    p.new = p %*% Q
    p.new = p.new/sum(p.new)
    pre = p
    p = p.new
    if((max(abs(pre - p)) < epsilon) || t > max.iter){
      break
    }
  }
  write(as.vector(p), file = paste0("data/wh_affinity/",v,".txt"))
}
