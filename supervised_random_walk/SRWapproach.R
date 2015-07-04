# Script will calculate affinity scores of nodes with supervised random walk
# This script should be copied and executed on many R session for faster computation
# Computation is very slow
#   Input: graph with edges with attribute vector
#   Output: vcount txt files with affinity score from SRW

# source("graphs/feature.graph.R")
library(igraph)
library(data.table)
source("supervised_random_walk/supervised.random.walk.R")

sample.graph = read.graph("data/sample.graph.ml", format = "graphml")

random.walk.restart = function(graph, weight.attr, restart.prob = 0.3, max.iter = 100) {
  epsilon = 1e-12
  A = get.adjacency(graph, attr = weight.attr)
  affinity.matrix = matrix(nrow = vcount(graph), ncol = vcount(graph))
  
  for(v in V(graph)) {
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
    affinity.matrix[, v] = as.vector(p)
  }
  affinity.matrix
}

calculate.D.nodes = function(affinity.matrix, graph, source.node, type = c("union", "intersect", "none"), part = 884) {
  whichpart <- function(x, n=30) {
    nx <- length(x)
    p <- nx-n
    xp <- sort(x, partial=p)[p]
    which(x > xp)
  }
  
  D = neighbors(graph, source.node)
  if(type != "none") {
    candidates = list()
    k = 1
    for(v in D) {
      candidates[[k]] = whichpart(affinity.matrix[, v], n = part)
      k = k + 1
    }
    if(type == "union") {
      candidates = Reduce(union, candidates)
      D = union(D, candidates)
      D = setdiff(D, source.node)
    } 
    if(type == "intersect") {
      candidates = Reduce(intersect, candidates)
      D = union(D, candidates)
      D = setdiff(D, source.node)
    }
  }
  D
}

t1 = Sys.time()
affinity.matrix = random.walk.restart(sample.graph, "wang_hybrid")
t2 = Sys.time()
print(paste("Calculating RWR affinity score in", t2 - t1))

current = 50:1000
t1 = Sys.time()
for(v in current){
  print(v)
  D = calculate.D.nodes(affinity.matrix, sample.graph, v, "intersect", part = round(vcount(sample.graph)*0.12))
  write(supervised.rw(sample.graph, v, D, quiet = T), file = paste0("data/SRW_wh_intersect/sample/", v, ".txt"))
}
t2 = Sys.time()
print(t2 - t1)




