# Sampling graph
#   Input: graph for sampling
#   Output: Smaller graph
source("graphs/feature.graph.R")

#Random Jump (RJ) Sample of ncount nodes from network graph
randomJumpSample<-function(graph,ncount,teleport=0.15) {
  
  node <- sample.int(vcount(graph), 1)
  selected <- rep(NA,ncount)
  selected[[1]]<-node
  i<-2 
  
  while(i<=ncount) {
    neigh<-neighbors(graph,node)
    if(length(neigh)==0 | runif(1)<=teleport){
      node <- sample.int(vcount(graph), 1)
    } else {
      node <- sample(neigh,1)
    }
    if (sum(node==selected,na.rm=TRUE)==0) {
      selected[[i]]<-node
      i<-i+1
      #print(paste0("We now have ",i," nodes."))
    }
  }
  return(induced.subgraph(graph, selected))
}

sample.graph = randomJumpSample(graph, 3000)
sample.graph = decompose.graph(sample.graph)[[1]]

write.graph(sample.graph, file = "data/sample.graph.ml", format = "graphml")


