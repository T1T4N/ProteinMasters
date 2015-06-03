source("create.interaction.graphs.BP.conf.R")

connected.graph.700 = decompose.graph(humanPPI700.graph)[[1]]
connected.graph.900 = decompose.graph(humanPPI900.graph)[[1]]

clustering.coef = function(x, graph) {
  print(graph)
  mean(transitivity(graph, type = "local", vids = names(x), isolates = "zero"))
}

degree.700 = degree.distribution(connected.graph.700)
degree.900 = degree.distribution(connected.graph.900)
plot(degree.700, log = "xy")
plot(degree.900, log = "xy")

pearson.cor.700 = assortativity.degree(connected.graph.700, FALSE)
pearson.cor.900 = assortativity.degree(connected.graph.900, FALSE)
print(pearson.cor.700)
print(pearson.cor.900)

clustering.spectum.700 = tapply(degree(connected.graph.700), as.factor(degree(connected.graph.700)), clustering.coef, connected.graph.700)
clustering.spectum.900 = tapply(degree(connected.graph.900), as.factor(degree(connected.graph.900)), clustering.coef, connected.graph.900)
plot(names(clustering.spectum.700), clustering.spectum.700 , xlab = "degree", ylab = "clustering coef")
plot(names(clustering.spectum.900), clustering.spectum.900 , xlab = "degree", ylab = "clustering coef")

shortest.paths.700 = shortest.paths(connected.graph.700, weights = NA)
shortest.paths.900 = shortest.paths(connected.graph.900, weights = NA)
max.shortest.path.700 = max(shortest.paths.700)
max.shortest.path.900 = max(shortest.paths.900)
percent.sp.pairs.700 = vector(length = max.shortest.path.700)
percent.sp.pairs.900 = vector(length = max.shortest.path.900)
for(i in 1:max(max.shortest.path.700, max.shortest.path.900)) {
  if(i <= max.shortest.path.700) {
    percent.sp.pairs.700[i] = sum(shortest.paths.700 == i)*100/length(shortest.paths.700)  
  }
  if(i <= max.shortest.path.900) {
    percent.sp.pairs.900[i] = sum(shortest.paths.900 == i)*100/length(shortest.paths.900)
  }
}
plot(percent.sp.pairs.700)
plot(percent.sp.pairs.900)

closeness.700 = -log(closeness(connected.graph.700, weights = NULL))
closeness.900 = -log(closeness(connected.graph.900, weights = NULL))
plot(table(round(closeness.700, 2))*100/length(closeness.700))
plot(table(round(closeness.900, 2))*100/length(closeness.900))

betweenness.700 = -log(betweenness(connected.graph.700, directed = FALSE, weights = NULL, nobigint = FALSE, normalized = TRUE))
betweenness.900 = -log(betweenness(connected.graph.900, directed = FALSE, weights = NULL, nobigint = FALSE, normalized = TRUE))
plot(table(round(betweenness.700))*100/length(betweenness.700))
plot(table(round(betweenness.900))*100/length(betweenness.900))

# Процент на функционални анотации кои се појавуват 
# за прв пат на k растојание
bp.conf.700.graph = decompose.graph(bp.conf.700.graph)[[1]]
bp.conf.900.graph = decompose.graph(bp.conf.900.graph)[[1]]
sp.700 = shortest.paths(bp.conf.700.graph)
sp.900 = shortest.paths(bp.conf.900.graph)
diameter.700 = max(sp.700) #13
diameter.900 = max(sp.900) #14
v.700 = vector(length = diameter.700)
v.900 = vector(length = diameter.900)
for(i in 1:vcount(bp.conf.700.graph)) {
  tmp = V(bp.conf.700.graph)[i]$GO[[1]]
  distance = 1
  while(length(tmp) > 0 && distance <= diameter.700) {
    x = unique(unlist(V(bp.conf.700.graph)[which(sp.700[i,] == distance, useNames = FALSE)]$GO))
    v.700[distance] = v.700[distance] + length(intersect(tmp, x))
    tmp = setdiff(tmp, x)
    distance = distance + 1
  }
}
for(i in 1:vcount(bp.conf.900.graph)) {
  tmp = V(bp.conf.900.graph)[i]$GO[[1]]
  distance = 1
  while(length(tmp) > 0 && distance <= diameter.900) {
    x = unique(unlist(V(bp.conf.900.graph)[which(sp.900[i,] == distance, useNames = FALSE)]$GO))
    v.900[distance] = v.900[distance] + length(intersect(tmp, x))
    tmp = setdiff(tmp, x)
    distance = distance + 1
  }
}
v.700 = v.700/sum(v.700) * 100
v.900 = v.900/sum(v.900) * 100
barplot(v.700, names.arg = c(1:diameter.700))
barplot(v.900, names.arg = c(1:diameter.900))
