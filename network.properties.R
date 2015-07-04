# Network properties for bp.700.conf and bp.900.conf graph
#   Input: bp.700.conf and bp.900.conf graph
#   Output: properties 
source("graphs/graphs.BP.conf.R")

connected.graph.700 = decompose.graph(humanPPI700.graph)[[1]]
connected.graph.900 = decompose.graph(humanPPI900.graph)[[1]]

clustering.coef = function(x, graph) {
  mean(transitivity(graph, type = "local", vids = names(x), isolates = "zero"))
}

degree.700 = degree.distribution(connected.graph.700)*100
degree.900 = degree.distribution(connected.graph.900)*100
plot(degree.700, log = "xy", xlab = "degree", ylab = "percent of proteins")
plot(degree.900, log = "xy", xlab = "degree", ylab = "percent of proteins")

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
plot(percent.sp.pairs.700, xlab = "shortest paths length", ylab = "percent of proteins pairs")
plot(percent.sp.pairs.900, xlab = "shortest paths length", ylab = "percent of proteins pairs")

closeness.700 = -log(closeness(connected.graph.700, weights = NULL))
closeness.900 = -log(closeness(connected.graph.900, weights = NULL))
plot(table(round(closeness.700, 2))*100/length(closeness.700), xlab = "closeness", ylab = "percent of proteins")
plot(table(round(closeness.900, 2))*100/length(closeness.900), xlab = "closeness", ylab = "percent of proteins")

betweenness.700 = -log(betweenness(connected.graph.700, directed = FALSE, weights = NULL, nobigint = FALSE, normalized = TRUE))
betweenness.900 = -log(betweenness(connected.graph.900, directed = FALSE, weights = NULL, nobigint = FALSE, normalized = TRUE))
plot(table(round(betweenness.700))*100/length(betweenness.700), xlab = "degree", ylab = "betweeness")
plot(table(round(betweenness.900))*100/length(betweenness.900), xlab = "degree", ylab = "betweeness")

# Percent of functional annotation which appear for first time in k distance
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
v.700 = v.700/sum(v.700)
v.900 = v.900/sum(v.900)
barplot(v.700, names.arg = c(1:diameter.700))
barplot(v.900, names.arg = c(1:diameter.900))

library(doParallel)
registerDoParallel(cores = 2)
# Percent of proteins who share one common function
# Dont do calculation again, 5-10 days for calculation
# distance is changed and that part of code is runned again
#   sp.700 = shortest.paths(bp.conf.700.graph)
#   sp.700[lower.tri(sp.700)] = 0
#   diameter.700 = max(sp.700) #13
#   distance = 1
#   pairs = which(sp.700 == distance, arr.ind = TRUE, useNames = FALSE)
#   current = 3500001:nrow(pairs)
#   sum = foreach(i = current, .combine = "+", .inorder = FALSE) %dopar% {
#     length(intersect(V(bp.conf.700.graph)[pairs[i, 1]]$GO[[1]], V(bp.conf.700.graph)[pairs[i, 2]]$GO[[1]])) > 0
#   }
#
#   sp.900 = shortest.paths(bp.conf.900.graph)
#   sp.900[lower.tri(sp.900)] = 0
#   diameter.900 = max(sp.900) #14
#   distance = 1
#   pairs = which(sp.900 == distance, arr.ind = TRUE, useNames = FALSE)
#   current = 1:nrow(pairs)
#   sum = foreach(i = current, .combine = "+", .inorder = FALSE) %dopar% {
#     length(intersect(V(bp.conf.900.graph)[pairs[i, 1]]$GO[[1]], V(bp.conf.900.graph)[pairs[i, 2]]$GO[[1]])) > 0
#   }
freq.700 = c(83569, 701200, 1890523, 1168566, 301229, 48944, 6285, 754, 149, 37, 12, 0, 0)
all.700 = c(135621, 2818306, 16213538, 15002615, 4140212, 645261, 84357, 26880, 23111, 4561, 399, 38, 4)
percents.700 = freq.700/all.700
barplot(percents.700, names.arg = 1:13)

freq.900 = c(46579, 346649, 1080329, 1085728, 452676, 123705, 29123, 6856, 2206, 524, 55, 9, 0, 0)
all.900 = c(63372, 977734, 6409388, 10889748, 6275461, 1995285, 527921, 128191, 29150, 5301, 705, 87, 11, 1)
percents.900 = freq.900/all.900
barplot(percents.900, names.arg = 1:14)
