source("create_interaction_graphs.R")

clustering_coef = function(x) {
  mean(transitivity(humanPPI700_graph, type = "local", vids = names(x), isolates = "zero"))
}

degree_700 = degree.distribution(humanPPI700_graph)
degree_900 = degree.distribution(humanPPI900_graph)
plot(degree_700, log = "xy")
plot(degree_900, log = "xy")

pearson_cor_coef = assortativity.degree(humanPPI700_graph, FALSE)

clustering_spectum = tapply(degree(humanPPI700_graph), as.factor(degree(humanPPI700_graph)), clustering_coef)
plot(names(clustering_spectum), clustering_spectum)

all_shortest_paths = shortest.paths(humanPPI700_graph, weights = NA)

#Премногу долго, заглавува ми R, успеав едншка но има Inf поради неповрзаност
#затоа го декомпозирам на подграфовите, a лесно ќе се одреде spectum
t1 = Sys.time()
matrices = list
for(g in decompose.graph(humanPPI700_graph)) {
  c(matrices, shortest.paths(g, weights = NA))
}
t2 = Sys.time()
print(t2 - t1)

#Премногу различни вредности, не знам како spectum да одредам, 14000 одредено
#со table(as.factor(closeness_vector)), 3min
t1 = Sys.time()
closeness_vector = closeness(humanPPI700_graph, weights = NULL)
t2 = Sys.time()
print(t2 - t1)

#исто како за closeness, 4min
t1 = Sys.time()
betweenness_vector = betweenness(humanPPI700_graph, directed = FALSE, weights = NULL)
t2 = Sys.time()
print(t2 - t1)

