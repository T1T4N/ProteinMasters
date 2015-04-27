source("create_interaction_graphs.R")

clustering_coef = function(x) {
  mean(transitivity(humanPPI700_graph, type = "local", vids = names(x), isolates = "zero"))
}

degree_700 = degree.distribution(humanPPI700_graph)
degree_900 = degree.distribution(humanPPI900_graph)
plot(degree_700, log = "xy")
plot(degree_900, log = "xy")

pearson_cor_coef_700 = assortativity.degree(humanPPI700_graph, FALSE)
pearson_cor_coef_900 = assortativity.degree(humanPPI900_graph, FALSE)
print(pearson_cor_coef_700)
print(pearson_cor_coef_900)

clustering_spectum_700 = tapply(degree(humanPPI700_graph), as.factor(degree(humanPPI700_graph)), clustering_coef)
clustering_spectum_900 = tapply(degree(humanPPI900_graph), as.factor(degree(humanPPI900_graph)), clustering_coef)
plot(names(clustering_spectum_700), clustering_spectum, xlab = "degree", ylab = "clustering coef")
plot(names(clustering_spectum_900), clustering_spectum, xlab = "degree", ylab = "clustering coef")

# large_component = decompose.graph(humanPPI700_graph)[[1]]
# shortest_paths = shortest.paths(large_component, weights = NA)
# longest_shortes_path = max(shortest_paths)
# large_component = decompose.graph(humanPPI900_graph)[[1]]
# shortest_paths = shortest.paths(large_component, weights = NA)
# longest_shortes_path = max(shortest_paths)
max_shortest_path_700 = 16
max_shortest_path_900 = 15
all_shortest_paths_700 = shortest.paths(humanPPI700_graph, weights = NA)
all_shortest_paths_900 = shortest_paths(humanPPI900_graph, weights = NA)
percent_sp_pairs_700 = vector(length = max_shortest_path_700)
percent_sp_pairs_900 = vector(length = max_shortest_path_900)
for(i in 1:max_shortest_path) {
  percent_sp_pairs_700[i] = sum(all_shortest_paths_700 == i)*100/length(all_shortest_paths_700)
}
for(i in 1:max_shortest_path) {
  percent_sp_pairs_900[i] = sum(all_shortest_paths_900 == i)*100/length(all_shortest_paths_900)
}
plot(percent_sp_pairs_700)
plot(percent_sp_pairs_900)

large_component_700 = decompose.graph(humanPPI700_graph)[[1]]
large_component_900 = decompose.graph(humanPPI900_graph)[[1]]
closeness_vector_700 = -log(closeness(large_component_700, weights = NULL))
closeness_vector_900 = -log(closeness(large_component_900, weights = NULL))
plot(table(round(closeness_vector_700, 2))*100/length(closeness_vector_700))
plot(table(round(closeness_vector_900, 2))*100/length(closeness_vector_900))

betweenness_vector_700 = -log(betweenness(large_component_700, directed = FALSE, weights = NULL, nobigint = FALSE, normalized = TRUE))
betweenness_vector_900 = -log(betweenness(large_component_900, directed = FALSE, weights = NULL, nobigint = FALSE, normalized = TRUE))
plot(table(round(betweenness_vector_700))*100/length(betweenness_vector_700))
plot(table(round(betweenness_vector_900))*100/length(betweenness_vector_900))


