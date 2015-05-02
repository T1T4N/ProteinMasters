source("create_interaction_graphs_for_BP_conf.R")
library(Rcpp)
library(GOSemSim)

jacard = function(T1, T2) {
  T1 = unique(T1)
  T2 = unique(T2)
  intersection_length = length(intersect(T1, T2))
  (intersection_length/length(T1)+intersection_length/length(T2))/2
}

connected_graph = decompose.graph(graph = bp_conf_700_graph)[[1]]

# Веќе пресметавме за директните соседи content similarity
jacard_content = as.data.table(read.table("data/jacard_content"))
jacard_content_graph = graph.data.frame(jacard_content, directed = FALSE)
jacard_content_graph = set.edge.attribute(jacard_content_graph, "weight", value = E(jacard_content_graph)$V3)
jacard_content_graph = remove.edge.attribute(jacard_content_graph, "V3")

resnik_content = as.data.table(read.table("data/resnik_content"))
resnik_content_graph = graph.data.frame(resnik_content, directed = FALSE)
resnik_content_graph = set.edge.attribute(resnik_content_graph, "weight", value = E(resnik_content_graph)$V3)
resnik_content_graph = remove.edge.attribute(resnik_content_graph, "V3")

wang_content = as.data.table(read.table("data/wang_content"))
wang_content_graph = graph.data.frame(wang_content, directed = FALSE)
wang_content_graph = set.edge.attribute(wang_content_graph, "weight", value = E(wang_content_graph)$V3)
wang_content_graph = remove.edge.attribute(wang_content_graph, "V3")

calculate_sim = function(p1, p2) {
  terms1 = V(connected_graph)[p1]$GO[[1]]
  terms2 = V(connected_graph)[p2]$GO[[1]]
  if(connected_graph[p1, p2] == 1) {
    j = jacard_content_graph[V(connected_graph)[p1]$name, V(connected_graph)[p2]$name]
    r = resnik_content_graph[V(connected_graph)[p1]$name, V(connected_graph)[p2]$name]
    w = wang_content_graph[V(connected_graph)[p1]$name, V(connected_graph)[p2]$name]
  }
  else {
    j = jacard(terms1, terms2)
    r = mgoSim(terms1, terms2, ont = "BP", organism = "human", measure = "Resnik", combine = "rcmax")
    w = mgoSim(terms1, terms2, ont = "BP", organism = "human", measure = "Wang", combine = "rcmax")
  }
  c(j, r, w)
}

get_names_from_index = function(p1, p2) {
  c(V(connected_graph)[p1]$name, V(connected_graph)[p2]$name)
}

cppFunction('DataFrame findAllSims(NumericVector p, int N, Function goSim, Function names) {
  int sumSet = sum(p);
  int sizeSet = p.size();
  StringVector proteins1(sizeSet*N-sumSet);
  StringVector proteins2(sizeSet*N-sumSet);
  DoubleVector jacard_sim(sizeSet*N-sumSet);
  DoubleVector resnik_sim(sizeSet*N-sumSet);
  DoubleVector wang_sim(sizeSet*N-sumSet);
  int count = 0;
  for(int i = p[0]; i <= p[sizeSet-1]; i++) {
    Rcout << i << std::endl;
    for(int j = i+1; j <= N ; j++) {
      NumericVector sim = as<NumericVector>(goSim(i, j));
      StringVector proteins = as<StringVector>(names(i, j));
      proteins1[count] = proteins[0];
      proteins2[count] = proteins[1];
      jacard_sim[count] = sim[0];
      resnik_sim[count] = sim[1];
      wang_sim[count] = sim[2];
      count++;
    }
  }
  return DataFrame::create(Named("protein1") = proteins1, 
                           Named("protein2") = proteins2, 
                           Named("jacard") = jacard_sim,
                           Named("resnik") = resnik_sim,
                           Named("wang") = wang_sim);
}', showOutput = TRUE)


# Целта е да се добие граф на јазли кои е целосно поврзан
# Вкупно врски = (8843*8843-8843)/2 = 39094903
# На три лаптопа = x13М
# Горно триаголна матрица па се рамномерно распределени
nproteins = vcount(connected_graph)
daniel = 1:1643
goran = 1644:3743
robert = 3743:8843
your_set = daniel
# your_set = goran
# your_set = robert

# working set вие колку вие сакате да си поделете
working_set = 1:10
t1 = Sys.time()
result = findAllSims(your_set[working_set], nproteins, calculate_sim, get_names_from_index)
t2 = Sys.time()
print(t2 - t1)

# Ако досега не сте, креирајте partial во data
filename = paste(paste("W", your_set[working_set][1], your_set[working_set][length(working_set)], sep = "_"), ".txt", sep = "")
write.table(result, paste("data/partial/", filename, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)