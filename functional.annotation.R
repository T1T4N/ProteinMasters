# Script will evaluate direct and clustering approach for protein interaction network
#   Input: graph for direct approach, and cluster members for clustering approach for all graphs
#   Output: results
source("graphs/graphs.BP.conf.R")
library(flux)

whichpart <- function(x, n=30) {
  nx <- length(x)
  p <- nx-n
  xp <- sort(x, partial=p)[p]
  which(x > xp)
}

evaluate = function(graph, method = c("direct", "clustering"), percent = 0.03, path) {
  new.sets = lapply(V(graph)$GO, unique)
  graph = set.vertex.attribute(graph, "GO", value = new.sets)
  
  confusion.matrix = matrix(0, nrow = 4, ncol = 10)
  rownames(confusion.matrix) = c("TP", "TN", "FP", "FN")
  colnames(confusion.matrix) = seq(0, 0.9, by = 0.1)
  
  all.terms = table(unlist(V(graph)$GO))/vcount(graph)
  global.freq = as.vector(all.terms)
  names(global.freq) = names(all.terms)
    
  for(v in V(graph)) {  
    print(v)
    terms = V(graph)[v]$GO[[1]]
    
    # chi square or frequency
    m = vector(length = length(all.terms))
    
    if(method == "direct") {
      aff.scores = scan(file = paste0(path,v,".txt"), quiet = T)
      aff.scores = aff.scores[-v]
      functional.neighborhood = whichpart(aff.scores, n=round(vcount(graph) * percent))
      freq = table(unlist(V(graph)[functional.neighborhood]$GO))
      e = global.freq[names(freq)] * length(functional.neighborhood)
      m = (freq - e)^2/e
      m = (m - min(m))/(max(m) - min(m))
    }
    if(method == "clustering") {
      functional.neighborhood = scan(file = paste0(path,v), quiet = T)
      freq = table(unlist(V(graph)[functional.neighborhood]$GO))
      m = freq / length(functional.neighborhood)
    }
    
    for(w in seq(0, 0.9, by = 0.1)) {
      predicted = names(m[m > w])
      confusion.matrix["TP", as.character(w)] = confusion.matrix["TP", as.character(w)] + length(intersect(terms, predicted))
      confusion.matrix["TN", as.character(w)] = confusion.matrix["TN", as.character(w)] + length(setdiff(names(all.terms), union(terms, predicted)))
      confusion.matrix["FP", as.character(w)] = confusion.matrix["FP", as.character(w)] + length(setdiff(predicted, terms))
      confusion.matrix["FN", as.character(w)] = confusion.matrix["FN", as.character(w)] + length(setdiff(terms, predicted))
    }
  }
  confusion.matrix
}



graph = decompose.graph(bp.conf.700.graph)[[1]]
remove(bp.conf.700.graph, bp.conf.900.graph, humanPPI700.graph, humanPPI900.graph)

# direct method 
sample.graph = read.graph("data/sample.graph.ml", format = "graphml")

for(v in V(sample.graph)) {
  sample.graph = set.vertex.attribute(sample.graph, "GO", V(sample.graph)[v], 
                                      value = V(graph)[V(sample.graph)[v]$name]$GO)
}

path = "data/SRW_wh_intersect/sample/"
results = evaluate(sample.graph, "direct", 0.03, path)

sensitivity = apply(results, 2, function(x) x[1]/(x[1]+x[4]))
fpr = apply(results, 2, function(x) x[3]/(x[3]+x[2]))
sensitivity = c(sensitivity, 1)
names(sensitivity) = seq(0, 1, 0.1)
fpr = c(fpr, 1)
names(fpr) = seq(0, 1, 0.1)
auc(fpr, sensitivity)


#clustering
#JC
path = "data/similarProteins/similarProteinsJC/protein_"
results = evaluate(graph, "clustering", path = path)

sensitivity = apply(results, 2, function(x) x[1]/(x[1]+x[4]))
fpr = apply(results, 2, function(x) x[3]/(x[3]+x[2]))
sensitivity = c(sensitivity, 1)
names(sensitivity) = seq(0, 1, 0.1)
fpr = c(fpr, 1)
names(fpr) = seq(0, 1, 0.1)
auc(fpr, sensitivity)

#JS
path = "data/similarProteins/similarProteinsJS/protein_"
results = evaluate(graph, "clustering", path = path)

sensitivity = apply(results, 2, function(x) x[1]/(x[1]+x[4]))
fpr = apply(results, 2, function(x) x[3]/(x[3]+x[2]))
sensitivity = c(sensitivity, 1)
names(sensitivity) = seq(0, 1, 0.1)
fpr = c(fpr, 1)
names(fpr) = seq(0, 1, 0.1)
auc(fpr, sensitivity)

#JH
path = "data/similarProteins/similarProteinsJH/protein_"
results = evaluate(graph, "clustering", path = path)

sensitivity = apply(results, 2, function(x) x[1]/(x[1]+x[4]))
fpr = apply(results, 2, function(x) x[3]/(x[3]+x[2]))
sensitivity = c(sensitivity, 1)
names(sensitivity) = seq(0, 1, 0.1)
fpr = c(fpr, 1)
names(fpr) = seq(0, 1, 0.1)
auc(fpr, sensitivity)

#RC
path = "data/similarProteins/similarProteinsRC/protein_"
results = evaluate(graph, "clustering", path = path)

sensitivity = apply(results, 2, function(x) x[1]/(x[1]+x[4]))
fpr = apply(results, 2, function(x) x[3]/(x[3]+x[2]))
sensitivity = c(sensitivity, 1)
names(sensitivity) = seq(0, 1, 0.1)
fpr = c(fpr, 1)
names(fpr) = seq(0, 1, 0.1)
auc(fpr, sensitivity)

#RS
path = "data/similarProteins/similarProteinsRS/protein_"
results = evaluate(graph, "clustering", path = path)

sensitivity = apply(results, 2, function(x) x[1]/(x[1]+x[4]))
fpr = apply(results, 2, function(x) x[3]/(x[3]+x[2]))
sensitivity = c(sensitivity, 1)
names(sensitivity) = seq(0, 1, 0.1)
fpr = c(fpr, 1)
names(fpr) = seq(0, 1, 0.1)
auc(fpr, sensitivity)

#RH
path = "data/similarProteins/similarProteinsRH/protein_"
results = evaluate(graph, "clustering", path = path)

sensitivity = apply(results, 2, function(x) x[1]/(x[1]+x[4]))
fpr = apply(results, 2, function(x) x[3]/(x[3]+x[2]))
sensitivity = c(sensitivity, 1)
names(sensitivity) = seq(0, 1, 0.1)
fpr = c(fpr, 1)
names(fpr) = seq(0, 1, 0.1)
auc(fpr, sensitivity)


#WC
path = "data/similarProteins/similarProteinsWC/protein_"
results = evaluate(graph, "clustering", path = path)

sensitivity = apply(results, 2, function(x) x[1]/(x[1]+x[4]))
fpr = apply(results, 2, function(x) x[3]/(x[3]+x[2]))
sensitivity = c(sensitivity, 1)
names(sensitivity) = seq(0, 1, 0.1)
fpr = c(fpr, 1)
names(fpr) = seq(0, 1, 0.1)
auc(fpr, sensitivity)

#WS
path = "data/similarProteins/similarProteinsWS/protein_"
results = evaluate(graph, "clustering", path = path)

sensitivity = apply(results, 2, function(x) x[1]/(x[1]+x[4]))
fpr = apply(results, 2, function(x) x[3]/(x[3]+x[2]))
sensitivity = c(sensitivity, 1)
names(sensitivity) = seq(0, 1, 0.1)
fpr = c(fpr, 1)
names(fpr) = seq(0, 1, 0.1)
auc(fpr, sensitivity)

#WH
path = "data/similarProteins/similarProteinsWH/protein_"
results = evaluate(graph, "clustering", path = path)

sensitivity = apply(results, 2, function(x) x[1]/(x[1]+x[4]))
fpr = apply(results, 2, function(x) x[3]/(x[3]+x[2]))
sensitivity = c(sensitivity, 1)
names(sensitivity) = seq(0, 1, 0.1)
fpr = c(fpr, 1)
names(fpr) = seq(0, 1, 0.1)
auc(fpr, sensitivity)
