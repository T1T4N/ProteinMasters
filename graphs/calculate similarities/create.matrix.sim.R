# It will calculate similarities of all pairs of protein in bp.conf.700
# This script is shared between computers to finish earlier
# Sets for all 3 computers are our names
#   Input: bp.conf.700
#   Output: partial W_sim files
# *After all similarities are computed partial files are integrated in 
#  one file W_sim
source("graphs/graphs.BP.conf.R")
library(Rcpp)
library(GOSemSim)

connected.graph = decompose.graph(graph = bp.conf.700.graph)[[1]]

jacard = function(T1, T2) {
  T1 = unique(T1)
  T2 = unique(T2)
  intersection.length = length(intersect(T1, T2))
  (intersection.length/length(T1)+intersection.length/length(T2))/2
}

calculate.sim = function(p1, p2) {
  terms1 = V(connected.graph)[p1]$GO[[1]]
  terms2 = V(connected.graph)[p2]$GO[[1]]
  j = jacard(terms1, terms2)
  r = mgoSim(terms1, terms2, ont = "BP", organism = "human", measure = "Resnik", combine = "rcmax")
  w = mgoSim(terms1, terms2, ont = "BP", organism = "human", measure = "Wang", combine = "rcmax")
  c(j, r, w)
}

get.names.from.index = function(p1, p2) {
  c(V(connected.graph)[p1]$name, V(connected.graph)[p2]$name)
}

cppFunction('DataFrame findAllSims(NumericVector p, int N, Function goSim, Function names) {
  int sumSet = sum(p);
  int sizeSet = p.size();
  StringVector proteins1(sizeSet*N-sumSet);
  StringVector proteins2(sizeSet*N-sumSet);
  DoubleVector jacardSim(sizeSet*N-sumSet);
  DoubleVector resnikSim(sizeSet*N-sumSet);
  DoubleVector wangSim(sizeSet*N-sumSet);
  int count = 0;
  for(int i = p[0]; i <= p[sizeSet-1]; i++) {
    Rcout << i << std::endl;
    for(int j = i+1; j <= N ; j++) {
      NumericVector sim = as<NumericVector>(goSim(i, j));
      StringVector proteins = as<StringVector>(names(i, j));
      proteins1[count] = proteins[0];
      proteins2[count] = proteins[1];
      jacardSim[count] = sim[0];
      resnikSim[count] = sim[1];
      wangSim[count] = sim[2];
      count++;
    }
  }
  return DataFrame::create(Named("protein1") = proteins1, 
                           Named("protein2") = proteins2, 
                           Named("jacard") = jacardSim,
                           Named("resnik") = resnikSim,
                           Named("wang") = wangSim);
}', showOutput = TRUE)

nproteins = vcount(connected.graph)

daniel = 1:1643
goran = 1644:3742
robert = 3743:8843
set = daniel
# set = goran
# set = robert

working.set = 1:10
t1 = Sys.time()
result = findAllSims(set[working.set], nproteins, calculate.sim, get.names.from.index)
t2 = Sys.time()
print(t2 - t1)
result = as.data.table(result)

filename = paste(paste("W", set[working.set][1], set[working.set][length(working.set)], sep = "_"), ".txt", sep = "")
write.table(result, paste("data/partial/", filename, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

remove(nproteins, daniel, goran, robert, 
       set, working.set, t1, t2, result, filename)