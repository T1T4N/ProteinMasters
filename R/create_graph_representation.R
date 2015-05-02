source("create_interaction_graphs_for_BP_conf.R")
library(GOSemSim)

# НЕ ГО ГЛЕДАЈТЕ ОВАЈ ФАЈЛ, ТУКА СЕ СМЕТАШЕ CONTENT И ПОЧНАВ СО СТРУКТУРНО. НО НЕ ИСПРАВНО ТАКА
jacard = function(T1, T2) {
  T1 = unique(T1)
  T2 = unique(T2)
  intersection_length = length(intersect(T1, T2))
  (intersection_length/length(T1)+intersection_length/length(T2))/2
}

connected_graph = decompose.graph(graph = bp_conf_700_graph)[[1]]

# Со наредните линии се добива матрица во која секој ред
# претставуват едно ребро, а колоните неговите јазли протеини
# Од коментирајте го вашето множество
edges45000 = get.edges(connected_graph, E(connected_graph, directed = FALSE)[1:45000])
edges90000 = get.edges(connected_graph, E(connected_graph, directed = FALSE)[45001:90000])
edgesEND = get.edges(connected_graph, E(connected_graph, directed = FALSE)[90001:ecount(connected_graph)])

# set_starting = 0
# currentSet = edges45000
# set_starting = 45000
# currentSet = edges90000
# set_starting = 90000
# currentSet = edgesEND

# currentWorkingEdges = c(30001:45000) 
# resnik_dt = data.frame(protein1 = character(), protein2 = character(), score = character(), stringsAsFactors = FALSE)
# wang_dt = data.frame(protein1 = character(), protein2 = character(), score = character(), stringsAsFactors = FALSE)
# jacard_dt = data.frame(protein1 = character(), protein2 = character(), score = character(), stringsAsFactors = FALSE)
# 
# t1 = Sys.time()
# for(i in currentWorkingEdges) {
#   name1 = V(connected_graph)[currentSet[i, 1]]$name
#   name2 = V(connected_graph)[currentSet[i, 2]]$name
#   terms1 = V(connected_graph)[currentSet[i, 1]]$GO[[1]]
#   terms2 = V(connected_graph)[currentSet[i, 2]]$GO[[1]]
#   
#   resnik_sim = mgoSim(terms1, terms2, ont = "BP", organism = "human", measure = "Resnik", combine = "rcmax")
#   wang_sim = mgoSim(terms1, terms2, ont = "BP", organism = "human", measure = "Wang", combine = "rcmax")
#   jacard_sim = jacard(terms1, terms2)
#   
#   resnik_dt = rbind(resnik_dt, data.frame(protein1 = name1, protein2 = name2, score = resnik_sim, stringsAsFactors = FALSE))
#   wang_dt = rbind(wang_dt, data.frame(protein1 = name1, protein2 = name2, score = wang_sim, stringsAsFactors = FALSE))
#   jacard_dt = rbind(jacard_dt, data.frame(protein1 = name1, protein2 = name2, score = jacard_sim, stringsAsFactors = FALSE))
# }
# t2 = Sys.time()
# print(t2 - t1)

# Креирајте folder partial во data
# resnik_dt = as.data.table(resnik_dt)
# wang_dt = as.data.table(wang_dt)
# jacard_dt = as.data.table(jacard_dt)
# resnik_content_file = paste(paste("resnik", "content", (currentWorkingEdges[1]+set_starting), (currentWorkingEdges[length(currentWorkingEdges)]+set_starting), sep = "_"), ".txt", sep = "")
# wang_content_file =  paste(paste("wang", "content", (currentWorkingEdges[1]+set_starting), (currentWorkingEdges[length(currentWorkingEdges)]+set_starting), sep = "_"), ".txt", sep = "")
# jacard_content_file =  paste(paste("jacard", "content", (currentWorkingEdges[1]+set_starting), (currentWorkingEdges[length(currentWorkingEdges)]+set_starting), sep = "_"), ".txt", sep = "")
# write.table(resnik_dt, paste("data/partial/", resnik_content_file, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
# write.table(wang_dt, paste("data/partial/", wang_content_file, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
# write.table(jacard_dt, paste("data/partial/", jacard_content_file, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)

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

# set_starting = 0
# currentSet = edges45000
# set_starting = 45000
# currentSet = edges90000
# set_starting = 90000
# currentSet = edgesEND

currentWorkingEdges = c(101:15000) 
resnik_dt = data.frame(protein1 = character(), protein2 = character(), score = character(), stringsAsFactors = FALSE)
wang_dt = data.frame(protein1 = character(), protein2 = character(), score = character(), stringsAsFactors = FALSE)
jacard_dt = data.frame(protein1 = character(), protein2 = character(), score = character(), stringsAsFactors = FALSE)
t1 = Sys.time()
for(e in currentWorkingEdges) {
  name1 = V(connected_graph)[currentSet[e, 1]]$name
  name2 = V(connected_graph)[currentSet[e, 2]]$name
  terms1 = V(connected_graph)[currentSet[e, 1]]$GO[[1]]
  terms2 = V(connected_graph)[currentSet[e, 2]]$GO[[1]]
  neighbors1 = neighbors(connected_graph, currentSet[e, 1])
  neighbors2 = neighbors(connected_graph, currentSet[e, 2])
  
  jacard_firstSum = 0
  resnik_firstSum = 0
  wang_firstSum = 0
  for(i in neighbors2) {
    nameNeighbor = V(connected_graph)[i]$name
    termsNeighbor = V(connected_graph)[i]$GO[[1]]
    if(i != currentSet[e, 1]) {
      if(connected_graph[name1, nameNeighbor] == 1) {
        jacard_firstSum = jacard_firstSum + jacard_content_graph[name1, nameNeighbor]
        resnik_firstSum = resnik_firstSum + resnik_content_graph[name1, nameNeighbor]
        wang_firstSum = wang_firstSum + wang_content_graph[name1, nameNeighbor]
      }
      else {
        jacard_firstSum = jacard_firstSum + jacard(terms1, termsNeighbor)
        resnik_firstSum = resnik_firstSum + mgoSim(terms1, termsNeighbor, ont = "BP", organism = "human", measure = "Resnik", combine = "rcmax")
        wang_firstSum = wang_firstSum + mgoSim(terms1, termsNeighbor, ont = "BP", organism = "human", measure = "Wang", combine = "rcmax")
      }
    }
  }
  jacard_firstSum = jacard_firstSum/(length(neighbors2)-1)
  resnik_firstSum = resnik_firstSum/(length(neighbors2)-1)
  wang_firstSum = wang_firstSum/(length(neighbors2)-1)
  

  jacard_secondSum = 0
  resnik_secondSum = 0
  wang_secondSum = 0
  for(i in neighbors1) {
    nameNeighbor = V(connected_graph)[i]$name
    termsNeighbor = V(connected_graph)[i]$GO[[1]]
    if(i != currentSet[e, 2]) {
      if(connected_graph[name2, nameNeighbor] == 1) {
        jacard_secondSum = jacard_secondSum + jacard_content_graph[name2, nameNeighbor]
        resnik_secondSum = resnik_secondSum + resnik_content_graph[name2, nameNeighbor]
        wang_secondSum = wang_secondSum + wang_content_graph[name2, nameNeighbor]
      }
      else {
        jacard_secondSum = jacard_secondSum + jacard(terms2, termsNeighbor)
        resnik_secondSum = resnik_secondSum + mgoSim(terms2, termsNeighbor, ont = "BP", organism = "human", measure = "Resnik", combine = "rcmax")
        wang_secondSum = wang_secondSum + mgoSim(terms2, termsNeighbor, ont = "BP", organism = "human", measure = "Wang", combine = "rcmax")
      }
    }
  }
  jacard_secondSum = jacard_secondSum/(length(neighbors1)-1)
  resnik_secondSum = resnik_secondSum/(length(neighbors1)-1)
  wang_secondSum = wang_secondSum/(length(neighbors1)-1)

  jacard_sim = (jacard_firstSum+jacard_secondSum)/2
  resnik_sim = (resnik_firstSum+resnik_secondSum)/2
  wang_sim = (wang_firstSum+wang_secondSum)/2
  
  resnik_dt = rbind(resnik_dt, data.frame(protein1 = name1, protein2 = name2, score = resnik_sim, stringsAsFactors = FALSE))
  wang_dt = rbind(wang_dt, data.frame(protein1 = name1, protein2 = name2, score = wang_sim, stringsAsFactors = FALSE))
  jacard_dt = rbind(jacard_dt, data.frame(protein1 = name1, protein2 = name2, score = jacard_sim, stringsAsFactors = FALSE))
}
t2 = Sys.time()
print(t2-t1)

resnik_dt = as.data.table(resnik_dt)
wang_dt = as.data.table(wang_dt)
jacard_dt = as.data.table(jacard_dt)
resnik_structure_file = paste(paste("resnik", "structure", (currentWorkingEdges[1]+set_starting), (currentWorkingEdges[length(currentWorkingEdges)]+set_starting), sep = "_"), ".txt", sep = "")
wang_structure_file =  paste(paste("wang", "structure", (currentWorkingEdges[1]+set_starting), (currentWorkingEdges[length(currentWorkingEdges)]+set_starting), sep = "_"), ".txt", sep = "")
jacard_structure_file =  paste(paste("jacard", "structure", (currentWorkingEdges[1]+set_starting), (currentWorkingEdges[length(currentWorkingEdges)]+set_starting), sep = "_"), ".txt", sep = "")
write.table(resnik_dt, paste("data/partial/", resnik_structure_file, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
write.table(wang_dt, paste("data/partial/", wang_structure_file, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
write.table(jacard_dt, paste("data/partial/", jacard_structure_file, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
