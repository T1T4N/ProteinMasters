source("create_interaction_graphs_for_BP_conf.R")
library(doParallel)
registerDoParallel(cores = 2)
# registerDoSEQ()

jacard = function(T1, T2) {
  T1 = unique(T1)
  T2 = unique(T2)
  intersection_length = length(intersect(T1, T2))
  (intersection_length/length(T1)+intersection_length/length(T2))/2
}

# За алгоритмити е потребно поврзана компонента, 
# ние ја зимаме нормално големата
connected_graph = decompose.graph(graph = bp_conf_700_graph)[[1]]

# Оваа скрипта ќе ја добиеме сите, и ќе поделеме вака работата
# Вкупно има 135621 ребро
#   - Горан 1:45000
#   - Даниел 45001:90000
#   - Роберт 90001:135621

# Со наредните линии се добива матрица во која секој ред
# претставуват едно ребро, а колоните неговите јазли протеини
# Од коментирајте го вашето множество
edges45000 = get.edges(connected_graph, E(connected_graph, directed = FALSE)[1:45000])
edges90000 = get.edges(connected_graph, E(connected_graph, directed = FALSE)[45001:90000])
edgesEND = get.edges(connected_graph, E(connected_graph, directed = FALSE)[90001:ecount(connected_graph)])

# Значи не е потребно да пуштате за целото ваше множество наредните функции,
# можете да го поделете за да можете да гасете лаптопот
# најдобро мислам е да се запишуват резултатите на следниот начин:
#   jacard_content_XXXXX_XXXXXX
#   jacard_structure_XXXX_XXXXX
#   resnik_content_XXXX_XXXX
#   resnik_structure_XXXX_XXXX
#   wang_content_XXXXX_XXXXX
#   wang_structure_XXXXX_XXXXX
# каде XXXXX_XXXX се од кое до кое ребро се наоѓа во тој фајл
# Структурата на фајлот:
#   nameProtein1 nameProtein2 similarity
# Така најлесно ќе се интегрираат на крај

# поставете го ова зависно од вашето множество 

set_starting = 0
currentSet = edges45000

# set_starting = 45000
# currentSet = edges90000

# set_starting = 90000
# # currentSet = edgesEND

# Во currentWorkingEdges ставете на кои јазли ќе пуштете 
# во зависност од вашето множество
# примет ако сум за 45000 поставувам c(1, 2) за 45001 и 45002
currentWorkingEdges = c(15001:30000) 
resnik_dt = data.frame(protein1 = character(), protein2 = character(), score = character(), stringsAsFactors = FALSE)
wang_dt = data.frame(protein1 = character(), protein2 = character(), score = character(), stringsAsFactors = FALSE)
jacard_dt = data.frame(protein1 = character(), protein2 = character(), score = character(), stringsAsFactors = FALSE)

t1 = Sys.time()
for(i in currentWorkingEdges) {
  name1 = V(connected_graph)[currentSet[i, 1]]$name
  name2 = V(connected_graph)[currentSet[i, 2]]$name
  terms1 = V(connected_graph)[currentSet[i, 1]]$GO[[1]]
  terms2 = V(connected_graph)[currentSet[i, 2]]$GO[[1]]
  
  resnik_sim = mgoSim(terms1, terms2, ont = "BP", organism = "human", measure = "Resnik", combine = "rcmax")
  wang_sim = mgoSim(terms1, terms2, ont = "BP", organism = "human", measure = "Wang", combine = "rcmax")
  jacard_sim = jacard(terms1, terms2)
  
  resnik_dt = rbind(resnik_dt, data.frame(protein1 = name1, protein2 = name2, score = resnik_sim, stringsAsFactors = FALSE))
  wang_dt = rbind(wang_dt, data.frame(protein1 = name1, protein2 = name2, score = wang_sim, stringsAsFactors = FALSE))
  jacard_dt = rbind(jacard_dt, data.frame(protein1 = name1, protein2 = name2, score = jacard_sim, stringsAsFactors = FALSE))
}
t2 = Sys.time()
print(t2 - t1)

# t1 = Sys.time()
# resnik_dt = foreach(i = currentWorkingEdges, .combine = rbind) %dopar% {
#   name1 = V(connected_graph)[currentSet[i, 1]]$name
#   name2 = V(connected_graph)[currentSet[i, 2]]$name
#   terms1 = V(connected_graph)[currentSet[i, 1]]$GO[[1]]
#   terms2 = V(connected_graph)[currentSet[i, 2]]$GO[[1]]
#   resnik_sim = mgoSim(terms1, terms2, ont = "BP", organism = "human", measure = "Resnik", combine = "rcmax")
#   data.frame(protein1 = name1, protein2 = name2, score = resnik_sim, stringsAsFactors = FALSE)
# }
# wang_dt = foreach(i = currentWorkingEdges, .combine = rbind) %dopar% {
#   name1 = V(connected_graph)[currentSet[i, 1]]$name
#   name2 = V(connected_graph)[currentSet[i, 2]]$name
#   terms1 = V(connected_graph)[currentSet[i, 1]]$GO[[1]]
#   terms2 = V(connected_graph)[currentSet[i, 2]]$GO[[1]]
#   wang_sim = mgoSim(terms1, terms2, ont = "BP", organism = "human", measure = "Wang", combine = "rcmax")
#   data.frame(protein1 = name1, protein2 = name2, score = wang_sim, stringsAsFactors = FALSE)
# }
# jacard_dt = foreach(i = currentWorkingEdges, .combine = rbind) %dopar% {
#   name1 = V(connected_graph)[currentSet[i, 1]]$name
#   name2 = V(connected_graph)[currentSet[i, 2]]$name
#   terms1 = V(connected_graph)[currentSet[i, 1]]$GO[[1]]
#   terms2 = V(connected_graph)[currentSet[i, 2]]$GO[[1]]
#   jacard_sim = jacard(terms1, terms2)
#   data.frame(protein1 = name1, protein2 = name2, score = jacard_sim, stringsAsFactors = FALSE)
# }
# t2 = Sys.time()
# print(t2 - t1)

# Креирајте folder partial во data
resnik_dt = as.data.table(resnik_dt)
wang_dt = as.data.table(wang_dt)
jacard_dt = as.data.table(jacard_dt)
resnik_content_file = paste(paste("resnik", "content", (currentWorkingEdges[1]+set_starting), (currentWorkingEdges[length(currentWorkingEdges)]+set_starting), sep = "_"), ".txt", sep = "")
wang_content_file =  paste(paste("wang", "content", (currentWorkingEdges[1]+set_starting), (currentWorkingEdges[length(currentWorkingEdges)]+set_starting), sep = "_"), ".txt", sep = "")
jacard_content_file =  paste(paste("jacard", "content", (currentWorkingEdges[1]+set_starting), (currentWorkingEdges[length(currentWorkingEdges)]+set_starting), sep = "_"), ".txt", sep = "")
write.table(resnik_dt, paste("data/partial/", resnik_content_file, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
write.table(wang_dt, paste("data/partial/", wang_content_file, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
write.table(jacard_dt, paste("data/partial/", jacard_content_file, sep = ""), sep = " ", row.names = FALSE, col.names = FALSE)
