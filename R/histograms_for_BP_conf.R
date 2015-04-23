source("create_interaction_graphs_for_BP_conf.R")

bp_conf_700_largest_component = decompose.graph(bp_conf_700_graph)[[1]]
matrix_shortest_paths = shortest.paths(bp_conf_700_largest_component)
maximum_distance_700 = max(matrix_shortest_paths)

# Ги содржи сите најкратки патеки, треба да се скокне во јава читањето на 
# првата колона и првиот ред, таму се имињата на протеините
# Бидејќи не се битни нема потреба од нив по индекси ќе се раководеме
# оти на наредниот фајл исто така ќе ти ги пратам по ред како што се
# во матрицата
# maximum_distance e 13, што значи ќе ми вратиш вектори со 13 елементи
t1 = Sys.time()
write.table(matrix_shortest_paths, file = "data/matrix_sp_bp_conf_700.txt")
t2 = Sys.time()
print(t2 - t1)

# Едншка се пушта, ако сакаш двапати избришеш фајл поради append
# Отвори го и сега не се рипа колони и редови
t1 = Sys.time()
for(i in V(bp_conf_700_largest_component)) {
  write.table(as.list(V(bp_conf_700_largest_component)[i]$GO[[1]]), file = "data/protein_term_700.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
}
t2 = Sys.time()
print(t2 - t1)

bp_conf_900_largest_component = decompose.graph(bp_conf_900_graph)[[1]]
matrix_shortest_paths = shortest.paths(bp_conf_900_largest_component)
maximum_distance_900 = max(matrix_shortest_paths)

# За 900 maximum_distance е 14
t1 = Sys.time()
write.table(matrix_shortest_paths, file = "data/matrix_sp_bp_conf_900.txt")
t2 = Sys.time()
print(t2 - t1)

t1 = Sys.time()
for(i in V(bp_conf_900_largest_component)) {
  write.table(as.list(V(bp_conf_900_largest_component)[i]$GO[[1]]), file = "data/protein_term_900.txt", append = TRUE, row.names = FALSE, col.names = FALSE)
}
t2 = Sys.time()
print(t2 - t1)

############################################################################
# Очекува се од тебе 2 вектора за првото барање во формат
# c(f(d1), f(d2), f(d3), ... f(d13)) и за 900 до f(d14)
# Исто така 2 вектора за второто барање
# Waiting JAVA processing
############################################################################



###########################################################################
# Ова ми стое како се користе мултикоре
# Ќе се пушта за секоја различна далечина (d)
# pairs - листа која содрже листа за секој протеин
# внатрешната листа e листа на протеини на растојание d и ги
# содржи GO термините
# За 3, 4, 5 пак е бавно, а има и bug и се откажувам
# t1 = Sys.time()
# pairs = foreach(protein = 1:(length(V(bp_conf_700_largest_component))-1)) %dopar% {
#   lists = V(bp_conf_700_largest_component)[which(matrix_shortest_paths[protein, (protein+1):nrow(matrix_shortest_paths)] == 5)]$GO
# }
# t2 = Sys.time()
# print(t2 - t1)
# t1 = Sys.time()
# result = 0
# for(protein in 1:100)  {
#   count = foreach(list = pairs[[protein]], .combine = "+") %dopar% {
#     as.numeric(length(intersect(V(bp_conf_700_largest_component)[protein]$GO[[1]], list))>0)
#   }
#   if(!is.null(count)){
#     result = result + count
#   }
# }
# t2 = Sys.time()
# print(t2 - t1)
# d1 17486
# d2 392221
# d6 55714
# d7 6807
# d8 1974
# d9 1574
# d10 292
# d11 25
# d12 1
# d13 1
