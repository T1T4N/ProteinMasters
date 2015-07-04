import igraph as ig
import louvain
import time

"""Change paths"""
graph_path = ""
output_path = ""

t1 = time.time()
graph = ig.Graph.Read_Ncol(graph_path)
t2 = time.time()
print "Reading graph finish in %d" % (t2-t1, )

part = louvain.find_partition(graph, method="Modularity", weight="weight")
t3 = time.time()
print "Finding communities finish in %d" % (t3-t2, )

for idx, cluster in enumerate(part):
    print idx, cluster
    graph.vs[cluster]["cluster"] = idx

print ig.VertexSeq(ig).attributes()

f = open(output_path, "w")
for cluster in part:
    f.writelines(" ".join(graph.vs[x]["name"] for x in cluster)+"\n")
