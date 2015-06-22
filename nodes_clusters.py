__author__ = 'goran'

# to read data from
edgeNodesPath = "edges.txt"
edgeClusteringPath = "jc.ncol" # vaka go imenuvav fajlot vo koj se dobivat klasterite dobienie od istoimeniot fajl jc.ncol

# to write data to
nodesClustersPath = "nodesClusters"

# map with key = edgeId and value = node1_id, node2_id
edgeNodes = {}

# loading node ids for each edge id
with open(edgeNodesPath) as f:
    for line in f:
        spl = line.replace('"', '').strip().split(" ")
        edgeNodes[int(spl[0])] = (int(spl[1]), int(spl[2].strip()))

# print edgeNodes

# map of clusters for each node
nodeClusters = {}

# assigning edge clusters ids to each node belonging to the corresponding edge
currCluster = 1 # 1 based indexing of edge clusters
with open(edgeClusteringPath) as f:
    for line in f:
        clustId = str(currCluster)
        spl = line.split(" ")
        # print spl
        for i in spl:
            print i
            (first, second) = edgeNodes[int(float(i))]# maybe should be handled in find_communities.py ???
            nodeClusters.setdefault(first, set())
            nodeClusters.setdefault(second, set())
            nodeClusters[first].add(clustId)
            nodeClusters[second].add(clustId)

        currCluster += 1

with open(nodesClustersPath, 'w') as f:
    for i in nodeClusters.keys():
        f.write(str(i) + ": "  + " ".join(nodeClusters[i]) + "\n")
# print nodeClusters


