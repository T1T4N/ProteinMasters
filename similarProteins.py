__author__ = 'goran'

clusterNodes = {}
nodeClusters = {}

nodesClustersPath = "data/nodesClustersList/nodesClustersJC"
similarProteinsPath = "data/similarProteinsList"

# loading node ids for each edge id
currCluster = 1
with open(nodesClustersPath) as f:
    for line in f:
        spl = line.replace('"', '').strip().split(" ")
        protein = int(spl[0][0:len(spl[0]) - 1])
        nodeClusters.setdefault(protein, set())
        for i in range(1, len(spl)):
            clusterI = int(spl[i])
            clusterNodes.setdefault(clusterI, set())
            clusterNodes[clusterI].add(protein)
            nodeClusters[protein].add(clusterI)

# print clusterNodes
# print nodeClusters

proteinsWithMutualCluster = {}

# print nodeClusters[1]
# print clusterNodes[2]
# print clusterNodes[12]


for i in nodeClusters: # we are looking for all the neighbours of node i
    proteinsWithMutualCluster.setdefault(i, set())
    clusters = nodeClusters[i] # all clusters in which node i belongs
    for j in clusters:#
        for k in clusterNodes[j]:
            if k != i:
                proteinsWithMutualCluster[i].add(str(k))

# writing the results
with open(similarProteinsPath, 'w') as f:
    for protein in proteinsWithMutualCluster:
        f.write(str(protein) + ": "  + " ".join(proteinsWithMutualCluster[protein]) + "\n")


