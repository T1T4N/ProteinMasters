__author__ = 'goran'

clusterNodes = {}
nodeClusters = {}

graph = "WH"

nodesClustersPath = "data/nodesClustersList/nodesClusters" + graph
similarProteinsFolderPath = "data/similarProteins" + graph

# initializing clusterNodes to contains all nodes which belong to the specific cluster
# key = cluster, value = (set of nodes which belong to cluster)
# and nodeClusters to contain all the clusters in which belongs the specific node
# key = node, value = (set of clusters in which node belongs)
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

# used to store all proteins which share at least one cluster with a given protein
# key = protein, value = (set of proteins which share cluster with protein)
proteinsWithMutualCluster = {}

for i in nodeClusters: # we are looking for all the neighbours of node i
    proteinsWithMutualCluster.setdefault(i, set())
    clusters = nodeClusters[i] # all clusters in which node i belongs
    for j in clusters:#
        for k in clusterNodes[j]:
            if k != i:
                proteinsWithMutualCluster[i].add(str(k))

# writing the results
# each line format is:
# protienId: similarProtein1 similarProtein2 ...(list of all similar proteins separated by one space)
for protein in proteinsWithMutualCluster:
    with open(similarProteinsFolderPath + "/protein_" + str(protein), 'w') as f:
        f.write(" ".join(proteinsWithMutualCluster[protein]))


