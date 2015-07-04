# Set GO terms on HumanPPI data
# Subset data with confidence score for terms above 3
# Split generated sets on 3 parts, BP, CC, MF
#   - keep only informative ties between terms
#   Input: humanPPI700.graph, humanPPI900.graph, 9606_go_knowledge_full, 
#          BPGOterms, CCGOterms, MFGOterms, BPGOfull, CCGOfull, MFGOfull
#   Output: HumanPPI700_GO, HumanPPI900_GO, 
#           HumanPPI700_GO_XX, HumanPPI700_GO_XX_conf
#           HumanPPI900_GO_XX, HumanPPI900_GO_XX_conf | XX in BP, CC, MF
source("graphs/graphs.humanPPI.R")

filter.child = function(term, graph, data) {
  tmp = TRUE
  vertices = V(graph)[neighbors(graph, as.character(term), mode = "in")]$name
  for(vertex in vertices) {
    if(!any(is.na(data[vertex]))) {
      if(data[vertex][1]$count > 30) {
        tmp = FALSE
        break
      } 
    }  
  }
  tmp
}

filter.informative.ties = function(data, graph) {
  data[, count := .N, by = GO]
  tmp = data[count >= 30]
  data = tmp[, if(filter.child(GO, graph, data)) .SD, by = GO]
  data
}

print("Reading GO...")
t1 = Sys.time()
go_knowledge_full = as.data.table(read.table("data/inputGO/9606_go_knowledge_full.tsv", sep="\t", quote = ""))
go_knowledge_full[, V1 := paste(V1, V2, sep = ".")]
go_knowledge_full[, V2 := NULL]
setnames(go_knowledge_full, c(1, 3, 7), c("protein", "GO", "conf"))
setkey(go_knowledge_full, "protein")
t2 = Sys.time()
print(paste("Reading finish in ", t2 - t1))

print("Filtering data...")
t1 = Sys.time()
HumanPPI700_GO = go_knowledge_full[protein %in% V(humanPPI700.graph)$name]
HumanPPI900_GO = go_knowledge_full[protein %in% V(humanPPI900.graph)$name]
HumanPPI700_GO_conf = HumanPPI700_GO[conf >= 3]
HumanPPI900_GO_conf = HumanPPI900_GO[conf >= 3]
setkey(HumanPPI700_GO, GO)
setkey(HumanPPI700_GO_conf, GO)
setkey(HumanPPI900_GO, GO)
setkey(HumanPPI900_GO_conf, GO)
t2 = Sys.time()
print(paste("Filtering finish in ", t2 - t1))

print("Writing to txt...")
t1 = Sys.time()
write.table(HumanPPI700_GO, file = "data/HumanPPI700_GO.txt")
write.table(HumanPPI700_GO_conf, file = "data/HumanPPI700_GO_conf.txt")
write.table(HumanPPI900_GO, file = "data/HumanPPI900_GO.txt")
write.table(HumanPPI900_GO_conf, file = "data/HumanPPI900_GO_conf.txt")
t2 = Sys.time()
print(paste("Writing finish in ", t2 - t1))

print("Generating new sets...")
t1 = Sys.time()
BPGOterms = as.data.table(read.table("data/inputGO/BPGOterms.txt"))
CCGOterms = as.data.table(read.table("data/inputGO/CCGOterms.txt"))
MFGOterms = as.data.table(read.table("data/inputGO/MFGOterms.txt"))
HumanPPI700_GO_CC = HumanPPI700_GO[GO %in% CCGOterms$V1]
HumanPPI700_GO_MF = HumanPPI700_GO[GO %in% MFGOterms$V1]
HumanPPI700_GO_BP = HumanPPI700_GO[GO %in% BPGOterms$V1]
HumanPPI700_GO_CC_conf = HumanPPI700_GO_conf[GO %in% CCGOterms$V1]
HumanPPI700_GO_MF_conf = HumanPPI700_GO_conf[GO %in% MFGOterms$V1]
HumanPPI700_GO_BP_conf = HumanPPI700_GO_conf[GO %in% BPGOterms$V1]
HumanPPI900_GO_CC = HumanPPI900_GO[GO %in% CCGOterms$V1]
HumanPPI900_GO_MF = HumanPPI900_GO[GO %in% MFGOterms$V1]
HumanPPI900_GO_BP = HumanPPI900_GO[GO %in% BPGOterms$V1]
HumanPPI900_GO_CC_conf = HumanPPI900_GO_conf[GO %in% CCGOterms$V1]
HumanPPI900_GO_MF_conf = HumanPPI900_GO_conf[GO %in% MFGOterms$V1]
HumanPPI900_GO_BP_conf = HumanPPI900_GO_conf[GO %in% BPGOterms$V1]
t2 = Sys.time()
print(paste("Generating finish in ", t2 - t1))

print("Processing new datasets...")
t1 = Sys.time()
CCGOfull = as.data.table(read.table("data/inputGO/CCGOfull.txt"))
CCGOfull = CCGOfull[V2 == "is_a"]
neworder = c(colnames(CCGOfull)[1], colnames(CCGOfull)[3], colnames(CCGOfull)[2])
setcolorder(CCGOfull, neworder)
MFGOfull = as.data.table(read.table("data/inputGO/MFGOfull.txt"))
MFGOfull = MFGOfull[V2 == "is_a"]
neworder = c(colnames(MFGOfull)[1], colnames(MFGOfull)[3], colnames(MFGOfull)[2])
setcolorder(MFGOfull, neworder)
BPGOfull = as.data.table(read.table("data/inputGO/BPGOfull.txt"))
BPGOfull = BPGOfull[V2 == "is_a"]
neworder = c(colnames(BPGOfull)[1], colnames(BPGOfull)[3], colnames(BPGOfull)[2])
setcolorder(BPGOfull, neworder)
CC_rel = graph.data.frame(CCGOfull, directed = TRUE)
MF_rel = graph.data.frame(MFGOfull, directed = TRUE)
BP_rel = graph.data.frame(BPGOfull, directed = TRUE)
HumanPPI700_GO_CC = filter.informative.ties(HumanPPI700_GO_CC, CC_rel)
HumanPPI700_GO_CC_conf = filter.informative.ties(HumanPPI700_GO_CC_conf, CC_rel)
HumanPPI900_GO_CC = filter.informative.ties(HumanPPI900_GO_CC, CC_rel)
HumanPPI900_GO_CC_conf = filter.informative.ties(HumanPPI900_GO_CC_conf, CC_rel)
HumanPPI700_GO_MF = filter.informative.ties(HumanPPI700_GO_MF, MF_rel)
HumanPPI700_GO_MF_conf = filter.informative.ties(HumanPPI700_GO_MF_conf, MF_rel)
HumanPPI900_GO_MF = filter.informative.ties(HumanPPI900_GO_MF, MF_rel)
HumanPPI900_GO_MF_conf = filter.informative.ties(HumanPPI900_GO_MF_conf, MF_rel)
HumanPPI700_GO_BP = filter.informative.ties(HumanPPI700_GO_BP, BP_rel)
HumanPPI700_GO_BP_conf = filter.informative.ties(HumanPPI700_GO_BP_conf, BP_rel)
HumanPPI900_GO_BP = filter.informative.ties(HumanPPI900_GO_BP, BP_rel)
HumanPPI900_GO_BP_conf = filter.informative.ties(HumanPPI900_GO_BP_conf, BP_rel)
t2 = Sys.time()
print(paste("Processing finish in ", t2 - t1))

print("Writing sets...")
t1 = Sys.time()
write.table(HumanPPI700_GO_CC, file = "data/HumanPPI700_GO_CC.txt")
write.table(HumanPPI700_GO_CC_conf, file = "data/HumanPPI700_GO_CC_conf.txt")
write.table(HumanPPI900_GO_CC, file = "data/HumanPPI900_GO_CC.txt")
write.table(HumanPPI900_GO_CC_conf, file = "data/HumanPPI900_GO_CC_conf.txt")
write.table(HumanPPI700_GO_MF, file = "data/HumanPPI700_GO_MF.txt")
write.table(HumanPPI700_GO_MF_conf, file = "data/HumanPPI700_GO_MF_conf.txt")
write.table(HumanPPI900_GO_MF, file = "data/HumanPPI900_GO_MF.txt")
write.table(HumanPPI900_GO_MF_conf, file = "data/HumanPPI900_GO_MF_conf.txt")
write.table(HumanPPI700_GO_BP, file = "data/HumanPPI700_GO_BP.txt")
write.table(HumanPPI700_GO_BP_conf, file = "data/HumanPPI700_GO_BP_conf.txt")
write.table(HumanPPI900_GO_BP, file = "data/HumanPPI900_GO_BP.txt")
write.table(HumanPPI900_GO_BP_conf, file = "data/HumanPPI900_GO_BP_conf.txt")
t2 = Sys.time()
print(paste("Writing finish in", t2 - t1))
