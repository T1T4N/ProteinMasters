library(data.table)

read.db = function(x, sep = "") {
  initial = read.table(x, sep, header = TRUE, nrow = 10)
  classes = sapply(initial, class)
  tmp = read.table(x, sep, header = TRUE, colClasses = classes, comment.char = "")
  tmp = as.data.table(tmp)
  tmp
}

make_interaction_dbs = function(data) {
  setnames(data, 1, colnames(entrez_mapping_db)[1])
  data = merge(data, entrez_mapping_db, by = colnames(entrez_mapping_db)[1])
  setnames(data, c(1, 2), c("First", colnames(entrez_mapping_db)[1]))
  data = merge(data, entrez_mapping_db, by = colnames(entrez_mapping_db)[1])
  data[, X.Entrez_Gene_ID := NULL]
  data[, First := NULL]
  setnames(data, c(1, 2), c("protein1", "protein2"))
  data
}

calculate_score = function(x, dataname) {
  neigh = 1 - x$neighborhood/1000
  fus = 1 - x$fusion/1000
  cooc = 1 - x$cooccurence/1000
  coex = 1 - x$coexpression/1000
  exp = 1 - x$experimental/1000
  db = 1 - x$database/1000
  txt = 1 - x$textmining/1000
    
  value = 0; hi = 1; lit = 1; ven = 1; yu = 1;
  
  if(dataname == "hi") {
    value = 950
    hi = 1 - value/1000
  } else if(dataname == "lit") {
    value = 900
    hi = 1 - x$hi/1000
    lit = 1 - value/1000
  } else if(dataname == "venkatesan") {
    value = 850
    hi = 1 - x$hi/1000
    lit = 1 - x$lit/1000
    ven = 1 - value/1000
  } else if(dataname == "yu") {
    value = 850
    hi = 1 - x$hi/1000
    lit = 1 - x$lit/1000
    ven = 1 - x$venkatesan/1000
    yu = 1 - value/1000
  }
  res = (1 - neigh * fus * cooc * coex * exp * db * txt * ven * yu * lit * hi) * 1000 
  list(value, as.integer(res))
}

integrate_dbs = function(data, name) {
  t4 = Sys.time()
  print(paste("Modifying", name))
  for(i in 1:nrow(data)) {
    string_db[list(as.character(data[i]$protein1), as.character(data[i]$protein2)), c(eval(name), "combined_score") := calculate_score(.SD, eval(name)), by = .EACHI]
    string_db[list(as.character(data[i]$protein2), as.character(data[i]$protein1)), c(eval(name), "combined_score") := calculate_score(.SD, eval(name)), by = .EACHI]
  }
  string_db[is.na(string_db[[eval(name)]]), eval(name) := 0]
  t5 = Sys.time()
  print(paste(name, "finish in", t5 - t4))
}

t1 = Sys.time()
print("Reading from txt...")
string = read.db("data/input/9606.protein.links.detailed.v9.1.txt.gz") 
entrez = read.db("data/input/entrez_gene_id.vs.string.v9.05.28122012.txt")
hi = read.db("data/input/HI-II-14.tsv", sep = "\t")
lit = read.db("data/input/Lit-BM-13.tsv", sep = "\t")
venkatesan = read.db("data/input/Venkatesan-09.tsv", sep = "\t")
yu = read.db("data/input/Yu-11.tsv", sep = "\t")
hi[, c(Symbol.A, Symbol.B) := NULL]
lit[, c(symbol_a, symbol_b) := NULL]
venkatesan[, c(DB.CCSB.ORF.ID, DB.Accession, AD.CCSB.ORF.ID, AD.Accession, Trial) := NULL]
yu[, c(Symbol.A, Symbol.B) := NULL]
t2 = Sys.time()
print(paste("Reading finish in ", t2 - t1))

print("Filtering interaction dbs...")
hi_db = unique(make_interaction_dbs(hi_db))
lit_db = unique(make_interaction_dbs(lit_db))
venkatesan_db = unique(make_interaction_dbs(venkatesan_db))
yu_db = unique(make_interaction_dbs(yu_db))
setkey(string_db, "protein1", "protein2")
setkey(hi_db, "protein1", "protein2")
setkey(lit_db, "protein1", "protein2")
setkey(venkatesan_db, "protein1", "protein2")
setkey(yu_db, "protein1", "protein2")
t3 = Sys.time()
print(paste("Filtering finish in ", t3 - t2))

print("Calculate scores...")
integrate_dbs(hi_db, "hi")
integrate_dbs(lit_db, "lit")
integrate_dbs(venkatesan_db, "venkatesan")
integrate_dbs(yu_db, "yu")
t4 = Sys.time()
print(paste("Modifying finish in ", t4 - t3))

print("Reorder columns...")
neworder = c(colnames(string_db)[1:9], colnames(string_db)[11:14], colnames(string_db)[10])
setcolorder(string_db, neworder)
t5 = Sys.time()
print(paste("Reordering finish in ", t5 - t4))

print("Writing to txt...")
t6 = Sys.time()
write.table(string_db, file = "data/HumanPPI.txt")
write.table(string_db[combined_score > 700], file = "data/HumanPPI700.txt")
write.table(string_db[combined_score > 900], file = "data/HumanPPI900.txt")
t7 = Sys.time()
print(paste("Writing finish in ", t7 - t6))




