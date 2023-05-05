#!/usr/bin/env Rscript
library(stringr)
d = commandArgs(trailingOnly = T)
input = d[1]
cotrv = d[2]
endir = d[3]
refs = d[4]
refbg = d[5]
tab = read.table(input, stringsAsFactors = F, header = T, sep = "\t", quote = "")
refs = read.table(refs, stringsAsFactors = F, header = F)
if(toupper(cotrv) == "TRUE" || toupper(cotrv) == "T"){
  tab2 = NULL
  for(i in 1:nrow(tab)){
    x = paste(tab$ref[i], tab$alt[i], sep = "")
    if(x != "CT" && x != "TC" && x != "AG" && x != "GA"){
      tab2 = rbind(tab2, tab[i,])
    }
  }
  tab = tab2
}
if(toupper(refbg) == "HG19"){
  tab = tab[,-which(names(tab) %in% c("GRCh38Chromosome", "GRCh38Location"))]
  names(tab)[names(tab) == 'GRCh37Location'] <- 'refpos'
  names(tab)[names(tab) == 'GRCh37Chromosome'] <- 'refchr'
} else if (toupper(refbg) == "HG38"){
  tab = tab[,-which(names(tab) %in% c("GRCh37Chromosome", "GRCh37Location"))]
  names(tab)[names(tab) == 'GRCh38Location'] <- 'refpos'
  names(tab)[names(tab) == 'GRCh38Chromosome'] <- 'refchr'
}
tab = tab[,which(names(tab) %in% c("refchr", "refpos"))]
if(refs[1,1] == "chr1"){
  for(i in 1:22){
    tab[tab == i] = paste("chr", i, sep = "")
  }
}
x = c("x", "X", 23)
for(i in 1:length(x)){
  ez = grep(x[i], refs$V1, value = T)
  if(length(ez) == 1){
    tab[tab == "X"] = ez
  }
}
y = c("y", "Y", 24)
for(i in 1:length(y)){
  ez = grep(y[i], refs$V1, value = T)
  if(length(ez) == 1){
    tab[tab == "Y"] = ez
  }
}
m = c("MT", "M", "MITO", "chrM", "chrMT", "chrMITO", "mitochondrial", "NC_012920.1", "rCRS", "rcrs", "CRS", "crs", "RSRS", "rsrs", "mt", "m", "chrmt", "chrm")
for(i in 1:length(m)){
  ez = grep(m[i], refs$V1, value = T)
  if(length(ez) == 1){
    tab[tab == "MT"] = ez
  }
}
tab$bedin = paste(tab$refchr, tab$refpos, sep = ":")
write.table(tab$bedin, file = paste(endir, "snps.bed", sep = ""), col.names = F, row.names = F, quote = F)





