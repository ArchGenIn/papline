#!/usr/bin/env Rscript
d = commandArgs(trailingOnly = T)
endir = d[1]
pos = read.table(d[2], stringsAsFactors = F, header = T)
base = read.table(d[3], stringsAsFactors = F, header = T)
basename = d[4]
coverage = d[5]
if(0 == coverage){
  coverage = 1
}
tab = cbind(pos, base)
colnames(tab) = c("chr", "pos", "depth", "A", "C", "G", "T")
header = paste(">", basename, sep = "")
blank = NULL
#if coverage does not start with position 1, then create Ns for the missing pos
if(as.numeric(tab$pos[1]) > 1){
  beginning = rep("N", (as.numeric(tab$pos[1]) - 1))
} else {
  beginning = NULL
}
blank = c(blank, beginning)
for(i in 1:nrow(tab)){
  if(i < nrow(tab)){
    ns = rep("N", (as.numeric(tab$pos[i + 1]) - as.numeric(tab$pos[i]) - 1))
  }
  dp = as.numeric(tab$depth[i])
  vec = c(as.numeric(tab$A[i]), as.numeric(tab$C[i]), as.numeric(tab$G[i]), as.numeric(tab$T[i]))
  names(vec) = c("A", "C", "G", "T")
  if(dp >= as.numeric(coverage)){
    actpos = names(which.max(vec))
  } else {
    actpos = "N"
  }
  blank = c(blank, ns, actpos)
}
finfa = paste(blank, collapse = "")
out = matrix(NA, 2, 1)
out[1,1] = header
out[2,1] = finfa
write.table(out, file = paste(endir, basename, ".mito.fasta", sep = ""), col.names = F, row.names = F, quote = F)
