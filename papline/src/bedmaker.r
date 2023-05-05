#!/usr/bin/env Rscript
d = commandArgs(trailingOnly = T)
input = d[1]
endir = d[2]
opt = d[3]
cotrv = d[4]
refs = d[5]
tab = as.data.frame(read.table(input, stringsAsFactors = F, header = F, sep = "\t"))
refs = read.table(refs, stringsAsFactors = F, header = F)
if(toupper(cotrv) == "TRUE" || toupper(cotrv) == "T"){
  tab2 = NULL
  for(i in 1:nrow(tab)){
    x = paste(tab$V5[i], tab$V6[i], sep = "")
    if(x != "CT" && x != "TC" && x != "AG" && x != "GA"){
      tab2 = rbind(tab2, tab[i,])
    }
  }
  tab = tab2
}
x = c("x", "X", 23)
for(i in 1:length(x)){
  ez = grep(x[i], refs$V1, value = T)
  if(length(ez) == 1){
    tab$V2[which(tab$V2 == 23)] = ez
  }
}
y = c("y", "Y", 24)
for(i in 1:length(y)){
  ez = grep(y[i], refs$V1, value = T)
  if(length(ez) == 1){
    tab$V2[which(tab$V2 == 24)] = ez
  }
}
if(refs[1,1] == "chr1"){
  tab$V2 = paste("chr", tab$V2, sep = "")
}
if(opt != "Diploid_GL"){
  write.table(tab[,c(2,4)], file = paste(endir, "genocall.bed", sep =""), row.names = F, col.names = F, quote = F, sep = "\t")
} else {
  write.table(tab[,c(2,4)], file = paste(endir, "genocall.bed", sep =""), row.names = F, col.names = F, quote = F, sep = ":")
}