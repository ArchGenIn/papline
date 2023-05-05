#!/usr/bin/env Rscript
d = commandArgs(trailingOnly = T)
input = d[1]
seq = toupper(rev(strsplit(input, "")[[1]]))
for(i in 1:length(seq)){
  if(seq[i] == "A"){
    seq[i] = "T"
  } else if(seq[i] == "T"){
    seq[i] = "A"
  } else if(seq[i] == "G"){
    seq[i] = "C"
  } else if(seq[i] == "C"){
    seq[i] = "G"
  }
}
cat(paste(seq, collapse = ""))
