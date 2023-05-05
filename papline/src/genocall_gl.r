#!/usr/bin/env Rscript
library(stringr)
d = commandArgs(trailingOnly = T)
beagle = d[1]
bedfile = d[2]
samples = d[3]
snpfile = d[4]
eigen = d[5]
runam = d[6]
endir = d[7]
gl = read.table(gzfile(beagle), stringsAsFactors = F, header = T)
tab = read.table(bedfile, stringsAsFactors = F, header = F)
ind = read.table(samples, stringsAsFactors = F, header = F)
snp = read.table(eigen, stringsAsFactors = F, header = F, sep = "\t")
snp2 = snp[,c(5,6)]
ref = cbind(tab, snp2)
colnames(ref) = c("pos", "alt", "ref")
gl[which(gl[,2] == 0),2] = "A"
gl[which(gl[,3] == 0),3] = "A"
gl[which(gl[,2] == 1),2] = "C"
gl[which(gl[,3] == 1),3] = "C"
gl[which(gl[,2] == 2),2] = "G"
gl[which(gl[,3] == 2),3] = "G"
gl[which(gl[,2] == 3),2] = "T"
gl[which(gl[,3] == 3),3] = "T"
tab$V1 = gsub(":", "_", tab$V1)
newcol1 = c("Position", "Allele1", "Allele2")
newcol2 = NULL
for(i in 1:nrow(ind)){
  ez = unlist(strsplit(ind$V1[i], "\\."))[1]
  newcol2 = c(newcol2, ez)
  vec = c(paste(ez, "Major", sep = "_"), paste(ez, "MaMi", sep = "_"), paste(ez, "Minor", sep = "_"))
  newcol1 = c(newcol1, vec)
}
colnames(gl) = newcol1
fin = gl[,c(1:3)]
for(i in 1:nrow(ind)){
  sample = paste(ind$V1[i], "_", sep = "")
  temptab = gl[,grep(sample, newcol1)]
  temptab2 = as.data.frame(matrix(NA, nrow(temptab), 1))
  colnames(temptab2) = ind$V1[i]
  for(j in 1:nrow(temptab)){
    if(temptab[j,1] == temptab[j,2] & temptab[j,1] == temptab[j,3]){
      temptab2[j,1] = 9
    } else {
      gt = which.max(temptab[j,])
      if(gt == 1){
        temptab2[j,1] = fin$Allele1[j]
        if(ref[which(fin$Position[j] == ref$pos),2] == fin$Allele1[j]){
          temptab2[j,1] = 0
        } else {
          temptab2[j,1] = 2
          }
      } else if (gt == 2){
        temptab2[j,1] = 1
      } else if (gt == 3){
        temptab2[j,1] = fin$Allele2[j]
        if(ref[which(fin$Position[j] == ref$pos),2] == fin$Allele2[j]){
          temptab2[j,1] = 0
        } else {
          temptab2[j,1] = 2
          }
      }
    }
  }
  fin = cbind(fin, temptab2)
}
fintab = as.data.frame(matrix(NA, nrow(ref), ncol(fin)))
colnames(fintab) = colnames(fin)
fintab$Position = ref$pos
fintab$Allele1 = ref$alt
fintab$Allele2 = ref$ref
fintab[fintab$Position %in% fin$Position,c(4:ncol(fintab))] = fin[,c(4:ncol(fin))]
fintab[is.na(fintab)] = 9
write.table(fintab[,c(4:ncol(fintab))], file = paste(endir, runam, "_geno.geno", sep = ""), row.names = F, col.names = F, quote = F, sep = "")
write.table(snp, file = paste(endir, runam, "_geno.snp", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
ind$v2 = "U"
ind$v3 = runam
write.table(ind, file = paste(endir, runam, "_geno.ind", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")