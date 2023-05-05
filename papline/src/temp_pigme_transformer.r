d = ("/media/archeogen/Elements/gerberetal2021/clinical/new20220518/")
setwd(d)
c = 0
refbg = "hg19"
fc = NULL
for(file in dir(path = d, pattern = "varcall.csv$")){
  tab = read.table(file, stringsAsFactors = F, header = T, sep = "\t", quote = "")
  pigme = tab[which(tab$label == "PIGME"),]
  name = unlist(strsplit(file, "\\."))[1]
  if(c == 0){
    fintab = as.data.frame(matrix(NA, nrow(pigme), 7))
    fintab$V1 = pigme$temp
    fintab$V2 = pigme$dbSNP.ID
    fintab$V3 = pigme$Gene.s.
    fintab$V4 = pigme$Condition.ref
    fintab$V5 = pigme$Condition.s.
    fintab$V6 = pigme$ref
    fintab$V7 = pigme$alt
    colnames(fintab) = c(paste("Location", refbg, sep = "_"), "SNP_ID", "Gene", "Condition_Ref", "Condition_Alt", "Ref", "Alt")
  }
  c = 1
  fintab$g = pigme$Genotype
  fintab$l = pigme$GL
  names(fintab)[names(fintab) == "g"] = paste(name, "GT", sep = "_")
  names(fintab)[names(fintab) == "l"] = paste(name, "GL", sep = "_")
  clintab = tab[which(tab$altCarrier == "Y"),]
  clintab = clintab[which(clintab$label == "CLINICAL"),]
  names(clintab)[names(clintab) == "Condition.ref"] = "Sample"
  if(nrow(clintab) > 0){clintab$Sample = name}
  fc = rbind(fc, clintab)
}
write.table(fintab, file = "Pigmentation_results.csv", row.names = F, quote = F, sep = "\t")
write.table(fc, file = "Clinical_results.csv", row.names = F, quote = F, sep = "\t")
