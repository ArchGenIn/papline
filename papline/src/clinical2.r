#!/usr/bin/env Rscript
library(stringr)
d = commandArgs(trailingOnly = T)
counts = d[1]
pos = d[2]
disis = d[3]
pigme = d[4]
panel = d[5]
refbg = d[6]
refs = d[7]
beagle = d[8]
endir = d[9]
name = d[10]
co = read.table(gzfile(counts), stringsAsFactors = F, header = T)
bases = c("A", "C", "G", "T")
colnames(co) = bases
po = read.table(gzfile(pos), stringsAsFactors = F, header = T)
gl = read.table(gzfile(beagle), stringsAsFactors = F, header = T)
gl[which(gl[,2] == 0),2] = "A"
gl[which(gl[,3] == 0),3] = "A"
gl[which(gl[,2] == 1),2] = "C"
gl[which(gl[,3] == 1),3] = "C"
gl[which(gl[,2] == 2),2] = "G"
gl[which(gl[,3] == 2),3] = "G"
gl[which(gl[,2] == 3),2] = "T"
gl[which(gl[,3] == 3),3] = "T"
colnames(gl) = c("marker", "a1", "a2", "r", "h", "a")
sample = cbind(po, co, gl)
tab = read.table(panel, stringsAsFactors = F, header = T, sep = "\t", quote = "")
refs = read.table(refs, stringsAsFactors = F, header = F)
if(toupper(refbg) == "HG19"){
  fintab = tab[,which(names(tab) %in% c("Gene.s.", "Condition.ref", "Condition.s.", "dbSNP.ID", "GRCh37Chromosome", "GRCh37Location", "ref", "alt", "label"))]
} else if(toupper(refbg) == "HG38"){
  fintab = tab[,which(names(tab) %in% c("Gene.s.", "Condition.ref", "Condition.s.", "dbSNP.ID", "GRCh38Chromosome", "GRCh38Location", "ref", "alt", "label"))]
}
fintab$allCov = NA
fintab$refCov = NA
fintab$altCov = NA
fintab$refCarrier = NA
fintab$altCarrier = NA
fintab$Genotype = NA
fintab$GL = NA
fintab$tmpa1 = NA
fintab$tmpa2 = NA
fintab$tmpr = NA
fintab$tmph = NA
if(refs[1,1] == "chr1"){
  for(i in 1:22){
    fintab[fintab == i] = paste("chr", i, sep = "")
  }
}
x = c("x", "X", 23)
for(i in 1:length(x)){
  ez = grep(x[i], refs$V1, value = T)
  if(length(ez) == 1){
    fintab[fintab == "X"] = ez
  }
}
y = c("y", "Y", 24)
for(i in 1:length(y)){
  ez = grep(y[i], refs$V1, value = T)
  if(length(ez) == 1){
    fintab[fintab == "Y"] = ez
  }
}
m = c("MT", "M", "MITO", "chrM", "chrMT", "chrMITO", "mitochondrial", "NC_012920.1", "rCRS", "rcrs", "CRS", "crs", "RSRS", "rsrs", "mt", "m", "chrmt", "chrm")
for(i in 1:length(m)){
  ez = grep(m[i], refs$V1, value = T)
  if(length(ez) == 1){
    fintab[fintab == "MT"] = ez
  }
}
if(toupper(refbg) == "HG19"){
  fintab$temp = paste(fintab$GRCh37Chromosome, fintab$GRCh37Location, sep = "_")
} else if(toupper(refbg) == "HG38"){
  fintab$temp = paste(fintab$GRCh38Chromosome, fintab$GRCh38Location, sep = "_")
}


for(i in 1:nrow(fintab)){
  x = sample$totDepth[which(sample$marker == fintab$temp[i])]
  if(length(x) != 0){
    fintab$allCov[i] = x
    fintab$refCov[i] = sample[which(sample$marker == fintab$temp[i]), which(colnames(sample) == toupper(fintab$ref[i]))]
    if(any(toupper(fintab$alt[i]) == bases)){
      fintab$altCov[i] = sample[which(sample$marker == fintab$temp[i]), which(colnames(sample) == toupper(fintab$alt[i]))]
    } else {
      if(toupper(fintab$alt[i]) == "R"){
        alts = c("A", "G")
      } else if(toupper(fintab$alt[i]) == "Y"){
        alts = c("C", "T")
      } else if(toupper(fintab$alt[i]) == "S"){
        alts = c("C", "G")
      } else if(toupper(fintab$alt[i]) == "W"){
        alts = c("A", "T")
      } else if(toupper(fintab$alt[i]) == "K"){
        alts = c("G", "T")
      } else if(toupper(fintab$alt[i]) == "M"){
        alts = c("C", "A")
      } else if(toupper(fintab$alt[i]) == "B"){
        alts = c("C", "T", "G")
      } else if(toupper(fintab$alt[i]) == "D"){
        alts = c("A", "T", "G")
      } else if(toupper(fintab$alt[i]) == "H"){
        alts = c("C", "T", "A")
      } else if(toupper(fintab$alt[i]) == "V"){
        alts = c("C", "G", "A")
      }
      for(j in 1:length(alts)){
        z = sample[which(sample$marker == fintab$temp[i]), which(colnames(sample) == alts[j])]
        if(as.numeric(z) > 0){
          fintab$altCov[i] = z
        }
      }
    }
    if(is.na(fintab$altCov[i])){
      fintab$altCov[i] = 0
    } else if(is.na(fintab$refCov[i])){
      fintab$refCov[i] = 0
    }
    if(fintab$label[i] == "CLINICAL"){
      if(as.numeric(fintab$refCov[i]) >= as.numeric(disis)){
        fintab$refCarrier[i] = "Y"
      }
      if(as.numeric(fintab$altCov[i]) >= as.numeric(disis)){
        fintab$altCarrier[i] = "Y"
      }
    } else if(fintab$label[i] == "PIGME"){
      if(as.numeric(fintab$refCov[i]) >= as.numeric(pigme)){
        fintab$refCarrier[i] = "Y"
      }
      if(as.numeric(fintab$altCov[i]) >= as.numeric(pigme)){
        fintab$altCarrier[i] = "Y"
      }
    }
    fintab$tmpa1[i] = sample$a1[which(sample$marker == fintab$temp[i])]
    fintab$tmpa2[i] = sample$a2[which(sample$marker == fintab$temp[i])]
    fintab$tmpr[i] = sample$r[which(sample$marker == fintab$temp[i])]
    fintab$tmph[i] = sample$h[which(sample$marker == fintab$temp[i])]
    if(fintab$tmpr[i] > fintab$tmph[i]){
      fintab$GL[i] = fintab$tmpr[i]
      fintab$Genotype[i] = fintab$tmpa1[i]
    } else {
      fintab$GL[i] = fintab$tmph[i]
      fintab$Genotype[i] = "H"
    }
  }
}
fintab = fintab[,-which(names(fintab) %in% c("tmpa1", "tmpa2", "tmpr", "tmph"))]
write.table(fintab, file = paste(endir, name, ".varcall.csv", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
