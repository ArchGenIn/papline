#!/usr/bin/env Rscript
d = commandArgs(trailingOnly = T)
hist = d[1]
endir = d[2]
chrx = d[3]
chry = d[4]
name = d[5]
xref = d[6]
yref = d[7]
base = d[8]
mtdna = d[9]
#reading data
chrvec = c("genome", 1:22, xref, yref)
tab = as.data.frame(matrix(NA, 1, 26))
colnames(tab) = c("Ind", chrvec)
tab$Ind = name
histtab = read.table(hist, stringsAsFactors = F)
if(is.na(mtdna)){
  mtdna = 0
  write.table(mtdna, file = paste(endir, name, ".mtdna_cov.txt", sep = ""), quote = F, row.names = F, col.names = F)
} else {
  sub = histtab[which(mtdna == histtab[,1]),]
  mtdna = (sum(as.numeric(sub[,2]) * as.numeric(sub[,3])))/as.numeric(sub[1,4])
  write.table(mtdna, file = paste(endir, name, ".mtdna_cov.txt", sep = ""), quote = F, row.names = F, col.names = F)
}
#transforming data
for(i in 1:length(chrvec)){
  sub = histtab[which(grep(chrvec[i], histtab[,1], value = T)[1] == histtab[,1]),]
  tab[1,(chrvec[i] == colnames(tab))] = (sum(as.numeric(sub[,2]) * as.numeric(sub[,3])))/as.numeric(sub[1,4])
}
avgcov = tab$genome
write.table(tab$X, file = paste(endir, name, ".chrx_cov.txt", sep = ""), quote = F, row.names = F, col.names = F)
refs = read.table(base, stringsAsFactors = F, header = T, sep = ",")
colnames(refs) = c("Ind", chrvec)
#if sample coverage is larger than minimum coverage from refs, use ZAC, else AR
if(as.numeric(min(refs$genome)) < as.numeric(tab$genome[1]) ){
  finz = rbind(refs, tab)
  for(i in 2:ncol(finz)){
    minz = min(as.numeric(finz[,i]))
    maxz = max(as.numeric(finz[,i]))
    for(j in 1:nrow(finz)){
      finz[j,i] = (as.numeric(finz[j,i]) - minz)/(maxz - minz)
    }
  }
  for(i in 1:nrow(finz)){
    finz[i,c(3:ncol(finz))] = as.numeric(finz[i,c(3:ncol(finz))]) / as.numeric(finz$genome[i])
  }
  finz = finz[which(finz$genome != 0),]
  tab = finz[(finz$Ind == name),]
  finy = finz[!(finz$Ind == name),]
  males = finz[(finz$Ind == "MALE"),]
  females = finz[(finz$Ind == "FEMALE"),]
  #presence of chrY
  yxrat = as.numeric(tab$Y) / as.numeric(tab$X)
  if(yxrat > 1){
    chrypres = "YES"
  } else {
    chrypres = "NO"
  }
  #trisomies 13 patau, 18 edwards, 21 down, X 
  mean13 = mean(as.numeric(finy$`13`))
  sd13 = sd(as.numeric(finy$`13`))
  mean18 = mean(as.numeric(finy$`18`))
  sd18 = sd(as.numeric(finy$`18`))
  mean21 = mean(as.numeric(finy$`21`))
  sd21 = sd(as.numeric(finy$`21`))
  meanX = mean(as.numeric(females$X))
  sdX = sd(as.numeric(females$X))
  meanY = mean(as.numeric(males$Y))
  sdY = mean(as.numeric(males$Y))
  meanXM =  mean(as.numeric(males$X))
  sdXM = sd(as.numeric(males$X))
  if(as.numeric(tab$`13`) <= (mean13 + sd13)){
    patau = "NO"
  } else if(as.numeric(tab$`13`) >= (1.5 * mean13)){
    patau = "YES"
  } else {
    if((as.numeric(tab$`13`) - (mean13 + sd13)) <= ((1.5 * mean13) - as.numeric(tab$`13`))){
      patau = "NO"
    } else {
      patau = "YES"
    }
  }
  if(as.numeric(tab$`18`) <= (mean18 + sd18)){
    edwards = "NO"
  } else if(as.numeric(tab$`18`) >= (1.5 * mean18)){
    edwards = "YES"
  } else {
    if((as.numeric(tab$`18`) - (mean18 + sd18)) <= ((1.5 * mean18) - as.numeric(tab$`18`))){
      edwards = "NO"
    } else {
      edwards = "YES"
    }
  }
  if(as.numeric(tab$`21`) <= (mean21 + sd21)){
    down = "NO"
  } else if(as.numeric(tab$`21`) >= (1.5 * mean21)){
    down = "YES"
  } else {
    if((as.numeric(tab$`21`) - (mean21 + sd21)) <= ((1.5 * mean21) - as.numeric(tab$`21`))){
      down = "NO"
    } else {
      down = "YES"
    }
  }
  if(chrypres == "NO"){
    if(as.numeric(tab$X) <= (meanX + sdX)){
      xtri = "NO"
    } else if(as.numeric(tab$X) >= (1.5 * meanX)){
      xtri = "YES"
    } else {
      if((as.numeric(tab$X) - (meanX + sdX)) <= ((1.5 * meanX) - as.numeric(tab$X))){
        xtri = "NO"
      } else {
        xtri = "YES"
      }
    }
  } else {
    xtri = "M"
  }
  #monosomy X female turner
  if(chrypres == "NO"){
    if(as.numeric(tab$X) >= (meanX - sdX)){
      xmon = "NO"
    } else if(as.numeric(tab$X) <= (0.5 * meanX)){
      xmon = "YES"
    } else {
      if(((meanX - sdX) - as.numeric(tab$X)) <= (as.numeric(tab$X) - (0.5 * meanX))){
        xmon = "NO"
      } else {
        xmon = "YES"
      }
    }
  } else {
    xmon = "M"
  }
  #disomy Y male jacobs
  if(chrypres == "YES"){
    if(as.numeric(tab$Y) <= (meanY + sdY)){
      jacobs = "NO"
    } else if(as.numeric(tab$Y) >= (2 * meanY)){
      jacobs = "YES"
    } else {
      if((as.numeric(tab$Y) - (meanY + sdY)) <= ((2 * meanY) - as.numeric(tab$Y))){
        jacobs = "NO"
      } else {
        jacobs = "YES"
      }
    }
  } else {
    jacobs = "F"
  }
  #XXY kleinefelter
  if(chrypres == "YES"){
    if(as.numeric(tab$X) <= (meanXM + sdXM)){
      klein = "NO"
    } else if(as.numeric(tab$X) >= (2 * meanXM)){
      klein = "YES"
    } else {
      if((as.numeric(tab$X) - (meanXM + sdXM)) <= ((2 * meanXM) - as.numeric(tab$X))){
        klein = "NO"
      } else {
        klein = "YES"
      }
    }
  } else {
    klein = "F"
  }
  tab$Average_Coverage = avgcov
  tab$mtDNA_Coverage = mtdna
  tab$chrY_present = chrypres
  tab$Patau = patau
  tab$Edwards = edwards
  tab$Down = down
  tab$X_trisomy = xtri
  tab$X_monosomy = xmon
  tab$Jacobs = jacobs
  tab$Kleinefelter = klein
  tab$AR = "ZAC"
  write.table(tab, file = paste(endir, name, ".coverages.txt", sep = ""))
  if((chrypres == "YES") & (klein == "NO")){
    sex = "MALE"
    write.table(sex, file = paste(endir, name, ".sex.txt", sep = ""), quote = F, row.names = F, col.names = F)
  } else {
    sex = "NOT_MALE"
    write.table(sex, file = paste(endir, name, ".sex.txt", sep = ""), quote = F, row.names = F, col.names = F)
  }
} else {
  if(as.numeric(chrx) < 100){
    ar = "LowCov"
    message("Genomic coverages are too low to infer genetic sex")
  } else {
    ar = (log10((as.numeric(chry)+1)*2.6)) / log10(((as.numeric(chry)+1)*2.6)*(as.numeric(chrx)+2.6))
    message("Genomic coverages are too low to infer aneuploidies, AR method is used to infer the presence of chrY")
  }
  tab$Average_Coverage = avgcov
  tab$mtDNA_Coverage = mtdna
  tab$chrY_present = "NA"
  tab$Patau = "NA"
  tab$Edwards = "NA"
  tab$Down = "NA"
  tab$X_trisomy = "NA"
  tab$X_monosomy = "NA"
  tab$Jacobs = "NA"
  tab$Kleinefelter = "NA"
  tab$AR = ar
  write.table(tab, file = paste(endir, name, ".coverages.txt", sep = ""))
}
