library(data.table)
library(gridExtra)
library(ggpubr)
library(stringr)
library(stringi)

crDNA <- function(dna)
{
  return(stri_reverse(chartr("acgtACGT","tgcaTGCA",dna)))
}

ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/MOTIFS_RT/"
MOTIFS_LIST <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/all_motifs.txt"

motifs <- read.csv(MOTIFS_LIST, header = FALSE)
motifs <- data.table(motifs)
setnames(motifs,c("Motif"))

cellCancer <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR/refCellTypesCancers.txt", sep = '\t')
cellCancer <- data.table(cellCancer)
setnames(cellCancer,c("Cancer","CellLine"))

data <- data.frame()
celllines <- data.frame()
for(i in 1:nrow(motifs))
{
  m <- as.character(motifs[i]$Motif)
  dt <- read.csv(paste0(ROOT_DIR,m,"_3.txt"), sep='\t')
  dt <- data.table(dt)
  dt[, Motif:=m]
  data <- rbind(data, dt)
  
  for(cl in c("IMR90","MCF7","NHEK"))
  {
    path <- paste0(ROOT_DIR,m,"_in_RTbins_",cl,".txt")
    dt <- read.csv(path, sep='\t')
    dt <- data.table(dt)
    dt[, Motif:=m]
    dt[, CellLine:=cl]
#    dt[ReplicationBin >=0, ReplicationBin := ReplicationBin - 1]
    celllines <- rbind(celllines,dt)
  }
  
}

data <- merge(data, cellCancer, by = "Cancer")
data <- merge(data, celllines, by = c("Motif","CellLine","ReplicationBin"))

for(s in unique(data$Sample))
{
  print(s)
  dt <- data[Sample == s]
  dtCnt <- dt[,.("Cnt"=sum(MutationCnt)),by=Motif]
  
  dtCoef <- data.frame()
  for(t in unique(dt$Motif))
  {
    dt2 <- dt[Motif==t]
    dt2 <- dt2[ReplicationBin %in% c(0,1,2,3,4,5,6)]
    cnt <- dt2[,sum(MutationCnt)]
    dt2[, MutationDensity := MutationCnt/TargetCnt]
    dt2[, DensityTotal := sum(MutationDensity)]
    dt2[, LeadingTotal := sum(LeadingCnt)]
    dt2[, LaggingTotal := sum(LaggingCnt)]
    dt2[, NormalizedDensity := MutationDensity/DensityTotal]
    StrandTotalRatio = dt2[, sum(LaggingCnt)/(sum(LaggingCnt)+sum(LeadingCnt))]
    res <- lm(NormalizedDensity ~ ReplicationBin,dt2)
    
    dtCoef <- rbind(dtCoef, data.frame("Motif"=paste(t,"/",crDNA(t)),
                                       "Cnt"=cnt,
                                       "Coef"=res$coefficients[[2]],
                                       "Pvalue"=summary(res)$coefficients[2,4],
                                       "StrandTotalRatio"=StrandTotalRatio,
                                       "logPvalue"=log(summary(res)$coefficients[2,4],10)))
  }
  
  dtCoef <- data.table(dtCoef)
  dtCoef[,CoefColor := as.factor(ifelse(Coef >= 0, 2, 1))]
  
  p1 <- ggbarplot(dtCoef, x = "Motif", y = "Coef", palette = "npg", fill = "CoefColor", 
                  color = "white", x.text.angle = 90, 
                  font.x = 10, font.xtickslab = 8, font.y = 10, font.ytickslab = 8)
  p1 <- ggpar(p1,legend = "none")
  p2 <- ggbarplot(dtCoef, x = "Motif", y = "logPvalue", palette = "npg", fill = "CoefColor",
                  color = "white", x.text.angle = 90,
                  font.x = 10, font.xtickslab = 8, font.y = 10, font.ytickslab = 8)
  p2 <- ggpar(p1,legend = "none")
  p3 <- ggbarplot(dtCoef, x = "Motif", y = "Cnt", palette = "npg", fill = "CoefColor", 
                  color = "white", x.text.angle = 90,
                  font.x = 10, font.xtickslab = 8, font.y = 10, font.ytickslab = 8) 
  p3 <- ggpar(p1,legend = "none")
  p <- grid.arrange(p1,p2,p3, nrow = 3)
  ggsave(filename=paste0(ROOT_DIR,"plots/",s,".jpg"), p, height = 8)

  
}

