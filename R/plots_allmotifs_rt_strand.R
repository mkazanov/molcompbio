library(data.table)
library(gridExtra)
library(ggplot2)

ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/MOTIFS/"
MOTIFS_LIST <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/all_motifs_compl.txt"

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
    dt2[, MutationDensity := MutationCnt/TargetCnt]
    dt2[, DensityTotal := sum(MutationDensity)]
    dt2[, LeadingTotal := sum(LeadingCnt)]
    dt2[, LaggingTotal := sum(LaggingCnt)]
    dt2[, NormalizedDensity := MutationDensity/DensityTotal]
    StrandTotalRatio = dt2[, sum(LaggingCnt)/(sum(LaggingCnt)+sum(LeadingCnt))]
    res <- lm(NormalizedDensity ~ ReplicationBin,dt2)
    
    dtCoef <- rbind(dtCoef, data.frame("Motif"=t,
                                       "Coef"=res$coefficients[[2]],
                                       "Pvalue"=summary(res)$coefficients[2,4],
                                       "StrandTotalRatio"=StrandTotalRatio))
  }
  
  coefCnt <- merge(dtCnt, dtCoef, by = "Motif")
  motifs[, rownum := .I]
  final <- merge(motifs, coefCnt, by = "Motif")
  final <- final[order(rownum)]
  final$Motif <- factor(final$Motif,levels=unique(final$Motif))
  
  p1 <- ggplot(final, aes(x = Motif, y = StrandTotalRatio)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- ggplot(final, aes(x = Motif, y = log(Pvalue,10))) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p3 <- ggplot(final, aes(x = Motif, y = Cnt)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p <- grid.arrange(p1,p2,p3, nrow = 3)
  ggsave(paste0(ROOT_DIR,"plotsStrand/",s,".jpg"), plot = p, device = "jpeg", width = 10)
  
}

