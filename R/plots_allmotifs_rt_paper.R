library(data.table)
library(gridExtra)
library(ggpubr)
library(stringr)
library(stringi)
library(svglite)

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
  if(!(s %in% c("TCGA-BT-A3PH-01A","TCGA-E2-A1LG-01A","TCGA-CV-6961-01A","TCGA-05-4422-01A","TCGA-60-2698-01A"))){
    next
  }
     
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
  
    if(s == "TCGA-BT-A3PH-01A" && t == "TGA"){
      p <- ggplot(dt2, aes(x=ReplicationBin, y=NormalizedDensity)) + 
        geom_bar(stat="identity", width = 0.7, color=rgb(212,42,47,maxColorValue = 255), 
                 fill=rgb(233,148,151,maxColorValue = 255)) +
        xlab("Replication timing") +
        ylab("Normalized mutation density") +
        geom_smooth(method='lm',formula=y~x) +
        theme(panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.position = "none",
              axis.text = element_text(size=16),
              axis.title = element_text(size=16)) +
        scale_x_continuous(
        breaks=c(0,1,2,3,4,5,6),
        labels=c("late","","","","","", "early")
      )
      ggsave(filename=paste0("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/paper/pics/RT_positive_slope.svg"), p)
      
    }
      
    if(s == "TCGA-BT-A3PH-01A" && t == "TAA"){
      p <- ggplot(dt2, aes(x=ReplicationBin, y=NormalizedDensity)) + 
        geom_bar(stat="identity", width = 0.7, color=rgb(212,42,47,maxColorValue = 255), 
                 fill=rgb(233,148,151,maxColorValue = 255)) +
        xlab("Replication timing") +
        ylab("Normalized mutation density") +
        geom_smooth(method='lm',formula=y~x) +
        theme(panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.position = "none",
              axis.text = element_text(size=16),
              axis.title = element_text(size=16)) +
        scale_x_continuous(
          breaks=c(0,1,2,3,4,5,6),
          labels=c("late","","","","","", "early")
        )
      ggsave(filename=paste0("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/paper/pics/RT_negative_slope.svg"), p)
      
    }    
    
    dtCoef <- rbind(dtCoef, data.frame("Motif"=paste0(t,"/",crDNA(t)),
                                       "Cnt"=cnt,
                                       "Coef"=res$coefficients[[2]],
                                       "Pvalue"=summary(res)$coefficients[2,4],
                                       "StrandTotalRatio"=StrandTotalRatio,
                                       "logPvalue"=log(summary(res)$coefficients[2,4],10)))
  }
  
  dtCoef <- data.table(dtCoef)
  dtCoef[,CoefColor := as.factor(ifelse(Coef >= 0, 2, 1))]
  
  p1 <- ggplot(dtCoef, aes(x = Motif, y = Coef, fill=CoefColor) ) + geom_bar(stat="identity", width = 0.8) +
    xlab("Cancer type") +
    ylab("Slope of mutation distribution along replication timing")
  p1 <- p1 + theme(panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.text.x = element_text(angle = 90),
                   legend.position = "none")
  ggsave(filename=paste0("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/paper/pics/",s,".svg"), p1, height = 3)

}

