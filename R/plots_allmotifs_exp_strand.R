library(data.table)
library(gridExtra)
library(ggplot2)
library(stringr)

ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/MOTIFS_EXP/"
MOTIFS_LIST <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/all_motifs.txt"

motifs <- read.csv(MOTIFS_LIST, header = FALSE)
motifs <- data.table(motifs)
setnames(motifs,c("Motif"))

dir <- list.files(ROOT_DIR,pattern = "ALLMOTIFS")

for(file in dir)
{
  if(str_detect(file,"onestrand")){
    next
  }
  targets <- read.csv(paste0(ROOT_DIR,file), sep='\t')
  if(exists("alltargets")){
    alltargets <- rbind(alltargets,targets)
  } else {
    alltargets <- targets
  }
}

dir <- list.files(ROOT_DIR,pattern = "ALLMUTS")

for(file in dir)
{
  if(str_detect(file,"onestrand")){
    next
  }
  muts <- read.csv(paste0(ROOT_DIR,file), sep='\t')
  if(exists("allmuts")){
    allmuts <- rbind(allmuts,muts)
  } else {
    allmuts <- muts
  }
}

data <- merge(allmuts, alltargets, by = c("Motif","Cancer","Sample","ExpressionBin"))
data <- data.table(data)

for(s in unique(data$Sample))
{
  print(s)
  dt <- data[Sample == s]
  
  dtCoef <- data.frame()
  for(t in unique(dt$Motif))
  {
    dt2 <- dt[Motif==t]
    cnt <- dt2[,sum(MutationCnt.x)]
    dt2[, MutationDensity := MutationCnt.x/MutationCnt.y]
    dt2[, DensityTotal := sum(MutationDensity)]
    dt2[, NormalizedDensity := MutationDensity/DensityTotal]
    StrandTotalRatio = dt2[, sum(PlusStrandConsistent.x)/(sum(PlusStrandConsistent.x)+sum(MinusStrandConsistent.x))]
    res <- lm(NormalizedDensity ~ ExpressionBin,dt2)
    
    dtCoef <- rbind(dtCoef, data.frame("Motif"=t,
                                       "Cnt"=cnt,
                                       "Coef"=res$coefficients[[2]],
                                       "Pvalue"=summary(res)$coefficients[2,4],
                                       "StrandTotalRatio"=StrandTotalRatio))
  }
  
 
  p1 <- ggplot(dtCoef, aes(x = Motif, y = StrandTotalRatio)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- ggplot(dtCoef, aes(x = Motif, y = log(Pvalue,10))) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p3 <- ggplot(dtCoef, aes(x = Motif, y = Cnt)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p <- grid.arrange(p1,p2,p3, nrow = 3)
  ggsave(paste0(ROOT_DIR,"plotsStrand/",s,".jpg"), plot = p, device = "jpeg", width = 10)
  
}

