library(data.table)
library(ggpubr)
library(ggplot2)
library(rlist)
library(reshape2)
library(MASS)
OUTPUT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalRall/RTEXP/"
STAT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/paper/pics/Fig5/"

samplesEnrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results/sample_enrichment.txt", sep='\t')
samplesEnrichment <- data.table(samplesEnrichment)
samplesEnrichment$SAMPLE <- as.character(samplesEnrichment$SAMPLE)
cancers <- unique(samplesEnrichment[, CANCER_TYPE])
cancers <- data.table(cancers)
setnames(cancers,c("CANCER_TYPE"))
cancers$CANCER_TYPE <- as.character(cancers$CANCER_TYPE)

ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL/"
MOTIFS_LIST <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/all/all_motifs.txt"
motifs <- read.csv(MOTIFS_LIST, header = FALSE)
motifs <- data.table(motifs)
setnames(motifs,c("Motif"))

# library(data.table)
# files <- list.files(ROOT_DIR)
# files <- data.table(files)
# setnames(files,c("sample"))
# samples <- files[grepl("_MUT_",sample) == TRUE]
# samples[, sample := tstrsplit(sample,"_MUT_")[2]]

results <- data.table()

coefplots <- list()

ci <- 1

samplesMutStat <- read.csv(paste0(STAT_DIR,"Fig5samplesStat.csv"))
samplesMutStat <- data.table(samplesMutStat)
filteredSamples <- samplesMutStat[APOBECcnt >= 3000]

for(c in 1:nrow(cancers)){
cancer <- cancers[c,CANCER_TYPE]
samples <- samplesEnrichment[CANCER_TYPE == cancer]
samples <- samples[order(APOBEC_ENRICHMENT)]

plots <- list()
dtl <- list()  

samples <- samples[SAMPLE %in% filteredSamples[,Sample]]

for(i in 1:nrow(samples)){
  sample <- samples[i, SAMPLE]
  aenrich <- samples[i, APOBEC_ENRICHMENT]
  aenrich <- round(aenrich,2)
  muts <- read.csv(paste0(ROOT_DIR,"RTEXP_MUT_",sample,".txt"), sep='\t', header=FALSE)
  muts <- data.table(muts)
  setnames(muts,c("Motif","ReplicationBin","ReplicationStrand","ExpressionBin","ExpressionStrand","MutateAllele","Cnt"))
  trgs <- read.csv(paste0(ROOT_DIR,"RTEXP_TRG_",sample,".txt"), sep='\t', header=FALSE)
  trgs <- data.table(trgs)
  setnames(trgs,c("Motif","ReplicationBin","ReplicationStrand","ExpressionBin","ExpressionStrand","TrgCnt"))
  
  muts <- muts[ReplicationBin != -2]
  trgs <- trgs[ReplicationBin != -2]
  muts <- muts[ReplicationBin != -1]
  trgs <- trgs[ReplicationBin != -1]
  
  muts <- muts[ReplicationStrand != -1]
  trgs <- trgs[ReplicationStrand != -1]
  muts[, LeadingStrand := 0]
  muts[, LaggingStrand := 0]
  muts[ReplicationStrand == 1, LeadingStrand := 1]
  muts[ReplicationStrand == 0, LaggingStrand := 1]
  trgs[, LeadingStrand := 0]
  trgs[, LaggingStrand := 0]
  trgs[ReplicationStrand == 1, LeadingStrand := 1]
  trgs[ReplicationStrand == 0, LaggingStrand := 1]

  muts <- muts[ExpressionBin != -1]
  trgs <- trgs[ExpressionBin != -1]
  muts <- muts[ExpressionBin != 0]
  trgs <- trgs[ExpressionBin != 0]
  muts[ExpressionBin == -2, ExpressionBin := 0]
  trgs[ExpressionBin == -2, ExpressionBin := 0]
  
  muts[, CodingStrand := 0]
  muts[, TemplateStrand := 0]
  muts[ExpressionStrand == 1, CodingStrand := 1]
  muts[ExpressionStrand == 0, TemplateStrand := 1]
  trgs[, CodingStrand := 0]
  trgs[, TemplateStrand := 0]
  trgs[ExpressionStrand == 1, CodingStrand := 1]
  trgs[ExpressionStrand == 0, TemplateStrand := 1]
    
  muts[,ReplicationStrand := NULL]
  muts[,ExpressionStrand := NULL]
  trgs[,ReplicationStrand := NULL]
  trgs[,ExpressionStrand := NULL]
  
  mutsAPOBEC <- muts[Motif %in% c("TCA","TCT","TCG","TCC") & MutateAllele %in% c("T","G")]
  mutsXCX <- muts[Motif %in% c("ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT")]
  mutsAPOBECrtexp <- mutsAPOBEC[, .("Cnt"=sum(Cnt)), by=.(ReplicationBin, ExpressionBin)]
  mutsXCXrtexp <- mutsXCX[, .("Cnt"=sum(Cnt)), by=.(ReplicationBin,ExpressionBin)]
  
  trgsAPOBEC <- trgs[Motif %in% c("TCA","TCT","TCG","TCC")]
  trgsXCX <- trgs[Motif %in% c("ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT")]
  trgsAPOBECrtexp <- trgsAPOBEC[, .("TrgCnt"=sum(TrgCnt)), by=.(ReplicationBin,ExpressionBin)]
  trgsXCXrtexp <- trgsXCX[, .("TrgCnt"=sum(TrgCnt)), by=.(ReplicationBin,ExpressionBin)]
  
  dataAPOBEC <- merge(trgsAPOBECrtexp,mutsAPOBECrtexp,by=c("ReplicationBin","ExpressionBin"),all.x = TRUE)
  dataAPOBEC[is.na(Cnt), Cnt:=0]
  dataAPOBEC[TrgCnt != 0, Density := Cnt/TrgCnt]
  dataAPOBEC[is.na(Density), Density := 0]
  DensSum <- dataAPOBEC[, sum(Density)]
  dataAPOBEC[, NormDensityAPOBEC := Density/DensSum]
  dataAPOBEC <- dataAPOBEC[ExpressionBin != 0]
  dataAPOBEC[, ReplicationBin := ReplicationBin + 1]
  
  lmDataAPOBEC <- dataAPOBEC[, .(NormDensityAPOBEC,ReplicationBin,ExpressionBin)]
  modelAPOBEC <- rlm(NormDensityAPOBEC ~ ., data=lmDataAPOBEC)

  dataXCX <- merge(trgsXCXrtexp,mutsXCXrtexp,by=c("ReplicationBin","ExpressionBin"),all.x = TRUE)
  dataXCX[is.na(Cnt), Cnt:=0]
  dataXCX[TrgCnt != 0, Density := Cnt/TrgCnt]
  dataXCX[is.na(Density), Density := 0]
  DensSum <- dataXCX[, sum(Density)]
  dataXCX[, NormDensityXCX := Density/DensSum]
  dataXCX <- dataXCX[ExpressionBin != 0]
  dataXCX[, ReplicationBin := ReplicationBin + 1]
  
  lmDataXCX <- dataXCX[, .(NormDensityXCX,ReplicationBin,ExpressionBin)]
  modelXCX <- rlm(NormDensityXCX ~ ., data=lmDataXCX)
  
  APOBECXCXrtexp <- merge(dataAPOBEC,dataXCX,by=c("ReplicationBin","ExpressionBin"))
  APOBECXCXrtexp[, diff := (NormDensityAPOBEC-NormDensityXCX)]
  
  lmDataAPOBECXCX <- APOBECXCXrtexp[, .(diff,ReplicationBin,ExpressionBin)]
  modelAPOBECXCX <- rlm(diff ~ ., data=lmDataAPOBECXCX)
  
  dtl[[(i-1)*3+1]] <- copy(dataAPOBEC)
  plots[[(i-1)*3+1]] <- local ({
    i <- i
    p1 <- ggplot(dtl[[(i-1)*3+1]], aes(x=factor(dtl[[(i-1)*3+1]]$ReplicationBin), 
                         y=factor(dtl[[(i-1)*3+1]]$ExpressionBin))) +
    geom_tile(aes(fill = dtl[[(i-1)*3+1]]$NormDensityAPOBEC), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    xlab("Late <- Replication timing -> Early") + ylab("Low <- Expression -> High")  +
    ggtitle(paste(cancer,", ",sample,", APOBEC_ENRICHMENT=",aenrich,sep="")) +
    geom_text(aes(x=dtl[[(i-1)*3+1]]$ReplicationBin, 
                  y=(dtl[[(i-1)*3+1]]$ExpressionBin), 
                  label=paste(round(dtl[[(i-1)*3+1]]$NormDensityAPOBEC,digits=3),"\n",dtl[[(i-1)*3+1]]$Cnt,"\n",dtl[[(i-1)*3+1]]$TrgCnt,sep="")), size = 1) +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(size=5),
          axis.text.y = element_text(size=5),
          axis.title.x = element_text(size=5),
          axis.title.y = element_text(size=5),
          legend.text = element_text(size = 5),
          plot.title = element_text(size = 5),
          legend.key.size = unit(2,"mm"))
    return(p1)
    })
  
    dtl[[(i-1)*3+2]] <- copy(dataXCX)
    plots[[(i-1)*3+2]] <- local ({
      i <- i
      p2 <- ggplot(dtl[[(i-1)*3+2]], aes(x=factor(dtl[[(i-1)*3+2]]$ReplicationBin), 
                      y=factor(dtl[[(i-1)*3+2]]$ExpressionBin))) +
    geom_tile(aes(fill = dtl[[(i-1)*3+2]]$NormDensityXCX), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    xlab("Late <- Replication timing -> Early") + ylab("Low <- Expression -> High")  +
    ggtitle(paste(cancer,", ",sample,", APOBEC_ENRICHMENT=",aenrich,sep="")) +
    geom_text(aes(x=dtl[[(i-1)*3+2]]$ReplicationBin, 
                  y=(dtl[[(i-1)*3+2]]$ExpressionBin), 
                  label=paste(format(dtl[[(i-1)*3+2]]$NormDensityXCX,digits=3),"\n",dtl[[(i-1)*3+2]]$Cnt,"\n",dtl[[(i-1)*3+2]]$TrgCnt,sep="")), size = 1)+
    theme(legend.title = element_blank(),
          axis.text.x = element_text(size=5),
          axis.text.y = element_text(size=5),
          axis.title.x = element_text(size=5),
          axis.title.y = element_text(size=5),
          legend.text = element_text(size = 5),
          legend.key.size = unit(2,"mm"),
          plot.title = element_text(size = 5))
    return(p2)
    })
    
    dtl[[(i-1)*3+3]] <- copy(APOBECXCXrtexp)
    plots[[(i-1)*3+3]] <- local({
    i <- i
    p3 <- ggplot(dtl[[(i-1)*3+3]], aes(x=factor(dtl[[(i-1)*3+3]]$ReplicationBin), 
                      y=factor(dtl[[(i-1)*3+3]]$ExpressionBin))) +
    geom_tile(aes(fill = dtl[[(i-1)*3+3]]$diff), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    xlab("Late <- Replication timing -> Early") + ylab("Low <- Expression -> High")  +
    ggtitle(paste(cancer,", ",sample,", APOBEC_ENRICHMENT=",aenrich,sep="")) +
    geom_text(aes(x=dtl[[(i-1)*3+3]]$ReplicationBin, 
                  y=(dtl[[(i-1)*3+3]]$ExpressionBin), 
                  label=round(dtl[[(i-1)*3+3]]$diff,3)), size = 1) +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(size=5),
          axis.text.y = element_text(size=5),
          axis.title.x = element_text(size=5),
          axis.title.y = element_text(size=5),
          legend.text = element_text(size = 5),
          legend.key.size = unit(2,"mm"),
          plot.title = element_text(size = 5))
    return(p3)
    })

  results <- rbind(results, data.table("Cancer"=cancer,
                                                   "Sample"=sample, 
                                                   "ApobecEnrichment"=aenrich,
                                                   "coefAPOBECReplicationBin"=summary(modelAPOBEC)$coefficients[[2]],
                                                   "coefAPOBECExpressionBin"=summary(modelAPOBEC)$coefficients[[3]],
                                                   #"coefAPOBECReplicationBin"=summary(modelAPOBEC)$coefficients[2,1],
                                                   #"coefAPOBECExpressionBin"=summary(modelAPOBEC)$coefficients[3,1],
                                                   #"stderrAPOBECReplicationBin"=summary(modelAPOBEC)$coefficients[2,2],
                                                   #"stderrAPOBECExpressionBin"==summary(modelAPOBEC)$coefficients[3,2],
                                                   #"pvalueAPOBECReplicationBin"=summary(modelAPOBEC)$coefficients[2,4],
                                                   #"pvalueAPOBECExpressionBin"=summary(modelAPOBEC)$coefficients[3,4],
                                                   "coefXCXReplicationBin"=summary(modelXCX)$coefficients[[2]],
                                                   "coefXCXExpressionBin"=summary(modelXCX)$coefficients[[3]],
                                                   #"coefXCXReplicationBin"=summary(modelXCX)$coefficients[2,1],
                                                   #"coefXCXExpressionBin"=summary(modelXCX)$coefficients[3,1],
                                                   #"stderrXCXReplicationBin"=summary(modelXCX)$coefficients[2,2],
                                                   #"stderrXCXExpressionBin"=summary(modelXCX)$coefficients[3,2],
                                                   #"pvalueXCXReplicationBin"=summary(modelXCX)$coefficients[2,4],
                                                   #"pvalueXCXExpressionBin"=summary(modelXCX)$coefficients[3,4],
                                                   "coefAPOXCXReplicationBin"=summary(modelAPOBECXCX)$coefficients[[2]],
                                                   "coefAPOXCXExpressionBin"=summary(modelAPOBECXCX)$coefficients[[3]]
                                                   #"coefAPOXCXReplicationBin"=summary(modelAPOBECXCX)$coefficients[2,1],
                                                   #"coefAPOXCXExpressionBin"=summary(modelAPOBECXCX)$coefficients[3,1],
                                                   #"stderrAPOXCXReplicationBin"=summary(modelAPOBECXCX)$coefficients[2,2],
                                                   #"stderrAPOXCXExpressionBin"=summary(modelAPOBECXCX)$coefficients[3,2],
                                                   #"pvalueAPOXCXReplicationBin"=summary(modelAPOBECXCX)$coefficients[2,4],
                                                   #"pvalueAPOXCXExpressionBin"=summary(modelAPOBECXCX)$coefficients[3,4]
                                                   ))
  
}

pl <- ggarrange(plotlist=plots, ncol=3, nrow=4)
ggexport(pl, filename=paste0(OUTPUT_DIR,"rtexp_",cancer,".pdf"))

# Coefficients for relative rep/exp surface

resCancer <- results[Cancer == cancer]
resCancer[, SampleEnrich := paste0(Sample,"_",ApobecEnrichment)]
dt4melt <- resCancer[, .(SampleEnrich,coefAPOXCXReplicationBin,coefAPOXCXExpressionBin)]
dtmelt <- melt(dt4melt)

dtmelt$SampleEnrich <- factor(dtmelt$SampleEnrich, levels=unique(dtmelt$SampleEnrich))

pp1 <- ggplot(data=dtmelt, aes(x=SampleEnrich,y=value,fill=variable)) +
  geom_bar(stat="identity",position=position_dodge()) +
  ggtitle(cancer) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5),
        axis.text.y = element_text(size=5),
        panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 5),
        legend.key.size = unit(2,"mm"))

coefplots[[(ci-1)*2+1]] <- pp1

resCancer[, diffAPOXCXReplicationBin := coefAPOBECReplicationBin - coefXCXReplicationBin]
resCancer[, diffAPOXCXExpressionBin := coefAPOBECExpressionBin - coefXCXExpressionBin]

dt4melt <- resCancer[, .(SampleEnrich,diffAPOXCXReplicationBin,diffAPOXCXExpressionBin)]
dtmelt <- melt(dt4melt)

dtmelt$SampleEnrich <- factor(dtmelt$SampleEnrich, levels=unique(dtmelt$SampleEnrich))

pp2 <- ggplot(data=dtmelt, aes(x=SampleEnrich,y=value,fill=variable)) +
  geom_bar(stat="identity",position=position_dodge()) +
  ggtitle(cancer) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5),
        axis.text.y = element_text(size=5),
        panel.background = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 5),
        legend.key.size = unit(2,"mm"))
coefplots[[(ci-1)*2+2]] <- pp2

 ci <- ci + 1
}

pl <- ggarrange(plotlist=coefplots, ncol=1, nrow=2)
ggexport(pl, filename=paste0(OUTPUT_DIR,"rtexp_coefs.pdf"))

write.csv(results,paste0(OUTPUT_DIR,"rtexp_coefs_table.csv"))
