library(data.table)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(viridis)

ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL/"
MOTIFS_LIST <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/all/all_motifs.txt"

samplesEnrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results/sample_enrichment.txt", sep='\t')
samplesEnrichment <- data.table(samplesEnrichment)
samplesEnrichment$SAMPLE <- as.character(samplesEnrichment$SAMPLE)
cancers <- unique(samplesEnrichment[, CANCER_TYPE])
cancers <- data.table(cancers)
setnames(cancers,c("CANCER_TYPE"))
cancers$CANCER_TYPE <- as.character(cancers$CANCER_TYPE)
motifs <- read.csv(MOTIFS_LIST, header = FALSE)
motifs <- data.table(motifs)
setnames(motifs,c("Motif"))

results <- data.table()

for(c in 1:nrow(cancers)){
  cancer <- cancers[c,CANCER_TYPE]
  samples <- samplesEnrichment[CANCER_TYPE == cancer]
  samples <- samples[order(APOBEC_ENRICHMENT)]
  
  for(i in 1:nrow(samples)){
    
    sample <- samples[i, SAMPLE]
    aenrich <- samples[i, APOBEC_ENRICHMENT]
    enrichChar <- paste0(gsub("\\.","_",round(aenrich,2)),"_")

    muts <- read.csv(paste0(ROOT_DIR,"RTEXP_MUT_",sample,".txt"), sep='\t', header=FALSE)
    muts <- data.table(muts)
    setnames(muts,c("Motif","ReplicationBin","ReplicationStrand","ExpressionBin","ExpressionStrand","MutateAllele","Cnt"))
    muts <- muts[Motif %in% motifs[,Motif]]
    
    trgs <- read.csv(paste0(ROOT_DIR,"RTEXP_TRG_",sample,".txt"), sep='\t', header=FALSE)
    trgs <- data.table(trgs)
    setnames(trgs,c("Motif","ReplicationBin","ReplicationStrand","ExpressionBin","ExpressionStrand","Cnt"))
    trgs <- trgs[Motif %in% motifs[,Motif]]
    
    muts[, inGenes := ifelse(ExpressionBin == -2, 0, 1)]
    mutsNotGenes <- muts[inGenes == 0, sum(Cnt)]
    mutsInGenes <- muts[inGenes == 1, sum(Cnt)]
    
    trgs[, inGenes := ifelse(ExpressionBin == -2, 0, 1)]
    trgsNotGenes <- trgs[inGenes == 0, sum(Cnt)]
    trgsInGenes <- trgs[inGenes == 1, sum(Cnt)]
    
    results <- rbind(results, data.table("Cancer"=cancer,
                                         "Sample"=sample,
                                         "MutsInGens"=mutsInGenes,
                                         "TrgsInGenes"=trgsInGenes,
                                         "MutsNotGenes"=mutsNotGenes,
                                         "TrgsNotGenes"=trgsNotGenes,
                                         "Aenrich"=aenrich))
    
    
  }
}

results[, DensityInGenes := MutsInGens/TrgsInGenes]
results[, DensityNotGenes := MutsNotGenes/TrgsNotGenes]
results[, DensityRatioLog2 := log2(DensityInGenes/DensityNotGenes)]
results[, TotalRatioLog2 := log2(MutsInGens/MutsNotGenes)]
results[, apobecMuts := MutsInGens + MutsNotGenes]

resultsTrim <- results[DensityRatioLog2 >= -1]

for(c in 1:nrow(cancers)){
  cancer <- cancers[c,CANCER_TYPE]
  dt <- resultsTrim[Cancer == cancer]

  ggscatter(dt, x = "Aenrich", y = "DensityRatioLog2", color="TotalRatioLog2",
           add = "reg.line",alpha = 0.75, size="apobecMuts", xlab = FALSE, ylab = FALSE, xticks.by = 0.5, ylim=c(-1,0.05))
    scale_size(range = c(0.5, 10))  
  
 ggsave(paste0("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/paper/pics/Fig6a/",cancer,".tiff"),units="mm",width=140,height=140,dpi=300)
}

results[, fractionlog2 := log2(MutsInGens/MutsNotGenes)]
results$Cancer <- as.factor(results$Cancer)

ggplot(data=results,aes(x=Cancer,y=fractionlog2,color=Aenrich)) +
  geom_beeswarm(size=3, cex=1.7) +
  scale_color_viridis(option = "plasma") +
  ylim(-1.75,0) +
  #stat_summary(aes(ymin = ..y.., max = ..y..), fun.y = 'mean', geom = 'crossbar', color = 'grey50', fatten=1) +
  theme(panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.position = "bottom")

ggsave("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/paper/pics/Fig6a/fig6a.tiff",units="mm",width=150,height=100,dpi=300)

