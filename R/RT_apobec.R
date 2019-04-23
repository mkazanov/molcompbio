library(data.table)
library(ggplot2)
library(reshape2)

### Input data

OUTPUT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR"
apobec <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/results_RT_APOBEC.txt", sep = '\t',header = TRUE)
apobec <- data.table(apobec)
other <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/results_RT_OTHER.txt", sep = '\t',header = TRUE)
other <- data.table(other)

# TCW targets

tcwIMR90 <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/TCW_in_RTbins_IMR90.txt", sep = '\t')
tcwIMR90 <- data.table(tcwIMR90)
tcwMCF7 <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/TCW_in_RTbins_MCF7.txt", sep = '\t')
tcwMCF7 <- data.table(tcwMCF7)
tcwNHEK <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/TCW_in_RTbins_NHEK.txt", sep = '\t')
tcwNHEK <- data.table(tcwNHEK)

tcwIMR90[, CellType := 'IMR90']
tcwMCF7[, CellType := 'MCF7']
tcwNHEK[, CellType := 'NHEK']

tcwNum <- rbind(tcwIMR90,tcwMCF7,tcwNHEK)

# All targets

allIMR90 <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL_in_RTbins_IMR90.txt", sep = '\t')
allIMR90 <- data.table(allIMR90)
allMCF7 <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL_in_RTbins_MCF7.txt", sep = '\t')
allMCF7 <- data.table(allMCF7)
allNHEK <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL_in_RTbins_NHEK.txt", sep = '\t')
allNHEK <- data.table(allNHEK)

allIMR90[, CellType := 'IMR90']
allMCF7[, CellType := 'MCF7']
allNHEK[, CellType := 'NHEK']

allNum <- rbind(allIMR90,allMCF7,allNHEK)


cellCancer <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR/refCellTypesCancers.txt", sep = '\t')
cellCancer <- data.table(cellCancer)

enrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/enrichment.txt", sep = '\t')
enrichment <- data.table(enrichment)

gordeninEnrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results/sample_enrichment.txt", sep = '\t')
gordeninEnrichment <- data.table(gordeninEnrichment)

# Processing data

cancers <- data.table("cancer"=unique(apobec[,Cancer]))

res <- data.frame()

# Create subdirs for cancer types

dir.create(file.path(OUTPUT_DIR,'RT'))
dir.create(file.path(OUTPUT_DIR,'RTstrand'))
dir.create(file.path(paste0(OUTPUT_DIR,'/RTstrand'),'RAW'))
dir.create(file.path(paste0(OUTPUT_DIR,'/RTstrand'),'RATIO'))
                      
for(i in 1:nrow(cancers))
{
  
  dir.create(file.path(paste0(OUTPUT_DIR,'/RT'), cancers[i]$cancer))
  dir.create(file.path(paste0(OUTPUT_DIR,'/RT/',cancers[i]$cancer), "APOBEC"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/RT/',cancers[i]$cancer), "OTHER"))
  
  dir.create(file.path(paste0(OUTPUT_DIR,'/RTstrand/RAW'), cancers[i]$cancer))
  dir.create(file.path(paste0(OUTPUT_DIR,'/RTstrand/RATIO'), cancers[i]$cancer))
  dir.create(file.path(paste0(OUTPUT_DIR,'/RTstrand/RAW/',cancers[i]$cancer), "APOBEC"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/RTstrand/RAW/',cancers[i]$cancer), "OTHER"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/RTstrand/RATIO/',cancers[i]$cancer), "APOBEC"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/RTstrand/RATIO/',cancers[i]$cancer), "OTHER"))
}
 
# Replication domains
 
for(i in 1:nrow(cancers))
{
 samples <- data.table("sample"=unique(apobec[Cancer == cancers[i]$cancer, Sample]))
 cellLine <- cellCancer[Cancer == cancers[i]$cancer, CellType]
 tcwNumCancer <- tcwNum[CellType == cellLine] 
 allNumCancer <- allNum[CellType == cellLine]
 
 for (j in 1:nrow(samples))
 {
  enrich <- enrichment[cancer == cancers[i]$cancer & sample == samples[j]$sample]
  enrichChar <- paste0(gsub("\\.","_",round(enrich$Enrichment_exclude_TCW,2)),"_")
  enrichGordenin <- gordeninEnrichment[CANCER_TYPE == cancers[i]$cancer & SAMPLE == samples[j]$sample]  
  
  # APOBEC 
  dt <- apobec[Sample == samples[j]$sample & Cancer == cancers[i]$cancer]
  dt <- dt[ReplicationBin %in% c(0,1,2,3,4,5,6)]
  dt1 <- merge(dt, tcwNumCancer, by = "ReplicationBin")
  dt1 <- dt1[, MutationDensity := MutationCnt/TargetCnt]
  dt1[, DensityTotal := sum(MutationDensity)]
  dt1[, NormalizedDensity := MutationDensity/DensityTotal]
  
  retAPOBEC <- lm(NormalizedDensity ~ ReplicationBin,dt1)
  
  p <- ggplot(dt1, aes(x=ReplicationBin, y=NormalizedDensity)) + 
    geom_bar(stat="identity") +
    geom_smooth(method='lm',formula=y~x) +
    ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
                   ", Enrich: ", round(enrich$Enrichment_exclude_TCW,2), ", Coef = ",
                   formatC(retAPOBEC$coefficients[[2]], format = "e", digits = 2))) +
    theme(plot.title = element_text(size = 10))
  
  ggsave(paste0(OUTPUT_DIR,"/RT/",cancers[i]$cancer, "/APOBEC/", enrichChar, samples[j]$sample, ".jpg"), plot = p, device = "jpeg")

  # Strand
  dt1 <- dt1[, LeadingDensity := LeadingCnt/TargetCnt]
  dt1 <- dt1[, LaggingDensity := LaggingCnt/TargetCnt]
  dt1 <- dt1[, StrandRatio := LaggingCnt/LeadingCnt]
  dt2 <- dt1[, .(ReplicationBin,LaggingDensity,LeadingDensity)]
  dt2 <- melt(dt2, id.vars = "ReplicationBin")
  
  retLagging <- lm(LaggingDensity ~ ReplicationBin,dt1)
  retLeading <- lm(LeadingDensity ~ ReplicationBin,dt1)
  
  p <- ggplot(dt2, aes(x=ReplicationBin, y=value, fill=variable)) + 
    geom_bar(stat="identity", position="dodge") +
    geom_smooth(method='lm') +
    ylab("Mutation density") +
    ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
                   ", Enrich: ", round(enrich$Enrichment_exclude_TCW,2),
                   ", CoefLagging = ", formatC(retLagging$coefficients[[2]], format = "e", digits = 2),
                   ", CoefLeading = ", formatC(retLeading$coefficients[[2]], format = "e", digits = 2))) +
    theme(plot.title = element_text(size = 10))
  
  ggsave(paste0(OUTPUT_DIR,"/RTstrand/RAW/",cancers[i]$cancer, "/APOBEC/", enrichChar, samples[j]$sample, ".jpg"), plot = p, device = "jpeg")
  
  retRatio <- lm(StrandRatio ~ ReplicationBin,dt1)
  
  p <- ggplot(dt1, aes(x=ReplicationBin, y=StrandRatio)) + 
    geom_bar(stat="identity") +
    geom_smooth(method='lm',formula=y~x) +
    ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
                   ", Enrich: ", round(enrich$Enrichment_exclude_TCW,2), ", Coef = ",round(retRatio$coefficients[[2]],2))) +
    theme(plot.title = element_text(size = 10))
  
  ggsave(paste0(OUTPUT_DIR,"/RTstrand/RATIO/",cancers[i]$cancer, "/APOBEC/", enrichChar, samples[j]$sample, ".jpg"), plot = p, device = "jpeg")
  
  # Other
  dt <- other[Sample == samples[j]$sample & Cancer == cancers[i]$cancer]
  dt <- dt[ReplicationBin %in% c(0,1,2,3,4,5,6)]
  dt1 <- merge(dt, allNumCancer, by = "ReplicationBin")
  dt1 <- dt1[, MutationDensity := MutationCnt/TargetCnt]
  dt1[, DensityTotal := sum(MutationDensity)]
  dt1[, NormalizedDensity := MutationDensity/DensityTotal]
  
  retOTHER <- lm(NormalizedDensity ~ ReplicationBin,dt1)
  
  p <- ggplot(dt1, aes(x=ReplicationBin, y=NormalizedDensity)) + 
    geom_bar(stat="identity") +
    geom_smooth(method='lm',formula=y~x) +
    ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
                   ", Enrich: ", enrich$Enrichment_exclude_TCW, ", Coef = ",retOTHER$coefficients[[2]]))
  
  ggsave(paste0(OUTPUT_DIR,"/RT/",cancers[i]$cancer, "/OTHER/", enrichChar, samples[j]$sample, ".jpg"), plot = p, device = "jpeg")
  
  res <- rbind(res, data.frame("Cancer"=cancers[i]$cancer,
                               "Sample"=samples[j]$sample,
                               "ApobecSlope"=retAPOBEC$coefficients[[2]],
                               "OtherSlope"=retOTHER$coefficients[[2]],
                               "LaggingSlope"=retLagging$coefficients[[2]],
                               "LeadingSlope"=retLeading$coefficients[[2]],
                               "StrandRatioSlope"=retRatio$coefficients[[2]],
                               "GordeninEnrichment"=enrichGordenin$APOBEC_ENRICHMENT,
                               "Enrichment"=enrich$Enrichment,
                               "Enrichment_exclude_TCW"=enrich$Enrichment_exclude_TCW))  

 }
}

write.csv(res,"/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR/RT/coefs.csv")

