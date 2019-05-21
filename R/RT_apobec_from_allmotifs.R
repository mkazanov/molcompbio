library(data.table)
library(ggplot2)
library(reshape2)

OUTPUT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalR"
MOTIFS_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/MOTIFS_RT/"

cancerSamples <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/samples.txt",sep='\t')
cancerSamples <- data.table(cancerSamples)

cellCancer <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR/refCellTypesCancers.txt", sep = '\t')
cellCancer <- data.table(cellCancer)

gordeninEnrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results/sample_enrichment.txt", sep = '\t')
gordeninEnrichment <- data.table(gordeninEnrichment)

cancers <- data.table("cancer"=unique(cancerSamples[,cancer]))

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

res <- data.frame()

for(i in 1:nrow(cancers))
{
  samples <- data.table("sample"=unique(cancerSamples[cancer == cancers[i]$cancer, sample]))
  cellLine <- cellCancer[Cancer == cancers[i]$cancer, CellLine]
  
  for (j in 1:nrow(samples))
  {
    enrichGordenin <- gordeninEnrichment[CANCER_TYPE == cancers[i]$cancer & SAMPLE == samples[j]$sample]  
    enrichChar <- paste0(gsub("\\.","_",round(enrichGordenin$APOBEC_ENRICHMENT,2)),"_")

  coefs <- data.frame()  
    
  for(m in c("TCT","TCA","TCC","TCG"))
  {
    tcxt <- read.csv(paste0(MOTIFS_DIR,m,"_T.txt"),sep='\t')
    tcxt <- data.table(tcxt)
    tcxt <- tcxt[Sample == samples[j]$sample & Cancer == cancers[i]$cancer]

    tcxg <- read.csv(paste0(MOTIFS_DIR,m,"_G.txt"),sep='\t')
    tcxg <- data.table(tcxg)
    tcxg <- tcxg[Sample == samples[j]$sample & Cancer == cancers[i]$cancer]
    
    tcx3 <- read.csv(paste0(MOTIFS_DIR,m,"_3.txt"), sep='\t')
    tcx3 <- data.table(tcx3)
    tcx3 <- tcx3[Sample == samples[j]$sample & Cancer == cancers[i]$cancer]
    
    tmp <- rbind(tcxt,tcxg)
    
    dt <- tmp[,.("MutationCnt"=sum(MutationCnt)),by=.(Cancer,Sample,ReplicationBin)]
    dt <- dt[ReplicationBin %in% c(0,1,2,3,4,5,6)]
    tcxNum <- read.csv(paste0(MOTIFS_DIR,m,"_in_RTbins_",cellLine,".txt"),sep="\t")
    tcxNum <- data.table(tcxNum)
    tcxNum <- tcxNum[ReplicationBin %in% c(0,1,2,3,4,5,6)]
    dt1 <- merge(dt, tcxNum, by = "ReplicationBin")
    dt1 <- dt1[, MutationDensity := MutationCnt/TargetCnt]
    dt1[, DensityTotal := sum(MutationDensity)]
    dt1[, NormalizedDensity := MutationDensity/DensityTotal]
    
    ret <- lm(NormalizedDensity ~ ReplicationBin,dt1)
    
    dt2 <- tcx3[ReplicationBin %in% c(0,1,2,3,4,5,6)]
    dt2 <- merge(dt2, tcxNum, by = "ReplicationBin")
    dt2 <- dt2[, MutationDensity := MutationCnt/TargetCnt]
    dt2[, DensityTotal := sum(MutationDensity)]
    dt2[, NormalizedDensity := MutationDensity/DensityTotal]
    
    ret2 <- lm(NormalizedDensity ~ ReplicationBin,dt2)
    
   # p <- ggplot(dt1, aes(x=ReplicationBin, y=NormalizedDensity)) + 
    #  geom_bar(stat="identity") +
    #  geom_smooth(method='lm',formula=y~x) +
    #  ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
    #                 ", Enrich: ", enrichChar, ", Coef = ",
    #                 formatC(ret$coefficients[[2]], format = "e", digits = 2))) +
    #  theme(plot.title = element_text(size = 10))

     # ggsave(paste0(OUTPUT_DIR,"/RT/",cancers[i]$cancer, "/APOBEC/", enrichChar, samples[j]$sample, "_",m,".jpg"), plot = p, device = "jpeg")
    
      mutcnt <- dt1[, sum(MutationCnt)]
      coefs <- rbind(coefs, data.frame("coef"=ret$coefficients[[2]], "mutcnt"=mutcnt, "coef3"=ret2$coefficients[[2]]))  
    }
 
  
  if(samples[j]$sample == "TCGA-05-4396-01A"){
    stop("ququ")
  }
  
   coefsAll <- data.frame()  
     
   for(m in c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","AGC","AGG","AGT",
              "ATA","ATC","ATG","CAA","CAG","CCA","CGG","CTA","GAA","GAC",
              "GAG","GCA","GCG","GGC","GGG","GTA","GTG","TAA")) 
   {
     dt <- read.csv(paste0(MOTIFS_DIR,m,"_3.txt"),sep='\t')
     dt <- data.table(dt)
     dt <- dt[Sample == samples[j]$sample & Cancer == cancers[i]$cancer]
     dt <- dt[ReplicationBin %in% c(0,1,2,3,4,5,6)]
     trgNum <- read.csv(paste0(MOTIFS_DIR,m,"_in_RTbins_",cellLine,".txt"),sep="\t")
     trgNum <- data.table(trgNum)
     trgNum <- trgNum[ReplicationBin %in% c(0,1,2,3,4,5,6)]
     dt1 <- merge(dt, trgNum, by = "ReplicationBin")
     dt1 <- dt1[, MutationDensity := MutationCnt/TargetCnt]
     dt1[, DensityTotal := sum(MutationDensity)]
     dt1[, NormalizedDensity := MutationDensity/DensityTotal]
     
     ret <- lm(NormalizedDensity ~ ReplicationBin,dt1)

     #p <- ggplot(dt1, aes(x=ReplicationBin, y=NormalizedDensity)) + 
     # geom_bar(stat="identity") +
     #  geom_smooth(method='lm',formula=y~x) +
     # ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
      #                ", Enrich: ", enrichChar, ", Coef = ",ret$coefficients[[2]]))
     
     #ggsave(paste0(OUTPUT_DIR,"/RT/",cancers[i]$cancer, "/OTHER/", enrichChar, samples[j]$sample, "_",m,".jpg"), plot = p, device = "jpeg")
     coefsAll <- rbind(coefsAll, data.frame("coef"=ret$coefficients[[2]], "mutcnt"=mutcnt,"motif"=m))  
     
   }
  
  
    coefs <- data.table(coefs)
    mutTotal <- coefs[1]$mutcnt + coefs[2]$mutcnt + coefs[3]$mutcnt + coefs[4]$mutcnt
    coefWeighted <- coefs[1]$coef*(coefs[1]$mutcnt/mutTotal) + 
                    coefs[2]$coef*(coefs[2]$mutcnt/mutTotal) +
                    coefs[3]$coef*(coefs[3]$mutcnt/mutTotal) + 
                    coefs[4]$coef*(coefs[4]$mutcnt/mutTotal)
    
    coefsAll <- data.table(coefsAll)
    coefOther <- coefsAll[,mean(coef)]
    
    res <- rbind(res, data.frame("Cancer"=cancers[i]$cancer,
                                 "Sample"=samples[j]$sample,
                                 "ApobecSlope"=coefWeighted,
                                 "coefTCT"=coefs[1]$coef,
                                 "coef3TCT"=coefs[1]$coef3,
                                 "mutcntTCT"=coefs[1]$mutcnt,
                                 "coefTCA"=coefs[2]$coef,
                                 "coef3TCA"=coefs[2]$coef3,
                                 "mutcntTCA"=coefs[2]$mutcnt,
                                 "coefTCC"=coefs[3]$coef,
                                 "coef3TCC"=coefs[3]$coef3,
                                 "mutcntTCC"=coefs[3]$mutcnt,
                                 "coefTCG"=coefs[4]$coef,
                                 "coef3TCG"=coefs[4]$coef3,
                                 "mutcntTCG"=coefs[4]$mutcnt,
                                 "OtherSlope"=coefOther,
                                 "GordeninEnrichment"=enrichGordenin$APOBEC_ENRICHMENT))  
    
  }
}

write.csv(res,paste0(OUTPUT_DIR,"/RT/coefsFinal.csv"))

