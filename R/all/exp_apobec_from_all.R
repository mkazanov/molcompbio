library(data.table)
library(ggplot2)
library(reshape2)

OUTPUT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalRall"
INPUT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL/"

cancerSamples <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/samples.txt",sep='\t')
cancerSamples <- data.table(cancerSamples)

cellCancer <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR/refCellTypesCancers.txt", sep = '\t')
cellCancer <- data.table(cellCancer)

gordeninEnrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results/sample_enrichment.txt", sep = '\t')
gordeninEnrichment <- data.table(gordeninEnrichment)

cancers <- data.table("cancer"=unique(cancerSamples[,cancer]))

# Create subdirs for cancer types

dir.create(file.path(OUTPUT_DIR,'EXP'))

for(i in 1:nrow(cancers))
{
  
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXP'), cancers[i]$cancer))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXP/',cancers[i]$cancer), "APOBEC"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXP/',cancers[i]$cancer), "OTHER"))

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
    
  dataMut <- read.csv(paste0(INPUT_DIR,"RTEXP_MUT_",samples[j]$sample,".txt"), sep='\t', header=FALSE)
  dataMut <- data.table(dataMut)
  setnames(dataMut,c("Motif","RTbin","RTstrand","expbin","senseStrand","mutAllele","MutationCnt"))
  
  dataTrg <- read.csv(paste0(INPUT_DIR,"RTEXP_TRG_",samples[j]$sample,".txt"),sep="\t", header=FALSE)
  dataTrg <- data.table(dataTrg)
  setnames(dataTrg,c("Motif","RTbin","RTstrand","expbin","senseStrand","TargetCnt"))
  
  for(m in c("TCT","TCA","TCC","TCG"))
  {
    tcx <- dataMut[Motif == m & mutAllele %in% c("T","G")]
    dt <- tcx[,.("MutationCnt"=sum(MutationCnt)),by=.(expbin)]
    dt <- dt[expbin %in% c(1,2,3,4,5,6,7)]
    
    tcxNum <- dataTrg[Motif == m]
    tcxNum <- tcxNum[expbin %in% c(1,2,3,4,5,6,7)]
    tcxNum <- tcxNum[,.("TargetCnt"=sum(TargetCnt)), by=.(expbin)]
    dt1 <- merge(dt, tcxNum, by = "expbin", all.y=TRUE)
    dt1[is.na(MutationCnt), MutationCnt:=0]
    dt1 <- dt1[, MutationDensity := MutationCnt/TargetCnt]
    dt1[, DensityTotal := sum(MutationDensity)]
    dt1[, NormalizedDensity := ifelse(DensityTotal==0,0,MutationDensity/DensityTotal)]
    
    ret <- lm(NormalizedDensity ~ expbin,dt1)
    
   # p <- ggplot(dt1, aes(x=ReplicationBin, y=NormalizedDensity)) + 
    #  geom_bar(stat="identity") +
    #  geom_smooth(method='lm',formula=y~x) +
    #  ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
    #                 ", Enrich: ", enrichChar, ", Coef = ",
    #                 formatC(ret$coefficients[[2]], format = "e", digits = 2))) +
    #  theme(plot.title = element_text(size = 10))

     # ggsave(paste0(OUTPUT_DIR,"/RT/",cancers[i]$cancer, "/APOBEC/", enrichChar, samples[j]$sample, "_",m,".jpg"), plot = p, device = "jpeg")
    
      mutcnt <- dt1[, sum(MutationCnt)]
      coefs <- rbind(coefs, data.frame("coef"=ret$coefficients[[2]], 
                                       "stderr"=summary(ret)$coefficients[2,2],
                                       "pvalue"=summary(ret)$coefficients[2,4],
                                       "mutcnt"=mutcnt))  
    }
  
   coefsAll <- data.frame()  
     
   # All other triplets with middle C
   C2motifList <- c("ACA","ACC","ACG","GCT","CCT","ACT","CCA","CCG","GCA","GCG","GCC","CCC")
   for(m in C2motifList) 
   {
     xt <- dataMut[Motif == m & mutAllele %in% c("T","G")]
     dt <- xt[,.("MutationCnt"=sum(MutationCnt)),by=.(expbin)]
     dt <- dt[expbin %in% c(1,2,3,4,5,6,7)]
     
     trgNum <- dataTrg[Motif == m]     
     trgNum <- trgNum[expbin %in% c(1,2,3,4,5,6,7)]
     trgNum <- trgNum[,.("TargetCnt"=sum(TargetCnt)), by=.(expbin)]
     dt1 <- merge(dt, trgNum, by = "expbin",all.y=TRUE)
     dt1[is.na(MutationCnt), MutationCnt:=0]
     dt1 <- dt1[, MutationDensity := MutationCnt/TargetCnt]
     dt1[, DensityTotal := sum(MutationDensity)]
     dt1[, NormalizedDensity := ifelse(DensityTotal==0,0,MutationDensity/DensityTotal)]
     
     ret <- lm(NormalizedDensity ~ expbin,dt1)

     #p <- ggplot(dt1, aes(x=ReplicationBin, y=NormalizedDensity)) + 
     # geom_bar(stat="identity") +
     #  geom_smooth(method='lm',formula=y~x) +
     # ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
      #                ", Enrich: ", enrichChar, ", Coef = ",ret$coefficients[[2]]))
     
     #ggsave(paste0(OUTPUT_DIR,"/RT/",cancers[i]$cancer, "/OTHER/", enrichChar, samples[j]$sample, "_",m,".jpg"), plot = p, device = "jpeg")
     mutcnt <- dt1[, sum(MutationCnt)]
     coefsAll <- rbind(coefsAll, data.frame("coef"=ret$coefficients[[2]], 
                                            "stderr"=summary(ret)$coefficients[2,2],
                                            "pvalue"=summary(ret)$coefficients[2,4],
                                            "mutcnt"=mutcnt,
                                            "motif"=m))  
     
   }
  
    coefs <- data.table(coefs)
    mutTotal <- coefs[1]$mutcnt + coefs[2]$mutcnt + coefs[3]$mutcnt + coefs[4]$mutcnt
    coefWeighted <- coefs[1]$coef*(coefs[1]$mutcnt/mutTotal) + 
                    coefs[2]$coef*(coefs[2]$mutcnt/mutTotal) +
                    coefs[3]$coef*(coefs[3]$mutcnt/mutTotal) + 
                    coefs[4]$coef*(coefs[4]$mutcnt/mutTotal)
    
    # All 4 APOBEC motifs
    
    tcx <- dataMut[Motif %in% c("TCT","TCA","TCC","TCG") & mutAllele %in% c("T","G")]
    dt <- tcx[,.("MutationCnt"=sum(MutationCnt)),by=.(expbin)]
    dt <- dt[expbin %in% c(1,2,3,4,5,6,7)]
    
    trgNum <- dataTrg[Motif %in% c("TCT","TCA","TCC","TCG")]     
    trgNum <- trgNum[expbin %in% c(1,2,3,4,5,6,7)]
    trgNum <- trgNum[,.("TargetCnt"=sum(TargetCnt)), by=.(expbin)]
    dt1 <- merge(dt, trgNum, by = "expbin",all.y=TRUE)
    dt1[is.na(MutationCnt), MutationCnt:=0]
    dt1 <- dt1[, MutationDensity := MutationCnt/TargetCnt]
    dt1[, DensityTotal := sum(MutationDensity)]
    dt1[, NormalizedDensity := MutationDensity/DensityTotal]
    
    ret <- lm(NormalizedDensity ~ expbin,dt1)
    
    # All C2-based motifs together
    
    xt <- dataMut[Motif %in% C2motifList & mutAllele %in% c("T","G")]
    dt <- xt[,.("MutationCnt"=sum(MutationCnt)),by=.(expbin)]
    dt <- dt[expbin %in% c(1,2,3,4,5,6,7)]
    
    trgNum <- dataTrg[Motif %in% C2motifList]     
    trgNum <- trgNum[expbin %in% c(1,2,3,4,5,6,7)]
    trgNum <- trgNum[,.("TargetCnt"=sum(TargetCnt)), by=.(expbin)]
    dt1 <- merge(dt, trgNum, by = "expbin",all.y=TRUE)
    dt1[is.na(MutationCnt), MutationCnt:=0]
    dt1 <- dt1[, MutationDensity := MutationCnt/TargetCnt]
    dt1[, DensityTotal := sum(MutationDensity)]
    dt1[, NormalizedDensity := MutationDensity/DensityTotal]
    
    с2oth <- lm(NormalizedDensity ~ expbin,dt1)
    
    # All other motifs
    notAPOBECmotifList <- c("AAA","AAC","AAG","AAT","ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT",
                            "CTA","CTC","CTG","CTT","GAA","GAC","GTA","TAA",
                            "ACA","ACC","ACG","GCT","CCT","ACT","CCA","CCG","GCA","GCG","GCC","CCC")
    
    xt <- dataMut[Motif %in% notAPOBECmotifList]
    dt <- xt[,.("MutationCnt"=sum(MutationCnt)),by=.(expbin)]
    dt <- dt[expbin %in% c(1,2,3,4,5,6,7)]
    
    trgNum <- dataTrg[Motif %in% notAPOBECmotifList]     
    trgNum <- trgNum[expbin %in% c(1,2,3,4,5,6,7)]
    trgNum <- trgNum[,.("TargetCnt"=sum(TargetCnt)), by=.(expbin)]
    dt1 <- merge(dt, trgNum, by = "expbin",all.y=TRUE)
    dt1[is.na(MutationCnt), MutationCnt:=0]
    dt1 <- dt1[, MutationDensity := MutationCnt/TargetCnt]
    dt1[, DensityTotal := sum(MutationDensity)]
    dt1[, NormalizedDensity := MutationDensity/DensityTotal] 
    
    oth <- lm(NormalizedDensity ~ expbin,dt1)
    
    coefsAll <- data.table(coefsAll)
    coefOther <- coefsAll[,mean(coef)]
    
    res <- rbind(res, data.frame("Cancer"=cancers[i]$cancer,
                                 "Sample"=samples[j]$sample,
                                 "coefWeighted"=coefWeighted,
                                 "coef"=ret$coefficients[[2]],
                                 "stderr"=summary(ret)$coefficients[2,2],
                                 "pvalue"=summary(ret)$coefficients[2,4],
                                 "coefTCT"=coefs[1]$coef,
                                 "stderrTCT"=coefs[1]$stderr,
                                 "pvalueTCT"=coefs[1]$pvalue,
                                 "mutcntTCT"=coefs[1]$mutcnt,
                                 "coefTCA"=coefs[2]$coef,
                                 "stderrTCA"=coefs[2]$stderr,
                                 "pvalueTCA"=coefs[2]$pvalue,
                                 "mutcntTCA"=coefs[2]$mutcnt,
                                 "coefTCC"=coefs[3]$coef,
                                 "stderrTCC"=coefs[3]$stderr,
                                 "pvalueTCC"=coefs[3]$pvalue,
                                 "mutcntTCC"=coefs[3]$mutcnt,
                                 "coefTCG"=coefs[4]$coef,
                                 "stderrTCG"=coefs[4]$stderr,
                                 "pvalueTCG"=coefs[4]$pvalue,
                                 "mutcntTCG"=coefs[4]$mutcnt,
                                 "meanС2OtherSlope"=coefOther,
                                 "coefC2Other"=с2oth$coefficients[[2]],
                                 "stderrC2Other"=summary(с2oth)$coefficients[2,2],
                                 "pvalueC2Other"=summary(с2oth)$coefficients[2,4],
                                 "coefOther"=oth$coefficients[[2]],
                                 "stderrOther"=summary(oth)$coefficients[2,2],
                                 "pvalueOther"=summary(oth)$coefficients[2,4],
                                 "GordeninEnrichment"=enrichGordenin$APOBEC_ENRICHMENT))  
    
  }
}

write.csv(res,paste0(OUTPUT_DIR,"/EXP/coefsFinal.csv"))

