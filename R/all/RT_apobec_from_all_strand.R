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

dir.create(file.path(OUTPUT_DIR,'RTStrand'))

for(i in 1:nrow(cancers))
{
  
  dir.create(file.path(paste0(OUTPUT_DIR,'/RTStrand'), cancers[i]$cancer))
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

    results <- data.frame()  
    
    dataMut <- read.csv(paste0(INPUT_DIR,"RTEXP_MUT_",samples[j]$sample,".txt"), sep='\t', header=FALSE)
    dataMut <- data.table(dataMut)
    setnames(dataMut,c("Motif","RTbin","RTstrand","expbin","senseStrand","mutAllele","MutationCnt"))
  
    dataTrg <- read.csv(paste0(INPUT_DIR,"RTEXP_TRG_",samples[j]$sample,".txt"),sep="\t", header=FALSE)
    dataTrg <- data.table(dataTrg)
    setnames(dataTrg,c("Motif","RTbin","RTstrand","expbin","senseStrand","TargetCnt"))
  
    for(m in c("TCT","TCA","TCC","TCG"))
    {
      tcx <- dataMut[Motif == m & mutAllele %in% c("T","G")]
      tcx[,LeadingCnt:=ifelse(RTstrand==1,MutationCnt,0)]
      tcx[,LaggingCnt:=ifelse(RTstrand==0,MutationCnt,0)]
    
      dt <- tcx[,.("MutationCnt"=sum(MutationCnt),
                  "LeadingCnt"=sum(LeadingCnt),
                 "LaggingCnt"=sum(LaggingCnt)),
              by=.(RTbin)]
      dt <- dt[RTbin %in% c(0,1,2,3,4,5,6)]
     leading <- dt[,sum(LeadingCnt)]
      lagging <- dt[,sum(LaggingCnt)]
      ratio <- lagging/leading
    
     leading0 <- dt[RTbin == 0,LeadingCnt]
     if(length(leading0) == 0) {leading0 <- NA}
      lagging0 <- dt[RTbin == 0,LaggingCnt]
      if(length(lagging0) == 0) {lagging0 <- NA}
      ratio0 <- lagging0/leading0
      if(length(ratio0) == 0){ratio0 <- NA}
    
      leading1 <- dt[RTbin == 1,LeadingCnt]
      if(length(leading1) == 0) {leading1 <- NA}
      lagging1 <- dt[RTbin == 1,LaggingCnt]
      if(length(lagging1) == 0) {lagging1 <- NA}
      ratio1 <- lagging1/leading1
      if(length(ratio1) == 0){ratio1 <- NA}
    
      leading2 <- dt[RTbin == 2,LeadingCnt]
      if(length(leading2) == 0) {leading2 <- NA}
      lagging2 <- dt[RTbin == 2,LaggingCnt]
      if(length(lagging2) == 0) {lagging2 <- NA}
      ratio2 <- lagging2/leading2
      if(length(ratio2) == 0){ratio2 <- NA}
    
      leading3 <- dt[RTbin == 3,LeadingCnt]
      if(length(leading3) == 0) {leading3 <- NA}
      lagging3 <- dt[RTbin == 3,LaggingCnt]
      if(length(lagging3) == 0) {lagging3 <- NA}
      ratio3 <- lagging3/leading3
      if(length(ratio3) == 0){ratio3 <- NA}
    
      leading4 <- dt[RTbin == 4,LeadingCnt]
      if(length(leading4) == 0) {leading4 <- NA}
      lagging4 <- dt[RTbin == 4,LaggingCnt]
      if(length(lagging4) == 0) {lagging4 <- NA}
      ratio4 <- lagging4/leading4
      if(length(ratio4) == 0){ratio4 <- NA}
    
      leading5 <- dt[RTbin == 5,LeadingCnt]
      if(length(leading5) == 0) {leading5 <- NA}
      lagging5 <- dt[RTbin == 5,LaggingCnt]
      if(length(lagging5) == 0) {lagging5 <- NA}
      ratio5 <- lagging5/leading5
      if(length(ratio5) == 0){ratio5 <- NA}
    
      leading6 <- dt[RTbin == 6,LeadingCnt]
      if(length(leading6) == 0) {leading6 <- NA}
      lagging6 <- dt[RTbin == 6,LaggingCnt]
      if(length(lagging6) == 0) {lagging6 <- NA}
      ratio6 <- lagging6/leading6
      if(length(ratio6) == 0){ ratio6 <- NA}
    
   # p <- ggplot(dt1, aes(x=ReplicationBin, y=NormalizedDensity)) + 
    #  geom_bar(stat="identity") +
    #  geom_smooth(method='lm',formula=y~x) +
    #  ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
    #                 ", Enrich: ", enrichChar, ", Coef = ",
    #                 formatC(ret$coefficients[[2]], format = "e", digits = 2))) +
    #  theme(plot.title = element_text(size = 10))

     # ggsave(paste0(OUTPUT_DIR,"/RT/",cancers[i]$cancer, "/APOBEC/", enrichChar, samples[j]$sample, "_",m,".jpg"), plot = p, device = "jpeg")
    
      mutcnt <- dt[, sum(MutationCnt)]
      results <- rbind(results, data.frame("motif"=m,
                                           "ratio"=ratio, 
                                           "leading"=leading,
                                           "lagging"=lagging,
                                           "ratio0"=ratio0,
                                           "leading0"=leading0,  
                                           "lagging0"=lagging0,  
                                           "ratio1"=ratio1,
                                           "leading1"=leading1,  
                                           "lagging1"=lagging1,  
                                           "ratio2"=ratio2,
                                           "leading2"=leading2,  
                                           "lagging2"=lagging2,  
                                           "ratio3"=ratio3,
                                           "leading3"=leading3,  
                                           "lagging3"=lagging3,  
                                           "ratio4"=ratio4,
                                           "leading4"=leading4,  
                                           "lagging4"=lagging4,  
                                           "ratio5"=ratio5,
                                           "leading5"=leading5,  
                                           "lagging5"=lagging5,  
                                           "ratio6"=ratio6,
                                           "leading6"=leading6,  
                                           "lagging6"=lagging6,  
                                           "mutcnt"=mutcnt))
    }
  
   resultsAll <- data.frame()  
     
   C2motifList <- c("ACA","ACC","ACG","GCT","CCT","ACT","CCA","CCG","GCA","GCG","GCC","CCC")
   
   # All other triplets with middle C
   for(m in C2motifList) 
   {
     xt <- dataMut[Motif == m & mutAllele %in% c("T","G")]
     xt[,LeadingCnt:=ifelse(RTstrand==1,MutationCnt,0)]
     xt[,LaggingCnt:=ifelse(RTstrand==0,MutationCnt,0)]
     
     dt <- xt[,.("MutationCnt"=sum(MutationCnt),
                  "LeadingCnt"=sum(LeadingCnt),
                  "LaggingCnt"=sum(LaggingCnt)),
               by=.(RTbin)]
     dt <- dt[RTbin %in% c(0,1,2,3,4,5,6)]
     
     leading <- dt[,sum(LeadingCnt)]
     lagging <- dt[,sum(LaggingCnt)]
     ratio <- lagging/leading

     leading0 <- dt[RTbin == 0,LeadingCnt]
     if(length(leading0) == 0) {leading0 <- NA}
     lagging0 <- dt[RTbin == 0,LaggingCnt]
     if(length(lagging0) == 0) {lagging0 <- NA}
     ratio0 <- lagging0/leading0
     if(length(ratio0) == 0){ratio0 <- NA}
     
     leading1 <- dt[RTbin == 1,LeadingCnt]
     if(length(leading1) == 0) {leading1 <- NA}
     lagging1 <- dt[RTbin == 1,LaggingCnt]
     if(length(lagging1) == 0) {lagging1 <- NA}
     ratio1 <- lagging1/leading1
     if(length(ratio1) == 0){ratio1 <- NA}
     
     leading2 <- dt[RTbin == 2,LeadingCnt]
     if(length(leading2) == 0) {leading2 <- NA}
     lagging2 <- dt[RTbin == 2,LaggingCnt]
     if(length(lagging2) == 0) {lagging2 <- NA}
     ratio2 <- lagging2/leading2
     if(length(ratio2) == 0){ratio2 <- NA}
     
     leading3 <- dt[RTbin == 3,LeadingCnt]
     if(length(leading3) == 0) {leading3 <- NA}
     lagging3 <- dt[RTbin == 3,LaggingCnt]
     if(length(lagging3) == 0) {lagging3 <- NA}
     ratio3 <- lagging3/leading3
     if(length(ratio3) == 0){ratio3 <- NA}
     
     leading4 <- dt[RTbin == 4,LeadingCnt]
     if(length(leading4) == 0) {leading4 <- NA}
     lagging4 <- dt[RTbin == 4,LaggingCnt]
     if(length(lagging4) == 0) {lagging4 <- NA}
     ratio4 <- lagging4/leading4
     if(length(ratio4) == 0){ratio4 <- NA}
     
     leading5 <- dt[RTbin == 5,LeadingCnt]
     if(length(leading5) == 0) {leading5 <- NA}
     lagging5 <- dt[RTbin == 5,LaggingCnt]
     if(length(lagging5) == 0) {lagging5 <- NA}
     ratio5 <- lagging5/leading5
     if(length(ratio5) == 0){ratio5 <- NA}
     
     leading6 <- dt[RTbin == 6,LeadingCnt]
     if(length(leading6) == 0) {leading6 <- NA}
     lagging6 <- dt[RTbin == 6,LaggingCnt]
     if(length(lagging6) == 0) {lagging6 <- NA}
     ratio6 <- lagging6/leading6
     if(length(ratio6) == 0){ ratio6 <- NA}
     
     #p <- ggplot(dt1, aes(x=ReplicationBin, y=NormalizedDensity)) + 
     # geom_bar(stat="identity") +
     #  geom_smooth(method='lm',formula=y~x) +
     # ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
      #                ", Enrich: ", enrichChar, ", Coef = ",ret$coefficients[[2]]))
     
     #ggsave(paste0(OUTPUT_DIR,"/RT/",cancers[i]$cancer, "/OTHER/", enrichChar, samples[j]$sample, "_",m,".jpg"), plot = p, device = "jpeg")
     mutcnt <- dt[, sum(MutationCnt)]
     resultsAll <- rbind(resultsAll, data.frame("motif"=m,
                                                "ratio"=ratio, 
                                                "leading"=leading,
                                                "lagging"=lagging,
                                                "ratio0"=ratio0,
                                                "leading0"=leading0,  
                                                "lagging0"=lagging0,  
                                                "ratio1"=ratio1,
                                                "leading1"=leading1,  
                                                "lagging1"=lagging1,  
                                                "ratio2"=ratio2,
                                                "leading2"=leading2,  
                                                "lagging2"=lagging2,  
                                                "ratio3"=ratio3,
                                                "leading3"=leading3,  
                                                "lagging3"=lagging3,  
                                                "ratio4"=ratio4,
                                                "leading4"=leading4,  
                                                "lagging4"=lagging4,  
                                                "ratio5"=ratio5,
                                                "leading5"=leading5,  
                                                "lagging5"=lagging5,  
                                                "ratio6"=ratio6,
                                                "leading6"=leading6,  
                                                "lagging6"=lagging6,  
                                                "mutcnt"=mutcnt))  
     
   }
  
   results <- data.table(results)
   resultsAll <- data.table(resultsAll)
   
   results <- results[, lapply(.SD, function(x) replace(x, is.infinite(x), NA))]
   resultsAll <- resultsAll[, lapply(.SD, function(x) replace(x, is.infinite(x), NA))]
   
   # results[,mutTotal := sum(mutcnt)]
   # resultsAll[, mutTotal := sum(mutcnt)]
   
   # results[, weight := mutcnt/mutTotal]
   # resultsAll[, weight := mutcnt/mutTotal]
   
   resultsLeading <- results[,sum(leading)]
   resultsLagging <- results[,sum(lagging)]
   resultsRatio <- resultsLagging/resultsLeading
   resultsAllLeading <- resultsAll[,sum(leading)]
   resultsAllLagging <- resultsAll[,sum(lagging)]
   resultsAllRatio <- resultsAllLagging/resultsAllLeading
   
   resultsLeading0 <- results[,sum(leading0, na.rm = TRUE)]
   resultsLagging0 <- results[,sum(lagging0, na.rm = TRUE)]
   resultsRatio0 <- resultsLagging0/resultsLeading0
   resultsAllLeading0 <- resultsAll[,sum(leading0, na.rm = TRUE)]
   resultsAllLagging0 <- resultsAll[,sum(lagging0, na.rm = TRUE)]
   resultsAllRatio0 <- resultsAllLagging0/resultsAllLeading0
   
   resultsLeading1 <- results[,sum(leading1, na.rm = TRUE)]
   resultsLagging1 <- results[,sum(lagging1, na.rm = TRUE)]
   resultsRatio1 <- resultsLagging1/resultsLeading1
   resultsAllLeading1 <- resultsAll[,sum(leading1, na.rm = TRUE)]
   resultsAllLagging1 <- resultsAll[,sum(lagging1, na.rm = TRUE)]
   resultsAllRatio1 <- resultsAllLagging1/resultsAllLeading1
   
   resultsLeading2 <- results[,sum(leading2, na.rm = TRUE)]
   resultsLagging2 <- results[,sum(lagging2, na.rm = TRUE)]
   resultsRatio2 <- resultsLagging2/resultsLeading2
   resultsAllLeading2 <- resultsAll[,sum(leading2, na.rm = TRUE)]
   resultsAllLagging2 <- resultsAll[,sum(lagging2, na.rm = TRUE)]
   resultsAllRatio2 <- resultsAllLagging2/resultsAllLeading2
   
   resultsLeading3 <- results[,sum(leading3, na.rm = TRUE)]
   resultsLagging3 <- results[,sum(lagging3, na.rm = TRUE)]
   resultsRatio3 <- resultsLagging3/resultsLeading3
   resultsAllLeading3 <- resultsAll[,sum(leading3, na.rm = TRUE)]
   resultsAllLagging3 <- resultsAll[,sum(lagging3, na.rm = TRUE)]
   resultsAllRatio3 <- resultsAllLagging3/resultsAllLeading3
   
   resultsLeading4 <- results[,sum(leading4, na.rm = TRUE)]
   resultsLagging4 <- results[,sum(lagging4, na.rm = TRUE)]
   resultsRatio4 <- resultsLagging4/resultsLeading4
   resultsAllLeading4 <- resultsAll[,sum(leading4, na.rm = TRUE)]
   resultsAllLagging4 <- resultsAll[,sum(lagging4, na.rm = TRUE)]
   resultsAllRatio4 <- resultsAllLagging4/resultsAllLeading4
   
   resultsLeading5 <- results[,sum(leading5, na.rm = TRUE)]
   resultsLagging5 <- results[,sum(lagging5, na.rm = TRUE)]
   resultsRatio5 <- resultsLagging5/resultsLeading5
   resultsAllLeading5 <- resultsAll[,sum(leading5, na.rm = TRUE)]
   resultsAllLagging5 <- resultsAll[,sum(lagging5, na.rm = TRUE)]
   resultsAllRatio5 <- resultsAllLagging5/resultsAllLeading5
   
   resultsLeading6 <- results[,sum(leading6, na.rm = TRUE)]
   resultsLagging6 <- results[,sum(lagging6, na.rm = TRUE)]
   resultsRatio6 <- resultsLagging6/resultsLeading6
   resultsAllLeading6 <- resultsAll[,sum(leading6, na.rm = TRUE)]
   resultsAllLagging6 <- resultsAll[,sum(lagging6, na.rm = TRUE)]
   resultsAllRatio6 <- resultsAllLagging6/resultsAllLeading6
   
   res <- rbind(res, data.frame("Cancer"=cancers[i]$cancer,
                                 "Sample"=samples[j]$sample,
                                 "Ratio"=resultsRatio,
                                 "Leading"=resultsLeading,
                                 "Lagging"=resultsLagging,
                                 "RatioAll"=resultsAllRatio,
                                 "LeadingAll"=resultsAllLeading,
                                 "LaggingAll"=resultsAllLagging,
                                 "RatioDiff"=resultsRatio/resultsAllRatio, 
                                 "Ratio0"=resultsRatio0,
                                 "Leading0"=resultsLeading0,
                                 "Lagging0"=resultsLagging0,
                                 "RatioAll0"=resultsAllRatio0,
                                 "LeadingAll0"=resultsAllLeading0,
                                 "LaggingAll0"=resultsAllLagging0,
                                 "RatioDiff0"=resultsRatio0/resultsAllRatio0, 
                                 "Ratio1"=resultsRatio1,
                                 "Leading1"=resultsLeading1,
                                 "Lagging1"=resultsLagging1,
                                 "RatioAll1"=resultsAllRatio1,
                                 "LeadingAll1"=resultsAllLeading1,
                                 "LaggingAll1"=resultsAllLagging1,
                                 "RatioDiff1"=resultsRatio1/resultsAllRatio1, 
                                 "Ratio2"=resultsRatio2,
                                 "Leading2"=resultsLeading2,
                                 "Lagging2"=resultsLagging2,
                                 "RatioAll2"=resultsAllRatio2,
                                 "LeadingAll2"=resultsAllLeading2,
                                 "LaggingAll2"=resultsAllLagging2,
                                 "RatioDiff2"=resultsRatio2/resultsAllRatio2, 
                                 "Ratio3"=resultsRatio3,
                                 "Leading3"=resultsLeading3,
                                 "Lagging3"=resultsLagging3,
                                 "RatioAll3"=resultsAllRatio3,
                                 "LeadingAll3"=resultsAllLeading3,
                                 "LaggingAll3"=resultsAllLagging3,
                                 "RatioDiff3"=resultsRatio3/resultsAllRatio3, 
                                 "Ratio4"=resultsRatio4,
                                 "Leading4"=resultsLeading4,
                                 "Lagging4"=resultsLagging4,
                                 "RatioAll4"=resultsAllRatio4,
                                 "LeadingAll4"=resultsAllLeading4,
                                 "LaggingAll4"=resultsAllLagging4,
                                 "RatioDiff4"=resultsRatio4/resultsAllRatio4, 
                                 "Ratio5"=resultsRatio5,
                                 "Leading5"=resultsLeading5,
                                 "Lagging5"=resultsLagging5,
                                 "RatioAll5"=resultsAllRatio5,
                                 "LeadingAll5"=resultsAllLeading5,
                                 "LaggingAll5"=resultsAllLagging5,
                                 "RatioDiff5"=resultsRatio5/resultsAllRatio5, 
                                 "Ratio6"=resultsRatio6,
                                 "Leading6"=resultsLeading6,
                                 "Lagging6"=resultsLagging6,
                                 "RatioAll6"=resultsAllRatio6,
                                 "LeadingAll6"=resultsAllLeading6,
                                 "LaggingAll6"=resultsAllLagging6,
                                 "RatioDiff6"=resultsRatio6/resultsAllRatio6, 
                                 "GordeninEnrichment"=enrichGordenin$APOBEC_ENRICHMENT
                                 ))  
    
  }
}

write.csv(res,paste0(OUTPUT_DIR,"/RTStrand/coefsFinalStrand.csv"))

