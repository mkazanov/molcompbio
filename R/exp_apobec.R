library(data.table)
library(ggplot2)
library(reshape2)

### Input data

OUTPUT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR"
apobec <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/results_exp_APOBEC.txt", sep = '\t',header = TRUE)
apobec <- data.table(apobec)
other <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/results_exp_OTHER.txt", sep = '\t',header = TRUE)
other <- data.table(other)

tcwNum <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/TCW_in_EXPbins_cluster.csv")
tcwNum <- data.table(tcwNum)

allNum <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL_in_EXPbins_cluster.csv")
allNum <- data.table(allNum)

enrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/enrichment.txt", sep = '\t')
enrichment <- data.table(enrichment)

gordeninEnrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results/sample_enrichment.txt", sep = '\t')
gordeninEnrichment <- data.table(gordeninEnrichment)

# Processing data

cancers <- data.table("cancer"=unique(apobec[,Cancer]))

res <- data.frame()

# Create subdirs for cancer types

dir.create(file.path(OUTPUT_DIR,'EXP'))
dir.create(file.path(OUTPUT_DIR,'EXPstrand'))
dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand'),'RAWall'))
dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand'),'RATIOall'))
dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand'),'RAWconsist'))
dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand'),'RATIOconsist'))

for(i in 1:nrow(cancers))
{
  
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXP'), cancers[i]$cancer))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXP/',cancers[i]$cancer), "APOBEC"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXP/',cancers[i]$cancer), "OTHER"))
  
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RAWall'), cancers[i]$cancer))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RATIOall'), cancers[i]$cancer))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RAWconsist'), cancers[i]$cancer))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RATIOconsist'), cancers[i]$cancer))
  
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RAWall/',cancers[i]$cancer), "APOBEC"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RAWall/',cancers[i]$cancer), "OTHER"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RATIOall/',cancers[i]$cancer), "APOBEC"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RATIOall/',cancers[i]$cancer), "OTHER"))

  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RAWconsist/',cancers[i]$cancer), "APOBEC"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RAWconsist/',cancers[i]$cancer), "OTHER"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RATIOconsist/',cancers[i]$cancer), "APOBEC"))
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand/RATIOconsist/',cancers[i]$cancer), "OTHER"))
}

res <- data.frame()

for(i in 1:nrow(cancers))
{
  samples <- data.table("sample"=unique(apobec[Cancer == cancers[i]$cancer, Sample]))

for (j in 1:nrow(samples))
{
  tcwNumSample <- tcwNum[Cancer == cancers[i]$cancer & Sample == samples[j]$sample] 
  allNumSample <- allNum[Cancer == cancers[i]$cancer & Sample == samples[j]$sample]
  enrich <- enrichment[cancer == cancers[i]$cancer & sample == samples[j]$sample]
  enrichChar <- paste0(gsub("\\.","_",round(enrich$Enrichment_exclude_TCW,2)),"_")
  enrichGordenin <- gordeninEnrichment[CANCER_TYPE == cancers[i]$cancer & SAMPLE == samples[j]$sample]  
  
  # APOBEC 
  dt <- apobec[Sample == samples[j]$sample & Cancer == cancers[i]$cancer]
  dt <- dt[ExpressionBin %in% c(1,2,3,4,5,6,7)]
  dt1 <- merge(dt, tcwNumSample, by = "ExpressionBin")
  dt1 <- dt1[, MutationDensity := MutationCnt/TargetCnt]
  dt1[, DensityTotal := sum(MutationDensity)]
  dt1[, NormalizedDensity := MutationDensity/DensityTotal]
  
  retAPOBEC <- lm(NormalizedDensity ~ ExpressionBin,dt1)
  
  p <- ggplot(dt1, aes(x=ExpressionBin, y=NormalizedDensity)) + 
    geom_bar(stat="identity") +
    geom_smooth(method='lm',formula=y~x) +
    ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
                   ", Enrich: ", round(enrich$Enrichment_exclude_TCW,2), ", Coef = ",
                   formatC(retAPOBEC$coefficients[[2]], format = "e", digits = 2))) +
    theme(plot.title = element_text(size = 10))
  
  ggsave(paste0(OUTPUT_DIR,"/EXP/",cancers[i]$cancer, "/APOBEC/", enrichChar, samples[j]$sample, ".jpg"), plot = p, device = "jpeg")
 
  # Strand
  dt1 <- dt1[, PlusConsistDensity := PlusStrandConsistent/TargetCnt]
  dt1 <- dt1[, MinusConsistDensity := MinusStrandConsistent/TargetCnt]
  dt1 <- dt1[, StrandRatio := PlusStrandConsistent/MinusStrandConsistent]
  dt1[!is.finite(StrandRatio), StrandRatio := 0]
  dt2 <- dt1[, .(ExpressionBin,PlusConsistDensity,MinusConsistDensity)]
  dt2 <- melt(dt2, id.vars = "ExpressionBin")
  
  retPlusConsist <- lm(PlusConsistDensity ~ ExpressionBin,dt1)
  retMinusConsist <- lm(MinusConsistDensity ~ ExpressionBin,dt1)
  
  p <- ggplot(dt2, aes(x=ExpressionBin, y=value, fill=variable)) + 
    geom_bar(stat="identity", position="dodge") +
    geom_smooth(method='lm') +
    ylab("Mutation density") +
    ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
                   ", Enrich: ", round(enrich$Enrichment_exclude_TCW,2),
                   ", CoefPlusConsist = ", formatC(retPlusConsist$coefficients[[2]], format = "e", digits = 2),
                   ", CoefMinusConsist = ", formatC(retMinusConsist$coefficients[[2]], format = "e", digits = 2))) +
    theme(plot.title = element_text(size = 10))
  
  ggsave(paste0(OUTPUT_DIR,"/EXPstrand/RAWconsist/",cancers[i]$cancer, "/APOBEC/", enrichChar, samples[j]$sample, ".jpg"), plot = p, device = "jpeg")
  
  retRatio <- lm(StrandRatio ~ ExpressionBin,dt1)
  
  p <- ggplot(dt1, aes(x=ExpressionBin, y=StrandRatio)) + 
    geom_bar(stat="identity") +
    geom_smooth(method='lm',formula=y~x) +
    ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
                   ", Enrich: ", round(enrich$Enrichment_exclude_TCW,2), ", Coef = ",round(retRatio$coefficients[[2]],2))) +
    theme(plot.title = element_text(size = 10))
  
  ggsave(paste0(OUTPUT_DIR,"/EXPstrand/RATIOconsist/",cancers[i]$cancer, "/APOBEC/", enrichChar, samples[j]$sample, ".jpg"), plot = p, device = "jpeg")
  

  dt1 <- dt1[, PlusAllDensity := PlusStrandAll/TargetCnt]
  dt1 <- dt1[, MinusAllDensity := MinusStrandAll/TargetCnt]
  dt1 <- dt1[, StrandRatio := PlusStrandAll/MinusStrandAll]
  dt1[!is.finite(StrandRatio), StrandRatio := 0]
  dt2 <- dt1[, .(ExpressionBin,PlusAllDensity,MinusAllDensity)]
  dt2 <- melt(dt2, id.vars = "ExpressionBin")
  
  retPlusAll <- lm(PlusAllDensity ~ ExpressionBin,dt1)
  retMinusAll <- lm(MinusAllDensity ~ ExpressionBin,dt1)
  
  p <- ggplot(dt2, aes(x=ExpressionBin, y=value, fill=variable)) + 
    geom_bar(stat="identity", position="dodge") +
    geom_smooth(method='lm') +
    ylab("Mutation density") +
    ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
                   ", Enrich: ", round(enrich$Enrichment_exclude_TCW,2),
                   ", CoefPlusAll = ", formatC(retPlusConsist$coefficients[[2]], format = "e", digits = 2),
                   ", CoefMinusAll = ", formatC(retMinusConsist$coefficients[[2]], format = "e", digits = 2))) +
    theme(plot.title = element_text(size = 10))
  
  ggsave(paste0(OUTPUT_DIR,"/EXPstrand/RAWall/",cancers[i]$cancer, "/APOBEC/", enrichChar, samples[j]$sample, ".jpg"), plot = p, device = "jpeg")
  
  retRatio <- lm(StrandRatio ~ ExpressionBin,dt1)
  
  p <- ggplot(dt1, aes(x=ExpressionBin, y=StrandRatio)) + 
    geom_bar(stat="identity") +
    geom_smooth(method='lm',formula=y~x) +
    ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
                   ", Enrich: ", round(enrich$Enrichment_exclude_TCW,2), ", Coef = ",round(retRatio$coefficients[[2]],2))) +
    theme(plot.title = element_text(size = 10))
  
  ggsave(paste0(OUTPUT_DIR,"/EXPstrand/RATIOall/",cancers[i]$cancer, "/APOBEC/", enrichChar, samples[j]$sample, ".jpg"), plot = p, device = "jpeg")
  
  
    
  # Other
  dt <- other[Sample == samples[j]$sample & Cancer == cancers[i]$cancer]
  dt <- dt[ExpressionBin %in% c(1,2,3,4,5,6,7)]
  dt1 <- merge(dt, allNumSample, by = "ExpressionBin")
  dt1 <- dt1[, MutationDensity := MutationCnt/TargetCnt]
  dt1[, DensityTotal := sum(MutationDensity)]
  dt1[, NormalizedDensity := MutationDensity/DensityTotal]

  retOTHER <- lm(NormalizedDensity ~ ExpressionBin,dt1)
  
  p <- ggplot(dt1, aes(x=ExpressionBin, y=NormalizedDensity)) + 
    geom_bar(stat="identity") +
    geom_smooth(method='lm',formula=y~x) +
    ggtitle(paste0("Cancer: ", cancers[i]$cancer, ", Sample: ", samples[j]$sample, 
                   ", Enrich: ", enrich$Enrichment_exclude_TCW, ", Coef = ",retOTHER$coefficients[[2]]))
  
  ggsave(paste0(OUTPUT_DIR,"/EXP/",cancers[i]$cancer, "/OTHER/", enrichChar, samples[j]$sample, ".jpg"), plot = p, device = "jpeg")
   
  res <- rbind(res, data.frame("Cancer"=cancers[i]$cancer,
                               "Sample"=samples[j]$sample,
                               "ApobecSlope"=retAPOBEC$coefficients[[2]],
                               "OtherSlope"=retOTHER$coefficients[[2]],
                               "PlusConsistSlope"=retPlusConsist$coefficients[[2]],
                               "MinusConsistSlope"=retMinusConsist$coefficients[[2]],
                               "PlusAllSlope"=retPlusConsist$coefficients[[2]],
                               "MinusAllSlope"=retMinusConsist$coefficients[[2]],
                               "StrandRatioSlope"=retRatio$coefficients[[2]],
                               "GordeninEnrichment"=enrichGordenin$APOBEC_ENRICHMENT,
                               "Enrichment"=enrich$Enrichment,
                               "Enrichment_exclude_TCW"=enrich$Enrichment_exclude_TCW))  
  
}
}

write.csv(res,"/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR/EXP/coefs.csv")
