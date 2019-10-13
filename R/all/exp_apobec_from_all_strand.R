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

dir.create(file.path(OUTPUT_DIR,'EXPstrand'))

for(i in 1:nrow(cancers))
{
  
  dir.create(file.path(paste0(OUTPUT_DIR,'/EXPstrand'), cancers[i]$cancer))
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
      tcx[,CodingCnt:=ifelse(senseStrand==1,MutationCnt,0)]
      tcx[,TemplateCnt:=ifelse(senseStrand==0,MutationCnt,0)]
      
      dt <- tcx[,.("MutationCnt"=sum(MutationCnt),
                   "CodingCnt"=sum(CodingCnt),
                   "TemplateCnt"=sum(TemplateCnt)),
                by=.(expbin)]
      dt <- dt[expbin %in% c(1,2,3,4,5,6,7)]
      
      tcxNum <- dataTrg[Motif == m]
      tcxNum <- tcxNum[expbin %in% c(1,2,3,4,5,6,7)]
      tcxNum[,CodingCnt:=ifelse(senseStrand==1,TargetCnt,0)]
      tcxNum[,TemplateCnt:=ifelse(senseStrand==0,TargetCnt,0)]
      tcxNum <- tcxNum[,.("TargetCnt"=sum(TargetCnt),
                          "TargetCodingCnt"=sum(CodingCnt),
                          "TargetTemplateCnt"=sum(TemplateCnt)), by=.(expbin)]
   
      coding <- dt[,sum(CodingCnt)]
      template <- dt[,sum(TemplateCnt)]
      ratio <- coding/template
      targetsCoding <- tcxNum[, sum(TargetCodingCnt)]
      targetsTemplate <- tcxNum[, sum(TargetTemplateCnt)]
      codingDensity <- coding/targetsCoding
      templateDensity <- template/targetsTemplate
      ratioDensity <- codingDensity/templateDensity
      
      coding1 <- dt[expbin == 1,CodingCnt]
      if(length(coding1) == 0) {coding1 <- NA}
      template1 <- dt[expbin == 1,TemplateCnt]
      if(length(template1) == 0) {template1 <- NA}
      ratio1 <- coding1/template1
      if(length(ratio1) == 0){ratio1 <- NA}
      targetsCoding1 <- tcxNum[expbin == 1, TargetCodingCnt]
      targetsTemplate1 <- tcxNum[expbin == 1, TargetTemplateCnt]
      codingDensity1 <- coding1/targetsCoding1
      templateDensity1 <- template1/targetsTemplate1
      ratioDensity1 <- codingDensity1/templateDensity1
      
      coding2 <- dt[expbin == 2,CodingCnt]
      if(length(coding2) == 0) {coding2 <- NA}
      template2 <- dt[expbin == 2,TemplateCnt]
      if(length(template2) == 0) {template2 <- NA}
      ratio2 <- coding2/template2
      if(length(ratio2) == 0){ratio2 <- NA}
      targetsCoding2 <- tcxNum[expbin == 2, TargetCodingCnt]
      targetsTemplate2 <- tcxNum[expbin == 2, TargetTemplateCnt]
      codingDensity2 <- coding2/targetsCoding2
      templateDensity2 <- template2/targetsTemplate2
      ratioDensity2 <- codingDensity2/templateDensity2
      
      coding3 <- dt[expbin == 3,CodingCnt]
      if(length(coding3) == 0) {coding3 <- NA}
      template3 <- dt[expbin == 3,TemplateCnt]
      if(length(template3) == 0) {template3 <- NA}
      ratio3 <- coding3/template3
      if(length(ratio3) == 0){ratio3 <- NA}
      targetsCoding3 <- tcxNum[expbin == 3, TargetCodingCnt]
      targetsTemplate3 <- tcxNum[expbin == 3, TargetTemplateCnt]
      codingDensity3 <- coding3/targetsCoding3
      templateDensity3 <- template3/targetsTemplate3
      ratioDensity3 <- codingDensity3/templateDensity3
      
      coding4 <- dt[expbin == 4,CodingCnt]
      if(length(coding4) == 0) {coding4 <- NA}
      template4 <- dt[expbin == 4,TemplateCnt]
      if(length(template4) == 0) {template4 <- NA}
      ratio4 <- coding4/template4
      if(length(ratio4) == 0){ratio4 <- NA}
      targetsCoding4 <- tcxNum[expbin == 4, TargetCodingCnt]
      targetsTemplate4 <- tcxNum[expbin == 4, TargetTemplateCnt]
      codingDensity4 <- coding4/targetsCoding4
      templateDensity4 <- template4/targetsTemplate4
      ratioDensity4 <- codingDensity4/templateDensity4
      
      coding5 <- dt[expbin == 5,CodingCnt]
      if(length(coding5) == 0) {coding5 <- NA}
      template5 <- dt[expbin == 5,TemplateCnt]
      if(length(template5) == 0) {template5 <- NA}
      ratio5 <- coding5/template5
      if(length(ratio5) == 0){ratio5 <- NA}
      targetsCoding5 <- tcxNum[expbin == 5, TargetCodingCnt]
      targetsTemplate5 <- tcxNum[expbin == 5, TargetTemplateCnt]
      codingDensity5 <- coding5/targetsCoding5
      templateDensity5 <- template5/targetsTemplate5
      ratioDensity5 <- codingDensity5/templateDensity5
      
      coding6 <- dt[expbin == 6,CodingCnt]
      if(length(coding6) == 0) {coding6 <- NA}
      template6 <- dt[expbin == 6,TemplateCnt]
      if(length(template6) == 0) {template6 <- NA}
      ratio6 <- coding6/template6
      if(length(ratio6) == 0){ratio6 <- NA}
      targetsCoding6 <- tcxNum[expbin == 6, TargetCodingCnt]
      targetsTemplate6 <- tcxNum[expbin == 6, TargetTemplateCnt]
      codingDensity6 <- coding6/targetsCoding6
      templateDensity6 <- template6/targetsTemplate6
      ratioDensity6 <- codingDensity6/templateDensity6
      
      coding7 <- dt[expbin == 7,CodingCnt]
      if(length(coding7) == 0) {coding7 <- NA}
      template7 <- dt[expbin == 7,TemplateCnt]
      if(length(template7) == 0) {template7 <- NA}
      ratio7 <- coding7/template7
      if(length(ratio7) == 0){ratio7 <- NA}
      targetsCoding7 <- tcxNum[expbin == 7, TargetCodingCnt]
      targetsTemplate7 <- tcxNum[expbin == 7, TargetTemplateCnt]
      codingDensity7 <- coding7/targetsCoding7
      templateDensity7 <- template7/targetsTemplate7
      ratioDensity7 <- codingDensity7/templateDensity7
      
      mutcnt <- dt[, sum(MutationCnt)]
      results <- rbind(results, data.frame("motif"=m,
                                           "ratio"=ratio, 
                                           "coding"=coding,
                                           "template"=template,
                                           "targetsCoding"=targetsCoding,
                                           "targetsTemplate"=targetsTemplate,
                                           "codingDensity"=codingDensity,
                                           "templateDensity"=templateDensity,
                                           "ratioDensity"=ratioDensity,
                                           "ratio1"=ratio1,
                                           "coding1"=coding1,  
                                           "template1"=template1,  
                                           "targetsCoding1"=targetsCoding1,
                                           "targetsTemplate1"=targetsTemplate1,
                                           "codingDensity1"=codingDensity1,
                                           "templateDensity1"=templateDensity1,
                                           "ratioDensity1"=ratioDensity1,
                                           "ratio2"=ratio2,
                                           "coding2"=coding2,  
                                           "template2"=template2, 
                                           "targetsCoding2"=targetsCoding2,
                                           "targetsTemplate2"=targetsTemplate2,
                                           "codingDensity2"=codingDensity2,
                                           "templateDensity2"=templateDensity2,
                                           "ratioDensity2"=ratioDensity2,
                                           "ratio3"=ratio3,
                                           "coding3"=coding3,  
                                           "template3"=template3,  
                                           "targetsCoding3"=targetsCoding3,
                                           "targetsTemplate3"=targetsTemplate3,
                                           "codingDensity3"=codingDensity3,
                                           "templateDensity3"=templateDensity3,
                                           "ratioDensity3"=ratioDensity3,
                                           "ratio4"=ratio4,
                                           "coding4"=coding4,  
                                           "template4"=template4, 
                                           "targetsCoding4"=targetsCoding4,
                                           "targetsTemplate4"=targetsTemplate4,
                                           "codingDensity4"=codingDensity4,
                                           "templateDensity4"=templateDensity4,
                                           "ratioDensity4"=ratioDensity4,
                                           "ratio5"=ratio5,
                                           "coding5"=coding5,  
                                           "template5"=template5,  
                                           "targetsCoding5"=targetsCoding5,
                                           "targetsTemplate5"=targetsTemplate5,
                                           "codingDensity5"=codingDensity5,
                                           "templateDensity5"=templateDensity5,
                                           "ratioDensity5"=ratioDensity5,
                                           "ratio6"=ratio6,
                                           "coding6"=coding6,  
                                           "template6"=template6,  
                                           "targetsCoding6"=targetsCoding6,
                                           "targetsTemplate6"=targetsTemplate6,
                                           "codingDensity6"=codingDensity6,
                                           "templateDensity6"=templateDensity6,
                                           "ratioDensity6"=ratioDensity6,
                                           "ratio7"=ratio7,
                                           "coding7"=coding7,  
                                           "template7"=template7,  
                                           "targetsCoding7"=targetsCoding7,
                                           "targetsTemplate7"=targetsTemplate7,
                                           "codingDensity7"=codingDensity7,
                                           "templateDensity7"=templateDensity7,
                                           "ratioDensity7"=ratioDensity7,
                                           "mutcnt"=mutcnt))
    }
    
    resultsAll <- data.frame()  
      
    C2motifList <- c("ACA","ACC","ACG","GCT","CCT","ACT","CCA","CCG","GCA","GCG","GCC","CCC")
    
    # All other triplets with middle C
    for(m in C2motifList) 
    {
      xt <- dataMut[Motif == m & mutAllele %in% c("T","G")]
      xt[,CodingCnt:=ifelse(senseStrand==1,MutationCnt,0)]
      xt[,TemplateCnt:=ifelse(senseStrand==0,MutationCnt,0)]
      
      dt <- xt[,.("MutationCnt"=sum(MutationCnt),
                  "CodingCnt"=sum(CodingCnt),
                  "TemplateCnt"=sum(TemplateCnt)),
               by=.(expbin)]
      dt <- dt[expbin %in% c(1,2,3,4,5,6,7)]
      
      coding <- dt[,sum(CodingCnt)]
      template <- dt[,sum(TemplateCnt)]
      ratio <- coding/template
      
      coding1 <- dt[expbin == 1,CodingCnt]
      if(length(coding1) == 0) {coding1 <- NA}
      template1 <- dt[expbin == 1,TemplateCnt]
      if(length(template1) == 0) {template1 <- NA}
      ratio1 <- coding1/template1
      if(length(ratio1) == 0){ratio1 <- NA}
      
      coding2 <- dt[expbin == 2,CodingCnt]
      if(length(coding2) == 0) {coding2 <- NA}
      template2 <- dt[expbin == 2,TemplateCnt]
      if(length(template2) == 0) {template2 <- NA}
      ratio2 <- coding2/template2
      if(length(ratio2) == 0){ratio2 <- NA}
      
      coding3 <- dt[expbin == 3,CodingCnt]
      if(length(coding3) == 0) {coding3 <- NA}
      template3 <- dt[expbin == 3,TemplateCnt]
      if(length(template3) == 0) {template3 <- NA}
      ratio3 <- coding3/template3
      if(length(ratio3) == 0){ratio3 <- NA}
      
      coding4 <- dt[expbin == 4,CodingCnt]
      if(length(coding4) == 0) {coding4 <- NA}
      template4 <- dt[expbin == 4,TemplateCnt]
      if(length(template4) == 0) {template4 <- NA}
      ratio4 <- coding4/template4
      if(length(ratio4) == 0){ratio4 <- NA}
      
      coding5 <- dt[expbin == 5,CodingCnt]
      if(length(coding5) == 0) {coding5 <- NA}
      template5 <- dt[expbin == 5,TemplateCnt]
      if(length(template5) == 0) {template5 <- NA}
      ratio5 <- coding5/template5
      if(length(ratio5) == 0){ratio5 <- NA}
      
      coding6 <- dt[expbin == 6,CodingCnt]
      if(length(coding6) == 0) {coding6 <- NA}
      template6 <- dt[expbin == 6,TemplateCnt]
      if(length(template6) == 0) {template6 <- NA}
      ratio6 <- coding6/template6
      if(length(ratio6) == 0){ratio6 <- NA}
      
      coding7 <- dt[expbin == 7,CodingCnt]
      if(length(coding7) == 0) {coding7 <- NA}
      template7 <- dt[expbin == 7,TemplateCnt]
      if(length(template7) == 0) {template7 <- NA}
      ratio7 <- coding7/template7
      if(length(ratio7) == 0){ratio7 <- NA}

      mutcnt <- dt[, sum(MutationCnt)]
      
      resultsAll <- rbind(resultsAll, data.frame("motif"=m,
                                                 "ratio"=ratio, 
                                                 "coding"=coding,
                                                 "template"=template,
                                                 "ratio1"=ratio1,
                                                 "coding1"=coding1,  
                                                 "template1"=template1,  
                                                 "ratio2"=ratio2,
                                                 "coding2"=coding2,  
                                                 "template2"=template2,  
                                                 "ratio3"=ratio3,
                                                 "coding3"=coding3,  
                                                 "template3"=template3,  
                                                 "ratio4"=ratio4,
                                                 "coding4"=coding4,  
                                                 "template4"=template4,  
                                                 "ratio5"=ratio5,
                                                 "coding5"=coding5,  
                                                 "template5"=template5,  
                                                 "ratio6"=ratio6,
                                                 "coding6"=coding6,  
                                                 "template6"=template6,  
                                                 "ratio7"=ratio7,
                                                 "coding7"=coding7,  
                                                 "template7"=template7,  
                                                 "mutcnt"=mutcnt))  
      
    }      
      
    results <- data.table(results)
    resultsAll <- data.table(resultsAll)
    
    results <- results[, lapply(.SD, function(x) replace(x, is.infinite(x), NA))]
    resultsAll <- resultsAll[, lapply(.SD, function(x) replace(x, is.infinite(x), NA))]
    
    resultsCoding <- results[,sum(coding)]
    resultsTemplate <- results[,sum(template)]
    resultsRatio <- resultsCoding/resultsTemplate
    resultsTargetsCoding <- results[, sum(targetsCoding)]
    resultsTargetsTemplate <- results[, sum(targetsTemplate)]
    resultsCodingDensity <- resultsCoding/resultsTargetsCoding
    resultsTemplateDensity <- resultsTemplate/resultsTargetsTemplate
    resultsRatioDensity <- resultsCodingDensity/resultsTemplateDensity
    resultsAllCoding <- resultsAll[,sum(coding)]
    resultsAllTemplate <- resultsAll[,sum(template)]
    resultsAllRatio <- resultsAllCoding/resultsAllTemplate
    
    resultsCoding1 <- results[,sum(coding1, na.rm = TRUE)]
    resultsTemplate1 <- results[,sum(template1, na.rm = TRUE)]
    resultsRatio1 <- resultsCoding1/resultsTemplate1
    resultsTargetsCoding1 <- results[, sum(targetsCoding1)]
    resultsTargetsTemplate1 <- results[, sum(targetsTemplate1)]
    resultsCodingDensity1 <- resultsCoding1/resultsTargetsCoding1
    resultsTemplateDensity1 <- resultsTemplate1/resultsTargetsTemplate1
    resultsRatioDensity1 <- resultsCodingDensity1/resultsTemplateDensity1
    resultsAllCoding1 <- resultsAll[,sum(coding1, na.rm = TRUE)]
    resultsAllTemplate1 <- resultsAll[,sum(template1, na.rm = TRUE)]
    resultsAllRatio1 <- resultsAllCoding1/resultsAllTemplate1
    
    resultsCoding2 <- results[,sum(coding2, na.rm = TRUE)]
    resultsTemplate2 <- results[,sum(template2, na.rm = TRUE)]
    resultsRatio2 <- resultsCoding2/resultsTemplate2
    resultsTargetsCoding2 <- results[, sum(targetsCoding2)]
    resultsTargetsTemplate2 <- results[, sum(targetsTemplate2)]
    resultsCodingDensity2 <- resultsCoding2/resultsTargetsCoding2
    resultsTemplateDensity2 <- resultsTemplate2/resultsTargetsTemplate2
    resultsRatioDensity2 <- resultsCodingDensity2/resultsTemplateDensity2
    resultsAllCoding2 <- resultsAll[,sum(coding2, na.rm = TRUE)]
    resultsAllTemplate2 <- resultsAll[,sum(template2, na.rm = TRUE)]
    resultsAllRatio2 <- resultsAllCoding2/resultsAllTemplate2
    
    resultsCoding3 <- results[,sum(coding3, na.rm = TRUE)]
    resultsTemplate3 <- results[,sum(template3, na.rm = TRUE)]
    resultsRatio3 <- resultsCoding3/resultsTemplate3
    resultsTargetsCoding3 <- results[, sum(targetsCoding3)]
    resultsTargetsTemplate3 <- results[, sum(targetsTemplate3)]
    resultsCodingDensity3 <- resultsCoding3/resultsTargetsCoding3
    resultsTemplateDensity3 <- resultsTemplate3/resultsTargetsTemplate3
    resultsRatioDensity3 <- resultsCodingDensity3/resultsTemplateDensity3
    resultsAllCoding3 <- resultsAll[,sum(coding3, na.rm = TRUE)]
    resultsAllTemplate3 <- resultsAll[,sum(template3, na.rm = TRUE)]
    resultsAllRatio3 <- resultsAllCoding3/resultsAllTemplate3
    
    resultsCoding4 <- results[,sum(coding4, na.rm = TRUE)]
    resultsTemplate4 <- results[,sum(template4, na.rm = TRUE)]
    resultsRatio4 <- resultsCoding4/resultsTemplate4
    resultsTargetsCoding4 <- results[, sum(targetsCoding4)]
    resultsTargetsTemplate4 <- results[, sum(targetsTemplate4)]
    resultsCodingDensity4 <- resultsCoding4/resultsTargetsCoding4
    resultsTemplateDensity4 <- resultsTemplate4/resultsTargetsTemplate4
    resultsRatioDensity4 <- resultsCodingDensity4/resultsTemplateDensity4
    resultsAllCoding4 <- resultsAll[,sum(coding4, na.rm = TRUE)]
    resultsAllTemplate4 <- resultsAll[,sum(template4, na.rm = TRUE)]
    resultsAllRatio4 <- resultsAllCoding4/resultsAllTemplate4

    resultsCoding5 <- results[,sum(coding5, na.rm = TRUE)]
    resultsTemplate5 <- results[,sum(template5, na.rm = TRUE)]
    resultsRatio5 <- resultsCoding5/resultsTemplate5
    resultsTargetsCoding5 <- results[, sum(targetsCoding5)]
    resultsTargetsTemplate5 <- results[, sum(targetsTemplate5)]
    resultsCodingDensity5 <- resultsCoding5/resultsTargetsCoding5
    resultsTemplateDensity5 <- resultsTemplate5/resultsTargetsTemplate5
    resultsRatioDensity5 <- resultsCodingDensity5/resultsTemplateDensity5
    resultsAllCoding5 <- resultsAll[,sum(coding5, na.rm = TRUE)]
    resultsAllTemplate5 <- resultsAll[,sum(template5, na.rm = TRUE)]
    resultsAllRatio5 <- resultsAllCoding5/resultsAllTemplate5
    
    resultsCoding6 <- results[,sum(coding6, na.rm = TRUE)]
    resultsTemplate6 <- results[,sum(template6, na.rm = TRUE)]
    resultsRatio6 <- resultsCoding6/resultsTemplate6
    resultsTargetsCoding6 <- results[, sum(targetsCoding6)]
    resultsTargetsTemplate6 <- results[, sum(targetsTemplate6)]
    resultsCodingDensity6 <- resultsCoding6/resultsTargetsCoding6
    resultsTemplateDensity6 <- resultsTemplate6/resultsTargetsTemplate6
    resultsRatioDensity6 <- resultsCodingDensity6/resultsTemplateDensity6
    resultsAllCoding6 <- resultsAll[,sum(coding6, na.rm = TRUE)]
    resultsAllTemplate6 <- resultsAll[,sum(template6, na.rm = TRUE)]
    resultsAllRatio6 <- resultsAllCoding6/resultsAllTemplate6

    resultsCoding7 <- results[,sum(coding7, na.rm = TRUE)]
    resultsTemplate7 <- results[,sum(template7, na.rm = TRUE)]
    resultsRatio7 <- resultsCoding7/resultsTemplate7
    resultsTargetsCoding7 <- results[, sum(targetsCoding7)]
    resultsTargetsTemplate7 <- results[, sum(targetsTemplate7)]
    resultsCodingDensity7 <- resultsCoding7/resultsTargetsCoding7
    resultsTemplateDensity7 <- resultsTemplate7/resultsTargetsTemplate7
    resultsRatioDensity7 <- resultsCodingDensity7/resultsTemplateDensity7
    resultsAllCoding7 <- resultsAll[,sum(coding7, na.rm = TRUE)]
    resultsAllTemplate7 <- resultsAll[,sum(template7, na.rm = TRUE)]
    resultsAllRatio7 <- resultsAllCoding7/resultsAllTemplate7
    
    res <- rbind(res, data.frame("Cancer"=cancers[i]$cancer,
                                 "Sample"=samples[j]$sample,
                                 "Ratio"=resultsRatio,
                                 "Coding"=resultsCoding,
                                 "Template"=resultsTemplate,
                                 "TargetsCoding"=resultsTargetsCoding,
                                 "TargetsTemplate"=resultsTargetsTemplate,
                                 "CodingDensity"=resultsCodingDensity,
                                 "TemplateDensity"=resultsTemplateDensity,
                                 "RatioDensity"=resultsRatioDensity,
                                 "RatioAll"=resultsAllRatio,
                                 "CodingAll"=resultsAllCoding,
                                 "TemplateAll"=resultsAllTemplate,
                                 "RatioDiff"=resultsRatio/resultsAllRatio, 
                                 "Ratio1"=resultsRatio1,
                                 "Coding1"=resultsCoding1,
                                 "Template1"=resultsTemplate1,
                                 "TargetsCoding1"=resultsTargetsCoding1,
                                 "TargetsTemplate1"=resultsTargetsTemplate1,
                                 "CodingDensity1"=resultsCodingDensity1,
                                 "TemplateDensity1"=resultsTemplateDensity1,
                                 "RatioDensity1"=resultsRatioDensity1,
                                 "RatioAll1"=resultsAllRatio1,
                                 "CodingAll1"=resultsAllCoding1,
                                 "TemplateAll1"=resultsAllTemplate1,
                                 "RatioDiff1"=resultsRatio1/resultsAllRatio1, 
                                 "Ratio2"=resultsRatio2,
                                 "Coding2"=resultsCoding2,
                                 "Template2"=resultsTemplate2,
                                 "TargetsCoding2"=resultsTargetsCoding2,
                                 "TargetsTemplate2"=resultsTargetsTemplate2,
                                 "CodingDensity2"=resultsCodingDensity2,
                                 "TemplateDensity2"=resultsTemplateDensity2,
                                 "RatioDensity2"=resultsRatioDensity2,
                                 "RatioAll2"=resultsAllRatio2,
                                 "CodingAll2"=resultsAllCoding2,
                                 "TemplateAll2"=resultsAllTemplate2,
                                 "RatioDiff2"=resultsRatio2/resultsAllRatio2, 
                                 "Ratio3"=resultsRatio3,
                                 "Coding3"=resultsCoding3,
                                 "Template3"=resultsTemplate3,
                                 "TargetsCoding3"=resultsTargetsCoding3,
                                 "TargetsTemplate3"=resultsTargetsTemplate3,
                                 "CodingDensity3"=resultsCodingDensity3,
                                 "TemplateDensity3"=resultsTemplateDensity3,
                                 "RatioDensity3"=resultsRatioDensity3,
                                 "RatioAll3"=resultsAllRatio3,
                                 "CodingAll3"=resultsAllCoding3,
                                 "TemplateAll3"=resultsAllTemplate3,
                                 "RatioDiff3"=resultsRatio3/resultsAllRatio3, 
                                 "Ratio4"=resultsRatio4,
                                 "Coding4"=resultsCoding4,
                                 "Template4"=resultsTemplate4,
                                 "TargetsCoding4"=resultsTargetsCoding4,
                                 "TargetsTemplate4"=resultsTargetsTemplate4,
                                 "CodingDensity4"=resultsCodingDensity4,
                                 "TemplateDensity4"=resultsTemplateDensity4,
                                 "RatioDensity4"=resultsRatioDensity4,
                                 "RatioAll4"=resultsAllRatio4,
                                 "CodingAll4"=resultsAllCoding4,
                                 "TemplateAll4"=resultsAllTemplate4,
                                 "RatioDiff4"=resultsRatio4/resultsAllRatio4, 
                                 "Ratio5"=resultsRatio5,
                                 "Coding5"=resultsCoding5,
                                 "Template5"=resultsTemplate5,
                                 "TargetsCoding5"=resultsTargetsCoding5,
                                 "TargetsTemplate5"=resultsTargetsTemplate5,
                                 "CodingDensity5"=resultsCodingDensity5,
                                 "TemplateDensity5"=resultsTemplateDensity5,
                                 "RatioDensity5"=resultsRatioDensity5,
                                 "RatioAll5"=resultsAllRatio5,
                                 "CodingAll5"=resultsAllCoding5,
                                 "TemplateAll5"=resultsAllTemplate5,
                                 "RatioDiff5"=resultsRatio5/resultsAllRatio5, 
                                 "Ratio6"=resultsRatio6,
                                 "Coding6"=resultsCoding6,
                                 "Template6"=resultsTemplate6,
                                 "TargetsCoding6"=resultsTargetsCoding6,
                                 "TargetsTemplate6"=resultsTargetsTemplate6,
                                 "CodingDensity6"=resultsCodingDensity6,
                                 "TemplateDensity6"=resultsTemplateDensity6,
                                 "RatioDensity6"=resultsRatioDensity6,
                                 "RatioAll6"=resultsAllRatio6,
                                 "CodingAll6"=resultsAllCoding6,
                                 "TemplateAll6"=resultsAllTemplate6,
                                 "RatioDiff6"=resultsRatio6/resultsAllRatio6, 
                                 "Ratio7"=resultsRatio7,
                                 "Coding7"=resultsCoding7,
                                 "Template7"=resultsTemplate7,
                                 "TargetsCoding7"=resultsTargetsCoding7,
                                 "TargetsTemplate7"=resultsTargetsTemplate7,
                                 "CodingDensity7"=resultsCodingDensity7,
                                 "TemplateDensity7"=resultsTemplateDensity7,
                                 "RatioDensity7"=resultsRatioDensity1,
                                 "RatioAll7"=resultsAllRatio7,
                                 "CodingAll7"=resultsAllCoding7,
                                 "TemplateAll7"=resultsAllTemplate7,
                                 "RatioDiff7"=resultsRatio7/resultsAllRatio7, 
                                 "GordeninEnrichment"=enrichGordenin$APOBEC_ENRICHMENT
    ))  
    
  }
}

write.csv(res,paste0(OUTPUT_DIR,"/EXPStrand/coefsFinalExpStrand.csv"))


      
