ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL/"
MOTIFS_LIST <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/all/all_motifs.txt"
OUTPUT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalRall/RTmotifsTCW"

crDNA <- function(dna)
{
  return(stri_reverse(chartr("acgtACGT","tgcaTGCA",dna)))
}


library(data.table)
library(stringr)
library(stringi)
library(ggpubr)

motifs <- read.csv(MOTIFS_LIST, header = FALSE)
motifs <- data.table(motifs)
setnames(motifs,c("Motif"))

samplesEnrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results/sample_enrichment.txt", sep='\t')
samplesEnrichment <- data.table(samplesEnrichment)
samplesEnrichment$SAMPLE <- as.character(samplesEnrichment$SAMPLE)
cancers <- unique(samplesEnrichment[, CANCER_TYPE])
cancers <- data.table(cancers)
setnames(cancers,c("CANCER_TYPE"))
cancers$CANCER_TYPE <- as.character(cancers$CANCER_TYPE)

results <- data.table()
resultsMotif <- data.table()

for(c in 1:nrow(cancers))
{
  dir.create(file.path(OUTPUT_DIR, cancers[c]$CANCER_TYPE))
}


mutInitialTemplate <- data.table()
trgInitialTemplate <- data.table()
for(j in 1:nrow(motifs)){
  m <- as.character(motifs[j]$Motif)
  for(rt in c(0,1,2,3,4,5,6)){
    trgInitialTemplate <- rbind(trgInitialTemplate, data.table("Motif"=m,"ReplicationBin"=rt))
    for(nt in c("A","C","G","T")){
      if(substr(m,2,2) != nt){
        mutInitialTemplate <- rbind(mutInitialTemplate, data.table("Motif"=m,
                                                     "ReplicationBin"=rt, 
                                                     "MutateAllele"=nt))
      }
    }
  }
}

for(c in 1:nrow(cancers)){
  cancer <- cancers[c,CANCER_TYPE]
  samples <- samplesEnrichment[CANCER_TYPE == cancer]
  samples <- samples[order(APOBEC_ENRICHMENT)]
  
  print(cancer)
  
  for(i in 1:nrow(samples)){
  
    sample <- samples[i, SAMPLE]
    aenrich <- samples[i, APOBEC_ENRICHMENT]
    enrichChar <- paste0(gsub("\\.","_",round(aenrich,2)),"_")
    
    print(sample)
  
    muts <- read.csv(paste0(ROOT_DIR,"RTEXP_MUT_",sample,".txt"), sep='\t', header=FALSE)
    muts <- data.table(muts)
    setnames(muts,c("Motif","ReplicationBin","ReplicationStrand","ExpressionBin","ExpressionStrand","MutateAllele","Cnt"))
    trgs <- read.csv(paste0(ROOT_DIR,"RTEXP_TRG_",sample,".txt"), sep='\t', header=FALSE)
    trgs <- data.table(trgs)
    setnames(trgs,c("Motif","ReplicationBin","ReplicationStrand","ExpressionBin","ExpressionStrand","Cnt"))
    
    muts <- muts[ReplicationBin != -2]
    trgs <- trgs[ReplicationBin != -2]
    muts <- muts[ReplicationBin != -1]
    trgs <- trgs[ReplicationBin != -1]

    mutsG <- muts[,.("Cnt"=sum(Cnt)), by = .(Motif,ReplicationBin,MutateAllele)]
    trgsG <- trgs[,.("TrgCnt"=sum(Cnt)), by = .(Motif,ReplicationBin)]
    
    mutTemplate <- copy(mutInitialTemplate)
    trgTemplate <- copy(trgInitialTemplate)
    
    mutsG <- merge(mutTemplate,mutsG,by=c("Motif","ReplicationBin","MutateAllele"), all.x=TRUE)
    mutsG[is.na(Cnt), Cnt:= 0] 
    trgsG <- merge(trgTemplate,trgsG,by=c("Motif","ReplicationBin"), all.x=TRUE)
    
    data <- merge(trgsG,mutsG,by=c("Motif","ReplicationBin"),all.x = TRUE)
    data[, ReplicationBin := 6 - ReplicationBin]
    
    data4max <- data[,.("Cnt"=sum(Cnt)),by=.(Motif)]
    maxMutCnt <- data4max[, max(Cnt,na.rm = TRUE)]
    
    plots <- list()
    
    totalCnt <- 0
    for(j in 1:nrow(motifs))
    {
      m <- as.character(motifs[j]$Motif)
      dt <- data[Motif == m]
      mname <- paste0(m," / ",crDNA(m))
      dt[, MotifName := mname]
      dt[, MutDensity := Cnt/TrgCnt]
      
      mutcnt <- dt[,sum(Cnt)]
      plotdata <- data.table("motif"=factor(mname),"MutCnt"=mutcnt)
      totalCnt <- totalCnt + mutcnt
      
      p <- ggplot(plotdata, aes(x=motif,y=MutCnt)) + geom_bar(stat="identity",fill="steelblue") + ylim(0,maxMutCnt) + coord_flip() +
        theme(panel.background = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(size=6),
              axis.text.y = element_text(size=6),
              axis.line.x = element_line(color="black"),
              axis.line.y = element_line(color="black"))
      plots[[(j-1)*3+1]] <- p
      
      p <- ggplot(dt, aes(x=ReplicationBin,y=Cnt,fill=MutateAllele)) + geom_bar(stat="identity") +
        theme(panel.background = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(size=6),
              axis.text.y = element_text(size=6),
              axis.line.x = element_line(color="black"),
              axis.line.y = element_line(color="black"),
              legend.title = element_blank(),
              legend.text = element_text(size = 6),
              legend.key.size = unit(2,"mm"))
      plots[[(j-1)*3+2]] <- p

      p <- ggplot(dt, aes(x=ReplicationBin,y=MutDensity,fill=MutateAllele)) + geom_bar(stat="identity") +
        theme(panel.background = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(size=6),
              axis.text.y = element_text(size=6),
              axis.line.x = element_line(color="black"),
              axis.line.y = element_line(color="black"),
              legend.title = element_blank(),
              legend.text = element_text(size = 6),
              legend.key.size = unit(2,"mm"))
      plots[[(j-1)*3+3]] <- p
      
      dtlm <- dt[,.("Cnt"=sum(Cnt),"TrgCnt"=max(TrgCnt)),by=.(ReplicationBin)]
      dtlm[, MutDensity := Cnt/TrgCnt]
      totalDensity <- dtlm[,sum(MutDensity)]
      dtlm[, RelativeDensity := MutDensity/totalDensity]
      
      motiflm <- lm(RelativeDensity ~ ReplicationBin, dtlm)
      
      resultsMotif <- rbind(resultsMotif,data.table("Cancer"=cancer,
                                                    "Sample"=sample,
                                                    "Motif"=mname,
                                                    "Cnt"=mutcnt,
                                                    "Coef"=motiflm$coefficients[[2]],
                                                    "Stderr"=summary(motiflm)$coefficients[2,2],
                                                    "Pvalue"=summary(motiflm)$coefficients[2,4]))
    }
    
    # APOBEC signature plots
    dataAPOBEC <- data[substr(Motif,1,3) %in% c("TCA","TCT") & MutateAllele %in% c("T","G")]
    dataAPOBECgrp <- dataAPOBEC[,.("TrgCnt"=sum(TrgCnt),"Cnt"=sum(Cnt)), by=.(ReplicationBin,MutateAllele)]
    dataAPOBECgrp[, MutDensity := Cnt/TrgCnt]
    dataAPOBECgrp[, MotifName := "APOBEC"]
 
    p <- ggplot(dataAPOBECgrp, aes(x=ReplicationBin,y=Cnt,fill=MutateAllele)) + geom_bar(stat="identity") +
      theme(panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            legend.key.size = unit(2,"mm"))
    plots[[98]] <- p
    
    p <- ggplot(dataAPOBECgrp, aes(x=ReplicationBin,y=MutDensity,fill=MutateAllele)) + geom_bar(stat="identity") +
      theme(panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            legend.key.size = unit(2,"mm"))
    plots[[99]] <- p
           
    dataAPOBEClm <- dataAPOBECgrp[,.("Cnt"=sum(Cnt),"TrgCnt"=max(TrgCnt)),by=.(ReplicationBin)]
    dataAPOBEClm[, MutDensity := Cnt/TrgCnt]
    totalDensity <- dataAPOBEClm[,sum(MutDensity)]
    dataAPOBEClm[, RelativeDensity := MutDensity/totalDensity]
    
    APOBEClm <- lm(RelativeDensity ~ ReplicationBin, dataAPOBEClm)
    
    APOBECcnt <- dataAPOBEClm[, sum(Cnt)]
    
    # *C* excuding APOBEC plot
    dataXCX <- data[Motif %in% c("ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCC","TCG")]
    dataXCXgrp <- dataXCX[,.("TrgCnt"=sum(TrgCnt),"Cnt"=sum(Cnt)), by=.(ReplicationBin,MutateAllele)]
    dataXCXgrp[, MutDensity := Cnt/TrgCnt]
    dataXCXgrp[, MotifName := "XCX"]
 
    p <- ggplot(dataXCXgrp, aes(x=ReplicationBin,y=Cnt,fill=MutateAllele)) + geom_bar(stat="identity") +
      theme(panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            legend.key.size = unit(2,"mm"))
    plots[[101]] <- p
    
    p <- ggplot(dataXCXgrp, aes(x=ReplicationBin,y=MutDensity,fill=MutateAllele)) + geom_bar(stat="identity") +
      theme(panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            legend.key.size = unit(2,"mm"))
    plots[[102]] <- p
    
    dataXCXlm <- dataXCXgrp[,.("Cnt"=sum(Cnt),"TrgCnt"=max(TrgCnt)),by=.(ReplicationBin)]
    dataXCXlm[, MutDensity := Cnt/TrgCnt]
    totalDensity <- dataXCXlm[,sum(MutDensity)]
    dataXCXlm[, RelativeDensity := MutDensity/totalDensity]
    
    XCXlm <- lm(RelativeDensity ~ ReplicationBin, dataXCXlm)
    
    XCXcnt <- dataXCXlm[, sum(Cnt)]
    
    pl <- ggarrange(plotlist=plots, ncol=3, nrow=10)
    ggexport(pl, filename=paste0(OUTPUT_DIR,"/",cancer,"/",enrichChar,sample,".pdf"))
    
    results <- rbind(results,data.table("Cancer"=cancer,
                                        "Sample"=sample,
                                        "APOBECcnt"=APOBECcnt,
                                        "APOBECcoef"=APOBEClm$coefficients[[2]],
                                        "APOBECstderr"=summary(APOBEClm)$coefficients[2,2],
                                        "APOBECpvalue"=summary(APOBEClm)$coefficients[2,4],
                                        "XCXcnt"=XCXcnt,
                                        "XCXcoef"=XCXlm$coefficients[[2]],
                                        "XCXstderr"=summary(XCXlm)$coefficients[2,2],
                                        "XCXpvalue"=summary(XCXlm)$coefficients[2,4],
                                        "TotalCnt"=totalCnt))
    
    resultsMotif <- rbind(resultsMotif,data.table("Cancer"=cancer,
                                                  "Sample"=sample,
                                                  "Motif"="APOBEC",
                                                  "Cnt"=APOBECcnt,
                                                  "Coef"=APOBEClm$coefficients[[2]],
                                                  "Stderr"=summary(APOBEClm)$coefficients[2,2],
                                                  "Pvalue"=summary(APOBEClm)$coefficients[2,4]))
    
    resultsMotif <- rbind(resultsMotif,data.table("Cancer"=cancer,
                                                  "Sample"=sample,
                                                  "Motif"="XCX",
                                                  "Cnt"=XCXcnt,
                                                  "Coef"=XCXlm$coefficients[[2]],
                                                  "Stderr"=summary(XCXlm)$coefficients[2,2],
                                                  "Pvalue"=summary(XCXlm)$coefficients[2,4]))
    
  }
  
}

write.csv(results,paste0(OUTPUT_DIR,"/coefs.csv"))
write.csv(resultsMotif,paste0(OUTPUT_DIR,"/coefsMotif.csv"))

