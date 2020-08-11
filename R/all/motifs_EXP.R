ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL/"
MOTIFS_LIST <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/all/all_motifs.txt"
OUTPUT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalRall/EXPmotifs"

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
  for(rt in c(0,1,2,3,4,5,6,7)){
    trgInitialTemplate <- rbind(trgInitialTemplate, data.table("Motif"=m,"ExpressionBin"=rt))
    for(nt in c("A","C","G","T")){
      if(substr(m,2,2) != nt){
        mutInitialTemplate <- rbind(mutInitialTemplate, data.table("Motif"=m,
                                                     "ExpressionBin"=rt, 
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
    
    muts <- muts[ExpressionBin != -1]
    trgs <- trgs[ExpressionBin != -1]
    muts <- muts[ExpressionBin != 0]
    trgs <- trgs[ExpressionBin != 0]
    muts[ExpressionBin == -2, ExpressionBin := 0]
    trgs[ExpressionBin == -2, ExpressionBin := 0]
    
    mutsG <- muts[,.("Cnt"=sum(Cnt)), by = .(Motif,ExpressionBin,MutateAllele)]
    trgsG <- trgs[,.("TrgCnt"=sum(Cnt)), by = .(Motif,ExpressionBin)]
    
    mutTemplate <- copy(mutInitialTemplate)
    trgTemplate <- copy(trgInitialTemplate)
    
    mutsG <- merge(mutTemplate,mutsG,by=c("Motif","ExpressionBin","MutateAllele"), all.x=TRUE)
    mutsG[is.na(Cnt), Cnt:= 0] 
    trgsG <- merge(trgTemplate,trgsG,by=c("Motif","ExpressionBin"), all.x=TRUE)
    
    data <- merge(trgsG,mutsG,by=c("Motif","ExpressionBin"),all.x = TRUE)

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
            
      p <- ggplot(plotdata, aes(x=motif,y=MutCnt)) + 
        geom_segment(aes(x=motif,xend=motif,y=0,yend=MutCnt), color="skyblue") + 
        geom_point(color="blue",size=4,alpha=0.6) +
        scale_y_continuous(expand = c(0, 0), limits=c(0,maxMutCnt+(0.1*maxMutCnt))) +
        ylab("SBS count") +
        coord_flip() +
        theme(panel.background = element_blank(),
              axis.title.x = element_text(size=5),
              axis.title.y = element_blank(),
              axis.text.x = element_text(size=5),
              axis.text.y = element_text(size=6),
              axis.line.x = element_line(color="black"),
              axis.line.y = element_line(color="black"))
      plots[[(j-1)*3+1]] <- p
      
      p <- ggplot(dt, aes(x=ExpressionBin,y=Cnt,fill=MutateAllele)) + geom_bar(stat="identity") +
        xlab("Low                    Gene expression level                    High") +
        ylab("SBS count") +
        theme(panel.background = element_blank(),
              axis.title.x = element_text(size=5),
              axis.title.y = element_text(size=5),
              axis.text.x = element_text(size=5),
              axis.text.y = element_text(size=5),
              axis.line.x = element_line(color="black"),
              axis.line.y = element_line(color="black"),
              legend.title = element_blank(),
              legend.text = element_text(size = 6),
              legend.key.size = unit(2,"mm"))
      plots[[(j-1)*3+2]] <- p
      
      p <- ggplot(dt, aes(x=ExpressionBin,y=MutDensity,fill=MutateAllele)) + geom_bar(stat="identity") +
        xlab("Low                    Gene expression level                    High") +
        ylab("SBS density") +
        theme(panel.background = element_blank(),
              axis.title.x = element_text(size=5),
              axis.title.y = element_text(size=5),
              axis.text.x = element_text(size=5),
              axis.text.y = element_text(size=5),
              axis.line.x = element_line(color="black"),
              axis.line.y = element_line(color="black"),
              legend.title = element_blank(),
              legend.text = element_text(size = 6),
              legend.key.size = unit(2,"mm"))
      plots[[(j-1)*3+3]] <- p
    
      dtlm <- dt[,.("Cnt"=sum(Cnt),"TrgCnt"=max(TrgCnt)),by=.(ExpressionBin)]
      dtlm[, MutDensity := Cnt/TrgCnt]
      totalDensity <- dtlm[,sum(MutDensity)]
      dtlm[, RelativeDensity := MutDensity/totalDensity]
      
      motiflm <- lm(RelativeDensity ~ ExpressionBin, dtlm)
      
      resultsMotif <- rbind(resultsMotif,data.table("Cancer"=cancer,
                                                    "Sample"=sample,
                                                    "Motif"=mname,
                                                    "Cnt"=mutcnt,
                                                    "Coef"=motiflm$coefficients[[2]],
                                                    "Stderr"=summary(motiflm)$coefficients[2,2],
                                                    "Pvalue"=summary(motiflm)$coefficients[2,4]))
      
    }
    
    # APOBEC signature plots
    dataAPOBEC <- data[substr(Motif,1,2) == "TC" & MutateAllele %in% c("T","G")]
    dataAPOBECgrp <- dataAPOBEC[,.("TrgCnt"=sum(TrgCnt),"Cnt"=sum(Cnt)), by=.(ExpressionBin,MutateAllele)]
    dataAPOBECgrp[, MutDensity := Cnt/TrgCnt]
    dataAPOBECgrp[, MotifName := "APOBEC"]

        sumApobecCnt <- dataAPOBECgrp[, sum(Cnt)]
        p <- ggplot(dataAPOBECgrp, aes(x="TCX / XGA",y=sumApobecCnt)) + 
          geom_segment(aes(x="TCX / XGA",xend="TCX / XGA",y=0,yend=sumApobecCnt), color="skyblue") + 
          geom_point(color="blue",size=4,alpha=0.6) +
          scale_y_continuous(expand = c(0, 0), limits=c(0,sumApobecCnt+(0.1*sumApobecCnt))) +
          ylab("APOBEC-induced SBS count") +
          coord_flip() +
          theme(panel.background = element_blank(),
                axis.title.x = element_text(size=5),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size=5),
                axis.text.y = element_text(size=6),
                axis.line.x = element_line(color="black"),
                axis.line.y = element_line(color="black"))
        plots[[97]] <- as.grob(p)
        
      p <- ggplot(dataAPOBECgrp, aes(x=ExpressionBin,y=Cnt,fill=MutateAllele)) + geom_bar(stat="identity") +
      xlab("Low                    Gene expression level                    High") +
      ylab("SBS count") +
      theme(panel.background = element_blank(),
            axis.title.x = element_text(size=5),
            axis.title.y = element_text(size=5),
            axis.text.x = element_text(size=5),
            axis.text.y = element_text(size=5),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            legend.key.size = unit(2,"mm"))
    plots[[98]] <- p

    p <- ggplot(dataAPOBECgrp, aes(x=ExpressionBin,y=MutDensity,fill=MutateAllele)) + geom_bar(stat="identity") +
      xlab("Low                    Gene expression level                    High") +
      ylab("SBS density") +
      theme(panel.background = element_blank(),
            axis.title.x = element_text(size=5),
            axis.title.y = element_text(size=5),
            axis.text.x = element_text(size=5),
            axis.text.y = element_text(size=5),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            legend.key.size = unit(2,"mm"))
    plots[[99]] <- p
    
    dataAPOBEClm <- dataAPOBECgrp[,.("Cnt"=sum(Cnt),"TrgCnt"=max(TrgCnt)),by=.(ExpressionBin)]
    dataAPOBEClm[, MutDensity := Cnt/TrgCnt]
    totalDensity <- dataAPOBEClm[,sum(MutDensity)]
    dataAPOBEClm[, RelativeDensity := MutDensity/totalDensity]
    
    APOBEClm <- lm(RelativeDensity ~ ExpressionBin, dataAPOBEClm)
    
    APOBECcnt <- dataAPOBEClm[, sum(Cnt)]
    
    
    # *C* excuding APOBEC plot
    dataXCX <- data[Motif %in% c("ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT")]
    dataXCXgrp <- dataXCX[,.("TrgCnt"=sum(TrgCnt),"Cnt"=sum(Cnt)), by=.(ExpressionBin,MutateAllele)]
    dataXCXgrp[, MutDensity := Cnt/TrgCnt]
    dataXCXgrp[, MotifName := "XCX"]
    
    sumXCXCnt <- dataXCXgrp[, sum(Cnt)]
    p <- ggplot(dataXCXgrp, aes(x="XCX / XGX",y=sumXCXCnt)) + 
      geom_segment(aes(x="XCX / XGX",xend="XCX / XGX",y=0,yend=sumXCXCnt), color="skyblue") + 
      geom_point(color="blue",size=4,alpha=0.6) +
      scale_y_continuous(expand = c(0, 0), limits=c(0,sumXCXCnt+(0.1*sumXCXCnt))) +
      ylab("non-APOBEC-induced SBS in cytosines count") +
      coord_flip() +
      theme(panel.background = element_blank(),
            axis.title.x = element_text(size=5),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=5),
            axis.text.y = element_text(size=6),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"))
    plots[[100]] <- as.grob(p)
    
    p <- ggplot(dataXCXgrp, aes(x=ExpressionBin,y=Cnt,fill=MutateAllele)) + geom_bar(stat="identity") +
      xlab("Low                    Gene expression level                    High") +
      ylab("SBS count") +
      theme(panel.background = element_blank(),
            axis.title.x = element_text(size=5),
            axis.title.y = element_text(size=5),
            axis.text.x = element_text(size=5),
            axis.text.y = element_text(size=5),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            legend.key.size = unit(2,"mm"))
    plots[[101]] <- p
    
    p <- ggplot(dataXCXgrp, aes(x=ExpressionBin,y=MutDensity,fill=MutateAllele)) + geom_bar(stat="identity") +
      xlab("Low                    Gene expression level                    High") +
      ylab("SBS density") +
      theme(panel.background = element_blank(),
            axis.title.x = element_text(size=5),
            axis.title.y = element_text(size=5),
            axis.text.x = element_text(size=5),
            axis.text.y = element_text(size=5),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            legend.key.size = unit(2,"mm"))
    plots[[102]] <- p
    
    dataXCXlm <- dataXCXgrp[,.("Cnt"=sum(Cnt),"TrgCnt"=max(TrgCnt)),by=.(ExpressionBin)]
    dataXCXlm[, MutDensity := Cnt/TrgCnt]
    totalDensity <- dataXCXlm[,sum(MutDensity)]
    dataXCXlm[, RelativeDensity := MutDensity/totalDensity]
    
    XCXlm <- lm(RelativeDensity ~ ExpressionBin, dataXCXlm)
    
    XCXcnt <- dataXCXlm[, sum(Cnt)]
    
    
    dtSampleMotifs <- resultsMotif[Sample == sample & Cancer == cancer]
    
    p <- ggplot(dtSampleMotifs, aes(x=Motif,y=Coef)) + geom_bar(stat="identity",
                                                                colour=rgb(38,120,178,maxColorValue = 255),
                                                                fill=rgb(146,187,216,maxColorValue = 255)) +
      ylab("Slope of mutation density distribution over gene expression levels") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5),
            #axis.text.x = element_blank(),
            axis.text.y = element_text(size=5),
            axis.ticks.x=element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=5))
    
    plots[[103]] <- p
  
    
    laymat <- c()
    for(i in 1:11){
      laymat <- rbind(laymat,c((i-1)*3+1,(i-1)*3+2,(i-1)*3+3))
    }
    
    pl <- marrangeGrob(grobs=plots[1:96], nrow=11, ncol=3, layout_matrix = laymat, 
                       top=paste0("Cancer: ",cancer,", Sample: ",sample,", APOBEC enrichment: ",round(aenrich,2)))
    ggexport(pl, filename=paste0(OUTPUT_DIR,"/",cancer,"/",enrichChar,sample,"_1.pdf"))
    
    laymat <- c()
    laymat <- rbind(c(97,98,99),c(100,101,102),c(103,103,103),c(103,103,103),c(104,105,106))  
    pl <- marrangeGrob(grobs=plots[97:103], nrow=11, ncol=3, layout_matrix = laymat, 
                       top=paste0("Cancer: ",cancer,", Sample: ",sample,", APOBEC enrichment: ",round(aenrich,2)))
    ggexport(pl, filename=paste0(OUTPUT_DIR,"/",cancer,"/",enrichChar,sample,"_2.pdf"))
    
    
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

