library(beeswarm)
library(data.table)
library(ggplot2)


gordeninEnrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results/sample_enrichment.txt", sep = '\t')
gordeninEnrichment <- data.table(gordeninEnrichment)
gordeninEnrichment[, color := ifelse(APOBEC_ENRICHMENT>=2, 
                                     rgb(212,42,47,maxColorValue = 255), 
                                     rgb(51,159,52,maxColorValue = 255))]
gordeninEnrichment[, fcolor := ifelse(APOBEC_ENRICHMENT>=2, 
                                      rgb(233,148,151,maxColorValue = 255), 
                                      rgb(152,206,153,maxColorValue = 255))]

tiff("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/paper/pics/enrichment_beeswarm.tiff")
beeswarm(APOBEC_ENRICHMENT ~ CANCER_TYPE, data=gordeninEnrichment, pch=21,
         pwcol=gordeninEnrichment[,color], pwbg=gordeninEnrichment[,fcolor], spacing=0.55, cex=1.5,
         xlab = "Cancer", ylab = "APOBEC-mutagenesis enrichment")
dev.off()

gordeninEnrichment[, grp := ifelse(APOBEC_ENRICHMENT>=2, 
                                     "enriched", 
                                     "notenriched")]

pctdata <- gordeninEnrichment[, .(cnt=.N), by=.(CANCER_TYPE,grp)]
pctdata[, total:=sum(cnt), by=CANCER_TYPE]
pctdata[, pct := round(cnt/total,2)]
