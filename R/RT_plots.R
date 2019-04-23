library(data.table)
library(ggplot2)
library(reshape2)

ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR/RTplots"

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR/RT/coefs.csv", sep = ',',header = TRUE)
data <- data.table(data)
data[, SlopeDifference := ApobecSlope - OtherSlope]

cancers <- data.table("cancer"=unique(data[,Cancer]))

dir.create(file.path(ROOT_DIR))
for(i in 1:nrow(cancers))
{
  dir.create(file.path(ROOT_DIR, cancers[i]$cancer))
} 

for(i in 1:nrow(cancers))
{
  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=ApobecSlope)) + 
    geom_point() #+
    #coord_cartesian(ylim = c(-10e-6, 10e-6)) 
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTApobecSlope.jpg"), plot = p, device = "jpeg")

  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=OtherSlope)) + 
    geom_point() #+
    #coord_cartesian(ylim = c(-10e-6, 10e-6)) 

  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTOtherSlope.jpg"), plot = p, device = "jpeg")
  
  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=LaggingSlope)) + 
    geom_point()
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTLaggingSlope.jpg"), plot = p, device = "jpeg")
  
  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=LeadingSlope)) + 
    geom_point()
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTLeadingSlope.jpg"), plot = p, device = "jpeg")

  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=StrandRatioSlope)) + 
    geom_point()
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTStrandRatioSlope.jpg"), plot = p, device = "jpeg")

  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=SlopeDifference)) + 
    geom_point() #+
    #coord_cartesian(ylim = c(-10e-6, 10e-6)) 
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTSlopeDifference.jpg"), plot = p, device = "jpeg")
  
}
