library(data.table)
library(ggplot2)
library(reshape2)

ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR/EXPplots"

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsR/EXP/coefs.csv", sep = ',',header = TRUE)
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
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPApobecSlope.jpg"), plot = p, device = "jpeg")

  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=OtherSlope)) + 
    geom_point() #+
    #coord_cartesian(ylim = c(-10e-6, 10e-6)) 

  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPOtherSlope.jpg"), plot = p, device = "jpeg")
  
  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=PlusConsistSlope)) + 
    geom_point()
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPPlusConsistSlope.jpg"), plot = p, device = "jpeg")
  
  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=MinusConsistSlope)) + 
    geom_point()
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPMinusConsistSlope.jpg"), plot = p, device = "jpeg")

  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=PlusAllSlope)) + 
    geom_point()
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPPlusAllSlope.jpg"), plot = p, device = "jpeg")
  
  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=MinusAllSlope)) + 
    geom_point()
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPMinusAllSlope.jpg"), plot = p, device = "jpeg")
  
  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=StrandRatioSlope)) + 
    geom_point()
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPStrandRatioSlope.jpg"), plot = p, device = "jpeg")

  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=SlopeDifference)) + 
    geom_point() #+
    #coord_cartesian(ylim = c(-10e-6, 10e-6)) 
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPSlopeDifference.jpg"), plot = p, device = "jpeg")
  
}
