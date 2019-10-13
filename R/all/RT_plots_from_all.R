library(data.table)
library(ggplot2)
library(reshape2)

ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalRall/RTplots"

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalRall/RT/coefsFinal.csv", sep = ',',header = TRUE)
data <- data.table(data)
data[, ApobecSlope := coef]
data[, ymin := coef - stderr]
data[, ymax := coef + stderr]
#data[, OtherSlope := meanС2OtherSlope]
data[, OtherSlope := coefOther]
data[, SlopeDifference := ApobecSlope - OtherSlope]
data[, ymindiff := SlopeDifference - stderr - stderrOther]
data[, ymaxdiff := SlopeDifference + stderr + stderrOther]

data <- data[(stderr + stderrOther) < 0.01]

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
    geom_point() +
    geom_errorbar(aes(ymin=ymin,ymax=ymax))
    #+
    #coord_cartesian(ylim = c(-10e-6, 10e-6)) 
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTApobecSlope.jpg"), plot = p, device = "jpeg")

  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=OtherSlope)) + 
    geom_point() #+
    #coord_cartesian(ylim = c(-10e-6, 10e-6)) 

  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTOtherSlope.jpg"), plot = p, device = "jpeg")
  
  dt <- data[Cancer == cancers[i]$cancer]
  p <- ggplot(dt, aes(x=GordeninEnrichment, y=SlopeDifference)) + 
    geom_point() +
    geom_errorbar(aes(ymin=ymindiff,ymax=ymaxdiff))
    #+
    #coord_cartesian(ylim = c(-10e-6, 10e-6)) 
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTSlopeDifference.jpg"), plot = p, device = "jpeg")
  
}
