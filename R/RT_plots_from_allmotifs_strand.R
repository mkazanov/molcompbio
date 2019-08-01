library(data.table)
library(ggplot2)
library(reshape2)

ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalR/RTplotsStrand"

data <- read.csv(paste0(ROOT_DIR,"/coefsFinalStrand.csv"), sep = ',',header = TRUE)
data <- data.table(data)

cancers <- data.table("cancer"=unique(data[,Cancer]))

dir.create(file.path(ROOT_DIR))
for(i in 1:nrow(cancers))
{
  dir.create(file.path(ROOT_DIR, cancers[i]$cancer))
} 

dtranges <- data.table()
dts <- data.table()

for(i in 1:nrow(cancers))
{
  
  statSmoothVals <- data.table("GordeninEnrichment"=seq(0.85,3.25,0.01))
  
  dt <- data[Cancer == cancers[i]$cancer]
  
  dtcur <- dt[RatioDiff < quantile(RatioDiff,probs=0.95)]
  dtcur <- dtcur[Leading > 20 & Lagging > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrand.jpg"), plot = p, device = "jpeg")

  dtcur <- dt[RatioDiff0 < quantile(RatioDiff0,probs=0.95)]
  dtcur <- dtcur[Leading0 > 20 & Lagging0 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff0)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 

  model <- loess(RatioDiff0 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff0"=pred[["fit"]]))

  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin0.jpg"), plot = p, device = "jpeg")

  dtcur <- dt[RatioDiff1 < quantile(RatioDiff1,probs=0.95)]
  dtcur <- dtcur[Leading1 > 20 & Lagging1 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff1)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff1 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff1"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin1.jpg"), plot = p, device = "jpeg")

  dtcur <- dt[RatioDiff2 < quantile(RatioDiff2,probs=0.95)]
  dtcur <- dtcur[Leading2 > 20 & Lagging2 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff2)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff2 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff2"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin2.jpg"), plot = p, device = "jpeg")
  
  dtcur <- dt[RatioDiff3 < quantile(RatioDiff3,probs=0.95)]
  dtcur <- dtcur[Leading3 > 20 & Lagging3 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff3)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff3 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff3"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin3.jpg"), plot = p, device = "jpeg")

  dtcur <- dt[RatioDiff4 < quantile(RatioDiff4,probs=0.95)]
  dtcur <- dtcur[Leading4 > 20 & Lagging4 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff4)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff4 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff4"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin4.jpg"), plot = p, device = "jpeg")
  
  dtcur <- dt[RatioDiff5 < quantile(RatioDiff5,probs=0.95)]
  dtcur <- dtcur[Leading5 > 20 & Lagging5 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff5)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff5 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff5"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin5.jpg"), plot = p, device = "jpeg")

  dtcur <- dt[RatioDiff6 < quantile(RatioDiff6,probs=0.95)]
  dtcur <- dtcur[Leading6 > 20 & Lagging6 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff6)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 

  model <- loess(RatioDiff6 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff6"=pred[["fit"]]))
    
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin6.jpg"), plot = p, device = "jpeg")
  
  dt1 <- dt[,.(GordeninEnrichment,RatioDiff0,RatioDiff1,RatioDiff2,RatioDiff3,RatioDiff4,RatioDiff5,RatioDiff6)]
  dt2 <- melt(dt1,id.vars=c("GordeninEnrichment"))
  
  p <- ggplot(dt2, aes(x=GordeninEnrichment, y=value, color=variable)) + 
    geom_point() #+
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBins.jpg"), plot = p, device = "jpeg")
  
  statMelt <- melt(statSmoothVals,id.vars = "GordeninEnrichment")
  p <- ggplot(statMelt,aes(x=GordeninEnrichment, y=variable)) + geom_tile(aes(fill = (value-1))) + 
      scale_fill_gradient2(low="darkblue", high="darkgreen")
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandColor.jpg"), plot = p, device = "jpeg")
 
#  dt[,EnrichmentRange := cut(GordeninEnrichment,
#                             breaks=c(0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,Inf),
#                             labels=c("0.75-1","1-1.25","1.25-1.5","1.5-1.75","1.75-2","2-2.25","2.25-2.5","2.5-2.75","2.75-3","3-3.25","3.25-3.5"))]
    dt[,EnrichmentRange := cut(GordeninEnrichment,
                               breaks=c(0,1,1.5,2,2.5,3,Inf),
                               labels=c("0-1","<1.5","1.5-2","2-2.5","2.5-3",">3"))]
  
    dt <- dt[!(EnrichmentRange == "0-1")]
    
  dtrange <- dt[,.("cnt"=.N,
         "RatioDiff"=mean(RatioDiff),
        "RatioDiff0"=mean(RatioDiff0),
        "RatioDiff1"=mean(RatioDiff1),
        "RatioDiff2"=mean(RatioDiff2),
        "RatioDiff3"=mean(RatioDiff3),
        "RatioDiff4"=mean(RatioDiff4),
        "RatioDiff5"=mean(RatioDiff5),
        "RatioDiff6"=mean(RatioDiff6)), by=EnrichmentRange]
  
  dtrange <- dtrange[,RatioDiff := NULL]
  
  statMelt <- melt(dtrange,id.vars = c("EnrichmentRange","cnt"))
  p <- ggplot(statMelt,aes(x=EnrichmentRange, y=variable)) + geom_tile(aes(fill = log(value,2))) + 
    scale_fill_gradient2(low="#dd1c77", high="#3182bd",limits=c(-1,1)) +
    xlab("APOBEC enrichment") +
    ylab("Replication timing") +
    scale_y_discrete(labels = c("Late","","","","","","Early")) +
    labs(fill = "") +
   # labs(fill = "Lagging/leading strand mutations log ratio")
    theme_bw() +
    theme(axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandColorRange.jpg"), plot = p, device = "jpeg")
  
  dtrange[, cancer := cancers[i]$cancer]
  dtranges <- rbind(dtranges,dtrange)
  
  dt[, cancer := cancers[i]$cancer]
  dts <- rbind(dts,dt)
}
