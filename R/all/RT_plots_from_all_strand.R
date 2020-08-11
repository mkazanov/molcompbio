library(data.table)
library(ggplot2)
library(reshape2)

INPUT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalRall/RTStrand"
ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalRall/RTplotsStrand"

data <- read.csv(paste0(INPUT_DIR,"/coefsFinalStrand.csv"), sep = ',',header = TRUE)
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
  
  # ratioDiff and ratio  
  
  dtcur <- dt[Leading > 20 & Lagging > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrand.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffRTstrand.jpg"), plot = p, device = "jpeg")

  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densRTstrand.jpg"), plot = p, device = "jpeg")
  
    
  # ratio0 and ratio diff0
  
  dtcur <- dt[Leading0 > 20 & Lagging0 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff0)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 

  model <- loess(RatioDiff0 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff0"=pred[["fit"]]))

  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin0.jpg"), plot = p, device = "jpeg")

  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio0)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio0 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio0"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffRTstrandBin0.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity0)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity0 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity0"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densRTstrandBin0.jpg"), plot = p, device = "jpeg")
  
  
    
  # ratio1 and ratio diff1
  
  dtcur <- dt[Leading1 > 20 & Lagging1 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff1)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff1 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff1"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin1.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio1)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio1 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio1"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffRTstrandBin1.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity1)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity1 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity1"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densRTstrandBin1.jpg"), plot = p, device = "jpeg")
  
  # ratio2 and ratio diff2

  dtcur <- dt[Leading2 > 20 & Lagging2 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff2)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff2 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff2"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin2.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio2)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio2 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio2"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffRTstrandBin2.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity2)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity2 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity2"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densRTstrandBin2.jpg"), plot = p, device = "jpeg")
  
  # ratio3 and ratio diff3
  
  dtcur <- dt[Leading3 > 20 & Lagging3 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff3)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff3 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff3"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin3.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio3)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio3 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio3"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffRTstrandBin3.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity3)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity3 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity3"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densRTstrandBin3.jpg"), plot = p, device = "jpeg")
  
  
  # ratio4 and ratio diff4
  
  dtcur <- dt[Leading4 > 20 & Lagging4 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff4)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff4 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff4"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin4.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio4)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio4 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio4"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffRTstrandBin4.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity4)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity4 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity4"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densRTstrandBin4.jpg"), plot = p, device = "jpeg")
  
  # ratio5 and ratio diff5
  
  dtcur <- dt[Leading5 > 20 & Lagging5 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff5)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff5 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff5"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin5.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio5)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio5 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio5"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffRTstrandBin5.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity5)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity5 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity5"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densRTstrandBin5.jpg"), plot = p, device = "jpeg")
  
  # ratio6 and ratio diff6
  
  dtcur <- dt[Leading6 > 20 & Lagging6 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff6)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 

  model <- loess(RatioDiff6 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff6"=pred[["fit"]]))
    
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_RTstrandBin6.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio6)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio6 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio6"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffRTstrandBin6.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity6)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity6 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity6"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densRTstrandBin6.jpg"), plot = p, device = "jpeg")
  
  #
    
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
                               breaks=c(0,1,1.5,2,2.5,Inf),
                               labels=c("0-1","<1.5","1.5-2","2-2.5",">2.5"))]
  
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
  p <- ggplot(statMelt,aes(x=EnrichmentRange, y=variable)) + 
    geom_tile(aes(fill = log(value,2))) + 
    geom_text(aes(label = round(log(value,2),3)), size = 2) +
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
  
  dtrange <- dt[,.("cnt"=.N,
                   "Ratio0"=mean(Ratio0),
                   "Ratio1"=mean(Ratio1),
                   "Ratio2"=mean(Ratio2),
                   "Ratio3"=mean(Ratio3),
                   "Ratio4"=mean(Ratio4),
                   "Ratio5"=mean(Ratio5),
                   "Ratio6"=mean(Ratio6)), by=EnrichmentRange]
  
  statMelt <- melt(dtrange,id.vars = c("EnrichmentRange","cnt"))
  
  p <- ggplot(statMelt,aes(x=EnrichmentRange, y=variable)) + 
    geom_tile(aes(fill = log(value,2))) + 
    geom_text(aes(label = round(log(value,2),3)), size = 2) +
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
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffRTstrandColorRange.jpg"), plot = p, device = "jpeg")
  
  dtrange <- dt[,.("cnt"=.N,
                   "RatioDensity0"=mean(RatioDensity0),
                   "RatioDensity1"=mean(RatioDensity1),
                   "RatioDensity2"=mean(RatioDensity2),
                   "RatioDensity3"=mean(RatioDensity3),
                   "RatioDensity4"=mean(RatioDensity4),
                   "RatioDensity5"=mean(RatioDensity5),
                   "RatioDensity6"=mean(RatioDensity6)), by=EnrichmentRange]
  
  statMelt <- melt(dtrange,id.vars = c("EnrichmentRange","cnt"))
  
  p <- ggplot(statMelt,aes(x=EnrichmentRange, y=variable)) + 
    geom_tile(aes(fill = log(value,2))) + 
    geom_text(aes(label = round(log(value,2),3)), size = 2) +
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
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densRTstrandColorRange.jpg"), plot = p, device = "jpeg")
  
  
  
  #dtrange[, cancer := cancers[i]$cancer]
  #dtranges <- rbind(dtranges,dtrange)
  
  #dt[, cancer := cancers[i]$cancer]
  #dts <- rbind(dts,dt)
}
