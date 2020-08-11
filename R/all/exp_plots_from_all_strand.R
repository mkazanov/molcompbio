library(data.table)
library(ggplot2)
library(reshape2)

INPUT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalRall/EXPStrand"
ROOT_DIR <- "/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/ResultsFinalRall/EXPplotsStrand"

data <- read.csv(paste0(INPUT_DIR,"/coefsFinalExpStrand.csv"), sep = ',',header = TRUE)
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
  
  dtcur <- dt[Coding > 20 & Template > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPstrand.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffEXPstrand.jpg"), plot = p, device = "jpeg")
 
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densEXPstrand.jpg"), plot = p, device = "jpeg")
  
  # ratio1 and ratio diff1
  
  dtcur <- dt[Coding1 > 20 & Template1 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff1)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff1 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff1"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPstrandBin1.jpg"), plot = p, device = "jpeg")
  
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio1)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio1 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio1"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffEXPstrandBin1.jpg"), plot = p, device = "jpeg")
 
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity1)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity1 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity1"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densEXPstrandBin1.jpg"), plot = p, device = "jpeg")
  
  # ratio2 and ratio diff2
  
  dtcur <- dt[Coding1 > 20 & Template1 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff2)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff2 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff2"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPstrandBin2.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio2)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio2 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio2"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffEXPstrandBin2.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity2)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity2 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity2"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densEXPstrandBin2.jpg"), plot = p, device = "jpeg")
  
  # ratio3 and ratio diff3
  
  dtcur <- dt[Coding3 > 20 & Template3 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff3)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff3 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff3"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPstrandBin3.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio3)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio3 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio3"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffEXPstrandBin3.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity3)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity3 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity3"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densEXPstrandBin3.jpg"), plot = p, device = "jpeg")
  
  # ratio4 and ratio diff4
  
  dtcur <- dt[Coding4 > 20 & Template4 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff4)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff4 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff4"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPstrandBin4.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio4)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio4 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio4"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffEXPstrandBin4.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity4)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity4 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity4"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densEXPstrandBin4.jpg"), plot = p, device = "jpeg")
  
  # ratio5 and ratio diff5
  
  dtcur <- dt[Coding5 > 20 & Template5 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff5)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff5 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff5"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPstrandBin5.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio5)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio5 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio5"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffEXPstrandBin5.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity5)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity5 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity5"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densEXPstrandBin5.jpg"), plot = p, device = "jpeg")
  
  # ratio6 and ratio diff6
  
  dtcur <- dt[Coding6 > 20 & Template6 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff6)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff6 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff6"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPstrandBin6.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio6)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio6 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio6"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffEXPstrandBin6.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity6)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity6 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity6"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densEXPstrandBin6.jpg"), plot = p, device = "jpeg")
  
  # ratio7 and ratio diff7
  
  dtcur <- dt[Coding7 > 20 & Template7 > 20]
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDiff7)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDiff7 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDiff7"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPstrandBin7.jpg"), plot = p, device = "jpeg")
  
  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=Ratio7)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(Ratio7 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("Ratio7"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffEXPstrandBin7.jpg"), plot = p, device = "jpeg")

  p <- ggplot(dtcur, aes(x=GordeninEnrichment, y=RatioDensity7)) + 
    geom_point() + stat_smooth() +
    coord_cartesian(ylim = c(0, 1.5)) 
  
  model <- loess(RatioDensity7 ~ GordeninEnrichment, data = dtcur)
  pred <- predict(model, newdata=seq(0.85,3.25,0.01) , se=TRUE)
  statSmoothVals <- cbind(statSmoothVals,data.table("RatioDensity7"=pred[["fit"]]))
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densEXPstrandBin7.jpg"), plot = p, device = "jpeg")
    
  #
  
  dt[,EnrichmentRange := cut(GordeninEnrichment,
                             breaks=c(0,1,1.5,2,2.5,Inf),
                             labels=c("0-1","<1.5","1.5-2","2-2.5",">2.5"))]
  
  dt <- dt[!(EnrichmentRange == "0-1")]
  
  dtrange <- dt[,.("cnt"=.N,
                   "RatioDiff"=mean(RatioDiff),
                   "RatioDiff1"=mean(RatioDiff1),
                   "RatioDiff2"=mean(RatioDiff2),
                   "RatioDiff3"=mean(RatioDiff3),
                   "RatioDiff4"=mean(RatioDiff4),
                   "RatioDiff5"=mean(RatioDiff5),
                   "RatioDiff6"=mean(RatioDiff6),
                   "RatioDiff7"=mean(RatioDiff7)), by=EnrichmentRange]
  
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
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_EXPstrandColorRange.jpg"), plot = p, device = "jpeg")
  
  dtrange <- dt[,.("cnt"=.N,
                   "Ratio1"=mean(Ratio1),
                   "Ratio2"=mean(Ratio2),
                   "Ratio3"=mean(Ratio3),
                   "Ratio4"=mean(Ratio4),
                   "Ratio5"=mean(Ratio5),
                   "Ratio6"=mean(Ratio6),
                   "Ratio7"=mean(Ratio7)), by=EnrichmentRange]
  
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
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_nodiffEXPstrandColorRange.jpg"), plot = p, device = "jpeg")

  
  dtrange <- dt[,.("cnt"=.N,
                   "RatioDensity1"=mean(RatioDensity1),
                   "RatioDensity2"=mean(RatioDensity2),
                   "RatioDensity3"=mean(RatioDensity3),
                   "RatioDensity4"=mean(RatioDensity4),
                   "RatioDensity5"=mean(RatioDensity5),
                   "RatioDensity6"=mean(RatioDensity6),
                   "RatioDensity7"=mean(RatioDensity7)), by=EnrichmentRange]
  
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
  
  ggsave(paste0(ROOT_DIR,"/",cancers[i]$cancer,"/",cancers[i]$cancer,"_densEXPstrandColorRange.jpg"), plot = p, device = "jpeg")
  
    
}
