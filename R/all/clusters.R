library(regioneR)
library(karyoploteR)

gdBLCA <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/mutations/my/2014_Fredriksson_BLCA_21_WGS_mutations_adjusted_anz1_NOrepeats_sorted_anz5.txt",sep='\t')
gdBLCA <- data.table(gdBLCA)
gdBRCA <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/mutations/my/2014_Fredriksson_BRCA_96_WGS_mutations_adjusted_anz1_NOrepeats_sorted_anz5.txt",sep='\t')
gdBRCA <- data.table(gdBRCA) 
gdHNSC <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/mutations/my/2014_Fredriksson_HNSC_27_WGS_mutations_adjusted_anz1_NOrepeats_sorted_anz5.txt",sep='\t')
gdHNSC <- data.table(gdHNSC)
gdLUAD <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/mutations/my/2014_Fredriksson_LUAD_46_WGS_mutations_adjusted_anz1_NOrepeats_sorted_anz5.txt",sep='\t')
gdLUAD <- data.table(gdLUAD)
gdLUSC <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/mutations/my/2014_Fredriksson_LUSC_45_WGS_mutations_adjusted_anz1_NOrepeats_sorted_anz5.txt",sep='\t')
gdLUSC <- data.table(gdLUSC)

clusters <- data.table()
dt <- gdBLCA[!is.na(StrainCluster_ID), .(Tumor_Sample_Barcode,cancer,Chromosome,Start_position,Reference_Allele,Tumor_Seq_Allele2,APOBEC_mutation,
                                         Variant_Type,Complex_ID,Complex_Size,StrainCluster_ID,Dataset_Cluster_ID,Cluster_Size_Mutations,Cluster_Size_Complexes)]
clusters <- rbind(clusters,dt)
dt <- gdBRCA[!is.na(StrainCluster_ID), .(Tumor_Sample_Barcode,cancer,Chromosome,Start_position,Reference_Allele,Tumor_Seq_Allele2,APOBEC_mutation,
                                         Variant_Type,Complex_ID,Complex_Size,StrainCluster_ID,Dataset_Cluster_ID,Cluster_Size_Mutations,Cluster_Size_Complexes)]
clusters <- rbind(clusters,dt)
dt <- gdHNSC[!is.na(StrainCluster_ID), .(Tumor_Sample_Barcode,cancer,Chromosome,Start_position,Reference_Allele,Tumor_Seq_Allele2,APOBEC_mutation,
                                         Variant_Type,Complex_ID,Complex_Size,StrainCluster_ID,Dataset_Cluster_ID,Cluster_Size_Mutations,Cluster_Size_Complexes)]
clusters <- rbind(clusters,dt)
dt <- gdLUAD[!is.na(StrainCluster_ID), .(Tumor_Sample_Barcode,cancer,Chromosome,Start_position,Reference_Allele,Tumor_Seq_Allele2,APOBEC_mutation,
                                         Variant_Type,Complex_ID,Complex_Size,StrainCluster_ID,Dataset_Cluster_ID,Cluster_Size_Mutations,Cluster_Size_Complexes)]
clusters <- rbind(clusters,dt)
dt <- gdLUSC[!is.na(StrainCluster_ID), .(Tumor_Sample_Barcode,cancer,Chromosome,Start_position,Reference_Allele,Tumor_Seq_Allele2,APOBEC_mutation,
                                         Variant_Type,Complex_ID,Complex_Size,StrainCluster_ID,Dataset_Cluster_ID,Cluster_Size_Mutations,Cluster_Size_Complexes)]
clusters <- rbind(clusters,dt)

write.csv(clusters,"/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/cluster/cluster_muts.csv", row.names = FALSE)





gdLUSC <- data.table(gdLUSC)
gmuts <- gdLUSC[Tumor_Sample_Barcode == "TCGA-60-2698-01A"]

 

somatic.mutations <- read.csv("/Users/mar/BIO/BIODATA/CancerMutations/Fredriksson_et_al_2014/mutations_apobec.tsv", sep='\t', header=FALSE)
setnames(somatic.mutations,c("cancer","sample","chr","start","ref","alt","isForwardStrand"))
somatic.mutations$end <- somatic.mutations$start
somatic.mutations <- split(somatic.mutations, somatic.mutations$sample)
sm <- somatic.mutations[["TCGA-60-2698-01A"]]
sm.gr <- toGRanges(sm[,c("chr", "start", "end")])
seqlevelsStyle(sm.gr) <- "UCSC"

kp <- plotKaryotype(plot.type=4, chromosomes="chr1")
kpAddBaseNumbers(kp)
kpPlotDensity(kp, data = sm.gr, window.size = 10e5, r0=0.62, r1=0.8)
kpPlotRainfall(kp, data = sm.gr, r0=0, r1=0.6)

