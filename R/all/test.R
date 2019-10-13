library(data.table)

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL/RTEXP_MUT_TCGA-BA-4076-01A.txt", sep='\t', header=FALSE)
data <- data.table(data)

d <- data[, sum(V7), by=.(V1,V2,V3,V6)]
d <- d[V1=="AAA"]
d <- d[V6=="C"]

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/ALL/RTEXP_TRG_TCGA-BA-4076-01A.txt", sep='\t', header=FALSE)
data <- data.table(data)

d <- data[, .(s=sum(V6)), by=.(V1,V2,V3)]
d <- d[V1=="AAA"]
d <- d[, sum(s), by=.(V1,V2)]
