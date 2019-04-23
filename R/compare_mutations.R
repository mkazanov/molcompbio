library(data.table)

cpp <- read.csv("/Users/mar/BIO/BIODATA/CancerMutations/Fredriksson_et_al_2014/mutations_apobec.tsv", sep = '\t', header = FALSE)
almira <- read.csv("/Users/mar/Dropbox/Almira/mutations/BLCA_short.txt", sep = '\t')

cpp <- data.table(cpp)
almira <- data.table(almira)

cpp <- cpp[V1 == "BLCA"]
almira <- almira[APOBEC_mutation == 1]

cpp[, V3 := paste0("chr",V3)]

cppKey <- cpp[,.(V3,V4)]
almiraKey <- almira[,.(Chromosome, Start_position)]

setnames(cppKey,c("chr","pos"))
setnames(almiraKey,c("chr","pos"))
