library(data.table)

allmuts <- read.csv("/Users/mar/BIO/BIODATA/CancerMutations/Fredriksson_et_al_2014/mutations.tsv",sep='\t')
allmuts <- data.table(allmuts)
muts <- allmuts[cancer %in% c("BLCA","BRCA","HNSC","LUAD","LUSC")]

m <- muts[,.N,by=.(cancer,barcode)]
