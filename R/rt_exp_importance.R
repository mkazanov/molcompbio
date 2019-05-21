library(data.table)

muts <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/results_RTexp_APOBEC.txt", sep='\t')
muts <- data.table(muts)

muts[ReplicationBin == 6, ReplicationBin:=7]
muts[ReplicationBin == 5, ReplicationBin:=6]
muts[ReplicationBin == 4, ReplicationBin:=5]
muts[ReplicationBin == 3, ReplicationBin:=4]
muts[ReplicationBin == 2, ReplicationBin:=3]
muts[ReplicationBin == 1, ReplicationBin:=2]
muts[ReplicationBin == 0, ReplicationBin:=1]
muts[ReplicationBin < 0, ReplicationBin:=0]

muts[ExpressionBin < 0, ExpressionBin:=0]

mt <- muts[,.("MutationCnt"=sum(MutationCnt)), by = .(Cancer,Sample,ReplicationBin,ExpressionBin)]

targets <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/TCW_in_RTEXPbins_cluster.csv")
targets <- data.table(targets)

targets[ReplicationBin == 6, ReplicationBin:=7]
targets[ReplicationBin == 5, ReplicationBin:=6]
targets[ReplicationBin == 4, ReplicationBin:=5]
targets[ReplicationBin == 3, ReplicationBin:=4]
targets[ReplicationBin == 2, ReplicationBin:=3]
targets[ReplicationBin == 1, ReplicationBin:=2]
targets[ReplicationBin == 0, ReplicationBin:=1]
targets[ReplicationBin < 0, ReplicationBin:=0]

targets[ExpressionBin < 0, ExpressionBin:=0]

trgs <- targets[,.("TargetCnt"=sum(MutationCnt)), by = .(Cancer,Sample,ReplicationBin,ExpressionBin)]

data <- merge(mt, trgs, by = c("Cancer","Sample","ReplicationBin","ExpressionBin"))
data[, SampleMutTotal := sum(MutationCnt), by = c("Cancer","Sample")]
data[, y:=as.numeric(MutationCnt/SampleMutTotal)]
data[, y:=as.numeric(y/TargetCnt)]

enrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/enrichment.txt", sep = '\t')
enrichment <- data.table(enrichment)

gordeninEnrichment <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results/sample_enrichment.txt", sep = '\t')
gordeninEnrichment <- data.table(gordeninEnrichment)

data <- merge(data, gordeninEnrichment, by.x = c("Cancer","Sample"), by.y = c("CANCER_TYPE","SAMPLE"))
data <- merge(data, enrichment, by.x = c("Cancer","Sample"), by.y = c("cancer","sample"))

dt <- data[APOBEC_ENRICHMENT > 2.0]


null <- lm(y ~ 1, data)
full <- lm(y ~ .^2, data)
step(null, scope = list(lower = null, upper = full), direction = 'forward') 











