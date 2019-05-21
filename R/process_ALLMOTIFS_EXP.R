library(stringi)
library(data.table)

crDNA <- function(dna)
{
  return(stri_reverse(chartr("acgtACGT","tgcaTGCA",dna)))
}

dir <- list.files("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/MOTIFS_EXP/",pattern = ".*onestrand.txt")

for(file in dir)
{
  data <- read.csv(paste0("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/MOTIFS_EXP/",file),sep='\t')
  dt <- data.table(data)
  
  motifs <- unique(dt$Motif)
 
  newdt <- copy(dt)
  for(m in motifs)
  {
    cm <- crDNA(m)
    mut <- dt[Motif == cm, MutationCnt]
    psc <- dt[Motif == cm, PlusStrandConsistent]
    msc <- dt[Motif == cm, MinusStrandConsistent]
    psa <- dt[Motif == cm, PlusStrandAll]
    msa <- dt[Motif == cm, MinusStrandAll]
    
    newdt[Motif == m, MutationCnt := MutationCnt + mut]
    newdt[Motif == m, PlusStrandConsistent := PlusStrandConsistent + msc]
    newdt[Motif == m, MinusStrandConsistent := MinusStrandConsistent + psc]
    newdt[Motif == m, PlusStrandAll := PlusStrandAll + msa]
    newdt[Motif == m, MinusStrandAll := MinusStrandAll + psa]
  }
    
  fname <- gsub("_onestrand","",file)
  write.table(newdt,paste0("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/MOTIFS_EXP/",fname), sep="\t",row.names=FALSE,quote=FALSE)
}
