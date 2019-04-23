library(ggplot2)
library(data.table)

data <- read.csv("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/Results_CPP/rt.txt", sep = '\t', header = FALSE)
data <- data.table(data)

d1 <- data[V1 == 0] 
  
p <- ggplot(d1, aes(x = V2, y = V4, color = as.factor(V5))) + geom_point(size=0.01)
ggsave("/Users/mar/BIO/PROJECTS/APOBEC/Project1_TranscriptionLevel/R/1.jpeg", device="jpeg", width=49.9)
