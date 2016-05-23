library(plyr)
library(readr)

#read
genes <- read_csv("data/GeneExpressionData.csv")
samples <- read_csv("data/SampleMap.csv")

#append sample type
genes$diseased <- samples[,4]
revalue(genes$diseased, c("GREEN"="Control", "RED"="Lupus")) -> genes$diseased

#id different sample sources 
controlTypes <- rep('C',53)
lupusTypes <- c(rep("L1", 151), rep('L2',83),rep('L3',22),rep('L4',106),rep('L5',5))
bothTypes <- factor(c(controlTypes, lupusTypes))
genes$Source <-  bothTypes

#for convenience
controls <- genes[genes$diseased=='Control',]
lupus <- genes[genes$diseased=='Lupus',]

annot <- read_csv('data/annot.csv')
