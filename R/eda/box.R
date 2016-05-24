#boxplots of all genes separated into patients and controls
library(ggplot2)
library(reshape2)

#need to reshape the data first for ggplot
#this creates a new 3 column data.frame. The columns are: 
#Gene= the name of the gene; 
#Expression = the gene expression level; 
#status= a categorical variable, Patient or Control
all <- melt(genes[,c(-1,-2,-ncol(genes))], id.vars = 'diseased')
colnames(all) <- c("status","Gene","Expression")

#draw the boxplots of all genes, there's a lot of information in these plots which can make 
#interpretation difficult but I think it gives a reasonable overview of the data
bp <- ggplot(all, aes(x=Gene, y=Expression)) + geom_boxplot() + ylab("Expression") + 
  facet_grid(status ~ .) 
bp + theme(axis.text.x=element_text(angle=270, vjust = 0.5))
#lupus patients have very high expression levels for many of these genes. Additionally many of genes have unusually large observations 
# as defined by the IQR rule

library(andrews)
andrews(genes[,-1],clr=95)
par(mfrow=c(2,1))
andrews(controls[,-1],type=4)
andrews(lupus[,-1],type=4)

allCM <- cor(as.matrix(genes[,c(-1,-ncol(genes))]))

corrplot(cor(lupus[,-c(1,2,ncol(lupus))]),tl.pos='d', order = "original", method="color", type='upper',tl.cex=.5)
corrplot(cor(controls[,-c(1,2,ncol(controls))]),tl.pos='d', order = "original", method="color", type='upper',tl.cex=.5)

controlCM <- cor(as.matrix(controls[,c(-1,-ncol(genes))]))
lupusCM <- cor(as.matrix(lupus[,c(-1,-ncol(genes))]))
