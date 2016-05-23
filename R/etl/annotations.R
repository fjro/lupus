bd = read.csv("Big Data.csv")
colnames(bd)
dim(bd)
setdiff(samples$Sample, colnames(bd))
?setdiff
summary(bd[,1:5])
summary(e[,1:5])

?ExpressionSet
ge = ExpressionSet(t(genes))
library(devtools)
#install_github("genomicsclass/GSE5859")
library(GSE5859)
data(GSE5859)
class(ge)
dat = exprs(ge)
dim(dat)

colnames(genes)
head(genes)
colnames(genes) = sub('X','', colnames(genes) )

require(hgfocus.db)
citation("hgfocus.db") 
annot = select(hgfocus.db, 
               keys=featureNames(ge), 
               keytype="PROBEID", 
               columns=c("CHR", "CHRLOC", "SYMBOL"))
## here we pick one column from the annotation
annot = annot[ match(featureNames(ge),annot$PROBEID), ]

length(which(is.na(annot$CHR)))
write.csv(annot, "annot.csv")
head(annot)
dim(annot)

sampleInfo = pData(ge)
dim(sampleInfo)
head(sampleInfo)

library(Homo)

?select
select(hgfocus.db, 
       keys=featureNames(ge)[3], 
       keytype="PROBEID", 
       columns=c("CHR", "CHRLOC", "SYMBOL"))

#source("http://www.bioconductor.org/biocLite.R")
#biocLite("limma", destdir="C:/Users/jroche1x/Documents")
require(limma)
?lmFit
?eBayes
?topTable
gm = t(as.matrix(genes[,-c(1,96)]))
dim(gm)
design = c(rep(1,53), rep(0,367))
length(design)
fit = lmFit(gm,design)
fit <- eBayes(fit)
topTable(fit)

sd <- 0.3*sqrt(4/rchisq(100,df=4))
y <- matrix(rnorm(100*6,sd=sd),100,6)
rownames(y) <- paste("Gene",1:100)
y[1:2,4:6] <- y[1:2,4:6] + 2
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
options(digits=3)

# Ordinary fit
fit <- lmFit(y,design)
fit <- eBayes(fit)
topTable(fit)
dim(fit)
colnames(fit)
rownames(fit)[1:10]
names(fit)

# Fold-change thresholding
fit2 <- treat(fit,lfc=0.1)
topTreat(fit2,coef=2)

# Volcano plot
volcanoplot(fit,highlight=2)
?volcanoplot

ge

# MA plot
plot(fit)

# Q-Q plot of moderated t-statistics
qqt(fit$t[,2],df=fit$df.residual+fit$df.prior)
abline(0,1)


####big
colnames(genes) = sub('X','', colnames(genes) )

require(hgfocus.db)
dim(ge)
annotBig = select(hgfocus.db, 
               keys=featureNames(ge), 
               keytype="PROBEID", 
               columns=c("CHR", "CHRLOC", "SYMBOL"))
## here we pick one column from the annotation
annotBig = annotBig[ match(featureNames(ge),annot$PROBEID), ]

?ExpressionSet
ge = ExpressionSet(t(subBD))
library(devtools)
#install_github("genomicsclass/GSE5859")
library(GSE5859)
data(GSE5859)
class(ge)
dat = exprs(ge)
dim(dat)

colnames(genes)
head(genes)
colnames(genes) = sub('X','', colnames(genes) )
