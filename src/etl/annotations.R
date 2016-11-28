library(devtools)
#install_github("genomicsclass/GSE5859")
source("http://www.bioconductor.org/biocLite.R")
#biocLite("hgfocus.db")
library(GSE5859)
library(hgfocus.db)
library(Homo)

#make an expression set from the training data
ge <- ExpressionSet(t(genes))

data(GSE5859)
class(ge)
dat <- exprs(ge)

citation("hgfocus.db") 

#annotate Affymetrix exprssion data with additional info such as chromosome location
#gene name etc.
select(hgfocus.db, 
       keys=featureNames(ge)[3], 
       keytype="PROBEID", 
       columns=c("CHR", "CHRLOC", "SYMBOL"))

annot <- select(hgfocus.db, 
               keys=featureNames(ge), 
               keytype="PROBEID", 
               columns=c("CHR", "CHRLOC", "SYMBOL"))
## here we pick one column from the annotation
annot = annot[ match(featureNames(ge),annot$PROBEID), ]



sampleInfo <- pData(ge)



write.csv(annot, "data/annot.csv")
