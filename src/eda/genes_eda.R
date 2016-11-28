#install.packages("andrews")

require(matie)
require(minerva)
require(corrplot)
require(ggplot2)
require(reshape)
library(plyr)
require('gtools')
require(grid)

getwd()


#look at distributions
gM <- melt(genes[,c(-1,-2)])
colnames(gM)[c(2,3)] <- c("Gene","Expression")
dp = ggplot(gM, aes(x=Expression)) + geom_line(aes(group=Gene),stat="density", alpha=0.15)
dp + facet_grid(Type ~ .)

#boxplots of all genes separated into patients and controls
patientsM <- melt(lupus[,c(-1,-2,-ncol(genes))])
colnames(patientsM) <- c("Gene","Expression")
patientsM$status <- rep('Lupus',length(patientsM[,1]))
controlsM <- melt(controls[,c(-1,-2,-ncol(genes))])
colnames(controlsM) <- c("Gene","Expression")
controlsM$status <- rep('Control',length(controlsM[,1]))
all <- rbind(patientsM, controlsM)
all$Gene = sub('X','', all$Gene)

geneAnnotations = read.csv('annot.csv')
ff = merge(x=all, y=geneAnnotations, by.x="Gene", by.y="PROBEID")
ff$CHR = factor(ff$CHR,levels = mixedsort(levels(ff$CHR)),ordered = TRUE)

factor(ff$PROBENO,levels(x)[c(4,5,1:3)])
ff = arrange(ff, ff$CHR, ff$CHRLOC, ff$Gene)
unique(ff$PROBENO)
ff$PROBENO = factor(ff$PROBENO,levels = unique(ff$PROBENO),ordered = TRUE)
levels(ff$Gene)
hmcol <- colorRampPalette(brewer.pal(9, "RdYlBu"))(19)
display.brewer.all()

#draw the boxplots of all genes, there's a lot of information in these plots which can make 
#interpretation difficult but I think it gives a reasonable overview of the data
bp <-ggplot(ff, aes(x=as.factor(PROBENO), y=Expression)) + geom_boxplot(aes(fill = CHR))+
  xlab("Probe No") + ylab("Expression")+ scale_fill_manual(values = hmcol)+
  facet_grid(status~.)  
bp + theme_bw() + theme(axis.text.x=element_text(angle=270, vjust = 0.4, size = 6), legend.position="bottom", legend.key.size = unit(.4, "cm"))  

#lupus patients have very high expression levels for many of these genes. Additionally many of genes have unusually large observations 
# as defined by the IQR rule






#plot pearson vs spearman for all pairwise combinations
to.upper<-function(X) X[upper.tri(X,diag=TRUE)]
pc <- to.upper(allCM)
sp <- to.upper(cor(as.matrix(genes[,c(-1,-ncol(genes))]), method="spearman"))
ken <- to.upper(cor(as.matrix(genes[,c(-1,-ncol(genes))]), method="kendall"))




#clean up the a scores first, set the diagonals to 1
?diag
asubset <- aScore[-c(1,96),-c(1,96)]
diag(asubset) <- 1
mas <- to.upper(asubset)

corrs <- as.data.frame(pc)
corrs <- cbind(corrs, sp)
corrs <- cbind(corrs, ken)
corrs <- cbind(corrs, mas)
corrs <- cbind(corrs, pc^2)
corrs <- cbind(corrs, sp^2)
corrs <- cbind(corrs, ken^2)

#all patients, A scores vs person, strong correlation between the two
ggplot(corrs, aes(x=pc, y=mas)) +
  geom_point(alpha=.3)
+
  geom_smooth(alpha=.2, size=1)

lupusA <- read.csv("lupusA.csv")
als <- lupusA[-c(1,97),-c(1,2,97)]
diag(als) <- 1
masLup <- to.upper(als)


corrsL <- as.data.frame(to.upper(cor(as.matrix(genes[genes$diseased=="GREEN",c(-1,-ncol(genes))]), method="pearson")))
corrsL <- cbind(corrsL, to.upper(cor(as.matrix(genes[genes$diseased=="GREEN",c(-1,-ncol(genes))]), method="spearman")))
corrsL <- cbind(corrsL, to.upper(cor(as.matrix(genes[genes$diseased=="GREEN",c(-1,-ncol(genes))]), method="kendall")))
corrsL <- cbind(corrsL, masLup)
corrsL <- cbind(corrsL, to.upper(micL$MIC))

colnames(corrsL) <- c('pc', 'sp', 'ken', 'a', 'mic')
pairs(corrsL)
#lupus patients, clearly A score idenifies some relationships missed by pearson
ggplot(corrsL, aes(x=pc, y=a)) +
  geom_point(alpha=.3)

ggplot(corrsL, aes(x=pc^2, y=a)) +
  geom_point(alpha=.3)

#plot all of the mic outputs pairwise
micM <- as.data.frame(to.upper(micL$MIC))
micM <- cbind(micM, to.upper(micL$MAS))
micM <- cbind(micM, to.upper(micL$MEV))
micM <- cbind(micM, to.upper(micL$MCN))
micM <- cbind(micM, to.upper(micL$MICR2))
colnames(micM) <- c('mic', 'mas', 'mev', 'mcn', 'mic-r2')

pairs(micM)
controlA <- read.csv("controlA.csv")
ctls <- controlA[-c(1,97),-c(1,2,97)]
colnames(ctls)
diag(ctls) <- 1
masCtl <- to.upper(ctls)
micC <- mine(as.matrix(controls[,c(-1,-ncol(genes))]))

#plot differences
diffs <- as.data.frame(masLup)
diffs <- cbind(diffs, to.upper(micL$MIC))
diffs <- cbind(diffs, masLup - to.upper(micL$MIC))
diffs <- cbind(diffs, to.upper(micL$MICR2))
colnames(diffs) <- c('A', 'mic', 'A-mic', 'mic-r2')
pairs(diffs)
plot(masLup - to.upper(micL$MIC))

par(mfrow=c(2,1))
plot(to.upper(micL$MICR2))
plot(to.upper(micC$MICR2))

corrsC <- as.data.frame(to.upper(cor(as.matrix(genes[genes$diseased=="RED",c(-1,-ncol(genes))]), method="pearson")))
corrsC <- cbind(corrsC, to.upper(cor(as.matrix(genes[genes$diseased=="RED",c(-1,-ncol(genes))]), method="spearman")))
corrsC <- cbind(corrsC, to.upper(cor(as.matrix(genes[genes$diseased=="RED",c(-1,-ncol(genes))]), method="kendall")))
corrsC <- cbind(corrsC, masCtl)
corrsC <- cbind(corrsC, to.upper(micC$MIC))

colnames(corrsC) <- c('pc', 'sp', 'ken', 'mas', 'mic')
pairs(corrsC)
#lupus patients, clearly A score idenifies some relationships missed by pearson
ggplot(corrsC, aes(x=pc, y=mas)) +
  geom_point(alpha=.3)


length(mas)



require(andrews)
andrews(genes[,-1],clr=95)
par(mfrow=c(2,1))
andrews(controls[,-1],type=4)
andrews(lupus[,-1],type=4)



?andrews
par(mfrow=c(2,1))
andrews(controls.pca$x[,1:15],type=2, ymax=5)
andrews(lupus.pca$x[,1:15],type=2, ymax=5)
#not very useful

#tails off around 6

dim(bd)




aScore$A
pValue <- ma.test(d,aScore)
hmAll <- amap(genes)
hmControls <- amap(genes[genes$diseased=='GREEN',])
hmLupus <- amap(genes[genes$diseased=='RED',])
#very strong gene clustering in lupus patients

hmControls
par(mfrow=c(1,2))
fdg(controls,dataName="Controls",method="A",cutoff=0.5,dim=2)
fdg(lupus,dataName="Lupus",method="A",cutoff=0.5,dim=2)


trainAllSD =  apply(genes[,-c(1,96)], 2, sd)
trainControlSD =  apply(controls[,-c(1,96)], 2, sd)
trainLupusSD =  apply(lupus[,-c(1,96)], 2, sd)

trainAllIQR=  apply(genes[,-c(1,96)], 2, IQR)
trainControlIQR =  apply(controls[,-c(1,96)], 2, IQR)
trainLupusIQR =  apply(lupus[,-c(1,96)], 2, IQR)


############################################exaxime variability
par(mfrow=c(2,3))
hist(trainAllSD)
hist(trainControlSD)
hist(trainLupusSD)

hist(trainAllIQR)
hist(trainControlIQR)
hist(trainLupusIQR)

twoWay = union(lll$X, lll$Y)
setdiff(twoWay, which(trainAllSD >= quantile(trainAllSD, 72/94)))
setdiff(twoWay, which(trainAllIQR>= quantile(trainAllIQR, 74/94)))

setdiff(de, which(trainAllSD >= quantile(trainAllSD, 70/94)))
setdiff(de, which(trainAllIQR>= quantile(trainAllIQR, 69/94)))
twoWay+1
twoWayNames = colnames(genes[,twoWay+1])
dim(controlBD)
bdAllSD =  apply(allBD, 2, sd)
bdControlSD =  apply(controlBD, 2, sd)
bdLupusSD =  apply(lupusBD, 2, sd)

bdAllIQR=  apply(allBD, 2, IQR)
bdControlIQR =  apply(controlBD, 2, IQR)
bdLupusIQR =  apply(lupusBD, 2, IQR)

par(mfrow=c(2,3))
hist(bdAllSD)
hist(bdControlSD)
hist(bdLupusSD)

hist(bdAllIQR)
hist(bdControlIQR)
hist(bdLupusIQR)

probs = (1:54675)/54675

quantile(bdAllIQR, probs=(53675/54675))
#probs[9657]
?quantile
#subBD = allBD[,bdAllIQR >= quantile(bdAllIQR, probs[9675])]
length(which(bdAllSD >= quantile(bdAllSD, probs)))
subBD = allBD[,bdAllIQR >= quantile(bdAllIQR, probs=(53675/54675))]
subBD2 = allBD[,bdAllSD >= quantile(bdAllSD, probs=(53675/54675))]

setdiff(twoWayNames, colnames(subBD2))
setdiff(colnames(genes[,-c(1,96)]), colnames(subBD2))

(length(which(controlCM >= 0.7 & controlCM < 1))/length(controlCM))/2

(length(which(lupusCM >= 0.7 & lupusCM < 1))/length(lupusCM))/2
)
140/8836
