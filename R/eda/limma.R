#source("http://www.bioconductor.org/biocLite.R")
#biocLite("limma", destdir="C:/Users/jroche1x/Documents")
require(reshape)
require(qvalue)
#plot differential expression
cm = mapply(mean,genes[genes$diseased=='Control',-c(1,96,97)])
lmm = mapply(mean,genes[genes$diseased=='Lupus',-c(1,96,97)])
hist(cm)
hist(lmm)
plot(lmm,cm)
#MA plot
badf = data.frame(cm, lmm)
badf = cbind(badf, (cm+lmm)/2)
badf = cbind(badf, lmm-cm)
colnames(badf) = c('Control', 'Lupus', 'Mean', 'Difference')
ggplot(badf,aes(x=Mean, y=Difference)) + geom_point()

length(cm)
plot((cm+lmm)/2, lmm-cm)

plot(lmm-cm,-log10(qvalue(fit$p.value[,2])$qvalues))
require(limma)
require(ggplot2)

#compare controls to lupus
#prepare the model
g = factor(genes$diseased)
y = t(as.matrix(genes[,-c(1,96,97)]))

#fit it
fit = lmFit(y, design=model.matrix(~ g))
fit = eBayes(fit)

#volcano plot based on fold cachange of 1
results = topTable(fit, number=94)
plot(results$logFC, -log10(results$adj.P.Val))
vdf = data.frame(results$logFC, -log10(results$adj.P.Val))
vdf = cbind(vdf, results$adj.P.Val < 0.05 & results$logFC > 1)
colnames(vdf) = c('LFC', 'pvalue', 'DE')
p = ggplot(vdf,aes(x=LFC, y=pvalue)) + geom_point() + geom_vline(xintercept = 1) + geom_hline(yintercept = -log10(0.05))
p + xlab('Log Fold Change') + ylab('Adjusted p-value')

which('202411_at' == colnames(genes[,-1]))
which('202411_at' == rownames(results))

diffExp = rownames(results[results$logFC > 1, ])

length(colnames(genes)[-c(1,96)])
which(colnames(genes) == '202411_at')
?match
match(colnames(genes)[-c(1,96)],rownames(results))
match(rownames(results),colnames(genes)[-c(1,96)])
which(colnames(genes)[-1] %in% diffExp)
results = cbind(results, match(rownames(results), colnames(genes)[-c(1,96)]))
results = cbind(results, rownames(results))

colnames(genes[,-1])[58]
colnames(results)[7:8] = c('ProbeNo', 'ProbeID')
#print results
require(xtable)
fit2 <- treat(fit,lfc=1)
dim(results)
pr = data.frame(results)
require(plyr)
pr = arrange(pr, desc(logFC))

?arrange

x.small <- xtable(pr[1:length(which(results$logFC > 1)),c(7,8,1)], label = 'tabsmall', caption = 'Differentially expressed genes')
print(x.small,latex.environments = "",table.placement = 'h')





fit$adj
# Valcano plot
plot(fit$Amean, -log10(fit$p.value[,2]))

#

genes$Type
mm = model.matrix(~ genes$Type)
colnames(mm) = c('C','L1','L2','L3','L4', 'L5')
mfit = lmFit(y, design=mm)
mfit = eBayes(mfit)
topTable(mfit)
?makeContrasts
lupusTypesFromControlsContrasts= c('L1-C','L2-C','L3-C','L4-C', 'L5-C')
controlContrasts= c('C2-C3','C2-C4','C3-C4')
lupusContrasts= c('L1-L2','L1-L3','L1-L4', 'L1-L5','L2-L3','L2-L4', 'L2-L5', 'L3-L4','L3-L5', 'L4-L5')
#contrast.matrix <- makeContrasts(contrasts=lupusContrasts, levels= c('C','L1','L2','L3','L4','L5'))
contrast.matrix <- makeContrasts(contrasts=lupusTypesFromControlsContrasts, levels= c('C','L1','L2','L3','L4','L5'))
fit2 <- contrasts.fit(mfit, contrast.matrix)

fit2 <- eBayes(fit2)

results2 <- decideTests(fit2, adjust.method = "BH",p.value=0.01)
vennDiagram(results2)
vennCounts(results2)

results2 = topTableF(fit2, lfc=1, number=94)
results2 = cbind(results2, match(rownames(results2), colnames(genes)[-c(1,96)]))
results2 = cbind(results2, rownames(results2))

dim(results2)
colnames(results2)[10:11] = c('ProbeNo', 'ProbeID')
pr = data.frame(results2)
pr = arrange(pr, desc(results2$F))
results2$F

betweenTypes = topTreat(fit2,lfc=1)
betweenTypes
#print results

topTable(fit2, coef=1, adjust="BH")
topTable(fit2, coef=2, adjust="BH")
topTable(fit2, coef=3, adjust="BH")
topTable(fit2, coef=4, adjust="BH")
topTable(fit2, coef=5, adjust="BH")

ggg = decideTests(fit2, lfc=1)
treat(ggg)
?decideTests
vennDiagram(ggg)
ttt = topTableF(fit2, lfc=1, number=94)

ttt$F


fit2 <- treat(fit,lfc=1)
pr2 = data.frame(betweenTypes[,c(1,2)])
x.small <- xtable(pr[1:20,10:11], label = 'tabsmall', caption = 'Differentially expressed genes between all types')
print(x.small,latex.environments = "",table.placement = 'h')

?treat
?topTreat

citation("limma")
