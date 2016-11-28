#source("http://www.bioconductor.org/biocLite.R")
#biocLite("limma", destdir="C:/Users/jroche1x/Documents")
require(reshape)
require(qvalue)
require(ggplot2)
#plot differential expression
safeControlBD = as.numeric(safeControlBD)

typeof(bd[1,1])
safeControlBD[1,1]
cm = apply(safeControlBD, 2, mean)
lmm = apply(safeLupusBD, 2, mean)
length(lmm)
hist(cm)
hist(lmm)
plot(lmm,cm)


#MA plot
badf = data.frame(cm, lmm)
badf = cbind(badf, (cm+lmm)/2)
badf = cbind(badf, lmm-cm)
colnames(badf) = c('Control', 'Lupus', 'Mean', 'Difference')
ggplot(badf,aes(x=Mean, y=Difference)) + geom_point(alpha=0.4)




require(limma)

#compare controls to lupus
#prepare the model

bigg = factor(c(rep('Control',nrow(controlBD)), rep('Lupus', nrow(lupusBD))))
bigg = factor(c(rep('Control',29), rep('Lupus',396-29)))
#g = sample(c(T, F), 599, replace=T)
bigy = t(as.matrix(safeBDAll))
dim(bigy)
#fit it
bigfit = lmFit(bigy, design=model.matrix(~ bigg))

bigfit = eBayes(bigfit)
dim(bigy)

resultsBig = topTable(bigfit, number=54675)

plot(lmm-cm,-log10(qvalue(bigfit$p.value[,2])$qvalues))
#volcano plot based on fold cachange of 1
vdf = data.frame(resultsBig$logFC, -log10(resultsBig$adj.P.Val))
vdf = cbind(vdf, resultsBig$adj.P.Val < 0.05 & resultsBig$logFC > 1)
colnames(vdf) = c('LFC', 'pvalue', 'DE')
p = ggplot(vdf,aes(x=LFC, y=pvalue)) + geom_point() + geom_vline(xintercept = c(1,-1)) + geom_hline(yintercept = -log10(0.05))
p + xlab('Log Fold Change') + ylab('Adjusted p-value')



colnames(allBD)[1:10]
ncol(allBD)

resultsBig = cbind(resultsBig, match(rownames(resultsBig), colnames(safeBDAll)))
resultsBig = cbind(resultsBig, rownames(resultsBig))
colnames(resultsBig)[7:8] = c('ProbeNo', 'ProbeID')


intersect(results[results$logFC >= 1, 8], resultsBig[resultsBig$logFC >= 1, 8])
intersect(results[results$logFC >= 1, 8], colnames(genes))
resultsBig$`Differential Expressed in Test` =  rep(0, nrow(resultsBig))
for(i in 1:35)
{
  if(resultsBig[i, 8] %in% intersect(results[results$logFC >= 1, 8], resultsBig[resultsBig$logFC >= 1, 8]))
  {
    resultsBig[i, 9] = 1
  }
  else
  {
    resultsBig[i, 9] = 0
  }
}
colnames(resultsBig)[9] = 'Differential Expressed in Test'

require(plyr)
resultsBig = arrange(resultsBig, desc(logFC))
count(resultsBig$logFC >= 1)
count(resultsBig$adj.P.Val <= 0.01)
count(resultsBig$P.Value <= 0.01)

#main results
x.small <- xtable(resultsBig[1:119,c(7,8,1)], label = 'tabsmall', caption = 'Differentially expressed genes')
print(x.small,latex.environments = "",table.placement = 'h')

#difference from test
diffRes = resultsBig[resultsBig$ProbeID %in%  setdiff(results[results$logFC >= 1, 8], resultsBig[resultsBig$logFC >= 1, 8]),]
x.small <- xtable(diffRes[,c(7,8,1)], label = 'tabsmall', caption = 'Differentially expressed genes')
print(x.small,latex.environments = "",table.placement = 'h')
