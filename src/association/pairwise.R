#results = read.csv('results.csv')
library(gplots)
require(rafalib)
require(ggplot2)
mypar()
res = unique(results[results$Sample == 'Control' & results$Value >= 0.90,c(1,2,4,5)])
res[,1]
sc = controls[,-c(1,96)]
par(mfrow=c(2,3))

plot(sc[,unlist(res[1,1])], sc[,unlist(res[1,2])],xlab = res[1,1], ylab = res[1,2])
plot(sc[,unlist(res[2,1])], sc[,unlist(res[2,2])],xlab = res[2,1], ylab = res[2,2])
plot(sc[,unlist(res[3,1])], sc[,unlist(res[3,2])],xlab = res[3,1], ylab = res[3,2])

sc = lupus[,-c(1,96)]
plot(sc[,unlist(res[1,1])], sc[,unlist(res[1,2])],xlab = res[1,1], ylab = res[1,2])
plot(sc[,unlist(res[2,1])], sc[,unlist(res[2,2])],xlab = res[2,1], ylab = res[2,2])
plot(sc[,unlist(res[3,1])], sc[,unlist(res[3,2])],xlab = res[3,1], ylab = res[3,2])

results[results$Sample == 'Lupus' & results$From > 0.85 & results$X %in% res[,1] & results$Y %in% res[,2],]

results[results]
allRes = unique(results[results$Sample == 'Lupus' & results$Value >= 0.90,3])

r2ResC = results[results$Sample == 'Control' & results$Association == "R2" & results$Value >= 0.90,3]
aResC = results[results$Sample == 'Control' & results$Association == "A" & results$Value >= 0.90,3]
dcorResC = results[results$Sample == 'Control' & results$Association == "dcor" & results$Value >= 0.90,3]
micResC = results[results$Sample == 'Control' & results$Association == "MIC" & results$Value >= 0.90,3]

r2ResL = results[results$Sample == 'Lupus' & results$Association == "R2" & results$Value >= 0.90,3]
aResL = results[results$Sample == 'Lupus' & results$Association == "A" & results$Value >= 0.90,3]
dcorResL = results[results$Sample == 'Lupus' & results$Association == "dcor" & results$Value >= 0.90,3]
micResL = results[results$Sample == 'Lupus' & results$Association == "MIC" & results$Value >= 0.90,3]

venn( list('R2'=r2ResC,'A'=aResC,'dcor'=dcorResC, 'MIC'=micResC) )
venn( list('R2'=r2ResL,'A'=aResL,'dcor'=dcorResL, 'MIC'=micResL) )

allRes = results[results$Sample == 'Lupus' & results$Value >= 0.95,c(1:3)]
r2Res = results[results$Sample == 'Lupus' & results$Association == "R2" & results$From >=0.95,c(1,2,3)]
aRes = results[results$Sample == 'Lupus' & results$Association == "A" & results$From >= 0.95,c(1,2,3)]
dcorRes = results[results$Sample == 'Lupus' & results$Association == "dcor" & results$From >= 0.95,c(1,2,3)]
micRes = results[results$Sample == 'Lupus' & results$Association == "MIC" & results$From >= 0.95,c(1,2,3)]

xgenes = results

nlRes = setdiff(allRes, r2Res)
rownames(nlRes)
dcorNL = results[results$Pair %in% nlRes$Pair & results$Sample=='Lupus' & results$Association == 'dcor',]

r2Respar(mfrow= c(4,4))
 
#mic
res = micRes
for(i in 1:4)
{
  plot(lupus[,res[i,1]], lupus[,res[i,2]],xlab = res[i,1], ylab = res[i,2])
}

#A
res = aRes
for(i in 1:4)
{
  plot(lupus[,res[i,1]], lupus[,res[i,2]],xlab = res[i,1], ylab = res[i,2])
}

#A
par(mfrow= c(4,4))
res = dcorNL
for(i in 33:48)
{
  plot(lupus[,res[i,1]], lupus[,res[i,2]],xlab = res[i,1], ylab = res[i,2])
}

setdiff(allRes, aRes)
setdiff(allRes, dcorRes)
setdiff(allRes, micRes)

#ooks like outliers
results[results$Sample == 'Lupus' & results$Association == "A" & results$From >= 0.95,]
dcor(lupus$'209835_x_at', lupus$'212014_x_a')
bdcor(lupus$'209835_x_at', lupus$'212014_x_a')

dcor(lupus$'213797_at', lupus$'219863_at')
mean(bdcor(lupus$'213797_at', lupus$'219863_at'))

?venn

nlRes
nlRes = results[results$Sample == 'Lupus' & results$Association == "R2" & results$From > 0.90,3]
dim(res)

results[results$Sample == 'Lupus' & results$Association == "MIC" & results$Value >= 0.99,]
results2 = arrange(results, X)
fit1 = glm(genes$diseased~genes$'213797_at'*genes$'219863_at', binomial(probit))
summary(fit1)

plot(glm(genes$diseased~genes$'210349_at'*genes$'229029_at', binomial(probit)))


ggplot(lupus, aes(x=lupus$'213797_at', y=lupus$'229029_at')) + geom_point(alpha=.9) + geom_smooth(method="gam", alpha=.8, size=1) + theme_bw() + theme(legend.position="none")

#calculate AUC
res = results[results$Sample == 'Lupus' & results$Association == 'dcor' & results$From >=0.95,]

auc = function(res, linpreds=F)
{
  y = genes$diseased == 'Control'
  rr = genes[,-c(1,96)]
  aucs = nrow(res)
  for(i in 1:nrow(res))
  {
    #preds=predict(earth(y = y, x = rr[,c(res$X[i],res$Y[i])], degree=2, linpreds=linpreds))
    #roc1=roc(y ~ preds)
    #aucs[i] = roc1$auc
    aucs[i] = earth(y = y, x = rr[,c(res$X[i],res$Y[i])], degree=2, linpreds=linpreds)$rsq
  }
  aucs
}
require(earth)
require(pROC)
max(aucDcor)
aucDcor = auc(results, T)
results = cbind(results, aucDcor)
colnames(results)[11] = "AUC"
results = cbind(results, auc(results, T))
colnames(results)[12] = "Linear AUC"

plot(res)
any(results$AUC < results$`Linear AUC`)
results[which(results$AUC < results$`Linear AUC`),]

plot(results$AUC, results$`Linear AUC`)
mean(results$AUC- results$`Linear AUC`)
sd(results$AUC- results$`Linear AUC`)

plot(roc1)

plot(aucDcor, res$Value)
require(ggplot2)
dp = ggplot(results, aes(x=Value, y=AUC)) + geom_point(alpha=.7) 
dp + facet_grid(Association ~ .)

results[which(results$Value >= 0.99),]

require(rafalib)
mypar()
par(mfrow=c(5,4))
lll = unique(results[results$AUC >= 0.91 & results$AUC > results$`Linear AUC` & results$Association != 'R2', c(1,2)])
for(i in 1:nrow(lll))
{
  plot(lupus[,lll[i,1]+1], lupus[,lll[i,2]+1])
}

ll = results[results$Association == "R2", 3]
nl = unique(results[results$Value >= 0.98, c(1,2)])

setdiff(ll$X, nl$X)


lll = arrange(results[!results$Pair %in% ll & results$Association =='dcor' & results$Value >= 0.925,],desc(Value))
?setdiff
par(mfrow=c(4,3))
require(ggplot2)
require(gridExtra)

#how many probes
su = union(results[results$Association =='dcor' & results$Value >= 0.95,1:2],results[results$Association =='dcor' & results$Value >= 0.95,2])
#probes correlated in both controls and lupus
clProbes = c(31,30,76,86,41,42)

length(su)

#how many genes)
gr = unique(annot[annot$PROBENO %in% su,6])
length(gr)

#how many chromosomes
cr = unique(annot[annot$PROBENO %in% su,3])
length(cr) 

#how many differentially expressed
de = annot[annot$PROBEID %in% diffExp,1]

length(de)
intersect(su, de)


ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[1,1]+1]), y=ggname(colnames(lupus)[lll[1,2]+1]), colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")

# pp1 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[1,1]+1]), y=ggname(colnames(lupus)[lll[1,2]+1]))) +  geom_point(size=1,alpha=.5) + geom_smooth() + theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# pp2 = ggplot(lupus, aes_string(x=ggname(colnames(lupus)[lll[2,1]+1]), y=ggname(colnames(lupus)[lll[2,2]+1]))) +  geom_point(size=1,alpha=.5) + geom_smooth() + theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# pp3 = ggplot(lupus, aes_string(x=ggname(colnames(lupus)[lll[3,1]+1]), y=ggname(colnames(lupus)[lll[3,2]+1]))) +  geom_point(size=1,alpha=.5) + geom_smooth()+ theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# pp4 = ggplot(lupus, aes_string(x=ggname(colnames(lupus)[lll[4,1]+1]), y=ggname(colnames(lupus)[lll[4,2]+1]))) +  geom_point(size=1,alpha=.5) + geom_smooth()+ theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# pp5 = ggplot(lupus, aes_string(x=ggname(colnames(lupus)[lll[5,1]+1]), y=ggname(colnames(lupus)[lll[5,2]+1]))) +  geom_point(size=1,alpha=.5) + geom_smooth()+ theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# pp6 = ggplot(lupus, aes_string(x=ggname(colnames(lupus)[lll[6,1]+1]), y=ggname(colnames(lupus)[lll[6,2]+1]))) +  geom_point(size=1,alpha=.5) + geom_smooth()+ theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# pp7 = ggplot(lupus, aes_string(x=ggname(colnames(lupus)[lll[7,1]+1]), y=ggname(colnames(lupus)[lll[7,2]+1]))) +  geom_point(size=1,alpha=.5) + geom_smooth()+ theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# pp8 = ggplot(lupus, aes_string(x=ggname(colnames(lupus)[lll[8,1]+1]), y=ggname(colnames(lupus)[lll[8,2]+1]))) +  geom_point(size=1,alpha=.5) + geom_smooth()+ theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# pp9 = ggplot(lupus, aes_string(x=ggname(colnames(lupus)[lll[9,1]+1]), y=ggname(colnames(lupus)[lll[9,2]+1]))) +  geom_point(size=1,alpha=.5) + geom_smooth()+ theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# pp10 = ggplot(lupus, aes_string(x=ggname(colnames(lupus)[lll[10,1]+1]), y=ggname(colnames(lupus)[lll[10,2]+1]))) +  geom_point(size=1,alpha=.5) + geom_smooth()+ theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# pp11 = ggplot(lupus, aes_string(x=ggname(colnames(lupus)[lll[11,1]+1]), y=ggname(colnames(lupus)[lll[11,2]+1]))) +  geom_point(size=1,alpha=.5) + geom_smooth()+ theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# pp12 = ggplot(lupus, aes_string(x=ggname(colnames(lupus)[lll[12,1]+1]), y=ggname(colnames(lupus)[lll[12,2]+1]))) + geom_point(size=1,alpha=.5) + geom_smooth()+ theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pp1 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[1,1]+1]), y=ggname(colnames(lupus)[lll[1,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth() + theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pp2 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[2,1]+1]), y=ggname(colnames(lupus)[lll[2,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth() + theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pp3 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[3,1]+1]), y=ggname(colnames(lupus)[lll[3,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pp4 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[4,1]+1]), y=ggname(colnames(lupus)[lll[4,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pp5 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[5,1]+1]), y=ggname(colnames(lupus)[lll[5,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pp6 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[6,1]+1]), y=ggname(colnames(lupus)[lll[6,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pp7 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[7,1]+1]), y=ggname(colnames(lupus)[lll[7,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pp8 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[8,1]+1]), y=ggname(colnames(lupus)[lll[8,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pp9 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[9,1]+1]), y=ggname(colnames(lupus)[lll[9,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pp10 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[10,1]+1]), y=ggname(colnames(lupus)[lll[10,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pp11 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[11,1]+1]), y=ggname(colnames(lupus)[lll[11,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
pp12 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[12,1]+1]), y=ggname(colnames(lupus)[lll[12,2]+1]), colour='diseased')) + geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


grid.arrange(pp1, pp2, pp3, pp4,pp5, pp6, pp7, pp8,pp9, pp10, pp11,pp11, ncol=3, nrow=4,  main="", plot=T)

colnames(lupus)[lll[i,1]+1]



plot(lupus[,lll[i,1]+1], lupus[,lll[i,2]+1])
curve(loess(lupus[,lll[i,1]+1]~lupus[,lll[i,2]+1]))
?geom_smooth

ggname <- function(x) {
  if (class(x) != "character") {
    return(x)
  }
  y <- sapply(x, function(s) {
    if (!grepl("^`", s)) {
      s <- paste("`", s, sep="", collapse="")
    }
    if (!grepl("`$", s)) {
      s <- paste(s, "`", sep="", collapse="")
    }
  }
  )
  y 
}


