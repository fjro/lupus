hX = results[results$Association == 'dcor' & results$Value > 0.95,1]
hY = results[results$Association == 'dcor' & results$Value > 0.95,2]
hb2 = c(hX, hY)
plot(table(hb2))



rR = lupusDcor - controlDcor
vv = as.matrix(rR)
vv[lower.tri(vv, diag = T)] <- 0
max(vv)
d2 = as.matrix(lupusDcor)
d2[lower.tri(d2, diag = T)] <- 0
r2 = as.matrix(lupusCM^2)
r2[lower.tri(r2, diag = T)] <- 0
res = which(vv >0.45 & d2 >0.90 & (d2 - r2) > 0.1, arr.ind = T)


hb = c(res[,1], res[,2])

setdiff(hb, de)
diffExp
unique(annot[annot$PROBENO %in% hb,6])
plot(table(hb))
res

lupusDcor.NL[74,81]
lupusA.NL[74,81]
lupusMIC.NL[74,81]

for (i in 1:nrow(res))
{
  print(lupusDcor.NL[res[i,1],res[i,2]])
  #print(lupusDcor[res[i,1],res[i,2]] - lupusCM[res[i,1],res[i,2]]^2)
}



lll = res
lll = cbind(lll, rep(0, nrow(lll)))
lll = cbind(lll, rep(0, nrow(lll)))
lll = cbind(lll, rep(0, nrow(lll)))
lll = cbind(lll, rep(0, nrow(lll)))
colnames(lll)[3] = 'diff'
colnames(lll)[4] = 'rA'
colnames(lll)[5] = 'rA/A'
colnames(lll)[6] = 'dcor'
require(plyr)
diff2 = lupusDcor - lupusCM^2
for (i in 1:nrow(res))
{
  #print(lupusDcor.NL[res[i,1],res[i,2]])
  lll[i,3] = diff2[res[i,1], res[i,2]]
  #lll[i,4] = lupusDcor.NL[res[i,1], res[i,2]]
  lll[i,5] = diff2[res[i,1], res[i,2]]/lupusDcor[res[i,1], res[i,2]]
  lll[i,6] = lupusDcor[res[i,1], res[i,2]]
}

dim(lll)
lll = arrange(as.data.frame(lll),desc(diff))
require(ggplot2)
require(gridExtra)


pp1 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[1,1]+1]), y=ggname(colnames(lupus)[lll[1,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth() + theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + xlab(lll[1,1]) +ylab(lll[1,2])+ scale_colour_brewer(palette="Set1")
pp2 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[2,1]+1]), y=ggname(colnames(lupus)[lll[2,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth() + theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[2,1]) +ylab(lll[2,2])+ scale_colour_brewer(palette="Set1")
pp3 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[3,1]+1]), y=ggname(colnames(lupus)[lll[3,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[3,1]) +ylab(lll[3,2])+ scale_colour_brewer(palette="Set1")
pp4 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[4,1]+1]), y=ggname(colnames(lupus)[lll[4,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[4,1]) +ylab(lll[4,2])+ scale_colour_brewer(palette="Set1")
pp5 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[5,1]+1]), y=ggname(colnames(lupus)[lll[5,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[5,1]) +ylab(lll[5,2])+ scale_colour_brewer(palette="Set1")
pp6 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[6,1]+1]), y=ggname(colnames(lupus)[lll[6,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[6,1]) +ylab(lll[6,2])+ scale_colour_brewer(palette="Set1")
pp7 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[7,1]+1]), y=ggname(colnames(lupus)[lll[7,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[7,1]) +ylab(lll[7,2])+ scale_colour_brewer(palette="Set1")
pp8 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[8,1]+1]), y=ggname(colnames(lupus)[lll[8,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[8,1]) +ylab(lll[8,2])+ scale_colour_brewer(palette="Set1")
pp9 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[9,1]+1]), y=ggname(colnames(lupus)[lll[9,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[9,1]) +ylab(lll[9,2])+ scale_colour_brewer(palette="Set1")
pp10 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[10,1]+1]), y=ggname(colnames(lupus)[lll[10,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[10,1]) +ylab(lll[10,2])+ scale_colour_brewer(palette="Set1")
pp11 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[11,1]+1]), y=ggname(colnames(lupus)[lll[11,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[11,1]) +ylab(lll[11,2])+ scale_colour_brewer(palette="Set1")
pp12 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[12,1]+1]), y=ggname(colnames(lupus)[lll[12,2]+1]), colour='diseased')) + geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[12,1]) +ylab(lll[12,2])+ scale_colour_brewer(palette="Set1")

pp13 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[13,1]+1]), y=ggname(colnames(lupus)[lll[13,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[13,1]) +ylab(lll[13,2])+ scale_colour_brewer(palette="Set1")
pp14 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[14,1]+1]), y=ggname(colnames(lupus)[lll[14,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[14,1]) +ylab(lll[14,2])+ scale_colour_brewer(palette="Set1")
pp15 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[15,1]+1]), y=ggname(colnames(lupus)[lll[15,2]+1]), colour='diseased')) + geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[15,1]) +ylab(lll[15,2])+ scale_colour_brewer(palette="Set1")
pp16 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[16,1]+1]), y=ggname(colnames(lupus)[lll[16,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[16,1]) +ylab(lll[16,2])+ scale_colour_brewer(palette="Set1")

pp17 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[17,1]+1]), y=ggname(colnames(lupus)[lll[17,2]+1]), colour='diseased')) + geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[17,1]) +ylab(lll[17,2])+ scale_colour_brewer(palette="Set1")
pp18 = ggplot(genes, aes_string(x=ggname(colnames(lupus)[lll[18,1]+1]), y=ggname(colnames(lupus)[lll[18,2]+1]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(lll[18,1]) +ylab(lll[18,2])+ scale_colour_brewer(palette="Set1")

?arrangeGrob
grid.arrange(arrangeGrob(pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10,pp11,pp12,pp13,pp14,pp15,pp16,pp17,pp18, ncol=4, nrow=5))

wilcox.test(to.upper(controlCM^2),to.upper(lupusCM^2),paired=T) 
wilcox.test(to.upper(controlA),to.upper(lupusA),paired=T) 
wilcox.test(to.upper(controlDcor),to.upper(lupusDcor),paired=T) 
wilcox.test(to.upper(controlMIC),to.upper(lupusMIC),paired=T) 

wilcox.test(to.upper(lupusCM^2),to.upper(lupusA),paired=T) 
wilcox.test(to.upper(lupusCM^2),to.upper(lupusDcor),paired=T) 
wilcox.test(to.upper(lupusCM^2),to.upper(lupusMIC),paired=T) 
wilcox.test(to.upper(lupusA),to.upper(lupusDcor),paired=T) 
wilcox.test(to.upper(lupusA),to.upper(lupusMIC),paired=T) 
wilcox.test(to.upper(lupusMIC),to.upper(lupusDcor),paired=) 

n = length(to.upper(controlCM^2))
binom.test(sum(to.upper(controlCM^2) >= to.upper(lupusCM^2)), n=n, alternative = 't')
binom.test(sum(to.upper(controlA) >= to.upper(lupusA)), n=n, alternative = 't')
binom.test(sum(to.upper(controlCM^2) >= to.upper(lupusCM^2)), n=n, alternative = 't')
binom.test(sum(to.upper(controlCM^2) >= to.upper(lupusCM^2)), n=n, alternative = 't')

binom.test(sum(to.upper(controlCM^2) >= to.upper(lupusCM^2)), n=n, alternative = 't')
binom.test(sum(to.upper(controlCM^2) >= to.upper(lupusCM^2)), n=n, alternative = 't')
binom.test(sum(to.upper(controlCM^2) >= to.upper(lupusCM^2)), n=n, alternative = 't')
binom.test(sum(to.upper(controlCM^2) >= to.upper(lupusCM^2)), n=n, alternative = 't')
binom.test(sum(to.upper(controlCM^2) >= to.upper(lupusCM^2)), n=n, alternative = 't')
binom.test(sum(to.upper(controlCM^2) >= to.upper(lupusCM^2)), n=n, alternative = 't')

median(to.upper(lupusCM^2))
median(to.upper(lupusA))
median(to.upper(lupusDcor))
median(to.upper(lupusMIC))


?wilcox.test

#copare rA v diff
cor(to.upper(lupusA.NL), to.upper(lupusA) - to.upper(lupusCM^2), method='spearman')
cor(to.upper(lupusDcor.NL), to.upper(lupusDcor) - to.upper(lupusCM^2), method='spearman')
cor(to.upper(lupusMIC.NL), to.upper(lupusMIC) - to.upper(lupusCM^2), method='spearman')

cor(round(to.upper(lupusA.NL), 2), round(to.upper(lupusA) - to.upper(lupusCM^2),2), method='kendall')
cor(round(to.upper(lupusDcor.NL), 2), round(to.upper(lupusDcor) - to.upper(lupusCM^2),2), method='kendall')
cor(to.upper(lupusMIC.NL), to.upper(lupusMIC) - to.upper(lupusCM^2), method='kendall')

cor(to.upper(lupusA) - to.upper(lupusCM^2), (to.upper(lupusA) - to.upper(lupusCM^2))/to.upper(lupusA), method='kendall')
cor(to.upper(lupusA) - to.upper(lupusCM^2), (to.upper(lupusA) - to.upper(lupusCM^2))/(1-to.upper(lupusCM^2)), method='kendall')

cor(to.upper(lupusDcor) - to.upper(lupusCM^2), (to.upper(lupusDcor) - to.upper(lupusCM^2))/to.upper(lupusDcor), method='kendall')
cor(to.upper(lupusDcor) - to.upper(lupusCM^2), (to.upper(lupusDcor) - to.upper(lupusCM^2))/(1-to.upper(lupusCM^2)), method='kendall')
cor(round(to.upper(lupusDcor) - to.upper(lupusCM^2),2), round((to.upper(lupusDcor) - to.upper(lupusCM^2))/to.upper(lupusDcor),2), method='kendall')
cor(to.upper(lupusDcor) - to.upper(lupusCM^2), (to.upper(lupusDcor) - to.upper(lupusCM^2))/(1-to.upper(lupusCM^2)), method='kendall')


cor(to.upper(lupusA) - to.upper(lupusCM^2), (to.upper(lupusDcor) - to.upper(lupusCM^2)), method='kendall')
cor(to.upper(lupusMIC) - to.upper(lupusCM^2), (to.upper(lupusDcor) - to.upper(lupusCM^2)), method='kendall')
cor(to.upper(lupusA) - to.upper(lupusCM^2), (to.upper(lupusMIC) - to.upper(lupusCM^2)), method='kendall')

cor(to.upper(lupusA) - to.upper(lupusCM^2), (to.upper(lupusDcor) - to.upper(lupusCM^2)), method='spearman')
cor(to.upper(lupusMIC) - to.upper(lupusCM^2), (to.upper(lupusDcor) - to.upper(lupusCM^2)), method='spearman')
cor(to.upper(lupusA) - to.upper(lupusCM^2), (to.upper(lupusMIC) - to.upper(lupusCM^2)), method='spearman')

#MARS
auc = function(res, linpreds=F)
{
  y = genes$diseased == 'Control'
  rr = genes[,-c(1,96)]
  aucs = nrow(res)
  for(i in 1:nrow(res))
  {
    aucs[i] = earth(y = y, x = rr[,c(res[i,1],res[i,2])], degree=2, linpreds=linpreds)$grsq
  }
  aucs
}
res
require(earth)
require(pROC)
max(aucDcor)
aucDcor = auc(lll, F)
grsq1 = data.frame(lll$dcor, auc(lll, F))
grsq1 = cbind(grsq1, rep('NonLinear', 18))
colnames(grsq1) = c('dcor', 'GRSq', 'Type')

grsq2 = data.frame(lll$dcor, auc(lll, T))
grsq2 = cbind(grsq2, rep('Linear', 18))
colnames(grsq2) = c('dcor', 'GRSq', 'Type')
grsq = rbind(grsq1, grsq2)

ggplot(grsq, aes(x=dcor, y=GRSq, colour=Type)) + geom_point(alpha=1, size=3) + scale_colour_brewer(palette="Set1")

par(mfrow=c(3,1))
hist(resultsNL$Value[which(resultsNL$Association == 'A')] - resultsNL$rValue[which(resultsNL$Association == 'A')])
hist(resultsNL$Value[which(resultsNL$Association == 'dcor')] - resultsNL$rValue[which(resultsNL$Association == 'dcor')])
hist(resultsNL$Value[which(resultsNL$Association == 'MIC')] - resultsNL$rValue[which(resultsNL$Association == 'MIC')])

length(resultsNL$Value[which(resultsNL$Association == 'A' & (resultsNL$Value - resultsNL$rValue <0)) ])/length(resultsNL$Value[which(resultsNL$Association == 'A')])
length(resultsNL$Value[which(resultsNL$Association == 'dcor' & (resultsNL$Value - resultsNL$rValue <0)) ])/length(resultsNL$Value[which(resultsNL$Association == 'dcor')])
length(resultsNL$Value[which(resultsNL$Association == 'MIC' & (resultsNL$Value - resultsNL$rValue <0)) ])/length(resultsNL$Value[which(resultsNL$Association == 'MIC')])

require(plyr)
arrange(resultsNL,desc(rValue))[1:50,]
cor(resultsNL$rValue[resultsNL$Association == 'A'], resultsNL$Difference[resultsNL$Association == 'A'], method='kendall')
cor(resultsNL$rValue[resultsNL$Association == 'dcor'], resultsNL$Difference[resultsNL$Association == 'dcor'], method='kendall')
cor(resultsNL$rValue[resultsNL$Association == 'MIC'], resultsNL$Difference[resultsNL$Association == 'MIC'], method='kendall')

cor(to.upper(controlA) - to.upper(controlCM^2), to.upper(controlA.NL))
cor(to.upper(controlDcor) - to.upper(controlCM^2), to.upper(controlDcor.NL))
cor(to.upper(controlMIC) - to.upper(controlCM^2), to.upper(controlMIC.NL))

cor(to.upper(lupusA) - to.upper(lupusCM^2), to.upper(lupusA.NL))
cor(to.upper(lupusDcor) - to.upper(lupusCM^2), to.upper(lupusDcor.NL))
cor(to.upper(lupusMIC) - to.upper(lupusCM^2), to.upper(lupusMIC.NL))

cor( to.upper(lupusCM^2), to.upper(lupusA.NL))
cor( to.upper(lupusCM^2), to.upper(lupusDcor.NL))
cor( to.upper(lupusCM^2), to.upper(lupusMIC.NL))

cor( to.upper(lupusA), to.upper(lupusA.NL))
cor(to.upper(lupusDcor) - to.upper(lupusCM^2), to.upper(lupusDcor.NL))
cor(to.upper(lupusMIC) - to.upper(lupusCM^2), to.upper(lupusMIC.NL))

plot(to.upper(lupusA) - to.upper(lupusCM^2), to.upper(lupusA.NL))
plot(to.upper(lupusDcor) - to.upper(lupusCM^2), to.upper(lupusDcor.NL))
plot(to.upper(lupusMIC) - to.upper(lupusCM^2), to.upper(lupusMIC.NL))

plot(to.upper(lupusCM^2), to.upper(lupusA.NL))
plot(to.upper(lupusCM^2), to.upper(lupusDcor.NL))
plot(to.upper(lupusCM^2), to.upper(lupusMIC.NL))

plot((to.upper(lupusCM^2) +to.upper(lupusA.NL))/2 , to.upper(lupusCM^2)- to.upper(lupusA.NL))
plot((to.upper(lupusCM^2) +to.upper(lupusDcor.NL))/2 , to.upper(lupusCM^2)- to.upper(lupusDcor.NL))
plot((to.upper(lupusCM^2) +to.upper(lupusMIC.NL))/2 , to.upper(lupusCM^2)- to.upper(lupusMIC.NL))

plot(((to.upper(lupusA) -to.upper(lupusCM^2)) +  /2 , to.upper(lupusCM^2)- to.upper(lupusA.NL))
plot((to.upper(lupusCM^2) +to.upper(lupusDcor.NL))/2 , to.upper(lupusCM^2)- to.upper(lupusDcor.NL))
plot((to.upper(lupusCM^2) +to.upper(lupusMIC.NL))/2 , to.upper(lupusCM^2)- to.upper(lupusMIC.NL))

which(lupusDcor.NL >= 0.3 & lupusDcor.NL < 1, arr.ind = T)

length(which(to.upper(lupusDcor) - to.upper(lupusDcor.NL) < 0))/4371
length(which(to.upper(lupusA) - to.upper(lupusA.NL) < 0))/4371
length(which(to.upper(lupusMIC) - to.upper(lupusMIC.NL) < 0))/4371
