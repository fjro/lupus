require(ggplot2)
require(reshape)
require(energy)
require(gridExtra)
library(plyr)

mc = melt(control3way)
colnames(mc) = c('Pair', 'Sample', 'Association', 'Value')
mc = mc[mc$Association %in% c('dcor', 'A', 'R2', 'nldcor', 'nlA'),]

lmc = melt(lupus3way)
colnames(lmc) = c('Pair', 'Sample', 'Association', 'Value')
lmc = lmc[lmc$Association %in% c('dcor', 'A', 'R2', 'nldcor', 'nlA'),]

mc = rbind(mc, lmc)
levels(mc$Association)
#density plots
dp = ggplot(mc, aes(x=Value)) + geom_area(aes(group=Sample),stat="density")
dp + facet_grid(Association ~ Sample)
#dp + facet_grid(Sample ~ Association)

levels(mc$Association)<- c('dcor' , "A","R2","rDcor","rA", "X" , "Y",  "Z")
#box
dp = ggplot(mc, aes(x=Sample, y=Value)) + geom_boxplot()
dp + facet_grid(. ~ Association)


#mean difference plots comparing lupus with contrls
#MA plot
mdFrame = function(x, y, association)
{
  badf = data.frame(x, y)
  badf = cbind(badf, (x+y)/2)
  badf = cbind(badf, x-y)
  badf = cbind(badf, rep(association, nrow(badf)))
  colnames(badf) = c('Control', 'Lupus', 'Mean', 'Difference', "Association")
  return(badf)
}

mdf = mdFrame(control3way$R2, lupus3way$R2, "R2")
mdf = rbind(mdf, mdFrame(control3way$A, lupus3way$A, "A"))
mdf = rbind(mdf, mdFrame(control3way$dcor, lupus3way$dcor, "dcor"))
mdf = rbind(mdf, mdFrame(control3way$nlA, lupus3way$nlA, "nlA"))
mdf = rbind(mdf, mdFrame(control3way$nldcor, lupus3way$nldcor, "nldcor"))
ggplot(mdf,aes(x=Mean, y=Difference)) + geom_point(alpha=0.5)+ facet_grid(Association ~ .)


mdFrame2 = function(x, y, association, sampleType)
{
  diff = x-y
  ll = mean(diff)-2*sd(diff)
  ul = mean(diff)+2*sd(diff)
  ml = mean(diff)
  badf = data.frame((x+y)/2, x-y)
  badf = cbind(rep(association, nrow(badf)), badf)
  badf = cbind(badf, rep(sampleType, nrow(badf)))
  badf = cbind(badf, rep(ll, nrow(badf)))
  badf = cbind(badf, rep(ul, nrow(badf)))
  badf = cbind(badf, rep(ml, nrow(badf)))
  colnames(badf) = c("Association", 'Mean', 'Difference', "Sample", "Lower", "Upper", "Middle")
  return(badf)
}
#compare measurement
mdf = mdFrame2(control3way$R2, control3way$A, "R2 v A", "Control")
mdf = rbind(mdf, mdFrame2(control3way$R2, control3way$dcor, "R2 v dcor", "Control"))
mdf = rbind(mdf, mdFrame2(control3way$A, control3way$dcor, "A v Dcor", "Control"))
mdf = rbind(mdf, mdFrame2(control3way$nlA, control3way$nldcor, "rA v rDcor", "Control"))

mdf = rbind(mdf, mdFrame2(lupus3way$R2, lupus3way$A, "R2 v A", "Lupus"))
mdf = rbind(mdf, mdFrame2(lupus3way$R2, lupus3way$dcor, "R2 v dcor", "Lupus"))
mdf = rbind(mdf, mdFrame2(lupus3way$A, lupus3way$dcor, "A v Dcor", "Lupus"))
mdf = rbind(mdf, mdFrame2(lupus3way$nlA, lupus3way$nldcor, "rA v rDcor", "Lupus"))

dim(mdf[mdf$Association == "nlA v nldcor" & mdf$Sample == "Control",])
#Bland ALtman plots
require(ggplot2)
ggplot(mdf,aes(x=Mean, y=Difference)) + geom_point(alpha=0.4, size=1) + 
  geom_hline(aes(yintercept=Lower, colour='Blue'))+
  geom_hline(aes(yintercept=Upper,colour='Blue'))+
  geom_hline(aes(yintercept=Middle,colour='Red'))+
  facet_grid(Association ~ Sample)

?cor

allCors = function(x, y)
{
  res = list('Pearson',  cor(to.upper(x), to.upper(y), method = 'pearson'),
             'Kendall', cor(to.upper(x), to.upper(y), method = 'kendall'),
             'Spearman',cor(to.upper(x), to.upper(y), method = 'spearman'), 
             'docr', dcor(to.upper(x), to.upper(y)))
  
  res
}

unlist(allCors(controlCM^2, controlA))
unlist(allCors(controlCM^2, controlDcor))
unlist(allCors(controlCM^2, controlMIC))
unlist(allCors(controlA, controlDcor))
unlist(allCors(controlA, controlMIC))
unlist(allCors(controlDcor, controlMIC))

unlist(allCors(lupusCM^2, lupusA))
unlist(allCors(lupusCM^2, lupusDcor))
unlist(allCors(lupusCM^2, lupusMIC))
unlist(allCors(lupusA, lupusDcor))
unlist(allCors(lupusA, lupusMIC))
unlist(allCors(lupusDcor, lupusMIC))

kcor = function(x, y)
{
  cor(x, y, method = 'kendall')
}

spear = function(x, y)
{
  cor(x, y, method = 'spearman')
}

#plot keddal results
kend = data.frame(kcor(control3way$R2, control3way$A))
kend = rbind(kend, kcor(control3way$R2, control3way$dcor))
kend = rbind(kend, kcor(control3way$A, control3way$dcor))
kend = rbind(kend, kcor(control3way$nlA, control3way$nldcor))

kend = rbind(kend, kcor(lupus3way$R2, lupus3way$A))
kend = rbind(kend, kcor(lupus3way$R2, lupus3way$dcor))
kend = rbind(kend, kcor(lupus3way$A, lupus3way$dcor))
kend = rbind(kend, kcor(lupus3way$nlA, lupus3way$nldcor))

kend = cbind(kend, c(rep("Control", 4), rep("Lupus", 4)))
pn = c("R2 v A", "R2 v dcor", "A v dcor", "rA v rDcor")
kend = cbind(kend, rep(pn, 2))
#kend$Association= rep(pn, 2)

colnames(kend) = c("Tau", "Sample", "Association")

require(ggplot2)
ggplot(kend,aes(x=Association, y=Tau)) + geom_bar(stat = "identity") + facet_grid(. ~ Sample) + ylab("Kendall's Tau") 

#plot keddal results
spearman = data.frame(spear(control3way$R2, control3way$A))
spearman = rbind(spearman, spear(control3way$R2, control3way$dcor))
spearman = rbind(spearman, spear(control3way$A, control3way$dcor))
spearman = rbind(spearman, spear(control3way$nlA, control3way$nldcor))

spearman = rbind(spearman, spear(lupus3way$R2, lupus3way$A))
spearman = rbind(spearman, spear(lupus3way$R2, lupus3way$dcor))
spearman = rbind(spearman, spear(lupus3way$A, lupus3way$dcor))
spearman = rbind(spearman, spear(lupus3way$nlA, lupus3way$nldcor))


#plot results 0.99, 0.95, 0.9, 0.75
dim(control3way)
hc = data.frame(length(which(control3way$R2 >= 0.99)))
hc = rbind(hc, length(which(control3way$R2 >= 0.95 & control3way$R2 < 0.99)))
hc = rbind(hc, length(which(control3way$R2 >= 0.90 & control3way$R2 < 0.95)))
hc = rbind(hc, length(which(control3way$R2 >= 0.85 & control3way$R2 < 0.90)))

hc = rbind(hc, length(which(control3way$A >= 0.99 )))
hc = rbind(hc, length(which(control3way$A >= 0.95 & control3way$A < 0.99)))
hc = rbind(hc, length(which(control3way$A >= 0.90 & control3way$A < 0.95)))
hc = rbind(hc, length(which(control3way$A >= 0.85 & control3way$A < 0.90)))

hc = rbind(hc, length(which(control3way$dcor >= 0.99 )))
hc = rbind(hc, length(which(control3way$dcor >= 0.95 & control3way$dcor < 0.99)))
hc = rbind(hc, length(which(control3way$dcor >= 0.90 & control3way$dcor < 0.95)))
hc = rbind(hc, length(which(control3way$dcor >= 0.85 & control3way$dcor < 0.90)))


hc = rbind(hc, length(which(lupus3way$R2 >= 0.99 )))
hc = rbind(hc, length(which(lupus3way$R2 >= 0.95 & lupus3way$R2 < 0.99)))
hc = rbind(hc, length(which(lupus3way$R2 >= 0.90 & lupus3way$R2 < 0.95)))
hc = rbind(hc, length(which(lupus3way$R2 >= 0.85 & lupus3way$R2 < 0.90)))

hc = rbind(hc, length(which(lupus3way$A >= 0.99 )))
hc = rbind(hc, length(which(lupus3way$A >= 0.95 & lupus3way$A < 0.99)))
hc = rbind(hc, length(which(lupus3way$A >= 0.90 & lupus3way$A < 0.95)))
hc = rbind(hc, length(which(lupus3way$A >= 0.85 & lupus3way$A < 0.90)))

hc = rbind(hc, length(which(lupus3way$dcor >= 0.99 )))
hc = rbind(hc, length(which(lupus3way$dcor >= 0.95 & lupus3way$dcor < 0.99)))
hc = rbind(hc, length(which(lupus3way$dcor >= 0.90 & lupus3way$dcor < 0.95)))
hc = rbind(hc, length(which(lupus3way$dcor >= 0.85 & lupus3way$dcor < 0.90)))

hc = cbind(hc, rep(c(rep("R2", 4), rep("A", 4), rep("dcor", 4)), 2))
hc = cbind(hc, c(rep("Control", 12), rep("Lupus",12)))
hc = cbind(hc, rep(c("0.99", "0.95", "0.90", "0.85"), 3))
colnames(hc) = c('Count', "Association", "Sample", "Level")

ggplot(hc[hc$Level != c("0.85"),],aes(x=Association, y=Count, fill=Level)) + geom_bar(stat = "identity") + 
  scale_fill_brewer(type="seq") + facet_grid(Sample ~ .) + theme_bw() + 
  coord_flip() + theme(legend.position='bottom',legend.key.size = unit(.4, "cm")) 


p1 = ggplot(hc[hc$Level != c("0.90","0.85"),],aes(x=Association, y=Count, fill=Level)) + geom_bar(stat = "identity") + facet_grid(Sample ~ .) + coord_flip() + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
p2 = ggplot(hc,aes(x=Association, y=Count, fill=Level)) + geom_bar(stat = "identity") + facet_grid(Sample ~ .) + coord_flip() + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
grid.arrange(p1,p2, ncol=2, nrow=1,   main="", plot=T)

#set diff and union #########################
r2ResC = control3way[control3way$R2 >= 0.90,9]
aResC = control3way[control3way$A >= 0.90,9]
dcorResC = control3way[control3way$dcor >= 0.90,9]

r2ResL = lupus3way[lupus3way$R2 >= 0.90,9]
aResL = lupus3way[lupus3way$A >= 0.90,9]
dcorResL = lupus3way[lupus3way$dcor >= 0.90,9]

require(gplots)
venn( list('R2'=r2ResC,'A'=aResC,'dcor'=dcorResC) )
venn( list('R2'=r2ResL,'A'=aResL,'dcor'=dcorResL) )






#controls first
vv = as.matrix(lupusA)
vv[lower.tri(vv, diag = T)] <- 0
res = which(vv >= 0.95 & vv < 0.99, arr.ind = T)
which(vv >= 0.95)
dim(res)
vv[which(as.list(vv) >= 0.95)]

typeof(controlA)
?cor.test
rownames(vv)[res[,1]]
colnames(vv)[res[,2]]
ddd = vv[res[,1],res[,2]]

idPairs = function(x, from, to, diseased, measure)
{
  vv = as.matrix(x)
  vv[lower.tri(vv, diag = T)] <- 0
  res = which(vv >= from & vv < to, arr.ind = T)
  df = data.frame(res[,1])
  df = cbind(df, res[,2])
  df = cbind(df, paste(res[,1], res[,2], sep=":"))
  df = cbind(df, rownames(vv)[res[,1]])
  df = cbind(df, colnames(vv)[res[,2]])
  n = nrow(df)
  df = cbind(df, rep(from, n))
  df = cbind(df, rep(to, n))
  df = cbind(df, rep(diseased, n))
  df = cbind(df, vv[which(vv >= from & vv < to)])
  df = cbind(df, rep(measure,n))
  
  colnames(df) = c('X','Y','Pair','XProbe','YProbe','From','To', 'Sample', 'Value', 'Association')
  
  arrange(df,desc(Value))
}

results = idPairs(controlCM^2, 0.99, 1, "Control", "R2")
results = rbind(results, idPairs(controlCM^2, 0.95, 0.99, "Control", "R2"))
results = rbind(results, idPairs(controlCM^2, 0.90, 0.95, "Control", "R2"))
results = rbind(results, idPairs(controlCM^2, 0.85, 0.90, "Control", "R2"))

results = rbind(results, idPairs(controlA, 0.99, 1, "Control", "A"))
results = rbind(results, idPairs(controlA, 0.95, 0.99, "Control", "A"))
results = rbind(results, idPairs(controlA, 0.90, 0.95, "Control", "A"))
results = rbind(results, idPairs(controlA, 0.85, 0.90, "Control", "A"))

results = rbind(results, idPairs(controlDcor, 0.99, 1, "Control", "dcor"))
results = rbind(results, idPairs(controlDcor, 0.95, 0.99, "Control", "dcor"))
results = rbind(results, idPairs(controlDcor, 0.90, 0.95, "Control", "dcor"))
results = rbind(results, idPairs(controlDcor, 0.85, 0.90, "Control", "dcor"))

results = rbind(results, idPairs(controlMIC, 0.99, 1, "Control", "MIC"))
results = rbind(results, idPairs(controlMIC, 0.95, 0.99, "Control", "MIC"))
results = rbind(results, idPairs(controlMIC, 0.90, 0.95, "Control", "MIC"))
results = rbind(results, idPairs(controlMIC, 0.85, 0.90, "Control", "MIC"))

results = rbind(results, idPairs(lupusCM^2, 0.99, 1, "Lupus", "R2"))
results = rbind(results, idPairs(lupusCM^2, 0.95, 0.99, "Lupus", "R2"))
results = rbind(results, idPairs(lupusCM^2, 0.90, 0.95, "Lupus", "R2"))
results = rbind(results, idPairs(lupusCM^2, 0.85, 0.90, "Lupus", "R2"))

results = rbind(results, idPairs(lupusA, 0.99, 1, "Lupus", "A"))
results = rbind(results, idPairs(lupusA, 0.95, 0.99, "Lupus", "A"))
results = rbind(results, idPairs(lupusA, 0.90, 0.95, "Lupus", "A"))
results = rbind(results, idPairs(lupusA, 0.85, 0.90, "Lupus", "A"))

results = rbind(results, idPairs(lupusDcor, 0.99, 1, "Lupus", "dcor"))
results = rbind(results, idPairs(lupusDcor, 0.95, 0.99, "Lupus", "dcor"))
results = rbind(results, idPairs(lupusDcor, 0.90, 0.95, "Lupus", "dcor"))
results = rbind(results, idPairs(lupusDcor, 0.85, 0.90, "Lupus", "dcor"))

results = rbind(results, idPairs(lupusMIC, 0.99, 1, "Lupus", "MIC"))
results = rbind(results, idPairs(lupusMIC, 0.95, 0.99, "Lupus", "MIC"))
results = rbind(results, idPairs(lupusMIC, 0.90, 0.95, "Lupus", "MIC"))
results = rbind(results, idPairs(lupusMIC, 0.85, 0.90, "Lupus", "MIC"))

write.csv(results, "results.csv")
table(results$Association)

results = results[results$Value >= 0.90,]
table(results$Sample, results$Association)

#calculate AUC
res = results[results$Sample == 'Lupus' & results$Association == 'dcor' & results$From >=0.95,]

auc3 = function(res, linpreds=F)
{
  y = genes$diseased == 'Control'
  rr = genes[,-c(1,96)]
  aucs = rep(0,nrow(res))
  for(i in 1:nrow(res))
  {
    aucs[i] = earth(y = y, x = rr[,c(res$X[i],res$Y[i], res$Z[i])], degree=3, linpreds=linpreds)$grsq
  }

  aucs
}
require(earth)
require(pROC)
require(ggplot2)
dim()


results3R2 = lupus3way[lupus3way$R2 >= 0.90,3]
results3R2 = cbind(results3R2, auc3(lupus3way[lupus3way$R2 >= 0.90,]))
results3R2 = cbind(results3R2, auc3(lupus3way[lupus3way$R2 >= 0.90,], T))
results3R2 = cbind(results3R2, rep('R2', nrow(results3R2)))
colnames(results3R2) = c("Value", "AUC", "Linear AUC", "Association")

results3A = lupus3way[lupus3way$A >= 0.90,2]
results3A = cbind(results3A, auc3(lupus3way[lupus3way$A >= 0.90,]))
results3A = cbind(results3A, auc3(lupus3way[lupus3way$A >= 0.90,], T))
results3A = cbind(results3A, rep('A', nrow(results3A)))
colnames(results3A) = c("Value", "AUC", "Linear AUC", "Association")

results3Dcor = lll3
results3Dcor = cbind(results3Dcor, auc3(lll3))
results3Dcor = cbind(results3Dcor, auc3(lll3, T))
results3Dcor = cbind(results3Dcor, rep('dcor', nrow(results3Dcor)))
colnames(results3Dcor) = c("Value", "AUC", "Linear AUC", "Association")

results3 = rbind(results3R2, results3A)
results3 = rbind(results3, results3Dcor)
results3 = as.data.frame(results3)

typeof(results3$AUC)
typeof(results3$Value)
results3$AUC = c(as.double(results3R2[,2]), as.double(results3A[,2]), as.double(results3Dcor[,2]))
results3$`Linear AUC` = c(as.double(results3R2[,3]), as.double(results3A[,3]), as.double(results3Dcor[,3]))
results3$Value = c(as.double(results3R2[,1]), as.double(results3A[,1]), as.double(results3Dcor[,1]))
plot(results3[results3$Association == 'R2',1], results3[results3$Association == 'R2',2])

dp = ggplot(results3, aes(x=Value, y=AUC)) + geom_point(alpha=0.5, size=1) + facet_grid(Association ~ .)

bp = ggplot(results3, aes(x=Association, y=AUC)) + geom_boxplot()  + facet_grid(Association ~ .) + ylab("") +theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())
bp





require(gridExtra)

grid.arrange(dp, bp, ncol=2, nrow=1, widths=3:1)

write.csv(results3, "auc3way.csv")

hist(results3$AUC - results3$`Linear AUC`)
mean(results3$AUC - results3$`Linear AUC`)
sd(results3$AUC - results3$`Linear AUC`)
sum(results3$AUC - results3$`Linear AUC` < 0)


grsq1 = data.frame(lll3$dcor, auc3(lll3, F))
grsq1 = cbind(grsq1, rep('NonLinear', nrow(grsq1)))
colnames(grsq1) = c('dcor', 'GRSq', 'Type')

grsq2 = data.frame(lll3$dcor, auc3(lll3, T))
grsq2 = cbind(grsq2, rep('Linear', nrow(grsq1)))
colnames(grsq2) = c('dcor', 'GRSq', 'Type')
grsq = rbind(grsq1, grsq2)

ggplot(grsq, aes(x=dcor, y=GRSq, colour=Type)) + geom_point(alpha=1, size=1) + scale_colour_brewer(palette="Set1")

mean(grsq1$GRSq)
mean(grsq1$GRSq - grsq2$GRSq)
sd(grsq1$GRSq - grsq2$GRSq)
sum(grsq1$GRSq < grsq2$GRSq)

#look at non-linear effects £##################################
library(plyr)

#look at top results
colnames(lupus3way)

lupus3way = cbind(lupus3way, lupus3way$dcor - lupus3way$R2)
lupus3way = cbind(lupus3way, lupus3way$A - lupus3way$R2)

lupus3way = cbind(lupus3way, (lupus3way$dcor - lupus3way$R2)/lupus3way$dcor)
lupus3way = cbind(lupus3way, (lupus3way$A - lupus3way$R2)/lupus3way$A)


lupus3way = cbind(lupus3way, (lupus3way$dcor - lupus3way$R2)/(1 - lupus3way$R2))
lupus3way = cbind(lupus3way, (lupus3way$A - lupus3way$R2)/(1 - lupus3way$R2))


max(lupus3way$nldcor/lupus3way$dcor)

colnames(lupus3way)[11:16] = c('diffDcor', 'diffA', 'nlpDcor', 'nlpA', 'ptvDcor', 'ptvA')

lupus3way = cbind(lupus3way, lupus3way$nldcor/lupus3way$dcor)
lupus3way = cbind(lupus3way, lupus3way$nlA/lupus3way$A)

colnames(lupus3way)[17:18] = c('rrDcor', 'rrA')

head(arrange(lupus3way,desc(nlA)), n = 20)
head(arrange(lupus3way,desc(diffA)), n = 20)
head(arrange(lupus3way,desc(nlpA)), n = 20)
head(arrange(lupus3way,desc(ptvA)), n = 20)

head(arrange(lupus3way,desc(nldcor)), n = 20)
head(arrange(lupus3way,desc(diffDcor)), n = 20)
head(arrange(lupus3way,desc(nlpDcor)), n = 20)
head(arrange(lupus3way,desc(ptvDcor)), n = 20)

head(arrange(lupus3way,desc(rrDcor)), n = 20)
head(arrange(lupus3way,desc(rrA)), n = 20)

plot(lupus3way$nlA, lupus3way$nlDcor)

head(arrange(lupus3way,desc(diffDcor)), n = 20)

head(arrange(resultsNL[resultsNL$Association == 'MIC',],desc(nlp)), n = 20)
head(arrange(resultsNL[resultsNL$Association == 'MIC',],desc(ptv)), n = 20)

#no obvious 3 ways
length(which(lupus3way$dcor >= 0.95 & lupus3way $R2 < 0.95))

suX = lupus3way[lupus3way$dcor >= 0.99 & lupus3way$R2 < 0.90,6]
suY = lupus3way[lupus3way$dcor >= 0.99 & lupus3way$R2 < 0.90,7]
suZ = lupus3way[lupus3way$dcor >= 0.99 & lupus3way$R2 < 0.90,8]

suX = lupus3way[lupus3way$dcor >= 0.99,6]
suY = lupus3way[lupus3way$dcor >= 0.99,7]
suZ = lupus3way[lupus3way$dcor >= 0.99,8]
su3 = union(suX, suY)
su3 = union(su3, suZ)
length(su3)


head(arrange(lupus3way[lupus3way$dcor >= 0.95,],desc(diffDcor)), n = 50)
mean(lupus3way[lupus3way$dcor >= 0.95,11])
sd(lupus3way[lupus3way$dcor >= 0.95,11])

res = lupus3way[lupus3way$dcor >= 0.989,c(6,7,8)]
res$X = as.factor(res$X)
res$Y = as.factor(res$Y)
res$Z = as.factor(res$Z)
rownames(res) = 1:nrow(res)
mosaicplot(res, main="")
hist(lupus3way[lupus3way$dcor >= 0.95 & lupus3way$diffDcor >=0.05 & !lupus3way$X %in% c(90,78,38),])

table(lupus3way[lupus3way$dcor >= 0.95 & lupus3way$diffDcor >=0.05,6])
table(lupus3way[lupus3way$dcor >= 0.95 & lupus3way$diffDcor >=0.05,7])
table(lupus3way[lupus3way$dcor >= 0.95 & lupus3way$diffDcor >=0.05,8])


lupus3way[lupus3way$dcor >= 0.989,]
require(scatterplot3d)
require(rafalib)

par(mfrow=c(2,2))
scatterplot3d(lupus[,91],   # x axis
              lupus[,20],     # y axis
              lupus[,78],    # z axis
              xlab = "90",
              ylab = "19",
              zlab = "78",
              main="")
scatterplot3d(lupus[,79],   # x axis
              lupus[,39],     # y axis
              lupus[,91],    # z axis
              xlab = "78",
              ylab = "38",
              zlab = "90",
              main="")
scatterplot3d(lupus[,91],   # x axis
              lupus[,44],     # y axis
              lupus[,79],    # z axis
              xlab = "90",
              ylab = "43",
              zlab = "79",
              main="")
scatterplot3d(lupus[,91],   # x axis
              lupus[,51],     # y axis
              lupus[,79],    # z axis
              xlab = "90",
              ylab = "19",
              zlab = "78",
              main="")
l2 = lupus[,-c(1,96)]
plot(l2[,90], lupus[,19])
plot(l2[,90], lupus[,78])
plot(l2[,78], lupus[,19])
cor(l2[,90], lupus[,19])
cor(l2[,90], lupus[,78])
cor(l2[,78], lupus[,19])



require(GGally)


l2 = l2[,c(90,19,78,38,43,50)]
colnames(l2) = c("A","B","C","D","E","F")
ggpairs(l2)
plot(l2[,90], l2[,19])
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

require(energy)
dcor(lupus[,4], lupus[,5:6])
dcor(lupus[,6], lupus[,4:5])
dcor(lupus[,5], lupus[,c(4,6)])

lupus3way = cbind(lupus3way, lupus3way$dcor - lupus3way$R2)
lupus3way = cbind(lupus3way, lupus3way$dcor - control3way$dcor)
colnames(lupus3way)[11:12] = c('betaDcor','gammaDcor')

#which(vv >0.45 & d2 >0.90 & (d2 - r2) > 0.1, arr.ind = T)
require(plyr)
bothProbes = c(30,31,86,76,42,41)
lll3 = arrange(lupus3way[lupus3way$dcor >=0.9 & lupus3way$betaDcor >= 0.1 & 
                           lupus3way$gammaDcor >=0.45 &
                           !(lupus3way$X %in% bothProbes) &
                           !(lupus3way$Y %in% bothProbes) &
                           !(lupus3way$Z %in% bothProbes),] ,desc(betaDcor))
su3 = union(lll3$X, lll3$Z)
yy = unique(lll3$Y)
intersect(yy, su3)



lll4 = arrange(lupus3way[lupus3way$dcor >=0.9 & lupus3way$betaDcor >= 0.1 & 
                           lupus3way$gammaDcor >=0.45 & !(lupus3way$X %in% lll$X) & 
                           !(lupus3way$Z %in% lll$Y) &
                           !(lupus3way$X %in% bothProbes) &
                           !(lupus3way$Y %in% bothProbes) &
                           !(lupus3way$Z %in% bothProbes),] ,desc(betaDcor))


lll5 = arrange(lupus3way[lupus3way$dcor >=0.9 & lupus3way$betaDcor >= 0.1 & 
                           lupus3way$gammaDcor >=0.45 & 
                          ( lupus3way$X %in% de) & 
                           (lupus3way$Z %in% de) &
                           !(lupus3way$X %in% bothProbes) &
                           !(lupus3way$Y %in% bothProbes) &
                           !(lupus3way$Z %in% bothProbes) &
                           lupus3way$Y %in% de,] ,desc(betaDcor))

lll6 = arrange(lupus3way[lupus3way$dcor >=0.9 & lupus3way$betaDcor >= 0.1 & 
                           lupus3way$gammaDcor >=0.45 & !(lupus3way$X %in% lll$X) & 
                           !(lupus3way$Z %in% lll$Y) &
                           !(lupus3way$X %in% bothProbes) &
                           !(lupus3way$Y %in% bothProbes) &
                           !(lupus3way$Z %in% bothProbes) &
                           lupus3way$Y %in% de,] ,desc(betaDcor))

par(mfrow=c(2,2))
plot(table(c(lll3[,6],lll3[,7],lll3[,8])), xlab = 'Probe No', ylab='Count',main='A', xlim=c(0,96))
plot(table(c(lll4[,6],lll4[,7],lll4[,8])), xlab = 'Probe No', ylab='Count',main='B', xlim=c(0,96))
plot(table(c(lll5[,6],lll5[,7],lll5[,8])), xlab = 'Probe No', ylab='Count',main='C', xlim=c(0,96))
plot(table(c(lll6[,6],lll6[,7],lll6[,8])), xlab = 'Probe No', ylab='Count',main='D', xlim=c(0,96))

de
unique(lll4$X)
unique(lll4$Z)
unique(lll4$Y)
su4 = union(lll4$X, lll4$Z)
yy = unique(lll4$Y)
length(yy)
union(su3, yy)
plot(lupus[,lll5$X[1] + 1], lupus[,lll5$Y[1] + 1])
plot(lupus[,lll5$X[1] + 1], lupus[,lll5$Z[1] + 1])
plot(lupus[,lll5$Y[1] + 1], lupus[,lll5$Z[1] + 1])

scatterplot3d(lupus[,lll5$X[1] + 1],   # x axis
              lupus[,lll5$Y[1] + 1],     # y axis
              lupus[,lll5$Z[1] + 1],    # z axis
              xlab = "81",
              ylab = "33",
              zlab = "74",
              main="")

my.col <- colorRampPalette(brewer.pal(11, "RdBu"))(diff(range(dat$z)))
xyplot(Y~ X, data=ff, col=my.col, pch=19, alpha=.5)

scatterplot3d(lupus[,lll4$X[2] + 1],   # x axis
              lupus[,lll4$Y[2] + 1],     # y axis
              lupus[,lll4$Z[2] + 1],    # z axis
              xlab = "90",
              ylab = "19",
              zlab = "78",
              main="")

scatterplot3d(lupus[,lll4$X[3] + 1],   # x axis
              lupus[,lll4$Y[3] + 1],     # y axis
              lupus[,lll4$Z[3] + 1],    # z axis
              xlab = "90",
              ylab = "19",
              zlab = "78",
              main="")

scatterplot3d(lupus[,lll4$X[168] + 1],   # x axis
              lupus[,lll4$Y[168] + 1],     # y axis
              lupus[,lll4$Z[168] + 1],    # z axis
              xlab = "90",
              ylab = "19",
              zlab = "78",
              main="")

require(lattice)
require(RColorBrewer)
?contourplot
ff = data.frame(lupus[,lll5$X[2] +1], lupus[,lll5$Y[1] +1], lupus[,lll5$Z[1] +1])
colnames(ff) = c('X', 'Y', 'Z')
my.col <- colorRampPalette(brewer.pal(11, "RdBu"))(diff(range(ff$X)))
levelplot(X ~ Y * Z, data = ff,region = TRUE, 
          col.regions = terrain.colors(100),
            xlab = "94",
            ylab = "Temperature (F)",
            main = "Cube Root Ozone (cube root ppb)")
?level.colors

xy = rep(0, nrow(lll5))
xz = rep(0, nrow(lll5))
zy = rep(0, nrow(lll5))
for(i in 1:nrow(lll5))
{
  xy[i] = lupusDcor[lll5$X[i],lll5$Y[i]]
  xz[i] = lupusDcor[lll5$X[i],lll5$Z[i]]
  zy[i] = lupusDcor[lll5$Z[i],lll5$Y[i]]
}

lll5 = cbind(lll5, xy)
lll5 = cbind(lll5, xz)
lll5 = cbind(lll5, zy)


require(energy)
length(lupus[,75])
length(runif(367))

x = 
  dcor(lupus[,82], lupus[,75])
dcor(lupus[,82], cbind(lupus[,75], runif(367)))
dcor(runif(367), cbind(lupus[,75], lupus[,82]))
res = rep(F, 1000)
for(i in 1:1000)
{
  res[i] = dcor(lupus[,82], cbind(lupus[,75], runif(367))) > x
}
ff = data.frame(lupus[,lll5$X[3] +1], lupus[,lll5$Y[3] +1], lupus[,lll5$Z[3] +1])
colnames(ff) = c('X', 'Y', 'Z')
elevation.loess = loess(Z ~ X*Y, data = ff, degree = 1, span = 0.25)
elevation.fit = expand.grid(list(X = seq(min(ff$X), max(ff$X), 0.2), Y =  seq(min(ff$Y), max(ff$Y), 0.2)))
Z = predict(elevation.loess, newdata = elevation.fit)
elevation.fit$Height = as.numeric(Z)
min(elevation.fit$Height)

contourplot(Height ~ X*Y, data = elevation.fit,
         xlab = "X Coordinate (feet)", ylab = "Y Coordinate (feet)",
          main = "Surface elevation dacta",
          col.regions = terrain.colors(100))
par(mfrow=c(2,2))

require(ggplot2)
ggplot(elevation.fit, aes(X, Y, fill = Height)) + geom_tile() +
  xlab("X Coordinate (feet)") + ylab("Y Coordinate (feet)") +
  scale_fill_gradient(limits = c(min(ff$Z), max(ff$Z)), low = "black", high = "white") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))
plot(lupus[,54], lupus[,44])
par(mfrow=c(4,3))
for(i in 1:12)
{
  plot(lupus[,lll5$X[i] + 1], lupus[,lll5$Y[i] + 1])
  plot(lupus[,lll5$X[i] + 1], lupus[,lll5$Z[i] + 1])
  plot(lupus[,lll5$Y[i] + 1], lupus[,lll5$Z[i] + 1])
}

median(lupus3way$dcor)
median(lupus3way$R2)
median(lupus3way$A)
