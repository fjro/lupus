require(ggplot2)
require(reshape)
require(energy)
require(gridExtra)
library(plyr)
#plot pearson vs spearman for all pairwise combinations
to.upper<-function(X) X[upper.tri(X,diag=F)]
?upper.tri

ml = length(to.upper(controlCM))
corrs <- as.data.frame(to.upper(allCM^2))
corrs <- cbind(corrs, to.upper(controlCM^2))
corrs <- cbind(corrs, to.upper(lupusCM^2))
colnames(corrs) = c("All", "Control", "Lupus")
mc = melt(corrs)

corrs <- as.data.frame(to.upper(allA))
corrs <- cbind(corrs, to.upper(controlA))
corrs <- cbind(corrs, to.upper(lupusA))
colnames(corrs) = c("All", "Control", "Lupus")
mc = rbind(mc, melt(corrs))

corrs <- as.data.frame(to.upper(allDcor))
corrs <- cbind(corrs, to.upper(controlDcor))
corrs <- cbind(corrs, to.upper(lupusDcor))
colnames(corrs) = c("All", "Control", "Lupus")
mc = rbind(mc, melt(corrs))

corrs <- as.data.frame(to.upper(allMIC))
corrs <- cbind(corrs, to.upper(controlMIC))
corrs <- cbind(corrs, to.upper(lupusMIC))
colnames(corrs) = c("All", "Control", "Lupus")
mc = rbind(mc, melt(corrs))

mc = cbind(mc, c(rep('R2', ml*3), rep('A', ml*3), rep('dcor', ml*3), rep('MIC', ml*3)))
colnames(mc) <- c("Sample","Value", "Association")

#density plots
dp = ggplot(mc, aes(x=Value)) + geom_area(aes(group=Sample),stat="density")
dp + facet_grid(Association ~ Sample)
#dp + facet_grid(Sample ~ Association)

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

mdf = mdFrame(to.upper(controlCM^2), to.upper(lupusCM^2), "Pearson")
mdf = rbind(mdf, mdFrame(to.upper(controlA), to.upper(lupusA), "A"))
mdf = rbind(mdf, mdFrame(to.upper(controlDcor), to.upper(lupusDcor), "dcor"))
mdf = rbind(mdf, mdFrame(to.upper(controlMIC), to.upper(lupusMIC), "MIC"))
ggplot(mdf,aes(x=Mean, y=Difference)) + geom_point()+ facet_grid(. ~Association)


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
mdf = mdFrame2(to.upper(controlCM^2), to.upper(controlA), "R2 v A", "Control")
mdf = rbind(mdf, mdFrame2(to.upper(controlCM^2), to.upper(controlDcor), "R2 v dcor", "Control"))
mdf = rbind(mdf, mdFrame2(to.upper(controlCM^2), to.upper(controlMIC), "R2 v MIC", "Control"))
mdf = rbind(mdf, mdFrame2(to.upper(controlA), to.upper(controlDcor), "A v Dcor", "Control"))
mdf = rbind(mdf, mdFrame2(to.upper(controlA), to.upper(controlMIC), "A v MIC", "Control"))
mdf = rbind(mdf, mdFrame2(to.upper(controlDcor), to.upper(controlMIC), "dcor v MIC", "Control"))
mdf = rbind(mdf, mdFrame2(to.upper(controlA.NL), to.upper(controlDcor.NL), "rA v rDcor", "Control"))
mdf = rbind(mdf, mdFrame2(to.upper(controlA.NL), to.upper(controlMIC.NL), "rA v rMIC", "Control"))
mdf = rbind(mdf, mdFrame2(to.upper(controlDcor.NL), to.upper(controlMIC.NL), "rDcor v rMIC", "Control"))

mdf = rbind(mdf, mdFrame2(to.upper(lupusCM^2), to.upper(lupusA), "R2 v A", "Lupus"))
mdf = rbind(mdf, mdFrame2(to.upper(lupusCM^2), to.upper(lupusDcor), "R2 v dcor", "Lupus"))
mdf = rbind(mdf, mdFrame2(to.upper(lupusCM^2), to.upper(lupusMIC), "R2 v MIC", "Lupus"))
mdf = rbind(mdf, mdFrame2(to.upper(lupusA), to.upper(lupusDcor), "A v Dcor", "Lupus"))
mdf = rbind(mdf, mdFrame2(to.upper(lupusA), to.upper(lupusMIC), "A v MIC", "Lupus"))
mdf = rbind(mdf, mdFrame2(to.upper(lupusDcor), to.upper(lupusMIC), "dcor v MIC", "Lupus"))
mdf = rbind(mdf, mdFrame2(to.upper(lupusA.NL), to.upper(lupusDcor.NL), "rA v rDcor", "Lupus"))
mdf = rbind(mdf, mdFrame2(to.upper(lupusA.NL), to.upper(lupusMIC.NL), "rA v rMIC", "Lupus"))
mdf = rbind(mdf, mdFrame2(to.upper(lupusDcor.NL), to.upper(lupusMIC.NL), "rDcor v rMIC", "Lupus"))

#Bland ALtman plots
ggplot(mdf,aes(x=Mean, y=Difference)) + geom_point(alpha=0.7, size=1) + 
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
  cor(to.upper(x), to.upper(y), method = 'kendall')
}

#plot keddal results
kend = data.frame(kcor(controlCM^2, controlA))
kend = rbind(kend, kcor(controlCM^2, controlDcor))
kend = rbind(kend, kcor(controlCM^2, controlMIC))
kend = rbind(kend, kcor(controlA, controlDcor))
kend = rbind(kend, kcor(controlA, controlMIC))
kend = rbind(kend, kcor(controlDcor, controlMIC))
kend = rbind(kend, kcor(controlA.NL, controlDcor.NL))
kend = rbind(kend, kcor(controlA.NL, controlMIC.NL))
kend = rbind(kend, kcor(controlDcor.NL, controlMIC.NL))

kend = rbind(kend, kcor(lupusCM^2, lupusA))
kend = rbind(kend, kcor(lupusCM^2, lupusDcor))
kend = rbind(kend, kcor(lupusCM^2, lupusMIC))
kend = rbind(kend, kcor(lupusA, lupusDcor))
kend = rbind(kend, kcor(lupusA, lupusMIC))
kend = rbind(kend, kcor(lupusDcor, lupusMIC))
kend = rbind(kend, kcor(lupusA.NL, lupusDcor.NL))
kend = rbind(kend, kcor(lupusA.NL, lupusMIC.NL))
kend = rbind(kend, kcor(lupusDcor.NL, lupusMIC.NL))

kend = cbind(kend, c(rep("Control", 9), rep("Lupus", 9)))
pn = c("R2 v A", "R2 v dcor", "R2 v MIC", "A v dcor", "A v MIC", "dcor v MIC", "rA v rDcor", "rA v rMIC", "rDcor vs rMIC")
kend = cbind(kend, rep(pn, 2))
colnames(kend) = c("Tau", "Sample", "Association")

ggplot(kend,aes(x=Association, y=Tau)) + geom_bar(stat = "identity") + facet_grid(Sample ~ .) + ylab("Kendall's Tau") 

 
#plot results 0.99, 0.95, 0.9, 0.75
hc = data.frame(length(which(to.upper(controlCM^2) >= 0.99)))
hc = rbind(hc, length(  which(to.upper(controlCM^2) >= 0.95 & to.upper(controlCM^2) < 0.99) ))
hc = rbind(hc, length(which(to.upper(controlCM^2) >= 0.90 & to.upper(controlCM^2) < 0.95)))
hc = rbind(hc, length(which(to.upper(controlCM^2) >= 0.85 & to.upper(controlCM^2) < 0.90)))

hc = rbind(hc, length(which(to.upper(controlA) >= 0.99)))
hc = rbind(hc, length(which(to.upper(controlA) >= 0.95 & to.upper(controlA) < 0.99)))
hc = rbind(hc, length(which(to.upper(controlA) >= 0.90 & to.upper(controlA) < 0.95)))
hc = rbind(hc, length(which(to.upper(controlA) >= 0.85 & to.upper(controlA) < 0.90)))

hc = rbind(hc, length(which(to.upper(controlDcor) >= 0.99)))
hc = rbind(hc, length(which(to.upper(controlDcor) >= 0.95 & to.upper(controlDcor) < 0.99)))
hc = rbind(hc, length(which(to.upper(controlDcor) >= 0.90 & to.upper(controlDcor) < 0.95)))
hc = rbind(hc, length(which(to.upper(controlDcor) >= 0.85 & to.upper(controlDcor) < 0.90)))

hc = rbind(hc, length(which(to.upper(controlMIC) >= 0.99 )))
hc = rbind(hc, length(which(to.upper(controlMIC) >= 0.95 & to.upper(controlMIC) < 0.99)))
hc = rbind(hc, length(which(to.upper(controlMIC) >= 0.90 & to.upper(controlMIC) < 0.95)))
hc = rbind(hc, length(which(to.upper(controlMIC) >= 0.85 & to.upper(controlMIC) < 0.90)))

hc = rbind(hc, length(which(to.upper(lupusCM^2) >= 0.99)))
hc = rbind(hc, length(which(to.upper(lupusCM^2) >= 0.95 & to.upper(lupusCM^2) < 0.99)))
hc = rbind(hc, length(which(to.upper(lupusCM^2) >= 0.90 & to.upper(lupusCM^2) < 0.95)))
hc = rbind(hc, length(which(to.upper(lupusCM^2) >= 0.85 & to.upper(lupusCM^2) < 0.90)))

hc = rbind(hc, length(which(to.upper(lupusA) >= 0.99)))
hc = rbind(hc, length(which(to.upper(lupusA) >= 0.95 & to.upper(lupusA) < 0.99)))
hc = rbind(hc, length(which(to.upper(lupusA) >= 0.90 & to.upper(lupusA) < 0.95)))
hc = rbind(hc, length(which(to.upper(lupusA) >= 0.85 & to.upper(lupusA) < 0.90)))

hc = rbind(hc, length(which(to.upper(lupusDcor) >= 0.99)))
hc = rbind(hc, length(which(to.upper(lupusDcor) >= 0.95  & to.upper(lupusDcor) < 0.99)))
hc = rbind(hc, length(which(to.upper(lupusDcor) >= 0.90 & to.upper(lupusDcor) < 0.95)))
hc = rbind(hc, length(which(to.upper(lupusDcor) >= 0.85 & to.upper(lupusDcor) < 0.90)))

hc = rbind(hc, length(which(to.upper(lupusMIC) >= 0.99)))
hc = rbind(hc, length(which(to.upper(lupusMIC) >= 0.95 & to.upper(lupusMIC) < 0.99)))
hc = rbind(hc, length(which(to.upper(lupusMIC) >= 0.90 & to.upper(lupusMIC) < 0.95)))
hc = rbind(hc, length(which(to.upper(lupusMIC) >= 0.85 & to.upper(lupusMIC) < 0.90)))

hc = cbind(hc, rep(c(rep("R2", 4), rep("A", 4), rep("dcor", 4), rep("MIC", 4)), 2))
hc = cbind(hc, c(rep("Control", 16), rep("Lupus",16)))
hc = cbind(hc, rep(c("0.99", "0.95", "0.90", "0.85"), 4))
colnames(hc) = c('Count', "Association", "Sample", "Level")

ggplot(hc[hc$Level != c("0.85"),],aes(x=Association, y=Count, fill=Level)) + geom_bar(stat = "identity") + 
  scale_fill_brewer(type="seq") + facet_grid(Sample ~ .) + theme_bw() + 
  coord_flip() + theme(legend.position='bottom',legend.key.size = unit(.4, "cm")) 


p1 = ggplot(hc[hc$Level != c("0.90","0.85"),],aes(x=Association, y=Count, fill=Level)) + geom_bar(stat = "identity") + facet_grid(Sample ~ .) + coord_flip() + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
p2 = ggplot(hc,aes(x=Association, y=Count, fill=Level)) + geom_bar(stat = "identity") + facet_grid(Sample ~ .) + coord_flip() + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
grid.arrange(p1,p2, ncol=2, nrow=1,   main="", plot=T)

#set diff and union #########################

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
require(plyr)
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

#residual association
par(mfrow=c(2,3))
hist(to.upper(controlA.NL))
hist(to.upper(controlDcor.NL))
hist(to.upper(controlMIC.NL))
hist(to.upper(lupusA.NL))
hist(to.upper(lupusDcor.NL))
hist(to.upper(lupusMIC.NL))

vv = as.matrix(lupusA.NL)
vv[lower.tri(vv, diag = T)] <- 0
res = which(vv >= 0.5 & vv < 0.99, arr.ind = T)
vv[which(vv >= 0.5 & vv < 0.99)]
res$row
vv[res[,1],res[,2]]
which(lupusDcor.NL >= 0.30 & lupusDcor.NL < 0.9, arr.ind = T)
which(lupusMIC.NL >= 0.450 & lupusMIC.NL < 0.9, arr.ind = T)

idPairsNL = function(nlX, x, from, to, diseased, measure)
{ 
#   x = lupusA
#   nlX = lupusA.NL
#   from = 0.5
#   to = 0.6
#   diseased = 'Lupus'
#   measure = 'A'
  x = as.matrix(x)
  vv = as.matrix(nlX)
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
  df = cbind(df, x[which(vv >= from & vv < to)])
  df = cbind(df, rep(measure,n))
  
  colnames(df) = c('X','Y','Pair','XProbe','YProbe','From','To', 'Sample', 'rValue', 'Value', 'Association')
  
  arrange(df,desc(Value))
}



resultsNL = idPairsNL(controlA.NL, controlA, 0.5, 0.6, "Control", "A")
resultsNL = rbind(resultsNL, idPairsNL(controlA.NL, controlA, 0.45, 0.5, "Control", "A"))
resultsNL = rbind(resultsNL, idPairsNL(controlA.NL, controlA, 0.4, 0.45, "Control", "A"))
resultsNL = rbind(resultsNL, idPairsNL(controlA.NL, controlA, 0.35, 0.4, "Control", "A"))

resultsNL = rbind(resultsNL, idPairsNL(controlDcor.NL, controlDcor, 0.5, 0.6, "Control", "dcor"))
resultsNL = rbind(resultsNL, idPairsNL(controlDcor.NL, controlDcor, 0.45, 0.5, "Control", "dcor"))
resultsNL = rbind(resultsNL, idPairsNL(controlDcor.NL, controlDcor, 0.4, 0.45, "Control", "dcor"))
resultsNL = rbind(resultsNL, idPairsNL(controlDcor.NL, controlDcor, 0.35, 0.4, "Control", "dcor"))

resultsNL = rbind(resultsNL, idPairsNL(controlMIC.NL, controlMIC, 0.5, 0.6, "Control", "MIC"))
resultsNL = rbind(resultsNL, idPairsNL(controlMIC.NL, controlMIC, 0.45, 0.5, "Control", "MIC"))
resultsNL = rbind(resultsNL, idPairsNL(controlMIC.NL, controlMIC, 0.4, 0.45, "Control", "MIC"))
resultsNL = rbind(resultsNL, idPairsNL(controlMIC.NL, controlMIC, 0.35, 0.4, "Control", "MIC"))

resultsNL = rbind(resultsNL, idPairsNL(lupusA.NL, lupusA, 0.7, 0.9, "Lupus", "A"))
resultsNL = rbind(resultsNL, idPairsNL(lupusA.NL, lupusA, 0.6, 0.7, "Lupus", "A"))
resultsNL = rbind(resultsNL, idPairsNL(lupusA.NL, lupusA, 0.5, 0.6, "Lupus", "A"))
resultsNL = rbind(resultsNL, idPairsNL(lupusA.NL, lupusA, 0.45, 0.5, "Lupus", "A"))
resultsNL = rbind(resultsNL, idPairsNL(lupusA.NL, lupusA, 0.4, 0.45, "Lupus", "A"))
resultsNL = rbind(resultsNL, idPairsNL(lupusA.NL, lupusA, 0.35, 0.4, "Lupus", "A"))

resultsNL = rbind(resultsNL, idPairsNL(lupusDcor.NL, lupusDcor, 0.7, 0.9, "Lupus", "dcor"))
resultsNL = rbind(resultsNL, idPairsNL(lupusDcor.NL, lupusDcor, 0.6, 0.7, "Lupus", "dcor"))
resultsNL = rbind(resultsNL, idPairsNL(lupusDcor.NL, lupusDcor, 0.5, 0.6, "Lupus", "dcor"))
resultsNL = rbind(resultsNL, idPairsNL(lupusDcor.NL, lupusDcor, 0.45, 0.5, "Lupus", "dcor"))
resultsNL = rbind(resultsNL, idPairsNL(lupusDcor.NL, lupusDcor, 0.4, 0.45, "Lupus", "dcor"))
resultsNL = rbind(resultsNL, idPairsNL(lupusDcor.NL, lupusDcor, 0.35, 0.4, "Lupus", "dcor"))

resultsNL = rbind(resultsNL, idPairsNL(lupusMIC.NL, lupusMIC, 0.7, 0.9, "Lupus", "MIC"))
resultsNL = rbind(resultsNL, idPairsNL(lupusMIC.NL, lupusMIC, 0.6, 0.7, "Lupus", "MIC"))
resultsNL = rbind(resultsNL, idPairsNL(lupusMIC.NL, lupusMIC, 0.5, 0.6, "Lupus", "MIC"))
resultsNL = rbind(resultsNL, idPairsNL(lupusMIC.NL, lupusMIC, 0.45, 0.5, "Lupus", "MIC"))
resultsNL = rbind(resultsNL, idPairsNL(lupusMIC.NL, lupusMIC, 0.4, 0.45, "Lupus", "MIC"))
resultsNL = rbind(resultsNL, idPairsNL(lupusMIC.NL, lupusMIC, 0.35, 0.4, "Lupus", "MIC"))

resultsNL[which(resultsNL$Value >0.80),]
cc = controls[,-c(1,96)]
ww = resultsNL[resultsNL$Sample == "Control",1:2]
res = rep(0,nrow(ww))
for(i in 1:nrow(ww))
{
  res[i] = cor(cc[,ww[i,1]], cc[,ww[i,2]])
}

ll = lupus[,-c(1,96)]
ww2 = resultsNL[resultsNL$Sample == "Lupus",1:2]
res2 = rep(0,nrow(ww2))
for(i in 1:nrow(ww2))
{
  res2[i] = cor(ll[,ww2[i,1]], ll[,ww2[i,2]])
}
length(res2) + length(res)
dim(resultsNL)
resultsNL = cbind(resultsNL, c(res, res2))
resultsNL = cbind(resultsNL, c(res^2, res2^2))
colnames(resultsNL)[12:13] = c('Pearson', 'R2')

resultsNL[which(resultsNL$Value >0.70 & resultsNL$R2 < 0.70),]

resultsNL = cbind(resultsNL, resultsNL$Value - resultsNL$R2)
resultsNL = cbind(resultsNL, (resultsNL$Value - resultsNL$R2)/resultsNL$Value)
resultsNL = cbind(resultsNL, (resultsNL$Value - resultsNL$R2)/(1-resultsNL$R2))

colnames(resultsNL)[14:16] = c('Difference', 'nlp', 'ptv')

write.csv(resultsNL, "resultsNL.csv")

library(plyr)

#look at top results
head(arrange(resultsNL[resultsNL$Association == 'MIC',],desc(rValue)), n = 20)
head(arrange(resultsNL[resultsNL$Association == 'MIC',],desc(nlp)), n = 20)
head(arrange(resultsNL[resultsNL$Association == 'MIC',],desc(ptv)), n = 20)

plot(lupus[,8],lupus[,22])
par(mfrow=c(2,3))
hist(to.upper(controlA.NL) - to.upper(controlA))
hist(to.upper(controlDcor.NL) - to.upper(controlDcor))
hist(to.upper(controlMIC.NL) - to.upper(controlMIC))

hist(to.upper(lupusA.NL) - to.upper(lupusA))
hist(to.upper(lupusDcor.NL) - to.upper(lupusDcor))
hist(to.upper(lupusMIC.NL) - to.upper(lupusMIC))

View(results)

