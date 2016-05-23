
require(ggplot2)
require(reshape)
require(energy)
require(gridExtra)
library(plyr)

alienBD2wayControl = controlBD2way
alienBD2wayLupus = lupusBD2way

mac1 = read.csv('controlBD2WayMAC.csv')
mac2 = read.csv('lupusBD2WayMAC.csv')
dim(mac1)
dim(mac2)
dim(controlBD2way)
controlBD2way = mac1[,-1]
lupusBD2way = mac2[,-1]

mcb = controlBD2way[,4]
mcb = cbind(mcb, rep('dcor', nrow(controlBD2way)))
colnames(mcb) = c('Value','Association')

mcb1 = controlBD2way[,5]
mcb1 = cbind(mcb1, rep('A', nrow(controlBD2way)))
colnames(mcb1) = c('Value','Association')

mcb2 = controlBD2way[,6]
mcb2 = cbind(mcb2, rep('R2', nrow(controlBD2way)))
colnames(mcb2) = c('Value','Association')

mcb3 = controlBD2way[,7]
mcb3 = cbind(mcb3, rep('MIC', nrow(controlBD2way)))
colnames(mcb3)= c('Value','Association')
mcb= rbind(mcb, mcb1)
mcb= rbind(mcb, mcb2)
mcb= rbind(mcb, mcb3)

##
mcb4 = lupusBD2way[,4]
mcb4 = cbind(mcb4, rep('dcor', nrow(controlBD2way)))
colnames(mcb4) = c('Value','Association')

mcb5 = lupusBD2way[,5]
mcb5 = cbind(mcb5, rep('A', nrow(controlBD2way)))
colnames(mcb5) = c('Value','Association')

mcb6 = lupusBD2way[,6]
mcb6 = cbind(mcb6, rep('R2', nrow(controlBD2way)))
colnames(mcb6) = c('Value','Association')

mcb7 = lupusBD2way[,7]
mcb7 = cbind(mcb7, rep('MIC', nrow(controlBD2way)))
colnames(mcb7) = c('Value','Association')

mcb= rbind(mcb, mcb4)
mcb= rbind(mcb, mcb5)
mcb= rbind(mcb, mcb6)
mcb = rbind(mcb, mcb7)
dim(mcb)
mcb = cbind(mcb, c(rep('Control', nrow(mcb)/2), rep('Lupus', nrow(mcb)/2)))
colnames(mcb)[3] = 'Sample'

mcb = data.frame(mcb)
#density plots
require(ggplot2)
#dp = ggplot(mcb, aes(x=Value)) + geom_area(aes(group=Sample),stat="density")
#dp + facet_grid(Association ~ Sample)
#dp + facet_grid(Sample ~ Association)


#box
dp = ggplot(mcb, aes(x=Sample, y=Value)) + geom_boxplot()
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

mdf = mdFrame(controlBD2way$R2, lupus3way$R2, "R2")
mdf = rbind(mdf, mdFrame(control3way$A, lupus3way$A, "A"))
mdf = rbind(mdf, mdFrame(control3way$dcor, lupus3way$dcor, "dcor"))
mdf = rbind(mdf, mdFrame(control3way$nlA, lupus3way$nlA, "nlA"))
mdf = rbind(mdf, mdFrame(control3way$nldcor, lupus3way$nldcor, "nldcor"))
ggplot(mdf,aes(x=Mean, y=Difference)) + geom_point(alpha=0.5)+ facet_grid(Association ~ .)



#compare measurement
mdf = mdFrame2(controlBD2way$R2, controlBD2way$A, "R2 v A", "Control")
mdf = rbind(mdf, mdFrame2(controlBD2way$R2, controlBD2way$dcor, "R2 v dcor", "Control"))
mdf = rbind(mdf, mdFrame2(controlBD2way$R2, controlBD2way$MIC, "R2 v MIC", "Control"))
mdf = rbind(mdf, mdFrame2(controlBD2way$A, controlBD2way$dcor, "A v dcor", "Control"))
mdf = rbind(mdf, mdFrame2(controlBD2way$A, controlBD2way$MIC, "A v MIC", "Control"))
mdf = rbind(mdf, mdFrame2(controlBD2way$dcor, controlBD2way$MIC, "dcor v MIC", "Control"))

mdf = rbind(mdf, mdFrame2(lupusBD2way$R2, lupusBD2way$A, "R2 v A", "Lupus"))
mdf = rbind(mdf, mdFrame2(lupusBD2way$R2, lupusBD2way$dcor, "R2 v dcor", "Lupus"))
mdf = rbind(mdf, mdFrame2(lupusBD2way$R2, lupusBD2way$MIC, "R2 v MIC", "Lupus"))
mdf = rbind(mdf, mdFrame2(lupusBD2way$A, lupusBD2way$dcor, "A v dcor", "Lupus"))
mdf = rbind(mdf, mdFrame2(lupusBD2way$A, lupusBD2way$MIC, "A v MIC", "Lupus"))
mdf = rbind(mdf, mdFrame2(lupusBD2way$dcor, lupusBD2way$MIC, "dcor v MIC", "Lupus"))


#Bland ALtman plots
require(ggplot2)
ggplot(mdf,aes(x=Mean, y=Difference)) + geom_point(alpha=0.4, size=1) + 
  geom_hline(aes(yintercept=Lower, colour='Blue'))+
  geom_hline(aes(yintercept=Upper,colour='Blue'))+
  geom_hline(aes(yintercept=Middle,colour='Red'))+
  facet_grid(Association ~ Sample)


kcor = function(x, y)
{
  cor(x, y, method = 'kendall')
}

spear = function(x, y)
{
  cor(x, y, method = 'spearman')
}

dim(lupus)

#plot keddal results
kend = data.frame(kcor(controlBD2way$R2, controlBD2way$A))
kend = rbind(kend, kcor(controlBD2way$R2, controlBD2way$dcor))
kend = rbind(kend, kcor(controlBD2way$R2, controlBD2way$MIC))
kend = rbind(kend, kcor(controlBD2way$A, controlBD2way$dcor))
kend = rbind(kend, kcor(controlBD2way$A, controlBD2way$MIC))
kend = rbind(kend, kcor(controlBD2way$dcor, controlBD2way$MIC))

kend = rbind(kend, kcor(lupusBD2way$R2, lupusBD2way$A))
kend = rbind(kend, kcor(lupusBD2way$R2, lupusBD2way$dcor))
kend = rbind(kend, kcor(lupusBD2way$R2, lupusBD2way$MIC))
kend = rbind(kend, kcor(lupusBD2way$A, lupusBD2way$dcor))
kend = rbind(kend, kcor(lupusBD2way$A, lupusBD2way$MIC))
kend = rbind(kend, kcor(lupusBD2way$dcor, lupusBD2way$MIC))

kend = cbind(kend, c(rep("Control", 6), rep("Lupus", 6)))
pn = c("R2 v A", "R2 v dcor", "R2 v MIC", "A v dcor", "A v MIC", "dcor v MIC")
kend = cbind(kend, rep(pn, 2))
#kend$Association= rep(pn, 2)

colnames(kend) = c("Tau", "Sample", "Association")

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

#performance
require(energy)
require(minerva)
require(matie)
x = rnorm(320)
y = rnorm(320)
system.time( replicate(10000, cor(x,y)^2 ) )
system.time( replicate(10000, dcor(x,y) ) )
system.time( replicate(10000, ma(data.frame(x,y)) ) )
system.time( replicate(10000, mine(x,y) ) )

#plot results 0.99, 0.95, 0.9, 0.75
dim(control3way)
hc = data.frame(length(which(controlBD2way$R2 >= 0.99)))
hc = rbind(hc, length(which(controlBD2way$R2 >= 0.95 & controlBD2way$R2 < 0.99)))
hc = rbind(hc, length(which(controlBD2way$R2 >= 0.90 & controlBD2way$R2 < 0.95)))
hc = rbind(hc, length(which(controlBD2way$R2 >= 0.85 & controlBD2way$R2 < 0.90)))

hc = rbind(hc, length(which(controlBD2way$A >= 0.99 )))
hc = rbind(hc, length(which(controlBD2way$A >= 0.95 & controlBD2way$A < 0.99)))
hc = rbind(hc, length(which(controlBD2way$A >= 0.90 & controlBD2way$A < 0.95)))
hc = rbind(hc, length(which(controlBD2way$A >= 0.85 & controlBD2way$A < 0.90)))

hc = rbind(hc, length(which(controlBD2way$dcor >= 0.99 )))
hc = rbind(hc, length(which(controlBD2way$dcor >= 0.95 & controlBD2way$dcor < 0.99)))
hc = rbind(hc, length(which(controlBD2way$dcor >= 0.90 & controlBD2way$dcor < 0.95)))
hc = rbind(hc, length(which(controlBD2way$dcor >= 0.85 & controlBD2way$dcor < 0.90)))

hc = rbind(hc, length(which(controlBD2way$MIC >= 0.99 )))
hc = rbind(hc, length(which(controlBD2way$MIC >= 0.95 & controlBD2way$MIC < 0.99)))
hc = rbind(hc, length(which(controlBD2way$MIC >= 0.90 & controlBD2way$MIC < 0.95)))
hc = rbind(hc, length(which(controlBD2way$MIC >= 0.85 & controlBD2way$MIC < 0.90)))


hc = rbind(hc, length(which(lupusBD2way$R2 >= 0.99 )))
hc = rbind(hc, length(which(lupusBD2way$R2 >= 0.95 & lupusBD2way$R2 < 0.99)))
hc = rbind(hc, length(which(lupusBD2way$R2 >= 0.90 & lupusBD2way$R2 < 0.95)))
hc = rbind(hc, length(which(lupusBD2way$R2 >= 0.85 & lupusBD2way$R2 < 0.90)))

hc = rbind(hc, length(which(lupusBD2way$A >= 0.99 )))
hc = rbind(hc, length(which(lupusBD2way$A >= 0.95 & lupusBD2way$A < 0.99)))
hc = rbind(hc, length(which(lupusBD2way$A >= 0.90 & lupusBD2way$A < 0.95)))
hc = rbind(hc, length(which(lupusBD2way$A >= 0.85 & lupusBD2way$A < 0.90)))

hc = rbind(hc, length(which(lupusBD2way$dcor >= 0.99 )))
hc = rbind(hc, length(which(lupusBD2way$dcor >= 0.95 & lupusBD2way$dcor < 0.99)))
hc = rbind(hc, length(which(lupusBD2way$dcor >= 0.90 & lupusBD2way$dcor < 0.95)))
hc = rbind(hc, length(which(lupusBD2way$dcor >= 0.85 & lupusBD2way$dcor < 0.90)))

hc = rbind(hc, length(which(lupusBD2way$MIC >= 0.99 )))
hc = rbind(hc, length(which(lupusBD2way$MIC >= 0.95 & lupusBD2way$MIC < 0.99)))
hc = rbind(hc, length(which(lupusBD2way$MIC >= 0.90 & lupusBD2way$MIC < 0.95)))
hc = rbind(hc, length(which(lupusBD2way$MIC >= 0.85 & lupusBD2way$MIC < 0.90)))

hc = cbind(hc, rep(c(rep("R2", 4), rep("A", 4), rep("dcor", 4), rep('MIC', 4)), 2))
hc = cbind(hc, c(rep("Control", 16), rep("Lupus",16)))
hc = cbind(hc, rep(c("0.99", "0.95", "0.90", "0.85"), 4))
colnames(hc) = c('Count', "Association", "Sample", "Level")
require(grid)
ggplot(hc[hc$Level != c("0.85"),],aes(x=Association, y=Count, fill=Level)) + geom_bar(stat = "identity") + 
  scale_fill_brewer(type="seq") + facet_grid(Sample ~ .) + theme_bw() + 
  coord_flip() + theme(legend.position='bottom',legend.key.size = unit(.4, "cm")) 



#set diff and union #########################
r2ResC = controlBD2way[controlBD2way$R2 >= 0.90,3]
aResC = controlBD2way[controlBD2way$A >= 0.90,3]
dcorResC = controlBD2way[controlBD2way$dcor >= 0.90,3]
micResC = controlBD2way[controlBD2way$MIC >= 0.90,3]

r2ResL = lupusBD2way[lupusBD2way$R2 >= 0.90,3]
aResL = lupusBD2way[lupusBD2way$A >= 0.90,3]
dcorResL = lupusBD2way[lupusBD2way$dcor >= 0.90,3]
micResL = lupusBD2way[lupusBD2way$dcor >= 0.90,3]

require(gplots)
venn( list('R2'=r2ResC,'A'=aResC,'dcor'=dcorResC, 'MIC'=micResC) )
venn( list('R2'=r2ResL,'A'=aResL,'dcor'=dcorResL, 'MIC'=micResL) )


lupusBD2way = cbind(lupusBD2way, lupusBD2way$dcor - lupusBD2way$R2)
lupusBD2way = cbind(lupusBD2way, lupusBD2way$dcor - controlBD2way$dcor)
colnames(lupusBD2way)[6:7] = c('betaDcor', 'gammaDcor')

controlBD2way = cbind(controlBD2way, controlBD2way$dcor - controlBD2way$R2)
controlBD2way = cbind(controlBD2way, controlBD2way$dcor - lupusBD2way$dcor)
colnames(controlBD2way)[6:7] = c('betaDcor', 'gammaDcor')

require(plyr)

which(controlBD2way$dcor >= 0.9 & lupusBD2way$dcor >= 0.9)
bothProbes = c(30,31,86,76,42,41)
lllBig = arrange(lupusBD2way[lupusBD2way$dcor >=0.90 & lupusBD2way$betaDcor >= 0.1 & 
                           lupusBD2way$gammaDcor >=0.45,] ,desc(betaDcor))
cccBig = arrange(controlBD2way[controlBD2way$dcor >=0.90 & controlBD2way$betaDcor >= 0.1 & 
                               controlBD2way$gammaDcor >=0.45,] ,desc(betaDcor))

lllBig = lllBig[,1:7]
lllBig = cbind(lllBig, resultsBig[match(lllBig$X, resultsBig$ProbeNo),8])
lllBig = cbind(lllBig, resultsBig[match(lllBig$Y, resultsBig$ProbeNo),8])

colnames(lllBig)[8:9] = c('XProbe','YProbe')
resultsBig[, 8]

dcor(subLupusBD[,49], subLupusBD[,50])
dcor(lupusBD[,8682], lupusBD[,19281])
dcor(lupusBD[,8682], lupusBD[,19282])
dcor(lupusBD[,8683], lupusBD[,19281])
dcor(lupusBD[,8684], lupusBD[,19282])

which(colnames(lupusBD) == '1552316_a_at')
which(colnames(lupusBD) == '1552485_at')

##outiers in lupus
lllBig$X[1]
lllBig$Y[1]
length(subLupusBD[,lllBig$Y[1]])
min(subLupusBD[,lllBig$Y[1]])
which(subLupusBD[,lllBig$Y[14]] < -2)
#MA43_1311 14, MA43_2222 61

cor(subLupusBD[,lllBig$X[1]], subLupusBD[,lllBig$Y[1]])^2
cor(subLupusBD[-c(14,61),lllBig$X[1]], subLupusBD[-c(14,61),lllBig$Y[1]])^2
plot(subLupusBD[,lllBig$X[1]], subLupusBD[,lllBig$Y[1]], xlab = lllBig$X[1], ylab = lllBig$Y[1], main="")

lupusBD2way = cbind(lupusBD2way, rep(0, nrow(lupusBD2way)))
colnames(lupusBD2way)[8] = 'R2NoOutliers'
for(i in 1:nrow(lupusBD2way))
{
  lupusBD2way[i, 8] = cor(subLupusBD[-c(14,61), lupusBD2way$X[i]], subLupusBD[-c(14,61), lupusBD2way$Y[i]])^2
}
lupusBD2way = cbind(lupusBD2way, lupusBD2way$dcor - lupusBD2way$R2NoOutliers)
colnames(lupusBD2way)[9] = 'betaDcorNoOutliers'


cor(subLupusBD[-270, lupusBD2way$X[1]], subLupusBD[-270, lupusBD2way$Y[1]])

#control outliers
which(subControlBD[,cccBig$Y[1]] < -1.9)
#MA43_M.107 MA43_M.184 
#3         19 
cor(subControlBD[,cccBig$X[1]], subControlBD[,cccBig$Y[1]])^2
cor(subControlBD[-3,cccBig$X[1]], subControlBD[-3,cccBig$Y[1]])^2
cor(subControlBD[-19,cccBig$X[1]], subControlBD[-19,cccBig$Y[1]])^2
cor(subControlBD[-c(3,19),cccBig$X[1]], subControlBD[-c(3,19),cccBig$Y[1]])^2

controlBD2way = cbind(controlBD2way, rep(0, nrow(controlBD2way)))
for(i in 1:nrow(controlBD2way))
{
  #controlBD2way[i, 10] = cor(subControlBD[-c(3,19), controlBD2way$X[i]], subControlBD[-c(3,19), controlBD2way$Y[i]])
  controlBD2way[i, 10] = controlBD2way[i,10]^2
}
colnames(controlBD2way)[10] = 'R2NoOutliers'
controlBD2way = cbind(controlBD2way, controlBD2way$dcor - controlBD2way$R2NoOutliers)
colnames(controlBD2way)[11] = 'betaDcorNoOutliers'
require(rafalib)
mypar(mfrow=c(4,4))
for(i in 1:16)
{
  plot(subBDAll[,lllBig$X[i]], subBDAll[,lllBig$Y[i]], xlab = lllBig$X[i], ylab = lllBig$Y[i], main="")
}

for(i in 17:32)
{
  plot(subLupusBD[,lllBig$X[i]], subLupusBD[,lllBig$Y[i]], xlab = lllBig$X[i], ylab = lllBig$Y[i], main="")
}

for(i in 33:48)
{
  plot(subLupusBD[,lllBig$X[i]], subLupusBD[,lllBig$Y[i]], xlab = lllBig$X[i], ylab = lllBig$Y[i], main="")
}

for(i in 49:64)
{
  plot(subLupusBD[,lllBig$X[i]], subLupusBD[,lllBig$Y[i]], xlab = lllBig$X[i], ylab = lllBig$Y[i], main="")
}
#####
for(i in 1:16)
{
  plot(subControlBD[,cccBig$X[i]], subControlBD[,cccBig$Y[i]], xlab = cccBig$X[i], ylab = cccBig$Y[i], main="")
}

for(i in 17:32)
{
  plot(subControlBD[,cccBig$X[i]], subControlBD[,cccBig$Y[i]], xlab = cccBig$X[i], ylab = cccBig$Y[i], main="")
}

for(i in 33:48)
{
  plot(subControlBD[,cccBig$X[i]], subControlBD[,cccBig$Y[i]], xlab = cccBig$X[i], ylab = cccBig$Y[i], main="")
}

for(i in 49:64)
{
  plot(subControlBD[,cccBig$X[i]], subControlBD[,cccBig$Y[i]], xlab = cccBig$X[i], ylab = cccBig$Y[i], main="")
}

dim(cccBig)
count(union(lllBig$X, lllBig$Y) < 119)
count(union(cccBig$X, cccBig$Y) < 119)

which(cccBig$X < 119 & cccBig$Y < 119)
cccBig$X[128]
cccBigDE = cccBig[c(438, 511),]
plot(subLupusBD[, 52], subLupusBD[, 99])
plot(subLupusBD[, 90], subLupusBD[, 99])

summary(subControlBD[,15226])
hist(subControlBD[,15226])
hist(subLupusBD[,15226])
plot(controlBD[,49], controlBD[,50])


su3 = union(lll3$X, lll3$Z)
yy = unique(lll3$Y)
intersect(yy, su3)

require(rafalib)
mypar(mfrow=c(2,1))
plot(table(c(lllBig$X,lllBig$Y)), xlab = 'Probe No', ylab='Count',main='Lupus')
plot(table(c(cccBig$X,cccBig$Y)), xlab = 'Probe No', ylab='Count',main='Controls')


#do associaiton s match in test and train
train2way = data.frame(rep(0, 15), rep(0, 15), rep(0, 15))
colnames(train2way) = c('X','Y','Pair')
for(i in 1:15)
{
  train2way[i,1] = colnames(genes[,-c(1,96)])[lll$X[i]]
  train2way[i,2] = colnames(genes[,-c(1,96)])[lll$Y[i]]
  train2way[i,3] = paste(train2way[i,1], train2way[i,2] , sep=':')
}

test2way = data.frame(rep(0, 57), rep(0, 57), rep(0, 57))
colnames(test2way) = c('X','Y','Pair')
for(i in 1:57)
{
  test2way[i,1] = colnames(safeBDAll)[lllBig$X[i]]
  test2way[i,2] = colnames(safeBDAll)[lllBig$Y[i]]
  test2way[i,3] = paste(test2way[i,1], test2way[i,2] , sep=':')
}

intersect(train2way$X, test2way$Y)
require(energy)
for(i in 1:15)
{
  print(cor(safeLupusBD[,train2way$X[i]], safeLupusBD[,train2way$Y[i]]))
}

for(i in 1:15)
{
  print(cor(safeControlBD[,train2way$X[i]], safeControlBD[,train2way$Y[i]]))
}



lupusBD2way = cbind(lupusBD2way, data.frame(rep(0, nrow(lupusBD2way)), rep(0, nrow(lupusBD2way))))
colnames(lupusBD2way)[8:9] = c('XProbe', 'YProbe')

lupusBD2way$XProbe = colnames(safeBDAll)[lupusBD2way$X]
lupusBD2way$YProbe = colnames(safeBDAll)[lupusBD2way$Y]

any(lupusBD2way$XProbe == train2way$X && lupusBD2way$YProbe == train2way$Y)

hist(subControlBD[,cccBig$Y[1]])
hist(subControlBD[,cccBig$Y[2]])
hist(subControlBD[,cccBig$Y[3]])
hist(subControlBD[,cccBig$Y[4]])
###which MIC in lupus?
lllBigMIC = arrange(lupusBD2way[lupusBD2way$MIC >=0.90 & lupusBD2way$dcor <= 0.9,] ,desc(betaDcor))

###which A in controls?
cccBigA = arrange(controlBD2way[controlBD2way$A >=0.90 & controlBD2way$dcor <= 0.9,] ,desc(betaDcor))
cccBigMIC = arrange(controlBD2way[controlBD2way$MIC >=0.90 & controlBD2way$dcor <= 0.9,] ,desc(betaDcor))
plot(subControlBD[,49994], subControlBD[,49995])


#now without outliers
lllBigNoOut = arrange(lupusBD2way[lupusBD2way$dcor >=0.90 & lupusBD2way$betaDcorNoOutliers >= 0.1 & 
                               lupusBD2way$gammaDcor >=0.45,] ,desc(betaDcorNoOutliers))
cccBigNoOut = arrange(controlBD2way[controlBD2way$dcor >=0.90 & controlBD2way$betaDcorNoOutliers >= 0.1 & 
                                 controlBD2way$gammaDcor >=0.45,] ,desc(betaDcorNoOutliers))

mypar(mfrow=c(4,4))
for(i in 1:16)
{
  plot(subLupusBD[,lllBigNoOut$X[i]], subLupusBD[,lllBigNoOut$Y[i]], xlab = lllBigNoOut$X[i], ylab = lllBigNoOut$Y[i], main="")
}
for(i in 17:32)
{
  plot(subLupusBD[,lllBigNoOut$X[i]], subLupusBD[,lllBigNoOut$Y[i]], xlab = lllBigNoOut$X[i], ylab = lllBigNoOut$Y[i], main="")
}
for(i in 33:48)
{
  plot(subLupusBD[,lllBigNoOut$X[i]], subLupusBD[,lllBigNoOut$Y[i]], xlab = lllBigNoOut$X[i], ylab = lllBigNoOut$Y[i], main="")
}
for(i in 49:64)
{
  plot(subLupusBD[,lllBigNoOut$X[i]], subLupusBD[,lllBigNoOut$Y[i]], xlab = lllBigNoOut$X[i], ylab = lllBigNoOut$Y[i], main="")
}

for(i in 1:16)
{
  plot(subControlBD[,cccBig$X[i]], subControlBD[,cccBig$Y[i]], xlab = cccBigNo$X[i], ylab = cccBigNoOut$Y[i], main="")
}

for(i in 17:32)
{
  plot(subControlBD[,cccBigNoOut$X[i]], subControlBD[,cccBigNoOut$Y[i]], xlab = cccBigNoOut$X[i], ylab = cccBigNoOut$Y[i], main="")
}
for(i in 33:48)
{
  plot(subControlBD[,cccBigNoOut$X[i]], subControlBD[,cccBigNoOut$Y[i]], xlab = cccBigNoOut$X[i], ylab = cccBigNoOut$Y[i], main="")
}
which(subControlBD[,cccBigNoOut$Y[1]] < -4)
dcor(subControlBD[,cccBigNoOut$Y[1]], subControlBD[,cccBigNoOut$X[1]])

require(gridExtra)
plot(subLupusBD[,lllBig$X[1]], subLupusBD[,lllBig$Y[1]], xlab = lllBig$X[1], ylab = lllBig$Y[1], main="")
plot(subBDAll[,lllBig$X[1]], subBDAll[,lllBig$Y[1]], xlab = lllBig$X[1], ylab = lllBig$Y[1], main="")
plot(forg[,lllBig$X[1]], forg[,lllBig$Y[1]], xlab = forg$X[1], ylab = forg$Y[1], main="")

forg = data.frame(subBDAll)
forg = cbind(forg, c(rep('Control', 29), rep('Lupus', nrow(safeLupusBD))))
colnames(forg)[ncol(forg)] = 'diseased'

lllBig[1,1:2]

pp1 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[1,1]]), y=ggname(colnames(forg)[lllBig[1,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth() + theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + xlab(resultsBig[lllBig[1,1],7]) +ylab(resultsBig[lllBig[1,2],7])+ scale_colour_brewer(palette="Set1")
pp2 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[2,1]]), y=ggname(colnames(forg)[lllBig[2,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth() + theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[2,1],7]) +ylab(resultsBig[lllBig[2,2],7])+ scale_colour_brewer(palette="Set1")
pp3 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[3,1]]), y=ggname(colnames(forg)[lllBig[3,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[3,1],7]) +ylab(resultsBig[lllBig[3,2],7])+ scale_colour_brewer(palette="Set1")
pp4 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[4,1]]), y=ggname(colnames(forg)[lllBig[4,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[4,1],7]) +ylab(resultsBig[lllBig[4,2],7])+ scale_colour_brewer(palette="Set1")
pp5 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[5,1]]), y=ggname(colnames(forg)[lllBig[5,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[5,1],7]) +ylab(resultsBig[lllBig[5,2],7])+ scale_colour_brewer(palette="Set1")
pp6 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[6,1]]), y=ggname(colnames(forg)[lllBig[6,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[6,1],7]) +ylab(resultsBig[lllBig[6,2],7])+ scale_colour_brewer(palette="Set1")
pp7 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[7,1]]), y=ggname(colnames(forg)[lllBig[7,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[7,1],7]) +ylab(resultsBig[lllBig[7,2],7])+ scale_colour_brewer(palette="Set1")
pp8 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[8,1]]), y=ggname(colnames(forg)[lllBig[8,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[8,1],7]) +ylab(resultsBig[lllBig[8,2],7])+ scale_colour_brewer(palette="Set1")
pp9 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[9,1]]), y=ggname(colnames(forg)[lllBig[9,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[9,1],7]) +ylab(resultsBig[lllBig[9,2],7])+ scale_colour_brewer(palette="Set1")
pp10 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[10,1]]), y=ggname(colnames(forg)[lllBig[10,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[10,1],7]) +ylab(resultsBig[lllBig[10,2],7])+ scale_colour_brewer(palette="Set1")
pp11 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[11,1]]), y=ggname(colnames(forg)[lllBig[11,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[11,1],7]) +ylab(resultsBig[lllBig[11,2],7])+ scale_colour_brewer(palette="Set1")
pp12 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[12,1]]), y=ggname(colnames(forg)[lllBig[12,2]]), colour='diseased')) + geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[12,1],7]) +ylab(resultsBig[lllBig[12,2],7])+ scale_colour_brewer(palette="Set1")

pp13 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[13,1]]), y=ggname(colnames(forg)[lllBig[13,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[13,1],7]) +ylab(resultsBig[lllBig[13,2],7])+ scale_colour_brewer(palette="Set1")
pp14 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[14,1]]), y=ggname(colnames(forg)[lllBig[14,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[14,1],7]) +ylab(resultsBig[lllBig[14,2],7])+ scale_colour_brewer(palette="Set1")
pp15 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[15,1]]), y=ggname(colnames(forg)[lllBig[15,2]]), colour='diseased')) + geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[15,1],7]) +ylab(resultsBig[lllBig[15,2],7])+ scale_colour_brewer(palette="Set1")
pp16 = ggplot(forg, aes_string(x=ggname(colnames(forg)[lllBig[16,1]]), y=ggname(colnames(forg)[lllBig[16,2]]), colour='diseased')) +  geom_point(size=1,alpha=.7) + geom_smooth()+ theme_bw()  + theme(legend.position="none", panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ xlab(resultsBig[lllBig[16,1],7]) +ylab(resultsBig[lllBig[16,2],7])+ scale_colour_brewer(palette="Set1")

grid.arrange(arrangeGrob(pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10,pp11,pp12,pp13,pp14,pp15,pp16, ncol=4, nrow=4))

