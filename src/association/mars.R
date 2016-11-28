library(earth)
library(pROC)

y <- genes$diseased
x <- genes[,-c(1,34,96)]
fit1 <- earth(y=y, x=x,glm=list(family=binomial), degree=2)
summary(fit1)
plot(fit1$glm.list[[1]], info=T) 
plotd(fit1, hist = T)
dd = fit1$bx
plot(dd[,4],dd[,8])
par(mfrow=c(2,2))
plot(fit1$glm.list[[1]])


#cross validated
fit= earth(y=y, x=x,glm=list(family=binomial), degree=2, nfold=10, keepxy=TRUE)
plot(fit, which=1, col.rsq=0) # which=1 for Model Selection plot only (optional)
plot(fit$glm.list[[1]])
plot.earth.models(fit$cv.list, which=1, ylim=c(0, .5))

#Using cross-validation to select the number of terms
fit <- earth(y=ytrain, x=xtrain, degree=2, ncross=10, nfold=2, stratify = T)
#plot(fit, which=1, col.rsq=0)
fs = summary(fit)

fit$glm.list
fs$
#variable importance
evimp(fit, trim=T)
predsT=predict(fit, newdata=xtest, type ='class')
preds=predict(fit, newdata=xtest)
roc1=roc(ytest ~ preds)
roc1
?predict.earth
sum(ytest != preds)

names(evimp(fit, trim=T))


preds=predict(fit1)
roc1=roc(y ~ preds)

preds2=predict(fit)
roc2=roc(y ~ preds2)

plot(roc1)
lines(roc2, col='red')

y = genes$diseased == 'Control'
x = genes[,-c(1,96)]
preds=predict(earth(y = y, x = genes[,c(31,32)], degree=2))
roc1=roc(y ~ preds)
roc1$auc
plot(roc1)

preds=predict(earth(y = y, x = genes[,c(85,77)], degree=2, linpreds=F))
roc1=roc(y ~ preds)
roc1$auc
lines(roc1)

preds=predict(earth(y = y, x = genes[,c(43,42)], degree=2))
roc1=roc(y ~ preds)
roc1$auc
lines(roc1)

preds2=predict(fit2)
roc2=roc(a ~ preds2)

plot(roc2, add=TRUE, col='red')
plot(earth(genes[,50],as.data.frame(genes[,60])))
fit.sp = smooth.spline(genes[,50] ~ genes[,50], nknots=15)
lines(fit1.sp)

which(colnames(genes)=="212418_at")

##########################################################################################
par(mfrow=c(2,2))
y = genes$diseased == 'Control'
xtrain = genes[testg,-c(1,96)]
ytrain = y[testg]

xtest = genes[-testg,-c(1,96)]
ytest = y[-testg]


preds=predict(earth(y = ytrain, x = xtrain, degree=1, ncross=10, nfold=2), newdata = xtest)
roc1=roc(ytest ~ preds)
roc1$auc
plot(roc1)

preds=predict(earth(y = ytrain, x = xtrain, degree=2, ncross=10, nfold=2), newdata = xtest)
roc2=roc(ytest ~ preds)
roc2$auc
lines(roc2, col='red')

preds=predict(earth(y = ytrain, x = xtrain, degree=3, ncross=10, nfold=2), newdata = xtest)
roc3=roc(ytest ~ preds)
roc3$auc
lines(roc3, col='blue')

preds=predict(earth(y = ytrain, x = xtrain, degree=1, linpreds=T, ncross=10, nfold=2), newdata = xtest)
roc4=roc(ytest ~ preds)
roc4$auc
plot(roc4)

preds=predict(earth(y = ytrain, x = xtrain, degree=2, linpreds=T, ncross=10, nfold=2), newdata = xtest)
roc5=roc(ytest ~ preds)
roc5$auc
lines(roc5, col='red')

preds=predict(earth(y = ytrain, x = xtrain, degree=3, linpreds=T, ncross=10, nfold=2), newdata = xtest)
roc6=roc(ytest ~ preds)
roc6$auc
lines(roc6, col='blue')

  
require(ggplot2)
n = length(roc1$sensitivities)
df = data.frame(roc1$sensitivities, roc1$specificities, rep('1', n), rep("Non-Linear", n))
colnames(df) = c('Sensitivity', 'Specificity', 'Degree', 'Type')

n = length(roc2$sensitivities)
df2 = data.frame(roc2$sensitivities, roc2$specificities, rep('2', n), rep("Non-Linear", n))
colnames(df2) = c('Sensitivity', 'Specificity', 'Degree', 'Type')

n = length(roc3$sensitivities)
df3 = data.frame(roc3$sensitivities, roc3$specificities, rep('3', n), rep("Non-Linear", n))
colnames(df3) = c('Sensitivity', 'Specificity', 'Degree', 'Type')

n = length(roc4$sensitivities)
df4 = data.frame(roc4$sensitivities, roc4$specificities, rep('1', n), rep("Linear", n))
colnames(df4) = c('Sensitivity', 'Specificity', 'Degree', 'Type')

n = length(roc5$sensitivities)
df5 = data.frame(roc5$sensitivities, roc5$specificities, rep('2', n), rep("Linear", n))
colnames(df5) = c('Sensitivity', 'Specificity', 'Degree', 'Type')

n = length(roc6$sensitivities)
df6 = data.frame(roc6$sensitivities, roc6$specificities, rep('3', n), rep("Linear", n))
colnames(df6) = c('Sensitivity', 'Specificity', 'Degree', 'Type')

df = rbind(df, df2)
df = rbind(df, df3)
df = rbind(df, df4)
df = rbind(df, df5)
df = rbind(df, df6)

ggplot(df, aes(x=Sensitivity, y=Specificity,group=Degree, colour=Degree)) + geom_line(size=1.2)+ theme(legend.position='bottom') + scale_x_reverse() + 
  facet_grid(. ~ Type) + scale_colour_brewer(palette="Set1")
  
#select variables
fit = earth(y = y, x = x, degree=2, ncross=10, nfold=2)
evimp(fit)


plot(fit3, which=1, col.rsq=0)










#need to use formula interface
ed = cbind(xtrain, ytrain)
colnames(ed)[95] = 'Y'
fit2 = earth(Y~., data=ed, degree = 2, ncross=10, nfold=2)
res = drop1(fit2)
res$AIC
which(res$AIC > 0.0485)
res$Df
summary(res$AIC)

fit3 = earth(y = y, x = x, degree=3, ncross=10, nfold=2)
res2 = evimp(fit2)
names(res[,1])
values(res[,1])

su = union(results[results$Association =='dcor' & results$Value >= 0.95,1:2],results[results$Association =='dcor' & results$Value >= 0.95,2])
#probes correlated in both controls and lupus
clProbes = c(31,30,76,86,41,42)
s2 = unique(unlist(su))
length(su)

s3 = union(lll$X,lll$Y)
s3 = union(lll5$X, s3)
s3 = union(lll5$Y, s3)
s3 = union(lll5$Z, s3)

#how many genes)
gr = unique(annot[annot$PROBENO %in% su,6])
length(gr)

unique(annot[annot$PROBEID %in% su,6])

#how many chromosomes
cr = unique(annot[annot$PROBENO %in% su,3])
length(cr) 

#how many differentially expressed
de = annot[annot$PROBEID %in% diffExp,1]
require(rafalib)
mypar()
venn( list('MARS'=marsProbes$`Probe No`,'DE'=de, 'dcor'=s3) )
DE
lll5$X
intersect(marsProbes$`Probe No`, s3)
require(gplots)
venn()

annot[annot$PROBEID %in% marsRes,1]

summary(fit)
res2 = evimp(fit)
marsRes = names(res2[,1])
marRes = gsub("`", "", marsRes)

marsProbes = data.frame(rep(0,length(marsRes)), rep(0,length(marsRes)), rep(0,length(marsRes)), rep(0,length(marsRes)))
for (i in 1:length(marsRes))
{
  marsProbes[i,1] = which(colnames(genes)[-c(1,96)] == marRes[i])
  marsProbes[i,2] = marRes[i]
  marsProbes[i,3] = sum(marsProbes[i,1] %in% lll$X) + sum(marsProbes[i,1] %in% lll$Y)
  marsProbes[i,4] = sum(marsProbes[i,1] %in% lll5$X) + sum(marsProbes[i,1] %in% lll5$Y)+ sum(marsProbes[i,1] %in% lll5$Z)
  marsProbes[i,5] = marsProbes[i,1] %in% de
}



colnames(marsProbes) = c('Probe No', 'ID', 'Two-way', 'Three-way', 'Differentially Expressed')
marsProbes = marsProbes[-11,]
marsProbes = marsProbes[-9,]
lll5

x.small <- xtable(marsProbes, label = 'tabsmall', caption = 'MARS', digits = 0)
print(x.small,latex.environments = "",table.placement = 'H')


marsProbes = which(colnames(genes)[-c(1,96)] %in% marRes)
marsProbes = data.frame(marsProbes, colnames(genes)[marsProbes +1])


as.character()marsRes = names(res2[,1])

s
### predict big

?predict.earth
colnames(allBD)[1:10]
dim[genes]
dim(genes)
cols = which(colnames(allBD )%in% intersect(colnames(allBD), colnames(genes[,-c(1,96)])))
nd = allBD[,cols]
dim(nd)


ta = c('219684_at', '205483_s_at')
which(colnames(allBD) %in% ta)

x2 = genes[,c(84,93)]
n2 = allBD[,c(14931, 28969)]

?glm
results

yBD = c(rep(T,232),c(rep(F,367)))
dim(x)
yL = genes$diseased == 'Lupus'
fit= earth(y = yL, x = x, degree=3, linpreds=F)
summary(fit)
fit$grsq
plot(fit)
preds = predict(fit,type='link')
rocBD = roc(y~preds)
plot(rocBD)
evimp(fit,trim=T)
?evimp
table(preds)

require(pROC)


roc5=roc(y ~ predict(fit))
roc5$auc
fit$rss
summary(fit)
fit$selected.terms
evimp(fit)

install.packages('caret')
require(caret)
?bagEarth
y = genes[,34]


#can't use bagging with a binary response
fit2 = earth(y = y, x = x, degree=2)
fit2 = bagEarth(y = y, x = x, degree=2, B=5,)
varImp(fit2)

testg = c(sample(1:53, 37), sample(54:420, 257))
xt = genes[testg,-c(1,96)]
fit = earth(y = yt, x = xt, degree=2, linpreds=F)
summary(fit)
pred = predict(fit, newdata = genes[-testg,])
plot(roc(y[-testg] ~ pred))



pred = predict(earth(y = yt, x = xt, degree=2, linpreds=T), newdata = genes[-testg,])
lines(roc(y[-testg] ~ pred), col= 'red')

pred = predict(earth(y = yt, x = xt, degree=3, linpreds=F), newdata = genes[-testg,])
lines(roc(y[-testg] ~ pred), col= 'green')

pred = predict(earth(y = yt, x = xt, degree=3, linpreds=T), newdata = genes[-testg,])
lines(roc(y[-testg] ~ pred), col= 'blue')

pred = predict(earth(y = yt, x = xt, degree=1, linpreds=F), newdata = genes[-testg,])
lines(roc(y[-testg] ~ pred), col= 'orange')

