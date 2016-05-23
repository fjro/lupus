#install.packages('mda')
library(mda)
require(xtable)

#first classify by disease
fit = fda(data=genes[,-c(1,97)],diseased~.,method=mars)
fit$confusion
#only a single dimension so no plot
predictions <- predict(fit, genes[,-c(1,97)])
# summarize accuracy
confusion(predictions, genes$diseased)

fit = fda(data=genes[,-c(1,97)],diseased~.,method=mars)
fit$confusion
#only a single dimension so no plot
predictions <- predict(fit, genes[,-c(1,97)])

#latex
x.small <- xtable(fit$confusion, label = 'tabsmall', caption = 'FDA confusion matrix for disease classification')
print(x.small,latex.environments = "",table.placement = 'h')



#now classify types
fit = fda(data=genes[,-c(1,96)],Type~.,method=mars)
fit$confusion
plot(fit)

#chaging degree made no difference
#marsfit2 <- update(fit, degree = 2)

x.small <- xtable(fit$confusion, label = 'tabsmall', caption = 'FDA confusion matrix for type classification')
print(x.small,latex.environments = "",table.placement = 'h')

coef(fit)

library(mda)
(fda <- fda(Species~., data=iris, keep.fitted=TRUE, method=earth, keepxy=TRUE))
summary(fda$fit) # examine earth model embedded in fda model
plot(fda) # right side of the figure

#Using plotmo we can plot the per-predictor dependence of the fda variates like this:
plotmo(fda, type="variates", nresponse=1, clip=F) # 1st disc var (Figure 5)
plotmo(fda, type="variates", nresponse=2, clip=F) # 2nd disc var (not shown)

#We can also look at the earth model embedded in the FDA model:
plotmo(fda$fit, nresponse=1, clip=F) # earth in FDA, 1st disc var (Figure 6)
plotmo(fda$fit, nresponse=2, clip=F) # earth in FDA, 2nd disc var (not shown)
#The graphs show the contribution of each predictor to the first discriminant variable.

##################
require(earth)
?earth
earthsimple = genes[,-c(1,97)]
colnames(simple)[1:94] = 1:94 
(fda <- fda(diseased~.,simple, keep.fitted=TRUE, method=earth, keepxy=TRUE, degree=1, linpreds=F))
fda$confusion
(summary(fda$fit)) # examine earth model embedded in fda model

#now allow hinge functions and interactions
(fda <- fda(diseased~.,simple, keep.fitted=TRUE, method=earth, keepxy=TRUE, degree=2, linpreds=F))
fda$confusion
(summary(fda$fit)) 

colnames(genes[,-c(1,97)])[c(33,75,80,94,14,24,25,68,80,84,47,52,61,68)]

plot(fda)
plot(fda$fit) # right side of the figure
plotd(fda$fit, hist = T)
fda$confusion
#Using plotmo we can plot the per-predictor dependence of the fda variates like this:
plotmo(fda, type="variates", clip=F) # 1st disc var (Figure 5)
plot(simple$'14', simple$'94')

#now try type



#We can also look at the earth model embedded in the FDA model:
plotmo(fda$fit, nresponse=1, clip=F,trace=1) # earth in FDA, 1st disc var (Figure 6)
#The graphs show the contribution of each predictor to the first discriminant variable.
plotmo(fda$fit, clip=F,caption="")
?plotmo






##try glm
# equivalent but using earth.default
a1a <- earth(x=simple[,-95], y=simple[,95],trace=1, glm=list(family=binomial), degree=2)
plotd(a1a, hist=T)
summary(a1a)
plotmo(a1a)
?earth
#now try fda with glm
fda <- fda(simple$diseased~.,simple, keep.fitted=TRUE, method=earth, keepxy=TRUE, degree=2, glm=list(family=binomial))

probe_nos = data.frame(1:94)
dim(probe_nos)
probe_nos = cbind(probe_nos, colnames(genes[-c(1,96,97)]))
colnames(probe_nos) = c('NO', 'ID')
colnames(genes)
sort(probe_nos[probe_nos$ID %in% c('210916_s_at', '210916_s_at',
                             '12418_at',
                             '1565868_at',
                             '200695_at',
                             '202086_at',
                             '202145_at',
                             '202411_at',
                             '202869_at',
                             '202886_s_at',
                             '203153_at',
                             '204415_at',
                             '204439_at',
                             '204489_s_at',
                             '204490_s_at',
                             '204613_at',
                             '204747_at',
                             '204972_at',
                             '205483_s_at',
                             '205569_at',
                             '205828_at',
                             '207540_s_at',
                             '207630_s_at',
                             '207892_at',
                             '209835_x_at',
                             '210171_s_at',
                             '212418_at',
                             '213797_at',
                             '214032_at',
                             '214059_at',
                             '216551_x_at',
                             '218400_at',
                             '219211_at',
                             '219684_at',
                             '219863_at',
                             '227609_at',
                             '228092_at',
                             '230511_at',
                             '241740_at',
                             '241812_at',
                             '241916_at',
                             '44673_at'),])
