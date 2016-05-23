install.packages('e1071')
require(randomForest)
?randomForest
head(genes)
subgene = genes[,-c(1)]
subgene = subgene[,-c(94,95,96,97)]
subgene = 2^subgene
subgene = cbind(subgene, genes$diseased)
colnames(subgene)[94] = 'diseased'
dim(subgene)
fit <- randomForest(x = genes[,-c(1,96)], y = genes$diseased,  importance=TRUE, ntree=2000)
varImpPlot(fit)
fit$confusion






fit <- randomForest(x = genes[,-c(1,96,97)], y = genes$Type,  importance=TRUE, ntree=2000)
varImpPlot(fit)
fit$confusion

#unsupervised, fails to discriminate
iris.urf <- randomForest(genes[,-c(1,96,97)])
MDSplot(iris.urf, genes$Type)
varImpPlot(iris.urf)




# load the package
library(e1071)
?naiveBayes
# fit model
fit <- naiveBayes(x = genes[,-c(1,96,97)], y = genes$Type)
# summarize the fit
summary(fit)
# make predictions
predictions <- predict(fit, genes[,-c(1,96,97)])
# summarize accuracy
table(predictions, genes$Type)



install.packages('kernlab')
# load the package
library(kernlab)
?ksvm
# fit model
fit <- ksvm(x = genes[,-c(1,96,97)], y = genes$Diseased)
# summarize the fit
summary(fit)
# make predictions
predictions <- predict(fit, genes[,-c(1,96,97)])
# summarize accuracy
table(predictions, iris$Species)

## create test and training set
ind <- sample(1:nrow(genes),60)
genetrain <- genes[-ind, -c(1,96,97)]
genetest <- genes[ind, -c(1,96,97)]

## train a support vector machine
length(genes$diseased[-ind])
nrow(genetrain)
gene <-  ksvm(y=genes$diseased[-ind],x=genetrain,kernel="rbfdot",kpar=list(sigma=0.015),C=70,cross=4,prob.model=TRUE)

## predict gene type probabilities on the test set
genetype <- predict(gene,genetest,type="probabilities")

(myrf <- randomForest(genetrain, genes$diseased[-ind]))
plot(predict(myrf, genetest),genes$diseased[ind])



