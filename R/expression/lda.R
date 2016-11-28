library(caret)
Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
                   Sp = rep(c("s","c","v"), rep(50,3)))
train <- sample(1:150, 75)
table(Iris$Sp[train])

z <- MASS::lda(Sp ~ ., Iris, prior = c(1,1,1)/3, subset = train)
predict(z, Iris[-train, ])$class

 




levels(j2m2$Status)
#priors <- c(length(which(normed$Status == 'Bad'))/nrow(normed),length(which(normed$Status == 'Good'))/nrow(normed))
train <- 1:5000
table(normed$Status[train])
lda_1 <- MASS::lda(normed[,1000:2500], grouping = normed$Status, subset = train)
confusionMatrix(predict(lda_1, normed[-train,1000:2500])$class, normed[-train,]$Status)
length(normed[-train,]$Status)

lda.cc <- (predict(lda_1, normed[-train,1000:2500])$x)


#build up the frame for plotting
lddf <- as.data.frame(lda.cc)
lddf$Sample <- 1:length(lda.cc)
lddf$Timestamp <- normed$Timestamp[-train]
lddf$Status <- normed$Status[-train]
lddf$Plasma <- normed$Plasma[-train]

ggplot(lddf, aes(y=LD1,x=Timestamp,shape=Plasma, colour=Status, alpha = 0.6)) +
  geom_point(size=3) + 
  theme(legend.position='bottom') + 
  xlab('Timestamp')


res <- apply(normed[,-c(1:2)], 2, function(x) any(is.na(x)))
which(res)
lda_1$

colnames(normed)[7273]



which(is.na(normed[,7252]))
colnames(normed[,7250:7273])
confusionMatrix()

complete.cases(normed[,100])
?complete.cases
