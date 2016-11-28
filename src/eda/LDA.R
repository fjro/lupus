library(MASS)
library(ggplot2)
library(gdata)
library(gridExtra)
library(reshape2)
library(grid)

#discriminate by disease status
lda.1 <- lda(data=genes[,-c(1,97)],diseased~.)
summary(lda.1)
lda.cc <- predict(lda.1)$x
summary(predict(lda.1))

#build up the frame for plotting
lddf <- as.data.frame(lda.cc)
lddf$Sample <- 1:length(lda.cc)
lddf$Patient <- genes$diseased
lddf$Source <- genes$Source
colnames(lddf)[2:4] = c('Sample',"Patient","Source")

ggplot(lddf, aes(y=LD1,x=Sample,shape=Patient, colour=Source)) +
  geom_point(size=3) + 
  theme(legend.position='bottom') + 
  xlab('Sample ID')
#reasonable discrimination based on Fisher's ldf

#now look at the scalings to see which genes are contributing most
scalings <- as.data.frame(lda.1$scaling)
scalings$Probe <- rownames(scalings)
scalings$ProbeNo <- 1:94

ggplot(scalings, aes(y=LD1,x=ProbeNo)) + 
  scale_x_discrete(limits = 1:94) + 
  xlab("Probe No.") +
  geom_bar(stat = "identity",position="dodge") + 
  theme(axis.text.x=element_text(angle=270, vjust = 0.5, size=6))

#now refine by type
lda.2 <- lda(data=genes[,-c(1,96)],Source~.)
lda.cc2 <- predict(lda.2)$x

#build up the frame for plotting
lddf2 <- as.data.frame(lda.cc2)
lddf2$Patient <- genes$diseased
lddf2$Source <- genes$Source

p1 <- ggplot(lddf2, aes(y=LD1,x=LD2,shape=Patient, colour = Source)) +geom_point(size=2) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
p2 <- ggplot(lddf2, aes(y=LD3,x=LD4,shape=Patient, colour = Source)) +geom_point(size=2)+ theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))


grid.arrange(arrangeGrob(p1,p2, ncol=2, nrow=1))
#no discrimination of interest but look at scalings anyway

scalings2 <- as.data.frame(lda.2$scaling)
scalings2$Probe <- rownames(scalings2)
scalings2$ProbeNo <- as.factor(1:94)

mm <- melt(scalings2)
colnames(mm)[3] <- 'Scaling'
p <-  ggplot(mm, aes(y=value, x=ProbeNo))+xlab("Probe No") + 
  scale_x_discrete(limits = 1:94) + 
  geom_bar(stat = "identity", position="dodge") + 
  theme(axis.text.x=element_text(angle=270, vjust = 0.5,size=6))
p + facet_grid(Scaling~.)  

