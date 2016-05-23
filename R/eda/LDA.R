require(MASS)
require(ggplot2)
require(gdata)
require(gridExtra)
require(reshape)

#discriminate by disease status
?lda
lda.1 = lda(data=genes[,-c(1,97)],diseased~.)
summary(lda.1)
lda.cc <- predict(lda.1)$x
summary(predict(lda.1))

#build up the frame for plotting
lddf = as.data.frame(lda.cc)
lddf = cbind(lddf, 1:length(lda.cc))
lddf = cbind(lddf, genes$diseased)
lddf = cbind(lddf, genes$Source)
colnames(lddf)[2:4] = c('Sample',"Patient","Source")

ggplot(lddf, aes(y=LD1,x=Sample,shape=Patient, colour=Source)) +geom_point(size=3) + theme(legend.position='bottom') + xlab('Sample ID')

scalings = as.data.frame(lda.1$scaling)
scalings = cbind(rownames(scalings),scalings)
colnames(scalings)[1] = 'Probes'
scalings = cbind(1:94, scalings)
colnames(scalings)[1] = 'ProbeNo'

ggplot(scalings, aes(y=LD1,x=ProbeNo))+ scale_x_discrete(limits = 1:94)+ xlab("Probe No.")+
  geom_bar(stat = "identity",position="dodge")+ theme(axis.text.x=element_text(angle=270, vjust = 0.5, size=6))


#now refine by type
lda.2 = lda(data=genes[,-c(1,96)],Source~.)
lda.cc2 <- predict(lda.2)$x

#build up the frame for plotting
lddf2 = as.data.frame(lda.cc2)
lddf2 = cbind(lddf2, genes$diseased)
lddf2 = cbind(lddf2, genes$Source)
colnames(lddf2)[6:7] = c("Patient","Source")
require(grid)
p1 = ggplot(lddf2, aes(y=LD1,x=LD2,shape=Patient, colour = Source)) +geom_point(size=2) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
p2 = ggplot(lddf2, aes(y=LD3,x=LD4,shape=Patient, colour = Source)) +geom_point(size=2)+ theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))

require(gridExtra)
grid.arrange(arrangeGrob(p1,p2, ncol=2, nrow=1))
grid.arrange(p1,p2, ncol=2, nrow=1, main="", plot=T)
#no discrimination of interest
#ggplot(lddf2, aes(y=LD4,x=LD5,shape=Disease, colour = Type)) +geom_point(size=4)

?melt
scalings2 = as.data.frame(lda.2$scaling)
scalings2 = cbind(rownames(scalings2),scalings2)
colnames(scalings2)[1] = 'Probe'
scalings2 = cbind(as.factor(1:94), scalings2)
colnames(scalings2)[1] = 'ProbeNo'
mm = melt(scalings2, id=)
colnames(mm)[3] = 'Scaling'
p= ggplot(mm, aes(y=value,x=ProbeNo))+xlab("Probe No")+ scale_x_discrete(limits = 1:94)+ geom_bar(stat = "identity",position="dodge")+ theme(axis.text.x=element_text(angle=270, vjust = 0.5,size=6))
p+facet_grid(Scaling~.)  

