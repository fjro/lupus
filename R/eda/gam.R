fromListToMatrix = function(x)
{
  m = matrix(0,94,94)
  index = 1
  for(r in 1:94)
  {
    for (c in 1:94)
    {
      if(c==r)
      {
        m[r,c] = 1
      }
      else if(c >r)
      {
        m[r,c] = x[index]
        index = index + 1
      } 
    }
  }
  m
}

test = fromListToMatrix(corrsL[,1])
rAMatL = fromListToMatrix(corrsL[,6])
rMICMatL = fromListToMatrix(corrsL[,7])
rDcorMatL = fromListToMatrix(corrsL[,8])
rAMatC = fromListToMatrix(corrsC[,6])
rMICMatC = fromListToMatrix(corrsC[,7])
rDcorMatC = fromListToMatrix(corrsC[,8])
diag(rAMatL) = -1
diag(rMICMatL) = -1
diag(rDcorMatL) = -1
diag(rAMatC) = -1
diag(rMICMatC) = -1
diag(rDcorMatC) = -1

cor(corrsL[1:50,6],corrsL[1:50,7], method='kendall')
cor(corrsL[,6],corrsL[,8], method='kendall')
cor(corrsL[1:50,7],corrsL[1:50,8], method='kendall')

?which

colnames(genes)
mainGAM = cbind(corrsL[1:50,8])


plotCurves = function(rMat, rVec)
{
  maxDF = colnames(genes)[which(rMat == max(rVec), arr.ind = TRUE)+1]
  so = sort(rVec,decreasing=T)[1:20]
  for(i in 2:20)
  {
    maxDF = rbind(maxDF,colnames(genes)[which(rMat == so[i], arr.ind = TRUE)+1])
  }
  
  
  colnames(genes)[which(rMat == max(rVec), arr.ind = TRUE)+1]
  plots = rep(p2,20)
  for (i in 1:20)
  {
    pt <- 
      ggplot(genes, aes_string(x=maxDF[i,1], y=maxDF[i,2], colour='diseased')) +
      geom_point(alpha=.5) +
      geom_smooth(alpha=.2, size=1)
    plots[i] = pt
  }
  p1 = ggplot(genes, aes_string(x=maxDF[1,1], y=maxDF[1,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p2 = ggplot(genes, aes_string(x=maxDF[2,1], y=maxDF[2,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p3 = ggplot(genes, aes_string(x=maxDF[3,1], y=maxDF[3,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p4 = ggplot(genes, aes_string(x=maxDF[4,1], y=maxDF[4,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p5 = ggplot(genes, aes_string(x=maxDF[5,1], y=maxDF[5,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p6 = ggplot(genes, aes_string(x=maxDF[6,1], y=maxDF[6,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p7 = ggplot(genes, aes_string(x=maxDF[7,1], y=maxDF[7,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p8 = ggplot(genes, aes_string(x=maxDF[8,1], y=maxDF[8,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p9 = ggplot(genes, aes_string(x=maxDF[9,1], y=maxDF[9,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p10 = ggplot(genes, aes_string(x=maxDF[10,1], y=maxDF[10,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p11 = ggplot(genes, aes_string(x=maxDF[11,1], y=maxDF[11,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p12 = ggplot(genes, aes_string(x=maxDF[12,1], y=maxDF[12,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p13 = ggplot(genes, aes_string(x=maxDF[13,1], y=maxDF[13,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p14 = ggplot(genes, aes_string(x=maxDF[14,1], y=maxDF[14,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p15 = ggplot(genes, aes_string(x=maxDF[15,1], y=maxDF[15,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  p16 = ggplot(genes, aes_string(x=maxDF[16,1], y=maxDF[16,2], colour='diseased')) + geom_point(alpha=.3) + geom_smooth(alpha=.2, size=1) + theme_bw() + theme(legend.position="none")
  
  #require(gridExtra)
  grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16, nrow=4, ncol=4) 
}

plotCurves(rAMatL, corrsL[,6])
plotCurves(rDcorMatC, corrsC[,8])
sort(corrsC[,8],decreasing=T)[1:20]
# what sort of non-linear relationships exist

?geom_smooth
which(rAMatC == max(rAMatC), arr.ind = TRUE)
colnames(genes)[c(4,37)]
p2 <- 
  ggplot(genes, aes(x=X241871_at, y=X209189_at, colour=diseased)) +
  geom_point(alpha=.5) +
  geom_smooth(alpha=.2, size=1) +
  ggtitle("GAM for Max Residual A in Controls")

#MIC, max same as A
which(rMICMatC == max(rMICMatC), arr.ind = TRUE) 
colnames(genes)[c(13,88)]
p3 <- 
  ggplot(genes, aes(x=X212014_x_at, y=X202086_at, colour=diseased)) +
  geom_point(alpha=.5) +
  geom_smooth(alpha=.2, size=1) +
  ggtitle("GAM for Max Residual MIC in Controls")
#looks like outliers

which(rDcorMatL == max(rDcorMatL), arr.ind = TRUE)
colnames(genes)[c(75,93)]
p4 <- 
  ggplot(genes, aes(x=X202411_at, y=X219684_at, colour=diseased)) +
  geom_point(alpha=.5) +
  geom_smooth(alpha=.2, size=1) +
  ggtitle("GAM for Max Residual dcor in Lupus")

which(rDcorMatC == max(rDcorMatC), arr.ind = TRUE) #same as A

pmax(1:5, 2)
?head
?sort
order(rAMatL)
which(rAMatL == max(rAMatL), arr.ind = TRUE) 
ctp = 50
l1 = which(rAMatL == max(corrsL[,6]), arr.ind = TRUE)+1 
for(i in 2:ctp)
{
  l1 = rbind(l1,which(rAMatL == sort(corrsL[,6],decreasing=T)[1:ctp][i], arr.ind = TRUE)+1)
}
l2 = which(rMICMatL == max(corrsL[,7]), arr.ind = TRUE)+1 
for(i in 2:ctp)
{
  l2 = rbind(l2,which(rMICMatL == sort(corrsL[,7],decreasing=T)[1:ctp][i], arr.ind = TRUE)+1)
}
l3 = which(rDcorMatL == max(corrsL[,8]), arr.ind = TRUE)+1 
for(i in 2:ctp)
{
  l3 = rbind(l3,which(rDcorMatL == sort(corrsL[,8],decreasing=T)[1:ctp][i], arr.ind = TRUE)+1)
}
#l1 = which(matrix(rAMatL %in% head(sort(rAMatL, T), ctp),nr = nrow(rAMatL)), arr.ind = TRUE)
#l2 = which(matrix(rMICMatL %in% head(sort(rMICMatL, T), ctp),nr = nrow(rMICMatL)), arr.ind = TRUE)
#l3 = which(matrix(rDcorMatL %in% head(sort(rDcorMatL, T), ctp),nr = nrow(rDcorMatL)), arr.ind = TRUE)

#paste(colnames(genes)[l1[,1]], colnames(genes)[l1[,2]], sep=" v ")

#make a ranked list of pairwise associations and compare pearson to A
a <- data.frame(genes = paste(colnames(genes)[l1[,1]], colnames(genes)[l1[,2]], sep=" v "), rank = 1:ctp)
m <- data.frame(genes = paste(colnames(genes)[l2[,1]], colnames(genes)[l2[,2]], sep=" v "), rank = 1:ctp)
d <- data.frame(genes = paste(colnames(genes)[l3[,1]], colnames(genes)[l3[,2]], sep=" v "), rank = 1:ctp)

a <- cbind(a, Association ='A')
a = cbind(a,colnames(genes)[l1[,1]])
colnames(a)[4] = 'gene1'
m <- cbind(m, Association = 'MIC')
m = cbind(m,colnames(genes)[l2[,1]])
colnames(m)[4] = 'gene1'
d <- cbind(d, Association = 'dcor')
d = cbind(d,colnames(genes)[l3[,1]])
colnames(d)[4] = 'gene1'

am.melt <- rbind(a, d)
am.melt = rbind(am.melt, m)
am.melt <- transform(am.melt, offset = c(rep(1.1, times = ctp), rep(15, times = ctp),rep(-.1, times = ctp)))
ggplot(data = am.melt, aes(x = Association, y = rank)) + geom_line(aes(group = genes, colour = gene1)) + geom_text(aes(colour=gene1,label = genes, hjust = offset,size=0.01))+ theme_bw() +
  
  #eliminates background, gridlines, and chart border
  theme(legend.position="none",
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank(),
    axis.line=element_blank(),
    #axis.text.y=element_blank(),
    axis.ticks=element_blank()
    #axis.title.y=element_blank()
  ) 
require(ggplot2)
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#install.packages('OrderedList')
#biocLite(c("OrderedList"))
require(OrderedList)
OL.data
data(OL.data)
t1 = which(matrix(rAMatL %in% head(sort(rAMatL, T), 4371),nr = nrow(rAMatL)), arr.ind = TRUE)
l2 = which(matrix(rMICMatL %in% head(sort(rMICMatL, T), ctp),nr = nrow(rMICMatL)), arr.ind = TRUE)
l3 = which(matrix(rDcorMatL %in% head(sort(rDcorMatL, T), ctp),nr = nrow(rDcorMatL)), arr.ind = TRUE)

t2=colnames(genes)[l2[,1]]
compareLists(1:5,3:7,two.sided=F,no.reverse=T)
?compareLists
list1 <- t1
list2 <- sample(t1[1:20])
xx <- compareLists(t1[-c(1,28,44)],t2)
xx
getOverlap(x)
setdiff(t1,t2)
which(t1 == "X240805_at" , arr.ind=T)
which(t1 =="X200695_at" , arr.ind=T)
which(t1 =="X209189_at", arr.ind=T)
