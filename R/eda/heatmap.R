library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e)
rownames(genes) = genes[,1]
colnames(genes)
d <- dist(genes[,-c(1,96)])

require(devtools)
#install_github("ririzarr/rafalib")
library(rafalib)
mypar2(1,1)
hc <- hclust(d)

plot(hc,labels=genes[,96])
as.character()
myplclust(hc, labels=genes[,96], lab.col=as.fumeric(as.character(genes[,96])))

abline(h=12)
hclusters <- cutree(hc, h=10)
table(true=genes[,96], cluster=hclusters)

#k-means
plot(e[1,], e[2,])
set.seed(1)
km <- kmeans(t(e[1:2,]), centers=7)
names(km)
plot(e[1,], e[2,], col=km$cluster, pch=16)
plot(e[1,], e[2,], col=as.fumeric(tissue), pch=16)
table(true=tissue,cluster=km$cluster)

#kmeans with cmd
mds <- cmdscale(d)
plot(mds[,1], mds[,2]) 
km <- kmeans(genes[,-c(1)], centers=7)
plot(mds[,1], mds[,2], col=km$cluster, pch=16)
table(true=tissue,cluster=km$cluster)

# install.packages("RColorBrewer")
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

library(genefilter)
rv <- rowVars(e)
idx <- order(-rv)[1:40]
heatmap(as.matrix(genes[,-c(1,96,97)]), col=hmcol)

# install.packages("gplots")
library(gplots)
#cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(as.character(genes[,96]))]
cols <- palette(brewer.pal(9, "RdYlBu"))[as.fumeric(as.character(genes[,96]))]

#cols = colorRampPalette(brewer.pal(9, "RdYlBu"))(19)

head(cbind(rownames(genes),cols))
heatmap.2(as.matrix(t(genes[,de_genes])), labCol=genes[,96], 
          trace="none", 
          ColSideColors=cols, 
          col=hmcol)

samps = substr(samples$Sample,0,4)
#cols <- palette(brewer.pal(8, "Dark2"))[as.fumeric(samps)]
head(cbind(rownames(genes),cols))
heatmap.2(as.matrix(t(genes[,-c(1,96)])), labCol=genes[,96], 
          trace="none", 
          ColSideColors=cols, 
          col=hmcol)
table(samps)

plot(hclust(as.dist(aa)))
annot = read.csv('annot.csv')

colors <- brewer.pal(11, "Spectral")

pal <- colorRampPalette(colors)
display.brewer.all()

#rows = palette(pal(20))[as.fumeric(as.character(annot$CHR))]
rows = colorRampPalette(brewer.pal(9, "RdYlBu"))(20)[as.fumeric(as.character(annot$CHR))]

?palette

heatmap.2(t(as.matrix(t(genes[,-c(1,96)]))), trace="none", density="none", col=hmcol,
          scale="row",
          labRow="", 
          ColSideColors = rows,
        RowSideColors=cols, 
        labCol = annot$CHR,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2))
require(energy)
heatmap.2(t(as.matrix(t(genes[,-c(1,96)]))), trace="none", density="none", col=hmcol, dendrogram='col',
          scale="row",
          labRow="", 
          ColSideColors = rows,
          RowSideColors=cols, 
          labCol = annot$CHR,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2))

# par(lend = 1)           # square line ends for the color legend
# legend("bottomleft",      # location of the legend on the heatmap plot
#        legend = unique(as.fumeric(as.character(annot$CHR))), # category labels
#        col = rows,  # color key
#        lty= 1,             # line style
#        lwd = 10            # line width
# )



require(matie)
?tap
length(rows)
length(cols)

heatmap.2(t(as.matrix(t(genes[,-c(1,96)]))), trace="none", density="none", col=hmcol, dendrogram='col',
          scale="row",
          labRow="", 
          ColSideColors = rows,
          RowSideColors=cols, 
          labCol = annot$CHR,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-tap(data.frame(t(x))))))

heatmap.2(t(as.matrix(t(genes[,-c(1,96,97)]))), trace="none", density="none", col=hmcol,
          scale="row",
          labRow="", 
          ColSideColors = rows,
          RowSideColors=cols, 
          labCol = annot$CHR,
          hclust=function(x) hclust(x,method="complete"))
tap(data.frame(as.matrix(genes[,-c(1,96)])))
