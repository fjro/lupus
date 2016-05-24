library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e)
rownames(genes) = genes[,1]
colnames(genes)
d <- dist(genes[,-c(1,96)])


library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

library(genefilter)
rv <- rowVars(e)
idx <- order(-rv)[1:40]
heatmap(as.matrix(genes[,-c(1,96,97)]), col=hmcol)

library(gplots)
cols <- palette(brewer.pal(9, "RdYlBu"))[as.fumeric(as.character(genes[,96]))]

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





require(matie)
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
