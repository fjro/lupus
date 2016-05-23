
library(RColorBrewer)
require(devtools)
library(rafalib)

mypar2(1,1)


dim(subBD)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)


library(gplots)

#cols = colorRampPalette(brewer.pal(9, "RdYlBu"))(19)
colors <- brewer.pal(11, "Spectral")

pal <- colorRampPalette(colors)
display.brewer.all()

dim(subBD)

bigStatus = c(rep('Control', nrow(safeControlBD)), rep('Lupus', nrow(safeLupusBD)))
cols <- palette(brewer.pal(3, "RdYlBu"))[as.fumeric(as.character(bigStatus))]


heatmap.2(as.matrix((subBD)), labRow=bigStatus, 
          trace="none", 
          RowSideColors=cols, 
          col=hmcol)
subBD = safeBD
require(genefilter)
?varFilter
heatmap.2(t(as.matrix(t(subBD))), trace="none", density="none", col=hmcol,
          scale="row",
          labCol = "",
          labRow="", 
          #ColSideColors = rows,
          RowSideColors=cols, 
          #labCol = annot$CHR,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-(cor(t(x))^2) )))


