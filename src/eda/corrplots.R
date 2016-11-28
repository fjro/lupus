#draw corrplots to visualise differences in various measures of gene association 
#between lupus patients and controls
library(corrplot)
corrplot(allCM, order <- "hclust", method="color")
corrplot(as.matrix(aa), method="color", type='upper')
par(mfrow=c(3,2))
colnames(controlCM) <- 1:94
corrplot(controlCM, title="", order = "original", method="color", tl.cex = .3, type='upper', tl.pos = 'd')

colnames(lupusCM) <- 1:94
corrplot(lupusCM, title="", order = "original", method="color",  tl.cex = .3, type='upper', tl.pos = 'd')
#controls vs lupus appear different, clusters of high positive and negative correlations

#'A' based corrplot
vv <- controlA
colnames(vv) <- 1:94
vv <- as.matrix(vv) # sign adjust using Pearson correlation to indicate direction of association
dP <- vv * sign(controlCM)
corrplot(dP, title="", order = "original", method="color", tl.cex = .3, type='upper', tl.pos = 'd')

bb <- lupusA
colnames(bb) = 1:94
bb <- as.matrix(bb)
dP <- bb * sign(lupusCM)
corrplot(dP, title="", order = "original", method="color", tl.cex = .3, type='upper', tl.pos = 'd')

#dcor plots
dvv <- tadcor(genes[genes$diseased=='Control',-c(1,96,97)])
which(dP > 1)
colnames(dvv) <- 1:94
#correlation adjust
dP <- dvv * sign(controlCM)
corrplot(dP, title="", order = "original", method="color", tl.cex = .3, type='upper', tl.pos = 'd')

dbb = tadcor(genes[genes$diseased=='Lupus',-c(1,96,97)])
colnames(dbb) = 1:94
dP <- dbb * sign(lupusCM)
corrplot(dP, title="", order = "original", method="color", tl.cex = .3, type='upper', tl.pos = 'd')

controlDcor <- tadcor(genes[genes$diseased=='Control',-c(1,96,97)])
lupusDcor <- tadcor(genes[genes$diseased=='Lupus',-c(1,96,97)])
