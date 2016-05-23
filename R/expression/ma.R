require(energy)
require(matie)
require(minerva)

#calculate associations
#pearson
allCM <- cor(as.matrix(genes[,-c(1,96,97)]))
controlCM <- cor(as.matrix(controls[,-c(1,96,97)]))
lupusCM <- cor(as.matrix(lupus[,-c(1,96,97)]))

#dcor
allDcor <- tadcor(genes[,-c(1,96,97)])
controlDcor <- tadcor(controls[,-c(1,96,97)])
lupusDcor <- tadcor(lupus[,-c(1,96,97)])

#A
allA <- tap(genes[,-c(1,96,97)])
diag(allA) <- 1
controlA <- tap(genes[genes$diseased=="Control",-c(1,96,97)])
diag(controlA) <- 1
lupusA <- tap(genes[genes$diseased=="Lupus",-c(1,96,97)])
diag(lupusA) <- 1

#MIC
allMIC <- mine(as.matrix(genes[,-c(1,96,97)]))$MIC
diag(allMIC) <- 1
controlMIC <- mine(as.matrix(controls[,-c(1,96,97)]))$MIC
diag(controlMIC) <- 1
lupusMIC <- mine(as.matrix(lupus[,-c(1,96,97)]))$MIC
diag(lupusMIC) <- 1

#Write to file
write.csv(allCM, "allCM.csv")
write.csv(lupusCM, "lupusCM.csv")
write.csv(controlCM, "controlCM.csv")

write.csv(allDcor, "allDcor.csv")
write.csv(lupusDcor, "lupusDcor.csv")
write.csv(controlDcor, "controlDcor.csv")

write.csv(allA, "allA.csv")
write.csv(lupusA, "lupusA.csv")
write.csv(controlA, "controlA.csv")

write.csv(allMIC, "allMIC.csv")
write.csv(lupusMIC, "lupusMIC.csv")
write.csv(controlMIC, "controlMIC.csv")

#read prepreped
allCM = read.csv("allCM.csv")
lupusCM = read.csv("lupusCM.csv")
controlCM = read.csv("controlCM.csv")

allDcor = read.csv("allDcor.csv")
lupusDcor = read.csv("lupusDcor.csv")
controlDcor = read.csv("controlDcor.csv")

allA = read.csv("allA.csv")
lupusA = read.csv("lupusA.csv")
controlA = read.csv("controlA.csv")

allMIC = read.csv("allMIC.csv")
lupusMIC = read.csv("lupusMIC.csv")
controlMIC = read.csv("controlMIC.csv")

require(plyr)
allCM = allCM[,-c(1,2)]
colnames(allCM) = sub('X','', colnames(allCM) )
rownames(allCM) = colnames(allCM)



fixCSV = function(x)
{
  x = x[,-1]
  colnames(x) = sub('X','', colnames(x) )
  rownames(x) = colnames(x)
  x
}

allA = fixCSV(allA)
controlA = fixCSV(controlA)
lupusA = fixCSV(lupusA)

allDcor = fixCSV(allDcor)
controlDcor = fixCSV(controlDcor)
lupusDcor = fixCSV(lupusDcor)

allMIC= fixCSV(allMIC)
controlMIC = fixCSV(controlMIC)
lupusMIC = fixCSV(lupusMIC)



#calculate residual associations
#pearson

#dcor
allDcor.NL <- nltadcor(genes[,-c(1,96,97)])
controlDcor.NL <- nltadcor(controls[,-c(1,96,97)])
lupusDcor.NL <- nltadcor(lupus[,-c(1,96,97)])

#A
allA.NL <- nltap(genes[,-c(1,96,97)])
diag(allA.NL) <- 1
controlA.NL <- nltap(genes[genes$diseased=="Control",-c(1,96,97)])
diag(controlA.NL) <- 1
lupusA.NL <- nltap(genes[genes$diseased=="Lupus",-c(1,96,97)])
diag(lupusA.NL) <- 1

#MIC
allMIC.NL <- nltaMIC(genes[,-c(1,96,97)])
diag(allMIC.NL) <- 1
controlMIC.NL <- nltaMIC(controls[,-c(1,96,97)])
diag(controlMIC.NL) <- 1
lupusMIC.NL <- nltaMIC(lupus[,-c(1,96,97)])
diag(lupusMIC.NL) <- 1

#Write to file
rownames(controlA.NL) = rownames(controlA)
rownames(controlDcor.NL) = rownames(controlDcor)
rownames(controlMIC.NL) = rownames(controlMIC)

rownames(lupusA.NL) = rownames(controlA)
rownames(lupusDcor.NL) = rownames(controlDcor)
rownames(lupusMIC.NL) = rownames(controlMIC)

write.csv(allDcor, "allDcor.csv")
write.csv(lupusDcor, "lupusDcor.csv")
write.csv(controlDcor, "controlDcor.csv")

write.csv(allA, "allA.csv")
write.csv(lupusA, "lupusA.csv")
write.csv(controlA, "controlA.csv")

write.csv(allMIC, "allMIC.csv")
write.csv(lupusMIC, "lupusMIC.csv")
write.csv(controlMIC, "controlMIC.csv")
