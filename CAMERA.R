### R code from vignette source 'CAMERA.Rnw'

###################################################
### code chunk number 1: EICPspec1
###################################################
library(CAMERA)
file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
xs   <- xcmsSet(file, method="centWave",ppm=30, peakwidth=c(5,10))
an   <- xsAnnotate(xs)
an   <- groupFWHM(an)
an   <- findAdducts(an, polarity="positive")
plotEICs(an, pspec=2, maxlabel=5)


###################################################
### code chunk number 2: Pspec1
###################################################
plotPsSpectrum(an, pspec=2, maxlabel=5)


