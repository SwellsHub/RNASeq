### R code from vignette source 'Glimma.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: basicanalysis
###################################################
library(edgeR)
library(limma)
library(Glimma)
data(lymphomaRNAseq)
x <- lymphomaRNAseq
class(x)

## Filter out genes with low counts
sel = rowSums(cpm(x$counts)>0.5)>=3
x = x[sel,]

## Make MDS plot
genotype = relevel(x$samples$group, "Smchd1-null")
par(mfrow=c(1,2))
plotMDS(x, label=1:ncol(x), main="MDS plot", col=as.numeric(genotype))
legend("topright",legend=c("Smchd1-null","Wild Type"), pch=20, col=1:2, text.col=1:2)

## Normalize the data using TMM
x = calcNormFactors(x, method="TMM")

## Set up design matrix
des = model.matrix(~genotype)
des

## Apply voom with sample quality weights and fit linear model
v=voomWithQualityWeights(x, design=des, normalization="none", plot=FALSE)
vfit = lmFit(v,des)

## Apply treat relative to a fold-change of 1.5
vtfit=treat(vfit,lfc=log2(1.5))
vfit= eBayes(vfit)
results = decideTests(vfit,p.value=0.01)
summary(results)

## Make a mean-difference (MD) plot of the results
plotMD(vfit,col=2, status=results[,2], hl.col=c("red", "blue"), 
       legend="topright", main="MD plot: Wild-type vs Smchd1")


###################################################
### code chunk number 3: interactiveMDS
###################################################
glMDSPlot(x, labels=1:7, groups=genotype, folder="Smchd1-Lymphoma", launch=FALSE)


###################################################
### code chunk number 4: interactiveMD
###################################################
glMDPlot(vfit, counts=x$counts, anno=x$genes, groups=genotype, samples=1:7, status=results[,2],
         display.columns=c("Symbols", "GeneID", "GeneName"), folder="Smchd1-Lymphoma",
         main="MD plot: Wild-type vs Smchd1", launch=FALSE)


###################################################
### code chunk number 5: interactivevolcano
###################################################
topt = topTable(vfit, coef=2, number=Inf, sort="none")
glMDPlot(topt, xval="logFC", yval="B", counts=x$counts, anno=x$genes, groups=genotype,
         samples=1:7, status=results[,2], display.columns=c("Symbols", "GeneID", "GeneName"), 
         folder="Smchd1-Lymphoma", html="volcano", launch=FALSE)


###################################################
### code chunk number 6: sessioninfo
###################################################
sessionInfo()


