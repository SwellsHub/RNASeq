#Pathway Analysis with CAMERA

#import libraries
library(edgeR)
library(limma)
library(biomaRt)
library(CAMERA)
library(limma)


#Set Working directory
setwd("C:/Users/wellsse/Documents/CorrectDirectionDE/CameraPathwaysRightDirection")

#Load count and phenotype label files to tables
x<-read.table("hEAT RNAseq 2018.txt", sep = "\t", header = TRUE, row.names = 1)
target<-read.table("hEAT RNAseq key file.csv",header=TRUE,sep=",", row.names = 1)


#Make matrix of labels for design of experiment and intended pairwise comparisons
design <- model.matrix(~target$Group + 0 )
colnames(design) <- c("FM", "FshL1", "FshScr", "FT", "FY", "Hypox", "IM", "IshL1", "IshScr", "IT", "IY", "Normox")
contrast.matrix <- makeContrasts(IM-FM, FT-FM, FY-FM,FshL1-FshScr, IshL1-IM, IY-IM, IT-IM, IshL1-FshL1, IT-FT, IY-FY, IshL1-IshScr, levels = design)

#Filter bad Data
genes <- as.matrix(rownames(x))
colnames(genes) = "ensembl_gene_id"
x <- x[order(rownames(x)), , drop = FALSE]
x <- x[unique(sub("\\..*","", genes)), , drop = FALSE]
genes <- unique(sub("\\..*","", genes))
rownames(x) <- sub("\\..*","", genes)

genes <- unique(sub("\\..*","", genes))
ensembl <- useEnsembl(biomart = ("ensembl"), dataset = "hsapiens_gene_ensembl" )
Conversion <- getBM(attributes=c('ensembl_gene_id', 'entrezgene' ),filters = 'ensembl_gene_id', values = genes, mart = ensembl)
Conversion <- Conversion[order(Conversion$ensembl_gene_id), , drop = FALSE]

#Map Ensembl to HGNC
x <- cbind(x, ensembl_gene_id = rownames(x))
finalsheet <- merge(x, Conversion, by = "ensembl_gene_id", drop = FALSE)
finalsheet <- finalsheet[!(finalsheet$entrezgene ==""),]
finalsheet <- na.omit(finalsheet)
finalsheet <- finalsheet[!duplicated(finalsheet$entrezgene), , drop = FALSE]
rownames(finalsheet) <- finalsheet$entrezgene
finalsheet$entrezgene <- NULL
finalsheet$ensembl_gene_id <- NULL
x <- finalsheet


#Make DGEList from tables
dgeFull<-DGEList(counts=x, target$Group)

#logCPM out low count genes
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)
normCounts <- cpm(dgeFull)
eff.lib.size <- dgeFull$samples$lib.size*dgeFull$samples$norm.factors

#Estimate Normalization factors
dgeFull <- calcNormFactors(dgeFull, "TMM")

#Normalize Data with voom
logCPM <- voom(dgeFull, design, plot=TRUE, header = FALSE, row.names = FALSE)

#Downloard and load desired gene sets
download.file("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata", 
              "human_c2_v5p2.rdata", mode = "wb")
load("human_c2_v5p2.rdata")

#prep camera test with gene set and index of probe identifiers

cIndex <- ids2indices(Hs.c2, rownames(logCPM))
for(i in seq(from = 1, to = ncol(contrast.matrix) )) {
cameraTest <- camera(logCPM, index = cIndex, design = design, contrast = contrast.matrix[,i], inter.gene.cor=0.01)
name <-colnames(contrast.matrix)[i]
write.table(cameraTest, file = paste(name, "Camera.01.csv"), sep = ",", quote = FALSE, col.names = NA, row.names = TRUE)
}
cameraTest[1:5,]
