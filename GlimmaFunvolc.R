


#Import Libraries
library(edgeR)
library(limma)
library(stringi)
library(biomaRt)
library(Glimma)


#Set Working directory

setwd("C:/Users/spenc/downloads")



#This is where key file MUST be upoaded
target <-
  read.table(
    "hEAT RNAseq key file.csv",
    header = TRUE,
    sep = ",",
    row.names = 1
  )


#This is where counts file MUST be uploaded
filter <-
  read.table(
    "hEAT RNAseq 2018.txt",
    sep = "\t",
    header = TRUE,
    row.names = 1
  )


#Below this line, no files are being imported anymore... consider it the inside of the machine
#------------------------------------------------------------------------------------------------------------



filter <- filter[order(rownames(filter)), , drop = FALSE]


#For Ensembl, trim version suffixes from IDs
species <- as.character(target[1, 4])
annotation <-  as.character(target[1, 5])

genes <- as.matrix(rownames(filter))
colnames(genes) = annotation
genes <- unique(sub("\\..*", "", genes))
filter <-
  filter[unique(sub("\\..*", "", rownames(filter))), , drop = FALSE]
rownames(filter) <- sub("\\..*", "", genes)
filter <- filter[complete.cases(filter), ]


#Convert Ensembl IDs to HGNC symbols

ensembl <- useEnsembl(biomart = ("ensembl"), dataset = species)
Conversion <-
  getBM(
    attributes = c(as.character(annotation), 'hgnc_symbol'),
    filters = as.character(annotation),
    values = rownames(filter),
    mart = ensembl
  )
Conversion <-
  Conversion[order(Conversion[annotation]), , drop = FALSE]
filter <- data.frame(filter, row.names = rownames(filter))
filter$annotation <- rownames(filter)
names(filter)[ncol(filter)] <- annotation
finalsheet <- merge(filter, Conversion, by = annotation)
finalsheet <- finalsheet[!(finalsheet$hgnc_symbol == ""),]
finalsheet <- finalsheet[complete.cases(finalsheet), ]
finalsheet <-
  finalsheet[!duplicated(finalsheet$hgnc_symbol), , drop = FALSE]
rownames(finalsheet) <- finalsheet$hgnc_symbol
finalsheet <- finalsheet[2:(ncol(finalsheet) - 1)]



#Load count and phenotype label files to tables
x <- finalsheet



#Make matrix of labels for design of experiment and intended pairwise comparisons
design <- model.matrix(~ target$Group + 0)
contrasts <- as.character(unlist(target[3]))
contrasts <- noquote(contrasts)
contrasts <- contrasts[contrasts != ""]
colnames(design) <- sort(as.character(levels(target$Group)))
contrast.matrix <-
  makeContrasts(contrasts = contrasts,  levels = design)



#Make DGEList from tables
dgeFull <- DGEList(counts = x, target$Group)



#Filter out low count genes
dgeFull <-
  DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
          group = dgeFull$samples$group)
normCounts <- cpm(dgeFull)
eff.lib.size <-
  dgeFull$samples$lib.size * dgeFull$samples$norm.factors


#Estimate Normalization factors
dgeFull <- calcNormFactors(dgeFull, "TMM")



if (target[4, 2] == 0) {
  #Get Log counts per million
  logCPM <- cpm(dgeFull, log = TRUE, prior.count = 3)
}
#Or apply Voom Normalization
if (target[4, 2] == 1) {
  logCPM <- voomWithQualityWeights(dgeFull, design, plot = TRUE)
}
#Apply and eBayes Linear trend
fit <- lmFit(logCPM, design)

##Apply to pairwise analysis
fit2 <- contrasts.fit(fit, contrast.matrix)




#Standard Method of fit calculation
if (target[5, 2] == 0) {
  fit2 <- eBayes(fit, trend = TRUE)
}

#Method for increased importance of logFC
if (target[5, 2] == 1) {
  fit2 <- treat(fit2, lfc = log2(1.5))
}

dt <- decideTests(fit2)


#Optional grouping if no design matrix is present
#groups <- as.integer((factor(target$Group)))

#Interactive MA Plot
#glMDPlot(fit2, coef = 4, status=dt, counts=logCPM, groups=as.integer(target$Group), side.main="Symbols")

#for (i in seq(from = 1, to = ncol(contrast.matrix))) {
  name <- colnames(contrast.matrix)[i]
  write.table(
    topTreat(fit2, coef = i, number = 50000),
    paste(name, "treat.csv"),
    sep = ",",
    col.names = NA,
    row.names = TRUE
  )
  
  
  #Filter data to GSEA Pre-ranked format
  filter <-
    read.table(
      paste(name, "treat.csv"),
      sep = ",",
      header = TRUE,
      row.names = 1
    )
  
  res2 <- filter
  
  status <- data.frame(row.names = rownames(filter))
  
  pos <- res2[(res2$adj.P.Val < .05 & res2$logFC > 0), ]
  if(length(pos$logFC != 0)) {
    pos$result <- 1
    pos <- pos[6]
  }
  
  neg <- res2[(res2$adj.P.Val < .05 & res2$logFC < 0), ]
  if(length(neg$logFC != 0)) {
    neg$result <- -1
    neg <- neg[6]
  }
  
  zeroes <- res2[res2$adj.P.Val >= .05, ]
  if(length(zeroes$logFC != 0)) {
    zeroes$result <- 0
    zeroes <- zeroes[6]
  }
  
  rownames(status) <-
    append(rownames(pos), c(rownames(neg), rownames(zeroes)))
  status$stats <- append(pos$result, c(neg$result, zeroes$result))
  
  #Make sure data from all inputs are in the same order
  logCPM$E <- logCPM$E[order(rownames(logCPM$E)), , drop = FALSE]
  filter <- filter[order(rownames(filter)), , drop = FALSE]
  status <- status[order(rownames(status)), , drop = FALSE]
  
  glXYPlot(
    x = filter$logFC,
    y = -log10(filter$adj.P.Val),
    status = status$stats,
    counts = logCPM,
    groups = as.integer(target$Group),
    xlab = "Log2 Fold Change",
    ylab = "-log10PFDR",
    side.main = "Symbols",
    html = name,
    lwd = 3
  )
  
  glMDPlot(fit2, coef = 1, status=dt, counts=logCPM, groups=as.integer(target$Group), side.main="Symbols")
#}