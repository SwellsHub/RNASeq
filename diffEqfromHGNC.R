#Pathway Analysis with CAMERA

#import libraries
library(edgeR)
library(limma)
library(biomaRt)
library("CAMERA")
library(UpSetR)
library(plyr)



#Set Working directory
setwd("C:/Users/spenc/Downloads/Documents/Millie")

#Load count and phenotype label/parameter files to tables
x <-
  read.table(
    "millieGroup1Counts.txt",
    sep = "\t",
    header = TRUE,
    row.names = NULL,
  )
target <-
  read.table(
    "millieGroup1Key.csv",
    header = TRUE,
    sep = ",",
    row.names = 1
  )


#No files being imported below this point... Consider it the cogs of the machine
#-----------------------------------------------------------------------------------------------------

rownames(x) <- make.names(x$X, unique = TRUE)

x <- x[2:ncol(x)]

#x <- x[order(rownames(x)), , drop = FALSE]


#For Ensembl, trim version suffixes from IDs
#species <- as.character(target[1, 4])
#annotation <-  as.character(target[1, 5])

#genes <- as.matrix(rownames(x))
#colnames(genes) = "ensembl_gene_id"
#x <- x[order(rownames(x)), , drop = FALSE]
#x <- x[unique(sub("\\..*","", genes)), , drop = FALSE]
#genes <- unique(sub("\\..*","", genes))
#rownames(x) <- sub("\\..*","", genes)

#genes <- unique(sub("\\..*","", genes))

#ensembl <-
  #useEnsembl(biomart = ("ensembl"), dataset = "hsapiens_gene_ensembl")
#Conversion <-
 # getBM(
    #attributes = c('ensembl_gene_id', 'hgnc_symbol'),
    #filters = 'ensembl_gene_id',
    #values = genes,
    #mart = ensembl
  #)
#Conversion <-
  #Conversion[order(Conversion$ensembl_gene_id), , drop = FALSE]

#filter <- x

#Commit the replacement HGNC symbols and trim the data frame to prepare for differential expression
#filter <- data.frame(filter, row.names = rownames(filter))
#filter$annotation <- rownames(filter)
#names(filter)[ncol(filter)] <- annotation
#finalsheet <- merge(filter, Conversion, by = annotation, drop = FALSE)
#finalsheet <- finalsheet[!(finalsheet$hgnc_symbol == ""),]
#finalsheet <- na.omit(finalsheet)
#finalsheet <- finalsheet[order(finalsheet$hgnc_symbol), , drop = FALSE]
#hgnc_symbols <- unique(finalsheet$hgnc_symbol)
#finalsheet <- ddply(finalsheet,"hgnc_symbols",numcolwise(sum))

#rownames(finalsheet) <- hgnc_symbols
#finalsheet <- finalsheet[2:(ncol(finalsheet))]


#x <- finalsheet



#Make matrix of labels for design of experiment and intended pairwise comparisons

#Factor group column of target file to construct design matrix
design <- model.matrix(~target$Group + 0)

#Extract the relevant contrasts of the experiment from the third column of the target file
contrasts <- as.character(unlist(target[3]))
contrasts <- noquote(contrasts)
contrasts <- contrasts[contrasts != ""]
#Arrange columns of design matrix in alphabetical order
colnames(design) <- sort(as.character(levels(target$Group)))
contrast.matrix <-
  makeContrasts(contrasts = contrasts,  levels = design)


#Make DGEList from tables
x1 <- data.frame(x)
x1 <- x1[!(rownames(x1) == ""),]
x1 <- na.omit(x1)

dgeFull <- DGEList(counts = x1, target$Group)

#logCPM out low count genes
dgeFull <-
  DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
          group = dgeFull$samples$group)
normCounts <- cpm(dgeFull)
eff.lib.size <-
  dgeFull$samples$lib.size * dgeFull$samples$norm.factors

#Estimate Normalization factors
dgeFull <- calcNormFactors(dgeFull, "TMM")


#The Fourth parameter in the second column of the target file (first column is rownames and does not count as column) 
#establishes whether logCPM or Voom normalization is applied to data. A 0 is logCPM, a 1 is Voom
if (target[4, 2] == 0) {
  #Get Log counts per million
  logCPM <- cpm(dgeFull, log = TRUE, prior.count = 3)
}
#Or apply Voom Normalization
if (target[4, 2] == 1) {
  logCPM <- voomWithQualityWeights(dgeFull, design, plot = TRUE, normalization = "none")
}
#Apply and eBayes Linear trend
fit <- lmFit(logCPM, design)

##Apply to pairwise analysis
fit2 <- contrasts.fit(fit, contrast.matrix)




#Standard Method of fit calculation
# The next cell in the parameters column of the key file determines wheter eBayes or treat is applied to the linear model.
#A 0 is eBayes, a 1 is treat (treat incorporates LogFC into the calculation of adj.p.Val).
if (target[5, 2] == 0) {
  fit2 <- eBayes(fit2, trend = TRUE)
}

#Method for increased importance of logFC
if (target[5, 2] == 1) {
  fit2 <- treat(fit2, lfc = .5)
}

#Loop through writing for every pairwise comparison
for (i in seq(from = 1, to = ncol(contrast.matrix))) {
  #Write table
  name <- colnames(contrast.matrix)[i]
  
  #Standard Output
  #The cell that determines if the treat normalization is used also determines what the naming convention of the output file is.
  
  if (target[5, 2] == 0) {
    write.table(topTable(fit2, coef = i, number = 50000),
                paste(name, "LogFCs.csv"),
                sep = ",")
  }
  
  #Output assigning more important to LogFC
  if (target[5, 2] == 1) {
    write.table(
      topTreat(fit2, coef = i, number = 50000),
      paste(name, "LogFCs.csv"),
      sep = ",",
      col.names = NA,
      row.names = TRUE
    )
  }
  

#Filter data to GSEA Pre-ranked format
filter <-
  read.table(
    paste(name, "LogFCs.csv"),
    sep = ",",
    header = TRUE,
    row.names = 1
  )



#final formatting
filter <- filter[1]
colnames(filter) <- NULL

#write table
write.table(
  filter,
  file = paste(name, ".rnk", sep = ""),
  sep = "\t",
  col.names = NA,
  row.names = TRUE,
  quote = FALSE
)
}
