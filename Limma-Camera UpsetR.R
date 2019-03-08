#Pathway Analysis with CAMERA, followed by Overlap Analysis with UpSetR

#Pathway Analysis with CAMERA

#import libraries
library(edgeR)
library(limma)
library(biomaRt)
library("CAMERA")
library(UpSetR)



#Set Working directory
setwd("C:/Users/spenc/Downloads")

#Load count and phenotype label/parameter files to tables
x <-
  read.table(
    "hEAT RNAseq 2018.txt",
    sep = "\t",
    header = TRUE,
    row.names = 1
  )
target <-
  read.table(
    "hEAT RNAseq key file.csv",
    header = TRUE,
    sep = ",",
    row.names = 1
  )


#No files being imported below this point... Consider it the cogs of the machine
#-----------------------------------------------------------------------------------------------------

x <- x[order(rownames(x)), , drop = FALSE]


#For Ensembl, trim version suffixes from IDs
species <- as.character(target[1, 4])
annotation <-  as.character(target[1, 5])

genes <- as.matrix(rownames(x))
colnames(genes) = "ensembl_gene_id"
x <- x[order(rownames(x)), , drop = FALSE]
x <- x[unique(sub("\\..*","", genes)), , drop = FALSE]
genes <- unique(sub("\\..*","", genes))
rownames(x) <- sub("\\..*","", genes)

genes <- unique(sub("\\..*","", genes))

ensembl <-
  useEnsembl(biomart = ("ensembl"), dataset = species)
Conversion <-
  getBM(
    attributes = c('ensembl_gene_id', 'entrezgene'),
    filters = 'ensembl_gene_id',
    values = genes,
    mart = ensembl
  )
Conversion <-
  Conversion[order(Conversion$ensembl_gene_id), , drop = FALSE]

filter <- x

#Commit the replacement HGNC symbols and trim the data frame to prepare for differential expression
filter <- data.frame(filter, row.names = rownames(filter))
filter$annotation <- rownames(filter)
names(filter)[ncol(filter)] <- annotation
finalsheet <- merge(filter, Conversion, by = annotation, drop = FALSE)
finalsheet <- finalsheet[!(finalsheet$entrezgene == ""),]
finalsheet <- na.omit(finalsheet)
finalsheet <- finalsheet[order(finalsheet$entrezgene), , drop = FALSE]
entrezgenes <- unique(finalsheet$entrezgene)
finalsheet <- ddply(finalsheet,"entrezgene",numcolwise(sum))

rownames(finalsheet) <- entrezgenes
finalsheet <- finalsheet[2:(ncol(finalsheet))]


x <- finalsheet



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
dgeFull <- DGEList(counts = x, target$Group)

#logCPM out low count genes
dgeFull <-
  DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
          group = dgeFull$samples$group)
normCounts <- cpm(dgeFull)
eff.lib.size <-
  dgeFull$samples$lib.size * dgeFull$samples$norm.factors

#Estimate Normalization factors
dgeFull <- calcNormFactors(dgeFull, "TMM")

#Normalize Data with voom
logCPM <-
  voom(
    dgeFull,
    design,
    plot = TRUE,
    header = FALSE,
    row.names = FALSE
  )

#Download and load desired gene sets
download.file(
  "http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata",
  "human_c2_v5p2.rdata",
  mode = "wb"
)
load("human_c2_v5p2.rdata")

#prep camera test with gene set and index of probe identifiers

cIndex <- ids2indices(Hs.c2, rownames(logCPM))



#Create vectors to coerce significant sets to usable list
xUp <- vector("list", ncol(contrast.matrix))
xDown <- vector("list", ncol(contrast.matrix))
LUp <- vector("list", ncol(contrast.matrix))
LDown <- vector("list", ncol(contrast.matrix))

#Conduct Camera test for each comparison
for (i in seq(from = 1, to = ncol(contrast.matrix))) {
  name <- colnames(contrast.matrix)[i]
  cameraTest <-
    camera(
      logCPM,
      index = cIndex,
      design = design,
      contrast = contrast.matrix[, i],
      inter.gene.cor = 0.01
    )
  
  
  cameraTest$names <- rownames(cameraTest)
  
  #Add upregulated sets to up array
  Pu <-
    cameraTest[(cameraTest$Direction == "Up" &
                  cameraTest$FDR < target[6, 2]), ]
  xUp[[i]] <- Pu
  xNamesU <- as.matrix(xUp[[i]]$names)
  LUp[[i]] <- xNamesU
  names(LUp)[[i]] <- paste(colnames(contrast.matrix)[i], "Up")
  
  #Add downregulated sets to down array
  Pd <-
    cameraTest[(cameraTest$Direction == "Down" &
                  cameraTest$FDR < target[6, 2]), ]
  xDown[[i]] <- Pd
  xNamesD <- as.matrix(xDown[[i]]$names)
  LDown[[i]] <- xNamesD
  names(LDown)[[i]] <- paste(colnames(contrast.matrix)[i], "Down")
  
  #Write Camera results to excel
  if (target[7, 2] == 1) {
    write.table(
      cameraTest,
      file = paste(
        "C:/Users/wellsse/Documents/CorrectDirectionDE/CameraPathwaysRightDirection/CameraChecks/Camera.04sh",
        name,
        ".csv"
      ),
      sep = ",",
      quote = FALSE,
      col.names = NA,
      row.names = TRUE
    )
  }

}
#Read which overlap comparisons are to be graphically represented
userInput <-
  read.table(file = "SomeCoolComps.csv",  sep = ",", fill = TRUE)
x <- 0


#Prepare the comparisons specified by user input for Overlap Analysis
for (l in seq(from = 1, to = ncol(userInput))) {
  InputTable <- vector("list", sum(userInput[, l] != ""))
  k <- 1
  for (h in seq(from = 1, to = sum(userInput[, l] != ""))) {
    for (j in seq(from = 1, to = ncol(contrast.matrix))) {
      if (identical(names(LUp[j]), toString(userInput[h, l]))) {
        InputTable[k] <- LUp[j]
        names(InputTable)[[k]] <- names(LUp[j])
        k <- k + 1
      }
      if (identical(names(LDown[j]), toString(userInput[h, l]))) {
        InputTable[k] <- LDown[j]
        names(InputTable)[[k]] <- names(LDown[j])
        k <- k + 1
      }
    }
  }
  #Construct graph to represent overlaps
  sets <- names(InputTable)
  setsd <- data.frame("names" = sets)
  
  #Subset to last two characters of each element
  setsdSub <-
    substr(setsd$names, nchar(as.character(setsd$names)) - 1, nchar(as.character(setsd$names)))
  
  
  #Coerce Data frames to greppable metadata for UpsetR (distinguish based on Up/Down)
  setlevels <- factor(setsdSub)
  setlevelsd <- data.frame("names" = setlevels)
  metadata <- as.data.frame(cbind(sets, setlevelsd))
  
  
  
  #Extract lists of Up and Down Regulated Sets to prepare calculation of overlaps
  grepgdown <- InputTable[grep("Down", setsd$names, value = TRUE)]
  grepgup <-  InputTable[grep("Up", setsd$names, value = TRUE)]
  n <- length(grepgdown)
  picName = paste(
    "C:/Users/spenc/Downloads/Documents/OverlapCharts/",
    toString(x),
    "Overlaps.png"
  )
  
  
  #Designate Overlaps for Excel sheet Generation
  Overlaps <-
    list(
      Down_Down = intersect(grepgdown[[1]], grepgdown[[2]]),
      Up_Up = intersect(grepgup[[2]], grepgup[[1]]),
      Down_Up = intersect(grepgdown[[1]], grepgup[[2]]),
      Up_Down = intersect(grepgdown[[2]], grepgup[[1]])
    )
  names(Overlaps) <-
    list(
      paste(names(InputTable[2]), "_", names(InputTable[4])),
      paste(names(InputTable[1]), "_", names(InputTable[3])),
      paste(names(InputTable[2]), "_", names(InputTable[3])),
      paste(names(InputTable[1]), "_", names(InputTable[4]))
    )
  noOverlaps <-
    list(
      setdiff(grepgup[[1]], c(grepgdown[[2]], grepgup[[2]])),
      setdiff(grepgdown[[1]], c(grepgdown[[2]], grepgup[[2]])),
      setdiff(grepgup[[2]], c(grepgdown[[1]], grepgup[[1]])),
      setdiff(grepgdown[[2]], c(grepgdown[[1]], grepgup[[1]]))
    )
  if (length(names(Overlaps)) == length(names(InputTable))) {
    names(noOverlaps) <- names(InputTable)
    breakdown <- append(Overlaps, noOverlaps)
    Overlapdf <-
      data.frame(lapply(breakdown, 'length<-', max(lengths(Overlaps))))
    
    #Write an Excel specifying which pathways correspond to which overlap
    write.csv(
      Overlapdf,
      paste(
        "C:/Users/spenc/Downloads/Documents/OverlapExcels/",
        toString(x),
        "dissection.csv"
      ),
      na = ''
    )
  }
  
  
  #Determine how to generate the plot based on what comparisons exist in data
  if ((length(intersect(grepgdown[[1]], grepgdown[[2]])) != 0) &
      (length(intersect(grepgup[[1]], grepgup[[2]])) != 0 &
       length(grepgup) <= 2)) {
    upset(
      fromList((InputTable)),
      sets = names(InputTable),
      keep.order = T,
      set.metadata = list(data = metadata, plots = list(
        list(
          type = "matrix_rows",
          column = "names",
          colors = c(Up = "navy", wn = "black"),
          alpha = .4
        )
      )),
      queries = list(
        list(
          query = intersects,
          params = grep("Up", setsd$names, value = TRUE),
          color = "blue",
          active = T
        ),
        list(
          query = intersects,
          params = grep("Down", setsd$names, value = TRUE),
          color = "blue",
          active = T
        )
      ),
      nsets = length(InputTable),
      text.scale = 2,
      main.bar.color = "black",
      mainbar.y.label = "number of sets"
    )
    
  } else if (length(intersect(grepgup[[1]], grepgup[[2]])) != 0 &
             length(intersect(grepgdown[[1]], grepgdown[[2]])) == 0 &
             length(grepgup) <= 2) {
    upset(
      fromList((InputTable)),
      sets = names(InputTable),
      keep.order = T,
      set.metadata = list(data = metadata, plots = list(
        list(
          type = "matrix_rows",
          column = "names",
          colors = c(Up = "navy", wn = "black"),
          alpha = .4
        )
      )),
      queries = list(
        list(
          query = intersects,
          params = grep("Up", setsd$names, value = TRUE),
          color = "blue",
          active = T
        )
      ),
      nsets = length(InputTable),
      text.scale = 2,
      main.bar.color = "black",
      mainbar.y.label = "number of sets"
    )
  } else if (length(intersect(grepgdown[[1]], grepgdown[[2]])) != 0 &
             length(intersect(grepgup[[1]], grepgup[[2]])) == 0 &
             length(grepgup) <= 2) {
    upset(
      fromList((InputTable)),
      sets = names(InputTable),
      keep.order = T,
      set.metadata = list(data = metadata, plots = list(
        list(
          type = "matrix_rows",
          column = "names",
          colors = c(Up = "navy", wn = "black"),
          alpha = .4
        )
      )),
      queries = list(
        list(
          query = intersects,
          params = grep("Down", setsd$names, value = TRUE),
          color = "blue",
          active = T
        )
      ),
      nsets = length(InputTable),
      text.scale = 2,
      main.bar.color = "black",
      mainbar.y.label = "number of sets"
    )
    
  } else if (length(grepgup) > 2) {
    upset(
      fromList((InputTable)),
      sets = names(InputTable),
      keep.order = T,
      set.metadata = list(data = metadata, plots = list(
        list(
          type = "matrix_rows",
          column = "names",
          colors = c(Up = "navy", wn = "black"),
          alpha = .4
        )
      )),
      nsets = length(InputTable),
      text.scale = 2,
      main.bar.color = "black",
      mainbar.y.label = "number of sets"
    )
    
  }
  
  
  #Ouput graphs
  dev.print(png,
            filename = picName,
            width = 862*2,
            height = 622*2)
  x <- x + 1
  
  
}
