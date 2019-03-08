
library(SnowballC)
library(wordcloud)
library(RColorBrewer)
library(RCurl)
library(XML)


setwd(
  "C:/Users/spenc/Downloads/Documents/Millie/pathwayAnalysis/c2"
)
path <-
  (
    "C:/Users/spenc/Downloads/Documents/Millie/PathwayAnalysis/c2"
  )
file.names <- dir(path, pattern = ".csv")
blackList <-
  read.table(
    "C:/Users/spenc/Downloads/Blacklist.csv",
    sep = ","
  )




#Masterlist <- list()

for (p in 1:length(file.names)) {
  file <- read.table(file.names[p], header = TRUE, sep = ",")
  file <- file[file$FDR < .05,]
  file1 <- file[order(file$X), , drop = FALSE]
  
  
  Modifiers <-
    read.csv(
      "C:/Users/spenc/Downloads/realModifiers.csv",
      stringsAsFactors = FALSE
    )
  Modifiers$Original <- tolower(Modifiers$Original)
  Modifiers$Original <- gsub("(^|_)([[:alpha:]])", "\\1\\U\\2", Modifiers$Original, perl=TRUE)
  
  #Get Positive List
  Positive <- file1[file1$Direction == "Up", ]

  Posnames1 <- Positive$X
  
  Posnames1 <- tolower(Posnames1)
  Posnames1 <- gsub("(^|_)([[:alpha:]])", "\\1\\U\\2", Posnames1, perl=TRUE)
  Posnames1sub <- Posnames1
  
  for (i in 1:nrow(Modifiers)) {
    Posnames1sub <-
      gsub(
        paste0("_", Modifiers[i, 1], "_"),
        paste0("_", Modifiers[i, 2], "_"),
        Posnames1sub
        
      )
  }
    
    for (i in 1:nrow(Modifiers)) {
      Posnames1sub <-
        gsub(
          paste0("_", Modifiers[i, 1], "$"),
          paste0("_", Modifiers[i, 2], ""),
          Posnames1sub
          
        )
    
    
  }
  
  Posnames1sub <- sub('^.*?_', '', Posnames1sub)
  Posnames1sub <- sub('_UP$', '', Posnames1sub)
  Posnames1sub <- sub('_DN$', '', Posnames1sub)
  
  Posnames1sub <- paste(Posnames1sub, collapse = " ")
  
  
  
  Posnames2 <- strsplit(Posnames1sub, split = '[_ ]+')
  PosVnames <- as.vector(Posnames2[[1]])
  VBLnames <- as.vector(tolower(blackList[[1]]))
  VBLnames <- gsub("(^|_)([[:alpha:]])", "\\1\\U\\2", VBLnames, perl=TRUE)
  
  Posnames3 <- PosVnames[!PosVnames %in% VBLnames]
  PosVnames2 <- unique(Posnames3)
  PosVnames2 <- sort(PosVnames2)
  
  
  
  
  #Get Negative List
  Negative <- file1[file1$Direction == "Down", ]

  Negnames1 <- Negative$X
  
  Negnames1 <- tolower(Negnames1)
  Negnames1 <- gsub("(^|_)([[:alpha:]])", "\\1\\U\\2", Negnames1, perl=TRUE)
  Negnames1sub <- Negnames1
  
  
#Carry out substitutions
  for (i in 1:nrow(Modifiers)) {
    Negnames1sub <-
      gsub(
        paste0("_", Modifiers[i, 1], "_"),
        paste0("_", Modifiers[i, 2], "_"),
        Negnames1sub
        
      )
    
    
  }
  
#Catch substitutions missed by the first passthrough
  for (i in 1:nrow(Modifiers)) {
    Negnames1sub <-
      gsub(
        paste0("_", Modifiers[i, 1], "$"),
        paste0("_", Modifiers[i, 2], ""),
        Negnames1sub
      
      )
    
    
  }
  
  #Filter out the first word (Author/source) of each pathway, as well as any 'Up's or 'DN's
  Negnames1sub <- sub('^.*?_', '', Negnames1sub)
  Negnames1sub <- sub('_UP$', '', Negnames1sub)
  Negnames1sub <- sub('_DN$', '', Negnames1sub)
  
  
  
  
  
  Negnames1sub <- paste(Negnames1sub, collapse = " ")
  
  
  Negnames2 <- strsplit(Negnames1sub, split = '[_ ]+')
  NegVnames <- as.vector(Negnames2[[1]])
  Negnames3 <- NegVnames[!NegVnames %in% VBLnames]
  NegVnames2 <- unique(Negnames3)
  NegVnames2 <- sort(NegVnames2)
  
  Colors = c(rep("red", length(PosVnames2)), rep("blue", length(NegVnames2)))
  
  wordList <- append(PosVnames2, NegVnames2)
  
  wordList <- data.frame(wordList)
  
  wordList <- cbind(wordList, Colors)
  
  wordListFreq <- append(Posnames3, Negnames3)
  wordListFreq <- wordListFreq[order(wordListFreq)]
  
  #wordListFreq <- gsub("(^|-)([[:alpha:]])", "\\1\\U\\2", wordListFreq, perl=TRUE)
  
  #Construct frequencies matrix
  
  FreqMatPos <- data.frame(ST = Posnames3,
                        
                        row.names = NULL)
  FreqMatPos <- FreqMatPos[order(FreqMatPos$ST), , drop = FALSE]
  
  #FreqMat$ST <- gsub("(^|-)([[:alpha:]])", "\\1\\U\\2", FreqMat$ST, perl=TRUE)
  
  FreqPos <- vector(length = length(unique(Posnames3)))
  j = 1
  
  if(!is.na(nrow(FreqMatPos))) {
    FreqPos[j] <- 1
    
    if(nrow(FreqMatPos) > 1) {
  for (i in 1:(nrow(FreqMatPos) - 1)) {
    check <- FreqMatPos[i, 1]
    
    if (check != FreqMatPos[(i + 1), 1]) {
      j <- j + 1
      FreqPos[j] <- 1
      
    } else {
      FreqPos[j] <- FreqPos[j] + 1
      
    }
    
  }
    }
  }
  
  freqMatPos <-
    noquote(data.frame(
      ST = unique(FreqMatPos$ST),
      freq = FreqPos,
      row.names = NULL
    ))
  
  
  
  FreqMatNeg <- data.frame(ST = Negnames3,
                           
                           row.names = NULL)
  FreqMatNeg <- FreqMatNeg[order(FreqMatNeg$ST), , drop = FALSE]
  
  #FreqMat$ST <- gsub("(^|-)([[:alpha:]])", "\\1\\U\\2", FreqMat$ST, perl=TRUE)
  
  FreqNeg <- vector(length = length(unique(Negnames3)))
  j = 1
  
  if(!is.na(nrow(FreqMatNeg))) {
  FreqNeg[j] <- 1
  
  if(nrow(FreqMatNeg) > 1) {
    
  for (i in 1:(nrow(FreqMatNeg) - 1)) {
    check <- FreqMatNeg[i, 1]
    
    if (check != FreqMatNeg[(i + 1), 1]) {
      j <- j + 1
      FreqNeg[j] <- 1
      
    } else {
      FreqNeg[j] <- FreqNeg[j] + 1
      
    }
  }
  }
  }
  
  freqMatNeg <-
    noquote(data.frame(
      ST = unique(FreqMatNeg$ST),
      freq = FreqNeg,
      row.names = NULL
    ))
  
  freqMat <- rbind(freqMatPos, freqMatNeg)
  
  
  
  
  #commented out bc some thuings need capitilization in wordcloud (TGFB1)
  #wordList$wordList <- tolower(wordList$wordList)
  
  
  #Final Modifications to Lists
  
  #wordList <- wordList[order(wordList$wordList),]
  
  #wordList$wordList <- gsub("(^|-)([[:alpha:]])", "\\1\\U\\2", wordList$wordList, perl=TRUE)
  
  #wordList <- wordList[!duplicated(wordList$wordList), ]
  
  wordList$wordList <- as.character(wordList$wordList)
  
  wordList$Colors <- as.character(wordList$Colors)
  
  wordList <- wordList[nchar(wordList$wordList) > 2, ]
  
  
  
  
  
  #commented out bc some thuings need capitilization in wordcloud (TGFB1)
  #freqMat$ST <- tolower(freqMat$ST)
  
  #freqMat <- freqMat[!duplicated(freqMat$ST), ]
  
  freqMat$ST <- as.character(freqMat$ST)
  
  freqMat <- freqMat[nchar(freqMat$ST) > 2, ]
  
  
  
  
  
  #Masterlist[[i]] <- freqMat$ST
  
  
  
  
  #data <- data.frame(ST = wordList$wordList, freq = freqMat$freq)
  #data <- data[data$freq > 1, , drop = FALSE]
  
  png( file = paste0(sub('.csv$', '', file.names[p]), "wordCloud.png"),
       width = 822,
       height = 622)
  
  wordcloud(
    words = wordList$wordList,
    freq = freqMat$freq,
    scale = c(4, .5),
    colors = wordList$Colors,
    ordered.colors = TRUE,
    min.freq = 2,
    random.order = FALSE,
    random.color = FALSE,
    max.words = 50,
    rot.per = 0,
    family = "sans",
    font = 2
  )
  
  dev.off()
  
 # wordcloud2(
    #data = data,
   # minSize = .5,
    #fontFamily = "sans",
    #color = wordList$Colors,
   # backgroundColor = "white",
   # gridSize = 10,
    #shape = "circle",
    #ellipticity = .65
  #)
  

  
  
}

#Masterlistfinal <- unique(unlist(Masterlist[1:11]))
#write.table(Masterlistfinal, "uniqueTags.csv", sep = ",")
