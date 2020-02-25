library("stringr")
library("protr")
library("randomForest")
library("compiler")
library("bio3d")
library("kernlab")

source("Blomap.R")

options(java.parameters="-Xmx80G")

#### DATA PREPARATION ####
data <- read.csv("data/Immunogenicity_Dataset.csv")

peptides <- data.frame(
  "MutatedPeptide" = data$Peptide,
  "Immunogenicity" = data$Immunogenicity
)
peptides$MutatedPeptide <- as.character(peptides$MutatedPeptide)

#######################
#### RANDOM FOREST ####
#######################

#### FACTOR ####
peptidesFactor <- do.call(rbind.data.frame, lapply(peptides$MutatedPeptide, function(x){
  strsplit(str_pad(x, 11, side = "right", pad = 'X'), split="")[[1]]
}))
colnames(peptidesFactor) <- paste("P", 1:11, sep="_")
peptidesFactor$Immunogenicity <- peptides$Immunogenicity

#### HYDROPHOBICITY / BULKINESS / POLARITY ####
pepPhys <- do.call(rbind.data.frame, lapply(peptides$MutatedPeptide, function(x){
  aas <- strsplit(str_pad(x, 11, side = "right", pad = 'X'), split="")[[1]]
  as.vector(rbind(
    "hydrophobicity" = aa2index(aas, index = "KYTJ820101", window=1),
    "bulkiness" = aa2index(aas, index = "ZIMJ680102", window = 1),
    "polarity" = aa2index(aas, index = "GRAR740102", window = 1)
  ))
}))
pepPhys[is.na(pepPhys)] <- 0

colnames(pepPhys) <- as.vector(t(outer(
  c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11"),
  c("Hydro", "Bulk", "Polarity"),
  paste, sep="_")))
pepPhys$Immunogenicity <- peptides$Immunogenicity

#### BLOMAP #####
pepBlomap <- do.call(rbind.data.frame, lapply(peptides$MutatedPeptide, function(x){
  aas <- strsplit(str_pad(x, 11, side = "right", pad = 'X'), split="")[[1]]
  as.vector(sapply(aas, blomap))
}))
colnames(pepBlomap) <- as.vector(t(outer(
  c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11"),
  c("X1", "X2", "X3", "X4", "X5"),
  paste, sep="_")))
pepBlomap$Immunogenicity <- peptides$Immunogenicity

#### AAC ####
aac <- do.call(rbind.data.frame, lapply(peptides$MutatedPeptide, extractAAC))
colnames(aac) <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P",  "S", "T", "W", "Y", "V")
aac$Immunogenicity <- peptides$Immunogenicity

#### QSO ####
qso <- do.call(rbind.data.frame, lapply(peptides$MutatedPeptide, extractQSO, nlag=7))
colnames(qso) <- names(extractQSO(peptides$MutatedPeptide[1], nlag=7))
qso$Immunogenicity <- peptides$Immunogenicity

#### PseAAC ####
paac <- do.call(rbind.data.frame, lapply(peptides$MutatedPeptide, extractPAAC, lambda = 7))
colnames(paac) <- names(extractPAAC(peptides$MutatedPeptide[1], lambda = 7))
paac$Immunogenicity <- peptides$Immunogenicity

### APAAC ###
apaac <- do.call(rbind.data.frame, lapply(peptides$MutatedPeptide, extractAPAAC, lambda = 7))
colnames(apaac) <- names(extractAPAAC(peptides$MutatedPeptide[1], lambda = 7))
apaac$Immunogenicity <- peptides$Immunogenicity

##############
#### KSVM ####
##############

### BLOSUM 62 ###
blosum62 <- as.kernelMatrix(data.matrix(read.table("data/gBlosum62.txt")))

### BLOSUM 35 ####
blosum35 <- as.kernelMatrix(data.matrix(read.table("data/gBlosum35.txt")))

### WEN ###
wen <- as.kernelMatrix(data.matrix(read.table("data/gWen.txt")))

#save.image(file="Bootstrap.Rdata")
#load(file="Bootstrap.Rdata")

###################
#### BOOTSTRAP ####
###################

set.seed(01015742)

do.Bootstrap <- function(){
  nIter <- 100
  set.seed(01015742)
  
  lapply(1:nIter, function(x){
    inBag <- sample(1:nrow(peptides), nrow(peptides), replace=TRUE)
    ooBag <- (1:nrow(peptides))[-inBag]
    
    print(paste("Iteration", x, Sys.time(), sep=" "))
    
    yOoBag <- as.numeric(peptides$Immunogenicity[ooBag])-1
    ooBag.data <- data.frame(
      "peptide" = peptides$MutatedPeptide[ooBag],
      "Immunogenicity" = yOoBag
    )
    
    ### RANDOM FORESTS
    rfF <- randomForest(Immunogenicity ~ ., data=peptidesFactor[inBag,])
    predF <- predict(rfF, peptidesFactor[ooBag,], type="prob")[,2]
    ooBag.data$Factor <- predF
    
    rfPh <- randomForest(Immunogenicity ~ ., data=pepPhys[inBag,])
    predPh <- predict(rfPh, pepPhys[ooBag,], type="prob")[,2]
    ooBag.data$Hydro <- predPh
    
    rfBm <- randomForest(Immunogenicity ~ ., data=pepBlomap[inBag,])
    predBm <- predict(rfBm, pepBlomap[ooBag,], type="prob")[,2]
    ooBag.data$Blomap <- predBm
    
    rfAAC <- randomForest(Immunogenicity ~ ., data=aac[inBag,])
    predAAC <- predict(rfAAC, aac[ooBag,], type="prob")[,2]
    ooBag.data$AAC <- predAAC
    
    rfQSO <- randomForest(Immunogenicity ~ ., data=qso[inBag,])
    predQSO <- predict(rfQSO, qso[ooBag,], type="prob")[,2]
    ooBag.data$QSO <- predQSO
    
    rfPAAC <- randomForest(Immunogenicity ~ ., data=paac[inBag,])
    predPAAC <- predict(rfPAAC, paac[ooBag,], type="prob")[,2]
    ooBag.data$PAAC <- predPAAC
    
    rfAPAAC <- randomForest(Immunogenicity ~ ., data=apaac[inBag,])
    predAPAAC <- predict(rfAPAAC, apaac[ooBag,], type="prob")[,2]
    ooBag.data$APAAC <- predAPAAC
    
    
    ### KSVM ###
    model62 <- ksvm(blosum62[inBag, inBag], peptides$Immunogenicity[inBag], type = "C-svc", kernel = "matrix", C=1, cache=40000)
    pred62 <- predict(model62, as.kernelMatrix((blosum62[ooBag, inBag])[,SVindex(model62)]), type="decision")
    ooBag.data$Blosum62 <- pred62
    
    model35 <- ksvm(blosum35[inBag, inBag], peptides$Immunogenicity[inBag], type = "C-svc", kernel = "matrix", C=1, cache=40000)
    pred35 <- predict(model35, as.kernelMatrix((blosum35[ooBag, inBag])[,SVindex(model35)]), type="decision")
    ooBag.data$Blosum35 <- pred35
    
    modelWen <- ksvm(wen[inBag, inBag], peptides$Immunogenicity[inBag], type = "C-svc", kernel = "matrix", C=1, cache=40000)
    predWen <- predict(modelWen, as.kernelMatrix((wen[ooBag, inBag])[,SVindex(modelWen)]), type="decision")
    ooBag.data$Wen <- predWen
    
    tableName <- paste(paste("iteration", x, sep="_"), ".txt", sep="")
    write.table(ooBag.data, file=tableName, row.names = FALSE)
    
  })
}

bootstrap <- cmpfun(do.Bootstrap)
bootstrap()
