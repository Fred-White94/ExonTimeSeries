#Fred White 01/10/2020 
#Notes on using exon level data to predict presence or absence of intronic SVA

library(GenomicRanges)
library(homologene)
library(ggplot2)
library(scales)
library(plyr)
library(gridExtra)



brainEx <- read.table("brainEx.txt", stringsAsFactors = F)

CanExon2 <- read.table("CanExonKG38.txt", stringsAsFactors = F)


SVA_cleanup <- read.table("SVAhg38.txt", stringsAsFactors = F)
Length <- abs(SVA_cleanup[,4] - SVA_cleanup[,5])
SVA_cleanup <- cbind(SVA_cleanup,Length)

CanExon2 <- CanExon2[-which(CanExon2$inout == 0),]

CanExon2 <- CanExon2[,c(1,4,5,7,9,10)]
SVA_cleanup <- SVA_cleanup[,c(1,4,5,7,10,15)]

colnames(CanExon2) <- c("chr","start","end","strand","Name","Exon_ID")
colnames(SVA_cleanup) <- c("chr","start","end","strand","Name","Length")


brainExCan <- brainEx[which(brainEx$Name %in% CanExon2$Exon_ID),]


a <- rle(CanExon2$Name)
i <- NULL
exons<- c()
for(i in 1:length(a$lengths)){
  
  exons <- c(exons,1:a$lengths[i])
  
  
  
}


CanExon2$Exon_ID <- paste(gsub("_.*","",CanExon2$Exon_ID),exons, sep = "_")




a <- rle(brainExCan$Description)
i <- NULL
exons<- c()
for(i in 1:length(a$lengths)){
  
  exons <- c(exons,1:a$lengths[i])
  
  
  
}


brainExCan$Name <- paste(gsub("_.*","",brainExCan$Name),exons, sep = "_")








introns <- function(genetable){
  
  uniquetrs <- unique(genetable[,5])
  
  i <- NULL
  j <- 1
  genelist <- list()
  intronlist <- list()
  for(i in uniquetrs){
    
    gene <- genetable[genetable[,5]==i,]
    genelist[[j]] <- gene
    
    
    k <- NULL
    intron <- data.frame()
    
    for(k in 1:(dim(gene)[1]-1)){
      
      if(dim(gene)[1] > 1){
        
        intron[k,1] <- gene[1,1]
        intron[k,2] <- gene[k,3]
        intron[k,3] <- gene[k+1,2]
        intron[k,4] <- gene[k,4]
        intron[k,5] <- gene[1,5]
        intron[k,6] <- gene[k,6]
        #intron[k,4] <- gene[1,10]
        
        intronlist[[j]] <- intron
        
      }
      
      else{next}
      
    }
    j <- j + 1
  }
  return(intronlist)
}

Introns <- introns(CanExon2)


Introns2 <- do.call(rbind,Introns)

colnames(Introns2) <- c("chr","start","end","strand","Name","Exon_ID")

Introns3 <- Introns2

Introns3[which(Introns3$start > Introns3$end),2:3] <- Introns3[which(Introns3$start > Introns3$end),3:2]

Introns3 <- cbind(Introns3, gsub(".*_","",Introns3$Exon_ID))

HT <- makeGRangesFromDataFrame(Introns3, keep.extra.columns=T)
HSVA <- makeGRangesFromDataFrame(SVA_cleanup, keep.extra.columns = T)

names(HT) <- paste(HT$Name,HT$Exon_ID, sep = "_")

b <- findOverlaps(HT, HSVA, ignore.strand = T)

b <- as.data.frame(b)

#names of genes with SVAs and the length of the SVAs 

hHits <- cbind.data.frame(names(HT[b$queryHits,]), HSVA$Length[b$subjectHits], stringsAsFactors = F)


split <- do.call(rbind,strsplit(hHits[,1], split = "_"))

hHits <- cbind(split[,1],hHits, split[,3])

colnames(hHits) <- c("Gene","ID","SVA_length","Exon_Num")






hHits.agg <- aggregate(SVA_length ~ ID, data = hHits, FUN = sum)


hHits.agg <- cbind(gsub("_.*","",hHits.agg$ID), hHits.agg, gsub(".*_","",hHits.agg$ID))
colnames(hHits.agg) <- c("Gene","ID","SVA_length","Intron")
exonnums <- as.data.frame(table(as.character(CanExon2$Name)))
colnames(exonnums) <- c("Gene","Exon_Num")
hHits.agg <- merge(hHits.agg,exonnums,by = "Gene")


hHits.agg[,4] <- as.numeric(as.character(hHits.agg[,4]))

hHits.aggALL <- hHits.agg

exonnumsNOSVA <- exonnums[-which(exonnums$Gene %in% hHits.aggALL$Gene),]
#exonnumsNOSVA <- exonnumsNOSVA[-which(exonnumsNOSVA$Gene %in% CandList$Gene),]

More <- as.character(hHits.agg$Gene[duplicated(hHits.agg$Gene)])

hHits.agg <- hHits.agg[-which(hHits.agg$Gene %in% More),]

SVAG3UD <- hHits.agg[which(hHits.agg$Intron >= 3 & hHits.agg$Intron <= (hHits.agg$Exon_Num - 3)),]

SVAG3UD <- SVAG3UD[-which(SVAG3UD$SVA_length < 500),]

#######


#TPM Calc
length <- as.data.frame(cbind(as.character(CanExon2$Exon_ID),abs(CanExon2$start - CanExon2$end)))
colnames(length)[1] <- "Name"
length <- join(brainExCan[,c(1,2)],length)
length <- as.numeric(as.character(length[,3]))
####
RPK <- brainExCan[,-c(1,2)]/length
RPK[is.na(RPK)] <- 0


TPM <- t( t(RPK) * 1e6 / colSums(RPK) )
TPM <- cbind(brainExCan[,c(1,2)],TPM)




tissue <- read.table("GTExTissueDat.txt", header = T, stringsAsFactors = F)
#brainTis <- tissue[which(tissue$Comment.histological.type. == "Brain"),]
#BT <- brainTis[,c(1,6,7)]
#BT <- BT[!duplicated(BT$Source.Name),]


##Filter for cell type
Cortex <- tissue[which(tissue$Comment.original.body.site.annotation. == "Brain - Cortex"),c(1,3,6,7)]
Cortex <- Cortex[!duplicated(Cortex$Characteristics.individual.),]
TPM <- cbind(TPM[,c(1,2)], TPM[,which(colnames(TPM) %in% Cortex$Source.Name)])




########
##filter for median exon expression in TPM
i <- NULL
Genes <- unique(TPM[,2])
AE <- data.frame()
for(i in 1:length(Genes)){
  
  temp1 <- TPM[which(TPM[,2] == Genes[i]),]
  temp <- median(as.matrix(temp1[,-c(1,2)]))
  max <- max(as.matrix(temp1[,-c(1,2)]))
  
  AE[i,1] <- temp
  AE[i,2] <- max
  print(i)
}


AE <- cbind(Genes,AE)
AE3AVG <- AE[which(AE[,2] > 3),]

AE30Max <- AE3AVG[which(AE3AVG$V2 <30),]


AE40Max <- AE3AVG[which(AE3AVG$V2 <40),]



##
TPMFilt <- TPM[which(TPM$Description %in% AE40Max$Genes),]




SVAG3UD <- SVAG3UD[which(SVAG3UD$Gene %in% TPMFilt$Description),]

sampleex <- as.data.frame(table(SVAG3UD$Exon_Num))


exonnumsNOSVA <- exonnumsNOSVA[which(exonnumsNOSVA$Gene %in% TPMFilt$Description),]




#######get a list of size matched control genes
#change this to solve problem of trying to sample a larger than possible number using AF NOSVAF

i <- NULL
l <- 0
CS <- list()
temp <- data.frame()
samples <- data.frame()
for(i in 1:dim(sampleex)[1]){
  
  #AF <- dim(exonnums[which(exonnums$Freq == 59),])[1]
  #NOSVAF <- dim(exonnumsNOSVA[which(exonnumsNOSVA$Freq == 59),])[1]
  
  #if(AF/NOSVAF < 2) {
  
  temp <- exonnumsNOSVA[which(exonnumsNOSVA$Exon_Num == sampleex[i,1]),]
  samples <- sample(as.character(temp[,1]), sampleex[i,2])
  
  j <- NULL
  for(j in 1:length(samples)){
    l <- l + 1
    CS[[l]] <- CanExon2[which(CanExon2$Name == samples[j]),]
    
    # CS[[l]] <- KG38d[grep(samples[j],KG38d$Gids),]
    
  }
  
  
  
  #}
  
  #else{next}
  
  
  
}


#SVAG3UD <- SVAG3UD[which(SVAG3UD$Exon_Num < 27),]


SVAG3UD <- SVAG3UD[order(SVAG3UD$Exon_Num),]

i <- NULL
ControlGenes <- vector()
for(i in 1: length(CS)){
  
  ControlGenes[i] <- as.character(CS[[i]][1,5])
  
  
}




CCgenes2 <- c(ControlGenes,as.character(SVAG3UD$Gene))
CCs2 <- cbind(CCgenes2, c(rep("control",(length(CCgenes2))/2), rep("case",(length(CCgenes2))/2)))
##
CCs2 <- cbind(CCs2, rep(SVAG3UD$Intron,2))
CCs2 <- as.data.frame(CCs2)


##

TPMFilt2 <- TPMFilt[which(TPMFilt$Description %in% CCs2$CCgenes2),]


exons <- gsub(".*_","",TPMFilt2[,1])
exons <- as.numeric(exons)
TPMFilt2 <- cbind(TPMFilt2[,1],exons,TPMFilt2[,-1])
colnames(TPMFilt2)[3] <- "GeneId"
colnames(CCs2)[1] <- "GeneId"
TPMFilt2 <- merge(CCs2,TPMFilt2, by = "GeneId")
TPMFilt2 <- TPMFilt2[order(TPMFilt2[,2],TPMFilt2[,1],TPMFilt2[,5]),]
Case <- CCs2[which(CCs2$V2 == "case"),]
Control <- CCs2[which(CCs2$V2 == "control"),]
Case$V3 <- as.numeric(as.character(Case$V3))
Control$V3 <- as.numeric(as.character(Control$V3))
#jump to line 384 for no 3up 3down filter







######
i <- NULL
CUp3 <- list()
CDown3 <- list()
ContUp3 <- list()
ContDown3 <- list()

for(i in 1:dim(Case)[1]){
  
  #caseTemp <- gds[grep(Case[i,1],gds[,1]),]
  
  caseTemp <- TPMFilt2[which(as.character(TPMFilt2[,1]) == as.character(Case[i,1])),]
  controlTemp <- TPMFilt2[which(as.character(TPMFilt2[,1]) == as.character(Control[i,1])),]
 
  
  
  
  CUp3[[i]] <- caseTemp[(Case[i,3]-2):Case[i,3],]
  CDown3[[i]] <- caseTemp[(Case[i,3] + 1):(Case[i,3] + 3),]
  CUp3[[i]]$exons <- c(-3,-2,-1)
  CDown3[[i]]$exons <- c(1,2,3)
  ###
  ContUp3[[i]] <- controlTemp[(Control[i,3]-2):Control[i,3],]
  ContDown3[[i]] <- controlTemp[(Control[i,3] + 1):(Control[i,3] + 3),]
  ContUp3[[i]]$exons <- c(-3,-2,-1)
  ContDown3[[i]]$exons <- c(1,2,3)

}


###########
CUp4 <- do.call(rbind,CUp3)

CDown4 <- do.call(rbind,CDown3)

ContUp4 <- do.call(rbind,ContUp3)

ContDown4 <- do.call(rbind,ContDown3)






ThreeExon <- rbind(CUp4,CDown4,ContUp4,ContDown4)




#####jump to here from line 326 

TEmeltCortex <- melt(ThreeExon, id.vars = c("GeneId","V2","V3","TPMFilt2[, 1]","exons"))
colnames(TEmeltCortex)[5] <- "Exon_Position"
colnames(TEmeltCortex)[7] <- "TPM_Counts"
colnames(TEmeltCortex)[2] <- "Class"



#CandList <- read.table("20_07_23HumanSpecificSVAimprovedRsquared.txt")


library(tidyr)
library(dplyr)




phenotypes <- read.table("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",stringsAsFactors=  F, sep = "\t", header = T)
head(phenotypes)
as.data.frame(table(phenotypes$AGE))
#phenotypes$SUBJID <- gsub("-",".",phenotypes$SUBJID)


TEmeltCortex$variable <- gsub("\\.","-",TEmeltCortex$variable)

IDs <- do.call(rbind,strsplit(TEmeltCortex$variable, split = "-"))[,c(1,2)]
IDs <- paste0(IDs[,1],"-",IDs[,2])
TEmeltCortex <- cbind(IDs,TEmeltCortex)
colnames(TEmeltCortex)[1] <- "SUBJID"

TEmeltCortex <- merge(phenotypes, TEmeltCortex, by = "SUBJID")


IDs <- paste(TEmeltCortex$GeneId,TEmeltCortex$SUBJID, sep = "_") 
TEmeltCortex <- cbind(TEmeltCortex,IDs)
TE4RF <- TEmeltCortex[,c(12,5,6,9,11)]
Pdat <- pivot_wider(TE4RF, names_from = Exon_Position, values_from = TPM_Counts)

#######


Case <- Pdat[which(Pdat$Class == "case"),]
Control <- Pdat[which(Pdat$Class == "control"),]


set.seed(124)
CaseTrain <- sample(unique(Case$GeneId),30)
ControlTrain <- sample(unique(Control$GeneId),30)

TrainGenes <- c(as.character(CaseTrain),as.character(ControlTrain))


train <- Pdat[which(Pdat$GeneId %in% TrainGenes),]
test  <- Pdat[-which(Pdat$GeneId %in% TrainGenes),]

##change this scaling 
train[,-c(1:3)] <- t(apply(train[,-c(1:3)],2,FUN = rescale))
test[,-c(1:3)] <- t(apply(test[,-c(1:3)],2,FUN = rescale))

trainclass <- as.character(train$Class)


testclass <- as.character(test$Class)

library(caret)
library(matrixStats)
library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

cat("Training Random Forest")

fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10)

rffit <- caret::train(x = as.matrix(train[,-c(1:3)]),y = trainclass, method = "rf", trControl = fitControl)

#svmfit <- caret::train(x = as.matrix(train[,-1]),y = trainclass, method = "svmLinearWeights2", trControl = fitControl)

######
CVACC <- rffit$results[which(rffit$results$Accuracy == max(rffit$results$Accuracy)),][1,]

cat("Predicting Test Data")

PREDICT <- predict(rffit, newdata = test[,-c(1:3)])
TESTACC <- length(which(PREDICT == testclass))/length(PREDICT)



stopCluster(cl)


