#Fred White 01/10/2020
#Script to label exon number and also compare overlapping ranges
#For better optimised overlapping ranges see file OverlappingRanges

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(qlcMatrix)
library(caret)



##File used by GTEx for exon refrence
KG38 <- read.csv("~/RNA_Seq/gencode.v26.GRCh38.genes.gtf.csv", header = F, stringsAsFactors = F)
KG38a <- KG38[which(KG38$V3 == "exon"),]
Gids <- gsub(".*transcript_name ","",KG38a[,9])
Gids <- gsub("\\;.*","",Gids)
Gids <- gsub("\\\"","",Gids)
Eids <- gsub(".*exon_id ","",KG38a[,9])
Eids <- gsub("\\;.*","",Eids)
Eids <- gsub("\\\"","",Eids)
KG38a <- KG38a[,-9]
KG38a <- cbind(KG38a,Gids,Eids)


kgXref <- read.table("kgXref2.txt",stringsAsFactors=F,header =T)


#G4EXON2 <- read.table("G4EXON2.txt",stringsAsFactors=F)
######TAF1 not in gen26......
#"TAF1" %in% kgXref$geneSymbol
#"TAF1" %in% gencan$geneSymbol
#"TAF1" %in% gen26$V13


gen26 <- read.table("gencode26CDS",stringsAsFactors=F)
#remove non coding transcripts
gen26 <- gen26[-which(gen26$V14 == "none"),]
#gen26SAFE <- gen26
gen26 <- gen26[-which(gen26$V14 == "incmpl"),]
gen26 <- gen26[-which(gen26$V15 == "incmpl"),]
colnames(gen26)[2] <- "kgID"


gencan <- read.table("gencode32canonical.txt",stringsAsFactors=F)
colnames(gencan)[5] <- "kgID"
gencan <- merge(kgXref,gencan, by = "kgID")
a <- c(1:22,"X","Y")
b <- paste0("chr",a)
gencan <- gencan[which(gencan$V1 %in% b),]




###change from here for machine learning of cell type specific isoform predicitons
gen26 <- gen26[which(gsub("\\..*","",gen26$kgID) %in% gsub("\\..*","",gencan$kgID)),]
#dirty way to remove duplicate genes - could improve here
gen26 <- gen26[-which(duplicated(gen26$V13)),]
gen26 <- gen26[-which(gen26$V9 < 2),]

KG38a <- KG38a[which(KG38a$Gids %in% gen26$V13),]
gen26 <- gen26[which(gen26$V13 %in% KG38a$Gids),]


#########filter annotation file based on overlap with canonical form file

i <- NULL
CanExon <- list()
for(i in 1:dim(gen26)[1]){
  
  canon <- cbind(as.numeric(unlist(strsplit(gen26[i,10],","))),
                as.numeric(unlist(strsplit(gen26[i,11],","))))
  
  
  if(gen26$V4[i] == "-"){
    rownames(canon) <- c(1:dim(canon)[1])
    canon <- canon[rev(rownames(canon)),]
    }

#####  
  query <- KG38a[which(KG38a$Gids == gen26$V13[i]),]
  res <- matrix(nrow = dim(query)[1],ncol = dim(canon))
  j <- NULL
  canexon <- matrix(nrow = dim(query)[1])
  for(j in 1:dim(canon)[1]){
    k <- NULL
    for(k in 1:dim(query)[1]){
      
      if(sum(canon[j,1]:canon[j,2] %in% query[k,4]:query[k,5])/length(canon[j,1]:canon[j,2]) > 0.8){
        canexon[k,1] <- 1
      }
      
      else{canexon[k,1] <- 0}
 
    }
    res[,j] <- canexon

  }
  
  inout <- rowSums(res)
  inout[which(inout > 1)] <- 1
  query <- cbind(query,inout)
  CanExon[[i]] <- query
  print(paste0("i = ",i," length of CanExon =",length(CanExon)))
}


##very slow:
CanExon2 <- do.call(rbind, CanExon)

write.table(CanExon2, "CanExonKG38.txt")

SVA_cleanup <- read.table("~/RNA_Seq/SVA_CleanupHG38.bed.bed", header = F, stringsAsFactors=F)
SVA_length <- abs(SVA_cleanup$V2 - SVA_cleanup$V3)
SVA_cleanup <- cbind(SVA_cleanup,SVA_length)


CanExon2 <- CanExon2[-which(CanExon2$inout == 0),]
#
CanExon3 <- CanExon2[which(CanExon2$Gids %in% G4A$geneSymbol),]#################################
roilist3 <- introns(CanExon3)
roilist3 <- compact(roilist3)

#this command takes a long time if used on all transcripts here it has already been subset to a list genereated in a previous analysis to save days
SVAG <- inROI(roilist3,SVA_cleanup)

exonnums <- as.data.frame(table(as.character(CanExon2$Gids)))

colnames(exonnums)[1] <- "Gene"
colnames(SVAG)[1] <- "Gene"

SVAG <- merge(SVAG,exonnums,by = "Gene")
exonnumsNOSVA <- exonnums[-which(exonnums$Var1 %in% SVAG$V1),]

SVAG3UD <- SVAG[which(SVAG$V6 >= 3 & SVAG$V6 <= (SVAG$Freq - 3)),]



#######



brainExCan <- brainEx[which(brainEx$Name %in% CanExon2$Eids),]



#TPM Calc
length <- as.data.frame(cbind(as.character(CanExon2$Eids),abs(CanExon2$V4 - CanExon2$V5)))
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



Cortex <- tissue[which(tissue$Comment.original.body.site.annotation. == "Brain - Putamen (basal ganglia)"),c(1,3,6,7)]
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


TPMFilt <- TPM[which(TPM$Description %in% AE30Max$Genes),]




manysvas <- SVAG3UD$Gene[duplicated(SVAG3UD$Gene)]
SVAGF <- SVAG3UD[-which(SVAG3UD$Gene %in% manysvas),]

#SVAGF <- SVAGF[which(SVAGF$Gene %in% AE3AVG$Genes),]

SVAGF <- SVAGF[which(SVAGF$Gene %in% AE30Max$Genes),]
#dirtily remove some genes of interest because no size matched control
SVAGFsafe2 <- SVAGF
SVAGF <- SVAGF[order(SVAGF$Freq),]
#SVAGF <- SVAGF[-which(SVAGF$Freq == 40),]
#SVAGF <- SVAGF[-which(SVAGF$Freq == 55),]
#SVAGF <- SVAGF[-which(SVAGF$Freq == 60),]
#SVAGF <- SVAGF[-which(SVAGF$Freq == 95),]


#####
sampleex <- as.data.frame(table(SVAGF$Freq))

exonnumsNOSVA3AVG <- exonnumsNOSVA[which(exonnumsNOSVA$Var1 %in% AE30Max$Genes),]

exonnumsNOSVA3AVG2 <- exonnumsNOSVA3AVG[-which(exonnumsNOSVA3AVG$Var1 %in% ControlGenes),]


#allControl
#exonnumsNOSVA3AVG2 <- exonnumsNOSVA3AVG[-which(exonnumsNOSVA3AVG$Var1 %in% allControl),]


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
  
  temp <- exonnumsNOSVA3AVG2[which(exonnumsNOSVA3AVG2$Freq == sampleex[i,1]),]
  samples <- sample(as.character(temp[,1]), sampleex[i,2])
  
  j <- NULL
  for(j in 1:length(samples)){
    l <- l + 1
    CS[[l]] <- CanExon2[which(CanExon2$Gids == samples[j]),]
    
    # CS[[l]] <- KG38d[grep(samples[j],KG38d$Gids),]
    
  }
  
  
  
  #}
  
  #else{next}
  
  
  
}



#SVAGF <- SVAGF[order(SVAGF$Freq),]


#####ordering wrong....


i <- NULL
ControlGenes <- vector()
for(i in 1: length(CS)){
  
  ControlGenes[i] <- as.character(CS[[i]][1,9])
  
  
}




CCgenes2 <- c(ControlGenes,as.character(SVAGF$Gene))
CCs2 <- cbind(CCgenes2, c(rep("control",(length(CCgenes2))/2), rep("case",(length(CCgenes2))/2)))
##
CCs2 <- cbind(CCs2, rep(SVAGF$V6,2))
CCs2 <- as.data.frame(CCs2)


TPMFilt2 <- TPMFilt[which(TPMFilt$Description %in% CCs2$CCgenes2),]


TPMFilt2 <- TPM[which(TPM$Description %in% CCs2$GeneId),]


a <- rle(TPMFilt2$Description)
a <- a$lengths
b <- rep(1,length(a))
c <- cbind(b,a)

i <- NULL
exons <- c()
for(i in 1:dim(c)[1]){
  
  
  exons <- c(exons,c(c[i,1]:c[i,2]))
  
}

#############################


TPMFilt2 <- cbind(TPMFilt2[,1],exons,TPMFilt2[,-1])


colnames(TPMFilt2)[3] <- "GeneId"
colnames(CCs2)[1] <- "GeneId"
TPMFilt2 <- merge(CCs2,TPMFilt2, by = "GeneId")
TPMFilt2 <- TPMFilt2[order(TPMFilt2[,2],TPMFilt2[,1],TPMFilt2[,5]),]
Case <- CCs2[which(CCs2$V2 == "case"),]
Control <- CCs2[which(CCs2$V2 == "control"),]
Case$V3 <- as.numeric(as.character(Case$V3))
Control$V3 <- as.numeric(as.character(Control$V3))

#diagnosing route of NAs in Control up/down3
#a <- 1:96
#a[-which(Control$GeneId %in% ContUp4$GeneId)]
#a[-which(Control$GeneId %in% ContDown4$GeneId)]



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
  
  #controlTemp <- gds[grep(Control[i,1],gds[,1]),]
  
  
  
  
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

TEmeltCortex <- melt(ThreeExon, id.vars = c("GeneId","V2","V3","TPMFilt2[, 1]","exons"))
colnames(TEmeltCortex)[5] <- "Exon_Position"
colnames(TEmeltCortex)[7] <- "TPM_Counts"
colnames(TEmeltCortex)[2] <- "Intronic_SVA"


########

#dim(cast(TEmeltCortex, Exon_Position~GeneId))

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
TE4RF <- TEmeltCortex[,c(12,6,9,11)]
Pdat <- pivot_wider(TE4RF, names_from = Exon_Position, values_from = TPM_Counts)




#dim(cast(TEmeltCortex, TPM_Counts~IDs))

#Class <- c(rep("case",dim(CaseR)[2]),c(rep("control",dim(ControlR)[2])))



Class <- Pdat$Intronic_SVA
Pdat <- Pdat[,-2]
trainIndex <- createDataPartition(Class, p = .75, list = FALSE, times = 1)
train <- Pdat[trainIndex,]
test  <- Pdat[-trainIndex,]

trainclass <- Class[trainIndex]
testclass <- Class[-trainIndex]
######
library(matrixStats)
library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

cat("Training Random Forest")

fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
rffit <- caret::train(x = as.matrix(train[,-1]),y = trainclass, method = "rf", trControl = fitControl)

#svmfit <- caret::train(x = as.matrix(train[,-1]),y = trainclass, method = "svmLinearWeights2", trControl = fitControl)

######
CVACC <- rffit$results[which(rffit$results$Accuracy == max(rffit$results$Accuracy)),][1,]

cat("Predicting Test Data")

PREDICT <- predict(rffit, newdata = test[,-1])
TESTACC <- length(which(PREDICT == testclass))/length(PREDICT)



stopCluster(cl)


trainVotes <- cbind(train,trainclass,rffit$finalModel$votes)
trainVotes$Genes <- gsub("_.*","",trainVotes$IDs)

#pdf("02_04_20meanTPMbyvotes.pdf")
#plot(trainVotes$case,rowMeans(trainVotes[,c(2:7)]), xlab = "Vote %", ylab  = "Mean Expression")
#dev.off()




pdf("02_04_20TPMExonVotesPut.pdf")
par(mfrow=c(2,3)) 
plot(trainVotes$case*100,trainVotes$'-3', xlab = "Vote %", ylab  = "Exon -3 Expression", ylim = c(0,30))
plot(trainVotes$case*100,trainVotes$'-2', xlab = "Vote %", ylab  = "Exon -2 Expression", ylim = c(0,30))
plot(trainVotes$case*100,trainVotes$'-1', xlab = "Vote %", ylab  = "Exon -1 Expression", ylim = c(0,30))
plot(trainVotes$case*100,trainVotes$'1', xlab = "Vote %", ylab  = "Exon 1 Expression", ylim = c(0,30))
plot(trainVotes$case*100,trainVotes$'2', xlab = "Vote %", ylab  = "Exon 2 Expression", ylim = c(0,30))
plot(trainVotes$case*100,trainVotes$'3', xlab = "Vote %", ylab  = "Exon 3 Expression", ylim = c(0,30))
dev.off()





#VARS <- colVars(RFdat)
#RFdatVF <- RFdat[,which(VARS > quantile(VARS,  .90))]
#train <- RFdatVF[trainIndex,]
#test  <- RFdatVF[-trainIndex,]
#trainclass <- Class[trainIndex]
#testclass <- Class[-trainIndex]

#PCD <- prcomp(train, scale. = T)
#test <- predict(PCD,test)[,1:10]
#train <- PCD$x[,1:5]

#design <- cbind(rownames(train),trainclass)
#pdf("PCARF.pdf")
#autoplot(PCD, data = design, colour = 'trainclass')
#dev.off()

#ICD <- icafast(train, nc = 35, center = T)
#test <- tcrossprod(as.matrix(test),ICD$W)
#train <- ICD$S
#colnames(train) <- paste("IC",c(1:35))
#colnames(test) <- paste("IC",c(1:35))











######graphs

pdf("20_04_02_avgexoncountsboxplotPut.pdf")
ggplot(TEmeltCortex, aes(x= as.factor(Exon_Position), y = TPM_Counts, fill = Intronic_SVA)) +
  geom_boxplot( outlier.shape = NA)  +  #position = position_dodge(1),
  #ylim(NA, 50)  +
  xlab("Exon Position") +
  ylab("TPM Counts")#
dev.off()






library(plyr)
mu <- ddply(TEmeltCortex, "Intronic_SVA", summarise, grp.mean=mean(TPM_Counts))
head(mu)

pdf("02_04_20ExonTPMDistributionPut.pdf")
ggplot(TEmeltCortex, aes(x = TPM_Counts, colour = Intronic_SVA)) +
  geom_histogram(alpha=0.5, position="identity") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Intronic_SVA),linetype="dashed")
dev.off()


#pdf("20_03_20_avgexoncountsboxplotAGEoutliersCONTROL.pdf")
#ggplot(TEmeltCortex[which(TEmeltCortex$Intronic_SVA == "control"),], aes(x= as.factor(Exon_Position), y = TPM_Counts, fill = AGE)) +
#  geom_boxplot( outlier.shape = NA)  +  #position = position_dodge(1),
#  ylim(NA, 50)  +
#  xlab("Exon Position") +
#  ylab("TPM Counts") +
#  ggtitle("Non-Intronic SVA Genes") #
#dev.off()



###########################################
#PARALLEL DOES NOT WORK
#roilist <- introns(CanExon2)
#some NAs produced due to no introns use function below from plyr to get rid of these
#roilist <- compact(roilist)

#list <- inROIPAR(roilist,SVA_cleanup)



######
i <- NULL
C1 <- list()

Cont1 <- list()


for(i in 1:dim(Case)[1]){
  
  #caseTemp <- gds[grep(Case[i,1],gds[,1]),]
  
  caseTemp <- TPMFilt2[which(as.character(TPMFilt2[,1]) == as.character(Case[i,1])),]
  controlTemp <- TPMFilt2[which(as.character(TPMFilt2[,1]) == as.character(Control[i,1])),]
  

  C1[[i]] <- caseTemp

  ###
  Cont1[[i]] <- controlTemp

  
}


###########
C2 <- do.call(rbind,C1)
Cont2 <- do.call(rbind,Cont1)


DF <- rbind(C2,Cont2)

TEmeltCortex2 <- melt(DF, id.vars = c("GeneId","V2","V3","TPMFilt2[, 1]","exons"))
colnames(TEmeltCortex2)[5] <- "Exon_Position"
colnames(TEmeltCortex2)[7] <- "TPM_Counts"
colnames(TEmeltCortex2)[2] <- "Intronic_SVA"

library(plyr)
mu <- ddply(TEmeltCortex2, "Intronic_SVA", summarise, grp.mean=mean(TPM_Counts))
head(mu)

pdf("30_03_20ExonTPMDistributionWG.pdf")
ggplot(TEmeltCortex2, aes(x = TPM_Counts, colour = Intronic_SVA)) +
  geom_histogram(alpha=0.5, position="identity") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Intronic_SVA),linetype="dashed")
dev.off()








#####################################################
set.seed(123)
#Class <- Pdat$Intronic_SVA
#Pdat <- Pdat[,-2]#
SPdat <- cbind(Pdat[,1],scale(Pdat[,-1]))
trainIndex <- createDataPartition(Class, p = .75, list = FALSE, times = 1)
train <- SPdat[trainIndex,]
test  <- SPdat[-trainIndex,]

trainclass <- Class[trainIndex]
testclass <- Class[-trainIndex]
######
library(matrixStats)
library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

cat("Training Random Forest")

fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
#rffit2 <- caret::train(x = as.matrix(train[,-1]),y = trainclass, method = "rf", trControl = fitControl)

svmfit <- caret::train(x = as.matrix(train[,-1]),y = trainclass, method = "svmRadial", trControl = fitControl)

######
CVACC <- rffit$results[which(rffit$results$Accuracy == max(rffit$results$Accuracy)),][1,]

cat("Predicting Test Data")

PREDICT <- predict(rffit, newdata = test[,-1])
TESTACC <- length(which(PREDICT == testclass))/length(PREDICT)




 