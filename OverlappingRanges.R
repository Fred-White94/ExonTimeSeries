#Fred White 01/10/2020

library(GenomicRanges)
library(homologene)
library(ggplot2)
library(scales)

#do this only once per session
#homologeneData <- updateHomologene()


hg38rep <- read.table("Repeatshg38.txt",stringsAsFactors=F, sep = "\t")

panTrorep <- read.table("RepeatsPanTro6.txt", stringsAsFactors=F, sep = "\t")



hg38rep$V9 <- gsub(".* transcript_id ","",hg38rep$V9)
panTrorep$V9 <- gsub(".* transcript_id ","",panTrorep$V9)

hg38rep$V9 <- gsub("_dup.*","",hg38rep$V9)
panTrorep$V9 <- gsub("_dup.*","",panTrorep$V9)


hg38rep$V9 <- gsub("; ","",hg38rep$V9)
panTrorep$V9 <- gsub("; ","",panTrorep$V9)

##ind <- c(grep("L1.*",hg38rep$V9), grep("Alu.*",hg38rep$V9))
#cind <- c(grep("L1.*",panTrorep$V9), grep("Alu.*",panTrorep$V9))

ind <- grep("SVA.*",hg38rep$V9)
cind <- grep("SVA.*",panTrorep$V9)

hg38rep <- hg38rep[ind,]
panTrorep <- panTrorep[cind,]



#hg38rep <- hg38rep[-grep("HAL.*",hg38rep$V9),]
#hg38rep <- hg38rep[-grep("HERV.*",hg38rep$V9),]

#panTrorep <- panTrorep[-grep("HAL.*",panTrorep$V9),]
#panTrorep <- panTrorep[-grep("HERV.*", panTrorep$V9),]


hg38rep <- cbind(hg38rep, abs(hg38rep$V4 - hg38rep$V5))
panTrorep <- cbind(panTrorep, abs(panTrorep$V4 - panTrorep$V5))

#h <- unique(hg38rep$V9)
#c <- unique(panTrorep$V9)


colnames(hg38rep) <- c("chr","Genome","Transcript","start","end","whoknows","strand","nothing","Name","Length")
colnames(panTrorep) <- c("chr","Genome","Transcript","start","end","whoknows","strand","nothing","Name","Length")


hg38rep <- hg38rep[,-c(2,3,6,8)]
panTrorep <- panTrorep[,-c(2,3,6,8)]





#######
HG <- read.table("hg38.ncbiRefSeq(1).gtf", stringsAsFactors=F, sep = "\t")
#HG <- read.table("/zfs/sils/jacobs/SRP121791/sra/hg38.knownGene.gtf", stringsAsFactors=F, sep = "\t")
#kgXref <- read.table("kgXref2.txt",stringsAsFactors=F,header =T)

HT <- HG[which(HG$V3 == "transcript"),]

HT$V9 <- gsub("gene_id ","",HT$V9)
HT$V9 <- gsub(";.*","",HT$V9)
#HT$V9 <- gsub(" transcript.*","",HT$V9)
#HT$V9 <- gsub(".*_name ","",HT$V9)
#HT$V9 <- gsub(";","",HT$V9)

l <- abs(HT[,4] - HT[,5])
HT <- cbind(HT,l)

#colnames(HT)[9] <- "kgID"

#HT <- merge(HT, kgXref, by = "kgID")

#order on decreasing transcript length and then remove shorter duplicates
HT <- HT[order(HT$l, decreasing = T),]
HT <- HT[!duplicated(HT$V9),]

#HT <- HT[-grep("LOC\\d",HT$geneSymbol),]
HT <- HT[-grep("MIR\\d",HT$V9),]
HT <- HT[-grep("MIRLET\\d",HT$V9),]

HT <- HT[-which(HT$l < 500),]


ChG <- read.table("/zfs/sils/jacobs/SRP121791/sra/panTro6.ncbiRefSeq.gtf", stringsAsFactors=F, sep = "\t")
ChT <- ChG[which(ChG$V3 == "transcript"),]
ChT$V9 <- gsub(".*_name ","",ChT$V9)
ChT$V9 <- gsub(";","",ChT$V9)

l <- abs(ChT[,4] - ChT[,5])
ChT <- cbind(ChT,l)

#order on decreasing transcript length and then remove shorter duplicates
ChT <- ChT[order(ChT$l, decreasing = T),]
ChT <- ChT[!duplicated(ChT$V9),]



HT <- HT[,-c(2,3,6,8,10)]
ChT <- ChT[,-c(2,3,6,8,10)]

colnames(HT) <- c("chr","start","end","strand","Name")
colnames(ChT) <- colnames(HT)




####GOanalysis line 200


GOt <- read.table("GOIDtranslate.txt", sep = "\t", stringsAsFactors=F)
colnames(GOt)[1] <- "geneSymbol"
#INTSVAS <- read.table("IntSVAs.txt")


DESeq <- read.table("FieldDESeq.txt", sep = "\t", header = T, stringsAsFactors=F)
colnames(DESeq)[1] <- "geneSymbol"
DESeq <- DESeq[which(DESeq$geneSymbol %in% GOt$geneSymbol),]

colnames(DESeq) <- gsub(".*_","",colnames(DESeq))


##

Human <- DESeq[,c(1:2,3:8)]
Chimp <- DESeq[,c(1:2,9:14)]
Orang <- DESeq[,c(1:2,15:20)]
Rhesus <- DESeq[,c(1:2,21:26)]

list <-list(Human,Chimp,Orang,Rhesus)

i <- NULL
ZEROES <- data.frame(matrix(nrow = lapply(list,dim)[[1]][1],ncol = 4))
X <- data.frame(matrix(nrow = lapply(list,dim)[[1]][1],ncol = 4))
for(i in 1:length(list)){
  
  ZEROES[,i] <- rowSums(data.matrix(list[[i]][,-c(1:2)]) < 1)
  X[,i] <- rowSums(data.matrix(list[[i]][,-c(1:2)]) > 20000)
  
}


i <- NULL
for(i in 1:dim(ZEROES)[2]){
  
  list[[i]] <- cbind(list[[i]],ZEROES[,i],X[,i])
  
  
}


i <- NULL
for(i in 1:4){
  list[[i]] <- list[[i]][-which(list[[i]]$"ZEROES[, i]" > 0),]
  list[[i]] <- list[[i]][-which(list[[i]]$"X[, i]" > 0),]
  
  list[[i]] <- list[[i]][,-c(9:10)]
  
}



Genes <- Reduce(intersect, list(list[[1]]$geneSymbol,list[[2]]$geneSymbol,list[[3]]$geneSymbol,list[[4]]$geneSymbol))


i <- NULL
for(i in 1:4){
  
  list[[i]] <- list[[i]][which(list[[i]]$geneSymbol %in% Genes),]
  
}
#############






#adapt human2mouse function from homologene package
human2chimp <- function (genes, db = homologeneData)
{
  out = homologene(genes, 9606, 9598, db)
  names(out) = c("humanGene", "chimpGene", "humanID", "chimpID")
  return(out)
}




genes <- list[[1]]$geneSymbol
genes <- human2chimp(genes)


HT <- HT[which(HT$Name %in% genes$humanGene),]
ChT <- ChT[which(ChT$Name %in% genes$chimpGene),]





####
HGrange <- makeGRangesFromDataFrame(hg38rep, keep.extra.columns=T)
CGrange <- makeGRangesFromDataFrame(panTrorep, keep.extra.columns = T)


HT <- makeGRangesFromDataFrame(HT, keep.extra.columns = T)
ChT <- makeGRangesFromDataFrame(ChT, keep.extra.columns = T)

names(HT) <- HT$Name
names(ChT) <- ChT$Name



b <- findOverlaps(HT, HGrange, ignore.strand = T)
c <- findOverlaps(ChT, CGrange, ignore.strand = T)


b <- as.data.frame(b)
c <- as.data.frame(c)

#names of genes with SVAs and the length of the SVAs 

hHits <- cbind.data.frame(names(HT[b$queryHits,]), HGrange$Length[b$subjectHits], stringsAsFactors = F)
cHits <- cbind.data.frame(names(ChT[c$queryHits,]), CGrange$Length[c$subjectHits], stringsAsFactors = F)

hHits <- hHits[which(hHits[,1] %in% genes$humanGene),]
cHits <- cHits[which(cHits[,1] %in% genes$chimpGene),]


colnames(hHits) <- c("humanGene", "SVAlength")
colnames(cHits) <- c("chimpGene", "SVAlength")
  

HNPG <- as.data.frame(table(hHits$humanGene))
CNPG <- as.data.frame(table(cHits$chimpGene))

colnames(HNPG) <- c("humanGene", "Hits")
colnames(CNPG) <- c("chimpGene", "Hits")

HNPG <- merge(HNPG, genes, by = "humanGene")
CNPG <- merge(CNPG, genes, by = "chimpGene")


HSpec <- HNPG[-which(HNPG$humanGene %in% CNPG$humanGene),]
CSpec <- CNPG[-which(CNPG$chimpGene %in% HNPG$chimpGene),]
Common <- HNPG[which(HNPG$humanGene %in% CNPG$humanGene),]



i <- NULL
NUMSVA <- data.frame()
l <- list()
for(i in 1:dim(Common)[1]){
  
  nh <- HNPG[which(HNPG$humanGene == Common$humanGene[i]),2]
  nc <- CNPG[which(CNPG$humanGene == Common$humanGene[i]),2]
  
  numdiff <- abs(nh - nc)
  
  if(numdiff == 0){
    
    cl <- cHits[which(cHits[,1] == Common$humanGene[i]),2]
    
    
    hl <- hHits[which(hHits[,1] == Common$humanGene[i]),2]
    
    
    #l[[i]] <- (c(length(cl), length(hl)))
    ldiff <- abs(cl-hl)
    
    if(sum(ldiff <  101) == length(ldiff)){
      
      numsva <- data.frame(Common[i,1], Common[i,2], Common[i,3], stringsAsFactors = F)
      
      NUMSVA <- rbind(NUMSVA,numsva)
      
    }
    
  }
  
  
  
}




hHitsSmol <- hHits[which(hHits$SVAlength < 500),]
cHitsSmol <- cHits[which(cHits$SVAlength < 500),]

NUMSVA <- NUMSVA[-which(NUMSVA[,1] %in% hHitsSmol[,1]),]
NUMSVA <- NUMSVA[-which(NUMSVA[,3] %in% cHitsSmol[,1]),]


CHComm <- NUMSVA
HSpecC <- HSpec
CSpecH <- CSpec


hs <- as.character(unique(HSpec$humanGene))
cs <- unique(CSpec$humanGene)


hs <- hs[-which(hs %in% hHitsSmol[,1])]
cs <- cs[-which(cs %in% cHitsSmol[,1])]
cg <- as.character(NUMSVA[,1])



#####which statement dodgy
genes <- genes[-which(genes$humanGene %in% HNPG$humanGene),]
genes <- genes[-which(genes$chimpGene %in% CNPG$chimpGene),]

alg <- unique(c(hs,cs,cg,genes$humanGene))




####
x <- list[[1]]
y <- list[[2]]

#z <- cg[which(cg %in% list[[1]]$geneSymbol)]
y <- y[which(y[,1] %in% x[,1]),]
x <- x[-which(duplicated(x$geneSymbol)),]
y <- y[-which(duplicated(y$geneSymbol)),]

x <- x[order(x$geneSymbol),]
y <- y[order(y$geneSymbol),]

x <- x[which(x$geneSymbol %in% alg),]
y <- y[which(y$geneSymbol %in% alg),]
#calculate the absolute difference between SVA contatining genes and their ortholog
diff <- cbind(x[,1],abs(x[,3:8] 
                        - y[,3:8]))
#half the difference for the first and last time point
diff[,2] <- diff[,2]/2
diff[,7] <- diff[,7]/2



diff <- data.frame(diff[,1],rowSums(diff[,-1]), stringsAsFactors = F)
diff <- cbind(diff,c(rep(0,dim(diff)[1])))
colnames(diff) <- c("Gene","Difference","Class")



diff[which(diff[,1] %in% cg),3] <- "Common SVAs"
diff[which(diff[,1] %in% cs),3] <- "Chimp SVAs"
diff[which(diff[,1] %in% hs),3] <- "Human SVAs"


diff[which(diff[,3] == 0),3] <- "Non SVAs"





diff <- cbind(diff,log2(diff$Difference))
colnames(diff)[4] <- "Log2Difference"



pdf("HCdensityplot200529.pdf")
ggplot(diff, aes(x = Log2Difference, fill = Class, linetype = Class)) + geom_density(alpha = 1)+
  scale_y_continuous(labels = percent_format()) + theme_bw()+
  ggtitle("Human-Chimp Gene Expression Differences in Development")
dev.off()


#


x <- list[[1]]
y <- list[[2]]

#z <- cg[which(cg %in% list[[1]]$geneSymbol)]
y <- y[which(y[,1] %in% x[,1]),]
x <- x[-which(duplicated(x$geneSymbol)),]
y <- y[-which(duplicated(y$geneSymbol)),]

x <- x[order(x$geneSymbol),]
y <- y[order(y$geneSymbol),]

x <- x[which(x$geneSymbol %in% alg),]
y <- y[which(y$geneSymbol %in% alg),]
#calculate the absolute difference between SVA contatining genes and their ortholog
diffdir <- cbind(x[,1],(x[,3:8] 
                        - y[,3:8]))
#half the difference for the first and last time point
diffdir[,2] <- diffdir[,2]/2
diffdir[,7] <- diffdir[,7]/2



diffdir <- data.frame(diffdir[,1],rowSums(diffdir[,-1]), stringsAsFactors = F)
diffdir <- cbind(diffdir,c(rep(0,dim(diffdir)[1])))
colnames(diffdir) <- c("Gene","Difference","Class")



diffdir[which(diffdir[,1] %in% cg),3] <- "Common SVAs"
diffdir[which(diffdir[,1] %in% cs),3] <- "Chimp SVAs"
diffdir[which(diffdir[,1] %in% hs),3] <- "Human SVAs"


diffdir[which(diffdir[,3] == 0),3] <- "Non SVAs"
######



diff <- diff[order(diff$Difference, decreasing = T),]

diffdir <- diffdir[order(diffdir$Difference, decreasing = T),]


diff2 <- diff[-which(diff$Class == "Non SVAs"),]
diffdir2 <- diffdir[-which(diffdir$Class == "Non SVAs"),]


pdf("HCdirdensityplot200615.pdf")
ggplot(diffdir2, aes(x = Difference, fill = Class, linetype = Class)) + geom_density(alpha = 0.5)+
  scale_y_continuous(labels = percent_format()) + theme_bw()+
  ggtitle("Human-Chimp Gene Expression Differences in Development")
dev.off()


TOPLIST <- diffdir2[which(abs(diffdir2$Difference) > 10000),]
write.table(TOPLIST, "200615TOPLIST.txt")


i <- NULL
meltlist <- list()
for(i in 1:4){
  
  meltlist[[i]] <- melt(list[[i]])
  
}




meltDF <- do.call(rbind,meltlist)

species <- meltDF$variable
species <- gsub("w.*","",species)
species <- gsub(".*_","",species)
week <- gsub(".*w","",meltDF$variable)

colnames(meltDF)[1] <- "Gene"

meltDF <- cbind(meltDF,species,week)

meltDF <- merge(meltDF,diffdir[,-c(2)], by = "Gene")

mDF <- meltDF[-which(meltDF$Class == "Non SVAs"),]

colnames(mDF)[4] <- "baseMean"

i <- NULL
for(i in 1:dim(TOPLIST)[1]){
  
  #pdf(paste0("200615_",TOPLIST$Gene[i],"_Cross-Species_Expr.pdf"))
  plot <- ggplot(mDF[which(mDF$Gene == TOPLIST$Gene[i]),], aes(x=week, y=baseMean, group=species, color=species)) + 
    geom_line() +
    geom_point()+
    theme_classic() +
    ggtitle(paste(TOPLIST$Gene[i],TOPLIST$Class[i], sep = " ")) 
  
  ggsave(paste0("200615_",TOPLIST$Gene[i],"_Cross-Species_Expr.pdf"), plot = plot)
  #dev.off()
  
  
  
}


########

pdf("20_04_24_MedianGeneExprALLCONTROL.pdf")
i <- NULL
g <- list()
for(i in 1:4){
  
  g[[i]] <- ggplot(results[[i]], aes(x=week, y=median_expression, group=species, color=species)) + 
    geom_line() +
    geom_point()+
    ylim(lim[i,1],lim[i,2])+
    theme_classic() +
    ggtitle(CELLS[i]) 
  
  
  
  
  
  
}

grid.arrange(g[[1]], g[[3]], g[[2]],g[[4]], ncol=2, nrow = 2)



dev.off()





##test on random set of genes assigned to interest condition
set.seed(132)
FakeSet <- sample(alg,250)
alg1 <- alg[-which(alg %in% FakeSet)]

i <- NULL
l <- 0
#IGs <- list()
DIFF <- c()
for(i in 1:9999){
  
  set.seed(i*3)
  
  g <- sample(alg1,250)
  
  DIFF[i] <- TSDiff(list[[1]],list[[2]],g)
  
  #if(DIFF[i] > Common){
  #  l <- l + 1 
  #  IGs[[l]] <- g
    
  }
  
  
#}

FS <- TSDiff(list[[1]],list[[2]],FakeSet)/250



#GOIS <- do.call(c,IGs)
#SUMMARY <- as.data.frame(table(GOIS))

DIFF <- DIFF/length(g)

DIFF <- c(DIFF,FS)



d <- density(DIFF)
pdf("200610ALL5.pdf")
plot(d, main = "Difference Density Distribution")
abline(v = FS, col="red", lwd=3, lty=2)

legend("topright", 
       legend = c("Fake Set"), 
       col = c("red"), pch = 1)
#text(SVADiff + 100, 2000000, "SVA Group", srt = 0.2, col = "red")
dev.off()









#########