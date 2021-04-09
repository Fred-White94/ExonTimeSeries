#Fred White 01/10/2020
#Functions and code for cross species analysis cleaning

library(GenomicRanges)
library(homologene)
library(ggplot2)
library(scales)

#do this only once per session
#homologeneData <- updateHomologene()


hg38rep <- read.table("Repeatshg38.txt",stringsAsFactors=F, sep = "\t")
panTrorep <- read.table("RepeatsPanTro6.txt", stringsAsFactors=F, sep = "\t")
ponAberep <- read.table("RepeatsPonAbe3.txt", stringsAsFactors=F, sep = "\t")


SVAfilt <- function(x){
  
  x$V9 <- gsub(".* transcript_id ","",x$V9)
  x$V9 <- gsub("_dup.*","",x$V9)
  x$V9 <- gsub("; ","",x$V9)
  ind <- grep("SVA.*",x$V9)
  x <- x[ind,]
  x <- cbind(x, abs(x$V4 - x$V5))
  colnames(x) <- c("chr","Genome","Transcript","start","end","whoknows","strand","nothing","Name","Length")
  x <- x[,-c(2,3,6,8)]
}

hg38rep <- SVAfilt(hg38rep)
panTrorep <- SVAfilt(panTrorep)
ponAberep <- SVAfilt(ponAberep)





HG <- read.table("hg38.ncbiRefSeq(1).gtf", stringsAsFactors=F, sep = "\t")
ChG <- read.table("/zfs/sils/jacobs/SRP121791/sra/panTro6.ncbiRefSeq.gtf", stringsAsFactors=F, sep = "\t")
OG <- read.table("/zfs/sils/jacobs/SRP121791/sra/ponAbe3.ncbiRefSeq.gtf", stringsAsFactors=F, sep = "\t")




tgtfclean <- function(HG){
  
  HT <- HG[which(HG$V3 == "transcript"),]
  HT$V9 <- gsub("gene_id ","",HT$V9)
  HT$V9 <- gsub(";.*","",HT$V9)
  l <- abs(HT[,4] - HT[,5])
  HT <- cbind(HT,l)
  HT <- HT[order(HT$l, decreasing = T),]
  HT <- HT[!duplicated(HT$V9),]
  HT <- HT[-which(HT$l < 500),]
  
  HT <- HT[,-c(2,3,6,8,10)]
  colnames(HT) <- c("chr","start","end","strand","Name")
  
  return(HT)
  
  
}

HT <- tgtfclean(HG)
ChT <- tgtfclean(ChG)
OT <- tgtfclean(OG)




##take two data frames (i.e. processed gtf files) and return 
hits <- function(x,y){
#input must be ready for Grange conversion
  
  x$start <- x$start - 50
  x$end <- x$end + 50
  
  Grange <- makeGRangesFromDataFrame(x, keep.extra.columns=T)
  
  Grange <- reduce(Grange)
  start(Grange@ranges) <- start(Grange@ranges) + 50
  end(Grange@ranges) <- end(Grange@ranges) - 50
  
  Tr <- makeGRangesFromDataFrame(y, keep.extra.columns = T)
  names(Tr) <- Tr$Name
  b <- findOverlaps(Tr, Grange, ignore.strand = T)
  b <- as.data.frame(b)
  
  #names of genes with SVAs and the length of the SVAs 
  Hits <- cbind.data.frame(names(Tr[b$queryHits,]), width(Grange@ranges)[b$subjectHits], start(Grange@ranges)[b$subjectHits], end(Grange@ranges)[b$subjectHits], stringsAsFactors = F)
  colnames(Hits) <- c("Gene", "SVAlength","SVAstart","SVAend")

  return(Hits)
  

}


#HNPG <- as.data.frame(table(Hits$Gene))
#colnames(HNPG) <- c("Gene", "Hits")

hHits <- hits(hg38rep,HT)
cHits <- hits(panTrorep,ChT)
oHits <- hits(ponAberep,OT)

#hhs <- hHits
#chs <- cHits
#ohs <- oHits


hHits <- hHits[which(hHits$SVAlength > 999),]
cHits <- cHits[which(cHits$SVAlength > 999),]
oHits <- oHits[which(oHits$SVAlength > 999),]

H <- unique(hHits$Gene)
C <- unique(cHits$Gene)
O <- unique(oHits$Gene)


HCO <- Reduce(intersect, list(H,C,O))

HC <- Reduce(intersect, list(H, C))
HC <- HC[-which(HC %in% HCO)]

HO <- Reduce(intersect, list(H, O))
HO <- HO[-which(HO %in% HCO)]


CO <- Reduce(intersect, list(C, O))
CO <- CO[-which(CO %in% HCO)]


ALLSVA <- unique(c(H,C,O))

HSpec <- H[-which(H %in% C)]
HSPec <- HSpec[-which(HSpec %in% O)]


CSpec <- C[-which(C %in% H)]
CSpec <- CSpec[-which(CSpec %in% O)]


OSpec <- O[-which(O %in% H)]
OSpec <- OSpec[-which(OSpec %in% C)]

###cannot convert orangutan genes......

#adapt human2mouse function from homologene package
#human2chimp <- function (genes, db = homologeneData)
#{
#  out = homologene(genes, 9606, 9598, db)
#  names(out) = c("humanGene", "chimpGene", "humanID", "chimpID")
#  return(out)
#}


HNPG <- as.data.frame(table(hHits$Gene), stringsAsFactors = F)
colnames(HNPG) <- c("Gene", "Hits")

CNPG <- as.data.frame(table(cHits$Gene), stringsAsFactors = F)
colnames(CNPG) <- c("Gene", "Hits")

ONPG <- as.data.frame(table(oHits$Gene), stringsAsFactors = F)
colnames(ONPG) <- c("Gene", "Hits")

HNPG[which(HNPG$Gene %in% HCO),]
CNPG[which(CNPG$Gene %in% HCO),]
ONPG[which(ONPG$Gene %in% HCO),]



Cm <- merge(HNPG,CNPG,by = "Gene", all = T)
Am <- merge(HCm, ONPG, by = "Gene", all = T)
Am[is.na(Am)] <- 0
colnames(Am) <- c("Gene","h","c","o")


Am <- cbind(Am,apply(Am[,-1],1,var))
colnames(Am) <- c("Gene","h","c","o","var")
  

Am <- Am[which(Am$Gene %in% list[[1]]$geneSymbol),]


DESeq <- DESeq[which(DESeq$geneSymbol %in% list[[1]]$geneSymbol),]








difference <- function(x,y){
  
  
  
  y <- y[which(y[,1] %in% x[,1]),]
  x <- x[-which(duplicated(x$geneSymbol)),]
  y <- y[-which(duplicated(y$geneSymbol)),]
  
  x <- x[order(x$geneSymbol),]
  y <- y[order(y$geneSymbol),]
  
  #x <- x[which(x$geneSymbol %in% alg),]
  #y <- y[which(y$geneSymbol %in% alg),]
  #calculate the absolute difference between SVA contatining genes and their ortholog
  diff <- cbind(x[,1],abs(x[,3:8] 
                          - y[,3:8]))
  #half the difference for the first and last time point
  diff[,2] <- diff[,2]/2
  diff[,7] <- diff[,7]/2
  
  
  
  diff <- data.frame(diff[,1],rowSums(diff[,-1]), stringsAsFactors = F)
  diff <- cbind(diff,c(rep(0,dim(diff)[1])))
  
  
  
  colnames(diff) <- c("Gene","Difference","Class")
  
  return(diff)
  
  
}



HC <- difference(list[[1]],list[[2]])
HO <- difference(list[[1]],list[[3]])
HR <- difference(list[[1]],list[[4]])
CO <- difference(list[[2]],list[[3]])
CR <- difference(list[[2]],list[[4]])
OR <- difference(list[[3]],list[[4]])

Diffs <- cbind(HC[,c(1,2)], HO[,2], HR[,2], CO[,2], CR[,2], OR[,2])
colnames(Diffs) <- c("Gene","HC","HO","HR","CO","CR","OR")


Diffs <- merge(Diffs,Am, by = "Gene", all = T)

Diffs[is.na(Diffs)] <- 0

Diffs <- Diffs[order(Diffs$HR, decreasing = T),]






meltDF <- melt(DESeq)

species <- meltDF$variable
species <- gsub("w.*","",species)
species <- gsub(".*_","",species)
week <- gsub(".*w","",meltDF$variable)

colnames(meltDF)[1] <- "Gene"

meltDF <- cbind(meltDF,species,week)
mDF <- meltDF
#meltDF <- merge(meltDF,diffdir[,-c(2)], by = "Gene")

#mDF <- meltDF[-which(meltDF$Class == "Non SVAs"),]

colnames(mDF)[4] <- "baseMean"

x <- "THOC5"
  #pdf(paste0("200615_",TOPLIST$Gene[i],"_Cross-Species_Expr.pdf"))
  plot <- ggplot(mDF[which(mDF$Gene == x),], aes(x=week, y=baseMean, group=species, color=species)) + 
    geom_line() +
    geom_point()+
    theme_classic() +
    ggtitle(x) 
  
  ggsave(paste0("200617_",x,"_Cross-Species_Expr.pdf"), plot = plot)



rownames(Diffs) <- Diffs$Gene
Diffs <- Diffs[,-c(1)]
  





Diffs <- Diffs[order(Diffs$HC, decreasing = T),]
hs1 <- Diffs[which(Diffs$h == 1 & Diffs$c == 0 & Diffs$o == 0 & Diffs$HC > Diffs$CR),][1:6,]



Diffs <- Diffs[order(Diffs$CO, decreasing = T),]
hc1 <- Diffs[which(Diffs$h == 1 & Diffs$c == 1 & Diffs$o == 0 & Diffs$CO < Diffs$OR),][1:6,]



Diffs <- Diffs[order(Diffs$CO, decreasing = T),]
ho1 <- Diffs[which(Diffs$h == 1 & Diffs$c == 0 & Diffs$o == 1 & Diffs$HO < Diffs$HC),][1:6,]


Diffs <- Diffs[order(Diffs$HO, decreasing = T),]
co1 <- Diffs[which(Diffs$h == 0 & Diffs$c == 1 & Diffs$o == 1 & Diffs$CO < Diffs$HC),][1:6,]


Diffs <- Diffs[order(Diffs$HC, decreasing = T),]
cs1 <- Diffs[which(Diffs$h == 0 & Diffs$c == 1 & Diffs$o == 0 & Diffs$HC > Diffs$HO),][1:6,]


Diffs <- Diffs[order(Diffs$CO, decreasing = T),]
os1 <- Diffs[which(Diffs$h == 0 & Diffs$c == 0 & Diffs$o == 1 & Diffs$HO > Diffs$HC),][1:6,]


cands1 <- rbind(hs1,hc1,ho1,co1,cs1,os1)
cands1 <- cbind(cands1, rep(c("HS","HC","HO","CO","CS","OS"), each = 6))
colnames(cands1)[11] <- "SVA_Species"




i <- NULL

for(i in 1:dim(cands1)[1]){
  
  x <- rownames(cands1)[i]
  #pdf(paste0("200618_",x,"_Cross-Species_Expr.pdf"))
  plot <- ggplot(mDF[which(mDF$Gene == x),], aes(x=week, y=baseMean, group=species, color=species)) + 
    geom_line() +
    geom_point()+
    theme_classic() +
    ggtitle(paste(x,cands1$SVA_Species[i], sep = " ")) 
  
  ggsave(paste0("200618_",x,"_Cross-Species_Expr.pdf"), plot = plot)
  
  #dev.off()
  
  
  
  
  
}



Diffs <- Diffs[order(Diffs$HR, decreasing = T),]
hco1 <- Diffs[which(Diffs$h == 1 & Diffs$c == 1 & Diffs$o == 1 & Diffs$HR > Diffs$HC),][1:6,]




diffs2 <- Diffs




diffs2 <- cbind(diffs2,apply(diffs2[,c(1:6)],1,var))

colnames(diffs2)[11] <- c("Diff_var")



x <- "LRRTM4"

plot <- ggplot(mDF[which(mDF$Gene == x),], aes(x=week, y=baseMean, group=species, color=species)) + 
  geom_line() +
  geom_point()+
  theme_classic() +
  ggtitle(paste(x,cands1$SVA_Species[i], sep = " ")) 

ggsave(paste0("200623_",x,"_Cross-Species_Expr.pdf"), plot = plot)
















Diffs <- Diffs[order(Diffs$HC, decreasing = T),]
hs2 <- Diffs[which(Diffs$h == 1 & Diffs$c == 0 & Diffs$o == 0 & Diffs$HC > Diffs$CR),][1:6,]



Diffs <- Diffs[order(Diffs$CO, decreasing = T),]
hc2 <- Diffs[which(Diffs$h == 1 & Diffs$c == 1 & Diffs$o == 0 & Diffs$CO < Diffs$OR),][1:6,]



Diffs <- Diffs[order(Diffs$CO, decreasing = T),]
ho2 <- Diffs[which(Diffs$h == 1 & Diffs$c == 0 & Diffs$o == 1 & Diffs$HO < Diffs$HC),][1:6,]


Diffs <- Diffs[order(Diffs$HO, decreasing = T),]
co2 <- Diffs[which(Diffs$h == 0 & Diffs$c == 1 & Diffs$o == 1 & Diffs$CO < Diffs$HC),][1:6,]


Diffs <- Diffs[order(Diffs$HC, decreasing = T),]
cs2 <- Diffs[which(Diffs$h == 0 & Diffs$c == 1 & Diffs$o == 0 & Diffs$HC > Diffs$HO),][1:6,]


Diffs <- Diffs[order(Diffs$CO, decreasing = T),]
os2 <- Diffs[which(Diffs$h == 0 & Diffs$c == 0 & Diffs$o == 1 & Diffs$HO > Diffs$HC),][1:6,]


cands2 <- rbind(hs2,hc2,ho2,co2,cs2,os2)
cands2 <- cbind(cands2, rep(c("HS","HC","HO","CO","CS","OS"), each = 6))
colnames(cands2)[12] <- "SVA_Species"




i <- NULL

for(i in 1:dim(cands2)[1]){
  
  x <- cands2[i,1]
  #pdf(paste0("200618_",x,"_Cross-Species_Expr.pdf"))
  plot <- ggplot(mDF[which(mDF$Gene == x),], aes(x=week, y=baseMean, group=species, color=species)) + 
    geom_line() +
    geom_point()+
    theme_classic() +
    ggtitle(paste(x,cands1$SVA_Species[i], sep = " ")) 
  
  ggsave(paste0("200623_",x,"_Cross-Species_Expr.pdf"), plot = plot)
  
  #dev.off()
  
  
  
  
  
}

