#Fred White 01/10/2020
#Exon level RNAseq analysis - analyse canonical transcript expression across species 

library(homologene)
library(tidyr)
library(reshape2)
library(viridis)



CANDS <- read.table("20_08_17CandListHR.txt", stringsAsFactors=F)

#do this only once per session
#homologeneData <- updateHomologene()





kgXref <- read.table("kgXref2.txt",stringsAsFactors=F,header =T)
kgXref <- kgXref[-which(duplicated(kgXref$kgID)),]

gencan32 <- read.table("KnownCanonicalHG38Gencode32.txt",stringsAsFactors=F)


Gencode32CDS <- read.table("Gencode32CDS.txt",stringsAsFactors=F)
Gencode32CDS <- Gencode32CDS[which(Gencode32CDS$V1 %in% gencan32$V5),]

Gencode32CDS <- Gencode32CDS[-which(Gencode32CDS$V2 == "none"),]
#Gencode32CDS <- Gencode32CDS[-which(Gencode32CDS$V3 == "none"),]

#take only protein coding genes 
#GOt <- read.table("GOIDtranslate.txt", sep = "\t", stringsAsFactors=F)





hg38 <- read.table("/zfs/sils/jacobs/SRP121791/sra/hg38.knownGene.gtf", stringsAsFactors=F, header = F, sep = "\t")



egtfclean <-function(hg38){
  
  hg38a <- hg38[which(hg38$V3 == "exon"),]
  Gids <- gsub(".*transcript_id ","",hg38a[,9])
  Gids <- gsub("\\;.*","",Gids)
  Gids <- gsub("\\\"","",Gids)
  Eids <- gsub(".*exon_id ","",hg38a[,9])
  Eids <- gsub("\\;.*","",Eids)
  Eids <- gsub("\\\"","",Eids)
  hg38a <- hg38a[,-9]
  hg38a <- cbind(hg38a,Gids,Eids)
  
  
  return(hg38a)
  
  
}

hg38a <- egtfclean(hg38)

hg38a <- hg38a[which(hg38a$Gids %in% gencan32$V5),]
colnames(hg38a)[9] <- "kgID"
hg38a <- merge(hg38a,kgXref, by = "kgID")
hg38a <- merge(hg38a, Gencode32CDS[,c(1,4,5)], by.x = "kgID", by.y = "V1")

###change the columns here
hg38a <- hg38a[order(hg38a$kgID,hg38a$V1,hg38a$V4.x),]



#kgXref for rhesus

RheXref <- read.table("RheXRef.txt",stringsAsFactors=F,header =F)
RheXref <- RheXref[-which(RheXref$V14 == "none"),]

RheXref2 <- RheXref[,c(2,13,16)]

RG <- read.table("/zfs/sils/jacobs/SRP121791/sra/rheMac10.ncbiRefSeq.gtf", stringsAsFactors=F, sep = "\t")

RGa <- egtfclean(RG)

RGa1 <- merge(RGa, RheXref2, by.x = "Gids", by.y = "V2")
RGa1 <- RGa1[order(RGa1$V1,RGa1$V4,RGa1$Gids),]


colnames(hg38a)[13] <- "exonframe"
colnames(RGa1)[12] <- "exonframe"
colnames(RGa1)[1] <- "kgID"


hg38asave <- hg38a

hg38a <- hg38a[which(hg38a$geneSymbol %in% CANDS$Gene),]

transcripts <- unique(hg38a$kgID)

#needs to be optimised
i <- NULL
hg38a2 <- data.frame()
for(i in 1:length(transcripts)){
  
  x <- hg38a[which(hg38a$kgID == transcripts[i]),]
  ef <- unlist(strsplit(x$exonframe[1], split = ","))
  x$exonframe <- ef
  
  hg38a2 <- rbind(hg38a2,x)
  print(i)
  
}


#GTF2IRD1 - fails: 2 different transcripts with the same annotation
#LEPROTL1 - fails

RGa1save <- RGa1

RGa1 <- RGa1[which(RGa1$V13 %in% CANDS$Gene),]

transcripts <- unique(RGa1$kgID)

i <- NULL
RGa2 <- data.frame()
for(i in 1:length(transcripts)){
  
  x <- RGa1[which(RGa1$kgID == transcripts[i]),]
  ef <- unlist(strsplit(x$exonframe[1], split = ","))
  x$exonframe <- ef
  
  RGa2 <- rbind(RGa2,x)
  print(i)
  
}


h <- hg38a2
r <- RGa2

hg38a2 <- hg38a2[-which(hg38a2$exonframe == "-1"),]
RGa2 <- RGa2[-which(RGa2$exonframe == "-1"),]

hg38a2 <- hg38a2[which(hg38a2$geneSymbol %in% RGa2$V13),]


a <- c(1:22,"X","Y")
b <- paste0("chr",a)
hg38a2 <- hg38a2[which(hg38a2[,2] %in% b),]
RGa2 <- RGa2[which(RGa2[,2] %in% b),]

RGa2$ID <- paste(RGa2$kgID, RGa2$V13, sep = "_")
hg38a2$ID <- paste(hg38a2$kgID, hg38a2$geneSymbol, sep = "_")

ht <- as.data.frame(table(hg38a2$ID))
rt <- as.data.frame(table(RGa2$ID))

ht$gene <- gsub(".*_","",ht$Var1)
rt$gene <- gsub(".*_","",rt$Var1)

ht$GE <- paste(ht$gene, ht$Freq, sep = "_")
rt$GE <- paste(rt$gene, rt$Freq, sep = "_")


rt <- rt[which(rt$GE %in% ht$GE),]
RGa2 <- RGa2[which(RGa2$ID %in% rt$Var1),]






hcounts <- read.table("/zfs/sils/jacobs/SRP121791/sra/hcounts.txt", header = T, stringsAsFactors = F)

rcounts <- read.table("/zfs/sils/jacobs/SRP121791/sra/rcounts.txt", header = T, stringsAsFactors = F)




hcounts2 <- cbind(gsub("\\.[^\\.]*$","", hcounts$Geneid),hcounts)
colnames(hcounts2)[1] <- "kgID"
hcounts3 <- hcounts2[which(hcounts2$kgID %in% gsub("_.*","",ht$Var1)),]

hcounts4 <- merge(hcounts3, kgXref, by = "kgID")


rcounts2 <- cbind(gsub("\\.[^\\.]*$","", rcounts$Geneid),rcounts)
colnames(rcounts2)[1] <- "kgID"
rcounts3 <- rcounts2[which(rcounts2$kgID %in% gsub("_[^\\.]*$","",rt$Var1)),]


rcounts4 <- merge(rcounts3,RheXref2[,1:2], by.x = "kgID", by.y = "V2")

hcounts4 <- hcounts4[which(hcounts4$geneSymbol %in% rcounts4$V13),]




rcounts4 <- rcounts4[which(rcounts4$Geneid %in% RGa2$Eids),]
hcounts4 <- hcounts4[which(hcounts4$Geneid %in% hg38a2$Eids),]



A <- rle(paste(as.character(hcounts4$kgID),as.character(hcounts4$Strand), sep = "__"))
strand <- gsub(".*__","",A$values)
a <- A$lengths
b <- rep(1,length(a))
c <- data.frame(cbind(b,a,strand), stringsAsFactors = F)

i <- NULL
EXONS <- c()
for(i in 1:dim(c)[1]){
  
  
  if(c$strand[i] == "+"){
    EXONS <- c(EXONS,c(c[i,1]:c[i,2]))
    
  }
  else if (c$strand[i] == "-"){
    EXONS <- c(EXONS,c(c[i,2]:c[i,1]))
    
  }
  
  
  
  
}

hcounts4$Geneid <- paste(gsub("\\.[^\\.]*$","", hcounts4$Geneid),EXONS,sep = ".")
hcounts4$GE <- paste(hcounts4$geneSymbol, EXONS, sep = "_")



A <- rle(paste(as.character(rcounts4$kgID),as.character(rcounts4$Strand), sep = "__"))
strand <- gsub(".*__","",A$values)
a <- A$lengths
b <- rep(1,length(a))
c <- data.frame(cbind(b,a,strand), stringsAsFactors = F)

i <- NULL
EXONS <- c()
for(i in 1:dim(c)[1]){
  
  
  if(c$strand[i] == "+"){
    EXONS <- c(EXONS,c(c[i,1]:c[i,2]))
    
  }
  else if (c$strand[i] == "-"){
    EXONS <- c(EXONS,c(c[i,2]:c[i,1]))
    
  }
  
  
  
  
}


rcounts4$Geneid <- paste(gsub("\\.[^\\.]*$","", rcounts4$Geneid),EXONS,sep = ".")
rcounts4$GE <- paste(rcounts4$V13, EXONS, sep = "_")


rcounts4 <- rcounts4[order(rcounts4$V13, rcounts4$kgID),]


#a <- unique(hcounts4$geneSymbol[-which(hcounts4$GE %in% rcounts4$GE)])
#hcounts4 <- hcounts4[-which(hcounts4$geneSymbol %in% a),]

i <- NULL
hcounts5 <- data.frame()
rcounts5 <- data.frame()
GENES <- unique(hcounts4$geneSymbol)
for(i in 1:length(GENES)){
  
  htemp <- hcounts4[which(hcounts4$geneSymbol == GENES[i]),]
  rtemp <- rcounts4[which(rcounts4$V13 == GENES[i]),]
  
  rtemp <- rtemp[1:dim(htemp)[1],]
  
  hcounts5 <- rbind(hcounts5,htemp)
  rcounts5 <- rbind(rcounts5, rtemp)
  
}


hcounts5 <- hcounts5[order(hcounts5$GE),]
rcounts5 <- rcounts5[order(rcounts5$GE),]



##########!!!!!!!!!!!!!


#####move this to just after canonical exon sorting stage
#only here because couldnt get all canon exons for rhesus in time
length <- hcounts5$Length
####
RPK <- hcounts5[,-c(1:7,20,21)]/length
RPK[is.na(RPK)] <- 0


TPM <- t( t(RPK) * 1e6 / colSums(RPK) )
TPM <- cbind(hcounts5[,c(1:7,20,21)],TPM)

hcounts5 <- TPM


length <- rcounts5$Length
####
RPK <- rcounts5[,-c(1:7,20,21)]/length
RPK[is.na(RPK)] <- 0


TPM <- t( t(RPK) * 1e6 / colSums(RPK) )
TPM <- cbind(rcounts5[,c(1:7,20,21)],TPM)

rcounts5 <- TPM



MD <- merge(hcounts5,rcounts5[-c(1:8)], by = "GE")

MD2 <- MD[,c(1,9:33)]
MD2 <- melt(MD2)


dictionaryField <- read.table("Dictionary_Field_2014.csv",sep = ",", stringsAsFactors = F)
colnames(dictionaryField)[c(1,11,14)] <- c("Sample","Species","Time_Point")

rownames(dictionaryField) <- dictionaryField[,1]


MD2$variable <- as.character(MD2$variable)
MD2$variable <- gsub(".*ed.","",MD2$variable)
MD2$variable <- gsub(".bam","",MD2$variable)

colnames(MD2)[c(3,4)] <- c("Sample","TPM")

MD2 <- merge(MD2,dictionaryField[,c(1,11,14)], by = "Sample")

MD2$Exon <- gsub(".*_","",MD2$GE)
MD2$Exon <- as.numeric(MD2$Exon)

MD2$Species[which(MD2$Species == "Homo sapiens")] <- "Human"
MD2$Species[which(MD2$Species == "Macaca mulatta")] <- "Rhesus"
MD2$humanGene <- gsub("_.*","",MD2$GE)
colnames(MD2)[6]  <- "week"
MD2$Condition <- paste(MD2$Species,MD2$week, sep = " ")

MD2 <- MD2[-which(MD2$week %in% c("week 0", "week 5")),]
##########################




#G4EXON <- read.table("G4EXON2.txt")




CANDS <- CANDS[which(CANDS$Gene %in% MD2$humanGene),]
#CANDS <- CANDS[which(CANDS$Gene %in% G4EXON$ID),]
#psen1 exon 3 (if removing out of frame exons)
#svapos <- c(1,1,3,1,11,1,5,14,2,45,1,4,1,6,1,5,5,3,4,5)
svapos <- c(5,24,1,1,12,3,6)
svapos <- svapos + 0.5
GENES <- c("GPATCH2","C2CD3","STAG1","ARFGEF2","CDK5RAP2","NAV2","CAMTA1")
#GENES <- CANDS[1:20,1]
i <- NULL
for(i in 1:length(GENES)){
  
  
  gene <- GENES[i]
  
  data <- MD2[which(MD2$humanGene == gene),]
  
  #svapos <- G4EXON[which(G4EXON$ID == gene),5] + 0.5
  
  #X <- data$exon
  plot<- ggplot(data, aes(x = Exon, y = TPM, group = Condition)) +
    #geom_line(aes(color = Condition, linetype = Species, stat = "identity")) +
    stat_summary(geom = "line", fun = mean, aes(color = Species, linetype = Species, group = Species), size = 2) +#, bins = 12) +
    #scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
    scale_x_continuous(breaks = c(1:max(data$Exon))) +
    
    
    geom_vline(xintercept = svapos[i], linetype = "dotted", color = "black", size = 1.5) +
    
    
    #geom_line(aes(color = Condition, linetype = Species, group = ID)) +
    #geom_point(aes(color = Condition)) +
    #scale_fill_viridis(discrete = T, option = "E") +
    ggtitle(gene) +
    #ylim(0,800) +
    
    #facet_wrap(~exon) +
    #theme_ipsum() +
    #theme(legend.position="none") +
    #theme(axis.text.x=element_blank()) + 
    theme_classic() +
    theme(axis.title=element_text(size=14), axis.text=element_text(size=15)) +
    theme(legend.text=element_text(size=12)) +
    theme(legend.title=element_text(size=14))
  #xlim(1,12)
  
  ggsave(paste0("200930_",GENES[i],"_Exon_Expression.pdf"), plot = plot)
  
  
  
}









RGenes <- unique(RGa1$V13)

####


#adapt human2mouse function from homologene package
rhesus2human <- function (genes, db = homologene::homologeneData)
{
  out = homologene(genes, 9544, 9606, db)
  names(out) = c("rhesusGene", "humanGene", "rhesusID", "humanID")
  return(out)
}



RHO <- rhesus2human(RGenes)

RGa2 <- merge(RGa1, RHO[,c(1,2)], by.x = "V13", by.y = "rhesusGene")
RGa2 <- RGa2[order(RGa2$humanGene, RGa2$Gids),]

hg38a <- hg38a[order(hg38a$geneSymbol, hg38a$kgID),]






