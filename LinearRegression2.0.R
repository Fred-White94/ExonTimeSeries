#Fred White 01/10/20
#take DESeq results table, preprocess and perform bootstrapped linear regression for candidate selection
#######
library(nlme)
library(ggplot2)
library(reshape2)
library(tidyverse)

GOt <- read.table("GOIDtranslate.txt", sep = "\t", stringsAsFactors=F)
colnames(GOt)[1] <- "geneSymbol"
Diffs <- read.table("20_08_10_Diffs.txt", stringsAsFactors=F, header = T)
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


#split data into individual species to apply species specific expression filter
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


#get rid of genes with zeroes 
i <- NULL
for(i in 1:4){
  list[[i]] <- list[[i]][-which(list[[i]]$"ZEROES[, i]" > 0),]
  list[[i]] <- list[[i]][-which(list[[i]]$"X[, i]" > 0),]
  list[[i]] <- list[[i]][,-c(9:10)]
  
}


#get list of genes that are expressed in all 4 species
Genes <- Reduce(intersect, list(list[[1]]$geneSymbol,list[[2]]$geneSymbol,list[[3]]$geneSymbol,list[[4]]$geneSymbol))
i <- NULL
for(i in 1:4){
  list[[i]] <- list[[i]][which(list[[i]]$geneSymbol %in% Genes),]
}
DESeq <- DESeq[which(DESeq$geneSymbol %in% list[[1]]$geneSymbol),]



des <- DESeq
des <- cbind(des,apply(des[,-c(1,2)],1,var))
colnames(des)[27] <- "var"

des <- merge(des,Diffs[,c(1,8:10)], by.x = "geneSymbol", by.y = "Gene")
des <- des[which(des$geneSymbol %in% Diffs$Gene),]

des <- cbind(des,paste(des$h,des$c,des$o, sep = "_"))

colnames(des)[31] <- "SVA_Species"
des$SVA_Species <- as.character(des$SVA_Species)

des <- des[!duplicated(des$geneSymbol),]
des2 <- des



ne <- c("0_1_2","0_2_1","1_0_2","1_1_2","1_2_0","1_2_1","1_3_0","2_0_1","2_1_0","2_1_1","2_1_2","2_2_1",
        "3_1_0","3_2_0","3_2_1","3_4_0","4_2_0","4_3_0")

des2 <- des2[-which(des2[,31] %in% ne),]

##convert SVA number labelling of previous dataset to strings to be used for graph labelling
des2[which(des2$SVA_Species %in% c("0_0_1","0_0_2","0_0_3")),31] <- "Orang_spec"
des2[which(des2$SVA_Species %in% c("0_1_0","0_2_0")),31] <- "Chimp_spec"
des2[which(des2$SVA_Species %in% c("0_1_1","0_2_2")),31] <- "Chimp_Orang"
des2[which(des2$SVA_Species %in% c("1_0_0","2_0_0","3_0_0")),31] <- "Human_spec"
des2[which(des2$SVA_Species %in% c("0_0_0")),31] <- "No_SVA"
des2[which(des2$SVA_Species %in% c("1_1_0","2_2_0","3_3_0")),31] <- "Human_Chimp"
des2[which(des2$SVA_Species %in% c("1_1_1")),31] <- "Human_Chimp_Orang"
des2[which(des2$SVA_Species %in% c("1_0_1")),31] <- "Human_Orang"

#remove previous SVA species number columns
des2 <- des2[,-c(27:30)]


#get rid of these genes as the class sizes are too small (less powerful regression model if included)
des2 <- des2[-which(des2$SVA_Species %in% c("Chimp_Orang","Human_Chimp_Orang","Human_Orang")),]



des2 <- melt(des2)

species <- des2$variable
species <- gsub("w.*","",species)
species <- gsub(".*_","",species)
week <- gsub(".*w","",des2$variable)

colnames(des2)[1] <- "Gene"
des2 <- cbind(des2,species,week)
des2$id <- paste(des2$Gene, des2$variable, sep = "_")
colnames(des2)[3] <- "Class"
des2$Class <- as.factor(des2$Class)
des2$Class <- relevel(des2$Class, ref= "No_SVA")
des2$week <- relevel(des2$week, ref = "1")
des2$species <- relevel(des2$species, ref = "rES")
colnames(des2)[5] <- "baseMean"


DES2save <- des2


#get rid of week 0 as cell organoids made from different cell types
des2 <- des2[-which(des2$week == 0),]
#also get rid of week 5 as these were harvested slightly differently 
des2 <- des2[-which(des2$week == 5),]


#des2 <- des2[which(des2$species %in% c("hES","rES")),]
des2$Class[which(des2$Class == "Human_Chimp")] <- "Human_spec"
des2$Class[which(des2$Class %in% c("Orang_spec","Chimp_spec"))] <- "No_SVA"
des2$Class <- droplevels(des2$Class)

hr <- des2[which(des2$species %in% c("hES","rES")),]

d2 <- des2
des2 <- hr

###
nsamp <- min(as.data.frame(table(des2$Class))$Freq)


#mean centre data for lme function which is optimised for this type of input
des2$baseMean <- scale(des2$baseMean)
des2$cellType <- gsub("^.","",des2$species)


#sample data and make i number of regression models both with and without SVA type information 
i <- NULL
res <- list()
for(i in 1:1000){
  
  
  set.seed(i*2)

  tr <- des2 %>% 
    group_by(Class) %>% 
    sample_n(nsamp)
  
  
  
  model6 <- lme(baseMean ~ species + week + Class, data = tr, random = ~1|Gene, control = lmeControl(opt = "optim"))
  model7 <- lme(baseMean ~ species + week, data = tr, random = ~1|Gene, control = lmeControl(opt = "optim"))
  
  x <- cbind(residuals(model6), residuals(model7))
  colnames(x) <- c("W/Class","W/oClass")
  
  #get residuals for each gene 
  y <- cbind(as.data.frame(tr), x)
  
  res[[i]] <- y
  
}



a <- do.call(rbind, res)


#remove non SVA containing genes here for larger sample size for averaging of residuals
b <- a[-which(a$Class == "No_SVA"),]
nsamp <- min(table(b$id))

b <- b %>% 
  group_by(id) %>% 
  sample_n(nsamp)

b <- as.data.frame(b)




#####
#calculate R2 values
b$R2W <- b$'W/Class'^2
b$R2Wo <- b$'W/oClass'^2







#calculate gene wise R2 for the models
grs <- b %>% 
  group_by(Gene) %>%
  summarise(GRS = sum(R2W))

grsmod2 <- b %>% 
  group_by(Gene) %>%
  summarise(GRS = sum(R2Wo))

GRS <- cbind(grs,grsmod2[,2])
colnames(GRS)[2:3] <- c("With","Without")

GRS$diff <- abs(GRS[,3]) - abs(GRS[,2])
GRS <- GRS[which(GRS$diff > 0),]

GRS$fract <- abs(GRS$diff)/abs(GRS[,3])
GRS <- GRS[order(GRS$fract, decreasing = T),]

hspc <- b[which(b$Class == "Human_spec"),]
cspc <- b[which(b$Class == "Chimp_spec"),]


grshs <- GRS[which(GRS$Gene %in% hspc$Gene),]


######
#to plot real expression values
#des2 <- DES2save
#des2 <- des2[-which(des2$week == 0),]

grshs <- grshs[sample(rownames(grshs),180),]



i <- NULL

for(i in 1:dim(grshs)[1]){
  
  x <- grshs[i,1]
  plot <- ggplot(des2[which(des2$Gene == x),], aes(x=week, y=baseMean, group=species, color=species)) + 
    geom_line() +
    geom_point()+
    theme_classic() +
    ggtitle(paste(x,des2$Class[which(des2$Gene == x)][1], sep = " ")) 
  
  ggsave(paste0("20_08_13/",x,"_Cross-Species_Expr.pdf"), plot = plot)
  
  
  
}





########
library(GenomicRanges)
Cands <- HRsave[which(HRsave$Gene %in% grshs$Gene),]

peaks <- read.table("SVA_vs_H3K4me3_ZNF91ko13_peaks_all.txt", sep = "\t")
gencan32 <- read.table("KnownCanonicalHG38Gencode32.txt",stringsAsFactors=F)
kgx <- read.table("kgXref2.txt", sep = "\t", header = T)


colnames(peaks) <- c("chr","start","end","SVA","peak","strand")
peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns=T)

colnames(gencan32) <- c("chr","start","end","score","transcript","protein")

gencan32 <- merge(gencan32, kgx, by.x = "transcript", by.y = "kgID")
gencan32 <- makeGRangesFromDataFrame(gencan32, keep.extra.columns=T)

b <- findOverlaps(gencan32, peaks, ignore.strand = T)
b <- as.data.frame(b)

c <- as.character(unique(gencan32$geneSymbol[b$queryHits]))
#hHits <- cbind.data.frame(gencan32$geneSymbol[b$queryHits], peaks$Length[b$subjectHits], stringsAsFactors = F)








