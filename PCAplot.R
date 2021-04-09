#Fred White 01/10/2020
#PCA for DESeq data

library(ggfortify)


DESeq <- read.table("FieldDESeq.txt", sep = "\t", header = T, stringsAsFactors=F)



df <- DESeq[,-c(1,2)]
df <- t(df)

species <- rep(c("Human","Chimp","Orangutan","Rhesus"), 6, each = T)
time <- paste("week ",rep(c(0:5),4))



PCD <- prcomp(df, scale. = F)
df <- cbind(species,time, df)

pdf("test.pdf")
autoplot(PCD, data = df, colour = 'time')
dev.off()

#########################################################################







diffs2 <- Diffs

diffs2 <- cbind(diffs2,paste(diffs2$h, diffs2$c, diffs2$o, sep = "_"))

ne <- c("0_1_2","0_2_1","1_0_2","1_1_2","1_2_0","1_2_1","1_3_0","2_0_1","2_1_0","2_1_1","2_1_2","2_2_1",
        "3_1_0","3_2_0","3_2_1","3_4_0","4_2_0","4_3_0")

diffs2 <- diffs2[-which(diffs2[,13] %in% ne),]

colnames(diffs2)[13] <- "SVA_Species"

diffs2 <- diffs2[,-c(8:12)]




diffs2 <- cbind(diffs2,apply(diffs2[,-c(1,8)],1,var))

colnames(diffs2)[9] <- "var"

diffs2 <- diffs2[-which(diffs2$SVA_Species == "0_0_0"),]

diffs2 <- diffs2[which(diffs2$var > quantile(diffs2$var, 0.95)),]



###
df <- DESeq[which(DESeq$geneSymbol %in% diffs2$Gene),]
df <- merge(df,diffs2[,c(1,8)], by.x = "geneSymbol", by.y = "Gene")


PCD <- prcomp(df[,-c(1,2,27)], scale. = F)

pdf("test.pdf")
autoplot(PCD, data = df, colour = "SVA_Species")
dev.off()


####PCA looks like trash




##try training random forest on gene class
des <- DESeq

des <- cbind(des,apply(des[,-c(1,2)],1,var))

colnames(des)[27] <- "var"


des <- merge(des,Diffs[,c(1,8:10)], by.x = "geneSymbol", by.y = "Gene")
des <- des[which(des$geneSymbol %in% Diffs$Gene),]


des <- cbind(des,paste(des$h,des$c,des$o, sep = "_"))

colnames(des)[31] <- "SVA_Species"
des$SVA_Species <- as.character(des$SVA_Species)

des2 <- des


ne <- c("0_1_2","0_2_1","1_0_2","1_1_2","1_2_0","1_2_1","1_3_0","2_0_1","2_1_0","2_1_1","2_1_2","2_2_1",
        "3_1_0","3_2_0","3_2_1","3_4_0","4_2_0","4_3_0")

des2 <- des2[-which(des2[,31] %in% ne),]


des2[which(des2$SVA_Species %in% c("0_0_1","0_0_2","0_0_3")),8] <- "Orang_spec"
des2[which(des2$SVA_Species %in% c("0_1_0","0_2_0")),8] <- "Chimp_spec"
des2[which(des2$SVA_Species %in% c("0_1_1","0_2_2")),8] <- "Chimp_Orang"
des2[which(des2$SVA_Species %in% c("1_0_0","2_0_0","3_0_0")),8] <- "Human_spec"
des2[which(des2$SVA_Species %in% c("0_0_0")),8] <- "No_SVA"
des2[which(des2$SVA_Species %in% c("1_1_0","2_2_0","3_3_0")),8] <- "Human_Chimp"
des2[which(des2$SVA_Species %in% c("1_1_1")),8] <- "Human_Chimp_Orang"
des2[which(des2$SVA_Species %in% c("1_0_1")),8] <- "Human_Orang"


des2 <- des2[which(des2$var > quantile(des2$var, 0.9)),]


PCD <- prcomp(des2[,-c(1,2,27:31)], scale. = T)

pdf("test.pdf")
autoplot(PCD, data = des2, colour = "SVA_Species")
dev.off()





PCD <- prcomp(des2[,c(2:7)], scale. = T)

pdf("test.pdf")
autoplot(PCD, data = des2, colour = "SVA_Species", loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)
dev.off()


x <- "CD2AP"
#pdf(paste0("200618_",x,"_Cross-Species_Expr.pdf"))
plot <- ggplot(mDF[which(mDF$Gene == x),], aes(x=week, y=baseMean, group=species, color=species)) + 
  geom_line() +
  geom_point()+
  theme_classic() +
  ggtitle(paste(x,"variance filtered", sep = " ")) 

ggsave(paste0("200624_",x,"_Cross-Species_Expr.pdf"), plot = plot)
