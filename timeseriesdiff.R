#Fred White 01/10/2020
#Function and code for calculating difference in time series data

TSDiff <- function(x,y,z){
  
  
  x <- x[-which(duplicated(x$geneSymbol)),]
  y <- y[-which(duplicated(y$geneSymbol)),]
  
  x <- x[order(x$geneSymbol),]
  y <- y[order(y$geneSymbol),]
  #calculate the absolute difference between SVA containing genes and their ortholog
  diff <- abs(x[which(x$geneSymbol %in% z),3:8] 
                 - y[which(y$geneSymbol %in% z),3:8])
  #half the difference for the first and last time point
  diff[,1] <- diff[,1]/2
  diff[,6] <- diff[,6]/2
  diff <- rowSums(diff)
  diff <- sum(diff)
  return(diff)
  
  
  
}


SVADiff <- TSDiff(list[[1]],list[[4]],SVAg)


genes <- Genes[-which(Genes %in% SVAg)]

i <- NULL
l <- 0
IGs <- list()
DIFF <- c()
for(i in 1:9999){
  
  set.seed(i*3)
  
  g <- sample(genes,length(SVAg))

  DIFF[i] <- TSDiff(list[[1]],list[[4]],g)
  
  if(DIFF[i] > SVADiff){
    l <- l + 1 
    IGs[[l]] <- g
    
  }
  
  
}



rGOIS <- do.call(c,IGs)
rSUMMARY <- as.data.frame(table(rGOIS))

#1 id wildly different ZBED1



#SVADiff > quantile(DIFF,0.95)  #FALSE
#SVADiff > quantile(DIFF,0.86)  #TRUE
DIFF <- c(DIFF,SVADiff)

d <- density(DIFF)
pdf("2005050.pdf")
plot(d, main = "Difference Density Distribution")
abline(v = SVADiff, col="red", lwd=3, lty=2)
text(SVADiff + 100, 2000000, "SVA Group", srt = 0.2, col = "red")
dev.off()















x <- list[[1]]
y <- list[[4]]
#z <- hsr
#z <- z[which(z %in% list[[1]]$geneSymbol)]
x <- x[-which(duplicated(x$geneSymbol)),]
y <- y[-which(duplicated(y$geneSymbol)),]

x <- x[order(x$geneSymbol),]
y <- y[order(y$geneSymbol),]

#calculate the absolute difference between SVA contatining genes and their ortholog
sdiff <- cbind(x[which(x$geneSymbol %in% z),1],abs(x[which(x$geneSymbol %in% z),3:8] 
            - y[which(y$geneSymbol %in% z),3:8]))
#half the difference for the first and last time point
sdiff[,2] <- sdiff[,2]/2
sdiff[,7] <- sdiff[,7]/2
sdiff <- data.frame(sdiff[,1],rowSums(sdiff[,-1]), stringsAsFactors = F)
squant95 <- quantile(sdiff[,2],0.95)


srq95 <- sdiff[which(sdiff[,2] > squant95),]
srq95 <- srq95[order(srq95[,2]),]


x <- list[[1]]
y <- list[[4]]
#z <- hsr
#z <- z[which(z %in% list[[1]]$geneSymbol)]
x <- x[-which(duplicated(x$geneSymbol)),]
y <- y[-which(duplicated(y$geneSymbol)),]

x <- x[order(x$geneSymbol),]
y <- y[order(y$geneSymbol),]

#calculate the absolute difference between SVA contatining genes and their ortholog
sdiff <- cbind(x[which(x$geneSymbol %in% z),1],abs(x[which(x$geneSymbol %in% z),3:8] 
                                                   - y[which(y$geneSymbol %in% z),3:8]))
#half the difference for the first and last time point
sdiff[,2] <- sdiff[,2]/2
sdiff[,7] <- sdiff[,7]/2
sdiff <- data.frame(sdiff[,1],rowSums(sdiff[,-1]), stringsAsFactors = F)
squant95 <- quantile(sdiff[,2],0.95)


srq95 <- sdiff[which(sdiff[,2] > squant95),]
srq95 <- srq95[order(srq95[,2]),]


diff <- cbind(diff,c(rep("Non-SVA",dim(diff)[1])))
sdiff <- cbind(sdiff,c(rep("SVA",dim(sdiff)[1])))

colnames(diff) <- c("Gene","Difference","Class")
colnames(sdiff) <- c("Gene","Difference","Class")






adiff <- rbind(diff,sdiff)
adiff <- cbind(adiff,log2(adiff$Difference))
colnames(adiff)[4] <- "Log2Difference"

pdf("HRdensityplot.pdf")
ggplot(adiff, aes(x = Log2Difference, fill = Class)) + geom_density(alpha = 0.5)+
  scale_y_continuous(labels = percent_format()) + theme_bw()+
  ggtitle("Human-Rhesus Gene Expression Differences in Development")
dev.off()




sum(x[x$geneSymbol %in% z,1] == y[y$geneSymbol %in% z,1])



x <- list[[1]]
y <- list[[4]]
z <- hsr
z <- z[which(z %in% list[[1]]$geneSymbol)]
x <- x[which(x$geneSymbol %in% z),]
y <- y[which(y$geneSymbol %in% z),]
x <- x[order(x$geneSymbol),]
y <- y[order(y$geneSymbol),]

dim(y)
dim(x)
head(x)
sum(x[,1] == y[,1])
x

f <- x[which(x$geneSymbol %in% rq95$z),]
g <- y[which(y$geneSymbol %in% rq95$z),]
abs(f[,3:8] - g[,3:8])


sdiff <- cbind(x[,1],(abs(x[,3:8] - y[,3:8])))





###################################
x <- list[[1]]
y <- list[[2]]

#z <- cg[which(cg %in% list[[1]]$geneSymbol)]
y <- y[which(y[,1] %in% x[,1]),]
x <- x[-which(duplicated(x$geneSymbol)),]
y <- y[-which(duplicated(y$geneSymbol)),]

x <- x[order(x$geneSymbol),]
y <- y[order(y$geneSymbol),]

#calculate the absolute difference between SVA contatining genes and their ortholog
diff <- cbind(x[,1],abs(x[,3:8] 
            - y[,3:8]))
#half the difference for the first and last time point
diff[,2] <- diff[,2]/2
diff[,7] <- diff[,7]/2


diff <- data.frame(diff[,1],rowSums(diff[,-1]), stringsAsFactors = F)
diff <- cbind(diff,c(rep(0,dim(diff)[1])))
colnames(diff) <- c("Gene","Difference","Class")


#diff[which(diff[,1] %in% z),3] <- "SVAs"

diff[which(diff[,1] %in% cg),3] <- "Common SVAs"
diff[which(diff[,1] %in% cs),3] <- "Chimp SVAs"
diff[which(diff[,1] %in% hs),3] <- "Human SVAs"


diff[which(diff[,3] == 0),3] <- "Non SVAs"


diff <- cbind(diff,log2(diff$Difference))
colnames(diff)[4] <- "Log2Difference"


genes <- diff$Gene
genes <- human2chimp(genes)
genes <- genes[-which(genes$humanGene %in% HISC[,2]),]
genes <- genes[-which(genes$chimpGene %in% ISC[,2]),]
genes <- genes[!duplicated(genes$humanGene),]

SVAg <- c(hs,cs,cg)


####downsample to get balanced group sizes cg is smallest at 145 genes present in current dataset
#ns <- sample(genes[,1],145)
#hs1 <- sample(hs[which(hs %in% diff$Gene)], 145)
#cs1 <- sample(cs[which(cs %in% diff$Gene)], 145)

#genes <- c(ns,hs1,cs1,cg[which(cg %in% diff$Gene)])




genes <- c(SVAg, genes$humanGene)


xdiff <- diff[which(diff$Gene %in% genes),]





pdf("HCdensityplot.pdf")
ggplot(xdiff, aes(x = Log2Difference, fill = Class, linetype = Class)) + geom_density(alpha = 0.3)+
  scale_y_continuous(labels = percent_format()) + theme_bw()+
  ggtitle("Human-Chimp Gene Expression Differences in Development")
dev.off()







quant95 <- quantile(diff,0.95)
z[which(diff > quant95)]



hscq95 <- z[which(diff > quant95)]


#rerun block above change z
csq95 <- z[which(diff > quant95)]

#rerun block above change z
cgq95 <- z[which(diff > quant95)]










