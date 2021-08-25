library(dplyr)
library(DESeq2)
library(tximport)
library(readr)
library(stringr)
library(Rsubread)
library(DEGreport)
library(GenomicFeatures)
library(tximeta)
library(ggplot2)
library(pheatmap)
library(FactoMineR)
library(factoextra)
library(ggrepel)

# Format counts dataframe obtained from featurecounts
countData <- read.csv("sus_star_featurecounts.txt", sep = "\t", header = F, skip=1)
colnames(countData) <- countData[1,]
rownames(countData) <- countData[,1]
countData <- countData[-1,]
countData <- countData[,-c(1:6)]

# If star file
colnames(countData) <- gsub(".*Star_mapped/|Aligned.sortedByCoord.out.bam", "",
                            colnames(countData))

col_order <- c("G10-64_ctrl", "G10-68_ctrl", "G10-71_ctrl", "G10-73_ctrl", "G10-75_ctrl", "G10-78_ctrl",
               "G10-64_EVs", "G10-68_EVs", "G10-71_EVs", "G10-73_EVs", "G10-75_EVs", "G10-78_EVs")

countdata <- countData[, col_order]
countdata = data.frame(lapply(countdata, function(x) as.numeric(as.character(x))),
                       check.names=F, row.names = rownames(countdata))

condition <- factor(c(rep("ctrl", 6), rep("EVs", 6)))
id <- rep(factor(seq(1,6,1)),2)
coldata <- data.frame(row.names=colnames(countdata), condition, id)

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~id + condition)
#dds <- DESeq(dds)

# Filter low counts
keep <- rowSums(counts(dds) > 0) >= 4
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "ctrl")
dds.res <- DESeq(dds)
res <- results(dds.res)
res.ctrl.evs <- results(dds.res, contrast=c("condition","EVs","ctrl"))
hist(res.ctrl.evs$pvalue[res.ctrl.evs$baseMean > 1], breaks = 0:50/50,
  col = "sienna2", border = "white", main = NULL, xlab = "Adjusted p-value")


res.evs.sort <- res.ctrl.evs[order(res.ctrl.evs$pvalue),]
evs <- as.data.frame(res.evs.sort) %>% filter(padj < 0.05)
evs.fc <- evs %>% filter(abs(log2FoldChange) >= 0.58)

# Heatmap all
rld <- rlogTransformation(dds.res, blind = T)
topVarianceGenes <- head(order(rowVars(assay(rld)), decreasing=T),1000) # Select genes
matrix <- assay(rld)[topVarianceGenes,]
matrix <- matrix - rowMeans(matrix)

annotation_data <- as.data.frame(colData(rld)) # Phyloseq z-score
new_annotation_data <- annotation_data[,-3]
annotation_colors = list(
  'condition' = c('ctrl' = '#F2AD00', 'EVs' = '#B40F20'), 
  'id' = c('1' = '#9986A5', '2' = '#79402E', '3' = '#CCBA72', '4' = '#0F0D0E',
            '5' = '#D9D0D3', '6' = '#5BBCD6'))

pheatmap(matrix[1:40,], annotation_col = new_annotation_data, 
         scale = "row", annotation_colors = annotation_colors, 
         border_color = 'white') # Add cutree_cols=6 to add heatmap cutting

# Heatmap sig.gene
rld.small <- assay(rld[rownames(evs.fc),])
rld.mat <- rld.small - rowMeans(rld.small)
annotation_data <- as.data.frame(colData(rld)) # Phyloseq z-score
new_annotation_data <- annotation_data[,-3]
annotation_colors = list(
  'condition' = c('ctrl' = '#F2AD00', 'EVs' = '#B40F20'), 
  'id' = c('1' = '#9986A5', '2' = '#79402E', '3' = '#CCBA72', '4' = '#0F0D0E',
           '5' = '#D9D0D3', '6' = '#5BBCD6'))

pheatmap(rld.mat, annotation_col = new_annotation_data, 
         scale = "row", annotation_colors = annotation_colors, 
         border_color = 'white') # Add cutree_cols=6 to add heatmap cutting

# Volcano plot
evs.vol <- data.frame(log2FC = res.evs.sort$log2FoldChange, PVal = -log10(res.evs.sort$padj),
                         row.names = rownames(res.evs.sort))

plot(evs.vol, col=rgb(128, 128, 128, max = 255, alpha = 200), 
     bg=rgb(128, 128, 128, max = 255, alpha = 125), pch=21, cex = 0.7, xlim=c(-3.2,3.2),
     xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), ylim = c(0, 6))
sigGenes = (abs(evs.vol$log2FC) > 0.58 & evs.vol$PVal > -log10(0.05))
#sigGenes2 = (abs(pgrmc1.vol$log2FC) > 3 & pgrmc1.vol$PVal > -log10(0.00005))
points(evs.vol[sigGenes, ], col="black", bg="orange", pch=21, cex = 0.8)
#points(pgrmc1.vol[sigGenes2, ], col="black", bg="red", pch=21, cex = 0.8)
#text(pgrmc1.vol[sigGenes2, ], labels=rownames(pgrmc1.vol), cex=0.5, line = 0.5)

abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-0.58, 0.58), col = "black", lty = 2)
mtext(paste("pval =", 0.05), side = 4, at = -log10(0.05), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", 1.5, "fold"), paste("+", 1.5, "fold")), side = 3, at = c(-1, 1),
      cex = 0.8, line = 0.5)

# Boxplot of individual genes
dds.norm <- counts(dds.res, normalized = T)
sel.gene.name <- "ENSSSCG00000027525"
sel.gene <- as.data.frame(dds.norm) %>% filter(rownames(dds.norm) == sel.gene.name)
sel.gene.t <- data.frame(val=t(sel.gene)[,1], id=rownames(t(sel.gene)), group=condition)

sel.gene.t %>% 
  ggplot(aes(x=group, y=val, label=id)) +
  stat_boxplot(geom = 'errorbar', width = .2) +
  geom_boxplot(width=.4, aes(fill=group)) +
  #  geom_text(aes(label = id), na.rm = F, hjust = -0.3) +
  theme_classic() +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  #  geom_point(size = 1.5) +
  geom_jitter(width = 0.1, size = 1.5) +
  ylim(0,2500)

# PCA
pca.t <- t(assay(rld))
pca.t <- cbind(pca.t, group=c(rep("ctrl", 6), rep("EVs", 6)))
pca.t <- as.data.frame(pca.t)
pca.t[seq(1,(dim(pca.t)[2])-1)] <- sapply(pca.t[seq(1,(dim(pca.t)[2])-1)],as.numeric)
pca.t.class <- prcomp(pca.t[,-dim(pca.t)[2]])
plot(pca.t.class$x[,1], pca.t.class$x[,2])

pca.df <- as.data.frame(pca.t.class$x) # All gene set
pca.df$group <- condition

pca.df.pc1.500 <- order(abs(pca.t.class$rotation[,1]), decreasing=TRUE)[1:500]
pca.df.pc2.500 <- order(abs(pca.t.class$rotation[,2]), decreasing=TRUE)[1:500]
pca.t.class$sdev^2/sum(pca.t.class$sdev^2) # calculate variation

p<-ggplot(pca.df,
          aes(x=PC1,y=PC2, label=row.names(pca.df)))
p + geom_point(shape=21, size = 3, colour="black", aes(fill=group)) + 
  theme_classic() +
  geom_text_repel(size=4) +
  scale_fill_manual(values = wes_palette("Darjeeling1"))

# Samples correlation heatmap
sample.dist <- dist(t(assay(rld)))
sample.dist.mat <- as.matrix(sample.dist)
rownames(sample.dist.mat) <- colnames(countdata)
colnames(sample.dist.mat) <- colnames(countdata)
pheatmap(sample.dist.mat, clustering_distance_rows=sample.dist, 
         clustering_distance_cols = sample.dist, border_color = 'white')
