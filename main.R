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

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pheatmap")
BiocManager::install("DEGreport")
BiocManager::install("Rsubread")
BiocManager::install("tximport")
BiocManager::install("tximeta")

# Format counts dataframe obtained from featurecounts
countData <- read.csv("star_featurecounts.txt", sep = "\t", header = F, skip=1)
colnames(countData) <- countData[1,]
rownames(countData) <- countData[,1]
countData <- countData[-1,]
countData <- countData[,-c(1:6)]

# If hisat2 file
colnames(countData) <- gsub(".bam|Mapped/", "", colnames(countData))

# If star file
colnames(countData) <- gsub(".*Star_mapped/|Aligned.sortedByCoord.out.bam", "", colnames(countData))

col_order <- c("S1", "S6", "S11", "S16", "S21", "S26",
               "S2", "S7", "S12", "S17", "S22", "S27",
               "S3", "S8", "S13", "S18", "S23", "S28",
               "S4", "S9", "S14", "S19", "S24", "S29",
               "S5", "S10", "S15", "S20", "S25", "S30")
countdata <- countData[, col_order]
countdata = data.frame(lapply(countdata, function(x) as.numeric(as.character(x))),
                       check.names=F, row.names = rownames(countdata))

condition <- factor(c(rep("ctl", 6), rep("pgrmc1", 6), rep("pgrmc2", 6), 
                      rep("mPRa", 6), rep("mPRb", 6)))

cow <- rep(factor(seq(1,6,1)),5)

coldata <- data.frame(row.names=colnames(countdata), condition, cow)

## Read in salmon data #########################################################
# Prepare samples dataframe
samples <- data.frame(rownames(coldata), coldata$condition)
colnames(samples) <- c("sample", "condition")
files <- file.path("quants", samples$sample, "quant.sf")
names(files) <- paste0(samples$sample)

# Create gene counts from transcript counts using gtf annotation file
txdb <- makeTxDbFromGFF("Bos_taurus.ARS-UCD1.2.104.chr.gtf", organism = 'Bos taurus')
k <- keys(txdb, keytype = "GENEID")
tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
new.order <- c("TXNAME", "GENEID")
tx2gene <- tx2gene[, new.order]
txi.salmon <- tximport(files, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = T)
dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~condition)
################################################################################

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~cow + condition)
dds <- DESeq(dds)

pd.set_option('display.max_colwidth', None)

# Filter low counts
keep <- rowSums(counts(dds) > 0) >= 5
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "ctl")
dds.res <- DESeq(dds)
res <- results(dds.res)

res.pgrmc1.pgrmc2 <- results(dds, contrast=c("condition","pgrmc2","pgrmc1"))
res.mPRa.mPRb <- results(dds, contrast=c("condition","mPRb","mPRa"))
res.mPRa.pgrmc1 <- results(dds, contrast=c("condition","pgrmc1","mPRa"))


res.pgrmc1 <- results(dds, contrast=c("condition","pgrmc1","ctl"))
res.pgrmc2 <- results(dds, contrast=c("condition","pgrmc2","ctl"))
res.mPRa <- results(dds, contrast=c("condition","mPRa","ctl"))
res.mPRb <- results(dds, contrast=c("condition","mPRb","ctl"))
res.pgrmc1.sort <- res.pgrmc1[order(res.pgrmc1$pvalue),]
res.pgrmc2.sort <- res.pgrmc2[order(res.pgrmc2$pvalue),]
res.mPRa.sort <- res.mPRa[order(res.mPRa$pvalue),]
res.mPRb.sort <- res.mPRb[order(res.mPRb$pvalue),]

pgrmc1 <- as.data.frame(res.pgrmc1.sort) %>% filter(padj < 0.05)
pgrmc2 <- as.data.frame(res.pgrmc2.sort) %>% filter(padj < 0.05)
mPRa <- as.data.frame(res.mPRa.sort) %>% filter(padj < 0.05)
mPRb <- as.data.frame(res.mPRb.sort) %>% filter(padj < 0.05)

pgrmc1.fc <- pgrmc1 %>% filter(abs(log2FoldChange) >= 1)
pgrmc2.fc <- pgrmc2 %>% filter(abs(log2FoldChange) >= 1)
mPRa.fc <- mPRa %>% filter(abs(log2FoldChange) >= 1)
mPRb.fc <- mPRb %>% filter(abs(log2FoldChange) >= 1)
pgrmc1.pgrmc2.fc <- pgrmc1.pgrmc2 %>% filter(abs(log2FoldChange) >= 1)

# Data check; intersect of groups; read counts of specified gene
length(intersect(rownames(pgrmc1), rownames(pgrmc2)))
as.data.frame(countdata) %>% filter(rownames(countdata) == 'ENSBTAG00000019552')

rld <- rlogTransformation(dds, blind = T)
vst <- vst(dds)
# P-value histogram QC
hist(res.pgrmc1$pvalue[res.pgrmc1$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

head(assay(rld))
hist(assay(rld))

## Graphical results
# Volcano 100%
pgrmc1.vol <- data.frame(log2FC = res.pgrmc1.sort$log2FoldChange, PVal = -log10(res.pgrmc1.sort$padj),
                         row.names = rownames(res.pgrmc1.sort))
par(mar = c(5, 4, 4, 4))
plot(pgrmc1.vol, col=rgb(128, 128, 128, max = 255, alpha = 200), 
     bg=rgb(128, 128, 128, max = 255, alpha = 125), pch=21, cex = 0.7, 
     xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), ylim = c(0, 18))
sigGenes = (abs(pgrmc1.vol$log2FC) > 1 & pgrmc1.vol$PVal > -log10(0.05))
#sigGenes2 = (abs(pgrmc1.vol$log2FC) > 3 & pgrmc1.vol$PVal > -log10(0.00005))
points(pgrmc1.vol[sigGenes, ], col="black", bg="orange", pch=21, cex = 0.8)
#points(pgrmc1.vol[sigGenes2, ], col="black", bg="red", pch=21, cex = 0.8)
#text(pgrmc1.vol[sigGenes2, ], labels=rownames(pgrmc1.vol), cex=0.5, line = 0.5)

abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-1, 1), col = "black", lty = 2)
mtext(paste("pval =", 0.05), side = 4, at = -log10(0.05), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", 2, "fold"), paste("+", 2, "fold")), side = 3, at = c(-1, 1),
      cex = 0.8, line = 0.5)

# Normalized counts
dds.norm <- counts(dds, normalized = T)
one.g0 <- as.data.frame(dds.norm) %>% filter(rownames(dds.norm) == 'ENSBTAG00000019552') # PGRMC1 gene
one.g1 <- as.data.frame(dds.norm) %>% filter(rownames(dds.norm) == 'ENSBTAG00000010843') # PGRMC2 gene
one.g2 <- as.data.frame(dds.norm) %>% filter(rownames(dds.norm) == 'ENSBTAG00000021787') # MPRA gene (PAQR7)
one.g3 <- as.data.frame(dds.norm) %>% filter(rownames(dds.norm) == 'ENSBTAG00000025494') # MPRB gene (PAQR8)


sel.gene.name <- "ENSBTAG00000006731"
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
  ylim(0,1500)

boxplot(data.frame(t(one.g3[1:6]), t(one.g3[7:12]), t(one.g3[13:18]), t(one.g3[19:24]), t(one.g3[25:30])))
stripchart(data.frame(t(one.g3[1:6]), t(one.g3[7:12]), t(one.g3[13:18]), t(one.g3[19:24]), t(one.g3[25:30])),
           method = "jitter", pch = 19, cex = 0.7, col = 2:4, vertical = TRUE,
           add = TRUE)


one.g0.t <- data.frame(val=t(one.g0)[,1], id=rownames(t(one.g0)), group=condition)
one.g1.t <- data.frame(val=t(one.g1)[,1], id=rownames(t(one.g1)), group=condition)
one.g2.t <- data.frame(val=t(one.g2)[,1], id=rownames(t(one.g2)), group=condition)
one.g3.t <- data.frame(val=t(one.g3)[,1], id=rownames(t(one.g3)), group=condition)


g2 <- one.g1.t %>% 
  ggplot(aes(x=group, y=val, label=id)) +
  stat_boxplot(geom = 'errorbar', width = .2) +
  geom_boxplot(width=.4, aes(fill=group)) +
#  geom_text(aes(label = id), na.rm = F, hjust = -0.3) +
  theme_classic() +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  geom_point(size = 1.5) 
#  stat_summary(fun=mean, geom="point", size=2) + 
#  stat_summary(fun.data = mean_se, geom = "errorbar")


# Plotting custom PCA
rld.pca <- prcomp(assay(rld))


library(datasets)
pca.t <- t(assay(rld))
pca.t <- cbind(pca.t, group=c(rep("ctl", 6), rep("pgrmc1", 6), rep("pgrmc2", 6), 
                              rep("mPRa", 6), rep("mPRb", 6)))
pca.t <- as.data.frame(pca.t)
pca.t[seq(1,18703)] <- sapply(pca.t[seq(1,18703)],as.numeric)
pca.t.class <- prcomp(pca.t[,-18704], graph = F)
plot(pca.t.class$x[,1], pca.t.class$x[,2])

pca.df <- as.data.frame(pca.t.class$x) # All gene set
pca.df$group <- condition

pca.df.pc1.500 <- order(abs(pca.t.class$rotation[,1]), decreasing=TRUE)[1:500]
pca.df.pc2.500 <- order(abs(pca.t.class$rotation[,2]), decreasing=TRUE)[1:500]
pca.t.class$sdev^2/sum(pca.t.class$sdev^2) # calculate variation

p<-ggplot(pca.df,
          aes(x=PC1,y=PC3, label=row.names(pca.df)))
p + geom_point(shape=21, size = 3, colour="black", aes(fill=group)) + 
  theme_classic() +
  geom_text_repel(size=4) +
  scale_fill_manual(values = wes_palette("Darjeeling1"))

# Heatmaps
rld <- rlog(dds, blind=F)
topVarianceGenes <- head(order(rowVars(assay(rld)), decreasing=T),40) # Select genes
matrix <- assay(rld)[topVarianceGenes,]
matrix <- matrix - rowMeans(matrix)

annotation_data <- as.data.frame(colData(rld)) # Phyloseq z-score
new_annotation_data <- annotation_data[,-3]
annotation_colors = list(
  'condition' = c('ctl' = '#FF0000', 'mPRa' = '#00A08A', 'mPRb' = '#F2AD00', 
                  'pgrmc1' = '#F98400', 'pgrmc2' = '#B40F20'), 
  'cow' = c('1' = '#9986A5', '2' = '#79402E', '3' = '#CCBA72', '4' = '#0F0D0E',
            '5' = '#D9D0D3', '6' = '#5BBCD6'))

pheatmap(matrix[1:30,], annotation_col = new_annotation_data, 
         scale = "row", annotation_colors = annotation_colors, 
         border_color = 'white') # Add cutree_cols=6 to add heatmap cutting

var.genes <- apply(dds.norm, 1, var)
select.var <- names(sort(var.genes, decreasing=TRUE))[1:50]
highly.variable.genecounts <- dds.norm[select.var,]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
heatmap.2(highly.variable.genecounts,col=rev(morecols(50)),trace="none", 
          main="Top 500 most variable genes across samples")

# Samples correlation heatmap
sample.dist <- dist(t(assay(rld)))
sample.dist.mat <- as.matrix(sample.dist)
rownames(sample.dist.mat) <- paste(rld$condition, rld$cow, sep="-")
colnames(sample.dist.mat) <- paste(rld$condition, rld$cow, sep="-")
pheatmap(sample.dist.mat, clustering_distance_rows=sample.dist, clusterin_distance_cols=sample.dist)

# UpSetR
listInput <- list(pgrmc1=pgrmc1.names, pgmrc2=pgrmc2.names, mPRa=mPRa.names, mPRb=mPRb.names)
upset(fromList(listInput), order.by = "freq", point.size = 3.5, line.size = 1.5,
      mainbar.y.label = "Group Intersections", sets.x.label = "Genes per group")

# Convert Ensembl ID to gene symbol
library('biomaRt')
mart <- useDataset("btaurus_gene_ensembl", useMart("ensembl"))
genes <- rownames(pgrmc1)
genes_list <- getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id", "external_gene_name"), 
                    values = genes, mart = mart)

df.test <- pgrmc1
df.test.gene.symbols <- inner_join(df.test %>% 
  mutate(ensembl_gene_id = rownames(df.test)), genes_list, by = "ensembl_gene_id")
