library(dplyr)
library(DESeq2)
library(tximport)
library(readr)
library(stringr)
library(Rsubread)
library(DEGreport)
library(GenomicFeatures)
library(tximeta)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

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
colnames(countData) <- gsub("Star_mapped/|Aligned.sortedByCoord.out.bam", "", colnames(countData))

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

coldata <- data.frame(row.names=colnames(countdata), condition)

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

dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds)

pd.set_option('display.max_colwidth', None)

# Filter low counts
keep <- rowSums(counts(dds) > 0) >= 5
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "ctl")
dds.res <- DESeq(dds)
res <- results(dds.res)

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

# Data check; intersect of groups; read counts of specified gene
length(intersect(rownames(pgrmc1), rownames(pgrmc2)))
as.data.frame(countdata) %>% filter(rownames(countdata) == 'ENSBTAG00000019552')

rld <- rlogTransformation(dds)
hist(res.pgrmc1$pvalue[res.pgrmc1$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

head(assay(rld))
hist(assay(rld))


# Plotting custom PCA
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}

png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()

