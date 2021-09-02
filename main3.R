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
library(AnnotationHub)
BiocManager::install("org.Mm.eg.db")

hub <- AnnotationHub()
q <- query(hub, "Mus")

# Format counts dataframe obtained from featurecounts
countData <- read.csv("PKoszalka_star_featurecounts.txt", sep = "\t", header = F, skip=1)
colnames(countData) <- countData[1,]
rownames(countData) <- countData[,1]
countData <- countData[-1,]
countData <- countData[,-c(1:6)]

# If star file
colnames(countData) <- gsub(".*Star_mapped/|Aligned.sortedByCoord.out.bam", "",
                            colnames(countData))
col_order <- c("KO1-1", "KO1-2", "KO1-3", "KO1-4", "KO2-1", "KO2-3", "KO2-4",
               "KO3-1", "KO3-2", "KO3-3A", "KO3-3D", "KO3-4", "KO3-5", "KO6-1",
               "KO6-2", "KO6-3", "KO6-4", "KO6-6", "KO6-7", "KO9-5", "KO9-6",
               "KO13-1", "KO13-3", "KO14-1", "KO14-3", "KO15-4", "KO15-5", 
               "KO14-5", "KO14-6", "KO15-1", "KO15-2", "KO15-3",
               "WT1-1", "WT1-4", "WT2-1", "WT2-2", "WT2-4", "WT3-2", "WT4-2", 
               "WT4-5A", "WT4-5B", "WT4-6", "WT5-2", "WT5-4", "WT5-5", "WT5-6",
               "WT7-1", "WT7-2", "WT7-3", "WT7-4A", "WT7-4B", "WT7-5", "WT7-6",
               "WT8-3", "WT8-4", "WT8-5", "WT8-6", "WT11-1", "WT11-5", "WT11-6",
               "WT12-1", "WT12-2")

countdata <- countData[, col_order]
countdata = data.frame(lapply(countdata, function(x) as.numeric(as.character(x))),
                       check.names=F, row.names = rownames(countdata))

type <- factor(c(rep("KO", 32), rep("WT", 30)))
stage <- factor(c(rep("progression", 27), rep("initiation", 5), rep("progression", 26), rep("initiation", 4)))
bc <- factor(c("mbcs", "mbcs", "sbc", "sbc", "mbcs", "mbcs", "sbc", "mbcs", "mbcs", "mbcs", "mbcs",
               "sbc", "mbcs", "sbc", "mbcs", "sbc", "sbc", "sbc", "sbc", "mbcs", "mbcs", "sbc",
               "sbc", "sbc", "mbcs", "sbc", "mbcs", "none", "none", "none", "none", "none",
               "sbc", "mbcs", "sbc", "mbcs", "mbcs", "sbc", "sbc", "mbcs", "mbcs", "sbc", "sbc",
               "sbc", "sbc", "sbc", "sbc", "mbcs", "mbcs", "mbcs", "mbcs", "sbc", "sbc", "sbc",
               "sbc", "mbcs", "sbc", "mbcs", "none", "none", "none", "none"))
condition <- factor(c(rep("KO_prog", 27), rep("KO_init", 5), rep("WT_prog", 26), rep("WT_init", 4)))

coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~condition)

# Filter low counts
keep <- rowSums(counts(dds) > 0) >= 2
dds <- dds[keep,]
dds$condiiton <- relevel(dds$condition, ref = "WT_init")
dds.res <- DESeq(dds)
res <- results(dds.res)
res.wt_init.ko_init <- results(dds.res, contrast=c("condition","KO_init","WT_init"))
res.wt_prog.ko_prog <- results(dds.res, contrast=c("condition","KO_prog","WT_prog"))
res.wt_init.wt_prog <- results(dds.res, contrast=c("condition","WT_prog","WT_init"))
res.ko_init.ko_prog <- results(dds.res, contrast=c("condition","KO_prog","KO_init"))


# Hist + volcano
par(mfrow=c(1,2))

hist(res.wt_init.ko_init$pvalue[res.wt_init.ko_init$baseMean > 1], breaks = 0:50/50,
     col = "sienna2", border = "white", main = NULL, xlab = "Adjusted p-value WT initiation vs KO initiation")
# Volcano plot
wt_init.ko_init.vol <- data.frame(log2FC = res.wt_init.ko_init$log2FoldChange, PVal = -log10(res.wt_init.ko_init$padj),
                      row.names = rownames(res.wt_init.ko_init))

plot(wt_init.ko_init.vol, col=rgb(128, 128, 128, max = 255, alpha = 200), 
     bg=rgb(128, 128, 128, max = 255, alpha = 125), pch=21, cex = 0.7, xlim=c(-8,8),
     xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), ylim = c(0, 6))
sigGenes = (abs(wt_init.ko_init.vol$log2FC) > 0.58 & wt_init.ko_init.vol$PVal > -log10(0.05))
points(wt_init.ko_init.vol[sigGenes, ], col="black", bg="orange", pch=21, cex = 0.8)
abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-0.58, 0.58), col = "black", lty = 2)
mtext(paste("pval =", 0.05), side = 4, at = -log10(0.05), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", 1.5, "fold"), paste("+", 1.5, "fold")), side = 3, at = c(-1, 1),
      cex = 0.8, line = 0.5)

par(mfrow=c(1,2))
hist(res.wt_prog.ko_prog$pvalue[res.wt_prog.ko_prog$baseMean > 1], breaks = 0:50/50,
     col = "sienna2", border = "white", main = NULL, xlab = "Adjusted p-value WT progression vs KO progression")
# Volcano plot
wt_prog.ko_prog.vol <- data.frame(log2FC = res.wt_prog.ko_prog$log2FoldChange, PVal = -log10(res.wt_prog.ko_prog$padj),
                                  row.names = rownames(res.wt_prog.ko_prog))

plot(wt_prog.ko_prog.vol, col=rgb(128, 128, 128, max = 255, alpha = 200), 
     bg=rgb(128, 128, 128, max = 255, alpha = 125), pch=21, cex = 0.7, xlim=c(-8,8),
     xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), ylim = c(0, 9))
sigGenes = (abs(wt_prog.ko_prog.vol$log2FC) > 0.58 & wt_prog.ko_prog.vol$PVal > -log10(0.05))
points(wt_prog.ko_prog.vol[sigGenes, ], col="black", bg="orange", pch=21, cex = 0.8)
abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-0.58, 0.58), col = "black", lty = 2)
mtext(paste("pval =", 0.05), side = 4, at = -log10(0.05), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", 1.5, "fold"), paste("+", 1.5, "fold")), side = 3, at = c(-1, 1),
      cex = 0.8, line = 0.5)

par(mfrow=c(1,2))
hist(res.wt_init.wt_prog$pvalue[res.wt_init.wt_prog$baseMean > 1], breaks = 0:50/50,
     col = "sienna2", border = "white", main = NULL, xlab = "Adjusted p-value WT initiation vs WT progression")
# Volcano plot
wt_init.wt_prog.vol <- data.frame(log2FC = res.wt_init.wt_prog$log2FoldChange, PVal = -log10(res.wt_init.wt_prog$padj),
                                  row.names = rownames(res.wt_init.wt_prog))

plot(wt_init.wt_prog.vol, col=rgb(128, 128, 128, max = 255, alpha = 200), 
     bg=rgb(128, 128, 128, max = 255, alpha = 125), pch=21, cex = 0.7, xlim=c(-10,10),
     xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), ylim = c(0, 24))
sigGenes = (abs(wt_init.wt_prog.vol$log2FC) > 0.58 & wt_init.wt_prog.vol$PVal > -log10(0.05))
points(wt_init.wt_prog.vol[sigGenes, ], col="black", bg="orange", pch=21, cex = 0.8)
abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-0.58, 0.58), col = "black", lty = 2)
mtext(paste("pval =", 0.05), side = 4, at = -log10(0.05), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", 1.5, "fold"), paste("+", 1.5, "fold")), side = 3, at = c(-1, 1),
      cex = 0.8, line = 0.5)

par(mfrow=c(1,2))
hist(res.ko_init.ko_prog$pvalue[res.ko_init.ko_prog$baseMean > 1], breaks = 0:50/50,
     col = "sienna2", border = "white", main = NULL, xlab = "Adjusted p-value KO initiation vs KO progression")
# Volcano plot
ko_init.ko_prog.vol <- data.frame(log2FC = res.ko_init.ko_prog$log2FoldChange, PVal = -log10(res.ko_init.ko_prog$padj),
                                  row.names = rownames(res.ko_init.ko_prog))

plot(ko_init.ko_prog.vol, col=rgb(128, 128, 128, max = 255, alpha = 200), 
     bg=rgb(128, 128, 128, max = 255, alpha = 125), pch=21, cex = 0.7, xlim=c(-12,12),
     xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), ylim = c(0, 28))
sigGenes = (abs(ko_init.ko_prog.vol$log2FC) > 0.58 & ko_init.ko_prog.vol$PVal > -log10(0.05))
points(ko_init.ko_prog.vol[sigGenes, ], col="black", bg="orange", pch=21, cex = 0.8)
abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-0.58, 0.58), col = "black", lty = 2)
mtext(paste("pval =", 0.05), side = 4, at = -log10(0.05), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", 1.5, "fold"), paste("+", 1.5, "fold")), side = 3, at = c(-1, 1),
      cex = 0.8, line = 0.5)
#############################

res.wt_init.ko_init.sort <- res.wt_init.ko_init[order(res.wt_init.ko_init$pvalue),]
res.wt_prog.ko_prog.sort <- res.wt_prog.ko_prog[order(res.wt_prog.ko_prog$pvalue),]
res.wt_init.wt_prog.sort <- res.wt_init.wt_prog[order(res.wt_init.wt_prog$pvalue),]
res.ko_init.ko_prog.sort <- res.ko_init.ko_prog[order(res.ko_init.ko_prog$pvalue),]

res.wt_init.ko_init.pval <- as.data.frame(res.wt_init.ko_init) %>% filter(padj < 0.05)
res.wt_prog.ko_prog.pval <- as.data.frame(res.wt_prog.ko_prog) %>% filter(padj < 0.05)
res.wt_init.wt_prog.pval <- as.data.frame(res.wt_init.wt_prog) %>% filter(padj < 0.05)
res.ko_init.ko_prog.pval <- as.data.frame(res.ko_init.ko_prog) %>% filter(padj < 0.05)

res.wt_init.ko_init.pval.sort <- res.wt_init.ko_init.pval[order(res.wt_init.ko_init.pval$pvalue),]
res.wt_prog.ko_prog.pval.sort <- res.wt_prog.ko_prog.pval[order(res.wt_prog.ko_prog.pval$pvalue),]
res.wt_init.wt_prog.pval.sort <- res.wt_init.wt_prog.pval[order(res.wt_init.wt_prog.pval$pvalue),]
res.ko_init.ko_prog.pval.sort <- res.ko_init.ko_prog.pval[order(res.ko_init.ko_prog.pval$pvalue),]

res.wt_init.ko_init.fc <- as.data.frame(res.wt_init.ko_init.pval) %>% filter(abs(log2FoldChange) >= 1 & baseMean >= 5)
res.wt_prog.ko_prog.fc <- as.data.frame(res.wt_prog.ko_prog.pval) %>% filter(abs(log2FoldChange) >= 1 & baseMean >= 5)
res.wt_init.wt_prog.fc <- as.data.frame(res.wt_init.wt_prog.pval) %>% filter(abs(log2FoldChange) >= 1 & baseMean >= 5)
res.ko_init.ko_prog.fc <- as.data.frame(res.ko_init.ko_prog.pval) %>% filter(abs(log2FoldChange) >= 1 & baseMean >= 5)

res.wt_init.ko_init.names <- row.names(res.wt_init.ko_init.fc)
res.wt_prog.ko_prog.names <- row.names(res.wt_prog.ko_prog.fc)
res.wt_init.wt_prog.names <- row.names(res.wt_init.wt_prog.fc)
res.ko_init.ko_prog.names <- row.names(res.ko_init.ko_prog.fc)

# Create dataframe of DEGs with gene symbols and ensembl id
library(biomaRt)
add_gene_symbols <- function(df){
  symb <- numeric(0)
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  symb <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name"),
                   values=rownames(df),mart= mart)
  rownames(symb) <- symb[,1]
  temp.list <- merge(symb, df, by=0, all=TRUE)
  temp.list <- temp.list[-1]
  return(temp.list)
}

test.wt_init.wt_prog <- add_gene_symbols(as.data.frame(res.wt_init.ko_init.sort))
wt_init.wt_prog.filtred.DEGs <- add_gene_symbols(res.wt_init.wt_prog.fc)
ko_init.ko_prog.filtred.DEGs <- add_gene_symbols(res.ko_init.ko_prog.fc)
wt_prog.ko_prog.filtred.DEGs <- add_gene_symbols(res.wt_prog.ko_prog.fc)
wt_init.ko_init.filtred.DEGs <- add_gene_symbols(res.wt_init.ko_init.fc)

# Save DEGs
write_DEGs <- function(df, name){
  write.table(df[order(df$pvalue),], paste(deparse(substitute(name)), ".tsv", sep=""), sep = "\t", quote = F, row.names = F)
  write.xlsx(df[order(df$pvalue),], paste(deparse(substitute(name)), ".xls", sep=""), row.names = F)
}

# topGO function
topgo_ready <- function(df){ #use add_gene_symbols output
  df.vector <- df[,8]
  names(df.vector) <- df[,2]
  df.vector <- df.vector[!is.na(df.vector)]
  selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
  allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Mm.eg.db", ID="symbol")
  GOdata <- new("topGOdata",
                ontology="BP",
                allGenes=df.vector,
                annot=annFUN.GO2genes,
                GO2genes=allGO2genes,
                geneSel=selection,
                nodeSize=10)
  
  results.fisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
  goEnrichment <- GenTable(GOdata, fisher=results.fisher, orderBy="fisher", topNodes=500, numChar=1000)
  goEnrichment$fisher <- as.numeric(goEnrichment$fisher)
  goEnrichment$fisher[is.na(goEnrichment$fisher)] <- as.numeric(1.1e-30)
  goEnrichment <- goEnrichment[goEnrichment$fisher<0.05,]
  annotGenes = lapply(goEnrichment$GO.ID, function(x) as.character(unlist(genesInTerm(object = GOdata, whichGO = x))))
  z <- as.data.frame(df) %>% filter(padj < 0.05)
  sigGenes = lapply(annotGenes, function(x) intersect(x, z[,2])) # where INT.GENES$V1 is your list of interesting genes
  goEnrichment <- goEnrichment[,c("GO.ID","Term","fisher")]
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
  temp.df <- goEnrichment
  temp.df$Identified_Genes <- vapply(sigGenes, paste, collapse = ", ", character(1L))
  return(temp.df)
}

wt_init.ko_init.GO <- topgo_ready(add_gene_symbols(as.data.frame(res.wt_init.ko_init.sort)))
wt_prog.ko_prog.GO <- topgo_ready(add_gene_symbols(as.data.frame(res.wt_prog.ko_prog.sort)))
wt_init.wt_prog.GO <- topgo_ready(add_gene_symbols(as.data.frame(res.wt_init.wt_prog.sort)))
ko_init.ko_prog.GO <- topgo_ready(add_gene_symbols(as.data.frame(res.ko_init.ko_prog.sort)))

# Vis topGO results
ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment") +
  ggtitle("Title") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
  theme_bw() +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=9, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=9, face="bold", vjust=0.5),
    axis.title=element_text(size=12, face="bold"),
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=9),  #Text size
    title=element_text(size=9)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()

##################

# PCA
rld <- rlogTransformation(dds.res, blind = T)

pca.t <- t(assay(rld))
pca.t <- cbind(pca.t, group = c(rep("KO_prog", 27), rep("KO_init", 5), rep("WT_prog", 26), rep("WT_init", 4)))
pca.t <- as.data.frame(pca.t)
pca.t[seq(1,(dim(pca.t)[2])-1)] <- sapply(pca.t[seq(1,(dim(pca.t)[2])-1)],as.numeric)
init.pca.t <- pca.t[-c(28:32, 59:62),]
init2.pca.t <- pca.t[c(28:32, 59:62),]
pca.t <- init2.pca.t


pca.t.class <- prcomp(pca.t[,-dim(pca.t)[2]])
plot(pca.t.class$x[,1], pca.t.class$x[,2])

init.pca.t <- pca.t[-c(28:32, 59:62),]
pca.t <- init.pca.t

pca.df <- as.data.frame(pca.t.class$x) # All gene set
pca.df$group <- condition[c(28:32, 59:62)]

pca.df.pc1.500 <- order(abs(pca.t.class$rotation[,1]), decreasing=TRUE)[1:500]
pca.df.pc2.500 <- order(abs(pca.t.class$rotation[,2]), decreasing=TRUE)[1:500]
pca.t.class$sdev^2/sum(pca.t.class$sdev^2) # calculate variation

p<-ggplot(pca.df,
          aes(x=PC1,y=PC2, label=row.names(pca.df)))
p + geom_point(shape=21, size = 3, colour="black", aes(fill=group)) + 
  theme_classic() +
  geom_text_repel(size=4) +
  scale_fill_manual(values = wes_palette("Moonrise2")) #Darjeeling1

# UpSetR
library(UpSetR)
listInput <- list(Initiation=sort(wt_init.ko_init.filtred.DEGs[,1]), 
                  Progression=sort(wt_prog.ko_prog.filtred.DEGs[,1]), 
                  WT=sort(wt_init.wt_prog.filtred.DEGs[,1]), 
                  KO=sort(ko_init.ko_prog.filtred.DEGs[,1]))
upset(fromList(listInput), order.by = "freq", point.size = 3.5, line.size = 1.5,
      mainbar.y.label = "Group Intersections", sets.x.label = "Genes per group", text.scale = 1.5)

# CD73 boxplots
wt.prog <- c("WT1-1", "WT1-4", "WT2-1", "WT2-2", "WT2-4", "WT3-2", "WT4-2",
"WT4-5A", "WT4-5B", "WT4-6", "WT5-2", "WT5-4", "WT5-5", "WT5-6",
"WT7-1", "WT7-2", "WT7-3", "WT7-4A", "WT7-4B", "WT7-5", "WT7-6",
"WT8-3", "WT8-4", "WT8-5", "WT8-6", "WT11-1")
wt.init <- c("WT11-5", "WT11-6","WT12-1", "WT12-2")
wt.sbc <- c("WT1-1", "WT2-1", "WT3-2", "WT4-2", "WT4-6", "WT5-2", "WT5-4", "WT5-5",
            "WT5-6", "WT7-1", "WT7-5", "WT7-6", "WT8-3", "WT8-4", "WT8-6")
wt.mbcs <- c("WT1-4", "WT2-2", "WT2-4", "WT4-5A", "WT4-5B", "WT7-2", "WT7-3",
             "WT7-4A", "WT7-4B", "WT8-5", "WT11-1")
wt.mbcs.first <- c("WT1-4", "WT4-5A", "WT7-3", "WT7-4B", "WT11-1")
wt.mbcs.next <- c("WT2-2", "WT2-4", "WT4-5B", "WT7-2", "WT7-4A", "WT8-5")

wt.prog.col <- sel.gene[,wt.prog]
wt.init.col <- sel.gene[,wt.init]
wt.sbc.col <- sel.gene[,wt.sbc]
wt.mbcs.col <- sel.gene[,wt.mbcs]
wt.mbcs.first.col <- sel.gene[,wt.mbcs.first]
wt.mbcs.next.col <- sel.gene[,wt.mbcs.next]

wt.bind <- cbind(wt.prog.col, wt.init.col, wt.sbc.col, wt.mbcs.col, wt.mbcs.first.col, wt.mbcs.next.col)

# Boxplot of individual genes
library(wesanderson)
library(RColorBrewer)
dds.norm <- counts(dds.res, normalized = T)
sel.gene.name <- "ENSMUSG00000100685"
sel.gene <- as.data.frame(dds.norm) %>% filter(rownames(dds.norm) == sel.gene.name)
sel.gene.t <- data.frame(val=t(wt.bind)[,1], id=seq(1,67,1), group=c(rep("WT prog", 26), rep("WT init", 4), rep("sBC", 15), rep("mBCs", 11), rep("mBCs.one", 5), rep("mBCs.next", 6)))

sel.gene.t <- data.frame(val=t(sel.gene)[,1], id=rownames(t(sel.gene)), group=condition)
#sel.gene.t.excl <- sel.gene.t[!grepl("KO_init|WT_init", sel.gene.t$group),]

sel.gene.t %>% 
  ggplot(aes(x=group, y=val, label=id)) +
  ggtitle("CD73 (NT5E) gene expression") +
  stat_boxplot(geom = 'errorbar', width = .2) +
  geom_boxplot(width=.4, aes(fill=group), outlier.shape = NA) +
  #  geom_text(aes(label = id), na.rm = F, hjust = -0.3) +
  theme_classic() +
  scale_fill_brewer(palette = "Pastel2") +
  #  geom_point(size = 1.5) +
  geom_jitter(width = 0.1, size = 1.5) +
  ylim(-0.1,90) +
  ylab("Normalized read counts") + xlab("Group") +
  geom_hline(yintercept=0, linetype="dashed", color = "lightgrey") +
  theme(plot.title = element_text(hjust = 0.5))

####
sel.gene.name <- "ENSMUSG00000066693"
sel.gene <- as.data.frame(dds.norm) %>% filter(rownames(dds.norm) == sel.gene.name)
ko.prog <- c("KO1-1", "KO1-2", "KO1-3", "KO1-4", "KO2-1", "KO2-3", "KO2-4",
  "KO3-1", "KO3-2", "KO3-3A", "KO3-3D", "KO3-4", "KO3-5", "KO6-1",
  "KO6-2", "KO6-3", "KO6-4", "KO6-6", "KO6-7", "KO9-5", "KO9-6",
  "KO13-1", "KO13-3", "KO14-1", "KO14-3", "KO15-4", "KO15-5")
ko.prog.col <- sel.gene[,ko.prog]
wt.prog.col <- sel.gene[,wt.prog]
wt.ko.bind <- cbind(wt.prog.col, ko.prog.col)

sel.gene.t <- data.frame(val=t(wt.ko.bind)[,1], id=seq(1,53,1), group=c(rep("WT prog", 26), rep("KO prog", 27)))

#sel.gene.t <- data.frame(val=t(sel.gene)[,1], id=rownames(t(sel.gene)), group=condition)
#sel.gene.t.excl <- sel.gene.t[!grepl("KO_init|WT_init", sel.gene.t$group),]

sel.gene.t %>% 
  ggplot(aes(x=group, y=val, label=id)) +
  ggtitle("ErbB-2 (HER2) gene expression") +
  stat_boxplot(geom = 'errorbar', width = .2) +
  geom_boxplot(width=.4, aes(fill=group), outlier.shape = NA) +
  #  geom_text(aes(label = id), na.rm = F, hjust = -0.3) +
  theme_classic() +
  scale_fill_brewer(palette = "Pastel2") +
  #  geom_point(size = 1.5) +
  geom_jitter(width = 0.1, size = 1.5) +
  ylim(0,7500) +
  theme(plot.title = element_text(hjust = 0.5))
