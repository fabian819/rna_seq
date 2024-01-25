#Change for appropriate working directory
#setwd("path/to/working/directory")
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
library("AnnotationDbi")
library("org.Hs.eg.db")

#Load the data and create a DESeqDataSet
cts <- read.csv("CDS_counts_processed.txt",
                header = TRUE, skip=1,
                sep="\t", row.name="Geneid")

cts <- cts[, c(3,4,1,2)]

colnames(cts) <- c("WT_Rep1", "WT_Rep2", "KO_Rep1", "KO_Rep2")

coldata <- data.frame("condition"=c("WT", "WT", "KO", "KO"),
                      "run"= c("WT_1", "WT_2", "KO_1", "KO_2"))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# Perform pre-filtering
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#Run the DES function
dds <- DESeq(dds)

#Log data transform
rld <- rlog(dds, blind=FALSE)

# Plotting sample distances
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$run
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("QC_sample_distance_matrix_CDS.pdf", paper="a4")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

#Plotting heatmap of count matrix
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing = TRUE)
df <- as.data.frame(colData(dds)[,c("condition","run")])
pdf("QC_count_matrix_CDS.pdf", paper="a4")
pheatmap(assay(rld)[select,],
         cluster_rows=FALSE,
         show_rownames=FALSE,
         cluster_cols=FALSE,
         annotation_col=df)
dev.off()

#Plotting PCA
pdf("QC_PCA_CDS.pdf")
pcaData <- plotPCA(rld, intgroup=c("condition", "run"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=run)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

res <- results(dds, contrast = c("condition", "KO", "WT"), alpha=0.05)
summary(res)


convertIDs <- function( ids, from, to, db, ifMultiple = c("putNA", "useFirst")) {
  stopifnot(inherits(db, "AnnotationDb"))
  ifMultiple <- match.arg(ifMultiple)
  suppressWarnings(selRes <- AnnotationDbi::select(
    db, keys = ids, keytype = from, columns = c(from, to)))
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[duplicated(selRes[ , 1] ), 1]
    selRes <- selRes[!selRes[ , 1] %in% duplicatedIds, ]
  }
  return(selRes[match(ids, selRes[ , 1]), 2])
}

res$GeneID <- row.names(res)
res$hgnc_symbol <- convertIDs(row.names(res), "ENSEMBL", "SYMBOL", org.Hs.eg.db)

## Convert DESeq object to dataframe
res_df <- as.data.frame(res)

## Define differentially expressed genes as:
## Upregulated (LFC > 0.5 & padj < 0.05)
## Downregulated (LFC < - 0.5 & padj < 0.05)

res_df$regulation_level <- ifelse((res_df$log2FoldChange > 0.5 & res_df$padj < 0.05), "Upregulated",
                                         ifelse((res_df$log2FoldChange < - 0.5 & res_df$padj < 0.05), "Downregulated",
                                                "Unchanged"))

write.table (res_df,
             file = "DESeq2_results.csv",
             sep = ",",
             row.names = F,
             col.names = T,
             quote = F)

res_df$regulation_level <- factor(res_df$regulation_level, levels = c("Upregulated", "Downregulated", "Unchanged"))

## Remove rows with NA in padj column before plotting
res_df <- res_df[!is.na(res_df$padj), ]

pdf("DESeq2_Volcano_Plot.pdf", width = 4, height = 5)
ggplot(res_df,
       aes(x = -log10(padj),
           y = log2FoldChange,
           color = regulation_level)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("#23bf15", "#bf1545", "#545656")) +
  xlab("-log10(adjusted p-value)") +
  ylab("Log2 fold change") +
  labs(color = "Regulation level") +
  theme_bw()
dev.off()