#Change for appropriate working directory
#setwd("path/to/working/directory")
library(grid)
library(gridExtra)
library(pathview)
library(AnnotationDbi)
library(clusterProfiler)
library(topGO)
library(ggplot2)
library(tidyr)

#Load data from DESeq2 differential expression analysis
#Data requires adjusted p values, ENSEMBL GeneIDs and log2-fold changes
df <- read.csv("DESeq2_results.csv")


#Get factor for up- and downregulated genes
upregulated <- ifelse(df$padj < 0.05 & df$log2FoldChange > 0, 1, 0)
downregulated <- ifelse(df$padj < 0.05 & df$log2FoldChange < 0, 1, 0)

genelist_up <- factor(upregulated)
names(genelist_up) <- df$GeneID

genelist_down <- factor(downregulated)
names(genelist_down) <- df$GeneID

upregulated <- df$GeneID[upregulated==1]
downregulated <- df$GeneID[downregulated==1]

#Map IDs
allGO2genes <- annFUN.org(whichOnto = "ALL",
                          feasibleGenes = NULL,
                          mapping = "org.Hs.eg.db",
                          ID = "ensembl")


#Wrapper function for creating TopGO object

topgo_object <- function(ontology, gene_list) {
  
  object <- new("topGOdata",
                ontology = ontology,
                allGenes = gene_list,
                annot = annFUN.GO2genes,
                GO2genes= allGO2genes,
                nodeSize = 10)
}

#Upregulated genes
GOdata_up_bp <- topgo_object("BP", genelist_up)
GOdata_up_mf <- topgo_object("MF", genelist_up)
GOdata_up_cc <- topgo_object("CC", genelist_up)

#Downregulated genes
GOdata_down_bp <- topgo_object("BP", genelist_down)
GOdata_down_mf <- topgo_object("MF", genelist_down)
GOdata_down_cc <- topgo_object("CC", genelist_down) 

#Fischer exact test
resultFis_up_bp <- runTest(GOdata_up_bp, algorithm = "elim", statistic = "fisher")
resultFis_up_mf <- runTest(GOdata_up_mf, algorithm = "elim", statistic = "fisher")
resultFis_up_cc <- runTest(GOdata_up_cc, algorithm = "elim", statistic = "fisher")
resultFis_down_bp <- runTest(GOdata_down_bp, algorithm = "elim", statistic = "fisher")
resultFis_down_mf <- runTest(GOdata_down_mf, algorithm = "elim", statistic = "fisher")
resultFis_down_cc <- runTest(GOdata_down_cc, algorithm = "elim", statistic = "fisher")

#Function to extract data from S4 objects
parse_tables <- function(GO_data, statistics) {
  
  #' Parse TopGO data
  
  goEnrichment <- GenTable(GO_data, weightFisher = statistics, topNodes = 20)
  sub("< ", "", goEnrichment$weightFisher)
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels = rev(goEnrichment$Term))
  goEnrichment$weightFisher <- as.numeric(sub("< ", "", goEnrichment$weightFisher))  
  goEnrichment
  
}

GOres_up_bp <- parse_tables(GOdata_up_bp, resultFis_up_bp)
GOres_up_mf <- parse_tables(GOdata_up_mf, resultFis_up_mf)
GOres_up_cc <- parse_tables(GOdata_up_cc, resultFis_up_cc)

GOres_down_bp <- parse_tables(GOdata_down_bp, resultFis_down_bp)
GOres_down_mf <- parse_tables(GOdata_down_mf, resultFis_down_mf)
GOres_down_cc <- parse_tables(GOdata_down_cc, resultFis_down_cc)

#Function for creating lollipop plots from top 20 tables
plot_GO <- function(GO_data, Ontology, Regulation, use_color) {
  

  GO_data$log_weightFisher <- (- log10(as.numeric(GO_data$weightFisher)))
  
  ggplot(GO_data, 
         aes(x = log_weightFisher,
             y = Term)) +
    geom_segment(aes(x = 0,
                     xend = log_weightFisher,
                     y = Term,
                     yend = Term),
                 colour = use_color)  +
    geom_point(aes(size = Significant),
               colour = use_color) +
    scale_size_area(name = "Gene counts") +
    xlab("Enrichment (- log10 Pvalue)") +
    ylab(Ontology) +
    ggtitle(Regulation) +
    scale_x_continuous() +
    theme_bw() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"))
}


#Make Plots
plot_up_BP <- plot_GO(GOres_up_bp, "Biological Proccess", "TopGO Up (fisher's exact test)", "#23bf15")
plot_up_MF <- plot_GO(GOres_up_mf, "Molecular Function", "TopGO Up (fisher's exact test)", "#23bf15")
plot_up_CC <- plot_GO(GOres_up_cc, "Cellular Component", "TopGO Up (fisher's exact test)", "#23bf15")


plot_down_BP <- plot_GO(GOres_down_bp, "Biological Proccess", "TopGO Down (fisher's exact test)", "#bf1545")
plot_down_MF <- plot_GO(GOres_down_mf, "Molecular Function", "TopGO Down (fisher's exact test)", "#bf1545")
plot_down_CC <- plot_GO(GOres_down_cc, "Cellular Component", "TopGO Down (fisher's exact test)", "#bf1545")


#Prepare data for writing
#retrieve genes2GO list from the "expanded" annotation in GOdata

extract_genes <- function(GOdata, regulation_level) {
  
  #Extract our genes in all identified GO terms
  
  all_GO_dereg <- genesInTerm(GOdata)
  data_to_extract <- lapply(all_GO_dereg, function(x) x[x %in% regulation_level])
  no_of_observations <- sapply(data_to_extract, length)
  seq_max <- seq_len(max(no_of_observations))
  GO_df <- t(sapply(data_to_extract, "[", i = seq_max))
  GO_df <- data.frame(GO_df[complete.cases(GO_df[ , 1:2]), ])
  GO_df$GO_ID <- row.names(GO_df)
  GO_df <- pivot_longer(GO_df,
                        cols = - GO_ID,
                        names_to = "Column_no",
                        values_to = "GeneID")
  GO_df <- GO_df[ , c(1, 3)]
  GO_df <- na.omit(GO_df)
  GO_df
  
}

GO_df_up_bp <- extract_genes(GOdata_up_bp, upregulated)
GO_df_up_mf <- extract_genes(GOdata_up_mf, upregulated)
GO_df_up_cc <- extract_genes(GOdata_up_cc, upregulated)

GO_df_down_bp <- extract_genes(GOdata_down_bp, downregulated)
GO_df_down_mf <- extract_genes(GOdata_down_mf, downregulated)
GO_df_down_cc <- extract_genes(GOdata_down_cc, downregulated)

## Export plots as PDF

sample_name = "RPF_KO_vs_WT"

pdf(paste(sample_name, "_Biological_Proccess_TopGO_Up_fisher.pdf", sep = ""))
plot_up_BP
dev.off()

pdf(paste(sample_name, "_Molecular_Function_TopGO_Up_fisher.pdf", sep = ""))
plot_up_MF
dev.off()

pdf(paste(sample_name, "_Cellular_Component_TopGO_Up_fisher.pdf", sep = ""))
plot_up_CC
dev.off()

pdf(paste(sample_name, "_Biological_Proccess_TopGO_down_fisher.pdf", sep = ""))
plot_down_BP
dev.off()

pdf(paste(sample_name, "_Molecular_Function_TopGO_down_fisher.pdf", sep = ""))
plot_down_MF
dev.off()

pdf(paste(sample_name, "_Cellular_Component_TopGO_down_fisher.pdf", sep = ""))
plot_down_CC
dev.off()

# Tables
write.table(GO_df_up_bp,
            file = paste(sample_name, "_BP_Up_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GO_df_up_mf,
            file = paste(sample_name, "_MF_Up_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GO_df_up_cc,
            file = paste(sample_name, "_CC_Up_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GO_df_down_bp,
            file = paste(sample_name, "_BP_Down_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GO_df_down_mf,
            file = paste(sample_name, "_MF_Down_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GO_df_down_cc,
            file = paste(sample_name, "_CC_Down_fisher_Genes_in_GO.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

## GO data

write.table(GOres_up_bp,
            file = paste(sample_name, "_BP_Up_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GOres_up_mf,
            file = paste(sample_name, "_MF_Up_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GOres_up_cc,
            file = paste(sample_name, "_CC_Up_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GOres_down_bp,
            file = paste(sample_name, "_BP_Down_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GOres_down_mf,
            file = paste(sample_name, "_MF_Down_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")

write.table(GOres_down_cc,
            file = paste(sample_name, "_CC_Down_fisher_GO_terms.txt", sep = ""),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F,
            na = "")
