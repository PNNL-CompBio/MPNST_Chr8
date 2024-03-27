# Differential expression & enrichment analyses: PDX transcriptomics, global & phospho
# Chr8: amplified vs. not amplified
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-03-12

# now including WU-225 in addition to samples shared w RNA-seq
# overview:
#### 1. Import metadata & crosstabs
#### 2. PCA plots
#### 3. Run panSEA: only shared MPNST PDX samples w known Chr8 status
#### 4. Run panSEA: all MPNST PDX samples w known Chr8 status
#### 5. Run panSEA: re-do proteomics w all MPNST PDX samples w known Chr8 status except JH-2-002 because it's not as Chr8-amplified
#### 5a. Run panSEA: re-do proteomics w all MPNST PDX samples w known Chr8 status except JH-2-002 and MN-2 because they're not as Chr8-amplified
#### 6. Pathways of interest: expression, log2FC heatmaps; TYK2, pSTAT3 boxplot

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr); library(plyr)
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
source("panSEA_helper.R")

#### 1. Import metadata & crosstabs ####
setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
)
meta.df <- readxl::read_excel("Chr8-MetaDataSheet.xlsx")
global.df <- read.table(
  "global_data/Chr8_crosstab_global_gene_corrected.txt",
  sep = "\t")
global.w.mouse.df <- read.table(
  "global_with_mouse_data/Chr8_crosstab_global_old_database_gene_corrected.txt", 
  sep = "\t")
phospho.df <- read.table(
  "phospho_data/Chr8_crosstab_phospho_SiteID_corrected.txt",
  sep = "\t")

setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/"
)
synapser::synGet("syn49554376", 
                 downloadLocation = getwd())
RNA.df <- read.table("salmon.merged.gene_counts.tsv", sep = "\t", 
                     header = TRUE)
RNA.df$gene_id <- NULL
colnames(RNA.df)[1] <- "Gene"

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)
global.w.mouse.df$Gene <- rownames(global.w.mouse.df)
phospho.df$SUB_SITE <- rownames(phospho.df)

# venn diagram
library(dplyr)
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/"
setwd(base.path)
phospho.venn <- phospho.df  %>% tidyr::extract(SUB_SITE, "Gene",
                                               remove = FALSE)
venn.list <- list('Global old' = unique(global.w.mouse.df$Gene),
                  'Global' = unique(global.df$Gene),
                  'RNA-Seq' = unique(RNA.df$Gene),
                  'Phospho' = unique(phospho.venn$Gene))
ggvenn::ggvenn(venn.list, set_name_size = 3)
ggplot2::ggsave("Chr8_venn_diagram_20240312.pdf")
ggvenn::ggvenn(venn.list, set_name_size = 3, show_percentage = FALSE)
ggplot2::ggsave("Chr8_venn_diagram_wo_percent_20240312.pdf")

venn.list <- list('Proteomics' = unique(global.w.mouse.df$Gene),
                  'RNA-Seq' = unique(RNA.df$Gene),
                  'Phospho' = unique(phospho.venn$Gene))
ggvenn::ggvenn(venn.list, set_name_size = 4, show_percentage = FALSE)
ggplot2::ggsave("Chr8_venn_diagram_global_old_20240312.pdf")

venn.list <- list('Proteomics' = unique(global.df$Gene),
                  'RNA-Seq' = unique(RNA.df$Gene),
                  'Phospho' = unique(phospho.venn$Gene))
ggvenn::ggvenn(venn.list, set_name_size = 4, show_percentage = FALSE)
ggplot2::ggsave("Chr8_venn_diagram_global_new_20240312.pdf")

# add Tube number to metadata
meta.df$Tube <- NA
tube.num <- 0
for (i in 1:nrow(meta.df)) {
  if (meta.df$SampleID[i] != "reference") {
    tube.num <- tube.num + 1
    meta.df$Tube[i] <- tube.num
  }
}

# fix sample IDs
meta.df[meta.df$SampleID == "WU-22t Cshim81", ]$SampleID <- "WU-225"
meta.df[meta.df$SampleID == "WU-487CPP1", ]$SampleID <- "WU-487"

# add Chr8 status
meta.df$Chr8_amplified <- TRUE
non.chr8.amp <- c("JH-2-055", "WU-487")
meta.df[meta.df$SampleID %in% non.chr8.amp, ]$Chr8_amplified <- FALSE
meta.df$Chr8_strongly_amplified <- TRUE
meta.df[meta.df$SampleID %in% non.chr8.amp, ]$Chr8_strongly_amplified <- FALSE
meta.df[meta.df$SampleID == "JH-2-002", ]$Chr8_strongly_amplified <- NA

meta.df$id <- paste0("X", meta.df$Tube)

#### 2. PCA plots ####
library(MSnSet.utils)

### Transcriptomics PCA
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/")
RNA.df$Gene <- NULL

# create metadata for RNA-Seq
Tube <- colnames(RNA.df)
pca.meta.df <- as.data.frame(Tube)
pca.meta.df$SampleID <- c("JH-2-009", "JH-2-055", "JH-2-079", "WU-487", 
                      "WU-536", "WU-561", "MN-2", "MN-3")
pca.meta.df$Sample <- "Sample"
pca.meta.df$SampleID <- pca.meta.df$Tube
rownames(meta.df) <- pca.meta.df$Tube

meta.df$SampleID <- c("JH-2-009", "JH-2-055", "JH-2-079", "WU-487", 
                      "WU-536", "WU-561", "MN-2", "MN-3")
m_RNA <- MSnSet(exprs = RNA.df %>% as.matrix(),
                pData = pca.meta.df)
plot_pca(m_RNA, phenotype = "SampleID") + ggtitle("RNA-Seq PCA")
ggsave("RNAseq_PDX_PCA_betterSampleID.pdf")

Chr8.RNA.df <- RNA.df[, c("X2.055_pdx", "X2.079_pdx", "WU.487_pdx", "WU.561_pdx", "mn2_pdx")]
# create metadata for RNA-Seq
Tube <- colnames(Chr8.RNA.df)
pca.meta.df <- as.data.frame(Tube)
pca.meta.df$SampleID <- c("JH-2-055", "JH-2-079", "WU-487", "WU-561", "MN-2")
pca.meta.df$Sample <- "Sample"
rownames(pca.meta.df) <- pca.meta.df$Tube
pca.meta.df$Chr8_status <- "Amplified"
non.chr8.amp.rna <- c("JH-2-055", "WU-487")
pca.meta.df[pca.meta.df$SampleID %in% non.chr8.amp.rna, ]$Chr8_status <- "Not amplified"

m_RNA <- MSnSet(exprs = Chr8.RNA.df %>% as.matrix(),
                pData = pca.meta.df)
plot_pca(m_RNA, phenotype = "SampleID") + ggtitle("RNA-Seq PCA")
ggsave("RNAseq_PDX_PCA_betterSampleID_knownChr8StatusOnly.pdf")

m_RNA <- MSnSet(exprs = Chr8.RNA.df %>% as.matrix(),
                pData = pca.meta.df)
plot_pca(m_RNA, phenotype = "Chr8_status") + ggtitle("RNA-Seq PCA")
ggsave("RNAseq_PDX_PCA_byChr8Status_knownChr8StatusOnly.pdf")

overlap.RNA.df <- RNA.df[, c("X2.055_pdx", "X2.079_pdx", "WU.487_pdx", "mn2_pdx")]
# create metadata for RNA-Seq
Tube <- colnames(overlap.RNA.df)
pca.meta.df <- as.data.frame(Tube)
pca.meta.df$SampleID <- c("JH-2-055", "JH-2-079", "WU-487", "MN-2")
pca.meta.df$Sample <- "Sample"
rownames(pca.meta.df) <- pca.meta.df$Tube
pca.meta.df$Chr8_status <- "Amplified"
non.chr8.amp.rna <- c("JH-2-055", "WU-487")
pca.meta.df[pca.meta.df$SampleID %in% non.chr8.amp.rna, ]$Chr8_status <- "Not amplified"

m_RNA <- MSnSet(exprs = overlap.RNA.df %>% as.matrix(),
                pData = pca.meta.df)
plot_pca(m_RNA, phenotype = "SampleID") + ggtitle("RNA-Seq PCA")
ggsave("RNAseq_PDX_PCA_betterSampleID_overlappingKnownChr8StatusOnly.pdf")

plot_pca(m_RNA, phenotype = "Chr8_status") + ggtitle("RNA-Seq PCA")
ggsave("RNAseq_PDX_PCA_byChr8Status_overlappingKnownChr8StatusOnly.pdf")

#### Global & phospho - new database ####
setwd(base.path)
meta <- meta.df[meta.df$Sample == "Sample",]
omics <- list("global" = global.df,
              "phospho" = phospho.df)
contrasts <- list(c(TRUE, FALSE))

contrast.type <- "Chr8_strongly_amplified"
temp.path <- file.path(base.path, "Chr8_strongly_amp_vs_not", "new database")
run_contrasts_global_phospho_human(contrasts, contrast.type, "id", meta,
                                   omics, gmt.list1 = "chr8", EA.types = c("KEGG", "Hallmark", "Positional", "Positional_Chr8_cancer"),
                                   gmt.list2 = "chr8", file.exists("chr8_gmt2_run_contrasts_global_phospho_human.rds"), 
                                   temp.path = temp.path,
                                   subfolder = FALSE,
                                   synapse_id = "syn54042241")

contrast.type <- "Chr8_amplified"
temp.path <- file.path(base.path, "Chr8_amp_vs_not", "new database")
run_contrasts_global_phospho_human(contrasts, contrast.type, "id", meta,
                                   omics, gmt.list1 = "chr8", EA.types = c("KEGG", "Hallmark", "Positional", "Positional_Chr8_cancer"),
                                   gmt.list2 = "chr8", base.path = base.path, temp.path = temp.path,
                                   subfolder = FALSE,
                                   synapse_id = "syn54042236")

#### Global - old database ####
omics <- list("global" = global.w.mouse.df)
contrasts <- list(c(TRUE, FALSE))

contrast.type <- "Chr8_strongly_amplified"
temp.path <- file.path(base.path, "Chr8_strongly_amp_vs_not", "old database")
run_contrasts_global_human(contrasts, contrast.type, "id", meta,
                           omics, gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                                "msigdb_Homo sapiens_H"), EA.types = c("KEGG", "Hallmark"),
                           base.path = base.path, temp.path = temp.path,
                           subfolder = FALSE,
                           synapse_id = "syn54042244")

contrast.type <- "Chr8_amplified"
temp.path <- file.path(base.path, "Chr8_amp_vs_not", "old database")
run_contrasts_global_human(contrasts, contrast.type, "id", meta,
                           omics, gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                                "msigdb_Homo sapiens_H"), EA.types = c("KEGG", "Hallmark"),
                           base.path = base.path, temp.path = temp.path,
                           subfolder = FALSE,
                           synapse_id = "syn54042237")

#### 6. Pathways of interest: expression, log2FC heatmaps; TYK2, pSTAT3 boxplot ####
# look at individual genes: Chr8q, Chr8 cancer, Myc targets, JAK/STAT (esp. TYK2, pSTAT3), TGF beta
# Chr8q
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")
chr8 <- unique(msigdb.info[grepl("chr8", msigdb.info$gs_name, 
                                 ignore.case = TRUE), ]$gene_symbol) # 1497
chr8p <- unique(msigdb.info[grepl("chr8p", msigdb.info$gs_name, 
                                  ignore.case = TRUE), ]$gene_symbol) # 560
chr8q <- unique(msigdb.info[grepl("chr8q", msigdb.info$gs_name, 
                                  ignore.case = TRUE), ]$gene_symbol) # 937

# Chr8 cancer genes
Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                       "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                       "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")

# Myc targets
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "H")
myc.targets <- unique(msigdb.info[grepl("MYC_TARGETS", msigdb.info$gs_name, 
                                 ignore.case = TRUE), ]$gene_symbol) # 240
myc.targetsv1 <- unique(msigdb.info[grepl("MYC_TARGETS_V1", msigdb.info$gs_name, 
                                        ignore.case = TRUE), ]$gene_symbol) # 200
myc.targetsv2 <- unique(msigdb.info[grepl("MYC_TARGETS_V2", msigdb.info$gs_name, 
                                          ignore.case = TRUE), ]$gene_symbol) # 58

# JAK/STAT
jak.stat <- unique(msigdb.info[grepl("JAK_STAT", msigdb.info$gs_name, 
                                     ignore.case = TRUE), ]$gene_symbol) 

# KRAS
kras <- unique(msigdb.info[grepl("KRAS", msigdb.info$gs_name, 
                                     ignore.case = TRUE), ]$gene_symbol) 

# TGF beta
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C2", "CP:KEGG")
tgf <- unique(msigdb.info[grepl("TGF_BETA", msigdb.info$gs_name, 
                                     ignore.case = TRUE), ]$gene_symbol) 

genesets <- list('Chr8' = chr8, 
                 'Chr8p' = chr8p,
                 'Chr8q' = chr8q,
                 'Chr8_cancer_genes' = Chr8.cancer.genes,
                 'Myc_targets' = myc.targets,
                 'Myc_targets_v1' = myc.targetsv1,
                 'Myc_targets_v2' = myc.targetsv2,
                 "JAK-STAT" = jak.stat,
                 "KRAS" = kras,
                 "TGF_Beta" = tgf)

# compile differential expression results
# all samples w known Chr8 status
setwd(paste0("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/",
             "global_with_mouse_data/analysis/"))
pdx.prot <- read.csv("Differential expression/DEG_results.csv")
pdx.prot <- pdx.prot[pdx.prot$type == "PDX Proteomics", ] # 6329 genes
pdx.prot$type <- "Proteomics"
colnames(pdx.prot)[2] <- "Gene"
pdx.prot$feature_name <- NULL
pdx.prot[pdx.prot$P.Value < 0.05 & pdx.prot$adj.P.Val < 0.05, ]$sig <- TRUE
pdx.prot <- pdx.prot[ , c("type", "Gene", "Log2FC", "minusLogFDR", "sig")]

setwd(paste0("~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/",
             "analysis/"))
pdx.transcr <- read.csv("Differential expression/DEG_results.csv") # 60347
pdx.transcr$minusLogP <- -log(pdx.transcr$P.Value, 10)
pdx.transcr$minusLogFDR <- -log(pdx.transcr$adj.P.Val, 10)
pdx.transcr$sig <- FALSE
pdx.transcr[pdx.transcr$P.Value < 0.05 & pdx.transcr$adj.P.Val < 0.05, ]$sig <- TRUE
pdx.transcr$type <- "Transcriptomics"
pdx.transcr <- pdx.transcr[ , c("type", "Gene", "Log2FC", "minusLogFDR", "sig")]

pdx.prot.transcr <- rbind(pdx.prot, pdx.transcr)

# samples overlapping between prot and transcr
setwd(paste0("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/",
             "global_with_mouse_data/analysis/PDX_RNAseq_and_proteomics"))
pdx.overlap <- read.csv("Differential expression/DEG_results.csv")
colnames(pdx.overlap)[2] <- "Gene"
pdx.overlap$feature_name <- NULL
pdx.overlap[pdx.overlap$P.Value < 0.05 & pdx.overlap$adj.P.Val < 0.05, ]$sig <- TRUE
pdx.overlap <- pdx.overlap[ , c("type", "Gene", "Log2FC", "minusLogFDR", "sig")]

setwd(paste0("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/",
             "global_with_mouse_data/analysis/"))
phospho <- read.csv("GSEA/GSEA_results.csv")
phospho <- phospho[phospho$type == "Phospho-proteomics",]

omics <- list('PDX_samples' = pdx.prot.transcr,
              'PDX_overlap_samples' = pdx.overlap)
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/global_with_mouse_data/analysis/"
setwd(base.path)
abs.max.Log2FC <- 0
abs.max.minusLogFDR <- 0
for (i in 1:length(omics)) {
  for (j in 1:length(genesets)) {
    # filter for gene set
    DEG.df <- omics[[i]]
    DEG.df <- DEG.df[DEG.df$Gene %in% genesets[[j]], ]

    # generate data frames for heatmaps
    Log2FC.df <- reshape2::dcast(DEG.df, Gene ~ type,
                                 value.var = "Log2FC", fill = NA
    )
    minusLogFDR.df <- reshape2::dcast(DEG.df, Gene ~ type,
                                      value.var = "minusLogFDR", fill = NA
    )
    # sig.df <- reshape2::dcast(DEG.df, Gene ~ type,
    #                                   value.var = "sig", fill = NA
    # ) # none are sig
    write.csv(Log2FC.df, paste0(names(omics)[i], "_", names(genesets)[j], "_Log2FC.csv"), row.names = FALSE)
    write.csv(minusLogFDR.df, paste0(names(omics)[i], "_", names(genesets)[j], "_minusLogFDR.csv"), row.names = FALSE)
    #write.csv(sig.df, paste0(names(omics)[i], "_", names(genesets)[j], "_sig.csv"), row.names = FALSE)

    if (max(abs(DEG.df$Log2FC), na.rm = TRUE) > abs.max.Log2FC) {
      abs.max.Log2FC <- max(abs(DEG.df$Log2FC), na.rm = TRUE)
    } # 9.591479, so using range of -10, 10 for heatmap colors
    
    if (max(abs(DEG.df$minusLogFDR), na.rm = TRUE) > abs.max.minusLogFDR) {
      abs.max.minusLogFDR <- max(abs(DEG.df$minusLogFDR), na.rm = TRUE)
    } # 0.1577173, so using range of 0, 0.2 for heatmap sizes
  }
}

# repeat but with phospho KSEA results
omics <- list('Phospho' = phospho)
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/global_with_mouse_data/analysis/"
setwd(base.path)
abs.max.NES <- 0
abs.max.minusLogFDR <- 0
for (i in 1:length(omics)) {
  for (j in 1:length(genesets)) {
    # filter for gene set
    DEG.df <- omics[[i]]
    DEG.df <- DEG.df[DEG.df$Feature_set %in% genesets[[j]], ]
    
    # generate data frames for heatmaps
    NES.df <- reshape2::dcast(DEG.df, Feature_set ~ type,
                                 value.var = "NES", fill = NA
    )
    minusLogFDR.df <- reshape2::dcast(DEG.df, Feature_set ~ type,
                                      value.var = "minusLogFDR", fill = NA
    )
    sig.df <- reshape2::dcast(DEG.df, Feature_set ~ type,
                                      value.var = "sig", fill = NA
    )
    write.csv(NES.df, paste0(names(omics)[i], "_", names(genesets)[j], "_NES.csv"), row.names = FALSE)
    write.csv(minusLogFDR.df, paste0(names(omics)[i], "_", names(genesets)[j], "_minusLogFDR.csv"), row.names = FALSE)
    write.csv(sig.df, paste0(names(omics)[i], "_", names(genesets)[j], "_sig.csv"), row.names = FALSE)
    
    if (max(abs(DEG.df$NES), na.rm = TRUE) > abs.max.NES) {
      abs.max.NES <- max(abs(DEG.df$NES), na.rm = TRUE)
    } # 2.62, so using range of -3, 3 for heatmap colors
    
    if (max(abs(DEG.df$minusLogFDR), na.rm = TRUE) > abs.max.minusLogFDR) {
      abs.max.minusLogFDR <- max(abs(DEG.df$minusLogFDR), na.rm = TRUE)
    } # 2.72, so using range of 0, 3 for heatmap sizes
  }
}

# compile expression data
setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
)
meta.df <- readxl::read_excel("Chr8-MetaDataSheet.xlsx")
global.df <- read.table(
  "global_data/Chr8_crosstab_global_gene_corrected_sampleNames.txt",
  sep = "\t")
colnames(global.df) <- global.df[1, ]
global.df <- global.df[2:nrow(global.df), ]
rownames(global.df) <- global.df[ , 1]
colnames(global.df)[1] <- "Gene"
#global.df <- global.df[ , 2:ncol(global.df)]
global.overlap <- global.df[ , c("Gene", "WU-487_rep1", "MN-2_rep1", 
                                    "JH-2-055_rep1", "JH-2-079_rep1", 
                                    "JH-2-055_rep2", "MN-2_rep2",
                                    "WU-487_rep2", "JH-2-079_rep2")]

global.w.mouse.df <- read.table(
  "global_with_mouse_data/Chr8_crosstab_global_gene_corrected_sampleNames.txt", 
  sep = "\t")
colnames(global.w.mouse.df) <- global.w.mouse.df[1, ]
global.w.mouse.df <- global.w.mouse.df[2:nrow(global.w.mouse.df), ]
rownames(global.w.mouse.df) <- global.w.mouse.df[ , 1]
colnames(global.w.mouse.df)[1] <- "Gene"
#global.w.mouse.df <- global.w.mouse.df[ , 2:ncol(global.w.mouse.df)]
global.w.mouse.overlap <- global.w.mouse.df[ , c("Gene", "WU-487_rep1", "MN-2_rep1", 
                                    "JH-2-055_rep1", "JH-2-079_rep1", 
                                    "JH-2-055_rep2", "MN-2_rep2",
                                    "WU-487_rep2", "JH-2-079_rep2")]

phospho.df <- read.table(
  "phospho_data/Chr8_crosstab_phospho_SiteID_corrected_sampleNames.txt",
  sep = "\t")
colnames(phospho.df) <- phospho.df[1, ]
phospho.df <- phospho.df[2:nrow(phospho.df), ]
rownames(phospho.df) <- phospho.df[ , 1]
phospho.df$Gene <- rownames(phospho.df)
phospho.df <- phospho.df[ , 2:ncol(phospho.df)]
phospho.df <- phospho.df %>% tidyr::extract(Gene, "Gene",
                                              remove = FALSE)
phospho.df <- phospho.df[ , c("Gene", colnames(phospho.df)[1:(ncol(phospho.df)-1)])]
phospho.overlap <- phospho.df[ , c("Gene", "WU-487_rep1", "MN-2_rep1", 
                                    "JH-2-055_rep1", "JH-2-079_rep1", 
                                    "JH-2-055_rep2", "MN-2_rep2",
                                    "WU-487_rep2", "JH-2-079_rep2")]

setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/"
)
RNA.df <- read.table("salmon.merged.gene_counts.tsv", sep = "\t", 
                     header = TRUE)
RNA.df$gene_id <- NULL
colnames(RNA.df) <- c("Gene", "JH-2-009", "JH-2-055", "JH-2-079", "WU-487",
                      "WU-536", "WU-561", "MN-2", "MN-3")
RNA.overlap <- RNA.df[ , c("Gene", "JH-2-055", "JH-2-079", "WU-487", "MN-2")]

omics <- list('PDX_global' = global.w.mouse.df,
              'PDX_RNAseq' = RNA.df,
              'global' = global.df,
              'phospho' = phospho.df,
              'PDX_global_overlap' = global.w.mouse.overlap,
              'PDX_RNAseq_overlap' = RNA.overlap,
              'global_overlap' = global.overlap,
              'phospho_overlap' = phospho.overlap)
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/heatmaps/"
setwd(base.path)
#abs.max <- 0
for (i in 1:length(omics)) {
  for (j in 1:length(genesets)) {
    # filter for gene set
    expr.df <- omics[[i]]
    expr.df <- expr.df[expr.df$Gene %in% genesets[[j]], ]
    write.csv(expr.df, paste0(names(omics)[i], "_", names(genesets)[j], ".csv"), row.names = FALSE)
  }
}

# proteomics wo JH-2-002
# expression
omics <- list('global' = global.w.mouse.stricter,
              'global_new_database' = global.stricter,
              'phospho' = phospho.stricter)
meta.stricter$id.info <- c("WU-487_rep1", "MN-2_rep1", "JH-2-055_rep1", 
                           "JH-2-079_rep1", "JH-2-055_rep2", "WU-225_rep1",
                           "MN-2_rep2", "WU-487_rep2", "WU-225_rep2", "JH-2-079_rep2")
global.w.mouse.stricter.w.names <- global.w.mouse.stricter
global.stricter.w.names <- global.stricter
phospho.stricter.w.names <- phospho.stricter
colnames(global.w.mouse.stricter.w.names) <- c(meta.stricter$id.info, "Gene")
colnames(global.stricter.w.names) <- c(meta.stricter$id.info, "Gene")
colnames(phospho.stricter.w.names) <- c(meta.stricter$id.info, "SUB_SITE")

global.w.mouse.stricter.w.names <- global.w.mouse.stricter.w.names[ , c("Gene", meta.stricter$id.info)]
global.stricter.w.names <- global.stricter.w.names[ , c("Gene", meta.stricter$id.info)]
phospho.stricter.w.names <- phospho.stricter.w.names[ , c("SUB_SITE", meta.stricter$id.info)]

omics <- list('global' = global.w.mouse.stricter.w.names,
              'global_new_database' = global.stricter.w.names,
              'phospho' = phospho.stricter.w.names)

synapse_ids <- c("syn53632588", "syn53632634", "syn53632674")
setwd(file.path(base.path, "Chr8_amp_vs_not_amp"))
for (i in 1:length(omics)) {
  setwd(file.path(base.path, "Chr8_amp_vs_not_amp", names(omics)[i]))
  dir.create("pathways_of_interest")
  setwd("pathways_of_interest")
  resultsFolder <- 
    synapser::synStore(synapser::Folder("pathways_of_interest",
                                        parent = synapse_ids[i]))
  expr.df <- omics[[i]]
  for (j in 1:length(genesets)) {
    setwd(file.path(base.path, "Chr8_amp_vs_not_amp", names(omics)[i], "pathways_of_interest"))
    dir.create(names(genesets)[j])
    setwd(names(genesets)[j])
    gene.expr.df <- data.frame()
    # filter for gene set
    if (names(omics)[i] == "phospho") {
      for (k in 1:length(genesets[[j]])) {
        temp.expr.df <- expr.df[grepl(genesets[[j]][k], expr.df$SUB_SITE, ignore.case = TRUE), ]
        gene.expr.df <- rbind(gene.expr.df, temp.expr.df)
      }
    } else {
      gene.expr.df <- expr.df[tolower(expr.df$Gene) %in% tolower(genesets[[j]]), ]
    }
    if (nrow(gene.expr.df) > 0) {
      write.csv(gene.expr.df, paste0(names(omics)[i], "_", names(genesets)[j], "_expression.csv"), row.names = FALSE)
      geneFolder <- 
        synapser::synStore(synapser::Folder(names(genesets)[j],
                                            parent = resultsFolder))
      CSV.file <- synapser::File(paste0(names(omics)[i], "_", names(genesets)[j], "_expression.csv"),
                     parent = geneFolder)
      synapser::synStore(CSV.file)
    }
  }
}

# log2FC
omics <- c('global', 'global_new_database', 'phospho')
setwd(file.path(base.path, "Chr8_amp_vs_not_amp"))
abs.max.Log2FC <- 0
abs.max.minusLogFDR <- 0
for (i in 1:length(omics)) {
  setwd(file.path(base.path, "Chr8_amp_vs_not_amp", omics[i]))
  #dir.create("pathways_of_interest")
  expr.df <- read.csv("Differential_expression/Differential_expression_results.csv")
  setwd("pathways_of_interest")
  resultsFolder <- 
    synapser::synStore(synapser::Folder("pathways_of_interest",
                                        parent = synapse_ids[i]))
  
  # filter for gene set
  for (j in 1:length(genesets)) {
    setwd(file.path(base.path, "Chr8_amp_vs_not_amp", omics[i], "pathways_of_interest"))
    #dir.create(names(genesets)[j])
    setwd(names(genesets)[j])
    if (omics[i] == "phospho") {
      DEG.df <- data.frame()
      for (k in 1:length(genesets[[j]])) {
        temp.DEG.df <- expr.df[grepl(genesets[[j]][k], expr.df$SUB_SITE, ignore.case = TRUE), ]
        DEG.df <- rbind(DEG.df, temp.DEG.df)
      }
    } else {
      DEG.df <- expr.df[tolower(expr.df$Gene) %in% tolower(genesets[[j]]), ]
    }
    
    if (nrow(DEG.df) > 0) {
      write.csv(DEG.df, paste0(omics[i], "_", names(genesets)[j], "_differential_expression.csv"), row.names = FALSE)
      geneFolder <- 
        synapser::synStore(synapser::Folder(names(genesets)[j],
                                            parent = resultsFolder))
      CSV.file <- synapser::File(paste0(omics[i], "_", names(genesets)[j], "_differential_expression.csv"),
                                 parent = geneFolder)
      synapser::synStore(CSV.file)
      
      if (max(abs(DEG.df$Log2FC), na.rm = TRUE) > abs.max.Log2FC) {
        abs.max.Log2FC <- max(abs(DEG.df$Log2FC), na.rm = TRUE)
      } # 7.293, so using range of -8, 8 for heatmap colors
      
      if (max(-log(DEG.df$adj.P.Val, 10), na.rm = TRUE) > abs.max.minusLogFDR) {
        abs.max.minusLogFDR <- max(-log(DEG.df$adj.P.Val, 10), na.rm = TRUE)
      } # 4.223, so using range of 0, 5 for heatmap sizes 
    }
  }
}

# NES for KSEA, substrate enrichment
omics <- list("Substrate_enrichment" = GSEA3$mGSEA.results$all.results$Phospho$result,
              "KSEA" = GSEA3_hallmark$mGSEA.results$all.results$Phospho$result)
omics <- c("Substrate_enrichment", "KSEA")
abs.max.NES <- 0
abs.max.minusLogFDR <- 0
for (i in 1:length(omics)) {
  setwd(file.path(base.path, "Chr8_amp_vs_not_amp", "phospho"))
  expr.df <- read.csv(paste0(omics[i], "/", omics[i], "_results.csv"))
  setwd("pathways_of_interest")
  resultsFolder <- 
    synapser::synStore(synapser::Folder("pathways_of_interest",
                                        parent = "syn53632674"))
  for (j in 1:length(genesets)) {
    # filter for gene set
    DEG.df <- expr.df[tolower(expr.df$Feature_set) %in% tolower(genesets[[j]]), ]
    if (nrow(DEG.df) > 0) {
      setwd(file.path(base.path, "Chr8_amp_vs_not_amp", "phospho", names(genesets)[j]))
      write.csv(DEG.df, paste0("phospho", "_", names(genesets)[j], "_", omics[i], ".csv"), row.names = FALSE)
      geneFolder <- 
        synapser::synStore(synapser::Folder(names(genesets)[j],
                                            parent = resultsFolder))
      CSV.file <- synapser::File(paste0("phospho", "_", names(genesets)[j], "_", omics[i], ".csv"),
                                 parent = geneFolder)
      synapser::synStore(CSV.file)
      
      if (max(abs(DEG.df$NES), na.rm = TRUE) > abs.max.NES) {
        abs.max.NES <- max(abs(DEG.df$NES), na.rm = TRUE)
      } # 2.62, so using range of -3, 3 for heatmap colors
      
      if (max(-log(DEG.df$FDR_q_value, 10), na.rm = TRUE) > abs.max.minusLogFDR) {
        abs.max.minusLogFDR <- max(-log(DEG.df$FDR_q_value, 10), na.rm = TRUE)
      } # 2.72, so using range of 0, 3 for heatmap sizes 
    }
  }
}
