# Differential expression & enrichment analyses: RNAseq, global PDX; global, phospho cell lines
# Chr8: amplified vs. not amplified
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-01-19

#### first analysis: only shared MPNST PDX samples w known Chr8 status
#### second analysis: only shared MPNST samples w known Chr8 status
#### third analysis: all MPNST samples w known Chr8 status

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr); library(plyr)

#### 1. Import metadata & crosstabs ####
setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
)
meta.df <- readxl::read_excel("Chr8-MetaDataSheet.xlsx")
global.df <- read.table(
  "global_data/Chr8_crosstab_global_gene_corrected.txt",
  sep = "\t")
global.w.mouse.df <- read.table(
  "global_with_mouse_data/Chr8_crosstab_global_gene_corrected.txt", 
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
phospho.venn <- phospho.df  %>% tidyr::extract(SUB_SITE, "Gene",
                                               remove = FALSE)
venn.list <- list('Proteomics' = unique(global.df$Gene),
                  'PDX Proteomics' = unique(global.w.mouse.df$Gene),
                  'PDX RNA-Seq' = unique(RNA.df$Gene),
                  'Phospho' = unique(phospho.df$SUB_SITE))
ggvenn::ggvenn(venn.list, set_name_size = 4)
ggplot2::ggsave("Chr8_venn_diagram.pdf")

venn.list <- list('Proteomics' = unique(global.df$Gene),
                  'PDX Proteomics' = unique(global.w.mouse.df$Gene),
                  'PDX RNA-Seq' = unique(RNA.df$Gene),
                  'Phospho' = unique(phospho.venn$Gene))
ggvenn::ggvenn(venn.list, set_name_size = 4)
ggplot2::ggsave("Chr8_venn_diagram_phospho_genes_instead_of_subsites.pdf")

# add Tube number to metadata
meta.df$Tube <- NA
tube.num <- 0
for (i in 1:nrow(meta.df)) {
  if (meta.df$SampleID[i] != "reference") {
    tube.num <- tube.num + 1
    meta.df$Tube[i] <- tube.num
  }
}


# add Chr8 status to metadata
meta.df$Chr8_amp <- TRUE
meta.df[meta.df$Sample == "Reference" | 
          meta.df$SampleID %in% c("JH-2-055", "WU-487CPP1"), ]$Chr8_amp <- FALSE

#### first analysis: only shared MPNST PDX samples w known Chr8 status
# identify global samples in RNAseq
meta.df$RNAseq <- FALSE
prot.samples.in.RNAseq <- c("WU-487CPP1", "MN-2", "JH-2-055", "JH-2-079")
meta.df[meta.df$SampleID %in% prot.samples.in.RNAseq, ]$RNAseq <- TRUE
prot.tubes.in.RNAseq <- paste0("X", meta.df[meta.df$RNAseq, ]$Tube)
#global.w.mouse.df1 <- global.w.mouse.df[ , c("Gene", prot.tubes.in.RNAseq)]
global.w.mouse.df1 <- plyr::ddply(global.w.mouse.df, .(Gene), summarize,
                                 'WU-487' = mean(c(X1, X9), na.rm = TRUE),
                                 'MN-2' = mean(c(X2, X8), na.rm = TRUE),
                                 'JH-2-055' = mean(c(X3, X6), na.rm = TRUE),
                                 'JH-2-079' = mean(c(X4, X12), na.rm = TRUE))

# identify RNAseq samples in global
RNAseq.samples.in.prot <- c("WU.487_pdx", "mn2_pdx", "X2.055_pdx", "X2.079_pdx")
RNA.df1 <- RNA.df[ , c("Gene", RNAseq.samples.in.prot)]
RNA.df1 <- plyr::ddply(RNA.df1, .(Gene), summarize,
                       'WU-487' = mean(WU.487_pdx, na.rm = TRUE),
                       'MN-2' = mean(mn2_pdx, na.rm = TRUE),
                       'JH-2-055' = mean(X2.055_pdx, na.rm = TRUE),
                       'JH-2-079' = mean(X2.079_pdx, na.rm = TRUE))
# colnames(RNA.df1) <- colnames(global.w.mouse.df1)
# RNA.df1 <- plyr::ddply(RNA.df1, .(Gene), summarize,
#                                   'WU-487' = mean('WU-487', na.rm = TRUE),
#                                   'MN-2' = mean('MN-2', na.rm = TRUE),
#                                   'JH-2-055' = mean('JH-2-055', na.rm = TRUE),
#                                   'JH-2-079' = mean('JH-2-079', na.rm = TRUE))

# run panSEA
data.list <- list(RNA.df1, global.w.mouse.df1)
types <- c("Transcriptomics", "Proteomics")
library(Biobase)
panSEA1 <- panSEA::panSEA(data.list, types, 
                          group.names=c("Chr8_amp", "not_amp"), 
                          group.samples = list(c("JH-2-079", "MN-2"),
                                               c("JH-2-055", "WU-487")))

DMEA1_prot <- panSEA::panSEA(data.list, types, GSEA = FALSE,
                          group.names=c("Chr8_amp", "not_amp"), 
                          group.samples = list(c("JH-2-079", "MN-2"),
                                               c("JH-2-055", "WU-487")),
                          expression = list("adherent CCLE", prot.df.noNA))

DEG.files <- list("DEG_results.csv" =
                     panSEA1$mDEG.results$compiled.results$results,
                   "DEG_mean_results.csv" =
                     panSEA1$mDEG.results$compiled.results$mean.results,
                   "DEG_correlation_matrix.pdf" =
                     panSEA1$mDEG.results$compiled.results$corr.matrix,
                   "DEG_venn_diagram.pdf" =
                     panSEA1$mDEG.results$compiled.results$venn.diagram,
                   "DEG_dot_plot.pdf" =
                     panSEA1$mDEG.results$compiled.results$dot.plot)
DMEA.files <- list("DMEA_results.csv" =
                     panSEA1$mDMEA.results$compiled.results$results,
                   "DMEA_mean_results.csv" =
                     panSEA1$mDMEA.results$compiled.results$mean.results,
                   "DMEA_correlation_matrix.pdf" =
                     panSEA1$mDMEA.results$compiled.results$corr.matrix,
                   "DMEA_venn_diagram.pdf" =
                     panSEA1$mDMEA.results$compiled.results$venn.diagram,
                   "DMEA_dot_plot.pdf" =
                     panSEA1$mDMEA.results$compiled.results$dot.plot,
                   "DMEA_interactive_network.graph.html" =
                     panSEA1$mDMEA.network$interactive,
                   "DMEA_transcriptomics_volcano_plot.pdf" =
                     panSEA1$mDMEA.results$all.results$Transcriptomics$volcano.plot,
                   "DMEA_proteomics_volcano_plot.pdf" =
                     panSEA1$mDMEA.results$all.results$Proteomics$volcano.plot,
                   "DMEA_proteomics_mountain_plot_adenosine_receptor_antagonist.pdf" =
                     panSEA1$mDMEA.results$all.results$Proteomics$mtn.plot[["adenosine receptor antagonist"]],
                   "DMEA_proteomics_mountain_plot_mTOR_inhibitor.pdf" =
                     panSEA1$mDMEA.results$all.results$Proteomics$mtn.plot[["mTOR inhibitor"]],
                   "DMEA_proteomics_mountain_plot_TGF_beta_receptor_inhibitor.pdf" =
                     panSEA1$mDMEA.results$all.results$Proteomics$mtn.plot[["TGF beta receptor inhibitor"]],
                   "DMEA_transcriptomics_mountain_plot_CHK_inhibitor.pdf" =
                     panSEA1$mDMEA.results$all.results$Transcriptomics$mtn.plot[["CHK inhibitor"]],
                   "DMEA_proteomics_scatter_plots.pdf" =
                     panSEA1$mDMEA.results$all.results$Proteomics$corr.scatter.plots,
                   "DMEA_transcriptomics_scatter_plots.pdf" =
                     panSEA1$mDMEA.results$all.results$Transcriptomics$corr.scatter.plots,
                   "DMEA_transcriptomics_correlations.csv" =
                     panSEA1$mDMEA.results$all.results$Transcriptomics$corr.result,
                   "DMEA_proteomics_correlations.csv" =
                     panSEA1$mDMEA.results$all.results$Proteomics$corr.result)
GSEA.files <- list("GSEA_results.csv" =
                     panSEA1$mGSEA.results$compiled.results$results,
                   "GSEA_mean_results.csv" =
                     panSEA1$mGSEA.results$compiled.results$mean.results,
                   "GSEA_correlation_matrix.pdf" =
                     panSEA1$mGSEA.results$compiled.results$corr.matrix,
                   "GSEA_venn_diagram.pdf" =
                     panSEA1$mGSEA.results$compiled.results$venn.diagram,
                   "GSEA_dot_plot.pdf" =
                     panSEA1$mGSEA.results$compiled.results$dot.plot,
                   "GSEA_interactive_network.graph.html" =
                     panSEA1$mGSEA.network$interactive,
                   "GSEA_transcriptomics_volcano_plot.pdf" =
                     panSEA1$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                   "GSEA_proteomics_volcano_plot.pdf" =
                     panSEA1$mGSEA.results$all.results$Proteomics$volcano.plot,
                   "GSEA_transcriptomics_mountain_plot_KEGG_RIBOSOME.pdf" =
                     panSEA1$mGSEA.results$all.results$Transcriptomics$mtn.plot[["KEGG_RIBOSOME"]])

# generate gmt.features beforehand to save time
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "H")

# extract necessary info into data frame
msigdb.info <- as.data.frame(msigdb.info[, c(
  "gene_symbol",
  "gs_name",
  "gs_description"
)])

gmt <- DMEA::as_gmt(
  msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
  descriptions = "gs_description"
) 
gmt.features <- list(gmt, gmt)
panSEA1_hallmarkGSEA <- panSEA::panSEA(data.list, types, DMEA = FALSE,
                              group.names=c("Chr8_amp", "not_amp"), 
                              group.samples = list(c("JH-2-079", "MN-2"), 
                                                   c("JH-2-055", "WU-487")),
                              gmt.features = gmt.features,
                              scatter.plots = FALSE)
panSEA1_hallmarkGSEA <- panSEA1_C1GSEA
# generate gmt.features beforehand to save time
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")

# extract necessary info into data frame
msigdb.info <- as.data.frame(msigdb.info[, c(
  "gene_symbol",
  "gs_name",
  "gs_description"
)])

gmt <- DMEA::as_gmt(
  msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
  descriptions = "gs_description"
) 
gmt.features <- list(gmt, gmt)
panSEA1_C1GSEA <- panSEA::panSEA(data.list, types, DMEA = FALSE,
                                       group.names=c("Chr8_amp", "not_amp"), 
                                       group.samples = list(c("JH-2-079", "MN-2"), 
                                                            c("JH-2-055", "WU-487")),
                                       gmt.features = gmt.features, 
                                 scatter.plots = FALSE)
# proteomics didn't have coverage for 2+ gene sets
panSEA1_C1GSEA <- panSEA::panSEA(list(RNA.df1), "Transcriptomics", DMEA = FALSE,
                                 group.names=c("Chr8_amp", "not_amp"), 
                                 group.samples = list(c("JH-2-079", "MN-2"), 
                                                      c("JH-2-055", "WU-487")),
                                 gmt.features = list(gmt), 
                                 scatter.plots = FALSE)

# try w Chr8 cancer-associated genes
Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                       "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                       "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")
gmt$genesets[[length(gmt$genesets)+1]] <- Chr8.cancer.genes
gmt$geneset.names[[length(gmt$geneset.names)+1]] <- "Chr8 cancer-associated genes"
gmt$geneset.descriptions[[length(gmt$geneset.descriptions)+1]] <- "Chr8 cancer-associated genes"
gmt.features <- list(gmt, gmt)
panSEA1_C1GSEA_custom <- panSEA::panSEA(data.list, types, DMEA = FALSE,
                                 group.names=c("Chr8_amp", "not_amp"), 
                                 group.samples = list(c("JH-2-079", "MN-2"), 
                                                      c("JH-2-055", "WU-487")),
                                 gmt.features = gmt.features, 
                                 scatter.plots = FALSE)
panSEA1_C1GSEA_custom <- panSEA::panSEA(list(RNA.df1), "Transcriptomics", DMEA = FALSE,
                                        group.names=c("Chr8_amp", "not_amp"), 
                                        group.samples = list(c("JH-2-079", "MN-2"), 
                                                             c("JH-2-055", "WU-487")),
                                        gmt.features = list(gmt), 
                                        scatter.plots = FALSE)

GSEA.H.files <- list("GSEA_results.csv" =
                     panSEA1_hallmarkGSEA$mGSEA.results$compiled.results$results,
                   "GSEA_mean_results.csv" =
                     panSEA1_hallmarkGSEA$mGSEA.results$compiled.results$mean.results,
                   "GSEA_correlation_matrix.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$compiled.results$corr.matrix,
                   "GSEA_venn_diagram.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$compiled.results$venn.diagram,
                   "GSEA_dot_plot.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$compiled.results$dot.plot,
                   "GSEA_interactive_network.graph.html" =
                     panSEA1_hallmarkGSEA$mGSEA.network$interactive,
                   "GSEA_transcriptomics_volcano_plot.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                   "GSEA_proteomics_volcano_plot.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$all.results$Proteomics$volcano.plot,
                   "GSEA_transcriptomics_mountain_plot_HALLMARK_MYC_TARGETS_V2.pdf" =
                     panSEA1_hallmarkGSEA$mGSEA.results$all.results$Transcriptomics$mtn.plots[["HALLMARK_MYC_TARGETS_V2"]])

GSEA.C1.files <- list("GSEA_transcriptomics_results.csv" =
                       panSEA1_C1GSEA$mGSEA.results$all.results$Transcriptomics$result,
                      "GSEA_interactive_network.graph.html" =
                        panSEA1_C1GSEA$mGSEA.network$interactive,
                     "GSEA_transcriptomics_volcano_plot.pdf" =
                       panSEA1_C1GSEA$mGSEA.results$all.results$Transcriptomics$volcano.plot)

GSEA.C1.custom.files <- list("GSEA_transcriptomics_results.csv" =
                        panSEA1_C1GSEA_custom$mGSEA.results$all.results$Transcriptomics$result,
                      "GSEA_interactive_network.graph.html" =
                        panSEA1_C1GSEA_custom$mGSEA.network$interactive,
                      "GSEA_transcriptomics_volcano_plot.pdf" =
                        panSEA1_C1GSEA_custom$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                      "GSEA_transcriptomics_mountain_plot_Chr8_cancer-associated_genes.pdf" =
                        panSEA1_C1GSEA_custom$mGSEA.results$all.results$Transcriptomics$mtn.plots[["Chr8 cancer-associated genes"]])

all.files <- list('GSEA_KEGG_LEGACY' = GSEA.files,
                  'DMEA' = DMEA.files,
                  'GSEA_HALLMARK' = GSEA.H.files,
                  'GSEA_POSITIONAL' = GSEA.C1.files
                  #,
                  #'GSEA_POSITIONAL_CUSTOM' = GSEA.C1.custom.files
                  )
all.files <- list('GSEA_POSITIONAL' = GSEA.C1.files,
                  'GSEA_POSITIONAL_CUSTOM' = GSEA.C1.custom.files
)
all.files <- list('Differential expression' = DEG.files)
all.files <- list('DMEA' = DMEA.files)

DMEA.prot.files <- list("DMEA_results.csv" =
                     DMEA1_prot$mDMEA.results$compiled.results$results,
                   "DMEA_mean_results.csv" =
                     DMEA1_prot$mDMEA.results$compiled.results$mean.results,
                   "DMEA_correlation_matrix.pdf" =
                     DMEA1_prot$mDMEA.results$compiled.results$corr.matrix,
                   "DMEA_venn_diagram.pdf" =
                     DMEA1_prot$mDMEA.results$compiled.results$venn.diagram,
                   "DMEA_dot_plot.pdf" =
                     DMEA1_prot$mDMEA.results$compiled.results$dot.plot,
                   "DMEA_interactive_network.graph.html" =
                     DMEA1_prot$mDMEA.network$interactive,
                   "DMEA_transcriptomics_volcano_plot.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Transcriptomics$volcano.plot,
                   "DMEA_proteomics_volcano_plot.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$volcano.plot,
                   "DMEA_proteomics_mountain_plot_adenosine_receptor_antagonist.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$mtn.plot[["adenosine receptor antagonist"]],
                   "DMEA_proteomics_mountain_plot_FGFR_inhibitor.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$mtn.plot[["FGFR inhibitor"]],
                   "DMEA_proteomics_mountain_plot_RET_tyrosine_kinase_inhibitor.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$mtn.plot[["RET tyrosine kinase inhibitor"]],
                   "DMEA_proteomics_mountain_plot_ALK_tyrosine_kinase_receptor_inhibitor.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$mtn.plot[["ALK tyrosine kinase receptor inhibitor"]],
                   "DMEA_transcriptomics_scatter_plots.pdf" =
                     DMEA1_prot$mDMEA.results$all.results$Transcriptomics$corr.scatter.plots,
                   "DMEA_transcriptomics_correlations.csv" =
                     DMEA1_prot$mDMEA.results$all.results$Transcriptomics$corr.result,
                   "DMEA_proteomics_correlations.csv" =
                     DMEA1_prot$mDMEA.results$all.results$Proteomics$corr.result)

all.files <- list('DMEA_prot' = DMEA.prot.files)

setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/global_with_mouse_data/analysis/"
)
saveRDS(panSEA1_C1GSEA_custom, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_positional_custom.rds")) # 8.5 GB for 7 contrasts
saveRDS(panSEA1_C1GSEA, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_positional.rds")) # 8.5 GB for 7 contrasts
saveRDS(panSEA1_hallmarkGSEA, file=paste0("Chr8_", "RNAseq_and_global_PDX", "_GSEA_hallmark.rds")) # 8.5 GB for 7 contrasts
saveRDS(panSEA1, file=paste0("Chr8_", "RNAseq_and_global_PDX", "_panSEA_CCLE.rds")) # 8.5 GB for 7 contrasts
saveRDS(DMEA1_prot, file=paste0("Chr8_", "RNAseq_and_global_PDX", "_DMEA_CCLE_proteomics.rds")) # 8.5 GB for 7 contrasts
panSEA1 <- readRDS(file=paste0("Chr8_", "RNAseq_and_global_PDX", "_panSEA_CCLE.rds"))

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
synapse_id_map <- c("syn53354456" = "global_with_mouse_data/")
k <- 1
for (i in 1:length(all.files)) {
  # create local folder for subset of results
  setwd(paste0(base.path, synapse_id_map[k], "analysis/PDX_RNAseq_and_proteomics"))
  dir.create(names(all.files)[i])
  setwd(names(all.files)[i])
  
  # save results locally
  temp.files <- all.files[[i]]
  
  CSV.files <- names(temp.files)[grepl(".csv", names(temp.files))]
  for (j in 1:length(CSV.files)) {
    write.csv(temp.files[[CSV.files[j]]], CSV.files[j], row.names = FALSE)
  }
  
  PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files))]
  for (j in 1:length(PDF.files)) {
    # if (grepl("dot", PDF.files[j])) { # wide plot for dot plots
    #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
    #                   device = "pdf", width = 35)
    # } else {
    #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
    #                   device = "pdf")
    # }
    ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                    device = "pdf")
  }
  
  HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
  if (length(HTML.files) > 0) {
    for (j in 1:length(HTML.files)) {
      visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j]) 
    }
  }
  
  # create folder on Synpase for subset of results
  dataFolder <- 
    synapser::synStore(synapser::Folder(names(all.files)[i],
                                        parent = names(synapse_id_map)[k]))
  
  # upload results to Synapse
  CSVs <- lapply(as.list(CSV.files), synapser::File,
                 parent = dataFolder)
  lapply(CSVs, synapser::synStore)
  
  PDFs <- lapply(as.list(PDF.files), synapser::File,
                 parent = dataFolder)
  lapply(PDFs, synapser::synStore)
  
  if (length(HTML.files) > 0) {
    HTMLs <- lapply(HTML.files, synapser::File,
                    parent = dataFolder)
    lapply(HTMLs, synapser::synStore) 
  }
}
#panSEA.results <- NULL # make space to process next omics type

#### second analysis: only shared MPNST samples w known Chr8 status


#### third analysis: all MPNST samples w known Chr8 status
### DMEA
data.list <- list(global.df, global.w.mouse.df)
types <- c("Proteomics", "PDX Proteomics")
Chr8.prot.samples <- paste0("X", meta.df[meta.df$Chr8_amp & 
                                           meta.df$Sample == "Sample", ]$Tube)
non.Chr8.prot.samples <- paste0("X", meta.df[!meta.df$Chr8_amp & 
                                           meta.df$Sample == "Sample", ]$Tube)
DMEA3 <- panSEA::panSEA(data.list, types, GSEA = FALSE, 
                        group.names = c("Chr8_amp","not_amp"),
                        group.samples = list(Chr8.prot.samples,
                                             non.Chr8.prot.samples))

data.list <- list(global.w.mouse.df, phospho.df)
types <- c("PDX Proteomics", "Phospho-proteomics")
gmt.ksea <- readRDS(paste0(
  base.path,"phospho_data/gmt_PNNL_kinase-substrate_Chr8_5344.rds"))
GSEA3 <- panSEA::panSEA(data.list, types, DMEA = FALSE, 
                        feature.names = c("Gene", "SUB_SITE"),
                        group.names = c("Chr8_amp","not_amp"),
                        group.samples = list(Chr8.prot.samples,
                                             non.Chr8.prot.samples),
                        gmt.features = list("msigdb_Homo sapiens_C2_CP:KEGG",
                                            gmt.ksea))
GSEA3_hallmark <- panSEA::panSEA(list(global.w.mouse.df), "PDX Proteomics", 
                                 DMEA = FALSE,
                                 group.names = c("Chr8_amp","not_amp"),
                                 group.samples = list(Chr8.prot.samples,
                                                      non.Chr8.prot.samples),
                                 gmt.features = list("msigdb_Homo sapiens_H"))
# GSEA3_positional <- panSEA::panSEA(list(global.w.mouse.df), "PDX Proteomics", 
#                                  DMEA = FALSE,
#                                  group.names = c("Chr8_amp","not_amp"),
#                                  group.samples = list(Chr8.prot.samples,
#                                                       non.Chr8.prot.samples),
#                                  gmt.features = list("msigdb_Homo sapiens_C1"))

msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")

# extract necessary info into data frame
msigdb.info <- as.data.frame(msigdb.info[, c(
  "gene_symbol",
  "gs_name",
  "gs_description"
)])

gmt <- DMEA::as_gmt(
  msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
  descriptions = "gs_description"
) 

Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                       "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                       "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")
gmt$genesets[[length(gmt$genesets)+1]] <- Chr8.cancer.genes
gmt$geneset.names[[length(gmt$geneset.names)+1]] <- "Chr8 cancer-associated genes"
gmt$geneset.descriptions[[length(gmt$geneset.descriptions)+1]] <- "Chr8 cancer-associated genes"

GSEA3_positional_custom <- panSEA::panSEA(list(global.w.mouse.df), "PDX Proteomics", 
                                   DMEA = FALSE,
                                   group.names = c("Chr8_amp","not_amp"),
                                   group.samples = list(Chr8.prot.samples,
                                                        non.Chr8.prot.samples),
                                   gmt.features = list(gmt))

GSEA.files <- list("GSEA_results.csv" =
                       GSEA3$mGSEA.results$compiled.results$results,
                     "GSEA_mean_results.csv" =
                       GSEA3$mGSEA.results$compiled.results$mean.results,
                     "GSEA_correlation_matrix.pdf" =
                       GSEA3$mGSEA.results$compiled.results$corr.matrix,
                     "GSEA_venn_diagram.pdf" =
                       GSEA3$mGSEA.results$compiled.results$venn.diagram,
                     "GSEA_dot_plot.pdf" =
                       GSEA3$mGSEA.results$compiled.results$dot.plot,
                     "GSEA_interactive_network.graph.html" =
                       GSEA3$mGSEA.network$interactive,
                     "GSEA_PDX_proteomics_volcano_plot.pdf" =
                       GSEA3$mGSEA.results$all.results$'PDX Proteomics'$volcano.plot,
                     "GSEA_phospho-proteomics_volcano_plot.pdf" =
                       GSEA3$mGSEA.results$all.results$'Phospho-proteomics'$volcano.plot)

GSEA.H.files <- list("GSEA_PDX_proteomics_results.csv" =
                        GSEA3_hallmark$mGSEA.results$all.results$'PDX Proteomics'$result,
                      "GSEA_interactive_network.graph.html" =
                        GSEA3_hallmark$mGSEA.network$interactive,
                      "GSEA_PDX_proteomics_volcano_plot.pdf" =
                        GSEA3_hallmark$mGSEA.results$all.results$'PDX Proteomics'$volcano.plot)

# GSEA.C1.files <- list("GSEA_PDX_proteomics_results.csv" =
#                         GSEA3_positional$mGSEA.results$all.results$'PDX Proteomics'$result,
#                       "GSEA_interactive_network.graph.html" =
#                         GSEA3_positional$mGSEA.network$interactive,
#                       "GSEA_PDX_proteomics_volcano_plot.pdf" =
#                         GSEA3_positional$mGSEA.results$all.results$'PDX Proteomics'$volcano.plot)

# GSEA.C1.custom.files <- list("GSEA_PDX_proteomics_results.csv" =
#                                GSEA3_positional_custom$mGSEA.results$all.results$'PDX Proteomics'$result,
#                              "GSEA_interactive_network.graph.html" =
#                                GSEA3_positional_custom$mGSEA.network$interactive,
#                              "GSEA_PDX_proteomics_volcano_plot.pdf" =
#                                GSEA3_positional_custom$mGSEA.results$all.results$'PDX Proteomics'$volcano.plot,
#                              "GSEA_PDX_proteomics_mountain_plot_Chr8_cancer-associated_genes.pdf" =
#                                GSEA3_positional_custom$mGSEA.results$all.results$'PDX Proteomics'$mtn.plots[["Chr8 cancer-associated genes"]])

# try using CCLE prot for prot DMEA
#### proteomics (CCLE from IMPROVE)
download.file("https://figshare.com/ndownloader/files/41466702", "proteomics.csv.gz")
prot.df <- read.csv(gzfile("proteomics.csv.gz"),fileEncoding="UTF-16LE")

allgenes = readr::read_csv("https://figshare.com/ndownloader/files/40576109")
genes = allgenes|>
  dplyr::select(gene_symbol,entrez_id)|>
  dplyr::distinct()
#genes <- genes[genes$gene_symbol %in% colnames(RNA.df)[2:ncol(RNA.df)], ]

allsamples = readr::read_csv('https://figshare.com/ndownloader/files/40576103')
CCLE.samples <- dplyr::distinct(allsamples[allsamples$id_source == "CCLE",
                                           c("other_id","improve_sample_id")])

# merge prot.df with genes, samples to stop using improve IDs
prot.df <- merge(prot.df, CCLE.samples)
prot.df <- merge(prot.df, genes)
prot.df$entrez_id <- NULL
prot.df <- dplyr::distinct(prot.df)

# convert to wide format for DMEA
prot.df <- reshape2::dcast(prot.df, other_id ~ gene_symbol, mean,
                           value.var = "proteomics")
colnames(prot.df)[1] <- "CCLE_ID"
prot.df.noNA <- prot.df[, colSums(is.na(prot.df)) == 0] # 23304 gene names
# DMEA
data.list <- list(global.df, global.w.mouse.df)
types <- c("Proteomics", "PDX Proteomics")
Chr8.prot.samples <- paste0("X", meta.df[meta.df$Chr8_amp & 
                                           meta.df$Sample == "Sample", ]$Tube)
non.Chr8.prot.samples <- paste0("X", meta.df[!meta.df$Chr8_amp & 
                                               meta.df$Sample == "Sample", ]$Tube)
DMEA3_prot <- panSEA::panSEA(data.list, types, GSEA = FALSE, 
                        group.names = c("Chr8_amp","not_amp"),
                        group.samples = list(Chr8.prot.samples,
                                             non.Chr8.prot.samples),
                        expression = list(prot.df.noNA, prot.df.noNA))

DEG.files <- list("DEG_results.csv" =
                    DMEA3$mDEG.results$compiled.results$results,
                  "DEG_mean_results.csv" =
                    DMEA3$mDEG.results$compiled.results$mean.results,
                  "DEG_correlation_matrix.pdf" =
                    DMEA3$mDEG.results$compiled.results$corr.matrix,
                  "DEG_venn_diagram.pdf" =
                    DMEA3$mDEG.results$compiled.results$venn.diagram,
                  "DEG_dot_plot.pdf" =
                    DMEA3$mDEG.results$compiled.results$dot.plot,
                  "DEG_results_phospho-proteomics.csv" = 
                    GSEA3$mDEG.results$all.results$`Phospho-proteomics`,
                  "DEG_dot_plot_phospho-proteomics.pdf" = 
                    GSEA3$mDEG.results$compiled.results$dot.plot)
DMEA.files <- list("DMEA_results.csv" =
                     DMEA3$mDMEA.results$compiled.results$results,
                   "DMEA_mean_results.csv" =
                     DMEA3$mDMEA.results$compiled.results$mean.results,
                   "DMEA_correlation_matrix.pdf" =
                     DMEA3$mDMEA.results$compiled.results$corr.matrix,
                   "DMEA_venn_diagram.pdf" =
                     DMEA3$mDMEA.results$compiled.results$venn.diagram,
                   "DMEA_dot_plot.pdf" =
                     DMEA3$mDMEA.results$compiled.results$dot.plot,
                   "DMEA_interactive_network.graph.html" =
                     DMEA3$mDMEA.network$interactive,
                   "DMEA_proteomics_volcano_plot.pdf" =
                     DMEA3$mDMEA.results$all.results$Proteomics$volcano.plot,
                   "DMEA_PDX_proteomics_volcano_plot.pdf" =
                     DMEA3$mDMEA.results$all.results$'PDX Proteomics'$volcano.plot,
                   "DMEA_PDX_proteomics_mountain_plot_TGF_beta_receptor_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'PDX Proteomics'$mtn.plot[["TGF beta receptor inhibitor"]],
                   "DMEA_PDX_proteomics_mountain_plot_mTOR_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'PDX Proteomics'$mtn.plot[["mTOR inhibitor"]],
                   "DMEA_PDX_proteomics_mountain_plot_RAF_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'PDX Proteomics'$mtn.plot[["RAF inhibitor"]],
                   "DMEA_proteomics_mountain_plot_TGF_beta_receptor_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'Proteomics'$mtn.plot[["TGF beta receptor inhibitor"]],
                   "DMEA_proteomics_mountain_plot_mTOR_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'Proteomics'$mtn.plot[["mTOR inhibitor"]],
                   "DMEA_proteomics_mountain_plot_RAF_inhibitor.pdf" =
                     DMEA3$mDMEA.results$all.results$'Proteomics'$mtn.plot[["RAF inhibitor"]],
                   "DMEA_proteomics_correlations.csv" =
                     DMEA3$mDMEA.results$all.results$Proteomics$corr.result,
                   "DMEA_PDX_proteomics_correlations.csv" =
                     DMEA3$mDMEA.results$all.results$'PDX Proteomics'$corr.result)
DMEA.prot.files <- list("DMEA_results.csv" =
                     DMEA3_prot$mDMEA.results$compiled.results$results,
                   "DMEA_mean_results.csv" =
                     DMEA3_prot$mDMEA.results$compiled.results$mean.results,
                   "DMEA_correlation_matrix.pdf" =
                     DMEA3_prot$mDMEA.results$compiled.results$corr.matrix,
                   "DMEA_venn_diagram.pdf" =
                     DMEA3_prot$mDMEA.results$compiled.results$venn.diagram,
                   "DMEA_dot_plot.pdf" =
                     DMEA3_prot$mDMEA.results$compiled.results$dot.plot,
                   "DMEA_interactive_network.graph.html" =
                     DMEA3_prot$mDMEA.network$interactive,
                   "DMEA_proteomics_volcano_plot.pdf" =
                     DMEA3_prot$mDMEA.results$all.results$Proteomics$volcano.plot,
                   "DMEA_PDX_proteomics_volcano_plot.pdf" =
                     DMEA3_prot$mDMEA.results$all.results$'PDX Proteomics'$volcano.plot,
                   "DMEA_PDX_proteomics_mountain_plot_FGFR_inhibitor.pdf" =
                     DMEA3_prot$mDMEA.results$all.results$'PDX Proteomics'$mtn.plot[["FGFR inhibitor"]],
                   "DMEA_PDX_proteomics_mountain_plot_mTOR_inhibitor.pdf" =
                     DMEA3_prot$mDMEA.results$all.results$'PDX Proteomics'$mtn.plot[["mTOR inhibitor"]],
                   "DMEA_proteomics_mountain_plot_mTOR_inhibitor.pdf" =
                     DMEA3_prot$mDMEA.results$all.results$'Proteomics'$mtn.plot[["mTOR inhibitor"]],
                   "DMEA_proteomics_correlations.csv" =
                     DMEA3_prot$mDMEA.results$all.results$Proteomics$corr.result,
                   "DMEA_PDX_proteomics_correlations.csv" =
                     DMEA3_prot$mDMEA.results$all.results$'PDX Proteomics'$corr.result)

all.files <- list('DMEA' = DMEA.files,
                  'DMEA_proteomics' = DMEA.prot.files,
                  'GSEA' = GSEA.files,
                  'GSEA_hallmark' = GSEA.H.files,
                  'Differential expression' = DEG.files)

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
synapse_id_map <- c("syn53357947" = "global_with_mouse_data/")
k <- 1
for (i in 1:length(all.files)) {
  # create local folder for subset of results
  setwd(paste0(base.path, synapse_id_map[k], "analysis"))
  dir.create(names(all.files)[i])
  setwd(names(all.files)[i])
  
  # save results locally
  temp.files <- all.files[[i]]
  
  CSV.files <- names(temp.files)[grepl(".csv", names(temp.files))]
  for (j in 1:length(CSV.files)) {
    write.csv(temp.files[[CSV.files[j]]], CSV.files[j], row.names = FALSE)
  }
  
  PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files))]
  for (j in 1:length(PDF.files)) {
    # if (grepl("dot", PDF.files[j])) { # wide plot for dot plots
    #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
    #                   device = "pdf", width = 35)
    # } else {
    #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
    #                   device = "pdf")
    # }
    ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                    device = "pdf")
  }
  
  HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
  if (length(HTML.files) > 0) {
    for (j in 1:length(HTML.files)) {
      visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j]) 
    }
  }
  
  # create folder on Synpase for subset of results
  dataFolder <- 
    synapser::synStore(synapser::Folder(names(all.files)[i],
                                        parent = names(synapse_id_map)[k]))
  
  # upload results to Synapse
  CSVs <- lapply(as.list(CSV.files), synapser::File,
                 parent = dataFolder)
  lapply(CSVs, synapser::synStore)
  
  PDFs <- lapply(as.list(PDF.files), synapser::File,
                 parent = dataFolder)
  lapply(PDFs, synapser::synStore)
  
  if (length(HTML.files) > 0) {
    HTMLs <- lapply(HTML.files, synapser::File,
                    parent = dataFolder)
    lapply(HTMLs, synapser::synStore) 
  }
}

# RNAseq only
RNA.df3 <- plyr::ddply(RNA.df, .(Gene), summarize,
                       'WU-487' = mean(WU.487_pdx, na.rm = TRUE),
                       'MN-2' = mean(mn2_pdx, na.rm = TRUE),
                       'JH-2-055' = mean(X2.055_pdx, na.rm = TRUE),
                       'JH-2-079' = mean(X2.079_pdx, na.rm = TRUE),
                       'WU-561' = mean(WU.561_pdx, na.rm = TRUE))
RNA.chr8 <- c("JH-2-079", "WU-561", "MN-2")
RNA.non.chr8 <- c("JH-2-055", "WU-487")
library(Biobase)
RNApanSEA <- panSEA::panSEA(list(RNA.df3), "Transcriptomics", 
                            group.samples = list(RNA.chr8, RNA.non.chr8))

DEG.files <- list("DEG_results.csv" =
                    RNApanSEA$mDEG.results$all.results$Transcriptomics)
DMEA.files <- list("DMEA_results.csv" =
                     RNApanSEA$mDMEA.results$all.results$Transcriptomics$result,
                   "DMEA_interactive_network.graph.html" =
                     RNApanSEA$mDMEA.network$interactive,
                   "DMEA_transcriptomics_volcano_plot.pdf" =
                     RNApanSEA$mDMEA.results$all.results$Transcriptomics$volcano.plot,
                   "DMEA_transcriptomics_mountain_plot_HMGCR_inhibitor.pdf" =
                     RNApanSEA$mDMEA.results$all.results$'Transcriptomics'$mtn.plot[["HMGCR inhibitor"]],
                   "DMEA_transcriptomics_mountain_plot_protein_synthesis_inhibitor.pdf" =
                     RNApanSEA$mDMEA.results$all.results$'Transcriptomics'$mtn.plot[["protein synthesis inhibitor"]],
                   "DMEA_transcriptomics_scatter_plots.pdf" =
                     RNApanSEA$mDMEA.results$all.results$Transcriptomics$corr.scatter.plots,
                   "DMEA_transcriptomics_correlations.csv" =
                     RNApanSEA$mDMEA.results$all.results$Transcriptomics$corr.result)
GSEA.files <- list("GSEA_results.csv" =
                     RNApanSEA$mGSEA.results$all.results$Transcriptomics$result,
                   "GSEA_interactive_network.graph.html" =
                     RNApanSEA$mGSEA.network$interactive,
                   "GSEA_transcriptomics_volcano_plot.pdf" =
                     RNApanSEA$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                   "GSEA_transcriptomics_correlations.csv" =
                     RNApanSEA$mGSEA.results$all.results$Transcriptomics$corr.result)
RNApanSEA_H <- panSEA::panSEA(list(RNA.df3), "Transcriptomics", DMEA = FALSE,
                            group.samples = list(RNA.chr8, RNA.non.chr8),
                            gmt.features = list("msigdb_Homo sapiens_H"))

msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")

# extract necessary info into data frame
msigdb.info <- as.data.frame(msigdb.info[, c(
  "gene_symbol",
  "gs_name",
  "gs_description"
)])

gmt <- DMEA::as_gmt(
  msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
  descriptions = "gs_description"
) 

Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                       "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                       "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")
gmt$genesets[[length(gmt$genesets)+1]] <- Chr8.cancer.genes
gmt$geneset.names[[length(gmt$geneset.names)+1]] <- "Chr8 cancer-associated genes"
gmt$geneset.descriptions[[length(gmt$geneset.descriptions)+1]] <- "Chr8 cancer-associated genes"

RNApanSEA_pos_custom <- panSEA::panSEA(list(RNA.df3), "Transcriptomics", DMEA = FALSE,
                              group.samples = list(RNA.chr8, RNA.non.chr8),
                              gmt.features = list(gmt))

GSEA.H.files <- list("GSEA_results.csv" =
                     RNApanSEA_H$mGSEA.results$all.results$Transcriptomics$result,
                   "GSEA_interactive_network.graph.html" =
                     RNApanSEA_H$mGSEA.network$interactive,
                   "GSEA_transcriptomics_volcano_plot.pdf" =
                     RNApanSEA_H$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                   "GSEA_transcriptomics_correlations.csv" =
                     RNApanSEA_H$mGSEA.results$all.results$Transcriptomics$corr.result,
                   "GSEA_transcriptomics_mountain_plot_HALLMARK_MYC_TARGETS_V2.pdf" =
                     RNApanSEA_H$mGSEA.results$all.results$'Transcriptomics'$mtn.plot[["HALLMARK_MYC_TARGETS_V2"]],
                   "GSEA_transcriptomics_mountain_plot_HALLMARK_MYC_TARGETS.pdf" =
                     RNApanSEA_H$mGSEA.results$all.results$'Transcriptomics'$mtn.plot[["HALLMARK_MYC_TARGETS"]])
GSEA.pos.custom.files <- list("GSEA_results.csv" =
                     RNApanSEA_pos_custom$mGSEA.results$all.results$Transcriptomics$result,
                   "GSEA_interactive_network.graph.html" =
                     RNApanSEA_pos_custom$mGSEA.network$interactive,
                   "GSEA_transcriptomics_volcano_plot.pdf" =
                     RNApanSEA_pos_custom$mGSEA.results$all.results$Transcriptomics$volcano.plot,
                   "GSEA_transcriptomics_correlations.csv" =
                     RNApanSEA_pos_custom$mGSEA.results$all.results$Transcriptomics$corr.result)

all.files <- list('DMEA' = DMEA.files,
                  'GSEA' = GSEA.files,
                  'GSEA_HALLMARK' = GSEA.H.files,
                  'GSEA_POSITIONAL_CUSTOM' = GSEA.pos.custom.files,
                  'Differential expression' = DEG.files)

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/"
for (i in 1:length(all.files)) {
  # create local folder for subset of results
  setwd(paste0(base.path, "analysis"))
  dir.create(names(all.files)[i])
  setwd(names(all.files)[i])
  
  # save results locally
  temp.files <- all.files[[i]]
  
  CSV.files <- names(temp.files)[grepl(".csv", names(temp.files))]
  for (j in 1:length(CSV.files)) {
    write.csv(temp.files[[CSV.files[j]]], CSV.files[j], row.names = FALSE)
  }
  
  PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files))]
  for (j in 1:length(PDF.files)) {
    # if (grepl("dot", PDF.files[j])) { # wide plot for dot plots
    #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
    #                   device = "pdf", width = 35)
    # } else {
    #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
    #                   device = "pdf")
    # }
    ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                    device = "pdf")
  }
  
  HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
  if (length(HTML.files) > 0) {
    for (j in 1:length(HTML.files)) {
      visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j]) 
    }
  }
}
setwd(paste0(base.path, "analysis"))
saveRDS(RNApanSEA_pos_custom, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_positional_custom.rds")) # 8.5 GB for 7 contrasts
saveRDS(RNApanSEA_H, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_hallmark.rds")) # 8.5 GB for 7 contrasts
saveRDS(RNApanSEA, file=paste0("Chr8_", "RNAseq_PDX", "_panSEA_CCLE.rds")) # 8.5 GB for 7 contrasts

# look at individual genes
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

genesets <- list('Chr8' = chr8, 
                 'Chr8p' = chr8p,
                 'Chr8q' = chr8q,
                 'Chr8_cancer_genes' = Chr8.cancer.genes,
                 'Myc_targets' = myc.targets,
                 'Myc_targets_v1' = myc.targetsv1,
                 'Myc_targets_v2' = myc.targetsv2)

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

# run GSEA on KSEA results
ksea.gsea <- panSEA::mGSEA(list(phospho), types = "Phospho-proteomics", 
                           feature.names = "Feature_set", rank.var = "NES")
ksea.gsea.hallmark <- panSEA::mGSEA(list(phospho), 
                                    gmt = list("msigdb_Homo sapiens_H"),
                                    types = "Phospho-proteomics", 
                                    feature.names = "Feature_set", 
                                    rank.var = "NES")
ksea.gsea.pos <- panSEA::mGSEA(list(phospho), 
                                    gmt = list("msigdb_Homo sapiens_C1"),
                                    types = "Phospho-proteomics", 
                                    feature.names = "Feature_set", 
                                    rank.var = "NES")

msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")

# extract necessary info into data frame
msigdb.info <- as.data.frame(msigdb.info[, c(
  "gene_symbol",
  "gs_name",
  "gs_description"
)])

gmt <- DMEA::as_gmt(
  msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
  descriptions = "gs_description"
) 

Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                       "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                       "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")
gmt$genesets[[length(gmt$genesets)+1]] <- Chr8.cancer.genes
gmt$geneset.names[[length(gmt$geneset.names)+1]] <- "Chr8 cancer-associated genes"
gmt$geneset.descriptions[[length(gmt$geneset.descriptions)+1]] <- "Chr8 cancer-associated genes"

ksea.gsea.pos.custom <- panSEA::mGSEA(list(phospho), 
                               gmt = list(gmt),
                               types = "Phospho-proteomics", 
                               feature.names = "Feature_set", 
                               rank.var = "NES")

ggplot2::ggsave("Phospho_KSEA_GSEA_mtn_plot_KEGG_JAK_STAT_SIGNALING_PATHWAY.pdf", 
       ksea.gsea$all.results$`Phospho-proteomics`$mtn.plots$KEGG_JAK_STAT_SIGNALING_PATHWAY)
ggplot2::ggsave("Phospho_KSEA_GSEA_mtn_plot_KEGG_TGF_BETA_SIGNALING_PATHWAY.pdf", 
                ksea.gsea$all.results$`Phospho-proteomics`$mtn.plots$KEGG_TGF_BETA_SIGNALING_PATHWAY)
saveRDS(ksea.gsea, "Phospho_KSEA_GSEA.rds")

# Transcriptomics PCA
setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/"
)
RNA.df <- read.table("salmon.merged.gene_counts.tsv", sep = "\t", 
                     header = TRUE)
rownames(RNA.df) <- RNA.df$gene_id
RNA.df$gene_id <- NULL
RNA.df$gene_name <- NULL

# create metadata for RNA-Seq
Tube <- colnames(RNA.df)
meta.df <- as.data.frame(Tube)
meta.df$SampleID <- c("JH-2-009", "JH-2-055", "JH-2-079", "WU-487", 
                      "WU-536", "WU-561", "MN-2", "MN-3")
meta.df$Sample <- "Sample"
meta.df$SampleID <- meta.df$Tube
rownames(meta.df) <- meta.df$Tube

library(MSnSet.utils)

RNA.mat <- RNA.df %>% as.matrix()
m_RNA <- MSnSet(exprs = RNA.df %>% as.matrix(),
                pData = meta.df)
plot_pca(m_RNA) + ggtitle("RNA-Seq PCA")
plot_pca(m_RNA, phenotype = "SampleID") + ggtitle("RNA-Seq PCA")
ggsave("RNAseq_PDX_PCA.pdf")

meta.df$SampleID <- c("JH-2-009", "JH-2-055", "JH-2-079", "WU-487", 
                      "WU-536", "WU-561", "MN-2", "MN-3")
m_RNA <- MSnSet(exprs = RNA.df %>% as.matrix(),
                pData = meta.df)
plot_pca(m_RNA, phenotype = "SampleID") + ggtitle("RNA-Seq PCA")
ggsave("RNAseq_PDX_PCA_betterSampleID.pdf")

Chr8.RNA.df <- RNA.df[, c("X2.055_pdx", "X2.079_pdx", "WU.487_pdx", "WU.561_pdx", "mn2_pdx")]
# create metadata for RNA-Seq
Tube <- colnames(Chr8.RNA.df)
meta.df <- as.data.frame(Tube)
meta.df$SampleID <- c("JH-2-055", "JH-2-079", "WU-487", "WU-561", "MN-2")
meta.df$Sample <- "Sample"
rownames(meta.df) <- meta.df$Tube
meta.df$Chr8_status <- "Amplified"
non.chr8.amp.rna <- c("JH-2-055", "WU-487")
meta.df[meta.df$SampleID %in% non.chr8.amp.rna, ]$Chr8_status <- "Not amplified"

m_RNA <- MSnSet(exprs = Chr8.RNA.df %>% as.matrix(),
                pData = meta.df)
plot_pca(m_RNA, phenotype = "SampleID") + ggtitle("RNA-Seq PCA")
ggsave("RNAseq_PDX_PCA_betterSampleID_knownChr8StatusOnly.pdf")

m_RNA <- MSnSet(exprs = Chr8.RNA.df %>% as.matrix(),
                pData = meta.df)
plot_pca(m_RNA, phenotype = "Chr8_status") + ggtitle("RNA-Seq PCA")
ggsave("RNAseq_PDX_PCA_byChr8Status_knownChr8StatusOnly.pdf")

overlap.RNA.df <- RNA.df[, c("X2.055_pdx", "X2.079_pdx", "WU.487_pdx", "mn2_pdx")]
# create metadata for RNA-Seq
Tube <- colnames(overlap.RNA.df)
meta.df <- as.data.frame(Tube)
meta.df$SampleID <- c("JH-2-055", "JH-2-079", "WU-487", "MN-2")
meta.df$Sample <- "Sample"
rownames(meta.df) <- meta.df$Tube
meta.df$Chr8_status <- "Amplified"
non.chr8.amp.rna <- c("JH-2-055", "WU-487")
meta.df[meta.df$SampleID %in% non.chr8.amp.rna, ]$Chr8_status <- "Not amplified"

m_RNA <- MSnSet(exprs = overlap.RNA.df %>% as.matrix(),
                pData = meta.df)
plot_pca(m_RNA, phenotype = "SampleID") + ggtitle("RNA-Seq PCA")
ggsave("RNAseq_PDX_PCA_betterSampleID_overlappingKnownChr8StatusOnly.pdf")

plot_pca(m_RNA, phenotype = "Chr8_status") + ggtitle("RNA-Seq PCA")
ggsave("RNAseq_PDX_PCA_byChr8Status_overlappingKnownChr8StatusOnly.pdf")

# more proteomics PCA
setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
)
meta.df <- readxl::read_excel("Chr8-MetaDataSheet.xlsx")
global.df <- read.table(
  "global_data/Chr8_crosstab_global_gene_corrected.txt",
  sep = "\t")
global.w.mouse.df <- read.table(
  "global_with_mouse_data/Chr8_crosstab_global_gene_corrected.txt", 
  sep = "\t")
phospho.df <- read.table(
  "phospho_data/Chr8_crosstab_phospho_SiteID_corrected.txt",
  sep = "\t")

# fix column names for MSnSet.utils
colnames(global.df) <- seq(1, ncol(global.df))
colnames(global.w.mouse.df) <- seq(1, ncol(global.w.mouse.df))
colnames(phospho.df) <- seq(1, ncol(phospho.df))

# add Tube number to metadata
meta.df$Tube <- NA
tube.num <- 0
for (i in 1:nrow(meta.df)) {
  if (meta.df$SampleID[i] != "reference") {
    tube.num <- tube.num + 1
    meta.df$Tube[i] <- tube.num
  }
}

# add Chr8 status to metadata
meta.df$Chr8_status <- "Amplified"
meta.df[meta.df$Sample == "Reference" | 
          meta.df$SampleID %in% c("JH-2-055", "WU-487CPP1"), ]$Chr8_status <- "Not amplified"

# identify global samples in RNAseq
meta.df$RNAseq <- FALSE
prot.samples.in.RNAseq <- c("WU-487CPP1", "MN-2", "JH-2-055", "JH-2-079")
meta.df[meta.df$SampleID %in% prot.samples.in.RNAseq, ]$RNAseq <- TRUE

# fix sampleIDs
meta.df[meta.df$SampleID == "WU-487CPP1", ]$SampleID <- "WU-487"
meta.df[meta.df$SampleID == "WU-22t Cshim81", ]$SampleID <- "WU-225"
meta.df <- meta.df[meta.df$Sample == "Sample", ]
meta.df <- data.frame(meta.df)

library(MSnSet.utils); library(plyr)
library(dplyr); library(ggplot2)
## global
# all samples
global.set <- MSnSet(exprs = global.df %>% as.matrix(),
                pData = meta.df)
plot_pca(global.set, phenotype = "SampleID") + ggtitle("Global PCA")
ggsave("Global_PCA_betterSampleID_KnownChr8StatusOnly.pdf")

plot_pca(global.set, phenotype = "Chr8_status") + ggtitle("Global PCA")
ggsave("Global_PCA_byChr8Status_KnownChr8StatusOnly.pdf")

# overlapping samples w RNAseq
global.overlap <- global.df[ , meta.df[meta.df$RNAseq, ]$Tube]
global.set <- MSnSet(exprs = global.overlap %>% as.matrix(),
                     pData = meta.df[meta.df$RNAseq, ])
plot_pca(global.set, phenotype = "SampleID") + ggtitle("Global PCA")
ggsave("Global_PCA_betterSampleID_overlappingKnownChr8StatusOnly.pdf")

plot_pca(global.set, phenotype = "Chr8_status") + ggtitle("Global PCA")
ggsave("Global_PCA_byChr8Status_overlappingKnownChr8StatusOnly.pdf")

## global PDX
# all samples
global.w.mouse.set <- MSnSet(exprs = global.w.mouse.df %>% as.matrix(),
                     pData = meta.df)
plot_pca(global.w.mouse.set, phenotype = "SampleID") + ggtitle("Global PCA")
ggsave("Global_PDX_PCA_betterSampleID_KnownChr8StatusOnly.pdf")

plot_pca(global.w.mouse.set, phenotype = "Chr8_status") + ggtitle("Global PCA")
ggsave("Global_PDX_PCA_byChr8Status_KnownChr8StatusOnly.pdf")

# overlapping samples w RNAseq
global.w.mouse.overlap <- global.w.mouse.df[ , meta.df[meta.df$RNAseq, ]$Tube]
global.w.mouse.set <- MSnSet(exprs = global.w.mouse.overlap %>% as.matrix(),
                     pData = meta.df[meta.df$RNAseq, ])
plot_pca(global.w.mouse.set, phenotype = "SampleID") + ggtitle("Global PCA")
ggsave("Global_PDX_PCA_betterSampleID_overlappingKnownChr8StatusOnly.pdf")

plot_pca(global.w.mouse.set, phenotype = "Chr8_status") + ggtitle("Global PCA")
ggsave("Global_PDX_PCA_byChr8Status_overlappingKnownChr8StatusOnly.pdf")

## phospho
# all samples
phospho.set <- MSnSet(exprs = phospho.df %>% as.matrix(),
                     pData = meta.df)
plot_pca(phospho.set, phenotype = "SampleID") + ggtitle("Phospho PCA")
ggsave("Phospho_PDX_PCA_betterSampleID_KnownChr8StatusOnly.pdf")

plot_pca(phospho.set, phenotype = "Chr8_status") + ggtitle("Phospho PCA")
ggsave("Phospho_PDX_PCA_byChr8Status_KnownChr8StatusOnly.pdf")

# overlapping samples w RNAseq
phospho.overlap <- phospho.df[ , meta.df[meta.df$RNAseq, ]$Tube]
phospho.set <- MSnSet(exprs = phospho.overlap %>% as.matrix(),
                     pData = meta.df[meta.df$RNAseq, ])
plot_pca(phospho.set, phenotype = "SampleID") + ggtitle("Phospho PCA")
ggsave("Phospho_PDX_PCA_betterSampleID_overlappingKnownChr8StatusOnly.pdf")

plot_pca(phospho.set, phenotype = "Chr8_status") + ggtitle("Phospho PCA")
ggsave("Phospho_PDX_PCA_byChr8Status_overlappingKnownChr8StatusOnly.pdf")

# redo phospho w only samples overlapping w RNAseq
setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
)
meta.df <- readxl::read_excel("Chr8-MetaDataSheet.xlsx")
phospho.df <- read.table(
  "phospho_data/Chr8_crosstab_phospho_SiteID_corrected.txt",
  sep = "\t")

# add column for feature names and later make it the first column
phospho.df$SUB_SITE <- rownames(phospho.df)

# add Tube number to metadata
meta.df$Tube <- NA
tube.num <- 0
for (i in 1:nrow(meta.df)) {
  if (meta.df$SampleID[i] != "reference") {
    tube.num <- tube.num + 1
    meta.df$Tube[i] <- tube.num
  }
}

# add Chr8 status to metadata
meta.df$Chr8_amp <- TRUE
meta.df[meta.df$Sample == "Reference" | 
          meta.df$SampleID %in% c("JH-2-055", "WU-487CPP1"), ]$Chr8_amp <- FALSE

# identify global samples in RNAseq
meta.df$RNAseq <- FALSE
prot.samples.in.RNAseq <- c("WU-487CPP1", "MN-2", "JH-2-055", "JH-2-079")
meta.df[meta.df$SampleID %in% prot.samples.in.RNAseq, ]$RNAseq <- TRUE
prot.tubes.in.RNAseq <- paste0("X", meta.df[meta.df$RNAseq, ]$Tube)
phospho.overlap <- phospho.df[ , c("SUB_SITE", prot.tubes.in.RNAseq)]

Chr8.prot.samples <- paste0("X", meta.df[meta.df$Chr8_amp & meta.df$RNAseq &
                                           meta.df$Sample == "Sample", ]$Tube)
non.Chr8.prot.samples <- paste0("X", meta.df[!meta.df$Chr8_amp & meta.df$RNAseq &
                                               meta.df$Sample == "Sample", ]$Tube)

data.list <- list(phospho.overlap)
types <- c("Phospho-proteomics")
gmt.ksea <- readRDS(paste0(
  "phospho_data/gmt_PNNL_kinase-substrate_Chr8_5344.rds"))
KSEA <- panSEA::panSEA(data.list, types, DMEA = FALSE, 
                        feature.names = c("SUB_SITE"),
                        group.names = c("Chr8_amp","not_amp"),
                        group.samples = list(Chr8.prot.samples,
                                             non.Chr8.prot.samples),
                        gmt.features = list(gmt.ksea))

KSEA.result <- KSEA$mGSEA.results$all.results$'Phospho-proteomics'$result
data.list <- list('Phospho-proteomics' = KSEA.result[ , c("Feature_set", "NES")])
types <- "Phospho-proteomics"

KSEA.GSEA <- panSEA::mGSEA(data.list, types = types, 
                        feature.names = "Feature_set",
                        rank.var = "NES")
# KSEA.GSEA2 <- panSEA::panSEA(data.list, types = types, 
#                            feature.names = "Feature_set",
#                            GSEA.rank.var = "NES", group.names = "Chr8_amp_vs_not_amp",
#                            group.samples = list("NES"))
KSEA.GSEA2 <- panSEA::panSEA(data.list, types = types,
                           feature.names = "Feature_set",
                           GSEA.rank.var = "NES", group.names = "Chr8_amp_vs_not_amp",
                           group.samples = list(2))

KSEA.GSEA.H <- panSEA::mGSEA(data.list, types = types, 
                            feature.names = "Feature_set",
                            rank.var = "NES",
                            gmt.features = list("msigdb_Homo sapiens_H"))
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")

# extract necessary info into data frame
msigdb.info <- as.data.frame(msigdb.info[, c(
  "gene_symbol",
  "gs_name",
  "gs_description"
)])

gmt <- DMEA::as_gmt(
  msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
  descriptions = "gs_description"
) 

Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                       "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                       "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")
gmt$genesets[[length(gmt$genesets)+1]] <- Chr8.cancer.genes
gmt$geneset.names[[length(gmt$geneset.names)+1]] <- "Chr8 cancer-associated genes"
gmt$geneset.descriptions[[length(gmt$geneset.descriptions)+1]] <- "Chr8 cancer-associated genes"
KSEA.GSEA.C1.custom <- panSEA::mGSEA(data.list, types = types, 
                                     feature.names = "Feature_set",
                                     rank.var = "NES",
                                     gmt = list(gmt))
KSEA.GSEA.C1 <- panSEA::mGSEA(data.list, types = types, 
                             feature.names = "Feature_set",
                             rank.var = "NES",
                             gmt = list("msigdb_Homo sapiens_C1"))
DEG.files <- list("DEG_results.csv" =
                    KSEA$mDEG.results$all.results$'Phospho-proteomics')
KSEA.files <- list("KSEA_results.csv" =
                       KSEA$mGSEA.results$all.results$'Phospho-proteomics'$result,
                     "KSEA_interactive_network.graph.html" =
                       KSEA$mGSEA.network$interactive,
                     "KSEA_phospho_volcano_plot.pdf" =
                       KSEA$mGSEA.results$all.results$'Phospho-proteomics'$volcano.plot,
                     "KSEA_phospho_mountain_plot_CDK17.pdf" =
                       KSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot[["CDK17"]],
                     "KSEA_phospho_mountain_plot_CDK18.pdf" =
                       KSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot[["CDK18"]]
                   # ,
                   # "KSEA_phospho_mountain_plots.rds" =
                   #   KSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot
                   )

KSEA.mtn.plots <- KSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot

DMEA.files <- list("DMEA_results.csv" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$result,
                   "DMEA_interactive_network.graph.html" =
                     KSEA$mDMEA.network$interactive,
                   "DMEA_phospho_volcano_plot.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$volcano.plot,
                   "DMEA_phospho_mountain_plot_HMGCR_inhibitor.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$mtn.plot[["HMGCR inhibitor"]],
                   "DMEA_phospho_mountain_plot_protein_synthesis_inhibitor.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$mtn.plot[["protein synthesis inhibitor"]],
                   "DMEA_phospho_scatter_plots.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$corr.scatter.plots,
                   "DMEA_phospho_correlations.csv" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$corr.result)
GSEA.files <- list("GSEA_results.csv" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$result,
                   "GSEA_interactive_network.graph.html" =
                     KSEA.GSEA$mGSEA.network$interactive,
                   "GSEA_phospho_volcano_plot.pdf" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$volcano.plot,
                   "GSEA_phospho_mountain_plot_KEGG_JAK_STAT_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_JAK_STAT_SIGNALING_PATHWAY"]],
                   "GSEA_phospho_mountain_plot_KEGG_TGF_BETA_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_TGF_BETA_SIGNALING_PATHWAY"]],
                   "GSEA_results.rds" = 
                     KSEA.GSEA)

GSEA.files <- list("GSEA_results.csv" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$result,
                   "GSEA_phospho_volcano_plot.pdf" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$volcano.plot,
                   "GSEA_phospho_mountain_plot_KEGG_JAK_STAT_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_JAK_STAT_SIGNALING_PATHWAY"]],
                   "GSEA_phospho_mountain_plot_KEGG_TGF_BETA_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_TGF_BETA_SIGNALING_PATHWAY"]],
                   "GSEA_results.rds" = 
                     KSEA.GSEA)
GSEA.H.files <- list("GSEA_results.csv" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$result,
                     "GSEA_phospho_volcano_plot.pdf" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$volcano.plot,
                     "GSEA_phospho_mountain_plot_HALLMARK_IL2_STAT5_SIGNALING.pdf" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$mtn.plot[["HALLMARK_IL2_STAT5_SIGNALING"]],
                     "GSEA_phospho_mountain_plot_HALLMARK_HEDGEHOG_SIGNALING.pdf" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$mtn.plot[["HALLMARK_HEDGEHOG_SIGNALING"]])
GSEA.pos.custom.files <- list("GSEA_results.csv" =
                                KSEA.GSEA.C1.custom$all.results$'Phospho-proteomics'$result,
                              "GSEA_phospho_volcano_plot.pdf" =
                                KSEA.GSEA.C1.custom$all.results$'Phospho-proteomics'$volcano.plot)
GSEA.pos.files <- list("GSEA_results.csv" =
                                KSEA.GSEA.C1$all.results$'Phospho-proteomics'$result,
                              "GSEA_phospho_volcano_plot.pdf" =
                                KSEA.GSEA.C1$all.results$'Phospho-proteomics'$volcano.plot)

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/phospho_data/analysis/overlap_w_PDX_RNAseq/"
all.files <- list('Differential expression' = DEG.files,
                  'KSEA' = KSEA.files)
all.files <- list('GSEA' = GSEA.files,
                  'GSEA_HALLMARK' = GSEA.H.files,
                  'GSEA_POSITIONAL' = GSEA.pos.files,
                  'GSEA_POSITIONAL_CUSTOM' = GSEA.pos.custom.files)
for (i in 1:length(all.files)) {
  # create local folder for subset of results
  setwd(paste0(base.path))
  dir.create(names(all.files)[i])
  setwd(names(all.files)[i])
  
  # save results locally
  temp.files <- all.files[[i]]
  
  CSV.files <- names(temp.files)[grepl(".csv", names(temp.files))]
  for (j in 1:length(CSV.files)) {
    write.csv(temp.files[[CSV.files[j]]], CSV.files[j], row.names = FALSE)
  }
  
  PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files))]
  if (length(PDF.files) > 0) {
    for (j in 1:length(PDF.files)) {
      # if (grepl("dot", PDF.files[j])) { # wide plot for dot plots
      #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
      #                   device = "pdf", width = 35)
      # } else {
      #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
      #                   device = "pdf")
      # }
      ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                      device = "pdf")
    } 
  }
  
  HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
  if (length(HTML.files) > 0) {
    for (j in 1:length(HTML.files)) {
      visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j]) 
    }
  }
  
  RDS.files <- names(temp.files)[grepl(".rds", names(temp.files))]
  if (length(RDS.files) > 0) {
    for (j in 1:length(RDS.files)) {
      saveRDS(temp.files[[RDS.files[j]]], RDS.files[j]) 
    }
  }
}
setwd(paste0(base.path, "analysis"))
saveRDS(RNApanSEA_pos_custom, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_positional_custom.rds")) # 8.5 GB for 7 contrasts
saveRDS(RNApanSEA_H, file=paste0("Chr8_", "RNAseq_PDX", "_GSEA_hallmark.rds")) # 8.5 GB for 7 contrasts
saveRDS(RNApanSEA, file=paste0("Chr8_", "RNAseq_PDX", "_panSEA_CCLE.rds")) # 8.5 GB for 7 contrasts

# heatmaps of gene sets
# look at individual genes
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

genesets <- list('Chr8' = chr8, 
                 'Chr8p' = chr8p,
                 'Chr8q' = chr8q,
                 'Chr8_cancer_genes' = Chr8.cancer.genes,
                 'Myc_targets' = myc.targets,
                 'Myc_targets_v1' = myc.targetsv1,
                 'Myc_targets_v2' = myc.targetsv2)

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

# try true ksea
library(plyr); library(dplyr); library(tidyr)
ksdb <- read.csv(paste0("https://raw.githubusercontent.com/BelindaBGarana/",
                        "panSEA/shiny-app/data/ksdb_20231101.csv"))
ksdb.human <- ksdb[
  ksdb$KIN_ORGANISM == "human" & ksdb$SUB_ORGANISM == "human", ]
ksdb <- NULL # save space
# ksdb.human$SUB_SITE <- paste(ksdb.human$SUBSTRATE, ksdb.human$SUB_MOD_RSD,
#                               collapse = "_")
ksdb.human$SUB_MOD_RSD_PNNL <- ksdb.human$SUB_MOD_RSD
ksdb.human$SUB_MOD_RSD_PNNL <- paste0(ksdb.human$SUB_MOD_RSD_PNNL, tolower(substr(ksdb.human$SUB_MOD_RSD_PNNL,1,1)))
ksdb.human <- ksdb.human %>% unite("SUB_SITE",
                                   c("SUBSTRATE", "SUB_MOD_RSD_PNNL"),
                                   sep = "-", remove = FALSE)
write.csv(ksdb.human, "ksdb_human_PNNL.csv")
gmt <- DMEA::as_gmt(ksdb.human, "SUB_SITE", "KINASE", min.per.set = 6,
                    descriptions = "KIN_ACC_ID")
saveRDS(gmt, "gmt_ksdb_human_PNNL.rds")

setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
)
meta.df <- readxl::read_excel("Chr8-MetaDataSheet.xlsx")
phospho.df <- read.table(
  "phospho_data/Chr8_crosstab_phospho_SiteID_corrected.txt",
  sep = "\t")

# add column for feature names and later make it the first column
phospho.df$SUB_SITE <- rownames(phospho.df)

# add Tube number to metadata
meta.df$Tube <- NA
tube.num <- 0
for (i in 1:nrow(meta.df)) {
  if (meta.df$SampleID[i] != "reference") {
    tube.num <- tube.num + 1
    meta.df$Tube[i] <- tube.num
  }
}

# add Chr8 status to metadata
meta.df$Chr8_amp <- TRUE
meta.df[meta.df$Sample == "Reference" | 
          meta.df$SampleID %in% c("JH-2-055", "WU-487CPP1"), ]$Chr8_amp <- FALSE

# identify global samples in RNAseq
meta.df$RNAseq <- FALSE
prot.samples.in.RNAseq <- c("WU-487CPP1", "MN-2", "JH-2-055", "JH-2-079")
meta.df[meta.df$SampleID %in% prot.samples.in.RNAseq, ]$RNAseq <- TRUE
prot.tubes.in.RNAseq <- paste0("X", meta.df[meta.df$RNAseq, ]$Tube)
phospho.overlap <- phospho.df[ , c("SUB_SITE", prot.tubes.in.RNAseq)]

Chr8.prot.samples <- paste0("X", meta.df[meta.df$Chr8_amp & meta.df$RNAseq &
                                           meta.df$Sample == "Sample", ]$Tube)
non.Chr8.prot.samples <- paste0("X", meta.df[!meta.df$Chr8_amp & meta.df$RNAseq &
                                               meta.df$Sample == "Sample", ]$Tube)

data.list <- list(phospho.overlap)
types <- c("Phospho-proteomics")
KSEA <- panSEA::panSEA(data.list, types, DMEA = FALSE, 
                       feature.names = c("SUB_SITE"),
                       group.names = c("Chr8_amp","not_amp"),
                       group.samples = list(Chr8.prot.samples,
                                            non.Chr8.prot.samples),
                       gmt.features = list(gmt))

KSEA.result <- KSEA$mGSEA.results$all.results$'Phospho-proteomics'$result
data.list <- list('Phospho-proteomics' = KSEA.result[ , c("Feature_set", "NES")])
types <- "Phospho-proteomics"

KSEA.GSEA <- panSEA::mGSEA(data.list, types = types, 
                           feature.names = "Feature_set",
                           rank.var = "NES") # 2+ gene sets not covered
# KSEA.GSEA2 <- panSEA::panSEA(data.list, types = types, 
#                            feature.names = "Feature_set",
#                            GSEA.rank.var = "NES", group.names = "Chr8_amp_vs_not_amp",
#                            group.samples = list("NES"))
# KSEA.GSEA2 <- panSEA::panSEA(data.list, types = types,
#                              feature.names = "Feature_set",
#                              GSEA.rank.var = "NES", group.names = "Chr8_amp_vs_not_amp",
#                              group.samples = list(2))

KSEA.GSEA.H <- panSEA::mGSEA(data.list, types = types, 
                             feature.names = "Feature_set",
                             rank.var = "NES",
                             gmt = list("msigdb_Homo sapiens_H")) # 2+ gene sets not covered
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")

# extract necessary info into data frame
msigdb.info <- as.data.frame(msigdb.info[, c(
  "gene_symbol",
  "gs_name",
  "gs_description"
)])

gmt.C1 <- DMEA::as_gmt(
  msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
  descriptions = "gs_description"
) 

Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                       "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                       "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")
gmt.C1$genesets[[length(gmt.C1$genesets)+1]] <- Chr8.cancer.genes
gmt.C1$geneset.names[[length(gmt.C1$geneset.names)+1]] <- "Chr8 cancer-associated genes"
gmt.C1$geneset.descriptions[[length(gmt.C1$geneset.descriptions)+1]] <- "Chr8 cancer-associated genes"
KSEA.GSEA.C1.custom <- panSEA::mGSEA(data.list, types = types, 
                                     feature.names = "Feature_set",
                                     rank.var = "NES",
                                     gmt = list(gmt.C1)) # 2+ gene sets not covered
KSEA.GSEA.C1 <- panSEA::mGSEA(data.list, types = types, 
                              feature.names = "Feature_set",
                              rank.var = "NES",
                              gmt = list("msigdb_Homo sapiens_C1")) # 2+ gene sets not covered
DEG.files <- list("DEG_results.csv" =
                    KSEA$mDEG.results$all.results$'Phospho-proteomics')
KSEA.files <- list("KSEA_results.csv" =
                     KSEA$mGSEA.results$all.results$'Phospho-proteomics'$result,
                   "KSEA_interactive_network_graph.html" =
                     KSEA$mGSEA.network$interactive,
                   "KSEA_phospho_volcano_plot.pdf" =
                     KSEA$mGSEA.results$all.results$'Phospho-proteomics'$volcano.plot,
                   "KSEA_phospho_mountain_plot_MELK.pdf" =
                     KSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot[["MELK"]])

#KSEA.mtn.plots <- KSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot

DMEA.files <- list("DMEA_results.csv" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$result,
                   "DMEA_interactive_network.graph.html" =
                     KSEA$mDMEA.network$interactive,
                   "DMEA_phospho_volcano_plot.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$volcano.plot,
                   "DMEA_phospho_mountain_plot_HMGCR_inhibitor.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$mtn.plot[["HMGCR inhibitor"]],
                   "DMEA_phospho_mountain_plot_protein_synthesis_inhibitor.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$mtn.plot[["protein synthesis inhibitor"]],
                   "DMEA_phospho_scatter_plots.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$corr.scatter.plots,
                   "DMEA_phospho_correlations.csv" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$corr.result)
GSEA.files <- list("GSEA_results.csv" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$result,
                   "GSEA_interactive_network.graph.html" =
                     KSEA.GSEA$mGSEA.network$interactive,
                   "GSEA_phospho_volcano_plot.pdf" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$volcano.plot,
                   "GSEA_phospho_mountain_plot_KEGG_JAK_STAT_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_JAK_STAT_SIGNALING_PATHWAY"]],
                   "GSEA_phospho_mountain_plot_KEGG_TGF_BETA_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_TGF_BETA_SIGNALING_PATHWAY"]],
                   "GSEA_results.rds" = 
                     KSEA.GSEA)

GSEA.files <- list("GSEA_results.csv" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$result,
                   "GSEA_phospho_volcano_plot.pdf" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$volcano.plot,
                   "GSEA_phospho_mountain_plot_KEGG_JAK_STAT_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_JAK_STAT_SIGNALING_PATHWAY"]],
                   "GSEA_phospho_mountain_plot_KEGG_TGF_BETA_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_TGF_BETA_SIGNALING_PATHWAY"]],
                   "GSEA_results.rds" = 
                     KSEA.GSEA)
GSEA.H.files <- list("GSEA_results.csv" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$result,
                     "GSEA_phospho_volcano_plot.pdf" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$volcano.plot,
                     "GSEA_phospho_mountain_plot_HALLMARK_IL2_STAT5_SIGNALING.pdf" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$mtn.plot[["HALLMARK_IL2_STAT5_SIGNALING"]],
                     "GSEA_phospho_mountain_plot_HALLMARK_HEDGEHOG_SIGNALING.pdf" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$mtn.plot[["HALLMARK_HEDGEHOG_SIGNALING"]])
GSEA.pos.custom.files <- list("GSEA_results.csv" =
                                KSEA.GSEA.C1.custom$all.results$'Phospho-proteomics'$result,
                              "GSEA_phospho_volcano_plot.pdf" =
                                KSEA.GSEA.C1.custom$all.results$'Phospho-proteomics'$volcano.plot)
GSEA.pos.files <- list("GSEA_results.csv" =
                         KSEA.GSEA.C1$all.results$'Phospho-proteomics'$result,
                       "GSEA_phospho_volcano_plot.pdf" =
                         KSEA.GSEA.C1$all.results$'Phospho-proteomics'$volcano.plot)

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/phospho_data/analysis/overlap_w_PDX_RNAseq/"
setwd(base.path)
dir.create("ksdb")
all.files <- list('Differential expression' = DEG.files,
                  'KSEA' = KSEA.files)
all.files <- list('GSEA' = GSEA.files,
                  'GSEA_HALLMARK' = GSEA.H.files,
                  'GSEA_POSITIONAL' = GSEA.pos.files,
                  'GSEA_POSITIONAL_CUSTOM' = GSEA.pos.custom.files)
for (i in 1:length(all.files)) {
  # create local folder for subset of results
  setwd(paste0(base.path, "ksdb/"))
  dir.create(names(all.files)[i])
  setwd(names(all.files)[i])
  
  # save results locally
  temp.files <- all.files[[i]]
  
  CSV.files <- names(temp.files)[grepl(".csv", names(temp.files))]
  for (j in 1:length(CSV.files)) {
    write.csv(temp.files[[CSV.files[j]]], CSV.files[j], row.names = FALSE)
  }
  
  PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files))]
  if (length(PDF.files) > 0) {
    for (j in 1:length(PDF.files)) {
      # if (grepl("dot", PDF.files[j])) { # wide plot for dot plots
      #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
      #                   device = "pdf", width = 35)
      # } else {
      #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
      #                   device = "pdf")
      # }
      ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                      device = "pdf")
    } 
  }
  
  HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
  if (length(HTML.files) > 0) {
    for (j in 1:length(HTML.files)) {
      visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j]) 
    }
  }
  
  RDS.files <- names(temp.files)[grepl(".rds", names(temp.files))]
  if (length(RDS.files) > 0) {
    for (j in 1:length(RDS.files)) {
      saveRDS(temp.files[[RDS.files[j]]], RDS.files[j]) 
    }
  }
}

## repeat for all samples (not just overlap)
data.list <- list(phospho.df)
types <- c("Phospho-proteomics")

Chr8.prot.samples <- paste0("X", meta.df[meta.df$Chr8_amp &
                                           meta.df$Sample == "Sample", ]$Tube)
non.Chr8.prot.samples <- paste0("X", meta.df[!meta.df$Chr8_amp &
                                               meta.df$Sample == "Sample", ]$Tube)

KSEA <- panSEA::panSEA(data.list, types, DMEA = FALSE, 
                       feature.names = c("SUB_SITE"),
                       group.names = c("Chr8_amp","not_amp"),
                       group.samples = list(Chr8.prot.samples,
                                            non.Chr8.prot.samples),
                       gmt.features = list(gmt))

KSEA.result <- KSEA$mGSEA.results$all.results$'Phospho-proteomics'$result
data.list <- list('Phospho-proteomics' = KSEA.result[ , c("Feature_set", "NES")])
types <- "Phospho-proteomics"
#dev.print(pdf, "KSEA_static_network_graph.pdf") # no significant enrichments

KSEA.GSEA <- panSEA::mGSEA(data.list, types = types, 
                           feature.names = "Feature_set",
                           rank.var = "NES") # 2+ gene sets not covered
# KSEA.GSEA2 <- panSEA::panSEA(data.list, types = types, 
#                            feature.names = "Feature_set",
#                            GSEA.rank.var = "NES", group.names = "Chr8_amp_vs_not_amp",
#                            group.samples = list("NES"))
# KSEA.GSEA2 <- panSEA::panSEA(data.list, types = types,
#                              feature.names = "Feature_set",
#                              GSEA.rank.var = "NES", group.names = "Chr8_amp_vs_not_amp",
#                              group.samples = list(2))

KSEA.GSEA.H <- panSEA::mGSEA(data.list, types = types, 
                             feature.names = "Feature_set",
                             rank.var = "NES",
                             gmt = list("msigdb_Homo sapiens_H")) # 2+ gene sets not covered
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")

# extract necessary info into data frame
msigdb.info <- as.data.frame(msigdb.info[, c(
  "gene_symbol",
  "gs_name",
  "gs_description"
)])

gmt.C1 <- DMEA::as_gmt(
  msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
  descriptions = "gs_description"
) 

Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                       "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                       "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")
gmt.C1$genesets[[length(gmt.C1$genesets)+1]] <- Chr8.cancer.genes
gmt.C1$geneset.names[[length(gmt.C1$geneset.names)+1]] <- "Chr8 cancer-associated genes"
gmt.C1$geneset.descriptions[[length(gmt.C1$geneset.descriptions)+1]] <- "Chr8 cancer-associated genes"
KSEA.GSEA.C1.custom <- panSEA::mGSEA(data.list, types = types, 
                                     feature.names = "Feature_set",
                                     rank.var = "NES",
                                     gmt = list(gmt.C1)) # 2+ gene sets not covered
KSEA.GSEA.C1 <- panSEA::mGSEA(data.list, types = types, 
                              feature.names = "Feature_set",
                              rank.var = "NES",
                              gmt = list("msigdb_Homo sapiens_C1")) # 2+ gene sets not covered
DEG.files <- list("DEG_results.csv" =
                    KSEA$mDEG.results$all.results$'Phospho-proteomics')
KSEA.files <- list("KSEA_results.csv" =
                     KSEA$mGSEA.results$all.results$'Phospho-proteomics'$result,
                   "KSEA_phospho_volcano_plot.pdf" =
                     KSEA$mGSEA.results$all.results$'Phospho-proteomics'$volcano.plot)

#KSEA.mtn.plots <- KSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot

DMEA.files <- list("DMEA_results.csv" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$result,
                   "DMEA_interactive_network.graph.html" =
                     KSEA$mDMEA.network$interactive,
                   "DMEA_phospho_volcano_plot.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$volcano.plot,
                   "DMEA_phospho_mountain_plot_HMGCR_inhibitor.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$mtn.plot[["HMGCR inhibitor"]],
                   "DMEA_phospho_mountain_plot_protein_synthesis_inhibitor.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$mtn.plot[["protein synthesis inhibitor"]],
                   "DMEA_phospho_scatter_plots.pdf" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$corr.scatter.plots,
                   "DMEA_phospho_correlations.csv" =
                     KSEA$mDMEA.results$all.results$'Phospho-proteomics'$corr.result)
GSEA.files <- list("GSEA_results.csv" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$result,
                   "GSEA_interactive_network.graph.html" =
                     KSEA.GSEA$mGSEA.network$interactive,
                   "GSEA_phospho_volcano_plot.pdf" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$volcano.plot,
                   "GSEA_phospho_mountain_plot_KEGG_JAK_STAT_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_JAK_STAT_SIGNALING_PATHWAY"]],
                   "GSEA_phospho_mountain_plot_KEGG_TGF_BETA_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$mGSEA.results$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_TGF_BETA_SIGNALING_PATHWAY"]],
                   "GSEA_results.rds" = 
                     KSEA.GSEA)

GSEA.files <- list("GSEA_results.csv" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$result,
                   "GSEA_phospho_volcano_plot.pdf" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$volcano.plot,
                   "GSEA_phospho_mountain_plot_KEGG_JAK_STAT_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_JAK_STAT_SIGNALING_PATHWAY"]],
                   "GSEA_phospho_mountain_plot_KEGG_TGF_BETA_SIGNALING_PATHWAY.pdf" =
                     KSEA.GSEA$all.results$'Phospho-proteomics'$mtn.plot[["KEGG_TGF_BETA_SIGNALING_PATHWAY"]],
                   "GSEA_results.rds" = 
                     KSEA.GSEA)
GSEA.H.files <- list("GSEA_results.csv" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$result,
                     "GSEA_phospho_volcano_plot.pdf" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$volcano.plot,
                     "GSEA_phospho_mountain_plot_HALLMARK_IL2_STAT5_SIGNALING.pdf" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$mtn.plot[["HALLMARK_IL2_STAT5_SIGNALING"]],
                     "GSEA_phospho_mountain_plot_HALLMARK_HEDGEHOG_SIGNALING.pdf" =
                       KSEA.GSEA.H$all.results$'Phospho-proteomics'$mtn.plot[["HALLMARK_HEDGEHOG_SIGNALING"]])
GSEA.pos.custom.files <- list("GSEA_results.csv" =
                                KSEA.GSEA.C1.custom$all.results$'Phospho-proteomics'$result,
                              "GSEA_phospho_volcano_plot.pdf" =
                                KSEA.GSEA.C1.custom$all.results$'Phospho-proteomics'$volcano.plot)
GSEA.pos.files <- list("GSEA_results.csv" =
                         KSEA.GSEA.C1$all.results$'Phospho-proteomics'$result,
                       "GSEA_phospho_volcano_plot.pdf" =
                         KSEA.GSEA.C1$all.results$'Phospho-proteomics'$volcano.plot)

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/phospho_data/analysis/"
setwd(base.path)
dir.create("ksdb")
all.files <- list('Differential expression' = DEG.files,
                  'KSEA' = KSEA.files)
all.files <- list('GSEA' = GSEA.files,
                  'GSEA_HALLMARK' = GSEA.H.files,
                  'GSEA_POSITIONAL' = GSEA.pos.files,
                  'GSEA_POSITIONAL_CUSTOM' = GSEA.pos.custom.files)
for (i in 1:length(all.files)) {
  # create local folder for subset of results
  setwd(paste0(base.path, "ksdb/"))
  dir.create(names(all.files)[i])
  setwd(names(all.files)[i])
  
  # save results locally
  temp.files <- all.files[[i]]
  
  CSV.files <- names(temp.files)[grepl(".csv", names(temp.files))]
  for (j in 1:length(CSV.files)) {
    write.csv(temp.files[[CSV.files[j]]], CSV.files[j], row.names = FALSE)
  }
  
  PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files))]
  if (length(PDF.files) > 0) {
    for (j in 1:length(PDF.files)) {
      # if (grepl("dot", PDF.files[j])) { # wide plot for dot plots
      #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
      #                   device = "pdf", width = 35)
      # } else {
      #   ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
      #                   device = "pdf")
      # }
      ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                      device = "pdf")
    } 
  }
  
  HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
  if (length(HTML.files) > 0) {
    for (j in 1:length(HTML.files)) {
      visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j]) 
    }
  }
  
  RDS.files <- names(temp.files)[grepl(".rds", names(temp.files))]
  if (length(RDS.files) > 0) {
    for (j in 1:length(RDS.files)) {
      saveRDS(temp.files[[RDS.files[j]]], RDS.files[j]) 
    }
  }
}
