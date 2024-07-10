# Differential expression & enrichment analyses: PDX transcriptomics, global & phospho
# Chr8: amplified vs. not amplified
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-04-25

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
source("panSEA_helper_20240508.R")

#### 1. Import metadata & crosstabs ####
setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/"
)
meta.df <- readxl::read_excel("Chr8-MetaDataSheet.xlsx")

global.df <- read.table(
  "global_data/Chr8_crosstab_global_gene_corrected.txt",
  sep = "\t")
phospho.df <- read.table(
  "phospho_data/Chr8_crosstab_phospho_SiteID_corrected.txt",
  sep = "\t")

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)
phospho.df$SUB_SITE <- rownames(phospho.df)

setwd(
  "~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/"
)
# download sf files from Synapse based on 01_mpnst_get_omics.R from Sara JC Gosline
synapser::synLogin()

##now get the manifest from synapse
manifest<-synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
  dplyr::rename(common_name='Sample')

##for now we only have tumor and PDX data
##they each get their own sample identifier
pdx_data<-manifest|>dplyr::select(common_name,Sex,RNASeq='PDX_RNASeq',Mutations='PDX_Somatic_Mutations',CopyNumber='PDX_CNV')
pdx_data$sampleType <- "PDX"

tumor_data<- manifest|>dplyr::select(common_name,Sex,RNASeq='Tumor_RNASeq',Mutations='Tumor_Somatic_Mutations',CopyNumber='Tumor_CNV')
tumor_data$sampleType <- "Tumor"

combined<-rbind(pdx_data,tumor_data)|>distinct()

# from coderdata 00-buildGeneFile.R
if(!require('org.Hs.eg.db')){
  BiocManager::install('org.Hs.eg.db')
  library(org.Hs.eg.db)
}

library(dplyr)
##get entrez ids to symbol
entrez<-as.data.frame(org.Hs.egALIAS2EG)

##get entriz ids to ensembl
ens<-as.data.frame(org.Hs.egENSEMBL2EG)

##get transcript ids as well
enst<-as.data.frame(org.Hs.egENSEMBLTRANS)

joined.df<-entrez%>%full_join(ens)%>%
  dplyr::rename(entrez_id='gene_id',gene_symbol='alias_symbol',other_id='ensembl_id')%>%
  mutate(other_id_source='ensembl_gene')

tdf<-entrez|>
  full_join(enst)|>
  dplyr::rename(entrez_id='gene_id',gene_symbol='alias_symbol',other_id='trans_id')|>
  dplyr::mutate(other_id_source='ensembl_transcript')

joined.df<-rbind(joined.df,tdf)|>
  distinct()

rnaseq<-do.call('rbind',lapply(setdiff(combined$RNASeq,NA),function(x){
  # if(x!=""){
  #print(x)
  sample<-base::subset(combined,RNASeq==x)
  #print(sample)
  res<-data.table::fread(synGet(x)$path)|>
    tidyr::separate(Name,into=c('other_id','vers'),sep='\\.')|>
    left_join(joined.df)|>
    dplyr::select(gene_symbol,TPM)|>
    subset(!is.na(gene_symbol)|>
    subset(TPM!=0))
  res$sample <- sample$common_name
  res$sex <- sample$Sex
  res$sampleType <- sample$sampleType
  return(distinct(res))
  # }
}))
rnaseq <- na.omit(rnaseq)

valid.chr8 <- c("Not amplified", "Amplified", "Strongly amplified")
Chr8.samples <- manifest[manifest$PDX_Chr8_status_BG %in% valid.chr8 | 
                           manifest$Tumor_Chr8_status_BG %in% valid.chr8,]$common_name # 13
rnaseq <- rnaseq[rnaseq$sample %in% Chr8.samples,]
rna <- reshape2::dcast(rnaseq, sampleType + gene_symbol ~ sample, mean, value.var = "TPM")
colnames(rna)[2] <- "Gene"

pdxRNA <- rna[rna$sampleType=="PDX", colnames(rna)[2:ncol(rna)]]
tumorRNA <- rna[rna$sampleType=="Tumor", colnames(rna)[2:ncol(rna)]]

#### 2. Venn diagrams ####
library(dplyr)
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/"
setwd(base.path)
phospho.venn <- phospho.df  %>% tidyr::extract(SUB_SITE, "Gene",
                                               remove = FALSE)
venn.list <- list('Global old' = unique(global.w.mouse.df$Gene),
                  'Global' = unique(global.df$Gene),
                  'RNA-Seq' = unique(pdxRNA$Gene),
                  'Phospho' = unique(phospho.venn$Gene))
ggvenn::ggvenn(venn.list, set_name_size = 3)
ggplot2::ggsave("Chr8_venn_diagram_20240415.pdf")
ggvenn::ggvenn(venn.list, set_name_size = 3, show_percentage = FALSE)
ggplot2::ggsave("Chr8_venn_diagram_wo_percent_20240415.pdf")

venn.list <- list('Proteomics' = unique(global.w.mouse.df$Gene),
                  'RNA-Seq' = unique(pdxRNA$Gene),
                  'Phospho' = unique(phospho.venn$Gene))
ggvenn::ggvenn(venn.list, set_name_size = 4, show_percentage = FALSE)
ggplot2::ggsave("Chr8_venn_diagram_global_old_20240415.pdf")

venn.list <- list('Proteomics' = unique(global.df$Gene),
                  'RNA-Seq' = unique(pdxRNA$Gene),
                  'Phospho' = unique(phospho.venn$Gene))
ggvenn::ggvenn(venn.list, set_name_size = 4, show_percentage = FALSE)
ggplot2::ggsave("Chr8_venn_diagram_global_new_20240415.pdf")

#### 3. run panSEA ####
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/"
setwd(base.path)
annotations <- c("Sex")

# make sure meta data has right contrasts
global.meta <- manifest[!is.na(manifest$PDX_Proteomics),
                        c("common_name", annotations, "PDX_Chr8_status_BG", "Tumor_Chr8_status_BG")]
global.meta.df2 <- merge(global.meta, meta.df, by.x="common_name", by.y="SampleID_revised")
global.meta.df2$id <- paste0("X", global.meta.df2$id)

global.meta.df2$Chr8 <- "Amplified"
global.meta.df2[global.meta.df2$PDX_Chr8_status_BG == "Not amplified", ]$Chr8 <- "Not_amplified"
rownames(global.meta.df2) <- global.meta.df2$id
colnames(global.meta.df2)[1] <- "Sample"

# rename proteomics columns with common sample names
prot.samples <- colnames(global.df)[1:12]
global.meta.df2 <- global.meta.df2[prot.samples,]
prot.names <- global.meta.df2$SampleName
colnames(global.df)[1:12] <- prot.names
global.df <- global.df[,c("Gene", prot.names)]
colnames(phospho.df)[1:12] <- prot.names
phospho.df <- phospho.df[,c("SUB_SITE", prot.names)]

rna.meta.df <- dplyr::distinct(global.meta.df2[,c("Sample", "Chr8", annotations)])
rownames(rna.meta.df) <- rna.meta.df$Sample

rownames(global.meta.df2) <- global.meta.df2$SampleName
global.meta.df2$Chr8_v2 <- global.meta.df2$Chr8
global.meta.df2[global.meta.df2$Sample == "JH-2-002",]$Chr8_v2 <- "Not_amplified"
rna.meta.df$Chr8_v2 <- rna.meta.df$Chr8
rna.meta.df[rna.meta.df$Sample == "JH-2-002",]$Chr8_v2 <- "Not_amplified"

synapser::synLogin()

# look at chr8q amplification
setwd(base.path)
dir.create("Chr8_quant")
setwd("Chr8_quant")
cnv <- read.csv("mpnst_copy_number.csv")
cnv.samples <- read.csv("mpnst_samples.csv")
cnv.samples <- dplyr::distinct(cnv.samples[,c("common_name", "improve_sample_id")])
prot.cnv.samples <- cnv.samples[cnv.samples$common_name %in% global.meta.df2$Sample,]
prot.cnv <- merge(prot.cnv.samples, cnv, by="improve_sample_id")
cnv <- merge(cnv.samples, cnv, by="improve_sample_id")
genes <- read.csv("genes.csv")
genes <- dplyr::distinct(genes[,c("entrez_id", "gene_symbol")])
prot.cnv <- merge(genes, prot.cnv, by="entrez_id")
cnv <- merge(genes, cnv, by="entrez_id")
pdxCNV <- prot.cnv[prot.cnv$study == "MPNST PDX MT",]
pdxCNV <- cnv[cnv$study == "MPNST PDX MT",]
omics2 <- list("Copy Number" = pdxCNV,
               "RNA-Seq" = pdxRNA,
               "Proteomics" = global.df)
dir.create("Single_sample_positional_GSEA")
setwd("Single_sample_positional_GSEA")
library(plyr); library(dplyr)
library(ggplot2)

# get annotations as df
dir.create("Single_sample_positional_medians")
setwd("Single_sample_positional_medians")
dir.create("Single_sample_positional_medians_allSamples")
setwd("Single_sample_positional_medians_allSamples")
omics2 <- list("PDX Copy Number" = pdxCNV,
               "PDX RNA-Seq" = pdxRNA,
               "Tumor RNA-Seq" = tumorRNA)
library(plyr); library(dplyr)
library(ggplot2)
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")
reduced.msigdb <- dplyr::distinct(msigdb.info[,c("gs_name", "gene_symbol")])
for (i in 1:length(omics2)) {
  setwd(file.path(base.path, "Chr8_quant", "Single_sample_positional_medians_allSamples"))
  dir.create(names(omics2)[i])
  setwd(names(omics2)[i])
  pos.GSEA <- list()
  if (names(omics2)[i] == "PDX Copy Number") {
    temp.samples <- unique(omics2[[i]]$common_name)

    for (j in temp.samples) { # assumes first column is feature names
      # prep input data
      temp.input <- omics2[[i]][omics2[[i]]$common_name == j,c("gene_symbol","copy_number")]
      temp.input <- merge(reduced.msigdb, temp.input, by="gene_symbol")
      temp.med <- plyr::ddply(temp.input, .(gs_name), summarize,
                              median_copy_number = median(copy_number, na.rm = TRUE),
                              mean_copy_number = mean(copy_number, na.rm = TRUE),
                              sd_copy_number = sd(copy_number, na.rm = TRUE))
      pos.GSEA[[as.character(j)]] <- temp.med
    }
  } else {
    for (j in 2:ncol(omics2[[i]])) { # assumes first column is feature names
      temp.input <- omics2[[i]][,c(1, j)]
      colnames(temp.input) <- c("gene_symbol", "copy_number")
      temp.input <- merge(reduced.msigdb, temp.input, by="gene_symbol")
      temp.med <- plyr::ddply(temp.input, .(gs_name), summarize,
                              median_copy_number = median(copy_number, na.rm = TRUE),
                              mean_copy_number = mean(copy_number, na.rm = TRUE),
                              sd_copy_number = sd(copy_number, na.rm = TRUE))
      pos.GSEA[[as.character(j)]] <- temp.med
    }
  }
  all.pos.GSEA <- data.table::rbindlist(pos.GSEA, use.names = TRUE, idcol = "Sample")
  write.csv(all.pos.GSEA, paste0(names(omics2)[i], "_positional_median_mean.csv"), row.names = FALSE)
  #all.pos.GSEA <- read.csv(paste0(names(omics2)[i], "_positional_median_mean.csv"))
  
  # bar graph of median copy number for each position of chr8q 
  chr8q.df <- all.pos.GSEA[grepl("chr8q",all.pos.GSEA$gs_name),]
  chr8q.df$"Chr8q Position" <- as.factor(sub("chr8q", "", chr8q.df$gs_name))
  write.csv(chr8q.df, paste0(names(omics2)[i], "_Chr8q_expression.csv"), row.names = FALSE)
  chr8q.df$`Median Copy Number` <- chr8q.df$median_copy_number
  ggplot2::ggplot(chr8q.df, aes(fill=Sample, x=`Chr8q Position`, y=`Median Copy Number`)) + 
    ggplot2::geom_bar(stat="identity", position="dodge") + ggplot2::theme_minimal() +
    geom_errorbar(aes(ymin=`Median Copy Number` - sd_copy_number, 
                      ymax = `Median Copy Number` + sd_copy_number), width=0.2,
                  position=position_dodge(0.9))
  ggplot2::ggsave(paste0(names(omics2)[i], "_Chr8q_positional_median.pdf"))
  
  chr8q.df$`Mean Copy Number` <- chr8q.df$mean_copy_number
  ggplot2::ggplot(chr8q.df, aes(fill=Sample, x=`Chr8q Position`, y=`Mean Copy Number`)) + 
    ggplot2::geom_bar(stat="identity", position="dodge") + ggplot2::theme_minimal() +
    geom_errorbar(aes(ymin=`Mean Copy Number` - sd_copy_number, 
                      ymax = `Mean Copy Number` + sd_copy_number), width=0.2,
                  position=position_dodge(0.9))
  ggplot2::ggsave(paste0(names(omics2)[i], "_Chr8q_positional_mean.pdf"))
  
  # median chr8q
  med.chr8q <- plyr::ddply(chr8q.df, .(Sample), summarize,
                           `Chr8q Median` = median(`Median Copy Number`, na.rm = TRUE),
                           `Chr8q Mean` = median(`Mean Copy Number`, na.rm = TRUE),
                           sd_med_copy_number = sd(`Median Copy Number`, na.rm = TRUE),
                           sd_mean_copy_number = sd(`Mean Copy Number`, na.rm = TRUE))
  pdf(paste0(names(omics2)[i], "_Chr8q_median_histogram.pdf"))
  hist(med.chr8q$`Chr8q Median`, xlab = "Chr8q Median", main = NULL)
  dev.off()
  med.chr8q <- med.chr8q[order(med.chr8q$`Chr8q Median`),]
  med.chr8q$'Sample ID' <- med.chr8q$Sample
  for (j in 1:nrow(med.chr8q)) {
    med.chr8q$Sample[j] <- stringr::str_split(med.chr8q$`Sample ID`[j], "_")[[1]][1]
  }
  ggplot2::ggplot(med.chr8q, aes(x=`Sample ID`, y=`Chr8q Median`, fill = Sample)) + 
    ggplot2::geom_bar(stat="identity", position="dodge") + 
    scale_x_discrete(limits = med.chr8q$`Sample ID`) + ggplot2::theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
    geom_errorbar(aes(ymin=`Chr8q Median` - sd_med_copy_number, 
                      ymax = `Chr8q Median` + sd_med_copy_number), width=0.2,
                  position=position_dodge(0.9))
  ggplot2::ggsave(paste0(names(omics2)[i], "_Chr8q_median.pdf"))
  write.csv(med.chr8q, paste0(names(omics2)[i], "_Chr8q_median.csv"), row.names = FALSE)
  
  pdf(paste0(names(omics2)[i], "_Chr8q_mean_histogram.pdf"))
  hist(med.chr8q$`Chr8q Mean`, xlab = "Chr8q Mean", main = NULL)
  dev.off()
  med.chr8q <- med.chr8q[order(med.chr8q$`Chr8q Mean`),]
  med.chr8q$'Sample ID' <- med.chr8q$Sample
  for (j in 1:nrow(med.chr8q)) {
    med.chr8q$Sample[j] <- stringr::str_split(med.chr8q$`Sample ID`[j], "_")[[1]][1]
  }
  ggplot2::ggplot(med.chr8q, aes(x=`Sample ID`, y=`Chr8q Mean`, fill = Sample)) + 
    ggplot2::geom_bar(stat="identity", position="dodge") + 
    scale_x_discrete(limits = med.chr8q$`Sample ID`) + ggplot2::theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
    geom_errorbar(aes(ymin=`Chr8q Mean` - sd_mean_copy_number, 
                      ymax = `Chr8q Mean` + sd_mean_copy_number), width=0.2,
                  position=position_dodge(0.9))
  ggplot2::ggsave(paste0(names(omics2)[i], "_Chr8q_mean.pdf"))
  write.csv(med.chr8q, paste0(names(omics2)[i], "_Chr8q_mean.csv"), row.names = FALSE)
  
  # Myc expression
  if (nrow(omics2[[i]][omics2[[i]]$Gene == "MYC",]) == 1) {
    pdf(paste0(names(omics2)[i], "_MYC_histogram.pdf"))
    hist(omics2[[i]][omics2[[i]]$Gene == "MYC",], xlab = "MYC Expression", main = NULL)
    dev.off()
  }
  #med.chr8q <- read.csv(paste0(names(omics2)[i], "_Chr8q_mean.csv"))
  
  # # reshape long if necessary
  # if (names(omics2)[i] == "PDX Copy Number") {
  #   temp.input <- omics2[[i]][,c("common_name", "gene_symbol", "copy_number")]
  # } else {
  #   temp.input <- reshape2::dcast()
  # }
  # 
  # med.chr8q$'Sample ID' <- med.chr8q$Sample
  # for (j in 1:nrow(med.chr8q)) {
  #   med.chr8q$Sample[j] <- stringr::str_split(med.chr8q$`Sample ID`[j], "_")[[1]][1]
  # }
  # ggplot2::ggplot(med.chr8q, aes(x=`Sample ID`, y=`Chr8q Mean`, fill = Sample)) + 
  #   ggplot2::geom_bar(stat="identity", position="dodge") + 
  #   scale_x_discrete(limits = med.chr8q$`Sample ID`) + ggplot2::theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  #   geom_errorbar(aes(ymin=`Chr8q Mean` - sd_mean_copy_number, 
  #                     ymax = `Chr8q Mean` + sd_mean_copy_number), width=0.2,
  #                 position=position_dodge(0.9))
  # ggplot2::ggsave(paste0(names(omics2)[i], "_Chr8q_mean.pdf"))
  # write.csv(med.chr8q, paste0(names(omics2)[i], "_Chr8q_mean.csv"), row.names = FALSE)
}

setwd("Chr8_quant")
cnv.med.chr8q <- read.csv("Single_sample_positional_medians/Copy Number/Copy Number_Chr8q_median.csv")
global.meta.df2$Chr8q_median <- NA
for (i in 1:nrow(global.meta.df2)) {
  temp.sample <- global.meta.df2$Sample[i]
  global.meta.df2$Chr8q_median[i] <- cnv.med.chr8q[cnv.med.chr8q$Sample == temp.sample,]$Chr8q.Median
}
rna.meta.df$Chr8q_median <- NA
for (i in 1:nrow(rna.meta.df)) {
  temp.sample <- rna.meta.df$Sample[i]
  rna.meta.df$Chr8q_median[i] <- cnv.med.chr8q[cnv.med.chr8q$Sample == temp.sample,]$Chr8q.Median
}

setwd(base.path)
setwd("Chr8_quant")
cnv.med.chr8q <- read.csv("Single_sample_positional_medians_allSamples/PDX Copy Number/PDX Copy Number_Chr8q_median.csv")
global.meta.df2$Chr8q_median <- NA
for (i in 1:nrow(global.meta.df2)) {
  temp.sample <- global.meta.df2$Sample[i]
  global.meta.df2$Chr8q_median[i] <- cnv.med.chr8q[cnv.med.chr8q$Sample == temp.sample,]$Chr8q.Median
}
rna.meta.df <- manifest[!is.na(manifest$PDX_RNASeq),c("common_name", "Sex")]
colnames(rna.meta.df)[1] <- "Sample"
rownames(rna.meta.df) <- rna.meta.df$Sample
rna.meta.df$Chr8q_median <- NA
for (i in 1:nrow(rna.meta.df)) {
  temp.sample <- rna.meta.df$Sample[i]
  if (temp.sample %in% cnv.med.chr8q$Sample) {
    rna.meta.df$Chr8q_median[i] <- cnv.med.chr8q[cnv.med.chr8q$Sample == temp.sample,]$Chr8q.Median 
  }
}
pdx.rna.meta.df <- rna.meta.df[colnames(pdxRNA)[2:ncol(pdxRNA)],]
tumor.rna.meta.df <- rna.meta.df[colnames(tumorRNA)[2:ncol(tumorRNA)],]

# get gene sets



omics <- list("Proteomics" = list("Global" = global.df, "Phospho" = phospho.df),
              "RNA-Seq" = list("PDX" = pdxRNA,"Tumor" = tumorRNA))
prot.feat <- c("Gene", "SUB_SITE")
rna.feat <- c("Gene", "Gene")

meta.list <- list("Proteomics" = global.meta.df2,
                  "RNA-Seq" = pdx.rna.meta.df)
expr.list <- list("Proteomics" = list("CCLE proteomics"),
                  "RNA-Seq" = list("adherent CCLE", "adherent CCLE"))
feature.list <- list("Proteomics" = prot.feat,
                     "RNA-Seq" = rna.feat)
dir.create("all_RNA_samples")
setwd("all_RNA_samples")
my.syn <- "syn60219614"
panSEA_corr(omics, meta.list, feature.list, rank.col = "Gain of C8",
            other.annotations = c("Ave CEP8", "Ave CMYC"), expr.list,
            temp.path = file.path(base.path, "Chr8_quant"), syn.id = my.syn)

#### 4. Revise plot formatting & highlight proteins of interest ####
# protein correlation volcano plot
# shape by: on Chr8q or not
# color by: passes FDR threshold or not

# get protein correlations
setwd(base.path)
plot.data <- read.csv("Chr8_quant/Proteomics/Global/Differential_expression/Differential_expression_results.csv")

# annotate whether on Chr8q or not
pos.info <- msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")
chr8q.info <- pos.info[grepl("chr8q", pos.info$gs_name),]
plot.data$Chr8q <- FALSE
plot.data[toupper(plot.data$Gene) %in% toupper(chr8q.info$gene_symbol),]$Chr8q <- TRUE

if (nrow(plot.data[plot.data$Pearson.p == 0, ]) > 0) {
  plot.data[plot.data$Pearson.p == 0, ]$Pearson.p <- 0.00099
}
FDR <- 0.05
# define x, y limits
limit.x <- ceiling(max(abs(as.numeric(plot.data$Pearson.est)), na.rm = TRUE))
limit.y <- ceiling(max(-as.numeric(log(plot.data$Pearson.p, 10)), na.rm = TRUE))

if (nrow(plot.data[plot.data$Pearson.q < FDR, ]) > 0) {
  plot.data$Significance <- paste0("FDR > ", FDR)
  plot.data[plot.data$Pearson.q < FDR, ]$Significance <- paste0("FDR < ", FDR)
  plot.data$Significance <- factor(plot.data$Significance,
                                   levels = c(paste0("FDR < ", FDR),
                                              paste0("FDR > ", FDR)))
  volc <- ggplot2::ggplot(data = plot.data, aes(
    x = Pearson.est, y = -log(Pearson.p, 10),
    color = Significance, shape = Chr8q
  )) +
    ggplot2::geom_point(size = 4) +
    ggrepel::geom_text_repel(
      data = subset(plot.data, Significance == paste0("FDR < ", FDR)),
      mapping = aes(label = Gene, size = I(6)), nudge_y = 0.25
    ) +
    ggplot2::scale_color_manual(
      values = c("red", "azure4"), name = "Significance",
      breaks = c(paste0("FDR < ", FDR), paste0("FDR > ", FDR))
    ) +
    ggplot2::xlim(-limit.x, limit.x) +
    ggplot2::ylim(0, limit.y) +
    ggplot2::xlab("Pearson Correlation Estimate") +
    ggplot2::ylab("-Log(p-value)") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey",
                        linewidth = 0.5) +
    ggplot2::theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_line(colour = "black", linewidth = 0.65),
      legend.text = element_text(size = 20),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 26, face = "bold"),
      panel.background = element_rect(
        fill = "white", colour = "white", linewidth = 0.5,
        linetype = "solid", color = "black"
      ), text = element_text(size = 20),
      legend.position = "bottom", legend.key = element_blank()
    )
} else {
  volc <- ggplot2::ggplot(data = plot.data,
                          aes(x = Pearson.est, y = -log(Pearson.p, 10),
                              shape = Chr8q)) +
    ggplot2::geom_point(size = 4, color = "azure4") +
    ggplot2::xlim(-limit.x, limit.x) +
    ggplot2::ylim(0, limit.y) +
    ggplot2::xlab("Pearson Correlation Estimate") +
    ggplot2::ylab("-Log(p-value)") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey",
                        linewidth = 0.5) +
    ggplot2::theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_line(colour = "black", linewidth = 0.65),
      axis.text = element_text(size = 20), axis.title =
        element_text(size = 26, face = "bold"),
      panel.background = element_rect(
        fill = "white", colour = "white", linewidth = 0.5,
        linetype = "solid", color = "black"
      ), text = element_text(size = 20)
    )
}
ggsave("Chr8_quant/Proteomics/Global/Differential_expression/Differential_expression_volcano_width11.pdf", volc, height=7, width=11)

# just annotate MYCBP, BCL2L2-PABPN1
goi <- c("MYCBP", "BCL2L2-PABPN1")
plot.data$Label <- FALSE
plot.data[plot.data$Gene %in% goi,]$Label <- TRUE
if (nrow(plot.data[plot.data$Pearson.q < FDR, ]) > 0) {
  plot.data$Significance <- paste0("FDR > ", FDR)
  plot.data[plot.data$Pearson.q < FDR, ]$Significance <- paste0("FDR < ", FDR)
  plot.data$Significance <- factor(plot.data$Significance,
                                   levels = c(paste0("FDR < ", FDR),
                                              paste0("FDR > ", FDR)))
  volc <- ggplot2::ggplot(data = plot.data, aes(
    x = Pearson.est, y = -log(Pearson.p, 10),
    color = Significance, shape = Chr8q
  )) +
    ggplot2::geom_point(size = 4) +
    ggrepel::geom_text_repel(
      data = subset(plot.data, Label),
      mapping = aes(label = Gene, size = I(6)), nudge_y = 0.25, nudge_x = -0.5
    ) +
    ggplot2::scale_color_manual(
      values = c("red", "azure4"), name = "Significance",
      breaks = c(paste0("FDR < ", FDR), paste0("FDR > ", FDR))
    ) +
    ggplot2::xlim(-limit.x, limit.x) +
    ggplot2::ylim(0, limit.y) +
    ggplot2::xlab("Pearson Correlation Estimate") +
    ggplot2::ylab("-Log(p-value)") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey",
                        linewidth = 0.5) +
    ggplot2::theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_line(colour = "black", linewidth = 0.65),
      legend.text = element_text(size = 20),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 26, face = "bold"),
      panel.background = element_rect(
        fill = "white", colour = "white", linewidth = 0.5,
        linetype = "solid", color = "black"
      ), text = element_text(size = 20),
      legend.position = "bottom", legend.key = element_blank()
    )
} else {
  volc <- ggplot2::ggplot(data = plot.data,
                          aes(x = Pearson.est, y = -log(Pearson.p, 10),
                              shape = Chr8q)) +
    ggplot2::geom_point(size = 4, color = "azure4") +
    ggplot2::xlim(-limit.x, limit.x) +
    ggplot2::ylim(0, limit.y) +
    ggplot2::xlab("Pearson Correlation Estimate") +
    ggplot2::ylab("-Log(p-value)") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey",
                        linewidth = 0.5) +
    ggplot2::theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_line(colour = "black", linewidth = 0.65),
      axis.text = element_text(size = 20), axis.title =
        element_text(size = 26, face = "bold"),
      panel.background = element_rect(
        fill = "white", colour = "white", linewidth = 0.5,
        linetype = "solid", color = "black"
      ), text = element_text(size = 20)
    )
}
ggsave("Chr8_quant/Proteomics/Global/Differential_expression/Differential_expression_volcano_only_MYCBP_BCL2L2_labeled_v8.pdf", volc, height=7, width=11)

# redo volcano plots with bigger text
sub.dir <- c("", "GSEA_KEGG/Global", "GSEA_phospho_ksdb")
ea.results <- c("DMEA", "GSEA", "Phospho_enrichment")
temp.omics <- c("Global", "Global", "Phospho")
ea.sets <- c("Drug_set", "Feature_set", "Feature_set")
FDR <- 0.05
library(ggplot2)
for (i in 2:length(ea.results)) {
  setwd(file.path(base.path, "Chr8_quant", "Proteomics"))
  setwd(temp.omics[i])
  setwd(ea.results[i])
  if (sub.dir[i] != "") {setwd(sub.dir[i])}
  # get results
  plot.data <- read.csv(paste0(ea.results[i], "_results.csv"))

  if (nrow(plot.data[plot.data$p_value == 0, ]) > 0) {
    plot.data[plot.data$p_value == 0, ]$p_value <- 0.00099
  }
  colnames(plot.data)[2] <- "Drug_set"
  
  # define x, y limits
  limit.x <- ceiling(max(abs(as.numeric(plot.data$NES)), na.rm = TRUE))
  limit.y <- ceiling(max(-as.numeric(log(plot.data$p_value, 10)), na.rm = TRUE))
  
  # categorize data by significance level if there are significant hits
  if (nrow(plot.data[plot.data$FDR_q_value < FDR, ]) > 0) {
    plot.data$Significance <- paste0("FDR > ", FDR)
    plot.data[plot.data$FDR_q_value < FDR, ]$Significance <- paste0("FDR < ", FDR)
    plot.data$Significance <- factor(plot.data$Significance,
                                     levels = c(paste0("FDR < ", FDR),
                                                paste0("FDR > ", FDR)))
    volc <- ggplot2::ggplot(data = plot.data, aes(
      x = NES, y = -log(p_value, 10),
      color = Significance
    )) +
      ggplot2::geom_point(size = 4) +
      ggrepel::geom_text_repel(
        data = subset(plot.data, Significance == paste0("FDR < ", FDR)),
        mapping = aes(label = Drug_set, size = I(4)), nudge_y = 0.25
      ) +
      ggplot2::scale_color_manual(
        values = c("red", "azure4"), name = "Significance",
        breaks = c(paste0("FDR < ", FDR), paste0("FDR > ", FDR))
      ) +
      ggplot2::xlim(-limit.x, limit.x) +
      ggplot2::ylim(0, limit.y) +
      ggplot2::xlab("Normalized Enrichment Score") +
      ggplot2::ylab("-Log(p-value)") +
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey",
                          linewidth = 0.5) +
      ggplot2::theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_line(colour = "black", linewidth = 0.65),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 26, face = "bold"),
        panel.background = element_rect(
          fill = "white", colour = "white", linewidth = 0.5,
          linetype = "solid", color = "black"
        ), text = element_text(size = 20),
        legend.position = "bottom", legend.key = element_blank()
      )
  } else {
    volc <- ggplot2::ggplot(data = plot.data,
                            aes(x = NES, y = -log(p_value, 10))) +
      ggplot2::geom_point(size = 4, color = "azure4") +
      ggplot2::xlim(-limit.x, limit.x) +
      ggplot2::ylim(0, limit.y) +
      ggplot2::xlab("Normalized Enrichment Score") +
      ggplot2::ylab("-Log(p-value)") +
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey",
                          linewidth = 0.5) +
      ggplot2::theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_line(colour = "black", linewidth = 0.65),
        axis.text = element_text(size = 20), axis.title =
          element_text(size = 26, face = "bold"),
        panel.background = element_rect(
          fill = "white", colour = "white", linewidth = 0.5,
          linetype = "solid", color = "black"
        ), text = element_text(size = 20)
      )
  }
  ggsave(paste0(ea.results[i], "_volcano_plot_larger_font_v5.pdf"), volc, height=7, width=7) # v4: size 6 labels, v5: size 4
}
