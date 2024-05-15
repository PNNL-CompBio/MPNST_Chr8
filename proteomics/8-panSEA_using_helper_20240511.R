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

#### 2. PCA plots ####
library(MSnSet.utils)
library(ggplot2)

### Transcriptomics PCA
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/transcriptomics/")

# create metadata for RNA-Seq
phenos <- c("common_name", "Sex", 
            "PDX_Chr8_status_BG", "Tumor_Chr8_status_BG")
pca.meta.df <- manifest[manifest$common_name %in% Chr8.samples, phenos]
row.names(pca.meta.df) <- pca.meta.df$common_name
m_RNA <- MSnSet(exprs = pdxRNA[,Chr8.samples] %>% as.matrix(),
                pData = pca.meta.df)
m_RNA_tumor <- MSnSet(exprs = tumorRNA[,Chr8.samples] %>% as.matrix(),
                pData = pca.meta.df)
# WU-561 tumor is all NaN so removing
Chr8.tumor.samples <- Chr8.samples[Chr8.samples != "WU-561"]
tumorRNA <- tumorRNA[,c("Gene", Chr8.tumor.samples)]
pca.tumor.meta.df <- manifest[manifest$common_name %in% Chr8.tumor.samples, phenos]
row.names(pca.tumor.meta.df) <- pca.tumor.meta.df$common_name
m_RNA_tumor <- MSnSet(exprs = tumorRNA[,Chr8.tumor.samples] %>% as.matrix(),
                      pData = pca.tumor.meta.df)
for (i in 1:length(phenos)) {
  MSnSet.utils::plot_pca(m_RNA, phenotype = phenos[i]) + ggtitle("PDX RNA-Seq PCA") # 250 complete rows for PCA out of 28,803
  ggsave(paste0("RNAseq_PDX_PCA_", phenos[i], "_", Sys.Date(), ".pdf"))
  MSnSet.utils::plot_pca(m_RNA_tumor, phenotype = phenos[i]) + ggtitle("Tumor RNA-Seq PCA") # 321 complete rows for PCA out of 28,672
  ggsave(paste0("RNAseq_tumor_PCA_", phenos[i], "_", Sys.Date(), ".pdf"))
}
# outliers:
# tumor: JH-2-023 (not amp but near strongly amplified group on PCA); JH-2-002 is near strongly amplified on PCA despite just being amp
# PDX: WU-356 and JH-2-031 strong amp but closer to not amp on PCA

# try filtering genes for being measured in at least 50% of samples
pdxRNA50 <- as.data.frame(pdxRNA[which(rowMeans(!is.na(pdxRNA)) >= 0.5),]) # 9629 out of 28803 pass
#pdxRNA502 <- as.data.frame(pdxRNA50[,which(colMeans(!is.na(pdxRNA50)) >= 0.5)]) # all 13 samples pass
tumorRNA50 <- as.data.frame(tumorRNA[which(rowMeans(!is.na(tumorRNA)) >= 0.5),]) # 8474 out of 28672 pass
#tumorRNA502 <- as.data.frame(tumorRNA50[,which(colMeans(!is.na(tumorRNA50)) >= 0.5)]) # all 13 samples pass

m_RNA <- MSnSet(exprs = pdxRNA50[,Chr8.samples] %>% as.matrix(),
                pData = pca.meta.df)
m_RNA_tumor <- MSnSet(exprs = tumorRNA50[,Chr8.tumor.samples] %>% as.matrix(),
                      pData = pca.tumor.meta.df)
for (i in 1:length(phenos)) {
  MSnSet.utils::plot_pca(m_RNA, phenotype = phenos[i]) + ggtitle("PDX RNA-Seq PCA") # still 250 complete rows for PCA
  ggsave(paste0("RNAseq_PDX_PCA_afterMissingnessFilter_", phenos[i], "_", Sys.Date(), ".pdf"))
  MSnSet.utils::plot_pca(m_RNA_tumor, phenotype = phenos[i]) + ggtitle("Tumor RNA-Seq PCA") # still 321 complete rows for PCA
  ggsave(paste0("RNAseq_tumor_PCA_afterMissingnessFilter_", phenos[i], "_", Sys.Date(), ".pdf"))
}

# also try only looking at samples which match Chr8 status across tumor & PDX
Chr8.matching <- manifest[manifest$PDX_Chr8_status_BG == manifest$Tumor_Chr8_status_BG,]$common_name
Chr8.matching <- Chr8.matching[Chr8.matching%in% Chr8.samples]
pdxRNAmatch <- pdxRNA[,c("Gene",Chr8.matching)]
tumorRNAmatch <- tumorRNA[,c("Gene",Chr8.matching[Chr8.matching != "WU-561"])]
m_RNA <- MSnSet(exprs = pdxRNAmatch[,2:ncol(pdxRNAmatch)] %>% as.matrix(),
                pData = pca.meta.df[colnames(pdxRNAmatch)[2:ncol(pdxRNAmatch)],])
m_RNA_tumor <- MSnSet(exprs = tumorRNAmatch[,2:ncol(tumorRNAmatch)] %>% as.matrix(),
                pData = pca.tumor.meta.df[colnames(tumorRNAmatch)[2:ncol(tumorRNAmatch)],])
for (i in 1:length(phenos)) {
  MSnSet.utils::plot_pca(m_RNA, phenotype = phenos[i]) + ggtitle("PDX RNA-Seq PCA") # 392 complete rows for PCA
  ggsave(paste0("RNAseq_PDX_PCA_matchingTumor_", phenos[i], "_", Sys.Date(), ".pdf"))
  MSnSet.utils::plot_pca(m_RNA_tumor, phenotype = phenos[i]) + ggtitle("Tumor RNA-Seq PCA") # 596 complete rows for PCA
  ggsave(paste0("RNAseq_tumor_PCA_matchingPDX_", phenos[i], "_", Sys.Date(), ".pdf"))
}

pdxRNAmatch50 <- as.data.frame(pdxRNAmatch[which(rowMeans(!is.na(pdxRNAmatch)) >= 0.5),]) # 11840 out of 28803 pass
#pdxRNAmatch502 <- as.data.frame(pdxRNAmatch50[,which(colMeans(!is.na(pdxRNAmatch50)) >= 0.5)]) # all 9 samples pass
tumorRNAmatch50 <- as.data.frame(tumorRNAmatch[which(rowMeans(!is.na(tumorRNAmatch)) >= 0.5),]) # 10189 out of 28672 pass
#tumorRNAmatch502 <- as.data.frame(tumorRNAmatch50[,which(colMeans(!is.na(tumorRNAmatch50)) >= 0.5)]) # all 8 samples pass (8 instead of 9 because no WU-561)

m_RNA <- MSnSet(exprs = pdxRNAmatch50[,2:ncol(pdxRNAmatch50)] %>% as.matrix(),
                pData = pca.meta.df[colnames(pdxRNAmatch50)[2:ncol(pdxRNAmatch50)],])
m_RNA_tumor <- MSnSet(exprs = tumorRNAmatch50[,2:ncol(tumorRNAmatch50)] %>% as.matrix(),
                      pData = pca.tumor.meta.df[colnames(tumorRNAmatch50)[2:ncol(tumorRNAmatch50)],])
for (i in 1:length(phenos)) {
  MSnSet.utils::plot_pca(m_RNA, phenotype = phenos[i]) + ggtitle("PDX RNA-Seq PCA") # still 392 complete rows for PCA
  ggsave(paste0("RNAseq_PDX_PCA_matchingTumor_afterMissingnessFilter_", phenos[i], "_", Sys.Date(), ".pdf"))
  MSnSet.utils::plot_pca(m_RNA_tumor, phenotype = phenos[i]) + ggtitle("Tumor RNA-Seq PCA") # still 596 complete rows for PCA
  ggsave(paste0("RNAseq_tumor_PCA_matchingPDX_afterMissingnessFilter_", phenos[i], "_", Sys.Date(), ".pdf"))
}

### Proteomics PCA
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
dir.create("PCA")
setwd("PCA")


# create metadata for pca
phenos <- c("common_name", "Sex", 
            "PDX_Chr8_status_BG", "Tumor_Chr8_status_BG")
global.meta.df <- merge(pca.meta.df, meta.df, by.x="common_name", by.y="SampleID")
Chr8.sample.ids <- global.meta.df[global.meta.df$common_name %in% Chr8.samples,]$id
rownames(global.meta.df) <- global.meta.df$id
m_global <- MSnSet(exprs = global.df[,Chr8.sample.ids] %>% as.matrix(),
                pData = global.meta.df)
m_global_old <- MSnSet(exprs = global.w.mouse.df[,Chr8.sample.ids] %>% as.matrix(),
                   pData = global.meta.df)
m_phospho <- MSnSet(exprs = phospho.df[,Chr8.sample.ids] %>% as.matrix(),
                   pData = global.meta.df)
for (i in 1:length(phenos)) {
  MSnSet.utils::plot_pca(m_global, phenotype = phenos[i]) + ggtitle("Global Proteomics PCA") # 8609 complete rows for PCA
  ggsave(paste0("Global_proteomics_PCA_", phenos[i], "_", Sys.Date(), ".pdf"))
  
  MSnSet.utils::plot_pca(m_global_old, phenotype = phenos[i]) + ggtitle("Global Proteomics PCA") # 6196 complete rows for PCA
  ggsave(paste0("Global_proteomics_oldDB_PCA_", phenos[i], "_", Sys.Date(), ".pdf"))
  
  MSnSet.utils::plot_pca(m_phospho, phenotype = phenos[i]) + ggtitle("Phospho-proteomics PCA") # 28592 complete rows for PCA
  ggsave(paste0("Phospho_proteomics_PCA_", phenos[i], "_", Sys.Date(), ".pdf"))
}

Chr8.match.ids <- global.meta.df[global.meta.df$common_name %in% Chr8.matching,]$id
global.df.match <- global.df[,c("Gene",Chr8.match.ids)]
global.w.mouse.df.match <- global.w.mouse.df[,c("Gene",Chr8.match.ids)]
phospho.df.match <- phospho.df[,c("SUB_SITE",Chr8.match.ids)]
m_global <- MSnSet(exprs = global.df.match[,2:ncol(global.df.match)] %>% as.matrix(),
                pData = global.meta.df[colnames(global.df.match)[2:ncol(global.df.match)],])
m_global_old <- MSnSet(exprs = global.w.mouse.df.match[,2:ncol(global.w.mouse.df.match)] %>% as.matrix(),
                   pData = global.meta.df[colnames(global.w.mouse.df.match)[2:ncol(global.w.mouse.df.match)],])
m_phospho <- MSnSet(exprs = phospho.df.match[,2:ncol(phospho.df.match)] %>% as.matrix(),
                   pData = global.meta.df[colnames(phospho.df.match)[2:ncol(phospho.df.match)],])
for (i in 1:length(phenos)) {
  MSnSet.utils::plot_pca(m_global, phenotype = phenos[i]) + ggtitle("Global Proteomics PCA") # 8622 complete rows for PCA
  ggsave(paste0("Global_proteomics_PCA_matchingTumor_", phenos[i], "_", Sys.Date(), ".pdf"))
  
  MSnSet.utils::plot_pca(m_global_old, phenotype = phenos[i]) + ggtitle("Global Proteomics PCA") # 6200 complete rows for PCA
  ggsave(paste0("Global_proteomics_oldDB_PCA_matchingTumor_", phenos[i], "_", Sys.Date(), ".pdf"))
  
  MSnSet.utils::plot_pca(m_phospho, phenotype = phenos[i]) + ggtitle("Phospho-proteomics PCA") # 28652 complete rows for PCA
  ggsave(paste0("Phospho_proteomics_PCA_matchingTumor_", phenos[i], "_", Sys.Date(), ".pdf"))
}

# histograms
omics <- list("PDX RNA-Seq" = pdxRNA,
              "PDX RNA-Seq with Missingness Filter" = pdxRNA50,
              "PDX RNA-Seq Matching with Tumor" = pdxRNAmatch,
              "PDX RNA-Seq Matching with Tumor with Missingness Filter" = pdxRNAmatch50,
              "Tumor RNA-Seq" = tumorRNA,
              "Tumor RNA-Seq with Missingness Filter" = tumorRNA50,
              "Tumor RNA-Seq Matching with PDX" = tumorRNAmatch,
              "Tumor RNA-Seq Matching with PDX with Missingness Filter" = tumorRNAmatch50,
              "PDX Global Proteomics" = global.df,
              "PDX Global Proteomics Matching with Tumor" = global.df.match,
              "PDX Global Proteomics (Old Database)" = global.w.mouse.df,
              "PDX Global Proteomics (Old Database) Matching with Tumor" = global.w.mouse.df.match,
              "PDX Phospho-proteomics" = phospho.df,
              "PDX Phospho-proteomics Matching with Tumor" = phospho.df.match)
for (i in 1:length(omics)) {
  # histogram
  melted.df <- reshape2::melt(omics[[i]])
  xlab <- names(omics)[i]
  pdf(file.path(paste0(names(omics)[i], "_histogram_", Sys.Date(), ".pdf")))
  hist(melted.df$value, xlab = xlab)
  dev.off()
}

# run panSEA
setwd(base.path)

# make sure meta data has right contrasts
global.meta.df$Chr8_amplified <- TRUE
global.meta.df$Chr8_strongly_amplified <- TRUE
global.meta.df$Chr8_strongly_amplified_inclusive <- TRUE
global.meta.df[global.meta.df$PDX_Chr8_status_BG == "Not amplified", ]$Chr8_amplified <- FALSE
global.meta.df[global.meta.df$PDX_Chr8_status_BG == "Not amplified" | 
                 global.meta.df$PDX_Chr8_status_BG == "Amplified", ]$Chr8_strongly_amplified_inclusive <- FALSE
global.meta.df[global.meta.df$PDX_Chr8_status_BG == "Not amplified", ]$Chr8_strongly_amplified <- FALSE
global.meta.df[global.meta.df$PDX_Chr8_status_BG == "Amplified", ]$Chr8_strongly_amplified <- NA


pca.meta.df$Chr8_amplified <- TRUE
pca.meta.df$Chr8_strongly_amplified <- TRUE
pca.meta.df$Chr8_strongly_amplified_inclusive <- TRUE
pca.meta.df[pca.meta.df$PDX_Chr8_status_BG == "Not amplified", ]$Chr8_amplified <- FALSE
pca.meta.df[pca.meta.df$PDX_Chr8_status_BG == "Not amplified" | 
                 pca.meta.df$PDX_Chr8_status_BG == "Amplified", ]$Chr8_strongly_amplified_inclusive <- FALSE
pca.meta.df[pca.meta.df$PDX_Chr8_status_BG == "Not amplified", ]$Chr8_strongly_amplified <- FALSE
pca.meta.df[pca.meta.df$PDX_Chr8_status_BG == "Amplified", ]$Chr8_strongly_amplified <- NA

prot <- list("Global" = global.df,
             "Phospho" = phospho.df)
prot.feat <- c("Gene", "SUB_SITE")
prot <- list("Global" = global.df)
prot.feat <- c("Gene")
global.w.mouse.df$Gene <- toupper(global.w.mouse.df$Gene)
prot.old <- list("Global" = global.w.mouse.df)
prot.old.feat <- "Gene"
rna.list <- list("PDX" = pdxRNA,
                 "Tumor" = tumorRNA)
rna.list.50 <- list("PDX" = pdxRNA50,
                 "Tumor" = tumorRNA50)
rna.feat <- c("Gene", "Gene")
omics <- list("Proteomics" = prot,
              "RNA-Seq" = rna.list,
              "RNA-Seq (After Missingness Filter)" = rna.list.50,
              "Proteomics (Old Database)" = prot.old)
pca.meta.df$id <- pca.meta.df$common_name
meta.list <- list("Proteomics" = global.meta.df,
                  "RNA-Seq" = pca.meta.df,
                  "RNA-Seq (After Missingness Filter)" = pca.meta.df,
                  "Proteomics (Old Database)" = global.meta.df)
expr.list <- list("Proteomics" = list("CCLE proteomics"),
                  "RNA-Seq" = list("adherent CCLE", "adherent CCLE"),
                  "RNA-Seq (After Missingness Filter)" = list("adherent CCLE", "adherent CCLE"),
                  "Proteomics (Old Database)" = list("CCLE proteomics"))
#contrasts <- c("Chr8_amplified", "Chr8_strongly_amplified", "Chr8_strongly_amplified_inclusive")
contrasts <- c("Chr8_strongly_amplified", "Chr8_strongly_amplified_inclusive")
# outliers <- list("Proteomics" = "JH-2-002",
#                  "RNA-Seq" = list("PDX" = c("WU-356", "JH-2-031"),
#                                   "Tumor" = "JH-2-023"),
#                  "RNA-Seq (After Missingness Filter)" = list("PDX" = c("WU-356", "JH-2-031"),
#                                   "Tumor" = "JH-2-023"),
#                  "Proteomics (Old Database)" = "JH-2-002")
#Chr8.match.ids <- global.meta.df[global.meta.df$common_name %in% Chr8.matching,]$id
match <- list("Proteomics" = Chr8.match.ids,
                 "RNA-Seq" = Chr8.matching,
              "RNA-Seq (After Missingness Filter)" = Chr8.matching,
              "Proteomics (Old Database)" = Chr8.match.ids)
feature.list <- list("Proteomics" = prot.feat,
                     "RNA-Seq" = rna.feat,
                     "RNA-Seq (After Missingness Filter)" = rna.feat,
                     "Proteomics (Old Database)" = prot.old.feat)
annotations <- c("Sex", "PDX_Chr8_status_BG", "Tumor_Chr8_status_BG")
global.meta <- manifest[!is.na(manifest$PDX_Proteomics),c("common_name", annotations)]
global.meta.df2 <- merge(global.meta, meta.df, by.x="common_name", by.y="SampleID_revised")
global.meta.df2$id <- paste0("X", global.meta.df2$id)

global.meta.df2$Chr8 <- "Not_amplified"
global.meta.df2$Chr8_inclusive <- "Not_amplified"
global.meta.df2[global.meta.df2$PDX_Chr8_status_BG == "Amplified", ]$Chr8 <- NA
global.meta.df2[global.meta.df2$PDX_Chr8_status_BG == "Strongly amplified", ]$Chr8 <- "Strongly_amplified"
rownames(global.meta.df2) <- global.meta.df2$id

pca.meta.df$Chr8 <- "Not_amplified"
pca.meta.df$Chr8_inclusive <- "Not_amplified"
pca.meta.df[pca.meta.df$PDX_Chr8_status_BG == "Amplified", ]$Chr8 <- NA
pca.meta.df[pca.meta.df$PDX_Chr8_status_BG == "Strongly amplified", ]$Chr8 <- "Strongly_amplified"
rownames(pca.meta.df) <- pca.meta.df$id

contrasts <- c("Chr8", "Chr8_inclusive")

meta.list <- list("Proteomics" = global.meta.df2,
                  "RNA-Seq" = pca.meta.df,
                  "RNA-Seq (After Missingness Filter)" = pca.meta.df,
                  "Proteomics (Old Database)" = global.meta.df2)

# get Chr8 set annotations
gmt1 <- get_chr8_gmt1()
names(gmt1) <- c("KEGG", "Hallmark", "Positional", "Positional_custom")
gmt2 <- get_chr8_gmt2()
gmt2 <- gmt2[[2]] # fails at ksdb
# for some reason failed at ssea too

syn_base <- "syn58649770"
synapser::synLogin()
#debug(run_contrasts2)
for (j in 1:length(contrasts)){
  setwd(base.path)
  dir.create(contrasts[j])
  setwd(contrasts[j])
  contrastFolder <- synapser::synStore(synapser::Folder(contrasts[j], syn_base))
  
  for (i in 1:length(omics)) {
    setwd(file.path(base.path, contrasts[j]))
    dir.create(names(omics)[i])
    setwd(names(omics)[i])
    omicsFolder <- synapser::synStore(synapser::Folder(names(omics)[i], contrastFolder))
    contrast.meta <- meta.list[[i]]
    contrast.meta <- na.omit(contrast.meta[,c(contrasts[j], annotations)])
    
    # run contrast w/ sex factor
    dir.create("no_filter")
    setwd("no_filter")
    noFiltFolder <- synapser::synStore(synapser::Folder("no_filter", omicsFolder))
    temp.path <- file.path(base.path, contrasts[j], names(omics)[i], 
                           "no_filter")
    scale <- TRUE
    cluster <- TRUE
    if (names(omics)[i] == "RNA-Seq" & contrasts[j] == "Chr8") {
      scale <- FALSE
      cluster <- FALSE
    }
    if (!grepl("RNA-Seq", names(omics)[i])) {
      panSEA2(contrasts[j], "Sex", contrast.meta, omics = omics[[i]],
              annotations = annotations, 
              gmt.list1 = gmt1, gmt.list2=gmt2, 
              expr = expr.list[[i]], base.path = base.path, temp.path = temp.path, 
              subfolder = FALSE, synapse_id = noFiltFolder, n.net = 10, scale = scale, cluster = cluster) 
    }
    
    # re-run with only samples for which PDX Chr8 status matches tumor
    if (!grepl("Proteomics", names(omics)[i])) {
      setwd(file.path(base.path, contrasts[j], names(omics)[i]))
      dir.create("PDX_Chr8_status_matches_tumor")
      setwd("PDX_Chr8_status_matches_tumor")
      matchFolder <- synapser::synStore(synapser::Folder("PDX_Chr8_status_matches_tumor", omicsFolder))
      temp.path <- file.path(base.path, contrasts[j], names(omics)[i],
                             "PDX_Chr8_status_matches_tumor")
      omics.match <- omics[[i]]
      for (k in 1:length(omics.match)) {
        match.samples <- c(feature.list[[i]][k], match[[i]])
        omics.match[[k]] <- omics.match[[k]][,match.samples]
      }
      meta.match <- contrast.meta[match[[i]],]
      
      scale <- TRUE
      cluster <- TRUE
      if (names(omics)[i] == "RNA-Seq" & contrasts[j] == "Chr8") {
        scale <- FALSE
        cluster <- FALSE
        }
      panSEA2(contrasts[j], "Sex", meta.match, omics = omics.match,
              annotations = annotations, 
              gmt.list1 = gmt1, gmt.list2=gmt2, 
              expr = expr.list[[i]], base.path = base.path, temp.path = temp.path, 
              subfolder = FALSE, synapse_id = matchFolder, n.net = 10, scale = scale, cluster = cluster) 
      omics.match <- NULL
    } else if (contrasts[j] != "Chr8") {
      setwd(file.path(base.path, contrasts[j], names(omics)[i]))
      dir.create("PDX_Chr8_status_matches_tumor")
      setwd("PDX_Chr8_status_matches_tumor")
      matchFolder <- synapser::synStore(synapser::Folder("PDX_Chr8_status_matches_tumor", omicsFolder))
      temp.path <- file.path(base.path, contrasts[j], names(omics)[i],
                             "PDX_Chr8_status_matches_tumor")
      omics.match <- omics[[i]]
      for (k in 1:length(omics.match)) {
        match.samples <- c(feature.list[[i]][k], match[[i]])
        omics.match[[k]] <- omics.match[[k]][,match.samples]
      }
      meta.match <- na.omit(contrast.meta[match[[i]],])
      
      panSEA2(contrasts[j], "Sex", meta.match, omics = omics.match,
              annotations = annotations, 
              gmt.list1 = gmt1, gmt.list2=gmt2, 
              expr = expr.list[[i]], base.path = base.path, temp.path = temp.path, 
              subfolder = FALSE, synapse_id = matchFolder, n.net = 10) 
      omics.match <- NULL
    }
    # re-run without outliers and w/wo samples which don't match across PDX/tumor
    # setwd(file.path(base.path, contrasts[j], names(omics)[i]))
    # dir.create("no_outliers")
    # setwd("no_outliers")
    # if (is.list(outliers[[i]])) {
    #   for (k in 1:length(outliers[[i]])) {
    #     temp.out <- outliers[[i]][[k]]
    #     outlier.name <- paste0(temp.out, collapse="_")
    #     setwd(file.path(base.path, contrasts[j], names(omics)[i], "no_outliers"))
    #     dir.create(paste0("no_", outlier.name))
    #     setwd(paste0("no_", outlier.name))
    #     noOutFolder <- synapser::synStore(synapser::Folder(paste0("no_", outlier.name), omicsFolder))
    #     temp.path <- file.path(base.path, contrasts[j], names(omics)[i], 
    #                            "no_outliers", paste0("no_", outlier.name))
    #     omics.wo.out <- omics[[i]]
    #     for (m in 1:length(omics.wo.out)) {
    #       omics.wo.out[[m]] <- omics.wo.out[[m]][,!(colnames(omics.wo.out[[m]]) %in% temp.out)]
    #     }
    #     meta.wo.out <- meta.list[[i]][!(meta.list[[i]]$common_name %in% temp.out),]
    #     run_contrasts2(contrasts[j], gmt.list1 = "chr8", gmt.list2 = "chr8", 
    #                    meta.df = meta.wo.out, 
    #                    omics = omics.wo.out, types = names(omics[[i]]), 
    #                    feature.names = feature.list[[i]], 
    #                    EA.types = c("KEGG", "Hallmark", "Positional", "Positional_Custom"),
    #                    base.path = base.path, 
    #                    temp.path = temp.path, subfolder = FALSE, 
    #                    synapse_id = noOutFolder)
    #     
    #     # also run with only samples matching across PDX/tumor
    #     dir.create("PDX_Chr8_status_matches_tumor")
    #     setwd("PDX_Chr8_status_matches_tumor")
    #     noOutMatchFolder <- synapser::synStore(synapser::Folder("PDX_Chr8_status_matches_tumor", noOutFolder))
    #     temp.path <- file.path(base.path, contrasts[j], names(omics)[i], 
    #                            "no_outliers", paste0("no_", outlier.name), 
    #                            "PDX_Chr8_status_matches_tumor")
    #     omics.wo.out.match <- omics.wo.out
    #     omics.wo.out <- NULL
    #     for (m in 1:length(omics.wo.out.match)) {
    #       match.samples <- c(feature.list[[i]][m], match[[i]])
    #       omics.wo.out.match[[m]] <- omics.wo.out.match[[m]][,colnames(omics.wo.out.match[[m]]) %in% match.samples]
    #     }
    #     meta.wo.out.match <- meta.wo.out[meta.wo.out$common_name %in% match[[i]],]
    #     run_contrasts2(contrasts[j], gmt.list1 = "chr8", gmt.list2 = "chr8", 
    #                    meta.df = meta.wo.out.match, 
    #                    omics = omics.wo.out.match, types = names(omics[[i]]), 
    #                    feature.names = feature.list[[i]], 
    #                    EA.types = c("KEGG", "Hallmark", "Positional", "Positional_Custom"),
    #                    base.path = base.path, temp.path = temp.path, 
    #                    subfolder = FALSE, synapse_id = noOutMatchFolder)
    #     omics.wo.out.match <- NULL
    #   }
    # } else {
    #   temp.out <- outliers[[i]]
    #   outlier.name <- paste0(temp.out, collapse="_")
    #   dir.create(paste0("no_", outlier.name))
    #   setwd(paste0("no_", outlier.name))
    #   noOutFolder <- synapser::synStore(synapser::Folder(paste0("no_", outlier.name), omicsFolder))
    #   temp.path <- file.path(base.path, contrasts[j], names(omics)[i], 
    #                          "no_outliers", paste0("no_", outlier.name))
    #   omics.wo.out <- omics[[i]]
    #   for (k in 1:length(omics.wo.out)) {
    #     omics.wo.out[[k]] <- omics.wo.out[[k]][,!(colnames(omics.wo.out[[k]]) %in% temp.out)]
    #   }
    #   meta.wo.out <- meta.list[[i]][!(meta.list[[i]]$common_name %in% temp.out),]
    #   run_contrasts2(contrasts[j], gmt.list1 = "chr8", gmt.list2 = "chr8", 
    #                  meta.df = meta.wo.out, 
    #                  omics = omics.wo.out, types = names(omics[[i]]), 
    #                  feature.names = feature.list[[i]], 
    #                  EA.types = c("KEGG", "Hallmark", "Positional", "Positional_Custom"),
    #                  base.path = base.path, 
    #                  temp.path = temp.path, subfolder = FALSE, 
    #                  synapse_id = noOutFolder) 
    #   
    #   # also run with only samples matching across PDX/tumor
    #   dir.create("PDX_Chr8_status_matches_tumor")
    #   setwd("PDX_Chr8_status_matches_tumor")
    #   noOutMatchFolder <- synapser::synStore(synapser::Folder("PDX_Chr8_status_matches_tumor", noOutFolder))
    #   temp.path <- file.path(base.path, contrasts[j], names(omics)[i], 
    #                          "no_outliers", paste0("no_", outlier.name), 
    #                          "PDX_Chr8_status_matches_tumor")
    #   omics.wo.out.match <- omics.wo.out
    #   omics.wo.out <- NULL
    #   for (k in 1:length(omics.wo.out.match)) {
    #     match.samples <- c(feature.list[[i]][k], match[[i]])
    #     omics.wo.out.match[[k]] <- omics.wo.out.match[[k]][,match.samples]
    #   }
    #   meta.wo.out.match <- meta.wo.out[match[[i]],]
    #   run_contrasts2(contrasts[j], gmt.list1 = "chr8", gmt.list2 = "chr8", 
    #                  meta.df = meta.wo.out.match, 
    #                  omics = omics.wo.out.match, types = names(omics[[i]]), 
    #                  feature.names = feature.list[[i]], 
    #                  EA.types = c("KEGG", "Hallmark", "Positional", "Positional_Custom"),
    #                  base.path = base.path, 
    #                  temp.path = temp.path, subfolder = FALSE, 
    #                  synapse_id = noOutMatchFolder)
    #   omics.wo.out.match <- NULL
    # }
  }
}
# old prot didn't have 2+ sets for positional gene sets
# for RNA, had to include tumor WU-561 even though all NaN

# Loading PRISM drug sensitivity AUC scores
# [1] "Running Chr8_Strongly_amplified_vs_Not_amplified with no filter"
# Running ssGSEA using phospho_ksdb data
# Running enrichment analysis...
# Loading required namespace: snow
# Loading required namespace: doSNOW
# Error in `$<-.data.frame`(`*tmp*`, "N_drugs", value = NA) : 
#   replacement has 1 row, data has 0

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

cc.df.global <- global.meta.df2[,annotations]
cc.df.rna <- pca.meta.df2[,annotations]
global.heatmaps <- make_heatmaps(global.df, cc.df.global, top.gmt = genesets, fontsize=6)
rna.heatmaps <- make_heatmaps(pdxRNA, cc.df.rna, top.gmt = genesets, fontsize=6)
tumorRna.heatmaps <- make_heatmaps(tumorRNA, cc.df.rna, top.gmt = genesets, fontsize=6)
heatmaps <- list("PDX_Global_Proteomics" = global.heatmaps,
                 "PDX_RNA-Seq" = rna.heatmaps,
                 "Tumor_RNA-Seq" = tumorRna.heatmaps)
dir.create("heatmaps_of_interest")
setwd("heatmaps_of_interest")
dir.create("fontsize_6")
setwd("fontsize_6")
save_to_synapse(heatmaps)


# compile differential expression results
# Chr8_strongly_amp_vs_not: new database
temp.base <- paste0("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/",
             "Chr8_strongly_amp_vs_not/new\ database/")
setwd(temp.base)
dir.create("heatmaps")
setwd("heatmaps")
diffexp <- read.csv("global/Differential_expression/Differential_expression_results.csv")
phospho_diffexp <- read.csv("phospho/Differential_expression/Differential_expression_results.csv")
abs.max.Log2FC <- 0
abs.max.minusLogFDR <- 0
omics <- list("Global" = diffexp,
              "Phospho" = phospho_diffexp)
for (i in 1:length(omics)) {
  if (names(omics)[i] == "Phospho") {
    feature.names <- "SUB_SITE"
  } else {
    feature.names <- "Gene"
  }
  DEG.df <- omics[[i]]
  DEG.df$minusLogFDR <- -log(DEG.df$adj.P.Val, 10)
  for (j in 1:length(genesets)) {
    # filter for gene set
    temp.DEG.df <- DEG.df[DEG.df[,feature.names] %in% genesets[[j]], ]
    if (nrow(temp.DEG.df) > 0) {
      # generate data frames for heatmaps
      Log2FC.df <- temp.DEG.df[,c(feature.names, "Log2FC")]
      minusLogFDR.df <- temp.DEG.df[,c(feature.names, "minusLogFDR")]
      write.csv(Log2FC.df, paste0(names(omics)[i], "_", names(genesets)[j], "_Log2FC.csv"), row.names = FALSE)
      write.csv(minusLogFDR.df, paste0(names(omics)[i], "_", names(genesets)[j], "_minusLogFDR.csv"), row.names = FALSE)
      
      if (max(abs(temp.DEG.df$Log2FC), na.rm = TRUE) > abs.max.Log2FC) {
        abs.max.Log2FC <- max(abs(temp.DEG.df$Log2FC), na.rm = TRUE)
      } # 9.591479, so using range of -10, 10 for heatmap colors
      
      if (max(abs(temp.DEG.df$minusLogFDR), na.rm = TRUE) > abs.max.minusLogFDR) {
        abs.max.minusLogFDR <- max(abs(temp.DEG.df$minusLogFDR), na.rm = TRUE)
      } # 0.1577173, so using range of 0, 0.2 for heatmap sizes 
    }
  }
}

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
setwd(base.path)
dir.create("heatmaps")

global.samples <- meta.df[meta.df$id %in% colnames(global.df),]$SampleID
colnames(global.df) <- c(global.samples, "Gene")
global.df <- global.df[, c("Gene", global.samples)]

global.w.mouse.samples <- meta.df[meta.df$id %in% colnames(global.w.mouse.df),]$SampleID
colnames(global.w.mouse.df) <- c(global.w.mouse.samples, "Gene")
global.w.mouse.df <- global.w.mouse.df[, c("Gene", global.w.mouse.samples)]

phospho.samples <- meta.df[meta.df$id %in% colnames(phospho.df),]$SampleID
colnames(phospho.df) <- c(phospho.samples, "SUB_SITE")
phospho.df <- phospho.df[, c("SUB_SITE", phospho.samples)]

omics <- list('Global_old' = global.w.mouse.df,
              'RNA-Seq' = pdxRNA,
              'Global' = global.df,
              'phospho' = phospho.df)
all.heatmaps <- list()
dir.create("heatmaps")
setwd("heatmaps")
for (i in 1:(length(omics)-1)) {
  all.heatmaps[[names(omics)[i]]] <- make_heatmaps(omics[[i]], genesets)
}
save_to_synapse(all.heatmaps)



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
# pdxRNA <- read.table("salmon.merged.gene_counts.tsv", sep = "\t", 
#                      header = TRUE)
# pdxRNA$gene_id <- NULL
# colnames(pdxRNA) <- c("Gene", "JH-2-009", "JH-2-055", "JH-2-079", "WU-487",
#                       "WU-536", "WU-561", "MN-2", "MN-3")
RNA.overlap <- pdxRNA[ , c("Gene", "JH-2-055", "JH-2-079", "WU-487", "MN-2")]

omics <- list('PDX_global' = global.w.mouse.df,
              'PDX_RNAseq' = pdxRNA,
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
