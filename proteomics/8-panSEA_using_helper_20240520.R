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
# global.w.mouse.df <- read.table(
#   "global_with_mouse_data/Chr8_crosstab_global_old_database_gene_corrected.txt", 
#   sep = "\t")
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

# cnv<-do.call(rbind,lapply(setdiff(combined$CopyNumber,NA),function(x){
# 
#   x2=x#gsub('"','',gsub("[",'',gsub("]",'',x,fixed=T),fixed=T),fixed=T)
#   print(x)
#   sample<-subset(combined,CopyNumber==x)
#   print(sample$improve_sample_id)
#   res<-data.table::fread(synGet(x2)$path)
# 
#   long_df<- res|>
#     tidyr::separate_rows(gene,sep=',')|>
#     dplyr::rename(gene_symbol='gene')|>
#     dplyr::left_join(joined.df)|>
#     subset(!is.na(entrez_id))|>
#     dplyr::select(entrez_id,log2)|>
#     dplyr::distinct()|>
#     dplyr::mutate(copy_number=2^log2)
# 
#   res<-long_df|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
#     dplyr::mutate(copy_call=ifelse(copy_number<0.5210507,'deep del',
#                                    ifelse(copy_number<0.7311832,'het loss',
#                                           ifelse(copy_number<1.214125,'diploid',
#                                                  ifelse(copy_number<1.422233,'gain','amp')))))|>
#     mutate(study='MPNST PDX MT',source='NF Data Portal',improve_sample_id=sample$improve_sample_id[1])|>
#     dplyr::distinct()
# 
#   # long_df <- res[, strsplit(as.character(gene), ","), by = .(chromosome, start, end, depth, log2)]
#   # filtered_df <- long_df |>
#   #     subset(is.finite(log2))|>
#   #     filter(V1 %in% genes_df$gene) # get only protein coding genes and remove empty gene symbols
#   # filtered_df <- filtered_df[, .(gene_symbol = V1,
#   #                        improve_sample_id = sample$improve_sample_id[1],
#   #                        copy_number = 2^log2,
#   #                        source = "NF Data Portal",
#   #                        study = "MPNST PDX MT")]
#   # res<-filtered_df|> ##deep del < 0.5210507 < het loss < 0.7311832 < diploid < 1.214125 < gain < 1.422233 < amp
#   #     dplyr::mutate(copy_call=ifelse(copy_number<0.5210507,'deep del',
#   #                                    ifelse(copy_number<0.7311832,'het loss',
#   #                                           ifelse(copy_number<1.214125,'diploid',
#   #                                           ifelse(copy_number<1.422233,'gain','amp')))))|>
#   #     left_join(genes_df)|>
#   #     dplyr::select(entrez_id,improve_sample_id,copy_number,copy_call,study,source)|>
#   #     subset(!is.na(entrez_id))|>
#   #     distinct()
#   # res|>group_by(copy_call)|>summarize(n_distinct(entrez_id))
#   return(res)
#   # }
# }))

valid.chr8 <- c("Not amplified", "Amplified", "Strongly amplified")
Chr8.samples <- manifest[manifest$PDX_Chr8_status_BG %in% valid.chr8 | 
                           manifest$Tumor_Chr8_status_BG %in% valid.chr8,]$common_name # 13
rnaseq <- rnaseq[rnaseq$sample %in% Chr8.samples,]
rna <- reshape2::dcast(rnaseq, sampleType + gene_symbol ~ sample, mean, value.var = "TPM")
colnames(rna)[2] <- "Gene"

pdxRNA <- rna[rna$sampleType=="PDX", colnames(rna)[2:ncol(rna)]]
tumorRNA <- rna[rna$sampleType=="Tumor", colnames(rna)[2:ncol(rna)]]

#cnv <- cnv[cnv$sample %in% Chr8.samples,]

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)
#global.w.mouse.df$Gene <- rownames(global.w.mouse.df)
phospho.df$SUB_SITE <- rownames(phospho.df)

#### 1a. Check Chr8 copy number to assign labels ####


#### 2. Venn diagram & PCA plots ####
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

# only keep RNA-Seq samples also in proteomics
common.names <- unique(global.meta.df2$Sample)
pdxRNA <- pdxRNA[,c("Gene", common.names)]
tumorRNA <- tumorRNA[,c("Gene", common.names)]
pdxRNA50 <- pdxRNA50[,c("Gene", common.names)]
tumorRNA50 <- tumorRNA50[,c("Gene", common.names)]

omics <- list("Proteomics" = list("Global" = global.df, "Phospho" = phospho.df),
             "RNA-Seq" = list("PDX" = pdxRNA,"Tumor" = tumorRNA),
             "RNA-Seq (After Missingness Filter)" = list("PDX" = pdxRNA50,"Tumor" = tumorRNA50))
prot.feat <- c("Gene", "SUB_SITE")
rna.feat <- c("Gene", "Gene")

meta.list <- list("Proteomics" = global.meta.df2,
                  "RNA-Seq" = rna.meta.df,
                  "RNA-Seq (After Missingness Filter)" = rna.meta.df)
expr.list <- list("Proteomics" = list("CCLE proteomics"),
                  "RNA-Seq" = list("adherent CCLE", "adherent CCLE"),
                  "RNA-Seq (After Missingness Filter)" = list("adherent CCLE", "adherent CCLE"))
outliers <- list("Proteomics" = c("JH-2-002_rep1", "JH-2-002_rep2"),
                 "RNA-Seq" = c(),
                 "RNA-Seq (After Missingness Filter)" = c())
Chr8.match.ids <- global.meta.df2[global.meta.df2$Sample != "MN-2",]$SampleName
Chr8.matching <- global.meta.df2[global.meta.df2$Sample != "MN-2" &
                                   global.meta.df2$Sample != "JH-2-023" &
                                   global.meta.df2$Sample != "WU-386" & 
                                   global.meta.df2$Sample != "WU-486",]$Sample
match <- list("Proteomics" = Chr8.match.ids,
                 "RNA-Seq" = Chr8.matching,
              "RNA-Seq (After Missingness Filter)" = Chr8.matching)
feature.list <- list("Proteomics" = prot.feat,
                     "RNA-Seq" = rna.feat,
                     "RNA-Seq (After Missingness Filter)" = rna.feat)
contrasts <- c("Chr8")
annotations <- c("Sex", "Sample")
contrasts <- c("Chr8_v2")
outliers <- list("Proteomics" = c(),
                 "RNA-Seq" = c("JH-2-002"),
                 "RNA-Seq (After Missingness Filter)" = c("JH-2-002"))

# global.meta.df2$Chr8_inclusive <- global.meta.df2$PDX_Chr8_status_BG
# global.meta.df2[global.meta.df2$Sample != "WU-225",]$Chr8_inclusive <- "Not_amplified"
# rna.meta.df$Chr8_inclusive <- rna.meta.df$Chr8
# rna.meta.df[rna.meta.df$Sample != "WU-225",]$Chr8_inclusive <- "Not_amplified"
# contrasts <- c("Chr8_inclusive")
# outliers <- list("Proteomics" = c(),
#                  "RNA-Seq" = list("PDX" = c("JH-2-031", "WU-487"), "Tumor" = c("JH-2-002", "JH-2-023")),
#                  "RNA-Seq (After Missingness Filter)" = list("PDX" = c("JH-2-031", "WU-487"), "Tumor" = c("JH-2-002", "JH-2-023")))

Chr8.v3.amp <- c("JH-2-079c", "WU-225")
global.meta.df2$Chr8_v3 <- "Not_amplified"
global.meta.df2[global.meta.df2$Sample %in% Chr8.v3.amp,]$Chr8_v3 <- "Amplified"
rna.meta.df$Chr8_v3 <- "Not_amplified"
rna.meta.df[rna.meta.df$Sample %in% Chr8.v3.amp,]$Chr8_v3 <- "Amplified"
contrasts <- c("Chr8_v3")
outliers <- list("Proteomics" = c("MN-2"),
                 "RNA-Seq" = list("PDX" = c("MN-2"), "Tumor" = c("JH-2-079c")),
                 "RNA-Seq (After Missingness Filter)" = list("PDX" = c("MN-2"), "Tumor" = c("JH-2-079c")))

# get Chr8 set annotations
#gmt1 <- get_chr8_gmt1()
gmt1 <- readRDS("gmt1.rds")
names(gmt1) <- c("KEGG", "Hallmark", "Positional", "Positional_custom")
gmt2 <- readRDS("gmt2.rds")

syn_base <- "syn58649770"
synapser::synLogin()
#debug(run_contrasts2)
for (j in 1:length(contrasts)){
  setwd(base.path)
  dir.create(contrasts[j])
  setwd(contrasts[j])
  contrastFolder <- synapser::synStore(synapser::Folder(contrasts[j], syn_base))
  
  for (i in 2:2) {
    setwd(file.path(base.path, contrasts[j]))
    dir.create(names(omics)[i])
    setwd(names(omics)[i])
    omicsFolder <- synapser::synStore(synapser::Folder(names(omics)[i], contrastFolder))
    contrast.meta <- meta.list[[i]]
    contrast.meta <- na.omit(contrast.meta[,c(contrasts[j], annotations)])
    
    # # run contrast w/ sex factor
    # dir.create("no_filter")
    # setwd("no_filter")
    # noFiltFolder <- synapser::synStore(synapser::Folder("no_filter", omicsFolder))
    # temp.path <- file.path(base.path, contrasts[j], names(omics)[i],
    #                        "no_filter")
    # scale <- TRUE
    # cluster <- TRUE
    # if (grepl("RNA-Seq", names(omics)[i])) {
    #   scale <- FALSE
    #   cluster <- FALSE
    # }
    #   panSEA2(contrasts[j], "Sex", contrast.meta, omics = omics[[i]],
    #           annotations = annotations,
    #           gmt.list1 = gmt1, gmt.list2=gmt2,
    #           expr = expr.list[[i]], base.path = base.path, temp.path = temp.path,
    #           subfolder = FALSE, synapse_id = noFiltFolder, n.net = 10, scale = scale, cluster = cluster)
    
    # # re-run with only samples for which PDX Chr8 status matches tumor
    # if (!grepl("Proteomics", names(omics)[i])) {
    #   setwd(file.path(base.path, contrasts[j], names(omics)[i]))
    #   dir.create("PDX_Chr8_status_matches_tumor")
    #   setwd("PDX_Chr8_status_matches_tumor")
    #   matchFolder <- synapser::synStore(synapser::Folder("PDX_Chr8_status_matches_tumor", omicsFolder))
    #   temp.path <- file.path(base.path, contrasts[j], names(omics)[i],
    #                          "PDX_Chr8_status_matches_tumor")
    #   omics.match <- omics[[i]]
    #   for (k in 1:length(omics.match)) {
    #     match.samples <- c(feature.list[[i]][k], match[[i]])
    #     omics.match[[k]] <- omics.match[[k]][,match.samples]
    #   }
    #   meta.match <- contrast.meta[match[[i]],]
    #   
    #   scale <- TRUE
    #   cluster <- TRUE
    #   if (grepl("RNA-Seq", names(omics)[i])) {
    #     scale <- FALSE
    #     cluster <- FALSE
    #     }
    #   panSEA2(contrasts[j], "Sex", meta.match, omics = omics.match,
    #           annotations = annotations, 
    #           gmt.list1 = gmt1, gmt.list2=gmt2, 
    #           expr = expr.list[[i]], base.path = base.path, temp.path = temp.path, 
    #           subfolder = FALSE, synapse_id = matchFolder, n.net = 10, scale = scale, cluster = cluster) 
    #   omics.match <- NULL
    # } else if (contrasts[j] != "Chr8") {
    #   setwd(file.path(base.path, contrasts[j], names(omics)[i]))
    #   dir.create("PDX_Chr8_status_matches_tumor")
    #   setwd("PDX_Chr8_status_matches_tumor")
    #   matchFolder <- synapser::synStore(synapser::Folder("PDX_Chr8_status_matches_tumor", omicsFolder))
    #   temp.path <- file.path(base.path, contrasts[j], names(omics)[i],
    #                          "PDX_Chr8_status_matches_tumor")
    #   omics.match <- omics[[i]]
    #   for (k in 1:length(omics.match)) {
    #     match.samples <- c(feature.list[[i]][k], match[[i]])
    #     omics.match[[k]] <- omics.match[[k]][,match.samples]
    #   }
    #   meta.match <- na.omit(contrast.meta[match[[i]],])
    #   
    #   panSEA2(contrasts[j], "Sex", meta.match, omics = omics.match,
    #           annotations = annotations, 
    #           gmt.list1 = gmt1, gmt.list2=gmt2, 
    #           expr = expr.list[[i]], base.path = base.path, temp.path = temp.path, 
    #           subfolder = FALSE, synapse_id = matchFolder, n.net = 10) 
    #   omics.match <- NULL
    # }
    # re-run without outliers and w/wo samples which don't match across PDX/tumor
    setwd(file.path(base.path, contrasts[j], names(omics)[i]))
    dir.create("no_outliers")
    setwd("no_outliers")
    if (is.list(outliers[[i]])) {
      for (k in 1:length(outliers[[i]])) {
        temp.out <- outliers[[i]][[k]]
        outlier.name <- paste0(temp.out, collapse="_")
        setwd(file.path(base.path, contrasts[j], names(omics)[i], "no_outliers"))
        dir.create(paste0("no_", outlier.name))
        setwd(paste0("no_", outlier.name))
        noOutFolder <- synapser::synStore(synapser::Folder(paste0("no_", outlier.name), omicsFolder))
        temp.path <- file.path(base.path, contrasts[j], names(omics)[i],
                               "no_outliers", paste0("no_", outlier.name))
        omics.wo.out <- omics[[i]]
        for (m in 1:length(omics.wo.out)) {
          omics.wo.out[[m]] <- omics.wo.out[[m]][,!(colnames(omics.wo.out[[m]]) %in% temp.out)]
        }
        meta.wo.out <- meta.list[[i]][!(rownames(meta.list[[i]]) %in% temp.out),]
        if (grepl("RNA", names(omics)[i])) {
          scale <- FALSE
          cluster <- FALSE
        } else {
          scale <- TRUE
          cluster <- TRUE
        }
        panSEA2(contrasts[j], "Sex", meta.wo.out, omics = omics.wo.out,
                annotations = annotations, 
                gmt.list1 = gmt1, gmt.list2=gmt2, 
                expr = expr.list[[i]], base.path = base.path, temp.path = temp.path, 
                subfolder = FALSE, synapse_id = noOutFolder, n.net = 10, scale = scale, cluster = cluster)

        # also run with only samples matching across PDX/tumor
        # dir.create("PDX_Chr8_status_matches_tumor")
        # setwd("PDX_Chr8_status_matches_tumor")
        # noOutMatchFolder <- synapser::synStore(synapser::Folder("PDX_Chr8_status_matches_tumor", noOutFolder))
        # temp.path <- file.path(base.path, contrasts[j], names(omics)[i],
        #                        "no_outliers", paste0("no_", outlier.name),
        #                        "PDX_Chr8_status_matches_tumor")
        # omics.wo.out.match <- omics.wo.out
        # omics.wo.out <- NULL
        # for (m in 1:length(omics.wo.out.match)) {
        #   match.samples <- c(feature.list[[i]][m], match[[i]])
        #   omics.wo.out.match[[m]] <- omics.wo.out.match[[m]][,colnames(omics.wo.out.match[[m]]) %in% match.samples]
        # }
        # meta.wo.out.match <- meta.wo.out[meta.wo.out$id %in% match[[i]],]
        # panSEA2(contrasts[j], "Sex", meta.wo.out.match, omics = omics.wo.out.match,
        #         annotations = annotations, 
        #         gmt.list1 = gmt1, gmt.list2=gmt2, 
        #         expr = expr.list[[i]], base.path = base.path, temp.path = temp.path, 
        #         subfolder = FALSE, synapse_id = noOutMatchFolder, n.net = 10, scale = scale, cluster = cluster)
        # omics.wo.out.match <- NULL
      }
    } else {
      temp.out <- outliers[[i]]
      if (length(temp.out) > 0) {
        outlier.name <- paste0(temp.out, collapse="_")
        dir.create(paste0("no_", outlier.name))
        setwd(paste0("no_", outlier.name))
        noOutFolder <- synapser::synStore(synapser::Folder(paste0("no_", outlier.name), omicsFolder))
        temp.path <- file.path(base.path, contrasts[j], names(omics)[i],
                               "no_outliers", paste0("no_", outlier.name))
        omics.wo.out <- omics[[i]]
        for (k in 1:length(omics.wo.out)) {
          omics.wo.out[[k]] <- omics.wo.out[[k]][,!(colnames(omics.wo.out[[k]]) %in% temp.out)]
        }
        meta.wo.out <- meta.list[[i]][!(rownames(meta.list[[i]]) %in% temp.out),]
        if (grepl("RNA", names(omics)[i])) {
          scale <- FALSE
          cluster <- FALSE
        } else {
          scale <- TRUE
          cluster <- TRUE
        }
        panSEA2(contrasts[j], "Sex", meta.wo.out, omics = omics.wo.out,
                annotations = annotations, 
                gmt.list1 = gmt1, gmt.list2=gmt2, 
                expr = expr.list[[i]], base.path = base.path, temp.path = temp.path, 
                subfolder = FALSE, synapse_id = noOutFolder, n.net = 10, scale = scale, cluster = cluster)
        
        # also run with only samples matching across PDX/tumor
        # dir.create("PDX_Chr8_status_matches_tumor")
        # setwd("PDX_Chr8_status_matches_tumor")
        # noOutMatchFolder <- synapser::synStore(synapser::Folder("PDX_Chr8_status_matches_tumor", noOutFolder))
        # temp.path <- file.path(base.path, contrasts[j], names(omics)[i],
        #                        "no_outliers", paste0("no_", outlier.name),
        #                        "PDX_Chr8_status_matches_tumor")
        # omics.wo.out.match <- omics.wo.out
        # omics.wo.out <- NULL
        # for (k in 1:length(omics.wo.out.match)) {
        #   match.samples <- c(feature.list[[i]][k], match[[i]])
        #   temp.cols <- colnames(omics.wo.out.match[[k]])
        #   omics.wo.out.match[[k]] <- omics.wo.out.match[[k]][,temp.cols[temp.cols %in% match.samples]]
        # }
        # meta.wo.out.match <- dplyr::distinct(na.omit(meta.wo.out[match[[i]],]))
        # panSEA2(contrasts[j], "Sex", meta.wo.out.match, omics = omics.wo.out.match,
        #         annotations = annotations,
        #         gmt.list1 = gmt1, gmt.list2=gmt2,
        #         expr = expr.list[[i]], base.path = base.path, temp.path = temp.path,
        #         subfolder = FALSE, synapse_id = noOutMatchFolder, n.net = 10, scale = scale, cluster = cluster)
        # omics.wo.out.match <- NULL
      }
    }
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

#### 4. Pathways of interest: expression, log2FC heatmaps; TYK2, pSTAT3 boxplot ####
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

cc.df.global <- as.data.frame(global.meta.df2[,annotations])
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

# look at chr8q amplification
setwd(base.path)
dir.create("Chr8_quant")
setwd("Chr8_quant")
cnv <- read.csv("mpnst_copy_number.csv")
cnv.samples <- read.csv("mpnst_samples.csv")
cnv.samples <- dplyr::distinct(cnv.samples[,c("common_name", "improve_sample_id")])
prot.cnv.samples <- cnv.samples[cnv.samples$common_name %in% global.meta.df2$Sample,]
prot.cnv <- merge(prot.cnv.samples, cnv, by="improve_sample_id")
genes <- read.csv("genes.csv")
genes <- dplyr::distinct(genes[,c("entrez_id", "gene_symbol")])
prot.cnv <- merge(genes, prot.cnv, by="entrez_id")
pdxCNV <- prot.cnv[prot.cnv$study == "MPNST PDX MT",]
omics2 <- list("Copy Number" = pdxCNV,
               "RNA-Seq" = pdxRNA,
               "Proteomics" = global.df)
dir.create("Single_sample_positional_GSEA")
setwd("Single_sample_positional_GSEA")
library(plyr); library(dplyr)
library(ggplot2)
gmt <- get_gmt1(gmt.list1 = "msigdb_Homo sapiens_C1")
for (i in 1:length(omics2)) {
  # setwd(file.path(base.path, "Single_sample_positional_GSEA"))
  # dir.create(names(omics2)[i])
  # setwd(names(omics2)[i])
  # run single-sample GSEA for positional gene sets
  # pos.GSEA <- list()
  # if (names(omics2)[i] == "Copy Number") {
  #   temp.samples <- unique(omics2[[i]]$common_name)
  #   
  #   for (j in temp.samples) { # assumes first column is feature names
  #     # prep input data
  #     temp.input <- omics2[[i]][omics2[[i]]$common_name == j,c("gene_symbol","copy_number")]
  #     if (any(duplicated(temp.input$gene_symbol))) {
  #       temp.input <- plyr::ddply(temp.input, .(gene_symbol), summarize,
  #                                 copy_number = mean(copy_number, na.rm = TRUE))
  #     }
  #     pos.GSEA[[as.character(j)]] <- panSEA::ssGSEA(temp.input, gmt = gmt[[1]])$result
  #   }
  # } else {
  #   for (j in 2:ncol(omics2[[i]])) { # assumes first column is feature names
  #     pos.GSEA[[colnames(omics2[[i]])[j]]] <- panSEA::ssGSEA(omics2[[i]][,c(1, j)],
  #                                                            gmt = gmt[[1]])$result
  #   }
  # }
  # all.pos.GSEA <- data.table::rbindlist(pos.GSEA, use.names = TRUE, idcol = "Sample")
  # write.csv(all.pos.GSEA, paste0(names(omics2)[i], "_positional_GSEA_results.csv"), row.names = FALSE)
  all.pos.GSEA <- read.csv(paste0(names(omics2)[i], "_positional_GSEA_results.csv"))
  
  # bar graph of NES for each position of chr8q 
  chr8q.df <- all.pos.GSEA[grepl("chr8q",all.pos.GSEA$Feature_set),]
  chr8q.df$"Chr8q Position" <- as.factor(sub("chr8q", "", chr8q.df$Feature_set))
  ggplot2::ggplot(chr8q.df, aes(fill=Sample, x=`Chr8q Position`, y=NES)) + 
    ggplot2::geom_bar(stat="identity", position="dodge") + ggplot2::theme_minimal()
  ggplot2::ggsave(paste0(names(omics2)[i], "_Chr8q_positional_GSEA_results.pdf"))
  
  # median chr8q NES
  med.chr8q <- plyr::ddply(chr8q.df, .(Sample), summarize,
                           `Chr8q Median NES` = median(NES, na.rm = TRUE))
  pdf(paste0(names(omics2)[i], "_Chr8q_median_NES_histogram.pdf"))
  hist(med.chr8q$`Chr8q Median NES`, xlab = "Chr8q Median NES", main = NULL)
  dev.off()
  med.chr8q <- med.chr8q[order(med.chr8q$`Chr8q Median NES`),]
  med.chr8q$'Sample ID' <- med.chr8q$Sample
  for (j in 1:nrow(med.chr8q)) {
    med.chr8q$Sample[j] <- stringr::str_split(med.chr8q$`Sample ID`[j], "_")[[1]][1]
  }
  ggplot2::ggplot(med.chr8q, aes(x=`Sample ID`, y=`Chr8q Median NES`, fill = Sample)) + 
    ggplot2::geom_bar(stat="identity", position="dodge") + 
    scale_x_discrete(limits = med.chr8q$`Sample ID`) + ggplot2::theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
  ggplot2::ggsave(paste0(names(omics2)[i], "_Chr8q_median_NES.pdf"))
  write.csv(med.chr8q, paste0(names(omics2)[i], "_Chr8q_median_NES.csv"), row.names = FALSE)
}

# get annotations as df
dir.create("Single_sample_positional_medians")
setwd("Single_sample_positional_medians")
library(plyr); library(dplyr)
library(ggplot2)
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")
reduced.msigdb <- dplyr::distinct(msigdb.info[,c("gs_name", "gene_symbol")])
for (i in 1:length(omics2)) {
  setwd(file.path(base.path, "Chr8_quant", "Single_sample_positional_medians"))
  dir.create(names(omics2)[i])
  setwd(names(omics2)[i])
  pos.GSEA <- list()
  if (names(omics2)[i] == "Copy Number") {
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
  
  # bar graph of NES for each position of chr8q 
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
  pdf(paste0(names(omics2)[i], "_MYC_histogram.pdf"))
  hist(omics2[[i]][omics2[[i]]$Gene == "MYC",], xlab = "MYC Expression", main = NULL)
  dev.off()
  med.chr8q <- omics2[[i]][omics2[[i]]$Gene == "MYC",]
  
  # reshape long if necessary
  if (names(omics2)[i] == "Copy Number") {
    temp.input <- omics2[[i]][,c("common_name", "gene_symbol", "copy_number")]
  } else {
    temp.input <- reshape2::dcast()
  }
  
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
}



# rank genes based on correlation between expression & Chr8q NES (points = samples)
# cnv.pos.GSEA <- read.csv("Copy Number_positional_GSEA_results.csv")
# global.cnv.NES <- cnv.pos.GSEA[cnv.pos.GSEA$Sample %in% global.meta.df2$Sample,]
# rownames(global.cnv.NES) <- global.cnv.NES$Sample
# global.meta.df2$Chr8_NES <- NA
# for (i in 1:nrow(global.meta.df2)) {
#   global.meta.df2$Chr8_NES[i] <- global.cnv.NES[global.meta.df2$Sample[i] & global.cnv.NES$Feature_set == "chr8q13",]$NES
# }
cnv.med.chr8q <- read.csv("Copy Number_Chr8q_median_NES.csv")
global.meta.df2$Chr8q_median_NES <- NA
for (i in 1:nrow(global.meta.df2)) {
  temp.sample <- global.meta.df2$Sample[i]
  global.meta.df2$Chr8q_median_NES[i] <- cnv.med.chr8q[cnv.med.chr8q$Sample == temp.sample,]$Chr8q.Median.NES
}
rna.meta.df$Chr8q_median_NES <- NA
for (i in 1:nrow(rna.meta.df)) {
  temp.sample <- rna.meta.df$Sample[i]
  rna.meta.df$Chr8q_median_NES[i] <- cnv.med.chr8q[cnv.med.chr8q$Sample == temp.sample,]$Chr8q.Median.NES
}
#global.meta.df3 <- merge(cnv.med.chr8q, global.meta.df2, by= "Sample", all.y=TRUE)

for (i in 1:length(omics)) {
  setwd(file.path(base.path, "Chr8_quant"))
  dir.create(names(omics)[i])
  setwd(names(omics)[i])
  temp.features <- feature.list[[i]]
  temp.omics <- omics[[i]]
  meta.df <- meta.list[[i]]
  red.meta <- meta.df[,c("Sample", "Chr8q_median_NES")]
  red.meta$Sample <- rownames(red.meta)
  
  omics.files <- list()
  for (j in 1:length(temp.omics)) {
    temp.omics.files <- list()
    omics.t <- as.data.frame(t(temp.omics[[j]]))
    omics.t <- omics.t[2:nrow(omics.t),]
    omics.t$Sample <- rownames(omics.t)
    # need to make sample names first col, gene symbols colnames
    omics.t <- merge(red.meta, omics.t)
    
    deg <- DMEA::rank_corr(omics.t, variable=temp.features[j], value = names(temp.omics)[j], plots = FALSE)
    temp.omics.files[["Differential_expression"]] <- list("DEG_correlation_results.csv" = deg$result)
    
    # GSEA
    if (grepl("phospho", names(temp.omics)[j], ignore.case = TRUE)) {
      names(gmt2) <- c("phospho_ksdb", "phospho_sub")
      gsea2 <- panSEA::mGSEA(list(deg$result, deg$result), gmt2, 
                                    types = names(gmt2), 
                                    feature.names = rep(temp.features[j],2),
                                    rank.var = rep("Pearson.est",2))
      # store GSEA results
      phospho.gsea.files <- list()
      for (k in 1:length(gmt2)) {
        gsea.name <- paste0("GSEA_", names(gmt2)[k])
        phospho.mtn <- get_top_mtn_plots(
          gsea2$all.results[[k]], 
          EA.type = names(gmt2)[k])
        if (length(phospho.mtn) > 1) {
          phospho.net <- panSEA::netSEA(list(deg$result),
                                       list(gsea2$all.results[[k]]$result),
                                       element.names = temp.features[j],
                                       rank.var = "Pearson.est",
                                       n.network.sets = 10)
          phospho.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                   gsea2$all.results[[k]]$result,
                                                 "GSEA_volcano_plot.pdf" =
                                                   gsea2$all.results[[k]]$volcano.plot,
                                                 "GSEA_network_graph.html" = 
                                                   phospho.net$interactive,
                                                 "mtn_plots" = phospho.mtn)
        } else {
          phospho.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                   gsea2$all.results[[k]]$result,
                                                 "GSEA_volcano_plot.pdf" =
                                                   gsea2$all.results[[k]]$volcano.plot,
                                                 "mtn_plots" = phospho.mtn)
        } 
        if (length(phospho.mtn) > 0) {
          phospho.gsea.files[[gsea.name]][["Pathways_of_interest"]] <- 
            get_pathways_of_interest(omics[[names(temp.omics)[j]]], 
                                     gsea2$all.results[[k]]$result, 
                                     gmt2[[k]], cc.df.global, n = 10,
                                     show_colnames = FALSE, 
                                     fontsize = 10, scale = TRUE, cluster = TRUE) 
        }
      }
      temp.omics.files[["Phospho_enrichment"]] <- phospho.gsea.files
      
      # prep for GSEA
      gsea1.inputs <- list(gsea2$all.results[[2]]$result)
      names(gsea1.inputs) <- names(gmt2)[2]
      features1 <- "Feature_set"
      rank.var <- "NES"
    } else {
      # prep for GSEA
      gsea1.inputs <- list(deg$result)
      names(gsea1.inputs) <- names(temp.omics)[j]
      features1 <- temp.features[j]
      rank.var <- "Pearson.est"
    }
    
    # run GSEA
    gsea1 <- list()
    all.global.gsea.files <- list()
    for (k in 1:length(gmt1)) {
      gsea.name <- paste0("GSEA_", names(gmt1)[k])
      gsea1[[gsea.name]] <- panSEA::mGSEA(gsea1.inputs, list(gmt1[[k]]), 
                                               types = names(gsea1.inputs), 
                                               feature.names = rep(features1,1),
                                               rank.var = rep(rank.var,1)) 
      # cannot run GSEA on KSEA results
      global.gsea.files <- list()
      for (m in 1:length(gsea1.inputs)) {
        temp.gsea.files <- list()
        global.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[[m]], 
          EA.type = names(gmt1)[k])
        if (length(global.mtn) > 1) {
          global.net <- panSEA::netSEA(list(gsea1.inputs[[m]]),
                                       list(gsea1[[gsea.name]]$all.results[[m]]$result),
                                       element.names = features1,
                                       rank.var = rank.var,
                                       n.network.sets = 10)
          temp.gsea.files <- list("GSEA_results.csv" =
                                                   gsea1[[gsea.name]]$all.results[[m]]$result,
                                                 "GSEA_volcano_plot.pdf" =
                                                   gsea1[[gsea.name]]$all.results[[m]]$volcano.plot,
                                                 "GSEA_network_graph.html" = 
                                                   global.net$interactive,
                                                 "mtn_plots" = global.mtn)
        } else {
          temp.gsea.files <- list("GSEA_results.csv" =
                                                   gsea1[[gsea.name]]$all.results[[m]]$result,
                                                 "GSEA_volcano_plot.pdf" =
                                                   gsea1[[gsea.name]]$all.results[[m]]$volcano.plot,
                                                 "mtn_plots" = global.mtn)
        }
        if (!grepl("phospho", names(temp.omics)[j], ignore.case = TRUE) & length(global.mtn) > 0) {
          temp.gsea.files[["Pathways_of_interest"]] <- 
            get_pathways_of_interest(omics[[names(temp.omics)[j]]], 
                                     gsea1[[gsea.name]]$all.results[[m]]$result, 
                                     gmt1[[k]], cc.df.global, n = 10,
                                     show_colnames = FALSE, 
                                     fontsize = 10, scale = TRUE, cluster = TRUE) 
        }
        global.gsea.files[[names(gsea1.inputs)[m]]] <- temp.gsea.files
      }
      all.global.gsea.files[[gsea.name]] <- global.gsea.files
    }
    temp.omics.files[["GSEA"]] <- global.gsea.files
    
    # DMEA
    if (names(omics)[i] == "Proteomics") {
      temp.expr <- read.csv(file.path(base.path, "CCLE_proteomics.csv"))
      temp.expr <- list(temp.expr)
    } else {
      temp.expr <- list("adherent CCLE")
    }
    dmea.results <- panSEA::mDMEA(expression = temp.expr, weights = gsea1.inputs, 
                                  types = names(temp.omics)[j], 
                                  feature.names = features1,
                                  weight.values = rank.var)
    DMEA.global.mtn <- get_top_mtn_plots(dmea.results$all.results[[1]],
                                         sets = "Drug_set",
                                         EA.type = "DMEA")
    DMEA.global.net <- panSEA::netSEA(list(dmea.results$all.results[[1]]$corr.result),
                                      list(dmea.results$all.results[[1]]$result),
                                      "Drug", "Pearson.est",
                                      n.network.sets = 10)
    global.DMEA.files <- list("DMEA_results.csv" =
                                dmea.results$all.results[[1]]$result,
                              "DMEA_correlation_results.csv" = 
                                dmea.results$all.results[[1]]$corr.result,
                              "DMEA_correlation_scatter_plots.pdf" = 
                                dmea.results$all.results[[1]]$corr.scatter.plots,
                              "DMEA_volcano_plot.pdf" =
                                dmea.results$all.results[[1]]$volcano.plot,
                              "DMEA_network_graph.html" = 
                                DMEA.global.net$interactive,
                              "mtn_plots" = DMEA.global.mtn)
    temp.omics.files[["DMEA"]] <- global.DMEA.files
    omics.files[[names(temp.omics)[j]]] <- temp.omics.files
  }
  save_to_synapse(omics.files)
}

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
#global.meta.df3 <- merge(cnv.med.chr8q, global.meta.df2, by= "Sample", all.y=TRUE)

for (i in 1:length(omics)) {
  setwd(file.path(base.path, "Chr8_quant"))
  dir.create(names(omics)[i])
  setwd(names(omics)[i])
  temp.features <- feature.list[[i]]
  temp.omics <- omics[[i]]
  meta.df <- meta.list[[i]]
  red.meta <- meta.df[,c("Sample", "Chr8q_median")]
  red.meta$Sample <- rownames(red.meta)
  
  omics.files <- list()
  for (j in 1:length(temp.omics)) {
    temp.omics.files <- list()
    omics.t <- as.data.frame(t(temp.omics[[j]]))
    colnames(omics.t) <- omics.t[1,]
    omics.t <- omics.t[2:nrow(omics.t),]
    omics.t$Sample <- rownames(omics.t)
    # need to make sample names first col, gene symbols colnames
    omics.t <- merge(red.meta, omics.t)
    
    deg <- DMEA::rank_corr(omics.t, variable=temp.features[j], value = names(temp.omics)[j], plots = FALSE)
    sig.degs <- na.omit(deg$result[deg$result$Pearson.q <= 0.05, ])
    
    # select numeric data
    just.expr <- dplyr::select_if(temp.omics[[j]], is.numeric)
    
    # identify feature name
    feature.name <- temp.features[j]
    
    if (nrow(sig.degs) > 0) {
      top.top.sig.degs <- sig.degs %>% slice_max(Pearson.est, n = 50/2)
      top.bot.sig.degs <- sig.degs %>% slice_min(Pearson.est, n = 50/2)
      top.sig.degs <- rbind(top.top.sig.degs, top.bot.sig.degs)
      heatmap.df <- temp.omics[[j]][temp.omics[[j]][,feature.name] %in% top.sig.degs[,feature.name],
                               c(feature.name, colnames(just.expr))]
      rownames(heatmap.df) <- heatmap.df[,feature.name]
      feature.order <- top.sig.degs[order(top.sig.degs$Pearson.est, decreasing = TRUE), feature.name]
      heatmap.df <- heatmap.df[feature.order,]
      heatmap.mat <- heatmap.df[,2:ncol(heatmap.df)]
      heatmap.mat <- filter_for_hclust(heatmap.mat)
      
      # create heatmaps
      if (nrow(heatmap.mat) > 1) {
        heatmap.mat <- as.matrix(heatmap.mat)
        deg.heatmap.clust <- pheatmap::pheatmap(heatmap.mat, color = 
                                                  colorRampPalette(
                                                    c("navy", "white", "firebrick3"))(50),
                                                scale = "row", annotation_col = cc.df.global, 
                                                angle_col = "45", 
                                                show_colnames = FALSE, 
                                                fontsize = 10)
        deg.heatmap <- pheatmap::pheatmap(heatmap.mat, 
                                          color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                          cluster_row = FALSE, 
                                          scale = "row", annotation_col = cc.df.global, 
                                          angle_col = "45", 
                                          show_colnames = FALSE, 
                                          fontsize = 10)
        deg.heatmap.abs <- pheatmap::pheatmap(heatmap.mat, color = 
                                                colorRampPalette(
                                                  c("navy", "white", "firebrick3"))(50), 
                                              annotation_col = cc.df.global, 
                                              angle_col = "45", 
                                              show_colnames = FALSE, 
                                              fontsize = 10)
      } else {
        deg.heatmap.clust <- list()
        deg.heatmap <- list()
        deg.heatmap.abs <- list()
      }
      
      temp.DEG.files <- list("Differential_expression_results.csv" = 
                               deg$result,
                             "Differential_expression_results_max_5_percent_FDR.csv" = 
                               sig.degs,
                             "Differential_expression_for_heatmap.csv" =
                               heatmap.df,
                             "Differential_expression_heatmap_not_scaled.bp" =
                               deg.heatmap.abs,
                             "Differential_expression_heatmap_scaled_ordered_by_Log2FC.bp" =
                               deg.heatmap,
                             "Differential_expression_heatmap_scaled.bp" =
                               deg.heatmap.clust) 
    } else {
      temp.DEG.files <- list("Differential_expression_results.csv" = 
                               deg$result)
    }
    temp.omics.files[["Differential_expression"]] <- temp.DEG.files
    
    # GSEA
    if (grepl("phospho", names(temp.omics)[j], ignore.case = TRUE)) {
      names(gmt2) <- c("phospho_ksdb", "phospho_sub")
      gsea2 <- panSEA::mGSEA(list(deg$result, deg$result), gmt2, 
                             types = names(gmt2), 
                             feature.names = rep(temp.features[j],2),
                             rank.var = rep("Pearson.est",2))
      # store GSEA results
      phospho.gsea.files <- list()
      for (k in 1:length(gmt2)) {
        gsea.name <- paste0("GSEA_", names(gmt2)[k])
        phospho.mtn <- get_top_mtn_plots(
          gsea2$all.results[[k]], 
          EA.type = names(gmt2)[k])
        if (length(phospho.mtn) > 1) {
          phospho.net <- panSEA::netSEA(list(deg$result),
                                        list(gsea2$all.results[[k]]$result),
                                        element.names = temp.features[j],
                                        rank.var = "Pearson.est",
                                        n.network.sets = 10)
          phospho.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                    gsea2$all.results[[k]]$result,
                                                  "GSEA_volcano_plot.pdf" =
                                                    gsea2$all.results[[k]]$volcano.plot,
                                                  "GSEA_network_graph.html" = 
                                                    phospho.net$interactive,
                                                  "mtn_plots" = phospho.mtn)
        } else {
          phospho.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                    gsea2$all.results[[k]]$result,
                                                  "GSEA_volcano_plot.pdf" =
                                                    gsea2$all.results[[k]]$volcano.plot,
                                                  "mtn_plots" = phospho.mtn)
        } 
        if (length(phospho.mtn) > 0) {
          phospho.gsea.files[[gsea.name]][["Pathways_of_interest"]] <- 
            get_pathways_of_interest(temp.omics[[names(temp.omics)[j]]], 
                                     gsea2$all.results[[k]]$result, 
                                     gmt2[[k]], cc.df.global, n = 10,
                                     show_colnames = FALSE, 
                                     fontsize = 10, scale = TRUE, cluster = TRUE) 
        }
      }
      temp.omics.files[["Phospho_enrichment"]] <- phospho.gsea.files
      
      # prep for GSEA
      gsea1.inputs <- list(gsea2$all.results[[2]]$result)
      names(gsea1.inputs) <- names(gmt2)[2]
      features1 <- "Feature_set"
      rank.var <- "NES"
    } else {
      # prep for GSEA
      gsea1.inputs <- list(deg$result)
      names(gsea1.inputs) <- names(temp.omics)[j]
      features1 <- temp.features[j]
      rank.var <- "Pearson.est"
    }
    
    # run GSEA
    gsea1 <- list()
    all.global.gsea.files <- list()
    for (k in 1:length(gmt1)) {
      gsea.name <- paste0("GSEA_", names(gmt1)[k])
      gsea1[[gsea.name]] <- panSEA::mGSEA(gsea1.inputs, list(gmt1[[k]]), 
                                          types = names(gsea1.inputs), 
                                          feature.names = rep(features1,1),
                                          rank.var = rep(rank.var,1)) 
      # cannot run GSEA on KSEA results
      global.gsea.files <- list()
      for (m in 1:length(gsea1.inputs)) {
        temp.gsea.files <- list()
        global.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[[m]], 
          EA.type = names(gmt1)[k])
        if (length(global.mtn) > 1) {
          global.net <- panSEA::netSEA(list(gsea1.inputs[[m]]),
                                       list(gsea1[[gsea.name]]$all.results[[m]]$result),
                                       element.names = features1,
                                       rank.var = rank.var,
                                       n.network.sets = 10)
          temp.gsea.files <- list("GSEA_results.csv" =
                                    gsea1[[gsea.name]]$all.results[[m]]$result,
                                  "GSEA_volcano_plot.pdf" =
                                    gsea1[[gsea.name]]$all.results[[m]]$volcano.plot,
                                  "GSEA_network_graph.html" = 
                                    global.net$interactive,
                                  "mtn_plots" = global.mtn)
        } else {
          temp.gsea.files <- list("GSEA_results.csv" =
                                    gsea1[[gsea.name]]$all.results[[m]]$result,
                                  "GSEA_volcano_plot.pdf" =
                                    gsea1[[gsea.name]]$all.results[[m]]$volcano.plot,
                                  "mtn_plots" = global.mtn)
        }
        if (!grepl("phospho", names(temp.omics)[j], ignore.case = TRUE) & length(global.mtn) > 0) {
          temp.gsea.files[["Pathways_of_interest"]] <- 
            get_pathways_of_interest(temp.omics[[names(temp.omics)[j]]], 
                                     gsea1[[gsea.name]]$all.results[[m]]$result, 
                                     gmt1[[k]], cc.df.global, n = 10,
                                     show_colnames = FALSE, 
                                     fontsize = 10, scale = TRUE, cluster = TRUE) 
        }
        global.gsea.files[[names(gsea1.inputs)[m]]] <- temp.gsea.files
      }
      all.global.gsea.files[[gsea.name]] <- global.gsea.files
    }
    temp.omics.files[["GSEA"]] <- all.global.gsea.files
    
    # DMEA
    if (names(omics)[i] == "Proteomics") {
      temp.expr <- read.csv(file.path(base.path, "CCLE_proteomics.csv"))
      temp.expr <- list(temp.expr)
    } else {
      temp.expr <- list("adherent CCLE")
    }
    dmea.results <- panSEA::mDMEA(expression = temp.expr, weights = gsea1.inputs, 
                                  types = names(gsea1.inputs), 
                                  feature.names = features1,
                                  weight.values = rank.var)
    DMEA.global.mtn <- get_top_mtn_plots(dmea.results$all.results[[1]],
                                         sets = "Drug_set",
                                         EA.type = "DMEA")
    DMEA.global.net <- panSEA::netSEA(list(dmea.results$all.results[[1]]$corr.result),
                                      list(dmea.results$all.results[[1]]$result),
                                      "Drug", "Pearson.est",
                                      n.network.sets = 10)
    global.DMEA.files <- list("DMEA_results.csv" =
                                dmea.results$all.results[[1]]$result,
                              "DMEA_correlation_results.csv" = 
                                dmea.results$all.results[[1]]$corr.result,
                              "DMEA_correlation_scatter_plots.pdf" = 
                                dmea.results$all.results[[1]]$corr.scatter.plots,
                              "DMEA_volcano_plot.pdf" =
                                dmea.results$all.results[[1]]$volcano.plot,
                              "DMEA_network_graph.html" = 
                                DMEA.global.net$interactive,
                              "mtn_plots" = DMEA.global.mtn)
    temp.omics.files[["DMEA"]] <- global.DMEA.files
    omics.files[[names(temp.omics)[j]]] <- temp.omics.files
  }
  save_to_synapse(omics.files)
}

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
FDR <- 0.25
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
