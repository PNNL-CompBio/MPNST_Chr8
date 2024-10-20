# Differential expression & enrichment analyses: PDX transcriptomics, global & phospho
# Chr8: median copy number vs. protein/phospho expression
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-09-10

# overview:
#### 1. Import metadata & crosstabs
#### 2. Run panSEA: only MPNST PDX with proteomics
#### 3. Run panSEA: all MPNST PDX
#### 4. PCSF network analysis

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr); library(plyr)
library(htmlwidgets); library(webshot); library(scales); library(msigdbr)
library(plyr); library(dplyr); library(R.utils)
library(ggplot2)
#webshot::install_phantomjs()
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
#source("panSEA_helper_20240508.R")
source("panSEA_helper_20240913.R")
source("https://raw.githubusercontent.com/PNNL-CompBio/amlresistancenetworks/master/R/proteinNetworks.R")

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
loadRNA <- function() {
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
  Chr8.samples <- manifest[manifest$PDX_Chr8_status %in% valid.chr8 | 
                             manifest$Tumor_Chr8_status %in% valid.chr8,]$common_name # 13
  rnaseq <- rnaseq[rnaseq$sample %in% Chr8.samples,]
  rna <- reshape2::dcast(rnaseq, sampleType + gene_symbol ~ sample, mean, value.var = "TPM")
  colnames(rna)[2] <- "Gene"
  
  pdxRNA <- rna[rna$sampleType=="PDX", colnames(rna)[2:ncol(rna)]]
  tumorRNA <- rna[rna$sampleType=="Tumor", colnames(rna)[2:ncol(rna)]]
  return(list(pdxRNA = pdxRNA, tumorRNA = tumorRNA))
}
RNA <- loadRNA()
pdxRNA <- RNA$pdxRNA
tumorRNA <- RNA$tumorRNA

#### 2. determine median chr8q copy number ####
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/"
setwd(base.path)
annotations <- c("Sex")

# make sure meta data has right contrasts
global.meta <- manifest[!is.na(manifest$PDX_Proteomics),
                        c("common_name", annotations)]
global.meta.df2 <- merge(global.meta, meta.df, by.x="common_name", by.y="SampleID_revised")
global.meta.df2$id <- paste0("X", global.meta.df2$id)
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

rna.meta.df <- dplyr::distinct(global.meta.df2[,c("Sample", annotations)])
rownames(rna.meta.df) <- rna.meta.df$Sample
rownames(global.meta.df2) <- global.meta.df2$SampleName

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

# actually should also incorporate WES
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/"
setwd(base.path)
dir.create("Chr8_quant")
setwd("Chr8_quant")
cnv <- read.csv("mpnst_copy_number.csv")
cnv.samples <- read.csv("mpnst_samples.csv")
cnv.samples <- dplyr::distinct(cnv.samples[,c("common_name", "improve_sample_id")])
#prot.cnv.samples <- cnv.samples[cnv.samples$common_name %in% global.meta.df2$Sample,]
#prot.cnv <- merge(prot.cnv.samples, cnv, by="improve_sample_id")
cnv <- merge(cnv.samples, cnv, by="improve_sample_id")
genes <- read.csv("genes.csv")
genes <- dplyr::distinct(genes[,c("entrez_id", "gene_symbol")])
#prot.cnv <- merge(genes, prot.cnv, by="entrez_id")
cnv <- merge(genes, cnv, by="entrez_id")
#cnv2 <- dplyr::distinct(cnv)
#cnv2 <- NULL
#pdxCNV <- prot.cnv[prot.cnv$study == "MPNST PDX MT",]
cnv <- cnv[cnv$study == "MPNST PDX MT",]

# format wide with common_name ~ gene_symbol
cnv <- reshape2::dcast(cnv, gene_symbol ~ common_name, mean, value.var = "copy_number")
colnames(cnv)[1] <- "Gene"

# get median Chr8q copy number for metadata
chr8q.info <- read.csv(synapser::synGet("syn61811211")$path)

# get other metadata
manifest<-synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
  dplyr::rename(common_name='Sample')

##for now we only have tumor and PDX data
##they each get their own sample identifier
pdx_data<-manifest|>dplyr::select(common_name,Sex,RNASeq='PDX_RNASeq',Mutations='PDX_Somatic_Mutations',CopyNumber='PDX_CNV')
pdx_data$sampleType <- "PDX"
colnames(pdx_data)[1] <- "Sample"

chr8q.info <- merge(pdx_data, chr8q.info, all.y = TRUE)
colnames(chr8q.info)[7] <- "Median Chr8q Copy Number"

# run corr analysis for copy number
cnv <- cnv[!is.na(cnv$Gene),]


# get more metadata for heatmaps
chr8q.syn <- "syn61811211"
syn.path <- synapser::synGet(chr8q.syn)$path
chr8q.info <- read.csv(syn.path)
chr8q.info <- chr8q.info[,c("Sample", "Chr8q.Median")]

pdx.info <- synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
  dplyr::rename(common_name='Sample')
pdx.info <- pdx.info[pdx.info$common_name %in% chr8q.info$Sample,c("common_name", "Sex", "PRC2_Status")]
colnames(pdx.info)[1] <- "Sample"

pdx.info2 <- merge(chr8q.info, pdx.info)
rownames(pdx.info2) <- pdx.info2$Sample
colnames(pdx.info2) <- c("Sample", "Median Chr8q Copy Number", "Sex", "PRC2 Status")

pdx.overlap <- pdx.info2[pdx.info2$Sample %in% global.meta.df2$Sample,]

global.meta.df3 <- distinct(global.meta.df2[,c("Sample","SampleName")])
global.meta.df3 <- merge(global.meta.df3, pdx.info2)
rownames(global.meta.df3) <- global.meta.df3$SampleName
global.meta.df3$PDX <- global.meta.df3$Sample
global.meta.df3$Sample <- global.meta.df3$SampleName

msigdb.info <- msigdbr::msigdbr()
msigdb.genes <- unique(msigdb.info$gene_symbol)
cnv.red <- cnv[cnv$Gene %in% msigdb.genes,]
pdxRNA.red <- pdxRNA[pdxRNA$Gene %in% msigdb.genes,]

omics <- list("Copy_number" = cnv, 
              "Proteomics" = list("Global" = global.df, "Phospho" = phospho.df),
              "RNA-Seq" = list("PDX" = pdxRNA,"Tumor" = tumorRNA))
prot.feat <- c("Gene", "SUB_SITE")
rna.feat <- c("Gene", "Gene")
meta.list <- list("Copy_number" = pdx.info2,
                  "Proteomics" = global.meta.df3,
                  "RNA-Seq" = pdx.info2,
                  "RNA-Seq (Overlapping with Proteomics)" = pdx.overlap)
expr.list <- list("Copy_number" = "adherent CCLE",
                  "Proteomics" = list("CCLE proteomics"),
                  "RNA-Seq" = list("adherent CCLE", "adherent CCLE"),
                  "RNA-Seq (Overlapping with Proteomics)" = list("adherent CCLE", "adherent CCLE"))
feature.list <- list("Copy_number" = "Gene",
                     "Proteomics" = prot.feat,
                     "RNA-Seq" = rna.feat,
                     "RNA-Seq (Overlapping with Proteomics)" = rna.feat)
my.syn <- "syn63138110"
setwd(base.path)
setwd("Chr8_quant")
gmt1 <- get_gmt1_v2()
gmt2 <- get_gmt2()
synapser::synLogin()

#### Figure 1 ####
# from Spearman folder
path.map <- list("DNA" = "Copy_Number_overlapping",
                 "RNA" = "RNA-Seq_overlapping",
                 "Protein" = "Proteomics/Proteomics",
                 "Phospho" = "Proteomics/Phospho")

### number of diffexp features ###
## compile diffexp features across omics
all.degs <- data.frame()
all.deg.list <- list()
up.deg.list <- list()
dn.deg.list <- list()
for (i in 1:length(pathway.map)) {
  setwd(base.path, "Chr8q_quant/Spearman")
  setwd(path.map[[i]])
  setwd("Differential_expression")
  temp.degs <- read.csv("Differential_expression_results.csv")
  temp.degs$Omics <- names(path.map)[i]
 
  # dot plot of top genes
  temp.degs$Direction <- "Upregulated"
  if (nrow(temp.degs[temp.degs$Log2FC < 0,]) > 0) {
    temp.degs[temp.degs$Log2FC < 0,]$Direction <- "Downregulated"
  }
  temp.degs$Significant <- FALSE
  if (nrow(temp.degs[temp.degs$adj.P.Val <= 0.05,]) > 0) {
    temp.degs[temp.degs$adj.P.Val <= 0.05,]$Significant <- TRUE
  }
  all.degs <- rbind(all.degs, temp.degs)
  all.deg.list[[names(path.map)[i]]] <- temp.degs
  up.deg.list[[names(path.map)[i]]] <- temp.degs[temp.degs$Significant & temp.degs$Direction == "Upregulated",]
  dn.deg.list[[names(path.map)[i]]] <- temp.degs[temp.degs$Significant & temp.degs$Direction == "Downregulated",]
}
setwd(base.path, "Chr8q_quant/Spearman")
dir.create("Compiled_results")
setwd("Compiled_results")
dir.create("Differential_expression")
setwd("Differential_expression")
write.csv(all.degs, "All_differential_expression_results.csv", row.names = FALSE)
filt.degs <- all.degs[all.degs$adj.P.Val <= 0.05, ]
write.csv(filt.degs, "All_differential_expression_results_maxFDR0.05.csv", row.names = FALSE)

bar.plot <- ggplot2::ggplot(all.degs, aes(x=Omics)) + geom_bar(stat="count") + ggplot2::theme_bw()
ggplot2::ggsave("Chr8_numberOfFeatures_barPlot.pdf", bar.plot, width = 7, height = 7, device = "pdf")

all.degs$Significant <- FALSE
if (nrow(filt.degs) > 0) {
  all.degs[all.degs$adj.P.Val <= 0.05,]$Significant <- TRUE
}
bar.plot1 <- ggplot2::ggplot(all.degs, aes(x=Omics, group = Significant)) + geom_bar(stat="count") + ggplot2::theme_bw()
ggplot2::ggsave("Chr8_numberOfFeatures_sigGroup_barPlot.pdf", bar.plot1, width = 7, height = 7, device = "pdf")

filt.degs <- temp.degs[temp.degs$adj.P.Val <= 0.05, ]
bar.plot2 <- ggplot2::ggplot(filt.degs, aes(x=Omics, group = Direction)) + geom_bar(stat="count") + ggplot2::theme_bw()
ggplot2::ggsave("Chr8_numberOfSigFeatures_directionGroup_barPlot.pdf", bar.plot2, width = 7, height = 7, device = "pdf")

all.degs.compiled <- panSEA::compile_mDEG(all.deg.list)
up.degs.compiled <- panSEA::compile_mDEG(up.deg.list)
dn.degs.compiled <- panSEA::compile_mDEG(dn.deg.list)

combo.files <- list()
all.compiled.deg.results <- list("All_features" = all.degs.compiled,
                                 "Significantly_upregulated_features" = up.degs.compiled,
                                 "Significantly_downregulated_features" = dn.degs.compiled)
for (i in 1:length(all.compiled.deg.results)) {
  compiled.DEG.files <- list("Differential_expression_results.csv" = all.compiled.deg.results[[i]]$results,
                             "Compiled_differential_expression_results.csv" = all.compiled.deg.results[[i]]$mean.results,
                             "Differential_expression_venn_diagram.pdf" = all.compiled.deg.results[[i]]$venn.diagram,
                             "Differential_expression_dot_plot.pdf" = all.compiled.deg.results[[i]]$dot.plot,
                             "Differential_expression_correlations.csv" = all.compiled.deg.results[[i]]$corr,
                             "Differential_expression_correlation_matrix.pdf" = all.compiled.deg.results[[i]]$corr.matrix) 
  combo.files[[names(all.compiled.deg.results)[i]]] <- compiled.DEG.files
}
save_to_synapse(combo.files)

#### number of enriched pathways ####
## compile diffexp features across omics
all.gsea <- data.frame()
gs.names <- c("BioCarta", "GO_BP", "GO_CC", "GO_MF", "Hallmark", "KEGG", "MIR_MIRDB", "Oncogenic_signatures", "PID", "Positional", "Positional_custom", "Reactome","WikiPathways")
gs.names <- c("TFT_GTRD", "MIR_MIRDB", "GO_BP", "GO_CC", "GO_MF", 
              "Oncogenic_signatures", "BioCarta", "KEGG", "PID", "Reactome", 
              "Hallmark", "Positional","WikiPathways")
for (i in 1:length(pathway.map)) {
  omics.gsea <- list()
  for (j in 1:length(gs.names)) {
    setwd(base.path, "Chr8q_quant/Spearman")
    setwd(path.map[[i]])
    setwd("GSEA")
    if (file.exists(gs.names[j])) {
      setwd(gs.names[j])
      temp.gsea <- read.csv("GSEA_results.csv")
      temp.gsea$Omics <- names(path.map)[i]
      temp.gsea$Collection <- gs.names[j]
      
      temp.gsea$Direction <- "Upregulated"
      if (nrow(temp.gsea[temp.gsea$NES < 0,]) > 0) {
        temp.gsea[temp.gsea$NES < 0,]$Direction <- "Downregulated"
      }
      temp.gsea$Significant <- FALSE
      if (nrow(temp.gsea[temp.gsea$FDR_q_value <= 0.25 & temp.gsea$p_value <= 0.05,]) > 0) {
        temp.gsea[temp.gsea$FDR_q_value <= 0.25 & temp.gsea$p_value <= 0.05,]$Significant <- TRUE
      }
      all.gsea <- rbind(all.gsea, temp.gsea)
      omics.gsea[[gs.names[j]]] <- temp.gsea
    }
  }
  setwd(base.path, "Chr8q_quant/Spearman")
  setwd(path.map[[i]])
  setwd("GSEA")
  omics.gsea.compiled <- panSEA::compile_mgsea(omics.gsea)
  omics.gsea.files <- list("GSEA_results.csv" = omics.gsea.compiled$results,
                           "Compiled_GSEA_results.csv" = omics.gsea.compiled$mean.results,
                           "GSEA_venn_diagram.pdf" = omics.gsea.compiled$venn.diagram,
                           "GSEA_dot_plot.pdf" = omics.gsea.compiled$dot.plot,
                           "GSEA_correlations.csv" = omics.gsea.compiled$corr,
                           "GSEA_correlation_matrix.pdf" = omics.gsea.compiled$corr.matrix)
  save_to_synapse(omics.gsea.files)
}
setwd(base.path, "Chr8q_quant/Spearman")
dir.create("Compiled_results")
setwd("Compiled_results")
dir.create("GSEA")
setwd("GSEA")
write.csv(all.gsea, "All_differential_expression_results.csv", row.names = FALSE)
filt.gsea <- all.gsea[all.gsea$Significant, ]
write.csv(filt.gsea, "All_differential_expression_results_maxFDR0.05.csv", row.names = FALSE)

bar.plot <- ggplot2::ggplot(all.gsea, aes(x=Omics)) + geom_bar(stat="count") + ggplot2::theme_bw()
ggplot2::ggsave("Chr8_numberOfFeatures_barPlot.pdf", bar.plot, width = 7, height = 7, device = "pdf")

all.gsea$Significant <- FALSE
if (nrow(filt.gsea) > 0) {
  all.gsea[all.gsea$adj.P.Val <= 0.05,]$Significant <- TRUE
}
bar.plot1 <- ggplot2::ggplot(all.gsea, aes(x=Omics, group = Significant)) + geom_bar(stat="count") + ggplot2::theme_bw()
ggplot2::ggsave("Chr8_numberOfFeatures_sigGroup_barPlot.pdf", bar.plot1, width = 7, height = 7, device = "pdf")

bar.plot2 <- ggplot2::ggplot(filt.gsea, aes(x=Omics, group = Direction)) + geom_bar(stat="count") + ggplot2::theme_bw()
ggplot2::ggsave("Chr8_numberOfFeatures_directionGroup_barPlot.pdf", bar.plot2, width = 7, height = 7, device = "pdf")

# organize results by gene set collection
for (i in 1:length(gs.names)) {
  if (gs.names[i] %in% all.gsea$Collection) {
    setwd(base.path, "Chr8q_quant/Spearman")
    setwd("Compiled_results")
    setwd("GSEA")
    dir.create(gs.names[i])
    setwd(gs.names[i])
    temp.gsea <- all.gsea[all.gsea$Collection == gs.names[i],]
    temp.omics <- unique(temp.gsea$Omics)
    all.gsea.list <- list()
    up.gsea.list <- list()
    dn.gsea.list <- list()
    for (j in 1:length(temp.omics)) {
      all.gsea.list[[temp.omics[j]]] <- temp.gsea[temp.gsea$Omics == temp.omics[j],]
      up.gsea.list[[temp.omics[j]]] <- temp.gsea[temp.gsea$Omics == temp.omics[j] & 
                                                   temp.gsea$Significant & 
                                                   temp.gsea$Direction == "Upregulated",]
      dn.gsea.list[[temp.omics[j]]] <- temp.gsea[temp.gsea$Omics == temp.omics[j] & 
                                                   temp.gsea$Significant & 
                                                   temp.gsea$Direction == "Downregulated",]
    }
    all.gsea.compiled <- panSEA::compile_mgsea(all.gsea.list)
    up.gsea.compiled <- panSEA::compile_mgsea(up.gsea.list)
    dn.gsea.compiled <- panSEA::compile_mgsea(dn.gsea.list)
    
    combo.files <- list()
    all.compiled.gsea.results <- list("All_features" = all.gsea.compiled,
                                      "Upregulated_features" = up.gsea.compiled,
                                      "Downregulated_features" = dn.gsea.compiled)
    for (i in 1:length(all.compiled.gsea.results)) {
      compiled.gsea.files <- list("GSEA_results.csv" = all.compiled.gsea.results[[i]]$results,
                                  "Compiled_GSEA_results.csv" = all.compiled.gsea.results[[i]]$mean.results,
                                  "GSEA_venn_diagram.pdf" = all.compiled.gsea.results[[i]]$venn.diagram,
                                  "GSEA_dot_plot.pdf" = all.compiled.gsea.results[[i]]$dot.plot,
                                  "GSEA_correlations.csv" = all.compiled.gsea.results[[i]]$corr,
                                  "GSEA_correlation_matrix.pdf" = all.compiled.gsea.results[[i]]$corr.matrix) 
      combo.files[[names(all.compiled.gsea.results)[i]]] <- compiled.gsea.files
    }
    save_to_synapse(combo.files)
  }
}

## actually need to re-do GSEA for RNA, copy number for just samples in proteomics
## also compile results across all gene sets to extract global top

#### number of enriched DMEA pathways ####
### bar plot, splitting up/down
### venn diagram of pathways
### dot plot of top pathways
## compile results across omics
## also just global proteomics
### dot plot of top PLKi
## compile results across omics
## also just global proteomics

#### Figure 2 ####
### heatmap of central upregulated nodes across omics ###

### dot plot of MYC enrichment across gene sets, omics ###

### dot plot of leading enriched MYC pathway genes across omics ###

### heatmap of phospho leading enriched pathway genes across samples ###

#### Figure 3 ####
### heatmap of central downregulated nodes across omics ###

### dot plot of WNT enrichment across gene sets, omics? ###

### dot plot of leading enriched WNT pathway genes across omics? ###

### heatmap of phospho leading enriched WNT pathway genes across samples? ###

#### FAK figure ####
### heatmap of central nodes across omics ###

### dot plot of WNT enrichment across gene sets, omics? ###

### dot plot of leading enriched MYC pathway genes across omics ###

### heatmap of phospho leading enriched pathway genes across samples ###


