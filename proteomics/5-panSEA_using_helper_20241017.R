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

# get annotations as df
dir.create("Single_sample_positional_medians")
setwd("Single_sample_positional_medians")
dir.create("Single_sample_positional_medians_allSamples")
setwd("Single_sample_positional_medians_allSamples")
omics2 <- list("PDX Copy Number" = pdxCNV,
               "PDX RNA-Seq" = pdxRNA,
               "Tumor RNA-Seq" = tumorRNA)

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
}

#### 3. run panSEA ####
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
              "RNA-Seq" = list("PDX" = pdxRNA,"Tumor" = tumorRNA),
              "RNA-Seq (Overlapping with Proteomics)" = list("PDX" = pdxRNA[,c("Gene",rownames(pdx.overlap))],"Tumor" = tumorRNA[,c("Gene",rownames(pdx.overlap))]))
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
panSEA_corr3(omics, meta.list, feature.list, rank.col = "Median Chr8q Copy Number",
             other.annotations = c("Sex", "PRC2 Status"), expr.list = expr.list, gmt1=gmt1, gmt2=gmt2,
             temp.path = file.path(base.path, "Chr8_quant", "Spearman"), syn.id = my.syn)

omics <- list("Copy_number_overlapping" = cnv[,c("Gene",rownames(pdx.overlap))], 
              "RNA-Seq_overlapping" = pdxRNA[,c("Gene",rownames(pdx.overlap))])
meta.list <- list("Copy_number_overlapping" = pdx.overlap,
                  "RNA-Seq_overlapping" = pdx.overlap)
expr.list <- list("Copy_number_overlapping" = "adherent CCLE",
                  "RNA-Seq_overlapping" = "adherent CCLE")
feature.list <- list("Copy_number_overlapping" = "Gene",
                     "RNA-Seq_overlapping" = "Gene")
my.syn <- "syn63138110"
setwd(base.path)
setwd("Chr8_quant")
gmt1 <- get_gmt1_v2()
gmt2 <- get_gmt2()
synapser::synLogin()
panSEA_corr3(omics, meta.list, feature.list, rank.col = "Median Chr8q Copy Number",
             other.annotations = c("Sex", "PRC2 Status"), expr.list = expr.list, gmt1=gmt1,
             temp.path = file.path(base.path, "Chr8_quant", "Spearman"), syn.id = my.syn)

omics <- list("Copy_number" = cnv.red, 
              "Proteomics" = global.df,
              "RNA-Seq" = pdxRNA)
meta.list <- list("Copy_number" = pdx.info2,
                  "Proteomics" = global.meta.df3,
                  "RNA-Seq" = pdx.info2)
expr.list <- list("Copy_number" = "adherent CCLE",
                  "Proteomics" = "CCLE proteomics",
                  "RNA-Seq" = "adherent CCLE")
feature.list <- list("Copy_number" = "Gene",
                     "Proteomics" = "Gene",
                     "RNA-Seq" = "Gene")
my.syn <- "syn63138110"
setwd(base.path)
setwd("Chr8_quant")
gmt1 <- get_gmt1_v2()
#gmt1[[11]] <- NULL
#gmt[[13]] <- NULL
# names(gmt1) <- c("TFT_GTRD", "MIR_MIRDB", "GO_BP", "GO_CC", "GO_MF", 
#                  "Oncogenic_signatures", "BioCarta", "KEGG", "PID", "Reactome", 
#                  "Hallmark", "Positional")
# also make sure to include WikiPathways
# considering excluding Reactome (1736 gene sets), GO (1-7k+), and MIR_MIRDB (2377 gene sets) for now

synapser::synLogin()
panSEA_corr3(omics, meta.list, feature.list, rank.col = "Median Chr8q Copy Number",
             other.annotations = c("Sex", "PRC2 Status"), expr.list = expr.list, gmt1=gmt1, gmt2=gmt2,
             temp.path = file.path(base.path, "Chr8_quant", "Spearman"), syn.id = my.syn)

omics <- list("Proteomics" = list("Phospho" = phospho.df))
prot.feat <- c("SUB_SITE")
meta.list <- list("Proteomics" = global.meta.df3)
expr.list <- list("Proteomics" = list("CCLE proteomics"))
feature.list <- list("Proteomics" = prot.feat)

panSEA_corr3(omics, meta.list, feature.list, rank.col = "Median Chr8q Copy Number",
             other.annotations = c("PDX", "Sex", "PRC2 Status"), expr.list = expr.list, gmt1=gmt1, gmt2=gmt2,
             temp.path = file.path(base.path, "Chr8_quant", "Spearman"), syn.id = my.syn)

#### PCSF sig v3 ####
# try using top 50 sig prot, phospho, WES - directional this time
# load feature weights
synapser::synLogin()
#rna.allSamples.result <- na.omit(read.csv(synapser::synGet("syn60939044")$path))
#rna.result <- na.omit(read.csv(synapser::synGet("syn61905522")$path)) # restricted to samples also in proteomics
global.result <- na.omit(read.csv(synapser::synGet("syn60219661")$path))
phospho.result <- na.omit(read.csv(synapser::synGet("syn61786190")$path))
gene.result <- read.csv(synapser::synGet("syn61924365")$path) # copy number
gene.allSamples.result <- read.csv(synapser::synGet("syn61920734")$path)

hist(global.result$Pearson.est)
hist(global.result$Spearman.est)
library(ggplot2)
ggplot(global.result, aes(x=Pearson.est, y=Spearman.est))+geom_point()

global.result <- global.result[global.result$Spearman.q <= 0.05, ]
phospho.result <- phospho.result[phospho.result$Spearman.q <= 0.05, ]
gene.result <- gene.result[gene.result$Spearman.q <= 0.05, ]
gene.allSamples.result <- gene.allSamples.result[gene.allSamples.result$Spearman.q <= 0.05, ]

# remove last characters from PNNL subsite names (e.g., change ABC-S325s to ABC-S325)
# don't need to worry about multiple sites because KSDB only has info on individual sites (e.g., ABC-S325 not ABC-S325T291)
phospho.result$SUB_SITE <- substr(phospho.result$SUB_SITE, 1, nchar(phospho.result$SUB_SITE)-1)

# try just looking at pathways of interest now that we noticed YWHAZ from WNT signaling is highly connected
msigdb.info <- msigdbr::msigdbr()
poi <- msigdb.info[grepl("wnt", msigdb.info$gs_name, ignore.case = TRUE) | 
                     grepl("yap", msigdb.info$gs_name, ignore.case = TRUE) |
                     grepl("hippo", msigdb.info$gs_name, ignore.case = TRUE) |
                     grepl("hedgehog", msigdb.info$gs_name, ignore.case = TRUE) |
                     grepl("chr8q", msigdb.info$gs_name, ignore.case = TRUE),]
poi.names <- unique(msigdb.info$gs_name)
poi.genes <- unique(msigdb.info$gene_symbol)

# filter for genes of interest
global.goi <- global.result[global.result$Gene %in% poi.genes,] # 1041
phospho.result$Gene <- stringr::str_split_fixed(phospho.result$SUB_SITE[1], "-", 2)[,1]
phospho.goi <- phospho.result[phospho.result$Gene %in% poi.genes,] # 3267
gene.goi <- gene.result[gene.result$Gene %in% poi.genes,] # 238
gene.allSamples.goi <- gene.allSamples.result[gene.allSamples.result$Gene %in% poi.genes,] # 385

chr8.poi <- poi[grepl("chr8q", poi$gs_name, ignore.case = TRUE),]
chr8.genes <- gene.result[gene.result$Gene %in% unique(chr8.poi$gene_symbol),] # 201
chr8.genes.allSamples <- gene.allSamples.result[gene.allSamples.result$Gene %in% unique(chr8.poi$gene_symbol),] # 354

# identify top 50 features
n.top <- c(5, 10, 25, 50, 100)
shared.inputs <- list("Phospho" = phospho.goi, "Prot" = global.goi, "Gene" = gene.goi)
all.inputs <- list("Phospho" = phospho.goi, "Prot" = global.goi, "Gene" = gene.allSamples.goi)
inputs <- list("sharedSamples_SigChr8qWntYapHippoHedgehogGenes" = shared.inputs,
               "allSamples_SigChr8qWntYapHippoHedgehogGenes" = all.inputs)
#directions <- c("absCorr","positivelyCorr", "negativelyCorr")

shared.inputs <- list("Phospho" = phospho.goi, "Prot" = global.goi, "Gene" = chr8.genes)
all.inputs <- list("Phospho" = phospho.goi, "Prot" = global.goi, "Gene" = chr8.genes.allSamples)
inputs <- list("sharedSamples_SigChr8qWntYapHippoHedgehogGenes_onlyChr8qCNV" = shared.inputs,
               "allSamples_SigChr8qWntYapHippoHedgehogGenes_onlyChr8qCNV" = all.inputs)

shared.inputs <- list("Phospho" = phospho.result, "Prot" = global.result, "Gene" = chr8.genes)
all.inputs <- list("Phospho" = phospho.result, "Prot" = global.result, "Gene" = chr8.genes.allSamples)
inputs <- list("sharedSamples_onlyChr8qCNV" = shared.inputs,
               "allSamples_onlyChr8qCNV" = all.inputs)

shared.inputs <- list("Phospho" = phospho.result, "Prot" = global.result, "Gene" = gene.result)
all.inputs <- list("Phospho" = phospho.result, "Prot" = global.result, "Gene" = gene.allSamples.result)
inputs <- list("sharedSamples" = shared.inputs,
               "allSamples" = all.inputs)

msigdb.info <- msigdbr::msigdbr()
poi2 <- msigdb.info[msigdb.info$gs_name == "PID_MYC_ACTIV_PATHWAY", ]
poi2.names <- unique(msigdb.info$gs_name)
poi2.genes <- unique(msigdb.info$gene_symbol)

# filter for genes of interest
global.goi2 <- global.result[global.result$Gene %in% poi.genes2,] # 1041
phospho.result$Gene <- stringr::str_split_fixed(phospho.result$SUB_SITE[1], "-", 2)[,1]
phospho.goi2 <- phospho.result[phospho.result$Gene %in% poi.genes2,] # 3267
gene.goi2 <- gene.result[gene.result$Gene %in% poi.genes2,] # 238
gene.allSamples.goi2 <- gene.allSamples.result[gene.allSamples.result$Gene %in% poi.genes2,] # 385


