# Differential expression & enrichment analyses: PDX transcriptomics, global & phospho
# Chr8: median copy number vs. protein/phospho expression
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-09-10
remove(list=ls())
# overview:
#### 1. Import metadata & crosstabs
#### 2. Run panSEA: only MPNST PDX with proteomics
#### 3. Run panSEA: all MPNST PDX
#### 4. PCSF network analysis

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr); library(plyr)
library(htmlwidgets); library(webshot); library(scales); library(msigdbr)
library(plyr); library(dplyr); library(R.utils); library(ggplot2)
#webshot::install_phantomjs()
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis"
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
#source("https://github.com/PNNL-CompBio/MPNST_Chr8/blob/main/panSEA_helper_20240913.R")
synapser::synLogin()

#### 1. Import metadata & crosstabs ####
### proteomics
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/")
meta.df <- readxl::read_excel(synapser::synGet('syn65986564')$path)

global.df <- read.table(synapser::synGet("syn65986566")$path, sep = "\t")
phospho.df <- read.table(synapser::synGet("syn65986573")$path, sep = "\t")
phospho.pep.df <- read.table(synapser::synGet("syn65986575")$path, sep = "\t")

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)
phospho.df$SUB_SITE <- rownames(phospho.df)
phospho.pep.df$SUB_SITE <- rownames(phospho.pep.df)

### other omics
manifest<-synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
  dplyr::rename(common_name='Sample')
pdx_data<-manifest|>dplyr::select(common_name,Sex,RNASeq='PDX_RNASeq',
                                  CopyNumber='PDX_CNV')

# from coderdata 00-buildGeneFile.R
# check gene symbols because some don't make sense
if(!require('org.Hs.eg.db')){
  BiocManager::install('org.Hs.eg.db')
  library(org.Hs.eg.db)
}

library(dplyr)
##get entrez ids to symbol
loadRNAandCN <- function(pdx_data) {
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
  
  rnaseq<-do.call('rbind',lapply(setdiff(pdx_data$RNASeq,NA),function(x){
    sample<-base::subset(pdx_data,RNASeq==x)
    res<-data.table::fread(synGet(x)$path)|>
      tidyr::separate(Name,into=c('other_id','vers'),sep='\\.')|>
      left_join(joined.df)|>
      dplyr::select(gene_symbol,TPM)|>
      subset(!is.na(gene_symbol)|>
               subset(TPM!=0))
    res$sample <- sample$common_name
    return(distinct(res))
  }))
  rnaseq <- na.omit(rnaseq)
  rnaseq <- reshape2::dcast(rnaseq, gene_symbol ~ sample, mean, value.var = "TPM")
  colnames(rnaseq)[1] <- "Gene"
  
  cn<-do.call(rbind,lapply(setdiff(pdx_data$CopyNumber,c(NA,"NA")),function(x){
    sample<-subset(pdx_data,CopyNumber==x)
    print(sample$improve_sample_id)
    res<-data.table::fread(synGet(x)$path)
    
    long_df<- res|>
      tidyr::separate_rows(gene,sep=',')|>
      dplyr::rename(gene_symbol='gene')|>
      dplyr::left_join(joined.df)|>
      subset(!is.na(gene_symbol))|>
      dplyr::select(gene_symbol,log2)|>
      dplyr::distinct()|>
      dplyr::mutate(copy_number=2*(2^log2))|>
      dplyr::select(-log2)

    long_df$sample <- sample$common_name
    return(long_df)
  }))
  colnames(cn)[1] <- "Gene"
  return(list(rna = rnaseq, cn = cn))
}
pdxOmics <- loadRNAandCN(pdx_data[pdx_data$common_name != "JH-2-009",]) # remove JH-2-009 since Ava said it is contaminated
pdxRNA <- pdxOmics$rna
pdxCN <- pdxOmics$cn

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
colnames(phospho.pep.df)[1:12] <- prot.names
phospho.pep.df <- phospho.pep.df[,c("SUB_SITE", prot.names)]

rna.meta.df <- dplyr::distinct(global.meta.df2[,c("Sample", annotations)])
rownames(rna.meta.df) <- rna.meta.df$Sample
rownames(global.meta.df2) <- global.meta.df2$SampleName

# look at chr8q amplification
setwd(base.path)
dir.create("Chr8_quant_20250409")
setwd("Chr8_quant_20250409")
dir.create("positional_medians")
setwd("positional_medians")
omics2 <- list("Copy Number" = pdxCN,
               "RNA-Seq" = pdxRNA,
               "Proteomics" = global.df)

msigdb.info <- msigdbr::msigdbr(category="C1")
reduced.msigdb <- dplyr::distinct(msigdb.info[,c("gs_name", "gene_symbol")])
for (i in 1:length(omics2)) {
  setwd(file.path(base.path, "Chr8_quant_20250409", "positional_medians"))
  dir.create(names(omics2)[i])
  setwd(names(omics2)[i])
  pos.GSEA <- list()
  if (names(omics2)[i] == "Copy Number") {
    temp.samples <- unique(omics2[[i]]$sample)

    for (j in temp.samples) { # assumes first column is feature names
      # prep input data
      temp.input <- omics2[[i]][omics2[[i]]$sample == j,c("Gene","copy_number")]
      colnames(temp.input)[1] <- "gene_symbol"
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
  ggplot2::ggsave(paste0(names(omics2)[i], "_Chr8q_positional_median.pdf"), width=5, height=5)
  
  chr8q.df$`Mean Copy Number` <- chr8q.df$mean_copy_number
  ggplot2::ggplot(chr8q.df, aes(fill=Sample, x=`Chr8q Position`, y=`Mean Copy Number`)) + 
    ggplot2::geom_bar(stat="identity", position="dodge") + ggplot2::theme_minimal() +
    geom_errorbar(aes(ymin=`Mean Copy Number` - sd_copy_number, 
                      ymax = `Mean Copy Number` + sd_copy_number), width=0.2,
                  position=position_dodge(0.9))
  ggplot2::ggsave(paste0(names(omics2)[i], "_Chr8q_positional_mean.pdf"), width=5, height=5)
  
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
  ggplot2::ggsave(paste0(names(omics2)[i], "_Chr8q_median.pdf"), width=5, height=5)
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
  ggplot2::ggsave(paste0(names(omics2)[i], "_Chr8q_mean.pdf"), width=5, height=5)
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
setwd("Chr8_quant_20250409")
cnv.med.chr8q <- read.csv("positional_medians/Copy Number/Copy Number_Chr8q_median.csv")
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

# correlate with peptide level phospho
phospho.pep.forCorr <- as.data.frame(t(phospho.pep.df[,prot.names]))
phospho.pep.forCorr$SampleName <- rownames(phospho.pep.forCorr)
reduced.meta <- global.meta.df2[,c("SampleName", "Chr8q_median")]
phospho.pep.corr.input <- merge(reduced.meta, phospho.pep.forCorr)
phospho.pep.corr <- DMEA::rank_corr(phospho.pep.corr.input, plots=FALSE)
colnames(phospho.pep.corr$result)[1] <- "Feature"
write.csv(phospho.pep.corr$result, "Phospho_peptide_correlations_Chr8q.csv", row.names = FALSE)

# get median Chr8q copy number for metadata
chr8q.info <- read.csv(file.path(base.path,"Chr8_quant_20250409/positional_medians/Copy Number/Copy Number_Chr8q_median.csv"))

# get other metadata
colnames(pdx_data)[1] <- "Sample"
chr8q.info <- merge(pdx_data, chr8q.info, all.y = TRUE)
colnames(chr8q.info)[5] <- "Median Chr8q Copy Number"

# get more metadata for heatmaps
chr8q.info <- chr8q.info[,c("Sample", "Sex", "Median Chr8q Copy Number")]

pdx.info <- synapser::synTableQuery("select * from syn53503360")$asDataFrame()|>
  dplyr::rename(common_name='Sample')
pdx.info <- pdx.info[pdx.info$common_name %in% chr8q.info$Sample,c("common_name", "Sex", "PRC2_Status")]
colnames(pdx.info)[1] <- "Sample"

pdx.info2 <- merge(chr8q.info, pdx.info, by=c("Sample","Sex"))
rownames(pdx.info2) <- pdx.info2$Sample
colnames(pdx.info2) <- c("Sample", "Sex", "Median Chr8q Copy Number", "PRC2 Status")
pdx.info2 <- pdx.info2[,c("Sample", "Median Chr8q Copy Number", "Sex", "PRC2 Status")]

global.meta.df3 <- distinct(global.meta.df2[,c("Sample","SampleName")])
global.meta.df3 <- merge(global.meta.df3, pdx.info2)
rownames(global.meta.df3) <- global.meta.df3$SampleName
global.meta.df3$PDX <- global.meta.df3$Sample
global.meta.df3$Sample <- global.meta.df3$SampleName

# format wide with common_name ~ gene_symbol
cn <- reshape2::dcast(pdxCN, Gene ~ sample, mean, value.var = "copy_number")

my.syn <- "syn65988130"
setwd(base.path)
setwd("Chr8_quant_20250409")
#gmt1 <- get_gmt1_v2()
gmt1 <- readRDS("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/proteomics/analysis/gmt1_more.rds")
#gmt2 <- get_gmt2()
gmt2 <- readRDS("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/proteomics/analysis/gmt2.rds")
synapser::synLogin()

# first, check positional enrichment on copy number
omics <- list("Copy_number" = cn)
meta.list <- list("Copy_number" = pdx.info2[pdx.info2$Sample!="JH-2-009",])
expr.list <- list("Copy_number" = "adherent CCLE")
feature.list <- list("Copy_number" = "Gene")
gmt1.cn <- list("Positional" = gmt1$Positional)
panSEA_corr3(omics, meta.list, feature.list, rank.col = "Median Chr8q Copy Number",
             other.annotations = c("Sex", "PRC2 Status"), expr.list = expr.list, gmt1=gmt1.cn, gmt2=gmt2,
             temp.path = file.path(base.path, "Chr8_quant_20250409", "Spearman"), syn.id = my.syn)

# next, proteomics and RNA
setwd(base.path)
setwd("Chr8_quant_20250409")
omics <- list("Proteomics" = list("Global" = global.df, "Phospho" = phospho.df),
              "RNA-Seq" = pdxRNA)
meta.list <- list("Proteomics" = global.meta.df3,
                  "RNA-Seq" = pdx.info2[pdx.info2$Sample!="JH-2-009",])
expr.list <- list("Proteomics" = "CCLE proteomics",
                  "RNA-Seq" = "adherent CCLE")
feature.list <- list("Proteomics" = c("Gene", "SUB_SITE"),
                     "RNA-Seq" = "Gene")
gmt1.rest <- list("Hallmark" = gmt1$Hallmark,
                  "KEGG" = gmt1$KEGG,
                  "Oncogenic" = gmt1$Oncogenic,
                  "PID" = gmt1$PID,
                  "TFT_GTRD" = gmt1$TFT_GTRD)
panSEA_corr3(omics, meta.list, feature.list, rank.col = "Median Chr8q Copy Number",
             other.annotations = c("Sex", "PRC2 Status"), expr.list = expr.list, gmt1=gmt1.rest, gmt2=gmt2,
             temp.path = file.path(base.path, "Chr8_quant_20250409", "Spearman"), syn.id = my.syn)


