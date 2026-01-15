

#### 3. run panSEA ####
#base.path <- file.path(getwd(),Sys.Date())
source("panSEAFunctions.R")
synLogin()
setwd(paste0(base.path,'/analysis'))
#dir.create("Chr8_quant_20250409")
#setwd("Chr8_quant_20250409")
#cnv.med.chr8q <- read.csv(synapser::synGet("syn66047330")$path)
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

# get median Chr8q copy number for metadata
chr8q.info <- read.csv(file.path(base.path,"analysis/positional_medians/Copy Number/Copy Number_Chr8q_median.csv"))
#chr8q.info <- read.csv(synapser::synGet("syn66047330")$path)

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
saveRDS(cn, "cn.rds")
saveRDS(pdxRNA, "pdxRNA.rds")

#my.syn <- "syn60219614"
setwd(paste0(base.path,'/analysis'))
#dir.create("Chr8_quant_20250409")
#setwd("Chr8_quant_20250409")

setwd('..')
# SG: this broke at some point, i'm not sure how, but now we have to load directly
#
gmt1 <-readRDS(synapser::synGet('syn72245876')$path)
# gmt1 <- get_gmt1_v2(gmt.list1=c("msigdb_Homo sapiens_HS_H","msigdb_Homo sapiens_HS_C1","msigdb_Homo sapiens_HS_C3_TFT:GTRD"),
#                     names1=c("Hallmark","Positional", "TFT_GTRD"))
# gmt2 <- get_gmt2(gmt.list2="ksdb_human", phospho=phospho.df)
# synapser::synLogin()
# 
# # # first, check positional enrichment on copy number
  omics <- list("Copy_number" = cn)
  meta.list <- list("Copy_number" = pdx.info2[pdx.info2$Sample!="JH-2-009",])
  expr.list <- list("Copy_number" = "adherent CCLE")
  feature.list <- list("Copy_number" = "Gene")
  gmt1.cn <- list("Positional" = gmt1$Positional)
  panSEA_corr3(omics, meta.list, feature.list, rank.col = "Median Chr8q Copy Number",
               other.annotations = c("Sex", "PRC2 Status"), expr.list = expr.list, gmt1=gmt1.cn, gmt2= NULL, #gmt2,
               temp.path = file.path(base.path, 'analysis'), syn.id = NULL)#my.syn)
 # 
 # ##how to stop the graphics?
  graphics.off()
  setwd(base.path)
  setwd('..')
  