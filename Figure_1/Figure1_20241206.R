# chr8 MPNST: figure 1
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(ggplot2)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_1")

# to do - change inputs to pull from Synapse

#### 1. Chr8 median copy number bar plot ####
# load data & format
med.chr8q <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/Github copy/Chr8/analysis/Chr8_quant/Single_sample_positional_medians_allSamples/PDX Copy Number/PDX Copy Number_Chr8q_median.csv")
colnames(med.chr8q)[2] <- "Median Chr8q Copy Number"
med.chr8q$Proteomics <- "Not Collected"
prot.samples <- c("JH-2-002", "JH-2-055", "JH-2-079c", "MN-2", "WU-225", "WU-487")
med.chr8q[med.chr8q$Sample %in% prot.samples,]$Proteomics <- "Collected"

# order lowest to highest
med.chr8q <- med.chr8q[order(med.chr8q$`Median Chr8q Copy Number`),]
sample.order <- med.chr8q$Sample

# create bar plot
ggplot2::ggplot(med.chr8q, aes(x=Sample, y=`Median Chr8q Copy Number`, fill = Proteomics)) + 
  ggplot2::geom_bar(stat="identity", position="dodge") + 
  scale_x_discrete(limits = sample.order) + ggplot2::theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  geom_errorbar(aes(ymin=`Median Chr8q Copy Number` - sd_med_copy_number, 
                    ymax = `Median Chr8q Copy Number` + sd_med_copy_number), width=0.2,
                position=position_dodge(0.9)) + xlab("MPNST PDX")
ggplot2::ggsave(paste0("PDX Copy Number_Chr8q_median_", Sys.Date(), ".pdf"), width = 7, height = 7)

#### 2. Chr8 enrichment along chromosome ####
# load enrichment results
chr8.enr <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/proteomics/analysis/Chr8_quant/Spearman/Copy_number/Copy_number/GSEA/Positional/GSEA_results.csv")

# load segment info for gene symbols
chr8.msigdb <- msigdbr::msigdbr(species = "Homo sapiens", category = "C1")
chr8.msigdb <- chr8.msigdb[grepl("chr8", chr8.msigdb$gs_name),]

# get chr8 segment positions
segments <- unique(chr8.enr[chr8.enr$FDR_q_value <= 0.25 & chr8.enr$p_value <= 0.05,]$Feature_set)
segments <- segments[grepl("chr8",segments)]

library(biomarRt)
ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl", host="https://useast.ensembl.org")
segment.pos <- data.frame(segments)
segment.pos[,c("start", "end")] <- NA
for (i in segments) {
  # get gene ids
  temp.genes <- unique(chr8.msigdb[chr8.msigdb$gs_name == i,]$ensembl_gene)
  
  # find start/end of genes
  if (length(temp.genes) > 0) {
    gene.pos <- biomaRt::getBM(attributes=c("ensembl_gene_id","start_position","end_position"),
                               filters="ensembl_gene_id",values = temp.genes,mart=ensembl) 
    segment.pos[segment.pos$segments == i,]$start <- min(gene.pos$start_position)
    segment.pos[segment.pos$segments == i,]$end <- max(gene.pos$end_position)
  }
}

# add enrichment scores
colnames(chr8.enr)[2] <- "segments"
segment.pos <- merge(segment.pos, chr8.enr[,c("segments","NES")], by="segments")
segment.pos$chr <- 8
write.csv(segment.pos, "Chr8_segment_enrichment.csv", row.names = FALSE)

library(RIdeogram)
# following vignette: https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html
data("human_karyotype", package="RIdeogram")
human8 <- as.data.frame(human_karyotype[8,])
ideo.df <- segment.pos[,c("chr", "start", "end", "NES")]
colnames(ideo.df) <- c("Chr", "Start", "End", "Value")
label.df <- segment.pos[,c("segments", "chr", "start", "end")]
label.df$Type <- label.df$segments
label.df$Shape <- "box"
label.df$color <- "black"
label.df[label.df$segments=="chr8q24",]$color <- "red"
label.df <- label.df[,c("Type", "Shape", "chr", "start", "end", "color")]
colnames(label.df)[3:5] <- c("Chr", "Start", "End")
RIdeogram::ideogram(karyotype = human_karyotype, overlaid = ideo.df, label = label.df, label_type = "marker", colorset1 = c("blue","white", "red"))
convertSVG("chromosome.svg", device="pdf")

#### 3. bar plot of # of diffexp features ####
all.degs <- data.frame()
all.deg.list <- list()
for (i in 1:length(path.map)) {
  setwd(file.path(base.path, "Chr8_quant/Spearman"))
  setwd(path.map[[i]])
  setwd("Differential_expression")
  temp.degs <- read.csv("Differential_expression_results.csv")
  temp.degs$Omics <- names(path.map)[i]
  
  # dot plot of top genes
  temp.degs$Direction <- "Upregulated"
  if (nrow(temp.degs[!is.na(temp.degs$Spearman.est) & 
                     temp.degs$Spearman.est < 0,]) > 0) {
    temp.degs[!is.na(temp.degs$Spearman.est) &
                temp.degs$Spearman.est < 0,]$Direction <- "Downregulated"
  }
  temp.degs$Significant <- FALSE
  if (nrow(temp.degs[!is.na(temp.degs$Spearman.q) &
                     temp.degs$Spearman.q <= 0.05,]) > 0) {
    temp.degs[!is.na(temp.degs$Spearman.q) &
                temp.degs$Spearman.q <= 0.05,]$Significant <- TRUE
  }
  temp.degs$Significance <- "Not Significant"
  if (any(temp.degs$Significant)) {
    temp.degs[temp.degs$Significant,]$Significance <- temp.degs[temp.degs$Significant,]$Direction
  }
  temp.degs$Feature <- temp.degs[,1]
  temp.degs$Feature_type <- colnames(temp.degs)[1]
  all.degs <- rbind(all.degs, temp.degs[,2:ncol(temp.degs)])
  nonFeatureNames <- colnames(temp.degs)[2:ncol(temp.degs)][colnames(temp.degs)[2:ncol(temp.degs)] != "Feature"]
  temp.degs <- temp.degs[,c("Feature", nonFeatureNames)]
  all.deg.list[[names(path.map)[i]]] <- temp.degs
}
setwd(file.path(base.path, "Chr8_quant/Spearman"))
dir.create(compil.path)
setwd(compil.path)
dir.create("Differential_expression")
setwd("Differential_expression")
write.csv(all.degs, "All_differential_expression_results.csv", row.names = FALSE)
filt.degs <- all.degs[!is.na(all.degs$Spearman.q) & all.degs$Spearman.q <= 0.05, ]
write.csv(filt.degs, "All_differential_expression_results_maxFDR0.05.csv", row.names = FALSE)

all.degs$Significance <- factor(all.degs$Significance, levels=c("Not Significant", "Downregulated", "Upregulated"))
bar.plot3log <- ggplot2::ggplot(all.degs, aes(fill = Significance, x=forcats::fct_infreq(Omics))) + 
  geom_bar(position = "dodge", stat="count") + ggplot2::theme_classic() + 
  scale_y_continuous(trans='log10') +
  ggplot2::xlab("Omics Type") +
  ggplot2::ylab("Number of Quantified Features")
ggplot2::ggsave("Chr8_numberOfFeatures_significanceGroup_barPlot_logScale.pdf", bar.plot3log, width = 7, height = 7, device = "pdf")


#### 4. venn diagram of diffexp features ####

#### 5. top up/dn diffexp features ####