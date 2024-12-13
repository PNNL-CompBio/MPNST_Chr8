# chr8 MPNST: figure 1
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(biomaRt);library(RIdeogram);library(viridis)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_1")

synapser::synLogin()
synID <- synapser::synStore(synapser::Folder("Figure_1", parent = "syn64367769"))
#### 1. Chr8 median copy number bar plot ####
# load data & format
#med.chr8q <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/Github copy/Chr8/analysis/Chr8_quant/Single_sample_positional_medians_allSamples/PDX Copy Number/PDX Copy Number_Chr8q_median.csv")
med.chr8q <- read.csv(synapser::synGet("syn61811211")$path)
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
                position=position_dodge(0.9)) + xlab("MPNST PDX") +
  scale_fill_manual(values = c("black", "grey"))
ggplot2::ggsave(paste0("PDX Copy Number_Chr8q_median_", Sys.Date(), ".pdf"), width = 7, height = 7)
synapser::synStore(synapser::File(paste0("PDX Copy Number_Chr8q_median_", Sys.Date(), ".pdf"), parent = synID))
#### 2. Chr8 enrichment along chromosome ####
# load enrichment results
#chr8.enr <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/proteomics/analysis/Chr8_quant/Spearman/Copy_number/Copy_number/GSEA/Positional/GSEA_results.csv")
chr8.enr <- read.csv(synapser::synGet("syn63434853")$path)

# load segment info for gene symbols
chr8.msigdb <- msigdbr::msigdbr(species = "Homo sapiens", category = "C1")
chr8.msigdb <- chr8.msigdb[grepl("chr8", chr8.msigdb$gs_name),]

# get chr8 segment positions
segments <- unique(chr8.enr[chr8.enr$FDR_q_value <= 0.25 & chr8.enr$p_value <= 0.05,]$Feature_set)
segments <- segments[grepl("chr8",segments)]

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
#synapser::synStore(synapser::File("Chr8_segment_enrichment.csv", parent = synID))

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
#synapser::synStore(synapser::File("chromosome.svg", parent = synID))
#synapser::synStore(synapser::File("chromosome.pdf", parent = synID))

#### 3. bar plot of # of diffexp Features ####
#source("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/compile_mCorr.R")
path.map <- list("RNA" = "syn63602832",
                 "Protein" = "syn63435101",
                 "Phospho" = "syn63425909")
compil.path <- paste0("Compiled_results_allSamples_",Sys.Date())

all.degs <- data.frame()
all.deg.list <- list()
sig.deg.list <- list()
up.deg.list <- list()
dn.deg.list <- list()
for (i in 1:length(path.map)) {
  temp.degs <- read.csv(synapser::synGet(path.map[[i]])$path)
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
  sig.deg.list[[names(path.map)[i]]] <- temp.degs[temp.degs$Significant,]
  up.deg.list[[names(path.map)[i]]] <- temp.degs[temp.degs$Significance == "Upregulated",]
  dn.deg.list[[names(path.map)[i]]] <- temp.degs[temp.degs$Significance == "Downregulated",]
}

setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_1")
# 
# all.degs$Significance <- factor(all.degs$Significance, levels=c("Not Significant", "Upregulated", "Downregulated"))
# all.degs$Direction <- factor(all.degs$Direction, levels=c("Upregulated", "Downregulated"))
# bar.plot3log <- ggplot2::ggplot(all.degs[all.degs$Significant,], aes(fill = Direction, x=forcats::fct_infreq(Omics))) + 
#   geom_bar(position = "dodge", stat="count") + ggplot2::theme_classic() + 
#   scale_y_continuous(trans='log10') + ggplot2::xlab("Omics Type") +
#   scale_fill_manual(values=c("red", "blue")) + 
#   ggplot2::ylab("Number of Differentially Expressed Features") + coord_flip()
# ggplot2::ggsave("Chr8_numberOfFeatures_directionGroup_barPlot_logScale.pdf", bar.plot3log, width = 7, height = 7, device = "pdf")
# 
# all.degs$Direction <- factor(all.degs$Direction, levels=c("Upregulated", "Downregulated"))
# all.degs$Omics <- factor(all.degs$Omics, levels=c("RNA", "Protein", "Phospho"))
# bar.plot3log <- ggplot2::ggplot(all.degs[all.degs$Significant,], aes(fill = Direction, x=Omics)) + 
#   geom_bar(position = "dodge", stat="count") + ggplot2::theme_classic() + 
#   scale_y_continuous(trans='log10') + ggplot2::xlab("Omics Type") +
#   scale_fill_manual(values=c("red", "blue")) + 
#   ggplot2::ylab("Number of Differentially Expressed Features") + coord_flip()
# ggplot2::ggsave("Chr8_numberOfFeatures_directionGroup_barPlot_logScale_manualOrder.pdf", bar.plot3log, width = 7, height = 7, device = "pdf")

all.degs$Omics <- factor(all.degs$Omics, levels=c("Phospho", "Protein", "RNA"))
all.degs$Direction <- factor(all.degs$Direction, levels=c("Upregulated", "Downregulated"))
bar.plot3log <- ggplot2::ggplot(all.degs[all.degs$Significant,], aes(fill = Direction, x=Omics)) + 
  geom_bar(position = "dodge", stat="count") + ggplot2::theme_classic() + 
  scale_y_continuous(trans='log10') + ggplot2::xlab("Omics Type") +
  scale_fill_manual(values=c("red", "blue")) + 
  ggplot2::ylab("Number of Differentially Expressed Features") + coord_flip() #+
  #theme(axis.text.y = element_text(angle = 45, vjust=0.75, hjust=0.75))
ggplot2::ggsave("Chr8_numberOfFeatures_directionGroup_barPlot_logScale_manualOrderV2.pdf", bar.plot3log, width = 7, height = 7, device = "pdf")

synapser::synStore(synapser::File("Chr8_numberOfFeatures_directionGroup_barPlot_logScale_manualOrderV2.pdf", parent = synID))

bar.plot3log <- ggplot2::ggplot(all.degs[all.degs$Significant,], aes(fill = Direction, x=Omics)) + 
  geom_bar(position = "dodge", stat="count") + ggplot2::theme_classic() + 
  ggplot2::xlab("Omics Type") + scale_fill_manual(values=c("red", "blue")) + 
  ggplot2::ylab("Number of Differentially Expressed Features") + coord_flip() #+
  #theme(axis.text.y = element_text(angle = 45, vjust=0.75, hjust=0.75))
ggplot2::ggsave("Chr8_numberOfFeatures_directionGroup_barPlot_manualOrderV2.pdf", bar.plot3log, width = 7, height = 7, device = "pdf")
synapser::synStore(synapser::File("Chr8_numberOfFeatures_directionGroup_barPlot_manualOrderV2.pdf", parent = synID))

#### 4. venn diagram of diffexp Features ####
# all.degs.compiled <- compile_mCorr(all.deg.list)
# sig.degs.compiled <- compile_mCorr(sig.deg.list)
# up.degs.compiled <- compile_mCorr(up.deg.list)
# dn.degs.compiled <- compile_mCorr(dn.deg.list)
# 
# combo.files <- list()
# all.compiled.deg.results <- list("All_Features" = all.degs.compiled,
#                                  "Significantly_correlated_Features" = sig.degs.compiled,
#                                  "Significantly_upregulated_Features" = up.degs.compiled,
#                                  "Significantly_downregulated_Features" = dn.degs.compiled)
# for (i in 1:length(all.compiled.deg.results)) {
#   compiled.DEG.files <- list("Differential_expression_results.csv" = all.compiled.deg.results[[i]]$results,
#                              "Compiled_differential_expression_results.csv" = all.compiled.deg.results[[i]]$mean.results,
#                              "Differential_expression_venn_diagram.pdf" = all.compiled.deg.results[[i]]$venn.diagram,
#                              "Differential_expression_dot_plot.pdf" = all.compiled.deg.results[[i]]$dot.plot,
#                              "Differential_expression_correlations.csv" = all.compiled.deg.results[[i]]$corr,
#                              "Differential_expression_correlation_matrix.pdf" = all.compiled.deg.results[[i]]$corr.matrix) 
#   combo.files[[names(all.compiled.deg.results)[i]]] <- compiled.DEG.files
# }
# synapser::synLogin()
# compilFolder <- synapser::synStore(synapser::Folder(compil.path, parent = "syn63138110"))
# diffexpFolder <- synapser::synStore(synapser::Folder("Differential_expression", parent = compilFolder))
# save_to_synapse(combo.files, resultsFolder = diffexpFolder, dot.scale=1, width=11)

venn.list <- list()
for (i in 1:length(sig.deg.list)) {
  temp.name <- names(sig.deg.list)[i]
  if (unique(sig.deg.list[[i]]$Feature_type) == "SUB_SITE") {
    venn.list[[temp.name]] <- unique(sub("-.*", "", sig.deg.list[[i]]$Feature))
  } else if (unique(sig.deg.list[[i]]$Feature_type) == "Gene") {
    venn.list[[temp.name]] <- unique(sig.deg.list[[i]]$Feature)
  } else {
    warning("No Gene or SUB_SITE column found")
  }
}
ggvenn::ggvenn(venn.list)
ggplot2::ggsave(paste0("Chr8_diffExp_vennDiagram_", Sys.Date(), ".pdf"), width=7, height=7)

ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_alpha = 0.25, text_size = 10, set_name_size = 8)
ggplot2::ggsave(paste0("Chr8_diffExp_vennDiagram_noPercentages_", Sys.Date(), ".pdf"), width=7, height=7)

ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_color = c("#F8766D", "#00BFC4", "#7CAE00"), 
               fill_alpha = 0.5, text_size = 10, set_name_size = 8)
ggplot2::ggsave(paste0("Chr8_diffExp_vennDiagram_noPercentages_", Sys.Date(), ".pdf"), width=7, height=7)
synapser::synStore(synapser::File(paste0("Chr8_diffExp_vennDiagram_", Sys.Date(), ".pdf"), parent = synID))

#### 5. top up/dn diffexp Features ####
# #all.degs.compiled <- compile_mCorr(all.deg.list, n.dot.Features = 100)
# createDotPlot <- function (DEGs, estCol = "Spearman.est", pCol = "Spearman.p",
#                            qCol = "Spearman.q", n.dot.Features = 100, p=0.05, FDR.Features=0.05,
#                            fname=paste0("Chr8_dotPlot_top_", n.dot.Features,"_Features.pdf")) {
#   DEG.df <- as.data.frame(data.table::rbindlist(DEGs, use.names = TRUE, idcol = "type"))
#   DEG.df$minusLogP <- -log(DEG.df[,pCol], base = 10)
#   DEG.df$minusLogFDR <- -log(DEG.df[,qCol], base = 10)
#   DEG.df$sig <- FALSE
#   DEG.df[!is.na(DEG.df[,pCol]) & !is.na(DEG.df[,qCol]) & 
#            DEG.df[,pCol] < p & DEG.df[,qCol] < FDR.Features, ]$sig <- TRUE
#   
#   # summarize results for each Feature
#   DEG.df$rankVal <- DEG.df[,estCol]
#   DEG.df$P.Value <- DEG.df[,pCol]
#   mean.DEG.df <- plyr::ddply(na.omit(DEG.df), .(Feature), summarize,
#                              mean_rankVal = mean(rankVal),
#                              sd_rankVal = sd(rankVal),
#                              Fisher_p = metap::sumlog(na.omit(P.Value))$p,
#                              types = paste0(type, collapse = ", "),
#                              N_types = length(unique(type)),
#                              N_sig = length(sig[sig]),
#                              sig_types = paste0(type[sig], collapse = ", "))
#   if (length(unique(mean.DEG.df$Fisher_p)) > 1) {
#     mean.DEG.df$adj_Fisher_p <- 
#       qvalue::qvalue(mean.DEG.df$Fisher_p, pi0=1)$qvalues
#   } else {
#     mean.DEG.df$adj_Fisher_p <- NA
#   }
#   
#   # order results by rankVal for dot plot
#   mean.DEG.df <- dplyr::arrange(mean.DEG.df, desc(mean_rankVal))
#   
#   # reduce plot data down to top results
#   sig.DEG.df <- mean.DEG.df[mean.DEG.df$N_sig > 0, ]
#   
#   if (nrow(sig.DEG.df) >= n.dot.Features) {
#     top.DEG.df <- 
#       sig.DEG.df %>% dplyr::slice_max(abs(mean_rankVal), n = n.dot.Features)
#   } else {
#     top.DEG.df <- 
#       mean.DEG.df %>% dplyr::slice_max(abs(mean_rankVal), n = n.dot.Features)
#   }
#   dot.df <- na.omit(DEG.df[DEG.df$Feature %in% top.DEG.df$Feature, ])
#   dot.df$qVal <- dot.df[,qCol]
#   if (any(dot.df$qVal == 0)) {
#     dot.df[dot.df$qVal == 0,]$qVal <- 0.0001
#   }
#   
#   dot.plot <- ggplot2::ggplot(
#     na.omit(dot.df),
#     ggplot2::aes(
#       x = type, y = Feature, color = rankVal,
#       size = -log10(qVal)
#     )
#   ) +
#     ggplot2::geom_point() +
#     ggplot2::scale_y_discrete(limits = mean.DEG.df[
#       mean.DEG.df$Feature %in% top.DEG.df$Feature, ]$Feature) +
#     viridis::scale_color_viridis() +
#     theme_classic() +
#     ggplot2::labs(
#       x = "Input Data",
#       y = "Feature Set",
#       color = estCol, size = "-log(FDR)"
#     )
#   ggplot2::ggsave(fname, dot.plot, width=7, height=7)
#   return(list(df = dot.df, mean.df = mean.DEG.df, plot = dot.plot))
# }
# dotPlot <- createDotPlot(all.deg.list)
# saveRDS(dotPlot, "dotPlot.rds")
# dotPlot <- readRDS("dotPlot.rds")

# try making dot plot of Features sig in 2+ omics types
dot.df <- read.csv(synapser::synGet("syn64024301")$path)
dot.df <- dot.df[grepl("RNA, Protein", dot.df$type),]
dot.df <- dot.df[dot.df$N_sig>1,] # 6 between RNA & protein
dot.df <- dotPlot$df # 789

dot.df <- na.omit(all.degs[all.degs$Significant & all.degs$N >= 9 & all.degs$Feature_type=="Gene",]) # 722; 757 if N>=6; 1622 if don't filter for N
nrow(dot.df[dot.df$Omics=="RNA",]) # 19; 52 if N>=6
nrow(dot.df[dot.df$Omics=="Protein",]) # 703; 705 if N>=6
hist(dot.df$Spearman.est,breaks=20)
min(abs(dot.df$Spearman.est)) # 0.7632938 if N>=6 or N>=9
nrow(dot.df[abs(dot.df$Spearman.est)>=0.9,]) # 117; 151 if N>=6
nrow(dot.df[abs(dot.df$Spearman.est)>=0.95,]) # 14; 31 if N>=6
nrow(dot.df[abs(dot.df$Spearman.est)>=0.92,]) # 58; 91 if N>=6
top.dot.df <- dot.df %>% dplyr::slice_max(abs(Spearman.est), n = 100) # 115; 111 if N>=6
min(abs(top.dot.df$Spearman.est)) # 0.9046445; 0.9187795 if N>=6

topGenes <- unique(dot.df[abs(dot.df$Spearman.est)>=0.9,]$Feature)
mean.dot.df <- read.csv(synapser::synGet("syn64024301")$path)
mean.dot.df <- mean.dot.df[mean.dot.df$feature %in% topGenes,]
geneOrder <- mean.dot.df[order(mean.dot.df$mean_rankVal),]$feature
dot.df <- all.degs[all.degs$Feature %in% topGenes,]
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.plot <- ggplot2::ggplot(
  na.omit(dot.df),
  ggplot2::aes(
    x = Omics, y = Feature, color = Spearman.est,
    size = minusLogFDR
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    #x = "Input Data",
    #y = "featureSet",
    color = "Spearman Correlation", size = "-log(FDR)"
  )
# most are only in prot and not in RNA
mean.dot.df <- read.csv(synapser::synGet("syn64024301")$path)
topGenes <- mean.dot.df[mean.dot.df$N_sig> 0,]$feature # 178 
topGenes <- na.omit(mean.dot.df[(grepl("RNA", mean.dot.df$sig_types) | grepl("Protein", mean.dot.df$sig_types)) &
                                  grepl("RNA, Protein", mean.dot.df$types) & mean.dot.df$adj_Fisher_p <= 0.05,]$feature) # 50
top.mean.df <- mean.dot.df[mean.dot.df$feature %in% topGenes,]
geneOrder <- top.mean.df[order(top.mean.df$mean_rankVal),]$feature
dot.df <- all.degs[all.degs$Feature %in% topGenes,]
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.plot <- ggplot2::ggplot(
  na.omit(dot.df),
  ggplot2::aes(
    x = Omics, y = Feature, color = Spearman.est,
    size = minusLogFDR
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_x_discrete(limits=c("RNA","Protein")) +
  viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Gene",
    color = "Spearman Correlation", size = "-log(FDR)"
  )
ggplot2::ggsave("Correlated_genes_adjFisherP0.05_dotPlot.pdf", dot.plot, width=4, height=9)

# another style
geneOrder <- mean.dot.df[order(mean.dot.df$mean_rankVal) & mean.dot.df$feature %in% topGenes,]$feature
geneFace <- rep("plain", length(geneOrder))
names(geneFace) <- geneOrder
sigInRNA <- mean.dot.df[grepl("RNA", mean.dot.df$sig_types) & mean.dot.df$feature %in% topGenes,]$feature
geneFace[sigInRNA] <- "italic"
sigInBoth <- mean.dot.df[grepl("RNA, Protein", mean.dot.df$sig_types) & mean.dot.df$feature %in% topGenes,]$feature
geneFace[sigInBoth] <- "bold"
# if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style

dot.df <- na.omit(all.degs[all.degs$Feature%in% topGenes &
                        all.degs$Omics %in% c("RNA", "Protein"),])
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
maxLogFDR <- max(dot.df$minusLogFDR)
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = Omics, y = Feature, color = Spearman.est,
    size = minusLogFDR
  )
) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,5)) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="blue",high="red") +
  #viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Gene",
    color = "Spearman Correlation", size = "-log(adjusted p-value)"
  ) + 
  #theme(axis.text.y=element_text(face=geneFace)) +
  geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Correlated_genes_adjFisherP0.05_dotPlot_v2.pdf", dot.plot, width=4, height=9)

