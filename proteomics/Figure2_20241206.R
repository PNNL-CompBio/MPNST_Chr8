# chr8 MPNST: figure 2
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_2")

synapser::synLogin()
synID <- synapser::synStore(synapser::Folder("Figure_2", parent = "syn64367769"))

#### 1. bar plot of # of enriched gene sets stacked by collection ####
all.gsea <- read.csv(synapser::synGet("syn64025221")$path)
gsea.rna.prot <- all.gsea[all.gsea$type %in% c("RNA", "Protein") &
                            all.gsea$Collection %in% c("Hallmark", "PID", "Oncogenic_signatures"),]
gsea.rna.prot[gsea.rna.prot$Collection=="Oncogenic_signatures",]$Collection <- "Oncogenic"
gsea.rna.prot$Omics <- factor(gsea.rna.prot$Omics, levels=c("RNA", "Protein"))
gsea.rna.prot$Direction <- factor(gsea.rna.prot$Direction, levels=c("Upregulated", "Downregulated"))
bar.plot3log <- ggplot2::ggplot(gsea.rna.prot[gsea.rna.prot$Significant,], aes(fill = Direction, x=Omics)) + 
  geom_bar(position = "dodge", stat="count") + ggplot2::theme_classic() + 
  facet_wrap(~Collection, scales = "free") +
  ggplot2::xlab("Omics Type") +
  scale_fill_manual(values=c("red", "blue")) + 
  ylim(c(0,15)) +
  ggplot2::ylab("Number of Enriched Gene Sets") #+
  #theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
ggplot2::ggsave("Chr8_numberOfGeneSets_directionGroup_barPlot_manualOrderV3.pdf", bar.plot3log, width = 7, height = 7, device = "pdf")

#synapser::synStore(synapser::File("Chr8_numberOfFeatures_directionGroup_barPlot_logScale_manualOrderV2.pdf", parent = synID))

#### 2. venn diagram of diffexp features ####
venn.list <- list()
for (i in 1:length(c("RNA", "Protein"))) {
  temp.name <- c("RNA", "Protein")[i]
  venn.list[[temp.name]] <- unique(gsea.rna.prot[gsea.rna.prot$Omics == temp.name &
                                                   gsea.rna.prot$Significant,]$Feature_set)
}
ggvenn::ggvenn(venn.list)
ggplot2::ggsave(paste0("Chr8_GSEA_vennDiagram_", Sys.Date(), ".pdf"), width=7, height=7)

ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_alpha = 0.25, text_size = 14, set_name_size = 7.5)
ggplot2::ggsave(paste0("Chr8_GSEA_vennDiagram_noPercentages_", Sys.Date(), ".pdf"), width=7, height=7)

ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_alpha = 0.25, text_size = 10, set_name_size = 7.5)
ggplot2::ggsave(paste0("Chr8_GSEA_vennDiagram_noPercentages_", Sys.Date(), "_v2.pdf"), width=7, height=7)

#### 3. top up/dn diffexp gene sets ####
mean.gsea <- read.csv(synapser::synGet("syn64025222")$path)

dot.df <- na.omit(gsea.rna.prot[gsea.rna.prot$Significant,]) # 48
nrow(dot.df[dot.df$Omics=="RNA",]) # 9
nrow(dot.df[dot.df$Omics=="Protein",]) # 39
hist(dot.df$NES,breaks=20)
min(abs(dot.df$NES)) # 1.36
topGenes <- unique(dot.df$Feature_set)

mean.gsea <- mean.gsea[mean.gsea$Feature_set %in% topGenes,]
geneOrder <- mean.gsea[order(mean.gsea$mean_NES),]$Feature_set
geneFace <- rep("plain", length(geneOrder))
names(geneFace) <- geneOrder
sigInRNA <- mean.gsea[grepl("RNA", mean.gsea$sig_types),]$Feature_set
geneFace[sigInRNA] <- "italic"
sigInBoth <- mean.gsea[grepl("RNA, Protein", mean.gsea$sig_types),]$Feature_set
geneFace[sigInBoth] <- "bold"
# if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style

dot.df <- all.gsea[all.gsea$Feature_set %in% topGenes &
                     all.gsea$Omics %in% c("RNA", "Protein"),]
if (any(dot.df$FDR_q_value == 0)) {
  dot.df[dot.df$FDR_q_value == 0,]$FDR_q_value <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"p_value"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"FDR_q_value"], base = 10)
dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
dot.df[dot.df$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
dot.plot <- ggplot2::ggplot(
  na.omit(dot.df),
  ggplot2::aes(
    x = Omics, y = Feature_set, color = NES,
    size = minusLogFDR
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="blue",high="red") +
  #viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Gene Set",
    color = "NES", size = "-log(FDR)"
  ) + 
  theme(axis.text.y=element_text(face=geneFace)) +
  facet_wrap(vars(cols = Collection))
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Enriched_geneSets_dotPlot_facetCollection.pdf", dot.plot, width=9, height=9)

# maybe do for each collection and then combine
collections <- c("Hallmark", "Oncogenic", "PID")
mean.gsea[mean.gsea$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
gsea.rna.prot[gsea.rna.prot$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
maxAbsNES <- max(abs(gsea.rna.prot[gsea.rna.prot$Collection %in% collections,]$NES))
maxLogFDR <- max(gsea.rna.prot[gsea.rna.prot$Collection %in% collections,]$minusLogFDR)
library(patchwork); library(ggplot2)
for (i in collections) {
  nGenes <- length(unique(na.omit(gsea.rna.prot[gsea.rna.prot$Collection == i,])$Feature_set))
  topGenes <- unique(na.omit(gsea.rna.prot[gsea.rna.prot$Significant &
                                             gsea.rna.prot$Collection == i,])$Feature_set) # 48
  nTopGenes <- length(topGenes)
  temp.mean.gsea <- mean.gsea[mean.gsea$Feature_set %in% topGenes,]
  geneOrder <- temp.mean.gsea[order(temp.mean.gsea$mean_NES),]$Feature_set
  geneFace <- rep("plain", length(geneOrder))
  names(geneFace) <- geneOrder
  sigInRNA <- temp.mean.gsea[grepl("RNA", temp.mean.gsea$sig_types),]$Feature_set
  geneFace[sigInRNA] <- "italic"
  sigInBoth <- temp.mean.gsea[grepl("RNA, Protein", temp.mean.gsea$sig_types),]$Feature_set
  geneFace[sigInBoth] <- "bold"
  # if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
  
  temp.dot.df <- all.gsea[all.gsea$Feature_set %in% topGenes &
                       all.gsea$Omics %in% c("RNA", "Protein"),]
  if (any(temp.dot.df$FDR_q_value == 0)) {
    temp.dot.df[temp.dot.df$FDR_q_value == 0,]$FDR_q_value <- 0.0001
  }
  temp.dot.df$minusLogP <- -log(temp.dot.df[,"p_value"], base = 10)
  temp.dot.df$minusLogFDR <- -log(temp.dot.df[,"FDR_q_value"], base = 10)
  temp.dot.df$Omics <- factor(temp.dot.df$Omics, levels=c("RNA", "Protein"))
  #temp.dot.df[temp.dot.df$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
  dot.plot <- ggplot2::ggplot(
    na.omit(temp.dot.df),
    ggplot2::aes(
      x = Omics, y = Feature_set, color = NES,
      size = minusLogFDR
    )
  ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = geneOrder) +
    scale_color_gradient2(low="blue",high="red", limits=c(-maxAbsNES, maxAbsNES)) +
    #viridis::scale_color_viridis() +
    theme_classic() +
    ggplot2::labs(
      x = "Omics Type",
      y = "Gene Set",
      color = "NES", size = "-log(FDR)"
    ) + 
    theme(axis.text.y=element_text(face=geneFace))
  plot.annot <- paste0(i, " (", nTopGenes, " / ", nGenes, " gene sets enriched)")
  dot.plot <- dot.plot + ggtitle(plot.annot)
  dot.plot
  if (i == collections[1]) {
    gsea.dot.plots <- (dot.plot + theme(legend.position = "none"))
  } else if (i == collections[2]) {
    gsea.dot.plots <- gsea.dot.plots / dot.plot
  } else {
    gsea.dot.plots <- gsea.dot.plots / (dot.plot + theme(legend.position = "none"))
  }
  # most are only in prot and not in RNA
  
}
gsea.dot.plots
ggplot2::ggsave("Enriched_geneSets_dotPlot_patchworkCollection.pdf", gsea.dot.plots, width=9, height=9)
gsea.dot.plots2 <- gsea.dot.plots + plot_annotation(caption = "Key\nBold: enriched in both RNA & protein\n Italics: enriched in RNA\n Plain: enriched in Protein")
gsea.dot.plots2
ggplot2::ggsave("Enriched_geneSets_dotPlot_patchworkCollection_withKey.pdf", gsea.dot.plots2, width=9, height=9)

#### Look closer at ESC_V6.5_UP_EARLY.V1_DN leading edge genes ####
corr.result <- read.csv(synapser::synGet("syn64024300")$path)
esc.genes <- msigdbr::msigdbr(species="Homo sapiens",category="C6")
esc.genes <- unique(esc.genes[esc.genes$gs_name == "ESC_V6.5_UP_EARLY.V1_DN",]$gene_symbol) # 161
esc.corr <- corr.result[corr.result$feature %in% esc.genes,] # 283
esc.corr <- esc.corr[esc.corr$type %in% c("RNA", "Protein"),] # 124
topGenes <- unique(esc.corr[esc.corr$sig,]$feature) # 11

mean.corr <- plyr::ddply(esc.corr[esc.corr$feature %in% topGenes,], .(feature), summarize,
                         mean_rankVal = mean(Spearman.est, na.rm = TRUE),
                         sig_types = paste(type, collapse = ", "))
geneOrder <- mean.corr[order(mean.corr$mean_rankVal),]$feature
geneFace <- rep("plain", length(geneOrder))
names(geneFace) <- geneOrder
sigInRNA <- mean.corr[grepl("RNA", mean.corr$sig_types),]$feature
geneFace[sigInRNA] <- "italic"
sigInBoth <- mean.corr[grepl("RNA, Protein", mean.corr$sig_types),]$feature
geneFace[sigInBoth] <- "bold"
# if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style

dot.df <- corr.result[corr.result$feature %in% topGenes &
                     corr.result$Omics %in% c("RNA", "Protein"),]
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
maxLogFDR.esc <- max(dot.df$minusLogFDR)
absMaxEst.esc <- max(abs(dot.df$Spearman.est))
esc.dot.plot <- ggplot2::ggplot(
  na.omit(dot.df),
  ggplot2::aes(
    x = Omics, y = feature, color = Spearman.est,
    size = minusLogFDR
  )
) + scale_size(limits=c(0,maxLogFDR.esc), range = c(0.5,5)) +
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
  theme(axis.text.y=element_text(face=geneFace))
esc.dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Correlated_ESC_genes_dotPlot.pdf", esc.dot.plot, width=4, height=4)

#### Look closer at HALLMARK_COAGULATION leading edge genes ####
corr.result <- read.csv(synapser::synGet("syn64024300")$path)
coag.genes <- msigdbr::msigdbr(species="Homo sapiens",category="H")
coag.genes <- unique(coag.genes[coag.genes$gs_name == "HALLMARK_COAGULATION",]$gene_symbol) # 138
coag.corr <- corr.result[corr.result$feature %in% coag.genes,] # 265
coag.corr <- coag.corr[coag.corr$type %in% c("RNA", "Protein"),] # 127
topGenes <- unique(coag.corr[coag.corr$sig,]$feature) # 12

mean.corr <- plyr::ddply(coag.corr[coag.corr$feature %in% topGenes,], .(feature), summarize,
                         mean_rankVal = mean(Spearman.est, na.rm = TRUE),
                         sig_types = paste(type, collapse = ", "))
geneOrder <- mean.corr[order(mean.corr$mean_rankVal),]$feature
geneFace <- rep("plain", length(geneOrder))
names(geneFace) <- geneOrder
sigInRNA <- mean.corr[grepl("RNA", mean.corr$sig_types),]$feature
geneFace[sigInRNA] <- "italic"
sigInBoth <- mean.corr[grepl("RNA, Protein", mean.corr$sig_types),]$feature
geneFace[sigInBoth] <- "bold"
# if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style

dot.df <- corr.result[corr.result$feature %in% topGenes &
                        corr.result$Omics %in% c("RNA", "Protein"),]
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
maxLogFDR.coag <- max(dot.df$minusLogFDR)
absMaxEst.coag <- max(abs(dot.df$Spearman.est))
coag.dot.plot <- ggplot2::ggplot(
  na.omit(dot.df),
  ggplot2::aes(
    x = Omics, y = feature, color = Spearman.est,
    size = minusLogFDR
  )
) + scale_size(limits=c(0,maxLogFDR.coag), range = c(0.5,5)) +
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
  theme(axis.text.y=element_text(face=geneFace))
coag.dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Correlated_coagulation_genes_dotPlot.pdf", coag.dot.plot, width=4, height=4)

#### combine esc and coag dot plots ####
absMaxEst <- max(c(absMaxEst.coag, absMaxEst.esc))
maxLogFDR <- max(c(maxLogFDR.coag, maxLogFDR.esc))
esc.dot.plot2 <- esc.dot.plot + ggtitle("ESC_V6.5_UP_EARLY.V1_DN") + 
  scale_color_gradient2(low="blue",high="red", limits = c(-absMaxEst, absMaxEst)) +
  scale_size(limits=c(0,maxLogFDR), range = c(0.5,5))
coag.dot.plot2 <- coag.dot.plot + ggtitle("HALLMARK_COAGULATION") + 
  scale_color_gradient2(low="blue",high="red", limits = c(-absMaxEst, absMaxEst)) +
  scale_size(limits=c(0,maxLogFDR), range = c(0.5,5))
esc.coag.dot.plots <-  (coag.dot.plot2 + theme(legend.position = "none"))/ esc.dot.plot2
esc.coag.dot.plots <- esc.coag.dot.plots + plot_annotation(caption="Key\nBold: correlated in both RNA & protein\n Italics: correlated in RNA\n Plain: correlated in Protein")
esc.coag.dot.plots
ggplot2::ggsave("Correlated_ESC_and_Hallmark_genes_dotPlot.pdf", esc.coag.dot.plots, width=4, height=9)
