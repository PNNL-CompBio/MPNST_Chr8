# chr8 MPNST: figure 3
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_TF")

synapser::synLogin()
synID <- synapser::synStore(synapser::Folder("Figure_3", parent = "syn64367769"))

#### 1. volcano plot of # of GSEA_TFT_GTRD results ####
all.gsea <- read.csv(synapser::synGet("syn64374317")$path)
plot.data <- all.gsea[all.gsea$Omics == "RNA" &
                            all.gsea$Collection == "TFT_GTRD",] # 489
plot.data[plot.data$Collection=="TFT_GTRD",]$Collection <- "Transcription Factor Targets"
FDR <- 0.25
plot.data$Significance <- paste0("FDR > ", FDR)
plot.data[plot.data$FDR_q_value < FDR, ]$Significance <- paste0("FDR < ", FDR)
plot.data$Significance <- factor(plot.data$Significance,
                                 levels = c(paste0("FDR < ", FDR),
                                            paste0("FDR > ", FDR)))
plot.data$TF <- sub("_TARGET_GENES.*","",plot.data$Feature_set)
if (any(plot.data$p_value == 0)) {
  plot.data[plot.data$p_value == 0,]$p_value <- 0.0001
}
limit.x <- ceiling(max(abs(plot.data$NES)))
limit.y <- ceiling(max(-log(plot.data$p_value, base=10)))
volc <- ggplot2::ggplot(data = plot.data, aes(
  x = NES, y = -log(p_value, 10),
  color = Significance
)) +
  ggplot2::geom_point(size = 4) +
  # ggrepel::geom_text_repel(
  #   data = subset(plot.data, Significant),
  #   mapping = aes(label = TF, size = I(6)), nudge_y = 0.25
  # ) +
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
volc

#synapser::synStore(synapser::File("Chr8_numberOfFeatures_directionGroup_barPlot_logScale_manualOrderV2.pdf", parent = synID))

#### 3. top up/dn diffexp gene sets ####
dot.df <- na.omit(plot.data[plot.data$Significant,]) # 12
topGenes <- unique(dot.df$Feature_set)
geneOrder <- dot.df[order(dot.df$NES),]$TF
dot.df$minusLogP <- -log(dot.df[,"p_value"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"FDR_q_value"], base = 10)
dot.plot <- ggplot2::ggplot(
  na.omit(dot.df),
  ggplot2::aes(
    x = Omics, y = TF, color = NES,
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
    y = "Transcription Factor",
    color = "NES", size = "-log(FDR)"
  ) +
  geom_point(data = subset(dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Enriched_TFs_dotPlot_v3.pdf", dot.plot, width=3, height=4)

#### correlated TF targets for each TF ####
tf.genes <- msigdbr::msigdbr(species="Homo sapiens",category="C3", subcategory="GTRD")
tf.genes.list <- unique(tf.genes[tf.genes$gs_name %in% topGenes,]$gene_symbol) # 1731
#topGenes <- unique(dot.df[order(-dot.df$NES),]$TF)
topGenes <- rev(geneOrder)
topGeneSets <- rev(dot.df[order(dot.df$NES),]$Feature_set)
synapser::synLogin()
corr.result <- read.csv(synapser::synGet("syn64024300")$path)
corr.result <- na.omit(corr.result[corr.result$type=="RNA" & corr.result$feature %in% tf.genes.list,]) # 205
if (any(corr.result$Spearman.q == 0)) {
  corr.result[corr.result$Spearman.q == 0,]$Spearman.q <- 0.0001
}
corr.result$minusLogFDR <- -log(corr.result$Spearman.q, base = 10)
maxAbsEst <- max(abs(corr.result$Spearman.est))
maxLogFDR <- max(abs(corr.result$minusLogFDR))
library(patchwork); library(ggplot2)
gsea.dot.plots <- NULL
for (i in 1:length(topGeneSets)) {
  temp.genes <- unique(tf.genes[tf.genes$gs_name == topGeneSets[i],]$gene_symbol)
  temp.gene.corr <- na.omit(corr.result[corr.result$feature %in% temp.genes &
                                          corr.result$type=="RNA" &
                                          corr.result$N>=6,])
  if (nrow(temp.gene.corr) > 0) {
    nGenes <- length(unique(temp.gene.corr$feature))
    temp.topGenes <- unique(temp.gene.corr[temp.gene.corr$Significant,]$feature) # 48
    nTopGenes <- length(temp.topGenes)
    geneOrder <- temp.gene.corr[order(temp.gene.corr$Spearman.est),]$feature

    dot.plot <- ggplot2::ggplot(
      na.omit(temp.gene.corr),
      ggplot2::aes(
        x = Omics, y = feature, color = Spearman.est,
        size = minusLogFDR
      )
    ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
      ggplot2::geom_point() +
      ggplot2::scale_y_discrete(limits = geneOrder) +
      scale_color_gradient2(low="blue",high="red", mid = "grey", limits=c(-maxAbsEst, maxAbsEst)) +
      #viridis::scale_color_viridis() +
      theme_classic() +
      ggplot2::labs(
        x = "Omics Type",
        y = "Transcription Factor Targets",
        color = "Spearman rho", size = "-log(adj. p)"
      ) + theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(), axis.title.y=element_blank()) +
      geom_point(data = subset(temp.gene.corr, Significant), col = "black", stroke = 1.5, shape = 21)
    #plot.annot <- paste0(topGenes[i], " (", nTopGenes, " / ", nGenes, " genes correlated)")
    dot.plot <- dot.plot + ggtitle(topGenes[i])
    # if (i == 1) {
    #   gsea.dot.plots <- (dot.plot + theme(legend.position = "none"))
    # } else if (i < 6) {
    #   gsea.dot.plots <- gsea.dot.plots + (dot.plot + theme(legend.position = "none"))
    # } else if (i == 6) {
    #   gsea.dot.plots <- gsea.dot.plots + dot.plot
    # } else if (i == 7) {
    #   gsea.dot.plots2 <- (dot.plot + theme(legend.position = "none"))
    # } else if (i > 7) {
    #   gsea.dot.plots2 <- gsea.dot.plots2 + (dot.plot + theme(legend.position = "none"))
    # }
    if (is.list(dot.plot) & length(dot.plot) > 1) {
      if (is.null(gsea.dot.plots)) {
        gsea.dot.plots <- (dot.plot + theme(legend.position = "none"))
      } else if (i == ceiling(length(topGenes)/2)) {
        gsea.dot.plots <- gsea.dot.plots + dot.plot #+ plot_layout(guides = 'collect')
      } else {
        gsea.dot.plots <- gsea.dot.plots + (dot.plot + theme(legend.position = "none"))
      }
    }
  }
}
#gsea.dot.plots3 <- gsea.dot.plots/gsea.dot.plots2
#gsea.dot.plots3 <- gsea.dot.plots3 + plot_layout(guides = 'collect')
#gsea.dot.plots
#ggplot2::ggsave("Correlated_TFTs_dotPlot_patchworkCollection_minN6_tall_v2.pdf", gsea.dot.plots, width=9, height=16)
#ggplot2::ggsave("Correlated_TFTs_dotPlot_patchworkCollection_minN6_v2.pdf", gsea.dot.plots, width=9, height=9)

source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
gsea.dot.plot2 <- gsea.dot.plots + plot_layout(nrow = 2, guides='collect')
ggplot2::ggsave(paste0("TFT_correlations_dotPlot_patchworkCollection_minN6_v4_",Sys.Date(),".pdf"), gsea.dot.plot2, width=11, height=11)
gsea.dot.plot2 <- gsea.dot.plots + plot_layout(nrow = 1, guides='collect')
ggplot2::ggsave(paste0("TFT_correlations_dotPlot_patchworkCollection_minN6_v5_",Sys.Date(),".pdf"), gsea.dot.plot2, width=14, height=7)
