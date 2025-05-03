# chr8 MPNST: figure 3
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_TF")

synapser::synLogin()

#### 1. volcano plot of # of GSEA_TFT_GTRD results ####
plot.data <- read.csv(synapser::synGet("syn66226952")$path)
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
sig.plot.data <- na.omit(plot.data[plot.data$FDR_q_value <= 0.25 & plot.data$p_value <= 0.05,])
dot.df <- na.omit(plot.data[plot.data$FDR_q_value <= 0.25 & plot.data$p_value <= 0.05,]) %>% slice_max(abs(NES), n=10)
topGenes <- unique(dot.df$Feature_set)
geneOrder <- dot.df[order(dot.df$NES),]$TF
dot.df$minusLogP <- -log(dot.df[,"p_value"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"FDR_q_value"], base = 10)
maxLogFDR <- max(abs(dot.df$minusLogFDR)) # 2.74...
absMaxNES <- max(abs(dot.df$NES)) # 2.34
# for same scales as kinase plot:
maxLogFDR <- 3.0128
absMaxNES <- 2.34
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = "RNA", y = TF, color = NES,
    size = minusLogFDR
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
  scale_color_gradient2(low="blue",high="red", mid="grey", limits = c(-absMaxNES, absMaxNES)) +
  #viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Transcription Factor",
    color = "NES", size = "-log(FDR)"
  ) +
  geom_point(data = dot.df, col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Enriched_TFs_dotPlot_top10_sameScaleAsKinase.pdf", dot.plot, width=2.1, height=2.1)

dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = "RNA", y = TF, color = NES,
    size = minusLogFDR
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
  scale_color_gradient2(low="blue",high="red") +
  #viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Transcription Factor",
    color = "NES", size = "-log(FDR)"
  ) +
  geom_point(data = dot.df, col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())
dot.plot2 <- dot.plot + coord_flip() + ggtitle(paste0("Transcription Factors (",nrow(sig.plot.data)," / ",nrow(plot.data)," Enriched)")) +
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        axis.title.y = element_blank(), legend.box = "vertical",
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1))
ggplot2::ggsave("Enriched_TFs_dotPlot_horizontal_2025-04-20.pdf", dot.plot2, width=2.5, height=2.5)

#### correlated TF targets for each TF ####
tf.genes <- msigdbr::msigdbr(collection="C3", subcollection="GTRD")
tf.genes.list <- unique(tf.genes[tf.genes$gs_name %in% topGenes,]$gene_symbol)
#topGenes <- unique(dot.df[order(-dot.df$NES),]$TF)
topGenes <- rev(geneOrder)
topGeneSets <- rev(dot.df[order(dot.df$NES),]$Feature_set)
synapser::synLogin()
corr.result <- read.csv(synapser::synGet("syn66226866")$path)
corr.result <- na.omit(corr.result[corr.result$Gene %in% tf.genes.list &
                                     corr.result$N>=6,]) # 691
corr.result$minusLogFDR <- -log(corr.result$Spearman.q, base = 10)
if (any(corr.result$Spearman.q == 0)) {
  corr.result[corr.result$Spearman.q == 0,]$minusLogFDR <- ceiling(max(corr.result[corr.result$Spearman.q != 0,]$minusLogFDR))
}
maxAbsEst <- max(abs(corr.result$Spearman.est))
maxLogFDR <- max(abs(corr.result$minusLogFDR))
library(patchwork); library(ggplot2)
corr.result$Significant <- FALSE
corr.result[corr.result$Spearman.q <= 0.05,]$Significant <- TRUE
gsea.dot.plots <- NULL
for (i in 1:length(topGeneSets)) {
  temp.genes <- unique(tf.genes[tf.genes$gs_name == topGeneSets[i],]$gene_symbol)
  temp.gene.corr <- na.omit(corr.result[corr.result$Gene %in% temp.genes,]) %>% slice_max(abs(Spearman.est),n=10)
  if (nrow(temp.gene.corr) > 0) {
    nGenes <- length(unique(temp.gene.corr$Gene))
    temp.topGenes <- unique(temp.gene.corr[temp.gene.corr$Significant,]$Gene) # 48
    nTopGenes <- length(temp.topGenes)
    geneOrder <- temp.gene.corr[order(temp.gene.corr$Spearman.est),]$Gene

    dot.plot <- ggplot2::ggplot(
      temp.gene.corr,
      ggplot2::aes(
        x = "RNA", y = Gene, color = Spearman.est,
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
      ) + theme(axis.title.x=element_blank(), #legend.box = "horizontal",
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

#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
gsea.dot.plot2 <- gsea.dot.plots + guide_area() + plot_layout(nrow = 2, guides='collect')
ggplot2::ggsave(paste0("TFT_correlations_dotPlot_patchworkCollection_minN6_v4_",Sys.Date(),".pdf"), gsea.dot.plot2, width=7, height=3.5)
gsea.dot.plot2 <- gsea.dot.plots + plot_layout(nrow = 1, guides='collect')
ggplot2::ggsave(paste0("TFT_correlations_dotPlot_patchworkCollection_minN6_v5_",Sys.Date(),".pdf"), gsea.dot.plot2, width=14, height=2)

#### also do for all TFs ####
tf.genes.list <- unique(tf.genes[tf.genes$gs_name %in% 
                                   plot.data[plot.data$p_value <= 0.05 & 
                                               plot.data$FDR_q_value <= 0.25,]$Feature_set,]$gene_symbol)
topGeneSets <- rev(plot.data[order(plot.data$NES),]$Feature_set)
topGenes <- sub("_TARGET_GENES", "", topGeneSets)
synapser::synLogin()
corr.result <- read.csv(synapser::synGet("syn66226866")$path)
corr.result <- na.omit(corr.result[corr.result$Gene %in% tf.genes.list &
                                     corr.result$N>=6,]) # 691
corr.result$minusLogFDR <- -log(corr.result$Spearman.q, base = 10)
if (any(corr.result$Spearman.q == 0)) {
  corr.result[corr.result$Spearman.q == 0,]$minusLogFDR <- ceiling(max(corr.result[corr.result$Spearman.q != 0,]$minusLogFDR))
}
maxAbsEst <- max(abs(corr.result$Spearman.est))
maxLogFDR <- max(abs(corr.result$minusLogFDR))
library(patchwork); library(ggplot2)
corr.result$Significant <- FALSE
corr.result[corr.result$Spearman.q <= 0.05,]$Significant <- TRUE
gsea.dot.plots <- NULL
for (i in 1:length(topGeneSets)) {
  temp.genes <- unique(tf.genes[tf.genes$gs_name == topGeneSets[i],]$gene_symbol)
  temp.gene.corr <- na.omit(corr.result[corr.result$Gene %in% temp.genes,]) %>% slice_max(abs(Spearman.est),n=10)
  if (nrow(temp.gene.corr) > 0) {
    nGenes <- length(unique(temp.gene.corr$Gene))
    temp.topGenes <- unique(temp.gene.corr[temp.gene.corr$Significant,]$Gene) # 48
    nTopGenes <- length(temp.topGenes)
    geneOrder <- temp.gene.corr[order(temp.gene.corr$Spearman.est),]$Gene
    
    dot.plot <- ggplot2::ggplot(
      temp.gene.corr,
      ggplot2::aes(
        x = "RNA", y = Gene, color = Spearman.est,
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
      ) + theme(axis.title.x=element_blank(), #legend.box = "horizontal",
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

#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
gsea.dot.plot2 <- gsea.dot.plots + plot_layout(nrow = 20, guides='collect')
gsea.dot.plot2
ggplot2::ggsave(paste0("All_TFT_correlations_dotPlot_patchworkCollection_minN6_",Sys.Date(),".pdf"), gsea.dot.plot2, width=40, height=32)

tf.genes.list <- unique(tf.genes[tf.genes$gs_name == "ZNF22_TARGET_GENES",]$gene_symbol)
topGeneSets <- "ZNF22_TARGET_GENES"
topGenes <- sub("_TARGET_GENES", "", topGeneSets)
synapser::synLogin()
corr.result <- read.csv(synapser::synGet("syn66226866")$path)
corr.result <- na.omit(corr.result[corr.result$Gene %in% tf.genes.list &
                                     corr.result$N>=6,]) # 691
corr.result$minusLogFDR <- -log(corr.result$Spearman.q, base = 10)
if (any(corr.result$Spearman.q == 0)) {
  corr.result[corr.result$Spearman.q == 0,]$minusLogFDR <- ceiling(max(corr.result[corr.result$Spearman.q != 0,]$minusLogFDR))
}
maxAbsEst <- max(abs(corr.result$Spearman.est))
maxLogFDR <- max(abs(corr.result$minusLogFDR))
library(patchwork); library(ggplot2)
corr.result$Significant <- FALSE
corr.result[corr.result$Spearman.q <= 0.05,]$Significant <- TRUE
gsea.dot.plots <- NULL
for (i in 1:length(topGeneSets)) {
  temp.genes <- unique(tf.genes[tf.genes$gs_name == topGeneSets[i],]$gene_symbol)
  temp.gene.corr <- na.omit(corr.result[corr.result$Gene %in% temp.genes,]) #%>% slice_max(abs(Spearman.est),n=10)
  if (nrow(temp.gene.corr) > 0) {
    nGenes <- length(unique(temp.gene.corr$Gene))
    temp.topGenes <- unique(temp.gene.corr[temp.gene.corr$Significant,]$Gene) # 48
    nTopGenes <- length(temp.topGenes)
    geneOrder <- temp.gene.corr[order(temp.gene.corr$Spearman.est),]$Gene
    
    dot.plot <- ggplot2::ggplot(
      temp.gene.corr,
      ggplot2::aes(
        x = "RNA", y = Gene, color = Spearman.est,
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
      ) + theme(axis.title.x=element_blank(), #legend.box = "horizontal",
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

#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
#gsea.dot.plot2 <- gsea.dot.plots + plot_layout(nrow = 10, guides='collect')
#gsea.dot.plot2
ggplot2::ggsave(paste0("ZNF22_TFT_correlations_dotPlot_patchworkCollection_minN6_",Sys.Date(),".pdf"), dot.plot, width=2.75, height=12)

tf.genes.list <- unique(tf.genes[tf.genes$gs_name == "ZNF22_TARGET_GENES",]$gene_symbol)
topGeneSets <- "ZNF22_TARGET_GENES"
topGenes <- sub("_TARGET_GENES", "", topGeneSets)
synapser::synLogin()
corr.result <- read.csv(synapser::synGet("syn66226866")$path)
corr.result <- na.omit(corr.result[corr.result$Gene %in% tf.genes.list,]) # 691
corr.result$minusLogFDR <- -log(corr.result$Spearman.q, base = 10)
if (any(corr.result$Spearman.q == 0)) {
  corr.result[corr.result$Spearman.q == 0,]$minusLogFDR <- ceiling(max(corr.result[corr.result$Spearman.q != 0,]$minusLogFDR))
}
maxAbsEst <- max(abs(corr.result$Spearman.est))
maxLogFDR <- max(abs(corr.result$minusLogFDR))
library(patchwork); library(ggplot2)
corr.result$Significant <- FALSE
corr.result[corr.result$Spearman.q <= 0.05,]$Significant <- TRUE
gsea.dot.plots <- NULL
for (i in 1:length(topGeneSets)) {
  temp.genes <- unique(tf.genes[tf.genes$gs_name == topGeneSets[i],]$gene_symbol)
  temp.gene.corr <- na.omit(corr.result[corr.result$Gene %in% temp.genes,]) #%>% slice_max(abs(Spearman.est),n=10)
  if (nrow(temp.gene.corr) > 0) {
    nGenes <- length(unique(temp.gene.corr$Gene))
    temp.topGenes <- unique(temp.gene.corr[temp.gene.corr$Significant,]$Gene) # 48
    nTopGenes <- length(temp.topGenes)
    geneOrder <- temp.gene.corr[order(temp.gene.corr$Spearman.est),]$Gene
    
    dot.plot <- ggplot2::ggplot(
      temp.gene.corr,
      ggplot2::aes(
        x = "RNA", y = Gene, color = Spearman.est,
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
      ) + theme(axis.title.x=element_blank(), #legend.box = "horizontal",
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

#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
#gsea.dot.plot2 <- gsea.dot.plots + plot_layout(nrow = 10, guides='collect')
#gsea.dot.plot2
ggplot2::ggsave(paste0("ZNF22_TFT_correlations_dotPlot_patchworkCollection_",Sys.Date(),".pdf"), dot.plot, width=2.75, height=14)
