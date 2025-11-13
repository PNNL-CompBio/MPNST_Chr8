# chr8 MPNST: figure 3: kinase enrichment
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase")

synapser::synLogin()

#all.gsea <- read.table(synapser::synGet("syn66279670")$path, sep="\t", header=TRUE) # no sig kinases with MIT
plot.data <- read.csv(synapser::synGet("syn66279699")$path)

#### 3. top up/dn diffexp gene sets ####
dot.df <- na.omit(plot.data[plot.data$FDR_q_value <= 0.25 & plot.data$p_value <= 0.05,]) # 12
topGenes <- unique(dot.df$Feature_set)
geneOrder <- dot.df[order(dot.df$NES),]$Feature_set
dot.df$minusLogP <- -log10(dot.df$p_value)
dot.df$minusLogFDR <- -log10(dot.df$FDR_q_value)
absMaxNES <- max(abs(dot.df$NES)) # 2.08
maxLogFDR <- max(abs(dot.df$minusLogFDR)) # slightly over 3.012781 so rounding to 3.0128
# for same scales as TF plot:
maxLogFDR <- 3.0128
absMaxNES <- 2.34
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(x="",y = Feature_set, color = NES,
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
    y = "Kinase",
    color = "NES", size = "-log(FDR)"
  ) + theme(axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank()) +
  geom_point(data = dot.df, col = "black", stroke = 1.5, shape = 21)
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Enriched_kinases_dotPlot_2025-04-20_sameScaleAsTF.pdf", dot.plot, width=2.2, height=0.75)
dot.plot2 <- dot.plot + theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical")
ggplot2::ggsave("Enriched_kinases_dotPlot_horizontalLegend_2025-04-20.pdf", dot.plot2, width=1.2, height=2.1)

dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(x="",y = Feature_set, color = NES,
               size = minusLogFDR
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="blue",high="red", mid="grey", limits = c(-absMaxNES, absMaxNES)) +
  #viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Kinase",
    color = "NES", size = "-log(FDR)"
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank()) +
  geom_point(data = dot.df, col = "black", stroke = 1.5, shape = 21)
dot.plot2 <- dot.plot + coord_flip() + ggtitle(paste0("Kinases (",nrow(dot.df)," / ",nrow(plot.data)," Enriched)")) +
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        axis.title.y = element_blank(), legend.box = "vertical", axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1))
ggplot2::ggsave("Enriched_kinases_dotPlot_horizontal_2025-04-20.pdf", dot.plot2, width=0.75, height=2.5)

#### correlated kinase targets for each kinase ####
gmt2 <- readRDS("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/gmt2.rds")[[1]]

topGenes <- rev(geneOrder)
corr.result <- read.csv(synapser::synGet("syn66226338")$path)
corr.result$minusLogFDR <- -log(corr.result$Spearman.q, base = 10)
corr.result$sig <- FALSE
corr.result[corr.result$Spearman.q <= 0.05,]$sig <- TRUE
topSites <- unique(unlist(gmt2$genesets[topGenes]))
corr.result <- na.omit(corr.result[corr.result$SUB_SITE %in% topSites & corr.result$N>=6,])
maxAbsEst <- ceiling(max(abs(corr.result$Spearman.est)))
maxLogFDR <- max(abs(corr.result$minusLogFDR))
library(patchwork); library(ggplot2)
gsea.dot.plots <- NULL
for (i in 1:length(topGenes)) {
  temp.sites <- unlist(gmt2$genesets[[topGenes[i]]])
  temp.gene.corr <- na.omit(corr.result[corr.result$SUB_SITE %in% temp.sites,]) %>% slice_max(abs(Spearman.est), n=10)
  if (nrow(temp.gene.corr) > 0) {
    nGenes <- length(unique(temp.gene.corr$SUB_SITE))
    nDuplicateGenes <- nrow(temp.gene.corr) - nGenes
    temp.topGenes <- unique(temp.gene.corr[temp.gene.corr$sig,]$SUB_SITE) # 48
    nTopGenes <- length(temp.topGenes)
    cat(topGenes[i], "has", nTopGenes, "correlated out of", nGenes, "with", nDuplicateGenes, "duplicates\n")
    geneOrder <- temp.gene.corr[order(temp.gene.corr$Spearman.est),]$SUB_SITE
    
    dot.plot <- ggplot2::ggplot(
      temp.gene.corr,
      ggplot2::aes(
        x = "Phospho", y = SUB_SITE, color = Spearman.est,
        size = -log(Spearman.q, base=10)
      )
    ) + ggplot2::geom_point() +
      ggplot2::scale_y_discrete(limits = geneOrder) +
      scale_color_gradient2(low="blue",high="red", mid="grey",limits=c(-maxAbsEst, maxAbsEst)) +
      scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
      #viridis::scale_color_viridis() +
      theme_classic() +
      ggplot2::labs(
        x = "Omics Type",
        y = "Kinase Substrate Sites",
        color = "Spearman rho", size = "-log(adj. p)"
      ) + theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(), axis.title.y=element_blank()) +
      geom_point(data = subset(temp.gene.corr, sig), col = "black", stroke = 1.5, shape = 21)
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
#ggplot2::ggsave("Correlated_kinaseSubsites_dotPlot_patchworkCollection_minN6_tall_v2.pdf", gsea.dot.plots, width=10, height=16)
#ggplot2::ggsave("Correlated_kinaseSubsites_dotPlot_patchworkCollection_minN6.pdf", gsea.dot.plots, width=10, height=10)
# actual min N is 10 but only 1 has 10, 4 have 11, and rest (171 - 5) have 12
#source("guides_build_mod.R")
#gsea.dot.plot2 <- gsea.dot.plots + guide_area() + plot_layout(nrow = 2, guides='collect')
#ggplot2::ggsave("Correlated_kinaseSubsites_dotPlot_patchworkCollection_minN6_wide_v3.pdf", gsea.dot.plot2, width=16, height=10)
#ggplot2::ggsave(paste0("kinaseSubsite_correlations_dotPlot_patchworkCollection_minN6_v2_",Sys.Date(),".pdf"), gsea.dot.plot2, width=12, height=10)
gsea.dot.plot2 <- gsea.dot.plots + plot_layout(nrow = 1, guides='collect')
ggplot2::ggsave(paste0("kinaseSubsite_correlations_dotPlot_patchworkCollection_minN6_wide_v3_",Sys.Date(),".pdf"), gsea.dot.plot2, width=3, height=2)
#ggplot2::ggsave(paste0("kinaseSubsite_correlations_dotPlot_patchworkCollection_minN6_v3_",Sys.Date(),".pdf"), gsea.dot.plot2, width=12, height=10)
gsea.dot.plot2 <- (gsea.dot.plots + plot_layout(guides='collect')) + ggplot2::theme(legend.position="bottom", legend.direction = "horizontal")
ggplot2::ggsave(paste0("kinaseSubsite_correlations_dotPlot_patchworkCollection_minN6_wide_v3_horizontalLegend",Sys.Date(),".pdf"), gsea.dot.plot2, width=5.5, height=2)
