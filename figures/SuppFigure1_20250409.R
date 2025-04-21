# chr8 MPNST: figure 1
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(biomaRt);library(RIdeogram);library(viridis)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/SuppFig_1")

synapser::synLogin()

#### 1. histogram of feature correlations ####
path.map <- list("Copy Number" = "syn66227257",
                 "RNA-Seq" = "syn66226866",
                 "Protein" = "syn66224803",
                 "Phospho" = "syn66226338")

inputs <- list()
for (i in 1:length(path.map)) {
  temp.degs <- read.csv(synapser::synGet(path.map[[i]])$path)
  temp.degs$Omics <- names(path.map)[i]
  
  # dot plot of top genes
  temp.degs$Direction <- "Positive"
  if (nrow(temp.degs[!is.na(temp.degs$Spearman.est) & 
                     temp.degs$Spearman.est < 0,]) > 0) {
    temp.degs[!is.na(temp.degs$Spearman.est) &
                temp.degs$Spearman.est < 0,]$Direction <- "Negative"
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
  nonFeatureNames <- colnames(temp.degs)[2:ncol(temp.degs)][colnames(temp.degs)[2:ncol(temp.degs)] != "Feature"]
  temp.degs <- temp.degs[,c("Feature", nonFeatureNames)]
  inputs[[names(path.map)[i]]] <- temp.degs
}

library(ggplot2); library(patchwork)
hist.plots <- NULL
for (i in 1:length(inputs)) {
  temp.df <- na.omit(inputs[[i]])
  temp.df$Significance <- "Adjusted p > 0.05"
  if (any(temp.df$Spearman.q <= 0.05)) {
    temp.df[which(temp.df$Spearman.q <= 0.05),]$Significance <- "Adjusted p <= 0.05"
  }
  n.distinct <- length(unique(temp.df$Spearman.est))
  n.ties <- nrow(temp.df) - n.distinct
  perc.ties <- round(n.ties * 100 / nrow(temp.df),0)
  temp.plot <- ggplot2::ggplot(temp.df, aes(x=Spearman.est, 
                                            group=Significance, 
                                            fill=Significance)) +
    geom_histogram(position="identity", alpha=0.5) + theme_classic() +
    scale_fill_manual(values = c("#00BFC4", "#F8766D"), breaks = c("Adjusted p <= 0.05", "Adjusted p > 0.05")) +
    xlab("Spearman Correlation Estimates") + ggtitle(paste0(names(inputs)[i], "\n(", perc.ties, "% tied)")) +
    theme(axis.text = element_text(size=12), axis.title = element_text(size=16), legend.title = element_text(size=16),
          legend.text=element_text(size=12), title = element_text(size = 24, hjust = 0.5),
          legend.position = "bottom")
  if (is.null(hist.plots)) {
    hist.plots <- temp.plot
  } else {
    hist.plots <- hist.plots / temp.plot
  } 
  #ggsave(paste0("Chr8_",names(inputs)[i],"_Spearman_histogram_", Sys.Date(),".pdf"),temp.plot, width=7, height=7)
}
hist.plots2 <- hist.plots + plot_layout(nrow=2, ncol=2, guides = 'collect') & theme(legend.position='bottom')
hist.plots2 <- hist.plots2 + plot_annotation(tag_levels='A')
ggplot2::ggsave("FeatureCorrelationHistograms_patchworkOmics.pdf", hist.plots2, width=10, height=10)
