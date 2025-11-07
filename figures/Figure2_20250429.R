# chr8 MPNST: figure 2
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)
#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/circBar.R")
#setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_2")
synapser::synLogin()

#### 0. compile GSEA ####
gsea.cn <- read.csv(synapser::synGet("syn66227265")$path)
gsea.cn$Omics <- "Copy_number"
gsea.cn$Collection <- "Positional"

gsea.rna.tf <- read.csv(synapser::synGet("syn66226952")$path)
gsea.rna.tf$Omics <- "RNA"
gsea.rna.tf$Collection <- "TFT_GTRD"
gsea.rna.hall <- read.csv(synapser::synGet("syn66226874")$path)
gsea.rna.hall$Omics <- "RNA"
gsea.rna.hall$Collection <- "Hallmark"

gsea.prot.hall <- read.csv(synapser::synGet("syn66224811")$path)
gsea.prot.hall$Omics <- "Protein"
gsea.prot.hall$Collection <- "Hallmark"

ksea <- read.csv(synapser::synGet("syn66279699")$path)
ksea$Omics <- "Phospho"
ksea$Collection <- "PhosphoSitePlus"
all.gsea <- rbind(gsea.cn, gsea.rna.hall,
                  gsea.prot.hall, gsea.rna.tf, ksea)
write.csv(all.gsea, "SupplementaryTable2_GSEA_and_KSEA.csv", row.names=FALSE)
all.gsea <- read.csv("SupplementaryTable2_GSEA_and_KSEA.csv")
all.gsea$Significant <- FALSE
all.gsea[all.gsea$p_value <= 0.05 & all.gsea$FDR_q_value <= 0.25,]$Significant <- TRUE

#### 3. top up/dn diffexp gene sets ####
# maybe do for each collection and then combine
collections <- c("Hallmark", "TFT_GTRD", "PhosphoSitePlus")
maxAbsNES <- max(abs(all.gsea[all.gsea$Collection %in% collections,]$NES))
all.gsea$minusLogFDR <- 4
all.gsea[all.gsea$FDR_q_value != 0,]$minusLogFDR <- -log10(all.gsea[all.gsea$FDR_q_value!=0,]$FDR_q_value)
all.gsea$Feature_set <- sub("HALLMARK_", "", all.gsea$Feature_set)
all.gsea$Feature_set <- sub("_TARGET_GENES", "", all.gsea$Feature_set)
maxLogFDR <- max(all.gsea[all.gsea$Collection %in% collections,]$minusLogFDR)

mean.gsea <- plyr::ddply(all.gsea, .(Feature_set, Collection), summarize,
                         mean_NES = mean(NES, na.rm=TRUE),
                         maxAbsNES = max(abs(NES[Significant]), na.rm=TRUE),
                         maxNES = max(NES[Significant], na.rm=TRUE),
                         minNES = min(NES[Significant], na.rm=TRUE),
                         sig_types = paste0(Omics[Significant], collapse = ", "))

n.cutoff <- c(25)
all.dot.plots <- list()
for (n in n.cutoff) {
  dot.plots <- list()
  for (i in collections) {
    nGenes <- length(unique(all.gsea[all.gsea$Collection == i,]$Feature_set))
    topGenes <- unique(all.gsea[all.gsea$Significant &
                                       all.gsea$Collection == i,]$Feature_set)
    nTopGenes <- length(topGenes)

    # let's cap it at 20 pathways per plot
    temp.max <- mean.gsea[mean.gsea$Feature_set %in% topGenes & 
                            mean.gsea$Collection == i,] %>% slice_max(maxAbsNES, n = n)
    geneOrder <- temp.max[order(temp.max$maxNES),]$Feature_set
    
    temp.dot.df <- all.gsea[all.gsea$Collection == i & all.gsea$Feature_set %in% geneOrder,]
    temp.dot.df$Significant <- FALSE
    temp.dot.df[temp.dot.df$p_value <= 0.05 & temp.dot.df$FDR_q_value <= 0.25,]$Significant <- TRUE
    if (any(temp.dot.df$FDR_q_value == 0)) {
      temp.dot.df[temp.dot.df$FDR_q_value == 0,]$FDR_q_value <- 0.0001
    }
    temp.dot.df$minusLogP <- -log(temp.dot.df[,"p_value"], base = 10)
    temp.dot.df$minusLogFDR <- -log(temp.dot.df[,"FDR_q_value"], base = 10)
    temp.dot.df$Omics <- factor(temp.dot.df$Omics, levels=c("RNA", "Protein", "Phospho"))
    dot.plot <- ggplot2::ggplot(
      na.omit(temp.dot.df),
      ggplot2::aes(
        x = Omics, y = Feature_set, color = NES,
        size = minusLogFDR
      )
    ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
      ggplot2::geom_point() +
      ggplot2::scale_y_discrete(limits = geneOrder) +
      scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
      theme_classic() +
      ggplot2::labs(color = "NES", size = "-log(FDR)") + 
      geom_point(data = subset(temp.dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1, vjust=1))
    if (i == "Hallmark") {
      plot.annot <- paste0(i, " Gene Sets\n(", nTopGenes, " / ", nGenes, " enriched)") 
    } else if (i == "TFT_GTRD") {
      plot.annot <- paste0("Transcription Factors\n(", nTopGenes, " / ", nGenes, " enriched)")
    } else if (i == "PhosphoSitePlus") {
      plot.annot <- paste0("Kinases\n(", nTopGenes, " / ", nGenes, " enriched)")
    }
    dot.plots[[i]] <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
    ggplot2::ggsave(paste0(i,"_Enriched_geneSets_dotPlot_patchworkCollection_",25,"maxByNES_2025-04-18.pdf"), dot.plots[[i]], width=2, height=length(geneOrder)/5)
    ggplot2::ggsave(paste0(i,"_Enriched_geneSets_dotPlot_patchworkCollection_",25,"maxByNES_2025-04-18_wider.pdf"), dot.plots[[i]], width=4.2, height=length(geneOrder)/5)
  }
  all.dot.plots[[as.character(n)]] <- dot.plots
}
reg.dot.plots <- (all.dot.plots$`25`$TFT_GTRD / all.dot.plots$`25`$PhosphoSitePlus) + 
  plot_layout(guides="collect", heights=c(5,0.5))
reg.dot.plots
gsea.dot.plots <- all.dot.plots$`25`$Hallmark + reg.dot.plots + plot_layout(guides="collect")
gsea.dot.plots
ggplot2::ggsave(paste0("Enriched_geneSets_dotPlot_patchworkCollection_",25,"maxByNES_2025-04-18.pdf"), gsea.dot.plots, width=6, height=6)
