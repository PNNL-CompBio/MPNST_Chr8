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
gsea.rna.prot$Direction <- factor(gsea.rna.prot$Direction, levels=c("Positive", "Negative"))
bar.plot <- ggplot2::ggplot(gsea.rna.prot[gsea.rna.prot$Significant,], aes(fill = Direction, x=Omics)) + 
  geom_bar(position = "dodge", stat="count", alpha=0.75) + ggplot2::theme_classic() + 
  facet_wrap(~Collection, scales = "free") +
  scale_fill_manual(values=c("red", "blue")) + 
  ylim(c(0,15)) +
  ggplot2::labs(
    x = "Omics Type",
    y = "Number of Enriched Gene Sets",
    color = "Correlation") +
  theme(axis.title.x=element_blank()) #+
  #theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
bar.plot
ggplot2::ggsave("Chr8_numberOfGeneSets_directionGroup_barPlot_manualOrderV3.pdf", bar.plot, width = 7, height = 7, device = "pdf")

circ.bar.plot <- ggplot2::ggplot(gsea.rna.prot[gsea.rna.prot$Significant,], aes(fill = Direction, x=Omics)) + 
  geom_bar(position = "dodge", stat="count", alpha=0.75) + ggplot2::theme_classic() + 
  facet_wrap(~Collection) +
  scale_fill_manual(values=c("red", "blue")) + 
  ylim(c(0,15)) +
  ggplot2::labs(
    x = "Omics Type",
    y = "Number of Enriched Gene Sets",
    color = "Correlation") +
  theme(axis.title.x=element_blank()) + coord_polar() #+
#theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
circ.bar.plot
ggplot2::ggsave("Chr8_numberOfGeneSets_directionGroup_barPlot_manualOrderV3.pdf", bar.plot3log, width = 7, height = 7, device = "pdf")

# source: https://r-graph-gallery.com/web-circular-barplot-with-R-and-ggplot2.html
plot_df <- plyr::ddply(gsea.rna.prot[gsea.rna.prot$Significant,] , .(Collection, Omics), summarize,
                       nEnriched = n())
#plot_df[6,] <- c("PID", "RNA", 0)


plt <- ggplot(plot_df) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(5, 10, 15)),
    color = "lightgrey"
  ) + 
  # Add bars to represent the cumulative track lengths
  # str_wrap(region, 5) wraps the text so each line has at most 5 characters
  # (but it doesn't break long words!)
  geom_col(
    aes(
      x = reorder(Collection, nEnriched),
      #x = forcats::fct_infreq(stringr::str_wrap(Collection, 5)),
      y = nEnriched,
      fill = Omics
    ),
    position = "dodge", show.legend = TRUE, alpha=0.5
  ) +
  geom_segment(
    aes(
      x = reorder(Collection, nEnriched),
      y = 0,
      xend = reorder(Collection, nEnriched),
      yend = max(plot_df$nEnriched)
    ),
    linetype = "dashed",
    color = "#5A5A5A"
  ) + 
  # Make it circular!
  coord_polar() + theme_classic()

plt

plt <- plt +
  # Annotate custom scale inside plot
  annotate(
    x = 0, 
    y = 4, 
    label = "5", 
    geom = "text", 
    color = "gray12"
  ) +
  annotate(
    x = 0, 
    y = 9, 
    label = "10", 
    geom = "text", 
    color = "gray12"
  ) +
  annotate(
    x = 0, 
    y =14, 
    label = "15", 
    geom = "text", 
    color = "gray12"
  ) +
  # Scale y axis so bars don't start in the center
  scale_y_continuous(
    limits = c(-5, 15),
    expand = c(0, 0),
    breaks = c(0, 5, 10, 15)
  ) + #scale_fill_manual(values=c("lightslateblue", "lightgoldenrod1")) +
  theme(
    # Remove axis ticks and text
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    # Use gray text for the region names
    axis.text.x = element_text(color = "gray12", size = 12),
    # Move the legend to the bottom
    legend.position = "bottom",
  ) + ggtitle("Number of Enriched Gene Sets") + 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
        axis.line=element_blank())

plt
ggplot2::ggsave("Chr8_numberOfGeneSets_circularBarPlot_greenOrange_matchingVenn.pdf", plt, width = 7, height = 7)
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

ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_color=c("#F8766D", "#00BFC4"), fill_alpha = 0.5, text_size = 10, set_name_size = 7.5)
ggplot2::ggsave(paste0("Chr8_GSEA_vennDiagram_noPercentages_", Sys.Date(), "_orangeGreen.pdf"), width=7, height=7)


# # # try using BioVenn for proportional circles
# library(BioVenn)
# venn.info <- draw.venn(venn.list[[1]], venn.list[[2]], list(), title = "Enriched Gene Sets",
#               subtitle = "", xtitle = "RNA", ytitle = "Protein")
# plot(eulerr::venn(c(RNA = length(venn.info$x_only), 
#                "RNA & Protein" = length(venn.info$xy_only),
#                Protein = length(venn.info$y_only))))
#install.packages("nVennR")
# devtools::install_github("vqf/nVennR")
# library(nVennR)
# plotVenn(venn.list)

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
#mean.gsea[mean.gsea$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
#gsea.rna.prot[gsea.rna.prot$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
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
  # geneFace <- rep("plain", length(geneOrder))
  # names(geneFace) <- geneOrder
  # sigInRNA <- temp.mean.gsea[grepl("RNA", temp.mean.gsea$sig_types),]$Feature_set
  # geneFace[sigInRNA] <- "italic"
  # sigInBoth <- temp.mean.gsea[grepl("RNA, Protein", temp.mean.gsea$sig_types),]$Feature_set
  # geneFace[sigInBoth] <- "bold"
  # if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
  
  geneFace <- rep("plain", length(geneOrder))
  names(geneFace) <- geneOrder
  MYC.pathways <- geneOrder[grepl("MYC", geneOrder)]
  SHH.pathways <- geneOrder[grepl("SHH", geneOrder) | grepl("HEDGEHOG", geneOrder)]
  poi <- c(MYC.pathways, SHH.pathways)
  geneFace[poi] <- "bold"
  
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
    theme(axis.text.y=element_text(face=geneFace)) +
    geom_point(data = subset(temp.dot.df, sig), col = "black", stroke = 1.5, shape = 21) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  plot.annot <- paste0(i, " (", nTopGenes, " / ", nGenes, " gene sets enriched)")
  dot.plot <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
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
#ggplot2::ggsave("Enriched_geneSets_dotPlot_patchworkCollection_boldMYCandSHH.pdf", gsea.dot.plots, width=9, height=9)
ggplot2::ggsave("Enriched_geneSets_dotPlot_patchworkCollection_boldMYCandSHH_boldTitles.pdf", gsea.dot.plots, width=9, height=9)

# gsea.dot.plots2 <- gsea.dot.plots + plot_annotation(caption = "Key\nBold: enriched in both RNA & protein\n Italics: enriched in RNA\n Plain: enriched in Protein")
# gsea.dot.plots2
# ggplot2::ggsave("Enriched_geneSets_dotPlot_patchworkCollection_withKey.pdf", gsea.dot.plots2, width=9, height=9)

# #### Look closer at ESC_V6.5_UP_EARLY.V1_DN leading edge genes ####
# corr.result <- read.csv(synapser::synGet("syn64024300")$path)
# esc.genes <- msigdbr::msigdbr(species="Homo sapiens",category="C6")
# esc.genes <- unique(esc.genes[esc.genes$gs_name == "ESC_V6.5_UP_EARLY.V1_DN",]$gene_symbol) # 161
# esc.corr <- corr.result[corr.result$feature %in% esc.genes,] # 283
# esc.corr <- esc.corr[esc.corr$type %in% c("RNA", "Protein"),] # 124
# topGenes <- unique(esc.corr[esc.corr$sig,]$feature) # 11
# 
# mean.corr <- plyr::ddply(esc.corr[esc.corr$feature %in% topGenes,], .(feature), summarize,
#                          mean_rankVal = mean(Spearman.est, na.rm = TRUE),
#                          sig_types = paste(type, collapse = ", "))
# geneOrder <- mean.corr[order(mean.corr$mean_rankVal),]$feature
# geneFace <- rep("plain", length(geneOrder))
# names(geneFace) <- geneOrder
# sigInRNA <- mean.corr[grepl("RNA", mean.corr$sig_types),]$feature
# geneFace[sigInRNA] <- "italic"
# sigInBoth <- mean.corr[grepl("RNA, Protein", mean.corr$sig_types),]$feature
# geneFace[sigInBoth] <- "bold"
# # if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
# 
# dot.df <- corr.result[corr.result$feature %in% topGenes &
#                      corr.result$Omics %in% c("RNA", "Protein"),]
# if (any(dot.df$Spearman.q == 0)) {
#   dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
# }
# dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
# dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
# dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
# maxLogFDR.esc <- max(dot.df$minusLogFDR)
# absMaxEst.esc <- max(abs(dot.df$Spearman.est))
# esc.dot.plot <- ggplot2::ggplot(
#   na.omit(dot.df),
#   ggplot2::aes(
#     x = Omics, y = feature, color = Spearman.est,
#     size = minusLogFDR
#   )
# ) + scale_size(limits=c(0,maxLogFDR.esc), range = c(0.5,5)) +
#   ggplot2::geom_point() +
#   ggplot2::scale_y_discrete(limits = geneOrder) +
#   scale_color_gradient2(low="blue",high="red") +
#   #viridis::scale_color_viridis() +
#   theme_classic() +
#   ggplot2::labs(
#     x = "Omics Type",
#     y = "Gene",
#     color = "Spearman Correlation", size = "-log(adjusted p-value)"
#   ) + 
#   theme(axis.text.y=element_text(face=geneFace))
# esc.dot.plot
# # most are only in prot and not in RNA
# ggplot2::ggsave("Correlated_ESC_genes_dotPlot.pdf", esc.dot.plot, width=4, height=4)
# 
# #### Look closer at HALLMARK_COAGULATION leading edge genes ####
# corr.result <- read.csv(synapser::synGet("syn64024300")$path)
# coag.genes <- msigdbr::msigdbr(species="Homo sapiens",category="H")
# coag.genes <- unique(coag.genes[coag.genes$gs_name == "HALLMARK_COAGULATION",]$gene_symbol) # 138
# coag.corr <- corr.result[corr.result$feature %in% coag.genes,] # 265
# coag.corr <- coag.corr[coag.corr$type %in% c("RNA", "Protein"),] # 127
# topGenes <- unique(coag.corr[coag.corr$sig,]$feature) # 12
# 
# mean.corr <- plyr::ddply(coag.corr[coag.corr$feature %in% topGenes,], .(feature), summarize,
#                          mean_rankVal = mean(Spearman.est, na.rm = TRUE),
#                          sig_types = paste(type, collapse = ", "))
# geneOrder <- mean.corr[order(mean.corr$mean_rankVal),]$feature
# geneFace <- rep("plain", length(geneOrder))
# names(geneFace) <- geneOrder
# sigInRNA <- mean.corr[grepl("RNA", mean.corr$sig_types),]$feature
# geneFace[sigInRNA] <- "italic"
# sigInBoth <- mean.corr[grepl("RNA, Protein", mean.corr$sig_types),]$feature
# geneFace[sigInBoth] <- "bold"
# # if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
# 
# dot.df <- corr.result[corr.result$feature %in% topGenes &
#                         corr.result$Omics %in% c("RNA", "Protein"),]
# if (any(dot.df$Spearman.q == 0)) {
#   dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
# }
# dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
# dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
# dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
# maxLogFDR.coag <- max(dot.df$minusLogFDR)
# absMaxEst.coag <- max(abs(dot.df$Spearman.est))
# coag.dot.plot <- ggplot2::ggplot(
#   na.omit(dot.df),
#   ggplot2::aes(
#     x = Omics, y = feature, color = Spearman.est,
#     size = minusLogFDR
#   )
# ) + scale_size(limits=c(0,maxLogFDR.coag), range = c(0.5,5)) +
#   ggplot2::geom_point() +
#   ggplot2::scale_y_discrete(limits = geneOrder) +
#   scale_color_gradient2(low="blue",high="red") +
#   #viridis::scale_color_viridis() +
#   theme_classic() +
#   ggplot2::labs(
#     x = "Omics Type",
#     y = "Gene",
#     color = "Spearman Correlation", size = "-log(adjusted p-value)"
#   ) + 
#   theme(axis.text.y=element_text(face=geneFace))
# coag.dot.plot
# # most are only in prot and not in RNA
# ggplot2::ggsave("Correlated_coagulation_genes_dotPlot.pdf", coag.dot.plot, width=4, height=4)
# 
# #### combine esc and coag dot plots ####
# absMaxEst <- max(c(absMaxEst.coag, absMaxEst.esc))
# maxLogFDR <- max(c(maxLogFDR.coag, maxLogFDR.esc))
# esc.dot.plot2 <- esc.dot.plot + ggtitle("ESC_V6.5_UP_EARLY.V1_DN") + 
#   scale_color_gradient2(low="blue",high="red", limits = c(-absMaxEst, absMaxEst)) +
#   scale_size(limits=c(0,maxLogFDR), range = c(0.5,5))
# coag.dot.plot2 <- coag.dot.plot + ggtitle("HALLMARK_COAGULATION") + 
#   scale_color_gradient2(low="blue",high="red", limits = c(-absMaxEst, absMaxEst)) +
#   scale_size(limits=c(0,maxLogFDR), range = c(0.5,5))
# esc.coag.dot.plots <-  (coag.dot.plot2 + theme(legend.position = "none"))/ esc.dot.plot2
# esc.coag.dot.plots <- esc.coag.dot.plots + plot_annotation(caption="Key\nBold: correlated in both RNA & protein\n Italics: correlated in RNA\n Plain: correlated in Protein")
# esc.coag.dot.plots
# ggplot2::ggsave("Correlated_ESC_and_Hallmark_genes_dotPlot.pdf", esc.coag.dot.plots, width=4, height=9)

#### Look closer at MYC pathway genes ####
corr.result <- read.csv(synapser::synGet("syn64024300")$path)
hallmark <- msigdbr::msigdbr(species="Homo sapiens",category="H")
myc.hallmark <- unique(hallmark[hallmark$gs_name %in% 
                                  c("HALLMARK_MYC_TARGETS_V1",
                                    "HALLMARK_MYC_TARGETS_V2"),]$gene_symbol) # 240
pid <- msigdbr::msigdbr(species="Homo sapiens",category="C2", subcategory="PID")
myc.pid <- unique(pid[pid$gs_name == "PID_HEDGEHOG_2PATHWAY",]$gene_symbol) # 22
myc.genes <- unique(c(myc.hallmark, myc.pid)) # 262

myc.corr <- corr.result[corr.result$feature %in% myc.genes,] # 530
myc.corr <- myc.corr[myc.corr$type %in% c("RNA", "Protein"),] # 268
topGenes <- unique(myc.corr[myc.corr$sig,]$feature) # 13

mean.corr <- plyr::ddply(myc.corr[myc.corr$feature %in% topGenes,], .(feature), summarize,
                         mean_rankVal = mean(Spearman.est, na.rm = TRUE),
                         sig_types = paste(type, collapse = ", "))
geneOrder <- mean.corr[order(mean.corr$mean_rankVal),]$feature
# geneFace <- rep("plain", length(geneOrder))
# names(geneFace) <- geneOrder
# sigInRNA <- mean.corr[grepl("RNA", mean.corr$sig_types),]$feature
# geneFace[sigInRNA] <- "italic"
# sigInBoth <- mean.corr[grepl("RNA, Protein", mean.corr$sig_types),]$feature
# geneFace[sigInBoth] <- "bold"
# if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
geneFace <- rep("plain", length(geneOrder))
names(geneFace) <- geneOrder
geneFace["SHH"] <- "bold"

dot.df <- na.omit(corr.result[corr.result$feature %in% topGenes &
                        corr.result$Omics %in% c("RNA", "Protein"),])
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
maxLogFDR.myc <- max(dot.df$minusLogFDR)
absMaxEst.myc <- max(abs(dot.df$Spearman.est))
myc.dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = Omics, y = feature, color = Spearman.est,
    size = minusLogFDR
  )
) + scale_size(limits=c(0,maxLogFDR.myc), range = c(0.5,5)) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="blue",high="red") +
  #viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Gene",
    color = "Spearman rho", size = "-log(adj. p)"
  ) + geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + ggtitle("MYC") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
  theme(axis.text.y=element_text(face=geneFace))
myc.dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Correlated_MYC_genes_dotPlot_boldSHH.pdf", myc.dot.plot, width=4, height=4)

#### Look closer at HALLMARK_shhULATION leading edge genes ####
corr.result <- read.csv(synapser::synGet("syn64024300")$path)
onc.sigs <- msigdbr::msigdbr(species="Homo sapiens",category="C6")
shh.onc <- unique(onc.sigs[onc.sigs$gs_name == "GCNP_SHH_UP_LATE.V1_UP",]$gene_symbol) # 181
shh.pid <- unique(pid[pid$gs_name == "PID_HEDGEHOG_2PATHWAY",]$gene_symbol) # 22
shh.genes <- unique(c(shh.onc, shh.pid)) # 202

shh.corr <- corr.result[corr.result$feature %in% shh.genes,] # 374
shh.corr <- shh.corr[shh.corr$type %in% c("RNA", "Protein"),] # 173
topGenes <- unique(shh.corr[shh.corr$sig,]$feature) # 8

mean.corr <- plyr::ddply(shh.corr[shh.corr$feature %in% topGenes,], .(feature), summarize,
                         mean_rankVal = mean(Spearman.est, na.rm = TRUE),
                         sig_types = paste(type, collapse = ", "))
geneOrder <- mean.corr[order(mean.corr$mean_rankVal),]$feature
# geneFace <- rep("plain", length(geneOrder))
# names(geneFace) <- geneOrder
# sigInRNA <- mean.corr[grepl("RNA", mean.corr$sig_types),]$feature
# geneFace[sigInRNA] <- "italic"
# sigInBoth <- mean.corr[grepl("RNA, Protein", mean.corr$sig_types),]$feature
# geneFace[sigInBoth] <- "bold"
# if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
geneFace <- rep("plain", length(geneOrder))
names(geneFace) <- geneOrder
geneFace["SHH"] <- "bold"

dot.df <- na.omit(corr.result[corr.result$feature %in% topGenes &
                        corr.result$Omics %in% c("RNA", "Protein"),])
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
maxLogFDR.shh <- max(dot.df$minusLogFDR)
absMaxEst.shh <- max(abs(dot.df$Spearman.est))
shh.dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = Omics, y = feature, color = Spearman.est,
    size = minusLogFDR
  )
) + scale_size(limits=c(0,maxLogFDR.shh), range = c(0.5,5)) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="blue",high="red") +
  #viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Gene",
    color = "Spearman rho", size = "-log(adj. p)"
  ) + geom_point(data = subset(dot.df, sig), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + ggtitle("SHH") +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
  theme(axis.text.y=element_text(face=geneFace))
shh.dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Correlated_SHH_genes_dotPlot_boldSHH.pdf", shh.dot.plot, width=4, height=4)

#### combine esc and shh dot plots ####
absMaxEst <- max(c(absMaxEst.shh, absMaxEst.myc))
maxLogFDR <- max(c(maxLogFDR.shh, maxLogFDR.myc))
myc.dot.plot2 <- myc.dot.plot + 
  scale_color_gradient2(low="blue",high="red", limits = c(-absMaxEst, absMaxEst)) +
  scale_size(limits=c(0,maxLogFDR), range = c(0.5,5))
shh.dot.plot2 <- shh.dot.plot +  
  scale_color_gradient2(low="blue",high="red", limits = c(-absMaxEst, absMaxEst)) +
  scale_size(limits=c(0,maxLogFDR), range = c(0.5,5))
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4/guides_build_mod.R")
myc.shh.dot.plots <- (myc.dot.plot2 / shh.dot.plot2) + plot_layout(guides='collect')
#myc.shh.dot.plots <-  (shh.dot.plot2 + theme(legend.position = "none"))/ myc.dot.plot2
#myc.shh.dot.plots <- myc.shh.dot.plots + plot_annotation(caption="Key\nBold: correlated in both RNA & protein\n Italics: correlated in RNA\n Plain: correlated in Protein")
myc.shh.dot.plots
#ggplot2::ggsave("Correlated_MYC_and_SHH_genes_dotPlot.pdf", myc.shh.dot.plots, width=4, height=9)
ggplot2::ggsave("Correlated_MYC_and_SHH_genes_dotPlot_w3h6.pdf", myc.shh.dot.plots, width=3, height=6)

overlap.genes <- myc.genes[myc.genes %in% shh.genes] # 38
topSHHGenes <- unique(shh.corr[shh.corr$sig,]$feature) # 8
topMYCGenes <- unique(myc.corr[myc.corr$sig,]$feature) # 13
topOverlapGenes <- topMYCGenes[topMYCGenes %in% topSHHGenes] # SHH
# go back and bold SHH in dot plots

# all.dot.plots <- gsea.dot.plots | myc.shh.dot.plots
# all.dot.plots
# ggplot2::ggsave("Fig2_dotPlots.pdf", all.dot.plots, width=9, height=9)
