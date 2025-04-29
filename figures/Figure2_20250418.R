# chr8 MPNST: figure 2
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/circBar.R")
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_2")
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
gsea.rna.kegg <- read.csv(synapser::synGet("syn66226901")$path)
gsea.rna.kegg$Omics <- "RNA"
gsea.rna.kegg$Collection <- "KEGG"
gsea.rna.pid <- read.csv(synapser::synGet("syn66226944")$path) # none sig, maybe use KEGG instead
gsea.rna.pid$Omics <- "RNA"
gsea.rna.pid$Collection <- "PID"
gsea.rna.onc <- read.csv(synapser::synGet("syn66226922")$path)
gsea.rna.onc$Omics <- "RNA"
gsea.rna.onc$Collection <- "Oncogenic_signatures"
gsea.rna.wp <- read.csv(synapser::synGet("syn66279012")$path)
gsea.rna.wp$Omics <- "RNA"
gsea.rna.wp$Collection <- "WikiPathways"
gsea.prot.hall <- read.csv(synapser::synGet("syn66224811")$path)
gsea.prot.hall$Omics <- "Protein"
gsea.prot.hall$Collection <- "Hallmark"
gsea.prot.kegg <- read.csv(synapser::synGet("syn66224890")$path)
gsea.prot.kegg$Omics <- "Protein"
gsea.prot.kegg$Collection <- "KEGG"
gsea.prot.pid <- read.csv(synapser::synGet("syn66225013")$path) # only 1 sig
gsea.prot.pid$Omics <- "Protein"
gsea.prot.pid$Collection <- "PID"
gsea.prot.onc <- read.csv(synapser::synGet("syn66224966")$path)
gsea.prot.onc$Omics <- "Protein"
gsea.prot.onc$Collection <- "Oncogenic_signatures"
gsea.prot.wp <- read.csv(synapser::synGet("syn66278971")$path)
gsea.prot.wp$Omics <- "Protein"
gsea.prot.wp$Collection <- "WikiPathways"
all.gsea <- rbind(gsea.cn, gsea.rna.hall, #gsea.rna.wp, gsea.rna.onc,
                  gsea.prot.hall)#, gsea.prot.wp, gsea.prot.onc)
write.csv(all.gsea, "SupplementaryTable2_GSEA.csv", row.names=FALSE)

#### 1. bar plot of # of enriched gene sets stacked by collection ####
gsea.rna.prot <- all.gsea[all.gsea$Omics %in% c("RNA", "Protein") &
                            all.gsea$Collection %in% c("Hallmark", "WikiPathways", "Oncogenic_signatures"),]
gsea.rna.prot[gsea.rna.prot$Collection=="Oncogenic_signatures",]$Collection <- "Oncogenic"
gsea.rna.prot$Omics <- factor(gsea.rna.prot$Omics, levels=c("RNA", "Protein"))
gsea.rna.prot$Direction <- "Positive"
gsea.rna.prot[gsea.rna.prot$NES < 0,]$Direction <- "Negative"
gsea.rna.prot$Direction <- factor(gsea.rna.prot$Direction, levels=c("Positive", "Negative"))
gsea.rna.prot$Significant <- FALSE
gsea.rna.prot[gsea.rna.prot$p_value <= 0.05 & gsea.rna.prot$FDR_q_value <= 0.25,]$Significant <- TRUE

# source: https://r-graph-gallery.com/web-circular-barplot-with-R-and-ggplot2.html
plot_df <- plyr::ddply(gsea.rna.prot[gsea.rna.prot$Significant,] , .(Collection, Omics), dplyr::summarize,
                       nEnriched = dplyr::n())
plot_df$log2N <- log2(plot_df$nEnriched) # wikipathways has a lot more than the others

plt <- ggplot(plot_df) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(2,3,4,5,6)),
    color = "lightgrey"
  ) + 
  # Add bars to represent the cumulative track lengths
  # str_wrap(region, 5) wraps the text so each line has at most 5 characters
  # (but it doesn't break long words!)
  geom_col(
    aes(
      x = reorder(Collection, log2N),
      #x = forcats::fct_infreq(stringr::str_wrap(Collection, 5)),
      y = log2N,
      fill = Omics
    ),
    position = "dodge", show.legend = TRUE, alpha=0.5
  ) +
  # geom_segment(
  #   aes(
  #     x = reorder(Collection, nEnriched),
  #     y = 0,
  #     xend = reorder(Collection, nEnriched),
  #     yend = max(plot_df$nEnriched)
  #   ),
  #   linetype = "dashed",
  #   color = "#5A5A5A"
  # ) + 
  # Make it circular!
  coord_polar() + theme_classic()

plt

plt <- plt +
  # Annotate custom scale inside plot
  annotate(
    x = 0, 
    y = 1.5, 
    label = as.character(2^2), 
    geom = "text", 
    color = "gray12"
  ) +
  annotate(
    x = 0, 
    y = 2.5, 
    label = as.character(2^3), 
    geom = "text", 
    color = "gray12"
  ) +
  annotate(
    x = 0, 
    y = 3.5, 
    label = as.character(2^4), 
    geom = "text", 
    color = "gray12"
  ) +
  annotate(
    x = 0, 
    y = 4.5, 
    label = as.character(2^5), 
    geom = "text", 
    color = "gray12"
  ) +
  annotate(
    x = 0, 
    y = 5.5, 
    label = as.character(2^6), 
    geom = "text", 
    color = "gray12"
  ) +
  # Scale y axis so bars don't start in the center
  scale_y_continuous(
    limits = c(-3, 6),
    expand = c(0, 0)#,
    #breaks = seq(2,6,by=1)
  ) + #scale_fill_manual(values=c("lightslateblue", "lightgoldenrod1")) +
  theme(
    # Remove axis ticks and text
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    # Use gray text for the region names
    axis.text.x = element_text(color = "gray12", size = 16), # was size 12
    # Move the legend to the bottom
    legend.position = "bottom",
  ) + ggtitle("Number of Enriched Gene Sets") + 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=24), # was size 16
        axis.line=element_blank())

plt
ggplot2::ggsave("Chr8_numberOfGeneSets_circularBarPlot_greenOrange_matchingVenn.pdf", plt, width = 7, height = 7)

#### 2. venn diagram of diffexp features ####
venn.list <- list()
for (i in 1:length(c("RNA", "Protein"))) {
  temp.name <- c("RNA", "Protein")[i]
  venn.list[[temp.name]] <- unique(gsea.rna.prot[gsea.rna.prot$Omics == temp.name &
                                                   gsea.rna.prot$Significant,]$Feature_set)
}
ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_color=c("#F8766D", "#00BFC4"), fill_alpha = 0.5, text_size = 10, set_name_size = 7.5)
ggplot2::ggsave(paste0("Chr8_GSEA_vennDiagram_noPercentages_", Sys.Date(), "_orangeGreen.pdf"), width=7, height=7)

#### 3. top up/dn diffexp gene sets ####
mean.gsea <- plyr::ddply(gsea.rna.prot, .(Feature_set), summarize,
                         mean_NES = mean(NES, na.rm=TRUE),
                         sig_types = paste0(Omics[Significant], collapse = ", "))

# maybe do for each collection and then combine
collections <- c("Oncogenic", "Hallmark", "WikiPathways")
all.sig.onc <- stringr::str_split_i(all.gsea[all.gsea$Collection == "Oncogenic_signatures" & all.gsea$p_value <= 0.05 &
                          all.gsea$FDR_q_value <= 0.25,]$Feature_set, "_", 1)
dup.sig.onc <- unique(all.sig.onc[duplicated(all.sig.onc)]) # 7
dup.sig.onc <- c(dup.sig.onc, "KRAS") # missed KRAS because only looked at things followed by underscore
dup.sig.onc.genes <- c("DCAF", "ESR1", "MAP2K") # for any proteins where gene symbol != protein name; ESR1 because of LTE2: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/LTE2_UP.V1_DN.html
dup.sig.onc.info <- unique(c(dup.sig.onc, dup.sig.onc.genes))
hall.info <- msigdbr(collection="H")
#kegg.info <- msigdbr(collection="C2", subcollection = "CP:KEGG_LEGACY")
wp.info <- msigdbr(collection="C2", subcollection = "CP:WIKIPATHWAYS")
rel.gs <- list() # 68 at end
dup.sig.onc.full <- c()
for (onc in dup.sig.onc.info) {
  if (onc == "IL2") {
    rel.hall <- unique(hall.info[hall.info$gene_symbol == onc | 
                                   grepl(onc, hall.info$gs_name, ignore.case=TRUE),]$gs_name)
    # rel.kegg <- unique(kegg.info[kegg.info$gene_symbol == onc | 
    #                                grepl(onc, kegg.info$gs_name, ignore.case=TRUE),]$gs_name)
    rel.wp <- unique(wp.info[wp.info$gene_symbol == onc | 
                                   grepl(onc, wp.info$gs_name, ignore.case=TRUE),]$gs_name)
  } else if (onc == "ESC") {
    rel.hall <- unique(hall.info[grepl("_ESC_", hall.info$gs_name, ignore.case=TRUE),]$gs_name) # else pick up E coli pathway
    #rel.kegg <- unique(kegg.info[grepl("_ESC_", kegg.info$gs_name, ignore.case=TRUE),]$gs_name)
    rel.wp <- unique(wp.info[grepl("_ESC_", wp.info$gs_name, ignore.case=TRUE),]$gs_name)
  } else {
    rel.hall <- unique(hall.info[startsWith(hall.info$gene_symbol, onc) | 
                                   grepl(onc, hall.info$gs_name, ignore.case=TRUE),]$gs_name)
    # rel.kegg <- unique(kegg.info[startsWith(kegg.info$gene_symbol, onc) | 
    #                                grepl(onc, kegg.info$gs_name, ignore.case=TRUE),]$gs_name) 
    rel.wp <- unique(wp.info[startsWith(wp.info$gene_symbol, onc) | 
                                   grepl(onc, wp.info$gs_name, ignore.case=TRUE),]$gs_name)
  }
  rel.gs[[onc]] <- c(rel.hall, rel.wp)
  dup.sig.onc.full <- unique(c(dup.sig.onc.full, all.gsea[all.gsea$Collection == "Oncogenic_signatures" & all.gsea$p_value <= 0.05 &
                                 all.gsea$FDR_q_value <= 0.25 & grepl(onc, all.gsea$Feature_set),]$Feature_set))
}
dup.sig.onc.full <- dup.sig.onc.full[!grepl("IL21", dup.sig.onc.full)] # IL21 accidentally included
all.rel.gs <- unlist(rel.gs, use.names=TRUE) # 121 MAP2K pathways, 93 MTOR, 58 KRAS, 37 IL2, 21 ESR1, 6 each of DCA/DCAF, 2 RPS14
all.rel.gs2 <- unique(unlist(rel.gs))
sig.rel.gs <- unique(gsea.rna.prot[gsea.rna.prot$Feature_set %in% all.rel.gs2 & gsea.rna.prot$Significant,]$Feature_set) # 32
sig.rel.gs.df <- gsea.rna.prot[gsea.rna.prot$Feature_set %in% all.rel.gs2 & gsea.rna.prot$Significant,] # 32
#sig.rel.gs.short <- sub("KEGG_", "", sig.rel.gs)
sig.rel.gs.short <- sub("HALLMARK_", "", sig.rel.gs)
sig.rel.gs.short <- sub("WP_", "", sig.rel.gs.short)

#mean.gsea[mean.gsea$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
#gsea.rna.prot[gsea.rna.prot$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
maxAbsNES <- max(abs(gsea.rna.prot[gsea.rna.prot$Collection %in% collections,]$NES))
gsea.rna.prot$minusLogFDR <- 4
gsea.rna.prot[gsea.rna.prot$FDR_q_value != 0,]$minusLogFDR <- -log10(gsea.rna.prot[gsea.rna.prot$FDR_q_value!=0,]$FDR_q_value)
maxLogFDR <- max(gsea.rna.prot[gsea.rna.prot$Collection %in% collections,]$minusLogFDR)
# library(patchwork); library(ggplot2)
# for (i in collections) {
#   nGenes <- length(unique(gsea.rna.prot[gsea.rna.prot$Collection == i,]$Feature_set))
#   topGenes <- unique(gsea.rna.prot[gsea.rna.prot$Significant &
#                                              gsea.rna.prot$Collection == i,]$Feature_set)
#   nTopGenes <- length(topGenes)
#   temp.mean.gsea <- mean.gsea[mean.gsea$Feature_set %in% topGenes,]
#   
#   # let's cap it at 20 pathways per plot
#   temp.mean.gsea <- temp.mean.gsea %>% slice_max(abs(mean_NES), n = 25)
#   temp.mean.gsea$Feature_set <- sub("WP_", "", temp.mean.gsea$Feature_set)
#   temp.mean.gsea$Feature_set <- sub("HALLMARK_", "", temp.mean.gsea$Feature_set)
#   geneOrder <- temp.mean.gsea[order(temp.mean.gsea$mean_NES),]$Feature_set
#   # geneFace <- rep("plain", length(geneOrder))
#   # names(geneFace) <- geneOrder
#   # sigInRNA <- temp.mean.gsea[grepl("RNA", temp.mean.gsea$sig_types),]$Feature_set
#   # geneFace[sigInRNA] <- "italic"
#   # sigInBoth <- temp.mean.gsea[grepl("RNA, Protein", temp.mean.gsea$sig_types),]$Feature_set
#   # geneFace[sigInBoth] <- "bold"
#   # if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
#   
#   geneFace <- rep("plain", length(geneOrder))
#   names(geneFace) <- geneOrder
#   poi <- c(geneOrder[geneOrder %in% sig.rel.gs.short], dup.sig.onc.full)
#            
#   geneFace[poi] <- "bold"
#   
#   temp.dot.df <- all.gsea[all.gsea$Feature_set %in% topGenes &
#                        all.gsea$Omics %in% c("RNA", "Protein"),]
#   temp.dot.df$Significant <- FALSE
#   temp.dot.df[temp.dot.df$p_value <= 0.05 & temp.dot.df$FDR_q_value <= 0.25,]$Significant <- TRUE
#   if (any(temp.dot.df$FDR_q_value == 0)) {
#     temp.dot.df[temp.dot.df$FDR_q_value == 0,]$FDR_q_value <- 0.0001
#   }
#   temp.dot.df$minusLogP <- -log(temp.dot.df[,"p_value"], base = 10)
#   temp.dot.df$minusLogFDR <- -log(temp.dot.df[,"FDR_q_value"], base = 10)
#   temp.dot.df$Omics <- factor(temp.dot.df$Omics, levels=c("RNA", "Protein"))
#   temp.dot.df$Feature_set <- sub("WP_", "", temp.dot.df$Feature_set)
#   temp.dot.df$Feature_set <- sub("HALLMARK_", "", temp.dot.df$Feature_set)
#   #temp.dot.df[temp.dot.df$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
#   dot.plot <- ggplot2::ggplot(
#     na.omit(temp.dot.df),
#     ggplot2::aes(
#       x = Omics, y = Feature_set, color = NES,
#       size = minusLogFDR
#     )
#   ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
#     ggplot2::geom_point() +
#     ggplot2::scale_y_discrete(limits = geneOrder) +
#     scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
#     #viridis::scale_color_viridis() +
#     theme_classic() +
#     ggplot2::labs(
#       x = "Omics Type",
#       y = "Gene Set",
#       color = "NES", size = "-log(FDR)"
#     ) + 
#     theme(axis.text.y=element_text(face=geneFace)) +
#     geom_point(data = subset(temp.dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1, vjust=1))
#   plot.annot <- paste0(i, "\n(", nTopGenes, " / ", nGenes, " gene sets enriched)")
#   dot.plot <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
#   dot.plot
#   if (i == collections[1]) {
#     gsea.dot.plots <- (dot.plot + theme(legend.position = "none"))
#   } else if (i == collections[2]) {
#     gsea.dot.plots <- gsea.dot.plots / dot.plot
#   } else {
#     gsea.dot.plots <- gsea.dot.plots / (dot.plot + theme(legend.position = "none"))
#   }
#   # most are only in prot and not in RNA
# }
# gsea.dot.plots
# ggplot2::ggsave("Enriched_geneSets_dotPlot_patchworkCollection_boldOnc_boldTitles_25max_2025-04-18.pdf", gsea.dot.plots, width=9, height=13)
# 
# # what if we only look at things sig in both rna and protein?
# for (i in collections) {
#   nGenes <- length(unique(gsea.rna.prot[gsea.rna.prot$Collection == i,]$Feature_set))
#   topGenes <- unique(gsea.rna.prot[gsea.rna.prot$Significant &
#                                      gsea.rna.prot$Collection == i,]$Feature_set)
#   nTopGenes <- length(topGenes)
#   temp.mean.gsea <- mean.gsea[mean.gsea$Feature_set %in% topGenes,]
#   
#   temp.mean.gsea <- temp.mean.gsea[temp.mean.gsea$sig_types == "RNA, Protein",]
#   
#   # let's cap it at 20 pathways per plot
#   temp.mean.gsea <- temp.mean.gsea %>% slice_max(abs(mean_NES), n = 25)
#   temp.mean.gsea$Feature_set <- sub("WP_", "", temp.mean.gsea$Feature_set)
#   temp.mean.gsea$Feature_set <- sub("HALLMARK_", "", temp.mean.gsea$Feature_set)
#   geneOrder <- temp.mean.gsea[order(temp.mean.gsea$mean_NES),]$Feature_set
#   # geneFace <- rep("plain", length(geneOrder))
#   # names(geneFace) <- geneOrder
#   # sigInRNA <- temp.mean.gsea[grepl("RNA", temp.mean.gsea$sig_types),]$Feature_set
#   # geneFace[sigInRNA] <- "italic"
#   # sigInBoth <- temp.mean.gsea[grepl("RNA, Protein", temp.mean.gsea$sig_types),]$Feature_set
#   # geneFace[sigInBoth] <- "bold"
#   # if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
#   
#   geneFace <- rep("plain", length(geneOrder))
#   names(geneFace) <- geneOrder
#   poi <- c(geneOrder[geneOrder %in% sig.rel.gs.short], dup.sig.onc.full)
#   
#   geneFace[poi] <- "bold"
#   
#   temp.dot.df <- all.gsea[all.gsea$Feature_set %in% topGenes &
#                             all.gsea$Omics %in% c("RNA", "Protein"),]
#   temp.dot.df$Significant <- FALSE
#   temp.dot.df[temp.dot.df$p_value <= 0.05 & temp.dot.df$FDR_q_value <= 0.25,]$Significant <- TRUE
#   if (any(temp.dot.df$FDR_q_value == 0)) {
#     temp.dot.df[temp.dot.df$FDR_q_value == 0,]$FDR_q_value <- 0.0001
#   }
#   temp.dot.df$minusLogP <- -log(temp.dot.df[,"p_value"], base = 10)
#   temp.dot.df$minusLogFDR <- -log(temp.dot.df[,"FDR_q_value"], base = 10)
#   temp.dot.df$Omics <- factor(temp.dot.df$Omics, levels=c("RNA", "Protein"))
#   temp.dot.df$Feature_set <- sub("WP_", "", temp.dot.df$Feature_set)
#   temp.dot.df$Feature_set <- sub("HALLMARK_", "", temp.dot.df$Feature_set)
#   #temp.dot.df[temp.dot.df$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
#   dot.plot <- ggplot2::ggplot(
#     na.omit(temp.dot.df),
#     ggplot2::aes(
#       x = Omics, y = Feature_set, color = NES,
#       size = minusLogFDR
#     )
#   ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
#     ggplot2::geom_point() +
#     ggplot2::scale_y_discrete(limits = geneOrder) +
#     scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
#     #viridis::scale_color_viridis() +
#     theme_classic() +
#     ggplot2::labs(
#       x = "Omics Type",
#       y = "Gene Set",
#       color = "NES", size = "-log(FDR)"
#     ) + 
#     theme(axis.text.y=element_text(face=geneFace)) +
#     geom_point(data = subset(temp.dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1, vjust=1))
#   plot.annot <- paste0(i, "\n(", nTopGenes, " / ", nGenes, " gene sets enriched)")
#   dot.plot <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
#   dot.plot
#   if (i == collections[1]) {
#     gsea.dot.plots <- (dot.plot + theme(legend.position = "none"))
#   } else if (i == collections[2]) {
#     gsea.dot.plots <- gsea.dot.plots / dot.plot
#   } else {
#     gsea.dot.plots <- gsea.dot.plots / (dot.plot + theme(legend.position = "none"))
#   }
#   # most are only in prot and not in RNA
# }
# gsea.dot.plots
# ggplot2::ggsave("sharedEnriched_geneSets_dotPlot_patchworkCollection_boldOnc_boldTitles_25max_2025-04-18.pdf", gsea.dot.plots, width=6, height=6)
# 
# # what if we only look at things which involve the oncogenic hits
# for (i in collections) {
#   nGenes <- length(unique(gsea.rna.prot[gsea.rna.prot$Collection == i,]$Feature_set))
#   topGenes <- unique(gsea.rna.prot[gsea.rna.prot$Significant &
#                                      gsea.rna.prot$Collection == i,]$Feature_set)
#   nTopGenes <- length(topGenes)
#   temp.mean.gsea <- mean.gsea[mean.gsea$Feature_set %in% topGenes,]
#   
#   temp.mean.gsea <- temp.mean.gsea[temp.mean.gsea$Feature_set %in% sig.rel.gs,]
#   
#   # let's cap it at 20 pathways per plot
#   temp.mean.gsea <- temp.mean.gsea %>% slice_max(abs(mean_NES), n = 25)
#   temp.mean.gsea$Feature_set <- sub("WP_", "", temp.mean.gsea$Feature_set)
#   temp.mean.gsea$Feature_set <- sub("HALLMARK_", "", temp.mean.gsea$Feature_set)
#   geneOrder <- temp.mean.gsea[order(temp.mean.gsea$mean_NES),]$Feature_set
#   # geneFace <- rep("plain", length(geneOrder))
#   # names(geneFace) <- geneOrder
#   # sigInRNA <- temp.mean.gsea[grepl("RNA", temp.mean.gsea$sig_types),]$Feature_set
#   # geneFace[sigInRNA] <- "italic"
#   # sigInBoth <- temp.mean.gsea[grepl("RNA, Protein", temp.mean.gsea$sig_types),]$Feature_set
#   # geneFace[sigInBoth] <- "bold"
#   # if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
#   
#   geneFace <- rep("plain", length(geneOrder))
#   names(geneFace) <- geneOrder
#   poi <- c(geneOrder[geneOrder %in% sig.rel.gs.short], dup.sig.onc.full)
#   
#   geneFace[poi] <- "bold"
#   
#   temp.dot.df <- all.gsea[all.gsea$Feature_set %in% topGenes &
#                             all.gsea$Omics %in% c("RNA", "Protein"),]
#   temp.dot.df$Significant <- FALSE
#   temp.dot.df[temp.dot.df$p_value <= 0.05 & temp.dot.df$FDR_q_value <= 0.25,]$Significant <- TRUE
#   if (any(temp.dot.df$FDR_q_value == 0)) {
#     temp.dot.df[temp.dot.df$FDR_q_value == 0,]$FDR_q_value <- 0.0001
#   }
#   temp.dot.df$minusLogP <- -log(temp.dot.df[,"p_value"], base = 10)
#   temp.dot.df$minusLogFDR <- -log(temp.dot.df[,"FDR_q_value"], base = 10)
#   temp.dot.df$Omics <- factor(temp.dot.df$Omics, levels=c("RNA", "Protein"))
#   temp.dot.df$Feature_set <- sub("WP_", "", temp.dot.df$Feature_set)
#   temp.dot.df$Feature_set <- sub("HALLMARK_", "", temp.dot.df$Feature_set)
#   #temp.dot.df[temp.dot.df$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
#   dot.plot <- ggplot2::ggplot(
#     na.omit(temp.dot.df),
#     ggplot2::aes(
#       x = Omics, y = Feature_set, color = NES,
#       size = minusLogFDR
#     )
#   ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
#     ggplot2::geom_point() +
#     ggplot2::scale_y_discrete(limits = geneOrder) +
#     scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
#     #viridis::scale_color_viridis() +
#     theme_classic() +
#     ggplot2::labs(
#       x = "Omics Type",
#       y = "Gene Set",
#       color = "NES", size = "-log(FDR)"
#     ) + 
#     theme(axis.text.y=element_text(face=geneFace)) +
#     geom_point(data = subset(temp.dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1, vjust=1))
#   plot.annot <- paste0(i, "\n(", nTopGenes, " / ", nGenes, " gene sets enriched)")
#   dot.plot <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
#   dot.plot
#   if (i == collections[1]) {
#     gsea.dot.plots <- (dot.plot + theme(legend.position = "none"))
#   } else if (i == collections[2]) {
#     gsea.dot.plots <- gsea.dot.plots / dot.plot
#   } else {
#     gsea.dot.plots <- gsea.dot.plots / (dot.plot + theme(legend.position = "none"))
#   }
#   # most are only in prot and not in RNA
# }
# gsea.dot.plots
# ggplot2::ggsave("oncEnriched_geneSets_dotPlot_patchworkCollection_boldOnc_boldTitles_25max_2025-04-18.pdf", gsea.dot.plots, width=8, height=12)
# 
# # what if we use NES instead of mean NES to sort
# for (i in collections) {
#   nGenes <- length(unique(gsea.rna.prot[gsea.rna.prot$Collection == i,]$Feature_set))
#   topGenes <- unique(gsea.rna.prot[gsea.rna.prot$Significant &
#                                      gsea.rna.prot$Collection == i,]$Feature_set)
#   nTopGenes <- length(topGenes)
#   temp.gsea <- gsea.rna.prot[gsea.rna.prot$Feature_set %in% topGenes,]
#   
#   # let's cap it at 20 pathways per plot
#   temp.gsea <- temp.gsea %>% slice_max(abs(NES), n = 20)
#   temp.gsea$Feature_set <- sub("WP_", "", temp.gsea$Feature_set)
#   temp.gsea$Feature_set <- sub("HALLMARK_", "", temp.gsea$Feature_set)
#   geneOrder <- unique(temp.gsea[order(temp.gsea$NES),]$Feature_set)
#   # geneFace <- rep("plain", length(geneOrder))
#   # names(geneFace) <- geneOrder
#   # sigInRNA <- temp.mean.gsea[grepl("RNA", temp.mean.gsea$sig_types),]$Feature_set
#   # geneFace[sigInRNA] <- "italic"
#   # sigInBoth <- temp.mean.gsea[grepl("RNA, Protein", temp.mean.gsea$sig_types),]$Feature_set
#   # geneFace[sigInBoth] <- "bold"
#   # if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
#   
#   geneFace <- rep("plain", length(geneOrder))
#   names(geneFace) <- geneOrder
#   poi <- c(geneOrder[geneOrder %in% sig.rel.gs.short], dup.sig.onc.full)
#   
#   geneFace[poi] <- "bold"
#   
#   temp.dot.df <- all.gsea[all.gsea$Feature_set %in% topGenes &
#                             all.gsea$Omics %in% c("RNA", "Protein"),]
#   temp.dot.df$Significant <- FALSE
#   temp.dot.df[temp.dot.df$p_value <= 0.05 & temp.dot.df$FDR_q_value <= 0.25,]$Significant <- TRUE
#   if (any(temp.dot.df$FDR_q_value == 0)) {
#     temp.dot.df[temp.dot.df$FDR_q_value == 0,]$FDR_q_value <- 0.0001
#   }
#   temp.dot.df$minusLogP <- -log(temp.dot.df[,"p_value"], base = 10)
#   temp.dot.df$minusLogFDR <- -log(temp.dot.df[,"FDR_q_value"], base = 10)
#   temp.dot.df$Omics <- factor(temp.dot.df$Omics, levels=c("RNA", "Protein"))
#   temp.dot.df$Feature_set <- sub("WP_", "", temp.dot.df$Feature_set)
#   temp.dot.df$Feature_set <- sub("HALLMARK_", "", temp.dot.df$Feature_set)
#   #temp.dot.df[temp.dot.df$Collection == "Oncogenic_signatures",]$Collection <- "Oncogenic"
#   dot.plot <- ggplot2::ggplot(
#     na.omit(temp.dot.df),
#     ggplot2::aes(
#       x = Omics, y = Feature_set, color = NES,
#       size = minusLogFDR
#     )
#   ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
#     ggplot2::geom_point() +
#     ggplot2::scale_y_discrete(limits = geneOrder) +
#     scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
#     #viridis::scale_color_viridis() +
#     theme_classic() +
#     ggplot2::labs(
#       x = "Omics Type",
#       y = "Gene Set",
#       color = "NES", size = "-log(FDR)"
#     ) + 
#     theme(axis.text.y=element_text(face=geneFace)) +
#     geom_point(data = subset(temp.dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1, vjust=1))
#   plot.annot <- paste0(i, "\n(", nTopGenes, " / ", nGenes, " gene sets enriched)")
#   dot.plot <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
#   dot.plot
#   if (i == collections[1]) {
#     gsea.dot.plots <- (dot.plot + theme(legend.position = "none"))
#   } else if (i == collections[2]) {
#     gsea.dot.plots <- gsea.dot.plots / dot.plot
#   } else {
#     gsea.dot.plots <- gsea.dot.plots / (dot.plot + theme(legend.position = "none"))
#   }
#   # most are only in prot and not in RNA
# }
# gsea.dot.plots
# ggplot2::ggsave("Enriched_geneSets_dotPlot_patchworkCollection_boldOnc_boldTitles_20maxByNES_2025-04-18.pdf", gsea.dot.plots, width=8, height=11)

# what if we use NES instead of mean NES to sort - no bold oncogenes, sort by n gene sets
collections <- c("Hallmark", "Oncogenic", "WikiPathways")
allTopSets <- c()
n.cutoff <- c(10)
for (n in n.cutoff) {
  for (i in collections) {
    nGenes <- length(unique(gsea.rna.prot[gsea.rna.prot$Collection == i,]$Feature_set))
    topGenes <- unique(gsea.rna.prot[gsea.rna.prot$Significant &
                                       gsea.rna.prot$Collection == i,]$Feature_set)
    nTopGenes <- length(topGenes)
    temp.gsea <- gsea.rna.prot[gsea.rna.prot$Feature_set %in% topGenes,]
    
    # let's cap it at 20 pathways per plot
    temp.gsea <- temp.gsea %>% slice_max(abs(NES), n = n)
    allTopSets <- unique(c(allTopSets, temp.gsea$Feature_set))
    temp.gsea$Feature_set <- sub("WP_", "", temp.gsea$Feature_set)
    temp.gsea$Feature_set <- sub("HALLMARK_", "", temp.gsea$Feature_set)
    geneOrder <- unique(temp.gsea[order(temp.gsea$NES),]$Feature_set)
    
    temp.dot.df <- all.gsea[all.gsea$Feature_set %in% topGenes &
                              all.gsea$Omics %in% c("RNA", "Protein"),]
    temp.dot.df$Significant <- FALSE
    temp.dot.df[temp.dot.df$p_value <= 0.05 & temp.dot.df$FDR_q_value <= 0.25,]$Significant <- TRUE
    if (any(temp.dot.df$FDR_q_value == 0)) {
      temp.dot.df[temp.dot.df$FDR_q_value == 0,]$FDR_q_value <- 0.0001
    }
    temp.dot.df$minusLogP <- -log(temp.dot.df[,"p_value"], base = 10)
    temp.dot.df$minusLogFDR <- -log(temp.dot.df[,"FDR_q_value"], base = 10)
    temp.dot.df$Omics <- factor(temp.dot.df$Omics, levels=c("RNA", "Protein"))
    temp.dot.df$Feature_set <- sub("WP_", "", temp.dot.df$Feature_set)
    temp.dot.df$Feature_set <- sub("HALLMARK_", "", temp.dot.df$Feature_set)
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
      scale_color_gradient2(low="blue",high="red", mid="grey", limits=c(-maxAbsNES, maxAbsNES)) +
      #viridis::scale_color_viridis() +
      theme_classic() +
      ggplot2::labs(
        x = "Omics Type",
        y = "Gene Set",
        color = "NES", size = "-log(FDR)"
      ) + 
      geom_point(data = subset(temp.dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1, vjust=1))
    plot.annot <- paste0(i, "\n(", nTopGenes, " / ", nGenes, " gene sets enriched)")
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
  ggplot2::ggsave(paste0("Enriched_geneSets_dotPlot_patchworkCollection_",n,"maxByNES_2025-04-18_shorter.pdf"), gsea.dot.plots, width=7, height=9)
  ggplot2::ggsave(paste0("Enriched_geneSets_dotPlot_patchworkCollection_",n,"maxByNES_2025-04-18_taller.pdf"), gsea.dot.plots, width=9, height=14)
  ggplot2::ggsave(paste0("Enriched_geneSets_dotPlot_patchworkCollection_",n,"maxByNES_2025-04-18.pdf"), gsea.dot.plots, width=8, height=11)
}

#### Look closer at top onc pathway genes ####
corr.result <- na.omit(read.csv(synapser::synGet("syn66227486")$path))
corr.result$minusLogFDR <- -log10(corr.result$Spearman.q) # max numeric is 37.327406
corr.result[corr.result$Spearman.q == 0,]$minusLogFDR <- 40
corr.result$Gene <- corr.result$Feature
corr.result[corr.result$Omics=="Phospho",]$Gene <- sub("-.*","",corr.result[corr.result$Omics=="Phospho",]$Feature)
absMaxEst <- 1
maxLogFDR <- max(corr.result$minusLogFDR)
onc.dot.plots <- NULL
all.topGenes <- c()
wnt.gs <- gsea.rna.prot[grepl("WNT",gsea.rna.prot$Feature_set) & gsea.rna.prot$Significant,]$Feature_set # 6
my.onc <- unique(onc.info[onc.info$gs_name %in% allTopSets,]$gene_symbol)
myc.hallmark <- unique(hall.info[hall.info$gs_name %in% allTopSets,]$gene_symbol)
myc.wp <- unique(wp.info[wp.info$gs_name %in% allTopSets,]$gene_symbol)
myc.genes <- unique(c(myc.hallmark, myc.wp))#, onc.gene))
myc.corr <- corr.result[corr.result$Gene %in% myc.genes,]

topGenes <- unique(myc.corr[myc.corr$Significant,]$Feature) # 95

dot.df <- myc.corr[myc.corr$Gene %in% topGenes,] # 112

dot.df <- na.omit(corr.result[corr.result$Feature %in% topGenes &
                                corr.result$Omics %in% c("RNA", "Protein"),])

dot.df <- na.omit(corr.result[corr.result$Omics %in% c("RNA", "Protein", "Phospho"),])
dot.df$Gene <- dot.df$Feature
dot.df[dot.df$Omics=="Phospho",]$Gene <- sub("-.*","",dot.df[dot.df$Omics=="Phospho",]$Feature)
dot.df <- dot.df[dot.df$Gene %in% topGenes,]
if (nrow(dot.df) > 0) {
  #dot.df <- dot.df %>% slice_max(abs(Spearman.est), n=10)
  dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
  
  geneOrder <- dot.df[order(dot.df$Spearman.est),]$Feature
  geneFace <- rep("plain", length(geneOrder))
  names(geneFace) <- geneOrder
  #geneFace[onc.gene] <- "bold"
  #temp.dup.topGenes <- geneOrder[geneOrder %in% dup.topGenes]
  #if (length(temp.dup.topGenes) > 0) {geneFace[temp.dup.topGenes] <- "bold"}
  
  myc.dot.plot <- ggplot2::ggplot(
    dot.df,
    ggplot2::aes(
      x = Omics, y = Feature, color = Spearman.est,
      size = minusLogFDR
    )
  ) + ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = geneOrder) +
    scale_color_gradient2(low="blue",high="red", mid="grey", limits = c(-absMaxEst, absMaxEst)) +
    scale_size(limits=c(0,maxLogFDR), range = c(0.5,5)) +
    theme_classic() +
    ggplot2::labs(
      x = "Omics Type",
      y = "Gene",
      color = "Spearman rho", size = "-log(adj. p)"
    ) + geom_point(data = subset(dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + ggtitle("Correlated") +
    theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
    theme(axis.text.y=element_text(face=geneFace), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
}
dup.topGenes <- all.topGenes[duplicated(all.topGenes)]
myc.dot.plot
#ggplot2::ggsave("Correlated_MAP2K-KRAS_genes_dotPlot.pdf", onc.dot.plots, width=3, height=10)
ggplot2::ggsave("Correlated_top10Pathway_genes_dotPlot_10max.pdf", onc.dot.plots, width=3, height=13)

wnt.gs <- gsea.rna.prot[grepl("WNT",gsea.rna.prot$Feature_set) & gsea.rna.prot$Significant,]$Feature_set # 6 - only oncogenic and WP, no hallmark
onc.info <- msigdbr::msigdbr(collection="C6")
myc.onc <- unique(onc.info[onc.info$gs_name %in% wnt.gs,]$gene_symbol)
myc.wp <- unique(wp.info[wp.info$gs_name %in% wnt.gs,]$gene_symbol)
myc.genes <- unique(c(myc.onc, myc.wp))#, onc.gene))

myc.corr <- corr.result[corr.result$Feature %in% myc.genes,] # 530
myc.corr <- myc.corr[myc.corr$Omics %in% c("RNA", "Protein"),] # 268
topGenes <- unique(myc.corr[myc.corr$Significant,]$Feature) # 13
all.topGenes <- c(all.topGenes, topGenes)

dot.df <- na.omit(corr.result[corr.result$Feature %in% topGenes &
                                corr.result$Omics %in% c("RNA", "Protein"),])
dot.df <- na.omit(corr.result[corr.result$Omics %in% c("RNA", "Protein", "Phospho"),])
dot.df$Gene <- dot.df$Feature
dot.df[dot.df$Omics=="Phospho",]$Gene <- sub("-.*","",dot.df[dot.df$Omics=="Phospho",]$Feature)
dot.df <- dot.df[dot.df$Gene %in% topGenes,]
dot.df$minusLogFDR <- -log10(dot.df$Spearman.q)
dot.df[dot.df$Spearman.q == 0,]$minusLogFDR <- ceiling(max(dot.df[dot.df$Spearman.q != 0,]$minusLogFDR))
if (nrow(dot.df) > 0) {
  dot.df <- dot.df %>% slice_max(abs(Spearman.est), n=10)
  dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
  
  mean.corr <- plyr::ddply(myc.corr[myc.corr$Feature %in% dot.df$Feature,], .(Feature), summarize,
                           mean_rankVal = mean(Spearman.est, na.rm = TRUE),
                           sig_types = paste(Omics, collapse = ", "))
  geneOrder <- mean.corr[order(mean.corr$mean_rankVal),]$Feature
  geneFace <- rep("plain", length(geneOrder))
  names(geneFace) <- geneOrder
  #geneFace[onc.gene] <- "bold"
  #temp.dup.topGenes <- geneOrder[geneOrder %in% dup.topGenes]
  #if (length(temp.dup.topGenes) > 0) {geneFace[temp.dup.topGenes] <- "bold"}
  
  myc.dot.plot <- ggplot2::ggplot(
    dot.df,
    ggplot2::aes(
      x = Omics, y = Feature, color = Spearman.est,
      size = minusLogFDR
    )
  ) + ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = geneOrder) +
    scale_color_gradient2(low="blue",high="red", mid="grey", limits = c(-absMaxEst, absMaxEst)) +
    scale_size(limits=c(0,maxLogFDR), range = c(0.5,5)) +
    theme_classic() +
    ggplot2::labs(
      x = "Omics Type",
      y = "Gene",
      color = "Spearman rho", size = "-log(adj. p)"
    ) + geom_point(data = subset(dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + ggtitle("WNT") +
    theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
    theme(axis.text.y=element_text(face=geneFace), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
}
#ggplot2::ggsave("Correlated_MAP2K-KRAS_genes_dotPlot.pdf", onc.dot.plots, width=3, height=10)
ggplot2::ggsave("Correlated_WNT_genes_dotPlot_10max.pdf", myc.dot.plot, width=2.5, height=2)
