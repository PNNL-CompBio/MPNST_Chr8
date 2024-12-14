# chr8 MPNST: figure 6
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_6")

synapser::synLogin()
synID <- synapser::synStore(synapser::Folder("Figure_6", parent = "syn64367769"))

#### 1. bar plot of # of corr drugs ####
drug.corr <- read.csv(synapser::synGet("syn64025324")$path)
mean.drugs <- read.csv(synapser::synGet("syn64025325")$path)
moa.results <- read.csv(synapser::synGet("syn64025318")$path)
mean.moa <- read.csv(synapser::synGet("syn64025319")$path)
gsea.rna.prot <- moa.results[moa.results$type %in% c("RNA", "Protein"),]
gsea.rna.prot$Omics <- factor(gsea.rna.prot$Omics, levels=c("RNA", "Protein"))
gsea.rna.prot$Direction <- factor(gsea.rna.prot$Direction, levels=c("Positive", "Negative"))
gsea.rna.prot$Toxicity <- "Chr8q-amplified"
gsea.rna.prot[gsea.rna.prot$NES > 0,]$Toxicity <- "Non-amplified"
bar.plot <- ggplot2::ggplot(gsea.rna.prot[gsea.rna.prot$Significant,], aes(fill = Toxicity, x=Omics)) + 
  geom_bar(position = "dodge", stat="count", alpha=0.75) + ggplot2::theme_classic() + 
  scale_fill_manual(values=c("red", "blue")) + 
  ylim(c(0,15)) +
  ggplot2::labs(
    x = "Omics Type",
    y = "Number of Enriched Drug Mechanisms of Action",
    fill = "Potential Toxicity",
    color = "Correlation") +
  theme(axis.title.x=element_blank()) #+
#theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
bar.plot
ggplot2::ggsave("Chr8_numberOfGeneSets_directionGroup_barPlot_manualOrderV3.pdf", bar.plot, width = 7, height = 7, device = "pdf")

# source: https://r-graph-gallery.com/web-circular-barplot-with-R-and-ggplot2.html
plot_df <- plyr::ddply(gsea.rna.prot[gsea.rna.prot$Significant,] , .(Toxicity, Omics), summarize,
                       nEnriched = n())
plot_df$Toxicity <- factor(plot_df$Toxicity, levels=c("Non-amplified", "Chr8q-amplified"))
plt <- ggplot(plot_df) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(5, 10)),
    color = "lightgrey"
  ) + 
  # Add bars to represent the cumulative track lengths
  # str_wrap(region, 5) wraps the text so each line has at most 5 characters
  # (but it doesn't break long words!)
  geom_col(
    aes(
      x = reorder(Omics, nEnriched),
      #x = forcats::fct_infreq(stringr::str_wrap(Collection, 5)),
      y = nEnriched,
      fill = Omics, alpha = Toxicity
    ),
    position = "dodge", show.legend = TRUE#, alpha=0.5
  ) +
  geom_segment(
    aes(
      x = reorder(Omics, nEnriched),
      y = 0,
      xend = reorder(Omics, nEnriched),
      yend = max(nEnriched)
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
  # Scale y axis so bars don't start in the center
  scale_y_continuous(
    limits = c(-5, 10),
    expand = c(0, 0),
    breaks = c(0, 5, 10)
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
  ) + ggtitle("Number of Enriched Drug Mechanisms") + 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
        axis.line=element_blank())

plt
ggplot2::ggsave("Chr8_numberOfDrugSets_circularBarPlot_greenOrange_matchingVenn.pdf", plt, width = 7, height = 7)
#synapser::synStore(synapser::File("Chr8_numberOfDrugs_directionGroup_barPlot_logScale_manualOrderV2.pdf", parent = synID))

#### 2. venn diagram of diffexp features ####
venn.list <- list()
for (i in 1:length(c("RNA", "Protein"))) {
  temp.name <- c("RNA", "Protein")[i]
  venn.list[[temp.name]] <- unique(gsea.rna.prot[gsea.rna.prot$Omics == temp.name &
                                                   gsea.rna.prot$Significant,]$Drug_set)
}

ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_color=c("#F8766D", "#00BFC4"), fill_alpha = 0.5, text_size = 10, set_name_size = 7.5)
ggplot2::ggsave(paste0("Chr8_GSEA_vennDiagram_noPercentages_", Sys.Date(), "_orangeGreen.pdf"), width=7, height=7)

#### 3. top up/dn diffexp gene sets ####
dot.df <- gsea.rna.prot[gsea.rna.prot$Significant,] # 19 - note glucocorticoid receptor agonist has FDR=0.25 in protein so use 'Significant' to include it
nrow(dot.df[dot.df$Omics=="RNA",]) # 7
nrow(dot.df[dot.df$Omics=="Protein",]) # 12
hist(dot.df$NES,breaks=20)
min(abs(dot.df$NES)) # 1.48
topGenes <- unique(dot.df$Drug_set) # 17

mean.moa <- mean.moa[mean.moa$Drug_set %in% topGenes,]
geneOrder <- mean.moa[order(mean.moa$mean_NES),]$Drug_set
# if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
moi <- na.omit(unique(mean.moa[grepl("RNA, Protein", mean.moa$sig_types),]$Drug_set)) # HDAC inhibitor
#moi <- c("HDAC inhibitor", "glucocorticoid receptor agonist")
geneFace <- rep("plain", length(geneOrder))
names(geneFace) <- geneOrder
geneFace[moi] <- "bold"
dot.df <- na.omit(moa.results[moa.results$Drug_set %in% topGenes &
                     moa.results$Omics %in% c("RNA", "Protein"),])
if (any(dot.df$FDR_q_value == 0)) {
  dot.df[dot.df$FDR_q_value == 0,]$FDR_q_value <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"p_value"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"FDR_q_value"], base = 10)
dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
absMaxNES <- max(abs(dot.df$NES))
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = Omics, y = Drug_set, color = NES,
    size = minusLogFDR
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="blue",high="red", limits=c(-absMaxNES, absMaxNES)) +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Drug Mechanism",
    color = "NES", size = "-log(FDR)"
  ) + geom_point(data = subset(dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + 
  ggtitle(paste0("Drug Mechanisms (", length(topGenes), " / ", length(unique(moa.results$Drug_set)), " enriched)")) +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
  theme(axis.text.y=element_text(face=geneFace))
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Enriched_drugSets_dotPlot_boldHDACi.pdf", dot.plot, width=5.5, height=4.5)

#### dot plot of corr drugs ####
hdac.drugs <- unique(drug.info[grepl("HDAC inhibitor", drug.info$moa, ignore.case = TRUE),]$name) # 25
hdac.drugs.noSpecialChar <- sub("-",".",hdac.drugs)

omics <- c("RNA", "Protein")
maxAbsEst <- max(na.omit(abs(drug.corr[drug.corr$Omics %in% omics,]$Pearson.est)))
#maxAbsEst <- 1
drug.corr <- na.omit(drug.corr)
drug.corr$`HDAC Inhibitor` <- FALSE
drug.corr[drug.corr$feature %in% hdac.drugs.noSpecialChar,]$`HDAC Inhibitor` <- TRUE
if (any(drug.corr$Pearson.q == 0)) {
  drug.corr[drug.corr$Pearson.q == 0,]$Pearson.q <- 0.0001
}
drug.corr$minusLogP <- -log(drug.corr[,"Pearson.p"], base = 10)
drug.corr$minusLogFDR <- -log(drug.corr[,"Pearson.q"], base = 10)
maxLogFDR <- max(na.omit(drug.corr[drug.corr$Omics %in% omics,]$minusLogFDR))
library(patchwork); library(ggplot2)
gsea.dot.plots <- NULL
for (i in omics) {
  nGenes <- length(unique(na.omit(drug.corr[drug.corr$Omics == i,])$feature))
  nCorrGenes <- length(unique(na.omit(drug.corr[drug.corr$Significant & 
                                                 drug.corr$Omics == i,])$feature))
  
  temp.dot.df <- drug.corr[drug.corr$Significant & drug.corr$Omics == i,] %>% slice_max(abs(Pearson.est), n = 50)
  topGenes <- temp.dot.df$feature
  nTopGenes <- length(topGenes)
  geneOrder <- temp.dot.df[order(temp.dot.df$Pearson.est),]$feature
  plot.annot <- paste0(i, "\n(", nCorrGenes, " / ", nGenes, " correlated)")
  
  geneFace <- rep("plain", length(geneOrder))
  names(geneFace) <- geneOrder
  geneFace[hdac.drugs.noSpecialChar] <- "bold"
  
  dot.plot <- ggplot2::ggplot(temp.dot.df,
                              ggplot2::aes(
                                x = Omics, y = feature, color = Pearson.est,
                                size = minusLogFDR
                              )
  ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = geneOrder) +
    scale_color_gradient2(low="blue",high="red", limits=c(-maxAbsEst, maxAbsEst)) +
    theme_classic() +
    ggplot2::labs(
      x = as.character(i),
      color = "Pearson r", size = "-log(adj. p)"
    ) + 
    theme(axis.text.x=element_text(face=geneFace)) +
    geom_point(data = subset(temp.dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
    #ylab(plot.annot)+
    theme(axis.title.y=element_text(face="bold", size=16#, angle = 90, vjust=1, hjust=1
    ), 
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x=element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
  
  cat(plot.annot, "\n")
  cat(i, "with", nTopGenes, "top features\n")
  dot.plot <- dot.plot + coord_flip() #+ ylab(plot.annot)#+ ggtitle(plot.annot) + 
  #theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
  dot.plot
  if (is.null(gsea.dot.plots)) {
    gsea.dot.plots <- dot.plot
  } else {
    gsea.dot.plots <- gsea.dot.plots / dot.plot
  } 
}
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4/guides_build_mod.R")
gsea.dot.plots <- gsea.dot.plots + plot_layout(guides = 'collect') + plot_annotation(caption = "Key\nBold: HDAC inhibitor")
gsea.dot.plots
#ggplot2::ggsave("Enriched_geneSets_dotPlot_patchworkOmics_boldMYCandSHH.pdf", gsea.dot.plots, width=9, height=9)
ggplot2::ggsave("CorrelatedDrugs_dotPlot_patchworkOmics_sliceMaxAbsPearson50_v2_boldHDACi.pdf", gsea.dot.plots, width=10, height=5)

# 316 sig drugs w Spearman vs 264 with Pearson

# color axis text by MOA?

#### waterfall plot of drug corr ####
drug.info <- read.csv("https://raw.githubusercontent.com/BelindaBGarana/DMEA/refs/heads/shiny-app/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv")
hdac.drugs <- unique(drug.info[grepl("HDAC inhibitor", drug.info$moa, ignore.case = TRUE),]$name) # 25
hdac.drugs.noSpecialChar <- sub("-",".",hdac.drugs)
omics <- c("RNA", "Protein")
maxAbsEst <- max(na.omit(abs(drug.corr[drug.corr$Omics %in% omics,]$Pearson.est)))
#maxAbsEst <- 1
drug.corr <- na.omit(drug.corr)
if (any(drug.corr$Pearson.q == 0)) {
  drug.corr[drug.corr$Pearson.q == 0,]$Pearson.q <- 0.0001
}
drug.corr$minusLogP <- -log(drug.corr[,"Pearson.p"], base = 10)
drug.corr$minusLogFDR <- -log(drug.corr[,"Pearson.q"], base = 10)
drug.info$Drug <- gsub("-",".",drug.info$name)
drug.info$Drug <- gsub("[(]",".",drug.info$Drug)
drug.info$Drug <- gsub("[)]",".",drug.info$Drug)
drug.info$Drug <- gsub("[+]",".",drug.info$Drug)
red.drug.info <- dplyr::distinct(drug.info[,c("Drug","name","moa","target")])
colnames(drug.corr)[2] <- "Drug"
drug.corr.wInfo <- merge(red.drug.info, drug.corr, by="Drug", all.y = TRUE)
drug.corr.wInfo$Mechanism <- "Other"
sig.moas <- unique(moa.results[moa.results$type %in% omics & moa.results$Significant,]$Drug_set)
drug.corr.wInfo[drug.corr.wInfo$moa %in% sig.moas,]$Mechanism <- drug.corr.wInfo[drug.corr.wInfo$moa %in% sig.moas,]$moa

maxAbsEst <- max(na.omit(abs(drug.corr.wInfo[drug.corr.wInfo$Omics %in% omics,]$Pearson.est)))
maxLogFDR <- max(na.omit(drug.corr.wInfo[drug.corr.wInfo$Omics %in% omics,]$minusLogFDR))
library(patchwork); library(ggplot2)
gsea.dot.plots <- NULL
moaOrder <- mean.moa[order(mean.moa$mean_NES),]$Drug_set
#drug.corr.wInfo$Mechanism <- factor(drug.corr.wInfo$Mechanism, levels = c(moaOrder, "Other"))
moaOrder2 <- moaOrder[moaOrder%in% unique(drug.corr.wInfo[drug.corr.wInfo$Significant,]$Mechanism)]
# HSP inhibitor (5th moa) isn't showing up on plots but others are
moaOrder3 <- moaOrder2[moaOrder2 != "HSP inhibitor"]

# take care of cases where there are multiple mechanisms
for (i in 1:length(moaOrder3)) {
  if (nrow(drug.corr.wInfo[drug.corr.wInfo$Mechanism == "Other" & 
                           grepl(moaOrder3[i],drug.corr.wInfo$moa),]) > 0) {
    drug.corr.wInfo[drug.corr.wInfo$Mechanism == "Other" & 
                      grepl(moaOrder3[i],drug.corr.wInfo$moa),]$Mechanism <- moaOrder3[i]
  }
}

library(plyr);library(dplyr)
for (i in omics) {
  nGenes <- length(unique(na.omit(drug.corr.wInfo[drug.corr.wInfo$Omics == i,])$Drug))
  nCorrGenes <- length(unique(na.omit(drug.corr.wInfo[drug.corr.wInfo$Significant & 
                                                  drug.corr.wInfo$Omics == i,])$Drug))
  
  temp.dot.df <- drug.corr.wInfo[drug.corr.wInfo$Significant & drug.corr.wInfo$Omics == i,]
  temp.dot.df <- temp.dot.df %>% slice_max(abs(Pearson.est),n=50)
  topGenes <- temp.dot.df$name
  nTopGenes <- length(topGenes)
  geneOrder <- temp.dot.df[order(temp.dot.df$Pearson.est),]$name
  
  geneFace <- rep("plain", length(geneOrder))
  names(geneFace) <- geneOrder
  geneFace[hdac.drugs.noSpecialChar] <- "bold"
  
  plot.annot <- paste0(i, "\n(", nCorrGenes, " / ", nGenes, " drugs correlated)")
  
  # geneFace <- rep("plain", length(geneOrder))
  # names(geneFace) <- geneOrder
  # MYC.pathways <- geneOrder[grepl("MYC", geneOrder)]
  # SHH.pathways <- geneOrder[grepl("SHH", geneOrder) | grepl("HEDGEHOG", geneOrder)]
  # poi <- c(MYC.pathways, SHH.pathways)
  # geneFace[poi] <- "bold"
  
  dot.plot <- ggplot2::ggplot(temp.dot.df,
                              ggplot2::aes(
                                x = name, y = Pearson.est, fill = Mechanism
                              )
  ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_discrete(limits = geneOrder) +
    theme_classic() + scale_fill_manual(breaks=c(moaOrder3,"Other"), 
                                          values = c(RColorBrewer::brewer.pal(12, "Set3"),"darkgrey")) +
    ggplot2::labs(
      #x = i,
      y = "Pearson r",
      fill = "Drug Mechanism"
    ) + theme(axis.title.x = element_blank()) +
    #theme(axis.text.y=element_text(face=geneFace)) +
    #theme(axis.title.y=element_text(face="bold", size=16), 
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    #axis.title.x=element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
  
  cat(plot.annot, "\n")
  cat(i, "with", nTopGenes, "top Drugs\n")
  dot.plot <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
  dot.plot
  if (is.null(gsea.dot.plots)) {
    gsea.dot.plots <- dot.plot
  } else {
    gsea.dot.plots <- gsea.dot.plots / dot.plot
  } 
}
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4/guides_build_mod.R")
gsea.dot.plots <- gsea.dot.plots + plot_layout(guides = 'collect')# + plot_annotation(caption = "Key\nBold: HDAC inhibitor")
gsea.dot.plots
#ggplot2::ggsave("Enriched_geneSets_dotPlot_patchworkOmics_boldMYCandSHH.pdf", gsea.dot.plots, width=9, height=9)
ggplot2::ggsave("CorrelatedDrugs_barPlot_patchworkOmics_moaFill_sliceMaxAbsPearson50_darkGreyOther.pdf", gsea.dot.plots, width=11, height=5)
