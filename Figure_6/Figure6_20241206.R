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

#### compile target ranks across all analyses ####
# get correlated genes (including from phospho)
corr.result <- read.csv(synapser::synGet("syn64024300")$path)
corr.result <- corr.result[corr.result$type != "DNA" & corr.result$Significant,]
corr.result$Gene <- corr.result$feature
corr.result[corr.result$type=="Phospho",]$Gene <- sub("-.*","",corr.result[corr.result$type=="Phospho",]$feature)
# corr.result <- plyr::ddply(corr.result, .(type, feature), summarize,
#                            Spearman.est = mean(Spearman.est, na.rm = TRUE),
#                            sd.Spearman.est = sd(Spearman.est, na.rm = TRUE),
#                            Spearman.p = mean(Spearman.p, na.rm = TRUE),
#                            sd.Spearman.p = sd(Spearman.p, na.rm = TRUE),
#                            Spearman.q = mean(Spearman.q, na.rm = TRUE),
#                            sd.Spearman.q = sd(Spearman.q, na.rm = TRUE),
#                            Significant = all(Significant)) # all SDs were NA so no duplicates
corr.result$Rank <- rank(corr.result$Spearman.est)
bestRank <- max(corr.result$Rank)
bestRankCorr <- max(corr.result$Rank)
corr.result$BestRank <- FALSE
corr.result[corr.result$Rank == bestRank,]$BestRank <- TRUE
corr.result$percentile <- percent_rank(corr.result$Spearman.est)

# get enriched pathways and their members
kin <- read.table(synapser::synGet("syn64374349")$path, sep="\t", header=TRUE)
kin$Gene <- kin$kinase
kin <- kin[kin$adjusted_p_value <= 0.25 & kin$p_value <= 0.05,] # 12
kin[kin$kinase == "P38G",]$Gene <- "MAPK12" # could also be ERK6; source: https://lincs.hms.harvard.edu/db/proteins/200777/
kin[kin$kinase == "ERK1",]$Gene <- "MAPK3" # could also be PRKM3; source: https://www.uniprot.org/uniprotkb/P27361/entry
# LTK should have LTK as gene name but could also be TYK1; source: https://www.uniprot.org/uniprotkb/P29376/entry
kin[kin$kinase == "P38B",]$Gene <- "MAPK11" # source: https://www.uniprot.org/uniprotkb/Q15759/entry
# CDK5 should be CDK5, same for CDK8, 10, 12, 19; MTOR should be MTOR
kin[kin$kinase == "P38A",]$Gene <- "MAPK14" # source: https://www.uniprot.org/uniprotkb/Q16539/entry
kin[kin$kinase == "JNK3",]$Gene <- "MAPK10" # source: https://www.uniprot.org/uniprotkb/P53779/entry
site.mapping <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4/site_mapping.csv")
site.mapping$Gene <- sub("-.*","",site.mapping$SUB_SITE)
site.mapping2 <- plyr::ddply(site.mapping, .(kinase), summarize,
                        Substrates = paste0(unique(Gene), collapse = " "))
kin2 <- merge(kin, site.mapping2, by="kinase", all.x = TRUE)
kin2$Rank <- rank(kin2$enrichment_value_log2)
bestRank <- max(kin2$Rank)
bestRankKin <- max(kin2$Rank)
kin2$BestRank <- FALSE
kin2[kin2$Rank == bestRank,]$BestRank <- TRUE
kin2$percentile <- percent_rank(kin2$enrichment_value_log2)

all.gsea <- read.csv(synapser::synGet("syn64025221")$path)
tf <- all.gsea[all.gsea$Collection=="TFT_GTRD" & all.gsea$type=="RNA" & all.gsea$Significant,] # 12
tf$Gene <- sub("_TARGET_GENES.*","",tf$Feature_set)
tf.genes <- msigdbr::msigdbr(species="Homo sapiens",category="C3", subcategory="GTRD")
tf.genes <- tf.genes[tf.genes$gs_name %in% tf$Feature_set,] # 1871 rows
tf.genes <- dplyr::distinct(tf.genes[tf.genes$gs_name %in% tf$Feature_set,c("gs_name","gene_symbol")]) # 1869 rows
colnames(tf.genes) <- c("Feature_set","Target")
tf.genes <- plyr::ddply(tf.genes, .(Feature_set), summarize,
                        Targets = paste0(unique(Target), collapse = " "))
tf$Rank <- rank(tf$NES)
bestRank <- max(tf$Rank)
bestRankTF <- max(tf2$Rank)
tf$BestRank <- FALSE
tf[tf$Rank == bestRank,]$BestRank <- TRUE
tf2 <- merge(tf, tf.genes, by="Feature_set", all.x = TRUE)
tf2$percentile <- percent_rank(tf2$NES)

# GSEA
coi <- c("Hallmark", "Oncogenic_signatures", "PID")
all.gsea <- all.gsea[all.gsea$type %in% c("RNA", "Protein") &
                       all.gsea$Collection %in% coi & all.gsea$Significant,] # 224
collections <- unique(all.gsea$Collection) # same as coi
hallmark <- msigdbr::msigdbr(species="Homo sapiens",category="H")
hallmark <- dplyr::distinct(hallmark[hallmark$gs_name %in% all.gsea$Feature_set,c("gs_name","gene_symbol")]) # 2358 rows
colnames(hallmark) <- c("Feature_set","Gene")
hallmark <- plyr::ddply(hallmark, .(Feature_set), summarize,
                        Gene = paste0(unique(Gene), collapse = " "))
all.gsea2 <- merge(hallmark, all.gsea, by="Feature_set", all.y = TRUE)

pid <- msigdbr::msigdbr(species="Homo sapiens",category="C2", subcategory="PID")
pid <- dplyr::distinct(pid[pid$gs_name %in% all.gsea$Feature_set,c("gs_name","gene_symbol")]) # 649 rows
colnames(pid) <- c("Feature_set","Gene")
pid <- plyr::ddply(pid, .(Feature_set), summarize,
                        Gene = paste0(unique(Gene), collapse = " "))
all.gsea2 <- merge(pid, all.gsea2, by="Feature_set", all.y = TRUE)

onc.sigs <- msigdbr::msigdbr(species="Homo sapiens",category="C6")
onc.sigs <- dplyr::distinct(onc.sigs[onc.sigs$gs_name %in% all.gsea$Feature_set,c("gs_name","gene_symbol")]) # 2824 rows
colnames(onc.sigs) <- c("Feature_set","Gene")
onc.sigs <- plyr::ddply(onc.sigs, .(Feature_set), summarize,
                   Gene = paste0(unique(Gene), collapse = " "))
all.gsea2 <- merge(onc.sigs, all.gsea2, by="Feature_set", all.y = TRUE)
all.gsea2[is.na(all.gsea2$Gene) & !is.na(all.gsea2$Gene.x),]$Gene <- 
  all.gsea2[is.na(all.gsea2$Gene) & !is.na(all.gsea2$Gene.x),]$Gene.x
all.gsea2[is.na(all.gsea2$Gene) & !is.na(all.gsea2$Gene.y),]$Gene <- 
  all.gsea2[is.na(all.gsea2$Gene) & !is.na(all.gsea2$Gene.y),]$Gene.y
all.gsea2$Gene.x <- NULL
all.gsea2$Gene.y <- NULL

all.gsea2$Rank <- rank(-all.gsea2$NES)
bestRank <- max(all.gsea2$Rank)
bestRankGSEA <- max(all.gsea2$Rank)
all.gsea2$BestRank <- FALSE
all.gsea2[all.gsea2$Rank == bestRank,]$BestRank <- TRUE
all.gsea2$percentile <- percent_rank(all.gsea2$NES)

# get correlated drugs and their targets
drug.corr <- read.csv(synapser::synGet("syn64025324")$path)
drug.corr <- drug.corr[drug.corr$type %in% c("RNA", "Protein") & drug.corr$Significant,] # 264
drug.info <- read.csv("https://raw.githubusercontent.com/BelindaBGarana/DMEA/refs/heads/shiny-app/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv")
drug.info$Drug <- gsub("-",".",drug.info$name)
drug.info$Drug <- gsub("[(]",".",drug.info$Drug)
drug.info$Drug <- gsub("[)]",".",drug.info$Drug)
drug.info$Drug <- gsub("[+]",".",drug.info$Drug)
drug.info <- dplyr::distinct(drug.info[drug.info$Drug %in% drug.corr$feature,c("name","Drug","moa","target")]) # 253
colnames(drug.corr)[2] <- "Drug"
drug.corr.wInfo <- merge(drug.info, drug.corr, by="Drug") # 256
drug.corr.wInfo <- drug.corr.wInfo[!is.na(drug.corr.wInfo$target),] # 207
drug.targets <- unique(unlist(strsplit(drug.corr.wInfo$target, ", "))) # 291
drug.corr.wInfo$Rank <- rank(-drug.corr.wInfo$Pearson.est)
bestRank <- max(drug.corr.wInfo$Rank)
bestRankDrug <- max(drug.corr.wInfo$Rank)
drug.corr.wInfo$BestRank <- FALSE
drug.corr.wInfo[drug.corr.wInfo$Rank == bestRank,]$BestRank <- TRUE
drug.corr.wInfo$percentile <- percent_rank(-drug.corr.wInfo$Pearson.est)

# # get enriched drug MOAs and their members 
moa.results <- read.csv(synapser::synGet("syn64025318")$path)
moa.results <- moa.results[moa.results$type %in% c("RNA", "Protein") & moa.results$Significant,] # 19
colnames(moa.results)[3] <- "moa"
moa.info <- dplyr::distinct(drug.info[drug.info$moa %in% moa.results$moa,c("moa","target")]) # 60 MOAs
moa.info <- plyr::ddply(moa.info, .(moa), summarize,
                        Gene = paste0(unique(target), collapse = " "))
moa.results2 <- merge(moa.results, moa.info, by="moa", all.x = TRUE)
moa.results2$Rank <- rank(-moa.results2$NES)
bestRank <- max(moa.results2$Rank)
bestRankMOA <- max(moa.results2$Rank)
moa.results2$BestRank <- FALSE
moa.results2[moa.results2$Rank == bestRank,]$BestRank <- TRUE
moa.results2$percentile <- percent_rank(-moa.results2$NES)

# get selected network nodes and their centrality
net.centr <- read.csv(file.path("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/",
                                paste0("PCSF_Protein_", "2024-12-12"),
                                "centrality_optimalRun.csv"))
net.centr$Rank <- rank(net.centr$eigen_centrality)
bestRank <- max(net.centr$Rank)
bestRankNet <- max(net.centr$Rank)
net.centr$BestRank <- FALSE
net.centr[net.centr$Rank == bestRank,]$BestRank <- TRUE
net.centr$percentile <- percent_rank(net.centr$eigen_centrality)

all.features <- unique(c(corr.result$Gene, kin2$Gene, tf2$Gene,
                         unlist(strsplit(tf2$Targets, " ")),
                         unlist(strsplit(kin2$Substrates, " ")),
                         unlist(strsplit(all.gsea2$Gene, " ")), 
                         unlist(strsplit(moa.results2$Gene, " ")), 
                         drug.targets, net.centr$name)) # 9585
other.features <- unique(c(corr.result$Gene, kin2$Gene, tf$Gene,
                           unlist(strsplit(tf2$Targets, " ")),
                           unlist(strsplit(kin2$Substrates, " ")),
                           unlist(strsplit(all.gsea2$Gene, " ")), 
                           net.centr$name)) # 9431

# all.features <- unique(c(kin$Gene, kin$Gene, tf$Gene,
#                          drug.targets, net.centr$name)) # 2063
# other.features <- unique(c(kin$Gene, kin$Gene, tf$Gene, net.centr$name)) # 1852
Target <- drug.targets[drug.targets %in% other.features] # 199; identify drug targets also implicated in at least 1 other analysis
excluded.targets <- drug.targets[!(drug.targets %in% other.features)] # 92; for example, ADRA1A is first and uniprot protein name is ADA1A
# so I think these are gene symbols

scores <- data.frame(Target)
scores[,c("Correlation", "GSEA", "TF", "TFT", "Kinase", "Substrate", "Network", 
          "Drug", "DMEA","SumRanks","SumBestRanks",
          "nRanks","nBestRanks")] <- NA
scores$BestRanks <- ""
sumScores <- data.frame(Target)
sumScores[,c("Correlation", "GSEA", "TF", "TFT", "Kinase", "Substrate", "Network", 
             "Drug", "DMEA","SumRanks","SumBestRanks",
             "nRanks","nBestRanks")] <- NA
sumScores$BestRanks <- ""
percScores <- data.frame(Target)
percScores[,c("Correlation", "GSEA", "TF", "TFT", "Kinase", "Substrate", "Network", 
             "Drug", "DMEA","MeanPerc","nPerc","nBestPerc")] <- NA
percScores$BestPercs <- ""
#scores[,c("Correlation", "TF", "Kinase", "Network")] <- NA
for (i in 1:length(Target)) {
  runningScore <- 0
  nRankVals <- 0
  sumBestRanks <- 0
  nBestRanks <- 0
  
  # correlations
  if (Target[i] %in% unique(corr.result$Gene)) {
    tempCorr <- corr.result[corr.result$Gene == Target[i],]
    tempTypes <- paste0(unique(tempCorr$type), collapse = " ")
    tempScores <- round(mean(tempCorr$Spearman.est), digits=2)
    tempSumRanks <- sum(tempCorr$Rank)
    tempRanks <- round(mean(tempCorr$Rank), digits=2)
    tempPerc <- mean(tempCorr$percentile)
    
    # keep score
    nRankVals <- nRankVals + nrow(tempCorr)
    runningScore <- runningScore + sum(tempCorr$Rank, na.rm = TRUE)
    if (any(tempCorr$Rank == bestRankCorr)) {
      nBestRanks <- nBestRanks + length(which(tempCorr$Rank == bestRankCorr))
      sumBestRanks <- sumBestRanks + bestRankCorr*length(which(tempCorr$Rank == bestRankCorr))
      scores$BestRanks[i] <- paste(scores$BestRanks[i], "Correlation")
    }
    
    # update data frame
    scores$Correlation[i] <- paste0(tempTypes, " | ", tempScores, " | ", tempRanks)
    sumScores$Correlation[i] <- tempSumRanks
    percScores$Correlation[i] <- tempPerc
  }
  
  # GSEA
  if (Target[i] %in% unique(unlist(strsplit(all.gsea2$Gene, " ")))) {
    tempGSEA <- all.gsea2[grepl(Target[i], all.gsea2$Gene),]
    tempTypes <- paste0(unique(tempGSEA$type), collapse = " ")
    tempSets <- paste0(tempGSEA$Feature_set, collapse = " ")
    tempScores <- round(mean(tempGSEA$NES), digits=2)
    tempRanks <- round(mean(tempGSEA$Rank), digits=2)
    tempSumRanks <- sum(tempGSEA$Rank)
    tempPerc <- mean(tempGSEA$percentile)
    
    # keep score
    nRankVals <- nRankVals + nrow(tempGSEA)
    runningScore <- runningScore + sum(tempGSEA$Rank, na.rm = TRUE)
    if (any(tempGSEA$Rank == bestRankGSEA)) {
      nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankGSEA))
      sumBestRanks <- sumBestRanks + bestRankGSEA*length(which(tempGSEA$Rank == bestRankGSEA))
      scores$BestRanks[i] <- paste(scores$BestRanks[i], "GSEA")
    }
    
    # update data frame
    scores$GSEA[i] <- paste0(tempTypes, " | ", tempSets, " | ", tempScores, " | ", tempRanks)
    sumScores$GSEA[i] <- tempSumRanks
    percScores$GSEA[i] <- tempPerc
  }
  
  # TF
  if (Target[i] %in% tf2$Gene) {
    tempGSEA <- tf2[tf2$Gene == Target[i],]
    tempScores <- round(mean(tempGSEA$NES), digits=2)
    tempRanks <- round(mean(tempGSEA$Rank), digits=2)
    tempSumRanks <- sum(tempGSEA$Rank)
    tempPerc <- mean(tempGSEA$percentile)
    
    # keep score
    nRankVals <- nRankVals + nrow(tempGSEA)
    runningScore <- runningScore + sum(tempGSEA$Rank, na.rm = TRUE)
    if (any(tempGSEA$Rank == bestRankTF)) {
      nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankTF))
      sumBestRanks <- sumBestRanks + bestRankTF*length(which(tempGSEA$Rank == bestRankTF))
      scores$BestRanks[i] <- paste(scores$BestRanks[i], "TF")
    }
    
    # update data frame
    scores$TF[i] <- paste0(tempScores, " | ", tempRanks)
    sumScores$TF[i] <- tempSumRanks
    percScores$TF[i] <- tempPerc
  }
  
  # TFT
  if (Target[i] %in% unique(unlist(strsplit(tf2$Targets, " ")))) {
    tempGSEA <- tf2[grepl(Target[i], tf2$Targets),]
    tempSets <- paste0(tempGSEA$Gene, collapse = " ")
    tempScores <- round(mean(tempGSEA$NES), digits=2)
    tempRanks <- round(mean(tempGSEA$Rank), digits=2)
    tempSumRanks <- sum(tempGSEA$Rank)
    tempPerc <- mean(tempGSEA$percentile)
    
    # keep score
    nRankVals <- nRankVals + nrow(tempGSEA)
    runningScore <- runningScore + sum(tempGSEA$Rank, na.rm = TRUE)
    if (any(tempGSEA$Rank == bestRankTF)) {
      nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankTF))
      sumBestRanks <- sumBestRanks + bestRankTF*length(which(tempGSEA$Rank == bestRankTF))
      scores$BestRanks[i] <- paste(scores$BestRanks[i], "TFT")
      sumScores$BestRanks[i] <- paste(scores$BestRanks[i], "TFT")
    }
    
    # update data frame
    scores$TFT[i] <- paste0(tempSets, " | ", tempScores, " | ", tempRanks)
    sumScores$TFT[i] <- tempSumRanks
    percScores$TFT[i] <- tempPerc
  }
  
  # Kinase
  if (Target[i] %in% kin2$Gene) {
    tempGSEA <- kin2[kin2$Gene == Target[i],]
    tempScores <- round(mean(tempGSEA$enrichment_value_log2),2)
    tempRanks <- round(mean(tempGSEA$Rank),2)
    tempSumRanks <- sum(tempGSEA$Rank)
    tempPerc <- mean(tempGSEA$percentile)
    
    # keep score
    nRankVals <- nRankVals + nrow(tempGSEA)
    runningScore <- runningScore + sum(tempGSEA$Rank, na.rm = TRUE)
    if (any(tempGSEA$Rank == bestRankKin)) {
      nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankKin))
      sumBestRanks <- sumBestRanks + bestRankKin*length(which(tempGSEA$Rank == bestRankKin))
      scores$BestRanks[i] <- paste(scores$BestRanks[i], "Kinase")
    }
    
    # update data frame
    scores$Kinase[i] <- paste0(tempScores, " | ", tempRanks)
    sumScores$Kinase[i] <- tempSumRanks
    percScores$Kinase[i] <- tempPerc
  }
  
  # Substrate
  if (Target[i] %in% unique(unlist(strsplit(kin2$Substrates, " ")))) {
    tempGSEA <- kin2[grepl(Target[i], kin2$Substrates),]
    tempSets <- paste0(tempGSEA$Gene, collapse = " ")
    tempScores <- round(mean(tempGSEA$enrichment_value_log2),2)
    tempRanks <- round(mean(tempGSEA$Rank),2)
    tempSumRanks <- sum(tempGSEA$Rank)
    tempPerc <- mean(tempGSEA$percentile)
    
    # keep score
    nRankVals <- nRankVals + nrow(tempGSEA)
    runningScore <- runningScore + sum(tempGSEA$Rank, na.rm = TRUE)
    if (any(tempGSEA$Rank == bestRankKin)) {
      nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankKin))
      sumBestRanks <- sumBestRanks + bestRankKin*length(which(tempGSEA$Rank == bestRankKin))
      scores$BestRanks[i] <- paste(scores$BestRanks[i], "Substrate")
    }
    
    # update data frame
    scores$Substrate[i] <- paste0(tempSets, " | ", tempScores, " | ", tempRanks)
    sumScores$Substrate[i] <- tempSumRanks
    percScores$Substrate[i] <- tempPerc
  }
  
  # Network
  if (Target[i] %in% net.centr$name) {
    tempGSEA <- net.centr[net.centr$name == Target[i],]
    tempScores <- round(mean(tempGSEA$eigen_centrality),2)
    tempRanks <- round(mean(tempGSEA$Rank),2)
    tempSumRanks <- sum(tempGSEA$Rank)
    tempPerc <- mean(tempGSEA$percentile)
    
    # keep score
    nRankVals <- nRankVals + nrow(tempGSEA)
    runningScore <- runningScore + sum(tempGSEA$Rank, na.rm = TRUE)
    if (any(tempGSEA$Rank == bestRankNet)) {
      nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankNet))
      sumBestRanks <- sumBestRanks + bestRankNet*length(which(tempGSEA$Rank == bestRankNet))
      scores$BestRanks[i] <- paste(scores$BestRanks[i], "Network")
    }
    
    # update data frame
    scores$Network[i] <- paste0(tempScores, " | ", tempRanks)
    sumScores$Network[i] <- tempSumRanks
    percScores$Network[i] <- tempPerc
  }
  
  # Drug correlation
  if (Target[i] %in% drug.corr.wInfo$target) {
    tempGSEA <- drug.corr.wInfo[grepl(Target[i], drug.corr.wInfo$target),]
    tempDrugs <- paste0(tempGSEA$name, collapse = " ")
    tempScores <- round(mean(tempGSEA$Pearson.est),2)
    tempRanks <- round(mean(tempGSEA$Rank),2)
    tempSumRanks <- sum(tempGSEA$Rank)
    tempPerc <- mean(tempGSEA$percentile)
    
    # keep score
    nRankVals <- nRankVals + nrow(tempGSEA)
    runningScore <- runningScore + sum(tempGSEA$Rank, na.rm = TRUE)
    if (any(tempGSEA$Rank == bestRankDrug)) {
      nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankDrug))
      sumBestRanks <- sumBestRanks + bestRankDrug*length(which(tempGSEA$Rank == bestRankDrug))
      scores$BestRanks[i] <- paste(scores$BestRanks[i], "Drug")
    }
    
    # update data frame
    scores$Drug[i] <- paste0(tempDrugs, " | ", tempScores, " | ", tempRanks)
    sumScores$Drug[i] <- tempSumRanks
    percScores$Drug[i] <- tempPerc
  }
  
  # DMEA
  if (Target[i] %in% unique(unlist(strsplit(moa.results2$Gene, " ")))) {
    tempGSEA <- moa.results2[grepl(Target[i], moa.results2$Gene),]
    tempTypes <- paste0(unique(tempGSEA$type), collapse = " ")
    tempSets <- paste0(tempGSEA$moa, collapse = " ")
    tempScores <- round(mean(tempGSEA$NES),2)
    tempRanks <- round(mean(tempGSEA$Rank),2)
    tempSumRanks <- sum(tempGSEA$Rank)
    tempPerc <- mean(tempGSEA$percentile)
    
    # keep score
    nRankVals <- nRankVals + nrow(tempGSEA)
    runningScore <- runningScore + sum(tempGSEA$Rank, na.rm = TRUE)
    if (any(tempGSEA$Rank == bestRankMOA)) {
      nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankMOA))
      sumBestRanks <- sumBestRanks + bestRankMOA*length(which(tempGSEA$Rank == bestRankMOA))
      scores$BestRanks[i] <- paste(scores$BestRanks[i], "DMEA")
    }
    
    # update data frame
    scores$DMEA[i] <- paste0(tempTypes, " | ", tempSets, " | ", tempScores, " | ", tempRanks)
    sumScores$DMEA[i] <- tempSumRanks
    percScores$DMEA[i] <- tempPerc
  }
  
  # update data frame
  scores$SumRanks[i] <- runningScore
  scores$SumBestRanks[i] <- sumBestRanks
  scores$nRanks[i] <- nRankVals
  scores$nBestRanks[i] <- nBestRanks
  percScores$MeanPerc[i] <- mean(as.numeric(percScores[i,]), na.rm = TRUE)
}
percScores$nBestPerc <- scores$nBestRanks
percScores$nPerc <- scores$nRanks
percScores$BestPercs <- scores$BestRanks
sumScores$BestRanks <- scores$BestRanks
sumScores$SumRanks <- scores$SumRanks
sumScores$SumBestRanks <- scores$SumBestRanks
sumScores$nBestRanks <- scores$nBestRanks
sumScores$nRanks <- scores$nRanks
sumScores$meanRank <- sumScores$SumRanks/sumScores$nRanks
percScores$MeanPercTimesN <- percScores$MeanPerc*percScores$nPerc
percScores$MeanPercTimesN <- percScores$MeanPerc*percScores$nPerc
write.csv(scores, "Rank_scores.csv", row.names = FALSE)
write.csv(sumScores, "Rank_score_sums.csv", row.names = FALSE)
write.csv(percScores, "Percentile_scores.csv", row.names = FALSE)
filtered.percScores <- percScores[percScores$nBestPerc > 0,] # 14
filtered.percScores <- filtered.percScores[order(-filtered.percScores$MeanPercTimesN),]
targetOrder <- filtered.percScores$Target
write.csv(filtered.percScores, "Percentile_scores_min1best.csv", row.names = FALSE)
filtered.percScores$nBestPerc <- NULL
filtered.percScores$BestPercs <- NULL
filtered.percScores <- reshape2::melt(filtered.percScores, id = "Target", variable.name = "Analysis")
filtered.percScores$`Top Score` <- FALSE
bestRanks <- percScores$BestPercs
names(bestRanks) <- percScores$Target
bestRanks <- bestRanks[bestRanks != ""]
bestTargets <- names(bestRanks)
for (i in bestTargets) {
  filtered.percScores[filtered.percScores$Target == i &
                        filtered.percScores$Analysis == sub(" ", "",bestRanks[i]),]$`Top Score` <- TRUE
}
# create dot plot to represent filtered.percScores
dot.plot <- ggplot2::ggplot(filtered.percScores[!(filtered.percScores$Analysis %in% c("MeanPerc","nPerc","MeanPercTimesN")),], 
                            aes(x = Target, y = Analysis, color = 100*as.numeric(value), size=4)) + 
                              ggplot2::geom_point() + scale_color_gradient2(low="white",high="red", limits=c(0, 100)) +
                              theme_classic() + ggplot2::labs(color = "Mean Percentile") + ggplot2::scale_x_discrete(limits = targetOrder) +
                              geom_point(data = subset(filtered.percScores, `Top Score`), col = "black", stroke = 1.5, shape = 21) +
                              theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
                                    #axis.text.x = element_text(angle=90,vjust=0.5),
                                    #axis.title.x=element_text(angle=180)
                                    )
dot.plot
saveRDS(dot.plot, "Target_dotPlot.rds")
ggsave("Target_dotPlot.pdf",dot.plot,width=6, height=3)
saveRDS(dot.plot, "Target_dotPlot_horizontal.rds")
ggsave("Target_dotPlot_horizontal.pdf",dot.plot,width=6, height=3)
