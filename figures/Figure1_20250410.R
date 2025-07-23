# chr8 MPNST: figure 1
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(biomaRt);library(RIdeogram);library(viridis)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_1")

synapser::synLogin()

#### 1. Chr8 median copy number bar plot ####
# load data & format
#med.chr8q <- read.csv(synapser::synGet("syn66047330")$path)
med.chr8q <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/proteomics/analysis/Chr8_quant_20250409/positional_medians/Copy Number/Copy Number_Chr8q_median.csv")
colnames(med.chr8q)[2] <- "Median Chr8q Copy Number"

# order lowest to highest
med.chr8q <- med.chr8q[order(med.chr8q$`Median Chr8q Copy Number`),]
sample.order <- med.chr8q$Sample

# create bar plot
med.chr8q$Omics <- "RNA-Seq"
prot.samples <- c("JH-2-002", "JH-2-055", "JH-2-079c", "MN-2", "WU-225", "WU-487")
med.chr8q[med.chr8q$Sample %in% prot.samples,]$Omics <- "RNA-Seq &\n Proteomics"
ggplot2::ggplot(med.chr8q, aes(x=Sample, y=`Median Chr8q Copy Number`, fill = Omics)) + 
  ggplot2::geom_bar(stat="identity", position="dodge", #alpha = 0.5
                    ) + 
  scale_x_discrete(limits = sample.order) + ggplot2::theme_classic(base_size=12) +
  scale_fill_manual(values=c("#C2A5CF","#FEE391"), breaks=c("RNA-Seq","RNA-Seq &\n Proteomics")) +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), legend.position=c(0.3,0.8)) + #, size = 16), 
        #axis.text.y = element_text(size = 16), axis.title=element_text(size=24),
        #legend.title=element_text(size=16),legend.text=element_text(size=12)) +
  geom_errorbar(aes(ymin=`Median Chr8q Copy Number` - sd_med_copy_number, 
                    ymax = `Median Chr8q Copy Number` + sd_med_copy_number), width=0.2,
                position=position_dodge(0.9)) + xlab("MPNST PDX")
ggplot2::ggsave(paste0("PDX Copy Number_Chr8q_median_", Sys.Date(), ".pdf"), width = 3.5, height = 3) # was width 5

#### 2. Chr8 enrichment along chromosome ####
# load enrichment results
chr8.enr <- read.csv(synapser::synGet("syn66227265")$path)

# load segment info for gene symbols
chr8.msigdb <- msigdbr::msigdbr(category = "C1")
chr8.msigdb <- chr8.msigdb[grepl("chr8", chr8.msigdb$gs_name),]

# get chr8 segment positions
segments <- unique(chr8.enr[chr8.enr$FDR_q_value <= 0.25 & chr8.enr$p_value <= 0.05,]$Feature_set)
segments <- segments[grepl("chr8",segments)]

ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")#, host="https://useast.ensembl.org")
segment.pos <- data.frame(segments)
segment.pos[,c("start", "end")] <- NA
for (i in segments) {
  # get gene ids
  temp.genes <- unique(chr8.msigdb[chr8.msigdb$gs_name == i,]$ensembl_gene)
  
  # find start/end of genes
  if (length(temp.genes) > 0) {
    gene.pos <- biomaRt::getBM(attributes=c("ensembl_gene_id","start_position","end_position"),
                               filters="ensembl_gene_id",values = temp.genes,mart=ensembl) 
    segment.pos[segment.pos$segments == i,]$start <- min(gene.pos$start_position)
    segment.pos[segment.pos$segments == i,]$end <- max(gene.pos$end_position)
  }
}

# add enrichment scores
colnames(chr8.enr)[2] <- "segments"
segment.pos <- merge(segment.pos, chr8.enr[,c("segments","NES")], by="segments")
segment.pos$chr <- 8
write.csv(segment.pos, "Chr8_segment_enrichment.csv", row.names = FALSE)
#segment.pos <- read.csv("Chr8_segment_enrichment.csv")

# following vignette: https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html
data("human_karyotype", package="RIdeogram")
human8 <- as.data.frame(human_karyotype[8,])
ideo.df <- segment.pos[,c("chr", "start", "end", "NES")]
colnames(ideo.df) <- c("Chr", "Start", "End", "Value")
label.df <- segment.pos[,c("segments", "chr", "start", "end")]
label.df$Type <- label.df$segments
label.df$Shape <- "box"
label.df$color <- "black"
label.df[label.df$segments=="chr8q24",]$color <- "red"
label.df <- label.df[,c("Type", "Shape", "chr", "start", "end", "color")]
colnames(label.df)[3:5] <- c("Chr", "Start", "End")
RIdeogram::ideogram(karyotype = human_karyotype, overlaid = ideo.df, label = label.df, label_type = "marker", colorset1 = c("blue","grey", "red"))
convertSVG("chromosome.svg", device="pdf")

#### 3. bar plot of # of diffexp Features ####
path.map <- list("RNA" = "syn66226866",
                 "Protein" = "syn66224803",
                 "Phospho" = "syn66226338")

all.degs <- data.frame()
all.deg.list <- list()
sig.deg.list <- list()
up.deg.list <- list()
dn.deg.list <- list()
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
  all.degs <- rbind(all.degs, temp.degs[,2:ncol(temp.degs)])
  nonFeatureNames <- colnames(temp.degs)[2:ncol(temp.degs)][colnames(temp.degs)[2:ncol(temp.degs)] != "Feature"]
  temp.degs <- temp.degs[,c("Feature", nonFeatureNames)]
  all.deg.list[[names(path.map)[i]]] <- temp.degs
  sig.deg.list[[names(path.map)[i]]] <- temp.degs[temp.degs$Significant,]
  up.deg.list[[names(path.map)[i]]] <- temp.degs[temp.degs$Significance == "Positive",]
  dn.deg.list[[names(path.map)[i]]] <- temp.degs[temp.degs$Significance == "Negative",]
}
all.degs$Omics <- factor(all.degs$Omics, levels=c("Phospho", "Protein", "RNA"))
all.degs$Direction <- factor(all.degs$Direction, levels=c("Positive", "Negative"))
write.csv(all.degs, "Differential_expression_results.csv", row.names = FALSE)
synapser::synStore(synapser::File("Differential_expression_results.csv", parent="syn65988130"))
all.degs <- read.csv("Differential_expression_results.csv")
all.degs$Omics <- factor(all.degs$Omics, levels=c("RNA", "Protein", "Phospho"))
chr8q.genes <- readRDS("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/proteomics/chr8qGenes.rds")

plot_df <- plyr::ddply(all.degs[all.degs$Significant,] , .(Direction, Omics), dplyr::summarize,
                       nCorr = dplyr::n(),
                       nChr8q = length(unique(Feature[tolower(sub("-.*", "",Feature)) %in% tolower(chr8q.genes)])))
plot_df$log10nCorr <- log(plot_df$nCorr, 10)

plot_dfMin6 <- plyr::ddply(all.degs[all.degs$Significant & all.degs$N>=6,] , .(Direction, Omics), dplyr::summarize,
                       nCorr = dplyr::n(),
                       nChr8q = length(unique(Feature[tolower(sub("-.*", "",Feature)) %in% tolower(chr8q.genes)])))
plot_dfMin6$log10nCorr <- log(plot_dfMin6$nCorr, 10)

plt <- ggplot(plot_df) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(1, 2, 3, 4)),
    color = "lightgrey"
  ) +
  # Add bars to represent the cumulative track lengths
  # str_wrap(region, 5) wraps the text so each line has at most 5 characters
  # (but it doesn't break long words!)
  geom_col(
    aes(
      x = reorder(Omics, log10nCorr),
      #x = forcats::fct_infreq(stringr::str_wrap(Collection, 5)),
      y = log10nCorr,
      fill = Omics, alpha = Direction # first version
      #fill = Direction
    ),
    position = "dodge", show.legend = TRUE#, alpha=0.5
  ) +
  geom_segment(
    aes(
      x = reorder(Omics, log10nCorr),
      y = 0,
      xend = reorder(Omics, log10nCorr),
      yend = max(log10nCorr)
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
    y =0.5, 
    label = "10", 
    geom = "text", 
    color = "gray12"
  ) +
  annotate(
    x = 0, 
    y = 1.5, 
    label = "100", 
    geom = "text", 
    color = "gray12"
  ) +
  annotate(
    x = 0, 
    y =2.5, 
    label = "1000", 
    geom = "text", 
    color = "gray12"
  )+
  annotate(
    x = 0, 
    y =3.5, 
    label = "10000", 
    geom = "text", 
    color = "gray12"
  ) +
  # Scale y axis so bars don't start in the center
  scale_y_continuous(
    limits = c(-2, 4),
    expand = c(0, 0),
    #breaks = c(0, 10, 100, 1000)#, trans = 'log10'
  ) + #scale_fill_manual(values=c("red", "blue")) +
  scale_fill_manual(values=c("#C2A5CF", "#FEE391", "#B2DF8A")) + # first version; colors were: "#C2A5CF", "#FEE391", "#B2DF8A"
  scale_alpha_discrete(range=c(0.5,1)) + 
  theme(
    # Remove axis ticks and text
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    # Use gray text for the region names
    axis.text.x = element_text(color = "gray12", size = 16), # was size 12
    # Move the legend to the bottom
    legend.position = "bottom",
  ) + ggtitle("Number of Correlated Features") + 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=24), # was size 16
        axis.line=element_blank())

plt
ggplot2::ggsave("Chr8_numberOfFeatures_circularBarPlot.pdf", plt, width = 7, height = 7)


plt <- ggplot(plot_dfMin6) +
  # Make custom panel grid
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(1, 2, 3)),
    color = "lightgrey"
  ) +
  # Add bars to represent the cumulative track lengths
  # str_wrap(region, 5) wraps the text so each line has at most 5 characters
  # (but it doesn't break long words!)
  geom_col(
    aes(
      x = reorder(Omics, log10nCorr),
      #x = forcats::fct_infreq(stringr::str_wrap(Collection, 5)),
      y = log10nCorr,
      fill = Omics, alpha = Direction # first version
      #fill = Direction
    ),
    position = "dodge", show.legend = TRUE#, alpha=0.5
  ) +
  geom_segment(
    aes(
      x = reorder(Omics, log10nCorr),
      y = 0,
      xend = reorder(Omics, log10nCorr),
      yend = max(log10nCorr)
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
    y =0.5, 
    label = "10", 
    geom = "text", 
    color = "gray12"
  ) +
  annotate(
    x = 0, 
    y = 1.5, 
    label = "100", 
    geom = "text", 
    color = "gray12"
  ) +
  annotate(
    x = 0, 
    y =2.5, 
    label = "1000", 
    geom = "text", 
    color = "gray12"
  )+
  # Scale y axis so bars don't start in the center
  scale_y_continuous(
    limits = c(-2, 3),
    expand = c(0, 0),
    #breaks = c(0, 10, 100, 1000)#, trans = 'log10'
  ) + #scale_fill_manual(values=c("red", "blue")) +
  scale_fill_manual(values=c("#C2A5CF", "#FEE391", "#B2DF8A")) + # first version; colors were: "#C2A5CF", "#FEE391", "#B2DF8A"
  scale_alpha_discrete(range=c(0.5,1)) + 
  theme(
    # Remove axis ticks and text
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    # Use gray text for the region names
    axis.text.x = element_text(color = "gray12", size = 16), # was size 12
    # Move the legend to the bottom
    legend.position = "bottom",
  ) + ggtitle("Number of Correlated Features") + 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=24), # was size 16
        axis.line=element_blank())

plt
ggplot2::ggsave("Chr8_numberOfFeatures_circularBarPlot_minN6.pdf", plt, width = 7, height = 7)


#### 4. venn diagram of diffexp Features ####
venn.list <- list()
for (i in 1:length(unique(all.degs$Omics))) {
  temp.name <- as.character(unique(all.degs$Omics))[i]
  if (unique(all.degs[all.degs$Omics==temp.name,]$Feature_type) == "SUB_SITE") {
    venn.list[[temp.name]] <- unique(sub("-.*", "", all.degs[all.degs$Omics==temp.name & all.degs$Significant,]$Feature))
  } else if (unique(all.degs[all.degs$Omics==temp.name,]$Feature_type) == "Gene") {
    venn.list[[temp.name]] <- unique(all.degs[all.degs$Omics==temp.name & all.degs$Significant,]$Feature)
  } else {
    warning("No Gene or SUB_SITE column found")
  }
}
ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_color = c("#C2A5CF", "#FEE391", "#B2DF8A"), 
               fill_alpha = 0.5, text_size = 10, set_name_size = 8)
ggplot2::ggsave(paste0("Chr8_diffExp_vennDiagram_noPercentages_matchingColors_", Sys.Date(), ".pdf"), width=7, height=7)

# what is shared between RNA & protein?
venn.list[["RNA"]][venn.list[["RNA"]] %in% venn.list[["Protein"]]] # "ERG28" - down in RNA with N=4, up in protein with N=12

# what is shared between RNA & phospho?
venn.list[["RNA"]][venn.list[["RNA"]] %in% venn.list[["Phospho"]]] # "GSPT2"  "RBM14"  "RNF169"
# GSPT2 up in RNA with N=3, GSPT2-S191s down in phospho with N=12, not significant in protein
# RBM14 up in RNA with N=14, RBM14-S663s up in phospho with N=12
# RNF169 up in RNA with N=3, RNF169-S292s down in phospho with N=12

# what is shared between protein & phospho?
venn.list[["Protein"]][venn.list[["Protein"]] %in% venn.list[["Phospho"]]] 
# "ARHGEF2" "CALD1"   "FLNB"    "H2AFY"   "MICAL3"  "MON1B"   "PDE4DIP" "PRG4"    "SPTAN1" "SPTBN1"  "WWTR1"  

## what if we require 6+ samples
venn.list <- list()
for (i in 1:length(unique(all.degs$Omics))) {
  temp.name <- as.character(unique(all.degs$Omics))[i]
  if (unique(all.degs[all.degs$Omics==temp.name,]$Feature_type) == "SUB_SITE") {
    venn.list[[temp.name]] <- unique(sub("-.*", "", all.degs[all.degs$Omics==temp.name & all.degs$Significant & all.degs$N>=6,]$Feature))
  } else if (unique(all.degs[all.degs$Omics==temp.name,]$Feature_type) == "Gene") {
    venn.list[[temp.name]] <- unique(all.degs[all.degs$Omics==temp.name & all.degs$Significant & all.degs$N>=6,]$Feature)
  } else {
    warning("No Gene or SUB_SITE column found")
  }
}
ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_color = c("#C2A5CF", "#FEE391", "#B2DF8A"), 
               fill_alpha = 0.5, text_size = 10, set_name_size = 8)
ggplot2::ggsave(paste0("Chr8_diffExp_vennDiagram_noPercentages_matchingColors_minN6_", Sys.Date(), ".pdf"), width=7, height=7)

# what is shared between RNA & protein?
venn.list[["RNA"]][venn.list[["RNA"]] %in% venn.list[["Protein"]]] # no overlap

# what is shared between RNA & phospho?
venn.list[["RNA"]][venn.list[["RNA"]] %in% venn.list[["Phospho"]]] # "RBM14"

# what is shared between protein & phospho?
venn.list[["Protein"]][venn.list[["Protein"]] %in% venn.list[["Phospho"]]] 
# "ARHGEF2" "CALD1"   "FLNB"    "H2AFY"   "MICAL3"  "MON1B"   "PDE4DIP" "PRG4"    "SPTAN1" "SPTBN1"  "WWTR1"  


#### 5. top up/dn diffexp Features ####
# try making dot plot of Features sig in 2+ omics types
mean.DEG.df <- plyr::ddply(all.degs, .(Feature), summarize,
                           mean_est = mean(Spearman.est),
                           sd_est = sd(Spearman.est),
                           Fisher_p = metap::sumlog(na.omit(Spearman.p))$p,
                           types = paste0(Omics, collapse = ", "),
                           N_types = length(unique(Omics)),
                           N_sig = length(Significant[Significant]),
                           sig_types = paste0(Omics[Significant], collapse = ", "))
if (length(unique(mean.DEG.df$Fisher_p)) > 1) {
  mean.DEG.df$adj_Fisher_p <- 
    qvalue::qvalue(mean.DEG.df$Fisher_p, pi0=1)$qvalues
} else {
  mean.DEG.df$adj_Fisher_p <- NA
}
write.csv(mean.DEG.df, "Compiled_differential_expression_results.csv", row.names=FALSE)
synapser::synStore(synapser::File("Compiled_differential_expression_results.csv", parent="syn65988130"))
mean.DEG.df <- read.csv(synapser::synGet("syn66227736")$path)

#### do bar plots instead and make sure each omics has the same number of features ####
omics <- c("RNA", "Protein", "Phospho")
maxAbsEst <- max(na.omit(abs(all.degs[all.degs$Omics %in% omics,]$Spearman.est)))
maxAbsEst <- 1
all.degs <- na.omit(all.degs)
if (any(all.degs$Spearman.q == 0)) {
  all.degs[all.degs$Spearman.q == 0,]$Spearman.q <- 0.0001
}
all.degs$minusLogP <- -log(all.degs[,"Spearman.p"], base = 10)
all.degs$minusLogFDR <- -log(all.degs[,"Spearman.q"], base = 10)
maxN <- max(all.degs[all.degs$Significant & 
                        all.degs$Omics %in% omics,]$N) # 15
library(patchwork); library(ggplot2)
chr8q.genes <- readRDS("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/proteomics/chr8qGenes.rds")
gsea.dot.plots <- NULL
for (i in omics) {
  nGenes <- length(unique(na.omit(all.degs[all.degs$Omics == i & 
                                             all.degs$N>=6,])$Feature))
  nCorrGenes <- length(unique(na.omit(all.degs[all.degs$Significant & 
                                                 all.degs$Omics == i & 
                                                 all.degs$N>=6,])$Feature))
  
  nChr8qGenes <- length(unique(na.omit(all.degs[all.degs$Omics == i & 
                                                  all.degs$N>=6 & 
                                                  tolower(sub("-.*", "",all.degs$Feature)) %in% 
                                                  tolower(chr8q.genes),])$Feature))
  nChr8qCorrGenes <- length(unique(na.omit(all.degs[all.degs$Significant & 
                                                 all.degs$Omics == i & 
                                                   all.degs$N>=6 & 
                                                   tolower(sub("-.*", "",all.degs$Feature)) %in% 
                                                   tolower(chr8q.genes),])$Feature))
  cat("corr chr8q genes:", paste0(unique(na.omit(all.degs[all.degs$Significant & 
                                                            all.degs$N>=6 & 
                                                     all.degs$Omics == i & 
                                                     tolower(sub("-.*", "",all.degs$Feature)) %in% 
                                                     tolower(chr8q.genes),])$Feature)), collapse=", ","\n")
  
  nGenes <- length(unique(na.omit(all.degs[all.degs$Omics == i,])$Feature))
  nCorrGenes <- length(unique(na.omit(all.degs[all.degs$Significant & 
                                                 all.degs$Omics == i,])$Feature))
  
  # topGenes <- unique(na.omit(all.degs[all.degs$Significant & all.degs$N >= 9 &
  #                                       abs(all.degs$Spearman.est) >= 0.9 &
  #                                            all.degs$Omics == i,])$Feature)
  temp.dot.df <- all.degs[all.degs$Significant & all.degs$N >= 6 &
                            all.degs$Omics == i,] %>% slice_max(abs(Spearman.est), n = 50)
  topGenes <- temp.dot.df$Feature
  nTopGenes <- length(topGenes)
  geneOrder <- temp.dot.df[order(-temp.dot.df$Spearman.est),]$Feature
  # plot.annot <- paste0(i, "\n(", nCorrGenes, " / ", nGenes, " correlated in 6+ samples, ",
  #                      nChr8qCorrGenes," of which belong to chr8q out of ", nChr8qGenes," quantified)")
  plot.annot <- paste0(i, "\n(", nCorrGenes, " / ", nGenes, " correlated)")
  
  
  # geneFace <- rep("plain", length(geneOrder))
  # names(geneFace) <- geneOrder
  # MYC.pathways <- geneOrder[grepl("MYC", geneOrder)]
  # SHH.pathways <- geneOrder[grepl("SHH", geneOrder) | grepl("HEDGEHOG", geneOrder)]
  # poi <- c(MYC.pathways, SHH.pathways)
  # geneFace[poi] <- "bold"
  
  geneFace <- rep("plain", length(geneOrder))
  names(geneFace) <- geneOrder
  if (i == "Phospho") {
    poi <- geneOrder[tolower(sub("-.*","",geneOrder)) %in% tolower(chr8q.genes)]
    if (length(poi) > 0) {
     cat("bolding chr8q genes: ", paste0(geneOrder[tolower(sub("-.*","",geneOrder)) %in% tolower(chr8q.genes)],collapse=", "))
    }
  } else {
    poi <- geneOrder[tolower(geneOrder) %in% tolower(chr8q.genes)]
    if (length(poi) > 0) {
      cat("bolding chr8q genes: ", paste0(geneOrder[tolower(geneOrder) %in% tolower(chr8q.genes)],collapse=", "))
    }
  }

  geneFace[poi] <- "bold"
  
  dot.plot <- ggplot2::ggplot(temp.dot.df,
                              ggplot2::aes(
                                y = Feature, x = N, fill = Spearman.est
                              )
  ) + #scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
    ggplot2::geom_col() + 
    ggplot2::scale_y_discrete(limits = geneOrder) +
    scale_fill_gradient2(low="blue",high="red", mid="gray",
      limits=c(-1,1)) +
    scale_x_continuous(breaks=c(0,3,6,9,12))+
    #viridis::scale_color_viridis() +
    theme_classic() +
    ggplot2::labs(
      x = "N",
      fill = "Spearman rho"
    ) + 
    theme(axis.text.y=element_text(face=geneFace)) +
    theme(axis.title.x=element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
  
  cat(plot.annot, "\n")
  cat(i, "with", nTopGenes, "top features\n")
  dot.plot <- dot.plot + coord_flip() + ggtitle(plot.annot) + 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
  dot.plot
  if (is.null(gsea.dot.plots)) {
    gsea.dot.plots <- dot.plot
  } else {
    gsea.dot.plots <- gsea.dot.plots / dot.plot
  } 
}
#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4/guides_build_mod.R")
gsea.dot.plots2 <- gsea.dot.plots + plot_layout(guides = 'collect')
gsea.dot.plots2
#ggplot2::ggsave("CorrelatedFeatures_dotPlot_patchworkOmics_minN6_sliceMaxAbsSpearman50_colorN.pdf", gsea.dot.plots, width=14, height=7)
#ggplot2::ggsave("CorrelatedFeatures_dotPlot_patchworkOmics_minN6_sliceMaxAbsSpearman50_yN_chr8q_minN6.pdf", gsea.dot.plots, width=14, height=7)
#ggplot2::ggsave("CorrelatedFeatures_dotPlot_patchworkOmics_minN6_sliceMaxAbsSpearman50_yNv2_chr8q_minN6.pdf", gsea.dot.plots2, width=14, height=7)
ggplot2::ggsave("CorrelatedFeatures_dotPlot_patchworkOmics_minN6_sliceMaxAbsSpearman50_yN_chr8qBold.pdf", gsea.dot.plots, width=14, height=7)
ggplot2::ggsave("CorrelatedFeatures_dotPlot_patchworkOmics_minN6_sliceMaxAbsSpearman50_yNv2_chr8qBold.pdf", gsea.dot.plots2, width=14, height=7)

# "EXT1" and "UBE2V2" are chr8q genes in protein - not sure why not bolded

# what % of up- or down-regulated genes are on chr8q? chr8p?
# chr8q corr in RNA: DSTNP3 GPR20 LRRC24 SCX TMEM276-ZFTRAF1 VCPIP1 ZHX1-C8orf76 - none have 6+ samples; 1 was down, 6 were up; DSTNP3 was down but N=3
# ZHX1 not quantified alone, TMEM276 not significant alone 
# chr8q corr in prot: DSCC1 ENY2 EXT1 LYPLA1 MRPL13 NUDCD1 PLEKHF2 RB1CC1 RNF19A UBE2V2 - all have 6+ samples, all up with chr8q
# chr8q corr in phospho: MTFR1-S119s NBN-S604s RIMS2-Y317y SHARPIN-S165s - all have 6+ samples; 3 down, 1 up: MTFR1-S119s is up with chr8q

# any overlapping genes?
temp.dot.df <- all.degs[all.degs$Significant & all.degs$N >= 6 &
                          all.degs$Omics == "RNA",] %>% slice_max(abs(Spearman.est), n = 50)
topGenesRNA <- temp.dot.df$Feature # 50
temp.dot.df <- all.degs[all.degs$Significant & all.degs$N >= 6 &
                          all.degs$Omics == "Protein",] %>% slice_max(abs(Spearman.est), n = 50)
topGenesProt <- temp.dot.df$Feature # 55
topGenesRNAProt <- topGenesProt[topGenesProt %in% topGenesRNA] # 0

temp.dot.df <- all.degs[all.degs$Significant & all.degs$N >= 6 &
                          all.degs$Omics == "Phospho",] %>% slice_max(abs(Spearman.est), n = 50)
topGenesPhos <- unique(sub("-.*","",temp.dot.df$Feature)) # 61
topGenesRNAPhos <- topGenesPhos[topGenesPhos %in% topGenesRNA] # 0
topGenesProtPhos <- topGenesPhos[topGenesPhos %in% topGenesProt] # 0
