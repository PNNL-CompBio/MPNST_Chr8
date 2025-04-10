# chr8 MPNST: figure 1
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(biomaRt);library(RIdeogram);library(viridis)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_1")

synapser::synLogin()

#### 1. Chr8 median copy number bar plot ####
# load data & format
med.chr8q <- read.csv(synapser::synGet("syn66047330")$path)
colnames(med.chr8q)[2] <- "Median Chr8q Copy Number"

# order lowest to highest
med.chr8q <- med.chr8q[order(med.chr8q$`Median Chr8q Copy Number`),]
sample.order <- med.chr8q$Sample

# create bar plot
med.chr8q$Omics <- "RNA-Seq"
prot.samples <- c("JH-2-002", "JH-2-055", "JH-2-079c", "MN-2", "WU-225", "WU-487")
med.chr8q[med.chr8q$Sample %in% prot.samples,]$Omics <- "RNA-Seq &\n Proteomics"
ggplot2::ggplot(med.chr8q, aes(x=Sample, y=`Median Chr8q Copy Number`, fill = Omics)) + 
  ggplot2::geom_bar(stat="identity", position="dodge", alpha = 0.5) + 
  scale_x_discrete(limits = sample.order) + ggplot2::theme_classic(base_size=12) +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) + #, size = 16), 
        #axis.text.y = element_text(size = 16), axis.title=element_text(size=24),
        #legend.title=element_text(size=16),legend.text=element_text(size=12)) +
  geom_errorbar(aes(ymin=`Median Chr8q Copy Number` - sd_med_copy_number, 
                    ymax = `Median Chr8q Copy Number` + sd_med_copy_number), width=0.2,
                position=position_dodge(0.9)) + xlab("MPNST PDX")
ggplot2::ggsave(paste0("PDX Copy Number_Chr8q_median_", Sys.Date(), ".pdf"), width = 5, height = 3)

#### 2. Chr8 enrichment along chromosome ####
# load enrichment results
#chr8.enr <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/proteomics/analysis/Chr8_quant/Spearman/Copy_number/Copy_number/GSEA/Positional/GSEA_results.csv")
#chr8.enr <- read.csv(synapser::synGet("syn63434853")$path) - update synID

# load segment info for gene symbols
chr8.msigdb <- msigdbr::msigdbr(species = "Homo sapiens", category = "C1")
chr8.msigdb <- chr8.msigdb[grepl("chr8", chr8.msigdb$gs_name),]

# get chr8 segment positions
segments <- unique(chr8.enr[chr8.enr$FDR_q_value <= 0.25 & chr8.enr$p_value <= 0.05,]$Feature_set)
segments <- segments[grepl("chr8",segments)]

ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl", host="https://useast.ensembl.org")
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
segment.pos <- read.csv("Chr8_segment_enrichment.csv")

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
path.map <- list("RNA" = "syn63602832",
                 "Protein" = "syn63435101",
                 "Phospho" = "syn63425909")
compil.path <- paste0("Compiled_results_allSamples_",Sys.Date())

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
bar.plot3log <- ggplot2::ggplot(all.degs[all.degs$Significant,], aes(fill = Direction, x=Omics)) + 
  geom_bar(position = "dodge", stat="count") + ggplot2::theme_classic() + 
  scale_y_continuous(trans='log10') + ggplot2::xlab("Omics Type") +
  scale_fill_manual(values=c("red", "blue")) + 
  ggplot2::ylab("Number of Differentially Expressed Features") + coord_flip() #+
ggplot2::ggsave("Chr8_numberOfFeatures_directionGroup_barPlot_logScale_manualOrderV2.pdf", bar.plot3log, width = 7, height = 7, device = "pdf")

all.degs$Omics <- factor(all.degs$Omics, levels=c("RNA", "Protein", "Phospho"))
plot_df <- plyr::ddply(all.degs[all.degs$Significant,] , .(Direction, Omics), dplyr::summarize,
                       nCorr = dplyr::n())
plot_df$log10nCorr <- log(plot_df$nCorr, 10)

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
  scale_fill_manual(values=c("#F8766D", "#00BFC4", "#7CAE00")) + # first version
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
      #fill = Omics, alpha = Direction # first version
      fill = Direction
    ),
    position = "dodge", show.legend = TRUE, alpha=0.5
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
  ) + scale_fill_manual(values=c("red", "blue")) + scale_alpha_discrete(range=c(0.5,1)) +
  #scale_fill_manual(values=c("#F8766D", "#00BFC4", "#7CAE00")) + # first version
  theme(
    # Remove axis ticks and text
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    # Use gray text for the region names
    axis.text.x = element_text(color = "gray12", size = 12),
    # Move the legend to the bottom
    legend.position = "bottom",
  ) + ggtitle("Number of Correlated Features") + 
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
        axis.line=element_blank())
plt
ggplot2::ggsave("Chr8_numberOfFeatures_circularBarPlot_v2.pdf", plt, width = 7, height = 7)

#### 4. venn diagram of diffexp Features ####
venn.list <- list()
for (i in 1:length(sig.deg.list)) {
  temp.name <- names(sig.deg.list)[i]
  if (unique(sig.deg.list[[i]]$Feature_type) == "SUB_SITE") {
    venn.list[[temp.name]] <- unique(sub("-.*", "", sig.deg.list[[i]]$Feature))
  } else if (unique(sig.deg.list[[i]]$Feature_type) == "Gene") {
    venn.list[[temp.name]] <- unique(sig.deg.list[[i]]$Feature)
  } else {
    warning("No Gene or SUB_SITE column found")
  }
}
ggvenn::ggvenn(venn.list)
ggplot2::ggsave(paste0("Chr8_diffExp_vennDiagram_", Sys.Date(), ".pdf"), width=7, height=7)

ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_alpha = 0.25, text_size = 10, set_name_size = 8)
ggplot2::ggsave(paste0("Chr8_diffExp_vennDiagram_noPercentages_", Sys.Date(), ".pdf"), width=7, height=7)

ggvenn::ggvenn(venn.list, show_percentage = FALSE, fill_color = c("#F8766D", "#00BFC4", "#7CAE00"), 
               fill_alpha = 0.5, text_size = 10, set_name_size = 8)
ggplot2::ggsave(paste0("Chr8_diffExp_vennDiagram_noPercentages_matchingColors_", Sys.Date(), ".pdf"), width=7, height=7)
synapser::synStore(synapser::File(paste0("Chr8_diffExp_vennDiagram_", Sys.Date(), ".pdf"), parent = synID))

#### 5. top up/dn diffexp Features ####
# try making dot plot of Features sig in 2+ omics types
dot.df <- read.csv(synapser::synGet("syn64024301")$path)
dot.df <- dot.df[grepl("RNA, Protein", dot.df$type),]
dot.df <- dot.df[dot.df$N_sig>1,] # 6 between RNA & protein
dot.df <- dotPlot$df # 789

dot.df <- na.omit(all.degs[all.degs$Significant & all.degs$N >= 9 & all.degs$Feature_type=="Gene",]) # 722; 757 if N>=6; 1622 if don't filter for N
nrow(dot.df[dot.df$Omics=="RNA",]) # 19; 52 if N>=6
nrow(dot.df[dot.df$Omics=="Protein",]) # 703; 705 if N>=6
hist(dot.df$Spearman.est,breaks=20)
min(abs(dot.df$Spearman.est)) # 0.7632938 if N>=6 or N>=9
nrow(dot.df[abs(dot.df$Spearman.est)>=0.9,]) # 117; 151 if N>=6
nrow(dot.df[abs(dot.df$Spearman.est)>=0.95,]) # 14; 31 if N>=6
nrow(dot.df[abs(dot.df$Spearman.est)>=0.92,]) # 58; 91 if N>=6
top.dot.df <- dot.df %>% dplyr::slice_max(abs(Spearman.est), n = 100) # 115; 111 if N>=6
min(abs(top.dot.df$Spearman.est)) # 0.9046445; 0.9187795 if N>=6

topGenes <- unique(dot.df[abs(dot.df$Spearman.est)>=0.9,]$Feature)
mean.dot.df <- read.csv(synapser::synGet("syn64024301")$path)
mean.dot.df <- mean.dot.df[mean.dot.df$feature %in% topGenes,]
geneOrder <- mean.dot.df[order(mean.dot.df$mean_rankVal),]$feature
dot.df <- all.degs[all.degs$Feature %in% topGenes,]
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.plot <- ggplot2::ggplot(
  na.omit(dot.df),
  ggplot2::aes(
    x = Omics, y = Feature, color = Spearman.est,
    size = minusLogFDR
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(color = "Spearman Correlation", size = "-log(FDR)")
dot.plot
# most are only in prot and not in RNA
mean.dot.df <- read.csv(synapser::synGet("syn64024301")$path)
topGenes <- mean.dot.df[mean.dot.df$N_sig> 0,]$feature # 178 
topGenes <- na.omit(mean.dot.df[(grepl("RNA", mean.dot.df$sig_types) | grepl("Protein", mean.dot.df$sig_types)) &
                                  grepl("RNA, Protein", mean.dot.df$types) & mean.dot.df$adj_Fisher_p <= 0.05,]$feature) # 50
top.mean.df <- mean.dot.df[mean.dot.df$feature %in% topGenes,]
geneOrder <- top.mean.df[order(top.mean.df$mean_rankVal),]$feature
dot.df <- all.degs[all.degs$Feature %in% topGenes,]
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.plot <- ggplot2::ggplot(
  na.omit(dot.df),
  ggplot2::aes(
    x = Omics, y = Feature, color = Spearman.est,
    size = minusLogFDR
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_x_discrete(limits=c("RNA","Protein")) +
  viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Gene",
    color = "Spearman Correlation", size = "-log(FDR)"
  )
dot.plot
ggplot2::ggsave("Correlated_genes_adjFisherP0.05_dotPlot.pdf", dot.plot, width=4, height=9)

# another style
geneOrder <- mean.dot.df[order(mean.dot.df$mean_rankVal) & mean.dot.df$feature %in% topGenes,]$feature
geneFace <- rep("plain", length(geneOrder))
names(geneFace) <- geneOrder
sigInRNA <- mean.dot.df[grepl("RNA", mean.dot.df$sig_types) & mean.dot.df$feature %in% topGenes,]$feature
geneFace[sigInRNA] <- "italic"
sigInBoth <- mean.dot.df[grepl("RNA, Protein", mean.dot.df$sig_types) & mean.dot.df$feature %in% topGenes,]$feature
geneFace[sigInBoth] <- "bold"
# if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style

dot.df <- na.omit(all.degs[all.degs$Feature%in% topGenes &
                        all.degs$Omics %in% c("RNA", "Protein"),])
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
maxLogFDR <- max(dot.df$minusLogFDR)
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = Omics, y = Feature, color = Spearman.est,
    size = minusLogFDR
  )
) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,5)) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="blue",high="red") +
  #viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Gene",
    color = "Spearman Correlation", size = "-log(adjusted p-value)"
  ) + 
  geom_point(data = subset(dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Correlated_genes_adjFisherP0.05_dotPlot_v2.pdf", dot.plot, width=4, height=9)

mean.dot.df <- read.csv(synapser::synGet("syn64024301")$path)
dot.df <- na.omit(all.degs[all.degs$Omics %in% c("RNA", "Protein"),])
# SHH is 0.9139077
topGenes <- unique(dot.df[dot.df$Significant & dot.df$N>=6 & dot.df$Spearman.est >= 0.9,]$Feature) # 83
topGenes <- unique(dot.df[dot.df$Significant & dot.df$N>=6 & dot.df$Spearman.est >= 0.91,]$Feature) # 70
#topGenes <- unique(dot.df[dot.df$Significant & dot.df$N>=6 & dot.df$Spearman.est >= 0.92,]$Feature) # 58
top.mean.df <- mean.dot.df[mean.dot.df$feature %in% topGenes,]
geneOrder <- top.mean.df[order(top.mean.df$mean_rankVal),]$feature
dot.df <- all.degs[all.degs$Feature %in% topGenes,]
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
maxLogFDR <- max(dot.df$minusLogFDR)
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = Omics, y = Feature, color = Spearman.est,
    size = minusLogFDR
  )
) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,5)) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="blue",high="red") +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Gene",
    color = "Spearman Correlation", size = "-log(adjusted p-value)"
  ) + 
  geom_point(data = subset(dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Correlated_genes_sig_minN6_Spearman91_dotPlot.pdf", dot.plot, width=4, height=12)

dot.df <- na.omit(all.degs[all.degs$Omics %in% c("RNA", "Protein"),])
topGenes <- unique(dot.df[dot.df$Significant & dot.df$N>=6 & dot.df$Spearman.est <= -0.9,]$Feature) # 68
topGenes <- unique(dot.df[dot.df$Significant & dot.df$N>=6 & dot.df$Spearman.est <= -0.91,]$Feature) # 45
#topGenes <- unique(dot.df[dot.df$Significant & dot.df$N>=6 & dot.df$Spearman.est <= -0.891,]$Feature) # 69
#topGenes <- unique(dot.df[dot.df$Significant & dot.df$N>=6 & dot.df$Spearman.est >= 0.92,]$Feature) # 58
top.mean.df <- mean.dot.df[mean.dot.df$feature %in% topGenes,]
geneOrder <- top.mean.df[order(top.mean.df$mean_rankVal),]$feature
dot.df <- all.degs[all.degs$Feature %in% topGenes,]
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.df$Omics <- factor(dot.df$Omics, levels=c("RNA", "Protein"))
maxLogFDR <- max(dot.df$minusLogFDR)
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = Omics, y = Feature, color = Spearman.est,
    size = minusLogFDR
  )
) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,5)) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="blue",high="red") +
  #viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Gene",
    color = "Spearman Correlation", size = "-log(adjusted p-value)"
  ) + 
  geom_point(data = subset(dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Correlated_genes_sig_minN6_negSpearman91_dotPlot.pdf", dot.plot, width=4, height=12)

9/12
dot.df <- na.omit(all.degs[all.degs$Omics %in% c("RNA", "Protein"),])
topGenes <- unique(dot.df[dot.df$Significant & dot.df$N>=9,]$Feature) # 722
topGenes <- unique(dot.df[dot.df$Significant & dot.df$N>=6 & abs(dot.df$Spearman.est) >= 0.9,]$Feature) # 151
topGenes <- unique(dot.df[dot.df$Significant & dot.df$N>=9 & abs(dot.df$Spearman.est) >= 0.9,]$Feature) # 117
top.mean.df <- mean.dot.df[mean.dot.df$feature %in% topGenes,]
geneOrder <- top.mean.df[order(top.mean.df$mean_rankVal),]$feature
dot.df <- all.degs[all.degs$Feature %in% topGenes,]
if (any(dot.df$Spearman.q == 0)) {
  dot.df[dot.df$Spearman.q == 0,]$Spearman.q <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"Spearman.p"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"Spearman.q"], base = 10)
dot.df$Omics <- factor(dot.df$Omics, levels=c("Protein", "RNA"))
maxLogFDR <- max(dot.df$minusLogFDR)
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = Omics, y = Feature, color = Spearman.est,
    size = minusLogFDR
  )
) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,5)) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="blue",high="red", limits=c(-1, 1)) +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Gene",
    color = "Spearman rho", size = "-log(adj. p)"
  ) + coord_flip() +
  geom_point(data = subset(dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Correlated_genes_sig_minN9_absSpearman90_dotPlot.pdf", dot.plot, width=18, height=4)
# no overlap

#### omics-specific dot plots ####
omics <- c("RNA", "Protein", "Phospho")
maxAbsEst <- max(na.omit(abs(all.degs[all.degs$Omics %in% omics,]$Spearman.est)))
maxAbsEst <- 1
all.degs <- na.omit(all.degs)
if (any(all.degs$Spearman.q == 0)) {
  all.degs[all.degs$Spearman.q == 0,]$Spearman.q <- 0.0001
}
all.degs$minusLogP <- -log(all.degs[,"Spearman.p"], base = 10)
all.degs$minusLogFDR <- -log(all.degs[,"Spearman.q"], base = 10)
maxLogFDR <- max(na.omit(all.degs[all.degs$Omics %in% omics,]$minusLogFDR))
library(patchwork); library(ggplot2)
gsea.dot.plots <- NULL
for (i in omics) {
  nGenes <- length(unique(na.omit(all.degs[all.degs$Omics == i,])$Feature))
  nCorrGenes <- length(unique(na.omit(all.degs[all.degs$Significant & 
                                        all.degs$Omics == i,])$Feature))
  
  # topGenes <- unique(na.omit(all.degs[all.degs$Significant & all.degs$N >= 9 &
  #                                       abs(all.degs$Spearman.est) >= 0.9 &
  #                                            all.degs$Omics == i,])$Feature)
  temp.dot.df <- all.degs[all.degs$Significant & all.degs$N >= 9 &
                         all.degs$Omics == i,] %>% slice_max(abs(Spearman.est), n = 50)
  topGenes <- temp.dot.df$Feature
  nTopGenes <- length(topGenes)
  geneOrder <- temp.dot.df[order(temp.dot.df$Spearman.est),]$Feature
  plot.annot <- paste0(i, "\n(", nCorrGenes, " / ", nGenes, " correlated)")

  # geneFace <- rep("plain", length(geneOrder))
  # names(geneFace) <- geneOrder
  # MYC.pathways <- geneOrder[grepl("MYC", geneOrder)]
  # SHH.pathways <- geneOrder[grepl("SHH", geneOrder) | grepl("HEDGEHOG", geneOrder)]
  # poi <- c(MYC.pathways, SHH.pathways)
  # geneFace[poi] <- "bold"
  
  dot.plot <- ggplot2::ggplot(temp.dot.df,
    ggplot2::aes(
      x = Omics, y = Feature, color = Spearman.est,
      size = minusLogFDR
    )
  ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = geneOrder) +
    scale_color_gradient2(low="blue",high="red", limits=c(-maxAbsEst, maxAbsEst)) +
    #viridis::scale_color_viridis() +
    theme_classic() +
    ggplot2::labs(
      x = as.character(i),
      color = "Spearman rho", size = "-log(adj. p)"
    ) + 
    #theme(axis.text.y=element_text(face=geneFace)) +
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
gsea.dot.plots <- gsea.dot.plots + plot_layout(guides = 'collect')
gsea.dot.plots
#ggplot2::ggsave("Enriched_geneSets_dotPlot_patchworkOmics_boldMYCandSHH.pdf", gsea.dot.plots, width=9, height=9)
ggplot2::ggsave("CorrelatedFeatures_dotPlot_patchworkOmics_minN9_sliceMaxAbsSpearman50.pdf", gsea.dot.plots, width=14, height=5)

# any overlapping genes?
temp.dot.df <- all.degs[all.degs$Significant & all.degs$N >= 9 &
                          all.degs$Omics == "RNA",] %>% slice_max(abs(Spearman.est), n = 50)
topGenesRNA <- temp.dot.df$Feature
temp.dot.df <- all.degs[all.degs$Significant & all.degs$N >= 9 &
                          all.degs$Omics == "Protein",] %>% slice_max(abs(Spearman.est), n = 50)
topGenesProt <- temp.dot.df$Feature
topGenesRNAProt <- topGenesProt[topGenesProt %in% topGenesRNA] # 0

temp.dot.df <- all.degs[all.degs$Significant & all.degs$N >= 9 &
                          all.degs$Omics == "Phospho",] %>% slice_max(abs(Spearman.est), n = 50)
topGenesPhos <- unique(sub("-.*","",temp.dot.df$Feature)) # 60
topGenesRNAPhos <- topGenesPhos[topGenesPhos %in% topGenesRNA] # 0
topGenesRNAPhos <- topGenesPhos[topGenesPhos %in% topGenesProt] # 3: MON1B, SPTBN1, WWTR1

# if just sig filter per omics:
# RNA (917 / 15773 correlated) 
# Protein (705 / 9013 correlated) 
# Phospho (3358 / 29545 correlated) 

# if also N9 filter:
# RNA (19 / 15773 correlated) 
# Protein (703 / 9013 correlated) 
# Phospho (3357 / 29545 correlated) 

# if sig, N9min, absSpearman.est0.9 min:
# RNA with 7 top features
# Protein with 110 top features
# Phospho with 448 top features

# if sig, N9min and pick top 50 by abs Spearman est:
# RNA (917 / 15773 correlated) 
# RNA with 19 top features
# Protein (705 / 9013 correlated) 
# Protein with 56 top features
# Phospho (3358 / 29545 correlated) 
# Phospho with 69 top features

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
                        all.degs$Omics %in% omics,]$N) # 12
library(patchwork); library(ggplot2)
gsea.dot.plots <- NULL
for (i in omics) {
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
  geneOrder <- temp.dot.df[order(temp.dot.df$Spearman.est),]$Feature
  plot.annot <- paste0(i, "\n(", nCorrGenes, " / ", nGenes, " correlated)")
  
  # geneFace <- rep("plain", length(geneOrder))
  # names(geneFace) <- geneOrder
  # MYC.pathways <- geneOrder[grepl("MYC", geneOrder)]
  # SHH.pathways <- geneOrder[grepl("SHH", geneOrder) | grepl("HEDGEHOG", geneOrder)]
  # poi <- c(MYC.pathways, SHH.pathways)
  # geneFace[poi] <- "bold"
  
  dot.plot <- ggplot2::ggplot(temp.dot.df,
                              ggplot2::aes(
                                y = Feature, x = Spearman.est, fill = N
                              )
  ) + #scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
    ggplot2::geom_col() + 
    ggplot2::scale_y_discrete(limits = geneOrder) +
    scale_fill_gradient2(#low="blue",high="red", 
      limits=c(6, maxN)) +
    #viridis::scale_color_viridis() +
    theme_classic() +
    ggplot2::labs(
      x = "Spearman rho",
      color = "Spearman rho"
    ) + 
    #theme(axis.text.y=element_text(face=geneFace)) +
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
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4/guides_build_mod.R")
gsea.dot.plots <- gsea.dot.plots + plot_layout(guides = 'collect'#, axis_titles = 'collect_y'
                                               )
gsea.dot.plots
#ggplot2::ggsave("Enriched_geneSets_dotPlot_patchworkOmics_boldMYCandSHH.pdf", gsea.dot.plots, width=9, height=9)
ggplot2::ggsave("CorrelatedFeatures_dotPlot_patchworkOmics_minN6_sliceMaxAbsSpearman50.pdf", gsea.dot.plots, width=14, height=7)
ggplot2::ggsave("CorrelatedFeatures_dotPlot_patchworkOmics_minN6_sliceMaxAbsSpearman50_colorN.pdf", gsea.dot.plots, width=14, height=7)

# any overlapping genes?
temp.dot.df <- all.degs[all.degs$Significant & all.degs$N >= 6 &
                          all.degs$Omics == "RNA",] %>% slice_max(abs(Spearman.est), n = 50)
topGenesRNA <- temp.dot.df$Feature
temp.dot.df <- all.degs[all.degs$Significant & all.degs$N >= 6 &
                          all.degs$Omics == "Protein",] %>% slice_max(abs(Spearman.est), n = 50)
topGenesProt <- temp.dot.df$Feature
topGenesRNAProt <- topGenesProt[topGenesProt %in% topGenesRNA] # 0

temp.dot.df <- all.degs[all.degs$Significant & all.degs$N >= 6 &
                          all.degs$Omics == "Phospho",] %>% slice_max(abs(Spearman.est), n = 50)
topGenesPhos <- unique(sub("-.*","",temp.dot.df$Feature)) # 60
topGenesRNAPhos <- topGenesPhos[topGenesPhos %in% topGenesRNA] # 0
topGenesProtPhos <- topGenesPhos[topGenesPhos %in% topGenesProt] # 3: MON1B, SPTBN1, WWTR1

# if sig, N6min:
# RNA with 50 top features
# Protein with 58 top features
# Phospho with 70 top features
