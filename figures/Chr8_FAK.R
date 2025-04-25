# FAK in Chr8q
# Author: Belinda B. Garana

# make bar plot of survival p-values
survival.p <- readxl::read_excel("~/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/chr8_crispr_target_survival.xlsx", sheet=3)
survival.p <- survival.p[survival.p$Groups=="Amp vs diploid",]
survival.p$Significance <- "p > 0.05"
survival.p[survival.p$p < 0.05,]$Significance <- "p < 0.05"
survival.p$Significance <- factor(survival.p$Significance, levels=c("p < 0.05", "p > 0.05"))
geneOrder <- survival.p[order(survival.p$p),]$Gene
ggplot2::ggplot(survival.p, aes(x=Gene, y=-log10(p), fill=Significance)) + geom_bar(stat="identity") +
  theme_classic(base_size=12) + scale_fill_manual(values=c("red","grey"), 
                                                  breaks=c("p < 0.05", "p > 0.05")) +
  ylab("-Log(p-value)") + ggtitle("Association with\nOverall Survival (N=236)") +
  scale_x_discrete(limits=geneOrder) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        plot.title=element_text(hjust=0.5))
ggsave("Association_w_OS_fillSig_2025-04-24.pdf",width=2.7, height=2.1)

# load human kinase-substrate database
ksdb <- read.csv("~/OneDrive - PNNL/Documents/ksdb_human_20231101.csv")

# identify kinases which phosphorylate PTK2
kin <- unique(ksdb[ksdb$SUB_GENE == "PTK2",]$GENE)

# identify substrates & sites phosphorylated by PTK2
sub <- unique(ksdb[ksdb$GENE == "PTK2",]$SUB_GENE)
sites <- ksdb[ksdb$GENE == "PTK2",]
sites$site_PNNL <- paste0(sites$SUB_GENE, "-", sites$SUB_MOD_RSD, tolower(substr(sites$SUB_MOD_RSD,1,1)))
subsites <- unique(sites$site_PNNL)

# load differential expression data
diffexp <- read.csv("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/Chr8_quant/Spearman/Compiled_results/Differential_expression/All_differential_expression_results.csv")
diffexp$Gene <- diffexp$Feature
diffexp[diffexp$Omics == "Phospho",]$Gene <- sub("\\-.*", "", diffexp[diffexp$Omics == "Phospho",]$Gene)

# heatmap of kinases which phosphorylate PTK2
diffexp.kin <- diffexp[diffexp$Gene %in% kin,] # 0; now 117?
sig.diffexp.kin <- diffexp.kin[diffexp.kin$Spearman.q <= 0.05,] # 7

# heatmap of substrates which are phosphorylated by PTK2
diffexp.sub <- diffexp[diffexp$Gene %in% sub,] # 127
sig.diffexp.sub <- diffexp.sub[diffexp.sub$Spearman.q <= 0.05,] # 21
sig.list <- list()
temp.omics <- unique(sig.diffexp.sub$Omics)
temp.omics <- temp.omics[temp.omics != "DNA"] # removing copy number because the only result there was PTK2 itself
for (i in 1:length(temp.omics)) {
  sig.list[[temp.omics[i]]] <- sig.diffexp.sub[sig.diffexp.sub$Omics == temp.omics[i],]
}
library(plyr);library(dplyr); library(ggplot2)
sig.compiled <- compile_mCorr(sig.list)
compiled.DEG.files <- list("Differential_expression_results_PTK2_substrates.csv" = sig.compiled$results,
                           "Compiled_differential_expression_results_PTK2_substrates.csv" = sig.compiled$mean.results,
                           "Differential_expression_venn_diagram_PTK2_substrates.pdf" = sig.compiled$venn.diagram,
                           "Differential_expression_dot_plot_PTK2_substrates.pdf" = sig.compiled$dot.plot,
                           "Differential_expression_correlations_PTK2_substrates.csv" = sig.compiled$corr,
                           "Differential_expression_correlation_matrix_PTK2_substrates.pdf" = sig.compiled$corr.matrix)
dir.create("Chr8_FAK")
setwd("Chr8_FAK")
save_to_synapse(compiled.DEG.files, dot.scale=1)
mean.DEG.df <- plyr::ddply(na.omit(sig.diffexp.sub), .(Feature), summarize,
                                           mean_rankVal = mean(Spearman.est),
                                           sd_rankVal = sd(Spearman.est),
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
mean.DEG.df <- dplyr::arrange(mean.DEG.df, desc(mean_rankVal))
sig.diffexp.sub$qVal <- sig.diffexp.sub$Spearman.q
if (any(sig.diffexp.sub$qVal == 0)) {
  sig.diffexp.sub[sig.diffexp.sub$qVal == 0,]$qVal <- 0.0001
}

bg.theme <- ggplot2::theme(
  legend.background = element_rect(), legend.position = "top",
  legend.text = element_text(size = 14),
  legend.key = element_blank(),
  legend.title = element_text(size = 16),
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16, colour = "black"),
  axis.title.y = element_text(size = 20),
  axis.text.y = element_text(size = 16, colour = "black"),
  plot.title = element_text(
    lineheight = .8, face = "bold", size = 36
  ),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(fill = NA),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks.x = element_line(colour = "black"),
  axis.ticks.y = element_line(colour = "black")
)


dot.plot <- ggplot2::ggplot(
  na.omit(sig.diffexp.sub),
  ggplot2::aes(
    x = Omics, y = Feature, color = Spearman.est,
    size = -log10(qVal)
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = mean.DEG.df$Feature) +
  viridis::scale_color_viridis() +
  bg.theme +
  ggplot2::labs(
    x = "Input Data",
    y = "Feature Set",
    color = "Correlation Estimate", size = "-log(FDR)"
  )
ggsave("PTK2_substrate_dot_plot.pdf", dot.plot, width = 7, height = 7)

# heatmap of substrate sites which are phosphorylated by PTK2
diffexp.sub <- diffexp[diffexp$Feature %in% subsites,] # 1
sig.diffexp.sub <- diffexp.sub[diffexp.sub$Spearman.q <= 0.05,] # 0
sig.list <- list()
temp.omics <- unique(sig.diffexp.sub$Omics)
temp.omics <- temp.omics[temp.omics != "DNA"] # removing copy number because the only result there was PTK2 itself
for (i in 1:length(temp.omics)) {
  sig.list[[temp.omics[i]]] <- sig.diffexp.sub[sig.diffexp.sub$Omics == temp.omics[i],]
}
library(plyr);library(dplyr); library(ggplot2)
sig.compiled <- compile_mCorr(sig.list)
compiled.DEG.files <- list("Differential_expression_results_PTK2_substrates.csv" = sig.compiled$results,
                           "Compiled_differential_expression_results_PTK2_substrates.csv" = sig.compiled$mean.results,
                           "Differential_expression_venn_diagram_PTK2_substrates.pdf" = sig.compiled$venn.diagram,
                           "Differential_expression_dot_plot_PTK2_substrates.pdf" = sig.compiled$dot.plot,
                           "Differential_expression_correlations_PTK2_substrates.csv" = sig.compiled$corr,
                           "Differential_expression_correlation_matrix_PTK2_substrates.pdf" = sig.compiled$corr.matrix)
dir.create("Chr8_FAK_sites")
setwd("Chr8_FAK_sites")
save_to_synapse(compiled.DEG.files, dot.scale=1)
mean.DEG.df <- plyr::ddply(na.omit(sig.diffexp.sub), .(Feature), summarize,
                           mean_rankVal = mean(Spearman.est),
                           sd_rankVal = sd(Spearman.est),
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
mean.DEG.df <- dplyr::arrange(mean.DEG.df, desc(mean_rankVal))
sig.diffexp.sub$qVal <- sig.diffexp.sub$Spearman.q
if (any(sig.diffexp.sub$qVal == 0)) {
  sig.diffexp.sub[sig.diffexp.sub$qVal == 0,]$qVal <- 0.0001
}

bg.theme <- ggplot2::theme(
  legend.background = element_rect(), legend.position = "top",
  legend.text = element_text(size = 14),
  legend.key = element_blank(),
  legend.title = element_text(size = 16),
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16, colour = "black"),
  axis.title.y = element_text(size = 20),
  axis.text.y = element_text(size = 16, colour = "black"),
  plot.title = element_text(
    lineheight = .8, face = "bold", size = 36
  ),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(fill = NA),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks.x = element_line(colour = "black"),
  axis.ticks.y = element_line(colour = "black")
)


dot.plot <- ggplot2::ggplot(
  na.omit(sig.diffexp.sub),
  ggplot2::aes(
    x = Omics, y = Feature, color = Spearman.est,
    size = -log10(qVal)
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = mean.DEG.df$Feature) +
  viridis::scale_color_viridis() +
  bg.theme +
  ggplot2::labs(
    x = "Input Data",
    y = "Feature Set",
    color = "Correlation Estimate", size = "-log(FDR)"
  )
ggsave("PTK2_substrate_site_dot_plot.pdf", dot.plot, width = 7, height = 7)

#### Chr8 plot of CRISPR pVals ####
# BiocManager::install("ggbio")
# library(ggbio)
# BiocManager::install("GenomicRanges")
# #library(GenomicRanges)
# # following example from: http://www.sthda.com/english/wiki/ggbio-visualize-genomic-data
# ideo <- ggbio::Ideogram(genome = "hg19")
# ideo
# # highlight each gene with p < 0.05 in crispr screen
# positions <- readxl::read_xlsx("58_genes_chr8 locations.xlsx")
# goi <- positions$symbol
# for (i in 1:length(goi)) {
#   ideo + xlim(GenomicRanges::GRanges("chr8", )) 
# }

# prep data
positions <- readxl::read_xlsx("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/data/58_genes_chr8 locations.xlsx")
crispr <- read.table("~/OneDrive - PNNL/Documents/MPNST/Chr8/Chr8_fitness.gene_summary.txt", sep="\t", header = TRUE)
crispr$minusLogP <- as.numeric(-log(crispr$neg.p.value, base=10))
crispr$Significant <- FALSE
crispr[crispr$neg.p.value<= 0.05,]$Significant <- TRUE
colnames(crispr)[1] <- "symbol"
df <- merge(crispr, positions)

# using RIdeogram
library(RIdeogram)
# following vignette: https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html
data("human_karyotype", package="RIdeogram")
human8 <- as.data.frame(human_karyotype[8,])
ideo.df <- df[,c("chr", "start", "end", "minusLogP")]
colnames(ideo.df) <- c("Chr", "Start", "End", "Value")
label.df <- df[,c("symbol", "chr", "start", "end")]
#label.df$Type <- ""
#label.df[label.df$symbol=="PTK2",]$Type <- "PTK2"
label.df$Type <- label.df$symbol
label.df$Shape <- "box"
label.df$color <- "black"
label.df[label.df$symbol=="PTK2",]$color <- "red"
label.df <- label.df[,c("Type", "Shape", "chr", "start", "end", "color")]
#label.df <- label.df[,c("symbol", "Shape", "chr", "start", "end", "color")]
colnames(label.df)[3:5] <- c("Chr", "Start", "End")
#label.df2 <- as.data.frame(label.df[label.df$Type == "PTK2",])
#RIdeogram::ideogram(karyotype = human_karyotype, overlaid = ideo.df, label = label.df2, label_type = "marker", colorset1 = c("pink", "red"))
RIdeogram::ideogram(karyotype = human_karyotype, overlaid = ideo.df, label = label.df, label_type = "marker", colorset1 = c("pink", "red"))
convertSVG("chromosome.svg", device="pdf")
# 
# # create plot just to scale text position
# plot(y = label.df$Start, x = rep(5, nrow(label.df)))
# text(x= rep(5, nrow(label.df)), y = label.df$Start, label = label.df$Type)

#RIdeogram::ideogram(karyotype = human_karyotype, overlaid = ideo.df, label = label.df2, label_type = "marker")
#convertSVG("chromosome.svg", device="pdf")

# using karyoploteR
df$y <- df$minusLogP/5.5 # scale to get y < 1
#df <- df[,c("chr", "start", "end","y")]
df$y0 <- 0
df$y1 <- df$y
df$chr <- "chr8"
df$x <- 0
df$x_jitter <- df$x
for (i in 1:nrow(df)) {
  df$x[i] <- mean(df$start[i], df$end[i])
}
# labels too close together after 1.41E8
df[df$x > 141000000,]$x_jitter <- jitter(df[df$x > 141000000,]$x, factor = 50)
df$y_jitter <- df$y
df[df$x > 141000000,]$y_jitter <- jitter(df[df$x > 141000000,]$y, factor = 1000)
#BiocManager::install("ChIPpeakAnno")
library(ChIPpeakAnno)
dd <- ChIPpeakAnno::toGRanges(df)
# duplicated or NA names found. Rename all the names by numbers.

#BiocManager::install("karyoploteR")
library(karyoploteR)
# kp <- karyoploteR::plotKaryotype(genome="hg38", chromosomes="chr8")
# karyoploteR::kpAxis(kp, labels = c(0, 0.5*5.5, 1*5.5))
# karyoploteR::kpPoints(kp, data=dd, col="red", border = "red")
# karyoploteR::kpPlotMarkers(kp, data=dd, labels=dd$symbol, y=dd$y, line.color=NULL, text.orientation = "horizontal")
# resource: https://bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html
kp <- karyoploteR::plotKaryotype(genome="hg38", chromosomes="chr8")
karyoploteR::kpAxis(kp, labels = c(0, 0.5*5.5, 1*5.5))
karyoploteR::kpPlotMarkers(kp, data=dd, labels=dd$symbol, y=dd$y, line.color="grey")
#karyoploteR::kpPlotMarkers(kp, data=dd, labels=dd$symbol, y=dd$y_jitter, line.color="grey", r1 = 1, x=dd$x_jitter)
karyoploteR::kpBars(kp, data=dd, col="red", border = "red")
#karyoploteR::kpPlotMarkers(kp, data=dd, labels=dd$symbol, y=dd$y_jitter, line.color="grey", r1 = 1, x=dd$x_jitter)
#karyoploteR::kpPlotMarkers(kp, data=dd, labels=dd$symbol, y=dd$y, line.color="grey")

# chr8 heatmap plot for Fig1
#BiocManager::install("ChromHeatMap")
library(ChromHeatMap)
crispr <- read.table("~/OneDrive - PNNL/Documents/MPNST/Chr8/Chr8_fitness.gene_summary.txt", sep="\t", header = TRUE)
crispr$minusLogP <- as.numeric(-log(crispr$neg.p.value, base=10))
crispr$minusLogFDR <- as.numeric(-log(crispr$neg.fdr, base=10) )
rownames(crispr) <- crispr$id
crispr$Significant <- FALSE
crispr[crispr$neg.p.value<= 0.05,]$Significant <- TRUE
sig.crispr <- crispr[crispr$neg.p.value <= 0.05,c("id", "minusLogP")]

heatmap.input <- as.matrix(crispr[,c("id","minusLogP","minusLogFDR","Significant")])
colnames(heatmap.input)[1] <- "SYMBOL"

# change from symbols to ids
#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
db.info <- AnnotationDbi::select(org.Hs.eg.db, columns = c("SYMBOL","ENTREZID"), keys=rownames(heatmap.input), keytype="SYMBOL")
heatmap.input <- na.omit(merge(heatmap.input, db.info, by="SYMBOL"))
rownames(heatmap.input) <- heatmap.input$ENTREZID
heatmap.input <- heatmap.input[,c("minusLogP", "minusLogFDR", "Significant")]
heatmap.input1 <- as.matrix(heatmap.input[,c("minusLogP", "minusLogFDR")])
heatmap.input1[,1] <- as.numeric(heatmap.input1[,1])
heatmap.input1[,2] <- as.numeric(heatmap.input1[,2])
#rownames(heatmap.input) <- rownames(crispr)
pos <- ChromHeatMap::makeChrStrandData(heatmap.input1, lib="org.Hs.eg.db")
ChromHeatMap::plotChrMap(pos, 8, strands="both")
## when input cols were minusLogP and Significant:
# Using interval argument set to 290278
# Error in plot.new() : figure margins too large
# In addition: There were 50 or more warnings (use warnings() to see the first 50)

## when using input cols minusLogP and minusLogFDR OR minusLogP and minusLogFDR:
# Using interval argument set to 290278
# Error in .External.graphics(C_layout, num.rows, num.cols, mat, as.integer(num.figures),  : 
#                               invalid graphics state
#                             In addition: There were 50 or more warnings (use warnings() to see the first 50)

heatmap.input2 <- as.matrix(heatmap.input[heatmap.input$Significant == "TRUE",c("minusLogP", "minusLogFDR")])
heatmap.input2[,1] <- as.numeric(heatmap.input2[,1])
heatmap.input2[,2] <- as.numeric(heatmap.input2[,2])
pos2 <- ChromHeatMap::makeChrStrandData(heatmap.input2, lib="org.Hs.eg.db")
chrPlot2 <- ChromHeatMap::plotChrMap(pos2, 8, strands="both")
# Using interval argument set to 290278
# Error in plot.new() : figure margins too large
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
debug(ChromHeatMap::plotChrMap)
# before running line 19 in debug, replace values for minusLogP and minusLogFDR
# actually idk which values should go where

# try doing minusLogFDR instead of minusLogFDR because minusLogFDR has negative values
# still only have grey bars. maybe it's because inputting matrix instead of expressionSet
heatmap.set2 <- ExpressionSet(assayData=heatmap.input2)
pos2 <- ChromHeatMap::makeChrStrandData(heatmap.set2, lib="org.Hs.eg.db")
chrPlot2 <- ChromHeatMap::plotChrMap(pos2, 8, strands="both")
# no change