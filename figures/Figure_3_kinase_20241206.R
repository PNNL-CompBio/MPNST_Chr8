# chr8 MPNST: figure 3: kinase enrichment
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase")

synapser::synLogin()
synID <- synapser::synStore(synapser::Folder("Figure_4", parent = "syn64367769"))

#### 1. volcano plot of # of kinase enrichment results ####
all.gsea <- read.table(synapser::synGet("syn64374349")$path, sep="\t", header=TRUE)
plot.data <- all.gsea
FDR <- 0.25
plot.data$Significance <- paste0("FDR > ", FDR)
plot.data[plot.data$adjusted_p_value < FDR, ]$Significance <- paste0("FDR < ", FDR)
plot.data$Significance <- factor(plot.data$Significance,
                                 levels = c(paste0("FDR < ", FDR),
                                            paste0("FDR > ", FDR)))
if (any(plot.data$p_value == 0)) {
  plot.data[plot.data$p_value == 0,]$p_value <- 0.0001
}
limit.x <- ceiling(max(abs(plot.data$enrichment_value_log2)))
limit.y <- ceiling(max(-log(plot.data$p_value, base=10)))
volc <- ggplot2::ggplot(data = plot.data, aes(
  x = enrichment_value_log2, y = -log(p_value, 10),
  color = Significance
)) +
  ggplot2::geom_point(size = 4) +
  # ggrepel::geom_text_repel(
  #   data = subset(plot.data, Significance == paste0("FDR > ", FDR)),
  #   mapping = aes(label = kinase, size = I(6)), nudge_y = 0.25
  # ) +
  ggplot2::scale_color_manual(
    values = c("red", "azure4"), name = "Significance",
    breaks = c(paste0("FDR < ", FDR), paste0("FDR > ", FDR))
  ) +
  ggplot2::xlim(-limit.x, limit.x) +
  ggplot2::ylim(0, limit.y) +
  ggplot2::xlab("log2(Enrichment Score)") +
  ggplot2::ylab("-log(p-value)") +
  ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey",
                      linewidth = 0.5) +
  ggplot2::theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_line(colour = "black", linewidth = 0.65),
    legend.text = element_text(size = 20),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 26, face = "bold"),
    panel.background = element_rect(
      fill = "white", colour = "white", linewidth = 0.5,
      linetype = "solid", color = "black"
    ), text = element_text(size = 20),
    legend.position = "bottom", legend.key = element_blank()
  )
volc

#synapser::synStore(synapser::File("Chr8_numberOfFeatures_directionGroup_barPlot_logScale_manualOrderV2.pdf", parent = synID))

#### 3. top up/dn diffexp gene sets ####
dot.df <- na.omit(plot.data[plot.data$Significance == paste0("FDR < ", FDR),]) # 12
topGenes <- unique(dot.df$kinase)
geneOrder <- dot.df[order(dot.df$enrichment_value_log2),]$kinase
dot.df$minusLogP <- -log(dot.df[,"p_value"], base = 10)
dot.df$adjusted_p_value_log10_abs <- -log(dot.df[,"adjusted_p_value"], base = 10)
dot.plot <- ggplot2::ggplot(
  na.omit(dot.df),
  ggplot2::aes(x="",y = kinase, color = enrichment_value_log2,
    size = adjusted_p_value_log10_abs
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="blue",high="red") +
  #viridis::scale_color_viridis() +
  theme_classic() +
  ggplot2::labs(
    x = "Omics Type",
    y = "Kinase",
    color = "log2(enrichment)", size = "-log(adj. p)"
  ) + theme(axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank()) +
  geom_point(data = subset(dot.df, Significance == paste0("FDR < ", FDR)), col = "black", stroke = 1.5, shape = 21)
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Enriched_kinases_dotPlot_v2.pdf", dot.plot, width=3, height=4)

#### correlated kinase targets for each kinase ####
site.files <- paste0(topGenes[topGenes != "LTK"], "_serineThreonineKinase_score-phosphoproteome-result-table.tsv")
site.files <- c(site.files, "LTK_tyrosineKinase_score-phosphoproteome-result-table.tsv")
site.scores <- data.frame()
for (i in topGenes) {
  # read site scores
  temp.file.name <- site.files[grepl(i, site.files)]
  if (length(temp.file.name) > 1) {
    stop("more than 1 filename detected for 1 kinase")
  } else {
    temp.file <- read.table(file.path("sites",temp.file.name), sep="\t", header = TRUE) 
  }
  temp.file$kinase <- i
  
  # combine with other site scores
  site.scores <- rbind(site.scores, temp.file)
}
# get gene symbols from UniProt
site.scores$Gene <- sub("_HUMAN.*", "", site.scores$protein)
site.scores$siteType <- tolower(substr(site.scores$position,1,1))
site.mapping <- data.frame()
for (i in 1:nrow(site.scores)) {
  SUB_SITE <- c()
  temp.prot <- site.scores$protein[i]
  download.file(paste0("https://rest.uniprot.org/uniprotkb/stream?fields=",
                       "accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2C",
                       "organism_name%2Clength&format=xlsx&query=%28%28id%3A",
                       temp.prot,"%29%29"), "uniprot_temp.xlsx")
  temp.results <- readxl::read_excel("uniprot_temp.xlsx")
  if (nrow(temp.results) > 0) {
    temp.genes <- strsplit(temp.results$`Gene Names`, " ")[[1]] 
  } else {
    temp.genes <- sub("_HUMAN.*", "", temp.prot)
  }
  SUB_SITE <- paste0(temp.genes, "-", site.scores$position[i], site.scores$siteType[i])
  temp.gene.df <- data.frame(SUB_SITE)
  temp.gene.df$kinase <- site.scores$kinase[i]
  site.mapping <- rbind(site.mapping, temp.gene.df)
}
write.csv(site.mapping, "site_mapping.csv", row.names = FALSE)
site.mapping <- read.csv("site_mapping.csv")

topGenes <- rev(geneOrder)
corr.result <- read.csv(synapser::synGet("syn64024300")$path)
corr.result <- na.omit(corr.result[corr.result$type=="Phospho" & corr.result$feature %in% site.mapping$SUB_SITE,]) # 171
if (any(corr.result$Spearman.q == 0)) {
  corr.result[corr.result$Spearman.q == 0,]$Spearman.q <- 0.0001
}
corr.result$minusLogFDR <- -log(corr.result$Spearman.q, base = 10)
maxAbsEst <- ceiling(max(abs(corr.result$Spearman.est)))
maxLogFDR <- max(abs(corr.result$minusLogFDR))
library(patchwork); library(ggplot2)
gsea.dot.plots <- NULL
for (i in 1:length(topGenes)) {
  temp.genes <- unique(site.mapping[site.mapping$kinase == topGenes[i],]$SUB_SITE)
  temp.gene.corr <- na.omit(corr.result[corr.result$feature %in% temp.genes & 
                                          corr.result$N>=6
                                        ,])
  if (nrow(temp.gene.corr) > 0) {
    nGenes <- length(unique(temp.gene.corr$feature))
    nDuplicateGenes <- nrow(temp.gene.corr) - nGenes
    temp.topGenes <- unique(temp.gene.corr[temp.gene.corr$sig,]$feature) # 48
    nTopGenes <- length(temp.topGenes)
    cat(topGenes[i], "has", nTopGenes, "correlated out of", nGenes, "with", nDuplicateGenes, "duplicates\n")
    geneOrder <- temp.gene.corr[order(temp.gene.corr$Spearman.est),]$feature
    
    dot.plot <- ggplot2::ggplot(
      temp.gene.corr,
      ggplot2::aes(
        x = Omics, y = feature, color = Spearman.est,
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
source("guides_build_mod.R")
gsea.dot.plot2 <- gsea.dot.plots + guide_area() + plot_layout(nrow = 2, guides='collect')
#ggplot2::ggsave("Correlated_kinaseSubsites_dotPlot_patchworkCollection_minN6_wide_v3.pdf", gsea.dot.plot2, width=16, height=10)
ggplot2::ggsave(paste0("kinaseSubsite_correlations_dotPlot_patchworkCollection_minN6_v2_",Sys.Date(),".pdf"), gsea.dot.plot2, width=12, height=10)
gsea.dot.plot2 <- gsea.dot.plots + plot_layout(nrow = 1, guides='collect')
ggplot2::ggsave(paste0("kinaseSubsite_correlations_dotPlot_patchworkCollection_minN6_wide_v3_",Sys.Date(),".pdf"), gsea.dot.plot2, width=18, height=6)
#ggplot2::ggsave(paste0("kinaseSubsite_correlations_dotPlot_patchworkCollection_minN6_v3_",Sys.Date(),".pdf"), gsea.dot.plot2, width=12, height=10)

