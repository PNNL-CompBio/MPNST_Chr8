# chr8 MPNST: figure 6
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4")

synapser::synLogin()

getCCLEprot <- function(){
  # get CCLE proteomics
    download.file("https://figshare.com/ndownloader/files/41466702", "proteomics.csv.gz")
    prot.df <- read.csv(gzfile("proteomics.csv.gz"),fileEncoding="UTF-16LE")
    
    allgenes = readr::read_csv("https://figshare.com/ndownloader/files/40576109")
    genes = allgenes|>
      dplyr::select(gene_symbol,entrez_id)|>
      dplyr::distinct()
    #genes <- genes[genes$gene_symbol %in% colnames(RNA.df)[2:ncol(RNA.df)], ]
    
    allsamples = readr::read_csv('https://figshare.com/ndownloader/files/40576103')
    CCLE.samples <- dplyr::distinct(allsamples[allsamples$id_source == "CCLE",
                                               c("other_id","improve_sample_id")])
    
    # merge prot.df with genes, samples to stop using improve IDs
    prot.df <- merge(prot.df, CCLE.samples)
    prot.df <- merge(prot.df, genes)
    prot.df$entrez_id <- NULL
    prot.df <- dplyr::distinct(prot.df)
    
    # convert to wide format for DMEA
    prot.df <- reshape2::dcast(prot.df, other_id ~ gene_symbol, mean,
                               value.var = "proteomics")
    colnames(prot.df)[1] <- "CCLE_ID"
    prot.df.noNA <- prot.df[, colSums(is.na(prot.df)) == 0] # 23304 gene names 
  }
  return(prot.df.noNA)
}

### redo just using significantly corr features
source("https://raw.githubusercontent.com/PNNL-CompBio/MPNST_Chr8/refs/heads/main/proteomics/panSEA_helper_20240913.R")

# load data
rna.corr <- read.csv(synapser::synGet("syn66226866")$path)
prot.corr <- read.csv(synapser::synGet("syn66224803")$path)
inputs <- list("RNA" = na.omit(rna.corr[rna.corr$Spearman.q <= 0.05,]), # 655 genes or 50 with min N6 / 17717
               "Protein" = na.omit(prot.corr[prot.corr$Spearman.q <= 0.05,])) # 208 genes / 9013

adh.RNA <- get_CCLE_RNA()
adh.RNA2 <- adh.RNA[,c("CCLE_ID",colnames(adh.RNA)[colnames(adh.RNA) %in% rna.corr$Gene])] # 4448 genes in 327 cell lines or 3601 with min N 6

temp.expr <- getCCLEprot()
sample.info <- read.csv("https://raw.githubusercontent.com/BelindaBGarana/DMEA/refs/heads/shiny-app/Inputs/CCLE_sample_info.csv")
temp.expr.adherent <- temp.expr[temp.expr$CCLE_ID %in% sample.info[tolower(sample.info$culture_type)=="adherent",]$CCLE_Name,]
adh.prot2 <- temp.expr.adherent[,c("CCLE_ID",colnames(temp.expr.adherent)[colnames(temp.expr.adherent) %in% prot.corr$Gene])] # 5035 proteins in 215 cell lines
soft.sarc.info <- sample.info[sample.info$lineage == "soft_tissue" & grepl("sarcoma", sample.info$lineage_subtype, ignore.case=TRUE),] # 49

other.sample.info <- sample.info[!(sample.info$CCLE_Name %in% temp.expr.adherent$CCLE_ID) &
                                   sample.info$CCLE_Name %in% temp.expr$CCLE_ID,]
adh.DMEA <- panSEA::mDMEA(expression = list(adh.RNA2, adh.prot2),
                              weights = inputs, types = c("RNA", "Protein"),
                              weight.values = rep("Spearman.est", 2),
                              ties = TRUE)
adh.DMEA.files <- list("DMEA_results.csv" = adh.DMEA$compiled.results$results,
                       "DMEA_results_compiled.csv" = adh.DMEA$compiled.results$mean.results,
                       "DMEA_venn_diagram.pdf" = adh.DMEA$compiled.results$venn.diagram,
                       "DMEA_dot_plot.pdf" = adh.DMEA$compiled.results$dot.plot,
                       "DMEA_correlation_matrix.csv" = adh.DMEA$compiled.results$corr,
                       "DMEA_correlation_matrix.pdf" = adh.DMEA$compiled.results$corr.matrix)
for (i in 1:length(inputs)) {
  temp.results <- adh.DMEA$all.results[[names(inputs)[i]]]
  temp.files <- list("DMEA_results.csv" = temp.results$result,
                     "DMEA_results_woShufflingTies.csv" = temp.results$result.w.ties,
                     "WV_results.csv" = temp.results$WV.scores,
                     "Unused_weights.csv" = temp.results$unused.weights,
                     "DMEA_correlations.csv" = temp.results$corr.result,
                     "DMEA_correlation_scatterPlots.pdf" = temp.results$corr.scatter.plots,
                     #"DMEA_mountainPlots.pdf" = temp.results$mtn.plots,
                     "DMEA_volcanoPlot.pdf" = temp.results$volcano.plot,
                     "DMEA_removed_sets.csv" = temp.results$removed.sets,
                     "DMEA_unannotated_drugs.csv" = temp.results$unannotated.drugs)
  adh.DMEA.files[[names(inputs)[i]]] <- temp.files
}
dir.create("DMEA")
setwd("DMEA")
save_to_synapse(adh.DMEA.files, "syn66295230")
saveRDS(adh.DMEA, "DMEA.rds")

source("https://raw.githubusercontent.com/PNNL-CompBio/MPNST_Chr8/refs/heads/main/figures/compile_mCorr.R")
compiled.drugCorr <- compile_mCorr(list("RNA" = adh.DMEA$all.results$RNA$corr.result,
                                        "Protein" = adh.DMEA$all.results$Protein$corr.result))
compiled.files <- list("DMEA_correlation_results.csv" = compiled.drugCorr$results,
                       "DMEA_correlation_results_compiled.csv" = compiled.drugCorr$mean.results,
                       "DMEA_correlation_venn_diagram.pdf" = compiled.drugCorr$venn.diagram,
                       "DMEA_correlation_dot_plot.pdf" = compiled.drugCorr$dot.plot,
                       "DMEA_correlation_correlation_matrix.csv" = compiled.drugCorr$corr,
                       "DMEA_correlation_correlation_matrix.pdf" = compiled.drugCorr$corr.matrix)

setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4")
setwd("DMEA")
dir.create("correlations")
setwd("correlations")
save_to_synapse(compiled.files, "syn66295272")

# sarcoma DMEA?
soft.sarc.RNA <- adh.RNA2[adh.RNA2$CCLE_ID %in% soft.sarc.info$CCLE_Name,] # 7 cell lines
soft.sarc.prot <- adh.prot2[adh.prot2$CCLE_ID %in% soft.sarc.info$CCLE_Name,] # 3 cell lines

#### 1. top up/dn enriched MOA sets ####
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4")
# drug.corr <- read.csv(synapser::synGet("syn66295273")$path)
# mean.drugs <- read.csv(synapser::synGet("syn66295274")$path)
# moa.results <- read.csv(synapser::synGet("syn66295241")$path)
# mean.moa <- read.csv(synapser::synGet("syn66295242")$path)

drug.corr <- read.csv(synapser::synGet("syn68747773")$path)
mean.drugs <- read.csv(synapser::synGet("syn68747774")$path)
moa.results <- read.csv(synapser::synGet("syn68747739")$path)
mean.moa <- read.csv(synapser::synGet("syn68747740")$path)
gsea.rna.prot <- moa.results[moa.results$type %in% c("RNA", "Protein"),]
gsea.rna.prot$type <- factor(gsea.rna.prot$type, levels=c("RNA", "Protein"))
gsea.rna.prot$Significant <- FALSE
gsea.rna.prot[gsea.rna.prot$p_value <= 0.05 & gsea.rna.prot$FDR_q_value <= 0.25,]$Significant <- TRUE
gsea.rna.prot$Direction <- "Negative"
gsea.rna.prot[gsea.rna.prot$NES > 0,]$Direction <- "Positive"
gsea.rna.prot$Direction <- factor(gsea.rna.prot$Direction, levels=c("Positive", "Negative"))
gsea.rna.prot$Toxicity <- "Chr8q-amplified"
gsea.rna.prot[gsea.rna.prot$NES > 0,]$Toxicity <- "Non-amplified"

dot.df <- gsea.rna.prot[gsea.rna.prot$Significant,] # 19 - note glucocorticoid receptor agonist has FDR=0.25 in protein so use 'Significant' to include it
nrow(dot.df[dot.df$type=="RNA",]) # 7
nrow(dot.df[dot.df$type=="Protein",]) # 12
hist(dot.df$NES,breaks=20)
min(abs(dot.df$NES)) # 1.48
topGenes <- unique(dot.df$Drug_set) # 17

#mean.moa <- mean.moa[mean.moa$Drug_set %in% topGenes,]
#geneOrder <- mean.moa[order(-mean.moa$mean_NES),]$Drug_set
mean.moa <- plyr::ddply(gsea.rna.prot, .(Drug_set), summarize,
                        maxAbsNES = max(abs(NES[Significant]), na.rm=TRUE),
                        maxNES = max(NES[Significant], na.rm=TRUE),
                        minNES = min(NES[Significant], na.rm=TRUE))
mean.moa <- mean.moa[mean.moa$maxAbsNES != "-Inf",]
geneOrder <- mean.moa[order(-mean.moa$maxNES),]$Drug_set
# if not sig in both or sig only in RNA, then sig only in Protein; if needed, could use "bold.italic" face style
#moi <- na.omit(unique(mean.moa[grepl("RNA, Protein", mean.moa$sig_types),]$Drug_set)) # HDAC inhibitor
#moi <- c("HDAC inhibitor", "glucocorticoid receptor agonist")
#geneFace <- rep("plain", length(geneOrder))
#names(geneFace) <- geneOrder
#geneFace[moi] <- "bold"
dot.df <- na.omit(moa.results[moa.results$Drug_set %in% topGenes &
                     moa.results$type %in% c("RNA", "Protein"),])
if (any(dot.df$FDR_q_value == 0)) {
  dot.df[dot.df$FDR_q_value == 0,]$FDR_q_value <- 0.0001
}
dot.df$minusLogP <- -log(dot.df[,"p_value"], base = 10)
dot.df$minusLogFDR <- -log(dot.df[,"FDR_q_value"], base = 10)
dot.df$type <- factor(dot.df$type, levels=c("RNA", "Protein"))
dot.df$Significant <- FALSE
dot.df[dot.df$p_value <= 0.05 & dot.df$FDR_q_value <= 0.25,]$Significant <- TRUE
absMaxNES <- max(abs(dot.df$NES))
dot.plot <- ggplot2::ggplot(
  dot.df,
  ggplot2::aes(
    x = type, y = Drug_set, color = NES,
    size = minusLogFDR
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_y_discrete(limits = geneOrder) +
  scale_color_gradient2(low="red",high="blue", mid="grey", limits=c(-absMaxNES, absMaxNES)) +
  theme_classic(base_size=12) +
  ggplot2::labs(
    x = "Omics Type",
    y = "Drug Mechanism",
    color = "NES", size = "-log(FDR)"
  ) + geom_point(data = subset(dot.df, Significant), col = "black", stroke = 1.5, shape = 21) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()#, axis.text = element_text(size=16)
        ) + 
  ggtitle(paste0("Drug Mechanisms (", length(topGenes), " / ", length(unique(moa.results$Drug_set)), " enriched)")) +
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
  theme(#axis.text.y=element_text(face=geneFace), 
    axis.text.x=element_text(angle=45, vjust=1, hjust=1))
dot.plot
# most are only in prot and not in RNA
ggplot2::ggsave("Enriched_drugSets_dotPlot_v2.pdf", dot.plot, width=6, height=6)

#### 2. waterfall plot of drug corr ####
drug.info <- read.csv("https://raw.githubusercontent.com/BelindaBGarana/DMEA/refs/heads/shiny-app/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv")
omics <- c("RNA", "Protein")
maxAbsEst <- max(na.omit(abs(drug.corr[drug.corr$type %in% omics,]$Pearson.est)))
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
drug.corr$sig <- FALSE
drug.corr[drug.corr$Pearson.q <= 0.05,]$sig <- TRUE # Synapse "sig" was mistakenly based on Spearman
drug.corr.wInfo <- merge(red.drug.info, drug.corr[drug.corr$sig,], by="Drug", all.y = TRUE)

sig.moas <- unique(moa.results[moa.results$type %in% omics & moa.results$sig,]$Drug_set) # 24

maxAbsEst <- max(na.omit(abs(drug.corr.wInfo[drug.corr.wInfo$type %in% omics,]$Pearson.est)))
maxLogFDR <- max(na.omit(drug.corr.wInfo[drug.corr.wInfo$type %in% omics,]$minusLogFDR))
library(patchwork); library(ggplot2)
#moaOrder <- mean.moa[order(mean.moa$mean_NES),]$Drug_set # 24

# simplify drug moa annotations as needed
drug.corr.wInfo$Mechanism <- "Other"
drug.corr.wInfo[drug.corr.wInfo$Drug == "danusertib",]$moa <- "Aurora kinase inhibitor" # also labeled as growth factor receptor inhibitor
#drug.corr.wInfo[drug.corr.wInfo$Drug == "daunorubicin",]$moa <- "topoisomerase inhibitor" # also labeled as RNA synthesis inhibitor
#drug.corr.wInfo[drug.corr.wInfo$Drug == "dovitinib",]$moa <- "EGFR inhibitor" # also labeled as FGFR inhibitor, FLT3 inhibitor, PDGFR tyrosine kinase receptor inhibitor, VEGFR inhibitor
#drug.corr.wInfo[drug.corr.wInfo$Drug == "evodiamine",]$moa <- "ATPase inhibitor" # also labeled as TRPV agonist
#drug.corr.wInfo[drug.corr.wInfo$Drug == "fosbretabulin",]$moa <- "tubulin polymerization inhibitor" # also labeled as VE-cadherin antagonist
#drug.corr.wInfo[drug.corr.wInfo$Drug == "KW.2449",]$moa <- "Aurora kinase inhibitor" # also labeled as Abl kinase inhibitor and FTL3 inhibitor
drug.corr.wInfo[drug.corr.wInfo$Drug == "valrubicin",]$moa <- "topoisomerase inhibitor" # also labeled as DNA inhibitor
## vinblastine is labeled as both microtubule inhibitor and tubulin polymerization inhibitor
drug.corr.wInfo[drug.corr.wInfo$Drug == "dexrazoxane",]$moa <- "topoisomerase inhibitor" # also labeled as chelating agent

drug.corr.wInfo[drug.corr.wInfo$moa %in% sig.moas,]$Mechanism <- drug.corr.wInfo[drug.corr.wInfo$moa %in% sig.moas,]$moa
drug.corr.wInfo[drug.corr.wInfo$Drug == "vinblastine",]$Mechanism <- drug.corr.wInfo[drug.corr.wInfo$Drug == "vinblastine",]$moa

moaOrder <- unique(drug.corr.wInfo[order(drug.corr.wInfo$Pearson.est),]$Mechanism)
moaOrder3 <- c(moaOrder[moaOrder != "Other"], "Other")

myColorPal <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(12, "Set3"))(length(moaOrder3))
show_col(myColorPal)
myColorPal <- c(myColorPal[myColorPal != "#D9D9D9"], "#D9D9D9") # move grey to end for "Other" category

library(plyr);library(dplyr)

gsea.dot.plots <- NULL
for (i in omics) {
  nGenes <- length(unique(na.omit(drug.corr[drug.corr$type == i,])$Drug))
  
  temp.dot.df <- drug.corr.wInfo[drug.corr.wInfo$sig & drug.corr.wInfo$type == i,]
  topGenes <- temp.dot.df$name
  nTopGenes <- length(topGenes)
  temp.dot.df <- temp.dot.df %>% slice_max(abs(Pearson.est), n = 5)
  geneOrder <- temp.dot.df[order(temp.dot.df$Pearson.est),]$name
  
  plot.annot <- paste0(i, "\n(", nTopGenes, " / ", nGenes, " drugs correlated)")
  
  dot.plot <- ggplot2::ggplot(temp.dot.df,
                              ggplot2::aes(
                                x = name, y = Pearson.est, fill = Mechanism
                              )
  ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_discrete(limits = geneOrder) +
    scale_y_continuous(limits=c(-maxAbsEst, maxAbsEst)) +
    theme_classic() + scale_fill_manual(breaks=moaOrder3, values = myColorPal) +
    ggplot2::labs(
      #x = i,
      y = "Pearson r",
      fill = "Drug Mechanism"
    ) + theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))#,
    #legend.direction = "horizontal", legend.position = "bottom")
  
  cat(plot.annot, "\n")
  cat(i, "with", nTopGenes, "top Drugs\n")
  dot.plot <- dot.plot + ggtitle(plot.annot) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
  dot.plot
  #customWidth <- ifelse(nrow(temp.dot.df) > 14, 14, nrow(temp.dot.df))
  #ggplot2::ggsave(paste0(i, "_CorrelatedDrugs_barPlot_moaFill.pdf"), dot.plot, width=2, height=2.5)
  ggplot2::ggsave(paste0(i, "_CorrelatedDrugs_barPlot_moaFill.pdf"), dot.plot, width=3.5, height=2)
  
  if (is.null(gsea.dot.plots)) {
    gsea.dot.plots <- dot.plot
  } else {
    gsea.dot.plots <- gsea.dot.plots + dot.plot
  } 
}
#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
gsea.dot.plots2 <- gsea.dot.plots + plot_layout(guides = 'collect')
gsea.dot.plots2
ggplot2::ggsave("CorrelatedDrugs_barPlot_patchworkOmics_moaFill_sliceMaxAbsPearson5_oppositeMOAorder_height3.pdf", gsea.dot.plots, width=7, height=3) # was height 4

gsea.dot.plots <- NULL
for (i in rev(omics)) {
  nGenes <- length(unique(na.omit(drug.corr[drug.corr$type == i,])$Drug))
  
  temp.dot.df <- drug.corr.wInfo[drug.corr.wInfo$sig & drug.corr.wInfo$type == i,]
  topGenes <- temp.dot.df$name
  nTopGenes <- length(topGenes)
  temp.dot.df <- temp.dot.df %>% slice_max(abs(Pearson.est), n = 5)
  geneOrder <- temp.dot.df[order(-temp.dot.df$Pearson.est),]$name
  
  #plot.annot <- paste0(i, "\n(", nTopGenes, " / ", nGenes, " drugs correlated)")
  
  dot.plot <- ggplot2::ggplot(temp.dot.df,
                              ggplot2::aes(
                                x = name, y = Pearson.est, fill = Mechanism
                              )
  ) + scale_size(limits=c(0,maxLogFDR), range = c(0.5,4)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_discrete(limits = geneOrder) +
    scale_y_continuous(limits=c(-maxAbsEst, maxAbsEst)) +
    theme_classic() + scale_fill_manual(breaks=moaOrder3, values = myColorPal) +
    ggplot2::labs(
      #x = i,
      y = "Pearson r",
      fill = "Drug Mechanism"
    ) + theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) + coord_flip()#,
  #legend.direction = "horizontal", legend.position = "bottom")
  
  cat(plot.annot, "\n")
  cat(i, "with", nTopGenes, "top Drugs\n")
  dot.plot <- dot.plot + ggtitle(i) + theme(plot.title = element_text(hjust = 0.5, face="bold", size=16), axis.title.y=element_blank()) + xlab("Pearson r")
  dot.plot
  #customWidth <- ifelse(nrow(temp.dot.df) > 14, 14, nrow(temp.dot.df))
  #ggplot2::ggsave(paste0(i, "_CorrelatedDrugs_barPlot_moaFill.pdf"), dot.plot, width=2, height=2.5)
  ggplot2::ggsave(paste0(i, "_CorrelatedDrugs_barPlot_moaFill.pdf"), dot.plot, width=3.5, height=2)
  
  if (is.null(gsea.dot.plots)) {
    gsea.dot.plots <- dot.plot
  } else {
    gsea.dot.plots <- gsea.dot.plots + dot.plot
  } 
}
#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
gsea.dot.plots2 <- gsea.dot.plots + plot_layout(guides = 'collect')
gsea.dot.plots2
ggplot2::ggsave("CorrelatedDrugs_barPlot_patchworkOmics_moaFill_sliceMaxAbsPearson5_oppositeMOAorder_vertical_oppOrder.pdf", gsea.dot.plots, width=7, height=2) # was height 4
