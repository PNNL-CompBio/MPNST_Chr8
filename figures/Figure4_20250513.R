# chr8 MPNST: figure 6
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4")

synapser::synLogin()

### redo just using significantly corr features
source("~/OneDrive - PNNL/Documents/GitHub/Chr8/panSEA_helper_20240913.R")

# load data
rna.corr <- read.csv(synapser::synGet("syn66226866")$path)
prot.corr <- read.csv(synapser::synGet("syn66224803")$path)
inputs <- list("RNA" = na.omit(rna.corr[rna.corr$Spearman.q <= 0.05,]), # 655 genes
               "Protein" = na.omit(prot.corr[prot.corr$Spearman.q <= 0.05,])) # 208 genes

adh.RNA <- get_CCLE_RNA()
adh.RNA2 <- adh.RNA[,c("CCLE_ID",colnames(adh.RNA)[colnames(adh.RNA) %in% rna.corr$Gene])] # 4448 genes in 327 cell lines

temp.expr <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/proteomics/analysis/CCLE_proteomics.csv")
sample.info <- read.csv("/Users/gara093/Documents/misc/DMEA-shiny-app/Inputs/CCLE_sample_info.csv")
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

source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/compile_mCorr.R")
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
drug.corr <- read.csv(synapser::synGet("syn66295273")$path)
mean.drugs <- read.csv(synapser::synGet("syn66295274")$path)
moa.results <- read.csv(synapser::synGet("syn66295241")$path)
mean.moa <- read.csv(synapser::synGet("syn66295242")$path)
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
ggplot2::ggsave("CorrelatedDrugs_barPlot_patchworkOmics_moaFill_sliceMaxAbsPearson5_oppositeMOAorder.pdf", gsea.dot.plots, width=7, height=4)

#### compile target ranks across all analyses - manual network ####
# get correlated genes (including from phospho)
corr.result <- read.csv(synapser::synGet("syn66227486")$path)
corr.result <- corr.result[corr.result$Omics != "DNA" & corr.result$Significant,]
corr.result$Gene <- corr.result$Feature
corr.result[corr.result$Omics=="Phospho",]$Gene <- sub("-.*","",corr.result[corr.result$Omics=="Phospho",]$Feature)
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
kin <- read.csv(synapser::synGet("syn66279699")$path)
kin$Gene <- kin$Feature_set
kin <- kin[kin$FDR_q_value <= 0.25 & kin$p_value <= 0.05,] # 2
gmt2 <- readRDS("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp21_NRAS_ASO_treated_patients/proteomics/analysis/gmt2.rds")[[1]]
site.mapping <- data.frame(kin$Feature_set, Substrates = "")
colnames(site.mapping)[1] <- "Feature_set"
for (i in kin$Feature_set) {
  tempSites <- unlist(gmt2$genesets[[i]])
  Gene <- sub("-.*","",tempSites)
  site.mapping[site.mapping$Feature_set == i,]$Substrates <- paste0(unique(Gene), collapse = " ")
}
kin2 <- merge(kin, site.mapping, by="Feature_set", all.x = TRUE)
kin2$Rank <- rank(kin2$NES)
bestRank <- max(kin2$Rank)
bestRankKin <- max(kin2$Rank)
kin2$BestRank <- FALSE
kin2[kin2$Rank == bestRank,]$BestRank <- TRUE
kin2$percentile <- percent_rank(kin2$NES)

tf <- read.csv(synapser::synGet("syn66226952")$path)
tf <- tf[tf$p_value <= 0.05 & tf$FDR_q_value <= 0.25,] # 206
tf$Gene <- sub("_TARGET_GENES.*","",tf$Feature_set)
tf.genes <- msigdbr::msigdbr(species="Homo sapiens",category="C3", subcategory="GTRD")
tf.genes <- tf.genes[tf.genes$gs_name %in% tf$Feature_set,] # 1871 rows
tf.genes <- dplyr::distinct(tf.genes[tf.genes$gs_name %in% tf$Feature_set,c("gs_name","gene_symbol")]) # 173217 rows
colnames(tf.genes) <- c("Feature_set","Target")
tf.genes <- plyr::ddply(tf.genes, .(Feature_set), summarize,
                        Targets = paste0(unique(Target), collapse = " "))
tf$Rank <- rank(tf$NES)
bestRank <- max(tf$Rank)
bestRankTF <- max(tf$Rank)
tf$BestRank <- FALSE
tf[tf$Rank == bestRank,]$BestRank <- TRUE
tf2 <- merge(tf, tf.genes, by="Feature_set", all.x = TRUE)
tf2$percentile <- percent_rank(tf2$NES)

# GSEA
coi <- c("Hallmark")
all.gsea <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_2/SupplementaryTable2_GSEA_and_KSEA.csv")
all.gsea$Significant <- FALSE
all.gsea[all.gsea$p_value <= 0.05 & all.gsea$FDR_q_value <= 0.25,]$Significant <- TRUE
all.gsea <- all.gsea[all.gsea$Omics %in% c("RNA", "Protein") &
                       all.gsea$Collection %in% coi & all.gsea$Significant,] # 181
collections <- unique(all.gsea$Collection) # same as coi
hallmark <- msigdbr::msigdbr(species="Homo sapiens",category="H")
hallmark <- dplyr::distinct(hallmark[hallmark$gs_name %in% all.gsea$Feature_set,c("gs_name","gene_symbol")]) # 2358 rows
colnames(hallmark) <- c("Feature_set","Gene")
hallmark <- plyr::ddply(hallmark, .(Feature_set), summarize,
                        Gene = paste0(unique(Gene), collapse = " "))
all.gsea2 <- merge(hallmark, all.gsea, by="Feature_set", all.y = TRUE)

all.gsea2$Rank <- rank(-all.gsea2$NES)
bestRank <- max(all.gsea2$Rank)
bestRankGSEA <- max(all.gsea2$Rank)
all.gsea2$BestRank <- FALSE
all.gsea2[all.gsea2$Rank == bestRank,]$BestRank <- TRUE
all.gsea2$percentile <- percent_rank(all.gsea2$NES)

# get correlated drugs and their targets
drug.corr <- read.csv(synapser::synGet("syn66295273")$path)
drug.corr$Significant <- FALSE
drug.corr[drug.corr$Pearson.q <= 0.05,]$Significant <- TRUE
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
moa.results <- read.csv(synapser::synGet("syn66295241")$path)
moa.results$Significant <- FALSE
moa.results[moa.results$p_value <= 0.05 & moa.results$FDR_q_value <= 0.25,]$Significant <- TRUE
moa.results <- moa.results[moa.results$type %in% c("RNA", "Protein") & moa.results$Significant,] # 19
colnames(moa.results)[3] <- "moa"
moa.info <- dplyr::distinct(drug.info[drug.info$moa %in% moa.results$moa,c("moa","target")]) # 19 MOAs
moa.info <- plyr::ddply(moa.info, .(moa), summarize,
                        Gene = paste0(unique(target), collapse = ", ")) # 10 moas
moa.results2 <- merge(moa.results, moa.info, by="moa", all.x = TRUE)
moa.results2$Rank <- rank(-moa.results2$NES)
bestRank <- max(moa.results2$Rank)
bestRankMOA <- max(moa.results2$Rank)
moa.results2$BestRank <- FALSE
moa.results2[moa.results2$Rank == bestRank,]$BestRank <- TRUE
moa.results2$percentile <- percent_rank(-moa.results2$NES)

# get selected network nodes and their centrality
pos.net.centr <- read.csv(file.path("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics",
                                "positive_centrality.csv"))
pos.net.centr$Rank <- rank(pos.net.centr$eigen_centrality)
bestRank <- max(pos.net.centr$Rank)
bestRankPosNet <- max(pos.net.centr$Rank)
pos.net.centr$BestRank <- FALSE
pos.net.centr[pos.net.centr$Rank == bestRank,]$BestRank <- TRUE
pos.net.centr$percentile <- percent_rank(pos.net.centr$eigen_centrality)

neg.net.centr <- read.csv(file.path("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics",
                                    "negative_centrality.csv"))
neg.net.centr$Rank <- rank(neg.net.centr$eigen_centrality)
bestRank <- max(neg.net.centr$Rank)
bestRankNegNet <- max(neg.net.centr$Rank)
neg.net.centr$BestRank <- FALSE
neg.net.centr[neg.net.centr$Rank == bestRank,]$BestRank <- TRUE
neg.net.centr$percentile <- percent_rank(neg.net.centr$eigen_centrality)

names(myColorPal) <- MOAsInTop50

all.features <- unique(c(corr.result$Gene, kin2$Gene, tf2$Gene,
                         unlist(strsplit(tf2$Targets, " ")),
                         unlist(strsplit(kin2$Substrates, " ")),
                         unlist(strsplit(all.gsea2$Gene, " ")), 
                         unlist(strsplit(moa.results2$Gene, " ")), 
                         drug.targets, pos.net.centr$name, neg.net.centr$name)) # 23853
other.features <- unique(c(corr.result$Gene, kin2$Gene, tf$Gene,
                           unlist(strsplit(tf2$Targets, " ")),
                           unlist(strsplit(kin2$Substrates, " ")),
                           unlist(strsplit(all.gsea2$Gene, " ")), 
                           pos.net.centr$name, neg.net.centr$name)) # 23816

Target <- drug.targets[drug.targets %in% other.features] # 87; identify drug targets also implicated in at least 1 other analysis
excluded.targets <- drug.targets[!(drug.targets %in% other.features)] # 8: CA7, TUBA3C, TUBA3D, TUBA3E, TUBB1, KDR, HTR2A, HTR2B

toxic <- c("chr8qAmp","nonAmp")
#amp.nonAmp.plot <- NULL
for (j in toxic) {
  if (j == "chr8qAmp") {
    toxCorr <- corr.result[corr.result$Spearman.est > 0,]
    toxGSEA <- all.gsea2[all.gsea2$NES > 0,]
    toxTF <- tf2[tf2$NES > 0,]
    toxKin <- kin2[kin2$NES > 0,]
    toxDrug <- drug.corr.wInfo[drug.corr.wInfo$Pearson.est < 0,]
    toxMOA <- moa.results2[moa.results2$NES < 0,]
    toxNet <- pos.net.centr
    bestRankNet <- bestRankPosNet
  } else if (j == "nonAmp") {
    toxCorr <- corr.result[corr.result$Spearman.est < 0,]
    toxGSEA <- all.gsea2[all.gsea2$NES < 0,]
    toxTF <- tf2[tf2$NES < 0,]
    toxKin <- kin2[kin2$NES < 0,]
    toxDrug <- drug.corr.wInfo[drug.corr.wInfo$Pearson.est > 0,]
    toxMOA <- moa.results2[moa.results2$NES > 0,]
    toxNet <- neg.net.centr
    bestRankNet <- bestRankNegNet
  }
  percScoresInfo <- data.frame(Target)
  percScoresInfo[,c("Correlation", "GSEA", "TF", "TFT", "Kinase", "Substrate", "Network", 
                    "Drug", "DMEA", "MeanPerc", "N","N_best")] <- NA
  percScoresInfo$Best <- ""
  percScores <- percScoresInfo
  nScores <- percScoresInfo
  for (i in 1:length(Target)) {
    nRankVals <- 0
    nBestRanks <- 0
    
    # correlations
    if (Target[i] %in% unique(toxCorr$Gene)) {
      tempCorr <- toxCorr[toxCorr$Gene == Target[i],]
      tempTypes <- paste0(unique(tempCorr$Omics), collapse = " ")
      tempScores <- round(mean(tempCorr$Spearman.est), digits=2)
      tempPerc <- round(mean(tempCorr$percentile), digits=2)
      
      # keep score
      nRankVals <- nRankVals + nrow(tempCorr)
      if (any(tempCorr$Rank == bestRankCorr)) {
        nBestRanks <- nBestRanks + length(which(tempCorr$Rank == bestRankCorr))
        percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "Correlation")
      }
      
      # update data frame
      percScoresInfo$Correlation[i] <- paste0(tempTypes, " | ", tempScores, " | ", tempPerc)
      percScores$Correlation[i] <- mean(tempCorr$percentile)
      nScores$Correlation[i] <- nrow(tempCorr)
    }
    
    # GSEA
    if (Target[i] %in% unique(unlist(strsplit(toxGSEA$Gene, " ")))) {
      tempGSEA <- toxGSEA[grepl(Target[i], toxGSEA$Gene),]
      tempTypes <- paste0(unique(tempGSEA$Omics), collapse = " ")
      tempSets <- paste0(tempGSEA$Feature_set, collapse = " ")
      tempScores <- round(mean(tempGSEA$NES), digits=2)
      tempPerc <- round(mean(tempGSEA$percentile), digits=2)
      
      # keep score
      nRankVals <- nRankVals + nrow(tempGSEA)
      if (any(tempGSEA$Rank == bestRankGSEA)) {
        nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankGSEA))
        percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "GSEA")
      }
      
      # update data frame
      percScoresInfo$GSEA[i] <- paste0(tempTypes, " | ", tempSets, " | ", tempScores, " | ", tempPerc)
      percScores$GSEA[i] <- mean(tempGSEA$percentile)
      nScores$GSEA[i] <- nrow(tempGSEA)
    }
    
    if (nrow(toxTF) > 0) {
      # TF
      if (Target[i] %in% toxTF$Gene) {
        tempGSEA <- toxTF[toxTF$Gene == Target[i],]
        tempScores <- round(mean(tempGSEA$NES), digits=2)
        tempPerc <- round(mean(tempGSEA$percentile), digits=2)
        
        # keep score
        nRankVals <- nRankVals + nrow(tempGSEA)
        if (any(tempGSEA$Rank == bestRankTF)) {
          nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankTF))
          percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "TF")
        }
        
        # update data frame
        percScoresInfo$TF[i] <- paste0(tempScores, " | ", tempPerc)
        percScores$TF[i] <- mean(tempGSEA$percentile)
        nScores$TF[i] <- nrow(tempGSEA)
      }
      
      # TFT
      if (Target[i] %in% unique(unlist(strsplit(toxTF$Targets, " ")))) {
        tempGSEA <- toxTF[grepl(Target[i], toxTF$Targets),]
        tempSets <- paste0(tempGSEA$Gene, collapse = " ")
        tempScores <- round(mean(tempGSEA$NES), digits=2)
        tempPerc <- round(mean(tempGSEA$percentile), digits=2)
        
        # keep score
        nRankVals <- nRankVals + nrow(tempGSEA)
        if (any(tempGSEA$Rank == bestRankTF)) {
          nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankTF))
          percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "TFT")
        }
        
        # update data frame
        percScoresInfo$TFT[i] <- paste0(tempSets, " | ", tempScores, " | ", tempPerc)
        percScores$TFT[i] <- mean(tempGSEA$percentile)
        nScores$TFT[i] <- nrow(tempGSEA)
      }
    }
    
    if (nrow(toxKin) > 0) {
      # Kinase
      if (Target[i] %in% toxKin$Gene) {
        tempGSEA <- toxKin[toxKin$Gene == Target[i],]
        tempScores <- round(mean(tempGSEA$NES),2)
        tempPerc <- round(mean(tempGSEA$percentile),2)
        
        # keep score
        nRankVals <- nRankVals + nrow(tempGSEA)
        if (any(tempGSEA$Rank == bestRankKin)) {
          nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankKin))
          percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "Kinase")
        }
        
        # update data frame
        percScoresInfo$Kinase[i] <- paste0(tempScores, " | ", tempPerc)
        percScores$Kinase[i] <- mean(tempGSEA$percentile)
        nScores$Kinase[i] <- nrow(tempGSEA)
      }
      
      # Substrate
      if (Target[i] %in% unique(unlist(strsplit(toxKin$Substrates, " ")))) {
        tempGSEA <- toxKin[grepl(Target[i], toxKin$Substrates),]
        tempSets <- paste0(tempGSEA$Gene, collapse = " ")
        tempScores <- round(mean(tempGSEA$NES),2)
        tempPerc <- round(mean(tempGSEA$percentile),2)
        
        # keep score
        nRankVals <- nRankVals + nrow(tempGSEA)
        if (any(tempGSEA$Rank == bestRankKin)) {
          nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankKin))
          percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "Substrate")
        }
        
        # update data frame
        percScoresInfo$Substrate[i] <- paste0(tempSets, " | ", tempScores, " | ", tempPerc)
        percScores$Substrate[i] <- mean(tempGSEA$percentile)
        nScores$Substrate[i] <- nrow(tempGSEA)
      } 
    }
    
    # Network
    if (Target[i] %in% toxNet$name) {
      tempGSEA <- toxNet[toxNet$name == Target[i],]
      tempScores <- round(mean(tempGSEA$eigen_centrality),2)
      tempPerc <- round(mean(tempGSEA$percentile),2)
      
      # keep score
      nRankVals <- nRankVals + nrow(tempGSEA)
      if (any(tempGSEA$Rank == bestRankNet)) {
        nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankNet))
        percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "Network")
      }
      
      # update data frame
      percScoresInfo$Network[i] <- paste0(tempScores, " | ", tempPerc)
      percScores$Network[i] <- mean(tempGSEA$percentile)
      nScores$Network[i] <- nrow(tempGSEA)
    }
    
    # Drug correlation
    allToxDrugTargets <- unique(unlist(strsplit(toxDrug$target, ", ")))
    if (Target[i] %in% allToxDrugTargets) {
      tempGSEA <- toxDrug[grepl(Target[i], toxDrug$target),]
      tempDrugs <- paste0(tempGSEA$name, collapse = " ")
      tempScores <- round(mean(tempGSEA$Pearson.est),2)
      tempPerc <- round(mean(tempGSEA$percentile),2)
      
      # keep score
      nRankVals <- nRankVals + nrow(tempGSEA)
      if (any(tempGSEA$Rank == bestRankDrug)) {
        nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankDrug))
        percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "Drug")
      }
      
      # update data frame
      percScoresInfo$Drug[i] <- paste0(tempDrugs, " | ", tempScores, " | ", tempPerc)
      percScores$Drug[i] <- mean(tempGSEA$percentile)
      nScores$Drug[i] <- nrow(tempGSEA)
    }
    
    # DMEA
    if (Target[i] %in% unique(unlist(strsplit(toxMOA$Gene, ", ")))) {
      tempGSEA <- toxMOA[grepl(Target[i], toxMOA$Gene),]
      tempTypes <- paste0(unique(tempGSEA$type), collapse = " ")
      tempSets <- paste0(tempGSEA$moa, collapse = " ")
      tempScores <- round(mean(tempGSEA$NES),2)
      tempPerc <- round(mean(tempGSEA$percentile),2)
      
      # keep score
      nRankVals <- nRankVals + nrow(tempGSEA)
      if (any(tempGSEA$Rank == bestRankMOA)) {
        nBestRanks <- nBestRanks + length(which(tempGSEA$Rank == bestRankMOA))
        percScoresInfo$Best[i] <- paste(percScoresInfo$Best[i], "DMEA")
      }
      
      # update data frame
      percScoresInfo$DMEA[i] <- paste0(tempTypes, " | ", tempSets, " | ", tempScores, " | ", tempPerc)
      percScores$DMEA[i] <- mean(tempGSEA$percentile)
      nScores$DMEA[i] <- nrow(tempGSEA)
    }
    
    # update data frame
    percScores$N[i] <- nRankVals
    percScores$N_best[i] <- nBestRanks
    percScores$MeanPerc[i] <- mean(as.numeric(percScores[i,2:10]), na.rm = TRUE)
  }
  percScoresInfo$N_best <- percScores$N_best
  percScoresInfo$N <- percScores$N
  percScoresInfo$MeanPerc <- percScores$MeanPerc
  percScoresInfo$Best <- substring(percScoresInfo$Best, 2)
  nScores$N_best <- percScoresInfo$N_best
  nScores$N <- percScoresInfo$N
  nScores$MeanPerc <- percScoresInfo$MeanPerc
  nScores$Best <- percScoresInfo$Best
  percScores$Best <- percScoresInfo$Best
  percScores$MeanPercTimesN <- percScores$MeanPerc*percScores$N
  percScoresInfo$MeanPercTimesN <- percScores$MeanPercTimesN
  
  percScores$N_analyses <- rowSums(!is.na(percScores[,2:10]))
  percScores$MeanPercTimesN_analyses <- percScores$N_analyses * percScores$MeanPerc
  percScoresInfo$MeanPercTimesN_analyses <- percScores$MeanPercTimesN_analyses
  write.csv(nScores, paste0(j,"_N_Percentile_scores_",Sys.Date(),".csv"), row.names = FALSE)
  #nScores <- read.csv(paste0(j,"_N_Percentile_scores_",Sys.Date(),".csv"))
  write.csv(percScores, paste0(j,"_Percentile_scores_",Sys.Date(),".csv"), row.names = FALSE)
  write.csv(percScoresInfo, paste0(j, "_Percentile_scores_info_",Sys.Date(),".csv"), row.names = FALSE)
  
  filtered.percScores <- percScores[!is.na(percScores$Network) & !is.na(percScores$DMEA),]
  if (nrow(filtered.percScores) > 0) {
    filtered.percScores <- filtered.percScores[order(-filtered.percScores$Drug),]
    targetOrder <- filtered.percScores$Target
    filtered.percScoresInfo <- percScoresInfo[percScoresInfo$Target %in% targetOrder,]
    filtered.percScoresInfo$moa <- stringr::str_split_i(filtered.percScoresInfo$DMEA, "[|]", 2)
    filtered.percScoresInfo$moa <- substr(filtered.percScoresInfo$moa, 2, nchar(filtered.percScoresInfo$moa)-1)
    
    write.csv(filtered.percScores, paste0(j,"_Percentile_scores_DMEA-Network_",Sys.Date(),".csv"), row.names = FALSE)
    write.csv(filtered.percScoresInfo, paste0(j,"_Percentile_scores_info_DMEA-Network_",Sys.Date(),".csv"), row.names = FALSE)
    
    # create tile plot for target MOAs
    
    # moa.plot <- ggplot2::ggplot(filtered.percScoresInfo, aes(x=Target,y=1, fill=moa)) + 
    #   geom_tile() + scale_x_discrete(limits=targetOrder) + theme_classic(base_size=12) +
    #   scale_fill_manual(values=myColorPal, breaks=moaOrder3) +
    #   theme(axis.ticks=element_blank(), axis.text = element_blank(), axis.title=element_blank(),
    #         axis.line=element_blank()) + labs(fill="Drug Mechanism")
    # ggplot2::ggsave(paste0(j,"_MOA_annotations_",Sys.Date(),".pdf"), moa.plot, width=6,height=3)
    
    filtered.percScores$N <- NULL
    filtered.percScores$N_best <- NULL
    filtered.percScores$N_analyses <- NULL
    filtered.percScores$MeanPerc <- NULL
    filtered.percScores$MeanPercTimesN <- NULL
    filtered.percScores$MeanPercTimesN_analyses <- NULL
    nScores2 <- reshape2::melt(nScores[nScores$Target %in% filtered.percScores$Target,
                                       colnames(filtered.percScores)], 
                               id = c("Target","Best"), variable.name = "Analysis", value.name="N")
    filtered.percScores <- reshape2::melt(filtered.percScores, id = c("Target", "Best"), 
                                          variable.name = "Analysis", value.name = "Mean Percentile")
    
    plot.df <- merge(filtered.percScores, nScores2, by=c("Target","Best","Analysis"))
    plot.df$`Mean Percentile` <- as.numeric(plot.df$`Mean Percentile`)
    plot.df$N <- as.numeric(plot.df$N)
    plot.df <- na.omit(plot.df)
    plot.df$`Top Score` <- FALSE
    if (any(plot.df$Best == plot.df$Analysis)) {
      plot.df[plot.df$Best == plot.df$Analysis,]$`Top Score` <- TRUE 
    }
    # create dot plot to represent filtered.percScores
    dot.plot <- ggplot2::ggplot(plot.df, aes(y = Target, x = Analysis, color = 100*`Mean Percentile`, size=N)) + 
      ggplot2::geom_point() + scale_color_gradient2(low="grey",high="red", limits=c(0, 100)) +
      theme_classic() + ggplot2::labs(color = "Mean Percentile", size = "N") + 
      ggplot2::scale_y_discrete(limits = targetOrder) +
      geom_point(data = subset(plot.df, `Top Score`), col = "black", stroke = 1.5, shape = 21) +
      theme(axis.text.y = element_text(angle = 45, vjust=1, hjust=1),
            axis.text.x = element_text(angle=90,vjust=0.5),
            axis.title.x=element_text(angle=180)
      ) + theme(#legend.position="top", 
        legend.direction="horizontal") + 
      coord_flip()
    dot.plot
    saveRDS(dot.plot, paste0(j,"_Target_DMEA-Network_dotPlot_",Sys.Date(),".rds"))
    ggsave(paste0(j,"_Target_DMEA-Network_dotPlot_wider_",Sys.Date(),".pdf"),dot.plot,width=10, height=3)
  }
}
