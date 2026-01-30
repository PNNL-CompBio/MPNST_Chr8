# chr8 MPNST: figure 6
# Author: Belinda B. Garana, belinda.garana@pnnl.gov
#remove(list=ls())
library(plyr);library(dplyr);library(ggplot2);library(synapser)
library(patchwork);library(msigdbr)

#this wont' work!
#setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4")

source('panSEAFunctions.R')
source("compile_mCorr.R")


if(!exists('base.path'))
  base.path <- file.path(getwd(),Sys.Date())
 

##hasn't this analysis been run before? why are we running again?
getCCLEprot <- function(){
#get CCLE proteomics
    poss.file = file.path(base.path,'analysis','Proteomics','proteomics.csv.gz')
    if(!file.exists(poss.file)){
      download.file("https://api.figshare.com/v2/file/download/41466702", "proteomics.csv.gz")
      prot.df <- read.csv(gzfile("proteomics.csv.gz"))#,fileEncoding="UTF-16LE")
    }else{
      prot.df <- read.csv(poss.file)#,fileEncoding="UTF-16LE")
    }

    allgenes = readr::read_csv("https://api.figshare.com/v2/file/download/40576109")
    genes = allgenes |>
      dplyr::select(gene_symbol,entrez_id)|>
      dplyr::distinct()
    #genes <- genes[genes$gene_symbol %in% colnames(RNA.df)[2:ncol(RNA.df)], ]

    allsamples = readr::read_csv('https://api.figshare.com/v2/file/download/40576103')
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
 # }
  return(prot.df.noNA)
}

#updated with recalculated values
rna.corr <- read.csv(file.path(base.path,'analysis','RNA-seq','Differential_expression','Differential_expression_results.csv'))#synapser::synGet("syn72399335")$path)
prot.corr <- read.csv(file.path(base.path,'analysis','Proteomics','Differential_expression','Differential_expression_results.csv'))#synapser::synGet("syn72399210")$path)

inputs <- list("RNA" = na.omit(rna.corr[rna.corr$Spearman.q <= 0.05,]), # 655 genes or 50 with min N6 / 17717
               "Protein" = na.omit(prot.corr[prot.corr$Spearman.q <= 0.05,])) # 208 genes / 9013

adh.RNA <- get_CCLE_RNA()
adh.RNA2 <- adh.RNA[,c("CCLE_ID",colnames(adh.RNA)[colnames(adh.RNA) %in% rna.corr$Gene])] # 4448 genes in 327 cell lines or 3601 with min N 6

#why do we have two functions to do teh same thing?
temp.expr <- getCCLEprot()
#temp.expr <- get_CCLE_prot() ##how is this different from abve?
sample.info <- read.csv("https://raw.githubusercontent.com/BelindaBGarana/DMEA/refs/heads/shiny-app/Inputs/CCLE_sample_info.csv")
temp.expr.adherent <- temp.expr[temp.expr$CCLE_ID %in% sample.info[tolower(sample.info$culture_type)=="adherent",]$CCLE_Name,]
adh.prot2 <- temp.expr.adherent[,c("CCLE_ID",colnames(temp.expr.adherent)[colnames(temp.expr.adherent) %in% prot.corr$Gene])] # 5035 proteins in 215 cell lines
soft.sarc.info <- sample.info[sample.info$lineage == "soft_tissue" & grepl("sarcoma", sample.info$lineage_subtype, ignore.case=TRUE),] # 49

other.sample.info <- sample.info[!(sample.info$CCLE_Name %in% temp.expr.adherent$CCLE_ID) &
                                   sample.info$CCLE_Name %in% temp.expr$CCLE_ID,]
set.seed(2562)
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
#saveRDS(adh.DMEA, "DMEA.rds")

#source("https://raw.githubusercontent.com/PNNL-CompBio/MPNST_Chr8/refs/heads/main/figures/compile_mCorr.R")

setwd(file.path(base.path,'analysis'))
dir.create("DMEA_corr")
setwd("DMEA_corr")
save_to_synapse_v2(adh.DMEA.files, NULL)#"syn66295230")


compiled.drugCorr <- compile_mCorr(list("RNA" = adh.DMEA$all.results$RNA$corr.result,
                                        "Protein" = adh.DMEA$all.results$Protein$corr.result))
compiled.files <- list("DMEA_correlation_results.csv" = compiled.drugCorr$results,
                       "DMEA_correlation_results_compiled.csv" = compiled.drugCorr$mean.results,
                       "DMEA_correlation_venn_diagram.pdf" = compiled.drugCorr$venn.diagram,
                       "DMEA_correlation_dot_plot.pdf" = compiled.drugCorr$dot.plot,
                       "DMEA_correlation_correlation_matrix.csv" = compiled.drugCorr$corr,
                       "DMEA_correlation_correlation_matrix.pdf" = compiled.drugCorr$corr.matrix)

##SG getting rid of directories that dont exist
#setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_4")
#setwd("DMEA")
#this save to synapse was not used previously, only v2
save_to_synapse_v2(compiled.files, NULL)

setwd(base.path)
setwd('..')
#