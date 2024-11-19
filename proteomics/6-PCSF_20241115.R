# Differential expression & enrichment analyses: PDX transcriptomics, global & phospho
# Chr8: median copy number vs. protein/phospho expression
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-09-10
remove(list=ls())
# overview:
#### 1. Import metadata & crosstabs
#### 2. Run panSEA: only MPNST PDX with proteomics
#### 3. Run panSEA: all MPNST PDX
#### 4. PCSF network analysis

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr); library(plyr)
library(htmlwidgets); library(webshot); library(scales); library(msigdbr)
library(plyr); library(dplyr); library(R.utils)
library(ggplot2)
#webshot::install_phantomjs()
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
#source("panSEA_helper_20240508.R")
source("panSEA_helper_20240913.R")
source("https://raw.githubusercontent.com/PNNL-CompBio/amlresistancenetworks/master/R/proteinNetworks.R")

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/"
setwd(base.path)
setwd("Chr8_quant")
# gmt1 <- readRDS("gmt1_less.rds")
# 
# # run TFT GSEA for all RNAseq sample corr
# gene.allSamples.result <- read.csv(synapser::synGet("syn63394665")$path)
# gmt.pos <- gmt1$Positional
# gmt1TF <- gmt1$TFT_GTRD
# RNA.allSamples.TF <- panSEA::drugSEA_ties(gene.allSamples.result, gmt1TF, "Gene", "Spearman.est", ties = TRUE)
# write.csv(RNA.allSamples.TF$result, "Chr8_RNA_allSamples_GSEA_TFT_GTRD_results.csv", row.names = FALSE)
# saveRDS(RNA.allSamples.TF, "Chr8_RNA_allSamples_GSEA_TFT_GTRD_results.rds")

RNA.allSamples.TF <- list()
RNA.allSamples.TF$result <- read.csv("Chr8_RNA_allSamples_GSEA_TFT_GTRD_results.csv")

#### PCSF sig v3 ####
# try using top 50 sig prot, phospho, WES - directional this time
# load feature weights
synapser::synLogin()
global.result <- na.omit(read.csv(synapser::synGet("syn63435101")$path))
rna.allSamples.result <- RNA.allSamples.TF$result
rna.allSamples.result$Spearman.est <- rna.allSamples.result$NES
rna.allSamples.result$Gene <- sub("\\_.*","",rna.allSamples.result$Drug_set) # need feature names to be first column
mit.kin <- read.table("~/OneDrive - PNNL/Documents/MPNST/FAK/enrichment-analysis-result-table-belinda-correlations.txt", sep = "\t", header = TRUE)
mit.kin$Gene <- mit.kin$kinase
mit.kin$Spearman.est <- mit.kin$enrichment_value
mit.kin$Spearman.q <- mit.kin$adjusted_p_value

# get total N features
n.tf <- nrow(RNA.allSamples.TF$result) # 489
n.prot <- nrow(global.result) # 9013
n.kin <- nrow(mit.kin) # 381

# filter for significance
global.result <- global.result[global.result$Spearman.q <= 0.05, ] # 705 / 9013 (7.82%)
rna.allSamples.result <- rna.allSamples.result[rna.allSamples.result$FDR_q_value <= 0.25 & rna.allSamples.result$p_value <= 0.05, ] # 12 / 489 (2.45%) for FDR < 0.25; 4 / 489 (0.0818%) if FDR < 0.05
rna.allSamples.result <- rna.allSamples.result[,c("Gene","Spearman.est")]
mit.kin <- mit.kin[mit.kin$Spearman.q <= 0.25,c("Gene", "Spearman.est")] # 12/381 (3.15%) for q < 0.25; 2 / 381 (0.0525%) if q < 0.05
# when q < 0.05 for mit.kin, only P38G and ERK1 were enriched, neither were in STRING interactome
# 100% of your terminal nodes are included in the interactome
# 96.59% of your terminal nodes are included in the interactome
# 58.33% of your terminal nodes are included in the interactome
# Solving the PCSF by adding random noise to the edge costs...

# # replace kinases with gene symbols when possible - actually maybe MIT are using gene names?
# kinases <- unique(mit.kin$Gene)
# uniprot.kin <- data.frame()
# for (i in kinases) {
#   download.file(paste0("https://rest.uniprot.org/uniprotkb/stream?fields",
#                        "=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names",
#                        "%2Corganism_name%2Clength&format=xlsx&query=%28%28",
#                        "protein_name%3A",i,"%29+AND+%28organism_id%3A9606%29%29"),
#                 "uniprot_temp.xlsx")
#   temp.results <- readxl::read_excel("uniprot_temp.xlsx")
#   uniprot.kin <- rbind(uniprot.kin, temp.results)
# }
# write.csv(uniprot.kin, paste0("Kinase_uniprot_results_", Sys.Date(), ".csv"), row.names = FALSE)

inputs <- list("Transcription_factor" = rna.allSamples.result, "Protein" = global.result, "Kinase" = mit.kin)

directions <- c("positivelyCorr", "negativelyCorr")
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date())
dir.create(temp.path)
setwd(temp.path)
base.path2 <- getwd()
beta.vals <- c(2, 4, 6, 8, 10) # 2 was the default before
mu.vals <- c(1E-7, 1E-6, 1E-5, 1E-4, 1E-3) # 5E-4 was the default before
omega.vals <- c(1, 2, 3, 4) # 4 was default before
min.per.set <- c(3,4,5,6)
all.vertices <- data.frame()
all.edges <- data.frame()
all.enrichr <- data.frame()
all.gsea <- data.frame()
enr.8q24 <- data.frame()
timeout <- 1200
for (j in directions) {
  setwd(base.path2)
  dir.create(j)
  setwd(j)
  
  # prep named feature lists
  final.inputs <- list()
  for (m in 1:length(inputs)){ # order must be phospho, prot, gene (if any gene)
    if (grepl("positive", j)) {
      top.phospho <- inputs[[m]][inputs[[m]]$Spearman.est > 0,] 
    } else if (grepl("negative", j)) {
      top.phospho <- inputs[[m]][inputs[[m]]$Spearman.est < 0,]
    } else {
      top.phospho <- inputs[[m]]  
    } 
    
    # compile inputs in named lists
    final.inputs[[names(inputs)[m]]] <- as.list(top.phospho$Spearman.est)
    names(final.inputs[[names(inputs)[m]]]) <- top.phospho[,1] # assuming first column contains feature names
  }
  
  # parameter sweep: mu, beta.tf, beta.prot, beta.kin, omega
  for (mu in mu.vals) {
    setwd(file.path(base.path2, j))
    dir.create(paste0("mu_",mu))
    setwd(paste0("mu_",mu))
    
    for (omega in omega.vals) {
      setwd(file.path(base.path2, j, paste0("mu_",mu)))
      dir.create(paste0("omega_",omega))
      setwd(paste0("omega_",omega))
      
      for (beta.tf in beta.vals) {
        for (beta.prot in beta.vals) {
          for (beta.kin in beta.vals) {
            # run PCSF
            temp.b <- c(beta.tf, beta.prot, beta.kin)
            temp.fname <- paste0("beta_TF_", beta.tf, "_prot_", beta.prot, "_kin_", beta.kin)
            pcsf.result <- try(R.utils::withTimeout(computeProteinNetwork_BG(final.inputs, betas = temp.b,
                                                    mu = mu, w = omega, 
                                                    fname = temp.fname), timeout = timeout, onTimeout="error"))
            if (!inherits(pcsf.result, "try-error")) {
              my.plot <- PCSF::plot.PCSF(pcsf.result$graph)
              saveWidget(my.plot, paste0(temp.fname,".html"))
              webshot(paste0(temp.fname,".html"), paste0(temp.fname,".png"))
              webshot(paste0(temp.fname,".html"), paste0(temp.fname,".pdf"))
              webshot(paste0(temp.fname,".html"), paste0(temp.fname,".jpeg"))
              
              # enrichment analysis
              my.enr <- try(R.utils::withTimeout(PCSF::enrichment_analysis(pcsf.result$graph), timeout = timeout, onTimeout="error"), silent = TRUE)
              if (!inherits(my.enr, "try-error")) {
                write.csv(my.enr$enrichment, paste0(temp.fname,"_enrichmentAnalysis.csv"), row.names = FALSE)
                
                # compile enrichr
                temp.enrichr <- my.enr$enrichment
                temp.enrichr$Direction <- j
                temp.enrichr[,"mu"] <- mu
                temp.enrichr$betaTF <- beta.tf
                temp.enrichr$betaProt <- beta.prot
                temp.enrichr$betaKin <- beta.kin
                temp.enrichr[,"omega"] <- omega
                all.enrichr <- rbind(all.enrichr, temp.enrichr) 
              }
              
              # store vertex, edge data
              pcsf.tab <- igraph::as_data_frame(pcsf.result$graph, what="both")
              write.csv(pcsf.tab$vertices, paste0(temp.fname, "_vertices.csv"), row.names = FALSE)
              write.csv(pcsf.tab$edges, paste0(temp.fname, "_edges.csv"), row.names = FALSE)
              
              # compile vertex data
              temp.vert <- pcsf.tab$vertices
              temp.vert$Direction <- j
              temp.vert[,"mu"] <- mu
              temp.vert$betaTF <- beta.tf
              temp.vert$betaProt <- beta.prot
              temp.vert$betaKin <- beta.kin
              temp.vert[,"omega"] <- omega
              all.vertices <- rbind(all.vertices, temp.vert)
              
              # compile edge data
              temp.edges <- pcsf.tab$edgesices
              temp.edges$Direction <- j
              temp.edges[,"mu"] <- mu
              temp.edges$betaTF <- beta.tf
              temp.edges$betaProt <- beta.prot
              temp.edges$betaKin <- beta.kin
              temp.edges[,"omega"] <- omega
              all.edges <- rbind(all.edges, temp.edges)
              
              # get centrality; inspired by: https://stackoverflow.com/questions/62739592/how-to-extract-node-id-from-a-igraph-object-into-a-data-frame-object-in-r
              centrality.df <- data.frame(name = V(pcsf.result$graph)$NodeID,
                                          degree = pcsf.result$graph$degree,
                                          closeness = pcsf.result$closeness,
                                          betweenness = pcsf.result$betweenness,
                                          eigenvector = pcsf.result$eig,
                                          hub = pcsf.result$hubs,
                                          authority = pcsf.result$authorities)
              
              # positional enrichment analysis
              for (q in min.per.set) {
                pos.enr <- try(R.utils::withTimeout(panSEA::drugSEA_ties(pcsf.tab$vertices, gmt = gmt.pos,
                                                                         drug="name", rank.metric="logFoldChange",
                                                                         stat.type="Classic", ties = FALSE,
                                                                         min.per.set = q), timeout = timeout, onTimeout="error"), silent = TRUE)
                if (!inherits(pos.enr, "try-error")) {
                  colnames(pos.enr$result)[2] <- "Feature_set"
                  
                  # compile GSEA results
                  temp.gsea <- as.data.frame(pos.enr$result)
                  temp.gsea$Direction <- j
                  temp.gsea[,"mu"] <- mu
                  temp.gsea$betaTF <- beta.tf
                  temp.gsea$betaProt <- beta.prot
                  temp.gsea$betaKin <- beta.kin
                  temp.gsea[,"omega"] <- omega
                  temp.gsea$minPerSet <- q
                  all.gsea <- rbind(all.gsea, temp.gsea)
                  
                  # make plots
                  GSEA.global.mtn <- get_top_mtn_plots(pos.enr,
                                                       sets = "Drug_set",
                                                       EA.type = "GSEA")
                  GSEA.global.net <- try(R.utils::withTimeout(panSEA::netSEA(list(pcsf.tab$vertices),
                                                                             list(pos.enr$result),
                                                                             "name", "prize",
                                                                             n.network.sets = 5), 
                                                              timeout = timeout, onTimeout="error"), silent = TRUE)
                  
                  # store files
                  global.GSEA.files <- list("GSEA_results.csv" = pos.enr$result,
                                            "GSEA_volcano_plot.pdf" = pos.enr$volcano.plot)
                  if (!inherits(GSEA.global.net, "try-error")) {
                    global.GSEA.files[["GSEA_network_graph.html"]] <- GSEA.global.net$interactive
                  }
                  if (length(GSEA.global.mtn) > 0) {
                    global.GSEA.files[["mtn_plots"]] <- GSEA.global.mtn
                  }
                  GSEA.files <- list("temp" = global.GSEA.files)
                  names(GSEA.files) <- paste0("GSEA_", temp.fname, "_minPerSet_", q)
                  save_to_synapse_v2(GSEA.files, dot.scale = 1)
                  
                  # extract chr8q24 results
                  if ("chr8q24" %in% pos.enr$result$Feature_set) {
                    temp.8q24 <- as.data.frame(pos.enr$result[pos.enr$result$Feature_set == "chr8q24",])
                    temp.8q24$Direction <- j
                    temp.8q24[,"mu"] <- mu
                    temp.8q24$betaTF <- beta.tf
                    temp.8q24$betaProt <- beta.prot
                    temp.8q24$betaKin <- beta.kin
                    temp.8q24[,"omega"] <- omega
                    temp.8q24$minPerSet <- q
                    enr.8q24 <- rbind(enr.8q24, temp.8q24)
                  }
                }
              } 
            }
            # omega 4, mu 1E-3, b = 2:
            # Error in other.vals[[j]] : subscript out of bounds
          }
        }
      }
    }
  }
}
setwd(base.path2)
write.csv(all.vertices,"vertices.csv", row.names = FALSE)
write.csv(all.edges,"edges.csv", row.names = FALSE)
write.csv(all.enrichr,"enrichr.csv", row.names = FALSE)
write.csv(all.gsea,"GSEA_positional.csv", row.names = FALSE)
write.csv(enr.8q24,"chr8q24_enrichment.csv", row.names = FALSE)

#### plot results of parameter sweep ####
### node degrees (steiner vs. terminal) for each parameter setting

### heatmap of nodes in vs. out clustered with parameter setting information

### % of terminal nodes included in output for each parameter setting

### combo plot (no clustering so that parameter settings are in same order)