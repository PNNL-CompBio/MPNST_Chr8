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
source("~/Downloads/proteinNetworks_BG.R")
source("~/Downloads/PCSF_rand_BG.R")
#source("https://raw.githubusercontent.com/PNNL-CompBio/amlresistancenetworks/master/R/proteinNetworks.R")

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/"
setwd(base.path)
setwd("Chr8_quant")
gmt1 <- readRDS("gmt1_less.rds")
# 
# # run TFT GSEA for all RNAseq sample corr
# gene.allSamples.result <- read.csv(synapser::synGet("syn63394665")$path)
gmt.pos <- gmt1$Positional
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

directions <- c("positive", "negative")
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
#timeout <- 10*60 # unit of seconds
#temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date(), "_", timeout/60, "_min_timeout")
#temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date())
temp.path <- paste0("PCSF_TFProteinKinase_", "2024-11-22")
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
all.centrality <- data.frame()
runs.so.far <- c()
run.times <- c()

all.vertices <- read.csv("vertices.csv")
all.edges <- read.csv("edges.csv")
all.enrichr <- read.csv("enrichr.csv")
all.gsea <- read.csv("GSEA_positional.csv")
enr.8q24 <- read.csv("chr8q24_enrichment.csv")
all.centrality <- read.csv("centrality.csv")
runs.so.far <- unique(all.vertices$runID)
directions <- c("negative")

beta.vals <- c(2, 10) # next try adding 4
mu.vals <- c(1E-7, 1E-3) # next try adding 1E-5
omega.vals <- c(1, 4) # next try adding 2

# # Find out how many cores are available
# if (requireNamespace("parallel") &
#     requireNamespace("snow") &
#     requireNamespace("doSNOW")) {
#   cores <- parallel::detectCores()
#   if (cores[1] > 1) {
#     cl <- snow::makeCluster(cores[1] - 1) # cluster using all but 1 core
#     doSNOW::registerDoSNOW(cl) # register cluster
#   }
# }

for (j in directions) {
  setwd(base.path2)
  dir.create(j)
  setwd(j)
  
  # prep named feature lists
  cat("Preparing input\n")
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
    if (nrow(top.phospho) > 0) {
      final.inputs[[names(inputs)[m]]] <- as.list(top.phospho$Spearman.est)
      names(final.inputs[[names(inputs)[m]]]) <- top.phospho[,1] # assuming first column contains feature names
    }
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
      
      if (length(final.inputs) == 3) {
        for (beta.tf in beta.vals) {
          for (beta.prot in beta.vals) {
            for (beta.kin in beta.vals) {
              setwd(file.path(base.path2, j, paste0("mu_",mu), paste0("omega_",omega)))
              
              # run PCSF
              temp.b <- c(beta.tf, beta.prot, beta.kin)
              temp.fname <- paste0("beta_TF_", beta.tf, "_prot_", beta.prot, "_kin_", beta.kin)
              temp.runID <- paste0(j, "_omega_",omega,"_mu_",mu,"_",temp.fname)
              if (!(temp.runID %in% runs.so.far)) {
                cat("Running PCSF\n")
                time.start <- Sys.time()
                pcsf.result <- try(computeProteinNetwork_BG(final.inputs, 
                                                            betas = temp.b, 
                                                            mu = mu, w = omega, 
                                                            fname = temp.fname), 
                                   silent = TRUE)
                runs.so.far <- c(runs.so.far, temp.runID)
                if (!inherits(pcsf.result, "try-error")) {
                  run.times <- c(run.times, Sys.time() - time.start)
                  
                  cat("Storing graph and tabulated data\n")
                  my.plot <- PCSF::plot.PCSF(pcsf.result$graph)
                  saveWidget(my.plot, paste0(temp.fname,".html"))
                  #webshot(paste0(temp.fname,".html"), paste0(temp.fname,".png"))
                  #webshot(paste0(temp.fname,".html"), paste0(temp.fname,".pdf"))
                  #webshot(paste0(temp.fname,".html"), paste0(temp.fname,".jpeg"))
                  
                  # store vertex, edge data
                  pcsf.tab <- igraph::as_data_frame(pcsf.result$graph, what="both")
                  #write.csv(pcsf.tab$vertices, paste0(temp.fname, "_vertices.csv"), row.names = FALSE)
                  #write.csv(pcsf.tab$edges, paste0(temp.fname, "_edges.csv"), row.names = FALSE)
                  
                  # compile vertex data
                  temp.vert <- pcsf.tab$vertices
                  temp.vert$Direction <- j
                  temp.vert[,"mu"] <- mu
                  temp.vert$betaTF <- beta.tf
                  temp.vert$betaProt <- beta.prot
                  temp.vert$betaKin <- beta.kin
                  temp.vert[,"omega"] <- omega
                  temp.vert$runID <- temp.runID
                  all.vertices <- rbind(all.vertices, temp.vert)
                  
                  # compile edge data
                  temp.edges <- pcsf.tab$edges
                  temp.edges$Direction <- j
                  temp.edges[,"mu"] <- mu
                  temp.edges$betaTF <- beta.tf
                  temp.edges$betaProt <- beta.prot
                  temp.edges$betaKin <- beta.kin
                  temp.edges[,"omega"] <- omega
                  temp.edges$runID <- temp.runID
                  all.edges <- rbind(all.edges, temp.edges)
                  
                  # get centrality; inspired by: https://stackoverflow.com/questions/62739592/how-to-extract-node-id-from-a-igraph-object-into-a-data-frame-object-in-r
                  centrality.df <- data.frame(name = V(pcsf.result$graph)$name,
                                              degree = igraph::degree(pcsf.result$graph, mode="all"),
                                              closeness = igraph::closeness(pcsf.result$graph, mode="all"),
                                              betweenness = igraph::betweenness(pcsf.result$graph, directed = FALSE),
                                              eigen_centrality = igraph::eigen_centrality(pcsf.result$graph, directed = FALSE)$vector,
                                              hub_score = igraph::hub_score(pcsf.result$graph)$vector,
                                              authority_score = igraph::authority_score(pcsf.result$graph)$vector)
                  centrality.df$runID <- temp.runID
                  all.centrality <- rbind(all.centrality, centrality.df)
                  
                  # enrichment analysis
                  cat("Running enrichr\n")
                  my.enr <- try(PCSF::enrichment_analysis(pcsf.result$graph), 
                                silent = TRUE)
                  if (!inherits(my.enr, "try-error")) {
                    #write.csv(my.enr$enrichment, paste0(temp.fname,"_enrichmentAnalysis.csv"), row.names = FALSE)
                    
                    # compile enrichr
                    temp.enrichr <- my.enr$enrichment
                    temp.enrichr$Direction <- j
                    temp.enrichr[,"mu"] <- mu
                    temp.enrichr$betaTF <- beta.tf
                    temp.enrichr$betaProt <- beta.prot
                    temp.enrichr$betaKin <- beta.kin
                    temp.enrichr[,"omega"] <- omega
                    temp.enrichr$runID <- temp.runID
                    all.enrichr <- rbind(all.enrichr, temp.enrichr) 
                  }
                  
                  # positional enrichment analysis
                  cat("Running GSEA\n")
                  for (q in min.per.set) {
                    pos.enr <- try(panSEA::drugSEA_ties(pcsf.tab$vertices, 
                                                        gmt = gmt.pos, drug="name", 
                                                        rank.metric="prize", 
                                                        stat.type="Classic", ties = FALSE,
                                                        min.per.set = q), silent = TRUE)
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
                      temp.gsea$runID <- temp.runID
                      all.gsea <- rbind(all.gsea, temp.gsea)
                      
                      # make plots
                      GSEA.global.mtn <- get_top_mtn_plots(pos.enr,
                                                           sets = "Drug_set",
                                                           EA.type = "GSEA")
                      GSEA.global.net <- try(panSEA::netSEA(list(pcsf.tab$vertices),
                                                            list(pos.enr$result),
                                                            "name", "prize",
                                                            n.network.sets = 5),
                                             silent = TRUE)
                      
                      # store files
                      global.GSEA.files <- list("GSEA_volcano_plot.pdf" = pos.enr$volcano.plot)
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
                        temp.8q24$runID <- temp.runID
                        enr.8q24 <- rbind(enr.8q24, temp.8q24)
                      }
                    }
                  } 
                  
                  # store files after each run
                  setwd(base.path2)
                  write.csv(all.vertices,"vertices.csv", row.names = FALSE)
                  write.csv(all.edges,"edges.csv", row.names = FALSE)
                  write.csv(all.enrichr,"enrichr.csv", row.names = FALSE)
                  write.csv(all.gsea,"GSEA_positional.csv", row.names = FALSE)
                  write.csv(enr.8q24,"chr8q24_enrichment.csv", row.names = FALSE)
                  write.csv(all.centrality, "centrality.csv", row.names = FALSE)
                } 
              }
            }
          }
        } 
      } else { # negative network only uses global proteomics
        for (beta.prot in beta.vals) {
          setwd(file.path(base.path2, j, paste0("mu_",mu), paste0("omega_",omega)))
          
          # run PCSF
          temp.fname <- paste0("beta_prot_", beta.prot)
          temp.runID <- paste0(j, "_omega_",omega,"_mu_",mu,"_",temp.fname)
          if (!(temp.runID %in% runs.so.far)) {
            cat("Running PCSF\n")
            runs.so.far <- c(runs.so.far, temp.runID)
            time.start <- Sys.time()
            final.input <- list("Protein" = final.inputs[[1]])
            pcsf.result <- try(computeProteinNetwork_BG(final.input, 
                                                        betas = beta.prot, 
                                                        mu = mu, w = omega, 
                                                        fname = temp.fname), 
                               silent = TRUE)
            
            if (!inherits(pcsf.result, "try-error")) {
              run.times <- c(run.times, Sys.time() - time.start)
              
              cat("Storing graph and tabulated data\n")
              my.plot <- PCSF::plot.PCSF(pcsf.result$graph)
              saveWidget(my.plot, paste0(temp.fname,".html"))
              #webshot(paste0(temp.fname,".html"), paste0(temp.fname,".png"))
              #webshot(paste0(temp.fname,".html"), paste0(temp.fname,".pdf"))
              #webshot(paste0(temp.fname,".html"), paste0(temp.fname,".jpeg"))
              
              # store vertex, edge data
              pcsf.tab <- igraph::as_data_frame(pcsf.result$graph, what="both")
              #write.csv(pcsf.tab$vertices, paste0(temp.fname, "_vertices.csv"), row.names = FALSE)
              #write.csv(pcsf.tab$edges, paste0(temp.fname, "_edges.csv"), row.names = FALSE)
              
              # compile vertex data
              temp.vert <- pcsf.tab$vertices
              temp.vert$Direction <- j
              temp.vert[,"mu"] <- mu
              temp.vert$betaTF <- beta.tf
              temp.vert$betaProt <- beta.prot
              temp.vert$betaKin <- beta.kin
              temp.vert[,"omega"] <- omega
              temp.vert$runID <- temp.runID
              all.vertices <- rbind(all.vertices, temp.vert)
              
              # compile edge data
              temp.edges <- pcsf.tab$edges
              temp.edges$Direction <- j
              temp.edges[,"mu"] <- mu
              temp.edges$betaTF <- beta.tf
              temp.edges$betaProt <- beta.prot
              temp.edges$betaKin <- beta.kin
              temp.edges[,"omega"] <- omega
              temp.edges$runID <- temp.runID
              all.edges <- rbind(all.edges, temp.edges)
              
              # get centrality; inspired by: https://stackoverflow.com/questions/62739592/how-to-extract-node-id-from-a-igraph-object-into-a-data-frame-object-in-r
              centrality.df <- data.frame(name = V(pcsf.result$graph)$name,
                                          degree = igraph::degree(pcsf.result$graph, mode="all"),
                                          closeness = igraph::closeness(pcsf.result$graph, mode="all"),
                                          betweenness = igraph::betweenness(pcsf.result$graph, directed = FALSE),
                                          eigen_centrality = igraph::eigen_centrality(pcsf.result$graph, directed = FALSE)$vector,
                                          hub_score = igraph::hub_score(pcsf.result$graph)$vector,
                                          authority_score = igraph::authority_score(pcsf.result$graph)$vector)
              centrality.df$runID <- temp.runID
              all.centrality <- rbind(all.centrality, centrality.df)
              
              # enrichment analysis
              cat("Running enrichr\n")
              my.enr <- try(PCSF::enrichment_analysis(pcsf.result$graph), 
                            silent = TRUE)
              if (!inherits(my.enr, "try-error")) {
                #write.csv(my.enr$enrichment, paste0(temp.fname,"_enrichmentAnalysis.csv"), row.names = FALSE)
                
                # compile enrichr
                temp.enrichr <- my.enr$enrichment
                temp.enrichr$Direction <- j
                temp.enrichr[,"mu"] <- mu
                temp.enrichr$betaTF <- beta.tf
                temp.enrichr$betaProt <- beta.prot
                temp.enrichr$betaKin <- beta.kin
                temp.enrichr[,"omega"] <- omega
                temp.enrichr$runID <- temp.runID
                all.enrichr <- rbind(all.enrichr, temp.enrichr) 
              }
              
              # positional enrichment analysis
              cat("Running GSEA\n")
              for (q in min.per.set) {
                pos.enr <- try(panSEA::drugSEA_ties(pcsf.tab$vertices, 
                                                    gmt = gmt.pos, drug="name", 
                                                    rank.metric="prize", 
                                                    stat.type="Classic", ties = FALSE,
                                                    min.per.set = q), silent = TRUE)
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
                  temp.gsea$runID <- temp.runID
                  all.gsea <- rbind(all.gsea, temp.gsea)
                  
                  # make plots
                  GSEA.global.mtn <- get_top_mtn_plots(pos.enr,
                                                       sets = "Drug_set",
                                                       EA.type = "GSEA")
                  GSEA.global.net <- try(panSEA::netSEA(list(pcsf.tab$vertices),
                                                        list(pos.enr$result),
                                                        "name", "prize",
                                                        n.network.sets = 5),
                                         silent = TRUE)
                  
                  # store files
                  global.GSEA.files <- list("GSEA_volcano_plot.pdf" = pos.enr$volcano.plot)
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
                    temp.8q24$runID <- temp.runID
                    enr.8q24 <- rbind(enr.8q24, temp.8q24)
                  }
                }
              } 
              
              # store files after each run
              setwd(base.path2)
              write.csv(all.vertices,"vertices.csv", row.names = FALSE)
              write.csv(all.edges,"edges.csv", row.names = FALSE)
              write.csv(all.enrichr,"enrichr.csv", row.names = FALSE)
              write.csv(all.gsea,"GSEA_positional.csv", row.names = FALSE)
              write.csv(enr.8q24,"chr8q24_enrichment.csv", row.names = FALSE)
              write.csv(all.centrality, "centrality.csv", row.names = FALSE)
              saveRDS(runs.so.far, "runs_so_far.rds")
            } 
          }
        }
      }
    }
  }
}
# 
# if (requireNamespace("parallel") &
#     requireNamespace("snow") &
#     requireNamespace("doSNOW")) {
#   if (cores[1] > 1) {
#     snow::stopCluster(cl) # stop cluster
#     rm(cl)
#   }
# }

avg.run.time <- mean(run.times) # 18.61
max.run.time <- max(run.times) # 19.98 as of mu = 1E-7, omega = 1, beta TF 6, prot 10, kin 10
median.run.time <- median(run.times) # 18.94

# all.vertices$runID <- paste0("omega_",all.vertices$omega,"_mu_",all.vertices$mu,"_beta_TF_",all.vertices$betaTF, "_prot_",all.vertices$betaTF,"_kin_",all.vertices$betaKin)
# all.edges$runID <- paste0("omega_",all.edges$omega,"_mu_",all.edges$mu,"_beta_TF_",all.edges$betaTF, "_prot_",all.edges$betaTF,"_kin_",all.edges$betaKin)
# all.enrichr$runID <- paste0("omega_",all.enrichr$omega,"_mu_",all.enrichr$mu,"_beta_TF_",all.enrichr$betaTF, "_prot_",all.enrichr$betaTF,"_kin_",all.enrichr$betaKin)
# all.gsea$runID <- paste0("omega_",all.gsea$omega,"_mu_",all.gsea$mu,"_beta_TF_",all.gsea$betaTF, "_prot_",all.gsea$betaTF,"_kin_",all.gsea$betaKin)
# enr.8q24$runID <- paste0("omega_",enr.8q24$omega,"_mu_",enr.8q24$mu,"_beta_TF_",enr.8q24$betaTF, "_prot_",enr.8q24$betaTF,"_kin_",enr.8q24$betaKin)
# runs.so.far <- unique(all.vertices$runID)

# ran out of time with mu = 1E-7, omega = 1, betas 4, 4, 8:
# Error in makeRestartList(...) : reached elapsed time limit
# In addition: There were 12 warnings (use warnings() to see them)
# Storing graph and tabulated data
# Error in pcsf.result$graph : $ operator is invalid for atomic vectors

#### plot results of parameter sweep ####
# create heatmap of parameter sweep
param.info <- dplyr::distinct(all.vertices[,c("Direction", "omega", "mu", "betaTF", "betaProt", "betaKin","runID")])
param.info <- reshape2::melt(param.info, id = "runID", variable.name = "parameter")
param.info <- param.info[order(param.info$runID),] # sort by runID to match order of node degree plot
param.info$rowNum <- seq(1, nrow(param.info))
# inspired by: https://r-graph-gallery.com/79-levelplot-with-ggplot2.html
library(viridis)
param.info2 <- dplyr::distinct(param.info[param.info$parameter != "Direction",])
param.info2$runID2 <- paste0(strsplit(param.info2$runID, "_", fixed = 2)[[1]][-1], collapse = "_")
param.info2 <- param.info2[order(param.info2$runID2),] # sort by runID to match order of node degree plot
param.info2$rowNum2 <- seq(1, nrow(param.info2))
param.plot <- ggplot2::ggplot(param.info2, 
                              aes(x=rowNum2, y = parameter, fill = value)) +
  geom_tile() + scale_fill_gradient(low="white", high = "black") #+ theme_ipsum()

for (j in directions) {
  ### node degrees (steiner vs. terminal) for each parameter setting
  node.types <- all.vertices[all.vertices$Direction == j, 
                             c("name", "type", "runID")]
  degree.info <- merge(node.types, all.centrality, by="runID")
  library(plyr)
  sum.degrees <- plyr::ddply(degree.info, .(type, runID), summarize,
                             degree_sum = sum(degree))
  
  # create line plot of node degrees vs. runID
  deg.plot <- ggplot2::ggplot(sum.degrees, aes(x = runID, y = degree_sum, group = type)) +
    geom_line() + ylab("Sum of Node Degrees")
  
  ### heatmap of nodes in vs. out clustered with parameter setting information
  # convert to wide
  temp.nodes <- reshape2::dcast(all.vertices[all.vertices$Direction = j,], 
                                name ~ runID, value = "Included", fill = "Not Included")
  rownames(temp.nodes) <- temp.nodes$name
  #node.mat <- as.matrix(temp.nodes)
  
  temp.heatmap <- pheatmap::pheatmap(as.matrix(temp.nodes), color = c("white","black"),
                                      cluster_rows = TRUE, cluster_cols = FALSE,
                                      scale = "row", annotation_col = dplyr::distinct(node.types[,c("name","type")]), 
                                      angle_col = "45", show_colnames = TRUE,
                                      fontsize = 6) 
  
  ### % of terminal nodes included in output for each parameter setting
  perc.term <- plyr::ddply(node.types, .(runID), summarize,
                           N_terminal = length(unique(name[type == "Terminal"])))
  if (j == "positive") {
    pos.inputs <- list(mit.kin, rna.allSamples.result, global.result[global.result$Spearman.est > 0,])
    N.inputs <- length(unique(unlist(pos.inputs))) 
  } else {
    N.inputs <- nrow(global.result[global.result$Spearman.est < 0,]) 
  }
  perc.term$N_inputs <- N.inputs
  perc.term$percent_terminal <- perc.term$N_terminal * 100 / perc.term$N_inputs
  
  perc.plot <- ggplot2::ggplot(perc.term, aes(x = runID, y = percent_terminal)) +
    geom_line() + ylab("% Terminal Nodes")
  
  ### chr8q24 enrichment plot
  #temp.8q24 <- enr.8q24[enr.8q24$Direction == j,]
  chr8.plot <- ggplot2::ggplot(enr.8q24[enr.8q24$Direction == j & 
                                          enr.8q24$Feature_set == "chr8q24",],
                               aes(x = runID, y = Feature_set, 
                                   size = -log(p_value, base = 10),
                                   color = NES)) + ggplot2::geom_point() +
    viridis::scale_color_viridis() + 
    ggplot2::labs(color = "NES", size = "-Log(p-value)")
  
  ### combine plots
  library(patchwork)
  combo.plot <- chr8.plot / perc.plot / deg.plot / temp.heatmap / param.plot
  if (j == "positive") {
    pos.plot <- combo.plot
  } else {
    neg.plot <- combo.plot
  }
  ggplot2::ggsave(paste0(j, "_parameter_sweep.pdf"), combo.plot, width = 7, height = 11)
  
  # now do clustering for in/out plot and parameters
  # temp.heatmap2 <- pheatmap::pheatmap(as.matrix(temp.nodes), color = c("white","black"),
  #                                    cluster_rows = TRUE, cluster_cols = TRUE,
  #                                    scale = "row", annotation_col = dplyr::distinct(node.types[,c("name","type")]), 
  #                                    angle_col = "45", show_colnames = TRUE,
  #                                    fontsize = 6)
}

# combine positive and negative network optimization
pos.inputs <- list(mit.kin, rna.allSamples.result, global.result[global.result$Spearman.est > 0,])
n.pos.inputs <- length(unique(unlist(pos.inputs))) 
n.neg.inputs <- nrow(global.result[global.result$Spearman.est < 0,])
n.pos.nodes <- length(unique(all.vertices[all.vertices$Direction == "positive", ]$name))
n.neg.nodes <- length(unique(all.vertices[all.vertices$Direction == "negative", ]$name))
pos.plot <- pos.plot + ggtitle(paste0("Positive Network (", n.pos.nodes, " nodes from ", n.pos.inputs, " inputs)"))
neg.plot <- neg.plot + ggtitle(paste0("Negative Network (", n.neg.nodes, " nodes from ", n.neg.inputs, " inputs)"))
final.combo <- pos.plot + neg.plot
ggplot2::ggsave("parameter_sweep.pdf", final.combo, width = 7, height = 11)


### combo plot (no clustering so that parameter settings are in same order)


