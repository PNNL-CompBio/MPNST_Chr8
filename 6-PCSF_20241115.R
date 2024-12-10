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

library(readxl); library(panSEA); 
library(synapser)
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

processPCSF <- function(pcsf.result, all.dfs, temp.runID, betas, mu, omega, 
                        base.path2, min.per.set = c(3,4,5,6), deg.exp) {
  if (!inherits(pcsf.result, "try-error")) {
    # get data frames from list
    all.centrality <- all.dfs[[1]]
    all.edges <- all.dfs[[2]]
    all.enrichr <- all.dfs[[3]]
    all.gsea <- all.dfs[[4]]
    all.vertices <- all.dfs[[5]]
    enr.8q24 <- all.dfs[[6]]
    
    if (length(betas) == 1) {
      betas <- c(NA, betas, NA)
    } else if (length(betas) == 2) {
      betas <- c(betas, NA)
    }
    
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
    temp.vert$betaTF <- betas[1]
    temp.vert$betaProt <- betas[2]
    temp.vert$betaKin <- betas[3]
    temp.vert[,"omega"] <- omega
    temp.vert[,"degExp"] <- deg.exp
    temp.vert$runID <- temp.runID
    all.vertices <- rbind(all.vertices, temp.vert)
    
    # compile edge data
    temp.edges <- pcsf.tab$edges
    temp.edges$Direction <- j
    temp.edges[,"mu"] <- mu
    temp.edges$betaTF <- betas[1]
    temp.edges$betaProt <- betas[2]
    temp.edges$betaKin <- betas[3]
    temp.edges[,"omega"] <- omega
    temp.edges[,"degExp"] <- deg.exp
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
      temp.enrichr$betaTF <- betas[1]
      temp.enrichr$betaProt <- betas[2]
      temp.enrichr$betaKin <- betas[3]
      temp.enrichr[,"omega"] <- omega
      temp.enrichr[,"degExp"] <- deg.exp
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
        temp.gsea$betaTF <- betas[1]
        temp.gsea$betaProt <- betas[2]
        temp.gsea$betaKin <- betas[3]
        temp.gsea[,"omega"] <- omega
        temp.gsea[,"degExp"] <- deg.exp
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
        saveRDS(GSEA.files, paste0("GSEA_", temp.fname, "_minPerSet_", q, ".rds"))
        #save_to_synapse_v2(GSEA.files, dot.scale = 1)
        
        # extract chr8q24 results
        if ("chr8q24" %in% pos.enr$result$Feature_set) {
          temp.8q24 <- as.data.frame(pos.enr$result[pos.enr$result$Feature_set == "chr8q24",])
          temp.8q24$Direction <- j
          temp.8q24[,"mu"] <- mu
          temp.8q24$betaTF <- betas[1]
          temp.8q24$betaProt <- betas[2]
          temp.8q24$betaKin <- betas[3]
          temp.8q24[,"omega"] <- omega
          temp.8q24[,"degExp"] <- deg.exp
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
  return(list(centr = all.centrality, edge = all.edges, enr = all.enrichr, 
              gsea = all.gsea, vert = all.vertices, q24 = enr.8q24))
}

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
mit.kin <- mit.kin[mit.kin$Spearman.q <= 0.25,c("Gene", "kinase","Spearman.est","Spearman.q")] # 12/381 (3.15%) for q < 0.25; 2 / 381 (0.0525%) if q < 0.05
# when q < 0.05 for mit.kin, only P38G and ERK1 were enriched, neither were in STRING interactome
# 100% of your terminal nodes are included in the interactome
# 96.59% of your terminal nodes are included in the interactome
# 58.33% of your terminal nodes are included in the interactome
# Solving the PCSF by adding random noise to the edge costs...

# replace gene names for kinases
mit.kin[mit.kin$kinase == "P38G",]$Gene <- "MAPK12" # could also be ERK6; source: https://lincs.hms.harvard.edu/db/proteins/200777/
mit.kin[mit.kin$kinase == "ERK1",]$Gene <- "MAPK3" # could also be PRKM3; source: https://www.uniprot.org/uniprotkb/P27361/entry
# LTK should have LTK as gene name but could also be TYK1; source: https://www.uniprot.org/uniprotkb/P29376/entry
mit.kin[mit.kin$kinase == "P38B",]$Gene <- "MAPK11" # source: https://www.uniprot.org/uniprotkb/Q15759/entry
# CDK5 should be CDK5, same for CDK8, 10, 12, 19; MTOR should be MTOR
mit.kin[mit.kin$kinase == "P38A",]$Gene <- "MAPK14" # source: https://www.uniprot.org/uniprotkb/Q16539/entry
mit.kin[mit.kin$kinase == "JNK3",]$Gene <- "MAPK10" # source: https://www.uniprot.org/uniprotkb/P53779/entry
#mit.kin <- mit.kin[mit.kin$Spearman.q <= 0.05,]
mit.kin <- mit.kin[,c("Gene", "Spearman.est")]

# now:
# 100% of your terminal nodes are included in the interactome (TFs)
# 96.59% of your terminal nodes are included in the interactome (protein corrs)
# 100% of your terminal nodes are included in the interactome (kinases)

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

setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
#timeout <- 10*60 # unit of seconds
#temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date(), "_", timeout/60, "_min_timeout")
#temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date(),"_scaled_0.75-1_kinaseAdjP0.05")
#temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date(),"_scaled_0.75-1_kinasesRenamed")
# temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date(),"_scaled_0.75-1_kinasesRenamed_kinaseAdjP0.05")
#temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date(),"_kinasesRenamed_kinaseAdjP0.05")
#temp.path <- paste0("PCSF_TFProteinKinase_", "2024-12-01","_kinasesRenamed_kinaseAdjP0.05")
#temp.path <- paste0("PCSF_TFProteinKinase_", "2024-12-02","_kinasesRenamed_FDRadjP0.05")
#temp.path <- paste0("PCSF_TFProteinKinase_", "2024-12-02","_scaled_0.75-1_kinasesRenamed_FDRadjP0.05")
#temp.path <- paste0("PCSF_TFProteinKinase_", "2024-11-22")
temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date()-1, "_kinasesRenamed")
#temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date(), "_kinasesRenamed")
#temp.path <- paste0("PCSF_TFProtein_", Sys.Date())
#temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date(), "_degExp2")
#temp.path <- paste0("PCSF_TFProteinKinase_2024-12-04_kinasesRenamed")
if (file.exists(temp.path)) {
  setwd(temp.path)
  all.vertices <- read.csv("vertices.csv")
  all.edges <- read.csv("edges.csv")
  all.enrichr <- read.csv("enrichr.csv")
  all.gsea <- read.csv("GSEA_positional.csv")
  enr.8q24 <- try(read.csv("chr8q24_enrichment.csv"),silent=TRUE)
  if (inherits(enr.8q24, "try-error")) {
    enr.8q24 <- data.frame() 
  }
  all.centrality <- read.csv("centrality.csv")
  runs.so.far <- unique(all.vertices$runID)
} else {
  dir.create(temp.path)
  setwd(temp.path)
  all.vertices <- data.frame()
  all.edges <- data.frame()
  all.enrichr <- data.frame()
  all.gsea <- data.frame()
  enr.8q24 <- data.frame()
  all.centrality <- data.frame()
  runs.so.far <- c()
  run.times <- c()
}


base.path2 <- getwd()
#beta.vals <- c(2, 4, 6, 8, 10) # 2 was the default before
#mu.vals <- c(1E-7, 1E-6, 1E-5, 1E-4, 1E-3) # 5E-4 was the default before
#omega.vals <- c(1, 2, 3, 4) # 4 was default before
min.per.set <- c(3,4,5,6)
deg.exp.vals <- c(4,8,1,2)
deg.exp.vals <- 22.3677
deg.exp.vals <- 1
all.dfs <- list(all.centrality, all.edges, all.enrichr, 
                all.gsea, all.vertices, enr.8q24)
beta.vals <- c(5, 6, 11) # next try adding 5
#beta.vals <- c(1E-3, 1, 1E3)
mu.vals <- c(1E3) # next try adding 1, 1E-7, 1E-3, 1E-5
omega.vals <- c(3) # next try adding 1, 2
inputs <- list("Transcription_factor" = rna.allSamples.result, "Protein" = global.result, "Kinase" = mit.kin)
#inputs <- list("Transcription_factor" = rna.allSamples.result, "Protein" = global.result)
#directions <- c("positive", "negative")
directions <- "positive"

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
min.prot.corr <- min(global.result$Spearman.est)
max.prot.corr <- max(global.result$Spearman.est)
#rescale.all <- TRUE
rescale.all <- FALSE

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
      hist(top.phospho$Spearman.est, xlab = names(inputs)[m], main = "")
      print(paste(j, "range of", names(inputs)[m], "is", min(top.phospho$Spearman.est), "to", max(top.phospho$Spearman.est)))
      if (rescale.all) {
        top.phospho$Spearman.est <- scales::rescale(top.phospho$Spearman.est, to = c(0.75,1))
        hist(top.phospho$Spearman.est, xlab = names(inputs)[m], main = "")
        print(paste(j, "range of", names(inputs)[m], "is now", min(top.phospho$Spearman.est), "to", max(top.phospho$Spearman.est)))
      }
      final.inputs[[names(inputs)[m]]] <- as.list(top.phospho$Spearman.est)
      names(final.inputs[[names(inputs)[m]]]) <- top.phospho[,1] # assuming first column contains feature names
    }
  }
  hist(unlist(final.inputs), xlab = paste0(j, " inputs"), main = "")
  
  # parameter sweep: mu, beta.tf, beta.prot, beta.kin, omega
  for (deg.exp in deg.exp.vals) {
    setwd(file.path(base.path2, j))
    dir.create(paste0("degExp_",deg.exp))
    setwd(paste0("degExp_",deg.exp))
    for (mu in mu.vals) {
      setwd(file.path(base.path2, j, paste0("degExp_",deg.exp)))
      dir.create(paste0("mu_",mu))
      setwd(paste0("mu_",mu))
      
      for (omega in omega.vals) {
        setwd(file.path(base.path2, j, paste0("degExp_",deg.exp), paste0("mu_",mu)))
        dir.create(paste0("omega_",omega))
        setwd(paste0("omega_",omega))
        
        if (length(final.inputs) == 3) {
          for (beta.tf in 5) {
            for (beta.prot in 11) {
              for (beta.kin in 6) {
                setwd(file.path(base.path2, j, paste0("degExp_",deg.exp), paste0("mu_",mu), paste0("omega_",omega)))
                
                # run PCSF
                temp.b <- c(beta.tf, beta.prot, beta.kin)
                temp.fname <- paste0("beta_TF_", beta.tf, "_prot_", beta.prot, "_kin_", beta.kin)
                temp.runID <- paste0(j, "_degExp_", deg.exp, "_omega_",omega,"_mu_",mu,"_",temp.fname)
                if (!(temp.runID %in% runs.so.far)) {
                  cat("Running PCSF\n")
                  time.start <- Sys.time()
                  pcsf.result <- try(computeProteinNetwork_BG(final.inputs, 
                                                              betas = temp.b, 
                                                              mu = mu, w = omega,
                                                              deg.exp = deg.exp,
                                                              fname = temp.fname), 
                                     silent = TRUE)
                  runs.so.far <- c(runs.so.far, temp.runID)
                  all.dfs <- processPCSF(pcsf.result, all.dfs, temp.runID, temp.b, mu, omega, base.path2, min.per.set, deg.exp)
                }
              }
            }
          } 
        } else if (length(final.inputs) == 2) {
          for (beta.tf in beta.vals) {
            for (beta.prot in beta.vals) {
              setwd(file.path(base.path2, j, paste0("degExp_",deg.exp), paste0("mu_",mu), paste0("omega_",omega)))
              
              # run PCSF
              temp.b <- c(beta.tf, beta.prot)
              temp.fname <- paste0("beta_TF_", beta.tf, "_prot_", beta.prot)
              temp.runID <- paste0(j, "_degExp_", deg.exp, "_omega_",omega,"_mu_",mu,"_",temp.fname)
              if (!(temp.runID %in% runs.so.far)) {
                cat("Running PCSF\n")
                time.start <- Sys.time()
                pcsf.result <- try(computeProteinNetwork_BG(final.inputs, 
                                                            betas = temp.b, 
                                                            mu = mu, w = omega,
                                                            deg.exp = deg.exp,
                                                            fname = temp.fname), 
                                   silent = TRUE)
                runs.so.far <- c(runs.so.far, temp.runID)
                all.dfs <- processPCSF(pcsf.result, all.dfs, temp.runID, temp.b, mu, omega, base.path2, min.per.set, deg.exp)
              }
            }
          } 
        } else { # negative network only uses global proteomics
          for (beta.prot in beta.vals) {
            setwd(file.path(base.path2, j, paste0("degExp_",deg.exp), paste0("mu_",mu), paste0("omega_",omega)))
            
            # run PCSF
            temp.fname <- paste0("beta_prot_", beta.prot)
            temp.runID <- paste0(j, "_degExp_", deg.exp, "_omega_",omega,"_mu_",mu,"_",temp.fname)
            if (!(temp.runID %in% runs.so.far)) {
              cat("Running PCSF\n")
              time.start <- Sys.time()
              final.input <- list("Protein" = final.inputs[[1]])
              pcsf.result <- try(computeProteinNetwork_BG(final.input, 
                                                          betas = beta.prot, 
                                                          mu = mu, w = omega, 
                                                          deg.exp = deg.exp,
                                                          fname = temp.fname), 
                                 silent = TRUE)
              runs.so.far <- c(runs.so.far, temp.runID)
              all.dfs <- processPCSF(pcsf.result, all.dfs, temp.runID, beta.prot, mu, omega, base.path2, min.per.set, deg.exp)
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

avg.run.time <- mean(run.times) # 18.61; now 20.43
max.run.time <- max(run.times) # 19.98 as of mu = 1E-7, omega = 1, beta TF 6, prot 10, kin 10; now 30
median.run.time <- median(run.times) # 18.94; now 19.53

enr.8q24.v2 <- enr.8q24[enr.8q24$minPerSet == 6,]
enr.8q24.v2$score <- -log(enr.8q24.v2$FDR_q_value, base = 10) * enr.8q24.v2$NES
enr.8q24.v2$rank <- rank(enr.8q24.v2$score) # best is mu = 1E-7, omega = 4, beta TF & prot = 10, beta kin = 2
# all.vertices$runID <- paste0("omega_",all.vertices$omega,"_mu_",all.vertices$mu,"_beta_TF_",all.vertices$betaTF, "_prot_",all.vertices$betaTF,"_kin_",all.vertices$betaKin)
# all.edges$runID <- paste0("omega_",all.edges$omega,"_mu_",all.edges$mu,"_beta_TF_",all.edges$betaTF, "_prot_",all.edges$betaTF,"_kin_",all.edges$betaKin)
# all.enrichr$runID <- paste0("omega_",all.enrichr$omega,"_mu_",all.enrichr$mu,"_beta_TF_",all.enrichr$betaTF, "_prot_",all.enrichr$betaTF,"_kin_",all.enrichr$betaKin)
# all.gsea$runID <- paste0("omega_",all.gsea$omega,"_mu_",all.gsea$mu,"_beta_TF_",all.gsea$betaTF, "_prot_",all.gsea$betaTF,"_kin_",all.gsea$betaKin)
# enr.8q24$runID <- paste0("omega_",enr.8q24$omega,"_mu_",enr.8q24$mu,"_beta_TF_",enr.8q24$betaTF, "_prot_",enr.8q24$betaTF,"_kin_",enr.8q24$betaKin)
# runs.so.far <- unique(all.vertices$runID)
setwd(base.path2)
library(ggplot2)
#### plot results of parameter sweep ####
# sort by rank and then enforce that runID order
for (j in directions) {
  ### parameter info
  all.vertices <- all.dfs[[5]]
  if (!("betaKin" %in% colnames(all.vertices))) {
    all.vertices$betaKin <- NA
  }
  param.info <- dplyr::distinct(all.vertices[all.vertices$Direction == j,c("omega", "mu", "betaTF", "betaProt", "betaKin","degExp","runID")])
  param.info <- reshape2::melt(param.info, id = "runID", variable.name = "parameter")
  param.info <- param.info[order(param.info$runID),] # sort by runID to match order of node degree plot
  param.info$rowNum <- seq(1, nrow(param.info))
  # inspired by: https://r-graph-gallery.com/79-levelplot-with-ggplot2.html
  library(viridis)
  param.info$value <- as.numeric(param.info$value)
  params <- unique(param.info$parameter)
  param.plot <- ggplot2::ggplot(param.info[param.info$parameter == params[1],], 
                                aes(x=rowNum, y = parameter, fill = value)) +
    geom_tile() + scale_fill_gradient(low="white", high = "black") + theme_classic() + 
    ylab(element_blank()) + xlab(element_blank()) + theme(legend.position = "none")
  library(patchwork)
  for (i in params[2:length(params)]) {
    param.plot <- param.plot / (ggplot2::ggplot(param.info[param.info$parameter == i,], 
                                                aes(x=rowNum, y = parameter, fill = value)) +
                                  geom_tile() + scale_fill_gradient(low="white", high = "black") + theme_classic() + 
                                  ylab(element_blank()) + xlab(element_blank()) + theme(legend.position = "none"))
    
  }
  
  ### node degrees (steiner vs. terminal) for each parameter setting
  node.types <- dplyr::distinct(all.vertices[all.vertices$Direction == j, 
                             c("name", "type", "runID")])
  all.centrality <- all.dfs[[1]]
  degree.info <- merge(node.types, dplyr::distinct(all.centrality[,c("name","runID","degree")]), by=c("name","runID"))
  library(plyr)
  rank.summary <- plyr::ddply(degree.info, .(type, runID), summarize,
                             degree_sum = sum(degree))
  rank.summary <- rank.summary[order(rank.summary$runID),]
  
  # create line plot of node degrees vs. runID
  deg.plot <- ggplot2::ggplot(rank.summary, aes(x = runID, y = degree_sum, group = type, color = type)) +
    geom_line() + ylab("Sum of Node Degrees") + theme_classic()
  
  ### heatmap of nodes in vs. out clustered with parameter setting information
  # convert to wide
  temp.nodes <- all.vertices[all.vertices$Direction == j,]
  temp.nodes$Included <- 1
  temp.nodes <- reshape2::dcast(temp.nodes, 
                                name ~ runID, value.var = "Included", fill = 0)
  rownames(temp.nodes) <- temp.nodes$name
  temp.nodes <- temp.nodes[,order(colnames(temp.nodes))]
  #temp.nodes <- temp.nodes %>% mutate_if(is.character, as.numeric) %>% select_if(colSums(.) != 0)
  #temp.nodes <- temp.nodes[,colSums(temp.nodes) > 0]
  node.mat <- as.matrix(temp.nodes[,2:ncol(temp.nodes)])
  node.mat <- node.mat[,colSums(node.mat) > 0]
  
  # temp.heatmap <- pheatmap::pheatmap(node.mat, color = c("white","black"),
  #                                     cluster_rows = FALSE, cluster_cols = FALSE,
  #                                     scale = "row", annotation_col = dplyr::distinct(node.types[,c("name","type")]), 
  #                                     angle_col = "45", show_colnames = TRUE,
  #                                     fontsize = 6) 
  # Error in check.length("fill") : 
  #   'gpar' element 'fill' must not be length 0
  
  ### % of terminal nodes included in output for each parameter setting
  perc.term <- plyr::ddply(node.types, .(runID), summarize,
                           N_terminal = length(unique(name[type == "Terminal"])))
  if (j == "positive") {
    #pos.inputs <- list(mit.kin, rna.allSamples.result, global.result[global.result$Spearman.est > 0,])
    #N.inputs <- length(unique(unlist(pos.inputs))) 
    N.inputs <- 288
  } else {
    N.inputs <- nrow(global.result[global.result$Spearman.est < 0,]) 
  }
  perc.term$N_inputs <- N.inputs
  perc.term$percent_terminal <- perc.term$N_terminal * 100 / perc.term$N_inputs
  perc.term <- perc.term[order(perc.term$runID),]
  
  perc.plot <- ggplot2::ggplot(perc.term, aes(x = runID, y = percent_terminal, group = 1)) +
    geom_line() + ylab("% Terminal Nodes") + theme_classic()
  
  ### chr8q24 enrichment plot - expect negative NES because highest prize would be at bottom of list even though unweighted calculation
  # is the score still affected by input order in unweighted GSEA? yes, well mostly affected by prize (i.e., rank metric) then input order
  # this paper uses betweenness centrality to rank nodes: https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2021.577623/full
  enr.8q24 <- all.dfs[[6]]
  temp.8q24 <- enr.8q24[enr.8q24$Direction == j & 
                          enr.8q24$Feature_set == "chr8q24" & enr.8q24$minPerSet == 6,]
  if (nrow(temp.8q24) > 0) {
    temp.8q24 <- temp.8q24[order(temp.8q24$runID),]
    temp.8q24$chr8q24 <- rank(-log(temp.8q24$FDR_q_value, 10) * temp.8q24$NES)
    chr8.plot <- ggplot2::ggplot(temp.8q24,
                                 aes(x = runID, y = Feature_set, 
                                     size = -log(p_value, base = 10),
                                     color = NES)) + ggplot2::geom_point() +
      viridis::scale_color_viridis() + theme_classic() +
      ggplot2::labs(color = "NES", size = "-Log(p-value)")
  }
  
  ### mean rank (min NES for chr8q24, min diff in node degrees, max percent terminal) - could maybe specify SUMO1, other general nodes to exclude
  rank.df <- data.frame(runID = perc.term$runID,
                        #chr8q24 = rank(-log(temp.8q24$FDR_q_value, 10) * temp.8q24$NES),
                        degree_diff = rank(abs(rank.summary[rank.summary$type == "Steiner",]$degree_sum - rank.summary[rank.summary$type == "Terminal",]$degree_sum)),
                        perc_terminal = rank(perc.term$percent_terminal))
  if (nrow(temp.8q24) > 0) {
    rank.df <- merge(temp.8q24[,c("runID", "chr8q24")], rank.df, by="runID")
    rank.summary <- plyr::ddply(rank.df, .(runID), summarize,
                                chr8q24 = chr8q24,
                                degree_diff = degree_diff,
                                perc_terminal = perc_terminal,
                                rank = mean(c(chr8q24, degree_diff, perc_terminal)))
  } else {
    rank.summary <- plyr::ddply(rank.df, .(runID), summarize,
                                degree_diff = degree_diff,
                                perc_terminal = perc_terminal,
                                rank = mean(c(degree_diff, perc_terminal))) 
  }
  rank.summary <- reshape2::melt(rank.summary, id = "runID", variable.name = "Parameter")
  rank.summary <- rank.summary[order(rank.summary$runID),]
  rank.plot <- ggplot2::ggplot(rank.summary, aes(x = runID, y = value, group = Parameter, color = Parameter)) +
    geom_line() + ylab("Rank") + theme_classic()
  
  ### combine plots
  library(patchwork)
  if (nrow(temp.8q24) > 0) {
    combo.plot <- rank.plot / chr8.plot / perc.plot / deg.plot / param.plot
  } else {
    combo.plot <- rank.plot / perc.plot / deg.plot / param.plot
  }
  #combo.plot <- rank.plot / chr8.plot / perc.plot / deg.plot / temp.heatmap / param.plot
  if (j == "positive") {
    pos.plot <- combo.plot
  } else {
    neg.plot <- combo.plot
  }
  ggplot2::ggsave(paste0(j, "_parameter_sweep.pdf"), combo.plot, width = 11, height = 44)
  saveRDS(combo.plot, paste0(j, "_parameter_sweep.rds"))
  
  # now do clustering for in/out plot and parameters
  # temp.heatmap2 <- pheatmap::pheatmap(as.matrix(temp.nodes), color = c("white","black"),
  #                                    cluster_rows = TRUE, cluster_cols = TRUE,
  #                                    scale = "row", annotation_col = dplyr::distinct(node.types[,c("name","type")]), 
  #                                    angle_col = "45", show_colnames = TRUE,
  #                                    fontsize = 6)
}
# all the 12 kinases are kept; none are from the TFs or correlations:
all.vertices <- all.dfs$vert
kept.terminals <- unique(all.vertices[all.vertices$type == "Terminal",]$name)
kept.kin.terminals <- kept.terminals[kept.terminals %in% mit.kin$Gene]
kept.tf.terminals <- kept.terminals[kept.terminals %in% rna.allSamples.result$Gene]
kept.corr.terminals <- kept.terminals[kept.terminals %in% global.result$Gene]

# maybe this is because outlier LTK (enrichment score of ~5 has adj p 0.07 and low # of hits, though low set sizes)
# this should be fixed by filtering kinase enrichment for 0.05 adj p

# combine positive and negative network optimization
pos.inputs <- list(mit.kin, rna.allSamples.result, global.result[global.result$Spearman.est > 0,])
#n.pos.inputs <- length(unique(unlist(pos.inputs))) 
n.pos.inputs <- 288
n.neg.inputs <- nrow(global.result[global.result$Spearman.est < 0,])
n.pos.nodes <- length(unique(all.vertices[all.vertices$Direction == "positive", ]$name))
n.neg.nodes <- length(unique(all.vertices[all.vertices$Direction == "negative", ]$name))
pos.plot <- pos.plot + ggtitle(paste0("Positive Network (", n.pos.nodes, " nodes from ", n.pos.inputs, " inputs)"))
neg.plot <- neg.plot + ggtitle(paste0("Negative Network (", n.neg.nodes, " nodes from ", n.neg.inputs, " inputs)"))
final.combo <- pos.plot + neg.plot
ggplot2::ggsave("parameter_sweep.pdf", final.combo, width = 22, height = 44)

# look at string network to see why kinases always dominate other proteins
data("STRING")
ppi <- construct_interactome(STRING)
string.centrality.df <- data.frame(name = V(ppi)$name,
                            degree = igraph::degree(ppi, mode="all"),
                            closeness = igraph::closeness(ppi, mode="all"),
                            betweenness = igraph::betweenness(ppi, directed = FALSE),
                            eigen_centrality = igraph::eigen_centrality(ppi, directed = FALSE)$vector,
                            hub_score = igraph::hub_score(ppi)$vector,
                            authority_score = igraph::authority_score(ppi)$vector)

kin.centrality <- string.centrality.df[string.centrality.df$name %in% mit.kin[mit.kin$Spearman.est>0,]$Gene,]
kin.degree <- sum(kin.centrality$degree)
kin.mean.degree <- mean(kin.centrality$degree)
kin.median.degree <- median(kin.centrality$degree)
kin.centrality$Source <- "Kinase enrichment"


TF.centrality <- string.centrality.df[string.centrality.df$name %in% rna.allSamples.result[rna.allSamples.result$Spearman.est>0,]$Gene,]
TF.degree <- sum(TF.centrality$degree)
TF.mean.degree <- mean(TF.centrality$degree)
TF.median.degree <- median(TF.centrality$degree)
TF.centrality$Source <- "Transcription factor enrichment"

corr.centrality <- string.centrality.df[string.centrality.df$name %in% global.result[global.result$Spearman.est>0,]$Gene,]
corr.degree <- sum(corr.centrality$degree)
corr.mean.degree <- mean(corr.centrality$degree)
corr.median.degree <- median(corr.centrality$degree)
corr.centrality$Source <- "Protein correlation"

term.centrality <- rbind(kin.centrality, TF.centrality, corr.centrality)
temp.plot <- ggplot2::ggplot(term.centrality, aes(x=degree, 
                                          group=Source, 
                                          fill=Source)) +
  geom_histogram(position="identity", alpha=0.5) + theme_classic() +
  xlab("Degree in STRING Network")
ggsave(paste0("Chr8_STRING_network_positive_terminal_degree_histogram_", Sys.Date(),".pdf"),temp.plot,
       width=7, height=7)
ggsave(paste0("Chr8_STRING_network_positive_terminal_degree_histogram_wider_", Sys.Date(),".pdf"),temp.plot,
       width=11, height=7)
temp.plot

Source <- unique(term.centrality$Source)
median.degree.df <- data.frame(Source)
median.degree.df$`Median Degree` <- 0
median.degree.df[median.degree.df$Source == "Kinase enrichment",]$`Median Degree` <- kin.median.degree
median.degree.df[median.degree.df$Source == "Transcription factor enrichment",]$`Median Degree` <- TF.median.degree
median.degree.df[median.degree.df$Source == "Protein correlation",]$`Median Degree` <- corr.median.degree
median.degree.df$`SD Degree` <- 0
median.degree.df[median.degree.df$Source == "Kinase enrichment",]$`SD Degree` <- sd(kin.centrality$degree)
median.degree.df[median.degree.df$Source == "Transcription factor enrichment",]$`SD Degree` <- sd(TF.centrality$degree)
median.degree.df[median.degree.df$Source == "Protein correlation",]$`SD Degree` <- sd(corr.centrality$degree)
median.plot <- ggplot2::ggplot(median.degree.df, aes(fill=Source, x=Source, y=`Median Degree`)) + 
  ggplot2::geom_bar(stat="identity", position="dodge") + ggplot2::theme_minimal() +
  geom_errorbar(aes(ymin=`Median Degree` - `SD Degree`, 
                    ymax = `Median Degree` + `SD Degree`), width=0.2,
                position=position_dodge(0.9)) + theme(axis.text.x=element_text(angle=45, vjust=1, hjust = 1))
median.plot
ggsave(paste0("Chr8_STRING_network_positive_terminal_degree_medians_", Sys.Date(),".pdf"),median.plot,
       width=7, height=7)
median.test <- t.test(kin.centrality$degree, c(corr.centrality$degree, TF.centrality$degree),
                      alternative = "greater")
median.test$p.value
# p = 0.00722011 that kinase degrees are higher than other inputs

median.test <- t.test(kin.centrality$degree, corr.centrality$degree,
                      alternative = "greater")
median.test$p.value
# p = 0.006521913 that kinase degrees are higher than corr inputs

median.test <- t.test(kin.centrality$degree, TF.centrality$degree,
                      alternative = "greater")
median.test$p.value
# p = 0.1035931 that kinase degrees are higher than TF inputs

median.test <- t.test(TF.centrality$degree, corr.centrality$degree,
                      alternative = "greater")
median.test$p.value
# p = 0.08557735 that TF degrees are higher than corr inputs

median.test <- t.test(TF.centrality$degree, c(corr.centrality$degree, kin.centrality$degree),
                      alternative = "greater")
median.test$p.value
# p = 0.1089297 that TF degrees are higher than other inputs

median.test <- t.test(corr.centrality$degree, c(TF.centrality$degree, kin.centrality$degree),
                      alternative = "greater")
median.test$p.value
# p = 0.997663 that corr degrees are higher than other inputs

term.centrality$Included <- FALSE
term.centrality[term.centrality$name %in% all.vertices[all.vertices$Direction=="positive",]$name,]$Included <- TRUE
temp.plot2 <- ggplot2::ggplot(term.centrality[term.centrality$Included & term.centrality$Source != "Kinase enrichment",], aes(x=degree, 
                                                  group=Source, 
                                                  fill=Source)) +
  geom_histogram(position="identity", alpha=0.5) + theme_classic() +
  xlab("Degree in STRING Network")
ggsave(paste0("Chr8_STRING_network_positive_included_terminal_degree_histogram_", Sys.Date(),".pdf"),temp.plot2,
       width=7, height=7)
temp.plot2
### combo plot (no clustering so that parameter settings are in same order)

# output initial prizes
rna.allSamples.result$Source <- "TF"
mit.kin$Source <- "Kinase"
global.result$Source <- "Correlation"
global.result <- global.result[,c("Gene", "Spearman.est", "Source")]
final.inputs.pos <- rbind(rna.allSamples.result, mit.kin, global.result)
saveRDS(final.inputs, "Positive_inputs.rds")
write.csv(final.inputs.pos, "Positive_inputs.csv", row.names = FALSE)

optFunc <- function(params, mode = "default", betaTF, betaProt, betaKin,
                    mu, deg.exp) {
  terminals <- readRDS("Positive_inputs.rds")
  if (length(params) == 1) {
    if (mode == "betaKin") {
      betaKin <- params
    } else if (mode == "betaTF") {
      betaTF <- params
    } else if (mode == "betaProt") {
      betaProt <- params
    } else if (mode == "mu") {
      mu <- params
    } else if (mode == "degExp") {
      deg.exp <- params 
    } else {
      stop("invalid mode for just 1 parameter entry")
    }
    b <- c(betaTF,betaProt,betaKin)
  } else if (length(params) == 5) {
    b <- params[1:3]
    mu <- params[4]
    deg.exp <- params[5]
  } else {
    stop("invalid number of parameters")
  }
  cat("Calculating prizes with betas",b, "mu", mu, "degExp",deg.exp,"\n")
  
  # load string network
  require(dplyr)
  if(!require('PCSF')){
    remotes::install_github('sgosline/PCSF')
    require('PCSF')
  }
  data("STRING")
  ppi <- construct_interactome(STRING)
  
  # calculate prizes for terminal nodes
  terminal_types <- names(terminals)
  node_names = V(ppi)$name
  node_prz = vector(mode = "numeric", length = length(node_names)) # 0 vector
  for (j in 1:length(terminal_types)) {
    terminal_names = names(terminals[[j]])
    terminal_values = as.numeric(terminals[[j]])
    
    # Incorporate the node prizes
    index = match(terminal_names, node_names)
    percent = signif((length(index) - sum(is.na(index)))/length(index)*100, 4)
    if (percent < 5)
      stop("  Less than 1% of your terminal nodes are matched in the interactome, check your terminals!")
    cat(paste0("  ", percent, "% of your terminal nodes are included in the interactome\n"))
    terminal_names = terminal_names[!is.na(index)]
    terminal_values = abs(terminal_values[!is.na(index)])
    index = index[!is.na(index)]
    node_prz[index] =  node_prz[index] + b[j]*terminal_values # in case node is in more than 1 terminal input
  }
  
  # Calculate the hub penalization scores
  node_degrees = igraph::degree(ppi)
  hub_penalization = - mu*(node_degrees^deg.exp)
  
  # Update the node prizes
  node_prizes = node_prz
  index = which(node_prizes==0)
  node_prizes[index] = hub_penalization[index]
  
  ## assign node types
  terms=unlist(terminals)
  lfcs<-terms[match(names(V(ppi)),names(terms))]
  lfcs[is.na(lfcs)]<-0.0
  
  types<-rep('Steiner',length(names(V(ppi))))
  names(types)<-names(V(ppi))
  
  ##assign node types (e.g., kinase, protein, TF)
  val.types <- names(terminals)
  types.so.far <- c()
  for (i in 1:length(val.types)) {
    types[intersect(names(V(ppi)),names(terminals[[i]]))]<-val.types[i]
    
    types.so.far <- c(types.so.far, val.types[i])
    other.types <- val.types[!(val.types %in% types.so.far)]
    if (length(other.types) > 0) {
      for (j in 1:length(other.types)) {
        temp.other.vals <- names(terminals[other.types[j]])
        types[intersect(names(V(ppi)),
                        intersect(names(terminals[[i]]),temp.other.vals))] <- 
          paste0(types[intersect(names(V(ppi)),
                                 intersect(names(terminals[[i]]),temp.other.vals))],
                 "+",other.types[j])
      } 
    }
  }
  kin.indices <- which(grepl("Kinase", types))
  other.indices <- which(!grepl("Kinase", types))
  
  score <- 1-t.test(node_prizes[kin.indices], node_prizes[other.indices], alternative="greater")$p.value
  cat("p =", -as.numeric(1-score), "\n")
  return(score)
}
optim.test <- optim(c(1,10,1,1E-3,22.3677), optFunc, method = "L-BFGS-B", lower = 0, hessian = TRUE)
optim.test <- optim(c(1,10,1,1E-3,22.3677), optFunc, hessian = TRUE)
optim.test <- optim(1, optFunc, method = "Brent", lower = 0, upper = 2^5, hessian = TRUE) # last run was 22.3677 with score 0.9092 aka p = 0.09082033
optim.test <- optim(22.3677, optFunc, method = "Brent", lower = 0, upper = 2^5, hessian = TRUE)
testResult <- optFunc(22.3677) # 0.9092 aka p = 0.09082033

# now optimize betaKin
#top.limit <- optim.test$par # 200
mode.types <- c("degExp","mu","betaKin","betaTF", "betaProt")
mode.vals <- c(22.3677, 1000, 1, 1, 10) # initially started with degExp = 1, mu = 1E-3
for (i in 2:length(mode.types)) {
  cat("optimizing", mode.types[i],"\n")
  optim.test <- list("par" = mode.vals[i])
  top.limit <- mode.vals[i]
  while (ceiling(optim.test$par) == top.limit | 
         abs(top.limit - optim.test$par) < 2) {
    top.limit <- top.limit*10
    optim.test <- try(optim(optim.test$par, optFunc, method = "Brent", lower = 0, 
                        upper = top.limit, hessian = TRUE, mode = mode.types[i], 
                        betaTF = mode.vals[4], betaProt = mode.vals[5],
                        betaKin = mode.vals[3], mu = mode.vals[2],
                        deg.exp = mode.vals[1]))
    if (inherits(optim.test, "try-error")) {
      break
    }
  } 
  if (inherits(optim.test, "try-error")) {
    optim.test <- list("par" = top.limit)
  }
  mode.vals[i] <- optim.test$par
  cat("p =", as.numeric(1-optim.test$value), "after optimizing", mode.types[i], "to", mode.vals[i],"\n")
}
mode.vals
mode.vals <- c(22.3677, 10000.0000, 10.0000, 10.0000, 100.0000) # after 1 try of above loop
mode.vals <- c(22.3677, 1E5, 100, 100, 1000) # after 2nd try
mode.vals <- c(22.3677, 1E6, 100, 100, 1000) # after 3rd try but error after 1 round of mu

# for betaKin: first tried 32, then 100, then 200 but kept approaching limit
#optim.test <- optim(1E5, optFunc, method = "Brent", lower = 0, upper = 1E6, hessian = TRUE, mode = "betaKin") # first tried 32, then 100, then 200 but kept approaching limit
1-optim.test$value # p=0.09082033
optim.test <- optFunc(c(mode.vals[4], mode.vals[5], mode.vals[3], mode.vals[2], mode.vals[1])) # first tried 32, then 100, then 200 but kept approaching limit
# returns NA

# try Rcgmin
packageurl <- "http://cran.r-project.org/src/contrib/Archive/Rcgmin/Rcgmin_2022-4.30.tar.gz"
install.packages("optextras")
install.packages(packageurl, repos=NULL, type="source")

optF <- function(params) {
  terminals <- readRDS("Positive_inputs.rds")
  if (length(params) == 5) {
    b <- params[1:3]
    mu <- params[4]
    deg.exp <- params[5]
  } else {
    stop("invalid number of parameters")
  }
  cat("Calculating prizes with betas",b, "mu", mu, "degExp",deg.exp,"\n")
  
  # load string network
  require(dplyr)
  if(!require('PCSF')){
    remotes::install_github('sgosline/PCSF')
    require('PCSF')
  }
  data("STRING")
  ppi <- construct_interactome(STRING)
  
  # calculate prizes for terminal nodes
  terminal_types <- names(terminals)
  node_names = V(ppi)$name
  node_prz = vector(mode = "numeric", length = length(node_names)) # 0 vector
  for (j in 1:length(terminal_types)) {
    terminal_names = names(terminals[[j]])
    terminal_values = as.numeric(terminals[[j]])
    
    # Incorporate the node prizes
    index = match(terminal_names, node_names)
    percent = signif((length(index) - sum(is.na(index)))/length(index)*100, 4)
    if (percent < 5)
      stop("  Less than 1% of your terminal nodes are matched in the interactome, check your terminals!")
    cat(paste0("  ", percent, "% of your terminal nodes are included in the interactome\n"))
    terminal_names = terminal_names[!is.na(index)]
    terminal_values = abs(terminal_values[!is.na(index)])
    index = index[!is.na(index)]
    node_prz[index] =  node_prz[index] + b[j]*terminal_values # in case node is in more than 1 terminal input
  }
  
  # Calculate the hub penalization scores
  node_degrees = igraph::degree(ppi)
  hub_penalization = - mu*(node_degrees^deg.exp)
  
  # Update the node prizes
  node_prizes = node_prz
  index = which(node_prizes==0)
  node_prizes[index] = hub_penalization[index]
  
  ## assign node types
  terms=unlist(terminals)
  lfcs<-terms[match(names(V(ppi)),names(terms))]
  lfcs[is.na(lfcs)]<-0.0
  
  types<-rep('Steiner',length(names(V(ppi))))
  names(types)<-names(V(ppi))
  
  ##assign node types (e.g., kinase, protein, TF)
  val.types <- names(terminals)
  types.so.far <- c()
  for (i in 1:length(val.types)) {
    types[intersect(names(V(ppi)),names(terminals[[i]]))]<-val.types[i]
    
    types.so.far <- c(types.so.far, val.types[i])
    other.types <- val.types[!(val.types %in% types.so.far)]
    if (length(other.types) > 0) {
      for (j in 1:length(other.types)) {
        temp.other.vals <- names(terminals[other.types[j]])
        types[intersect(names(V(ppi)),
                        intersect(names(terminals[[i]]),temp.other.vals))] <- 
          paste0(types[intersect(names(V(ppi)),
                                 intersect(names(terminals[[i]]),temp.other.vals))],
                 "+",other.types[j])
      } 
    }
  }
  kin.indices <- which(grepl("Kinase", types))
  other.indices <- which(!grepl("Kinase", types))
  
  score <- abs(mean(node_prizes[kin.indices]) - mean(node_prizes[other.indices]))
  cat("absDiff =", score, "\n")
  return(score)
}

# calculates gradient of prize function instead of typical prizes
optG <- function(params) {
  terminals <- readRDS("Positive_inputs.rds")
  if (length(params) == 5) {
    b <- params[1:3]
    mu <- params[4]
    deg.exp <- params[5]
  } else {
    stop("invalid number of parameters")
  }
  cat("Calculating prizes with betas",b, "mu", mu, "degExp",deg.exp,"\n")
  
  # load string network
  require(dplyr)
  if(!require('PCSF')){
    remotes::install_github('sgosline/PCSF')
    require('PCSF')
  }
  data("STRING")
  ppi <- construct_interactome(STRING)
  
  # calculate prizes for terminal nodes
  terminal_types <- names(terminals)
  node_names = V(ppi)$name
  node_prz = vector(mode = "numeric", length = length(node_names)) # 0 vector
  for (j in 1:length(terminal_types)) {
    terminal_names = names(terminals[[j]])
    terminal_values = as.numeric(terminals[[j]])
    
    # Incorporate the node prizes
    index = match(terminal_names, node_names)
    percent = signif((length(index) - sum(is.na(index)))/length(index)*100, 4)
    if (percent < 5)
      stop("  Less than 1% of your terminal nodes are matched in the interactome, check your terminals!")
    cat(paste0("  ", percent, "% of your terminal nodes are included in the interactome\n"))
    terminal_names = terminal_names[!is.na(index)]
    terminal_values = abs(terminal_values[!is.na(index)])
    index = index[!is.na(index)]
    node_prz[index] =  node_prz[index] + terminal_values # in case node is in more than 1 terminal input
  }
  
  # Calculate the hub penalization scores
  node_degrees = igraph::degree(ppi)
  hub_penalization = - mu*(node_degrees^deg.exp)*log(node_degrees) - (node_degrees^deg.exp)
  
  # Update the node prizes
  node_prizes = node_prz
  index = which(node_prizes==0)
  node_prizes[index] = hub_penalization[index]
  
  ## assign node types
  terms=unlist(terminals)
  lfcs<-terms[match(names(V(ppi)),names(terms))]
  lfcs[is.na(lfcs)]<-0.0
  
  types<-rep('Steiner',length(names(V(ppi))))
  names(types)<-names(V(ppi))
  
  ##assign node types (e.g., kinase, protein, TF)
  val.types <- names(terminals)
  types.so.far <- c()
  for (i in 1:length(val.types)) {
    types[intersect(names(V(ppi)),names(terminals[[i]]))]<-val.types[i]
    
    types.so.far <- c(types.so.far, val.types[i])
    other.types <- val.types[!(val.types %in% types.so.far)]
    if (length(other.types) > 0) {
      for (j in 1:length(other.types)) {
        temp.other.vals <- names(terminals[other.types[j]])
        types[intersect(names(V(ppi)),
                        intersect(names(terminals[[i]]),temp.other.vals))] <- 
          paste0(types[intersect(names(V(ppi)),
                                 intersect(names(terminals[[i]]),temp.other.vals))],
                 "+",other.types[j])
      } 
    }
  }
  kin.indices <- which(grepl("Kinase", types))
  other.indices <- which(!grepl("Kinase", types))
  
  score <- abs(mean(node_prizes[kin.indices]) - mean(node_prizes[other.indices]))
  cat("absDiffGradient =", score, "\n")
  return(score)
}

mode.types <- c("betaTF", "betaProt","betaKin","mu", "degExp")
mode.vals <- c(5,10,1,10,10) # initially started with degExp = 1, mu = 1E-3
optim.test <- Rcgmin::Rcgmin(par=mode.vals, fn=optF, gr=optG, 
                             #lower=rep(0, length(mode.vals)),
                             #upper=rep(1E9,length(mode.vals)),
                             bdmsk = rep(1,length(mode.vals)))
# last try was: 
mode.vals <- c(-6.164362e+33, -6.164362e+33, -6.164362e+33, -6.164362e+33, -6.164362e+33)
optim.test <- Rcgmin::Rcgmin(par=mode.vals, fn=optF, gr=optG, 
                             bdmsk = rep(1,length(mode.vals)))

# Calculating prizes with betas -6.164362e+33 -6.164362e+33 -6.164362e+33 mu -6.164362e+33 degExp -6.164362e+33 
# 100% of your terminal nodes are included in the interactome
# 96.59% of your terminal nodes are included in the interactome
# 100% of your terminal nodes are included in the interactome
# absDiff = 9.140631e+33 
# Calculating prizes with betas -6.164362e+33 -6.164362e+33 -6.164362e+33 mu -6.164362e+33 degExp -6.164362e+33 
# 100% of your terminal nodes are included in the interactome
# 96.59% of your terminal nodes are included in the interactome
# 100% of your terminal nodes are included in the interactome
# absDiff = 1.482819 

mode.vals <- rep(1,5)
optim.test <- Rcgmin::Rcgmin(par=mode.vals, fn=optF, gr=optG, 
                             bdmsk = rep(1,length(mode.vals)))
# result: -19.28594 -19.28594 -19.28594 -19.28594 -19.28594 

# try nleqslv
install.packages("nleqslv")

optComparison <- function(params, comparison) {
  terminals <- readRDS("Positive_inputs.rds")
  if (length(params) == 5) {
    b <- params[1:3]
    mu <- params[4]
    deg.exp <- params[5]
  } else {
    stop("invalid number of parameters")
  }
  cat("Calculating prizes with betas",b, "mu", mu, "degExp",deg.exp,"\n")
  
  # load string network
  require(dplyr)
  if(!require('PCSF')){
    remotes::install_github('sgosline/PCSF')
    require('PCSF')
  }
  data("STRING")
  ppi <- construct_interactome(STRING)
  
  # calculate prizes for terminal nodes
  terminal_types <- names(terminals)
  node_names = V(ppi)$name
  node_prz = vector(mode = "numeric", length = length(node_names)) # 0 vector
  for (j in 1:length(terminal_types)) {
    terminal_names = names(terminals[[j]])
    terminal_values = as.numeric(terminals[[j]])
    
    # Incorporate the node prizes
    index = match(terminal_names, node_names)
    percent = signif((length(index) - sum(is.na(index)))/length(index)*100, 4)
    if (percent < 5)
      stop("  Less than 1% of your terminal nodes are matched in the interactome, check your terminals!")
    cat(paste0("  ", percent, "% of your terminal nodes are included in the interactome\n"))
    terminal_names = terminal_names[!is.na(index)]
    terminal_values = abs(terminal_values[!is.na(index)])
    index = index[!is.na(index)]
    node_prz[index] =  node_prz[index] + b[j]*terminal_values # in case node is in more than 1 terminal input
  }
  
  # Calculate the hub penalization scores
  node_degrees = igraph::degree(ppi)
  hub_penalization = - mu*(node_degrees^deg.exp)
  
  # Update the node prizes
  node_prizes = node_prz
  index = which(node_prizes==0)
  node_prizes[index] = hub_penalization[index]
  
  ## assign node types
  terms=unlist(terminals)
  lfcs<-terms[match(names(V(ppi)),names(terms))]
  lfcs[is.na(lfcs)]<-0.0
  
  types<-rep('Steiner',length(names(V(ppi))))
  names(types)<-names(V(ppi))
  
  ##assign node types (e.g., kinase, protein, TF)
  val.types <- names(terminals)
  types.so.far <- c()
  for (i in 1:length(val.types)) {
    types[intersect(names(V(ppi)),names(terminals[[i]]))]<-val.types[i]
    
    types.so.far <- c(types.so.far, val.types[i])
    other.types <- val.types[!(val.types %in% types.so.far)]
    if (length(other.types) > 0) {
      for (j in 1:length(other.types)) {
        temp.other.vals <- names(terminals[other.types[j]])
        types[intersect(names(V(ppi)),
                        intersect(names(terminals[[i]]),temp.other.vals))] <- 
          paste0(types[intersect(names(V(ppi)),
                                 intersect(names(terminals[[i]]),temp.other.vals))],
                 "+",other.types[j])
      } 
    }
  }
  kin.indices <- which(grepl(comparison[2], types))
  other.indices <- which(grepl(comparison[1], types))
  
  score <- abs(mean(node_prizes[kin.indices]) - mean(node_prizes[other.indices]))
  cat("absDiff between", comparison[1], "and", comparison[2], "=", score, "\n")
  return(score)
}

optSingle <- function(params, type) {
  terminals <- readRDS("Positive_inputs.rds")
  if (length(params) == 5) {
    b <- params[1:3]
    mu <- params[4]
    deg.exp <- params[5]
  } else {
    stop("invalid number of parameters")
  }
  cat("Calculating prizes with betas",b, "mu", mu, "degExp",deg.exp,"\n")
  
  # load string network
  require(dplyr)
  if(!require('PCSF')){
    remotes::install_github('sgosline/PCSF')
    require('PCSF')
  }
  data("STRING")
  ppi <- construct_interactome(STRING)
  
  # calculate prizes for terminal nodes
  terminal_types <- names(terminals)
  node_names = V(ppi)$name
  node_prz = vector(mode = "numeric", length = length(node_names)) # 0 vector
  for (j in 1:length(terminal_types)) {
    terminal_names = names(terminals[[j]])
    terminal_values = as.numeric(terminals[[j]])
    
    # Incorporate the node prizes
    index = match(terminal_names, node_names)
    percent = signif((length(index) - sum(is.na(index)))/length(index)*100, 4)
    if (percent < 5)
      stop("  Less than 1% of your terminal nodes are matched in the interactome, check your terminals!")
    cat(paste0("  ", percent, "% of your terminal nodes are included in the interactome\n"))
    terminal_names = terminal_names[!is.na(index)]
    terminal_values = abs(terminal_values[!is.na(index)])
    index = index[!is.na(index)]
    node_prz[index] =  node_prz[index] + terminal_values # in case node is in more than 1 terminal input
  }
  
  # Calculate the hub penalization scores
  node_degrees = igraph::degree(ppi)
  hub_penalization = - mu*(node_degrees^deg.exp)*log(node_degrees) - (node_degrees^deg.exp)
  
  # Update the node prizes
  node_prizes = node_prz
  index = which(node_prizes==0)
  node_prizes[index] = hub_penalization[index]
  
  ## assign node types
  terms=unlist(terminals)
  lfcs<-terms[match(names(V(ppi)),names(terms))]
  lfcs[is.na(lfcs)]<-0.0
  
  types<-rep('Steiner',length(names(V(ppi))))
  names(types)<-names(V(ppi))
  
  ##assign node types (e.g., kinase, protein, TF)
  val.types <- names(terminals)
  types.so.far <- c()
  for (i in 1:length(val.types)) {
    types[intersect(names(V(ppi)),names(terminals[[i]]))]<-val.types[i]
    
    types.so.far <- c(types.so.far, val.types[i])
    other.types <- val.types[!(val.types %in% types.so.far)]
    if (length(other.types) > 0) {
      for (j in 1:length(other.types)) {
        temp.other.vals <- names(terminals[other.types[j]])
        types[intersect(names(V(ppi)),
                        intersect(names(terminals[[i]]),temp.other.vals))] <- 
          paste0(types[intersect(names(V(ppi)),
                                 intersect(names(terminals[[i]]),temp.other.vals))],
                 "+",other.types[j])
      } 
    }
  }
  kin.indices <- which(grepl(type, types))
  
  score <- mean(node_prizes[kin.indices])
  cat("Gradient for", type, "=", score, "\n")
  return(score)
}

soe <- function(params) {
  diffEqTF <- optComparison(params, comparison = c("Transcription_factor", "Kinase"))
  diffEqCorr <- optComparison(params, comparison = c("Protein", "Kinase"))
  gradEqTF <- optSingle(params, type = "Transcription_factor")
  gradEqCorr <- optSingle(params, type = "Protein")
  gradEqKin <- optSingle(params, type = "Kinase")
  return(c(diffEqTF, diffEqCorr, gradEqTF, gradEqCorr, gradEqKin))
}

optim.test.soe <- nleqslv::nleqslv(mode.vals, soe, control=list("allowSingular"=TRUE))
# last run: Calculating prizes with betas 0.6063449 1.325895 0.7378821 mu 1 degExp 1
# absDiff between Transcription_factor and Kinase = 5.150322e-10 
# absDiff between Protein and Kinase = 4.133756e-10  
# absDiffGradient = 1.82 
# absDiffGradient = 0.8323039
# absDiffGradient = 1.495561 

mode.vals <- optim.test.soe$x
optim.test.soe2 <- nleqslv::nleqslv(mode.vals, soe, control=list("allowSingular"=TRUE))
saveRDS(optim.test.soe, "Chr8_PCSF_optimization_nleqslv.rds")
mode.vals # 0.6063449 1.3258952 0.7378821 1.0000000 1.0000000
optim.test.soe$x # 0.6063449 1.3258952 0.7378821 1.0000000 1.0000000
optim.test.soe2$x # 0.6063449 1.3258952 0.7378821 1.0000000 1.0000000
optim.test.soe$fvec # 5.150322e-10 4.133756e-10 1.820000e+00 8.323039e-01 1.495561e+00
optim.test.soe2$fven # 3.404592e-10 1.350362e-10 1.820000e+00 8.323039e-01 1.495561e+00
optim.test.soe$termcd # 3
optim.test.soe2$termcd # 3
optim.test.soe$message # "No better point found (algorithm has stalled)"
optim.test.soe2$message # "No better point found (algorithm has stalled)"
optim.test.soe$scalex # 1 1 1 1 1
optim.test.soe2$scalex # 1 1 1 1 1
optim.test.soe$nfcnt # 13
optim.test.soe2$nfcnt # 1
optim.test.soe$njcnt # 2
optim.test.soe2$njcnt # 1
optim.test.soe$iter # 3
optim.test.soe2$iter # 1

# same beta ratios would be TF = 5, prot = 11, kin = 6
mode.vals <- c(5, 11, 6, 1, 1)
optim.test.soe3 <- nleqslv::nleqslv(mode.vals, soe, control=list("allowSingular"=TRUE))
optim.test.soe3$x #  5.005833 10.946261  6.091771  1.000000  1.000000
saveRDS(optim.test.soe3, "Chr8_PCSF_optimization_nleqslv3.rds")

# try maintaining betas and only optimize mu, degExp for min diff
soe_mu_deg <- function(params) {
  diffEqTF <- optComparison(c(5,11,6,params), comparison = c("Transcription_factor", "Kinase"))
  diffEqCorr <- optComparison(c(5,11,6,params), comparison = c("Protein", "Kinase"))
  return(c(diffEqTF, diffEqCorr))
}
optim.test.soe3 <- nleqslv::nleqslv(mode.vals[4:5], soe_mu_deg, control=list("allowSingular"=TRUE))
optim.test.soe3$message # "Jacobian is completely unusable (all zero entries?)"

optim.test.soe3 <- nleqslv::nleqslv(c(1E-3,1), soe_mu_deg, control=list("allowSingular"=TRUE))
optim.test.soe3$message # "Jacobian is completely unusable (all zero entries?)"

# try maintaining mu, n and only optimize beta ratios for min diff
soe_betas <- function(params) {
  diffEqTF <- optComparison(c(params[1],1,params[2],1E-3,1), comparison = c("Transcription_factor", "Kinase"))
  diffEqCorr <- optComparison(c(params[1],1,params[2],1E-3,1), comparison = c("Protein", "Kinase"))
  return(c(diffEqTF, diffEqCorr))
}
mode.vals <- c(0.6,0.7)
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/PCSF_TFProteinKinase_2024-12-06_kinasesRenamed")
optim.test.soe.b <- nleqslv::nleqslv(mode.vals, soe_betas, control=list("allowSingular"=TRUE))
optim.test.soe.b$message # 'Function criterion near zero.' Convergence of function values has been achieved.
mode.vals <- optim.test.soe.b$x # 0.4573098 0.5565161
optim.test.soe.b$termcd # 1
optim.test.soe.b$scalex # 1 1
optim.test.soe.b$nfcnt # 1
optim.test.soe.b$njcnt # 1
optim.test.soe.b$iter # 1

optim.test.soe.b2 <- nleqslv::nleqslv(mode.vals, soe_betas, control=list("allowSingular"=TRUE))
optim.test.soe.b2$message # 'Function criterion near zero.' Convergence of function values has been achieved.
mode.vals <- optim.test.soe.b2$x # 0.4573098 0.5565161
optim.test.soe.b2$termcd # 1
optim.test.soe.b2$scalex # 1 1
optim.test.soe.b2$nfcnt # 0
optim.test.soe.b2$njcnt # 0
optim.test.soe.b2$iter # 0

mode.vals*11 # 5.030408 6.121678

# try to optimize steiner nodes via mu, degExp now
getPCSFResults <- function(all.dfs, j = "positive") {
  all.vertices <- all.dfs[[5]]
  if (!("betaKin" %in% colnames(all.vertices))) {
    all.vertices$betaKin <- NA
  }
  param.info <- dplyr::distinct(all.vertices[all.vertices$Direction == j,c("omega", "mu", "betaTF", "betaProt", "betaKin","degExp","runID")])
  param.info <- reshape2::melt(param.info, id = "runID", variable.name = "parameter")
  param.info <- param.info[order(param.info$runID),] # sort by runID to match order of node degree plot
  param.info$rowNum <- seq(1, nrow(param.info))
  # inspired by: https://r-graph-gallery.com/79-levelplot-with-ggplot2.html
  library(viridis)
  param.info$value <- as.numeric(param.info$value)
  params <- unique(param.info$parameter)
  param.plot <- ggplot2::ggplot(param.info[param.info$parameter == params[1],], 
                                aes(x=rowNum, y = parameter, fill = value)) +
    geom_tile() + scale_fill_gradient(low="white", high = "black") + theme_classic() + 
    ylab(element_blank()) + xlab(element_blank()) + theme(legend.position = "none")
  library(patchwork)
  for (i in params[2:length(params)]) {
    param.plot <- param.plot / (ggplot2::ggplot(param.info[param.info$parameter == i,], 
                                                aes(x=rowNum, y = parameter, fill = value)) +
                                  geom_tile() + scale_fill_gradient(low="white", high = "black") + theme_classic() + 
                                  ylab(element_blank()) + xlab(element_blank()) + theme(legend.position = "none"))
    
  }
  
  ### node degrees (steiner vs. terminal) for each parameter setting
  node.types <- dplyr::distinct(all.vertices[all.vertices$Direction == j, 
                                             c("name", "type", "runID")])
  all.centrality <- all.dfs[[1]]
  degree.info <- merge(node.types, dplyr::distinct(all.centrality[,c("name","runID","degree")]), by=c("name","runID"))
  library(plyr)
  rank.summary <- plyr::ddply(degree.info, .(type, runID), summarize,
                              degree_sum = sum(degree))
  rank.summary <- rank.summary[order(rank.summary$runID),]
  
  # create line plot of node degrees vs. runID
  deg.plot <- ggplot2::ggplot(rank.summary, aes(x = runID, y = degree_sum, group = type, color = type)) +
    geom_line() + ylab("Sum of Node Degrees") + theme_classic()
  
  ### heatmap of nodes in vs. out clustered with parameter setting information
  # convert to wide
  temp.nodes <- all.vertices[all.vertices$Direction == j,]
  temp.nodes$Included <- 1
  temp.nodes <- reshape2::dcast(temp.nodes, 
                                name ~ runID, value.var = "Included", fill = 0)
  rownames(temp.nodes) <- temp.nodes$name
  temp.nodes <- temp.nodes[,order(colnames(temp.nodes))]
  #temp.nodes <- temp.nodes %>% mutate_if(is.character, as.numeric) %>% select_if(colSums(.) != 0)
  #temp.nodes <- temp.nodes[,colSums(temp.nodes) > 0]
  node.mat <- as.matrix(temp.nodes[,2:ncol(temp.nodes)])
  node.mat <- node.mat[,colSums(node.mat) > 0]
  
  # temp.heatmap <- pheatmap::pheatmap(node.mat, color = c("white","black"),
  #                                     cluster_rows = FALSE, cluster_cols = FALSE,
  #                                     scale = "row", annotation_col = dplyr::distinct(node.types[,c("name","type")]), 
  #                                     angle_col = "45", show_colnames = TRUE,
  #                                     fontsize = 6) 
  # Error in check.length("fill") : 
  #   'gpar' element 'fill' must not be length 0
  
  ### % of terminal nodes included in output for each parameter setting
  perc.term <- plyr::ddply(node.types, .(runID), summarize,
                           N_terminal = length(unique(name[type == "Terminal"])))
  if (j == "positive") {
    #pos.inputs <- list(mit.kin, rna.allSamples.result, global.result[global.result$Spearman.est > 0,])
    #N.inputs <- length(unique(unlist(pos.inputs))) 
    N.inputs <- 288
  } else {
    N.inputs <- nrow(global.result[global.result$Spearman.est < 0,]) 
  }
  perc.term$N_inputs <- N.inputs
  perc.term$percent_terminal <- perc.term$N_terminal * 100 / perc.term$N_inputs
  perc.term <- perc.term[order(perc.term$runID),]
  
  perc.plot <- ggplot2::ggplot(perc.term, aes(x = runID, y = percent_terminal, group = 1)) +
    geom_line() + ylab("% Terminal Nodes") + theme_classic()
  
  ### chr8q24 enrichment plot - expect negative NES because highest prize would be at bottom of list even though unweighted calculation
  # is the score still affected by input order in unweighted GSEA? yes, well mostly affected by prize (i.e., rank metric) then input order
  # this paper uses betweenness centrality to rank nodes: https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2021.577623/full
  enr.8q24 <- all.dfs[[6]]
  temp.8q24 <- enr.8q24[enr.8q24$Direction == j & 
                          enr.8q24$Feature_set == "chr8q24" & enr.8q24$minPerSet == 6,]
  if (nrow(temp.8q24) > 0) {
    temp.8q24 <- temp.8q24[order(temp.8q24$runID),]
    temp.8q24$chr8q24 <- rank(-log(temp.8q24$FDR_q_value, 10) * temp.8q24$NES)
    chr8.plot <- ggplot2::ggplot(temp.8q24,
                                 aes(x = runID, y = Feature_set, 
                                     size = -log(p_value, base = 10),
                                     color = NES)) + ggplot2::geom_point() +
      viridis::scale_color_viridis() + theme_classic() +
      ggplot2::labs(color = "NES", size = "-Log(p-value)")
  }
  
  ### mean rank (min NES for chr8q24, min diff in node degrees, max percent terminal) - could maybe specify SUMO1, other general nodes to exclude
  rank.df <- data.frame(runID = perc.term$runID,
                        #chr8q24 = rank(-log(temp.8q24$FDR_q_value, 10) * temp.8q24$NES),
                        degree_diff = rank(abs(rank.summary[rank.summary$type == "Steiner",]$degree_sum - rank.summary[rank.summary$type == "Terminal",]$degree_sum)),
                        perc_terminal = rank(1-perc.term$percent_terminal))
  if (nrow(temp.8q24) > 0) {
    rank.df <- merge(temp.8q24[,c("runID", "chr8q24")], rank.df, by="runID")
    rank.summary2 <- plyr::ddply(rank.df, .(runID), summarize,
                                chr8q24 = chr8q24,
                                degree_diff = degree_diff,
                                perc_terminal = perc_terminal,
                                rank = mean(c(chr8q24, degree_diff, perc_terminal)))
  } else {
    rank.summary2 <- plyr::ddply(rank.df, .(runID), summarize,
                                degree_diff = degree_diff,
                                perc_terminal = perc_terminal,
                                rank = mean(c(degree_diff, perc_terminal))) 
  }
  rank.summary2 <- reshape2::melt(rank.summary2, id = "runID", variable.name = "Parameter")
  rank.summary2 <- rank.summary2[order(rank.summary2$runID),]
  rank.plot <- ggplot2::ggplot(rank.summary2, aes(x = runID, y = value, group = Parameter, color = Parameter)) +
    geom_line() + ylab("Rank") + theme_classic()
  
  ### combine plots
  library(patchwork)
  if (nrow(temp.8q24) > 0) {
    combo.plot <- rank.plot / chr8.plot / perc.plot / deg.plot / param.plot
  } else {
    combo.plot <- rank.plot / perc.plot / deg.plot / param.plot
  }

  ggplot2::ggsave(paste0(j, "_parameter_sweep.pdf"), combo.plot, width = 11, height = 44)
  saveRDS(combo.plot, paste0(j, "_parameter_sweep.rds"))
  return(list(q24 = temp.8q24$chr8q24,
              perc.term = (1-perc.term$percent_terminal),
              degree.diff = abs(rank.summary[rank.summary$type == "Steiner",]$degree_sum - rank.summary[rank.summary$type == "Terminal",]$degree_sum)))
}

soe_steiner <- function(params) {
  base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/PCSF_TFProteinKinase_2024-12-06_kinasesRenamed"
  all.dfs <- list(data.frame(), data.frame(), data.frame(), data.frame(),
                  data.frame())
  temp.fname <- paste0("beta_TF_", 5, "_prot_", 11, "_kin_", 6)
  temp.runID <- paste0("positive", "_degExp_", params[2], "_omega_",params[3],"_mu_",params[1],"_",temp.fname)
  final.inputs <- readRDS(file.path(base.path,"Positive_inputs.rds"))
  cat("Running PCSF with mu =", params[1], "degExp =", params[2], "omega =", params[3], "\n")
  pcsf.result <- computeProteinNetwork_BG(final.inputs, 
                                              betas = c(5,11,6), 
                                              mu = params[1], w = params[3],
                                              deg.exp = params[2],
                                              fname = temp.runID)
  all.dfs <- processPCSF(pcsf.result, all.dfs, temp.runID, c(5,11,6), params[1],
                         params[3], getwd(), deg.exp = params[2])
  rankSummary <- getPCSFResults(all.dfs)
  
  #diffEqTF <- optComparison(c(5,11,6,params), comparison = c("Transcription_factor", "Kinase"))
  #diffEqCorr <- optComparison(c(5,11,6,params), comparison = c("Protein", "Kinase"))
  cat("Results:", -(1-rankSummary$perc.term),"terminals used,\n",
      rankSummary$degree.diff,"absDiff in sum of degrees,\n",
      #diffEqTF, "absDiff in kinase & TF prizes,\n",
      #diffEqCorr,"absDiff in corr & TF prizes,\n",
      rankSummary$q24, "chr8q24 enrichment strength\n")
  return(c(rankSummary$perc.term, rankSummary$degree.diff, 
           #diffEqTF, diffEqCorr, 
           rankSummary$q24))
}
mode.vals <- c(1E-3,1,3)
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/PCSF_TFProteinKinase_2024-12-06_kinasesRenamed")
optim.soe <- nleqslv::nleqslv(mode.vals, soe_steiner, control=list("allowSingular"=TRUE))
optim.test.soe$message # 'Function criterion near zero.' Convergence of function values has been achieved.
mode.vals <- optim.test.soe$x # 0.4573098 0.5565161
optim.test.soe$termcd # 1
optim.test.soe$scalex # 1 1
optim.test.soe$nfcnt # 1
optim.test.soe$njcnt # 1
optim.test.soe.b$iter # 1

optim.test.soe2 <- nleqslv::nleqslv(mode.vals, soe_steiner, control=list("allowSingular"=TRUE))
optim.test.soe2$message # 'Function criterion near zero.' Convergence of function values has been achieved.
mode.vals <- optim.test.soe2$x # 0.4573098 0.5565161
optim.test.soe2$termcd # 1
optim.test.soe2$scalex # 1 1
optim.test.soe2$nfcnt # 0
optim.test.soe2$njcnt # 0
optim.test.soe2$iter # 0
