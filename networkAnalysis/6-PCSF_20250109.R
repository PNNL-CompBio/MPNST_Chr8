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
library(plyr); library(dplyr); library(R.utils); library(ggplot2)
#webshot::install_phantomjs()
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
source("https://github.com/PNNL-CompBio/MPNST_Chr8/blob/main/panSEA_helper_20240913.R")
source("https://github.com/PNNL-CompBio/MPNST_Chr8/blob/main/networkAnalysis/proteinNetworks_BG.R")
source("https://github.com/PNNL-CompBio/MPNST_Chr8/blob/main/networkAnalysis/PCSF_rand_BG.R")

posEnr <- function(all.vertices, gmt.pos) {
  all.runIDs <- unique(all.vertices$runID)
  all.gsea <- data.frame()
  enr.8q24 <- data.frame()
  for (i in all.runIDs) {
    # prep for GSEA
    temp.vertices <- all.vertices[all.vertices$runID == i,]
    temp.vertices <- temp.vertices[order(temp.vertices$frequency),] 
    
    # positional enrichment analysis
    cat("Running GSEA\n")
    for (q in 6) {
      pos.enr <- try(panSEA::drugSEA_ties(temp.vertices, 
                                          gmt = gmt.pos, drug="name", 
                                          rank.metric="frequency", 
                                          stat.type="Classic", 
                                          ties = TRUE,
                                          min.per.set = q), silent = TRUE)
      if (!inherits(pos.enr, "try-error")) {
        colnames(pos.enr$result)[2] <- "Feature_set"
        
        # compile GSEA results
        temp.gsea <- as.data.frame(pos.enr$result)
        temp.gsea$minPerSet <- q
        temp.gsea$runID <- i
        all.gsea <- rbind(all.gsea, temp.gsea)
        
        # extract chr8q24 results
        if ("chr8q24" %in% pos.enr$result$Feature_set) {
          temp.8q24 <- as.data.frame(pos.enr$result[pos.enr$result$Feature_set == "chr8q24",])
          temp.8q24$minPerSet <- q
          temp.8q24$runID <- i
          enr.8q24 <- rbind(enr.8q24, temp.8q24)
        }
      }
    } 
  }
  #write.csv(all.gsea,paste0("GSEA_positional_unweighted_",Sys.Date(),".csv"), row.names = FALSE)
  #write.csv(enr.8q24,paste0("chr8q24_enrichment_unweighted_",Sys.Date(),".csv"), row.names = FALSE)
  return(list(gsea = all.gsea, q24 = enr.8q24))
}
redidPosEnr <- posEnr(all.vertices, gmt.pos)

processPCSF <- function(pcsf.result, all.dfs, temp.runID, beta, mu, omega, 
                        base.path2, deg.exp, j = "positive", temp.fname) {
  if (!inherits(pcsf.result, "try-error")) {
    # get data frames from list
    all.centrality <- all.dfs[[1]]
    all.edges <- all.dfs[[2]]
    all.enrichr <- all.dfs[[3]]
    all.gsea <- all.dfs[[4]]
    all.vertices <- all.dfs[[5]]
    enr.8q24 <- all.dfs[[6]]

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
    temp.vert[,"beta"] <- beta
    temp.vert[,"omega"] <- omega
    temp.vert[,"degExp"] <- deg.exp
    temp.vert$runID <- temp.runID
    all.vertices <- rbind(all.vertices, temp.vert)
    
    # compile edge data
    temp.edges <- pcsf.tab$edges
    temp.edges$Direction <- j
    temp.edges[,"mu"] <- mu
    temp.edges[,"beta"] <- beta
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
      temp.enrichr[,"beta"] <- beta
      temp.enrichr[,"omega"] <- omega
      temp.enrichr[,"degExp"] <- deg.exp
      temp.enrichr$runID <- temp.runID
      all.enrichr <- rbind(all.enrichr, temp.enrichr) 
    }
    
    # positional enrichment analysis
    cat("Running GSEA\n")
    EAresults <- posEnr(pcsf.tab$vertices, gmt.pos)
    
    # store files after each run
    setwd(base.path2)
    write.csv(all.vertices,"vertices.csv", row.names = FALSE)
    write.csv(all.edges,"edges.csv", row.names = FALSE)
    write.csv(all.enrichr,"enrichr.csv", row.names = FALSE)
    write.csv(EAresults$gsea,"GSEA_positional.csv", row.names = FALSE)
    write.csv(EAresults$q24,"chr8q24_enrichment.csv", row.names = FALSE)
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
#temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date(), "_kinasesRenamed")

deg.exp.vals <- 1
beta.vals <- c(1,5,10)
mu.vals <- c(1E-7, 1E-5, 1E-3) # next try adding 1, 1E-7, 1E-3, 1E-5
omega.vals <- c(1,2,3) # next try adding 1, 2
#inputs <- list("Transcription_factor" = rna.allSamples.result, "Protein" = global.result, "Kinase" = mit.kin)
#inputs <- list("Protein" = global.result)
inputs <- list("Transcription_factor" = rna.allSamples.result, "Kinase" = mit.kin)
directions <- c("positive", "negative")

for (m in 1:length(inputs)) {
  setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
  temp.path <- paste0("PCSF_", names(inputs)[m], "_", Sys.Date())
  #temp.path <- paste0("PCSF_Protein_", "2024-12-12")
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
  all.dfs <- list(all.centrality, all.edges, all.enrichr, 
                  all.gsea, all.vertices, enr.8q24)
  
  for (j in directions) {
    setwd(base.path2)
    dir.create(j)
    setwd(j)
    
    # prep named feature lists
    cat("Preparing input\n")
    final.inputs <- list()
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
            
            for (beta in beta.vals) {
              setwd(file.path(base.path2, j, paste0("degExp_",deg.exp), paste0("mu_",mu), paste0("omega_",omega)))
              
              # run PCSF
              temp.fname <- paste0("beta_", beta)
              temp.runID <- paste0(j, "_degExp_", deg.exp, "_omega_",omega,"_mu_",mu,"_",temp.fname)
              if (!(temp.runID %in% runs.so.far)) {
                cat("Running PCSF\n")
                pcsf.result <- try(computeProteinNetwork_BG(final.inputs, 
                                                        betas = beta, 
                                                        mu = mu, w = omega,
                                                        deg.exp = deg.exp,
                                                        fname = temp.fname), silent=TRUE)
                runs.so.far <- c(runs.so.far, temp.runID)
                all.dfs <- processPCSF(pcsf.result, all.dfs, temp.runID, 
                                       beta, mu, omega, base.path2, deg.exp, j, 
                                       temp.fname)
              }
            }
          } 
        }
      }
    }
  } 
}

library(ggplot2)
#### see why kinases always dominate other proteins ####
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

#### plot all results ####
library(patchwork)

# sort by rank and then enforce that runID order
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
vert.prot <- read.csv("PCSF_Protein_2024-12-12/vertices.csv")
vert.prot$inputType <- "Protein"
vert.tf <- read.csv("PCSF_Transcription_factor_2025-01-10/vertices.csv")
vert.tf$inputType <- "Transcription Factor"
vert.tf$betaTF <- vert.tf$beta
vert.kin <- read.csv("PCSF_Kinase_2025-01-10/vertices.csv")
vert.kin$inputType <- "Kinase"
vert.kin$betaKin <- vert.kin$beta

keepCols <- c(colnames(vert.prot),"betaTF","betaKin") # should not include "beta" but instead "betaTF", "betaProt", "betaKin"
vert.prot <- vert.prot[vert.prot$Direction == "positive",keepCols[keepCols %in% colnames(vert.prot)]]
vert.tf <- vert.tf[vert.tf$Direction == "positive",keepCols[keepCols %in% colnames(vert.tf)]]
vert.kin <- vert.kin[vert.kin$Direction == "positive",keepCols[keepCols %in% colnames(vert.kin)]]
minVertCols <- c("name","type","mu","omega")
all.vertices <- merge(vert.tf[,c(minVertCols,"betaTF")], 
                      vert.prot[,c(minVertCols,"betaProt")], 
                      by=minVertCols, all = TRUE)
all.vertices <- merge(all.vertices[,c(minVertCols,"betaTF","betaProt")], 
                      vert.kin[,c(minVertCols,"betaKin")], 
                      by=minVertCols, all = TRUE)

param.info <- as.data.frame(tidyr::expand_grid(omega=omega.vals, mu=mu.vals, 
                                               betaTF=beta.vals, 
                                               betaProt=beta.vals, 
                                               betaKin=beta.vals))

centr.prot <- read.csv("PCSF_Protein_2024-12-12/centrality.csv")
centr.prot$inputType <- "Protein"
centr.prot$runID <- sub("_prot_","_",centr.prot$runID)
centr.tf <- read.csv("PCSF_Transcription_factor_2025-01-10/centrality.csv")
centr.tf$inputType <- "Transcription Factor"
centr.kin <- read.csv("PCSF_Kinase_2025-01-10/centrality.csv")
centr.kin$inputType <- "Kinase"
keepCols <- colnames(centr.prot)
centr.tf <- centr.tf[,keepCols]
centr.kin <- centr.kin[,keepCols]
all.centrality <- rbind(centr.tf, centr.prot, centr.kin)
all.centrality$runID2 <- sub(".*tive_degExp_1_","",all.centrality$runID)
all.centrality$Direction <- sub("\\_.*","",all.centrality$runID)

prot.8q24 <- read.csv("PCSF_Protein_2024-12-12/chr8q24_enrichment.csv")
tf.8q24 <- read.csv("PCSF_Transcription_factor_2025-01-10/chr8q24_enrichment.csv") # NA
kin.8q24 <- read.csv("PCSF_Kinase_2025-01-10/chr8q24_enrichment.csv") # NA
enr.8q24 <- prot.8q24[prot.8q24$Direction=="positive",]
enr.8q24$runID2 <- sub(".*tive_degExp_1_","",enr.8q24$runID)

# just focus on positive since there are no negative terminals for TFs or kinases
all.centrality <- all.centrality[all.centrality$Direction == "positive",
                                 c("name","degree","runID2","inputType")]|> 
  tidyr::separate_wider_delim(runID2, delim="_", names=c("var1","omega","var2", "mu","var3","beta"))
all.centrality[,c("var1","var2","var3")] <- NULL
all.centrality[,c("betaTF","betaProt","betaKin")] <- NA
all.centrality[all.centrality$inputType=="Transcription Factor",]$betaTF <- 
  all.centrality[all.centrality$inputType=="Transcription Factor",]$beta
all.centrality[all.centrality$inputType=="Protein",]$betaProt <- 
  all.centrality[all.centrality$inputType=="Protein",]$beta
all.centrality[all.centrality$inputType=="Kinase",]$betaKin <- 
  all.centrality[all.centrality$inputType=="Kinase",]$beta
all.centrality$beta <- NULL

### node degrees (steiner vs. terminal) for each parameter setting
degree.info <- dplyr::distinct(merge(all.vertices, all.centrality))
degree.summary <- plyr::ddply(degree.info, .(type, omega, mu, betaTF, betaProt, betaKin), summarize,
                              degree_sum = sum(degree, na.rm=TRUE))

### % of terminal nodes included in output for each parameter setting
perc.term <- plyr::ddply(all.vertices[all.vertices$type=="Terminal",], 
                         .(omega, mu, betaTF, betaProt, betaKin), summarize,
                         N_terminal = length(unique(name)))
#perc.term$N_inputs <- 441 # number of negatively correlated proteins with adj. p <= 0.05 (no negatively enriched TFs or kinases)
perc.term$N_inputs <- 264
perc.term[is.na(perc.term$betaProt),]$N_inputs <- 12
perc.term$percent_terminal <- perc.term$N_terminal * 100 / perc.term$N_inputs
perc.term$runID2 <- paste0("omega_", perc.term$omega, "_mu_", perc.term$mu,
                           "_beta_TF_",perc.term$betaTF, "_prot_",perc.term$betaProt, 
                           "_kin_",perc.term$betaKin)

### chr8q24 enrichment plot - expect negative NES because highest prize would be at bottom of list even though unweighted calculation
# is the score still affected by input order in unweighted GSEA? yes, well mostly affected by prize (i.e., rank metric) then input order
# this paper uses betweenness centrality to rank nodes: https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2021.577623/full
temp.8q24 <- enr.8q24[enr.8q24$Feature_set == "chr8q24" & enr.8q24$minPerSet == 6,]
if (nrow(temp.8q24) > 0) {
  if (any(temp.8q24$FDR_q_value == 0)) {
    temp.8q24[temp.8q24$FDR_q_value == 0,]$FDR_q_value <- 1E-4
  }
  temp.8q24$chr8q24 <- rank(-log(temp.8q24$FDR_q_value, 10) * temp.8q24$NES)
  temp.8q24$runID2 <- sub(".*tive_degExp_1_","",temp.8q24$runID)
}

### mean rank (min NES for chr8q24, min diff in node degrees, max percent terminal) - could maybe specify SUMO1, other general nodes to exclude
degree.summary2 <- reshape2::dcast(degree.summary, omega+mu+betaTF+betaProt+betaKin ~ type, value.var = "degree_sum")
degree.summary2$runID2 <- paste0("omega_", degree.summary2$omega, "_mu_", degree.summary2$mu,
                                 "_beta_TF_",degree.summary2$betaTF, "_prot_",degree.summary2$betaProt, 
                                 "_kin_",degree.summary2$betaKin)
perc.term2 <- plyr::ddply(perc.term, .(runID2), summarize,
                          perc_terminal = mean(percent_terminal, na.rm = TRUE))
rank.df <- dplyr::distinct(merge(perc.term2, degree.summary2, by="runID2"))
rank.df$degree_diff <- abs(rank.df$Steiner - rank.df$Terminal)
rank.df$rank_degDiff <- rank(rank.df$degree_diff)
rank.df$rank_percTerm <- rank(100-rank.df$perc_terminal)
rank.summary <- dplyr::distinct(rank.df[,c("runID2","rank_degDiff", "rank_percTerm")])
temp.8q24 <- data.frame()
if (nrow(temp.8q24) > 0) {
  rank.summary2 <- merge(temp.8q24[,c("runID2", "chr8q24")], rank.summary, by="runID2")
  rank.summary2 <- plyr::ddply(rank.summary2, .(runID2), summarize,
                               degree_diff = mean(rank_degDiff, na.rm = TRUE),
                               perc_terminal = mean(rank_percTerm, na.rm = TRUE),
                               chr8q24 = mean(chr8q24, na.rm = TRUE),
                               rank = mean(c(rank_degDiff, rank_percTerm, chr8q24), na.rm = TRUE))
} else {
  rank.summary2 <- plyr::ddply(rank.summary, .(runID2), summarize,
                               degree_diff = mean(rank_degDiff, na.rm = TRUE),
                               perc_terminal = mean(rank_percTerm, na.rm = TRUE),
                               rank = mean(c(rank_degDiff, rank_percTerm), na.rm = TRUE))
}
runOrder <- rank.summary2[order(rank.summary2$rank),]$runID2
rank.summary3 <- reshape2::melt(rank.summary2, id = "runID2", variable.name = "Parameter")
library(ggplot2)
rank.plot <- ggplot2::ggplot(rank.summary3, 
                             aes(x = runID2, y = value, group = Parameter, 
                                 color = Parameter)) +
  geom_line() + ylab("Rank") + theme_classic() +
  scale_color_discrete(labels=c("Chr8q24\nenrichment", "Difference\nin degree","% terminal\nnodes used","Overall"))+
  ggplot2::scale_x_discrete(limits = runOrder) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
rank.plot <- ggplot2::ggplot(rank.summary3, 
                             aes(x = runID2, y = value, group = Parameter, 
                                 color = Parameter)) +
  geom_line() + ylab("Rank") + theme_classic() +
  scale_color_discrete(labels=c("Difference\nin degree","% terminal\nnodes used","Overall"))+
  ggplot2::scale_x_discrete(limits = runOrder) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
rank.plot

### create plots now that runID overall ranks are known
# create line plot of node degrees vs. runID
#runOrder <- sort(unique(rank.summary$runID))
deg.plot <- ggplot2::ggplot(degree.summary, aes(x = runID2, y = degree_sum, group = type, color = type)) +
  geom_line() + ylab("Sum of Node Degrees") + theme_classic() + 
  ggplot2::scale_x_discrete(limits = runOrder) +
  labs(color="Node Type")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
deg.plot

# percent terminal nodes used
perc.plot <- ggplot2::ggplot(perc.term, aes(x = runID2, y = percent_terminal, group = 1)) +
  geom_line() + ylab("% Terminal Nodes") + theme_classic() + 
  ggplot2::scale_x_discrete(limits = runOrder)+
  scale_y_continuous(limits=c(94,97)) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
perc.plot

# chr8q24 enrichment
if (nrow(temp.8q24) > 0) {
  temp.8q24$sig <- FALSE
  temp.8q24$Direction <- "Negative"
  temp.8q24[grepl("positive",temp.8q24$runID),]$Direction <- "Positive"
  if (any(temp.8q24$p_value <= 0.05 & temp.8q24&FDR_q_value <= 0.25)) {
    temp.8q24[temp.8q24$p_value <= 0.05 & temp.8q24$FDR_q_value<=0.25,]$sig <- TRUE
  }
  maxAbsNES <- max(abs(temp.8q24$NES))
  chr8.plot <- ggplot2::ggplot(temp.8q24,
                               aes(x = runID2, y = Direction, 
                                   size = -log(FDR_q_value, base = 10),
                                   color = NES, na.rm=TRUE)) + ggplot2::geom_point() +
    scale_color_gradient2(low="blue",high="red", mid="grey",limits=c(-maxAbsNES, maxAbsNES)) + theme_classic() +
    ggplot2::labs(color = "chr8q24 NES", size = "-log(FDR)") + 
    ggplot2::scale_x_discrete(limits = runOrder) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    geom_point(data = subset(temp.8q24, sig), col = "black", stroke = 1.5, shape = 21)
  chr8.plot
}

### parameter info
param.info$runID2 <- paste0("omega_", param.info$omega, "_mu_", param.info$mu,
                            "_beta_TF_",param.info$betaTF, "_prot_",param.info$betaProt, 
                            "_kin_",param.info$betaKin)
param.info <- reshape2::melt(param.info, id = "runID2", variable.name = "parameter")
# inspired by: https://r-graph-gallery.com/79-levelplot-with-ggplot2.html
param.info$value <- as.numeric(param.info$value)
params <- unique(param.info$parameter)
temp.param.info <- param.info[param.info$parameter == params[1],]
maxParam <- max(temp.param.info$value)
minParam <- min(temp.param.info$value)
temp.param.info$Value <- ""
temp.param.info[temp.param.info$value == maxParam,]$Value <- "High"
temp.param.info[temp.param.info$value == minParam,]$Value <- "Low"
param.plot <- ggplot2::ggplot(temp.param.info, aes(x=runID2, y = parameter, fill = Value)) +
  geom_tile() + scale_fill_manual(breaks=c("High","","Low"),values=c("black","grey","white")) + theme_classic() + 
  ylab(element_blank()) + xlab(element_blank()) + 
  ggplot2::scale_x_discrete(limits = runOrder) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())#+ 
  #scale_fill_discrete(labels=c("Low", "High"))
param.plot
library(patchwork)
for (i in params[2:length(params)]) {
  temp.param.info <- na.omit(param.info[param.info$parameter == i,])
  if (nrow(temp.param.info) > 0) {
    maxParam <- max(temp.param.info$value)
    minParam <- min(temp.param.info$value)
    temp.param.info$Value <- ""
    temp.param.info[temp.param.info$value == maxParam,]$Value <- "High"
    temp.param.info[temp.param.info$value == minParam,]$Value <- "Low"
    
    param.plot <- param.plot / (ggplot2::ggplot(temp.param.info, aes(x=runID2, y = parameter, fill = Value)) +
                                  geom_tile() + scale_fill_manual(breaks=c("High","","Low"),values=c("black","grey","white")) +
                                  theme_classic() + ylab(element_blank()) + xlab(element_blank()) +
                                  theme(legend.position = "none") + 
                                  ggplot2::scale_x_discrete(limits = runOrder)) +
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()) 
  }
}
param.plot
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
param.plot <- param.plot + plot_layout(guides = 'collect')
param.plot

### combine plots
library(patchwork)
if (nrow(temp.8q24) > 0) {
  combo.plot <- rank.plot / deg.plot / perc.plot / chr8.plot / param.plot
} else {
  combo.plot <- rank.plot / perc.plot / deg.plot / param.plot
}
combo.plot
ggplot2::ggsave(paste0("parameter_sweep_v3_",Sys.Date(),".pdf"), combo.plot, width = 6, height = 16)
saveRDS(combo.plot, paste0("parameter_sweep_v3_",Sys.Date(),".rds"))

#### plot all results - wo chr8q24, heatmap ####
library(patchwork)

### node degrees (steiner vs. terminal) for each parameter setting
node.types <- dplyr::distinct(all.vertices[,c("name", "type", "Direction", "runID")])
node.types$runID2 <- sub(".*tive_degExp_1_","",node.types$runID)
degree.info <- merge(node.types, dplyr::distinct(all.centrality[,c("name","runID","degree")]), by=c("name","runID"))
library(plyr)
degree.summary <- plyr::ddply(degree.info, .(type, runID2), summarize,
                              degree_sum = sum(degree))

### % of terminal nodes included in output for each parameter setting
perc.term <- plyr::ddply(node.types, .(Direction, runID2), summarize,
                         N_terminal = length(unique(name[type == "Terminal"])))
N.pos.inputs <- nrow(global.result[global.result$Spearman.est > 0,]) # 264
N.neg.inputs <- nrow(global.result[global.result$Spearman.est < 0,]) # 441
perc.term$N_inputs <- N.neg.inputs
perc.term[perc.term$Direction == "positive",]$N_inputs <- N.pos.inputs
perc.term$percent_terminal <- perc.term$N_terminal * 100 / perc.term$N_inputs

### mean rank (min NES for chr8q24, min diff in node degrees, max percent terminal) - could maybe specify SUMO1, other general nodes to exclude
rank.df <- data.frame(runID2 = perc.term$runID2,
                      degree_diff = rank(abs(degree.summary[degree.summary$type == "Steiner",]$degree_sum - 
                                               degree.summary[degree.summary$type == "Terminal",]$degree_sum)),
                      perc_terminal = rank(100-perc.term$percent_terminal))
rank.summary <- plyr::ddply(rank.df, .(runID2), summarize,
                            degree_diff = mean(degree_diff, na.rm = TRUE),
                            perc_terminal = mean(perc_terminal, na.rm = TRUE),
                            rank = mean(c(degree_diff, perc_terminal), na.rm = TRUE)) 
runOrder <- rank.summary[order(rank.summary$rank),]$runID2
rank.summary <- reshape2::melt(rank.summary, id = "runID2", variable.name = "Metric")
rank.summary$Metric <- factor(rank.summary$Metric, levels=c("rank","degree_diff","perc_terminal"))
# make sure levels, labels for rank.plot, and order of plots in patchwork match
rank.plot <- ggplot2::ggplot(rank.summary, 
                             aes(x = runID2, y = value, group = Metric, 
                                 color = Metric)) +
  geom_line() + ylab("Rank") + theme_classic() +
  scale_color_discrete(labels=c("Overall","Difference\nin degree","% terminal\nnodes used"))+
  ggplot2::scale_x_discrete(limits = runOrder) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
rank.plot

### create plots now that runID overall ranks are known
# create line plot of node degrees vs. runID
#runOrder <- sort(unique(rank.summary$runID))
deg.plot <- ggplot2::ggplot(degree.summary, aes(x = runID2, y = degree_sum, group = type, color = type)) +
  geom_line() + ylab("Sum of Node Degrees") + theme_classic() + 
  ggplot2::scale_x_discrete(limits = runOrder) +
  labs(color="Node Type")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
deg.plot

# percent terminal nodes used
perc.plot <- ggplot2::ggplot(perc.term, aes(x = runID2, y = percent_terminal, group = 1)) +
  geom_line() + ylab("% Terminal Nodes") + theme_classic() + 
  ggplot2::scale_x_discrete(limits = runOrder)+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
perc.plot

### parameter info
param.info <- dplyr::distinct(all.vertices[,c("omega", "mu", "betaProt", "runID")])
param.info$runID2 <- sub(".*tive_degExp_1_","",param.info$runID)
param.info <- dplyr::distinct(param.info[,c("runID2","omega","mu","betaProt")])
colnames(param.info)[4] <- "beta"
param.info <- reshape2::melt(param.info, id = "runID2", variable.name = "parameter")
# inspired by: https://r-graph-gallery.com/79-levelplot-with-ggplot2.html
param.info$value <- as.numeric(param.info$value)
params <- unique(param.info$parameter)
temp.param.info <- param.info[param.info$parameter == params[1],]
maxParam <- max(temp.param.info$value)
minParam <- min(temp.param.info$value)
temp.param.info$Value <- ""
temp.param.info[temp.param.info$value == maxParam,]$Value <- "High"
temp.param.info[temp.param.info$value == minParam,]$Value <- "Low"
param.plot <- ggplot2::ggplot(temp.param.info, aes(x=runID2, y = parameter, fill = Value)) +
  geom_tile() + scale_fill_manual(breaks=c("High","","Low"),values=c("black","grey","white")) + theme_classic() + 
  ylab(element_blank()) + xlab(element_blank()) + 
  ggplot2::scale_x_discrete(limits = runOrder) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())#+ 
#scale_fill_discrete(labels=c("Low", "High"))
param.plot
library(patchwork)
for (i in params[2:length(params)]) {
  temp.param.info <- na.omit(param.info[param.info$parameter == i,])
  if (nrow(temp.param.info) > 0) {
    maxParam <- max(temp.param.info$value)
    minParam <- min(temp.param.info$value)
    temp.param.info$Value <- ""
    temp.param.info[temp.param.info$value == maxParam,]$Value <- "High"
    temp.param.info[temp.param.info$value == minParam,]$Value <- "Low"
    
    param.plot <- param.plot / (ggplot2::ggplot(temp.param.info, aes(x=runID2, y = parameter, fill = Value)) +
                                  geom_tile() + scale_fill_manual(breaks=c("High","","Low"),values=c("black","grey","white")) +
                                  theme_classic() + ylab(element_blank()) + xlab(element_blank()) +
                                  theme(legend.position = "none") + 
                                  ggplot2::scale_x_discrete(limits = runOrder)) +
      theme(axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()) 
  }
}
param.plot
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
param.plot <- param.plot + plot_layout(guides = 'collect')
param.plot

### combine plots
library(patchwork)
combo.plot <- rank.plot / deg.plot / perc.plot / param.plot
combo.plot
ggplot2::ggsave(paste0("parameter_sweep_woChr8q24Ranks_",Sys.Date(),".pdf"), combo.plot, width = 4, height = 7)
saveRDS(combo.plot, paste0("parameter_sweep_woChr8q24Ranks_",Sys.Date(),".rds"))
write.csv(perc.term,"PercentTerminalNodesUsed.csv", row.names = FALSE)
write.csv(degree.summary,"NodeTypeDegrees.csv", row.names = FALSE)
write.csv(rank.summary,"Ranks.csv", row.names = FALSE)
write.csv(rank.df,"Rank_info.csv", row.names = FALSE)

#### save optimal run ####
selectedRun <- runOrder[1] # "omega_2_mu_1e-05_beta_10"
selectedRun <- "omega_2_mu_1e-05_beta_10"

vert.prot <- read.csv("PCSF_Protein_2024-12-12/vertices.csv")
vert.prot$inputType <- "Protein"
vert.prot$beta <- vert.prot$betaProt
vert.prot[,c("betaTF","betaKin","betaProt")] <- NULL
vert.prot <- vert.prot[vert.prot$runID == "positive_degExp_1_omega_2_mu_1e-05_beta_prot_10",]
vert.tf <- read.csv("PCSF_Transcription_factor_2025-01-10/vertices.csv")
vert.tf$inputType <- "Transcription Factor"
vert.tf <- vert.tf[vert.tf$runID == "positive_degExp_1_omega_2_mu_1e-05_beta_10",]
vert.kin <- read.csv("PCSF_Kinase_2025-01-10/vertices.csv")
vert.kin$inputType <- "Kinase"
vert.kin <- vert.kin[vert.kin$runID == "positive_degExp_1_omega_2_mu_1e-05_beta_10",]
selectedVertices <- rbind(vert.prot,vert.tf,vert.kin)

centr.prot <- read.csv("PCSF_Protein_2024-12-12/centrality.csv")
centr.prot$inputType <- "Protein"
centr.prot$runID <- sub("_prot_","_",centr.prot$runID)
centr.tf <- read.csv("PCSF_Transcription_factor_2025-01-10/centrality.csv")
centr.tf$inputType <- "Transcription Factor"
centr.kin <- read.csv("PCSF_Kinase_2025-01-10/centrality.csv")
centr.kin$inputType <- "Kinase"
keepCols <- colnames(centr.prot)
centr.tf <- centr.tf[,keepCols]
centr.kin <- centr.kin[,keepCols]
all.centrality <- rbind(centr.tf, centr.prot, centr.kin)
all.centrality$runID2 <- sub(".*tive_degExp_1_","",all.centrality$runID)
all.centrality$Direction <- sub("\\_.*","",all.centrality$runID)
selectedCentrality <- all.centrality[all.centrality$runID == "positive_degExp_1_omega_2_mu_1e-05_beta_10",]

edge.prot <- read.csv("PCSF_Protein_2024-12-12/edges.csv")
edge.prot$inputType <- "Protein"
edge.prot$runID <- sub("_prot_","_",edge.prot$runID)
edge.prot$beta <- edge.prot$betaProt
edge.prot[,c("betaTF","betaKin","betaProt")] <- NULL
edge.tf <- read.csv("PCSF_Transcription_factor_2025-01-10/edges.csv")
edge.tf$inputType <- "Transcription Factor"
edge.kin <- read.csv("PCSF_Kinase_2025-01-10/edges.csv")
edge.kin$inputType <- "Kinase"
all.edges <- rbind(edge.prot, edge.tf, edge.kin)
selectedEdges <- all.edges[all.edges$runID == "positive_degExp_1_omega_2_mu_1e-05_beta_10",]
write.csv(selectedCentrality, "centrality_optimalRun.csv", row.names = FALSE)
write.csv(selectedEdges, "edges_optimalRun.csv", row.names = FALSE)
write.csv(selectedVertices, "vertices_optimalRun.csv", row.names = FALSE)
directions <- c("positive")
for (j in directions) {
  temp.centr <- selectedCentrality[grepl(j,selectedCentrality$runID),]
  temp.vert <- selectedVertices[grepl(j, selectedVertices$runID),]
  temp.edges <- selectedEdges[grepl(j, selectedEdges$runID),]
  write.csv(selectedCentrality, paste0("centrality_optimalRun_",j,"_", Sys.Date(),".csv"), row.names = FALSE)
  write.csv(selectedEdges, paste0("edges_optimalRun_",j,"_", Sys.Date(),".csv"), row.names = FALSE)
  write.csv(selectedVertices, paste0("vertices_optimalRun_",j,"_", Sys.Date(),".csv"), row.names = FALSE)
}

### save as xlsx
library(openxlsx)
Sheet <- c("Analysis", "Vertices", "Edges")

all.vertices <- read.csv("vertices.csv")
all.vertices$betaKin <- NULL
all.vertices$betaTF <- NULL
colnames(all.vertices)[3] <- "outputPrize"
colnames(all.vertices)[9] <- "beta"

all.edges$betaKin <- NULL
all.edges$betaTF <- NULL
colnames(all.edges)[7] <- "beta"

readme <- data.frame(Sheet)
readme$Description <- c("Network analysis including node degree, closeness, betweenness, eigen centrality, hub score, and authority score",
                        "Information about nodes including frequency in Prize Collecting Steiner Forest randomization, output prize, input prize, and node type",
                        "Information about edges including weight")
full.list <- list(readme, all.centrality, all.vertices, all.edges)
names(full.list) <- c("README", Sheet)
openxlsx::write.xlsx(full.list, file=paste0("SupplementaryTable5_NetworkOptimization_",Sys.Date(),".xlsx"), rowNames=FALSE)


selectedRun <- runOrder[1] # "omega_3_mu_0.001_beta_prot_10"
selectedRun <- "omega_2_mu_1e-05_beta_prot_10"
runIDs <- paste0(c("positive_degExp_1_", "negative_degExp_1_"), selectedRun)
selectedCentrality <- all.centrality[all.centrality$runID %in% runIDs,]
selectedEdges <- all.edges[all.edges$runID %in% runIDs,]
selectedVertices <- all.vertices[all.vertices$runID %in% runIDs,]
#selectedq24 <- temp.8q24[temp.8q24$runID %in% runIDs,] # 0 rows
sel.list <- list(readme, selectedCentrality, selectedVertices, selectedEdges)
names(sel.list) <- c("README", Sheet)
openxlsx::write.xlsx(sel.list, file=paste0("SupplementaryTable6_OptimizedNetworks_",Sys.Date(),".xlsx"), rowNames=FALSE)

### save full networks as PDF
library(igraph)
pos.graph <- igraph::read_graph("positive/degExp_1/mu_1e-05/omega_2/beta_prot_10.gml", format="gml")
plot(pos.graph)
temp.fname <- "positive/degExp_1/mu_1e-05/omega_2/beta_prot_10"
webshot(paste0(temp.fname,".html"), paste0(temp.fname,".pdf"))
temp.fname <- "negative/degExp_1/mu_1e-05/omega_2/beta_prot_10"
webshot(paste0(temp.fname,".html"), paste0(temp.fname,".pdf"))

### add Spearman correlation estimates to vertex information
# maybe networks are too big - reduce to top 10-15 nodes by eigen centrality
setwd("PCSF_Protein_2024-12-12")
selectedCentrality <- read.csv("centrality_optimalRun.csv")
selectedEdges <- read.csv("edges_optimalRun.csv")
selectedVertices <- read.csv("vertices_optimalRun.csv")
selectedVertices$Spearman.est <- NA
global.result <- na.omit(read.csv(synapser::synGet("syn63435101")$path))
for (i in 1:nrow(selectedVertices)) {
  if (selectedVertices$name[i] %in% global.result$Gene) {
    selectedVertices$Spearman.est[i] <- global.result[global.result$Gene == selectedVertices$name[i],]$Spearman.est
    selectedVertices$Spearman.q[i] <- global.result[global.result$Gene == selectedVertices$name[i],]$Spearman.q
  }
}
selectedVertices$abs.Spearman.est <- abs(as.numeric(selectedVertices$Spearman.est))
selectedVertices$sig.corr <- FALSE
selectedVertices[selectedVertices$Spearman.q <= 0.05,]$sig.corr <- TRUE
selectedVertices$Absolute.Spearman.rho <- NULL
write.csv(selectedVertices, paste0("vertices_optimalRun_withSpearmanRho_",Sys.Date(),".csv"), row.names = FALSE)

### plot in Cytoscape
#install.packages("BiocManager")
BiocManager::install("RCy3")
#setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/")
setwd("~/OneDrive - PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Network/unified/")
selectedCentrality <- read.csv("centrality_optimalRun_positive_2025-01-12.csv")
selectedEdges <- read.csv("edges_optimalRun_positive_2025-01-12.csv")
selectedVertices <- read.csv("vertices_optimalRun_withSpearmanRho_2025-01-21.csv")
directions <- "positive"
library(RCy3)
nTop.list <- c(50, 25, 20,15,10,5)
#nTop.list <- 5
library(plyr); library(dplyr)
for (nTop in nTop.list) {
  for (j in directions) {
    topGenes <- selectedCentrality[grepl(j, selectedCentrality$runID),] %>% 
      slice_max(eigen_centrality,n=nTop) # or maybe 2% but different # of nodes (862 positive, 1219 negative)
    # topVert <- selectedVertices[selectedVertices$Direction ==j & 
    #                               selectedVertices$name %in% topGenes$name,]
    tempEdges <- selectedEdges[selectedEdges$Direction==j & 
                                 (selectedEdges$to %in% topGenes$name | 
                                    selectedEdges$from %in% topGenes$name),]
    tempGenes <- unique(c(tempEdges$to, tempEdges$from)) # 74
    tempVert <- selectedVertices[selectedVertices$Direction ==j & 
                                   selectedVertices$name %in% tempGenes,]
    tempVert2 <- plyr::ddply(tempVert, .(name, Direction, mu, omega, beta, degExp, Spearman.est, abs.Spearman.est, Spearman.q, sig.corr), 
                            summarize,
                            frequency = mean(frequency),
                            prize = mean(prize),
                            type = ifelse(length(unique(type))==1, unique(type), "Terminal"),
                            inputPrize = mean(inputPrize),
                            nodeType = paste0(unique(nodeType), collapse=" and "),
                            inputType = paste0(unique(inputType), collapse=" and "))
    cat(j,"network has", nrow(tempVert2), "nodes and", nrow(tempEdges), "edges from top", nTop, "central nodes")
    topGraph <- igraph::graph_from_data_frame(tempEdges, directed=FALSE, vertices=tempVert2) # error: duplicate vertex names; e.g., DNMT1 is steiner in one but terminal in another
    plot(topGraph)
    tempTitle <- paste0(j, "_", nTop,"_topCentralNodes_wTFsKinases_",Sys.Date())
    RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
  } 
}
# set node colors using: RColorBrewer::brewer.pal(7, "Set2") ("#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494")
# where kinase: "#66C2A5", protein: "#FC8D62", TF: "#FFD92F"

### create venn diagram of nodes in each network
venn.list <- list()
inputTypes <- sort(unique(selectedVertices$inputType))
my.color.pal <- c("#66C2A5", "#FC8D62", "#FFD92F")
for (i in inputTypes) {
  venn.list[[i]] <- selectedVertices[grepl(i,selectedVertices$inputType),]$name
}
ggvenn::ggvenn(venn.list, fill_color = RColorBrewer::brewer.pal(7, "Set2"), show_percentage = FALSE, text_size=6)
ggplot2::ggsave("PCSF_ProteinTFKinase_unified_vennDiagram.pdf", width=4, height=4)

### also try expanding network based on PPI to connect vertices
library(PCSF)
data("STRING") # edges (1991761)
#ppi <- construct_interactome(STRING)

# keep string db edges if both "to" and "from" are in the vertices
stringEdges <- STRING[STRING$to %in% selectedVertices$name &
                        STRING$from %in% selectedVertices$name,] # 35910

# remake networks
library(RCy3)
nTop.list <- c(50, 25, 20,15,10,5)
#nTop.list <- 5
library(plyr); library(dplyr)
directions <- "positive"
for (nTop in nTop.list) {
  for (j in directions) {
    topGenes <- selectedCentrality[grepl(j, selectedCentrality$runID),] %>% 
      slice_max(eigen_centrality,n=nTop) # or maybe 2% but different # of nodes (862 positive, 1219 negative)
    # topVert <- selectedVertices[selectedVertices$Direction ==j & 
    #                               selectedVertices$name %in% topGenes$name,]
    tempEdges <- stringEdges[stringEdges$to %in% topGenes$name | 
                                    stringEdges$from %in% topGenes$name,]
    tempGenes <- unique(c(tempEdges$to, tempEdges$from)) # 74
    tempVert <- selectedVertices[selectedVertices$Direction ==j & 
                                   selectedVertices$name %in% tempGenes,]
    tempVert2 <- plyr::ddply(tempVert, .(name, Direction, mu, omega, beta, degExp, Spearman.est, abs.Spearman.est, Spearman.q, sig.corr), 
                             summarize,
                             frequency = mean(frequency),
                             prize = mean(prize),
                             type = ifelse(length(unique(type))==1, unique(type), "Terminal"),
                             inputPrize = mean(inputPrize),
                             nodeType = paste0(unique(nodeType), collapse=" and "),
                             inputType = paste0(unique(inputType), collapse=" and "))
    cat(j,"network has", nrow(tempVert2), "nodes and", nrow(tempEdges), "edges from top", nTop, "central nodes")
    topGraph <- igraph::graph_from_data_frame(tempEdges, directed=FALSE, vertices=tempVert2) # error: duplicate vertex names; e.g., DNMT1 is steiner in one but terminal in another
    plot(topGraph)
    tempTitle <- paste0(j, "_", nTop,"_topCentralNodes_wTFsKinases_expandedBySTRING_",Sys.Date())
    RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
    my.enr <- try(PCSF::enrichment_analysis(topGraph), 
                  silent = TRUE)
    if (!inherits(my.enr, "try-error")) {
     write.csv(my.enr, paste0("enrichr_",tempTitle,".csv"), row.names = FALSE)
    }
  } 
}

# Error in PCSF::enrichment_analysis(topGraph) : 
# The subnetwork must be a "PCSF" object derived from an "igraph" class.

# maybe just plot enrichr results from individual omics instead
setwd("PCSF_Protein_2024-12-12")
my.enr <- read.csv("enrichr.csv")
my.enr <- my.enr[my.enr$mu == "1e-05" & my.enr$betaProt == 10 &
                   my.enr$omega == 2,]
my.enr$N_included <- sub("/.*","",my.enr$Overlap)
my.enr$N_total <- sub(".*/","",my.enr$Overlap)
my.enr$perc_included <- 100*as.numeric(my.enr$N_included)/as.numeric(my.enr$N_total)
my.enr2 <- my.enr[my.enr$N_included >= 6 & my.enr$perc_included >= 10,] # 341

pos.enr <- my.enr2[my.enr2$Direction == "positive",] # 148
sig.pos.enr <- pos.enr[pos.enr$Adjusted.P.value <= 0.05,] # 148
sig.pos.enr.info <- plyr::ddply(sig.pos.enr, .(Term), summarize,
                                nClusters = length(unique(Cluster)),
                                meanOR = mean(Odds.Ratio)) # only 2 / 146 are in more than 1 cluster (and those are in 2)
clusters <- unique(sig.pos.enr$Cluster)
top.pos.enr <- data.frame()
for (i in clusters) {
  temp.pos.enr <- sig.pos.enr[sig.pos.enr$Cluster == i,] %>% slice_max(Odds.Ratio, n = 5) 
  top.pos.enr <- rbind(top.pos.enr, temp.pos.enr)
}

# make a circular bar plot for each omics' positive network
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/"
setwd(base.path)
path.map <- list("Transcription Factor" = "PCSF_Transcription_factor_2025-01-10",
                 "Protein" = "PCSF_Protein_2024-12-12",
                 "Kinase" = "PCSF_Kinase_2025-01-10")
path.map <- list("Protein" = "PCSF_Protein_2024-12-12")
n.top <- 3
enr.plots <- NULL
#source("~/OneDrive - PNNL/Documents/circBar.R")
directions <- c("positive","negative")
for (i in names(path.map)) {
  # load enrichr results
  setwd(path.map[[i]])
  my.enr <- read.csv("enrichr.csv")
  
  # filter for parameters, gene set coverage, significance
  my.enr <- my.enr[my.enr$mu == "1e-05" & my.enr$betaProt == 10 &
                     my.enr$omega == 2,]
  my.enr$N_included <- sub("/.*","",my.enr$Overlap)
  my.enr$N_total <- sub(".*/","",my.enr$Overlap)
  my.enr$perc_included <- 100*as.numeric(my.enr$N_included)/as.numeric(my.enr$N_total)
  my.enr2 <- my.enr[my.enr$N_included >= 6 & my.enr$perc_included >= 10,] # 341
  for (k in directions) {
    pos.enr <- my.enr2[my.enr2$Direction == k,] # 148
    sig.pos.enr <- pos.enr[pos.enr$Adjusted.P.value <= 0.05,] # 148
    
    # get top gene sets for each cluster
    clusters <- unique(sig.pos.enr$Cluster)
    top.pos.enr <- data.frame()
    for (j in clusters) {
      temp.pos.enr <- sig.pos.enr[sig.pos.enr$Cluster == j,] %>% slice_max(Odds.Ratio, n = 3) 
      top.pos.enr <- rbind(top.pos.enr, temp.pos.enr)
    }
    top.pos.enr$minusLogFDR <- -log(top.pos.enr$Adjusted.P.value, base=10)
    top.pos.enr$Odds.Ratio <- as.numeric(top.pos.enr$Odds.Ratio)
    top.pos.enr$Cluster2 <- paste("Cluster",top.pos.enr$Cluster)
    top.pos.enr$alpha <- 0.5 # 47 terms
    top.pos.enr$`Odds Ratio` <- top.pos.enr$Odds.Ratio
    top.pos.enr$`-Log(FDR)` <- top.pos.enr$minusLogFDR
    
    # generate and save plot
    setwd(base.path)
    # temp.enr.plot <- circBar(top.pos.enr, x="Term", y="Odds.Ratio", 
    #                          fill="Cluster2", alpha="alpha", ymin = 0,
    #                          ymax= 1200,
    #                          fillVals=grDevices::colorRampPalette(
    #                            RColorBrewer::brewer.pal(12, "Set3"))(length(clusters)),
    #                          title = i, fname=paste0(i,"_",k,"_enrichr_circBarPlot.pdf"))
    ggplot2::ggplot(top.pos.enr, aes(x=Term, y=`Odds Ratio`, fill=as.factor(Cluster), alpha=0.5)) +
      geom_col() + scale_fill_manual(values = c(grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(8, "Set2"))(length(clusters)))) +
      theme_classic(base_size=12)
    ggplot2::ggplot(top.pos.enr, aes(x=Term, y=`Odds Ratio`, fill=as.factor(Cluster), alpha=0.5)) +
      geom_col() + facet_wrap(~Cluster) +
      theme_minimal(base_size=12)
    absMax <- max(abs(top.pos.enr$Odds.Ratio))
    ggplot2::ggplot(top.pos.enr, aes(x=Cluster, y=Term, fill=`Odds Ratio`, 
                                     size=`-Log(FDR)`)) + geom_point() + 
      theme_minimal(base_size=12) +
      scale_color_gradient2(low="red",high="blue", mid="grey", limits=c(-absMax, absMax))
    
    top.top.pos.enr <- sig.pos.enr %>% slice_max(Odds.Ratio, n = 10) # 11 terms
    absMax <- max(abs(top.top.pos.enr$Odds.Ratio))
    minVal <- floor(min(top.top.pos.enr$Odds.Ratio))
    top.top.pos.enr$`-Log(FDR)` <- -log(top.top.pos.enr$Adjusted.P.value, base=10)
    top.top.pos.enr$`Odds Ratio` <- top.top.pos.enr$Odds.Ratio
    termOrder <- top.top.pos.enr[order(top.top.pos.enr$Combined.Score),]$Term
    ggplot2::ggplot(top.top.pos.enr, aes(x=as.factor(Cluster), y=Term, color=`Odds Ratio`, 
                                         size=`-Log(FDR)`)) + geom_point() + 
      theme_minimal(base_size=12) + ggplot2::scale_y_discrete(limits = termOrder) +
      scale_color_gradient2(high="red",low="grey", limits=c(minVal, absMax)) + xlab("Cluster")
    ggplot2::ggsave(paste0(i,"_",k,"_enrichr_dotPlot.pdf"), width=9, height=3)
  }
}

# don't have enrichr results for TF, kinase; need to check if it just wasn't run or if the network was too small
inputs <- list("Transcription_factor" = rna.allSamples.result,
               "Kinase" = mit.kin)
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/"
for (i in 1:length(inputs)) {
  final.inputs <- list()
  final.inputs[[names(inputs)[i]]] <- as.list(inputs[[i]]$Spearman.est)
  names(final.inputs[[names(inputs)[i]]]) <- inputs[[i]][,1]
  pcsf.result <- try(computeProteinNetwork_BG(final.inputs, 
                                              betas = 10, 
                                              mu = 1e-05, w = 2,
                                              deg.exp = 1,
                                              fname = names(inputs)[i]), silent=TRUE)
  my.enr <- try(PCSF::enrichment_analysis(pcsf.result$graph), 
                silent = TRUE) 
  temp.runID <- paste0(j, "_degExp_", deg.exp, "_omega_",omega,"_mu_",mu,"_beta_", beta)
  if (!inherits(my.enr, "try-error")) {
    temp.enrichr <- my.enr$enrichment
    temp.enrichr$Direction <- "positive"
    temp.enrichr[,"mu"] <- 1e-05
    temp.enrichr[,"beta"] <- 10
    temp.enrichr[,"omega"] <- 2
    temp.enrichr[,"degExp"] <- 1
    temp.enrichr$runID <- temp.runID 
  }
  temp.path <- paste0("PCSF_",names(inputs)[i],"_2025-01-10")
  setwd(temp.path)
  write.csv(temp.enrichr, "enrichr.csv", row.names = FALSE)
  setwd(base.path)
}

# Performing enrichment analysis...
# 
# Enrichment is being performed by EnrichR (http://amp.pharm.mssm.edu/Enrichr) API ...
# Error in `[<-.data.frame`(`*tmp*`, i, , value = c("", "<title>Enrichr</title>" : 
#                                                     replacement has 2 rows, data has 1
#                                                   In addition: Warning message:
#                                                     In cluster_edge_betweenness(subnet) :
#                                                     At vendor/cigraph/src/community/edge_betweenness.c:504 : Membership vector will be selected based on the highest modularity score.

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/"
setwd(base.path)
path.map <- list("Transcription Factor" = "PCSF_Transcription_factor_2025-01-10",
                 "Protein" = "PCSF_Protein_2024-12-12",
                 "Kinase" = "PCSF_Kinase_2025-01-10")
inputs <- list("Transcription_factor" = rna.allSamples.result, "Protein" = global.result, "Kinase" = mit.kin)
for (i in 1:length(path.map)) {
  setwd(path.map[[i]])
  temp.vert <- read.csv("vertices.csv")
  if (names(path.map)[i]=="Protein") {
    temp.vert <- temp.vert[temp.vert$mu==1e-05 & temp.vert$betaProt==10 &
                             temp.vert$omega==2,]
  } else {
    temp.vert <- temp.vert[temp.vert$mu==1e-05 & temp.vert$beta==10 &
                             temp.vert$omega==2,] 
  }
  nInputs <- length(na.omit(unique(inputs[[i]][,1])))
  nIncluded <- length(na.omit(unique(temp.vert[temp.vert$name %in% inputs[[i]][,1],]$name)))
  cat(names(path.map)[i], nIncluded, "/", nInputs, "included\n")
  setwd(base.path)
}
