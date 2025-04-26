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
#source("https://github.com/PNNL-CompBio/MPNST_Chr8/blob/main/panSEA_helper_20240913.R")
source("~/OneDrive - PNNL/Documents/GitHub/Chr8/panSEA_helper_20240913.R")
#source("https://github.com/PNNL-CompBio/MPNST_Chr8/blob/main/networkAnalysis/proteinNetworks_BG.R")
source("~/OneDrive - PNNL/Documents/GitHub/Chr8/networkAnalysis/proteinNetworks_BG.R")
#source("https://github.com/PNNL-CompBio/MPNST_Chr8/blob/main/networkAnalysis/PCSF_rand_BG.R")
source("~/OneDrive - PNNL/Documents/GitHub/Chr8/networkAnalysis/PCSF_rand_BG.R")

#BiocManager::install("topGO") # need for PCSF
#remotes::install_github("sgosline/amlresistancenetworks") # need for PCSF; same repo on PNNL-CompBio


if (file.exists("gmt_human_positional.rds")) {
  gmt.pos <- readRDS("gmt_human_positional.rds")
} else {
  pos.info <- msigdbr::msigdbr(collection="C1")
  gmt.pos <- DMEA::as_gmt(as.data.frame(pos.info), element.names="gene_symbol", set.names = "gs_name")
  saveRDS(gmt.pos,"gmt_human_positional.rds")
}

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
setwd("Chr8_quant_20250409")


#### prep inputs ####
# try using top 50 sig prot, phospho, WES - directional this time
# load feature weights
synapser::synLogin()
global.result <- na.omit(read.csv(synapser::synGet("syn66224803")$path))
tf.result <- na.omit(read.csv(synapser::synGet("syn66226952")$path))
kin.result <- read.csv(synapser::synGet("syn66279699")$path)

# get total N features
n.prot <- nrow(global.result) # 9013
n.tf <- nrow(tf.result) # 408
n.kin <- nrow(kin.result) # 184

# filter for significance
global.result <- global.result[global.result$Spearman.q <= 0.05, ] # 208 / 9013 (2.31%); 98.02% of 101 positive, 100% of 107 negative are in the interactome
tf.result <- tf.result[tf.result$p_value <= 0.05 & tf.result$FDR_q_value <= 0.25, ] # 206 / 408 (50.49%); 96.1% of 205 positive, 100% of 1 negative are in the interactome
kin.result <- kin.result[kin.result$p_value <= 0.05 & kin.result$FDR_q_value <= 0.25, ] # 2 / 184 (1.09%); 100% of 2 negative are in the interactome

tf.result$Gene <- sub("_.*","",tf.result$Feature_set)

setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/Chr8_quant_20250409")

#### use all terminals since PCSF only uses kinases ####
library(PCSF)
data("STRINGv12") # 1477610
pos.terms <- c(global.result[global.result$Spearman.est>0,]$Gene, 
                      tf.result[tf.result$NES>0,]$Gene) # 305; no pos kinases
# are any terms from both protein or TF or kinase?
any(duplicated(pos.terms)) # yes
dup.pos <- pos.terms[duplicated(pos.terms)] # ZNF22 up in both protein and TF

neg.terms <- c(global.result[global.result$Spearman.est<0,]$Gene, 
                      tf.result[tf.result$NES<0,]$Gene, kin.result[kin.result$NES<0,]$Feature_set) # 110
# are any terms from both protein or TF or kinase?
any(duplicated(neg.terms)) # no

pos.edges <- STRINGv12[STRINGv12$from %in% pos.terms | STRINGv12$to %in% pos.terms,] # 58206
pos.genes <- unique(c(pos.edges$from, pos.edges$to)) # 9133
pos.vert <- data.frame(pos.genes, type = "Inferred", Omics = NA)
pos.vert[pos.vert$pos.genes %in% pos.terms,]$type <- "Terminal"
nrow(pos.vert[pos.vert$type == "Inferred",]) # 8834
pos.vert[pos.vert$pos.genes %in% global.result[global.result$Spearman.est>0,]$Gene,]$Omics <- "Protein"
pos.vert[pos.vert$pos.genes %in% tf.result[tf.result$NES>0,]$Gene,]$Omics <- "TF"
pos.vert[pos.vert$pos.genes == dup.pos,]$Omics <- "Protein & TF"

neg.edges <- STRINGv12[STRINGv12$from %in% neg.terms | STRINGv12$to %in% neg.terms,] # 26390
neg.genes <- unique(c(neg.edges$from, neg.edges$to)) # 6374
neg.vert <- data.frame(neg.genes, type = "Inferred", Omics=NA)
neg.vert[neg.vert$neg.genes %in% neg.terms,]$type <- "Terminal"
nrow(neg.vert[neg.vert$type == "Inferred",]) # 6265
neg.vert[neg.vert$neg.genes %in% global.result[global.result$Spearman.est<0,]$Gene,]$Omics <- "Protein"
neg.vert[neg.vert$neg.genes %in% tf.result[tf.result$NES<0,]$Gene,]$Omics <- "TF"
neg.vert[neg.vert$neg.genes %in% kin.result[kin.result$NES<0,]$Feature_set,]$Omics <- "Kinase"

# narrow it down by filtering for nodes with most interactions
posN <- plyr::ddply(pos.edges, .(from), summarize,
                    N = n()) # only need 'from' since edges are symmetrical
colnames(pos.vert)[1] <- "from"
posN <- merge(posN, pos.vert, by="from")
hist(posN$N)
hist(posN[posN$N < 50,]$N)
hist(posN[posN$N < 10,]$N)
quantile(posN$N)
# 0%  25%  50%  75% 100% 
# 1    1    2    4  897 
nrow(posN[posN$N > 4,]) # 2057
pos500 <- posN %>% slice_max(N, n=500) # 525
min(pos500$N) # 13
pos1000 <- posN %>% slice_max(N, n=1000) # 1024
min(pos1000$N) # 8
pos1000edges <- pos.edges[pos.edges$from %in% pos1000$from,] # fine because symmetrical
pos1000vert <- pos.vert[pos.vert$from %in% c(pos1000edges$from, pos1000edges$to),] # back to 9133

negN <- plyr::ddply(neg.edges, .(from), summarize,
                    N = n()) # only need 'from' since edges are symmetrical
colnames(neg.vert)[1] <- "from"
negN <- merge(negN, neg.vert, by="from")
hist(negN$N)
hist(negN[negN$N < 50,]$N)
hist(negN[negN$N < 10,]$N)
quantile(negN$N)
# 0%  25%  50%  75% 100% 
# 1    1    1    2  592 
nrow(negN[negN$N > 2,]) # 1544
neg500 <- negN %>% slice_max(N, n=500) # 605
min(neg500$N) # 5
neg1000 <- negN %>% slice_max(N, n=1000) # 1544
min(neg1000$N) # 3

# make graphs
topGraph <- igraph::graph_from_data_frame(pos.edges, directed=FALSE, vertices=posN) 
plot(topGraph)
tempTitle <- paste0("positive", "_", nrow(pos1000),"_topConnectedNodes_manual_",Sys.Date())
RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)

topGraph <- igraph::graph_from_data_frame(neg.edges, directed=FALSE, vertices=negN) 
plot(topGraph)
tempTitle <- paste0("negative", "_", nrow(neg1000),"_topConnectedNodes_manual_",Sys.Date())
RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)

# try filtering for nodes which connect terminals
pos.con <- pos.edges[pos.edges$from %in% pos.terms & pos.edges$to %in% pos.terms,] # 752
pos.con.vert <- posN[posN$from %in% c(pos.con$from, pos.con$to),] # 206 - all terminal

pos.con.vert <- posN[posN$N > 1,] # 5630
pos.con <- pos.edges[pos.edges$from %in% pos.con.vert$from & pos.edges$to %in% pos.con.vert$from,] # 51200
topGraph <- igraph::graph_from_data_frame(pos.con, directed=FALSE, vertices=pos.con.vert) 
#plot(topGraph)
tempTitle <- paste0("positive", "_", nrow(pos.con.vert),"_topConnectedNodes_manual_1PlusConnection_",Sys.Date())
RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)

pos.con.centrality <- data.frame(name = V(topGraph)$name,
                                   degree = igraph::degree(topGraph, mode="all"),
                                   closeness = igraph::closeness(topGraph, mode="all"),
                                   betweenness = igraph::betweenness(topGraph, directed = FALSE),
                                   eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                                   hub_score = igraph::hub_score(topGraph)$vector,
                                   authority_score = igraph::authority_score(topGraph)$vector)

#### multiomics PCSF ####
#temp.path <- paste0("PCSF_TFProteinKinase_", Sys.Date(), "_kinasesRenamed")

deg.exp.vals <- 1
beta.vals <- c(1,5,10)
beta.combos <- tidyr::expand_grid("Protein" = c(1,5,10), "TF" = c(1,5,10), "Kinase" = c(1,5,10))
mu.vals <- c(1E-7, 1E-5, 1E-3) # next try adding 1, 1E-7, 1E-3, 1E-5
omega.vals <- c(1,2,3) # next try adding 1, 2
inputs <- list("Protein" = global.result[,c("Gene","Spearman.est")], 
               "TF" = tf.result[,c("Gene","NES")], 
               "Kinase" = kin.result[,c("Feature_set","NES")])
directions <- c("negative","positive")

setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/Chr8_quant_20250409")
temp.path <- paste0("PCSF_", Sys.Date())
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
  for (m in 1:length(inputs)) { # second column of each input must be rank value
    if (grepl("positive", j)) {
      top.phospho <- inputs[[m]][inputs[[m]][,2] > 0,] 
    } else if (grepl("negative", j)) {
      top.phospho <- inputs[[m]][inputs[[m]][,2] < 0,]
    } else {
      top.phospho <- inputs[[m]]  
    } 
    
    # compile inputs in named lists
    if (nrow(top.phospho) > 0) {
      final.inputs[[names(inputs)[m]]] <- as.list(top.phospho[,2])
      names(final.inputs[[names(inputs)[m]]]) <- top.phospho[,1] # assuming first column contains feature names
    }
  }
  temp.beta.combos <- dplyr::distinct(beta.combos[,names(final.inputs)])
  
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
        
        for (beta in 1:nrow(temp.beta.combos)) {
          setwd(file.path(base.path2, j, paste0("degExp_",deg.exp), paste0("mu_",mu), paste0("omega_",omega)))
          
          # run PCSF
          temp.fname <- paste0("beta_", paste0(unlist(temp.beta.combos[beta,]), collapse="-"))
          temp.runID <- paste0(j, "_degExp_", deg.exp, "_omega_",omega,"_mu_",mu,"_",temp.fname)
          if (!(temp.runID %in% runs.so.far)) {
            cat("Running PCSF\n")
            pcsf.result <- try(computeProteinNetwork_BG(final.inputs, 
                                                        betas = unlist(temp.beta.combos[beta,]), 
                                                        mu = mu, w = omega,
                                                        deg.exp = deg.exp,
                                                        fname = temp.fname))
            runs.so.far <- c(runs.so.far, temp.runID)
            all.dfs <- processPCSF(pcsf.result, all.dfs, temp.runID, 
                                   paste0(unlist(temp.beta.combos[beta,]), collapse="-"), 
                                   mu, omega, base.path2, deg.exp, j, temp.fname)
          }
        }
      } 
    }
  }
}

library(ggplot2)
#### plot all results ####
library(patchwork)

# sort by rank and then enforce that runID order
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/Chr8_quant_20250409")
all.vertices <- read.csv("PCSF_Protein_2025-04-16/vertices.csv")

param.info <- as.data.frame(tidyr::expand_grid(omega=omega.vals, mu=mu.vals,  
                                               beta=beta.vals))

all.centrality <- read.csv("PCSF_Protein_2025-04-16/centrality.csv")
all.centrality$runID2 <- sub(".*tive_degExp_1_","",all.centrality$runID)
all.centrality$Direction <- sub("\\_.*","",all.centrality$runID)

# enr.8q24 <- read.csv("PCSF_Protein_2025-04-16/chr8q24_enrichment.csv") # 0 rows
# enr.8q24$runID2 <- sub(".*tive_degExp_1_","",enr.8q24$runID)
enr.8q24 <- data.frame()

### node degrees (steiner vs. terminal) for each parameter setting
degree.info <- dplyr::distinct(merge(all.vertices, all.centrality))
degree.summary <- plyr::ddply(degree.info, .(type, omega, mu, beta, Direction), summarize,
                              degree_sum = sum(degree, na.rm=TRUE))

### % of terminal nodes included in output for each parameter setting
perc.term <- plyr::ddply(all.vertices[all.vertices$type=="Terminal",], 
                         .(omega, mu, beta, Direction), summarize,
                         N_terminal = length(unique(name)))
#perc.term$N_inputs <- 441 # number of negatively correlated proteins with adj. p <= 0.05 (no negatively enriched TFs or kinases)
perc.term$N_inputs <- nrow(global.result[global.result$Spearman.est < 0,])
perc.term[perc.term$Direction == "positive",]$N_inputs <- nrow(global.result[global.result$Spearman.est > 0,])
perc.term$percent_terminal <- perc.term$N_terminal * 100 / perc.term$N_inputs
perc.term$runID2 <- paste0("omega_", perc.term$omega, "_mu_", perc.term$mu,
                           "_beta_",perc.term$beta)

### chr8q24 enrichment plot - expect negative NES because highest prize would be at bottom of list even though unweighted calculation
# is the score still affected by input order in unweighted GSEA? yes, well mostly affected by prize (i.e., rank metric) then input order
# this paper uses betweenness centrality to rank nodes: https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2021.577623/full
if (nrow(enr.8q24) > 0) {
  temp.8q24 <- enr.8q24[enr.8q24$Feature_set == "chr8q24" & enr.8q24$minPerSet == 6,]
  if (nrow(temp.8q24) > 0) {
    if (any(temp.8q24$FDR_q_value == 0)) {
      temp.8q24[temp.8q24$FDR_q_value == 0,]$FDR_q_value <- 1E-4
    }
    temp.8q24$chr8q24 <- rank(-log10(temp.8q24$FDR_q_value) * temp.8q24$NES)
    temp.8q24$runID2 <- sub(".*tive_degExp_1_","",temp.8q24$runID)
  } 
}

### mean rank (min NES for chr8q24, min diff in node degrees, max percent terminal) - could maybe specify SUMO1, other general nodes to exclude
degree.summary2 <- reshape2::dcast(degree.summary, omega+mu+beta+Direction ~ type, value.var = "degree_sum")
degree.summary2$runID2 <- paste0("omega_", degree.summary2$omega, "_mu_", 
                                 degree.summary2$mu, "_beta_",degree.summary2$beta)
perc.term2 <- plyr::ddply(perc.term, .(runID2), summarize,
                          perc_terminal = mean(percent_terminal, na.rm = TRUE),
                          sd_perc_terminal = sd(percent_terminal, na.rm = TRUE))
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
} else {
  rank.summary2 <- plyr::ddply(rank.summary, .(runID2), summarize,
                               degree_diff = mean(rank_degDiff, na.rm = TRUE),
                               perc_terminal = mean(rank_percTerm, na.rm = TRUE),
                               rank = mean(c(rank_degDiff, rank_percTerm), na.rm = TRUE))
  runOrder <- rank.summary2[order(rank.summary2$rank),]$runID2
  rank.summary3 <- reshape2::melt(rank.summary2, id = "runID2", variable.name = "Parameter")
  library(ggplot2)
  rank.plot <- ggplot2::ggplot(rank.summary3, 
                               aes(x = runID2, y = value, group = Parameter, 
                                   color = Parameter)) +
    geom_line() + ylab("Rank") + theme_classic() +
    scale_color_discrete(labels=c("Difference\nin degree","% terminal\nnodes used","Overall"))+
    ggplot2::scale_x_discrete(limits = runOrder) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
}
rank.plot

### create plots now that runID overall ranks are known
# create line plot of node degrees vs. runID
#runOrder <- sort(unique(rank.summary$runID))
degree.summary$runID2 <- paste0("omega_", degree.summary$omega, "_mu_", 
                                 degree.summary$mu, "_beta_",degree.summary$beta)
degree.summary$DirectionType <- paste0(degree.summary$Direction, "-", degree.summary$type)
deg.plot <- ggplot2::ggplot(degree.summary, aes(x = runID2, y = degree_sum, 
                                                group = DirectionType, color = type)) +
  geom_line(aes(linetype=Direction)) + ylab("Sum of Node Degrees") + theme_classic() + 
  ggplot2::scale_x_discrete(limits = runOrder) + 
  scale_linetype_discrete(labels=c("Negative", "Positive")) +
  labs(color="Node Type")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
deg.plot

# percent terminal nodes used
perc.plot <- ggplot2::ggplot(perc.term, aes(x = runID2, y = percent_terminal, 
                                            group=Direction, linetype = Direction)) +
  geom_line() + ylab("% Terminal Nodes") + theme_classic() + 
  ggplot2::scale_x_discrete(limits = runOrder)+
  scale_linetype_discrete(labels=c("Negative", "Positive")) +
  scale_y_continuous(limits=c(0, 100)) +
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
                            "_beta_",param.info$beta)
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
#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
param.plot <- param.plot + plot_layout(guides = 'collect')
param.plot

### combine plots
if (nrow(temp.8q24) > 0) {
  combo.plot <- rank.plot / deg.plot / perc.plot / chr8.plot / param.plot
} else {
  combo.plot <- rank.plot / perc.plot / deg.plot / param.plot
}
combo.plot
ggplot2::ggsave(paste0("parameter_sweep_v3_",Sys.Date(),".pdf"), combo.plot, width = 6, height = 16)
ggplot2::ggsave(paste0("parameter_sweep_v3_shorter_",Sys.Date(),".pdf"), combo.plot, width = 4, height = 7)
saveRDS(combo.plot, paste0("parameter_sweep_v3_",Sys.Date(),".rds"))

#### save optimal run ####
selectedRun <- runOrder[1] # "omega_3_mu_0.001_beta_10"
selectedRun <- "omega_3_mu_0.001_beta_10"

vert.prot <- read.csv("PCSF_Protein_2025-04-16/vertices.csv")
selectedVertices <- vert.prot[grepl(selectedRun, vert.prot$runID),]

centr.prot <- read.csv("PCSF_Protein_2025-04-16/centrality.csv")
all.centrality <- centr.prot
all.centrality$runID2 <- sub(".*tive_degExp_1_","",all.centrality$runID)
all.centrality$Direction <- sub("\\_.*","",all.centrality$runID)
selectedCentrality <- all.centrality[grepl(selectedRun, all.centrality$runID),]

edge.prot <- read.csv("PCSF_Protein_2025-04-16/edges.csv")
all.edges <- edge.prot
selectedEdges <- all.edges[grepl(selectedRun, all.edges$runID),]
write.csv(selectedCentrality, "centrality_optimalRun.csv", row.names = FALSE)
write.csv(selectedEdges, "edges_optimalRun.csv", row.names = FALSE)
write.csv(selectedVertices, "vertices_optimalRun.csv", row.names = FALSE)
directions <- c("positive", "negative")
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

all.vertices <- read.csv("PCSF_Protein_2025-04-16/vertices.csv")
colnames(all.vertices)[3] <- "outputPrize"

readme <- data.frame(Sheet)
readme$Description <- c("Network analysis including node degree, closeness, betweenness, eigen centrality, hub score, and authority score",
                        "Information about nodes including frequency in Prize Collecting Steiner Forest randomization, output prize, input prize, and node type",
                        "Information about edges including weight")
full.list <- list(readme, all.centrality, all.vertices, all.edges)
names(full.list) <- c("README", Sheet)
openxlsx::write.xlsx(full.list, file=paste0("SupplementaryTable5_NetworkOptimization_",Sys.Date(),".xlsx"), rowNames=FALSE)


selectedRun <- runOrder[1] # "omega_3_mu_0.001_beta_10"
selectedRun <- "omega_3_mu_0.001_beta_10"
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
setwd("PCSF_Protein_2025-04-16")
pos.graph <- igraph::read_graph("positive/degExp_1/mu_0.001/omega_3/beta_10.gml", format="gml")
plot(pos.graph)
temp.fname <- "positive/degExp_1/mu_0.001/omega_3/beta_10"
webshot(paste0(temp.fname,".html"), paste0(temp.fname,".pdf"))
temp.fname <- "negative/degExp_1/mu_0.001/omega_3/beta_10"
webshot(paste0(temp.fname,".html"), paste0(temp.fname,".pdf"))

### add Spearman correlation estimates to vertex information
# maybe networks are too big - reduce to top 10-15 nodes by eigen centrality
selectedCentrality <- read.csv("centrality_optimalRun.csv")
selectedEdges <- read.csv("edges_optimalRun.csv")
selectedVertices <- read.csv("vertices_optimalRun.csv")
selectedVertices$Spearman.est <- NA
global.result <- na.omit(read.csv(synapser::synGet("syn66224803")$path))
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
#setwd("~/OneDrive - PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Network/unified/")
selectedCentrality <- read.csv("centrality_optimalRun.csv")
selectedEdges <- read.csv("edges_optimalRun.csv")
selectedVertices <- read.csv("vertices_optimalRun_withSpearmanRho_2025-04-17.csv")
directions <- c("positive", "negative")
library(RCy3)
nTop.list <- c(50, 25, 20,15,10,5, 100, "max") # 101 upregulated, 107 downregulated so "positive_107" network is really "positive_101"
#nTop.list <- 5
library(plyr); library(dplyr)
for (nTop in nTop.list) {
  for (j in directions) {
    if (nTop == "max") {
      topGenes <- selectedCentrality[grepl(j, selectedCentrality$runID),]
    } else {
      topGenes <- selectedCentrality[grepl(j, selectedCentrality$runID),] %>% 
        slice_max(eigen_centrality,n=nTop)
    }
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
                            nodeType = paste0(unique(nodeType), collapse=" and "))
    cat(j,"network has", nrow(tempVert2), "nodes and", nrow(tempEdges), "edges from top", nTop, "central nodes")
    topGraph <- igraph::graph_from_data_frame(tempEdges, directed=FALSE, vertices=tempVert2) # error: duplicate vertex names; e.g., DNMT1 is steiner in one but terminal in another
    plot(topGraph)
    tempTitle <- paste0(j, "_", nrow(topGenes),"_topCentralNodes",Sys.Date())
    RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
  } 
}
# set node colors using: RColorBrewer::brewer.pal(7, "Set2") ("#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494")
# where kinase: "#66C2A5", protein: "#FC8D62", TF: "#FFD92F"

### create venn diagram of nodes in each network
venn.list <- list()
#my.color.pal <- c("#66C2A5", "#FC8D62", "#FFD92F")
my.color.pal <- c("red", "blue")
for (i in directions) {
  if (i == "positive") {
    temp.name <- "Positive"
  } else { temp.name <- "Negative"}
  venn.list[[temp.name]] <- selectedVertices[selectedVertices$Direction==i,]$name
}
#ggvenn::ggvenn(venn.list, fill_color = RColorBrewer::brewer.pal(7, "Set2"), show_percentage = FALSE, text_size=6)
ggvenn::ggvenn(venn.list, fill_color = my.color.pal, show_percentage = FALSE, text_size=6, set_name_size=4.25)
ggplot2::ggsave("PCSF_Protein_Direction_unified_vennDiagram.pdf", width=4, height=4)
venn.list[["Positive"]][venn.list[["Positive"]] %in% venn.list[["Negative"]]]
# "ATG12"   "CHCHD4"  "CSK"     "DEPDC5"  "FBLIM1"  "FGFR2"   "GFER"    "ITGA1"   "ITGA4"   "MADCAM1" "MON1B"   "RAB7A"   "RRAGA"   "SDC2"    "SELL"    "VASP"  

# maybe just plot enrichr results from individual omics instead
setwd("PCSF_Protein_2025-04-16")
my.enr <- read.csv("enrichr.csv")
my.enr <- my.enr[my.enr$mu == 0.001 & my.enr$beta == 10 &
                   my.enr$omega == 3,] # none