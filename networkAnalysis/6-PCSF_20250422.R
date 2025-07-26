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
chr8q.genes <- unlist(gmt.pos$genesets[grepl("chr8q",names(gmt.pos$genesets))]) # 1010
saveRDS(chr8q.genes, "chr8qGenes.rds")
chr8q.genes <- readRDS("chr8qGenes.rds")

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
pos.vert$chr8q <- FALSE
pos.vert[pos.vert$pos.genes %in% chr8q.genes,]$chr8q <- TRUE
write.csv(pos.vert, "positive_vertices.csv", row.names=FALSE)
write.csv(pos.edges, "positive_edges.csv", row.names=FALSE)
topGraph <- igraph::graph_from_data_frame(pos.edges, directed=FALSE, vertices=pos.vert) 
pos.centrality <- data.frame(name = V(topGraph)$name,
                             degree = igraph::degree(topGraph, mode="all"),
                             closeness = igraph::closeness(topGraph, mode="all"),
                             betweenness = igraph::betweenness(topGraph, directed = FALSE),
                             eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                             hub_score = igraph::hub_score(topGraph)$vector,
                             authority_score = igraph::authority_score(topGraph)$vector)
write.csv(pos.centrality, "positive_centrality.csv", row.names=FALSE)
pos.vert <- read.csv("positive_vertices.csv")
pos.edges <- read.csv("positive_edges.csv")
pos.centrality <- read.csv("positive_centrality.csv")
chr8.pos.centr <- pos.centrality[pos.centrality$name %in% chr8q.genes,]

neg.edges <- STRINGv12[STRINGv12$from %in% neg.terms | STRINGv12$to %in% neg.terms,] # 26390
neg.genes <- unique(c(neg.edges$from, neg.edges$to)) # 6374
neg.vert <- data.frame(neg.genes, type = "Inferred", Omics=NA)
neg.vert[neg.vert$neg.genes %in% neg.terms,]$type <- "Terminal"
nrow(neg.vert[neg.vert$type == "Inferred",]) # 6265
neg.vert[neg.vert$neg.genes %in% global.result[global.result$Spearman.est<0,]$Gene,]$Omics <- "Protein"
neg.vert[neg.vert$neg.genes %in% tf.result[tf.result$NES<0,]$Gene,]$Omics <- "TF"
neg.vert[neg.vert$neg.genes %in% kin.result[kin.result$NES<0,]$Feature_set,]$Omics <- "Kinase"
neg.vert$chr8q <- FALSE
neg.vert[neg.vert$neg.genes %in% chr8q.genes,]$chr8q <- TRUE
write.csv(neg.vert, "negative_vertices.csv", row.names=FALSE)
write.csv(neg.edges, "negative_edges.csv", row.names=FALSE)
topGraph <- igraph::graph_from_data_frame(neg.edges, directed=FALSE, vertices=neg.vert) 
neg.centrality <- data.frame(name = V(topGraph)$name,
                                    degree = igraph::degree(topGraph, mode="all"),
                                    closeness = igraph::closeness(topGraph, mode="all"),
                                    betweenness = igraph::betweenness(topGraph, directed = FALSE),
                                    eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                                    hub_score = igraph::hub_score(topGraph)$vector,
                                    authority_score = igraph::authority_score(topGraph)$vector)
write.csv(neg.centrality, "negative_centrality.csv", row.names=FALSE)
neg.vert <- read.csv("negative_vertices.csv")
neg.edges <- read.csv("negative_edges.csv")
neg.centrality <- read.csv("negative_centrality.csv")

chr8.neg.centr <- neg.centrality[neg.centrality$name %in% chr8q.genes,]
chr8.neg.centr$Direction <- "negative"
chr8.pos.centr$Direction <- "positive"
chr8.centr <- rbind(chr8.neg.centr, chr8.pos.centr)
chr8.degree <- plyr::ddply(chr8.centr, .(name), summarize,
                                Directions = paste0(sort(unique(Direction)), collapse=", "),
                                meanDegree = mean(degree),
                                sumDegree = sum(degree),
                                meanCentrality = mean(eigen_centrality))
write.csv(chr8.degree, "chr8q_centrality.csv", row.names=FALSE)
chr8.degree$centrPerDegree <- chr8.degree$meanCentrality/chr8.degree$meanDegree
chr8.degree$centrPerSumDegree <- chr8.degree$meanCentrality/chr8.degree$sumDegree

pos.reg <- c(tf.result[tf.result$NES>0,]$Gene) # 205
pos.reg.edges <- STRINGv12[STRINGv12$from %in% pos.reg | STRINGv12$to %in% pos.reg,] # 38180
pos.reg.genes <- unique(c(pos.reg.edges$from, pos.reg.edges$to)) # 7109
pos.reg.vert <- data.frame(pos.reg.genes, type = "Inferred", Omics=NA)
pos.reg.vert[pos.reg.vert$pos.reg.genes %in% pos.reg,]$type <- "Terminal"
pos.reg.vert$chr8q <- FALSE
pos.reg.vert[pos.reg.vert$pos.reg.genes %in% chr8q.genes,]$chr8q <- TRUE
topGraph <- igraph::graph_from_data_frame(pos.reg.edges, directed=FALSE, vertices=pos.reg.vert) 
pos.reg.centrality <- data.frame(name = V(topGraph)$name,
                                 degree = igraph::degree(topGraph, mode="all"),
                                 closeness = igraph::closeness(topGraph, mode="all"),
                                 betweenness = igraph::betweenness(topGraph, directed = FALSE),
                                 eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                                 hub_score = igraph::hub_score(topGraph)$vector,
                                 authority_score = igraph::authority_score(topGraph)$vector)
pos.reg.centrality$chr8q <- FALSE
pos.reg.centrality[pos.reg.centrality$name %in% chr8q.genes,]$chr8q <- TRUE
write.csv(pos.reg.centrality, "positive_TFKin_centrality.csv", row.names=FALSE)
write.csv(pos.reg.vert, "positive_TFKin_vertices.csv", row.names=FALSE)
write.csv(pos.reg.edges, "positive_TFKin_edges.csv", row.names=FALSE)
pos.reg.centrality <- read.csv("positive_TFKin_centrality.csv")

neg.reg <- c(tf.result[tf.result$NES<0,]$Gene, kin.result[kin.result$NES<0,]$Feature_set) # 3
neg.reg.edges <- STRINGv12[STRINGv12$from %in% neg.reg | STRINGv12$to %in% neg.reg,] # 1172
neg.reg.genes <- unique(c(neg.reg.edges$from, neg.reg.edges$to)) # 574
neg.reg.vert <- data.frame(neg.reg.genes, type = "Inferred", Omics=NA)
neg.reg.vert[neg.reg.vert$neg.reg.genes %in% neg.reg,]$type <- "Terminal"
neg.reg.vert$chr8q <- FALSE
neg.reg.vert[neg.reg.vert$neg.reg.genes %in% chr8q.genes,]$chr8q <- TRUE
topGraph <- igraph::graph_from_data_frame(neg.reg.edges, directed=FALSE, vertices=neg.reg.vert) 
neg.reg.centrality <- data.frame(name = V(topGraph)$name,
                             degree = igraph::degree(topGraph, mode="all"),
                             closeness = igraph::closeness(topGraph, mode="all"),
                             betweenness = igraph::betweenness(topGraph, directed = FALSE),
                             eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                             hub_score = igraph::hub_score(topGraph)$vector,
                             authority_score = igraph::authority_score(topGraph)$vector)
neg.reg.centrality$chr8q <- FALSE
neg.reg.centrality[neg.reg.centrality$name %in% chr8q.genes,]$chr8q <- TRUE
write.csv(neg.reg.centrality, "negative_TFKin_centrality.csv", row.names=FALSE)
write.csv(neg.reg.vert, "negative_TFKin_vertices.csv", row.names=FALSE)
write.csv(neg.reg.edges, "negative_TFKin_edges.csv", row.names=FALSE)
neg.reg.centrality <- read.csv("negative_TFKin_centrality.csv")

neg.reg.centrality$Direction <- "negative"
pos.reg.centrality$Direction <- "positive"
chr8q.reg.centrality <- rbind(pos.reg.centrality[pos.reg.centrality$chr8q,], 
                              neg.reg.centrality[neg.reg.centrality$chr8q,]) # 173
write.csv(chr8q.reg.centrality, "chr8q_TFKin_centrality.csv", row.names=FALSE)
chr8q.reg.degree <- plyr::ddply(chr8q.reg.centrality, .(name), summarize,
                                Directions = paste0(sort(unique(Direction)), collapse=", "),
                                meanDegree = mean(degree),
                                sdDegree = sd(degree),
                                sumDegree = sum(degree),
                                meanCentrality = mean(eigen_centrality),
                                sdCentrality = sd(eigen_centrality))
write.csv(chr8q.reg.degree, "chr8q_TFKin_meanCentrality.csv", row.names=FALSE)
chr8q.reg.degree <- read.csv("chr8q_TFKin_meanCentrality.csv")

bar.df <- chr8q.reg.degree[chr8q.reg.degree$Directions=="negative, positive",]
ggplot2::ggplot(bar.df, aes(x=meanCentrality, y=name)) + geom_bar(stat="identity", alpha=0.7) +
  geom_errorbar(aes(xmin=meanCentrality-sdCentrality, xmax=meanCentrality+sdCentrality, y=name)) + 
  theme_classic(base_size=12) + labs(x="Mean Centrality", y=NULL) + 
  scale_y_discrete(limits=bar.df[order(bar.df$meanCentrality),]$name)

bar.df2 <- chr8q.reg.centrality[chr8q.reg.centrality$name %in% bar.df$name,]
ggplot2::ggplot(bar.df2, aes(x=eigen_centrality, y=name)) + #geom_violin() +
  geom_boxplot() + geom_point(aes(shape=Direction))+ 
  scale_shape_manual(values=c(2,6), labels=c("Upregulated\nNetwork","Downregulated\nNetwork"), breaks=c("positive","negative")) +
  theme_classic(base_size=12) + labs(x="Centrality", y=NULL) + 
  theme(legend.position="inside", legend.position.inside = c(0.7,0.5)) +
  scale_y_discrete(limits=bar.df[order(bar.df$meanCentrality),]$name)
ggsave("chr8qGeneCentralityTFKinNetworks.pdf",width=3.5,height=3)

ggplot2::ggplot(bar.df2, aes(x=eigen_centrality, y=name)) + #geom_violin() +
  geom_bar(stat="summary", fun="mean", fill="grey") + geom_point(aes(shape=Direction))+ 
  scale_shape_manual(values=c(2,6), labels=c("Upregulated\nNetwork","Downregulated\nNetwork"), breaks=c("positive","negative")) +
  theme_classic(base_size=12) + labs(x="Centrality", y=NULL) + 
  theme(legend.position="inside", legend.position.inside = c(0.7,0.5)) +
  scale_y_discrete(limits=bar.df[order(bar.df$meanCentrality),]$name)
ggsave("chr8qGeneCentralityTFKinNetworks_barPlotWithPoints.pdf",width=3.5,height=3)

ggplot2::ggplot(bar.df2, aes(x=eigen_centrality, y=name, fill=Direction)) +
  geom_bar(position="dodge",stat="identity", alpha=0.7) + 
  scale_fill_manual(values=c("red","blue"), labels=c("Upregulated\nNetwork","Downregulated\nNetwork"), breaks=c("positive","negative")) +
  theme_classic(base_size=12) + labs(x="Centrality", y=NULL) + 
  theme(legend.position="inside", legend.position.inside = c(0.7,0.5)) +
  scale_y_discrete(limits=bar.df[order(bar.df$meanCentrality),]$name)
ggsave("chr8qGeneCentralityTFKinNetworks_barPlot.pdf",width=3.5,height=3)

ggplot2::ggplot(bar.df2, aes(x=eigen_centrality, y=name, fill=Direction)) +
  geom_bar(stat="identity", alpha=0.7) + 
  scale_fill_manual(values=c("red","blue"), labels=c("Upregulated\nNetwork","Downregulated\nNetwork"), breaks=c("positive","negative")) +
  theme_classic(base_size=12) + labs(x="Centrality", y=NULL) + 
  theme(legend.position="inside", legend.position.inside = c(0.7,0.5)) +
  scale_y_discrete(limits=bar.df[order(bar.df$meanCentrality),]$name)
ggsave("chr8qGeneCentralityTFKinNetworks_barPlotStacked.pdf",width=3.5,height=3)
  
# look for connections between chr8q genes and altered TFs/kinases
colnames(neg.vert)[1] <- "from"
neg.tf.con <- neg.edges[(neg.edges$from %in% c(tf.result$Gene, kin.result$Feature_set) | 
                         neg.edges$to %in% c(tf.result$Gene, kin.result$Feature_set)) &
                        (neg.edges$from %in% chr8q.genes | neg.edges$to %in% chr8q.genes),] # 40
neg.tf.con.vert <- neg.vert[neg.vert$from %in% c(neg.tf.con$from, neg.tf.con$to),] # 23
topGraph <- igraph::graph_from_data_frame(neg.tf.con, directed=FALSE, vertices=neg.tf.con.vert) 
#plot(topGraph)
tempTitle <- paste0("negative", "_", nrow(neg.tf.con.vert),"_TFChr8Nodes_manual_",Sys.Date())
RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
write.csv(neg.tf.con.vert, "negative_TFChr8_vertices.csv", row.names=FALSE)
write.csv(neg.tf.con, "negative_TFChr8_edges.csv", row.names=FALSE)
neg.tf.con.centrality <- data.frame(name = V(topGraph)$name,
                                    degree = igraph::degree(topGraph, mode="all"),
                                    closeness = igraph::closeness(topGraph, mode="all"),
                                    betweenness = igraph::betweenness(topGraph, directed = FALSE),
                                    eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                                    hub_score = igraph::hub_score(topGraph)$vector,
                                    authority_score = igraph::authority_score(topGraph)$vector)
write.csv(neg.tf.con.centrality, "negative_TFChr8_centrality.csv", row.names=FALSE)

colnames(pos.vert)[1] <- "from"
pos.tf.con <- pos.edges[(pos.edges$from %in% c(tf.result$Gene, kin.result$Feature_set) | 
                           pos.edges$to %in% c(tf.result$Gene, kin.result$Feature_set)) &
                          (pos.edges$from %in% chr8q.genes | pos.edges$to %in% chr8q.genes),] # 2310
pos.tf.con.vert <- pos.vert[pos.vert$from %in% c(pos.tf.con$from, pos.tf.con$to),] # 908
topGraph <- igraph::graph_from_data_frame(pos.tf.con, directed=FALSE, vertices=pos.tf.con.vert) 
#plot(topGraph)
tempTitle <- paste0("positive", "_", nrow(pos.tf.con.vert),"_TFChr8Nodes_manual_",Sys.Date())
RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
write.csv(pos.tf.con.vert, "positive_TFChr8_vertices.csv", row.names=FALSE)
write.csv(pos.tf.con, "positive_TFChr8_edges.csv", row.names=FALSE)
pos.tf.con.centrality <- data.frame(name = V(topGraph)$name,
                                    degree = igraph::degree(topGraph, mode="all"),
                                    closeness = igraph::closeness(topGraph, mode="all"),
                                    betweenness = igraph::betweenness(topGraph, directed = FALSE),
                                    eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                                    hub_score = igraph::hub_score(topGraph)$vector,
                                    authority_score = igraph::authority_score(topGraph)$vector)
write.csv(pos.tf.con.centrality, "positive_TFChr8_centrality.csv", row.names=FALSE)

# just generally look for chr8q connections with terminals
neg.tf.con <- neg.edges[neg.edges$from %in% chr8q.genes | neg.edges$to %in% chr8q.genes,] # 5050
neg.tf.con.vert <- neg.vert[neg.vert$from %in% c(neg.tf.con$from, neg.tf.con$to),] # 815
topGraph <- igraph::graph_from_data_frame(neg.tf.con, directed=FALSE, vertices=neg.tf.con.vert) 
#plot(topGraph)
tempTitle <- paste0("negative", "_", nrow(neg.tf.con.vert),"_Chr8Nodes_manual_",Sys.Date())
RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
write.csv(neg.tf.con.vert, "negative_Chr8_vertices.csv", row.names=FALSE)
write.csv(neg.tf.con, "negative_Chr8_edges.csv", row.names=FALSE)
neg.tf.con.centrality <- data.frame(name = V(topGraph)$name,
                                    degree = igraph::degree(topGraph, mode="all"),
                                    closeness = igraph::closeness(topGraph, mode="all"),
                                    betweenness = igraph::betweenness(topGraph, directed = FALSE),
                                    eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                                    hub_score = igraph::hub_score(topGraph)$vector,
                                    authority_score = igraph::authority_score(topGraph)$vector)
write.csv(neg.tf.con.centrality, "negative_Chr8_centrality.csv", row.names=FALSE)

pos.tf.con <- pos.edges[pos.edges$from %in% chr8q.genes | pos.edges$to %in% chr8q.genes,] # 5050
pos.tf.con.vert <- pos.vert[pos.vert$from %in% c(pos.tf.con$from, pos.tf.con$to),] # 815
topGraph <- igraph::graph_from_data_frame(pos.tf.con, directed=FALSE, vertices=pos.tf.con.vert) 
#plot(topGraph)
tempTitle <- paste0("positive", "_", nrow(pos.tf.con.vert),"_Chr8Nodes_manual_",Sys.Date())
RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
write.csv(pos.tf.con.vert, "positive_Chr8_vertices.csv", row.names=FALSE)
write.csv(pos.tf.con, "positive_Chr8_edges.csv", row.names=FALSE)
pos.tf.con.centrality <- data.frame(name = V(topGraph)$name,
                                    degree = igraph::degree(topGraph, mode="all"),
                                    closeness = igraph::closeness(topGraph, mode="all"),
                                    betweenness = igraph::betweenness(topGraph, directed = FALSE),
                                    eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                                    hub_score = igraph::hub_score(topGraph)$vector,
                                    authority_score = igraph::authority_score(topGraph)$vector)
write.csv(pos.tf.con.centrality, "positive_Chr8_centrality.csv", row.names=FALSE)

# look for connections between MYC and PLK1/EGFR
neg.tf.con <- neg.edges[neg.edges$from %in% c("MYC", "EGFR") | neg.edges$to %in% c("MYC", "EGFR"),] # 68
neg.tf.con.vert <- neg.vert[neg.vert$from %in% c(neg.tf.con$from, neg.tf.con$to),] # 32
topGraph <- igraph::graph_from_data_frame(neg.tf.con, directed=FALSE, vertices=neg.tf.con.vert) 
#plot(topGraph)
tempTitle <- paste0("negative", "_", nrow(neg.tf.con.vert),"_MYC_EGFR_manual_",Sys.Date())
RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
write.csv(neg.tf.con.vert, "negative_MYC_EGFR_vertices.csv", row.names=FALSE)
write.csv(neg.tf.con, "negative_MYC_EGFR_edges.csv", row.names=FALSE)
neg.tf.con.centrality <- data.frame(name = V(topGraph)$name,
                                    degree = igraph::degree(topGraph, mode="all"),
                                    closeness = igraph::closeness(topGraph, mode="all"),
                                    betweenness = igraph::betweenness(topGraph, directed = FALSE),
                                    eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                                    hub_score = igraph::hub_score(topGraph)$vector,
                                    authority_score = igraph::authority_score(topGraph)$vector)
write.csv(neg.tf.con.centrality, "negative_MYC_EGFR_centrality.csv", row.names=FALSE)


neg.tf.con <- neg.edges[neg.edges$from %in% c("MYC", "HRH1") | neg.edges$to %in% c("MYC", "HRH1"),] # 68
neg.tf.con.vert <- neg.vert[neg.vert$from %in% c(neg.tf.con$from, neg.tf.con$to),] # 32
topGraph <- igraph::graph_from_data_frame(neg.tf.con, directed=FALSE, vertices=neg.tf.con.vert) 
#plot(topGraph)
tempTitle <- paste0("negative", "_", nrow(neg.tf.con.vert),"_MYC_HRH1_manual_",Sys.Date())
RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
write.csv(neg.tf.con.vert, "negative_MYC_HRH1_vertices.csv", row.names=FALSE)
write.csv(neg.tf.con, "negative_MYC_HRH1_edges.csv", row.names=FALSE)
neg.tf.con.centrality <- data.frame(name = V(topGraph)$name,
                                    degree = igraph::degree(topGraph, mode="all"),
                                    closeness = igraph::closeness(topGraph, mode="all"),
                                    betweenness = igraph::betweenness(topGraph, directed = FALSE),
                                    eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                                    hub_score = igraph::hub_score(topGraph)$vector,
                                    authority_score = igraph::authority_score(topGraph)$vector)
write.csv(neg.tf.con.centrality, "negative_MYC_HRH1_centrality.csv", row.names=FALSE)
# MYC and HRH1 not connected in this network
# try pulling from STRING again
neg.tf.con <- STRINGv12[STRINGv12$from %in% c("MYC", "HRH1") | STRINGv12$to %in% c("MYC", "HRH1"),]
from <- unique(c(neg.tf.con$from, neg.tf.con$to))
neg.tf.con.vert <- data.frame(from, type="Inferred", Omics=NA, chr8q=FALSE)
neg.tf.con.vert[neg.tf.con.vert$from %in% sig.global[sig.global$Spearman.est<0,]$Gene,]$Omics <- "Protein"
neg.tf.con.vert[neg.tf.con.vert$from %in% sig.tf[sig.tf$NES<0,]$Gene,]$Omics <- "TF"
neg.tf.con.vert[neg.tf.con.vert$from %in% sig.kin[sig.kin$NES<0,]$Feature_set,]$Omics <- "Kinase"
neg.tf.con.vert[!is.na(neg.tf.con.vert$Omics),]$type <- "Terminal"
neg.tf.con.vert[neg.tf.con.vert$from %in% chr8q.genes,]$chr8q <- TRUE
topGraph <- igraph::graph_from_data_frame(neg.tf.con, directed=FALSE, vertices=neg.tf.con.vert) 
#plot(topGraph)
tempTitle <- paste0("STRING_negative", "_", nrow(neg.tf.con.vert),"_MYC_HRH1_manual_",Sys.Date())
RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
write.csv(neg.tf.con.vert, "STRING_negative_MYC_HRH1_vertices.csv", row.names=FALSE)
write.csv(neg.tf.con, "STRING_negative_MYC_HRH1_edges.csv", row.names=FALSE)
neg.tf.con.centrality <- data.frame(name = V(topGraph)$name,
                                    degree = igraph::degree(topGraph, mode="all"),
                                    closeness = igraph::closeness(topGraph, mode="all"),
                                    betweenness = igraph::betweenness(topGraph, directed = FALSE),
                                    eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                                    hub_score = igraph::hub_score(topGraph)$vector,
                                    authority_score = igraph::authority_score(topGraph)$vector)
write.csv(neg.tf.con.centrality, "STRING_negative_MYC_HRH1_centrality.csv", row.names=FALSE)

pos.tf.con <- pos.edges[pos.edges$from %in% c("MYC", "PLK1") | pos.edges$to %in% c("MYC", "PLK1"),] # 118
pos.tf.con.vert <- pos.vert[pos.vert$from %in% c(pos.tf.con$from, pos.tf.con$to),] # 58
topGraph <- igraph::graph_from_data_frame(pos.tf.con, directed=FALSE, vertices=pos.tf.con.vert) 
#plot(topGraph)
tempTitle <- paste0("positive", "_", nrow(pos.tf.con.vert),"_MYC_PLK1_manual_",Sys.Date())
RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
write.csv(pos.tf.con.vert, "positive_MYC_PLK1_vertices.csv", row.names=FALSE)
write.csv(pos.tf.con, "positive_MYC_PLK1_edges.csv", row.names=FALSE)
pos.tf.con.centrality <- data.frame(name = V(topGraph)$name,
                                    degree = igraph::degree(topGraph, mode="all"),
                                    closeness = igraph::closeness(topGraph, mode="all"),
                                    betweenness = igraph::betweenness(topGraph, directed = FALSE),
                                    eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
                                    hub_score = igraph::hub_score(topGraph)$vector,
                                    authority_score = igraph::authority_score(topGraph)$vector)
write.csv(pos.tf.con.centrality, "positive_MYC_PLK1_centrality.csv", row.names=FALSE)

# how many connections does MYC have in these networks?
pos.myc.con <- pos.edges[pos.edges$from=="MYC"|pos.edges$to=="MYC",]
pos.myc.con <- unique(c(pos.myc.con$from,pos.myc.con$to))
pos.myc.con <- pos.myc.con[pos.myc.con!="MYC"] # 50

neg.myc.con <- neg.edges[neg.edges$from=="MYC"|neg.edges$to=="MYC",]
neg.myc.con <- unique(c(neg.myc.con$from,neg.myc.con$to))
neg.myc.con <- neg.myc.con[neg.myc.con!="MYC"] # 21

myc.con <- unique(c(pos.myc.con, neg.myc.con)) # 71

con <- unique(c(pos.vert$from, neg.vert$from)) # 11084
con <- con[con != "MYC"] # 11083

# how many of those are quantified in proteomics/TF/kinases?
synapser::synLogin()
global.result <- na.omit(read.csv(synapser::synGet("syn66224803")$path))
tf.result <- na.omit(read.csv(synapser::synGet("syn66226952")$path))
tf.result$Gene <- sub("_.*","",tf.result$Feature_set)
kin.result <- read.csv(synapser::synGet("syn66279699")$path)
quant.genes <- unique(c(global.result$Gene, tf.result$Gene, kin.result$Feature_set))
quant.con <- con[con %in% quant.genes] # 6969
quant.myc.con <- myc.con[myc.con %in% quant.genes] # all 71

# how many of those pass our significance thresholds?
sig.global <- global.result[global.result$Spearman.q <= 0.05, ] # 208 / 9013 (2.31%); 98.02% of 101 positive, 100% of 107 negative are in the interactome
sig.tf <- tf.result[tf.result$p_value <= 0.05 & tf.result$FDR_q_value <= 0.25, ] # 206 / 408 (50.49%); 96.1% of 205 positive, 100% of 1 negative are in the interactome
sig.kin <- kin.result[kin.result$p_value <= 0.05 & kin.result$FDR_q_value <= 0.25, ] # 2 / 184 (1.09%); 100% of 2 negative are in the interactome
sig.quant <- unique(c(sig.global$Gene, sig.tf$Gene, sig.kin$Feature_set)) # 415
sig.quant.con <- quant.con[quant.con %in% sig.quant] # 408
sig.quant.myc.con <- quant.myc.con[quant.myc.con %in% sig.quant] # all 71

# fisher's exact test: are MYC's connections which we quantified altered with Chr8q status more than by random chance?
non.myc.quant.con <- quant.con[!(quant.con %in% quant.myc.con)] # 6898
sig.non.myc.quant.con <- sig.quant.con[!(sig.quant.con %in% quant.myc.con)] # 337
contingency.matrix <- data.frame("MYC_interaction" = c(length(sig.quant.myc.con), # 71
                                                       length(quant.myc.con) - length(sig.quant.myc.con)), # 0
                                 "Non_MYC_interaction" = c(length(sig.non.myc.quant.con), # 337
                                                           length(non.myc.quant.con) - length(sig.non.myc.quant.con)), # 6561
                                 row.names = c("Chr8q_association", "Non_associated"))
# contingency.matrix <- data.frame("MYC_interaction" = c(length(sig.quant.myc.con), 
#                                                        length(quant.myc.con) - length(sig.quant.myc.con)),
#                                  "Non_MYC_interaction" = c(length(sig.quant.con) - length(sig.quant.myc.con), 
#                                                            length(quant.con) - length(sig.quant.con)),
#                                  row.names = c("Chr8q_association", "Non_associated"))
if (length(quant.con) == sum(c(contingency.matrix$MYC_interaction,contingency.matrix$Non_MYC_interaction))) {
  mosaicplot(contingency.matrix, color=TRUE)
  fisher.result <- fisher.test(contingency.matrix) # < 2.2e-16
  fisher.p <- fisher.result$p.value # 6.83817802826197E-91
} else {
  warning("number of genes in matrix does not match number quantified in network")
}

# how many connections does MYC have in STRING network?
myc.con <- STRINGv12[STRINGv12$from == "MYC" | STRINGv12$to == "MYC",]
myc.con <- unique(c(myc.con$from, myc.con$to))
myc.con <- myc.con[myc.con != "MYC"] # 2023

con <- unique(c(STRINGv12$from, STRINGv12$to)) # 18767
con <- con[con != "MYC"] # 18766

# how many of those are quantified in proteomics/TF/kinases?
synapser::synLogin()
global.result <- na.omit(read.csv(synapser::synGet("syn66224803")$path))
tf.result <- na.omit(read.csv(synapser::synGet("syn66226952")$path))
tf.result$Gene <- sub("_.*","",tf.result$Feature_set)
kin.result <- read.csv(synapser::synGet("syn66279699")$path)
quant.genes <- unique(c(global.result$Gene, tf.result$Gene, kin.result$Feature_set)) # 9281
quant.con <- con[con %in% quant.genes] # 9115
quant.myc.con <- myc.con[myc.con %in% quant.genes] # 1711

# how many of those pass our significance thresholds?
sig.global <- global.result[global.result$Spearman.q <= 0.05, ] # 208 / 9013 (2.31%); 98.02% of 101 positive, 100% of 107 negative are in the interactome
sig.tf <- tf.result[tf.result$p_value <= 0.05 & tf.result$FDR_q_value <= 0.25, ] # 206 / 408 (50.49%); 96.1% of 205 positive, 100% of 1 negative are in the interactome
sig.kin <- kin.result[kin.result$p_value <= 0.05 & kin.result$FDR_q_value <= 0.25, ] # 2 / 184 (1.09%); 100% of 2 negative are in the interactome
sig.quant <- unique(c(sig.global$Gene, sig.tf$Gene, sig.kin$Feature_set)) # 415
sig.pos.quant <- unique(c(sig.global[sig.global$Spearman.est>0,]$Gene, 
                          sig.tf[sig.tf$NES>0,]$Gene, 
                          sig.kin[sig.kin$NES>0,]$Feature_set)) # 305
sig.neg.quant <- unique(c(sig.global[sig.global$Spearman.est<0,]$Gene, 
                          sig.tf[sig.tf$NES<0,]$Gene, 
                          sig.kin[sig.kin$NES<0,]$Feature_set)) # 110
sig.quant.con <- quant.con[quant.con %in% sig.quant] # 408
sig.quant.myc.con <- quant.myc.con[quant.myc.con %in% sig.quant] # 71

sig.quant.pos.con <- quant.con[quant.con %in% sig.pos.quant] # 299
sig.quant.pos.myc.con <- quant.myc.con[quant.myc.con %in% sig.pos.quant] # 50

sig.quant.neg.con <- quant.con[quant.con %in% sig.neg.quant] # 109
sig.quant.neg.myc.con <- quant.myc.con[quant.myc.con %in% sig.neg.quant] # 20

# fisher's exact test: are MYC's connections which we quantified altered with Chr8q status more than by random chance?
non.myc.quant.con <- quant.con[!(quant.con %in% quant.myc.con)] # 7404
sig.non.myc.quant.con <- sig.quant.con[!(sig.quant.con %in% quant.myc.con)] # 337
contingency.matrix <- data.frame("MYC_interaction" = c(length(sig.quant.myc.con), # 71
                                                       length(quant.myc.con) - length(sig.quant.myc.con)), # 1640
                                 "Non_MYC_interaction" = c(length(sig.non.myc.quant.con), # 337
                                                           length(non.myc.quant.con) - length(sig.non.myc.quant.con)), # 7067
                                 row.names = c("Chr8q_association", "Non_associated"))
if (length(quant.con) == sum(c(contingency.matrix$MYC_interaction,contingency.matrix$Non_MYC_interaction))) {
  mosaicplot(contingency.matrix, color=TRUE)
  fisher.result <- fisher.test(contingency.matrix) # 0.5165
  fisher.p <- fisher.result$p.value # 0.516528134580918
} else {
  warning("number of genes in matrix does not match number quantified in network")
}

non.myc.quant.con <- quant.con[!(quant.con %in% quant.myc.con)] # 7404
sig.pos.non.myc.quant.con <- sig.quant.pos.con[!(sig.quant.pos.con %in% quant.myc.con)] # 337
contingency.matrix <- data.frame("MYC_interaction" = c(length(sig.quant.pos.myc.con), # 71
                                                       length(quant.myc.con) - length(sig.quant.pos.myc.con)), # 1640
                                 "Non_MYC_interaction" = c(length(sig.pos.non.myc.quant.con), # 337
                                                           length(non.myc.quant.con) - length(sig.pos.non.myc.quant.con)), # 7067
                                 row.names = c("Chr8q_association", "Non_associated"))
if (length(quant.con) == sum(c(contingency.matrix$MYC_interaction,contingency.matrix$Non_MYC_interaction))) {
  mosaicplot(contingency.matrix, color=TRUE)
  fisher.result <- fisher.test(contingency.matrix) # 0.4073
  fisher.p <- fisher.result$p.value # 0.407280880077656
} else {
  warning("number of genes in matrix does not match number quantified in network")
}

non.myc.quant.con <- quant.con[!(quant.con %in% quant.myc.con)] # 7404
sig.neg.non.myc.quant.con <- sig.quant.neg.con[!(sig.quant.neg.con %in% quant.myc.con)] # 337
contingency.matrix <- data.frame("MYC_interaction" = c(length(sig.quant.neg.myc.con), # 21
                                                       length(quant.myc.con) - length(sig.quant.neg.myc.con)), # 1690
                                 "Non_MYC_interaction" = c(length(sig.neg.non.myc.quant.con), # 88
                                                           length(non.myc.quant.con) - length(sig.neg.non.myc.quant.con)), # 7316
                                 row.names = c("Chr8q_association", "Non_associated"))
if (length(quant.con) == sum(c(contingency.matrix$MYC_interaction,contingency.matrix$Non_MYC_interaction))) {
  mosaicplot(contingency.matrix, color=TRUE)
  fisher.result <- fisher.test(contingency.matrix) # 0.9017634
  fisher.p <- fisher.result$p.value # 0.90176337076464
} else {
  warning("number of genes in matrix does not match number quantified in network")
}

#### Fisher's exact test for all non-terminal nodes ####
con <- unique(c(STRINGv12$from, STRINGv12$to)) # 18767
# how many of those are quantified in proteomics/TF/kinases?
synapser::synLogin()
global.result <- na.omit(read.csv(synapser::synGet("syn66224803")$path))
tf.result <- na.omit(read.csv(synapser::synGet("syn66226952")$path))
tf.result$Gene <- sub("_.*","",tf.result$Feature_set)
kin.result <- read.csv(synapser::synGet("syn66279699")$path)
quant.genes <- unique(c(global.result$Gene, tf.result$Gene, kin.result$Feature_set)) # 9281
quant.con <- con[con %in% quant.genes] # 9115

# how many of those pass our significance thresholds?
sig.global <- global.result[global.result$Spearman.q <= 0.05, ] # 208 / 9013 (2.31%); 98.02% of 101 positive, 100% of 107 negative are in the interactome
sig.tf <- tf.result[tf.result$p_value <= 0.05 & tf.result$FDR_q_value <= 0.25, ] # 206 / 408 (50.49%); 96.1% of 205 positive, 100% of 1 negative are in the interactome
sig.kin <- kin.result[kin.result$p_value <= 0.05 & kin.result$FDR_q_value <= 0.25, ] # 2 / 184 (1.09%); 100% of 2 negative are in the interactome
sig.quant <- unique(c(sig.global$Gene, sig.tf$Gene, sig.kin$Feature_set)) # 415
sig.pos.quant <- unique(c(sig.global[sig.global$Spearman.est>0,]$Gene, 
                          sig.tf[sig.tf$NES>0,]$Gene, 
                          sig.kin[sig.kin$NES>0,]$Feature_set)) # 305
sig.neg.quant <- unique(c(sig.global[sig.global$Spearman.est<0,]$Gene, 
                          sig.tf[sig.tf$NES<0,]$Gene, 
                          sig.kin[sig.kin$NES<0,]$Feature_set)) # 110
sig.quant.con <- quant.con[quant.con %in% sig.quant] # 408
sig.quant.pos.con <- quant.con[quant.con %in% sig.pos.quant] # 299
sig.quant.neg.con <- quant.con[quant.con %in% sig.neg.quant] # 109

directions <- list("positive" = list(#"edges" = pos.edges,
                                  "vert" = pos.vert),
                "negative" = list(#"edges" = neg.edges,
                                  "vert" = neg.vert))
mosaics <- list()
contingency.matrices <- list()
fisher.df <- data.frame()
for (dir in names(directions)) {
  #temp.edges <- directions[[dir]][["edges"]]
  temp.vert <- directions[[dir]][["vert"]]
  
  inferred <- unique(temp.vert[temp.vert$type == "Inferred" & is.na(temp.vert$Omics),]$from)
  titleSigInt <- ifelse(dir=="positive"," interactors tend to be Chr8q-upregulated proteins (p = ",
                   " interactors tend to be Chr8q-downregulated proteins (p = ")
  titleSig <- ifelse(dir=="positive"," interactors tend not to be Chr8q-upregulated proteins (p = ",
                     " interactors tend not to be Chr8q-downregulated proteins (p = ")
  titleInsig <- ifelse(dir=="positive"," interactors are equally likely to be Chr8q-upregulated proteins (p = ",
                     " interactors are equally likely to be Chr8q-downregulated proteins (p = ")
  dir.rowname <- ifelse(dir=="positive", "Up with Chr8q", "Down with Chr8q")
  dir.mosaics <- list()
  dir.matrices <- list()
  dir.fisher <- data.frame(inferred, Direction = dir, N_interactor_terminals = NA, 
                           N_noninteractor_terminals = NA, p = NA)
  for (node in inferred) {
    temp.con <- con[con != node]
   
    # how many connections does node have in STRING network?
    node.con <- STRINGv12[STRINGv12$from == node | STRINGv12$to == node,]
    node.con <- unique(c(node.con$from, node.con$to))
    node.con <- node.con[node.con != node] # 2023 
    
    # how many of those are quantified
    if (dir == "positive") {
      dir.quant.con <- quant.con[!(quant.con %in% sig.quant.neg.con)] # exclude negative terminals from positive analysis
    } else {
      dir.quant.con <- quant.con[!(quant.con %in% sig.quant.pos.con)] # exclude positive terminals from negative analysis
    }
    quant.node.con <- node.con[node.con %in% dir.quant.con] # 1711
    
    # how many are up- or down-regulated?
    sig.quant.node.con <- quant.node.con[quant.node.con %in% sig.quant]
    if (dir == "positive") {
      sig.quant.dir.node.con <- quant.node.con[quant.node.con %in% sig.pos.quant] 
      sig.quant.dir.con <- sig.quant.pos.con
    } else {
      sig.quant.dir.node.con <- quant.node.con[quant.node.con %in% sig.neg.quant]
      sig.quant.dir.con <- sig.quant.neg.con
    }
    
    # prepare contingency matrix
    non.node.quant.con <- dir.quant.con[!(dir.quant.con %in% quant.node.con)]
    sig.dir.non.node.quant.con <- sig.quant.dir.con[!(sig.quant.dir.con %in% quant.node.con)]
    contingency.matrix <- data.frame("node_interaction" = c(length(sig.quant.dir.node.con),
                                                           length(quant.node.con) - length(sig.quant.dir.node.con)),
                                     "Non_node_interaction" = c(length(sig.dir.non.node.quant.con),
                                                               length(non.node.quant.con) - length(sig.dir.non.node.quant.con)),
                                     row.names = c(dir.rowname,"Unaltered with Chr8q"))
    colnames(contingency.matrix) <- paste0(node, c(" Interactor", " Non-interactor"))
    if (length(dir.quant.con) == sum(c(contingency.matrix[,1],contingency.matrix[,2]))) {
      fisher.p <- fisher.test(contingency.matrix)$p.value
      title <- ifelse(fisher.p <= 0.05 & contingency.matrix[1,1] > contingency.matrix[2,1], 
                      paste0(node, titleSigInt, fisher.p, ")"), # if tends to be chr8q terminal interactor
                      ifelse(fisher.p <= 0.05, paste0(node, titleSig, fisher.p, ")"), # if tends to NOT be chr8q terminal interactor
                             paste0(node, titleInsig, fisher.p, ")"))) # if not significant difference
      dir.mosaics[[node]] <- mosaicplot(contingency.matrix, main = title, color=TRUE)
      dir.matrices[[node]] <- contingency.matrix
      dir.fisher[dir.fisher$inferred==node,]$N_interactor_terminals <- contingency.matrix[1,1]
      dir.fisher[dir.fisher$inferred==node,]$N_noninteractor_terminals <- contingency.matrix[1,2]
      dir.fisher[dir.fisher$inferred==node,]$N_interactor_nonterminals <- contingency.matrix[2,1]
      dir.fisher[dir.fisher$inferred==node,]$N_noninteractor_nonterminals <- contingency.matrix[2,2]
      dir.fisher[dir.fisher$inferred==node,]$p <- fisher.p
    } else {
      stop("number of genes in matrix does not match number quantified in network")
    }
  }
  mosaics[[dir]] <- dir.mosaics
  contingency.matrices[[dir]] <- dir.matrices
  fisher.df <- rbind(fisher.df, dir.fisher)
}
fisher.df$Terminal_interactor_ratio <- fisher.df$N_interactor_terminals/(fisher.df$N_noninteractor_terminals + fisher.df$N_interactor_terminals)
fisher.df$minusLogP <- -log10(fisher.df$p)
fisher.df$chr8q <- FALSE
fisher.df[fisher.df$inferred %in% chr8q.genes,]$chr8q <- TRUE
fisher.df$drugTarget <- FALSE
drug.info <- read.csv("https://raw.githubusercontent.com/BelindaBGarana/DMEA/refs/heads/shiny-app/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv")
drug.targets <- unique(unlist(strsplit(unlist(drug.info$target), ", "))) # 1027
fisher.df[fisher.df$inferred %in% drug.targets,]$drugTarget <- TRUE
write.csv(fisher.df, "FishersExactTest_inferredSTRINGNodes.csv",row.names=FALSE)
if (any(fisher.df$p > 1)) {
  fisher.df[fisher.df$p>1,]$p <- 1
}
fisher.df$q <- NA
fisher.df$q <- qvalue::qvalue(fisher.df$p,pi0=1)$qvalues
fisher.df$chr8q_q <- NA
fisher.df[fisher.df$chr8q==TRUE,]$chr8q_q <- qvalue::qvalue(fisher.df[fisher.df$chr8q==TRUE,]$p,pi0=1)$qvalues
fisher.df$drugTarget_q <- NA
fisher.df[fisher.df$drugTarget==TRUE,]$drugTarget_q <- qvalue::qvalue(fisher.df[fisher.df$drugTarget==TRUE,]$p,pi0=1)$qvalues
write.csv(fisher.df, "FishersExactTest_inferredSTRINGNodes_adjustedP.csv",row.names=FALSE)
saveRDS(contingency.matrices, "contingencyMatrices_inferredSTRINGNodes.rds")
saveRDS(mosaics, "mosaics_inferredSTRINGNodes.rds")

synapser::synLogin()
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
sig.pos.corr.targets <- na.omit(unique(unlist(strsplit(drug.corr.wInfo[drug.corr.wInfo$Pearson.q <= 0.05 & drug.corr.wInfo$Pearson.est>0,]$target, ", ")))) # 88
fisher.df$drugCorrPos <- FALSE
fisher.df[fisher.df$inferred %in% sig.pos.corr.targets,]$drugCorrPos <- TRUE
sig.neg.corr.targets <- na.omit(unique(unlist(strsplit(drug.corr.wInfo[drug.corr.wInfo$Pearson.q <= 0.05 & drug.corr.wInfo$Pearson.est>0,]$target, ", ")))) # 88
fisher.df$drugCorrNeg <- FALSE
fisher.df[fisher.df$inferred %in% sig.neg.corr.targets,]$drugCorrNeg <- TRUE
fisher.df$drugCorr <- FALSE
fisher.df[fisher.df$drugCorrNeg | fisher.df$drugCorrPos,]$drugCorr <- TRUE

posMOAs <- na.omit(unique(gsea.rna.prot[gsea.rna.prot$p_value <= 0.05 & gsea.rna.prot$FDR_q_value <= 0.25 &
                           gsea.rna.prot$NES > 0,]$Drug_set))
negMOAs <- na.omit(unique(gsea.rna.prot[gsea.rna.prot$p_value <= 0.05 & gsea.rna.prot$FDR_q_value <= 0.25 &
                                          gsea.rna.prot$NES < 0,]$Drug_set))
"VEGFR inhibitor" %in% c(posMOAs,negMOAs) # FALSE: make sure we don't get VEGFRi when looking for EGFRi
posMOAtargets <- c()
negMOAtargets <- c()
for (i in 1:nrow(red.drug.info)) {
  tempMOAs <- na.omit(unique(strsplit(red.drug.info$moa[i], ", ")))
  tempTargets <- na.omit(unique(strsplit(red.drug.info$target[i], ", ")))
  if (any(tempMOAs %in% posMOAs)) {
    posMOAtargets <- c(posMOAtargets, red.drug.info$target[i])
  }
  if (any(tempMOAs %in% negMOAs)) {
    negMOAtargets <- c(negMOAtargets, red.drug.info$target[i])
  }
}
fisher.df$drugMOATargetPos <- FALSE
fisher.df[fisher.df$inferred %in% posMOAtargets,]$drugMOATargetPos <- TRUE
fisher.df$drugMOATargetNeg <- FALSE
fisher.df[fisher.df$inferred %in% negMOAtargets,]$drugMOATargetNeg <- TRUE
write.csv(fisher.df, "FishersExactTest_inferredSTRINGNodes_adjustedP_wDrugInfo.csv",row.names=FALSE)
fisher.df <- read.csv("FishersExactTest_inferredSTRINGNodes_adjustedP_wDrugInfo.csv")
contingency.matrices <- readRDS("contingencyMatrices_inferredSTRINGNodes.rds")
chr8q.nodes <- fisher.df[fisher.df$chr8q,]$inferred
fisher.df[,c("N_interactor_nonterminals","N_noninteractor_nonterminals")] <- NA
for (i in chr8q.nodes) {
  directions <- na.omit(unique(fisher.df[fisher.df$inferred==i,]$Direction))
  for (dir in directions) {
    temp.mat <- contingency.matrices[[dir]][[i]]
    fisher.df[fisher.df$inferred == i & fisher.df$Direction == dir,]$N_interactor_nonterminals <- temp.mat[2,1]
    fisher.df[fisher.df$inferred == i & fisher.df$Direction == dir,]$N_noninteractor_nonterminals <- temp.mat[2,2]
  }
}
fisher.df$Nonterminal_interactor_ratio <- fisher.df$N_interactor_nonterminals/(fisher.df$N_noninteractor_nonterminals+fisher.df$N_interactor_nonterminals)
fisher.df$Terminal_minus_nonterminal_interactor_ratio <- fisher.df$Terminal_interactor_ratio - fisher.df$Nonterminal_interactor_ratio

target.nodes <- fisher.df[fisher.df$drugTarget,]$inferred
for (i in target.nodes) {
  directions <- na.omit(unique(fisher.df[fisher.df$inferred==i,]$Direction))
  for (dir in directions) {
    temp.mat <- contingency.matrices[[dir]][[i]]
    fisher.df[fisher.df$inferred == i & fisher.df$Direction == dir,]$N_interactor_nonterminals <- temp.mat[2,1]
    fisher.df[fisher.df$inferred == i & fisher.df$Direction == dir,]$N_noninteractor_nonterminals <- temp.mat[2,2]
  }
}
fisher.df$Nonterminal_interactor_ratio <- fisher.df$N_interactor_nonterminals/(fisher.df$N_noninteractor_nonterminals+fisher.df$N_interactor_nonterminals)
fisher.df$Terminal_minus_nonterminal_interactor_ratio <- fisher.df$Terminal_interactor_ratio - fisher.df$Nonterminal_interactor_ratio

other.nodes <- fisher.df[is.na(fisher.df$Nonterminal_interactor_ratio),]$inferred
for (i in other.nodes) {
  directions <- na.omit(unique(fisher.df[fisher.df$inferred==i,]$Direction))
  for (dir in directions) {
    temp.mat <- contingency.matrices[[dir]][[i]]
    fisher.df[fisher.df$inferred == i & fisher.df$Direction == dir,]$N_interactor_nonterminals <- temp.mat[2,1]
    fisher.df[fisher.df$inferred == i & fisher.df$Direction == dir,]$N_noninteractor_nonterminals <- temp.mat[2,2]
  }
}
fisher.df$Nonterminal_interactor_ratio <- fisher.df$N_interactor_nonterminals/(fisher.df$N_noninteractor_nonterminals+fisher.df$N_interactor_nonterminals)
fisher.df$Terminal_minus_nonterminal_interactor_ratio <- fisher.df$Terminal_interactor_ratio - fisher.df$Nonterminal_interactor_ratio
fisher.df$minusLogP_signed <- fisher.df$minusLogP * sign(fisher.df$Terminal_minus_nonterminal_interactor_ratio)
write.csv(fisher.df, "FishersExactTest_inferredSTRINGNodes_adjustedP_wDrugInfo_v2.csv",row.names=FALSE)

# pseudo volcano plot
dot.df <- fisher.df
dot.df$Score <- dot.df$Terminal_interactor_ratio
dot.df[dot.df$Direction=="negative",]$Score <- -dot.df[dot.df$Direction=="negative",]$Score
dot.df$chr8q <- factor(dot.df$chr8q, levels=c(TRUE, FALSE))
dot.df$drugTarget <- factor(dot.df$drugTarget, levels=c(TRUE, FALSE))
n.sig.pos <- nrow(dot.df[dot.df$Direction=="positive" & dot.df$p<=0.05,]) # 658
n.pos <- nrow(dot.df[dot.df$Direction=="positive",]) # 8834
perc.sig.pos <- round(100*n.sig.pos/n.pos,2) # 7.45%
n.sig.neg <- nrow(dot.df[dot.df$Direction=="negative" & dot.df$p<=0.05,]) # 803
n.neg <- nrow(dot.df[dot.df$Direction=="negative",]) # 6265
perc.sig.neg <- round(100*n.sig.neg/n.neg,2) # 12.82%
n.sig <- nrow(dot.df[dot.df$p<=0.05,]) # 1461
n.total <- nrow(dot.df) # 8834
perc.sig <- round(100*n.sig/n.total,2) # 9.68%
title <- paste0(n.sig," / ",n.total, " (", perc.sig, "%) Proteins More Likely to\nInteract with Chr8q-affected Proteins")
ggplot2::ggplot(dot.df, aes(x=Score, y = -log10(p), shape=chr8q, 
                            color=drugTarget, size=N_interactor_terminals)) + 
  geom_hline(yintercept=-log10(0.05),color="lightgrey", linetype="dashed") +
  geom_point() + theme_classic(base_size=12) + scale_color_manual(values=c("red", "darkgrey"), 
                                                                  breaks=levels(dot.df$drugTarget)) + 
  scale_shape_manual(values=c(17,19)) + ggrepel::geom_text_repel(data=subset(dot.df, p <= 0.05), 
                                                                 aes(label=inferred), show.legend=FALSE) + 
  labs(x="Fraction of Interactions with Chr8q-affected Proteins",y=expression(paste("-log"[10]," (P-value)")),
       size = "# of Interactions", color="Drug Target", shape = "Chr8q Gene") +
  geom_point(data = subset(dot.df, p <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 21) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
ggsave("inferredSTRINGNodes_volcanoPlot_p.pdf",width=8,height=6)

n.sig <- nrow(dot.df[dot.df$q<=0.05,]) # 79
n.total <- nrow(dot.df) # 8834
perc.sig <- round(100*n.sig/n.total,2) # 0.52%
title <- paste0(n.sig," / ",n.total, " (", perc.sig, "%) Proteins More Likely to\nInteract with Chr8q-affected Proteins")
ggplot2::ggplot(dot.df, aes(x=Score, y = -log10(p), shape=chr8q, 
                            color=drugTarget, size=N_interactor_terminals)) + 
  geom_hline(yintercept=-log10(0.05),color="lightgrey", linetype="dashed") +
  geom_point() + theme_classic(base_size=12) + scale_color_manual(values=c("red", "darkgrey"), 
                                                      breaks=levels(dot.df$drugTarget)) + 
  scale_shape_manual(values=c(17,19)) + ggrepel::geom_text_repel(data=subset(dot.df, p <= 0.05), 
                                                                 aes(label=inferred), show.legend=FALSE) + 
  labs(x="Fraction of Interactions with Chr8q-affected Proteins",y=expression(paste("-log"[10]," (P-value)")),
       size = "# of Interactions", color="Drug Target", shape = "Chr8q Gene") +
  geom_point(data = subset(dot.df, q <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 21) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
ggsave("inferredSTRINGNodes_volcanoPlot_q.pdf",width=8,height=6)

ggplot2::ggplot(dot.df, aes(x=Score, y = -log10(p), color=chr8q, size=N_interactor_terminals)) + 
  geom_hline(yintercept=-log10(0.05),color="lightgrey", linetype="dashed") +
  geom_point() + theme_classic(base_size=12) + scale_color_manual(values=c("red", "darkgrey"), 
                                                                  breaks=levels(dot.df$chr8q)) + 
  ggrepel::geom_text_repel(data=subset(dot.df, q <= 0.05), aes(label=inferred), show.legend=FALSE) + 
  labs(x="Fraction of Interactions with Chr8q-affected Proteins",y=expression(paste("-log"[10]," (P-value)")),
       size = "# of Interactions", color="Chr8q Gene") +
  geom_point(data = subset(dot.df, q <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 21) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
ggsave("inferredSTRINGNodes_volcanoPlot_q_v2.pdf",width=8,height=6)

ggplot2::ggplot(dot.df, aes(x=Score, y = -log10(p), color=chr8q, size=N_interactor_terminals)) + 
  geom_hline(yintercept=-log10(0.05),color="lightgrey", linetype="dashed") +
  geom_point() + theme_classic(base_size=12) + scale_color_manual(values=c("red", "darkgrey"), 
                                                                  breaks=levels(dot.df$chr8q)) + 
  ggrepel::geom_text_repel(data=subset(dot.df, p <= 0.05), aes(label=inferred), show.legend=FALSE) + 
  labs(x="Fraction of Interactions with Chr8q-affected Proteins",y=expression(paste("-log"[10]," (P-value)")),
       size = "# of Interactions", color="Chr8q Gene") +
  geom_point(data = subset(dot.df, q <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 21) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
ggsave("inferredSTRINGNodes_volcanoPlot_q_v3.pdf",width=8,height=6)

dot.df$minusLogP_signed <- dot.df$minusLogP * sign(dot.df$Terminal_minus_nonterminal_interactor_ratio)
ggplot2::ggplot(dot.df, aes(x=Score, y = minusLogP_signed, color=chr8q, size=N_interactor_terminals)) + 
  geom_hline(yintercept=-log10(0.05),color="lightgrey", linetype="dashed") +
  geom_hline(yintercept=log10(0.05),color="lightgrey", linetype="dashed") +
  geom_point() + theme_classic(base_size=12) + scale_color_manual(values=c("red", "darkgrey"), 
                                                                  breaks=levels(dot.df$chr8q)) + 
  ggrepel::geom_text_repel(data=subset(dot.df, p <= 0.05), aes(label=inferred), show.legend=FALSE) + 
  labs(x="Fraction of Interactions with Chr8q-affected Proteins",y=expression(paste("-log"[10]," (P-value)")),
       size = "# of Interactions", color="Chr8q Gene") +
  geom_point(data = subset(dot.df, q <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 21) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
ggsave("inferredSTRINGNodes_volcanoPlot_q_v4.pdf",width=8,height=8)

n.sig <- nrow(dot.df[dot.df$q<=0.05,]) # 79
n.total <- nrow(dot.df) # 8834
perc.sig <- round(100*n.sig/n.total,2) # 0.52%
title <- paste0(n.sig," / ",n.total, " (", perc.sig, "%) Proteins More Likely to\nInteract with Chr8q-affected Proteins")
ggplot2::ggplot(dot.df, aes(x=Score, y = minusLogP_signed, color=chr8q, shape = drugTarget, size=N_interactor_terminals)) +
  geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  geom_hline(yintercept=-log10(0.05),color="lightgrey", linetype="dashed") +
  geom_hline(yintercept=log10(0.05),color="lightgrey", linetype="dashed") +
  geom_point() + theme_minimal(base_size=12) + scale_color_manual(values=c("red", "darkgrey"), 
                                                                  breaks=levels(dot.df$chr8q)) + scale_shape_manual(values=c(17,19)) +
  ggrepel::geom_text_repel(data=subset(dot.df, p <= 0.05), aes(label=inferred), show.legend=FALSE) + 
  labs(x="Fraction of Interactions with Chr8q-affected Proteins",y=expression(paste("-log"[10]," (P-value)")),
       size = "# of Interactions", color="Chr8q Gene", shape = "Drug Target") +
  geom_point(data = subset(dot.df, chr8q == TRUE & q <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 21) +
  geom_point(data = subset(dot.df, drugTarget == TRUE & q <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 24) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
ggsave("inferredSTRINGNodes_volcanoPlot_q_v5.pdf",width=10,height=8)

dot.drug <- dot.df[dot.df$drugTarget==TRUE | dot.df$chr8q==TRUE,]
n.sig <- nrow(dot.drug[dot.drug$p<=0.05,]) # 141
n.total <- nrow(dot.drug) # 1372
perc.sig <- round(100*n.sig/n.total,2) # 10.28%
title <- paste0(n.sig," / ",n.total, " (", perc.sig, "%) Proteins More Likely to\nInteract with Chr8q-affected Proteins")
ggplot2::ggplot(dot.drug, aes(x=Score, y = -log10(p), shape=chr8q, 
                            color=drugTarget, size=N_interactor_terminals)) + 
  geom_hline(yintercept=-log10(0.05),color="lightgrey", linetype="dashed") +
  geom_point() + theme_classic(base_size=12) + scale_color_manual(values=c("red", "darkgrey"), 
                                                                  breaks=levels(dot.drug$drugTarget)) + 
  scale_shape_manual(values=c(17,19)) + ggrepel::geom_text_repel(data=subset(dot.drug, p <= 0.05), 
                                                                 aes(label=inferred), show.legend=FALSE) + 
  labs(x="Fraction of Interactions with Chr8q-affected Proteins",y=expression(paste("-log"[10]," (P-value)")),
       size = "# of Interactions", color="Drug Target", shape = "Chr8q Gene") +
  geom_point(data = subset(dot.drug, p <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 21) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
ggsave("inferredSTRINGNodes_drugTargetsAndChr8q_volcanoPlot_p.pdf",width=8,height=6)

n.sig <- nrow(dot.drug[dot.drug$q<=0.05,]) # 79
n.total <- nrow(dot.drug) # 8834
perc.sig <- round(100*n.sig/n.total,2) # 0.52%
title <- paste0(n.sig," / ",n.total, " (", perc.sig, "%) Proteins More Likely to\nInteract with Chr8q-affected Proteins")
ggplot2::ggplot(dot.drug, aes(x=Score, y = -log10(p), shape=chr8q, 
                            color=drugTarget, size=N_interactor_terminals)) + 
  geom_hline(yintercept=-log10(0.05),color="lightgrey", linetype="dashed") +
  geom_point() + theme_classic(base_size=12) + scale_color_manual(values=c("red", "darkgrey"), 
                                                                  breaks=levels(dot.drug$drugTarget)) + 
  scale_shape_manual(values=c(17,19)) + ggrepel::geom_text_repel(data=subset(dot.drug, p <= 0.05), 
                                                                 aes(label=inferred), show.legend=FALSE) + 
  labs(x="Fraction of Interactions with Chr8q-affected Proteins",y=expression(paste("-log"[10]," (P-value)")),
       size = "# of Interactions", color="Drug Target", shape = "Chr8q Gene") +
  geom_point(data = subset(dot.drug, q <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 21) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
ggsave("inferredSTRINGNodes_drugTargetsAndChr8q_volcanoPlot_q.pdf",width=8,height=6)

# repeat for terminals
directions <- list("positive" = list(#"edges" = pos.edges,
  "vert" = pos.vert),
  "negative" = list(#"edges" = neg.edges,
    "vert" = neg.vert))
con <- unique(c(STRINGv12$from, STRINGv12$to)) # 18767
# how many of those are quantified in proteomics/TF/kinases?
synapser::synLogin()
global.result <- na.omit(read.csv(synapser::synGet("syn66224803")$path))
tf.result <- na.omit(read.csv(synapser::synGet("syn66226952")$path))
tf.result$Gene <- sub("_.*","",tf.result$Feature_set)
kin.result <- read.csv(synapser::synGet("syn66279699")$path)
quant.genes <- unique(c(global.result$Gene, tf.result$Gene, kin.result$Feature_set)) # 9281
quant.con <- con[con %in% quant.genes] # 9115

# how many of those pass our significance thresholds?
sig.global <- global.result[global.result$Spearman.q <= 0.05, ] # 208 / 9013 (2.31%); 98.02% of 101 positive, 100% of 107 negative are in the interactome
sig.tf <- tf.result[tf.result$p_value <= 0.05 & tf.result$FDR_q_value <= 0.25, ] # 206 / 408 (50.49%); 96.1% of 205 positive, 100% of 1 negative are in the interactome
sig.kin <- kin.result[kin.result$p_value <= 0.05 & kin.result$FDR_q_value <= 0.25, ] # 2 / 184 (1.09%); 100% of 2 negative are in the interactome
sig.quant <- unique(c(sig.global$Gene, sig.tf$Gene, sig.kin$Feature_set)) # 415
sig.pos.quant <- unique(c(sig.global[sig.global$Spearman.est>0,]$Gene, 
                          sig.tf[sig.tf$NES>0,]$Gene, 
                          sig.kin[sig.kin$NES>0,]$Feature_set)) # 305
sig.neg.quant <- unique(c(sig.global[sig.global$Spearman.est<0,]$Gene, 
                          sig.tf[sig.tf$NES<0,]$Gene, 
                          sig.kin[sig.kin$NES<0,]$Feature_set)) # 110
sig.quant.con <- quant.con[quant.con %in% sig.quant] # 408
sig.quant.pos.con <- quant.con[quant.con %in% sig.pos.quant] # 299
sig.quant.neg.con <- quant.con[quant.con %in% sig.neg.quant] # 109

mosaics <- list()
contingency.matrices <- list()
fisher.df <- data.frame()
for (dir in names(directions)) {
  #temp.edges <- directions[[dir]][["edges"]]
  temp.vert <- directions[[dir]][["vert"]]
  
  terminal <- unique(temp.vert[temp.vert$type == "Terminal",]$from)
  titleSigInt <- ifelse(dir=="positive"," interactors tend to be Chr8q-upregulated proteins (p = ",
                        " interactors tend to be Chr8q-downregulated proteins (p = ")
  titleSig <- ifelse(dir=="positive"," interactors tend not to be Chr8q-upregulated proteins (p = ",
                     " interactors tend not to be Chr8q-downregulated proteins (p = ")
  titleInsig <- ifelse(dir=="positive"," interactors are equally likely to be Chr8q-upregulated proteins (p = ",
                       " interactors are equally likely to be Chr8q-downregulated proteins (p = ")
  dir.rowname <- ifelse(dir=="positive", "Up with Chr8q", "Down with Chr8q")
  dir.mosaics <- list()
  dir.matrices <- list()
  dir.fisher <- data.frame(terminal, Direction = dir, N_interactor_terminals = NA, 
                           N_noninteractor_terminals = NA, 
                           N_interactor_nonterminals = NA, 
                           N_noninteractor_nonterminals = NA, p = NA)
  for (node in terminal) {
    temp.con <- con[con != node]
    
    # how many connections does node have in STRING network?
    node.con <- STRINGv12[STRINGv12$from == node | STRINGv12$to == node,]
    node.con <- unique(c(node.con$from, node.con$to))
    node.con <- node.con[node.con != node] # 2023 
    
    # how many of those are quantified
    if (dir == "positive") {
      dir.quant.con <- quant.con[!(quant.con %in% sig.quant.neg.con)] # exclude negative terminals from positive analysis
    } else {
      dir.quant.con <- quant.con[!(quant.con %in% sig.quant.pos.con)] # exclude positive terminals from negative analysis
    }
    quant.node.con <- node.con[node.con %in% dir.quant.con] # 1711
    
    # how many are up- or down-regulated?
    sig.quant.node.con <- quant.node.con[quant.node.con %in% sig.quant]
    if (dir == "positive") {
      sig.quant.dir.node.con <- quant.node.con[quant.node.con %in% sig.pos.quant] 
      sig.quant.dir.con <- sig.quant.pos.con
    } else {
      sig.quant.dir.node.con <- quant.node.con[quant.node.con %in% sig.neg.quant]
      sig.quant.dir.con <- sig.quant.neg.con
    }
    
    # prepare contingency matrix
    non.node.quant.con <- dir.quant.con[!(dir.quant.con %in% quant.node.con)]
    sig.dir.non.node.quant.con <- sig.quant.dir.con[!(sig.quant.dir.con %in% quant.node.con)]
    contingency.matrix <- data.frame("node_interaction" = c(length(sig.quant.dir.node.con),
                                                            length(quant.node.con) - length(sig.quant.dir.node.con)),
                                     "Non_node_interaction" = c(length(sig.dir.non.node.quant.con),
                                                                length(non.node.quant.con) - length(sig.dir.non.node.quant.con)),
                                     row.names = c(dir.rowname,"Unaltered with Chr8q"))
    colnames(contingency.matrix) <- paste0(node, c(" Interactor", " Non-interactor"))
    if (length(dir.quant.con) == sum(c(contingency.matrix[,1],contingency.matrix[,2]))) {
      fisher.p <- fisher.test(contingency.matrix)$p.value
      title <- ifelse(fisher.p <= 0.05 & contingency.matrix[1,1] > contingency.matrix[2,1], 
                      paste0(node, titleSigInt, fisher.p, ")"), # if tends to be chr8q terminal interactor
                      ifelse(fisher.p <= 0.05, paste0(node, titleSig, fisher.p, ")"), # if tends to NOT be chr8q terminal interactor
                             paste0(node, titleInsig, fisher.p, ")"))) # if not significant difference
      dir.mosaics[[node]] <- mosaicplot(contingency.matrix, main = title, color=TRUE)
      dir.matrices[[node]] <- contingency.matrix
      dir.fisher[dir.fisher$terminal==node,]$N_interactor_terminals <- contingency.matrix[1,1]
      dir.fisher[dir.fisher$terminal==node,]$N_noninteractor_terminals <- contingency.matrix[1,2]
      dir.fisher[dir.fisher$terminal==node,]$N_interactor_nonterminals <- contingency.matrix[2,1]
      dir.fisher[dir.fisher$terminal==node,]$N_noninteractor_nonterminals <- contingency.matrix[2,2]
      dir.fisher[dir.fisher$terminal==node,]$p <- fisher.p
    } else {
      stop("number of genes in matrix does not match number quantified in network")
    }
  }
  mosaics[[dir]] <- dir.mosaics
  contingency.matrices[[dir]] <- dir.matrices
  fisher.df <- rbind(fisher.df, dir.fisher)
}
fisher.df$Terminal_interactor_ratio <- fisher.df$N_interactor_terminals/(fisher.df$N_noninteractor_terminals + fisher.df$N_interactor_terminals)
fisher.df$Nonterminal_interactor_ratio <- fisher.df$N_interactor_nonterminals/(fisher.df$N_noninteractor_nonterminals + fisher.df$N_interactor_nonterminals)
fisher.df$minusLogP <- -log10(fisher.df$p)
fisher.df$chr8q <- FALSE
fisher.df[fisher.df$terminal %in% chr8q.genes,]$chr8q <- TRUE
fisher.df$drugTarget <- FALSE
drug.info <- read.csv("https://raw.githubusercontent.com/BelindaBGarana/DMEA/refs/heads/shiny-app/Inputs/PRISM_secondary-screen-replicate-treatment-info.csv")
drug.targets <- unique(unlist(strsplit(unlist(drug.info$target), ", "))) # 1027
fisher.df[fisher.df$terminal %in% drug.targets,]$drugTarget <- TRUE
write.csv(fisher.df, "FishersExactTest_terminalSTRINGNodes.csv",row.names=FALSE)
if (any(fisher.df$p > 1)) {
  fisher.df[fisher.df$p>1,]$p <- 1
}
fisher.df$q <- NA
fisher.df$q <- qvalue::qvalue(fisher.df$p,pi0=1)$qvalues
fisher.df$chr8q_q <- NA
fisher.df[fisher.df$chr8q==TRUE,]$chr8q_q <- qvalue::qvalue(fisher.df[fisher.df$chr8q==TRUE,]$p,pi0=1)$qvalues
fisher.df$drugTarget_q <- NA
fisher.df[fisher.df$drugTarget==TRUE,]$drugTarget_q <- qvalue::qvalue(fisher.df[fisher.df$drugTarget==TRUE,]$p,pi0=1)$qvalues
fisher.df$Terminal_minus_nonterminal_interactor_ratio <- fisher.df$Terminal_interactor_ratio - fisher.df$Nonterminal_interactor_ratio
fisher.df$minusLogP_signed <- fisher.df$minusLogP * sign(fisher.df$Terminal_minus_nonterminal_interactor_ratio)
write.csv(fisher.df, "FishersExactTest_terminalSTRINGNodes_adjustedP.csv",row.names=FALSE)
saveRDS(contingency.matrices, "contingencyMatrices_terminalSTRINGNodes.rds")
saveRDS(mosaics, "mosaics_terminalSTRINGNodes.rds")

dot.df <- fisher.df
dot.df$Score <- dot.df$Terminal_interactor_ratio
dot.df[dot.df$Direction=="negative",]$Score <- -dot.df[dot.df$Direction=="negative",]$Score
dot.df$chr8q <- factor(dot.df$chr8q, levels=c(TRUE, FALSE))
dot.df$drugTarget <- factor(dot.df$drugTarget, levels=c(TRUE, FALSE))

n.sig <- nrow(dot.df[dot.df$q<=0.05 & dot.df$Terminal_minus_nonterminal_interactor_ratio>0,]) # 79
n.total <- nrow(dot.df) # 8834
perc.sig <- round(100*n.sig/n.total,2) # 0.52%
title <- paste0(n.sig," / ",n.total, " (", perc.sig, "%) Proteins More Likely to\nInteract with Chr8q-affected Proteins")
ggplot2::ggplot(dot.df, aes(x=Score, y = minusLogP_signed, color=chr8q, shape = drugTarget, size=N_interactor_terminals)) +
  geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  geom_hline(yintercept=-log10(0.05),color="lightgrey", linetype="dashed") +
  geom_hline(yintercept=log10(0.05),color="lightgrey", linetype="dashed") +
  geom_point() + theme_minimal(base_size=12) + scale_color_manual(values=c("red", "darkgrey"), 
                                                                  breaks=levels(dot.df$chr8q)) + scale_shape_manual(values=c(17,19)) +
  ggrepel::geom_text_repel(data=subset(dot.df, p <= 0.05), aes(label=terminal), show.legend=FALSE) + 
  labs(x="Fraction of Interactions with Chr8q-affected Proteins",y=expression(paste("-log"[10]," (P-value)")),
       size = "# of Interactions", color="Chr8q Gene", shape = "Drug Target") +
  geom_point(data = subset(dot.df, drugTarget == FALSE & q <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 21) +
  geom_point(data = subset(dot.df, drugTarget == TRUE & q <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 24) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
ggsave("terminalSTRINGNodes_volcanoPlot_q_v5.pdf",width=10,height=8)

# went back and revised inferred nodes
term.fisher <- read.csv("FishersExactTest_terminalSTRINGNodes_adjustedP.csv")
term.fisher$type <- "Terminal"
colnames(term.fisher)[1] <- "Gene"
fisher.df$type <- "Inferred"
colnames(fisher.df)[1] <- "Gene"
all.fisher <- rbind(term.fisher, fisher.df)

# actually need to add drug corr and MOA results to terminals
synapser::synLogin()
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
sig.pos.corr.targets <- na.omit(unique(unlist(strsplit(drug.corr.wInfo[drug.corr.wInfo$Pearson.q <= 0.05 & drug.corr.wInfo$Pearson.est>0,]$target, ", ")))) # 88
term.fisher$drugCorrPos <- FALSE
term.fisher[term.fisher$Gene %in% sig.pos.corr.targets,]$drugCorrPos <- TRUE
sig.neg.corr.targets <- na.omit(unique(unlist(strsplit(drug.corr.wInfo[drug.corr.wInfo$Pearson.q <= 0.05 & drug.corr.wInfo$Pearson.est>0,]$target, ", ")))) # 88
term.fisher$drugCorrNeg <- FALSE
term.fisher[term.fisher$Gene %in% sig.neg.corr.targets,]$drugCorrNeg <- TRUE

posMOAs <- na.omit(unique(gsea.rna.prot[gsea.rna.prot$p_value <= 0.05 & gsea.rna.prot$FDR_q_value <= 0.25 &
                                          gsea.rna.prot$NES > 0,]$Drug_set))
negMOAs <- na.omit(unique(gsea.rna.prot[gsea.rna.prot$p_value <= 0.05 & gsea.rna.prot$FDR_q_value <= 0.25 &
                                          gsea.rna.prot$NES < 0,]$Drug_set))
"VEGFR inhibitor" %in% c(posMOAs,negMOAs) # FALSE: make sure we don't get VEGFRi when looking for EGFRi
posMOAtargets <- c()
negMOAtargets <- c()
for (i in 1:nrow(red.drug.info)) {
  tempMOAs <- na.omit(unique(strsplit(red.drug.info$moa[i], ", ")))
  tempTargets <- na.omit(unique(strsplit(red.drug.info$target[i], ", ")))
  if (any(tempMOAs %in% posMOAs)) {
    posMOAtargets <- c(posMOAtargets, red.drug.info$target[i])
  }
  if (any(tempMOAs %in% negMOAs)) {
    negMOAtargets <- c(negMOAtargets, red.drug.info$target[i])
  }
}
term.fisher$drugMOATargetPos <- FALSE
term.fisher[term.fisher$Gene %in% posMOAtargets,]$drugMOATargetPos <- TRUE
term.fisher$drugMOATargetNeg <- FALSE
term.fisher[term.fisher$Gene %in% negMOAtargets,]$drugMOATargetNeg <- TRUE
term.fisher$drugCorr <- FALSE
term.fisher[term.fisher$drugCorrNeg | term.fisher$drugCorrPos,]$drugCorr <- TRUE
write.csv(term.fisher, "FishersExactTest_terminalSTRINGNodes_adjustedP_wDrugInfo.csv",row.names=FALSE)

all.fisher <- rbind(term.fisher[,colnames(fisher.df)], fisher.df)
all.fisher$ranking <- all.fisher$Terminal_interactor_ratio * all.fisher$minusLogP_signed
all.fisher$q <- NA
all.fisher$q <- qvalue::qvalue(all.fisher$p,pi0=1)$qvalues
write.csv(all.fisher, "FishersExactTest_STRINGNodes_adjustedP_wDrugInfo.csv")
nrow(all.fisher[all.fisher$chr8q == TRUE & all.fisher$drugTarget == TRUE,]) # 17
both.nets <- unique(all.fisher$Gene[all.fisher$Gene %in% all.fisher[all.fisher$Direction=="positive",]$Gene &
                               all.fisher$Gene %in% all.fisher[all.fisher$Direction=="negative",]$Gene]) # 4423
target.both.nets <- unique(all.fisher[all.fisher$drugTarget == TRUE,]$Gene[all.fisher[all.fisher$drugTarget == TRUE,]$Gene %in% both.nets]) # 354
target.both.nets.info <- all.fisher[all.fisher$Gene %in% target.both.nets,]

chr8q.both.nets <- unique(all.fisher[all.fisher$chr8q == TRUE,]$Gene[all.fisher[all.fisher$chr8q == TRUE,]$Gene %in% both.nets]) # 102
chr8q.both.nets.info <- all.fisher[all.fisher$Gene %in% chr8q.both.nets,]

all.fisher$drugMOATarget <- FALSE
all.fisher[all.fisher$drugMOATargetNeg == TRUE | all.fisher$drugMOATargetPos == TRUE,]$drugMOATarget <- TRUE

dot.df <- all.fisher
dot.df$Score <- dot.df$Terminal_interactor_ratio
dot.df[dot.df$Direction=="negative",]$Score <- -dot.df[dot.df$Direction=="negative",]$Score
dot.df$chr8q <- factor(dot.df$chr8q, levels=c(TRUE, FALSE))
dot.df$drugTarget <- factor(dot.df$drugTarget, levels=c(TRUE, FALSE))

n.sig <- nrow(dot.df[dot.df$q<=0.05 & dot.df$Terminal_minus_nonterminal_interactor_ratio>0,]) # 79
n.total <- nrow(dot.df) # 8834
perc.sig <- round(100*n.sig/n.total,2) # 0.52%
title <- paste0(n.sig," / ",n.total, " (", perc.sig, "%) Proteins More Likely to\nInteract with Chr8q-affected Proteins")
ggplot2::ggplot(dot.df, aes(x=Score, y = minusLogP_signed, color=chr8q, shape = drugTarget, size=N_interactor_terminals)) +
  geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  geom_hline(yintercept=-log10(0.05),color="lightgrey", linetype="dashed") +
  geom_hline(yintercept=log10(0.05),color="lightgrey", linetype="dashed") +
  geom_point() + theme_minimal(base_size=12) + scale_color_manual(values=c("red", "darkgrey"), 
                                                                  breaks=levels(dot.df$chr8q)) + scale_shape_manual(values=c(17,19)) +
  ggrepel::geom_text_repel(data=subset(dot.df, q <= 0.05), aes(label=Gene), show.legend=FALSE) + 
  labs(x="Fraction of Interactions with Chr8q-affected Proteins",y=expression(paste("-log"[10]," (P-value)")),
       size = "# of Interactions", color="Chr8q Gene", shape = "Drug Target") +
  geom_point(data = subset(dot.df, drugTarget == FALSE & q <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 21) +
  geom_point(data = subset(dot.df, drugTarget == TRUE & q <= 0.05), aes(size=N_interactor_terminals), col = "black", shape = 24) +
  ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
ggsave("STRINGNodes_volcanoPlot_q_v5.pdf",width=10,height=8)


# 
# maxLogP <- max(fisher.df$minusLogP)
# maxRatio <- max(fisher.df$Terminal_interactor_ratio)
# directions <- c("positive","negative")
# filter <- c("none", "chr8q", "drugTarget")
# fname <- "FishersExactTest_dotPlots_top5Ratio"
# for (f in filter) {
#   dot.plots <- NULL
#   for (dir in directions) {
#     dot.df <- fisher.df[fisher.df$Direction==dir & fisher.df$p <= 0.05,] 
#     if (f == "chr8q") {
#       dot.df <- dot.df[dot.df$chr8q,]
#       fname2 <- paste0(fname,"_",f)
#     } else if (f == "drugTarget") {
#       dot.df <- dot.df[dot.df$drugTarget,]
#       fname2 <- paste0(fname,"_",f)
#     } else {
#       fname2 <- fname
#     }
#     dot.df <- dot.df %>% slice_max(Terminal_interactor_ratio,n=5)
#     dot.plot <- ggplot2::ggplot(dot.df, aes(x=dir,y=reorder(inferred, -Terminal_interactor_ratio), 
#                                             fill=100*Terminal_interactor_ratio, size=minusLogP)) +
#       geom_point() + theme_classic(base_size=12) + scale_size_continuous(limits=c(0,maxLogP)) +
#       scale_fill_continuous(limits=c(0,maxRatio)) + ggtitle(dir) + 
#       labs(fill="% Interactors", size="-log(p-value)") +
#       geom_point(data = dot.df, col = "black", stroke = 1.5, shape = 21) +
#       theme(axis.title=element_blank(), plot.title=element_text(face="bold", hjust=0.5))
#     if (is.null(dot.plots)) {
#       dot.plots <- dot.plot
#     } else {
#       dot.plots <- dot.plots + dot.plot + plot_layout(guides="collect")
#     }
#   }
#   dot.plots
#   ggsave(paste0(fname2,".pdf"), width=5, height=5)
# }

# 
# # narrow it down by filtering for nodes with most interactions
# posN <- plyr::ddply(pos.edges, .(from), summarize,
#                     N = n()) # only need 'from' since edges are symmetrical
# colnames(pos.vert)[1] <- "from"
# posN <- merge(posN, pos.vert, by="from")
# hist(posN$N)
# hist(posN[posN$N < 50,]$N)
# hist(posN[posN$N < 10,]$N)
# quantile(posN$N)
# # 0%  25%  50%  75% 100% 
# # 1    1    2    4  897 
# nrow(posN[posN$N > 4,]) # 2057
# pos500 <- posN %>% slice_max(N, n=500) # 525
# min(pos500$N) # 13
# pos1000 <- posN %>% slice_max(N, n=1000) # 1024
# min(pos1000$N) # 8
# pos1000edges <- pos.edges[pos.edges$from %in% pos1000$from,] # fine because symmetrical
# pos1000vert <- pos.vert[pos.vert$from %in% c(pos1000edges$from, pos1000edges$to),] # back to 9133
# 
# negN <- plyr::ddply(neg.edges, .(from), summarize,
#                     N = n()) # only need 'from' since edges are symmetrical
# colnames(neg.vert)[1] <- "from"
# negN <- merge(negN, neg.vert, by="from")
# hist(negN$N)
# hist(negN[negN$N < 50,]$N)
# hist(negN[negN$N < 10,]$N)
# quantile(negN$N)
# # 0%  25%  50%  75% 100% 
# # 1    1    1    2  592 
# nrow(negN[negN$N > 2,]) # 1544
# neg500 <- negN %>% slice_max(N, n=500) # 605
# min(neg500$N) # 5
# neg1000 <- negN %>% slice_max(N, n=1000) # 1544
# min(neg1000$N) # 3
# 
# # # make graphs
# # topGraph <- igraph::graph_from_data_frame(pos.edges, directed=FALSE, vertices=posN) 
# # plot(topGraph)
# # tempTitle <- paste0("positive", "_", nrow(pos1000),"_topConnectedNodes_manual_",Sys.Date())
# # RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
# # 
# # topGraph <- igraph::graph_from_data_frame(neg.edges, directed=FALSE, vertices=negN) 
# # plot(topGraph)
# # tempTitle <- paste0("negative", "_", nrow(neg1000),"_topConnectedNodes_manual_",Sys.Date())
# # RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
# 
# # try filtering for nodes which connect terminals
# # pos.con <- pos.edges[pos.edges$from %in% pos.terms & pos.edges$to %in% pos.terms,] # 752
# # pos.con.vert <- posN[posN$from %in% c(pos.con$from, pos.con$to),] # 206 - all terminal
# 
# pos.con.vert <- posN[posN$N > 1,] # 5630
# pos.con.vert$chr8q <- FALSE
# pos.con.vert[pos.con.vert$from %in% chr8q.genes,]$chr8q <- TRUE
# length(pos.terms[pos.terms %in% pos.con.vert$from]) * 100 / length(pos.terms) # 97.71242% of terminals included
# pos.con <- pos.edges[pos.edges$from %in% pos.con.vert$from & pos.edges$to %in% pos.con.vert$from,] # 51200
# write.csv(pos.con.vert, "positive_1PlusConnection_vertices.csv", row.names=FALSE)
# write.csv(pos.con, "positive_1PlusConnection_edges.csv", row.names=FALSE)
# topGraph <- igraph::graph_from_data_frame(pos.con, directed=FALSE, vertices=pos.con.vert) 
# #plot(topGraph)
# tempTitle <- paste0("positive", "_", nrow(pos.con.vert),"_topConnectedNodes_manual_1PlusConnection_",Sys.Date())
# RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
# 
# pos.con.centrality <- data.frame(name = V(topGraph)$name,
#                                    degree = igraph::degree(topGraph, mode="all"),
#                                    closeness = igraph::closeness(topGraph, mode="all"),
#                                    betweenness = igraph::betweenness(topGraph, directed = FALSE),
#                                    eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
#                                    hub_score = igraph::hub_score(topGraph)$vector,
#                                    authority_score = igraph::authority_score(topGraph)$vector)
# write.csv(pos.con.centrality, "positive_1PlusConnection_centrality.csv", row.names=FALSE)
# 
# pos.tf.con <- pos.con[(pos.con$from %in% tf.result$Gene | pos.con$to %in% tf.result$Gene) &
#                         (pos.con$from %in% chr8q.genes | pos.con$to %in% chr8q.genes),] # 5050
# pos.tf.con.vert <- pos.con.vert[pos.con.vert$from %in% c(pos.tf.con$from, pos.tf.con$to),] # 815
# topGraph <- igraph::graph_from_data_frame(pos.tf.con, directed=FALSE, vertices=pos.tf.con.vert) 
# #plot(topGraph)
# tempTitle <- paste0("positive", "_", nrow(pos.tf.con.vert),"_connectedTFChr8Nodes_manual_1PlusConnection_",Sys.Date())
# RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
# write.csv(pos.tf.con.vert, "positive_TFChr8_1PlusConnection_vertices.csv", row.names=FALSE)
# write.csv(pos.tf.con, "positive_TFChr8_1PlusConnection_edges.csv", row.names=FALSE)
# pos.tf.con.centrality <- data.frame(name = V(topGraph)$name,
#                                  degree = igraph::degree(topGraph, mode="all"),
#                                  closeness = igraph::closeness(topGraph, mode="all"),
#                                  betweenness = igraph::betweenness(topGraph, directed = FALSE),
#                                  eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
#                                  hub_score = igraph::hub_score(topGraph)$vector,
#                                  authority_score = igraph::authority_score(topGraph)$vector)
# write.csv(pos.tf.con.centrality, "positive_TFChr8_1PlusConnection_centrality.csv", row.names=FALSE)
# 
# neg.con.vert <- negN[negN$N > 1,] # 2890
# neg.con.vert$chr8q <- FALSE
# neg.con.vert[neg.con.vert$from %in% chr8q.genes,]$chr8q <- TRUE
# length(neg.terms[neg.terms %in% neg.con.vert$from]) * 100 / length(neg.terms) # 99.09091% of terminals included
# neg.con <- neg.edges[neg.edges$from %in% neg.con.vert$from & neg.edges$to %in% neg.con.vert$from,] # 19422
# write.csv(neg.con.vert, "negative_1PlusConnection_vertices.csv", row.names=FALSE)
# write.csv(neg.con, "negative_1PlusConnection_edges.csv", row.names=FALSE)
# topGraph <- igraph::graph_from_data_frame(neg.con, directed=FALSE, vertices=neg.con.vert) 
# #plot(topGraph)
# tempTitle <- paste0("negative", "_", nrow(neg.con.vert),"_topConnectedNodes_manual_1PlusConnection_",Sys.Date())
# RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
# 
# neg.con.centrality <- data.frame(name = V(topGraph)$name,
#                                  degree = igraph::degree(topGraph, mode="all"),
#                                  closeness = igraph::closeness(topGraph, mode="all"),
#                                  betweenness = igraph::betweenness(topGraph, directed = FALSE),
#                                  eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
#                                  hub_score = igraph::hub_score(topGraph)$vector,
#                                  authority_score = igraph::authority_score(topGraph)$vector)
# write.csv(neg.con.centrality, "negative_1PlusConnection_centrality.csv", row.names=FALSE)
# 
# neg.tf.con <- neg.con[(neg.con$from %in% c(tf.result$Gene, kin.result$Feature_set) | 
#                          neg.con$to %in% c(tf.result$Gene, kin.result$Feature_set)) &
#                         (neg.con$from %in% chr8q.genes | neg.con$to %in% chr8q.genes),] # 5050
# neg.tf.con.vert <- neg.con.vert[neg.con.vert$from %in% c(neg.tf.con$from, neg.tf.con$to),] # 815
# topGraph <- igraph::graph_from_data_frame(neg.tf.con, directed=FALSE, vertices=neg.tf.con.vert) 
# #plot(topGraph)
# tempTitle <- paste0("negative", "_", nrow(neg.tf.con.vert),"_connectedTFChr8Nodes_manual_1PlusConnection_",Sys.Date())
# RCy3::createNetworkFromIgraph(topGraph, title=tempTitle)
# write.csv(neg.tf.con.vert, "negative_TFChr8_1PlusConnection_vertices.csv", row.names=FALSE)
# write.csv(neg.tf.con, "negative_TFChr8_1PlusConnection_edges.csv", row.names=FALSE)
# neg.tf.con.centrality <- data.frame(name = V(topGraph)$name,
#                                     degree = igraph::degree(topGraph, mode="all"),
#                                     closeness = igraph::closeness(topGraph, mode="all"),
#                                     betweenness = igraph::betweenness(topGraph, directed = FALSE),
#                                     eigen_centrality = igraph::eigen_centrality(topGraph, directed = FALSE)$vector,
#                                     hub_score = igraph::hub_score(topGraph)$vector,
#                                     authority_score = igraph::authority_score(topGraph)$vector)
# write.csv(neg.tf.con.centrality, "negative_TFChr8_1PlusConnection_centrality.csv", row.names=FALSE)

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