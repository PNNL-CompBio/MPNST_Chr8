# STRING protein networks
# full: all correlated proteins in addition to enriched TFs and kinases
# regulatory: enriched TFs and kinases
remove(list=ls())
library(synapser); library(msigdbr); library(PCSF)
library(plyr); library(ggplot2)
setwd("~/OneDrive - PNNL/Documents/GitHub/Chr8/proteomics/analysis/Chr8_quant_20250409")

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

# get chr8q genes
chr8q.genes <- msigdbr::msigdbr(collection="C1")
chr8q.genes <- unique(chr8q.genes[startsWith(chr8q.genes$gs_name,"chr8q"),]$gene_symbol)

# get latest STRING data (physical protein interactions)
library(PCSF)
data("STRINGv12") # 1477610

#### full networks ####
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
colnames(pos.vert)[1] <- "name"
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
pos.centrality$chr8q <- FALSE
pos.centrality[pos.centrality$name %in% chr8q.genes,]$chr8q <- TRUE
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
colnames(neg.vert)[1] <- "name"
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
neg.centrality$chr8q <- FALSE
neg.centrality[neg.centrality$name %in% chr8q.genes,]$chr8q <- TRUE
write.csv(neg.centrality, "negative_centrality.csv", row.names=FALSE)
neg.vert <- read.csv("negative_vertices.csv")
neg.edges <- read.csv("negative_edges.csv")
neg.centrality <- read.csv("negative_centrality.csv")

# look for connections between MYC, PLK1/EGFR, and other proteins in full network
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


#### regulatory networks (TFs and kinases) #### 
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
colnames(pos.reg.vert)[1] <- "name"
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
colnames(neg.reg.vert)[1] <- "name"
write.csv(neg.reg.vert, "negative_TFKin_vertices.csv", row.names=FALSE)
write.csv(neg.reg.edges, "negative_TFKin_edges.csv", row.names=FALSE)
neg.reg.centrality <- read.csv("negative_TFKin_centrality.csv")

neg.reg.centrality$Direction <- "negative"
pos.reg.centrality$Direction <- "positive"

# look at centrality of chr8q genes in regulatory networks
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
bar.df2 <- chr8q.reg.centrality[chr8q.reg.centrality$name %in% bar.df$name,]
ggplot2::ggplot(bar.df2, aes(x=eigen_centrality, y=name, fill=Direction)) +
  geom_bar(stat="identity", alpha=0.7) + 
  scale_fill_manual(values=c("red","blue"), labels=c("Upregulated\nNetwork","Downregulated\nNetwork"), breaks=c("positive","negative")) +
  theme_classic(base_size=12) + labs(x="Centrality", y=NULL) + 
  theme(legend.position="inside", legend.position.inside = c(0.7,0.5)) +
  scale_y_discrete(limits=bar.df[order(bar.df$meanCentrality),]$name)
ggsave("chr8qGeneCentralityTFKinNetworks_barPlotStacked.pdf",width=3.5,height=3)