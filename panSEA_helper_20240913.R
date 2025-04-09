# Differential expression & enrichment analyses: global & phospho
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-04-25

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr)
library(dplyr); library(pheatmap); library(grid)

# get gmt information for GSEA relevant to Chr8
get_chr8_gmt1 <- function(gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                        "msigdb_Homo sapiens_H",
                                        "msigdb_Homo sapiens_C1",
                                        "chr8_cancer_human"), min.per.set=6) {
  if (file.exists("chr8_gmt1_run_contrasts_global_phospho_human.rds")) {
    gmt1 <- readRDS("chr8_gmt1_run_contrasts_global_phospho_human.rds")
  } else {
    gmt1 <- list()
    for (i in 1:length(gmt.list1)) {
      if (is.character(gmt.list1[i])) {
        if (grepl("msigdb", gmt.list1[i], ignore.case = TRUE)) {
          gmt.info <- stringr::str_split(gmt.list1[i], "_")[[1]]
          if (length(gmt.info) > 1) {
            if (length(gmt.info) == 2) {
              msigdb.info <- msigdbr::msigdbr(gmt.info[2])
            } else if (length(gmt.info) == 3) {
              msigdb.info <- msigdbr::msigdbr(gmt.info[2], gmt.info[3])
            } else {
              msigdb.info <- msigdbr::msigdbr(gmt.info[2], gmt.info[3], gmt.info[4])
            }
            
            # extract necessary info into data frame
            msigdb.info <- as.data.frame(msigdb.info[, c(
              "gene_symbol",
              "gs_name",
              "gs_description"
            )])
            
            gmt1[[i]] <- DMEA::as_gmt(
              msigdb.info, "gene_symbol", "gs_name", min.per.set,
              descriptions = "gs_description"
            )
          }
        } else if (gmt.list1[i] == "chr8_cancer_human") {
          msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C1")
          msigdb.info <- as.data.frame(msigdb.info[, c(
            "gene_symbol",
            "gs_name",
            "gs_description"
          )])
          
          gmt <- DMEA::as_gmt(
            msigdb.info, "gene_symbol", "gs_name", min.per.set,
            descriptions = "gs_description"
          )
          Chr8.cancer.genes <- c("PLAG1", "CHCHD7", "SOX17", "TCEA1", "NCOA2", "TCEB1", 
                                 "HEY1", "RUNX1T1", "NBN", "CNBD1", "COX6C", "UBR5", 
                                 "RAD21", "EXT1", "MYC", "NDRG1", "EPPK1", "RECQL4")
          gmt$genesets[[length(gmt$genesets)+1]] <- Chr8.cancer.genes
          gmt$geneset.names[[length(gmt$geneset.names)+1]] <- "Chr8 cancer-associated genes"
          gmt$geneset.descriptions[[length(gmt$geneset.descriptions)+1]] <- "Chr8 cancer-associated genes"
          gmt1[[i]] <- gmt
        }
      } else {
        gmt1[[i]] <- gmt.list1[[i]]
      }
    } 
    saveRDS(gmt1, "gmt1.rds")
  }
  return(gmt1)
}

# get gmt information for KSEA/SSEA relevant to chr8
get_chr8_gmt2 <- function() {
  if (file.exists("chr8_gmt2_run_contrasts_global_phospho_human.rds")) {
    gmt2 <- readRDS("chr8_gmt2_run_contrasts_global_phospho_human.rds")
  } else {
    stop("Chr8 gmt2 file not found in current directory") 
  }
  saveRDS(gmt2, "gmt2.rds")
}

save_base_plot <- function(base.plot, filename, width = 7, height = 7) {
  pdf(filename, width, height)
  grid::grid.newpage()
  grid::grid.draw(base.plot$gtable)
  dev.off()
  return()
}

# check if data frame is valid for hclust
filter_for_hclust <- function(expr.mat) {
  if (nrow(expr.mat) > 1) { 
    validRows <- rowSums(is.na(expr.mat)) < (ncol(expr.mat)-1)
    validCols <- colSums(is.na(expr.mat)) < (nrow(expr.mat)-1)
    keep <- apply(dplyr::select_if(expr.mat, is.numeric), 1, function(x) length(unique(x[!is.na(x)])) != 1)
    while (sum(keep, na.rm = TRUE) < nrow(expr.mat) | 
           sum(validRows, na.rm = TRUE) < nrow(expr.mat) |
           sum(validCols, na.rm = TRUE) < ncol(expr.mat)) {
      # remove invalid rows/columns
      expr.mat <- expr.mat[validRows,] # require 2+ numeric values along rows
      expr.mat <- expr.mat[,validCols] # require 2+ numeric values along columns 
      if(is.data.frame(expr.mat) & length(keep) > 0){
        expr.mat <- expr.mat[keep,] # remove rows where all values are the same
        
        # re-evaluate if there are any invalid rows/columns
        validRows <- rowSums(is.na(expr.mat)) < (ncol(expr.mat)-1)
        validCols <- colSums(is.na(expr.mat)) < (nrow(expr.mat)-1)
        keep <- apply(dplyr::select_if(expr.mat, is.numeric), 1, function(x) length(unique(x[!is.na(x)])) != 1)
      } else {
        expr.mat <- data.frame()
        return(expr.mat)
      } 
    }
    expr.mat <- as.data.frame(expr.mat[which(rowMeans(!is.na(expr.mat)) >= 0.5),]) # require at least 50% coverage for each feature
  }
  return(expr.mat)
}

## get expression of features from pathways of interest and heatmaps
get_pathways_of_interest <- function(expr.df, gsea.result, gmt, cc.df, n=5, 
                                     FDR = 0.25, p = 0.05, show_colnames = FALSE, 
                                     fontsize = 10, scale = TRUE, cluster = TRUE) {
  # get top pathways of interest
  sig.result <- gsea.result[gsea.result$FDR_q_value <= FDR & 
                              gsea.result$p_value <= p,]
  if (nrow(sig.result) > 0) {
    top.pathways <- sig.result %>% slice_max(abs(NES), n = n)
    top.gmt <- gmt$genesets[top.pathways$Feature_set]
    
    poi.files <- try(R.utils::withTimeout(make_heatmaps(expr.df, cc.df, top.pathways, top.gmt, show_colnames,
                               fontsize, scale, cluster), timeout = 600, onTimeout="error"), silent = TRUE)
    if (inherits(poi.files, "try-error")) {
      poi.files <- list()
    }
  } else {
    poi.files <- list()
  }
  return(poi.files)
}

make_heatmaps <- function(expr.df, cc.df, top.pathways = NULL, top.gmt, show_colnames = FALSE,
                          fontsize = 10, scale = TRUE, cluster = TRUE) {
  expr.list <- list()
  lead.list <- list()
  
  my.clust.heatmaps <- list()
  my.abs.heatmaps <- list()
  
  my.clust.heatmaps.leads <- list()
  my.abs.heatmaps.leads <- list()
  for (i in 1:length(top.gmt)) {
    # create heatmap names
    clust.name <- file.path(paste0("Expression_of_", names(top.gmt)[i], "_heatmap_scaled.bp"))
    abs.name <- file.path(paste0("Expression_of_", names(top.gmt)[i], "_heatmap_not_scaled.bp"))
    clust.name.leads <- file.path(paste0("Expression_of_", names(top.gmt)[i], "_leading_edge_heatmap_scaled.bp"))
    abs.name.leads <- file.path(paste0("Expression_of_", names(top.gmt)[i], "_leading_edge_heatmap_not_scaled.bp"))
    
    # select numeric data
    just.expr <- dplyr::select_if(expr.df, is.numeric)
    
    # identify feature name
    feature.name <- colnames(expr.df)[!(colnames(expr.df) %in% colnames(just.expr))]
    
    # filter for gene set
    temp.expr <- expr.df[expr.df[,feature.name] %in% top.gmt[[i]], ]
    rownames(temp.expr) <- temp.expr[,feature.name]
    
    expr.list[[names(top.gmt)[i]]] <- temp.expr
    expr.mat <- dplyr::select_if(temp.expr, is.numeric)
    expr.mat <- filter_for_hclust(expr.mat)
    
    if (nrow(expr.mat) > 1) {
      # create heatmaps
      expr.mat <- as.matrix(expr.mat) 
      if (scale) {
        my.clust.heatmaps[[clust.name]] <- pheatmap::pheatmap(expr.mat, color = 
                                                                colorRampPalette(
                                                                  c("navy", "white", "firebrick3"))(50),
                                                              cluster_rows = cluster, cluster_cols = cluster,
                                                              scale = "row", annotation_col = cc.df, 
                                                              angle_col = "45", 
                                                              show_colnames = show_colnames,
                                                              fontsize = fontsize) 
      }
      my.abs.heatmaps[[abs.name]] <- pheatmap::pheatmap(expr.mat, color = 
                                                          colorRampPalette(
                                                            c("navy", "white", "firebrick3"))(50), 
                                                        cluster_rows = cluster, cluster_cols = cluster,
                                                        annotation_col = cc.df, 
                                                        angle_col = "45", 
                                                        show_colnames = show_colnames,
                                                        fontsize = fontsize)
      
      # filter to leading edge features
      if (!is.null(top.pathways)) {
        temp.set <- names(top.gmt)[i]
        temp.leads <- stringr::str_split(top.pathways[top.pathways$Feature_set == temp.set,]$Leading_edge, ", ")[[1]]
        temp.expr.leads <- temp.expr[temp.leads,]
        lead.list[[names(top.gmt)[i]]] <- temp.expr.leads
        lead.mat <- dplyr::select_if(temp.expr.leads, is.numeric)
        lead.mat <- filter_for_hclust(lead.mat)
        
        if (nrow(lead.mat) > 1) {
          # create heatmap
          lead.mat <- as.matrix(lead.mat)
          if (scale) {
            my.clust.heatmaps.leads[[clust.name.leads]] <- pheatmap::pheatmap(lead.mat, color = 
                                                                                colorRampPalette(
                                                                                  c("navy", "white", "firebrick3"))(50),
                                                                              cluster_rows = cluster, cluster_cols = cluster,
                                                                              scale = "row", annotation_col = cc.df, 
                                                                              angle_col = "45", 
                                                                              show_colnames = show_colnames,
                                                                              fontsize = fontsize) 
          }
          my.abs.heatmaps.leads[[abs.name.leads]] <- pheatmap::pheatmap(lead.mat, color = 
                                                                          colorRampPalette(
                                                                            c("navy", "white", "firebrick3"))(50), 
                                                                        cluster_rows = cluster, cluster_cols = cluster,
                                                                        annotation_col = cc.df, 
                                                                        angle_col = "45", 
                                                                        show_colnames = show_colnames,
                                                                        fontsize = fontsize) 
        }
      }
    }
  }
  all.expr.df <- data.table::rbindlist(expr.list, use.names = TRUE, idcol = "Feature_set")
  all.leads.df <- data.table::rbindlist(lead.list, use.names = TRUE, idcol = "Feature_set")
  
  lead.heatmaps <- list("Not_scaled" = my.abs.heatmaps.leads,
                        "Scaled" = my.clust.heatmaps.leads)
  lead.files <- list("Expression_of_leading_edge_features.csv" =
                       all.leads.df,
                     "Heatmaps" = lead.heatmaps)
  
  all.heatmaps <- list("Not_scaled" = my.abs.heatmaps,
                       "Scaled" = my.clust.heatmaps,
                       "Leading_edge_features_only" = lead.files)
  poi.files <- list("Expression_of_pathways_of_interest.csv" =
                      all.expr.df,
                    "Heatmaps" = all.heatmaps)
  return(poi.files)
}

## get gmt information for GSEA
get_gmt1 <- function(gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                   "msigdb_Homo sapiens_H",
                                   "msigdb_Homo sapiens_C1"), min.per.set=6) {
  if (file.exists("gmt1.rds")) {
    gmt1 <- readRDS("gmt1.rds")
  } else {
    if ("chr8" %in% gmt.list1) {
      gmt1 <- get_chr8_gmt1(gmt.list1)
    } else {
      gmt1 <- list()
      for (i in 1:length(gmt.list1)) {
        if (is.character(gmt.list1[i])) {
          if (grepl("msigdb", gmt.list1[i], ignore.case = TRUE)) {
            gmt.info <- stringr::str_split(gmt.list1[i], "_")[[1]]
            if (length(gmt.info) > 1) {
              if (length(gmt.info) == 2) {
                msigdb.info <- msigdbr::msigdbr(gmt.info[2])
              } else if (length(gmt.info) == 3) {
                msigdb.info <- msigdbr::msigdbr(gmt.info[2], gmt.info[3])
              } else {
                msigdb.info <- msigdbr::msigdbr(gmt.info[2], gmt.info[3], gmt.info[4])
              }
              
              # extract necessary info into data frame
              msigdb.info <- as.data.frame(msigdb.info[, c(
                "gene_symbol",
                "gs_name",
                "gs_description"
              )])
              
              gmt1[[i]] <- DMEA::as_gmt(
                msigdb.info, "gene_symbol", "gs_name", min.per.set,
                descriptions = "gs_description"
              )
            }
          }
        }else {
          gmt1[[i]] <- gmt.list1[[i]]
        } 
      }
    } 
  }
  saveRDS(gmt1, "gmt1.rds")
  return(gmt1)
}

# uses more pathways
get_gmt1_v2 <- function(gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                   "msigdb_Homo sapiens_H",
                                   "msigdb_Homo sapiens_C1",
                                   "msigdb_Homo sapiens_C3_TFT:GTRD",
                                   "msigdb_Homo sapiens_C3_MIR:MIRDB",
                                   "msigdb_Homo sapiens_C5_GO:BP",
                                   "msigdb_Homo sapiens_C5_GO:CC",
                                   "msigdb_Homo sapiens_C5_GO:MF",
                                   "msigdb_Homo sapiens_C6",
                                   "msigdb_Homo sapiens_C2_CP:BIOCARTA",
                                   "msigdb_Homo sapiens_C2_CP:PID",
                                   "msigdb_Homo sapiens_C2_CP:REACTOME",
                                   "msigdb_Homo sapiens_C2_CP:WIKIPATHWAYS"),
                        names1 = c("KEGG", "Hallmark", "Positional", "TFT_GTRD", 
                                  "MIR_MIRDB", "GO_BP", "GO_CC", "GO_MF", 
                                  "Oncogenic_signatures", "BioCarta", "KEGG", 
                                  "PID", "Reactome", "WikiPathways"),
                        min.per.set=6) {
  if (file.exists("gmt1.rds")) {
    gmt1 <- readRDS("gmt1_more.rds")
  } else {
    if ("chr8" %in% gmt.list1) {
      gmt1 <- get_chr8_gmt1(gmt.list1)
    } else {
      gmt1 <- list()
      for (i in 1:length(gmt.list1)) {
        if (is.character(gmt.list1[i])) {
          if (grepl("msigdb", gmt.list1[i], ignore.case = TRUE)) {
            gmt.info <- stringr::str_split(gmt.list1[i], "_")[[1]]
            if (length(gmt.info) > 1) {
              if (length(gmt.info) == 2) {
                msigdb.info <- msigdbr::msigdbr(gmt.info[2])
              } else if (length(gmt.info) == 3) {
                msigdb.info <- msigdbr::msigdbr(gmt.info[2], gmt.info[3])
              } else {
                msigdb.info <- msigdbr::msigdbr(gmt.info[2], gmt.info[3], gmt.info[4])
              }
              
              # extract necessary info into data frame
              msigdb.info <- as.data.frame(msigdb.info[, c(
                "gene_symbol",
                "gs_name",
                "gs_description"
              )])
              
              gmt1[[i]] <- DMEA::as_gmt(
                msigdb.info, "gene_symbol", "gs_name", min.per.set,
                descriptions = "gs_description"
              )
            }
          }
        }else {
          gmt1[[i]] <- gmt.list1[[i]]
        } 
      }
    } 
  }
  names(gmt1) <- names1
  saveRDS(gmt1, "gmt1_more.rds")
  return(gmt1)
}

# uses feature names like ABC-S325s instead of ABC-S325
get_gmt2 <- function(gmt.list2 = c("ksdb_human", "sub"), phospho) {
  if (file.exists("gmt2.rds")) {
    gmt2 <- readRDS("gmt2.rds")
  } else if ("chr8" %in% gmt.list2) {
    gmt2 <- get_chr8_gmt2(gmt.list2)
  } else {
      gmt2 <- list()
      for (i in 1:length(gmt.list2)) {
        if (is.character(gmt.list2[i])) {
          if (grepl("ksdb", gmt.list2[i], ignore.case = TRUE)) {
            org <- stringr::str_split(gmt.list2[i], "_")[[1]][2]
            gmt2[[i]] <- get_ksdb(organism = org)
          } else if (gmt.list2[i] == "sub") {
            SUB_SITE <- phospho$SUB_SITE
            phospho.ref <- data.frame(SUB_SITE)
            phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
                                                          remove = FALSE)
            SUB_SITE <- NULL
            gmt2[[i]] <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE")
          }
        } else {
          gmt2 <- gmt.list2
          break
        }
      } 
    saveRDS(gmt2, "gmt2.rds") 
  }
  return(gmt2)
}

# uses feature names like ABC-S325 instead of ABC-S325s
get_gmt2_v2 <- function(gmt.list2 = c("ksdb_human", "sub"), phospho) {
  if (file.exists("gmt2_v2.rds")) {
    gmt2 <- readRDS("gmt2_v2.rds")
  } else if ("chr8" %in% gmt.list2) {
    gmt2 <- get_chr8_gmt2(gmt.list2)
  } else {
    gmt2 <- list()
    for (i in 1:length(gmt.list2)) {
      if (is.character(gmt.list2[i])) {
        if (grepl("ksdb", gmt.list2[i], ignore.case = TRUE)) {
          org <- stringr::str_split(gmt.list2[i], "_")[[1]][2]
          gmt2[[i]] <- get_ksdb_v2(organism = org)
        } else if (gmt.list2[i] == "sub") {
          SUB_SITE <- phospho$SUB_SITE
          phospho.ref <- data.frame(SUB_SITE)
          phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
                                                        remove = FALSE)
          SUB_SITE <- NULL
          gmt2[[i]] <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE")
        }
      } else {
        gmt2 <- gmt.list2
        break
      }
    } 
    saveRDS(gmt2, "gmt2_v2.rds") 
  }
  return(gmt2)
}

# get kinase-substrate database information in gmt format for KSEA
get_ksdb <- function(organism="human"){
  if (!is.character(organism)) {
    stop("organism must be character entry")
  }
  
  if (file.exists(file.path(paste0("gmt_ksdb_", organism, ".rds")))) {
    gmt <- readRDS(file.path(paste0("gmt_ksdb_", organism, ".rds")))
  } else if (organism == "human") {
    url <- "https://raw.github.com/BelindaBGarana/panSEA/main/data/gmt_ksdb_human.rds"
    httr::GET(url, httr::write_disk("gmt_ksdb_human.rds")) 
    gmt <- readRDS("gmt_ksdb_human.rds")
  } else {
    ksdb <- read.csv(paste0("https://raw.githubusercontent.com/BelindaBGarana/",
                            "panSEA/shiny-app/data/ksdb_20231101.csv"))
    if (organism %in% unique(na.omit(ksdb$KIN_ORGANISM)) &
        organism %in% unique(na.omit(ksdb$SUB_ORGANISM))) {
      ksdb <- ksdb[ksdb$KIN_ORGANISM == organism &
                     ksdb$SUB_ORGANISM == organism, ]
    }
    ksdb$SUB_MOD_RSD_PNNL <- ksdb$SUB_MOD_RSD
    ksdb$SUB_MOD_RSD_PNNL <- paste0(ksdb$SUB_MOD_RSD_PNNL,
                                    tolower(substr(ksdb$SUB_MOD_RSD_PNNL,1,1)))
    ksdb <- ksdb %>% unite("SUB_SITE", c("SUBSTRATE", "SUB_MOD_RSD_PNNL"),
                           sep = "-", remove = FALSE)
    gmt <- DMEA::as_gmt(ksdb, "SUB_SITE", "KINASE",
                        descriptions = "KIN_ACC_ID")
    saveRDS(gmt, file.path(paste0("ksdb_", organism, ".rds")))
  } 
  return(gmt)
}

get_ksdb_v2 <- function(organism="human"){
  if (!is.character(organism)) {
    stop("organism must be character entry")
  }
  
  if (file.exists(file.path(paste0("gmt_ksdb_", organism, "_v2.rds")))) {
    gmt <- readRDS(file.path(paste0("gmt_ksdb_", organism, "_v2.rds")))
  } else {
    ksdb <- read.csv(paste0("https://raw.githubusercontent.com/BelindaBGarana/",
                            "panSEA/shiny-app/data/ksdb_20231101.csv"))
    if (organism %in% unique(na.omit(ksdb$KIN_ORGANISM)) &
        organism %in% unique(na.omit(ksdb$SUB_ORGANISM))) {
      ksdb <- ksdb[ksdb$KIN_ORGANISM == organism &
                     ksdb$SUB_ORGANISM == organism, ]
    }
    ksdb <- ksdb %>% unite("SUB_SITE", c("SUBSTRATE", "SUB_MOD_RSD"),
                           sep = "-", remove = FALSE)
    gmt <- DMEA::as_gmt(ksdb, "SUB_SITE", "KINASE",
                        descriptions = "KIN_ACC_ID")
    saveRDS(gmt, file.path(paste0("ksdb_", organism, "_v2.rds")))
  } 
  return(gmt)
}

# get CCLE global proteomics data
get_CCLE_prot <- function() {
  message("Loading CCLE proteomics")
  if (file.exists("CCLE_proteomics.csv")) {
    prot.df.noNA <- read.csv("CCLE_proteomics.csv")
  } else {
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
  write.csv(prot.df.noNA, "CCLE_proteomics.csv", row.names = FALSE)
  return(prot.df.noNA)
}

get_CCLE_RNA <- function() {
  message("Loading adherent CCLE RNA-seq data version 19Q4")
  if (!file.exists("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")) {
    download.file(
      paste0(
        "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"
      ),
      destfile =
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"
    )
  }
  
  if (!file.exists("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")) {
    download.file(
      paste0(
        "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"
      ),
      destfile =
        "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"
    )
  }
  
  load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
  load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
  RNA.df <- rbind(RNA.first200, RNA.rest)
  return(RNA.df)
}

## load BeatAML data formatted for DMEA
# input: file path where BeatAML data should be saved
# output: list of BeatAML meta, drug AUC, global, and phospho data frames
load_BeatAML_for_DMEA <- function(BeatAML.path = "BeatAML_DMEA_inputs") {
  message("Loading Beat AML data for DMEA")
  BeatAML_synapse_id <- list("drug_response.csv" = "syn51674470", 
                             "Ex10_metadata.txt" = "syn25807733",
                             "ptrc_ex10_crosstab_global_gene_corrected.txt" = "syn25714248",
                             "ptrc_ex10_crosstab_phospho_siteID_corrected(1).txt" = "syn25714921")
  
  gmt.drug <- readRDS(gzcon(url("https://raw.github.com/BelindaBGarana/panSEA/main/Examples/Inputs/gmt_BeatAML_drug_MOA.rds")))
  
  ### download files if any not already downloaded
  if (!file.exists(BeatAML.path)) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  } else if (!any(FALSE %in% lapply(names(BeatAML_synapse_id), file.exists))) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  }
  
  ### load files
  drug.BeatAML <- read.csv(file.path(BeatAML.path, names(BeatAML_synapse_id)[1]))
  meta.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[2]), 
                             sep = "\t", header = TRUE)
  global.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[3]),
                               sep = "\t", header = TRUE)
  phospho.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[4]),
                                sep = "\t", header = TRUE)
  
  ### format BeatAML data for DMEA
  sample.names <- "Barcode.ID"
  
  ## format drug sensitivity data frame
  # format drug.BeatAML wide (samples in first column, drug names for rest of columns)
  drug.BeatAML <- reshape2::dcast(drug.BeatAML, sample_id ~ inhibitor, 
                                  value.var = "auc", fill = NA)
  
  # # remove drugs without moa annotations and drug combos
  # valid.drugs <- 
  #   names(drug.BeatAML)[names(drug.BeatAML) %in% 
  #                         moa.BeatAML[!is.na(moa.BeatAML),]$Drug] # 167 drugs
  # drug.BeatAML <- drug.BeatAML[ , c("sample_id", valid.drugs)] # 167 drugs
  # moa.BeatAML <- 
  #   moa.BeatAML[moa.BeatAML$Drug %in% names(drug.BeatAML)[2:ncol(drug.BeatAML)], ]
  
  # change sample column name to match expression data
  names(drug.BeatAML)[1] <- sample.names
  
  ## format global proteomics data frame
  # change global.BeatAML column names from SampleID.abbrev to 
  # Barcode.ID to match drug.BeatAML
  global.ids <- names(global.BeatAML)
  
  # remove X and any 0's from start of each column name and then
  # replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
  for(i in seq_len(length(global.ids))){
    global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
    
    if(substring(global.ids[i], 1, 1) == 0){
      global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
    }
    
    if(global.ids[i] %in% meta.BeatAML$SampleID.abbrev){
      global.ids[i] <- meta.BeatAML[meta.BeatAML$SampleID.abbrev == global.ids[i], ]$Barcode.ID
    }
  }
  
  # replace global.BeatAML column names 
  names(global.BeatAML) <- global.ids
  
  # transpose global.BeatAML so that first column is Barcode.ID and 
  # rest of columns are gene symbols
  global.BeatAML <- as.data.frame(t(global.BeatAML))
  
  # make first column Barcode.ID
  global.BeatAML[, "Barcode.ID"] <- rownames(global.BeatAML)
  global.BeatAML <- 
    global.BeatAML[ , c("Barcode.ID", 
                        names(global.BeatAML[ , 1:(ncol(global.BeatAML)-1)]))]
  
  ## format phospho-proteomics data frame
  # change global.BeatAML column names from SampleID.abbrev to Barcode.ID to match drug.BeatAML
  phospho.ids <- names(phospho.BeatAML)
  
  # remove X and any 0's from start of each column name and then
  # replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
  for(i in seq_len(length(phospho.ids))){
    phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
    
    if(substring(phospho.ids[i], 1, 1) == 0){
      phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
    }
    
    if(phospho.ids[i] %in% meta.BeatAML$SampleID.abbrev){
      phospho.ids[i] <- meta.BeatAML[
        meta.BeatAML$SampleID.abbrev == phospho.ids[i], ]$Barcode.ID
    }
  }
  
  # replace phospho.BeatAML column names
  names(phospho.BeatAML) <- phospho.ids
  
  # transpose phospho.BeatAML so that first column is Barcode.ID and rest of columns are gene symbols
  phospho.BeatAML <- as.data.frame(t(phospho.BeatAML))
  
  # make first column Barcode.ID
  phospho.BeatAML[, "Barcode.ID"] <- rownames(phospho.BeatAML)
  phospho.BeatAML <- phospho.BeatAML[ , c("Barcode.ID", names(phospho.BeatAML[ , 1:(ncol(phospho.BeatAML)-1)]))]
  
  return(list(meta = meta.BeatAML, drug = drug.BeatAML, gmt = gmt.drug,
              global = global.BeatAML, phospho = phospho.BeatAML))
}

load_not_norm_BeatAML_for_DMEA <- function(BeatAML.path = "BeatAML_DMEA_inputs_not_normalized") {
  message("Loading Beat AML data for DMEA")
  BeatAML_synapse_id <- list("drug_response.csv" = "syn51674470", 
                             "Ex10_metadata.txt" = "syn25807733",
                             "ptrc_ex10_crosstab_global_gene_original.txt" = "syn25714254",
                             "ptrc_ex10_crosstab_phospho_siteID_original.txt" = "syn25714936")
  
  gmt.drug <- readRDS(gzcon(url("https://raw.github.com/BelindaBGarana/panSEA/main/Examples/Inputs/gmt_BeatAML_drug_MOA.rds")))
  
  ### download files if any not already downloaded
  if (!file.exists(BeatAML.path)) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  } else if (!any(FALSE %in% lapply(names(BeatAML_synapse_id), file.exists))) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  }
  
  ### load files
  drug.BeatAML <- read.csv(file.path(BeatAML.path, names(BeatAML_synapse_id)[1]))
  meta.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[2]), 
                             sep = "\t", header = TRUE)
  global.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[3]),
                               sep = "\t", header = TRUE)
  phospho.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[4]),
                                sep = "\t", header = TRUE)
  
  ### format BeatAML data for DMEA
  sample.names <- "Barcode.ID"
  
  ## format drug sensitivity data frame
  # format drug.BeatAML wide (samples in first column, drug names for rest of columns)
  drug.BeatAML <- reshape2::dcast(drug.BeatAML, sample_id ~ inhibitor, 
                                  value.var = "auc", fill = NA)
  
  # # remove drugs without moa annotations and drug combos
  # valid.drugs <- 
  #   names(drug.BeatAML)[names(drug.BeatAML) %in% 
  #                         moa.BeatAML[!is.na(moa.BeatAML),]$Drug] # 167 drugs
  # drug.BeatAML <- drug.BeatAML[ , c("sample_id", valid.drugs)] # 167 drugs
  # moa.BeatAML <- 
  #   moa.BeatAML[moa.BeatAML$Drug %in% names(drug.BeatAML)[2:ncol(drug.BeatAML)], ]
  
  # change sample column name to match expression data
  names(drug.BeatAML)[1] <- sample.names
  
  ## format global proteomics data frame
  # change global.BeatAML column names from SampleID.abbrev to 
  # Barcode.ID to match drug.BeatAML
  global.ids <- names(global.BeatAML)
  
  # remove X and any 0's from start of each column name and then
  # replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
  for(i in seq_len(length(global.ids))){
    global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
    
    if(substring(global.ids[i], 1, 1) == 0){
      global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
    }
    
    if(global.ids[i] %in% meta.BeatAML$SampleID.abbrev){
      global.ids[i] <- meta.BeatAML[meta.BeatAML$SampleID.abbrev == global.ids[i], ]$Barcode.ID
    }
  }
  
  # replace global.BeatAML column names 
  names(global.BeatAML) <- global.ids
  
  # subtract sample medians
  sample.names <- colnames(dplyr::select_if(global.BeatAML, is.numeric))
  #global.BeatAML[,sample.names] <- log(global.BeatAML[,sample.names], 2)
  global_sample_coef <- apply(global.BeatAML[,sample.names], 2, median, na.rm = T)
  global.BeatAML[,sample.names] <- sweep(global.BeatAML[,sample.names], 2, global_sample_coef, FUN = '-')
  
  # transpose global.BeatAML so that first column is Barcode.ID and 
  # rest of columns are gene symbols
  global.BeatAML <- as.data.frame(t(global.BeatAML))
  
  # make first column Barcode.ID
  global.BeatAML[,"Barcode.ID"] <- rownames(global.BeatAML)
  global.BeatAML <- 
    global.BeatAML[ , c("Barcode.ID", 
                        names(global.BeatAML[ , 1:(ncol(global.BeatAML)-1)]))]
  
  ## format phospho-proteomics data frame
  # change global.BeatAML column names from SampleID.abbrev to Barcode.ID to match drug.BeatAML
  phospho.ids <- names(phospho.BeatAML)
  
  # remove X and any 0's from start of each column name and then
  # replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
  for(i in seq_len(length(phospho.ids))){
    phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
    
    if(substring(phospho.ids[i], 1, 1) == 0){
      phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
    }
    
    if(phospho.ids[i] %in% meta.BeatAML$SampleID.abbrev){
      phospho.ids[i] <- meta.BeatAML[
        meta.BeatAML$SampleID.abbrev == phospho.ids[i], ]$Barcode.ID
    }
  }
  
  # replace phospho.BeatAML column names
  names(phospho.BeatAML) <- phospho.ids
  
  # subtract sample medians
  sample.names <- colnames(dplyr::select_if(phospho.BeatAML, is.numeric))
  #phospho.BeatAML[,sample.names] <- log(phospho.BeatAML[,sample.names], 2)
  phospho_sample_coef <- apply(phospho.BeatAML[,sample.names], 2, median, na.rm = T)
  phospho.BeatAML[,sample.names] <- sweep(phospho.BeatAML[,sample.names], 2, phospho_sample_coef, FUN = '-')
  
  # transpose phospho.BeatAML so that first column is Barcode.ID and rest of columns are gene symbols
  phospho.BeatAML <- as.data.frame(t(phospho.BeatAML))
  
  # make first column Barcode.ID
  phospho.BeatAML[, "Barcode.ID"] <- rownames(phospho.BeatAML)
  phospho.BeatAML <- phospho.BeatAML[ , c("Barcode.ID", names(phospho.BeatAML[ , 1:(ncol(phospho.BeatAML)-1)]))]
  
  return(list(meta = meta.BeatAML, drug = drug.BeatAML, gmt = gmt.drug,
              global = global.BeatAML, phospho = phospho.BeatAML))
}


## get list of top n mountain plots (if any)
# input: EA output (e.g., mGSEA.result or mDMEA.result)
# output: list of top n mountain plots
get_top_mtn_plots <- function(base.result, n.top = 10, EA.type = "GSEA", 
                              sets = "Feature_set") {
  all.mtn <- base.result$mtn.plots
  temp.results <- base.result$result
  if (length(all.mtn) > 0) {
    if (length(all.mtn) > n.top) {
      # identify top significant enrichments
      mtn.results <- temp.results[temp.results[ , sets] %in% names(all.mtn), ]
      top.mtn.results <- mtn.results %>% slice_max(abs(NES), n = n.top)
      all.top.mtn <- all.mtn[names(all.mtn) %in% top.mtn.results[ , sets]]
    } else {
      all.top.mtn <- all.mtn
    }
    
    # create list of mtn plots for files
    mtn.file.names <- c()
    for (i in 1:length(all.top.mtn)) {
      mtn.file.names <- c(mtn.file.names, 
                          paste0(EA.type, "_mtn_plot_", names(all.top.mtn)[i], ".pdf"))
    }
    names(all.top.mtn) <- mtn.file.names
  } else {
    all.top.mtn <- list()
  }
  return(all.top.mtn)
}

# save folder and nested files/subfolders to Synapse - more efficient
save_to_synapse_v2 <- function(temp.files, resultsFolder = NULL,
                            width = 7, height = 7, dot.scale = 4) {
  for (i in names(temp.files)) {
    # if file, save appropriately
    if (grepl("[.]", i) & is.list(temp.files[[i]]) & length(temp.files[[i]]) > 0) {
      # local save
      if (endsWith(tolower(i), ".csv")) {
        write.csv(temp.files[[i]], i, row.names = FALSE)
      } else if (endsWith(tolower(i), ".pdf")) {
          if (grepl("_bar_plot", i) | grepl("_dot_plot", i)) {
            ggplot2::ggsave(i, temp.files[[i]], device = "pdf", 
                            width = dot.scale*width, height = height) 
          } else {
            ggplot2::ggsave(i, temp.files[[i]], device = "pdf", 
                            width = width, height = height)
          } 
      } else if (endsWith(tolower(i), ".svg")) {
          if (grepl("_bar_plot", i) | grepl("_dot_plot", i)) {
            ggplot2::ggsave(i, temp.files[[i]], device = "svg", 
                            width = dot.scale*width, height = height) 
          } else {
            ggplot2::ggsave(i, temp.files[[i]], device = "svg", 
                            width = width, height = height)
          } 
      } else if (endsWith(tolower(i), ".html")) {
          visNetwork::visSave(temp.files[[i]], i) 
      } else if (endsWith(tolower(i), ".bp")) {
        temp.name <- paste0(substr(i, 1, nchar(i)-3), ".pdf")
        save_base_plot(temp.files[[i]], temp.name, 
                       width = width, height = height)
        i <- temp.name
      }
      
      # save to synapse if relevant
      if (!is.null(resultsFolder)) {
        # save to synapse
        mySynFile <- synapser::File(i, parent = resultsFolder)
        synapser::synStore(mySynFile) 
      }
    } else {
     # else create subfolder
      temp.base <- getwd()
      dir.create(i)
      setwd(i)
      if (!is.null(resultsFolder)) {
        subFolder <- 
          synapser::synStore(synapser::Folder(i, parent = resultsFolder))
      } else {
        subFolder <- NULL
      }
      
      sub.base <- file.path(temp.base, i)
      save_to_synapse_v2(temp.files[[i]], subFolder, width, height, dot.scale)
      setwd(temp.base)
    }
  }
}

# save folder and nested files/subfolders to Synapse
save_to_synapse <- function(temp.files, resultsFolder = NULL,
                            width = 7, height = 7, dot.scale = 4) {
  CSV.files <- names(temp.files)[grepl(".csv", names(temp.files), 
                                       ignore.case = TRUE)]
  if (length(CSV.files) > 0) {
    # save locally
    for (j in 1:length(CSV.files)) {
      write.csv(temp.files[[CSV.files[j]]], CSV.files[j], row.names = FALSE)
    }
    
    if (!is.null(resultsFolder)) {
      # save to synapse
      CSVs <- lapply(as.list(CSV.files), synapser::File,
                     parent = resultsFolder)
      lapply(CSVs, synapser::synStore) 
    }
  }
  
  PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files), 
                                       ignore.case = TRUE)]
  if (length(PDF.files) > 0) {
    # save locally
    for (j in 1:length(PDF.files)) {
      if (is.list(temp.files[[PDF.files[j]]])) {
        if (endsWith(PDF.files[j], "_bar_plot.pdf") | 
            endsWith(PDF.files[j], "_dot_plot.pdf") |
            endsWith(PDF.files[j], "_dot_plot_withSD.pdf")) {
          ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                          device = "pdf", width = dot.scale*width, height = height) 
        } else {
          ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                          device = "pdf", width = width, height = height)
        }
      }
    }
    
    if (!is.null(resultsFolder)) {
      # save to synapse
      PDF.files <- list.files(pattern = ".*.pdf", full.names = TRUE)
      PDFs <- lapply(as.list(PDF.files), synapser::File,
                     parent = resultsFolder)
      lapply(PDFs, synapser::synStore) 
    }
  }
  
  SVG.files <- names(temp.files)[grepl(".svg", names(temp.files), 
                                       ignore.case = TRUE)]
  if (length(SVG.files) > 0) {
    # save locally
    for (j in 1:length(SVG.files)) {
      ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                      device = "svg", width = width, height = height)
    }
    
    if (!is.null(resultsFolder)) {
      # save to synapse
      SVGs <- lapply(as.list(SVG.files), synapser::File,
                     parent = resultsFolder)
      lapply(SVGs, synapser::synStore) 
    }
  }
  
  HTML.files <- names(temp.files)[grepl(".html", names(temp.files), 
                                        ignore.case = TRUE)]
  if (length(HTML.files) > 0) {
    # save locally
    for (j in 1:length(HTML.files)) {
      if (is.list(temp.files[[HTML.files[j]]])) {
        if (length(temp.files[[HTML.files[j]]]) > 0) {
          visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j]) 
        }
      } else {
        HTML.files <- HTML.files[HTML.files != HTML.files[j]]
      }
    }
    
    if (!is.null(resultsFolder)) {
      # save to synapse
      HTMLs <- lapply(HTML.files, synapser::File,
                      parent = resultsFolder)
      lapply(HTMLs, synapser::synStore) 
    }
  }
  
  base.plot.files <- names(temp.files)[grepl(".bp", names(temp.files), 
                                             ignore.case = TRUE)]
  if (length(base.plot.files) > 0) {
    # save locally
    for (j in 1:length(base.plot.files)) {
      if (is.list(temp.files[[base.plot.files[j]]])) {
        if (length(temp.files[[base.plot.files[j]]]) > 0) {
          temp.name <- paste0(substr(base.plot.files[j], 1, nchar(base.plot.files[j])-3), ".pdf")
          save_base_plot(temp.files[[base.plot.files[j]]], temp.name, 
                         width = width, height = height)
          base.plot.files[j] <- temp.name
        }
      } else {
        base.plot.files <- base.plot.files[base.plot.files != base.plot.files[j]]
      }
    }
    
    if (!is.null(resultsFolder)) {
      # save to synapse
      base.plots <- lapply(base.plot.files, synapser::File,
                           parent = resultsFolder)
      lapply(base.plots, synapser::synStore) 
    }
  }
  
  # save subfolders if relevant
  subfolders <- names(temp.files)[!grepl("[.]", names(temp.files))]
  if (length(subfolders) > 0) {
    temp.base <- getwd()
    for (m in 1:length(subfolders)) {
      # create folder for mtn plots
      dir.create(subfolders[m])
      setwd(subfolders[m])
      if (!is.null(resultsFolder)) {
        subFolder <- 
          synapser::synStore(synapser::Folder(subfolders[m],
                                              parent = resultsFolder))
      } else {
        subFolder <- NULL
      }
      
      sub.base <- file.path(temp.base, subfolders[m])
      save_to_synapse(temp.files[[subfolders[m]]], subFolder)
      setwd(temp.base)
    } 
  }
}

get_expr_for_mDMEA <- function(expression) {
  if ("adherent CCLE" %in% expression) {
    for (i in which(expression == "adherent CCLE")) {
      expression[[i]] <- get_CCLE_RNA()
    }
  } else if ("CCLE proteomics" %in% expression) {
    for (i in which(expression == "CCLE proteomics")) {
      expression[[i]] <- get_CCLE_prot()
    }
  } else if ("Beat AML proteomics" %in% expression) {
    BeatAML <- load_BeatAML_for_DMEA()
    for (i in which(expression == "Beat AML proteomics")) {
      expression[[i]] <- BeatAML$global
    }
  } else if ("Beat AML phospho" %in% expression) {
    BeatAML <- load_BeatAML_for_DMEA()
    for (i in which(expression == "Beat AML phospho")) {
      expression[[i]] <- BeatAML$phospho
    }
  }
  return(expression)
}

get_gmt_for_DMEA <- function(gmt.drug) {
  if (is.character(gmt.drug)) {
    if (gmt.drug == "PRISM") {
      message("Loading PRISM drug mechanism of action annotations")
      gmt.drug <- GSA::GSA.read.gmt(
        file = paste0(
          "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/",
          "Inputs/MOA_gmt_file_n6_no_special_chars.gmt"
        )
      )
    } else if (gmt.drug == "Beat AML") {
      BeatAML <- load_BeatAML_for_DMEA()
      gmt.drug <- BeatAML$gmt
    } else {
      stop("gmt must be either 'PRISM', 'Beat AML', or list object per documentation")
    }
  }
  return(gmt.drug)
}

get_sensitivity_for_DMEA <- function(drug.sensitivity) {
  if (is.character(drug.sensitivity)) {
    if (drug.sensitivity == "PRISM") {
      message("Loading PRISM drug sensitivity AUC scores")
      drug.sensitivity <- read.csv(file = paste0(
        "https://raw.github.com/BelindaBGarana/",
        "DMEA/shiny-app/Inputs/PRISM_drug_mean_AUC_6-23-21.csv"
      )) # 481 cell lines
      drug.sensitivity$X <- NULL
    } else if (drug.sensitivity == "Beat AML") {
      BeatAML <- load_BeatAML_for_DMEA()
      drug.sensitivity <- BeatAML$drug
    } else {
      stop(paste(
        "drug.sensitivity must be either 'PRISM', 'Beat AML', or data frame per",
        "documentation"
      ))
    }
  }
  return(drug.sensitivity)
}

prep_for_panSEA2 <- function(meta.df, omics,
                             gmt.list1 = c("KEGG" = "msigdb_Homo sapiens_C2_CP:KEGG",
                                           "Hallmark" = "msigdb_Homo sapiens_H",
                                           "Positional" = "msigdb_Homo sapiens_C1"),
                             gmt.list2 = c("ksdb_human", "sub"),
                             expr = as.list(
                               rep("CCLE proteomics", 
                                   length(omics)-ifelse(any(grepl("phospho", names(omics), 
                                                                  ignore.case=TRUE)), 1, 0))),
                             gmt.drug = "PRISM", drug.sensitivity = "PRISM", 
                             base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/") {
  
  setwd(base.path)
  types <- names(omics)
  non.phospho.types <- types[!grepl("phospho", types, ignore.case = TRUE)]
  
  # make sure samples in omics all match and have annotations in meta.df
  #meta.df <- na.omit(meta.df)
  meta.samples <- rownames(meta.df)
  sample.names <- colnames(dplyr::select_if(omics[[1]], is.numeric))
  feature.names <- colnames(omics[[1]])[!(colnames(omics[[1]]) %in% sample.names)]
  
  sample.names <- sample.names[sample.names %in% meta.samples] # only keep samples with annotations
  omics[[1]] <- omics[[1]][, c(feature.names, sample.names)]
  if (length(types) > 1) {
    for (i in 2:length(omics)) {
      temp.samples <- colnames(dplyr::select_if(omics[[i]], is.numeric))
      temp.features <- colnames(omics[[i]])[!(colnames(omics[[i]]) %in% temp.samples)]
      if (all(sample.names %in% temp.samples)) {
        omics[[i]] <- omics[[i]][,c(temp.features, sample.names)]
        feature.names <- c(feature.names, temp.features)
      } else {
        stop("All numeric column names must be in all omics data frames")
      }
    } 
  }
  
  # load set annotations
  gmt1 <- get_gmt1(gmt.list1)
  
  if (any(grepl("phospho", types, ignore.case = TRUE))) {
    phospho.index <- grep("phospho", names(omics), ignore.case = TRUE)
    gmt2 <- get_gmt2(gmt.list2, omics[[phospho.index]])
  } else {
    gmt2 <- list()
  }
  
  # load drug data for DMEA
  gmt.drug <- get_gmt_for_DMEA(gmt.drug)
  drug.sensitivity <- get_sensitivity_for_DMEA(drug.sensitivity)
  
  # load expression data for DMEA
  expr <- get_expr_for_mDMEA(expr)
  
  return(list(meta = meta.df, features = feature.names, 
              gmt1 = gmt1, gmt2 = gmt2, gmt.drug = gmt.drug, 
              drug = drug.sensitivity, expr = expr))
}

extract_files_for_save <- function(omics, deg, mDEG.results, dmea.results, 
                                   all.gsea.files, combo.gsea.files,
                                   gsea2 = NULL, cc.df, n.degs = 50, n.net = 5,
                                   show_colnames = FALSE, fontsize = 10, 
                                   KSEA = FALSE, SSEA = FALSE, gmt2 = NULL, scale = TRUE, cluster = TRUE,
                                   ties = FALSE) {
  # create mountain plots and network graphs for KSEA, SSEA
  if (!is.null(gsea2)) {
    n.phospho <- grep("phospho", names(deg), ignore.case = TRUE)
    gsea2.inputs <- list(deg[[n.phospho]], deg[[n.phospho]])
    kin.mtn2 <- get_top_mtn_plots(gsea2$all.results[[1]],
                                  EA.type = "KSEA")
    if (ties) {
      if (KSEA) {
        kin.net2 <- panSEA::netSEA(list(gsea2.inputs[[1]]),
                                   list(gsea2$all.results[[1]]$result.w.ties),
                                   "SUB_SITE",
                                   n.network.sets = n.net) 
      } else {
        kin.net2 <- list()
      }
      
      if (SSEA) {
        sub.net2 <- panSEA::netSEA(list(gsea2.inputs[[2]]),
                                   list(gsea2$all.results[[2]]$result.w.ties),
                                   "SUB_SITE",
                                   n.network.sets = n.net) 
      } else {
        sub.net2 <- list()
      }
    } else {
      if (KSEA) {
        kin.net2 <- panSEA::netSEA(list(gsea2.inputs[[1]]),
                                   list(gsea2$all.results[[1]]$result),
                                   "SUB_SITE",
                                   n.network.sets = n.net) 
      } else {
        kin.net2 <- list()
      }
      
      if (SSEA) {
        sub.net2 <- panSEA::netSEA(list(gsea2.inputs[[2]]),
                                   list(gsea2$all.results[[2]]$result),
                                   "SUB_SITE",
                                   n.network.sets = n.net) 
      } else {
        sub.net2 <- list()
      }
    }
    sub.mtn2 <- get_top_mtn_plots(gsea2$all.results[[2]],
                                  EA.type = "Substrate_enrichment")
    if (length(kin.mtn2) > 1) {
      if (ties) {
        kin.poi <- get_pathways_of_interest(omics[[n.phospho]], 
                                            gsea2$all.results[[1]]$result.w.ties, 
                                            gmt2[[1]], cc.df, n = n.net, 
                                            show_colnames = show_colnames, 
                                            fontsize = fontsize, scale = scale, cluster = cluster) 
      } else {
        kin.poi <- get_pathways_of_interest(omics[[n.phospho]], 
                                            gsea2$all.results[[1]]$result, 
                                            gmt2[[1]], cc.df, n = n.net, 
                                            show_colnames = show_colnames, 
                                            fontsize = fontsize, scale = scale, cluster = cluster)
      }
    } else {kin.poi <- list()}
    
    if (length(sub.mtn2) > 1) {
      sub.poi <- get_pathways_of_interest(omics[[n.phospho]], 
                                          gsea2$all.results[[2]]$result, 
                                          gmt2[[2]], cc.df, n = n.net,
                                          show_colnames = show_colnames, 
                                          fontsize = fontsize, scale = scale, cluster = cluster) 
    } else {sub.poi <- list()} 
  }
  
  all.files <- list()
  types <- names(omics)
  for (i in 1:length(types)) {
    sig.degs <- na.omit(deg[[i]][deg[[i]]$adj.P.Val <= 0.05, ])
    
    # select numeric data
    just.expr <- dplyr::select_if(omics[[i]], is.numeric)
    
    # identify feature name
    feature.name <- colnames(omics[[i]])[!(colnames(omics[[i]]) %in% colnames(just.expr))]

    if (nrow(sig.degs) > 0) {
      top.top.sig.degs <- sig.degs %>% slice_max(Log2FC, n = n.degs/2)
      top.bot.sig.degs <- sig.degs %>% slice_min(Log2FC, n = n.degs/2)
      top.sig.degs <- rbind(top.top.sig.degs, top.bot.sig.degs)
      heatmap.df <- omics[[i]][omics[[i]][,feature.name] %in% top.sig.degs[,feature.name],
                               c(feature.name, colnames(just.expr))]
      rownames(heatmap.df) <- heatmap.df[,feature.name]
      feature.order <- top.sig.degs[order(top.sig.degs$Log2FC, decreasing = TRUE), feature.name]
      heatmap.df <- heatmap.df[feature.order,]
      heatmap.mat <- heatmap.df[,2:ncol(heatmap.df)]
      heatmap.mat <- filter_for_hclust(heatmap.mat)
      
      # create heatmaps
      if (nrow(heatmap.mat) > 1) {
        heatmap.mat <- as.matrix(heatmap.mat)
        deg.heatmap.clust <- pheatmap::pheatmap(heatmap.mat, color = 
                                                  colorRampPalette(
                                                    c("navy", "white", "firebrick3"))(50),
                                                scale = "row", annotation_col = cc.df, 
                                                angle_col = "45", 
                                                show_colnames = show_colnames, 
                                                fontsize = fontsize)
        deg.heatmap <- pheatmap::pheatmap(heatmap.mat, 
                                          color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                          cluster_row = FALSE, 
                                          scale = "row", annotation_col = cc.df, 
                                          angle_col = "45", 
                                          show_colnames = show_colnames, 
                                          fontsize = fontsize)
        deg.heatmap.abs <- pheatmap::pheatmap(heatmap.mat, color = 
                                                colorRampPalette(
                                                  c("navy", "white", "firebrick3"))(50), 
                                              annotation_col = cc.df, 
                                              angle_col = "45", 
                                              show_colnames = show_colnames, 
                                              fontsize = fontsize)
      } else {
        deg.heatmap.clust <- list()
        deg.heatmap <- list()
        deg.heatmap.abs <- list()
      }
      
      temp.DEG.files <- list("Differential_expression_results.csv" = 
                               deg[[i]],
                             "Differential_expression_results_max_5_percent_FDR.csv" = 
                               sig.degs,
                             "Differential_expression_for_heatmap.csv" =
                               heatmap.df,
                             "Differential_expression_heatmap_not_scaled.bp" =
                               deg.heatmap.abs,
                             "Differential_expression_heatmap_scaled_ordered_by_Log2FC.bp" =
                               deg.heatmap,
                             "Differential_expression_heatmap_scaled.bp" =
                               deg.heatmap.clust) 
    } else {
      temp.DEG.files <- list("Differential_expression_results.csv" = 
                               deg[[i]])
    }
    
    if (grepl("phospho", types[i], ignore.case = TRUE)) {
      if (KSEA) {
        kin.DMEA.mtn <- get_top_mtn_plots(dmea.results$all.results[["phospho_ksdb"]],
                                          sets = "Drug_set",
                                          EA.type = "DMEA")
        DMEA.kin.net <- panSEA::netSEA(list(dmea.results$all.results[["phospho_ksdb"]]$corr.result),
                                       list(dmea.results$all.results[["phospho_ksdb"]]$result),
                                       "Drug", "Pearson.est",
                                       n.network.sets = n.net)
        kin.DMEA.files <- list("DMEA_results.csv" =
                                 dmea.results$all.results[["phospho_ksdb"]]$result,
                               "DMEA_correlation_results.csv" = 
                                 dmea.results$all.results[["phospho_ksdb"]]$corr.result,
                               "DMEA_correlation_scatter_plots.pdf" = 
                                 dmea.results$all.results[["phospho_ksdb"]]$corr.scatter.plots,
                               "DMEA_volcano_plot.pdf" =
                                 dmea.results$all.results[["phospho_ksdb"]]$volcano.plot,
                               "DMEA_network_graph.html" = 
                                 DMEA.kin.net$interactive,
                               "mtn_plots" = kin.DMEA.mtn)
        kin.gsea.files <- all.gsea.files[["phospho_ksdb"]]
      } else {
        kin.DMEA.files <- list()
        kin.gsea.files <- list()
      }
      
      if (SSEA) {
        sub.DMEA.mtn <- get_top_mtn_plots(dmea.results$all.results[["phospho_sub"]],
                                          sets = "Drug_set",
                                          EA.type = "DMEA")
        DMEA.sub.net <- panSEA::netSEA(list(dmea.results$all.results[["phospho_sub"]]$corr.result),
                                       list(dmea.results$all.results[["phospho_sub"]]$result),
                                       "Drug", "Pearson.est",
                                       n.network.sets = n.net)
        sub.DMEA.files <- list("DMEA_results.csv" =
                                 dmea.results$all.results[["phospho_sub"]]$result,
                               "DMEA_correlation_results.csv" = 
                                 dmea.results$all.results[["phospho_sub"]]$corr.result,
                               "DMEA_correlation_scatter_plots.pdf" = 
                                 dmea.results$all.results[["phospho_sub"]]$corr.scatter.plots,
                               "DMEA_volcano_plot.pdf" =
                                 dmea.results$all.results[["phospho_sub"]]$volcano.plot,
                               "DMEA_network_graph.html" = 
                                 DMEA.sub.net$interactive,
                               "mtn_plots" = sub.DMEA.mtn)
        sub.gsea.files <- all.gsea.files[["phospho_sub"]]
      } else {
        sub.DMEA.files <- list()
        sub.gsea.files <- list()
      }
      
      phospho.kin.files <- list("KSEA_results.csv" =
                                  gsea2$all.results[["phospho_ksdb"]]$result,
                                "KSEA_volcano_plot.pdf" =
                                  gsea2$all.results[["phospho_ksdb"]]$volcano.plot,
                                "KSEA_network_graph.html" = 
                                  kin.net2$interactive,
                                "mtn_plots" = kin.mtn2,
                                "Pathways_of_interest" = kin.poi,
                                "GSEA" = kin.gsea.files,
                                "DMEA" = kin.DMEA.files) 
      phospho.sub.files <- list("Substrate_enrichment_results.csv" =
                                  gsea2$all.results[["phospho_sub"]]$result,
                                "Substrate_enrichment_volcano_plot.pdf" =
                                  gsea2$all.results[["phospho_sub"]]$volcano.plot,
                                "Substrate_enrichment_network_graph.html" = 
                                  sub.net2$interactive,
                                "mtn_plots" = sub.mtn2,
                                "Pathways_of_interest" = sub.poi,
                                "GSEA" = sub.gsea.files,
                                "DMEA" = sub.DMEA.files)
      phospho.files <- list('Differential_expression' = temp.DEG.files, 
                            'KSEA' = phospho.kin.files,
                            'Substrate_enrichment' = phospho.sub.files)
      all.files[[types[i]]] <- phospho.files
    } else {
      DMEA.global.mtn <- get_top_mtn_plots(dmea.results$all.results[[i]],
                                           sets = "Drug_set",
                                           EA.type = "DMEA")
      DMEA.global.net <- panSEA::netSEA(list(dmea.results$all.results[[i]]$corr.result),
                                        list(dmea.results$all.results[[i]]$result),
                                        "Drug", "Pearson.est",
                                        n.network.sets = n.net)
      global.DMEA.files <- list("DMEA_results.csv" =
                                  dmea.results$all.results[[i]]$result,
                                "DMEA_correlation_results.csv" = 
                                  dmea.results$all.results[[i]]$corr.result,
                                "DMEA_correlation_scatter_plots.pdf" = 
                                  dmea.results$all.results[[i]]$corr.scatter.plots,
                                "DMEA_volcano_plot.pdf" =
                                  dmea.results$all.results[[i]]$volcano.plot,
                                "DMEA_network_graph.html" = 
                                  DMEA.global.net$interactive,
                                "mtn_plots" = DMEA.global.mtn)
      global.files <- list('Differential_expression' = temp.DEG.files, 
                           'GSEA' = all.gsea.files[[types[i]]],
                           'DMEA' = global.DMEA.files)
      all.files[[types[i]]] <- global.files
    }
  }
  
  ## combo
  if (length(types) > 1) {
    combo.files <- list()
    if (!grepl("phospho", types[i], ignore.case = TRUE)) {
    combo.files[["Differential_expression"]] <- list("Differential_expression_results.csv" =
                                                       mDEG.results$compiled.results$results,
                                                     "Differential_expression_mean_results.csv" =
                                                       mDEG.results$compiled.results$mean.results,
                                                     "Differential_expression_correlation_matrix.pdf" =
                                                       mDEG.results$compiled.results$corr.matrix,
                                                     "Differential_expression_dot_plot.pdf" =
                                                       mDEG.results$compiled.results$dot.plot)
    }
    combo.files[["GSEA"]] <- combo.gsea.files
    if (length(dmea.results$compiled.results) > 1) {
      combo.files[["DMEA"]] <- list("DMEA_results.csv" =
                                      dmea.results$compiled.results$results,
                                    "DMEA_mean_results.csv" =
                                      dmea.results$compiled.results$mean.results,
                                    "DMEA_correlation_matrix.pdf" =
                                      dmea.results$compiled.results$corr.matrix,
                                    "DMEA_dot_plot.pdf" =
                                      dmea.results$compiled.results$dot.plot) 
    }
    combo.name <- paste0(types, collapse = "_and_")
    all.files[[combo.name]] <- combo.files
  }
  return(all.files)
}

extract_DMEA_files <- function(dmea.results, index, n.net=5) {
  kin.DMEA.mtn <- get_top_mtn_plots(dmea.results$all.results[[index]],
                                    sets = "Drug_set",
                                    EA.type = "DMEA")
  DMEA.kin.net <- panSEA::netSEA(list(dmea.results$all.results[[index]]$corr.result),
                                 list(dmea.results$all.results[[index]]$result),
                                 "Drug", "Pearson.est",
                                 n.network.sets = n.net)
  kin.DMEA.files <- list("DMEA_results.csv" =
                           dmea.results$all.results[[index]]$result,
                         "DMEA_correlation_results.csv" = 
                           dmea.results$all.results[[index]]$corr.result,
                         "DMEA_correlation_scatter_plots.pdf" = 
                           dmea.results$all.results[[index]]$corr.scatter.plots,
                         "DMEA_volcano_plot.pdf" =
                           dmea.results$all.results[[index]]$volcano.plot,
                         "DMEA_network_graph.html" = 
                           DMEA.kin.net$interactive,
                         "mtn_plots" = kin.DMEA.mtn)
  return(kin.DMEA.files)
}

## run panSEA across global & phospho types for multiple contrasts
# input: vector of contrasts, contrast.type column for meta.df extraction, 
#   omics data.list, beatAML expression list, gmt set information, 
#   BeatAML drug AUC, base.path, Synapse ID info
# output: panSEA result files are saved locally & uploaded to Synapse

# inputs:
# meta.df: data frame where sample IDs are rownames, contrasts and contrast2 are columns,
# other columns can be used for filtering
# meta.df: data frame where sample IDs are rownames, columns are sample annotations for heatmaps
# gmt.list1: all gene sets to evaluate 
panSEA2 <- function(contrasts, contrast2 = NULL, meta.df, omics, 
                    annotations = c("Sort Type", "Sample Type", "Patient"),
                    gmt.list1 = c("KEGG" = "msigdb_Homo sapiens_C2_CP:KEGG",
                                  "Hallmark" = "msigdb_Homo sapiens_H"),
                    gmt.list2 = c("ksdb_human", "sub"),
                    expr = as.list(
                      rep("CCLE proteomics", 
                          length(omics)-ifelse(any(grepl("phospho", names(omics), 
                                                         ignore.case=TRUE)), 1, 0))),
                    gmt.drug = "PRISM", drug.sens = "PRISM", 
                    base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                    temp.path = base.path, subfolder = TRUE, synapse_id = NULL, 
                    filter = NA, filterID = NULL, n.degs = 50, n.net = 5,
                    width = 7, height = 7, show_colnames = FALSE, fontsize = 10,
                    scale = TRUE, cluster = TRUE) {
  # prep to run contrasts
  types <- names(omics)
  EA.types <- names(gmt.list1)
  
  prep <- prep_for_panSEA2(meta.df, omics,
                           gmt.list1,
                           gmt.list2,
                           expr,
                           gmt.drug, drug.sens, 
                           base.path)
  feature.names <- prep$features
  gmt1 <- prep$gmt1
  names(gmt1) <- EA.types
  gmt2 <- prep$gmt2
  gmt.drug <- prep$gmt.drug
  drug.sensitivity <- prep$drug
  expr <- prep$expr
  
  # filter meta.df
  if (!is.null(filterID) & !is.na(filter)) {
    meta.df <- meta.df[meta.df[,filterID] == filter,] 
  }
  
  all.degs <- data.frame()
  
  ## for each contrast: 
  for (k in 1:length(contrasts)) {
    contrast.type <- contrasts[k]
    
    # require 2 levels for contrast
    if (length(unique(meta.df[,contrast.type])) != 2) {
      next
    }
    
    # prepare annotations
    if (!is.null(contrast2) & ncol(meta.df) > 1) {
      factor.info <- meta.df[,c(contrast.type, contrast2)]
      colnames(factor.info) <- c(contrast.type, contrast2)
      levels(factor.info[,2]) <- unique(factor.info[,2])
    } else {
      factor.info <- as.data.frame(meta.df[,contrast.type])
      rownames(factor.info) <- rownames(meta.df)
      colnames(factor.info) <- contrast.type
    }
    if ("TRUE" %in% factor.info[,1]) {
      factor.info[factor.info[,1],1] <- "True"
      factor.info[factor.info[,1]=="FALSE",1] <- "False"
      levels(factor.info[,1]) <- c("True", "False")
    } else if ("Sensitive" %in% factor.info[,1]) {
      levels(factor.info[,1]) <- c("Sensitive", "Resistant")
    } else if ("Pos" %in% factor.info[,1]) {
      levels(factor.info[,1]) <- c("Pos", "Neg")
    } else if ("MSC" %in% factor.info[,1]) {
      levels(factor.info[,1]) <- c("MSC", "Non_MSC")
    } else if ("Strongly_amplified" %in% factor.info[,1]) {
      levels(factor.info[,1]) <- c("Strongly_amplified", "Not_amplified")
    } else if ("Amplified" %in% factor.info[,1]) {
      levels(factor.info[,1]) <- c("Amplified", "Not_amplified")
    } else {
      levels(factor.info[,1]) <- unique(factor.info[,1])
    }
    contrast.name <- paste0(contrast.type, "_", levels(factor.info[,1])[1], 
                            "_vs_", levels(factor.info[,1])[2])
    
    setwd(temp.path)
    if (subfolder) {
      if (!dir.exists(file.path(contrast.name))) {
        dir.create(file.path(contrast.name))
        setwd(file.path(contrast.name))
        contrastFolder <- 
          synapser::synStore(synapser::Folder(file.path(contrast.name),
                                              parent = synapse_id)) 
      } else {
        contrastFolder <- 
          synapser::synStore(synapser::Folder(file.path(contrast.name),
                                              parent = synapse_id)) 
        contrastDEGs <- as.list(synapser::synGetChildren(contrastFolder, list("file"), sortBy = 'NAME'))
        if (length(contrastDEGs) > 0) {
          if (contrastDEGs[[1]]$name == "Differential_expression_results.csv") {
            contrastFile <- synapser::synGet(contrastDEGs[[1]]$id)
            contrastDEG <- read.csv(contrastFile$path)
            
            if (is.na(filter) & is.null(filterID)) {contrastDEG$Filter <- NA} else {
              contrastDEG$Filter <- paste0(filterID, "_", filter) 
            }
            contrastDEG$Contrast <- contrast.name
            
            all.degs <- rbind(all.degs, contrastDEG)
          }
        }
        next
      }
    } else {
      contrastFolder <- synapse_id
    }
    
    # prep annotations for heatmaps
    annotations.df <- as.data.frame(meta.df[,annotations])
    colnames(annotations.df) <- annotations
    cc.df <- cbind(annotations.df, factor.info)
    
    # run panSEA across omics types
    # CAUTION: this only works because the # of samples for each treatment type 
    # is equal; otherwise would have to run panSEA for each contrast separately 
    # and perhaps set the group.samples input parameter for panSEA
    
    # run global against KEGG sets, phospho against kinase sets from ksdb
    # run differential expression analysis
    if (is.null(filterID) & is.na(filter)) {
      print(paste("Running", contrast.name, "with no filter"))
    } else{
      print(paste("Running", contrast.name, "with", filterID, "==", filter))
    }
    
    mDEG.results <- panSEA::mDEG(omics, factor.info, feature.names)
    deg <- mDEG.results$all.results
    
    ## for phospho data:
    # run GSEA for each gmt in gmt2 and then check coverage of 2+ sets in gmt1
    gsea1.inputs <- deg
    features1 <- feature.names
    rank.var <- rep("Log2FC", length(feature.names))
    KSEA <- FALSE
    SSEA <- FALSE
    if (any(grepl("phospho", names(deg), ignore.case = TRUE))) {
      n.phospho <- grep("phospho", names(deg), ignore.case = TRUE)
      gsea1.inputs <- gsea1.inputs[[-n.phospho]]
      if (is.data.frame(gsea1.inputs)) {
        gsea1.inputs <- list(gsea1.inputs)
        names(gsea1.inputs) <- types[1]
      }
      features1 <- features1[-n.phospho]
      rank.var <- rank.var[-n.phospho]
      if (length(n.phospho) > 1) {
        stop("Currently not accepting more than 1 phospho input")
      }
      names(gmt2) <- c("phospho_ksdb", "phospho_sub")
      gsea2.inputs <- list()
      for (n in 1:length(gmt.list2)) {
        gsea2.inputs[[n]] <- deg[[n.phospho]]
      }
      names(gsea2.inputs) <- names(gmt2)
      gsea2 <- panSEA::mGSEA(gsea2.inputs, gmt2, types = names(gmt2),
                             feature.names = rep(feature.names[n.phospho], length(gmt2))) 
      if (grepl("ksdb", names(gmt2)[1], ignore.case = TRUE)){
        if (nrow(gsea2$all.results[[1]]$result) > 100) {
          KSEA <- TRUE
          gsea1.inputs[["phospho_ksdb"]] <- gsea2$all.results[[1]]$result
          expr[[length(expr)+1]] <- expr[[1]]
          features1 <- c(features1, "Feature_set")
          rank.var <- c(rank.var, "NES")
        } else {
          KSEA <- FALSE
        }
        
        if (grepl("sub", names(gmt2)[2], ignore.case = TRUE)){
          if (nrow(gsea2$all.results[[2]]$result) > 100) {
            SSEA <- TRUE
            gsea1.inputs[["phospho_sub"]] <- gsea2$all.results[[2]]$result
            expr[[length(expr)+1]] <- expr[[1]]
            features1 <- c(features1, "Feature_set")
            rank.var <- c(rank.var, "NES")
          } else {
            SSEA <- FALSE
          }
        }
      } else {
        if (grepl("sub", names(gmt2)[1], ignore.case = TRUE)){
          if (nrow(gsea2$all.results[[1]]$result) > 100) {
            SSEA <- TRUE
            gsea1.inputs[["phospho_sub"]] <- gsea2$all.results[[1]]$result
            expr[[length(expr)+1]] <- expr[[1]]
            features1 <- c(features1, "Feature_set")
            rank.var <- c(rank.var, "NES")
          } else {
            SSEA <- FALSE
          }
        }
      }
    } else {
      gsea2 <- NULL
      n.phospho <- NULL
    }
    
    # run GSEA for each gmt in gmt1
    gsea1 <- list()
    combo.gsea.files <- list()
    for (i in 1:length(gmt1)) {
      gsea.name <- paste0("GSEA_", EA.types[i])
      
      # prep set annotations
      temp.gmt1 <- list()
      for (j in 1:length(gsea1.inputs)) {
        temp.gmt1[[j]] <- gmt1[[i]]
      }
      
      # run GSEA
      gsea1[[gsea.name]] <- panSEA::mGSEA(gsea1.inputs, temp.gmt1, 
                                          types = names(gsea1.inputs), 
                                          feature.names = features1,
                                          rank.var = rank.var)
      if (length(gsea1.inputs) > 1) {
        combo.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                gsea1[[gsea.name]]$compiled.results$results,
                                              "GSEA_mean_results.csv" =
                                                gsea1[[gsea.name]]$compiled.results$mean.results,
                                              "GSEA_correlation_matrix.pdf" =
                                                gsea1[[gsea.name]]$compiled.results$corr.matrix,
                                              "GSEA_dot_plot.pdf" =
                                                gsea1[[gsea.name]]$compiled.results$dot.plot)
      }
    }
    
    # store results for each input type
    all.gsea.files <- list()
    for (j in 1:length(gsea1.inputs)) {
      global.gsea.files <- list()
      for (i in 1:length(EA.types)) {
        gsea.name <- paste0("GSEA_", EA.types[i])
        global.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[[j]], 
          EA.type = EA.types[i])
        if (length(global.mtn) > 1) {
          global.net <- panSEA::netSEA(list(gsea1.inputs[[j]]),
                                       list(gsea1[[gsea.name]]$all.results[[j]]$result),
                                       element.names = features1[j],
                                       rank.var = rank.var[j],
                                       n.network.sets = n.net)
          global.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                   gsea1[[gsea.name]]$all.results[[j]]$result,
                                                 "GSEA_volcano_plot.pdf" =
                                                   gsea1[[gsea.name]]$all.results[[j]]$volcano.plot,
                                                 "GSEA_network_graph.html" = 
                                                   global.net$interactive,
                                                 "mtn_plots" = global.mtn)
          
          if (!grepl("phospho", names(gsea1.inputs)[j], ignore.case = TRUE)) {
            global.gsea.files[[gsea.name]][["Pathways_of_interest"]] <- 
              get_pathways_of_interest(omics[[names(gsea1.inputs)[j]]], 
                                       gsea1[[gsea.name]]$all.results[[j]]$result, 
                                       gmt1[[i]], cc.df, n = n.net,
                                       show_colnames = show_colnames, 
                                       fontsize = fontsize, scale = scale, cluster = cluster) 
          }
        } else {
          global.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                   gsea1[[gsea.name]]$all.results[[j]]$result,
                                                 "GSEA_volcano_plot.pdf" =
                                                   gsea1[[gsea.name]]$all.results[[j]]$volcano.plot,
                                                 "mtn_plots" = global.mtn)
        }
      }
      all.gsea.files[[j]] <- global.gsea.files
    }
    names(all.gsea.files) <- names(gsea1.inputs)
    
    # run DMEA
    dmea.results <- panSEA::mDMEA(drug.sens, gmt.drug, expr, gsea1.inputs, 
                                  names(gsea1.inputs), 
                                  feature.names = features1,
                                  weight.values = rank.var)
    
    #### save results & upload to Synapse
    ### set file names
    ## global
    all.files <- extract_files_for_save(omics, deg, mDEG.results, dmea.results, 
                                        all.gsea.files, combo.gsea.files,
                                        gsea2, cc.df, 
                                        n.degs, n.net, show_colnames, 
                                        fontsize, KSEA, SSEA, gmt2)
    
    save_to_synapse(all.files, contrastFolder)
    
    # compile DEGs
    # add feature type for each omics
    for (i in 1:length(deg)) {
      deg[[i]]$Feature_type <- colnames(deg[[i]])[1]
      colnames(deg[[i]])[1] <- "Feature"
    }
    
    # rbind and add contrast information
    temp.degs <- na.omit(data.table::rbindlist(deg, use.names = TRUE, fill = TRUE))
    if (is.na(filter) & is.null(filterID)) {temp.degs$Filter <- NA} else {
      temp.degs$Filter <- paste0(filterID, "_", filter) 
    }
    temp.degs$Contrast <- contrast.name
    all.degs <- rbind(all.degs, temp.degs)
    
    # make space to process next contrast
    gsea1 <- NULL
    gsea2 <- NULL
    dmea.results <- NULL
    deg <- NULL
    gsea1.inputs <- NULL
  }
  
  if (nrow(all.degs) > 0) {
    all.DEG.files <- list("Differential_expression_results.csv" = 
                            all.degs,
                          "Differential_expression_results_max_5_percent_FDR.csv" = 
                            all.degs[all.degs$adj.P.Val <= 0.05, ])
    setwd(temp.path)
    save_to_synapse(all.DEG.files, synapse_id)
  }
}



panSEA2_combos <- function(contrasts, contrast2 = NULL, meta.df, omics, 
                           annotations = c("Sort Type", "Sample Type", "Patient"),
                           gmt.list1 = c("KEGG" = "msigdb_Homo sapiens_C2_CP:KEGG",
                                         "Hallmark" = "msigdb_Homo sapiens_H",
                                         "Positional" = "msigdb_Homo sapiens_C1"),
                           gmt.list2 = c("ksdb_human", "sub"),
                           expr = as.list(
                             rep("CCLE proteomics", 
                                 length(omics)-ifelse(any(grepl("phospho", names(omics), 
                                                                ignore.case=TRUE)), 
                                                      length(gmt.list2)-1, 0))),
                           gmt.drug = "PRISM", drug.sens = "PRISM", 
                           base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                           temp.path = base.path, subfolder = TRUE, synapse_id = NULL, 
                           filters = contrasts, n.degs = 50, n.net = 5,
                           width = 7, height = 7, show_colnames = FALSE, fontsize = 10) {
  # prep to run contrasts
  types <- names(omics)
  prep <- prep_for_panSEA2(meta.df, omics,
                           gmt.list1,
                           gmt.list2,
                           expr,
                           gmt.drug, drug.sens, 
                           base.path)
  meta.df <- prep$meta
  feature.names <- prep$features
  gmt1 <- prep$gmt1
  gmt2 <- prep$gmt2
  gmt.drug <- prep$gmt.drug
  drug.sensitivity <- prep$drug
  expr <- prep$expr
  
  all.degs <- data.frame()
  # run contrasts with no filters
  setwd(temp.path)
  dir.create("no_filter")
  setwd("no_filter")
  nullPath <- file.path(temp.path, "no_filter")
  nullFolder <- 
    synapser::synStore(synapser::Folder("no_filter",
                                        parent = synapse_id))
  panSEA2(contrasts, contrast2, meta.df, omics, annotations, gmt.list1, gmt.list2,
          expr, gmt.drug, drug.sens, base.path, nullPath, subfolder, 
          synapse_id = nullFolder, filter = NA, filterID = NULL, 
          n.degs = n.degs, n.net = n.net, width = width, height = height, 
          fontsize = fontsize)
  
  nullFiles <- as.list(synapser::synGetChildren(nullFolder, list("file"), sortBy = 'NAME'))
  if (length(nullFiles) > 0) {
    if (nullFiles[[1]]$name == "Differential_expression_results.csv") {
      nullFile <- synapser::synGet(nullFiles[[1]]$id)
      nullDEGs <- read.csv(nullFile$path)
      all.degs <- rbind(all.degs, nullDEGs)
    }
  }
  
  for (m in 1:length(filters)) {
    filter.types <- unique(meta.df[,filters[m]])
    
    # run contrasts with each filter
    for (n in 1:length(filter.types)) {
      setwd(temp.path)
      dir.create(file.path(paste0(filters[m], "_", filter.types[n])))
      setwd(file.path(paste0(filters[m], "_", filter.types[n])))
      truePath <- file.path(temp.path, paste0(filters[m], "_", filter.types[n]))
      trueFolder <- 
        synapser::synStore(synapser::Folder(file.path(paste0(filters[m], "_", 
                                                             filter.types[n])),
                                            parent = synapse_id))
      panSEA2(contrasts, contrast2, meta.df, omics, annotations, gmt.list1, gmt.list2,
              expr, gmt.drug, drug.sens, base.path, truePath, subfolder, 
              synapse_id = trueFolder, filter = filter.types[n], filterID = filters[m], 
              n.degs = n.degs, n.net = n.net, width = width, height = height, 
              fontsize = fontsize)
      trueDEGfiles <- as.list(synapser::synGetChildren(trueFolder, list("file"), sortBy = 'NAME'))
      if (length(trueDEGfiles) > 0) {
        if (trueDEGfiles[[1]]$name == "Differential_expression_results.csv") {
          trueDEGfile <- synapser::synGet(trueDEGfiles[[1]]$id)
          trueDEGs <- read.csv(trueDEGfile$path)
          all.degs <- rbind(all.degs, trueDEGs)
        }
      } 
    }
  }
  
  if (nrow(all.degs) > 0) {
    setwd(temp.path)
    all.DEG.files <- list("Differential_expression_results.csv" = 
                            all.degs,
                          "Differential_expression_results_max_5_percent_FDR.csv" = 
                            all.degs[all.degs$adj.P.Val <= 0.05, ])
    save_to_synapse(all.DEG.files, synapse_id)
  }
}

# inputs:
# omics: list of dataframes where first column is feature names, rest of columns are samples
# meta.list: list of dataframes with metadata matching omics dataframes where sample IDs are in first column
# feature.list: list of column names containing features in corresponding omics dataframes
# temp.path: where the files should be stored locally
# syn.id: synapse ID where files should be stored
panSEA_corr <- function(omics, meta.list, feature.list, rank.col = "Gain of C8",
                        other.annotations = c("Ave CEP8", "Ave CMYC"),
                        expr.list, temp.path, syn.id = NULL, 
                        n.top = 50, n.top.sets = 10) {
  for (i in 1:length(omics)) {
    setwd(file.path(temp.path))
    dir.create(names(omics)[i])
    setwd(names(omics)[i])
    temp.features <- feature.list[[i]]
    temp.omics <- omics[[i]]
    meta.df <- as.data.frame(meta.list[[i]])
    rownames(meta.df) <- meta.df[,1]
    meta.df$Sample <- rownames(meta.df)
    red.meta <- meta.df[,c("Sample", rank.col)]
    cc.df.global <- meta.df[,c("Sample", other.annotations, rank.col)]
    
    omics.files <- list()
    for (j in 1:length(temp.omics)) {
      temp.omics.files <- list()
      omics.t <- as.data.frame(t(temp.omics[[j]]))
      colnames(omics.t) <- omics.t[1,]
      omics.t <- omics.t[2:nrow(omics.t),]
      omics.t$Sample <- rownames(omics.t)
      # need to make sample names first col, gene symbols colnames
      omics.t <- merge(red.meta, omics.t)
      
      deg <- DMEA::rank_corr(omics.t, variable=temp.features[j], 
                             value = names(temp.omics)[j], plots = FALSE)
      sig.degs <- na.omit(deg$result[deg$result$Pearson.q <= 0.05, ])
      
      # select numeric data
      just.expr <- dplyr::select_if(temp.omics[[j]], is.numeric)
      
      # identify feature name
      feature.name <- temp.features[j]
      
      if (nrow(sig.degs) > 0) {
        top.top.sig.degs <- sig.degs %>% slice_max(Pearson.est, n = n.top/2)
        top.bot.sig.degs <- sig.degs %>% slice_min(Pearson.est, n = n.top/2)
        top.sig.degs <- rbind(top.top.sig.degs, top.bot.sig.degs)
        heatmap.df <- temp.omics[[j]][temp.omics[[j]][,feature.name] %in% top.sig.degs[,feature.name],
                                      c(feature.name, colnames(just.expr))]
        rownames(heatmap.df) <- heatmap.df[,feature.name]
        feature.order <- top.sig.degs[order(top.sig.degs$Pearson.est, decreasing = TRUE), feature.name]
        heatmap.df <- heatmap.df[feature.order,]
        heatmap.mat <- heatmap.df[,2:ncol(heatmap.df)]
        heatmap.mat <- filter_for_hclust(heatmap.mat)
        
        # create heatmaps
        if (nrow(heatmap.mat) > 1) {
          heatmap.mat <- as.matrix(heatmap.mat)
          deg.heatmap.clust <- pheatmap::pheatmap(heatmap.mat, color = 
                                                    colorRampPalette(
                                                      c("navy", "white", "firebrick3"))(50),
                                                  scale = "row", annotation_col = cc.df.global, 
                                                  angle_col = "45", 
                                                  show_colnames = FALSE, 
                                                  fontsize = 10)
          deg.heatmap <- pheatmap::pheatmap(heatmap.mat, 
                                            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                            cluster_row = FALSE, 
                                            scale = "row", annotation_col = cc.df.global, 
                                            angle_col = "45", 
                                            show_colnames = FALSE, 
                                            fontsize = 10)
          deg.heatmap.abs <- pheatmap::pheatmap(heatmap.mat, color = 
                                                  colorRampPalette(
                                                    c("navy", "white", "firebrick3"))(50), 
                                                annotation_col = cc.df.global, 
                                                angle_col = "45", 
                                                show_colnames = FALSE, 
                                                fontsize = 10)
        } else {
          deg.heatmap.clust <- list()
          deg.heatmap <- list()
          deg.heatmap.abs <- list()
        }
        
        temp.DEG.files <- list("Differential_expression_results.csv" = 
                                 deg$result,
                               "Differential_expression_results_max_5_percent_FDR.csv" = 
                                 sig.degs,
                               "Differential_expression_for_heatmap.csv" =
                                 heatmap.df,
                               "Differential_expression_heatmap_not_scaled.bp" =
                                 deg.heatmap.abs,
                               "Differential_expression_heatmap_scaled_ordered_by_Log2FC.bp" =
                                 deg.heatmap,
                               "Differential_expression_heatmap_scaled.bp" =
                                 deg.heatmap.clust) 
      } else {
        temp.DEG.files <- list("Differential_expression_results.csv" = 
                                 deg$result)
      }
      temp.omics.files[["Differential_expression"]] <- temp.DEG.files
      
      # GSEA
      if (grepl("phospho", names(temp.omics)[j], ignore.case = TRUE)) {
        names(gmt2) <- c("phospho_ksdb", "phospho_sub")
        gsea2 <- panSEA::mGSEA(list(deg$result, deg$result), gmt2, 
                               types = names(gmt2), 
                               feature.names = rep(temp.features[j],2),
                               rank.var = rep("Pearson.est",2))
        # store GSEA results
        phospho.gsea.files <- list()
        for (k in 1:length(gmt2)) {
          gsea.name <- paste0("GSEA_", names(gmt2)[k])
          phospho.mtn <- get_top_mtn_plots(
            gsea2$all.results[[k]], 
            EA.type = names(gmt2)[k])
          if (length(phospho.mtn) > 1) {
            phospho.net <- panSEA::netSEA(list(deg$result),
                                          list(gsea2$all.results[[k]]$result),
                                          element.names = temp.features[j],
                                          rank.var = "Pearson.est",
                                          n.network.sets = n.top.sets)
            phospho.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                      gsea2$all.results[[k]]$result,
                                                    "GSEA_volcano_plot.pdf" =
                                                      gsea2$all.results[[k]]$volcano.plot,
                                                    "GSEA_network_graph.html" = 
                                                      phospho.net$interactive,
                                                    "mtn_plots" = phospho.mtn)
          } else {
            phospho.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                      gsea2$all.results[[k]]$result,
                                                    "GSEA_volcano_plot.pdf" =
                                                      gsea2$all.results[[k]]$volcano.plot,
                                                    "mtn_plots" = phospho.mtn)
          } 
          if (length(phospho.mtn) > 0) {
            phospho.gsea.files[[gsea.name]][["Pathways_of_interest"]] <- 
              get_pathways_of_interest(temp.omics[[names(temp.omics)[j]]], 
                                       gsea2$all.results[[k]]$result, 
                                       gmt2[[k]], cc.df.global, n = n.top.sets,
                                       show_colnames = FALSE, 
                                       fontsize = 10, scale = TRUE, cluster = TRUE) 
          }
        }
        temp.omics.files[["Phospho_enrichment"]] <- phospho.gsea.files
        
        # prep for GSEA
        gsea1.inputs <- list(gsea2$all.results[[2]]$result)
        names(gsea1.inputs) <- names(gmt2)[2]
        features1 <- "Feature_set"
        rank.var <- "NES"
      } else {
        # prep for GSEA
        gsea1.inputs <- list(deg$result)
        names(gsea1.inputs) <- names(temp.omics)[j]
        features1 <- temp.features[j]
        rank.var <- "Pearson.est"
      }
      
      # run GSEA
      gsea1 <- list()
      all.global.gsea.files <- list()
      for (k in 1:length(gmt1)) {
        gsea.name <- paste0("GSEA_", names(gmt1)[k])
        gsea1[[gsea.name]] <- panSEA::mGSEA(gsea1.inputs, list(gmt1[[k]]), 
                                            types = names(gsea1.inputs), 
                                            feature.names = rep(features1,1),
                                            rank.var = rep(rank.var,1)) 
        # cannot run GSEA on KSEA results
        global.gsea.files <- list()
        for (m in 1:length(gsea1.inputs)) {
          temp.gsea.files <- list()
          global.mtn <- get_top_mtn_plots(
            gsea1[[gsea.name]]$all.results[[m]], 
            EA.type = names(gmt1)[k])
          temp.gsea.files <- list("GSEA_results.csv" =
                                    gsea1[[gsea.name]]$all.results[[m]]$result,
                                  "GSEA_volcano_plot.pdf" =
                                    gsea1[[gsea.name]]$all.results[[m]]$volcano.plot,
                                  "mtn_plots" = global.mtn)
          if (length(global.mtn) > 1) {
            temp.gsea.files[["GSEA_network_graph.html"]] <- panSEA::netSEA(list(gsea1.inputs[[m]]),
                                                                           list(gsea1[[gsea.name]]$all.results[[m]]$result),
                                                                           element.names = features1,
                                                                           rank.var = rank.var,
                                                                           n.network.sets = n.top.sets)$interactive
          }
          if (!grepl("phospho", names(temp.omics)[j], ignore.case = TRUE) & length(global.mtn) > 0) {
            temp.gsea.files[["Pathways_of_interest"]] <-
              get_pathways_of_interest(temp.omics[[names(temp.omics)[j]]],
                                       gsea1[[gsea.name]]$all.results[[m]]$result,
                                       gmt1[[k]], cc.df.global, n = n.top.sets,
                                       show_colnames = FALSE,
                                       fontsize = 10, scale = TRUE, cluster = TRUE)
          }
          global.gsea.files[[names(gsea1.inputs)[m]]] <- temp.gsea.files
        }
        all.global.gsea.files[[gsea.name]] <- global.gsea.files
      }
      temp.omics.files[["GSEA"]] <- all.global.gsea.files
      
      # DMEA
      if (names(omics)[i] == "Proteomics") {
        temp.expr <- read.csv(file.path(base.path, "CCLE_proteomics.csv"))
        temp.expr <- list(temp.expr)
      } else {
        temp.expr <- list("adherent CCLE")
      }
      dmea.results <- panSEA::mDMEA(expression = temp.expr, weights = gsea1.inputs, 
                                    types = names(gsea1.inputs), 
                                    feature.names = features1,
                                    weight.values = rank.var)
      DMEA.global.mtn <- get_top_mtn_plots(dmea.results$all.results[[1]],
                                           sets = "Drug_set",
                                           EA.type = "DMEA")
      DMEA.global.net <- panSEA::netSEA(list(dmea.results$all.results[[1]]$corr.result),
                                        list(dmea.results$all.results[[1]]$result),
                                        "Drug", "Pearson.est",
                                        n.network.sets = n.top.sets)
      global.DMEA.files <- list("DMEA_results.csv" =
                                  dmea.results$all.results[[1]]$result,
                                "DMEA_correlation_results.csv" = 
                                  dmea.results$all.results[[1]]$corr.result,
                                "DMEA_correlation_scatter_plots.pdf" = 
                                  dmea.results$all.results[[1]]$corr.scatter.plots,
                                "DMEA_volcano_plot.pdf" =
                                  dmea.results$all.results[[1]]$volcano.plot,
                                "DMEA_network_graph.html" = 
                                  DMEA.global.net$interactive,
                                "mtn_plots" = DMEA.global.mtn)
      temp.omics.files[["DMEA"]] <- global.DMEA.files
      omics.files[[names(temp.omics)[j]]] <- temp.omics.files
    }
    save_to_synapse(omics.files, syn.id)
  }
}

panSEA_corr2 <- function(omics, meta.list, feature.list, rank.col = "Gain of C8",
                        other.annotations = c("Ave CEP8", "Ave CMYC"),
                        expr.list, temp.path, syn.id = NULL, 
                        gmt1 = get_gmt1_v2(), gmt2 = NULL,
                        n.top = 50, n.top.sets = 10) {
  for (i in 1:length(omics)) {
    setwd(file.path(temp.path))
    dir.create(names(omics)[i])
    setwd(names(omics)[i])
    temp.features <- feature.list[[i]]
    temp.omics <- omics[[i]]
    meta.df <- as.data.frame(meta.list[[i]])
    rownames(meta.df) <- meta.df[,1]
    meta.df$Sample <- rownames(meta.df)
    red.meta <- meta.df[,c("Sample", rank.col)]
    cc.df.global <- meta.df[,c("Sample", other.annotations, rank.col)]
    
    temp.omics.files <- list()
    omics.t <- as.data.frame(t(temp.omics))
    colnames(omics.t) <- omics.t[1,]
    omics.t <- omics.t[2:nrow(omics.t),]
    omics.t$Sample <- rownames(omics.t)
    # need to make sample names first col, gene symbols colnames
    omics.t <- merge(red.meta, omics.t)
    
    deg <- DMEA::rank_corr(omics.t, variable=temp.features, 
                           value = names(omics)[i], plots = FALSE)
    sig.degs <- na.omit(deg$result[deg$result$Pearson.q <= 0.05, ])
    
    # select numeric data
    just.expr <- dplyr::select_if(temp.omics, is.numeric)
    
    # identify feature name
    feature.name <- temp.features
    
    if (nrow(sig.degs) > 0) {
      top.top.sig.degs <- sig.degs %>% slice_max(Pearson.est, n = n.top/2)
      top.bot.sig.degs <- sig.degs %>% slice_min(Pearson.est, n = n.top/2)
      top.sig.degs <- rbind(top.top.sig.degs, top.bot.sig.degs)
      heatmap.df <- temp.omics[temp.omics[,feature.name] %in% top.sig.degs[,feature.name],
                                    c(feature.name, colnames(just.expr))]
      rownames(heatmap.df) <- heatmap.df[,feature.name]
      feature.order <- top.sig.degs[order(top.sig.degs$Pearson.est, decreasing = TRUE), feature.name]
      heatmap.df <- heatmap.df[feature.order,]
      heatmap.mat <- heatmap.df[,2:ncol(heatmap.df)]
      heatmap.mat <- filter_for_hclust(heatmap.mat)
      
      # create heatmaps
      if (nrow(heatmap.mat) > 1) {
        heatmap.mat <- as.matrix(heatmap.mat)
        deg.heatmap.clust <- pheatmap::pheatmap(heatmap.mat, color = 
                                                  colorRampPalette(
                                                    c("navy", "white", "firebrick3"))(50),
                                                scale = "row", annotation_col = cc.df.global, 
                                                angle_col = "45", 
                                                show_colnames = FALSE, 
                                                fontsize = 10)
        deg.heatmap <- pheatmap::pheatmap(heatmap.mat, 
                                          color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                          cluster_row = FALSE, 
                                          scale = "row", annotation_col = cc.df.global, 
                                          angle_col = "45", 
                                          show_colnames = FALSE, 
                                          fontsize = 10)
        deg.heatmap.abs <- pheatmap::pheatmap(heatmap.mat, color = 
                                                colorRampPalette(
                                                  c("navy", "white", "firebrick3"))(50), 
                                              annotation_col = cc.df.global, 
                                              angle_col = "45", 
                                              show_colnames = FALSE, 
                                              fontsize = 10)
      } else {
        deg.heatmap.clust <- list()
        deg.heatmap <- list()
        deg.heatmap.abs <- list()
      }
      
      temp.DEG.files <- list("Differential_expression_results.csv" = 
                               deg$result,
                             "Differential_expression_results_max_5_percent_FDR.csv" = 
                               sig.degs,
                             "Differential_expression_for_heatmap.csv" =
                               heatmap.df,
                             "Differential_expression_heatmap_not_scaled.bp" =
                               deg.heatmap.abs,
                             "Differential_expression_heatmap_scaled_ordered_by_Log2FC.bp" =
                               deg.heatmap,
                             "Differential_expression_heatmap_scaled.bp" =
                               deg.heatmap.clust) 
    } else {
      temp.DEG.files <- list("Differential_expression_results.csv" = 
                               deg$result)
    }
    temp.omics.files[["Differential_expression"]] <- temp.DEG.files
    
    # GSEA
    if (grepl("phospho", names(omics)[i], ignore.case = TRUE)) {
      names(gmt2) <- c("phospho_ksdb", "phospho_sub")
      gsea2 <- panSEA::mGSEA(list(deg$result, deg$result), gmt2, 
                             types = names(gmt2), 
                             feature.names = rep(temp.features,2),
                             rank.var = rep("Pearson.est",2))
      # store GSEA results
      phospho.gsea.files <- list()
      for (k in 1:length(gmt2)) {
        gsea.name <- paste0("GSEA_", names(gmt2)[k])
        phospho.mtn <- get_top_mtn_plots(
          gsea2$all.results[[k]], 
          EA.type = names(gmt2)[k])
        if (length(phospho.mtn) > 1) {
          phospho.net <- panSEA::netSEA(list(deg$result),
                                        list(gsea2$all.results[[k]]$result),
                                        element.names = temp.features,
                                        rank.var = "Pearson.est",
                                        n.network.sets = n.top.sets)
          phospho.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                    gsea2$all.results[[k]]$result,
                                                  "GSEA_volcano_plot.pdf" =
                                                    gsea2$all.results[[k]]$volcano.plot,
                                                  "GSEA_network_graph.html" = 
                                                    phospho.net$interactive,
                                                  "mtn_plots" = phospho.mtn)
        } else {
          phospho.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                    gsea2$all.results[[k]]$result,
                                                  "GSEA_volcano_plot.pdf" =
                                                    gsea2$all.results[[k]]$volcano.plot,
                                                  "mtn_plots" = phospho.mtn)
        } 
        if (length(phospho.mtn) > 0) {
          phospho.gsea.files[[gsea.name]][["Pathways_of_interest"]] <- 
            get_pathways_of_interest(temp.omics, 
                                     gsea2$all.results[[k]]$result, 
                                     gmt2[[k]], cc.df.global, n = n.top.sets,
                                     show_colnames = FALSE, 
                                     fontsize = 10, scale = TRUE, cluster = TRUE) 
        }
      }
      temp.omics.files[["Phospho_enrichment"]] <- phospho.gsea.files
      
      # prep for GSEA
      gsea1.inputs <- list(gsea2$all.results[[2]]$result)
      names(gsea1.inputs) <- names(gmt2)[2]
      features1 <- "Feature_set"
      rank.var <- "NES"
    } else {
      # prep for GSEA
      gsea1.inputs <- list(deg$result)
      names(gsea1.inputs) <- names(omics)[i]
      features1 <- temp.features
      rank.var <- "Pearson.est"
    }
    
    # run GSEA
    gsea1 <- list()
    all.global.gsea.files <- list()
    for (k in 1:length(gmt1)) {
      gsea.name <- paste0("GSEA_", names(gmt1)[k])
      gsea1[[gsea.name]] <- panSEA::mGSEA(gsea1.inputs, list(gmt1[[k]]), 
                                          types = names(gsea1.inputs), 
                                          feature.names = rep(features1,1),
                                          rank.var = rep(rank.var,1)) 
      # cannot run GSEA on KSEA results
      global.gsea.files <- list()
      for (m in 1:length(gsea1.inputs)) {
        temp.gsea.files <- list()
        global.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[[m]], 
          EA.type = names(gmt1)[k])
        temp.gsea.files <- list("GSEA_results.csv" =
                                  gsea1[[gsea.name]]$all.results[[m]]$result,
                                "GSEA_volcano_plot.pdf" =
                                  gsea1[[gsea.name]]$all.results[[m]]$volcano.plot,
                                "mtn_plots" = global.mtn)
        if (length(global.mtn) > 1) {
          temp.gsea.files[["GSEA_network_graph.html"]] <- panSEA::netSEA(list(gsea1.inputs[[m]]),
                                                                         list(gsea1[[gsea.name]]$all.results[[m]]$result),
                                                                         element.names = features1,
                                                                         rank.var = rank.var,
                                                                         n.network.sets = n.top.sets)$interactive
        }
        if (!grepl("phospho", names(omics)[i], ignore.case = TRUE) & length(global.mtn) > 0) {
          temp.gsea.files[["Pathways_of_interest"]] <-
            get_pathways_of_interest(temp.omics,
                                     gsea1[[gsea.name]]$all.results[[m]]$result,
                                     gmt1[[k]], cc.df.global, n = n.top.sets,
                                     show_colnames = FALSE,
                                     fontsize = 10, scale = TRUE, cluster = TRUE)
        }
        global.gsea.files[[names(gsea1.inputs)[m]]] <- temp.gsea.files
      }
      all.global.gsea.files[[gsea.name]] <- global.gsea.files
    }
    temp.omics.files[["GSEA"]] <- all.global.gsea.files
    
    # DMEA
    if (names(omics)[i] == "Proteomics") {
      temp.expr <- read.csv(file.path(base.path, "CCLE_proteomics.csv"))
      temp.expr <- list(temp.expr)
    } else {
      temp.expr <- list("adherent CCLE")
    }
    dmea.results <- panSEA::mDMEA(expression = temp.expr, weights = gsea1.inputs, 
                                  types = names(gsea1.inputs), 
                                  feature.names = features1,
                                  weight.values = rank.var)
    DMEA.global.mtn <- get_top_mtn_plots(dmea.results$all.results[[1]],
                                         sets = "Drug_set",
                                         EA.type = "DMEA")
    DMEA.global.net <- panSEA::netSEA(list(dmea.results$all.results[[1]]$corr.result),
                                      list(dmea.results$all.results[[1]]$result),
                                      "Drug", "Pearson.est",
                                      n.network.sets = n.top.sets)
    global.DMEA.files <- list("DMEA_results.csv" =
                                dmea.results$all.results[[1]]$result,
                              "DMEA_correlation_results.csv" = 
                                dmea.results$all.results[[1]]$corr.result,
                              "DMEA_correlation_scatter_plots.pdf" = 
                                dmea.results$all.results[[1]]$corr.scatter.plots,
                              "DMEA_volcano_plot.pdf" =
                                dmea.results$all.results[[1]]$volcano.plot,
                              "DMEA_network_graph.html" = 
                                DMEA.global.net$interactive,
                              "mtn_plots" = DMEA.global.mtn)
    temp.omics.files[["DMEA"]] <- global.DMEA.files
    save_to_synapse(temp.omics.files, syn.id)
  }
}

# panSEA_corr2 but permuting ties in enrichment analyses and allowing option for rank metric
panSEA_corr3 <- function(omics, meta.list, feature.list, rank.col = "Gain of C8",
                         other.annotations = c("Ave CEP8", "Ave CMYC"), gene.rank.val = "Spearman.est",
                         gene.sig="Spearman.q",
                         expr.list, temp.path, syn.id = NULL, 
                         gmt1 = get_gmt1_v2(), gmt2 = NULL,
                         n.top = 50, n.top.sets = 10, ties=TRUE, timeout = 300) {
  DEG.forCompile <- list()
  GSEA.forCompile <- list()
  for (i in 1:length(gmt1)) {
    gsea.name <- paste0("GSEA_", names(gmt1)[i])
    GSEA.forCompile[[gsea.name]] <- list()
  }
  DMEA.forCompile <- list()
  files.for.Synapse <- list()
  for (i in 1:length(omics)) {
    setwd(file.path(temp.path))
    dir.create(names(omics)[i])
    setwd(names(omics)[i])
    
    temp.features <- feature.list[[i]]
    temp.omics <- omics[[i]]
    meta.df <- as.data.frame(meta.list[[i]])
    rownames(meta.df) <- meta.df[,1]
    meta.df$Sample <- rownames(meta.df)
    red.meta <- meta.df[,c("Sample", rank.col)]
    cc.df.global <- meta.df[,c("Sample", other.annotations, rank.col)]
    
    if (is.data.frame(temp.omics)) {
      temp.omics.files <- list()
      omics.t <- as.data.frame(t(temp.omics))
      colnames(omics.t) <- omics.t[1,]
      omics.t <- omics.t[2:nrow(omics.t),]
      omics.t$Sample <- rownames(omics.t)
      # need to make sample names first col, gene symbols colnames
      omics.t <- merge(red.meta, omics.t)
      
      deg <- DMEA::rank_corr(omics.t, variable=temp.features, 
                             value = names(omics)[i], plots = FALSE)
      DEG.forCompile[[names(omics)[i]]] <- deg
      sig.degs <- na.omit(deg$result[deg$result[,gene.sig] <= 0.05, ])
      sig.degs$rank.val <- sig.degs[,gene.rank.val]
      
      # select numeric data
      non.features <- colnames(temp.omics)[!(colnames(temp.omics) %in% temp.features)]
      just.expr <- temp.omics[,non.features]
      
      # identify feature name
      feature.name <- temp.features
      
      if (nrow(sig.degs) > 0) {
        top.top.sig.degs <- sig.degs %>% slice_max(rank.val, n = n.top/2)
        top.bot.sig.degs <- sig.degs %>% slice_min(rank.val, n = n.top/2)
        top.sig.degs <- rbind(top.top.sig.degs, top.bot.sig.degs)
        heatmap.df <- temp.omics[temp.omics[,feature.name] %in% top.sig.degs[,feature.name],
                                 c(feature.name, colnames(just.expr))]
        rownames(heatmap.df) <- heatmap.df[,feature.name]
        feature.order <- top.sig.degs[order(top.sig.degs$rank.val, decreasing = TRUE), feature.name]
        heatmap.df <- heatmap.df[feature.order,]
        heatmap.mat <- heatmap.df[,2:ncol(heatmap.df)]
        heatmap.mat <- filter_for_hclust(heatmap.mat)
        
        # create heatmaps
        if (nrow(heatmap.mat) > 1) {
          heatmap.mat <- as.matrix(heatmap.mat)
          deg.heatmap.clust <- try(R.utils::withTimeout(pheatmap::pheatmap(heatmap.mat, color = 
                                                    colorRampPalette(
                                                      c("navy", "white", "firebrick3"))(50),
                                                  scale = "row", annotation_col = cc.df.global, 
                                                  angle_col = "45", 
                                                  show_colnames = FALSE, 
                                                  fontsize = 10), timeout = timeout, onTimeout="error"), silent = TRUE)
          if (inherits(deg.heatmap.clust, "try-error")) {
            deg.heatmap.clust <- list()
          }
          deg.heatmap <- try(R.utils::withTimeout(pheatmap::pheatmap(heatmap.mat, 
                                            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                            cluster_row = FALSE, 
                                            scale = "row", annotation_col = cc.df.global, 
                                            angle_col = "45", 
                                            show_colnames = FALSE, 
                                            fontsize = 10), timeout = timeout, onTimeout="error"), silent = TRUE)
          if (inherits(deg.heatmap, "try-error")) {
            deg.heatmap <- list()
          }
          deg.heatmap.abs <- try(R.utils::withTimeout(pheatmap::pheatmap(heatmap.mat, color = 
                                                  colorRampPalette(
                                                    c("navy", "white", "firebrick3"))(50), 
                                                annotation_col = cc.df.global, 
                                                angle_col = "45", 
                                                show_colnames = FALSE, 
                                                fontsize = 10), timeout = timeout, onTimeout="error"), silent = TRUE)
          if (inherits(deg.heatmap.abs, "try-error")) {
            deg.heatmap.abs <- list()
          }
        } else {
          deg.heatmap.clust <- list()
          deg.heatmap <- list()
          deg.heatmap.abs <- list()
        }
        
        temp.DEG.files <- list("Differential_expression_results.csv" = 
                                 deg$result,
                               "Differential_expression_results_max_5_percent_FDR.csv" = 
                                 sig.degs,
                               "Differential_expression_for_heatmap.csv" =
                                 heatmap.df,
                               "Differential_expression_heatmap_not_scaled.bp" =
                                 deg.heatmap.abs,
                               "Differential_expression_heatmap_scaled_ordered_by_Log2FC.bp" =
                                 deg.heatmap,
                               "Differential_expression_heatmap_scaled.bp" =
                                 deg.heatmap.clust) 
      } else {
        temp.DEG.files <- list("Differential_expression_results.csv" = 
                                 deg$result)
      }
      temp.omics.files[["Differential_expression"]] <- temp.DEG.files
      
      # GSEA
      if (grepl("phospho", names(omics)[i], ignore.case = TRUE)) {
        names(gmt2) <- c("phospho_ksdb", "phospho_sub")
        if (check_coverage(deg$result, gmt2[[1]], temp.features, gene.rank.val) &
            check_coverage(deg$result, gmt2[[2]], temp.features, gene.rank.val)) {
        gsea2 <- panSEA::mGSEA(list(deg$result, deg$result), gmt2, 
                               types = names(gmt2), 
                               feature.names = rep(temp.features,2),
                               rank.var = rep(gene.rank.val,2), ties=ties)
        # store GSEA results
        phospho.gsea.files <- list()
        for (k in 1:length(gmt2)) {
          gsea.name <- paste0("GSEA_", names(gmt2)[k])
          phospho.mtn <- get_top_mtn_plots(
            gsea2$all.results[[k]], 
            EA.type = names(gmt2)[k])
          if (length(phospho.mtn) > 1) {
            if (ties) {
              phospho.net <- try(R.utils::withTimeout(panSEA::netSEA(list(deg$result),
                                            list(gsea2$all.results[[k]]$result.w.ties),
                                            element.names = temp.features,
                                            rank.var = gene.rank.val,
                                            n.network.sets = n.top.sets), timeout = timeout, onTimeout="error"), silent = TRUE)
              if (inherits(phospho.net, "try-error")) {
                phospho.net <- list()
              }
            } else {
              phospho.net <- try(R.utils::withTimeout(panSEA::netSEA(list(deg$result),
                                            list(gsea2$all.results[[k]]$result),
                                            element.names = temp.features,
                                            rank.var = gene.rank.val,
                                            n.network.sets = n.top.sets), timeout = timeout, onTimeout="error"), silent = TRUE)
              if (inherits(phospho.net, "try-error")) {
                phospho.net <- list()
              }
            }
            phospho.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                      gsea2$all.results[[k]]$result,
                                                    "GSEA_volcano_plot.pdf" =
                                                      gsea2$all.results[[k]]$volcano.plot,
                                                    "GSEA_network_graph.html" = 
                                                      phospho.net$interactive,
                                                    "mtn_plots" = phospho.mtn)
          } else {
            phospho.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                      gsea2$all.results[[k]]$result,
                                                    "GSEA_volcano_plot.pdf" =
                                                      gsea2$all.results[[k]]$volcano.plot,
                                                    "mtn_plots" = phospho.mtn)
          } 
          if (ties) {
            phospho.gsea.files[[gsea.name]][["GSEA_results_withoutShufflingTies.csv"]] <- gsea2$all.results[[k]]$result.w.ties
            phospho.gsea.files[[gsea.name]][["GSEA_bar_plot.pdf"]] <- gsea2$all.results[[k]]$bar.plot
            phospho.gsea.files[[gsea.name]][["GSEA_dot_plot.pdf"]] <- gsea2$all.results[[k]]$dot.plot
            phospho.gsea.files[[gsea.name]][["GSEA_dot_plot_withSD.pdf"]] <- gsea2$all.results[[k]]$dot.sd
          }
          if (length(phospho.mtn) > 0) {
            if (ties) {
              phospho.gsea.files[[gsea.name]][["Pathways_of_interest"]] <- 
                get_pathways_of_interest(temp.omics, 
                                         gsea2$all.results[[k]]$result.w.ties, 
                                         gmt2[[k]], cc.df.global, n = n.top.sets,
                                         show_colnames = FALSE, 
                                         fontsize = 10, scale = TRUE, cluster = TRUE)
            } else {
              phospho.gsea.files[[gsea.name]][["Pathways_of_interest"]] <- 
                get_pathways_of_interest(temp.omics, 
                                         gsea2$all.results[[k]]$result, 
                                         gmt2[[k]], cc.df.global, n = n.top.sets,
                                         show_colnames = FALSE, 
                                         fontsize = 10, scale = TRUE, cluster = TRUE) 
            }
          }
        }
        temp.omics.files[["Phospho_enrichment"]] <- phospho.gsea.files
        
        # prep for GSEA
        gsea1.inputs <- list(gsea2$all.results[[2]]$result)
        names(gsea1.inputs) <- names(gmt2)[2]
        features1 <- "Feature_set"
        rank.var <- "NES"
        } else {
          gsea2.result <- data.frame()
          gsea1.inputs <- list(gsea2.result)
          features1 <- "Feature_set"
          rank.var <- "NES"
        }
      } else {
        # prep for GSEA
        gsea1.inputs <- list(deg$result)
        names(gsea1.inputs) <- names(omics)[i]
        features1 <- temp.features
        rank.var <- gene.rank.val
      }
      
      # run GSEA
      gsea1 <- list()
      all.global.gsea.files <- list()
      for (k in 1:length(gmt1)) {
        gsea.name <- paste0("GSEA_", names(gmt1)[k])
        if (check_coverage(gsea1.inputs[[1]], gmt1[[k]], features1, rank.var)) {
          gsea1[[gsea.name]] <- panSEA::mGSEA(gsea1.inputs, list(gmt1[[k]]), 
                                              types = names(gsea1.inputs), 
                                              feature.names = rep(features1,1),
                                              rank.var = rep(rank.var,1), ties=ties)
          GSEA.forCompile[[gsea.name]][[names(omics)[i]]] <- gsea1[[gsea.name]]$all.results[[1]]
          global.gsea.files <- list()
          for (m in 1:length(gsea1.inputs)) {
            temp.gsea.files <- list()
            global.mtn <- get_top_mtn_plots(
              gsea1[[gsea.name]]$all.results[[m]], 
              EA.type = names(gmt1)[k])
            temp.gsea.files <- list("GSEA_results.csv" =
                                      gsea1[[gsea.name]]$all.results[[m]]$result,
                                    "GSEA_volcano_plot.pdf" =
                                      gsea1[[gsea.name]]$all.results[[m]]$volcano.plot,
                                    "mtn_plots" = global.mtn)
            if (ties) {
              temp.gsea.files[["GSEA_results_withoutShufflingTies.csv"]] <- gsea1[[gsea.name]]$all.results[[m]]$result.w.ties
              temp.gsea.files[["GSEA_bar_plot.pdf"]] <- gsea1[[gsea.name]]$all.results[[m]]$bar.plot
              temp.gsea.files[["GSEA_dot_plot.pdf"]] <- gsea1[[gsea.name]]$all.results[[m]]$dot.plot
              temp.gsea.files[["GSEA_dot_plot_withSD.pdf"]] <- gsea1[[gsea.name]]$all.results[[m]]$dot.sd
            }
            if (length(global.mtn) > 1) {
              if (ties) {
                temp.gsea.files[["GSEA_network_graph.html"]] <- try(R.utils::withTimeout(panSEA::netSEA(list(gsea1.inputs[[m]]),
                                                                               list(gsea1[[gsea.name]]$all.results[[m]]$result.w.ties),
                                                                               element.names = features1,
                                                                               rank.var = rank.var,
                                                                               n.network.sets = n.top.sets), timeout = timeout, onTimeout="error")$interactive, silent = TRUE)
                if (inherits(temp.gsea.files[["GSEA_network_graph.html"]], "try-error")) {
                  temp.gsea.files[["GSEA_network_graph.html"]] <- list()
                }
              } else {
                temp.gsea.files[["GSEA_network_graph.html"]] <- try(R.utils::withTimeout(panSEA::netSEA(list(gsea1.inputs[[m]]),
                                                                               list(gsea1[[gsea.name]]$all.results[[m]]$result),
                                                                               element.names = features1,
                                                                               rank.var = rank.var,
                                                                               n.network.sets = n.top.sets), timeout = timeout, onTimeout="error")$interactive, silent = TRUE)
                if (inherits(temp.gsea.files[["GSEA_network_graph.html"]], "try-error")) {
                  temp.gsea.files[["GSEA_network_graph.html"]] <- list()
                }
              }
            }
            if (!grepl("phospho", names(omics)[i], ignore.case = TRUE) & length(global.mtn) > 0) {
              if (ties) {
                temp.gsea.files[["Pathways_of_interest"]] <-
                  get_pathways_of_interest(temp.omics,
                                           gsea1[[gsea.name]]$all.results[[m]]$result.w.ties,
                                           gmt1[[k]], cc.df.global, n = n.top.sets,
                                           show_colnames = FALSE,
                                           fontsize = 10, scale = TRUE, cluster = TRUE)
              } else {
                temp.gsea.files[["Pathways_of_interest"]] <-
                  get_pathways_of_interest(temp.omics,
                                           gsea1[[gsea.name]]$all.results[[m]]$result,
                                           gmt1[[k]], cc.df.global, n = n.top.sets,
                                           show_colnames = FALSE,
                                           fontsize = 10, scale = TRUE, cluster = TRUE)
              }
            }
          }
        } else {
          temp.gsea.files <- list()
        }
        all.global.gsea.files[[gsea.name]] <- temp.gsea.files
      }
      if (length(all.global.gsea.files) > 1) {
        temp.compiled.GSEA <- try(R.utils::withTimeout(panSEA::compile_mGSEA(all.global.gsea.files), timeout = timeout, onTimeout="error"), silent = TRUE)
        if (!inherits(temp.compiled.GSEA, "try-error")) {
          temp.files <- list("GSEA_results.csv" = temp.compiled.GSEA$results,
                             "Compiled_GSEA_results.csv" = temp.compiled.GSEA$mean.results,
                             "GSEA_venn_diagram.pdf" = temp.compiled.GSEA$venn.diagram,
                             "GSEA_dot_plot.pdf" = temp.compiled.GSEA$dot.plot,
                             "GSEA_correlations.csv" = temp.compiled.GSEA$corr,
                             "GSEA_correlation_matrix.pdf" = temp.compiled.GSEA$corr.matrix)
          all.global.gsea.files[["Compiled_results"]] <- temp.files 
        } 
      }
      temp.omics.files[["GSEA"]] <- all.global.gsea.files
      
      # DMEA
      if (!grepl("phospho", names(omics)[i], ignore.case = TRUE)) {
        temp.expr <- expr.list[[i]]
        if (temp.expr == "CCLE proteomics") {
          temp.expr <- read.csv(file.path(base.path, "CCLE_proteomics.csv"))
          temp.expr <- list(temp.expr)
        } else {
          temp.expr <- list(temp.expr)
        }
        dmea.results <- panSEA::mDMEA(expression = temp.expr, weights = gsea1.inputs, 
                                      types = names(gsea1.inputs), 
                                      feature.names = features1,
                                      weight.values = rank.var, ties=ties)
        DMEA.global.mtn <- get_top_mtn_plots(dmea.results$all.results[[1]],
                                             sets = "Drug_set",
                                             EA.type = "DMEA")
        if (ties) {
          DMEA.global.net <- try(R.utils::withTimeout(panSEA::netSEA(list(dmea.results$all.results[[1]]$corr.result),
                                            list(dmea.results$all.results[[1]]$result.w.ties),
                                            "Drug", "Pearson.est",
                                            n.network.sets = n.top.sets), timeout = timeout, onTimeout="error"), silent = TRUE)
          if (inherits(DMEA.global.net, "try-error")) {
            DMEA.global.net <- list()
          }
          global.DMEA.files <- list("DMEA_results.csv" =
                                      dmea.results$all.results[[1]]$result,
                                    "DMEA_results_withoutShufflingTies.csv" =
                                      dmea.results$all.results[[1]]$result.w.ties,
                                    "DMEA_correlation_results.csv" = 
                                      dmea.results$all.results[[1]]$corr.result,
                                    "DMEA_correlation_scatter_plots.pdf" = 
                                      dmea.results$all.results[[1]]$corr.scatter.plots,
                                    "DMEA_volcano_plot.pdf" =
                                      dmea.results$all.results[[1]]$volcano.plot,
                                    "DMEA_network_graph.html" = 
                                      DMEA.global.net$interactive,
                                    "mtn_plots" = DMEA.global.mtn)
        } else {
          DMEA.global.net <- try(R.utils::withTimeout(panSEA::netSEA(list(dmea.results$all.results[[1]]$corr.result),
                                            list(dmea.results$all.results[[1]]$result),
                                            "Drug", "Pearson.est",
                                            n.network.sets = n.top.sets), timeout = timeout, onTimeout="error"), silent = TRUE)
          if (inherits(DMEA.global.net, "try-error")) {
            DMEA.global.net <- list()
          }
          global.DMEA.files <- list("DMEA_results.csv" =
                                      dmea.results$all.results[[1]]$result,
                                    "DMEA_correlation_results.csv" = 
                                      dmea.results$all.results[[1]]$corr.result,
                                    "DMEA_correlation_scatter_plots.pdf" = 
                                      dmea.results$all.results[[1]]$corr.scatter.plots,
                                    "DMEA_volcano_plot.pdf" =
                                      dmea.results$all.results[[1]]$volcano.plot,
                                    "DMEA_network_graph.html" = 
                                      DMEA.global.net$interactive,
                                    "mtn_plots" = DMEA.global.mtn)
        }
        temp.omics.files[["DMEA"]] <- global.DMEA.files
        DMEA.forCompile[[names(omics)[i]]] <- dmea.results$all.results[[1]]
      }
      files.for.Synapse[[names(omics)[i]]] <- temp.omics.files
      temp.files <- list(temp.omics.files)
      names(temp.files) <- names(omics)[i]
      save_to_synapse_v2(temp.files, syn.id)
    } else {
      temp.meta <- list()
      for (j in 1:length(temp.omics)) {
        temp.meta[[j]] <- cc.df.global
      }
      panSEA_corr3(omics=temp.omics, meta.list=temp.meta, 
                   feature.list=temp.features, rank.col=rank.col,
                   other.annotations=other.annotations, gene.rank.val=gene.rank.val,
                   expr.list=expr.list[[i]], temp.path = getwd(), syn.id=syn.id, 
                   gmt1=gmt1, gmt2=gmt2,
                   n.top=n.top, n.top.sets=n.top.sets, ties=ties)
    }
  }
  # compile results
  if (length(names(omics)) > 1) {
    combo.files <- list()
    compiled.DEGs <- try(R.utils::withTimeout(panSEA::compile_mDEG(DEG.forCompile), timeout = timeout, onTimeout="error"), silent = TRUE)
    if (!inherits(compiled.DEGs, "try-error")) {
      compiled.DEG.files <- list("Differential_expression_results.csv" = compiled.DEGs$results,
                                 "Compiled_differential_expression_results.csv" = compiled.DEGs$mean.results,
                                 "Differential_expression_venn_diagram.pdf" = compiled.DEGs$venn.diagram,
                                 "Differential_expression_dot_plot.pdf" = compiled.DEGs$dot.plot,
                                 "Differential_expression_correlations.csv" = compiled.DEGs$corr,
                                 "Differential_expression_correlation_matrix.pdf" = compiled.DEGs$corr.matrix)
      combo.files[["Differential_expression"]] <- compiled.DEG.files 
    }
    
    compiled.GSEA.files <- list()
    for (i in 1:length(gmt1)) {
      if (length(GSEA.forCompile[[names(gmt1)[i]]] > 1)) {
        compiled.GSEA <- try(R.utils::withTimeout(panSEA::compile_mGSEA(GSEA.forCompile[[names(gmt1)[i]]]), timeout = timeout, onTimeout="error"), silent = TRUE)
        if (!inherits(compiled.GSEA, "try-error")) {
          temp.files <- list("GSEA_results.csv" = compiled.GSEA$results,
                             "Compiled_GSEA_results.csv" = compiled.GSEA$mean.results,
                             "GSEA_venn_diagram.pdf" = compiled.GSEA$venn.diagram,
                             "GSEA_dot_plot.pdf" = compiled.GSEA$dot.plot,
                             "GSEA_correlations.csv" = compiled.GSEA$corr,
                             "GSEA_correlation_matrix.pdf" = compiled.GSEA$corr.matrix)
          compiled.GSEA.files[[paste0("GSEA_",names(gmt1)[i])]] <- temp.files 
        }
      }
    }
    if (length(compiled.GSEA.files) > 0) {
      combo.files[["GSEA"]] <- compiled.GSEA.files
    }
    
    if (length(DMEA.forCompile) > 1) {
      compiled.DMEA <- try(R.utils::withTimeout(panSEA::compile_mDMEA(DMEA.forCompile), timeout = timeout, onTimeout="error"), silent = TRUE)
      if (!inherits(compiled.DMEA, "try-error")) {
        compiled.DMEA.files <- list("DMEA_results.csv" = compiled.DMEA$results,
                                    "Compiled_DMEA_results.csv" = compiled.DMEA$mean.results,
                                    "DMEA_venn_diagram.pdf" = compiled.DMEA$venn.diagram,
                                    "DMEA_dot_plot.pdf" = compiled.DMEA$dot.plot,
                                    "DMEA_correlations.csv" = compiled.DMEA$corr,
                                    "DMEA_correlation_matrix.pdf" = compiled.DMEA$corr.matrix) 
        combo.files[["DMEA"]] <- compiled.DMEA.files 
      }
    }
    if (length(combo.files) > 0) {
      files.for.Synapse[["Compiled_results"]] <- combo.files
      save_to_synapse_v2(list("Compiled_results" = combo.files), syn.id) 
    }
  }
  #save_to_synapse(files.for.Synapse, syn.id)
}

# just use regular panSEA but add 
panSEA_corr4 <- function(omics, meta.list, feature.list, rank.col = "Gain of C8",
                         other.annotations = c("Ave CEP8", "Ave CMYC"), gene.rank.val = "Spearman.est",
                         gene.sig="Spearman.q",
                         expr.list, temp.path, syn.id = NULL, 
                         gmt1 = get_gmt1_v2(), gmt2 = NULL,
                         n.top = 50, n.top.sets = 10, ties=TRUE) {
  
}

# check set coverage
check_coverage_list <- function(gsea.inputs, gmt.list, 
                                feature.names = rep("Gene", length(gsea.inputs)), 
                                rank.var = rep("Log2FC", length(gsea.inputs)), 
                                min.per.set = 6) {
  # add check to see if lengths of lists are all equal
  input.names <- names(gsea.inputs)
  gmt.names <- names(gmt.list)
  names(gmt.list) <- input.names
  gmt.name.map <- as.list(gmt.names)
  names(gmt.name.map) <- input.names
  
  features <- as.list(feature.names)
  names(features) <- input.names
    
  ranks <- as.list(rank.var)
  names(ranks) <- input.names
  for (i in 1:length(input.names)) {
    if (!check_coverage(gsea.inputs[[input.names[i]]], gmt.list[[input.names[i]]], 
                        feature.names[i], rank.var[i], min.per.set)) {
      gsea.inputs <- gsea.inputs[-which(names(gsea.inputs) == input.names[i])]
      gmt.list <- gmt.list[-which(names(gmt.list) == input.names[i])]
      gmt.name.map <- gmt.name.map[-which(names(gmt.name.map) == input.names[i])]
      features <- features[-which(names(features) == input.names[i])]
      ranks <- ranks[-which(names(ranks) == input.names[i])]
    }
  }
  ranks2 <- unlist(ranks)
  features2 <- unlist(features)
  names(gmt.list) <- unlist(gmt.name.map)
  return(list(gsea.inputs = gsea.inputs, gmt.list = gmt.list, 
              feature.names = features2, rank.var = ranks2))
}

check_coverage <- function(gsea.input, gmt, feature.names = "Gene", 
                           rank.var = "Log2FC", min.per.set = 6) {
  covered <- FALSE
  n.sets <- 0
  for (i in 1:length(gmt$genesets)) {
    temp.genes <- gmt$genesets[[i]]
    temp.gsea.input <- gsea.input[gsea.input[,feature.names] %in% temp.genes, ]
    if (nrow(temp.gsea.input) >= min.per.set) {
      n.sets <- n.sets + 1
    }
  }
  if (n.sets >= 2) {
    covered <- TRUE
  }
  return(covered)
}
