# Differential expression & enrichment analyses: global & phospho
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-04-25

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr)
library(dplyr)

## get gmt information for GSEA
get_gmt1 <- function(gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                   "msigdb_Homo sapiens_H",
                                   "msigdb_Homo sapiens_C1")) {
  if (file.exists("gmt1_run_contrasts_global_phospho_human.rds")) {
    gmt1 <- readRDS("gmt1_run_contrasts_global_phospho_human.rds")
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
  return(gmt1)
}

# get gmt information for GSEA relevant to Chr8
get_chr8_gmt1 <- function(gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                   "msigdb_Homo sapiens_H",
                                   "msigdb_Homo sapiens_C1",
                                   "chr8_cancer_human")) {
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
  }
  return(gmt1)
}

# get kinase-substrate database information in gmt format for KSEA
get_ksdb <- function(organism="human"){
  if (is.character(organism)) {
    if (organism == "human") {
      if (!file.exists("gmt_ksdb_human.rds")) {
        url <- "https://raw.github.com/BelindaBGarana/panSEA/main/data/gmt_ksdb_human.rds"
        httr::GET(url, httr::write_disk("gmt_ksdb_human.rds")) 
      }
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
    } 
  }
  return(gmt)
}

# get gmt information for KSEA/SSEA relevant to chr8
get_chr8_gmt2 <- function() {
  if (file.exists("chr8_gmt2_run_contrasts_global_phospho_human.rds")) {
    gmt2 <- readRDS("chr8_gmt2_run_contrasts_global_phospho_human.rds")
  } else {
    stop("Chr8 gmt2 file not found in current directory") 
  }
}

# get CCLE global proteomics data
get_CCLE_prot <- function() {
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
  return(prot.df.noNA)
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

# save folder and nested files/subfolders to Synapse
save_to_synapse <- function(temp.files, resultsFolder = NULL) {
  CSV.files <- names(temp.files)[grepl(".csv", names(temp.files))]
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
  
  PDF.files <- names(temp.files)[grepl(".pdf", names(temp.files))]
  if (length(PDF.files) > 0) {
    # save locally
    for (j in 1:length(PDF.files)) {
      if (is.list(temp.files[[PDF.files[j]]])) {
        ggplot2::ggsave(PDF.files[j], temp.files[[PDF.files[j]]], 
                        device = "pdf") 
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
  
  HTML.files <- names(temp.files)[grepl(".html", names(temp.files))]
  if (length(HTML.files) > 0) {
    # save locally
    for (j in 1:length(HTML.files)) {
      if (is.list(temp.files[[HTML.files[j]]])) {
        visNetwork::visSave(temp.files[[HTML.files[j]]], HTML.files[j])
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

## run panSEA across global & phospho types for multiple contrasts
# input: vector of contrasts, contrast.type column for meta.df extraction, 
#   omics data.list, beatAML expression list, gmt set information, 
#   BeatAML drug AUC, base.path, Synapse ID info
# output: panSEA result files are saved locally & uploaded to Synapse
run_contrasts_global_phospho_human <- function(contrasts, contrast.type, id.type, meta.df, 
                          omics, gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                              "msigdb_Homo sapiens_H",
                                              "msigdb_Homo sapiens_C1"),
                          EA.types = c("KEGG", "Hallmark", "Positional"),
                          gmt.list2 = c("ksdb_human", "sub"),
                          expr = as.list(rep("adherent CCLE", 3)),
                          gmt.drug = "PRISM", drug.sens = "PRISM", 
                          base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                          temp.path, subfolder = TRUE,
                          synapse_id = NULL) {
  setwd(base.path)
  types = c("global", "phospho")
  feature.names = c("Gene", "SUB_SITE")
  
  # load set annotations
  if (gmt.list1[1] == "chr8") {
    gmt1 <- get_chr8_gmt1()
  } else {
    gmt1 <- get_gmt1(gmt.list1)
  }
  
  if (gmt.list2[1] == "chr8") {
    gmt2 <- get_chr8_gmt2()
  } else {
    gmt2 <- list()
    for (i in 1:length(gmt.list2)) {
      if (is.character(gmt.list2[i])) {
        if (grepl("ksdb", gmt.list2[i], ignore.case = TRUE)) {
          org <- stringr::str_split(gmt.list2[i], "_")[[1]][2]
          gmt2[[i]] <- get_ksdb(organism = org)
        } else if (gmt.list2[i] == "sub") {
          SUB_SITE <- omics[[2]]$SUB_SITE
          phospho.ref <- data.frame(SUB_SITE)
          phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
                                                        remove = FALSE)
          SUB_SITE <- NULL
          gmt2[[i]] <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE")
        }
      } else {
        gmt2[[i]] <- gmt.list2[i]
      }
    }
  }
  
  prot.df.noNA <- get_CCLE_prot()
  
  ## for each contrast: 
  for (k in 1:length(contrasts)) {
    # identify samples for each side of contrast
    c1 <- contrasts[[k]][1]
    c2 <- contrasts[[k]][2]
    group.names <- c(c1, c2)
    contrast.name <- paste0(contrast.type, "_", contrasts[[k]][1], "_vs_", contrasts[[k]][2])
    group.samples <- list(as.vector(na.omit(meta.df[meta.df[,contrast.type] == c1, id.type]))[[id.type]],
                          as.vector(na.omit(meta.df[meta.df[,contrast.type] == c2, id.type]))[[id.type]])
    
    if (length(group.samples[[1]]) > 0 & length(group.samples[[2]]) > 0) {
      # run panSEA across omics types
      # CAUTION: this only works because the # of samples for each treatment type 
      # is equal; otherwise would have to run panSEA for each contrast separately 
      # and perhaps set the group.samples input parameter for panSEA
      
      # run global against KEGG sets, phospho against kinase sets from ksdb
      # run differential expression analysis
      deg <- panSEA::mDEG(omics, types, group.names, group.samples, 
                                   feature.names)$all.results
      
      ## for phospho data:
      # run GSEA for each gmt in gmt2 and then check coverage of 2+ sets in gmt1
      gsea2.inputs <- list(deg[[2]], deg[[2]])
      gsea2 <- panSEA::mGSEA(gsea2.inputs, gmt2, types = c("phospho_ksdb", "phospho_sub"),
                             feature.names = c("SUB_SITE", "SUB_SITE"))
      
      # create network graphs for KSEA, SSEA
      kin.net2 <- panSEA::netSEA(list(gsea2.inputs[[1]]),
                                list(gsea2$all.results[[1]]$result),
                                "SUB_SITE",
                                n.network.sets = 5)
      sub.net2 <- panSEA::netSEA(list(gsea2.inputs[[2]]),
                                 list(gsea2$all.results[[2]]$result),
                                 "SUB_SITE",
                                 n.network.sets = 5)
      
      if (nrow(gsea2$all.results[[1]]$result) > 100 & 
          nrow(gsea2$all.results[[2]]$result) > 100) {
        KSEA <- TRUE
        SSEA <- TRUE
        gsea1.inputs <- list("global" = deg[[1]],
                             "phospho_ksdb" = gsea2$all.results[[1]]$result,
                             "phospho_sub" = gsea2$all.results[[2]]$result) 
        prot.expr <- list(prot.df.noNA, prot.df.noNA, prot.df.noNA)
        features1 <- c("Gene", "Feature_set", "Feature_set")
        rank.var <- c("Log2FC", "NES", "NES")
      } else if (nrow(gsea2$all.results[[1]]$result) > 100) {
        KSEA <- TRUE
        SSEA <- FALSE
        gsea1.inputs <- list("global" = deg[[1]],
                             "phospho_ksdb" = gsea2$all.results[[1]]$result)
        prot.expr <- list(prot.df.noNA, prot.df.noNA)
        features1 <- c("Gene", "Feature_set")
        rank.var <- c("Log2FC", "NES")
      } else if (nrow(gsea2$all.results[[2]]$result) > 100) {
        KSEA <- FALSE
        SSEA <- TRUE
        gsea1.inputs <- list("global" = deg[[1]],
                             "phospho_sub" = gsea2$all.results[[2]]$result) 
        prot.expr <- list(prot.df.noNA, prot.df.noNA)
        features1 <- c("Gene", "Feature_set")
        rank.var <- c("Log2FC", "NES")
      } else {
        KSEA <- FALSE
        SSEA <- FALSE
        gsea1.inputs <- list("global" = deg[[1]])
        prot.expr <- list(prot.df.noNA)
        features1 <- c("Gene")
        rank.var <- c("Log2FC")
      }
      
      # run GSEA for each gmt in gmt1
      gsea1 <- list()
      global.gsea.files <- list()
      kin.gsea.files <- list()
      sub.gsea.files <- list()
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
        
        # store results for each input type
        global.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[["global"]], 
          EA.type = EA.types[i])
        global.net <- panSEA::netSEA(list(gsea1.inputs[["global"]]),
                                                  list(gsea1[[gsea.name]]$all.results[["global"]]$result),
                                                  n.network.sets = 5)
        global.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                 gsea1[[gsea.name]]$all.results[["global"]]$result,
                                               "GSEA_volcano_plot.pdf" =
                                                 gsea1[[gsea.name]]$all.results[["global"]]$volcano.plot,
                                               "GSEA_network_graph.html" = 
                                                 global.net$interactive,
                                               "mtn_plots" = global.mtn)
        if (KSEA) {
          kin.mtn <- get_top_mtn_plots(
            gsea1[[gsea.name]]$all.results[["phospho_ksdb"]], 
            EA.type = EA.types[i])
          kin.net <- panSEA::netSEA(list(gsea1.inputs[["phospho_ksdb"]]),
                                                 list(gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$result),
                                                 "Feature_set", "NES",
                                                 n.network.sets = 5)
          kin.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$result,
                                              "GSEA_volcano_plot.pdf" =
                                                gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$volcano.plot,
                                              "GSEA_network_graph.html" = 
                                                kin.net$interactive,
                                              "mtn_plots" = kin.mtn) 
        }
        
        if (SSEA) {
          sub.mtn <- get_top_mtn_plots(
            gsea1[[gsea.name]]$all.results[["phospho_sub"]], 
            EA.type = EA.types[i])
          sub.net <- panSEA::netSEA(list(gsea1.inputs[["phospho_sub"]]),
                                                 list(gsea1[[gsea.name]]$all.results[["phospho_sub"]]$result),
                                                 "Feature_set", "NES",
                                                 n.network.sets = 5)
          sub.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                gsea1[[gsea.name]]$all.results[["phospho_sub"]]$result,
                                              "GSEA_volcano_plot.pdf" =
                                                gsea1[[gsea.name]]$all.results[["phospho_sub"]]$volcano.plot,
                                              "GSEA_network_graph.html" = 
                                                sub.net$interactive,
                                              "mtn_plots" = sub.mtn) 
        }
      }
      
      # run DMEA
      dmea.results <- panSEA::mDMEA(drug.sens, gmt.drug, expr, gsea1.inputs, 
                                    names(gsea1.inputs), 
                                    feature.names = features1,
                                    weight.values = rank.var)
      dmea.prot.results <- panSEA::mDMEA(drug.sens, gmt.drug, prot.expr, gsea1.inputs, 
                                         names(gsea1.inputs), 
                                         feature.names = features1,
                                         weight.values = rank.var)
      
      #### save results & upload to Synapse
      ### set file names
      ## global
      global.DEG.files <- list("Differential_expression_results.csv" = 
                                 deg[["global"]])
      DMEA.global.mtn <- get_top_mtn_plots(dmea.results$all.results[["global"]],
                                           sets = "Drug_set",
                                           EA.type = "DMEA")
      DMEA.global.net <- panSEA::netSEA(list(dmea.results$all.results[["global"]]$corr.result),
                                list(dmea.results$all.results[["global"]]$result),
                                "Drug", "Pearson.est",
                                n.network.sets = 5)
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
      prot.DMEA.global.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[[1]],
                                                sets = "Drug_set",
                                                EA.type = "DMEA")
      prot.DMEA.global.net <- panSEA::netSEA(list(dmea.prot.results$all.results[[1]]$corr.result),
                                        list(dmea.prot.results$all.results[[1]]$result),
                                        "Drug", "Pearson.est",
                                        n.network.sets = 5)
      prot.global.DMEA.files <- list("DMEA_results.csv" =
                                       dmea.prot.results$all.results[[1]]$result,
                                     "DMEA_correlation_results.csv" = 
                                       dmea.prot.results$all.results[[1]]$corr.result,
                                     "DMEA_correlation_scatter_plots.pdf" = 
                                       dmea.prot.results$all.results[[1]]$corr.scatter.plots,
                                     "DMEA_volcano_plot.pdf" =
                                       dmea.prot.results$all.results[[1]]$volcano.plot,
                                     "DMEA_network_graph.html" = 
                                       prot.DMEA.global.net$interactive,
                                     "mtn_plots" = prot.DMEA.global.mtn)
      global.files <- list('Differential_expression' = global.DEG.files, 
                           'DMEA' = global.DMEA.files,
                           'DMEA_proteomics' = prot.global.DMEA.files,
                           'GSEA' = global.gsea.files)
      
      ## phospho
      phospho.DEG.files <- list("Differential_expression_results.csv" = 
                                  deg[[2]])

      kin.mtn2 <- get_top_mtn_plots(gsea2$all.results[[1]],
                                    EA.type = "KSEA")
      if (KSEA) {
        kin.DMEA.mtn <- get_top_mtn_plots(dmea.results$all.results[[2]],
                                          sets = "Drug_set",
                                          EA.type = "DMEA")
        DMEA.kin.net <- panSEA::netSEA(list(dmea.results$all.results[["phospho_ksdb"]]$corr.result),
                                            list(dmea.results$all.results[["phospho_ksdb"]]$result),
                                            "Drug", "Pearson.est",
                                            n.network.sets = 5)
        kin.DMEA.files <- list("DMEA_results.csv" =
                                 dmea.results$all.results[[2]]$result,
                               "DMEA_correlation_results.csv" = 
                                 dmea.results$all.results[[2]]$corr.result,
                               "DMEA_correlation_scatter_plots.pdf" = 
                                 dmea.results$all.results[[2]]$corr.scatter.plots,
                               "DMEA_volcano_plot.pdf" =
                                 dmea.results$all.results[[2]]$volcano.plot,
                               "DMEA_network_graph.html" = 
                                 DMEA.kin.net$interactive,
                               "mtn_plots" = kin.DMEA.mtn)
        
        prot.kin.DMEA.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[[2]],
                                               sets = "Drug_set",
                                               EA.type = "DMEA")
        prot.DMEA.kin.net <- panSEA::netSEA(list(prot.dmea.results$all.results[["phospho_ksdb"]]$corr.result),
                                            list(prot.dmea.results$all.results[["phospho_ksdb"]]$result),
                                            "Drug", "Pearson.est",
                                            n.network.sets = 5)
        prot.kin.DMEA.files <- list("DMEA_results.csv" =
                                      dmea.prot.results$all.results[[2]]$result,
                                    "DMEA_correlation_results.csv" = 
                                      dmea.prot.results$all.results[[2]]$corr.result,
                                    "DMEA_correlation_scatter_plots.pdf" = 
                                      dmea.prot.results$all.results[[2]]$corr.scatter.plots,
                                    "DMEA_volcano_plot.pdf" =
                                      dmea.prot.results$all.results[[2]]$volcano.plot,
                                    "DMEA_network_graph.html" = 
                                      prot.DMEA.kin.net$interactive,
                                    "mtn_plots" = prot.kin.DMEA.mtn)
      } else {
        kin.DMEA.files <- list()
        prot.kin.DMEA.files <- list()
      }
      phospho.kin.files <- list("KSEA_results.csv" =
                                  gsea2$all.results[[1]]$result,
                                "KSEA_volcano_plot.pdf" =
                                  gsea2$all.results[[1]]$volcano.plot,
                                "KSEA_network_graph.html" = 
                                  kin.net2$interactive,
                                "mtn_plots" = kin.mtn2,
                                "GSEA" = kin.gsea.files,
                                "DMEA" = kin.DMEA.files,
                                "DMEA_proteomics" = prot.kin.DMEA.files) 
      
      
      sub.mtn2 <- get_top_mtn_plots(gsea2$all.results[[2]],
                                    EA.type = "Substrate_enrichment")
      if (SSEA) {
        sub.DMEA.mtn <- get_top_mtn_plots(dmea.results$all.results[["phospho_sub"]],
                                          sets = "Drug_set",
                                          EA.type = "DMEA")
        DMEA.sub.net <- panSEA::netSEA(list(dmea.results$all.results[["phospho_sub"]]$corr.result),
                                       list(dmea.results$all.results[["phospho_sub"]]$result),
                                       "Drug", "Pearson.est",
                                       n.network.sets = 5)
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
        prot.sub.DMEA.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[["phospho_sub"]],
                                               sets = "Drug_set",
                                               EA.type = "DMEA")
        prot.DMEA.sub.net <- panSEA::netSEA(list(dmea.prot.results$all.results[["phospho_sub"]]$corr.result),
                                       list(dmea.prot.results$all.results[["phospho_sub"]]$result),
                                       "Drug", "Pearson.est",
                                       n.network.sets = 5)
        prot.sub.DMEA.files <- list("DMEA_results.csv" =
                                      dmea.prot.results$all.results[["phospho_sub"]]$result,
                                    "DMEA_correlation_results.csv" = 
                                      dmea.prot.results$all.results[["phospho_sub"]]$corr.result,
                                    "DMEA_correlation_scatter_plots.pdf" = 
                                      dmea.prot.results$all.results[["phospho_sub"]]$corr.scatter.plots,
                                    "DMEA_volcano_plot.pdf" =
                                      dmea.prot.results$all.results[["phospho_sub"]]$volcano.plot,
                                    "DMEA_network_graph.html" = 
                                      prot.DMEA.sub.net$interactive,
                                    "mtn_plots" = prot.sub.DMEA.mtn)
      } else {
        sub.DMEA.files <- list()
        prot.sub.DMEA.files <- list()
      }
      phospho.sub.files <- list("Substrate_enrichment_results.csv" =
                                  gsea2$all.results[[2]]$result,
                                "Substrate_enrichment_volcano_plot.pdf" =
                                  gsea2$all.results[[2]]$volcano.plot,
                                "Substrate_enrichment_network_graph.html" = 
                                  sub.net2$interactive,
                                "mtn_plots" = sub.mtn2,
                                "GSEA" = sub.gsea.files,
                                "DMEA" = sub.DMEA.files,
                                "DMEA_proteomics" = prot.sub.DMEA.files) 
      
      phospho.files <- list('Differential_expression' = phospho.DEG.files, 
                            'KSEA' = phospho.kin.files,
                            'Substrate_enrichment' = phospho.sub.files)
      
      ## combo
      combo.DMEA.files <- list("DMEA_results.csv" =
                                 dmea.results$compiled.results$results,
                               "DMEA_mean_results.csv" =
                                 dmea.results$compiled.results$mean.results,
                               "DMEA_correlation_matrix.pdf" =
                                 dmea.results$compiled.results$corr.matrix,
                               "DMEA_dot_plot.pdf" =
                                 dmea.results$compiled.results$dot.plot)
      prot.combo.DMEA.files <- list("DMEA_results.csv" =
                                      dmea.prot.results$compiled.results$results,
                                    "DMEA_mean_results.csv" =
                                      dmea.prot.results$compiled.results$mean.results,
                                    "DMEA_correlation_matrix.pdf" =
                                      dmea.prot.results$compiled.results$corr.matrix,
                                    "DMEA_dot_plot.pdf" =
                                      dmea.prot.results$compiled.results$dot.plot)
      combo.files <- list('DMEA' = combo.DMEA.files,
                          'DMEA_proteomics' = prot.combo.DMEA.files)
      
      all.files <- list('global_and_phospho' = combo.files,
                        'global' = global.files,
                        'phospho' = phospho.files)
      
      # create folder for contrast
      setwd(temp.path)
      if (subfolder) {
        dir.create(contrast.name)
        setwd(contrast.name)
        contrastFolder <- 
          synapser::synStore(synapser::Folder(contrast.name,
                                              parent = synapse_id))
      } else {
        contrastFolder <- synapse_id
        contrast.name <- ""
      }
      
      save_to_synapse(all.files, contrastFolder)
      
      # make space to process next contrast
      gsea1 <- NULL
      gsea2 <- NULL
      dmea.results <- NULL
      deg <- NULL
      gsea1.inputs <- NULL
    }
  }
}

run_contrasts2 <- function(contrast.type, contrast.type2 = NULL, 
                          contrasts = as.list(rep(c(TRUE, FALSE), length(contrast.type))),
                          contrasts2 = NULL, id.type, meta.df, omics, 
                          types = c("global", "phospho"), 
                          feature.names = c("Gene", "SUB_SITE"),
                          gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                        "msigdb_Homo sapiens_H",
                                        "msigdb_Homo sapiens_C1"),
                          EA.types = c("KEGG", "Hallmark", "Positional"),
                          gmt.list2 = c("ksdb_human", "sub"),
                          expr = as.list(rep("adherent CCLE", length(types) + ifelse(grepl("phospho", types, ignore.case = TRUE),length(gmt.list2)-1,0))),
                          gmt.drug = "PRISM", drug.sens = "PRISM", 
                          base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                          temp.path = base.path, subfolder = TRUE, synapse_id = NULL,
                          compileDEGs = TRUE, filter = NA, filterID = NULL) {
  setwd(base.path)
  # load set annotations
  if (gmt.list1[1] == "chr8") {
    gmt1 <- get_chr8_gmt1()
  } else {
    gmt1 <- get_gmt1(gmt.list1)
  }
  
  if (grepl("phospho", types, ignore.case = TRUE)) {
    if (gmt.list2[1] == "chr8") {
      gmt2 <- get_chr8_gmt2()
    } else {
      gmt2 <- list()
      for (i in 1:length(gmt.list2)) {
        if (is.character(gmt.list2[i])) {
          if (grepl("ksdb", gmt.list2[i], ignore.case = TRUE)) {
            org <- stringr::str_split(gmt.list2[i], "_")[[1]][2]
            gmt2[[i]] <- get_ksdb(organism = org)
          } else if (gmt.list2[i] == "sub") {
            SUB_SITE <- omics[[2]]$SUB_SITE
            phospho.ref <- data.frame(SUB_SITE)
            phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
                                                          remove = FALSE)
            SUB_SITE <- NULL
            gmt2[[i]] <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE")
          }
        } else {
          gmt2[[i]] <- gmt.list2[i]
        }
      }
    } 
  }
  
  prot.df.noNA <- get_CCLE_prot()
  all.degs <- data.frame()
  
  ## for each contrast: 
  for (k in 1:length(contrasts)) {
    # identify samples for each side of contrast
    c1 <- contrasts[[k]][1]
    c2 <- contrasts[[k]][2]
    group.names <- c(c1, c2)
    contrast.name <- paste0(contrast.type, "_", contrasts[[k]][1], "_vs_", contrasts[[k]][2])
    group.samples <- list(as.vector(na.omit(meta.df[meta.df[,contrast.type] == c1, id.type]))[[id.type]],
                          as.vector(na.omit(meta.df[meta.df[,contrast.type] == c2, id.type]))[[id.type]])
    
    # identify samples for second contrast if relevant
    if (!is.null(contrast.type2) & !is.null(contrasts2)) {
      c1a <- contrasts2[[k]][1]
      c2a <- contrasts2[[k]][2]
      group.names2 <- c(c1a, c2a)
      contrast.name2 <- paste0(contrast.type2, "_", contrasts2[[k]][1], "_vs_", contrasts2[[k]][2])
      group.samples2 <- list(as.vector(na.omit(meta.df[meta.df[,contrast.type2] == c1a, id.type]))[[id.type]],
                            as.vector(na.omit(meta.df[meta.df[,contrast.type2] == c2a, id.type]))[[id.type]])
      if (length(group.samples2[[1]]) == 0 | length(group.samples2[[2]]) == 0) {
        group.names2 <- NULL
        group.samples2 <- NULL
      }
    } else {
      group.names2 <- NULL
      group.samples2 <- NULL
    }
    
    if (length(group.samples[[1]]) > 0 & length(group.samples[[2]]) > 0) {
      # run panSEA across omics types
      # CAUTION: this only works because the # of samples for each treatment type 
      # is equal; otherwise would have to run panSEA for each contrast separately 
      # and perhaps set the group.samples input parameter for panSEA
      
      # run global against KEGG sets, phospho against kinase sets from ksdb
      # run differential expression analysis
      deg <- panSEA::mDEG(omics, types, group.names, group.samples, 
                          group.names2, group.samples2,
                          feature.names)$all.results
      
      ## for phospho data:
      # run GSEA for each gmt in gmt2 and then check coverage of 2+ sets in gmt1
      gsea1.inputs <- deg
      KSEA <- FALSE
      SSEA <- FALSE
      if (grepl("phospho", names(deg), ignore.case = TRUE)) {
        n.phospho <- grep("phospho", names(deg), ignore.case = TRUE)
        gsea2.inputs <- list(deg[[n.phospho]], deg[[n.phospho]])
        gsea2 <- panSEA::mGSEA(gsea2.inputs, gmt2, types = c("phospho_ksdb", "phospho_sub"),
                               feature.names = rep(feature.names[n.phospho], 2)) 
        
        # create mountain plots and network graphs for KSEA, SSEA
        kin.mtn2 <- get_top_mtn_plots(gsea2$all.results[[1]],
                                      EA.type = "KSEA")
        kin.net2 <- panSEA::netSEA(list(gsea2.inputs[[1]]),
                                   list(gsea2$all.results[[1]]$result),
                                   "SUB_SITE",
                                   n.network.sets = 5)
        sub.mtn2 <- get_top_mtn_plots(gsea2$all.results[[2]],
                                      EA.type = "Substrate_enrichment")
        sub.net2 <- panSEA::netSEA(list(gsea2.inputs[[2]]),
                                   list(gsea2$all.results[[2]]$result),
                                   "SUB_SITE",
                                   n.network.sets = 5)
        
        if (nrow(gsea2$all.results[[1]]$result) > 100 & 
            nrow(gsea2$all.results[[2]]$result) > 100) {
          KSEA <- TRUE
          SSEA <- TRUE
          gsea1.inputs[["phospho_ksdb"]] <- gsea2$all.results[[1]]$result
          gsea1.inputs[["phospho_sub"]] <- gsea2$all.results[[2]]$result 
          prot.expr <- list(prot.df.noNA, prot.df.noNA, prot.df.noNA)
          features1 <- c("Gene", "Feature_set", "Feature_set")
          rank.var <- c("Log2FC", "NES", "NES")
        } else if (nrow(gsea2$all.results[[1]]$result) > 100) {
          KSEA <- TRUE
          gsea1.inputs[["phospho_ksdb"]] <- gsea2$all.results[[1]]$result
          prot.expr <- list(prot.df.noNA, prot.df.noNA)
          features1 <- c("Gene", "Feature_set")
          rank.var <- c("Log2FC", "NES")
        } else if (nrow(gsea2$all.results[[2]]$result) > 100) {
          SSEA <- TRUE
          gsea1.inputs[["phospho_sub"]] <- gsea2$all.results[[2]]$result 
          prot.expr <- list(prot.df.noNA, prot.df.noNA)
          features1 <- c("Gene", "Feature_set")
          rank.var <- c("Log2FC", "NES")
        } else {
          prot.expr <- list(prot.df.noNA)
          features1 <- c("Gene")
          rank.var <- c("Log2FC")
        }
      } else {
        gsea2 <- NULL
        prot.expr <- list(prot.df.noNA)
        features1 <- feature.names[feature.names != "SUB_SITE"]
        rank.var <- c("Log2FC")
      }
      
      if (KSEA) {
        kin.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[["phospho_ksdb"]], 
          EA.type = EA.types[i])
        kin.net <- panSEA::netSEA(list(gsea1.inputs[["phospho_ksdb"]]),
                                  list(gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$result),
                                  "Feature_set", "NES",
                                  n.network.sets = 5)
        kin.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                              gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$result,
                                            "GSEA_volcano_plot.pdf" =
                                              gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$volcano.plot,
                                            "GSEA_network_graph.html" = 
                                              kin.net$interactive,
                                            "mtn_plots" = kin.mtn) 
      }
      
      if (SSEA) {
        sub.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[["phospho_sub"]], 
          EA.type = EA.types[i])
        sub.net <- panSEA::netSEA(list(gsea1.inputs[["phospho_sub"]]),
                                  list(gsea1[[gsea.name]]$all.results[["phospho_sub"]]$result),
                                  "Feature_set", "NES",
                                  n.network.sets = 5)
        sub.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                              gsea1[[gsea.name]]$all.results[["phospho_sub"]]$result,
                                            "GSEA_volcano_plot.pdf" =
                                              gsea1[[gsea.name]]$all.results[["phospho_sub"]]$volcano.plot,
                                            "GSEA_network_graph.html" = 
                                              sub.net$interactive,
                                            "mtn_plots" = sub.mtn) 
      }
      
      # run GSEA for each gmt in gmt1
      gsea1 <- list()
      global.gsea.files <- list()
      kin.gsea.files <- list()
      sub.gsea.files <- list()
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
        
        # store results for each input type
        global.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[["global"]], 
          EA.type = EA.types[i])
        global.net <- panSEA::netSEA(list(gsea1.inputs[["global"]]),
                                     list(gsea1[[gsea.name]]$all.results[["global"]]$result),
                                     n.network.sets = 5)
        global.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                 gsea1[[gsea.name]]$all.results[["global"]]$result,
                                               "GSEA_volcano_plot.pdf" =
                                                 gsea1[[gsea.name]]$all.results[["global"]]$volcano.plot,
                                               "GSEA_network_graph.html" = 
                                                 global.net$interactive,
                                               "mtn_plots" = global.mtn)
        
      }
      
      # run DMEA
      dmea.results <- panSEA::mDMEA(drug.sens, gmt.drug, expr, gsea1.inputs, 
                                    names(gsea1.inputs), 
                                    feature.names = features1,
                                    weight.values = rank.var)
      dmea.prot.results <- panSEA::mDMEA(drug.sens, gmt.drug, prot.expr, gsea1.inputs, 
                                         names(gsea1.inputs), 
                                         feature.names = features1,
                                         weight.values = rank.var)
      
      #### save results & upload to Synapse
      ### set file names
      ## global
      for (i in 1:length(types)) {
        temp.DEG.files <- list("Differential_expression_results.csv" = 
                                   deg[[i]])
        
        if (grepl("phospho", types[i], ignore.case = TRUE)) {
          if (KSEA) {
            kin.DMEA.mtn <- get_top_mtn_plots(dmea.results$all.results[["phospho_ksdb"]],
                                              sets = "Drug_set",
                                              EA.type = "DMEA")
            DMEA.kin.net <- panSEA::netSEA(list(dmea.results$all.results[["phospho_ksdb"]]$corr.result),
                                           list(dmea.results$all.results[["phospho_ksdb"]]$result),
                                           "Drug", "Pearson.est",
                                           n.network.sets = 5)
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
            
            prot.kin.DMEA.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[["phospho_ksdb"]],
                                                   sets = "Drug_set",
                                                   EA.type = "DMEA")
            prot.DMEA.kin.net <- panSEA::netSEA(list(prot.dmea.results$all.results[["phospho_ksdb"]]$corr.result),
                                                list(prot.dmea.results$all.results[["phospho_ksdb"]]$result),
                                                "Drug", "Pearson.est",
                                                n.network.sets = 5)
            prot.kin.DMEA.files <- list("DMEA_results.csv" =
                                          dmea.prot.results$all.results[["phospho_ksdb"]]$result,
                                        "DMEA_correlation_results.csv" = 
                                          dmea.prot.results$all.results[["phospho_ksdb"]]$corr.result,
                                        "DMEA_correlation_scatter_plots.pdf" = 
                                          dmea.prot.results$all.results[["phospho_ksdb"]]$corr.scatter.plots,
                                        "DMEA_volcano_plot.pdf" =
                                          dmea.prot.results$all.results[["phospho_ksdb"]]$volcano.plot,
                                        "DMEA_network_graph.html" = 
                                          prot.DMEA.kin.net$interactive,
                                        "mtn_plots" = prot.kin.DMEA.mtn)
          } else {
            kin.DMEA.files <- list()
            prot.kin.DMEA.files <- list()
          }
          
          if (SSEA) {
            sub.DMEA.mtn <- get_top_mtn_plots(dmea.results$all.results[["phospho_sub"]],
                                              sets = "Drug_set",
                                              EA.type = "DMEA")
            DMEA.sub.net <- panSEA::netSEA(list(dmea.results$all.results[["phospho_sub"]]$corr.result),
                                           list(dmea.results$all.results[["phospho_sub"]]$result),
                                           "Drug", "Pearson.est",
                                           n.network.sets = 5)
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
            prot.sub.DMEA.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[["phospho_sub"]],
                                                   sets = "Drug_set",
                                                   EA.type = "DMEA")
            prot.DMEA.sub.net <- panSEA::netSEA(list(dmea.prot.results$all.results[["phospho_sub"]]$corr.result),
                                                list(dmea.prot.results$all.results[["phospho_sub"]]$result),
                                                "Drug", "Pearson.est",
                                                n.network.sets = 5)
            prot.sub.DMEA.files <- list("DMEA_results.csv" =
                                          dmea.prot.results$all.results[["phospho_sub"]]$result,
                                        "DMEA_correlation_results.csv" = 
                                          dmea.prot.results$all.results[["phospho_sub"]]$corr.result,
                                        "DMEA_correlation_scatter_plots.pdf" = 
                                          dmea.prot.results$all.results[["phospho_sub"]]$corr.scatter.plots,
                                        "DMEA_volcano_plot.pdf" =
                                          dmea.prot.results$all.results[["phospho_sub"]]$volcano.plot,
                                        "DMEA_network_graph.html" = 
                                          prot.DMEA.sub.net$interactive,
                                        "mtn_plots" = prot.sub.DMEA.mtn)
          } else {
            sub.DMEA.files <- list()
            prot.sub.DMEA.files <- list()
          }
          
          phospho.kin.files <- list("KSEA_results.csv" =
                                      gsea2$all.results[["phospho_ksdb"]]$result,
                                    "KSEA_volcano_plot.pdf" =
                                      gsea2$all.results[["phospho_ksdb"]]$volcano.plot,
                                    "KSEA_network_graph.html" = 
                                      kin.net2$interactive,
                                    "mtn_plots" = kin.mtn2,
                                    "GSEA" = kin.gsea.files,
                                    "DMEA" = kin.DMEA.files,
                                    "DMEA_proteomics" = prot.kin.DMEA.files) 
          phospho.sub.files <- list("Substrate_enrichment_results.csv" =
                                      gsea2$all.results[["phospho_sub"]]$result,
                                    "Substrate_enrichment_volcano_plot.pdf" =
                                      gsea2$all.results[["phospho_sub"]]$volcano.plot,
                                    "Substrate_enrichment_network_graph.html" = 
                                      sub.net2$interactive,
                                    "mtn_plots" = sub.mtn2,
                                    "GSEA" = sub.gsea.files,
                                    "DMEA" = sub.DMEA.files,
                                    "DMEA_proteomics" = prot.sub.DMEA.files) 
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
                                            n.network.sets = 5)
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
          prot.DMEA.global.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[[i]],
                                                    sets = "Drug_set",
                                                    EA.type = "DMEA")
          prot.DMEA.global.net <- panSEA::netSEA(list(dmea.prot.results$all.results[[i]]$corr.result),
                                                 list(dmea.prot.results$all.results[[i]]$result),
                                                 "Drug", "Pearson.est",
                                                 n.network.sets = 5)
          prot.global.DMEA.files <- list("DMEA_results.csv" =
                                           dmea.prot.results$all.results[[i]]$result,
                                         "DMEA_correlation_results.csv" = 
                                           dmea.prot.results$all.results[[i]]$corr.result,
                                         "DMEA_correlation_scatter_plots.pdf" = 
                                           dmea.prot.results$all.results[[i]]$corr.scatter.plots,
                                         "DMEA_volcano_plot.pdf" =
                                           dmea.prot.results$all.results[[i]]$volcano.plot,
                                         "DMEA_network_graph.html" = 
                                           prot.DMEA.global.net$interactive,
                                         "mtn_plots" = prot.DMEA.global.mtn)
          global.files <- list('Differential_expression' = temp.DEG.files, 
                               'DMEA' = global.DMEA.files,
                               'DMEA_proteomics' = prot.global.DMEA.files,
                               'GSEA' = global.gsea.files)
          all.files[[types[i]]] <- global.files
        }
      }

      ## combo
      if (length(types) > 1) {
        combo.DMEA.files <- list("DMEA_results.csv" =
                                   dmea.results$compiled.results$results,
                                 "DMEA_mean_results.csv" =
                                   dmea.results$compiled.results$mean.results,
                                 "DMEA_correlation_matrix.pdf" =
                                   dmea.results$compiled.results$corr.matrix,
                                 "DMEA_dot_plot.pdf" =
                                   dmea.results$compiled.results$dot.plot)
        prot.combo.DMEA.files <- list("DMEA_results.csv" =
                                        dmea.prot.results$compiled.results$results,
                                      "DMEA_mean_results.csv" =
                                        dmea.prot.results$compiled.results$mean.results,
                                      "DMEA_correlation_matrix.pdf" =
                                        dmea.prot.results$compiled.results$corr.matrix,
                                      "DMEA_dot_plot.pdf" =
                                        dmea.prot.results$compiled.results$dot.plot)
        combo.files <- list('DMEA' = combo.DMEA.files,
                            'DMEA_proteomics' = prot.combo.DMEA.files)
        combo.name <- paste0(types, collapse = "_and_")
        all.files[[combo.name]] <- combo.files
      }
      
      # create folder for contrast
      setwd(temp.path)
      if (subfolder) {
        dir.create(contrast.name)
        setwd(contrast.name)
        contrastFolder <- 
          synapser::synStore(synapser::Folder(contrast.name,
                                              parent = synapse_id))
      } else {
        contrastFolder <- synapse_id
        contrast.name <- ""
      }
      
      save_to_synapse(all.files, contrastFolder)
      
      # compile DEGs if relevant
      if (compileDEGs) {
        # add feature type for each omics
        for (i in 1:length(deg)) {
          deg[[i]]$Feature_type <- colnames(deg[[i]])[1]
          colnames(deg[[i]])[1] <- "Feature"
        }
        
        # rbind and add contrast information
        temp.degs <- na.omit(data.table::rbindlist(deg, use.names = TRUE, fill = TRUE))
        temp.degs$Filter <- paste0(filterID, "_", filter)
        temp.degs$Contrast <- contrast.name
        all.degs <- rbind(all.degs, temp.degs)
      }
      
      # make space to process next contrast
      gsea1 <- NULL
      gsea2 <- NULL
      dmea.results <- NULL
      deg <- NULL
      gsea1.inputs <- NULL
    }
  }
  
  if (compileDEGs) {
    all.DEG.files <- list("Differential_expression_results.csv" = 
                            all.degs,
                          "Differential_expression_results_max_5_percent_FDR.csv" = 
                            all.degs[all.degs$adj.P.Val <= 0.05, ])
    save_to_synapse(all.DEG.files, synapse_id)
  }
}

run_contrasts <- function(contrast.type, 
                          contrasts = as.list(rep(c(TRUE, FALSE), length(contrast.type))),
                          id.type, meta.df, omics, 
                          types = c("global", "phospho"), 
                          feature.names = c("Gene", "SUB_SITE"),
                          gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                        "msigdb_Homo sapiens_H",
                                        "msigdb_Homo sapiens_C1"),
                          EA.types = c("KEGG", "Hallmark", "Positional"),
                          gmt.list2 = c("ksdb_human", "sub"),
                          expr = as.list(rep("adherent CCLE", length(types) + ifelse(grepl("phospho", types, ignore.case = TRUE),length(gmt.list2)-1,0))),
                          gmt.drug = "PRISM", drug.sens = "PRISM", 
                          base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                          temp.path = base.path, subfolder = TRUE, synapse_id = NULL,
                          compileDEGs = TRUE, filter = NA, filterID = NULL) {
  setwd(base.path)
  # load set annotations
  if (gmt.list1[1] == "chr8") {
    gmt1 <- get_chr8_gmt1()
  } else {
    gmt1 <- get_gmt1(gmt.list1)
  }
  
  if (grepl("phospho", types, ignore.case = TRUE)) {
    if (gmt.list2[1] == "chr8") {
      gmt2 <- get_chr8_gmt2()
    } else {
      gmt2 <- list()
      for (i in 1:length(gmt.list2)) {
        if (is.character(gmt.list2[i])) {
          if (grepl("ksdb", gmt.list2[i], ignore.case = TRUE)) {
            org <- stringr::str_split(gmt.list2[i], "_")[[1]][2]
            gmt2[[i]] <- get_ksdb(organism = org)
          } else if (gmt.list2[i] == "sub") {
            SUB_SITE <- omics[[2]]$SUB_SITE
            phospho.ref <- data.frame(SUB_SITE)
            phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
                                                          remove = FALSE)
            SUB_SITE <- NULL
            gmt2[[i]] <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE")
          }
        } else {
          gmt2[[i]] <- gmt.list2[i]
        }
      }
    } 
  }
  
  prot.df.noNA <- get_CCLE_prot()
  all.degs <- data.frame()
  
  ## for each contrast: 
  for (k in 1:length(contrasts)) {
    # identify samples for each side of contrast
    c1 <- contrasts[[k]][1]
    c2 <- contrasts[[k]][2]
    group.names <- c(c1, c2)
    contrast.name <- paste0(contrast.type, "_", contrasts[[k]][1], "_vs_", contrasts[[k]][2])
    group.samples <- list(as.vector(na.omit(meta.df[meta.df[,contrast.type] == c1, id.type]))[[id.type]],
                          as.vector(na.omit(meta.df[meta.df[,contrast.type] == c2, id.type]))[[id.type]])
    
    if (length(group.samples[[1]]) > 0 & length(group.samples[[2]]) > 0) {
      # run panSEA across omics types
      # CAUTION: this only works because the # of samples for each treatment type 
      # is equal; otherwise would have to run panSEA for each contrast separately 
      # and perhaps set the group.samples input parameter for panSEA
      
      # run global against KEGG sets, phospho against kinase sets from ksdb
      # run differential expression analysis
      deg <- panSEA::mDEG(omics, types, group.names, group.samples, 
                          feature.names)$all.results
      
      ## for phospho data:
      # run GSEA for each gmt in gmt2 and then check coverage of 2+ sets in gmt1
      gsea1.inputs <- deg
      KSEA <- FALSE
      SSEA <- FALSE
      if (grepl("phospho", names(deg), ignore.case = TRUE)) {
        n.phospho <- grep("phospho", names(deg), ignore.case = TRUE)
        gsea2.inputs <- list(deg[[n.phospho]], deg[[n.phospho]])
        gsea2 <- panSEA::mGSEA(gsea2.inputs, gmt2, types = c("phospho_ksdb", "phospho_sub"),
                               feature.names = rep(feature.names[n.phospho], 2)) 
        
        # create mountain plots and network graphs for KSEA, SSEA
        kin.mtn2 <- get_top_mtn_plots(gsea2$all.results[[1]],
                                      EA.type = "KSEA")
        kin.net2 <- panSEA::netSEA(list(gsea2.inputs[[1]]),
                                   list(gsea2$all.results[[1]]$result),
                                   "SUB_SITE",
                                   n.network.sets = 5)
        sub.mtn2 <- get_top_mtn_plots(gsea2$all.results[[2]],
                                      EA.type = "Substrate_enrichment")
        sub.net2 <- panSEA::netSEA(list(gsea2.inputs[[2]]),
                                   list(gsea2$all.results[[2]]$result),
                                   "SUB_SITE",
                                   n.network.sets = 5)
        
        if (nrow(gsea2$all.results[[1]]$result) > 100 & 
            nrow(gsea2$all.results[[2]]$result) > 100) {
          KSEA <- TRUE
          SSEA <- TRUE
          gsea1.inputs[["phospho_ksdb"]] <- gsea2$all.results[[1]]$result
          gsea1.inputs[["phospho_sub"]] <- gsea2$all.results[[2]]$result 
          prot.expr <- list(prot.df.noNA, prot.df.noNA, prot.df.noNA)
          features1 <- c("Gene", "Feature_set", "Feature_set")
          rank.var <- c("Log2FC", "NES", "NES")
        } else if (nrow(gsea2$all.results[[1]]$result) > 100) {
          KSEA <- TRUE
          gsea1.inputs[["phospho_ksdb"]] <- gsea2$all.results[[1]]$result
          prot.expr <- list(prot.df.noNA, prot.df.noNA)
          features1 <- c("Gene", "Feature_set")
          rank.var <- c("Log2FC", "NES")
        } else if (nrow(gsea2$all.results[[2]]$result) > 100) {
          SSEA <- TRUE
          gsea1.inputs[["phospho_sub"]] <- gsea2$all.results[[2]]$result 
          prot.expr <- list(prot.df.noNA, prot.df.noNA)
          features1 <- c("Gene", "Feature_set")
          rank.var <- c("Log2FC", "NES")
        } else {
          prot.expr <- list(prot.df.noNA)
          features1 <- c("Gene")
          rank.var <- c("Log2FC")
        }
      } else {
        gsea2 <- NULL
        prot.expr <- list(prot.df.noNA)
        features1 <- feature.names[feature.names != "SUB_SITE"]
        rank.var <- c("Log2FC")
      }
      
      if (KSEA) {
        kin.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[["phospho_ksdb"]], 
          EA.type = EA.types[i])
        kin.net <- panSEA::netSEA(list(gsea1.inputs[["phospho_ksdb"]]),
                                  list(gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$result),
                                  "Feature_set", "NES",
                                  n.network.sets = 5)
        kin.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                              gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$result,
                                            "GSEA_volcano_plot.pdf" =
                                              gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$volcano.plot,
                                            "GSEA_network_graph.html" = 
                                              kin.net$interactive,
                                            "mtn_plots" = kin.mtn) 
      }
      
      if (SSEA) {
        sub.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[["phospho_sub"]], 
          EA.type = EA.types[i])
        sub.net <- panSEA::netSEA(list(gsea1.inputs[["phospho_sub"]]),
                                  list(gsea1[[gsea.name]]$all.results[["phospho_sub"]]$result),
                                  "Feature_set", "NES",
                                  n.network.sets = 5)
        sub.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                              gsea1[[gsea.name]]$all.results[["phospho_sub"]]$result,
                                            "GSEA_volcano_plot.pdf" =
                                              gsea1[[gsea.name]]$all.results[["phospho_sub"]]$volcano.plot,
                                            "GSEA_network_graph.html" = 
                                              sub.net$interactive,
                                            "mtn_plots" = sub.mtn) 
      }
      
      # run GSEA for each gmt in gmt1
      gsea1 <- list()
      global.gsea.files <- list()
      kin.gsea.files <- list()
      sub.gsea.files <- list()
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
        
        # store results for each input type
        global.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[["global"]], 
          EA.type = EA.types[i])
        global.net <- panSEA::netSEA(list(gsea1.inputs[["global"]]),
                                     list(gsea1[[gsea.name]]$all.results[["global"]]$result),
                                     n.network.sets = 5)
        global.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                 gsea1[[gsea.name]]$all.results[["global"]]$result,
                                               "GSEA_volcano_plot.pdf" =
                                                 gsea1[[gsea.name]]$all.results[["global"]]$volcano.plot,
                                               "GSEA_network_graph.html" = 
                                                 global.net$interactive,
                                               "mtn_plots" = global.mtn)
        
      }
      
      # run DMEA
      dmea.results <- panSEA::mDMEA(drug.sens, gmt.drug, expr, gsea1.inputs, 
                                    names(gsea1.inputs), 
                                    feature.names = features1,
                                    weight.values = rank.var)
      dmea.prot.results <- panSEA::mDMEA(drug.sens, gmt.drug, prot.expr, gsea1.inputs, 
                                         names(gsea1.inputs), 
                                         feature.names = features1,
                                         weight.values = rank.var)
      
      #### save results & upload to Synapse
      ### set file names
      ## global
      for (i in 1:length(types)) {
        temp.DEG.files <- list("Differential_expression_results.csv" = 
                                 deg[[i]])
        
        if (grepl("phospho", types[i], ignore.case = TRUE)) {
          if (KSEA) {
            kin.DMEA.mtn <- get_top_mtn_plots(dmea.results$all.results[["phospho_ksdb"]],
                                              sets = "Drug_set",
                                              EA.type = "DMEA")
            DMEA.kin.net <- panSEA::netSEA(list(dmea.results$all.results[["phospho_ksdb"]]$corr.result),
                                           list(dmea.results$all.results[["phospho_ksdb"]]$result),
                                           "Drug", "Pearson.est",
                                           n.network.sets = 5)
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
            
            prot.kin.DMEA.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[["phospho_ksdb"]],
                                                   sets = "Drug_set",
                                                   EA.type = "DMEA")
            prot.DMEA.kin.net <- panSEA::netSEA(list(prot.dmea.results$all.results[["phospho_ksdb"]]$corr.result),
                                                list(prot.dmea.results$all.results[["phospho_ksdb"]]$result),
                                                "Drug", "Pearson.est",
                                                n.network.sets = 5)
            prot.kin.DMEA.files <- list("DMEA_results.csv" =
                                          dmea.prot.results$all.results[["phospho_ksdb"]]$result,
                                        "DMEA_correlation_results.csv" = 
                                          dmea.prot.results$all.results[["phospho_ksdb"]]$corr.result,
                                        "DMEA_correlation_scatter_plots.pdf" = 
                                          dmea.prot.results$all.results[["phospho_ksdb"]]$corr.scatter.plots,
                                        "DMEA_volcano_plot.pdf" =
                                          dmea.prot.results$all.results[["phospho_ksdb"]]$volcano.plot,
                                        "DMEA_network_graph.html" = 
                                          prot.DMEA.kin.net$interactive,
                                        "mtn_plots" = prot.kin.DMEA.mtn)
          } else {
            kin.DMEA.files <- list()
            prot.kin.DMEA.files <- list()
          }
          
          if (SSEA) {
            sub.DMEA.mtn <- get_top_mtn_plots(dmea.results$all.results[["phospho_sub"]],
                                              sets = "Drug_set",
                                              EA.type = "DMEA")
            DMEA.sub.net <- panSEA::netSEA(list(dmea.results$all.results[["phospho_sub"]]$corr.result),
                                           list(dmea.results$all.results[["phospho_sub"]]$result),
                                           "Drug", "Pearson.est",
                                           n.network.sets = 5)
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
            prot.sub.DMEA.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[["phospho_sub"]],
                                                   sets = "Drug_set",
                                                   EA.type = "DMEA")
            prot.DMEA.sub.net <- panSEA::netSEA(list(dmea.prot.results$all.results[["phospho_sub"]]$corr.result),
                                                list(dmea.prot.results$all.results[["phospho_sub"]]$result),
                                                "Drug", "Pearson.est",
                                                n.network.sets = 5)
            prot.sub.DMEA.files <- list("DMEA_results.csv" =
                                          dmea.prot.results$all.results[["phospho_sub"]]$result,
                                        "DMEA_correlation_results.csv" = 
                                          dmea.prot.results$all.results[["phospho_sub"]]$corr.result,
                                        "DMEA_correlation_scatter_plots.pdf" = 
                                          dmea.prot.results$all.results[["phospho_sub"]]$corr.scatter.plots,
                                        "DMEA_volcano_plot.pdf" =
                                          dmea.prot.results$all.results[["phospho_sub"]]$volcano.plot,
                                        "DMEA_network_graph.html" = 
                                          prot.DMEA.sub.net$interactive,
                                        "mtn_plots" = prot.sub.DMEA.mtn)
          } else {
            sub.DMEA.files <- list()
            prot.sub.DMEA.files <- list()
          }
          
          phospho.kin.files <- list("KSEA_results.csv" =
                                      gsea2$all.results[["phospho_ksdb"]]$result,
                                    "KSEA_volcano_plot.pdf" =
                                      gsea2$all.results[["phospho_ksdb"]]$volcano.plot,
                                    "KSEA_network_graph.html" = 
                                      kin.net2$interactive,
                                    "mtn_plots" = kin.mtn2,
                                    "GSEA" = kin.gsea.files,
                                    "DMEA" = kin.DMEA.files,
                                    "DMEA_proteomics" = prot.kin.DMEA.files) 
          phospho.sub.files <- list("Substrate_enrichment_results.csv" =
                                      gsea2$all.results[["phospho_sub"]]$result,
                                    "Substrate_enrichment_volcano_plot.pdf" =
                                      gsea2$all.results[["phospho_sub"]]$volcano.plot,
                                    "Substrate_enrichment_network_graph.html" = 
                                      sub.net2$interactive,
                                    "mtn_plots" = sub.mtn2,
                                    "GSEA" = sub.gsea.files,
                                    "DMEA" = sub.DMEA.files,
                                    "DMEA_proteomics" = prot.sub.DMEA.files) 
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
                                            n.network.sets = 5)
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
          prot.DMEA.global.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[[i]],
                                                    sets = "Drug_set",
                                                    EA.type = "DMEA")
          prot.DMEA.global.net <- panSEA::netSEA(list(dmea.prot.results$all.results[[i]]$corr.result),
                                                 list(dmea.prot.results$all.results[[i]]$result),
                                                 "Drug", "Pearson.est",
                                                 n.network.sets = 5)
          prot.global.DMEA.files <- list("DMEA_results.csv" =
                                           dmea.prot.results$all.results[[i]]$result,
                                         "DMEA_correlation_results.csv" = 
                                           dmea.prot.results$all.results[[i]]$corr.result,
                                         "DMEA_correlation_scatter_plots.pdf" = 
                                           dmea.prot.results$all.results[[i]]$corr.scatter.plots,
                                         "DMEA_volcano_plot.pdf" =
                                           dmea.prot.results$all.results[[i]]$volcano.plot,
                                         "DMEA_network_graph.html" = 
                                           prot.DMEA.global.net$interactive,
                                         "mtn_plots" = prot.DMEA.global.mtn)
          global.files <- list('Differential_expression' = temp.DEG.files, 
                               'DMEA' = global.DMEA.files,
                               'DMEA_proteomics' = prot.global.DMEA.files,
                               'GSEA' = global.gsea.files)
          all.files[[types[i]]] <- global.files
        }
      }
      
      ## combo
      if (length(types) > 1) {
        combo.DMEA.files <- list("DMEA_results.csv" =
                                   dmea.results$compiled.results$results,
                                 "DMEA_mean_results.csv" =
                                   dmea.results$compiled.results$mean.results,
                                 "DMEA_correlation_matrix.pdf" =
                                   dmea.results$compiled.results$corr.matrix,
                                 "DMEA_dot_plot.pdf" =
                                   dmea.results$compiled.results$dot.plot)
        prot.combo.DMEA.files <- list("DMEA_results.csv" =
                                        dmea.prot.results$compiled.results$results,
                                      "DMEA_mean_results.csv" =
                                        dmea.prot.results$compiled.results$mean.results,
                                      "DMEA_correlation_matrix.pdf" =
                                        dmea.prot.results$compiled.results$corr.matrix,
                                      "DMEA_dot_plot.pdf" =
                                        dmea.prot.results$compiled.results$dot.plot)
        combo.files <- list('DMEA' = combo.DMEA.files,
                            'DMEA_proteomics' = prot.combo.DMEA.files)
        combo.name <- paste0(types, collapse = "_and_")
        all.files[[combo.name]] <- combo.files
      }
      
      # create folder for contrast
      setwd(temp.path)
      if (subfolder) {
        dir.create(contrast.name)
        setwd(contrast.name)
        contrastFolder <- 
          synapser::synStore(synapser::Folder(contrast.name,
                                              parent = synapse_id))
      } else {
        contrastFolder <- synapse_id
        contrast.name <- ""
      }
      
      save_to_synapse(all.files, contrastFolder)
      
      # compile DEGs if relevant
      if (compileDEGs) {
        # add feature type for each omics
        for (i in 1:length(deg)) {
          deg[[i]]$Feature_type <- colnames(deg[[i]])[1]
          colnames(deg[[i]])[1] <- "Feature"
        }
        
        # rbind and add contrast information
        temp.degs <- na.omit(data.table::rbindlist(deg, use.names = TRUE, fill = TRUE))
        temp.degs$Filter <- paste0(filterID, "_", filter)
        temp.degs$Contrast <- contrast.name
        all.degs <- rbind(all.degs, temp.degs)
      }
      
      # make space to process next contrast
      gsea1 <- NULL
      gsea2 <- NULL
      dmea.results <- NULL
      deg <- NULL
      gsea1.inputs <- NULL
    }
  }
  
  if (compileDEGs) {
    all.DEG.files <- list("Differential_expression_results.csv" = 
                            all.degs,
                          "Differential_expression_results_max_5_percent_FDR.csv" = 
                            all.degs[all.degs$adj.P.Val <= 0.05, ])
    save_to_synapse(all.DEG.files, synapse_id)
  }
}

run_contrast_combos <- function(contrast.type, contrasts = as.list(rep(c(TRUE, FALSE), length(contrast.type))),
                          id.type, meta.df, 
                          omics, types = c("global", "phospho"),
                          feature.names = c("Gene", "SUB_SITE"),
                          gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                        "msigdb_Homo sapiens_H",
                                        "msigdb_Homo sapiens_C1"),
                          EA.types = c("KEGG", "Hallmark", "Positional"),
                          gmt.list2 = c("ksdb_human", "sub"),
                          expr = as.list(rep("adherent CCLE", 3)),
                          gmt.drug = "PRISM", drug.sens = "PRISM", 
                          base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                          temp.path, subfolder = TRUE,
                          synapse_id = NULL, compileDEGs = TRUE, 
                          filter = NA, filterID = NULL) {
  setwd(base.path)
  # load set annotations
  if (gmt.list1[1] == "chr8") {
    gmt1 <- get_chr8_gmt1()
  } else {
    gmt1 <- get_gmt1(gmt.list1)
  }
  
  if (grepl("phospho", types, ignore.case = TRUE)) {
    if (gmt.list2[1] == "chr8") {
      gmt2 <- get_chr8_gmt2()
    } else {
      gmt2 <- list()
      for (i in 1:length(gmt.list2)) {
        if (is.character(gmt.list2[i])) {
          if (grepl("ksdb", gmt.list2[i], ignore.case = TRUE)) {
            org <- stringr::str_split(gmt.list2[i], "_")[[1]][2]
            gmt2[[i]] <- get_ksdb(organism = org)
          } else if (gmt.list2[i] == "sub") {
            SUB_SITE <- omics[[2]]$SUB_SITE
            phospho.ref <- data.frame(SUB_SITE)
            phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
                                                          remove = FALSE)
            SUB_SITE <- NULL
            gmt2[[i]] <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE")
          }
        } else {
          gmt2[[i]] <- gmt.list2[i]
        }
      }
    } 
  }
  
  all.degs <- data.frame()
  # run contrasts with no filters
  setwd(temp.path)
  dir.create("no_filter")
  setwd("no_filter")
  nullPath <- file.path(temp.path, "no_filter")
  nullFolder <- 
    synapser::synStore(synapser::Folder("no_filter",
                                        parent = synapse_id))
  run_contrasts(contrast.types, id.type, meta.df, 
                omics, gmt.list1 = gmt1,
                EA.types, gmt.list2 = gmt2,
                expr,
                gmt.drug, drug.sens, 
                base.path, temp.path = nullPath, subfolder,
                synapse_id = nullFolder, compileDEGs)
  if (compileDEGs) {
    nullFiles <- as.list(synapser::synGetChildren(nullFolder, list("file"), sortBy = 'NAME'))
    nullFile <- synapser::synGet(nullFiles[[1]]$id)
    nullDEGs <- read.csv(nullFile$path)
    all.degs <- rbind(all.degs, nullDEGs)
  }
  
  for (m in 1:length(contrast.types)) {
    # run contrasts with TRUE filters
    setwd(temp.path)
    dir.create(file.path(paste0(contrast.types[m], "_TRUE")))
    setwd(file.path(paste0(contrast.types[m], "_TRUE")))
    truePath <- file.path(temp.path, paste0(contrast.types[m], "_TRUE"))
    trueFolder <- 
      synapser::synStore(synapser::Folder(file.path(paste0(contrast.types[m], "_TRUE")),
                                          parent = synapse_id))
    run_contrasts(contrast.types, id.type, meta.df, 
                                          omics, gmt.list1 = gmt1,
                                          EA.types,
                                          gmt.list2 = gmt2,
                                          expr,
                                          gmt.drug, drug.sens, 
                                          base.path, temp.path = truePath, subfolder,
                                          synapse_id = trueFolder, compileDEGs, 
                                          filter = "TRUE", filterID = contrast.types[m])
    if (compileDEGs) {
      trueDEGfiles <- as.list(synapser::synGetChildren(trueFolder, list("file"), sortBy = 'NAME'))
      trueDEGfile <- synapser::synGet(trueDEGfiles[[1]]$id)
      trueDEGs <- read.csv(trueDEGfile$path)
      all.degs <- rbind(all.degs, trueDEGs)
    }
    
    # run contrasts with FALSE filters
    setwd(temp.path)
    dir.create(file.path(paste0(contrast.types[m], "_FALSE")))
    setwd(file.path(paste0(contrast.types[m], "_FALSE")))
    falsePath <- file.path(temp.path, paste0(contrast.types[m], "_FALSE"))
    falseFolder <- 
      synapser::synStore(synapser::Folder(file.path(paste0(contrast.types[m], "_FALSE")),
                                          parent = synapse_id))
    run_contrasts(contrast.types, id.type, meta.df, 
                                          omics, gmt.list1 = gmt1,
                                          EA.types,
                                          gmt.list2 = gmt2,
                                          expr,
                                          gmt.drug, drug.sens, 
                                          base.path, temp.path = falsePath, subfolder,
                                          synapse_id = falseFolder, compileDEGs, 
                                          filter = "FALSE", filterID = contrast.types[m])
    if (compileDEGs) {
      falseDEGfiles <- as.list(synapser::synGetChildren(falseFolder, list("file"), sortBy = 'NAME'))
      falseDEGfile <- synapser::synGet(falseDEGfiles[[1]]$id)
      falseDEGs <- read.csv(falseDEGfile$path)
      all.degs <- rbind(all.degs, falseDEGs)
    }
  }
  
  if (compileDEGs) {
    setwd(temp.path)
    all.DEG.files <- list("Differential_expression_results.csv" = 
                            all.degs,
                          "Differential_expression_results_max_5_percent_FDR.csv" = 
                            all.degs[all.degs$adj.P.Val <= 0.05, ])
    save_to_synapse(all.DEG.files, synapse_id)
  }
}

run_TF_contrasts_global_phospho_human <- function(contrast.types, id.type, meta.df, 
                                               omics, gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                                                    "msigdb_Homo sapiens_H"),
                                               EA.types = c("KEGG", "Hallmark"),
                                               gmt.list2 = c("ksdb_human", "sub"),
                                               expr = as.list(rep("adherent CCLE", 3)),
                                               gmt.drug = "PRISM", drug.sens = "PRISM", 
                                               base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                                               temp.path = base.path, subfolder = TRUE,
                                               synapse_id = NULL, compileDEGs = TRUE, 
                                               filter = NA, filterID = NULL) {
  setwd(base.path)
  types = c("global", "phospho")
  feature.names = c("Gene", "SUB_SITE")
  if (!is.null(filterID)) {
    meta.df <- meta.df[meta.df[,filterID] == filter,]
    contrast.types <- contrast.types[contrast.types != filterID]
  }
  
  # load set annotations
  if (gmt.list1[1] == "chr8") {
    gmt1 <- get_chr8_gmt1()
  } else {
    gmt1 <- get_gmt1(gmt.list1)
  }
  
  if (gmt.list2[1] == "chr8") {
    gmt2 <- get_chr8_gmt2()
  } else {
    gmt2 <- list()
    for (i in 1:length(gmt.list2)) {
      if (is.character(gmt.list2[i])) {
        if (grepl("ksdb", gmt.list2[i], ignore.case = TRUE)) {
          org <- stringr::str_split(gmt.list2[i], "_")[[1]][2]
          gmt2[[i]] <- get_ksdb(organism = org)
        } else if (gmt.list2[i] == "sub") {
          SUB_SITE <- omics[[2]]$SUB_SITE
          phospho.ref <- data.frame(SUB_SITE)
          phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
                                                        remove = FALSE)
          SUB_SITE <- NULL
          gmt2[[i]] <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE")
        }
      } else {
        gmt2[[i]] <- gmt.list2[i]
      }
    }
  }
  
  prot.df.noNA <- get_CCLE_prot()
  all.degs <- data.frame()
  
  ## for each contrast: 
  for (k in 1:length(contrast.types)) {
    # identify samples for each side of contrast
    c1 <- "TRUE"
    c2 <- "FALSE"
    group.names <- c(c1, c2)
    contrast.name <- paste0(contrast.types[k], "_", c1, "_vs_", c2)
    group.samples <- list(as.vector(na.omit(meta.df[meta.df[,contrast.types[k]] == c1, id.type])),
                          as.vector(na.omit(meta.df[meta.df[,contrast.types[k]] == c2, id.type])))
    
    if (length(group.samples[[1]]) > 0 & length(group.samples[[2]]) > 0 & 
        length(unique(c(group.samples[[1]], group.samples[[2]]))) > 2) {
      # run panSEA across omics types
      # CAUTION: this only works because the # of samples for each treatment type 
      # is equal; otherwise would have to run panSEA for each contrast separately 
      # and perhaps set the group.samples input parameter for panSEA
      
      # run global against KEGG sets, phospho against kinase sets from ksdb
      # run differential expression analysis
      deg <- panSEA::mDEG(omics, types, group.names, group.samples, 
                          feature.names)$all.results
      
      ## for phospho data:
      # run GSEA for each gmt in gmt2 and then check coverage of 2+ sets in gmt1
      gsea2.inputs <- list(deg[[2]], deg[[2]])
      gsea2 <- panSEA::mGSEA(gsea2.inputs, gmt2, types = c("phospho_ksdb", "phospho_sub"),
                             feature.names = c("SUB_SITE", "SUB_SITE"))
      
      # create network graphs for KSEA, SSEA
      kin.net2 <- panSEA::netSEA(list(gsea2.inputs[[1]]),
                                 list(gsea2$all.results[[1]]$result),
                                 "SUB_SITE",
                                 n.network.sets = 5)
      sub.net2 <- panSEA::netSEA(list(gsea2.inputs[[2]]),
                                 list(gsea2$all.results[[2]]$result),
                                 "SUB_SITE",
                                 n.network.sets = 5)
      
      if (nrow(gsea2$all.results[[1]]$result) > 100 & 
          nrow(gsea2$all.results[[2]]$result) > 100) {
        KSEA <- TRUE
        SSEA <- TRUE
        gsea1.inputs <- list("global" = deg[[1]],
                             "phospho_ksdb" = gsea2$all.results[[1]]$result,
                             "phospho_sub" = gsea2$all.results[[2]]$result) 
        prot.expr <- list(prot.df.noNA, prot.df.noNA, prot.df.noNA)
        features1 <- c("Gene", "Feature_set", "Feature_set")
        rank.var <- c("Log2FC", "NES", "NES")
      } else if (nrow(gsea2$all.results[[1]]$result) > 100) {
        KSEA <- TRUE
        SSEA <- FALSE
        gsea1.inputs <- list("global" = deg[[1]],
                             "phospho_ksdb" = gsea2$all.results[[1]]$result)
        prot.expr <- list(prot.df.noNA, prot.df.noNA)
        features1 <- c("Gene", "Feature_set")
        rank.var <- c("Log2FC", "NES")
      } else if (nrow(gsea2$all.results[[2]]$result) > 100) {
        KSEA <- FALSE
        SSEA <- TRUE
        gsea1.inputs <- list("global" = deg[[1]],
                             "phospho_sub" = gsea2$all.results[[2]]$result) 
        prot.expr <- list(prot.df.noNA, prot.df.noNA)
        features1 <- c("Gene", "Feature_set")
        rank.var <- c("Log2FC", "NES")
      } else {
        KSEA <- FALSE
        SSEA <- FALSE
        gsea1.inputs <- list("global" = deg[[1]])
        prot.expr <- list(prot.df.noNA)
        features1 <- c("Gene")
        rank.var <- c("Log2FC")
      }
      
      # run GSEA for each gmt in gmt1
      gsea1 <- list()
      global.gsea.files <- list()
      kin.gsea.files <- list()
      sub.gsea.files <- list()
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
        
        # store results for each input type
        global.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[["global"]], 
          EA.type = EA.types[i])
        global.net <- panSEA::netSEA(list(gsea1.inputs[["global"]]),
                                     list(gsea1[[gsea.name]]$all.results[["global"]]$result),
                                     n.network.sets = 5)
        global.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                 gsea1[[gsea.name]]$all.results[["global"]]$result,
                                               "GSEA_volcano_plot.pdf" =
                                                 gsea1[[gsea.name]]$all.results[["global"]]$volcano.plot,
                                               "GSEA_network_graph.html" = 
                                                 global.net$interactive,
                                               "mtn_plots" = global.mtn)
        if (KSEA) {
          kin.mtn <- get_top_mtn_plots(
            gsea1[[gsea.name]]$all.results[["phospho_ksdb"]], 
            EA.type = EA.types[i])
          kin.net <- panSEA::netSEA(list(gsea1.inputs[["phospho_ksdb"]]),
                                    list(gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$result),
                                    "Feature_set", "NES",
                                    n.network.sets = 5)
          kin.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$result,
                                              "GSEA_volcano_plot.pdf" =
                                                gsea1[[gsea.name]]$all.results[["phospho_ksdb"]]$volcano.plot,
                                              "GSEA_network_graph.html" = 
                                                kin.net$interactive,
                                              "mtn_plots" = kin.mtn) 
        }
        
        if (SSEA) {
          sub.mtn <- get_top_mtn_plots(
            gsea1[[gsea.name]]$all.results[["phospho_sub"]], 
            EA.type = EA.types[i])
          sub.net <- panSEA::netSEA(list(gsea1.inputs[["phospho_sub"]]),
                                    list(gsea1[[gsea.name]]$all.results[["phospho_sub"]]$result),
                                    "Feature_set", "NES",
                                    n.network.sets = 5)
          sub.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                gsea1[[gsea.name]]$all.results[["phospho_sub"]]$result,
                                              "GSEA_volcano_plot.pdf" =
                                                gsea1[[gsea.name]]$all.results[["phospho_sub"]]$volcano.plot,
                                              "GSEA_network_graph.html" = 
                                                sub.net$interactive,
                                              "mtn_plots" = sub.mtn) 
        }
      }
      
      # run DMEA
      dmea.results <- panSEA::mDMEA(drug.sens, gmt.drug, expr, gsea1.inputs, 
                                    names(gsea1.inputs), 
                                    feature.names = features1,
                                    weight.values = rank.var)
      dmea.prot.results <- panSEA::mDMEA(drug.sens, gmt.drug, prot.expr, gsea1.inputs, 
                                         names(gsea1.inputs), 
                                         feature.names = features1,
                                         weight.values = rank.var)
      
      #### save results & upload to Synapse
      ### set file names
      ## global
      global.DEG.files <- list("Differential_expression_results.csv" = 
                                 deg[["global"]])
      DMEA.global.mtn <- get_top_mtn_plots(dmea.results$all.results[["global"]],
                                           sets = "Drug_set",
                                           EA.type = "DMEA")
      DMEA.global.net <- panSEA::netSEA(list(dmea.results$all.results[["global"]]$corr.result),
                                        list(dmea.results$all.results[["global"]]$result),
                                        "Drug", "Pearson.est",
                                        n.network.sets = 5)
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
      prot.DMEA.global.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[[1]],
                                                sets = "Drug_set",
                                                EA.type = "DMEA")
      prot.DMEA.global.net <- panSEA::netSEA(list(dmea.prot.results$all.results[[1]]$corr.result),
                                             list(dmea.prot.results$all.results[[1]]$result),
                                             "Drug", "Pearson.est",
                                             n.network.sets = 5)
      prot.global.DMEA.files <- list("DMEA_results.csv" =
                                       dmea.prot.results$all.results[[1]]$result,
                                     "DMEA_correlation_results.csv" = 
                                       dmea.prot.results$all.results[[1]]$corr.result,
                                     "DMEA_correlation_scatter_plots.pdf" = 
                                       dmea.prot.results$all.results[[1]]$corr.scatter.plots,
                                     "DMEA_volcano_plot.pdf" =
                                       dmea.prot.results$all.results[[1]]$volcano.plot,
                                     "DMEA_network_graph.html" = 
                                       prot.DMEA.global.net$interactive,
                                     "mtn_plots" = prot.DMEA.global.mtn)
      global.files <- list('Differential_expression' = global.DEG.files, 
                           'DMEA' = global.DMEA.files,
                           'DMEA_proteomics' = prot.global.DMEA.files,
                           'GSEA' = global.gsea.files)
      
      ## phospho
      phospho.DEG.files <- list("Differential_expression_results.csv" = 
                                  deg[[2]])
      
      kin.mtn2 <- get_top_mtn_plots(gsea2$all.results[[1]],
                                    EA.type = "KSEA")
      if (KSEA) {
        kin.DMEA.mtn <- get_top_mtn_plots(dmea.results$all.results[[2]],
                                          sets = "Drug_set",
                                          EA.type = "DMEA")
        DMEA.kin.net <- panSEA::netSEA(list(dmea.results$all.results[["phospho_ksdb"]]$corr.result),
                                       list(dmea.results$all.results[["phospho_ksdb"]]$result),
                                       "Drug", "Pearson.est",
                                       n.network.sets = 5)
        kin.DMEA.files <- list("DMEA_results.csv" =
                                 dmea.results$all.results[[2]]$result,
                               "DMEA_correlation_results.csv" = 
                                 dmea.results$all.results[[2]]$corr.result,
                               "DMEA_correlation_scatter_plots.pdf" = 
                                 dmea.results$all.results[[2]]$corr.scatter.plots,
                               "DMEA_volcano_plot.pdf" =
                                 dmea.results$all.results[[2]]$volcano.plot,
                               "DMEA_network_graph.html" = 
                                 DMEA.kin.net$interactive,
                               "mtn_plots" = kin.DMEA.mtn)
        
        prot.kin.DMEA.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[[2]],
                                               sets = "Drug_set",
                                               EA.type = "DMEA")
        prot.DMEA.kin.net <- panSEA::netSEA(list(prot.dmea.results$all.results[["phospho_ksdb"]]$corr.result),
                                            list(prot.dmea.results$all.results[["phospho_ksdb"]]$result),
                                            "Drug", "Pearson.est",
                                            n.network.sets = 5)
        prot.kin.DMEA.files <- list("DMEA_results.csv" =
                                      dmea.prot.results$all.results[[2]]$result,
                                    "DMEA_correlation_results.csv" = 
                                      dmea.prot.results$all.results[[2]]$corr.result,
                                    "DMEA_correlation_scatter_plots.pdf" = 
                                      dmea.prot.results$all.results[[2]]$corr.scatter.plots,
                                    "DMEA_volcano_plot.pdf" =
                                      dmea.prot.results$all.results[[2]]$volcano.plot,
                                    "DMEA_network_graph.html" = 
                                      prot.DMEA.kin.net$interactive,
                                    "mtn_plots" = prot.kin.DMEA.mtn)
      } else {
        kin.DMEA.files <- list()
        prot.kin.DMEA.files <- list()
      }
      phospho.kin.files <- list("KSEA_results.csv" =
                                  gsea2$all.results[[1]]$result,
                                "KSEA_volcano_plot.pdf" =
                                  gsea2$all.results[[1]]$volcano.plot,
                                "KSEA_network_graph.html" = 
                                  kin.net2$interactive,
                                "mtn_plots" = kin.mtn2,
                                "GSEA" = kin.gsea.files,
                                "DMEA" = kin.DMEA.files,
                                "DMEA_proteomics" = prot.kin.DMEA.files) 
      
      
      sub.mtn2 <- get_top_mtn_plots(gsea2$all.results[[2]],
                                    EA.type = "Substrate_enrichment")
      if (SSEA) {
        sub.DMEA.mtn <- get_top_mtn_plots(dmea.results$all.results[["phospho_sub"]],
                                          sets = "Drug_set",
                                          EA.type = "DMEA")
        DMEA.sub.net <- panSEA::netSEA(list(dmea.results$all.results[["phospho_sub"]]$corr.result),
                                       list(dmea.results$all.results[["phospho_sub"]]$result),
                                       "Drug", "Pearson.est",
                                       n.network.sets = 5)
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
        prot.sub.DMEA.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[["phospho_sub"]],
                                               sets = "Drug_set",
                                               EA.type = "DMEA")
        prot.DMEA.sub.net <- panSEA::netSEA(list(dmea.prot.results$all.results[["phospho_sub"]]$corr.result),
                                            list(dmea.prot.results$all.results[["phospho_sub"]]$result),
                                            "Drug", "Pearson.est",
                                            n.network.sets = 5)
        prot.sub.DMEA.files <- list("DMEA_results.csv" =
                                      dmea.prot.results$all.results[["phospho_sub"]]$result,
                                    "DMEA_correlation_results.csv" = 
                                      dmea.prot.results$all.results[["phospho_sub"]]$corr.result,
                                    "DMEA_correlation_scatter_plots.pdf" = 
                                      dmea.prot.results$all.results[["phospho_sub"]]$corr.scatter.plots,
                                    "DMEA_volcano_plot.pdf" =
                                      dmea.prot.results$all.results[["phospho_sub"]]$volcano.plot,
                                    "DMEA_network_graph.html" = 
                                      prot.DMEA.sub.net$interactive,
                                    "mtn_plots" = prot.sub.DMEA.mtn)
      } else {
        sub.DMEA.files <- list()
        prot.sub.DMEA.files <- list()
      }
      phospho.sub.files <- list("Substrate_enrichment_results.csv" =
                                  gsea2$all.results[[2]]$result,
                                "Substrate_enrichment_volcano_plot.pdf" =
                                  gsea2$all.results[[2]]$volcano.plot,
                                "Substrate_enrichment_network_graph.html" = 
                                  sub.net2$interactive,
                                "mtn_plots" = sub.mtn2,
                                "GSEA" = sub.gsea.files,
                                "DMEA" = sub.DMEA.files,
                                "DMEA_proteomics" = prot.sub.DMEA.files) 
      
      phospho.files <- list('Differential_expression' = phospho.DEG.files, 
                            'KSEA' = phospho.kin.files,
                            'Substrate_enrichment' = phospho.sub.files)
      
      ## combo
      combo.DMEA.files <- list("DMEA_results.csv" =
                                 dmea.results$compiled.results$results,
                               "DMEA_mean_results.csv" =
                                 dmea.results$compiled.results$mean.results,
                               "DMEA_correlation_matrix.pdf" =
                                 dmea.results$compiled.results$corr.matrix,
                               "DMEA_dot_plot.pdf" =
                                 dmea.results$compiled.results$dot.plot)
      prot.combo.DMEA.files <- list("DMEA_results.csv" =
                                      dmea.prot.results$compiled.results$results,
                                    "DMEA_mean_results.csv" =
                                      dmea.prot.results$compiled.results$mean.results,
                                    "DMEA_correlation_matrix.pdf" =
                                      dmea.prot.results$compiled.results$corr.matrix,
                                    "DMEA_dot_plot.pdf" =
                                      dmea.prot.results$compiled.results$dot.plot)
      combo.files <- list('DMEA' = combo.DMEA.files,
                          'DMEA_proteomics' = prot.combo.DMEA.files)
      
      all.files <- list('global_and_phospho' = combo.files,
                        'global' = global.files,
                        'phospho' = phospho.files)
      
      # create folder for contrast
      setwd(temp.path)
      if (subfolder) {
        dir.create(contrast.name)
        setwd(contrast.name)
        contrastFolder <- 
          synapser::synStore(synapser::Folder(contrast.name,
                                              parent = synapse_id))
      } else {
        contrastFolder <- synapse_id
        contrast.name <- ""
      }
      
      save_to_synapse(all.files, contrastFolder)
      
      # compile DEGs if relevant
      if (compileDEGs) {
        # get global and phospho degs
        global.degs <- deg[["global"]]
        phospho.degs <- deg[[2]]
        
        # add feature_type
        global.degs$Feature_type <- colnames(global.degs)[1]
        phospho.degs$Feature_type <- colnames(phospho.degs)[1]
        colnames(global.degs)[1] <- "Feature"
        colnames(phospho.degs)[1] <- "Feature"
        
        # rbind and add contrast information
        temp.degs <- na.omit(rbind(global.degs, phospho.degs))
        temp.degs$Filter <- paste0(filterID, "_", filter)
        temp.degs$Contrast <- contrast.name
        all.degs <- rbind(all.degs, temp.degs)
      }
      
      # make space to process next contrast
      gsea1 <- NULL
      gsea2 <- NULL
      dmea.results <- NULL
      deg <- NULL
      gsea1.inputs <- NULL
    } 
  }
  
  if (compileDEGs) {
    all.DEG.files <- list("Differential_expression_results.csv" = 
                                all.degs,
                          "Differential_expression_results_max_5_percent_FDR.csv" = 
                            all.degs[all.degs$adj.P.Val <= 0.05, ])
    save_to_synapse(all.DEG.files, synapse_id)
  }
}

run_TF_contrast_combos_global_phospho_human <- function(contrast.types, id.type, meta.df, 
                                                  omics, gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                                                       "msigdb_Homo sapiens_H"),
                                                  EA.types = c("KEGG", "Hallmark"),
                                                  gmt.list2 = c("ksdb_human", "sub"),
                                                  expr = as.list(rep("adherent CCLE", 3)),
                                                  gmt.drug = "PRISM", drug.sens = "PRISM", 
                                                  base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                                                  temp.path = base.path, subfolder = TRUE,
                                                  synapse_id = NULL, compileDEGs = TRUE) {
  setwd(base.path)
  types = c("global", "phospho")
  feature.names = c("Gene", "SUB_SITE")
  
  # load set annotations
  if (gmt.list1[1] == "chr8") {
    gmt1 <- get_chr8_gmt1()
  } else {
    gmt1 <- get_gmt1(gmt.list1)
  }
  
  if (gmt.list2[1] == "chr8") {
    gmt2 <- get_chr8_gmt2()
  } else {
    gmt2 <- list()
    for (i in 1:length(gmt.list2)) {
      if (is.character(gmt.list2[i])) {
        if (grepl("ksdb", gmt.list2[i], ignore.case = TRUE)) {
          org <- stringr::str_split(gmt.list2[i], "_")[[1]][2]
          gmt2[[i]] <- get_ksdb(organism = org)
        } else if (gmt.list2[i] == "sub") {
          SUB_SITE <- omics[[2]]$SUB_SITE
          phospho.ref <- data.frame(SUB_SITE)
          phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
                                                        remove = FALSE)
          SUB_SITE <- NULL
          gmt2[[i]] <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE")
        } else {
          gmt2[[i]] <- gmt.list2[i]
        }
      }
    }
  }
  
  all.degs <- data.frame()
  # run contrasts with no filters
  setwd(temp.path)
  dir.create("no_filter")
  setwd("no_filter")
  nullPath <- file.path(temp.path, "no_filter")
  nullFolder <- 
    synapser::synStore(synapser::Folder("no_filter",
                                        parent = synapse_id))
  run_TF_contrasts_global_phospho_human(contrast.types, id.type, meta.df, 
                                omics, gmt.list1 = gmt1,
                                EA.types, gmt.list2 = gmt2,
                                expr,
                                gmt.drug, drug.sens, 
                                base.path, temp.path = nullPath, subfolder,
                                synapse_id = nullFolder, compileDEGs)
  if (compileDEGs) {
    nullFiles <- as.list(synapser::synGetChildren(nullFolder, list("file"), sortBy = 'NAME'))
    nullFile <- synapser::synGet(nullFiles[[1]]$id)
    nullDEGs <- read.csv(nullFile$path)
    all.degs <- rbind(all.degs, nullDEGs)
  }
  
  for (m in 1:length(contrast.types)) {
    # run contrasts with TRUE filters
    setwd(temp.path)
    dir.create(file.path(paste0(contrast.types[m], "_TRUE")))
    setwd(file.path(paste0(contrast.types[m], "_TRUE")))
    truePath <- file.path(temp.path, paste0(contrast.types[m], "_TRUE"))
    trueFolder <- 
      synapser::synStore(synapser::Folder(file.path(paste0(contrast.types[m], "_TRUE")),
                                          parent = synapse_id))
    run_TF_contrasts_global_phospho_human(contrast.types, id.type, meta.df, 
                                          omics, gmt.list1 = gmt1,
                                          EA.types,
                                          gmt.list2 = gmt2,
                                          expr,
                                          gmt.drug, drug.sens, 
                                          base.path, temp.path = truePath, subfolder,
                                          synapse_id = trueFolder, compileDEGs, 
                                          filter = "TRUE", filterID = contrast.types[m])
    if (compileDEGs) {
      trueDEGfiles <- as.list(synapser::synGetChildren(trueFolder, list("file"), sortBy = 'NAME'))
      trueDEGfile <- synapser::synGet(trueDEGfiles[[1]]$id)
      trueDEGs <- read.csv(trueDEGfile$path)
      all.degs <- rbind(all.degs, trueDEGs)
    }
    
    # run contrasts with FALSE filters
    setwd(temp.path)
    dir.create(file.path(paste0(contrast.types[m], "_FALSE")))
    setwd(file.path(paste0(contrast.types[m], "_FALSE")))
    falsePath <- file.path(temp.path, paste0(contrast.types[m], "_FALSE"))
    falseFolder <- 
      synapser::synStore(synapser::Folder(file.path(paste0(contrast.types[m], "_FALSE")),
                                          parent = synapse_id))
    run_TF_contrasts_global_phospho_human(contrast.types, id.type, meta.df, 
                                          omics, gmt.list1 = gmt1,
                                          EA.types,
                                          gmt.list2 = gmt2,
                                          expr,
                                          gmt.drug, drug.sens, 
                                          base.path, temp.path = falsePath, subfolder,
                                          synapse_id = falseFolder, compileDEGs, 
                                          filter = "FALSE", filterID = contrast.types[m])
    if (compileDEGs) {
      falseDEGfiles <- as.list(synapser::synGetChildren(falseFolder, list("file"), sortBy = 'NAME'))
      falseDEGfile <- synapser::synGet(falseDEGfiles[[1]]$id)
      falseDEGs <- read.csv(falseDEGfile$path)
      all.degs <- rbind(all.degs, falseDEGs)
    }
  }
  
  if (compileDEGs) {
    setwd(temp.path)
    all.DEG.files <- list("Differential_expression_results.csv" = 
                            all.degs,
                          "Differential_expression_results_max_5_percent_FDR.csv" = 
                            all.degs[all.degs$adj.P.Val <= 0.05, ])
    save_to_synapse(all.DEG.files, synapse_id)
  }
}

run_contrasts_global_human <- function(contrasts, contrast.type, id.type, meta.df, 
                                               omics, gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                                                    "msigdb_Homo sapiens_H",
                                                                    "msigdb_Homo sapiens_C1"),
                                               EA.types = c("KEGG", "Hallmark", "Positional"),
                                               expr = as.list(rep("adherent CCLE", 1)),
                                               gmt.drug = "PRISM", drug.sens = "PRISM", 
                                               base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                                               temp.path = base.path, subfolder = TRUE,
                                               synapse_id = NULL) {
  setwd(base.path)
  types = c("global")
  feature.names = c("Gene")
  
  # load set annotations
  if (gmt.list1[1] == "chr8") {
    gmt1 <- get_chr8_gmt1()
  } else {
    gmt1 <- get_gmt1(gmt.list1)
  }
  
  prot.df.noNA <- get_CCLE_prot()
  
  ## for each contrast: 
  for (k in 1:length(contrasts)) {
    # identify samples for each side of contrast
    c1 <- contrasts[[k]][1]
    c2 <- contrasts[[k]][2]
    group.names <- c(c1, c2)
    contrast.name <- paste0(contrast.type, "_", contrasts[[k]][1], "_vs_", contrasts[[k]][2])
    group.samples <- list(as.vector(na.omit(meta.df[meta.df[,contrast.type] == c1, id.type]))[[id.type]],
                          as.vector(na.omit(meta.df[meta.df[,contrast.type] == c2, id.type]))[[id.type]])
    
    if (length(group.samples[[1]]) > 0 & length(group.samples[[2]]) > 0) {
      # run panSEA across omics types
      # CAUTION: this only works because the # of samples for each treatment type 
      # is equal; otherwise would have to run panSEA for each contrast separately 
      # and perhaps set the group.samples input parameter for panSEA
      
      # run global against KEGG sets, phospho against kinase sets from ksdb
      # run differential expression analysis
      deg <- panSEA::mDEG(omics, types, group.names, group.samples, 
                          feature.names)$all.results
      
      # run GSEA for each gmt in gmt1
      gsea1 <- list()
      global.gsea.files <- list()
      for (i in 1:length(gmt1)) {
        gsea.name <- paste0("GSEA_", EA.types[i])
        
        # prep set annotations
        temp.gmt1 <- list()
        for (j in 1:length(deg)) {
          temp.gmt1[[j]] <- gmt1[[i]]
        }
        
        # run GSEA
        gsea1[[gsea.name]] <- panSEA::mGSEA(deg, temp.gmt1, 
                                            types = "global")
        
        # store results for each input type
        global.mtn <- get_top_mtn_plots(
          gsea1[[gsea.name]]$all.results[["global"]], 
          EA.type = EA.types[i])
        global.net <- panSEA::netSEA(list(deg[["global"]]),
                                     list(gsea1[[gsea.name]]$all.results[["global"]]$result),
                                     n.network.sets = 5)
        global.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                 gsea1[[gsea.name]]$all.results[["global"]]$result,
                                               "GSEA_volcano_plot.pdf" =
                                                 gsea1[[gsea.name]]$all.results[["global"]]$volcano.plot,
                                               "GSEA_network_graph.html" = 
                                                 global.net$interactive,
                                               "mtn_plots" = global.mtn)
      }
      
      # run DMEA
      dmea.results <- panSEA::mDMEA(drug.sens, gmt.drug, expr, deg, 
                                    names(deg))
      prot.expr <- list(prot.df.noNA)
      dmea.prot.results <- panSEA::mDMEA(drug.sens, gmt.drug, prot.expr, deg, 
                                         names(deg))
      
      #### save results & upload to Synapse
      ### set file names
      ## global
      global.DEG.files <- list("Differential_expression_results.csv" = 
                                 deg[["global"]])
      DMEA.global.mtn <- get_top_mtn_plots(dmea.results$all.results[["global"]],
                                           sets = "Drug_set",
                                           EA.type = "DMEA")
      DMEA.global.net <- panSEA::netSEA(list(dmea.results$all.results[["global"]]$corr.result),
                                        list(dmea.results$all.results[["global"]]$result),
                                        "Drug", "Pearson.est",
                                        n.network.sets = 5)
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
      prot.DMEA.global.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[[1]],
                                                sets = "Drug_set",
                                                EA.type = "DMEA")
      prot.DMEA.global.net <- panSEA::netSEA(list(dmea.prot.results$all.results[[1]]$corr.result),
                                             list(dmea.prot.results$all.results[[1]]$result),
                                             "Drug", "Pearson.est",
                                             n.network.sets = 5)
      prot.global.DMEA.files <- list("DMEA_results.csv" =
                                       dmea.prot.results$all.results[[1]]$result,
                                     "DMEA_correlation_results.csv" = 
                                       dmea.prot.results$all.results[[1]]$corr.result,
                                     "DMEA_correlation_scatter_plots.pdf" = 
                                       dmea.prot.results$all.results[[1]]$corr.scatter.plots,
                                     "DMEA_volcano_plot.pdf" =
                                       dmea.prot.results$all.results[[1]]$volcano.plot,
                                     "DMEA_network_graph.html" = 
                                       prot.DMEA.global.net$interactive,
                                     "mtn_plots" = prot.DMEA.global.mtn)
      global.files <- list('Differential_expression' = global.DEG.files, 
                           'DMEA' = global.DMEA.files,
                           'DMEA_proteomics' = prot.global.DMEA.files,
                           'GSEA' = global.gsea.files)
      all.files <- list('global' = global.files)
      
      # create folder for contrast
      setwd(temp.path)
      if (subfolder) {
        dir.create(contrast.name)
        setwd(contrast.name)
        contrastFolder <- 
          synapser::synStore(synapser::Folder(contrast.name,
                                              parent = synapse_id))
      } else {
        contrastFolder <- synapse_id
        contrast.name <- ""
      }
      
      save_to_synapse(all.files, contrastFolder)
      
      # make space to process next contrast
      gsea1 <- NULL
      dmea.results <- NULL
      deg <- NULL
    }
  }
}

run_TF_contrasts_global_human <- function(contrast.types, id.type, meta.df, 
                                                  omics, gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                                                       "msigdb_Homo sapiens_H"),
                                                  EA.types = c("KEGG", "Hallmark"),
                                                  expr = as.list(rep("adherent CCLE", 1)),
                                                  gmt.drug = "PRISM", drug.sens = "PRISM", 
                                                  base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                                                  temp.path = base.path, subfolder = TRUE,
                                                  synapse_id = NULL, compileDEGs = TRUE, 
                                                  filter = NA, filterID = NULL) {
  setwd(base.path)
  types = c("global")
  feature.names = c("Gene")
  if (!is.null(filterID)) {
    meta.df <- meta.df[meta.df[,filterID] == filter,]
    contrast.types <- contrast.types[contrast.types != filterID]
  }
  
  all.degs <- data.frame()
  if (nrow(meta.df) > 2) {
    # load set annotations
    if (gmt.list1[1] == "chr8") {
      gmt1 <- get_chr8_gmt1()
    } else {
      gmt1 <- get_gmt1(gmt.list1)
    }
    
    prot.df.noNA <- get_CCLE_prot()
    
    ## for each contrast: 
    for (k in 1:length(contrast.types)) {
      # identify samples for each side of contrast
      c1 <- "TRUE"
      c2 <- "FALSE"
      group.names <- c(c1, c2)
      contrast.name <- paste0(contrast.types[k], "_", c1, "_vs_", c2)
      group.samples <- list(as.vector(na.omit(meta.df[meta.df[,contrast.types[k]] == c1, id.type])),
                            as.vector(na.omit(meta.df[meta.df[,contrast.types[k]] == c2, id.type])))
      
      if (length(group.samples[[1]]) > 0 & length(group.samples[[2]]) > 0 & 
          length(unique(c(group.samples[[1]], group.samples[[2]]))) > 2) {
        # run panSEA across omics types
        # CAUTION: this only works because the # of samples for each treatment type 
        # is equal; otherwise would have to run panSEA for each contrast separately 
        # and perhaps set the group.samples input parameter for panSEA
        
        # run global against KEGG sets, phospho against kinase sets from ksdb
        # run differential expression analysis
        deg <- panSEA::mDEG(omics, types, group.names, group.samples, 
                            feature.names)$all.results
        gsea1.inputs <- list("global" = deg[[1]])
        prot.expr <- list(prot.df.noNA)
        features1 <- c("Gene")
        rank.var <- c("Log2FC")
        
        # run GSEA for each gmt in gmt1
        gsea1 <- list()
        global.gsea.files <- list()
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
          
          # store results for each input type
          global.mtn <- get_top_mtn_plots(
            gsea1[[gsea.name]]$all.results[["global"]], 
            EA.type = EA.types[i])
          global.net <- panSEA::netSEA(list(gsea1.inputs[["global"]]),
                                       list(gsea1[[gsea.name]]$all.results[["global"]]$result),
                                       n.network.sets = 5)
          global.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                   gsea1[[gsea.name]]$all.results[["global"]]$result,
                                                 "GSEA_volcano_plot.pdf" =
                                                   gsea1[[gsea.name]]$all.results[["global"]]$volcano.plot,
                                                 "GSEA_network_graph.html" = 
                                                   global.net$interactive,
                                                 "mtn_plots" = global.mtn)
        }
        
        # run DMEA
        dmea.results <- panSEA::mDMEA(drug.sens, gmt.drug, expr, gsea1.inputs, 
                                      names(gsea1.inputs), 
                                      feature.names = features1,
                                      weight.values = rank.var)
        dmea.prot.results <- panSEA::mDMEA(drug.sens, gmt.drug, prot.expr, gsea1.inputs, 
                                           names(gsea1.inputs), 
                                           feature.names = features1,
                                           weight.values = rank.var)
        
        #### save results & upload to Synapse
        ### set file names
        ## global
        global.DEG.files <- list("Differential_expression_results.csv" = 
                                   deg[["global"]])
        DMEA.global.mtn <- get_top_mtn_plots(dmea.results$all.results[["global"]],
                                             sets = "Drug_set",
                                             EA.type = "DMEA")
        DMEA.global.net <- panSEA::netSEA(list(dmea.results$all.results[["global"]]$corr.result),
                                          list(dmea.results$all.results[["global"]]$result),
                                          "Drug", "Pearson.est",
                                          n.network.sets = 5)
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
        prot.DMEA.global.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[[1]],
                                                  sets = "Drug_set",
                                                  EA.type = "DMEA")
        prot.DMEA.global.net <- panSEA::netSEA(list(dmea.prot.results$all.results[[1]]$corr.result),
                                               list(dmea.prot.results$all.results[[1]]$result),
                                               "Drug", "Pearson.est",
                                               n.network.sets = 5)
        prot.global.DMEA.files <- list("DMEA_results.csv" =
                                         dmea.prot.results$all.results[[1]]$result,
                                       "DMEA_correlation_results.csv" = 
                                         dmea.prot.results$all.results[[1]]$corr.result,
                                       "DMEA_correlation_scatter_plots.pdf" = 
                                         dmea.prot.results$all.results[[1]]$corr.scatter.plots,
                                       "DMEA_volcano_plot.pdf" =
                                         dmea.prot.results$all.results[[1]]$volcano.plot,
                                       "DMEA_network_graph.html" = 
                                         prot.DMEA.global.net$interactive,
                                       "mtn_plots" = prot.DMEA.global.mtn)
        global.files <- list('Differential_expression' = global.DEG.files, 
                             'DMEA' = global.DMEA.files,
                             'DMEA_proteomics' = prot.global.DMEA.files,
                             'GSEA' = global.gsea.files)
        
        all.files <- list('global' = global.files)
        
        # create folder for contrast
        setwd(temp.path)
        if (subfolder) {
          dir.create(contrast.name)
          setwd(contrast.name)
          contrastFolder <- 
            synapser::synStore(synapser::Folder(contrast.name,
                                                parent = synapse_id))
        } else {
          contrastFolder <- synapse_id
          contrast.name <- ""
        }
        
        save_to_synapse(all.files, contrastFolder)
        
        # compile DEGs if relevant
        if (compileDEGs) {
          # get global degs
          global.degs <- deg[["global"]]
          
          # add feature_type
          global.degs$Feature_type <- colnames(global.degs)[1]
          colnames(global.degs)[1] <- "Feature"
          
          # add contrast information
          global.degs$Filter <- paste0(filterID, "_", filter)
          global.degs$Contrast <- contrast.name
          all.degs <- rbind(all.degs, global.degs)
        }
        
        # make space to process next contrast
        gsea1 <- NULL
        dmea.results <- NULL
        deg <- NULL
        gsea1.inputs <- NULL
      } 
    }
  } 
  
  if (compileDEGs) {
    if (nrow(all.degs) > 0) {
      all.DEG.files <- list("Differential_expression_results.csv" = 
                              all.degs,
                            "Differential_expression_results_max_5_percent_FDR.csv" = 
                              all.degs[all.degs$adj.P.Val <= 0.05, ]) 
    } else {
      all.DEG.files <- list("Differential_expression_results.csv" = 
                              all.degs)
    }
    save_to_synapse(all.DEG.files, synapse_id)
  }
}

run_TF_contrasts_human <- function(contrast.types, id.type, meta.df, 
                                          omics, types = c("global", "phospho"),
                                   feature.names = c("Gene", "SUB_SITE"),
                                   gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                                               "msigdb_Homo sapiens_H"),
                                          EA.types = c("KEGG", "Hallmark"),
                                   gmt.list2 = c("ksdb_human",
                                                 "sub"),
                                          expr = as.list(rep("adherent CCLE", 3)),
                                          gmt.drug = "PRISM", drug.sens = "PRISM", 
                                          base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                                          temp.path = base.path, subfolder = TRUE,
                                          synapse_id = NULL, compileDEGs = TRUE, 
                                          filter = NA, filterID = NULL) {
  setwd(base.path)
  if (!is.null(filterID)) {
    meta.df <- meta.df[meta.df[,filterID] == filter,]
    contrast.types <- contrast.types[contrast.types != filterID]
  }
  
  all.degs <- data.frame()
  if (nrow(meta.df) > 2) {
    # load set annotations
    if (gmt.list1[1] == "chr8") {
      gmt1 <- get_chr8_gmt1()
    } else {
      gmt1 <- get_gmt1(gmt.list1)
    }
    
    prot.df.noNA <- get_CCLE_prot()
    
    ## for each contrast: 
    for (k in 1:length(contrast.types)) {
      # identify samples for each side of contrast
      c1 <- "TRUE"
      c2 <- "FALSE"
      group.names <- c(c1, c2)
      contrast.name <- paste0(contrast.types[k], "_", c1, "_vs_", c2)
      group.samples <- list(as.vector(na.omit(meta.df[meta.df[,contrast.types[k]] == c1, id.type])),
                            as.vector(na.omit(meta.df[meta.df[,contrast.types[k]] == c2, id.type])))
      
      if (length(group.samples[[1]]) > 0 & length(group.samples[[2]]) > 0 & 
          length(unique(c(group.samples[[1]], group.samples[[2]]))) > 2) {
        # run panSEA across omics types
        # CAUTION: this only works because the # of samples for each treatment type 
        # is equal; otherwise would have to run panSEA for each contrast separately 
        # and perhaps set the group.samples input parameter for panSEA
        
        # run global against KEGG sets, phospho against kinase sets from ksdb
        # run differential expression analysis
        deg <- panSEA::mDEG(omics, types, group.names, group.samples, 
                            feature.names)$all.results
        gsea1.inputs <- list("global" = deg[[1]])
        prot.expr <- list(prot.df.noNA)
        features1 <- c("Gene")
        rank.var <- c("Log2FC")
        
        # run GSEA for each gmt in gmt1
        gsea1 <- list()
        global.gsea.files <- list()
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
          
          # store results for each input type
          global.mtn <- get_top_mtn_plots(
            gsea1[[gsea.name]]$all.results[["global"]], 
            EA.type = EA.types[i])
          global.net <- panSEA::netSEA(list(gsea1.inputs[["global"]]),
                                       list(gsea1[[gsea.name]]$all.results[["global"]]$result),
                                       n.network.sets = 5)
          global.gsea.files[[gsea.name]] <- list("GSEA_results.csv" =
                                                   gsea1[[gsea.name]]$all.results[["global"]]$result,
                                                 "GSEA_volcano_plot.pdf" =
                                                   gsea1[[gsea.name]]$all.results[["global"]]$volcano.plot,
                                                 "GSEA_network_graph.html" = 
                                                   global.net$interactive,
                                                 "mtn_plots" = global.mtn)
        }
        
        # run DMEA
        dmea.results <- panSEA::mDMEA(drug.sens, gmt.drug, expr, gsea1.inputs, 
                                      names(gsea1.inputs), 
                                      feature.names = features1,
                                      weight.values = rank.var)
        dmea.prot.results <- panSEA::mDMEA(drug.sens, gmt.drug, prot.expr, gsea1.inputs, 
                                           names(gsea1.inputs), 
                                           feature.names = features1,
                                           weight.values = rank.var)
        
        #### save results & upload to Synapse
        ### set file names
        ## global
        global.DEG.files <- list("Differential_expression_results.csv" = 
                                   deg[["global"]])
        DMEA.global.mtn <- get_top_mtn_plots(dmea.results$all.results[["global"]],
                                             sets = "Drug_set",
                                             EA.type = "DMEA")
        DMEA.global.net <- panSEA::netSEA(list(dmea.results$all.results[["global"]]$corr.result),
                                          list(dmea.results$all.results[["global"]]$result),
                                          "Drug", "Pearson.est",
                                          n.network.sets = 5)
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
        prot.DMEA.global.mtn <- get_top_mtn_plots(dmea.prot.results$all.results[[1]],
                                                  sets = "Drug_set",
                                                  EA.type = "DMEA")
        prot.DMEA.global.net <- panSEA::netSEA(list(dmea.prot.results$all.results[[1]]$corr.result),
                                               list(dmea.prot.results$all.results[[1]]$result),
                                               "Drug", "Pearson.est",
                                               n.network.sets = 5)
        prot.global.DMEA.files <- list("DMEA_results.csv" =
                                         dmea.prot.results$all.results[[1]]$result,
                                       "DMEA_correlation_results.csv" = 
                                         dmea.prot.results$all.results[[1]]$corr.result,
                                       "DMEA_correlation_scatter_plots.pdf" = 
                                         dmea.prot.results$all.results[[1]]$corr.scatter.plots,
                                       "DMEA_volcano_plot.pdf" =
                                         dmea.prot.results$all.results[[1]]$volcano.plot,
                                       "DMEA_network_graph.html" = 
                                         prot.DMEA.global.net$interactive,
                                       "mtn_plots" = prot.DMEA.global.mtn)
        global.files <- list('Differential_expression' = global.DEG.files, 
                             'DMEA' = global.DMEA.files,
                             'DMEA_proteomics' = prot.global.DMEA.files,
                             'GSEA' = global.gsea.files)
        
        all.files <- list('global' = global.files)
        
        # create folder for contrast
        setwd(temp.path)
        if (subfolder) {
          dir.create(contrast.name)
          setwd(contrast.name)
          contrastFolder <- 
            synapser::synStore(synapser::Folder(contrast.name,
                                                parent = synapse_id))
        } else {
          contrastFolder <- synapse_id
          contrast.name <- ""
        }
        
        save_to_synapse(all.files, contrastFolder)
        
        # compile DEGs if relevant
        if (compileDEGs) {
          # get global degs
          global.degs <- deg[["global"]]
          
          # add feature_type
          global.degs$Feature_type <- colnames(global.degs)[1]
          colnames(global.degs)[1] <- "Feature"
          
          # add contrast information
          global.degs$Filter <- paste0(filterID, "_", filter)
          global.degs$Contrast <- contrast.name
          all.degs <- rbind(all.degs, global.degs)
        }
        
        # make space to process next contrast
        gsea1 <- NULL
        dmea.results <- NULL
        deg <- NULL
        gsea1.inputs <- NULL
      } 
    }
  } 
  
  if (compileDEGs) {
    if (nrow(all.degs) > 0) {
      all.DEG.files <- list("Differential_expression_results.csv" = 
                              all.degs,
                            "Differential_expression_results_max_5_percent_FDR.csv" = 
                              all.degs[all.degs$adj.P.Val <= 0.05, ]) 
    } else {
      all.DEG.files <- list("Differential_expression_results.csv" = 
                              all.degs)
    }
    save_to_synapse(all.DEG.files, synapse_id)
  }
}

run_TF_contrast_combos_global_human <- function(contrast.types, id.type, meta.df, 
                                                        omics, gmt.list1 = c("msigdb_Homo sapiens_C2_CP:KEGG",
                                                                             "msigdb_Homo sapiens_H"),
                                                        EA.types = c("KEGG", "Hallmark"),
                                                        expr = as.list(rep("adherent CCLE", 1)),
                                                        gmt.drug = "PRISM", drug.sens = "PRISM", 
                                                        base.path = "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/",
                                                        temp.path = base.path, subfolder = TRUE,
                                                        synapse_id = NULL, compileDEGs = TRUE) {
  setwd(base.path)
  types = c("global")
  feature.names = c("Gene")
  
  # load set annotations
  if (gmt.list1[1] == "chr8") {
    gmt1 <- get_chr8_gmt1()
  } else {
    gmt1 <- get_gmt1(gmt.list1)
  }
  
  all.degs <- data.frame()
  # run contrasts with no filters
  setwd(temp.path)
  dir.create("no_filter")
  setwd("no_filter")
  nullPath <- file.path(temp.path, "no_filter")
  nullFolder <- 
    synapser::synStore(synapser::Folder("no_filter",
                                        parent = synapse_id))
  run_TF_contrasts_global_human(contrast.types, id.type, meta.df, 
                                omics, gmt.list1 = gmt1,
                                EA.types,
                                expr,
                                gmt.drug, drug.sens, 
                                base.path, temp.path = nullPath, subfolder,
                                synapse_id = nullFolder, compileDEGs)
  if (compileDEGs) {
    nullFiles <- as.list(synapser::synGetChildren(nullFolder, list("file"), sortBy = 'NAME'))
    nullFile <- synapser::synGet(nullFiles[[1]]$id)
    nullDEGs <- read.csv(nullFile$path)
    all.degs <- rbind(all.degs, nullDEGs)
  }
  
  for (m in 1:length(contrast.types)) {
    # run contrasts with TRUE filters
    setwd(temp.path)
    dir.create(file.path(paste0(contrast.types[m], "_TRUE")))
    setwd(file.path(paste0(contrast.types[m], "_TRUE")))
    truePath <- file.path(temp.path, paste0(contrast.types[m], "_TRUE"))
    trueFolder <- 
      synapser::synStore(synapser::Folder(file.path(paste0(contrast.types[m], "_TRUE")),
                                          parent = synapse_id))
    run_TF_contrasts_global_human(contrast.types, id.type, meta.df, 
                                          omics, gmt.list1 = gmt1,
                                          EA.types,
                                          expr,
                                          gmt.drug, drug.sens, 
                                          base.path, temp.path = truePath, subfolder,
                                          synapse_id = trueFolder, compileDEGs, 
                                          filter = "TRUE", filterID = contrast.types[m])
    if (compileDEGs) {
      trueDEGfiles <- as.list(synapser::synGetChildren(trueFolder, list("file"), sortBy = 'NAME'))
      trueDEGfile <- synapser::synGet(trueDEGfiles[[1]]$id)
      trueDEGs <- read.csv(trueDEGfile$path)
      all.degs <- rbind(all.degs, trueDEGs)
    }
    
    # run contrasts with FALSE filters
    setwd(temp.path)
    dir.create(file.path(paste0(contrast.types[m], "_FALSE")))
    setwd(file.path(paste0(contrast.types[m], "_FALSE")))
    falsePath <- file.path(temp.path, paste0(contrast.types[m], "_FALSE"))
    falseFolder <- 
      synapser::synStore(synapser::Folder(file.path(paste0(contrast.types[m], "_FALSE")),
                                          parent = synapse_id))
    run_TF_contrasts_global_human(contrast.types, id.type, meta.df, 
                                          omics, gmt.list1 = gmt1,
                                          EA.types,
                                          expr,
                                          gmt.drug, drug.sens, 
                                          base.path, temp.path = falsePath, subfolder,
                                          synapse_id = falseFolder, compileDEGs, 
                                          filter = "FALSE", filterID = contrast.types[m])
    if (compileDEGs) {
      falseDEGfiles <- as.list(synapser::synGetChildren(falseFolder, list("file"), sortBy = 'NAME'))
      falseDEGfile <- synapser::synGet(falseDEGfiles[[1]]$id)
      falseDEGs <- read.csv(falseDEGfile$path)
      all.degs <- rbind(all.degs, falseDEGs)
    }
  }
  
  if (compileDEGs) {
    setwd(temp.path)
    all.DEG.files <- list("Differential_expression_results.csv" = 
                            all.degs,
                          "Differential_expression_results_max_5_percent_FDR.csv" = 
                            all.degs[all.degs$adj.P.Val <= 0.05, ])
    save_to_synapse(all.DEG.files, synapse_id)
  }
}

## load BeatAML data formatted for DMEA
# input: file path where BeatAML data should be saved
# output: list of BeatAML meta, drug AUC, global, and phospho data frames
load_BeatAML_for_DMEA <- function(BeatAML.path) {
  BeatAML_synapse_id <- list("drug_response.csv" = "syn51674470", 
                             "Ex10_metadata.txt" = "syn25807733",
                             "ptrc_ex10_crosstab_global_gene_corrected.txt" = "syn25714248",
                             "ptrc_ex10_crosstab_phospho_siteID_corrected(1).txt" = "syn25714921")
  
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
  
  # remove drugs without moa annotations and drug combos
  valid.drugs <- 
    names(drug.BeatAML)[names(drug.BeatAML) %in% 
                          moa.BeatAML[!is.na(moa.BeatAML),]$Drug] # 167 drugs
  drug.BeatAML <- drug.BeatAML[ , c("sample_id", valid.drugs)] # 167 drugs
  moa.BeatAML <- 
    moa.BeatAML[moa.BeatAML$Drug %in% names(drug.BeatAML)[2:ncol(drug.BeatAML)], ]
  
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
  global.BeatAML[, sample.names] <- rownames(global.BeatAML)
  global.BeatAML <- 
    global.BeatAML[ , c(sample.names, 
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
  phospho.BeatAML[, sample.names] <- rownames(phospho.BeatAML)
  phospho.BeatAML <- phospho.BeatAML[ , c(sample.names, names(phospho.BeatAML[ , 1:(ncol(phospho.BeatAML)-1)]))]
  
  return(list(meta = meta.BeatAML, drug = drug.BeatAML,
              global = global.BeatAML, phospho = phospho.BeatAML))
}
