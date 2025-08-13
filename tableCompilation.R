# compile supplementary tables
library(purrr); library(data.table); library(synapser)
synapser::synLogin()
de <- list("Copy_number" = read.csv(synapser::synGet("syn66227257")$path),
           "RNA" = read.csv(synapser::synGet("syn66226866")$path),
           "Protein" = read.csv(synapser::synGet("syn66224803")$path),
           "Phospho" = read.csv(synapser::synGet("syn66226338")$path))
Sheet <- names(de)
readme <- data.frame(Sheet)
readme$Description <- c("Correlations between copy number and median chr8q copy number",
                        "Correlations between RNA-seq TPM and median chr8q copy number",
                        "Correlations between relative protein abundance and median chr8q copy number",
                        "Correlations between relative phospho-protein abundance and median chr8q copy number")
de <- append(list("README" = readme), de)

gsea <- list("Copy_number" = data.table::rbindlist(list("Positional" = read.csv(synapser::synGet("syn66227265")$path)),
                                                   use.names=TRUE, idcol="Collection"),
           "RNA" = data.table::rbindlist(list("TFT_GTRD" = read.csv(synapser::synGet("syn66226952")$path), 
                        "Hallmark" = read.csv(synapser::synGet("syn66226874")$path)),
                        use.names=TRUE, idcol="Collection"),
           "Protein" = data.table::rbindlist(list("Hallmark" = read.csv(synapser::synGet("syn66224811")$path)),
                                             use.names=TRUE, idcol="Collection"))
Sheet <- names(gsea)
readme <- data.frame(Sheet)
readme$Description <- c("Gene set enrichment analysis of copy number Spearman correlations with median chr8q copy number",
                        "Gene set enrichment analysis of RNA-Seq TPM Spearman correlations with median chr8q copy number",
                        "Gene set enrichment analysis of relative protein abundance Spearman correlations with median chr8q copy number")
gsea <- append(list("README" = readme), gsea)

ksea <- list("KSEA" = read.csv(synapser::synGet("syn66279699")$path))
Sheet <- names(ksea)
readme <- data.frame(Sheet)
readme$Description <- c("Kinase-substrate enrichment analysis of relative phospho-protein abundance Spearman correlations with median chr8q copy number")
ksea <- append(list("README" = readme), ksea)

net.analysis <- list("Full" = data.table::rbindlist(list("Positive" = read.csv(synapser::synGet("syn68869131")$path), 
                                                         "Negative" = read.csv(synapser::synGet("syn68869133")$path)), 
                                                    use.names=TRUE, idcol="Direction"), 
                     "Regulatory" = data.table::rbindlist(list("Positive" = read.csv(synapser::synGet("syn68869720")$path), 
                                                               "Negative" = read.csv(synapser::synGet("syn68869690")$path)), 
                                                          use.names=TRUE, idcol="Direction"))
net.nodes <- list("Full" = data.table::rbindlist(list("Positive" = read.csv(synapser::synGet("syn68869127")$path), 
                                                         "Negative" = read.csv(synapser::synGet("syn68869130")$path)), 
                                                    use.names=TRUE, idcol="Direction"), 
                     "Regulatory" = data.table::rbindlist(list("Positive" = read.csv(synapser::synGet("syn68869697")$path), 
                                                               "Negative" = read.csv(synapser::synGet("syn68869692")$path)), 
                                                          use.names=TRUE, idcol="Direction"))
net.edges <- list("Full" = data.table::rbindlist(list("Positive" = read.csv(synapser::synGet("syn68869132")$path), 
                                                      "Negative" = read.csv(synapser::synGet("syn68869128")$path)), 
                                                 use.names=TRUE, idcol="Direction"), 
                  "Regulatory" = data.table::rbindlist(list("Positive" = read.csv(synapser::synGet("syn68869703")$path), 
                                                             "Negative" = read.csv(synapser::synGet("syn68869691")$path)), 
                                                       use.names=TRUE, idcol="Direction"))
net <- list("Analysis" = data.table::rbindlist(net.analysis, use.names=TRUE, idcol="Network"),
            "Nodes" = data.table::rbindlist(net.nodes, use.names=TRUE, idcol="Network"),
            "Edges" = data.table::rbindlist(net.edges, use.names=TRUE, idcol="Network"))
Sheet <- names(net)
readme <- data.frame(Sheet)
readme$Description <- c("Network analysis including node degree, closeness, betweenness, eigen centrality, hub score, and authority score",
               "Information about nodes",
               "Information about edges")
net <- append(list("README" = readme), net)

dmea <- list("DrugCorrelations" = read.csv(synapser::synGet("syn66295273")$path),
           "DMEA" = read.csv(synapser::synGet("syn66295241")$path))
colnames(dmea[["DrugCorrelations"]])[1:2] <- c("Omics","Drug")
dmea[["DrugCorrelations"]]$feature_name <- NULL
colnames(dmea[["DMEA"]])[1] <- "Omics"
Sheet <- names(dmea)
readme <- data.frame(Sheet)
readme$Description <- c("Drug correlations between drug sensitivity (i.e., AUC) and CCLE scores based on correlations with median chr8q copy number",
                        "Drug mechanism enrichment analysis results")
dmea <- append(list("README" = readme), dmea)

synIDs <- list("Correlations" = de,
               "GSEA" = gsea,
               "KSEA" = ksea,
               "Network_analyses" = net,
               "DMEA" = dmea)
for (i in 1:length(synIDs)) {
  openxlsx::write.xlsx(synIDs[[i]], file=paste0("SupplementaryTable",i,"_",names(synIDs)[i],".xlsx"), rowNames=FALSE)
}
