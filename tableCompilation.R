# compile supplementary tables
de <- list("Copy_number" = "",
           "RNA" = "",
           "Protein" = "",
           "Phospho" = "")
Sheet <- names(de)
readme <- data.frame(Sheet)
readme$Description <- c("Correlations between copy number and median chr8q copy number",
                        "Correlations between RNA-seq TPM and median chr8q copy number",
                        "Correlations between relative protein abundance and median chr8q copy number",
                        "Correlations between relative phospho-protein abundance and median chr8q copy number")
de <- append(list("README" = readme), de)

gsea <- list("Copy_number" = "",
           "RNA" = "",
           "Protein" = "",
           "Phospho" = "")
Sheet <- names(gsea)
readme <- data.frame(Sheet)
readme$Description <- c("Gene set enrichment analysis of copy number Spearman correlations with median chr8q copy number",
                        "Gene set enrichment analysis of RNA-Seq TPM Spearman correlations with median chr8q copy number",
                        "Gene set enrichment analysis of relative protein abundance Spearman correlations with median chr8q copy number",
                        "Kinase-substrate enrichment analysis of relative phospho-protein abundance Spearman correlations with median chr8q copy number")
gsea <- append(list("README" = readme), gsea)

net <- list("Analysis" = "",
            "Nodes" = "",
            "Edges" = "")
Sheet <- names(net)
readme <- data.frame(Sheet)
readme$Description <- c("Network analysis including node degree, closeness, betweenness, eigen centrality, hub score, and authority score",
               "Information about nodes including frequency in Prize Collecting Steiner Forest randomization, output prize, input prize, and node type",
               "Information about edges including weight")
net <- append(list("README" = readme), net)

dmea <- list("RNA" = "",
           "Protein" = "")
Sheet <- names(dmea)
readme <- data.frame(Sheet)
readme$Description <- c("Drug correlations between drug sensitivity (i.e., AUC) and CCLE scores based on correlations with median chr8q copy number",
                        "Drug mechanism enrichment analysis results")
dmea <- append(list("README" = readme), dmea)

synIDs <- list("Correlations" = de,
               "Enrichment_analyses" = gsea,
               "Network_analyses" = net,
               "Drug_sensitivity_predictions" = dmea)
for (i in 1:length(synIDs)) {
  openxlsx::write.xlsx(synIDs[[i]], file=paste0("SupplementaryTable",i,"_",names(synIDs)[i],".xlsx"), rowNames=FALSE)
}