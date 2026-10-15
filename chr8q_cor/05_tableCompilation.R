# compile supplementary tables
library(purrr); library(data.table); #library(synapser)
#synapser::synLogin()
#removing synapse integration since many files were not uploaded to synapse
base.dir=file.path(base.path,'analysis')
net.dir = file.path(base.dir,'network')
table.dir=file.path('2026-10-05','data')

de <- list("Copy_number" = read.csv(file.path(base.dir,'Copy_number','Differential_expression','Differential_expression_results.csv')),#synapser::synGet("syn66227257")$path),
           "RNA" = read.csv(file.path(base.dir,'RNA-Seq','Differential_expression','Differential_expression_results.csv')),#synapser::synGet("syn66226866")$path),
           "Protein" = read.csv(file.path(base.dir,'Proteomics','Differential_expression','Differential_expression_results.csv')),#synapser::synGet("syn66224803")$path),
           "Phospho" = read.csv(file.path(base.dir,'Phospho','Differential_expression','Differential_expression_results.csv')))#synapser::synGet("syn66226338")$path))


##Reviewer comment: find shared features
sig <- do.call(rbind,lapply(names(de), function(vn) {
  val <- de[[vn]]
  val <- val[abs(as.numeric(val$Spearman.est))>0.75 & val$Spearman.q<0.05,]
  val$chr8Cor = val$Spearman.est>0
  val$dataType = vn
  colnames(val)[1] <- 'Gene'
  if(vn=='Phospho')
    val$Gene <- sapply(val$Gene, function(x) gsub('-.+','',x))
  return(val)}))

##munging to do the counts
duplicates <- sig |> group_by(chr8Cor,Gene) |> summarize(nvals = n_distinct(dataType)) |> subset(nvals>1) |> arrange(desc(nvals)) |> unique()
sig_source <- sig |> as.data.frame() |> drop_na() |> dplyr::select(Gene,dataType,chr8Cor)|> pivot_wider(names_from=chr8Cor, values_from = dataType)
sig_source$`TRUE` <- sapply(sig_source$`TRUE`, function(x) paste(unique(x), collapse=','))
sig_source$`FALSE` <- sapply(sig_source$`FALSE`, function(x) paste(unique(x), collapse=','))
sig_source <- sig_source |> tidyr:: pivot_longer(c(`TRUE`,`FALSE`),names_to='chr8Cor',values_to='dataTypes')|>
  subset(dataTypes!="")|> 
  subset(!dataTypes %in% c('Copy_number','RNA','Protein','Phospho'))
sig_summary <-  sig_source |> group_by(dataTypes,chr8Cor) |> summarize(n_distinct(Gene))

de <- append(list('Overlapping features'= sig_source,
                  'Overlap statistics' = sig_summary), de)

Sheet <- names(de)
readme <- data.frame(Sheet)
readme$Description <- c("Correlations between copy number and median chr8q copy number",
                        "Correlations between RNA-seq TPM and median chr8q copy number",
                        "Correlations between relative protein abundance and median chr8q copy number",
                        "Correlations between relative phospho-protein abundance and median chr8q copy number",
                        "Summary of overlap between correlated (abs Rho >0.75 and q<0.05) genes",
                        "Total number of features agreeing between data types")

de <- append(list("README" = readme), de)


gsea <- list("Copy_number" = data.table::rbindlist(list("Positional" = read.csv(file.path(base.dir,'Copy_number','GSEA','GSEA_Positional','GSEA_results.csv'))),#read.csv(synapser::synGet("syn66227265")$path)),
                                                   use.names=TRUE, idcol="Collection"),
           "RNA" = data.table::rbindlist(list("TFT_GTRD" = read.csv(file.path(base.dir,'RNA-Seq','GSEA','GSEA_TFT_GTRD','GSEA_results.csv')),#read.csv(synapser::synGet("syn66226952")$path), 
                        "Hallmark" = read.csv(file.path(base.dir,'RNA-Seq','GSEA','GSEA_Hallmark','GSEA_results.csv'))),#read.csv(synapser::synGet("syn66226874")$path)),
                        use.names=TRUE, idcol="Collection"),
           "Protein" = data.table::rbindlist(list("Hallmark" = read.csv(file.path(base.dir,'Proteomics','GSEA','GSEA_Hallmark','GSEA_results.csv'))),#read.csv(synapser::synGet("syn66224811")$path)),
                                             use.names=TRUE, idcol="Collection"))
Sheet <- names(gsea)
readme <- data.frame(Sheet)
readme$Description <- c("Gene set enrichment analysis of copy number Spearman correlations with median chr8q copy number",
                        "Gene set enrichment analysis of RNA-Seq TPM Spearman correlations with median chr8q copy number",
                        "Gene set enrichment analysis of relative protein abundance Spearman correlations with median chr8q copy number")
gsea <- append(list("README" = readme), gsea)

ksea <- list("KSEA" = read.csv(file.path(base.dir,'Phospho','Phospho_enrichment','GSEA_ksdb_human','GSEA_results.csv')))#read.csv(synapser::synGet("syn66279699")$path))
Sheet <- names(ksea)
readme <- data.frame(Sheet)
readme$Description <- c("Kinase-substrate enrichment analysis of relative phospho-protein abundance Spearman correlations with median chr8q copy number")
ksea <- append(list("README" = readme), ksea)

net.analysis <- list("Full" = data.table::rbindlist(list("Positive" = read.csv(file.path(net.dir,'positive_centrality.csv')),#,synapser::synGet("syn68869131")$path), 
                                                         "Negative" = read.csv(file.path(net.dir,'negative_centrality.csv'))),#synapser::synGet("syn68869133")$path)), 
                                                    use.names=TRUE, idcol="Direction"), 
                     "Regulatory" = data.table::rbindlist(list("Positive" = read.csv(file.path(net.dir,'positive_TFKin_centrality.csv')),#synapser::synGet("syn68869720")$path), 
                                                               "Negative" = read.csv(file.path(net.dir,'negative_TFKin_centrality.csv'))),#synapser::synGet("syn68869690")$path)), 
                                                          use.names=TRUE, idcol="Direction"))
net.nodes <- list("Full" = data.table::rbindlist(list("Positive" = read.csv(file.path(net.dir,'positive_vertices.csv')),#synapser::synGet("syn68869127")$path), 
                                                         "Negative" = read.csv(file.path(net.dir,'negative_vertices.csv'))),#synapser::synGet("syn68869130")$path)), 
                                                    use.names=TRUE, idcol="Direction"), 
                     "Regulatory" = data.table::rbindlist(list("Positive" = read.csv(file.path(net.dir,'positive_TFKin_vertices.csv')),#synapser::synGet("syn68869697")$path), 
                                                               "Negative" = read.csv(file.path(net.dir,'negative_TFKin_vertices.csv'))),#synapser::synGet("syn68869692")$path)), 
                                                          use.names=TRUE, idcol="Direction"))
net.edges <- list("Full" = data.table::rbindlist(list("Positive" = read.csv(file.path(net.dir,'positive_edges.csv')),#synapser::synGet("syn68869132")$path), 
                                                      "Negative" = read.csv(file.path(net.dir,'negative_edges.csv'))),#synapser::synGet("syn68869128")$path)), 
                                                 use.names=TRUE, idcol="Direction"), 
                  "Regulatory" = data.table::rbindlist(list("Positive" = read.csv(file.path(net.dir,'positive_TFKin_edges.csv')),#synapser::synGet("syn68869703")$path), 
                                                             "Negative" = read.csv(file.path(net.dir,'negative_TFKin_edges.csv'))),#synapser::synGet("syn68869691")$path)), 
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


##reviewer requested a link between drug mechanism and the drug name
gmt.drug <- GSA::GSA.read.gmt(synGet('syn74955107')$path)
drug_moa <- do.call(rbind,
                    lapply(1:length(gmt.drug$geneset.names), function(i){
                      gsid <- gmt.drug$geneset.descriptions[i]
                      drugs <- setdiff(gmt.drug$genesets[[i]],"")
                      return(data.frame(Drug=drugs,moa = rep(gsid,length(drugs))))
                    }))

dmea <- list("DrugCorrelations" = read.csv(file.path(base.dir,'DMEA_corr','DMEA_correlation_results.csv')),#synapser::synGet("syn66295273")$path),
           "DMEA" = read.csv(file.path(base.dir,'DMEA_corr','DMEA_results.csv')))#synapser::synGet("syn66295241")$path))
colnames(dmea[["DrugCorrelations"]])[1:2] <- c("Omics","Drug")
dmea[["DrugCorrelations"]]$feature_name <- NULL
dmea$DrugCorrelations <- dmea$DrugCorrelations |> left_join(drug_moa)
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
  openxlsx::write.xlsx(synIDs[[i]], file=paste0(table.dir,"/SupplementaryTable",i+2,"_",names(synIDs)[i],".xlsx"), rowNames=FALSE)
}
