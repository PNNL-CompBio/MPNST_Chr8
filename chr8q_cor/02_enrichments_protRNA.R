
# next, proteomics and RNA
setwd(paste0(base.path,'/analysis'))
if(!exists('gmt1'))
  gmt1 <-readRDS(synapser::synGet('syn72245876')$path)
#gmt1 <- get_gmt1_v2(gmt.list1=c("msigdb_Homo sapiens_HS_H","msigdb_Homo sapiens_HS_C1","msigdb_Homo sapiens_HS_C3_TFT:GTRD"),
#                   names1=c("Hallmark","Positional", "TFT_GTRD"))
#updated the phospho list to have the correct length

#gmt2 <- get_gmt2(gmt.list2=c("ksdb_human",'sub'), phospho=phospho.df)
if(!exists('gmt2'))
  gmt2 <- list(`ksdb_human`=get_leapR_ksdb_gmt(), `phospho`=phospho.df)

#setwd("Chr8_quant_20250409")
omics <- list( "Phospho" = phospho.df,
               "Proteomics" = global.df,#list("Global" = global.df,"Phospho" = phospho.df),
              # "Proteomics" = list("Phospho" = phospho.df),
              "RNA-Seq" = pdxRNA)
meta.list <- list("Phospho" = global.meta.df3,
              "Proteomics" = global.meta.df3,
                  "RNA-Seq" = pdx.info2[pdx.info2$Sample != "JH-2-009",])
expr.list <- list("Phospho" = "CCLE proteomics",
                  "Proteomics" = "CCLE proteomics",
                  "RNA-Seq" = "adherent CCLE")
feature.list <- list( "Phospho" = "SUB_SITE",
                      "Proteomics" = "Gene",# "SUB_SITE"),
                     "RNA-Seq" = "Gene")
gmt1.rest <- list("Hallmark" = gmt1$Hallmark,
                  "KEGG" = gmt1$KEGG,
                  "Oncogenic" = gmt1$Oncogenic,
                  "PID" = gmt1$PID,
                  "TFT_GTRD" = gmt1$TFT_GTRD,
                  "WikiPathways" = gmt1$WikiPathways)
set.seed(12340)
ofiles <- panSEA_corr3(omics, meta.list, feature.list, 
             rank.col = "Median Chr8q Copy Number",
             other.annotations = c("Sex", "PRC2 Status"), expr.list = expr.list, gmt1=gmt1.rest, gmt2=gmt2,
             temp.path = file.path(base.path, "analysis"), syn.id = NULL)#my.syn)

graphics.off()
##SG: let's manually update files of interest to synapse instead of trying to do it all
setwd(paste0(base.path,'/analysis'))


# re-do KSEA
#gmt2 <- get_gmt2(gmt.list2="ksdb_human", phospho=phospho.df)
# load correlations
corr.result <- read.csv(file.path(base.path,'analysis','Phospho','Phospho','Differential_expression','Differential_expression_results.csv'))#read.csv(synapser::synGet("syn66226338")$path)
gsea.input <- corr.result[,c("SUB_SITE","Spearman.est")]

# run KSEA
set.seed(2361)
gsea.result <- panSEA::ssGSEA(gsea.input, gmt2[[1]], ties=TRUE) 

# save results
gsea.files <- list("KSEA_results.csv" = gsea.result$result,
                   "KSEA_results_withoutShufflingTies.csv" = gsea.result$result.w.ties,
                   "KSEA_volcano_plot.pdf" = gsea.result$volcano.plot,
                   "KSEA_dot_plot.pdf" = gsea.result$dot.plot,
                   "KSEA_dot_plot_withSD.pdf" = gsea.result$dot.sd,
                   "KSEA_bar_plot.pdf" = gsea.result$bar.plot,
                   "mtn_plots" = get_top_mtn_plots(gsea.result))

# setwd(file.path(base.path,"analysis","Spearman","Phospho", "Phospho"))
# dir.create("KSEA")
# setwd("KSEA")
# gseaFolder <- synapser::synStore(synapser::Folder("KSEA", parent = "syn66226336"))
save_to_synapse_v2(gsea.files, NULL)#gseaFolder)
