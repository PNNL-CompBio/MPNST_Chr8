#New Figure S2 that compares ties vs  no-ties across calculations




##define helper function to join tied and untied together, plot correlations,
## and identify differences
check_ties <- function(wties_file, noties_file,data_type=''){
  wties<-readr::read_csv(wties_file)
  colnames(wties)[c(3:11)] <- paste('Shuffled',colnames(wties)[c(3:11)])
  
  noties<-readr::read_csv(noties_file)
  colnames(noties)[c(3:11)] <- paste('Original',colnames(noties)[c(3:11)])
  
  allres<- wties |> left_join(noties)
  
  ##ES CORRELATION
  p1<-allres |> ggplot(aes(x=`Original ES`,y=`Shuffled ES`)) + 
    geom_point()+#aes(col=values.to.vary,size=as.numeric(SetSize)))  +
    geom_smooth(method = "lm", se = FALSE, col='darkgrey')  +# Adds the linear regression line
    stat_cor(method = "pearson", r.digits=8,label.x.npc = "left", label.y.npc = "top") + 
    theme_minimal()
  
  ##NES CORRELATION
  p2<-allres |> ggplot(aes(x=`Original NES`,y=`Shuffled NES`)) + 
    geom_point()  +
    geom_smooth(method = "lm", se = FALSE, col='darkgrey')  +# Adds the linear regression line
    stat_cor(method = "pearson", r.digits=8,label.x.npc = "left", label.y.npc = "top") + 
    theme_minimal()
  
  #pvalue correlation
  ##for p-values, we want to check if there are conflicts
  allres$Outliers = ""
  ol <- allres$`Shuffled p_value`>0.05&allres$`Original p_value`<0.05
  il <- allres$`Shuffled p_value`<0.05&allres$`Original p_value`>0.05
  
  allres$Outliers[ol|il]<-allres$Feature_set[ol|il]
  allres$Color = ''
  allres$Color[ol|il] = 'outlier'
  # print(allres$Outliers)
  p3<-allres |> ggplot(aes(x=`Original p_value`,y=`Shuffled p_value`)) + 
    geom_point(aes(color=Color))+#aes(col=values.to.vary,size=as.numeric(SetSize)))  +
    ggrepel::geom_text_repel(aes(label=Outliers),max.overlaps = Inf)+
    scale_color_manual(values=c('black','darkred'))+
    geom_hline(yintercept =.05,col='red',linetype='dashed')+
    geom_vline(xintercept=.05,col='red',linetype='dashed')+
    geom_smooth(method = "lm", se = FALSE, col='darkgrey')  +# Adds the linear regression line
    stat_cor(method = "pearson", r.digits=8,label.x.npc = "left", label.y.npc = "top") + 
    theme_minimal()+
    theme(plot.margin = margin(10, 30, 10, 10),legend.position = "none")
  
  ##FDR correlation
  
  allres$Outliers = ""
  ol <- allres$`Shuffled FDR_q_value`>0.25&allres$`Original FDR_q_value`<0.25
  il <- allres$`Shuffled FDR_q_value`<0.25&allres$`Original FDR_q_value`>0.25
  
  allres$Outliers[ol|il]<-allres$Feature_set[ol|il]
  allres$Color = 'black'
  allres$Color[ol|il] = 'red'
  # print(allres$Outliers)
  
  p4<-allres |> ggplot(aes(x=`Original FDR_q_value`,y=`Shuffled FDR_q_value`)) + 
    geom_point(aes(color=Color))+#aes(col=values.to.vary,size=as.numeric(SetSize)))  +
    ggrepel::geom_text_repel(aes(label=Outliers),max.overlaps = Inf)+
    scale_color_manual(values=c('black','darkred'))+
    geom_hline(yintercept =.25,col='red',linetype='dashed')+
    geom_vline(xintercept=.25,col='red',linetype='dashed')+
    geom_smooth(method = "lm", se = FALSE, col='darkgrey')  +# Adds the linear regression line
    stat_cor(group = 1,method = "pearson", r.digits=8,label.x.npc = "left", label.y.npc = "top") + 
    theme_minimal()+
    theme(plot.margin = margin(10, 30, 10, 10),legend.position = "none")
  
  cowplot::plot_grid(plotlist = c(p1,p2,p3,p4),nrow=2,rel_widths = c(1, 1.4))
  
}


##can we do the same on original gene sets?
##let's just do the chromosomal enrichment
cnv_plots <- check_ties(wties_file=file.path(base.path,"analysis/Copy_number/GSEA/GSEA_Positional/GSEA_results.csv"),
                        noties_file =file.path(base.path,"analysis/Copy_number/GSEA/GSEA_Positional/GSEA_results_withoutShufflingTies.csv"),
                        data_type='Copy number')

rna_plots <-check_ties(wties_file=file.path(base.path,"analysis/RNA-Seq/GSEA/GSEA_Hallmark/GSEA_results.csv"),
                       noties_file =file.path(base.path,"analysis/RNA-Seq/GSEA/GSEA_Hallmark/GSEA_results_withoutShufflingTies.csv"),
                       data_type='RNA')

tf_plots <-check_ties(wties_file=file.path(base.path,"analysis/RNA-Seq/GSEA/GSEA_TFT_GTRD/GSEA_results.csv"),
                      noties_file =file.path(base.path,"analysis/RNA-Seq/GSEA/GSEA_TFT_GTRD/GSEA_results_withoutShufflingTies.csv"),
                      data_type='RNA_TF') 
prot_plots <-check_ties(wties_file=file.path(base.path,"analysis/Proteomics/GSEA/GSEA_Hallmark/GSEA_results.csv"),
                        noties_file =file.path(base.path,"analysis/Proteomics/GSEA/GSEA_Hallmark/GSEA_results_withoutShufflingTies.csv"),
                        data_type='Protein')

kin_plots <-check_ties(wties_file=file.path(base.path,"analysis/Phospho/Phospho_enrichment/GSEA_ksdb_human/GSEA_results.csv"),
                       noties_file =file.path(base.path,"analysis/Phospho/Phospho_enrichment/GSEA_ksdb_human/GSEA_results_withoutShufflingTies.csv"),
                       data_type='Protein')

ggsave(cnv_plots,file=file.path(base.path,'figures/fig_s2_cnv_tie_shuffles.pdf'), width =8)
ggsave(rna_plots,file=file.path(base.path,'figures/fig_s2_rna_tie_shuffles.pdf'), width =8)
ggsave(prot_plots,file=file.path(base.path,'figures/fig_s2_prot_tie_shuffles.pdf'), width =8)
ggsave(tf_plots,file=file.path(base.path,'figures/fig_s2_tf_tie_shuffles.pdf'), width =8)

ggsave(kin_plots,file=file.path(base.path,'figures/fig_s2_kin_tie_shuffles.pdf'), width =8)
#ggsave(rna_drg_plots,file=file.path(base.path,'figures/fig_s2_rna_drug_tie_shuffles.pdf'), width =8)
#ggsave(prot_drg_plots,file=file.path(base.path,'figures/fig_s2_prot__drug_tie_shuffles.pdf'), width = 8)


