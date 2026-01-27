# Author: Belinda B. Garana
# can GSEA_ties detect enrichment of synthetic gene set?
### Step 1: Generate gene rank metrics that have a normal distribution, then perturb 1 set of X genes
### Step 2: Run GSEA_ties for each different variant (+/- Y gene rank)
### Step 3: Re-run for 50 replicates
### Step 4: Repeat with different values of X (5, 10, 20, 40)

rm(list=ls(all=TRUE))
if(!require(devtools)){install.packages("dev.tools")}
if(!require(DMEA)){devtools::install_github('BelindaBGarana/DMEA')}
library(DMEA);library(dplyr);library(GSA);library(reshape2);library(data.table);library(ggplot2);library(msigdbr)



### Step 0: Prep rows, columns, values
##where does this file come from!!??
plot.data <- read.csv("Proteomics/Global/Differential_expression/Differential_expression_results.csv")
plot.data <- na.omit(plot.data[,c("Gene", "Spearman.est")])
plot.data <- plot.data[plot.data$Spearman.est != 0,]

# shuffle gene names so that no gene set should be enriched unless synthetically
plot.data$Gene <- sample(plot.data$Gene)
Gene <- sample(plot.data$Gene)
gene.names <- Gene

# # model distribution based on real correlation
# hist(plot.data$Spearman.est)
# sd.spearman <- sd(plot.data$Spearman.est) # 0.4421582
# mean.spearman <- mean(plot.data$Spearman.est) # -0.06911027
# norm.dist <- rnorm(nrow(plot.data), mean=mean.spearman, sd = sd.spearman)
# plot(density(plot.data$Spearman.est))
# lines(density(norm.dist), col="red") # close enough

# median number of ties is 90, avg is 101.032287
# create normal distribution of 100 and then add 90 duplicates for total of 9000
Spearman.est <- rnorm(100, mean = 0, sd = 0.5)
Spearman.est <- Spearman.est[abs(Spearman.est) < 1]
Spearman.est <- rep(Spearman.est, length.out = length(Gene))
norm.df <- data.frame(Gene, Spearman.est)

# How far do we want to vary values here?
abs.vary <- 0.25
values.to.vary <- round(seq(from = -abs.vary, to = abs.vary, by = abs.vary/5), digits=2)
synthetic.cell.names <- "Spearman.est"

# Import MOA set information
gmt.info <- msigdbr::msigdbr(category="H") # hallmark gene sets for homo sapiens
gmt <- list("genesets" = list(), "geneset.names" = list(), "geneset.descriptions" = list())
gs <- unique(gmt.info$gs_name)
for (i in 1:length(gs)) {
  gmt$genesets[[i]] <- gmt.info[gmt.info$gs_name == gs[i],]$gene_symbol
  gmt$geneset.names[[i]] <- gs[i]
  gmt$geneset.descriptions[[i]] <- gs[i]
}

path.sim <- paste0("GSEA_ties_sim_results_",Sys.Date())
dir.create(path.sim)
# Change path to where you want files saved
setwd(path.sim)
setwd("valueAdd0.25")
og.all.synth.sim.results <- list()
og.all.sim.results <- list()
og.all.perc.sig <- list()

all.synth.sim.results <- list()
all.sim.results <- list()
all.perc.sig <- list()
synth.gene.sets <- data.frame(N_genes=integer(), synthetic.gene.set=character())
gene.set.sizes <- seq(30, 300, length.out=10)
for(i in 1:(length(gene.set.sizes))){
  synthetic.gene.set <- sample(gene.names,gene.set.sizes[i])
  synthetic.gene.df <- as.data.frame(synthetic.gene.set)
  synthetic.gene.df$N_genes <- gene.set.sizes[i]
  synth.gene.sets <- rbind(synth.gene.sets, synthetic.gene.df)
 
  # make new gmt with synthetic gene set
  new.gene.sets <- gmt$genesets
  new.gene.set.names <- gmt$geneset.names
  new.gene.set.descr <- gmt$geneset.descriptions
  new.gene.set.names[[length(gmt$genesets)+1]] <- paste0(gene.set.sizes[i],"_random_genes")
  new.gene.set.descr[[length(gmt$genesets)+1]] <- paste0(gene.set.sizes[i],"_random_genes")
  new.gene.sets[[length(gmt$genesets)+1]] <- synthetic.gene.set
  new.gmt <- list(genesets=new.gene.sets, geneset.names=new.gene.set.names, geneset.descriptions=new.gene.set.descr)
  
  og.synth.sim.results <- list()
  og.sim.results <- list()
  synth.sim.results <- list()
  sim.results <- list()
  for(N in 1:20){
    ### Step 3: Generate gene AUC that have a normal distribution for each cell line, then perturb synthetic gene set
    og.synthetic.results <- as.data.frame(values.to.vary)
    og.synthetic.results[,c("ES","NES","p","q")] <- NA
    og.results <- list()
    
    synthetic.results <- as.data.frame(values.to.vary)
    synthetic.results[,c("ES","NES","p","q")] <- NA
    results <- list()
    for(j in 1:length(values.to.vary)){
      rank.matrix <- norm.df
      row.names(rank.matrix) <- rank.matrix$Gene
      
      # perturb synthetic gene set
      curr.vals <- rank.matrix[synthetic.gene.set,synthetic.cell.names]
      rank.matrix[synthetic.gene.set,synthetic.cell.names] <- curr.vals + values.to.vary[j]
      
      # run GSEA accounting for ties
      DMEA.result <- drugSEA_ties(data = rank.matrix, drug = "Gene", gmt = new.gmt, rank.metric=synthetic.cell.names, plots=FALSE, min.per.set=5)
      # i = 1 (5 genes), j=7 (varied val +0.75), N=2: Error in if (temp.NES >= 0) { : missing value where TRUE/FALSE needed
      og.DMEA.result <- DMEA.result$result.w.ties
      og.synthetic.results[j, c("ES", "NES", "p", "q")] <- og.DMEA.result[
        og.DMEA.result$Drug_set==paste0(gene.set.sizes[i],"_random_genes"), c("ES", "NES", "p_value", "FDR_q_value")]
      og.results[[as.character(values.to.vary[j])]] <- og.DMEA.result
      
      DMEA.result <- DMEA.result$result
      synthetic.results[j, c("ES", "NES", "p", "q")] <- DMEA.result[
        DMEA.result$Drug_set==paste0(gene.set.sizes[i],"_random_genes"), c("ES", "NES", "p_value", "FDR_q_value")]
      results[[as.character(values.to.vary[j])]] <- DMEA.result
    }
    og.synth.sim.results[[N]] <- og.synthetic.results 
    og.sim.results[[N]] <- data.table::rbindlist(og.results, use.names=TRUE, idcol="Value Added", fill=TRUE)
    
    synth.sim.results[[N]] <- synthetic.results 
    sim.results[[N]] <- data.table::rbindlist(results, use.names=TRUE, idcol="Value Added", fill=TRUE)
  }
  og.all.synth.sim.results[[as.character(gene.set.sizes[i])]] <- data.table::rbindlist(og.synth.sim.results, use.names=TRUE, idcol="Simulation Number", fill=TRUE)
  og.all.sim.results[[as.character(gene.set.sizes[i])]] <- data.table::rbindlist(og.sim.results, use.names=TRUE, idcol="Simulation Number", fill=TRUE)
  write.csv(og.all.sim.results[[as.character(gene.set.sizes[i])]],file=paste0(gene.set.sizes[i],"_genes_",N,"_replicates_all_GSEA_ties_results_originalOrder_",Sys.Date(),".csv"))
  write.csv(og.all.synth.sim.results[[as.character(gene.set.sizes[i])]],file=paste0(gene.set.sizes[i],"_genes_",N,"_replicates_Synthetic_gene_Set_GSEA_ties_results_originalOrder_",Sys.Date(),".csv"))
  og.all.synth.sim.results.df <- og.all.synth.sim.results[[as.character(gene.set.sizes[i])]]
  
  all.synth.sim.results[[as.character(gene.set.sizes[i])]] <- data.table::rbindlist(synth.sim.results, use.names=TRUE, idcol="Simulation Number", fill=TRUE)
  all.sim.results[[as.character(gene.set.sizes[i])]] <- data.table::rbindlist(sim.results, use.names=TRUE, idcol="Simulation Number", fill=TRUE)
  write.csv(all.sim.results[[as.character(gene.set.sizes[i])]],file=paste0(gene.set.sizes[i],"_genes_",N,"_replicates_all_GSEA_ties_results_",Sys.Date(),".csv"))
  write.csv(all.synth.sim.results[[as.character(gene.set.sizes[i])]],file=paste0(gene.set.sizes[i],"_genes_",N,"_replicates_Synthetic_gene_Set_GSEA_ties_results_",Sys.Date(),".csv"))
  all.synth.sim.results.df <- all.synth.sim.results[[as.character(gene.set.sizes[i])]]
  
  perc.sig <- as.data.frame(values.to.vary)
  perc.sig[,c("N.replicates","avg.NES","sd.NES",
              "perc.sig","perc.sig.pos","perc.sig.neg",
              "perc.moreSig","perc.moreSig.pos","perc.moreSig.neg")] <- NA
  og.perc.sig <- perc.sig
  for(k in 1:length(values.to.vary)){
    temp.data <- all.synth.sim.results.df[all.synth.sim.results.df$values.to.vary==values.to.vary[k],]
    og.temp.data <- og.all.synth.sim.results.df[og.all.synth.sim.results.df$values.to.vary==values.to.vary[k],]
    
    # N and NES
    og.perc.sig$N.replicates[k] <- nrow(og.temp.data)
    og.perc.sig$avg.NES[k] <- mean(as.numeric(og.temp.data$NES))
    og.perc.sig$sd.NES[k] <- sd(as.numeric(og.temp.data$NES))
    
    perc.sig$N.replicates[k] <- nrow(temp.data)
    perc.sig$avg.NES[k] <- mean(as.numeric(temp.data$NES))
    perc.sig$sd.NES[k] <- sd(as.numeric(temp.data$NES))
    
    # % significant
    ## og
    # p < 0.05 & q < 0.25
    og.sig.temp.data <- og.temp.data[og.temp.data$p<0.05 & og.temp.data$q<0.25,]
    og.pos.sig.temp.data <- og.sig.temp.data[og.sig.temp.data$NES>0,]
    og.neg.sig.temp.data <- og.sig.temp.data[og.sig.temp.data$NES<0,]
    og.perc.sig$perc.sig[k] <- 100*nrow(og.sig.temp.data)/nrow(og.temp.data)
    og.perc.sig$perc.sig.pos[k] <- 100*nrow(og.pos.sig.temp.data)/nrow(og.temp.data)
    og.perc.sig$perc.sig.neg[k] <- 100*nrow(og.neg.sig.temp.data)/nrow(og.temp.data)
    
    # p < 0.05 & q < 0.25
    og.sig.temp.data <- og.temp.data[og.temp.data$p<0.05 & og.temp.data$q<0.05,]
    og.pos.sig.temp.data <- og.sig.temp.data[og.sig.temp.data$NES>0,]
    og.neg.sig.temp.data <- og.sig.temp.data[og.sig.temp.data$NES<0,]
    og.perc.sig$perc.moreSig[k] <- 100*nrow(og.sig.temp.data)/nrow(og.temp.data)
    og.perc.sig$perc.moreSig.pos[k] <- 100*nrow(og.pos.sig.temp.data)/nrow(og.temp.data)
    og.perc.sig$perc.moreSig.neg[k] <- 100*nrow(og.neg.sig.temp.data)/nrow(og.temp.data)
    
    ## shuffled ties
    # p < 0.05 & q < 0.25
    sig.temp.data <- temp.data[temp.data$p<0.05 & temp.data$q<0.25,]
    pos.sig.temp.data <- sig.temp.data[sig.temp.data$NES>0,]
    neg.sig.temp.data <- sig.temp.data[sig.temp.data$NES<0,]
    perc.sig$perc.sig[k] <- 100*nrow(sig.temp.data)/nrow(temp.data)
    perc.sig$perc.sig.pos[k] <- 100*nrow(pos.sig.temp.data)/nrow(temp.data)
    perc.sig$perc.sig.neg[k] <- 100*nrow(neg.sig.temp.data)/nrow(temp.data)
    
    # p < 0.05 & q < 0.25
    sig.temp.data <- temp.data[temp.data$p<0.05 & temp.data$q<0.05,]
    pos.sig.temp.data <- sig.temp.data[sig.temp.data$NES>0,]
    neg.sig.temp.data <- sig.temp.data[sig.temp.data$NES<0,]
    perc.sig$perc.moreSig[k] <- 100*nrow(sig.temp.data)/nrow(temp.data)
    perc.sig$perc.moreSig.pos[k] <- 100*nrow(pos.sig.temp.data)/nrow(temp.data)
    perc.sig$perc.moreSig.neg[k] <- 100*nrow(neg.sig.temp.data)/nrow(temp.data)
  }
  write.csv(og.perc.sig,file=paste0(gene.set.sizes[i],"_genes_",N,"_replicates_percent_significant_synthetic_gene_set_enrichments_originalOrder.csv"))
  og.all.perc.sig[[as.character(gene.set.sizes[i])]] <- og.perc.sig
  write.csv(perc.sig,file=paste0(gene.set.sizes[i],"_genes_",N,"_replicates_percent_significant_synthetic_gene_set_enrichments.csv"))
  all.perc.sig[[as.character(gene.set.sizes[i])]] <- perc.sig
}
og.complete.synth.sim.results <- data.table::rbindlist(og.all.synth.sim.results, use.names=TRUE, idcol="Number of genes in Set", fill=TRUE)
og.complete.sim.results <- data.table::rbindlist(og.all.sim.results, use.names=TRUE, idcol="Number of genes in Set", fill=TRUE)
write.csv(og.complete.sim.results,file=paste0(N,"_replicates_all_GSEA_ties_results_simulation_originalOrder_",Sys.Date(),".csv"), row.names=FALSE)
write.csv(og.complete.synth.sim.results,file=paste0(N,"_replicates_synthetic_gene_Set_GSEA_ties_results_originalOrder_",Sys.Date(),".csv"), row.names=FALSE)
og.complete.perc.sig <- data.table::rbindlist(og.all.perc.sig, use.names=TRUE, idcol="Number of genes in Set", fill=TRUE)
write.csv(og.complete.perc.sig,file=paste0(N,"_replicates_percent_significant_synthetic_gene_set_enrichments_originalOrder_",Sys.Date(),".csv"), row.names=FALSE)

complete.synth.sim.results <- data.table::rbindlist(all.synth.sim.results, use.names=TRUE, idcol="Number of genes in Set", fill=TRUE)
complete.sim.results <- data.table::rbindlist(all.sim.results, use.names=TRUE, idcol="Number of genes in Set", fill=TRUE)
write.csv(complete.sim.results,file=paste0(N,"_replicates_all_GSEA_ties_results_simulation_",Sys.Date(),".csv"), row.names=FALSE)
write.csv(complete.synth.sim.results,file=paste0(N,"_replicates_synthetic_gene_Set_GSEA_ties_results_",Sys.Date(),".csv"), row.names=FALSE)
complete.perc.sig <- data.table::rbindlist(all.perc.sig, use.names=TRUE, idcol="Number of genes in Set", fill=TRUE)
write.csv(complete.perc.sig,file=paste0(N,"_replicates_percent_significant_synthetic_gene_set_enrichments_",Sys.Date(),".csv"), row.names=FALSE)


## Generate tile plots
# prepare collapsed results
og.complete.perc.sig$N_genes <- as.numeric(og.complete.perc.sig$`Number of genes in Set`)
og.complete.perc.sig$values.to.vary <- as.numeric(og.complete.perc.sig$values.to.vary)

complete.perc.sig$N_genes <- as.numeric(complete.perc.sig$`Number of genes in Set`)
complete.perc.sig$values.to.vary <- as.numeric(complete.perc.sig$values.to.vary)

# load theme for plot
ng.theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.border = element_rect(fill=NA), panel.background = element_blank(),
                  axis.line = element_line(colour = "black"), axis.text.x = element_text(colour = "black"),
                  axis.text.y = element_text(colour = "black"), axis.ticks.x = element_line(colour="black"),
                  axis.ticks.y = element_line(colour="black"), legend.title = element_blank(),
                  axis.title.y = element_text(size=8, colour="black"))

# produce tile plots
cutoff <- c(0.05, 0.25)
for (i in cutoff) {
  temp.perc.col <- ifelse(i == 0.05, "perc.moreSig", "perc.sig")
  temp.cols <- c("N_genes", "values.to.vary", temp.perc.col)
  # original order
  temp.og <- og.complete.perc.sig[,..temp.cols]
  colnames(temp.og)[3] <- "perc"
  q.tile.plot <- ggplot(temp.og, aes(x = as.factor(N_genes), y = values.to.vary, fill = perc)) +
    geom_tile(height=0.25) + ng.theme + theme_light() + scale_fill_gradient2() +
    scale_x_discrete(breaks = gene.set.sizes) + 
    labs(x="Number of Genes in Set", y = "Value Added", fill = paste0("% q < ", i))
  ggsave(q.tile.plot, file = paste0("GSEA_ties_sim_p0.05_q",i,"_N",N,"_originalOrder_",Sys.Date(),".pdf"))
  
  # shuffled ties
  temp.df <- complete.perc.sig[,..temp.cols]
  colnames(temp.df)[3] <- "perc"
  q.tile.plot <- ggplot(temp.df, aes(x = as.factor(N_genes), y = values.to.vary, fill = perc)) +
    geom_tile(height=0.25) + ng.theme + theme_light() + scale_fill_gradient2() +
    scale_x_discrete(breaks = gene.set.sizes) + 
    labs(x="Number of Genes in Set", y = "Value Added", fill = paste0("% q < ", i))
  ggsave(q.tile.plot, file = paste0("GSEA_ties_sim_p0.05_q",i,"_N",N,"_",Sys.Date(),".pdf"))
}

NES.tile.plot <- ggplot(og.complete.perc.sig, aes(x = as.factor(N_genes), y = values.to.vary, fill=avg.NES)) +
  geom_tile(height=0.25) + ng.theme + theme_light() + scale_fill_gradient2() +
  scale_x_discrete(breaks = gene.set.sizes) +  
  labs(x="Number of Genes in Set", y = "Value Added", fill = "Mean NES")
ggsave(NES.tile.plot, file = paste0("GSEA_ties_sim_NES_N",N,"_originalOrder_",Sys.Date(),".pdf"))

NES.tile.plot <- ggplot(complete.perc.sig, aes(x = as.factor(N_genes), y = values.to.vary, fill=avg.NES)) +
  geom_tile(height=0.25) + ng.theme + theme_light() + scale_fill_gradient2() +
  scale_x_discrete(breaks = gene.set.sizes) +  
  labs(x="Number of Genes in Set", y = "Value Added", fill = "Mean NES")
ggsave(NES.tile.plot, file = paste0("GSEA_ties_sim_NES_N",N,"_",Sys.Date(),".pdf"))
