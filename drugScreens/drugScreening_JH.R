# analyzing drug viability results from Ava
library(plyr); library(dplyr); library(tidyr); library(ggplot2)
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/drugScreens/JH")
dataPath <- "~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/drugScreens/JH"

exp <- list.files(dataPath, ".xlsx")
meta <- data.frame(exp)
meta <- meta %>% tidyr::separate_wider_delim(exp, "_", 
                                             names=c("Drug","N"),
                                             cols_remove=FALSE)
all.df <- data.frame()
for (i in exp) {
  # read each file 
  temp.file <- readxl::read_excel(file.path(dataPath, i))
  
  # fix colnames since they were merged
  for (col in 2:ncol(temp.file)) {
    temp.colname <- colnames(temp.file)[col]
    if (!startsWith(temp.colname,"...")) {
      # restart replicate number and update cell line if appropriate
      repN <- 1
      cellLine <- temp.colname
    }
    colnames(temp.file)[col] <- paste0(cellLine,"_",repN)
    repN <- repN + 1
  }
  temp.file <- as.data.frame(temp.file)
  temp.file$concUM <- (10^(temp.file[,1]))/1000
  temp.file <- temp.file[,2:ncol(temp.file)]
  
  # reformat long
  file.long <- reshape2::melt(temp.file, id.var = "concUM", 
                              variable.name="sample")
  file.long$unit <- "Percent Relative Viability"
  
  # prep sample name format for extracting metadata
  file.long <- file.long %>% 
    tidyr::separate_wider_delim(sample, "_", names=c("MPNST", "replicate"),
                                cols_remove=FALSE)
  file.long$Drug <- strsplit(i, "_")[[1]][1]
  # concatenate data
  all.df <- rbind(all.df, file.long)
}
all.df$author <- "AW"
write.csv(all.df, paste0("AW_drugViability_data_", Sys.Date(),".csv"), row.names=FALSE)
all.df <- read.csv(paste0("AW_drugViability_data_2025-06-16.csv"))
all.df$time <- 120
all.df$time_unit <- "h"
all.df$concUnitUM <- "um"

drugPlots <- function(death.df, unit="Death") {
  # plot each drug
  if ("Drug" %in% colnames(death.df)) {
    drugs <- unique(death.df$Drug)
  } else {
    drugs <- unique(death.df$improve_drug_id)
  }
  for (d in drugs) {
    if ("DOSE" %in% colnames(death.df)) {
      drug.df <- death.df[death.df$Drug==d,]
      drug.df$DOSEum <- paste0(signif(drug.df$DOSE,1), "uM")
      drug.df$DOSEum <- factor(drug.df$DOSEum, levels=paste0(sort(unique(signif(drug.df$DOSE,1))), "uM"))
      drug.df$sample <- paste0(drug.df$DOSEum, " after ", drug.df$time,"h") 
    } else {
      drug.df <- death.df[death.df$improve_drug_id==d,]
      drug.df$sample <- paste0(drug.df$time,"h")
    }
    
    # reduce to overlapping parameters (dose, time) and filter for cell types
    gain <- drug.df[drug.df$improve_sample_id %in% c("JH-2-002","JH-2-079c"),]
    gain.params <- unique(gain$sample)
    if ("JH-2-055d" %in% drug.df$improve_sample_id & 
        nrow(drug.df[drug.df$improve_sample_id == "JH-2-055d",])>1) {
      jh <- drug.df[drug.df$improve_sample_id == "JH-2-055d",]
      jh.params <- unique(jh$sample)
      shared.params <- gain.params[gain.params %in% jh.params]
      gain <- gain[gain$sample %in% shared.params,]
      jh <- jh[jh$sample %in% shared.params,]
      
      # compare paired results
      if ("GROWTH" %in% colnames(gain)) {
        gain <- gain[order(gain$sample),]$GROWTH
        jh <- jh[order(jh$sample),]$GROWTH
      } else {
        gain <- gain[order(gain$sample),]$dose_response_value
        jh <- jh[order(jh$sample),]$dose_response_value
      }
      if (length(gain) == length(jh)) {
        paired <- TRUE
      } else {paired <- FALSE}
      if (grepl("Death",unit, ignore.case=TRUE)) {
        if (median(gain) < median(jh)) {
          drug.p <- t.test(jh, gain, alternative="greater", paired=paired)$p.value 
          if (paired) {
            title <- paste0("JH-2-055d (diploid) is more sensitive to ",d," (paired p=", signif(drug.p,2),")")
          } else {
            title <- paste0("JH-2-055d (diploid) is more sensitive to ",d," (p=", signif(drug.p,2),")") 
          }
        } else {
          drug.p <- t.test(gain, jh, alternative="greater", paired=paired)$p.value 
          if (paired) {
            title <- paste0("Cell lines with chr8q gain are more sensitive to ",d," (paired p=", signif(drug.p,2),")")
          } else {
            title <- paste0("Cell lines with chr8q gain are more sensitive to ",d," (p=", signif(drug.p,2),")") 
          }
        } 
      } else {
        if (median(gain) > median(jh)) {
          drug.p <- t.test(gain, jh, alternative="greater", paired=paired)$p.value 
          if (paired) {
            title <- paste0("JH-2-055d (diploid) is more sensitive to ",d," (paired p=", signif(drug.p,2),")")
          } else {
            title <- paste0("JH-2-055d (diploid) is more sensitive to ",d," (p=", signif(drug.p,2),")") 
          }
        } else {
          drug.p <- t.test(jh, gain, alternative="greater", paired=paired)$p.value 
          if (paired) {
            title <- paste0("Cell lines with chr8q gain are more sensitive to ",d," (paired p=", signif(drug.p,2),")")
          } else {
            title <- paste0("Cell lines with chr8q gain are more sensitive to ",d," (p=", signif(drug.p,2),")") 
          }
        }
      } 
    } else {
      title <- paste0("Cell lines treated with ",d)
    }
    if ("GROWTH" %in% colnames(drug.df)) {
      ggplot2::ggplot(drug.df, aes(x=time, y=GROWTH, color=improve_sample_id)) + 
        geom_point() + facet_grid(.~DOSEum) + ggtitle(title) + theme_classic() +
        labs(y=paste("Normalized",unit), x="Time (h)", color="MPNST") + 
        theme(plot.title=element_text(hjust=0.5))
      ggsave(paste0(tolower(unit),"_",d,".pdf"),width=12, height=6) 
    } else {
      ggplot2::ggplot(drug.df, aes(x=time, y=dose_response_value, color=improve_sample_id)) + 
        geom_point() + ggtitle(title) + theme_classic() +
        labs(y=paste("Normalized",unit), x="Time (h)", color="MPNST") + 
        theme(plot.title=element_text(hjust=0.5))
      ggsave(paste0(tolower(unit),"_",d,".pdf"),width=12, height=6)
    }
  }
}

# get curve fitting script
if (!file.exists("fit_curve.py")) {download.file('https://raw.githubusercontent.com/PNNL-CompBio/coderdata/refs/heads/main/build/utils/fit_curve.py',
              destfile="fit_curve.py")}
maxConc <- c(max(all.df$concUM)) # uM; all non-control concentrations: 0.0001000000 0.0002999991 0.0010000000 0.0029999982 0.0100000000 0.0299999824 0.1000000000 0.2999998240 1.0000000000
for (m in maxConc) {
  # change colnames to use existing curve fitting script
  confluence.df <- all.df[all.df$concUM <= m,] 
  confluence.df$study <- "Chr8"
  #oldCols <- c("concentration","value", "study", "author", "MPNST", "Drug", "Elapsed", "time_unit")
  oldCols <- c("concUM","value", "study", "author", "MPNST", "Drug", "time", "time_unit")
  newCols <- c('DOSE','GROWTH','study','source','improve_sample_id','Drug','time','time_unit')
  confluence.df <- confluence.df %>% rename_at(vars(oldCols), ~ newCols)
  write.table(confluence.df[,newCols], paste0("chr8_normConfluence_max",m,"um.tsv"), sep="\t", row.names=FALSE)
  
  # plot each drug
  drugPlots(confluence.df, "Confluence")
  
  # run curve fitting
  library(reticulate)
  #setwd("~/")
  system2("python3", paste0("fit_curve.py --input chr8_normConfluence_max",m,
                            "um.tsv --output chr8_normConfluenceCurves_max",m,"um"))
  # if you get an error "No module named 'pandas'":
  # open the terminal and run the following before trying again:
  # python3 -m venv env
  # source env/bin/activate
  # pip3 install pandas
  
  # from there, can just run it in the terminal: python3 fit_curve.py --input chr8_normConfluence_max1um.tsv --output chr8_normConfluenceCurves_max1um
  
  confluence <- read.table(paste0("chr8_normConfluenceCurves_max",m,"um.0"), header=TRUE, fill=NA)
  ## plot
  # AUC
  auc <- confluence[confluence$dose_response_metric=="fit_auc",]
  drugPlots(auc, "Confluence Area Under the Curve")
  mean.auc <- plyr::ddply(auc, .(improve_drug_id), summarize,
                          auc=mean(dose_response_value, na.rm=TRUE))
  
  # r2
  r2 <- confluence[confluence$dose_response_metric=="fit_r2",]
  
  # AUC with r2 fill
  auc.r2 <- merge(auc, r2, by=c("improve_drug_id","improve_sample_id", "time"))
  auc.r2[auc.r2$improve_sample_id=="JH-2-055d",]$improve_sample_id <- "JH-2-055d (diploid)"
  auc.r2[auc.r2$improve_sample_id=="JH-2-002",]$improve_sample_id <- "JH-2-002 (gain)"
  auc.r2[auc.r2$improve_sample_id=="JH-2-079c",]$improve_sample_id <- "JH-2-079c (gain)"
  drugOrder <- unique(mean.auc[order(mean.auc$auc),]$improve_drug_id)
  ggplot2::ggplot(auc.r2, 
                  aes(x=time, y=dose_response_value.x, alpha=dose_response_value.y, color=improve_sample_id)) + 
    geom_point() + theme_classic(base_size=12) + facet_grid(.~improve_drug_id)+
    labs(y="AUC: Normalized Confluence", x= "Time (h)", alpha=expression(paste("Fit R"^"2")), color="MPNST") + 
    theme(#axis.title.x=element_blank(), 
      axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0("fit_auc_r2_confluence_barPlot_max",m,"um.pdf"), width=10, height=3)
  
  quantiles <- quantile(r2$dose_response_value)
  auc.r2$Quantile <- NA
  for (q in names(quantiles)) {
    if (any(auc.r2$dose_response_value.y >= quantiles[[q]])){
      auc.r2[auc.r2$dose_response_value.y >= quantiles[[q]],]$Quantile <- q 
    }
  }
  ggplot2::ggplot(auc.r2, 
                  aes(x=time, y=dose_response_value.x, shape=Quantile, color=improve_sample_id)) + 
    geom_point() + theme_classic(base_size=12) + facet_grid(.~improve_drug_id)+
    labs(y="AUC: Normalized Confluence", x= "Time (h)", shape=expression(paste("R"^"2","Quantile")), color="MPNST") + 
    theme(#axis.title.x=element_blank(), 
      axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + 
    scale_shape_manual(values=seq(1,5),labels=paste(c("0%", "25%", "50%", "75%", "100%"),quantiles), 
                                breaks=c("0%", "25%", "50%", "75%", "100%"))
  ggsave(paste0("fit_auc_r2_confluence_barPlot_max",m,"um_v2.pdf"), width=10, height=3)
  
  ggplot2::ggplot(auc.r2, 
                  aes(x=time, y=dose_response_value.x, color=improve_sample_id)) + 
    geom_point() + theme_classic(base_size=12) + facet_grid(.~improve_drug_id)+
    labs(y="AUC: Normalized Confluence", x= "Time (h)", color="MPNST") + 
    theme(#axis.title.x=element_blank(), 
      axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0("fit_auc_confluence_barPlot_max",m,"um.pdf"), width=10, height=3)
}

#### association between AUC and chr8/chr8q/MYC ####
med.chr8q <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/proteomics/analysis/Chr8_quant_20250409/positional_medians/Copy Number/Copy Number_Chr8q_median.csv")
auc[auc$improve_sample_id=="JH-2-079",]$improve_sample_id <- "JH-2-079c" # confirm with Ava - discrepancy is with BI-2536
auc.chr8q <- merge(auc, med.chr8q, by.x="improve_sample_id",by.y="Sample") # only 2 samples: JH-2-002 and JH-2-079c
# does AUC correlate with chr8q copy number? MYC copy number? - will have to skip due to lack of copy number data; plus the copy number is PDX vs. cell lines screened with drugs
# actually, have copy number from FISH analysis
fish <- read.csv("FISH_results.csv")
colnames(fish)[1] <- "improve_sample_id"
auc.fish <- merge(auc, fish, by="improve_sample_id")
Drug <- unique(auc$improve_drug_id) # 5
t.df <- data.frame(Drug=rep(Drug,2), Pearson.est=NA, Pearson.p=NA, 
                   Spearman.est=NA, Spearman.p=NA, 
                   test=rep(c("MYC","chr8"),length(Drug)))

for (d in Drug) {
  temp.auc.fish <- auc.fish[auc.fish$improve_drug_id==d,]
  for (t in c("chr8","MYC")) {
    temp.auc.fish$rank <- as.numeric(temp.auc.fish[,t])
    temp.auc.fish <- na.omit(temp.auc.fish) # remove NF11.1 because not sure exact copy number or rank
    temp.test <- cor.test(temp.auc.fish$rank,temp.auc.fish$dose_response_value)
    temp.test.sp <- cor.test(temp.auc.fish$rank,temp.auc.fish$dose_response_value, method="spearman")
    
    t.df[t.df$Drug == d & t.df$test==t,]$Pearson.est <- temp.test$estimate
    t.df[t.df$Drug == d & t.df$test==t,]$Pearson.p <- temp.test$p.value
    t.df[t.df$Drug == d & t.df$test==t,]$Spearman.est <- temp.test.sp$estimate
    t.df[t.df$Drug == d & t.df$test==t,]$Spearman.p <- temp.test.sp$p.value
    
    if (is.numeric(temp.test$p.value) | is.numeric(temp.test.sp$p.value)) { # was checking for p < 0.05 but now want plots either way
      stats_pearson <- substitute(
        r == est * "," ~ ~"p" ~ "=" ~ p,
        list(
          est = format(as.numeric(temp.test$estimate), digits = 3),
          p = format(temp.test$p.value, digits = 3)
        )
      )
      ggplot(temp.auc.fish, aes(x=rank, y=dose_response_value)) + geom_point() + theme_minimal() + 
        labs(x=paste(t,"Copy Number"), y = paste(d,"Sensitivity (AUC)")) +
        geom_smooth(se=FALSE, linetype="dashed", method="lm") + 
        ggrepel::geom_label_repel(aes(label=improve_sample_id)) + 
        ggplot2::geom_text(
          x = Inf, y = -Inf, vjust = "inward", hjust = "inward",
          colour = "blue", parse = TRUE,
          label = as.character(as.expression(stats_pearson)), size = 5
        ) 
      ggsave(paste0(d,"_",t,"_scatterPlot_PearsonCorrelation.pdf"), width=5, height=5)
      
      stats_spearman <- substitute(
        rho == est * "," ~ ~"p" ~ "=" ~ p,
        list(
          est = format(as.numeric(temp.test.sp$estimate), digits = 3),
          p = format(temp.test.sp$p.value, digits = 3)
        )
      )
      ggplot(temp.auc.fish, aes(x=rank, y=dose_response_value)) + geom_point() + theme_minimal() + 
        labs(x=paste(t,"Copy Number"), y = paste(d,"Sensitivity (AUC)")) +
        geom_smooth(se=FALSE, linetype="dashed", method="lm") + 
        ggrepel::geom_label_repel(aes(label=improve_sample_id)) + 
        ggplot2::geom_text(
          x = Inf, y = -Inf, vjust = "inward", hjust = "inward",
          colour = "blue", parse = TRUE,
          label = as.character(as.expression(stats_spearman)), size = 5
        ) 
      ggsave(paste0(d,"_",t,"_scatterPlot_SpearmanCorrelation.pdf"), width=5, height=5)
    }
  }
}
write.csv(t.df, "correlations_AUC_chr8_MYC_2025-07-28.csv", row.names=FALSE) # none significant

# try plotting PLKi and EGFRi together and facet by chr8/MYC copy number
auc.fish$moa <- NA
auc.fish[auc.fish$improve_drug_id %in% c("BI-2536","Volasertib"),]$moa <- "PLK Inhibitor"
auc.fish[auc.fish$improve_drug_id %in% c("Erlotinib","Osimertib"),]$moa <- "EGFR Inhibitor"
auc.fish[auc.fish$improve_drug_id %in% c("Mirdametinib"),]$moa <- "MEK Inhibitor"
moa <- na.omit(unique(auc.fish$moa)) # 5
t.df <- data.frame(moa=rep(moa,2), Pearson.est=NA, Pearson.p=NA, 
                   Spearman.est=NA, Spearman.p=NA, 
                   test=rep(c("MYC","chr8"),length(moa)))
for (m in moa) {
  temp.auc.fish <- auc.fish[auc.fish$moa==m,]
  for (t in c("chr8","MYC")) {
    temp.auc.fish$rank <- as.numeric(temp.auc.fish[,t])
    temp.auc.fish <- na.omit(temp.auc.fish) # remove NF11.1 because not sure exact copy number or rank
    temp.test <- cor.test(temp.auc.fish$rank,temp.auc.fish$dose_response_value)
    temp.test.sp <- cor.test(temp.auc.fish$rank,temp.auc.fish$dose_response_value, method="spearman")
    
    t.df[t.df$moa == m & t.df$test==t,]$Pearson.est <- temp.test$estimate
    t.df[t.df$moa == m & t.df$test==t,]$Pearson.p <- temp.test$p.value
    t.df[t.df$moa == m & t.df$test==t,]$Spearman.est <- temp.test.sp$estimate
    t.df[t.df$moa == m & t.df$test==t,]$Spearman.p <- temp.test.sp$p.value
    
    if (is.numeric(temp.test$p.value) | is.numeric(temp.test.sp$p.value)) { # was checking for p < 0.05 but now want plots either way
      stats_pearson <- substitute(
        r == est * "," ~ ~"p" ~ "=" ~ p,
        list(
          est = format(as.numeric(temp.test$estimate), digits = 3),
          p = format(temp.test$p.value, digits = 3)
        )
      )
      ggplot(temp.auc.fish, aes(x=rank, y=dose_response_value)) + geom_point() + theme_minimal() + 
        labs(x=paste(t,"Copy Number"), y = paste(m,"Sensitivity (AUC)")) +
        geom_smooth(se=FALSE, linetype="dashed", method="lm") + 
        ggrepel::geom_label_repel(aes(label=improve_sample_id)) + 
        ggplot2::geom_text(
          x = Inf, y = -Inf, vjust = "inward", hjust = "inward",
          colour = "blue", parse = TRUE,
          label = as.character(as.expression(stats_pearson)), size = 5
        ) 
      ggsave(paste0(m,"_",t,"_scatterPlot_PearsonCorrelation.pdf"), width=5, height=5)
      
      stats_spearman <- substitute(
        rho == est * "," ~ ~"p" ~ "=" ~ p,
        list(
          est = format(as.numeric(temp.test.sp$estimate), digits = 3),
          p = format(temp.test.sp$p.value, digits = 3)
        )
      )
      ggplot(temp.auc.fish, aes(x=rank, y=dose_response_value)) + geom_point() + theme_minimal() + 
        labs(x=paste(t,"Copy Number"), y = paste(m,"Sensitivity (AUC)")) +
        geom_smooth(se=FALSE, linetype="dashed", method="lm") + 
        ggrepel::geom_label_repel(aes(label=improve_sample_id)) + 
        ggplot2::geom_text(
          x = Inf, y = -Inf, vjust = "inward", hjust = "inward",
          colour = "blue", parse = TRUE,
          label = as.character(as.expression(stats_spearman)), size = 5
        ) 
      ggsave(paste0(m,"_",t,"_scatterPlot_SpearmanCorrelation.pdf"), width=5, height=5)
    }
  }
}
write.csv(t.df, "moa_correlations_AUC_chr8_MYC_2025-07-28.csv", row.names=FALSE) # none significant

# is PLKi AUC lower with MYC gain? is EGFRi AUC lower with MYC diploid? Likewise with chr8
auc$chr8 <- "Gain"
auc[auc$improve_sample_id %in% c("JH-2-055d","ST8814","NF10.1"),]$chr8 <- "Diploid"
auc$chr8q <- "Gain"
auc[auc$improve_sample_id %in% c("JH-2-055d","ST8814","NF10.1"),]$chr8q <- "Diploid"
auc$MYC <- "Gain"
auc[auc$improve_sample_id == "JH-2-055d",]$MYC <- "Diploid"
Drug <- unique(auc$improve_drug_id) # 5
t.df <- data.frame(Drug=rep(Drug,2), p=NA, mean_gain = NA, median_gain=NA, sd_gain=NA, N_gain=NA,
                   mean_diploid = NA, median_diploid=NA, sd_diploid=NA, N_diploid=NA, test=rep(c("MYC","chr8"),length(Drug)),
                   gain_moreSensitive=FALSE)

for (d in Drug) {
  for (t in c("chr8","MYC")) {
    gain <- na.omit(auc[auc[,t]=="Gain" & auc$improve_drug_id==d,]$dose_response_value)
    diploid <- na.omit(auc[auc[,t]=="Diploid" & auc$improve_drug_id==d,]$dose_response_value)
    if (length(diploid) > 1) {
      if (mean(gain) > mean(diploid)) {
        temp.test <- t.test(gain, diploid, "greater")
      } else {
        temp.test <- t.test(gain, diploid, "less")
        t.df[t.df$Drug == d & t.df$test==t,]$gain_moreSensitive <- TRUE
      }
      t.df[t.df$Drug == d & t.df$test==t,]$p <- temp.test$p.value
      t.df[t.df$Drug == d & t.df$test==t,]$mean_gain <- mean(gain) 
      t.df[t.df$Drug == d & t.df$test==t,]$median_gain <- median(gain)
      t.df[t.df$Drug == d & t.df$test==t,]$sd_gain <- sd(gain)
      t.df[t.df$Drug == d & t.df$test==t,]$N_gain <- length(gain)
      
      t.df[t.df$Drug == d & t.df$test==t,]$mean_diploid <- mean(diploid) 
      t.df[t.df$Drug == d & t.df$test==t,]$median_diploid <- median(diploid)
      t.df[t.df$Drug == d & t.df$test==t,]$sd_diploid <- sd(diploid)
      t.df[t.df$Drug == d & t.df$test==t,]$N_diploid <- length(diploid)
    }
  }
}
write.csv(t.df, "tTests_AUC_chr8_MYC.csv", row.names=FALSE)
# lowest p-value is 0.0874 (all for chr8; can't do MYC because only 1 MYC diploid)

# repeat with moa instead of drug
auc$moa <- NA
auc[auc$improve_drug_id %in% c("BI-2536","Volasertib"),]$moa <- "PLK Inhibitor"
auc[auc$improve_drug_id %in% c("Erlotinib","Osimertinib"),]$moa <- "EGFR Inhibitor"
auc[auc$improve_drug_id %in% c("Mirdametinib"),]$moa <- "MEK Inhibitor"
moa <- na.omit(unique(auc$moa)) # 5
t.df <- data.frame(moa=rep(moa,2), p=NA, mean_gain = NA, median_gain=NA, sd_gain=NA, N_gain=NA,
                   mean_diploid = NA, median_diploid=NA, sd_diploid=NA, N_diploid=NA, test=rep(c("MYC","chr8"),length(moa)),
                   gain_moreSensitive=FALSE)

for (m in moa) {
  for (t in c("chr8","MYC")) {
    gain <- na.omit(auc[auc[,t]=="Gain" & auc$moa==m,]$dose_response_value)
    diploid <- na.omit(auc[auc[,t]=="Diploid" & auc$moa==m,]$dose_response_value)
    if (length(diploid) > 1) {
      if (mean(gain) > mean(diploid)) {
        temp.test <- t.test(gain, diploid, "greater")
      } else {
        temp.test <- t.test(gain, diploid, "less")
        t.df[t.df$moa==m & t.df$test==t,]$gain_moreSensitive <- TRUE
      }
      t.df[t.df$moa==m & t.df$test==t,]$p <- temp.test$p.value
      t.df[t.df$moa==m & t.df$test==t,]$mean_gain <- mean(gain) 
      t.df[t.df$moa==m & t.df$test==t,]$median_gain <- median(gain)
      t.df[t.df$moa==m & t.df$test==t,]$sd_gain <- sd(gain)
      t.df[t.df$moa==m & t.df$test==t,]$N_gain <- length(gain)
      
      t.df[t.df$moa==m & t.df$test==t,]$mean_diploid <- mean(diploid) 
      t.df[t.df$moa==m & t.df$test==t,]$median_diploid <- median(diploid)
      t.df[t.df$moa==m & t.df$test==t,]$sd_diploid <- sd(diploid)
      t.df[t.df$moa==m & t.df$test==t,]$N_diploid <- length(diploid)
    }
  }
}
write.csv(t.df, "moa_tTests_AUC_chr8_MYC.csv", row.names=FALSE)

temp.df <- reshape2::melt(auc, id.var=c("source", "improve_sample_id", 
                                        "improve_drug_id", "study", "time",
                                        "time_unit", "dose_response_metric",
                                        "dose_response_value","moa"))
temp.df <- temp.df[temp.df$variable != "chr8q" & temp.df$moa != "MEK Inhibitor",]

all.mpnst <- c("JH-2-055d","ST8814","NF10.1","JH-2-002","JH-2-079c","JH-2-103","JH-2-155b","NF11.1","NF90.8","S462")
ggplot(temp.df, aes(x=value, y=dose_response_value)) + geom_violin(alpha=0) + 
  geom_boxplot() + geom_point(aes(color=improve_drug_id, shape=improve_sample_id)) + 
  facet_grid(moa~variable) + theme_classic() + labs(color="Drug", shape="MPNST Cell Line", y ="Viability Area Under the Curve") +
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst) + theme(axis.title.x=element_blank())
ggsave("AUC_vs_moa_boxPlot.pdf",width=4, height=4.5)

temp.df$moa <- factor(temp.df$moa, levels=c("PLK Inhibitor","EGFR Inhibitor"))
ggplot(temp.df, aes(x=value, y=dose_response_value)) + geom_violin(alpha=0) + 
  geom_boxplot() + geom_point(aes(color=improve_drug_id, shape=improve_sample_id)) + 
  facet_grid(variable~moa, switch="y") + theme_classic() + 
  labs(color="Drug", shape="MPNST Cell Line", y ="Viability Area Under the Curve") +
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst) + 
  theme(axis.title.y=element_blank(), legend.position="bottom") + coord_flip() + 
  scale_x_discrete(position="top")
ggsave("AUC_vs_moa_boxPlot_horizontal.pdf",width=4.5, height=3)
ggsave("AUC_vs_moa_boxPlot_horizontal_wider.pdf",width=6, height=3)

temp.df$Drug <- factor(temp.df$improve_drug_id, levels=c("Erlotinib","Osimertinib","BI-2536","Volasertib"))
temp.df$Drug <- factor(temp.df$improve_drug_id, levels=c("BI-2536","Volasertib","Erlotinib","Osimertinib"))
ggplot(temp.df, aes(x=value, y=dose_response_value)) + geom_violin(alpha=0) + 
  geom_boxplot() + geom_point(aes(color=moa, shape=improve_sample_id)) + 
  facet_grid(Drug~variable) + theme_classic() + labs(color="MOA", shape="MPNST Cell Line", y ="Viability Area Under the Curve") +
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst) + theme(axis.title.x=element_blank())
ggsave("AUC_vs_drug_boxPlot_v2.pdf",width=4, height=6)
ggsave("AUC_vs_drug_boxPlot_v2_shorter.pdf",width=4, height=5)

temp.df$Drug <- factor(temp.df$improve_drug_id, levels=c("BI-2536","Volasertib","Erlotinib","Osimertinib"))
ggplot(temp.df, aes(x=value, y=dose_response_value)) + geom_violin(alpha=0) + 
  geom_boxplot() + geom_point(aes(color=moa, shape=improve_sample_id)) + 
  facet_grid(variable~Drug, switch="y") + theme_classic() + 
  labs(color="MOA", shape="MPNST Cell Line", y ="Viability Area Under the Curve") +
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst) + 
  theme(axis.title.y=element_blank(), legend.position="bottom") + coord_flip() + 
  scale_x_discrete(position="top")
ggsave("AUC_vs_drug_boxPlot_horizontal.pdf",width=6, height=3)
ggsave("AUC_vs_drug_boxPlot_horizontal_wider.pdf",width=8, height=3)
write.csv(temp.df,"AUC_vs_drug_or_moa.csv", row.names=FALSE)
#### manually checking for significant correlations between MYC copy number and PLK1 or EGFR RNAi ####
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure4_viability/depmap")
##### PLK1 #####
soft <- read.csv("MYC Copy Number log2(relative to ploidy + 1) vs PLK1 Gene Effect (DEMETER2) filtered by Soft Tissue.csv")
soft.corr <- cor.test(soft$MYC.Copy.Number.log2.relative.to.ploidy...1..Copy.Number.Public.25Q2..Log2.transformed.,
                      soft$PLK1.Gene.Effect..DEMETER2..RNAi..Achilles.DRIVE.Marcotte..DEMETER2., method="spearman")
# p = 0.2912 spearman, 0.6095 pearson

#all <- read.csv("MYC Copy Number log2(relative to ploidy + 1) vs PLK1 Gene Effect (DEMETER2).csv")


##### EGFR #####
soft <- read.csv("MYC Copy Number log2(relative to ploidy + 1) vs EGFR Gene Effect (DEMETER2) filtered by Soft Tissue.csv")
soft.corr <- cor.test(soft$MYC.Copy.Number.log2.relative.to.ploidy...1..Copy.Number.Public.25Q2..Log2.transformed.,
                      soft$EGFR.Gene.Effect..DEMETER2..RNAi..Achilles.DRIVE.Marcotte..DEMETER2.
                      #, method="spearman"
                      )
soft.corr
# p = 0.4366 spearman, 0.8074709 pearson