# analyzing drug viability results from Ava
library(plyr); library(dplyr); library(tidyr); library(ggplot2)
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/drugScreens/JH")
dataPath <- "~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/drugScreens/JH"

#### compile data ####
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
#all.df <- read.csv(paste0("AW_drugViability_data_2025-06-16.csv"))
all.df <- read.csv(paste0("AW_drugViability_data_2025-08-11.csv"))
all.df$time <- 120
all.df$time_unit <- "h"
all.df$concUnitUM <- "um"

#### fit curves ####
# get curve fitting script
download.file('https://raw.githubusercontent.com/PNNL-CompBio/coderdata/refs/heads/main/build/utils/fit_curve.py',
              destfile="fit_curve.py")
m <- c(max(all.df$concUM)) # uM; all non-control concentrations: 0.0001000000 0.0002999991 0.0010000000 0.0029999982 0.0100000000 0.0299999824 0.1000000000 0.2999998240 1.0000000000

# change colnames to use existing curve fitting script
confluence.df <- all.df[all.df$concUM <= m,] 
confluence.df$study <- "Chr8"
#oldCols <- c("concentration","value", "study", "author", "MPNST", "Drug", "Elapsed", "time_unit")
oldCols <- c("concUM","value", "study", "author", "MPNST", "Drug", "time", "time_unit")
newCols <- c('DOSE','GROWTH','study','source','improve_sample_id','Drug','time','time_unit')
confluence.df <- confluence.df %>% rename_at(vars(oldCols), ~ newCols)
write.table(confluence.df[,newCols], paste0("chr8_normConfluence_max",m,"um.tsv"), sep="\t", row.names=FALSE)

# plot each drug
#drugPlots(confluence.df, "Confluence")

# run curve fitting
library(reticulate)
#setwd("~/")
#reticulate::py_install("pandas")
#system2("brew","install","pandas")
system2("python3", paste0("fit_curve.py --input chr8_normConfluence_max",m,
                          "um.tsv --output chr8_normConfluenceCurves_max",m,"um"))
# if you get an error "No module named 'pandas'":
# open the terminal and run the following before trying again:
# python3 -m venv env
# source env/bin/activate
# pip3 install pandas

# from there, can just run it in the terminal: python3 fit_curve.py --input chr8_normConfluence_max10um.tsv --output chr8_normConfluenceCurves_max10um
# make sure to cd first

#### association between AUC and chr8/chr8q/MYC ####
# does AUC correlate with chr8q copy number? MYC copy number? from FISH analysis
confluence <- read.table(paste0("chr8_normConfluenceCurves_max",m,"um.0"), header=TRUE, fill=NA)

# AUC
auc <- confluence[confluence$dose_response_metric=="auc",] # was fit AUC before 20250811

# FISH
fish <- read.csv("FISH/FISH_results.csv")
colnames(fish)[1] <- "improve_sample_id"
fish[grepl("79c",fish$improve_sample_id),]$improve_sample_id <- "JH-2-079c"

# merge
auc.fish <- merge(auc, fish, by="improve_sample_id")
Drug <- unique(auc$improve_drug_id) # 5

##### correlations for each drug #####
t.df <- data.frame(Drug=sort(rep(Drug,2)), Pearson.est=NA, Pearson.p=NA, 
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

##### correlations for each moa #####
auc.fish$moa <- NA
auc.fish[auc.fish$improve_drug_id %in% c("BI-2536","Volasertib"),]$moa <- "PLK Inhibitor"
auc.fish[auc.fish$improve_drug_id %in% c("Erlotinib","Osimertib"),]$moa <- "EGFR Inhibitor"
moa <- na.omit(unique(auc.fish$moa))
t.df <- data.frame(moa=sort(rep(moa,2)), Pearson.est=NA, Pearson.p=NA, 
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

##### t-tests for each drug #####
# is PLKi AUC lower with MYC gain? is EGFRi AUC lower with MYC diploid? Likewise with chr8
auc$chr8 <- "Gain"
auc[auc$improve_sample_id %in% c("JH-2-055d","ST8814","NF10.1"),]$chr8 <- "Diploid"
auc$chr8q <- "Gain"
auc[auc$improve_sample_id %in% c("JH-2-055d","ST8814","NF10.1"),]$chr8q <- "Diploid"
auc$MYC <- "Gain"
auc[auc$improve_sample_id == "JH-2-055d",]$MYC <- "Diploid"
Drug <- unique(auc$improve_drug_id) # 5
t.df <- data.frame(Drug=sort(rep(Drug,2)), p=NA, mean_gain = NA, median_gain=NA, sd_gain=NA, N_gain=NA,
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

##### t-tests for each moa #####
auc$moa <- NA
auc[auc$improve_drug_id %in% c("BI-2536","Volasertib"),]$moa <- "PLK Inhibitor"
auc[auc$improve_drug_id %in% c("Erlotinib","Osimertinib"),]$moa <- "EGFR Inhibitor"
moa <- na.omit(unique(auc$moa))
t.df <- data.frame(moa=sort(rep(moa,2)), p=NA, mean_gain = NA, median_gain=NA, sd_gain=NA, N_gain=NA,
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

#### plot AUCs ####
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

#### plot dose-response curves ####
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/drugScreens")
chr8_tTest <- function(drug.df, unit, factors=c("chr8q","MYC")) {
  d <- unique(drug.df$Drug)
  if (length(d) > 1) {
    stop("more than 1 drug input")
  }
  titles <- list()
  for (f in factors) {
    if (any(drug.df[,f] == "Diploid") & any(drug.df[,f] == "Gain")) {
      diploid <- drug.df[drug.df[,f] == "Diploid",]
      gain <- drug.df[drug.df[,f] == "Gain",]
      # reduce to overlapping parameters (dose, time) and filter for cell types
      gain.params <- unique(gain$sample)
      diploid.params <- unique(diploid$sample)
      shared.params <- gain.params[gain.params %in% diploid.params]
      gain <- gain[gain$sample %in% shared.params,]
      diploid <- diploid[diploid$sample %in% shared.params,]
      
      # compare paired results
      if ("GROWTH" %in% colnames(gain)) {
        gain <- gain[order(gain$sample),]$GROWTH
        diploid <- diploid[order(diploid$sample),]$GROWTH
      } else {
        gain <- gain[order(gain$sample),]$dose_response_value
        diploid <- diploid[order(diploid$sample),]$dose_response_value
      }
      if (length(gain) == length(diploid)) {
        paired <- TRUE
      } else {paired <- FALSE}
      if (grepl("Death",unit, ignore.case=TRUE)) {
        if (median(gain) < median(diploid)) {
          drug.p <- t.test(diploid, gain, alternative="greater", paired=paired)$p.value 
          if (paired) {
            titles[[f]] <- paste0(f," diploid is more sensitive to ",d," (paired p=", signif(drug.p,2),")")
          } else {
            titles[[f]] <- paste0(f," diploid is more sensitive to ",d," (p=", signif(drug.p,2),")") 
          }
        } else {
          drug.p <- t.test(gain, diploid, alternative="greater", paired=paired)$p.value 
          if (paired) {
            titles[[f]] <- paste0(f," gain is more sensitive to ",d," (paired p=", signif(drug.p,2),")")
          } else {
            titles[[f]] <- paste0(f, " gain is more sensitive to ",d," (p=", signif(drug.p,2),")") 
          }
        } 
      } else {
        if (median(gain) > median(diploid)) {
          drug.p <- t.test(gain, diploid, alternative="greater", paired=paired)$p.value 
          if (paired) {
            titles[[f]] <- paste0(f," diploid is more sensitive to ",d," (paired p=", signif(drug.p,2),")")
          } else {
            titles[[f]] <- paste0(f," diploid is more sensitive to ",d," (p=", signif(drug.p,2),")") 
          }
        } else {
          drug.p <- t.test(diploid, gain, alternative="greater", paired=paired)$p.value 
          if (paired) {
            titles[[f]] <- paste0(f," gain is more sensitive to ",d," (paired p=", signif(drug.p,2),")")
          } else {
            titles[[f]] <- paste0(f," gain is more sensitive to ",d," (p=", signif(drug.p,2),")") 
          }
        }
      } 
    } else {
      titles[[f]] <- ""
    } 
  }
  return(titles)
}

# get data
rel.conf <- read.table("JH/chr8_normConfluence_max10um.tsv", sep="\t", header=TRUE)

# calculate mean and sd for each drug, dose, mpnst, time combo; also preserve sample and chr8q columns
library(plyr);library(dplyr)
mean.conf <- plyr::ddply(rel.conf, .(improve_sample_id, time, Drug, DOSE), summarize,
                         meanGROWTH = mean(GROWTH, na.rm=TRUE),
                         sdGROWTH = sd(GROWTH, na.rm=TRUE))

# annotate MYC and chr8q status
known.dip <- c("JH-2-055d","ST8814","NF10.1")
known.amp <- unique(mean.conf$improve_sample_id[!(mean.conf$improve_sample_id %in% known.dip)])
mean.conf$chr8q <- "Not yet known"
mean.conf[mean.conf$improve_sample_id %in% known.amp,]$chr8q <- "Gain"
mean.conf[mean.conf$improve_sample_id %in% known.dip,]$chr8q <- "Diploid"
mean.conf$chr8qv2 <- mean.conf$chr8q
mean.conf[mean.conf$improve_sample_id %in% c("ST8814","NF10.1"),]$chr8qv2 <- "Diploid with MYC gain"
mean.conf$MYC <- "Gain"
mean.conf[mean.conf$improve_sample_id == 'JH-2-055d',]$MYC <- "Diploid"

rel.conf$chr8q <- "Not yet known"
rel.conf[rel.conf$improve_sample_id %in% known.amp,]$chr8q <- "Gain"
rel.conf[rel.conf$improve_sample_id %in% known.dip,]$chr8q <- "Diploid"
rel.conf$MYC <- "Gain"
rel.conf[rel.conf$improve_sample_id == 'JH-2-055d',]$MYC <- "Diploid"
mean.conf$knownChr8q <- FALSE
mean.conf[mean.conf$improve_sample_id %in% c(known.dip,known.amp),]$knownChr8q <- TRUE
rel.conf$knownChr8q <- FALSE
rel.conf[rel.conf$improve_sample_id %in% c(known.dip,known.amp),]$knownChr8q <- TRUE
mean.conf$alpha <- 0.5
mean.conf[mean.conf$improve_sample_id %in% c(known.dip,known.amp),]$alpha <- 1

# add sample column with all parameters
mean.conf$sample <- paste0(mean.conf$Drug,
                           "_",mean.conf$DOSE,"uM")
rel.conf$sample <- paste0(rel.conf$Drug,
                          "_",rel.conf$DOSE,"uM")

all.mpnst <- unique(mean.conf$improve_sample_id)
all.mpnst <- c(known.dip,known.amp, all.mpnst[!(all.mpnst %in% c(known.dip,known.amp))])
dir.create("singleTimePlots_20250811")
setwd("singleTimePlots_20250811")
library(drc) # curve source: answer by greenjune: https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc
# Sara shared these 2 links: https://stackoverflow.com/questions/68209998/plot-drc-using-ggplot; https://forum.posit.co/t/extract-out-points-in-dose-response-curve-produced-by-drc/159433
Metric <- c("% Relative Viability")
library(ggplot2)
all.p.df <- data.frame()
drugs <- unique(mean.conf$Drug)
mpnst <- unique(mean.conf$improve_sample_id)

# plot relative confluence and relative death for each drug
for (d in drugs) {
  # filter for drug
  drug.conf <- rel.conf[rel.conf$Drug == d,]
  mean.drug.conf <- mean.conf[mean.conf$Drug == d,]
  
  p.df <- data.frame(Metric, chr8_moreSensitive = NA, MYC_moreSensitive = NA,
                     chr8_p = NA, MYC_p=NA, Drug=d, Time=t, maxConcUM=c, minValue=NA)
  # run t-tests between gain and diploid
  titles <- chr8_tTest(na.omit(drug.conf), "% Relative Viability")
  conf.more.sens <- strsplit(titles$chr8, " ")[[1]][2]
  conf.p <- strsplit(titles$chr8, "p=")[[1]][2]
  conf.p <- as.numeric(substr(conf.p,1,nchar(conf.p)-1))
  
  conf.myc <- strsplit(titles$MYC, " ")[[1]][2]
  conf.p.myc <- strsplit(titles$MYC, "p=")[[1]][2]
  conf.p.myc <- as.numeric(substr(conf.p.myc,1,nchar(conf.p.myc)-1))
  
  p.df[p.df$Metric=="% Relative Viability",]$chr8_moreSensitive <- conf.more.sens
  p.df[p.df$Metric=="% Relative Viability",]$MYC_moreSensitive <- conf.myc
  p.df[p.df$Metric=="% Relative Viability",]$chr8_p <- conf.p
  p.df[p.df$Metric=="% Relative Viability",]$MYC_p <- conf.p.myc
  p.df[p.df$Metric=="% Relative Viability",]$minValue <- min(drug.conf[drug.conf$DOSE <= c,]$GROWTH, na.rm=TRUE)
  all.p.df <- rbind(all.p.df, p.df)
  
  conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$improve_sample_id %in% c("JH-2-055d","JH-2-002","ST8814","NF10.1"),],
                      aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
    geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
    scale_color_manual(values=c("blue","black","red"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + 
    scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
    scale_x_continuous(transform="log10") + geom_point() +
    geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
    ggtitle(titles$chr8) +
    theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                       color = "Chr8q Status", shape = "MPNST Cell Line") + 
    theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8q_log10_v2_angled_JH5donlyMYC.pdf"),conf.plot,width=4,height=2.5)
  
  conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$improve_sample_id %in% c("JH-2-055d","JH-2-002","ST8814","NF10.1"),],
                      aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
    geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
    scale_color_manual(values=c("blue","black","red"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + 
    scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
    scale_x_continuous(transform="log10") + geom_point() +
    geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
    ggtitle(d) +
    theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                       color = "Chr8q Status", shape = "MPNST Cell Line") + 
    theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8q_log10_v2_angled_JH5donlyMYCv2.pdf"),conf.plot,width=4,height=2.5)
  
  conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$improve_sample_id %in% c("JH-2-055d","JH-2-002","ST8814","NF10.1"),],
                      aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
    geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
    scale_color_manual(values=c("blue","black","red"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + 
    scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
    scale_x_continuous(transform="log10") + geom_point() +
    geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
    ggtitle(d) +
    theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                       color = "Chr8q Status", shape = "MPNST Cell Line") + 
    theme(plot.title=element_text(face="bold",hjust=0.5, color="blue"), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8q_log10_v2_angled_JH5donlyMYCv2blue.pdf"),conf.plot,width=4,height=2.5)
  
  conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$improve_sample_id %in% c("JH-2-055d","JH-2-002","ST8814","NF10.1"),],
                      aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
    geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
    scale_color_manual(values=c("blue","black","red"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + 
    scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
    scale_x_continuous(transform="log10") + geom_point() +
    geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
    ggtitle(d) +
    theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                       color = "Chr8q Status", shape = "MPNST Cell Line") + 
    theme(plot.title=element_text(face="bold",hjust=0.5, color="red"), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8q_log10_v2_angled_JH5donlyMYCv2red.pdf"),conf.plot,width=4,height=2.5)
  
  conf.plot <- ggplot(mean.drug.conf,
                      aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
    geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
    scale_color_manual(values=c("blue","black","red"), breaks=c("Diploid","Diploid with MYC gain","Gain")) +  # was black, darkgrey, lightgrey
    scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
    scale_x_continuous(transform="log10") + geom_point() +
    geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
    ggtitle(titles$MYC) +
    theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                       color = "Chr8q Status", shape = "MPNST Cell Line") + 
    theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8qMYC_log10_v2_angled.pdf"),conf.plot,width=4,height=2.5)
  
  conf.plot <- ggplot(mean.drug.conf,
                      aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
    geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
    scale_color_manual(values=c("blue","black","red"), breaks=c("Diploid","Diploid with MYC gain","Gain")) +  # was black, darkgrey, lightgrey
    scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
    scale_x_continuous(transform="log10") + geom_point() +
    geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
    ggtitle(d) +
    theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                       color = "Chr8q Status", shape = "MPNST Cell Line") + 
    theme(plot.title=element_text(face="bold",hjust=0.5, color="blue"), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8qMYC_log10_v2_angledBlue.pdf"),conf.plot,width=4,height=2.5)
  ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8qMYC_log10_v2_angledBlue_height5width5.pdf"),conf.plot,width=5,height=5)
  
  conf.plot <- ggplot(mean.drug.conf,
                      aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
    geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
    scale_color_manual(values=c("blue","black","red"), breaks=c("Diploid","Diploid with MYC gain","Gain")) +  # was black, darkgrey, lightgrey
    scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
    scale_x_continuous(transform="log10") + geom_point() +
    geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
    ggtitle(d) +
    theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                       color = "Chr8q Status", shape = "MPNST Cell Line") + 
    theme(plot.title=element_text(face="bold",hjust=0.5, color="red"), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8qMYC_log10_v2_angledRed.pdf"),conf.plot,width=4,height=2.5)
  ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8qMYC_log10_v2_angledRed_height5width5.pdf"),conf.plot,width=5,height=5)
  
  conf.plot <- ggplot(mean.drug.conf,
                      aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
    geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
    scale_color_manual(values=c("blue","black","red"), breaks=c("Diploid","Diploid with MYC gain","Gain")) +  # was black, darkgrey, lightgrey
    scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
    scale_x_continuous(transform="log10") + geom_point() +
    geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
    ggtitle(d) +
    theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                       color = "Chr8q Status", shape = "MPNST Cell Line") + 
    theme(plot.title=element_text(face="bold",hjust=0.5, color="black"), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8qMYC_log10_v2_angledBlack.pdf"),conf.plot,width=4,height=2.5)
}
write.csv(all.p.df, "p_values_relPercViability_JHandWU.csv", row.names=FALSE)

#### combined IC50 pVals calculated in Excel ####
chr8.IC50.p <- data.frame(Target = c("PLK","PLK","EGFR","EGFR"), 
                          Drug=c("Volasertib","BI-2536","Erlotinib","Osimertinib"),
                     p=c(1.52159E-32, 0.013775419, 0.000146753, 4.35539E-05))
MYC.IC50.p <- data.frame(Target = c("PLK","PLK","EGFR","EGFR"), 
                          Drug=c("Volasertib","BI-2536","Erlotinib","Osimertinib"),
                          p=c(0.039410244, 0.050906472, 0.12883391, 0.038412181))
library(poolr)
chr8.plk <- poolr::fisher(chr8.IC50.p[chr8.IC50.p$Target=="PLK",]$p)$p
chr8.plk # 1.646405e-32
chr8.egfr <- poolr::fisher(chr8.IC50.p[chr8.IC50.p$Target=="EGFR",]$p)$p
chr8.egfr # 1.269913e-07
myc.plk <- poolr::fisher(MYC.IC50.p[MYC.IC50.p$Target=="PLK",]$p)$p
myc.plk # 0.01446796
myc.egfr <- poolr::fisher(MYC.IC50.p[MYC.IC50.p$Target=="EGFR",]$p)$p
myc.egfr # 0.03122
