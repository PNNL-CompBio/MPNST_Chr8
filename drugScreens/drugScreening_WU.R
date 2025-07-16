# analyzing IncuCyte results
library(plyr); library(dplyr); library(tidyr);library(ggplot2)
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/drugScreens/WU")
dataPath <- "~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/drugScreens/WU/Chr8 IncuCyte drug testing"
exp <- list.files(dataPath)
meta <- data.frame(exp)
meta <- meta %>% tidyr::separate_wider_delim(exp, " ", 
                                             names=c("author","MPNST",
                                                     "Drug1","Drug2","Date"),
                                             cols_remove=FALSE)
all.df <- data.frame()
for (i in exp) {
  # get names of files inside experiment folder
  exp.files <- list.files(file.path(dataPath, i))
  
  exp.df <- data.frame()
  for (j in exp.files) {
    # read each file 
    temp.file <- readxl::read_excel(file.path(dataPath, i, j), skip=7)
    
    # reformat long
    file.long <- reshape2::melt(temp.file, id.var = "Elapsed", 
                                variable.name="sample")
    file.long <- file.long[file.long$sample != "Date Time",]
    file.long$unit <- strsplit(strsplit(j, paste0(i, " "))[[1]][2],".xlsx")[[1]]
    
    # prep sample name format for extracting metadata
    file.long$sample <- sub(" / well", "", file.long$sample)
    file.long$sample <- sub("Con [(]vehicle, outer[)]", "Control", file.long$sample)
    file.long$sample <- sub("Con [(]vehicle, inner[)]", "Control", file.long$sample)
    file.long <- file.long %>% 
      tidyr::separate_wider_delim(sample, " ", names=c("MPNST", "passage",
                                                       "N_cells_plated", "Drug", 
                                                       "concentration", 
                                                       "concUnit", "well"),
                                  cols_remove=FALSE)
    exp.df <- rbind(exp.df, file.long)
  }
  
  # if norm confluence is missing, calculate it
  if (all(!grepl("Confluence Norm", exp.files))) {
    # get confluence
    temp.conf <- exp.df[grepl("confluence",exp.df$unit, ignore.case=TRUE),]
    conf0 <- temp.conf[temp.conf$Elapsed==0,]
    
    # divide by time 0 for each well
    temp.samples <- unique(conf0$sample)
    for (s in temp.samples) {
      temp.conf[temp.conf$sample==s,]$value <- temp.conf[temp.conf$sample==s,]$value / conf0[conf0$sample==s,]$value
    }
    
    # store result
    temp.conf$unit <- "Confluence Norm"
    exp.df <- rbind(exp.df, temp.conf)
  }
  
  exp.df$exp <- i
  
  # concatenate data
  all.df <- rbind(all.df, exp.df)
}
all.df$time_unit <- "h"
all.df$author <- "DB"
write.csv(all.df, paste0("DB_IncuCyte_data_", Sys.Date(),".csv"), row.names=FALSE)
all.df <- read.csv("DB_IncuCyte_data_2025-05-15.csv")

# add molecular weights to convert to micromolar
all.df$MW <- NA
all.df[all.df$Drug == "Volas",]$MW <- 618.8 # g/mol; source: https://pubchem.ncbi.nlm.nih.gov/compound/Volasertib
all.df[all.df$Drug == "Gefitinib",]$MW <- 446.9 # g/mol; source: https://pubchem.ncbi.nlm.nih.gov/compound/123631
all.df[all.df$Drug == "Control",]$MW <- 78.14 # g/mol for DMSO; source: https://pubchem.ncbi.nlm.nih.gov/compound/679
all.df[all.df$Drug == "AT7519",]$MW <- 382.2 # g/mol; source: https://pubchem.ncbi.nlm.nih.gov/compound/11338033
all.df[all.df$Drug == "Mirdametinib",]$MW <- 482.19 # g/mol; source: https://pubchem.ncbi.nlm.nih.gov/compound/9826528
all.df[all.df$Drug == "Ribociclib",]$MW <- 434.5 # g/mol; source: https://pubchem.ncbi.nlm.nih.gov/compound/44631912
all.df[all.df$Drug == "Alrizo",]$MW <- 642.6 # g/mol; source: https://pubchem.ncbi.nlm.nih.gov/compound/91972012
all.df[all.df$concUnit=="ÂµM",]$concUnit <- "um"

all.df$concUM <- NA
all.df[all.df$concUnit == "um",]$concUM <- all.df[all.df$concUnit == "um",]$concentration
all.df$concUM <- as.numeric(all.df$concUM)
all.df[all.df$concUnit == "mg/ml",]$concUM <- (10^6)*as.numeric(all.df[all.df$concUnit == "mg/ml",]$concentration) / as.numeric(all.df[all.df$concUnit == "mg/ml",]$MW)
all.df$concUnitUM <- "um"

chr8_tTest <- function(drug.df, unit) {
  # reduce to overlapping parameters (dose, time) and filter for cell types
  wu <- drug.df[drug.df$improve_sample_id == "WU-356",]
  wu.params <- unique(wu$sample)
  jh <- drug.df[drug.df$improve_sample_id == "JH-2-055d",]
  jh.params <- unique(jh$sample)
  shared.params <- wu.params[wu.params %in% jh.params]
  wu <- wu[wu$sample %in% shared.params,]
  jh <- jh[jh$sample %in% shared.params,]
  
  # compare paired results
  if ("GROWTH" %in% colnames(wu)) {
    wu <- wu[order(wu$sample),]$GROWTH
    jh <- jh[order(jh$sample),]$GROWTH
  } else {
    wu <- wu[order(wu$sample),]$dose_response_value
    jh <- jh[order(jh$sample),]$dose_response_value
  }
  if (length(wu) == length(jh)) {
    paired <- TRUE
  } else {paired <- FALSE}
  if (grepl("Death",unit, ignore.case=TRUE)) {
    if (median(wu) < median(jh)) {
      drug.p <- t.test(jh, wu, alternative="greater", paired=paired)$p.value 
      if (paired) {
        title <- paste0("JH-2-055d (diploid) is more sensitive to ",d," (paired p=", signif(drug.p,2),")")
      } else {
        title <- paste0("JH-2-055d (diploid) is more sensitive to ",d," (p=", signif(drug.p,2),")") 
      }
    } else {
      drug.p <- t.test(wu, jh, alternative="greater", paired=paired)$p.value 
      if (paired) {
        title <- paste0("WU-356 (gain) is more sensitive to ",d," (paired p=", signif(drug.p,2),")")
      } else {
        title <- paste0("WU-356 (gain) is more sensitive to ",d," (p=", signif(drug.p,2),")") 
      }
    } 
  } else {
    if (median(wu) > median(jh)) {
      drug.p <- t.test(wu, jh, alternative="greater", paired=paired)$p.value 
      if (paired) {
        title <- paste0("JH-2-055d (diploid) is more sensitive to ",d," (paired p=", signif(drug.p,2),")")
      } else {
        title <- paste0("JH-2-055d (diploid) is more sensitive to ",d," (p=", signif(drug.p,2),")") 
      }
    } else {
      drug.p <- t.test(jh, wu, alternative="greater", paired=paired)$p.value 
      if (paired) {
        title <- paste0("WU-356 (gain) is more sensitive to ",d," (paired p=", signif(drug.p,2),")")
      } else {
        title <- paste0("WU-356 (gain) is more sensitive to ",d," (p=", signif(drug.p,2),")") 
      }
    }
  }
  return(title)
}
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
      drug.df$DOSEum <- paste0(drug.df$DOSE, "uM")
      drug.df$DOSEum <- factor(drug.df$DOSEum, levels=paste0(sort(unique(drug.df$DOSE)), "uM"))
      drug.df$sample <- paste0(drug.df$DOSEum, "uM after ", drug.df$time,"h") 
    } else {
      drug.df <- death.df[death.df$improve_drug_id==d,]
      drug.df$sample <- paste0(drug.df$time,"h")
    }
    
    title <- chr8_tTest(drug.df, unit)
    if ("GROWTH" %in% colnames(drug.df)) {
      ggplot2::ggplot(drug.df, aes(x=time, y=GROWTH, color=improve_sample_id)) + 
        geom_point() + facet_grid(.~DOSEum) + ggtitle(title) + theme_classic() +
        labs(y=paste("Normalized",unit), x="Time (h)", color="MPNST") + 
        theme(plot.title=element_text(hjust=0.5))
      ggsave(paste0(tolower(unit),"_",d,"_no_4-23-25.pdf"),width=12, height=6) 
    } else {
      ggplot2::ggplot(drug.df, aes(x=time, y=dose_response_value, color=improve_sample_id)) + 
        geom_point() + ggtitle(title) + theme_classic() +
        labs(y=paste("Normalized",unit), x="Time (h)", color="MPNST") + 
        theme(plot.title=element_text(hjust=0.5))
      ggsave(paste0(tolower(unit),"_",d,"_no_4-23-25.pdf"),width=12, height=6)
    }
  }
}

# get curve fitting script
if (!file.exists("fit_curve.py")) {download.file('https://raw.githubusercontent.com/PNNL-CompBio/coderdata/refs/heads/main/build/utils/fit_curve.py',
              destfile="fit_curve.py")} # comment out line 268 (don't need to divide values by 100)
maxConc <- c(20) # uM; all non-control concentrations: 0.04  0.08  0.16  0.32  0.63  1.25  2.50  5.00 10.00 20.00
for (m in maxConc) {
  # change colnames to use existing curve fitting script
  confluence.df <- all.df[grepl("Confluence Norm",all.df$unit) & all.df$Drug != "Control" &
                            all.df$concUM <= m & 
                            all.df$exp != "DB JH-2-055d Gef Ribo 4-23-25" &
                            all.df$exp != "DB WU-356 Gef Ribo 4-23-25",] # 20 um may be too high for gefitinib
  confluence.df$study <- "Chr8"
  #oldCols <- c("concentration","value", "study", "author", "MPNST", "Drug", "Elapsed", "time_unit")
  oldCols <- c("concUM","value", "study", "author", "MPNST", "Drug", "Elapsed", "time_unit")
  newCols <- c('DOSE','GROWTH','study','source','improve_sample_id','Drug','time','time_unit')
  confluence.df <- confluence.df %>% rename_at(vars(oldCols), ~ newCols)
  write.table(confluence.df[,newCols], paste0("chr8_normConfluence_max",m,"um_no_4-23-25.tsv"), sep="\t", row.names=FALSE)
  
  # plot each drug
  drugPlots(confluence.df, "Confluence")
  
  # run curve fitting
  system2("python3", paste0("fit_curve.py --input chr8_normConfluence_max",m,
                            "um_no_4-23-25.tsv --output chr8_normConfluenceCurves_max",m,"um_no_4-23-25"))
  
  confluence <- read.table(paste0("chr8_normConfluenceCurves_max",m,"um_no_4-23-25.0"), header=TRUE, fill=NA)
  ## plot
  # AUC
  auc <- confluence[confluence$dose_response_metric=="fit_auc",]
  drugPlots(auc, "Confluence Area Under the Curve")
  quantiles <- quantile(auc$dose_response_value)
  auc$Quantile <- NA
  for (q in names(quantiles)) {
    if (any(auc$dose_response_value >= quantiles[[q]])){
      auc[auc$dose_response_value >= quantiles[[q]],]$Quantile <- q 
    }
  }
  auc$Quantile <- factor(auc$Quantile, names(quantiles))
  #auc$time <- paste0(auc$time, "h")
  #auc$time <- factor(auc$time, levels=c("48h", "120h"))
  mean.auc <- plyr::ddply(auc, .(improve_drug_id), summarize,
                          auc=mean(dose_response_value, na.rm=TRUE))
  
  # r2
  r2 <- confluence[confluence$dose_response_metric=="fit_r2",]
  quantiles <- quantile(r2$dose_response_value)
  r2$Quantile <- NA
  for (q in names(quantiles)) {
    if (any(r2$dose_response_value >= quantiles[[q]])){
      r2[r2$dose_response_value >= quantiles[[q]],]$Quantile <- q 
    }
  }
  r2$Quantile <- factor(r2$Quantile, names(quantiles))
  
  # AUC with r2 fill
  auc.r2 <- merge(auc, r2, by=c("improve_drug_id","improve_sample_id", "time"))
  auc.r2[auc.r2$improve_sample_id=="JH-2-055d",]$improve_sample_id <- "JH-2-055d (diploid)"
  auc.r2[auc.r2$improve_sample_id=="WU-356",]$improve_sample_id <- "WU-356 (gain)"
  drugOrder <- unique(mean.auc[order(mean.auc$auc),]$improve_drug_id)
  ggplot2::ggplot(auc.r2, 
                  aes(x=time, y=dose_response_value.x, color=dose_response_value.y, shape=improve_sample_id)) + 
    geom_point() + theme_classic(base_size=12) + facet_grid(.~improve_drug_id)+
    labs(y="AUC: Normalized Confluence", x= "Time (h)", color=expression(paste("Fit R"^"2")), shape="MPNST") + 
    theme(#axis.title.x=element_blank(), 
      axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0("fit_auc_r2_confluence_dotPlot_max",m,"um_no_4-23-25.pdf"), width=10, height=3)
  
  # also do relative confluence by dividing by control at each timepoint
  # change colnames to use existing curve fitting script
  confluence.df <- all.df[grepl("Confluence Norm",all.df$unit) & 
                            all.df$concUM <= m & 
                            all.df$exp != "DB JH-2-055d Gef Ribo 4-23-25" &
                            all.df$exp != "DB WU-356 Gef Ribo 4-23-25",] # 20 um may be too high for gefitinib
  #ctrl.confl <- confluence.df[confluence.df$Drug == "Control",] # 0 rows
  ctrl.confl <- all.df[grepl("Confluence Norm",all.df$unit) & 
                         #all.df$concUM <= m & # this was filtering out DMSO which matches max conc (20 uM)
                         all.df$Drug == "Control" &
                         all.df$exp != "DB JH-2-055d Gef Ribo 4-23-25" &
                         all.df$exp != "DB WU-356 Gef Ribo 4-23-25",]
  for (cellLine in unique(confluence.df$MPNST)) {
    for (expt in unique(confluence.df$exp)) {
      for (t in unique(confluence.df$Elapsed)) {
        temp.val <- mean(ctrl.confl[ctrl.confl$Elapsed == t & ctrl.confl$MPNST == cellLine &
                                      ctrl.confl$exp == expt,]$value, na.rm=TRUE)
        confluence.df[confluence.df$Elapsed == t & confluence.df$MPNST == cellLine & 
                        confluence.df$exp == expt,]$value <- 
          (confluence.df[confluence.df$Elapsed == t & confluence.df$MPNST == cellLine & 
                           confluence.df$exp == expt,]$value)/temp.val
      }  
    }
  }
  confluence.df <- confluence.df[confluence.df$Drug != "Control",]
  confluence.df$study <- "Chr8"
  #oldCols <- c("concentration","value", "study", "author", "MPNST", "Drug", "Elapsed", "time_unit")
  oldCols <- c("concUM","value", "study", "author", "MPNST", "Drug", "Elapsed", "time_unit")
  newCols <- c('DOSE','GROWTH','study','source','improve_sample_id','Drug','time','time_unit')
  confluence.df <- confluence.df %>% rename_at(vars(oldCols), ~ newCols)
  write.table(confluence.df[,newCols], paste0("chr8_relConfluence_max",m,"um_no_4-23-25.tsv"), sep="\t", row.names=FALSE)
  
  # plot each drug
  drugPlots(confluence.df, "Relative Confluence")
  
  # run curve fitting
  system2("python3", paste0("fit_curve.py --input chr8_relConfluence_max",m,
                            "um_no_4-23-25.tsv --output chr8_relConfluenceCurves_max",m,"um_no_4-23-25"))
  
  confluence <- read.table(paste0("chr8_relConfluenceCurves_max",m,"um_no_4-23-25.0"), header=TRUE, fill=NA)
  ## plot
  # AUC
  auc <- confluence[confluence$dose_response_metric=="fit_auc",]
  drugPlots(auc, "Relative Confluence Area Under the Curve")
  quantiles <- quantile(auc$dose_response_value)
  auc$Quantile <- NA
  for (q in names(quantiles)) {
    if (any(auc$dose_response_value >= quantiles[[q]])){
      auc[auc$dose_response_value >= quantiles[[q]],]$Quantile <- q 
    }
  }
  auc$Quantile <- factor(auc$Quantile, names(quantiles))
  #auc$time <- paste0(auc$time, "h")
  #auc$time <- factor(auc$time, levels=c("48h", "120h"))
  mean.auc <- plyr::ddply(auc, .(improve_drug_id), summarize,
                          auc=mean(dose_response_value, na.rm=TRUE))
  
  # r2
  r2 <- confluence[confluence$dose_response_metric=="fit_r2",]
  quantiles <- quantile(r2$dose_response_value)
  r2$Quantile <- NA
  for (q in names(quantiles)) {
    if (any(r2$dose_response_value >= quantiles[[q]])){
      r2[r2$dose_response_value >= quantiles[[q]],]$Quantile <- q 
    }
  }
  r2$Quantile <- factor(r2$Quantile, names(quantiles))
  
  # AUC with r2 fill
  auc.r2 <- merge(auc, r2, by=c("improve_drug_id","improve_sample_id", "time"))
  auc.r2[auc.r2$improve_sample_id=="JH-2-055d",]$improve_sample_id <- "JH-2-055d (diploid)"
  auc.r2[auc.r2$improve_sample_id=="WU-356",]$improve_sample_id <- "WU-356 (gain)"
  drugOrder <- unique(mean.auc[order(mean.auc$auc),]$improve_drug_id)
  ggplot2::ggplot(auc.r2, 
                  aes(x=time, y=dose_response_value.x, color=dose_response_value.y, shape=improve_sample_id)) + 
    geom_point() + theme_classic(base_size=12) + facet_grid(.~improve_drug_id)+
    labs(y="AUC: Normalized Confluence", x= "Time (h)", color=expression(paste("Fit R"^"2")), shape="MPNST") + 
    theme(#axis.title.x=element_blank(), 
      axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0("fit_auc_r2_relconfluence_dotPlot_max",m,"um_no_4-23-25.pdf"), width=10, height=3)
  hist(r2$dose_response_value)
  
  death.df <- all.df[grepl("norm",all.df$unit, ignore.case=TRUE) & 
                       !grepl("confluence", all.df$unit, ignore.case=TRUE) & 
                       all.df$Drug != "Control" & all.df$concUM <= 10 &
                       all.df$exp != "DB JH-2-055d Gef Ribo 4-23-25" &
                       all.df$exp != "DB WU-356 Gef Ribo 4-23-25",]
  death.df$study <- "Chr8"
  #oldCols <- c("concentration","value", "study", "author", "MPNST", "Drug", "Elapsed", "time_unit")
  oldCols <- c("concUM","value", "study", "author", "MPNST", "Drug", "Elapsed", "time_unit")
  newCols <- c('DOSE','GROWTH','study','source','improve_sample_id','Drug','time','time_unit')
  death.df <- death.df %>% rename_at(vars(oldCols), ~ newCols)
  write.table(death.df[,newCols], paste0("chr8_normDeath_max",m,"um_no_4-23-25.tsv"), sep="\t", row.names=FALSE)
  
  # plot each drug
  drugPlots(death.df, "Death")
  
  # run curve fitting
  system2("python3", paste0("fit_curve.py --input chr8_normDeath_max",m,
                            "um_no_4-23-25.tsv --output chr8_normDeathCurves_max",m,"um_no_4-23-25"))
  death <- read.table(paste0("chr8_normDeathCurves_max",m,"um_no_4-23-25.0"), header=TRUE, fill=NA)
  ## plot
  # AUC
  auc <- death[death$dose_response_metric=="fit_auc",]
  drugPlots(auc, "Death Area Under the Curve")
  quantiles <- quantile(auc$dose_response_value)
  auc$Quantile <- NA
  for (q in names(quantiles)) {
    if (any(auc$dose_response_value >= quantiles[[q]])){
      auc[auc$dose_response_value >= quantiles[[q]],]$Quantile <- q 
    }
  }
  auc$Quantile <- factor(auc$Quantile, names(quantiles))
  #auc$time <- paste0(auc$time, "h")
  #auc$time <- factor(auc$time, levels=c("48h", "120h"))
  mean.auc <- plyr::ddply(auc, .(improve_drug_id), summarize,
                          auc=mean(dose_response_value, na.rm=TRUE))
  
  # r2
  r2 <- death[death$dose_response_metric=="fit_r2",]
  quantiles <- quantile(r2$dose_response_value)
  r2$Quantile <- NA
  for (q in names(quantiles)) {
    if (any(r2$dose_response_value >= quantiles[[q]])){
      r2[r2$dose_response_value >= quantiles[[q]],]$Quantile <- q 
    }
  }
  r2$Quantile <- factor(r2$Quantile, names(quantiles))
  
  # AUC with r2 fill
  auc.r2 <- merge(auc, r2, by=c("improve_drug_id","improve_sample_id", "time"))
  auc.r2[auc.r2$improve_sample_id=="JH-2-055d",]$improve_sample_id <- "JH-2-055d (diploid)"
  auc.r2[auc.r2$improve_sample_id=="WU-356",]$improve_sample_id <- "WU-356 (gain)"
  drugOrder <- unique(mean.auc[order(mean.auc$auc),]$improve_drug_id)
  ggplot2::ggplot(auc.r2, 
                  aes(x=time, y=dose_response_value.x, color=dose_response_value.y, shape=improve_sample_id)) + 
    geom_point() + theme_classic(base_size=12) + facet_grid(.~improve_drug_id)+
    labs(y="AUC: Normalized Death", x= "Time (h)", color=expression(paste("Fit R"^"2")), shape="MPNST") + 
    theme(#axis.title.x=element_blank(), 
      axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0("fit_auc_r2_death_dotPlot_max",m,"um_no_4-23-25.pdf"), width=10, height=3)
  
  death.df <- all.df[grepl("norm",all.df$unit, ignore.case=TRUE) & 
                       !grepl("confluence", all.df$unit, ignore.case=TRUE) & 
                       all.df$Drug != "Control" & all.df$concUM <= 10 &
                       all.df$exp != "DB JH-2-055d Gef Ribo 4-23-25" &
                       all.df$exp != "DB WU-356 Gef Ribo 4-23-25",]
  ctrl.death <- all.df[grepl("norm",all.df$unit, ignore.case=TRUE) & 
                         !grepl("confluence", all.df$unit, ignore.case=TRUE) &
                         #all.df$concUM <= m & # this was filtering out DMSO which matches max conc (20 uM)
                         all.df$Drug == "Control" &
                         all.df$exp != "DB JH-2-055d Gef Ribo 4-23-25" &
                         all.df$exp != "DB WU-356 Gef Ribo 4-23-25",]
  for (cellLine in unique(death.df$MPNST)) {
    for (expt in unique(death.df$exp)) {
      for (t in unique(death.df$Elapsed)) {
        temp.val <- mean(ctrl.death[ctrl.death$Elapsed == t & ctrl.death$MPNST == cellLine &
                                      ctrl.death$exp == expt,]$value, na.rm=TRUE)
        death.df[death.df$Elapsed == t & death.df$MPNST == cellLine & 
                   death.df$exp == expt,]$value <- 
          (death.df[death.df$Elapsed == t & death.df$MPNST == cellLine & 
                      death.df$exp == expt,]$value)/temp.val
      }  
    }
  }
  death.df <- death.df[death.df$Drug != "Control",]
  death.df$study <- "Chr8"
  #oldCols <- c("concentration","value", "study", "author", "MPNST", "Drug", "Elapsed", "time_unit")
  oldCols <- c("concUM","value", "study", "author", "MPNST", "Drug", "Elapsed", "time_unit")
  newCols <- c('DOSE','GROWTH','study','source','improve_sample_id','Drug','time','time_unit')
  death.df <- death.df %>% rename_at(vars(oldCols), ~ newCols)
  write.table(death.df[,newCols], paste0("chr8_relDeath_max",m,"um_no_4-23-25.tsv"), sep="\t", row.names=FALSE)
  
  # plot each drug
  drugPlots(death.df, "Relative Death")
  
  # run curve fitting
  system2("python3", paste0("fit_curve.py --input chr8_relDeath_max",m,
                            "um_no_4-23-25.tsv --output chr8_relDeathCurves_max",m,"um_no_4-23-25"))
  death <- read.table(paste0("chr8_relDeathCurves_max",m,"um_no_4-23-25.0"), header=TRUE, fill=NA)
  ## plot
  # AUC
  auc <- death[death$dose_response_metric=="fit_auc",]
  drugPlots(auc, "Relative Death Area Under the Curve")
  quantiles <- quantile(auc$dose_response_value)
  auc$Quantile <- NA
  for (q in names(quantiles)) {
    if (any(auc$dose_response_value >= quantiles[[q]])){
      auc[auc$dose_response_value >= quantiles[[q]],]$Quantile <- q 
    }
  }
  auc$Quantile <- factor(auc$Quantile, names(quantiles))
  #auc$time <- paste0(auc$time, "h")
  #auc$time <- factor(auc$time, levels=c("48h", "120h"))
  mean.auc <- plyr::ddply(auc, .(improve_drug_id), summarize,
                          auc=mean(dose_response_value, na.rm=TRUE))
  
  # r2
  r2 <- death[death$dose_response_metric=="fit_r2",]
  quantiles <- quantile(r2$dose_response_value)
  r2$Quantile <- NA
  for (q in names(quantiles)) {
    if (any(r2$dose_response_value >= quantiles[[q]])){
      r2[r2$dose_response_value >= quantiles[[q]],]$Quantile <- q 
    }
  }
  r2$Quantile <- factor(r2$Quantile, names(quantiles))
  
  # AUC with r2 fill
  auc.r2 <- merge(auc, r2, by=c("improve_drug_id","improve_sample_id", "time"))
  auc.r2[auc.r2$improve_sample_id=="JH-2-055d",]$improve_sample_id <- "JH-2-055d (diploid)"
  auc.r2[auc.r2$improve_sample_id=="WU-356",]$improve_sample_id <- "WU-356 (gain)"
  drugOrder <- unique(mean.auc[order(mean.auc$auc),]$improve_drug_id)
  ggplot2::ggplot(auc.r2, 
                  aes(x=time, y=dose_response_value.x, color=dose_response_value.y, shape=improve_sample_id)) + 
    geom_point() + theme_classic(base_size=12) + facet_grid(.~improve_drug_id)+
    labs(y="AUC: Normalized Death", x= "Time (h)", color=expression(paste("Fit R"^"2")), shape="MPNST") + 
    theme(#axis.title.x=element_blank(), 
      axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  ggsave(paste0("fit_auc_r2_reldeath_dotPlot_max",m,"um_no_4-23-25.pdf"), width=10, height=3)
}

confluence <- read.table("chr8_normConfluenceCurves.0", header=TRUE, fill=NA)
## plot
# AUC
auc <- confluence[confluence$dose_response_metric=="fit_auc",]
quantiles <- quantile(auc$dose_response_value)
auc$Quantile <- NA
for (q in names(quantiles)) {
  if (any(auc$dose_response_value >= quantiles[[q]])){
    auc[auc$dose_response_value >= quantiles[[q]],]$Quantile <- q 
  }
}
auc$Quantile <- factor(auc$Quantile, names(quantiles))
#auc$time <- paste0(auc$time, "h")
#auc$time <- factor(auc$time, levels=c("48h", "120h"))
mean.auc <- plyr::ddply(auc, .(improve_drug_id), summarize,
                        auc=mean(dose_response_value, na.rm=TRUE))

# r2
r2 <- confluence[confluence$dose_response_metric=="fit_r2",]
quantiles <- quantile(r2$dose_response_value)
r2$Quantile <- NA
for (q in names(quantiles)) {
  if (any(r2$dose_response_value >= quantiles[[q]])){
    r2[r2$dose_response_value >= quantiles[[q]],]$Quantile <- q 
  }
}
r2$Quantile <- factor(r2$Quantile, names(quantiles))

# AUC with r2 fill
auc.r2 <- merge(auc, r2, by=c("improve_drug_id","improve_sample_id", "time"))
auc.r2[auc.r2$improve_sample_id=="JH-2-055d",]$improve_sample_id <- "JH-2-055d (diploid)"
auc.r2[auc.r2$improve_sample_id=="WU-356",]$improve_sample_id <- "WU-356 (gain)"
drugOrder <- unique(mean.auc[order(mean.auc$auc),]$improve_drug_id)
ggplot2::ggplot(auc.r2, 
                aes(x=time, y=dose_response_value.x, color=dose_response_value.y, shape=improve_sample_id)) + 
  geom_point() + theme_classic(base_size=12) + facet_grid(.~improve_drug_id)+
  labs(y="AUC: Normalized Confluence", x= "Time (h)", color=expression(paste("Fit R"^"2")), shape="MPNST") + 
  theme(#axis.title.x=element_blank(), 
    axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("fit_auc_r2_confluence_dotPlot.pdf", width=7, height=3)

ggplot2::ggplot(auc.r2[auc.r2$dose_response_value.y>=0.5,], 
                aes(x=time, y=dose_response_value.x, color=dose_response_value.y, shape=improve_sample_id)) + 
  geom_point() + theme_classic(base_size=12) + facet_grid(.~improve_drug_id)+
  labs(y="AUC: Normalized Confluence", x= "Time (h)", color=expression(paste("Fit R"^"2")), shape="MPNST") + 
  theme(#axis.title.x=element_blank(), 
    axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("fit_auc_r2Min0.5_confluence_dotPlot.pdf", width=7, height=3)

death <- read.table("chr8_normDeathCurves.0", header=TRUE, fill=NA)
## plot
# AUC
auc <- death[death$dose_response_metric=="fit_auc",]
quantiles <- quantile(auc$dose_response_value)
auc$Quantile <- NA
for (q in names(quantiles)) {
  if (any(auc$dose_response_value >= quantiles[[q]])){
    auc[auc$dose_response_value >= quantiles[[q]],]$Quantile <- q 
  }
}
auc$Quantile <- factor(auc$Quantile, names(quantiles))
#auc$time <- paste0(auc$time, "h")
#auc$time <- factor(auc$time, levels=c("48h", "120h"))
mean.auc <- plyr::ddply(auc, .(improve_drug_id), summarize,
                        auc=mean(dose_response_value, na.rm=TRUE))

# r2
r2 <- death[death$dose_response_metric=="fit_r2",]
quantiles <- quantile(r2$dose_response_value)
r2$Quantile <- NA
for (q in names(quantiles)) {
  if (any(r2$dose_response_value >= quantiles[[q]])){
    r2[r2$dose_response_value >= quantiles[[q]],]$Quantile <- q 
  }
}
r2$Quantile <- factor(r2$Quantile, names(quantiles))

# AUC with r2 fill
auc.r2 <- merge(auc, r2, by=c("improve_drug_id","improve_sample_id", "time"))
auc.r2[auc.r2$improve_sample_id=="JH-2-055d",]$improve_sample_id <- "JH-2-055d (diploid)"
auc.r2[auc.r2$improve_sample_id=="WU-356",]$improve_sample_id <- "WU-356 (gain)"
drugOrder <- unique(mean.auc[order(mean.auc$auc),]$improve_drug_id)
ggplot2::ggplot(auc.r2, 
                aes(x=time, y=dose_response_value.x, color=dose_response_value.y, shape=improve_sample_id)) + 
  geom_point() + theme_classic(base_size=12) + facet_grid(.~improve_drug_id)+
  labs(y="AUC: Normalized Death", x= "Time (h)", color=expression(paste("Fit R"^"2")), shape="MPNST") + 
  theme(#axis.title.x=element_blank(), 
    axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("fit_auc_r2_death_dotPlot.pdf", width=10, height=3)

ggplot2::ggplot(auc.r2[auc.r2$dose_response_value.y>=0.5,], 
                aes(x=time, y=dose_response_value.x, color=dose_response_value.y, shape=improve_sample_id)) + 
  geom_point() + theme_classic(base_size=12) + facet_grid(.~improve_drug_id)+
  labs(y="AUC: Normalized Death", x= "Time (h)", color=expression(paste("Fit R"^"2")), shape="MPNST") + 
  theme(#axis.title.x=element_blank(), 
    axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("fit_auc_r2Min0.5_death_dotPlot.pdf", width=10, height=3)

# what is going on with Gefitinib fits?
confluence <- read.table("chr8_normConfluence.tsv", header=TRUE)
confluence$unit <- "Normalized Confluence"
death <- read.table("chr8_normDeath.tsv", header=TRUE)
death$unit <- "Normalized Death"
gef <- rbind(confluence[confluence$Drug=="Gefitinib",], death[death$Drug == "Gefitinib",])
ggplot2::ggplot(gef, aes(x=time, y=GROWTH, shape=improve_sample_id, color=unit)) +geom_point()+ facet_grid(.~DOSE)
ggsave("confluenceAndDeath_Gefitinib.pdf",width=12, height=6)

# crowded plot, just look at JH-2-055d which had the worst fits
ggplot2::ggplot(gef[gef$improve_sample_id=="JH-2-055d",], aes(x=time, y=GROWTH, shape=improve_sample_id, color=unit)) +geom_point()+ facet_grid(.~DOSE)
ggsave("confluenceAndDeath_Gefitinib_JH-2-055d.pdf",width=12, height=6)

ggplot2::ggplot(gef[gef$improve_sample_id=="WU-356",], aes(x=time, y=GROWTH, shape=improve_sample_id, color=unit)) +geom_point()+ facet_grid(.~DOSE)
ggsave("confluenceAndDeath_Gefitinib_WU-356.pdf",width=12, height=6)

ggplot2::ggplot(gef[gef$improve_sample_id=="JH-2-055d" & gef$DOSE<=10,], aes(x=time, y=GROWTH, shape=improve_sample_id, color=unit)) +geom_point()+ facet_grid(.~DOSE)
ggsave("confluenceAndDeath_Gefitinib_JH-2-055d_max10um.pdf",width=12, height=6)

ggplot2::ggplot(gef[gef$improve_sample_id=="WU-356" & gef$DOSE<=10,], aes(x=time, y=GROWTH, shape=improve_sample_id, color=unit)) +geom_point()+ facet_grid(.~DOSE)
ggsave("confluenceAndDeath_Gefitinib_WU-356_max10um.pdf",width=12, height=6)


vol <- rbind(confluence[confluence$Drug=="Volas",], death[death$Drug == "Volas",])
ggplot2::ggplot(vol, aes(x=time, y=GROWTH, shape=improve_sample_id, color=unit)) +geom_point()+ facet_grid(.~DOSE)
ggsave("confluenceAndDeath_Volas.pdf",width=12, height=6)

# crowded plot, just look at JH-2-055d which had the worst fits
ggplot2::ggplot(vol[vol$improve_sample_id=="JH-2-055d",], aes(x=time, y=GROWTH, shape=improve_sample_id, color=unit)) +geom_point()+ facet_grid(.~DOSE)
ggsave("confluenceAndDeath_Volas_JH-2-055d.pdf",width=12, height=6)

ggplot2::ggplot(vol[vol$improve_sample_id=="WU-356",], aes(x=time, y=GROWTH, shape=improve_sample_id, color=unit)) +geom_point()+ facet_grid(.~DOSE)
ggsave("confluenceAndDeath_Volas_WU-356.pdf",width=12, height=6)

ggplot2::ggplot(vol[vol$improve_sample_id=="JH-2-055d" & vol$DOSE<=10,], aes(x=time, y=GROWTH, shape=improve_sample_id, color=unit)) +geom_point()+ facet_grid(.~DOSE)
ggsave("confluenceAndDeath_Volas_JH-2-055d_max10um.pdf",width=12, height=6)

ggplot2::ggplot(vol[vol$improve_sample_id=="WU-356" & vol$DOSE<=10,], aes(x=time, y=GROWTH, shape=improve_sample_id, color=unit)) +geom_point()+ facet_grid(.~DOSE)
ggsave("confluenceAndDeath_Volas_WU-356_max10um.pdf",width=12, height=6)

# gefitinib death peaks very early, then decreases
# 20 uM seems very strong for both gefitinib and volasertib

# repeat for all drugs
drugs <- unique(confluence$Drug)
for (d in drugs) {
  drug.df <- confluence[confluence$Drug==d,]
  ggplot2::ggplot(drug.df, aes(x=time, y=GROWTH, color=improve_sample_id)) +geom_point()+ facet_grid(.~DOSE)
  ggsave(paste0("confluence_",d,".pdf"),width=12, height=6)
  
  ggplot2::ggplot(drug.df[drug.df$DOSE <= 10,], aes(x=time, y=GROWTH, color=improve_sample_id)) +geom_point()+ facet_grid(.~DOSE)
  ggsave(paste0("confluence_",d,"_max10um.pdf"),width=12, height=6)
}

# is WU-356 (gain) more sensitive to volasertib than JH-2-055d (diploid)?
wu.vol.death <- death[death$Drug == "Volas" & death$improve_sample_id == "WU-356",]
wu.vol.death <- wu.vol.death[order(wu.vol.death$time),]$GROWTH
jh.vol.death <- death[death$Drug == "Volas" & death$improve_sample_id == "JH-2-055d",]
jh.vol.death <- jh.vol.death[order(jh.vol.death$time),]$GROWTH
vol.test <- t.test(wu.vol.death, jh.vol.death, alternative="greater", paired=TRUE)
# wu.vol.death is greater (p=4.0424E-21; WU-256 mean 3.04393, JH-2-055d mean 2.262132) - supports hypothesis

# is WU-356 (gain) more sensitive to gefitinib than JH-2-055d (diploid)?
wu.gef.death <- death[death$Drug == "Gefitinib" & death$improve_sample_id == "WU-356",]
wu.gef.death <- wu.gef.death[order(wu.gef.death$time),]$GROWTH
jh.gef.death <- death[death$Drug == "Gefitinib" & death$improve_sample_id == "JH-2-055d",]
jh.gef.death <- jh.gef.death[order(jh.gef.death$time),]$GROWTH
gef.test <- t.test(wu.gef.death, jh.gef.death, alternative="greater", paired=TRUE)
# wu.vol.death is greater (p=3.948041e-274; WU-256 mean 1.39417, JH-2-055d mean 0.7310665) - does not support hypothesis
mean(wu.gef.death)
mean(jh.gef.death)

#### also just plot 5d to match data from JH ####
# get data
rel.conf <- read.table("WU/chr8_relConfluence_max20um_no_4-23-25.tsv", sep="\t", header=TRUE)
rel.death <- read.table("WU/chr8_relDeath_max20um_no_4-23-25.tsv", sep="\t", header=TRUE)

# calculate mean and sd for each drug, dose, mpnst, time combo; also preserve sample and chr8q columns
mean.conf <- plyr::ddply(rel.conf, .(improve_sample_id, time, Drug, DOSE), summarize,
                         meanGROWTH = mean(GROWTH, na.rm=TRUE),
                         sdGROWTH = sd(GROWTH, na.rm=TRUE))
mean.death <- plyr::ddply(rel.death, .(improve_sample_id, time, Drug, DOSE), summarize,
                         meanGROWTH = mean(GROWTH, na.rm=TRUE),
                         sdGROWTH = sd(GROWTH, na.rm=TRUE))

# annotate chr8q status
mean.conf$chr8q <- "Amplified"
mean.conf[mean.conf$improve_sample_id == "JH-2-055d",]$chr8q <- "Diploid"
mean.death$chr8q <- "Amplified"
mean.death[mean.death$improve_sample_id == "JH-2-055d",]$chr8q <- "Diploid"
rel.conf$chr8q <- "Amplified"
rel.conf[rel.conf$improve_sample_id == "JH-2-055d",]$chr8q <- "Diploid"
rel.death$chr8q <- "Amplified"
rel.death[rel.death$improve_sample_id == "JH-2-055d",]$chr8q <- "Diploid"

# add sample column with all parameters
mean.conf$sample <- paste0(mean.conf$Drug,
                          "_",mean.conf$DOSE,"uM")
mean.death$sample <- paste0(mean.death$Drug,
                           "_",mean.death$DOSE,"uM")
rel.conf$sample <- paste0(rel.conf$Drug,
                           "_",rel.conf$DOSE,"uM")
rel.death$sample <- paste0(rel.death$Drug,
                            "_",rel.death$DOSE,"uM")

#times <- unique(c(rel.conf$time, rel.death$time)) # too many
times <- c(24,48,72,96,120)
setwd("WU")
dir.create("singleTimePlots")
setwd("singleTimePlots")
library(drc) # curve source: answer by greenjune: https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc
# Sara shared these 2 links: https://stackoverflow.com/questions/68209998/plot-drc-using-ggplot; https://forum.posit.co/t/extract-out-points-in-dose-response-curve-produced-by-drc/159433
Metric <- c("Relative Confluence", "Relative Death")
all.p.df <- data.frame()
for (t in times) {
  rel.conf.96 <- rel.conf[rel.conf$time==t,]
  rel.death.96 <- rel.death[rel.death$time==t,]
  mean.conf.96 <- mean.conf[mean.conf$time==t,]
  mean.death.96 <- mean.death[mean.death$time==t,]
  drugs <- unique(c(mean.conf.96$Drug, mean.death.96$Drug)) # all 6 drugs are now present: AT7519, Mirdametinib, Gefitinib, Ribociclib, Alrizo, Volas
  mpnst <- unique(c(mean.conf.96$improve_sample_id,mean.death.96$improve_sample_id)) # all 2 cell lines: JH-2-055d and WU-356
  
  # plot relative confluence and relative death for each drug
  for (d in drugs) {
    # filter for drug
    drug.conf <- rel.conf.96[rel.conf.96$Drug == d,]
    drug.death <- rel.death.96[rel.death.96$Drug == d,]
    mean.drug.conf <- mean.conf.96[mean.conf.96$Drug == d,]
    mean.drug.death <- mean.death.96[mean.death.96$Drug == d,]
    
    p.df <- data.frame(Metric, moreSensitive = NA, p = NA, Drug=d, Time=t, maxConcUM=20, minValue=NA)
    # run t-tests between WU-356 (chr8q-gain) and JH-2-055d (chr8q-diploid)
    conf.title <- chr8_tTest(drug.conf, "Relative Confluence")
    conf.more.sens <- strsplit(conf.title, " ")[[1]][1]
    conf.p <- strsplit(conf.title, "p=")[[1]][2]
    conf.p <- as.numeric(substr(conf.p,1,nchar(conf.p)-1))
    p.df[p.df$Metric=="Relative Confluence",]$moreSensitive <- conf.more.sens
    p.df[p.df$Metric=="Relative Confluence",]$p <- conf.p
    p.df[p.df$Metric=="Relative Confluence",]$minValue <- min(drug.conf$GROWTH, na.rm=TRUE)
    
    death.title <- chr8_tTest(drug.death, "Relative Death")
    death.more.sens <- strsplit(death.title, " ")[[1]][1]
    death.p <- strsplit(death.title, "p=")[[1]][2]
    death.p <- as.numeric(substr(death.p,1,nchar(death.p)-1))
    p.df[p.df$Metric=="Relative Death",]$moreSensitive <- death.more.sens
    p.df[p.df$Metric=="Relative Death",]$p <- death.p
    p.df[p.df$Metric=="Relative Death",]$minValue <- min(drug.death$GROWTH, na.rm=TRUE)
    all.p.df <- rbind(all.p.df, p.df)
    
    # generate plots
    conf.plot <- ggplot(mean.drug.conf, aes(x=DOSE, y=meanGROWTH, shape=chr8q, color=improve_sample_id)) +
      geom_smooth(linetype="dashed",se=FALSE, method=drc::drm, method.args=list(fct=LL.4())) + 
      geom_point() + geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      ggtitle(conf.title) +
      theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "Relative Confluence", 
                                         shape = "Chr8q Status", color = "MPNST Cell Line") + theme(plot.title=element_text(face="bold",hjust=0.5))
    ggsave(paste0(d,"_",t,"h_relConfluence_max20um_no_4-23-25_v2.pdf"),conf.plot,width=4,height=3)
    
    death.plot <- ggplot(mean.drug.death, aes(x=DOSE, y=meanGROWTH, shape=chr8q, color=improve_sample_id)) +
      geom_smooth(linetype="dashed",se=FALSE, method=drc::drm, method.args=list(fct=LL.4())) + 
      geom_point() + geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      ggtitle(death.title) +
      theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "Relative Death", 
                                         shape = "Chr8q Status", color = "MPNST Cell Line") + 
      theme(plot.title=element_text(face="bold",hjust=0.5))
    ggsave(paste0(d,"_",t,"h_relDeath_max20um_no_4-23-25_v2.pdf"),death.plot,width=4,height=3)
    
    conf.plot <- ggplot(mean.drug.conf, aes(x=DOSE, y=meanGROWTH, shape=chr8q, color=improve_sample_id)) +
      geom_smooth(linetype="dashed",se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) + 
      scale_x_continuous(trans="log10") + geom_point() + 
      geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      ggtitle(conf.title) + theme_classic(base_size=12) + 
      labs(x="Concentration (uM)", y = "Relative Confluence", shape = "Chr8q Status", color = "MPNST Cell Line") + 
      theme(plot.title=element_text(face="bold",hjust=0.5))
    ggsave(paste0(d,"_",t,"h_relConfluence_max20um_no_4-23-25_v2_log10.pdf"),conf.plot,width=4,height=3)
    
    death.plot <- ggplot(mean.drug.death, aes(x=DOSE, y=meanGROWTH, shape=chr8q, color=improve_sample_id)) +
      geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) + 
      scale_x_continuous(trans="log10") + geom_point() + geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      ggtitle(death.title) + theme_classic(base_size=12) + 
      labs(x="Concentration (uM)", y = "Relative Death", shape = "Chr8q Status", color = "MPNST Cell Line") + 
      theme(plot.title=element_text(face="bold",hjust=0.5))
    ggsave(paste0(d,"_",t,"h_relDeath_max20um_no_4-23-25_v2_log10.pdf"),death.plot,width=4,height=3)
    
    ## try filtering to remove very high concentrations
    conc.max <- c(0.2,0.5,1,1.25,2.5,5)
    for (c in conc.max) {
      p.df <- data.frame(Metric, moreSensitive = NA, p = NA, Drug=d, Time=t, maxConcUM=c, minValue=NA)
      # run t-tests between WU-356 (chr8q-gain) and JH-2-055d (chr8q-diploid)
      conf.title <- chr8_tTest(drug.conf[drug.conf$DOSE <= c,], "Relative Confluence")
      conf.more.sens <- strsplit(conf.title, " ")[[1]][1]
      conf.p <- strsplit(conf.title, "p=")[[1]][2]
      conf.p <- as.numeric(substr(conf.p,1,nchar(conf.p)-1))
      p.df[p.df$Metric=="Relative Confluence",]$moreSensitive <- conf.more.sens
      p.df[p.df$Metric=="Relative Confluence",]$p <- conf.p
      p.df[p.df$Metric=="Relative Confluence",]$minValue <- min(drug.conf[drug.conf$DOSE <= c,]$GROWTH, na.rm=TRUE)
      
      death.title <- chr8_tTest(drug.death[drug.death$DOSE <= c,], "Relative Death")
      death.more.sens <- strsplit(death.title, " ")[[1]][1]
      death.p <- strsplit(death.title, "p=")[[1]][2]
      death.p <- as.numeric(substr(death.p,1,nchar(death.p)-1))
      p.df[p.df$Metric=="Relative Death",]$moreSensitive <- death.more.sens
      p.df[p.df$Metric=="Relative Death",]$p <- death.p
      p.df[p.df$Metric=="Relative Death",]$minValue <- min(drug.death[drug.death$DOSE <= c,]$GROWTH, na.rm=TRUE)
      all.p.df <- rbind(all.p.df, p.df)
      
      # generate plots
      conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,], aes(x=DOSE, y=meanGROWTH, shape=chr8q, color=improve_sample_id)) +
        geom_smooth(linetype="dashed",se=FALSE, method=drc::drm, method.args=list(fct=LL.4())) + 
        geom_point() + geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
        ggtitle(conf.title) +
        theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "Relative Confluence", 
                                           shape = "Chr8q Status", color = "MPNST Cell Line") + theme(plot.title=element_text(face="bold",hjust=0.5))
      ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_v2.pdf"),conf.plot,width=4,height=3)
      
      death.plot <- ggplot(mean.drug.death[mean.drug.death$DOSE <= c,], aes(x=DOSE, y=meanGROWTH, shape=chr8q, color=improve_sample_id)) +
        geom_smooth(linetype="dashed",se=FALSE, method=drc::drm, method.args=list(fct=LL.4())) + 
        geom_point() + geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
        ggtitle(death.title) +
        theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "Relative Death", 
                                           shape = "Chr8q Status", color = "MPNST Cell Line") + 
        theme(plot.title=element_text(face="bold",hjust=0.5))
      ggsave(paste0(d,"_",t,"h_relDeath_max",c,"um_no_4-23-25_v2.pdf"),death.plot,width=4,height=3) 
      
      conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,], 
                          aes(x=DOSE, y=meanGROWTH, shape=chr8q, color=improve_sample_id)) +
        geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) + 
        scale_x_continuous(transform="log10") + geom_point() + 
        geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
        ggtitle(conf.title) +
        theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "Relative Confluence", 
                                           shape = "Chr8q Status", color = "MPNST Cell Line") + theme(plot.title=element_text(face="bold",hjust=0.5))
      ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_v2_log10.pdf"),conf.plot,width=4,height=3)
      
      death.plot <- ggplot(mean.drug.death[mean.drug.death$DOSE <= c,], 
                           aes(x=DOSE, y=meanGROWTH, shape=chr8q, color=improve_sample_id)) +
        geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
        scale_x_continuous(transform="log10") + geom_point() + 
        geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) + 
        ggtitle(death.title) + theme_classic(base_size=12) + 
        labs(x="Concentration (uM)", y = "Relative Death", shape = "Chr8q Status", color = "MPNST Cell Line") + 
        theme(plot.title=element_text(face="bold",hjust=0.5))
      ggsave(paste0(d,"_",t,"h_relDeath_max",c,"um_no_4-23-25_v2_log10.pdf"),death.plot,width=4,height=3) 
    }
  }
}
write.csv(all.p.df, "p_values_relConfluenceORDeath.csv", row.names=FALSE)

#### also plot with JH data ####
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
rel.conf <- read.table("WU/chr8_relConfluence_max20um_no_4-23-25.tsv", sep="\t", header=TRUE)
rel.conf$GROWTH <- 100*rel.conf$GROWTH
#jh.rel.conf <- read.table("JH/chr8_normConfluence_max1um.tsv", sep="\t", header=TRUE)
jh.rel.conf <- read.table("JH/chr8_normConfluence_max10um.tsv", sep="\t", header=TRUE)
rel.conf <- rbind(rel.conf, jh.rel.conf)

# calculate mean and sd for each drug, dose, mpnst, time combo; also preserve sample and chr8q columns
library(plyr);library(dplyr)
mean.conf <- plyr::ddply(rel.conf, .(improve_sample_id, time, Drug, DOSE), summarize,
                         meanGROWTH = mean(GROWTH, na.rm=TRUE),
                         sdGROWTH = sd(GROWTH, na.rm=TRUE))

# annotate chr8q status
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

#times <- unique(c(rel.conf$time, rel.death$time)) # too many
times <- c(24,48,72,96,120)
all.mpnst <- unique(mean.conf$improve_sample_id)
all.mpnst <- c(known.dip,known.amp, all.mpnst[!(all.mpnst %in% c(known.dip,known.amp))])
dir.create("singleTimePlots_20250714")
setwd("singleTimePlots_20250714")
library(drc) # curve source: answer by greenjune: https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc
# Sara shared these 2 links: https://stackoverflow.com/questions/68209998/plot-drc-using-ggplot; https://forum.posit.co/t/extract-out-points-in-dose-response-curve-produced-by-drc/159433
Metric <- c("% Relative Viability")
library(ggplot2)
all.p.df <- data.frame()
for (t in times) {
  rel.conf.96 <- rel.conf[rel.conf$time==t,]
  mean.conf.96 <- mean.conf[mean.conf$time==t,]
  drugs <- unique(mean.conf.96$Drug) # all 6 drugs are now present: AT7519, Mirdametinib, Gefitinib, Ribociclib, Alrizo, Volas
  mpnst <- unique(mean.conf.96$improve_sample_id) # all 2 cell lines: JH-2-055d and WU-356
  
  # plot relative confluence and relative death for each drug
  for (d in drugs) {
    # filter for drug
    drug.conf <- rel.conf.96[rel.conf.96$Drug == d,]
    mean.drug.conf <- mean.conf.96[mean.conf.96$Drug == d,]
    
    ## try filtering to remove very high concentrations
    #conc.max <- c(1,1.25,2.5,5,10,20)
    #conc.max <- conc.max[conc.max %in% drug.conf$DOSE]
    conc.max <- max(drug.conf$DOSE)
    for (c in conc.max) {
      p.df <- data.frame(Metric, chr8_moreSensitive = NA, MYC_moreSensitive = NA,
                         chr8_p = NA, MYC_p=NA, Drug=d, Time=t, maxConcUM=c, minValue=NA)
      # run t-tests between WU-356 (chr8q-gain) and JH-2-055d (chr8q-diploid)
      titles <- chr8_tTest(na.omit(drug.conf[drug.conf$DOSE <= c,]), "% Relative Viability")
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
      
      # generate plots
      # conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c & mean.drug.conf$chr8q != "Not yet known",], 
      #                     aes(x=DOSE, y=meanGROWTH, shape=chr8q, color=improve_sample_id)) +
      #   geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) + 
      #   scale_x_continuous(transform="log10") + geom_point() + 
      #   geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      #   ggtitle(conf.title) +
      #   theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability", 
      #                                      shape = "Chr8q Status", color = "MPNST Cell Line") + theme(plot.title=element_text(face="bold",hjust=0.5))
      # ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8q_log10.pdf"),conf.plot,width=5,height=4)
      # 
      # conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,],
      #                     aes(x=DOSE, y=meanGROWTH, color=chr8q, shape=improve_sample_id)) +
      #   geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
      #   scale_color_manual(values=c(scales::hue_pal()(2),"lightgrey")) + 
      #   scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
      #   scale_x_continuous(transform="log10") + geom_point() +
      #   geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      #   ggtitle(conf.title) +
      #   theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
      #                                      color = "Chr8q Status", shape = "MPNST Cell Line") + theme(plot.title=element_text(face="bold",hjust=0.5))
      # ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8q_log10_alpha.pdf"),conf.plot,width=5,height=4)
      # 
      # conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,],
      #                     aes(x=DOSE, y=meanGROWTH, color=chr8q, shape=improve_sample_id)) +
      #   geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
      #   scale_color_manual(values=c("black","grey"), breaks=c("Diploid","Gain")) + 
      #   scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
      #   scale_x_continuous(transform="log10") + geom_point() +
      #   geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      #   ggtitle(conf.title) +
      #   theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
      #                                      color = "Chr8q Status", shape = "MPNST Cell Line") + theme(plot.title=element_text(face="bold",hjust=0.5))
      # ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8q_log10_v2.pdf"),conf.plot,width=4,height=2.5)
      # 
      # conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,],
      #                     aes(x=DOSE, y=meanGROWTH, color=chr8q, shape=improve_sample_id)) +
      #   geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
      #   scale_color_manual(values=c("black","grey"), breaks=c("Diploid","Gain")) + 
      #   scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
      #   scale_x_continuous(transform="log10") + geom_point() +
      #   geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      #   ggtitle(titles$chr8) +
      #   theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
      #                                      color = "Chr8q Status", shape = "MPNST Cell Line") + 
      #   theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
      # ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8q_log10_v2_angled.pdf"),conf.plot,width=4,height=2.5)
      # 
      # conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c & mean.drug.conf$improve_sample_id %in% c("JH-2-055d","JH-2-002"),],
      #                     aes(x=DOSE, y=meanGROWTH, color=chr8q, shape=improve_sample_id)) +
      #   geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
      #   scale_color_manual(values=c("blue","red"), breaks=c("Diploid","Gain")) + 
      #   scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
      #   scale_x_continuous(transform="log10") + geom_point() +
      #   geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      #   ggtitle(titles$chr8) +
      #   theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
      #                                      color = "Chr8q Status", shape = "MPNST Cell Line") + 
      #   theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
      # ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8q_log10_v2_angled_JH5donly.pdf"),conf.plot,width=4,height=2.5)
      
      conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c & mean.drug.conf$improve_sample_id %in% c("JH-2-055d","JH-2-002","ST8814","NF10.1"),],
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
      
      conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c & mean.drug.conf$improve_sample_id %in% c("JH-2-055d","JH-2-002","ST8814","NF10.1"),],
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
      
      conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c & mean.drug.conf$improve_sample_id %in% c("JH-2-055d","JH-2-002","ST8814","NF10.1"),],
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
      
      conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c & mean.drug.conf$improve_sample_id %in% c("JH-2-055d","JH-2-002","ST8814","NF10.1"),],
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
      
      # conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,],
      #                     aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
      #   geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
      #   scale_color_manual(values=c("blue","black","red"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + # was black, darkgrey, lightgrey
      #   scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
      #   scale_x_continuous(transform="log10") + geom_point() +
      #   geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      #   ggtitle(titles$chr8) +
      #   theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
      #                                      color = "Chr8q Status", shape = "MPNST Cell Line") + theme(plot.title=element_text(face="bold",hjust=0.5))
      # ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8qMYC_log10_v2.pdf"),conf.plot,width=4,height=2.5)
      
      conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,],
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
      
      conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,],
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
      
      conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,],
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
      
      conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,],
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

      # conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,], 
      #                     aes(x=DOSE, y=meanGROWTH, shape=chr8q, color=improve_sample_id, alpha=alpha)) +
      #   geom_smooth(linetype="dashed", se=FALSE, alpha=mean.drug.conf[mean.drug.conf$DOSE <= c,]$alpha, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) + 
      #   scale_alpha_identity() + scale_x_continuous(transform="log10") + geom_point() + 
      #   geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      #   ggtitle(conf.title) +
      #   theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability", 
      #                                      shape = "Chr8q Status", color = "MPNST Cell Line") + theme(plot.title=element_text(face="bold",hjust=0.5))
      # ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_knownChr8q_log10_alpha_v2.pdf"),conf.plot,width=5,height=4)
      # 
      # 
      # conf.plot <- ggplot(mean.drug.conf[mean.drug.conf$DOSE <= c,], 
      #                     aes(x=DOSE, y=meanGROWTH, shape=chr8q, color=improve_sample_id)) +
      #   geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) + 
      #   scale_x_continuous(transform="log10") + geom_point() + 
      #   geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
      #   ggtitle(conf.title) +
      #   theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability", 
      #                                      shape = "Chr8q Status", color = "MPNST Cell Line") + theme(plot.title=element_text(face="bold",hjust=0.5))
      # ggsave(paste0(d,"_",t,"h_relConfluence_max",c,"um_no_4-23-25_v2_log10.pdf"),conf.plot,width=5,height=4)
    }
  }
}
write.csv(all.p.df, "p_values_relPercViability_JHandWU.csv", row.names=FALSE)
all.p.long <- reshape2::melt(all.p.df[all.p.df$Time==120,])
all.p.longer <- reshape2::melt(all.p.long, id.var=c("Metric","Drug","variable","value"))
colnames(all.p.longer) <- c("Metric", "Drug", "Test", "p", "whichSens", "moreSensitive")
all.p.longer <- all.p.longer[,c("Drug","Test","p","moreSensitive")]
all.p.longer$Test <- sub("_p","",all.p.longer$Test)
all.p.longer <- dplyr::distinct(na.omit(all.p.longer[all.p.longer$Test %in% c("chr8","MYC"),]))
ggplot(all.p.longer, aes(x=reorder(Drug,log10(p)), y=-log10(p), fill=moreSensitive)) + 
  geom_bar(stat="identity", position="dodge") + facet_wrap(.~Test) + theme_classic() + 
  scale_fill_manual(values=c("red","blue"), breaks=c("gain","diploid"), labels=c("Gain","Diploid")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray") + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
                                                  axis.title.x=element_blank(), legend.position="top") + 
  labs(fill="More Sensitive",y="-log(P-value)")
ggsave("pValues_v2.pdf", width=3.5, height=3) # was width 5

ggplot(all.p.longer[all.p.longer$Drug != "AT7519",], aes(x=reorder(Drug,log10(p)), y=-log10(p), fill=moreSensitive)) + 
  geom_bar(stat="identity", position="dodge") + facet_wrap(.~Test) + theme_classic() + 
  scale_fill_manual(values=c("red","blue"), breaks=c("gain","diploid"), labels=c("Gain","Diploid")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray") + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.title.x=element_blank(), legend.position="top") + 
  labs(fill="More Sensitive",y="-log(P-value)")
ggsave("pValues_noAT7519_v2.pdf", width=3, height=3) # was width 4

ggplot(all.p.longer[!(all.p.longer$Drug %in% c("AT7519","BI-2536")),], aes(x=reorder(Drug,log10(p)), y=-log10(p), fill=moreSensitive)) + 
  geom_bar(stat="identity", position="dodge") + facet_wrap(.~Test) + theme_classic() + 
  scale_fill_manual(values=c("red","blue"), breaks=c("gain","diploid"), labels=c("Gain","Diploid")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray") + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.title.x=element_blank()) + 
  labs(fill="More Sensitive",y="-log(P-value)")
ggsave("pValues_noAT7519-BI2536.pdf", width=3.5, height=3)

ggplot(all.p.longer[!(all.p.longer$Drug %in% c("AT7519","BI-2536")),], aes(x=reorder(Drug,log10(p)), y=-log10(p), fill=moreSensitive)) + 
  geom_bar(stat="identity", position="dodge") + facet_wrap(.~Test) + theme_classic() + 
  scale_fill_manual(values=c("red","blue"), breaks=c("gain","diploid"), labels=c("Gain","Diploid")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray") + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.title.x=element_blank(), legend.position="top") + 
  labs(fill="More Sensitive",y="-log(P-value)")
ggsave("pValues_noAT7519-BI2536_v2.pdf", width=3, height=3)

ggplot(all.p.longer[!(all.p.longer$Drug %in% c("AT7519","Mirdametinib")),], aes(x=reorder(Drug,log10(p)), y=-log10(p), fill=moreSensitive)) + 
  geom_bar(stat="identity", position="dodge") + facet_wrap(.~Test) + theme_classic() + 
  scale_fill_manual(values=c("red","blue"), breaks=c("gain","diploid"), labels=c("Gain","Diploid")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray") + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.title.x=element_blank()) + 
  labs(fill="More Sensitive",y="-log(P-value)")
ggsave("pValues_noAT7519-Mirdametinib.pdf", width=3.5, height=3)

ggplot(all.p.longer[!(all.p.longer$Drug %in% c("AT7519","Mirdametinib")),], aes(x=reorder(Drug,log10(p)), y=-log10(p), fill=moreSensitive)) + 
  geom_bar(stat="identity", position="dodge") + facet_wrap(.~Test) + theme_classic() + 
  scale_fill_manual(values=c("red","blue"), breaks=c("gain","diploid"), labels=c("Gain","Diploid")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray") + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.title.x=element_blank(), legend.position="top") + 
  labs(fill="More Sensitive",y="-log(P-value)")
ggsave("pValues_noAT7519-Mirdametinib_v2.pdf", width=3, height=3)

ggplot(all.p.longer[!(all.p.longer$Drug %in% c("AT7519","BI-2536","Mirdametinib")),], aes(x=reorder(Drug,log10(p)), y=-log10(p), fill=moreSensitive)) + 
  geom_bar(stat="identity", position="dodge") + facet_wrap(.~Test) + theme_classic() + 
  scale_fill_manual(values=c("red","blue"), breaks=c("gain","diploid"), labels=c("Gain","Diploid")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray") + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.title.x=element_blank(), legend.position="top") + 
  labs(fill="More Sensitive",y="-log(P-value)")
ggsave("pValues_EGFR-PLKonly.pdf", width=2, height=3)

ggplot(all.p.longer[all.p.longer$Drug %in% c("Erlotinib","Volasertib"),], aes(x=reorder(Drug,log10(p)), y=-log10(p), fill=moreSensitive)) + 
  geom_bar(stat="identity", position="dodge") + facet_wrap(.~Test) + theme_classic() + 
  scale_fill_manual(values=c("red","blue"), breaks=c("gain","diploid"), labels=c("Gain","Diploid")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray") + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.title.x=element_blank()) + 
  labs(fill="More Sensitive",y="-log(P-value)")
ggsave("pValues_erlotinibVolasertibOnly.pdf", width=3, height=3)

ggplot(all.p.longer[all.p.longer$Drug %in% c("Erlotinib","Volasertib"),], aes(x=reorder(Drug,log10(p)), y=-log10(p), fill=moreSensitive)) + 
  geom_bar(stat="identity", position="dodge") + facet_wrap(.~Test) + theme_classic() + 
  scale_fill_manual(values=c("red","blue"), breaks=c("gain","diploid"), labels=c("Gain","Diploid")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray") + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.title.x=element_blank(), legend.position="top") + 
  labs(fill="More Sensitive",y="-log(P-value)")
ggsave("pValues_erlotinibVolasertibOnly_v2.pdf", width=1.5, height=3)


#### compare auc to FISH copy numbers ####
auc <- read.table("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/drugScreens/JH/chr8_normConfluenceCurves_max1um.0", 
                  fill=TRUE, header=TRUE)
auc2 <- auc[auc$improve_drug_id %in% c("Erlotinib","Volasertib"),]



all.p.df[all.p.df$Drug == "Volas",]$Drug <- "Volasertib"
all.p.df$DrugTime <- paste0(all.p.df$Drug, " (", as.numeric(all.p.df$Time)/24,"d)")

# just volasertib, mirdametinib, gefitinib, erlotinib, osimertinib 4-5d
plot.df <- mean.conf[mean.conf$time >= 96 & mean.conf$Drug %in% c("Volas","Volasertib","Mirdametinib","Gefitinib","Erlotinib","Osimertinib"),]
plot.df[plot.df$Drug == "Volas",]$Drug <- "Volasertib"
plot.df$DrugTime <- paste0(plot.df$Drug, " (", as.numeric(plot.df$time)/24,"d)")
doi <- c("Volasertib (4d)", "Volasertib (5d)","Mirdametinib (5d)",
         "Gefitinib (4d)", "Erlotinib (5d)", "Osimertinib (5d)")
doi.p <- all.p.df[all.p.df$DrugTime %in% doi,]
p0.05 <- unique(doi.p[doi.p$p < 0.05,]$DrugTime)
p0.01 <- unique(doi.p[doi.p$p < 0.01,]$DrugTime)
p0.001 <- unique(doi.p[doi.p$p < 0.001,]$DrugTime)
p0.0001 <- unique(doi.p[doi.p$p < 0.0001,]$DrugTime)
plot.df <- plot.df[plot.df$DrugTime %in% doi,]
plot.df$DrugTimeSig <- plot.df$DrugTime
plot.df[plot.df$DrugTime %in% p0.05,]$DrugTimeSig <- paste0(plot.df[plot.df$DrugTime %in% p0.05,]$DrugTime, "*")
plot.df[plot.df$DrugTime %in% p0.01,]$DrugTimeSig <- paste0(plot.df[plot.df$DrugTime %in% p0.01,]$DrugTime, "**")
plot.df[plot.df$DrugTime %in% p0.001,]$DrugTimeSig <- paste0(plot.df[plot.df$DrugTime %in% p0.001,]$DrugTime, "***")

# plot.df$DrugTime <- factor(plot.df$DrugTime, levels=c("Volasertib (4d)", "Volasertib (5d)","Mirdametinib (5d)",
#                                                       "Gefitinib (4d)", "Erlotinib (5d)", "Osimertinib (5d)"))
plot.df$DrugTimeSig <- factor(plot.df$DrugTimeSig, levels=c("Volasertib (4d)***", "Volasertib (5d)**","Mirdametinib (5d)*",
                                                      "Gefitinib (4d)***", "Erlotinib (5d)**", "Osimertinib (5d)"))

conf.plot <- ggplot(plot.df[plot.df$DrugTime %in% c("Volasertib (4d)", "Volasertib (5d)","Mirdametinib (5d)",
                                                    "Gefitinib (4d)", "Erlotinib (5d)", "Osimertinib (5d)"),],
                    aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
  geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_color_manual(values=c("black","darkgrey","lightgrey"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + 
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
  scale_x_continuous(transform="log10") + geom_point() +
  geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
  facet_wrap(.~DrugTime) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                     color = "Chr8q Status", shape = "MPNST Cell Line") + 
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v2_angled.pdf",conf.plot,width=7,height=4)

conf.plot <- ggplot(plot.df,
                    aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
  geom_smooth(linewidth=0.5, linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_color_manual(values=c("black","gold2","lightgrey"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + 
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
  scale_x_continuous(transform="log10") + geom_point() +
  geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
  facet_wrap(.~DrugTimeSig) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                     color = "Chr8q Status", shape = "MPNST Cell Line") + 
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v2_angled_diffColors.pdf",conf.plot,width=7,height=4)
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v2_angled_diffColors_wider.pdf",conf.plot,width=8,height=4)

conf.plot <- ggplot(plot.df[plot.df$DrugTime %in% c("Volasertib (4d)", "Volasertib (5d)","Mirdametinib (5d)",
                                                    "Gefitinib (4d)", "Erlotinib (5d)", "Osimertinib (5d)"),],
                    aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id, linewidth=chr8qv2)) +
  geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_color_manual(values=c("black","darkgrey","lightgrey"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + 
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
  scale_x_continuous(transform="log10") + geom_point() +
  geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
  facet_wrap(.~DrugTime) + scale_linewidth_manual(values=c(1,0.5,0.25),breaks=c("Diploid","Diploid with MYC gain","Gain"))+
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                     color = "Chr8q Status", shape = "MPNST Cell Line") + 
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v3_angled.pdf",conf.plot,width=7,height=4)
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v2_angled_taller.pdf",conf.plot,width=7,height=6)

conf.plot <- ggplot(plot.df[plot.df$DrugTime %in% c("Volasertib (4d)", "Volasertib (5d)","Mirdametinib (5d)",
                                                    "Gefitinib (4d)", "Erlotinib (5d)", "Osimertinib (5d)"),],
                    aes(x=DOSE, y=meanGROWTH, color=chr8q, shape=improve_sample_id, linetype=MYC)) +
  geom_smooth(se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_color_manual(values=c("black","grey"), breaks=c("Diploid","Gain")) +
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
  scale_x_continuous(transform="log10") + geom_point() +
  geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
  facet_wrap(.~DrugTime) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                     color = "Chr8q Status", shape = "MPNST Cell Line") + 
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v4_angled.pdf",conf.plot,width=7,height=4)

conf.plot <- ggplot(plot.df[plot.df$DrugTime %in% c("Volasertib (4d)", "Volasertib (5d)","Mirdametinib (5d)",
                                                    "Gefitinib (4d)", "Erlotinib (5d)", "Osimertinib (5d)"),],
                    aes(x=DOSE, y=meanGROWTH, color=chr8q, shape=improve_sample_id, linetype=MYC)) +
  geom_smooth(linewidth=0.5, se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_color_manual(values=c("black","grey"), breaks=c("Diploid","Gain")) +
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
  scale_x_continuous(transform="log10") + geom_point() +
  geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
  facet_wrap(.~DrugTime) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                     color = "Chr8q Status", shape = "MPNST Cell Line") + 
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position="bottom", legend.direction="horizontal")
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v5_angled.pdf",conf.plot,width=7,height=4)

conf.plot <- ggplot(plot.df[plot.df$DrugTime %in% c("Volasertib (4d)", "Volasertib (5d)","Mirdametinib (5d)",
                                                    "Gefitinib (4d)", "Erlotinib (5d)", "Osimertinib (5d)"),],
                    aes(x=DOSE, y=meanGROWTH, color=chr8q, shape=improve_sample_id, linetype=MYC)) +
  geom_smooth(se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_color_manual(values=c("black","grey"), breaks=c("Diploid","Gain")) +
  scale_linetype_manual(values=c("solid","dashed"), breaks=c("FALSE","TRUE"), labels=c("Diploid","Gain")) +
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
  scale_x_continuous(transform="log10") + geom_point() +
  geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
  facet_wrap(.~DrugTime) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                     color = "Chr8q Status", shape = "MPNST Cell Line") + 
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        #legend.position="bottom", legend.direction="horizontal"
        )
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v5_angled.pdf",conf.plot,width=7,height=4)



conf.plot <- ggplot(plot.df[plot.df$DrugTime %in% c("Volasertib (4d)", "Volasertib (5d)","Mirdametinib (5d)",
                                                    "Gefitinib (4d)", "Erlotinib (5d)", "Osimertinib (5d)"),],
                    aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
  geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_color_manual(values=c("black","darkgrey","lightgrey"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + 
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
  scale_x_continuous(transform="log10") + geom_point() +
  geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
  facet_wrap(.~DrugTime) + scale_y_continuous(limits=c(0,100)) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                     color = "Chr8q Status", shape = "MPNST Cell Line") + 
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v2_angled_max100.pdf",conf.plot,width=7,height=4)

conf.plot <- ggplot(plot.df[plot.df$DrugTime %in% c("Volasertib (4d)", "Volasertib (5d)","Mirdametinib (5d)",
                                                    "Gefitinib (4d)", "Erlotinib (5d)", "Osimertinib (5d)"),],
                    aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
  geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_color_manual(values=c("black","darkgrey","lightgrey"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + 
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
  scale_x_continuous(transform="log10") + geom_point() +
  geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
  facet_wrap(.~DrugTime) + scale_y_continuous(limits=c(0,150)) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                     color = "Chr8q Status", shape = "MPNST Cell Line") + 
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v2_angled_max150.pdf",conf.plot,width=7,height=4)

conf.plot <- ggplot(plot.df[plot.df$DrugTime %in% c("Volasertib (4d)", "Volasertib (5d)","Mirdametinib (5d)",
                                                    "Gefitinib (4d)", "Erlotinib (5d)", "Osimertinib (5d)"),],
                    aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
  geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_color_manual(values=c("black","darkgrey","lightgrey"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + 
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
  scale_x_continuous(transform="log10") + geom_point() +
  geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
  facet_wrap(.~DrugTime) + scale_y_continuous(limits=c(0,120)) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                     color = "Chr8q Status", shape = "MPNST Cell Line") + 
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v2_angled_max120.pdf",conf.plot,width=7,height=4)

conf.plot <- ggplot(plot.df[plot.df$DrugTime %in% c("Volasertib (4d)", "Volasertib (5d)","Mirdametinib (5d)",
                                                    "Gefitinib (4d)", "Erlotinib (5d)", "Osimertinib (5d)"),],
                    aes(x=DOSE, y=meanGROWTH, color=chr8qv2, shape=improve_sample_id)) +
  geom_smooth(linetype="dashed", se=FALSE, method=drc::drm, method.args=list(fct=L.4(), se=FALSE)) +
  scale_color_manual(values=c("black","darkgrey","lightgrey"), breaks=c("Diploid","Diploid with MYC gain","Gain")) + 
  scale_shape_manual(values=1:length(all.mpnst), breaks=all.mpnst)+
  scale_x_continuous(transform="log10") + geom_point() +
  geom_errorbar(aes(ymin=meanGROWTH-sdGROWTH, ymax=meanGROWTH+sdGROWTH), width=0.2) +
  facet_wrap(.~DrugTime) + scale_y_continuous(limits=c(0,175)) +
  theme_classic(base_size=12) + labs(x="Concentration (uM)", y = "% Relative Viability",
                                     color = "Chr8q Status", shape = "MPNST Cell Line") + 
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("PLK-MEK-EGFRi_knownChr8qMYC_log10_v2_angled_max175.pdf",conf.plot,width=7,height=4)
