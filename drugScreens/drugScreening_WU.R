# analyzing IncuCyte results
library(plyr); library(dplyr); library(tidyr)
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
maxConc <- c(10,20) # uM; all non-control concentrations: 0.04  0.08  0.16  0.32  0.63  1.25  2.50  5.00 10.00 20.00
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
