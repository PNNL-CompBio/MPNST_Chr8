##compare the IC50 values from the drug test with FISH-measured copy number


library(ggplot2)
library(dplyr)
library(synapser)
library(ggpubr)

synLogin()

#synapser::synLogin()
#base.path = '2026-10-05'
if(!exists('base.path'))
  base.path <- file.path(getwd(),Sys.Date())

fig.path <- file.path(base.path,'figures')
if (!dir.exists(fig.path))
  dir.create(fig.path)

##get ic50 values
ic50s <- readxl::read_xlsx(synapser::synGet('syn70365118')$path, sheet = 'IC50')

drugs <- colnames(ic50s)[c(2,9,17,24)]

cvals <- ic50s[2:9,grep('Cell line',ic50s[1,])]
mvals <- ic50s[2:9,grep('Mean',ic50s[1,])]

full.res <- do.call(rbind,lapply(1:length(drugs), function(i)
  data.frame(drug = rep(drugs[i],nrow(mvals)),
             `Cell line` = unlist(cvals[,i]),
             ic50 = as.numeric(unlist(mvals[,i])))))
colnames(full.res)[2] <- 'Cell line'

##get FISH results
fish <- readxl::read_xlsx(synapser::synGet('syn70284125')$path, 
                          sheet = 'Raw data for plots')

fish[,2]<-as.numeric(unlist(fish[,2]))

tab <- full.res |> left_join(fish) |>
  tidyr::pivot_longer(c(4:5), names_to='Fish measurement', values_to='Copy number')

ggplot(tab, aes(x=log10(ic50), y = `Copy number`,col = drug,group=drug)) + 
  geom_point() + facet_grid(drug~`Fish measurement`) +   
  geom_smooth(method = "lm", se = FALSE)  +# Adds the linear regression line
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") + 
  theme_minimal()

ggsave(paste0(fig.path,'/figS3_ic50cor.pdf'), width =11)