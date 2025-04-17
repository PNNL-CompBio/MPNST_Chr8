# prep MIT kinase enrichment input

# load peptide correlations
pep.cor <- read.csv("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Chr8/proteomics/analysis/Chr8_quant_20250409/Phospho_peptide_correlations_Chr8q.csv")

# extract peptide sequences
pep.cor$seq <- sub(".*@", "", pep.cor$Feature)

# remove amino acids and periods before and after sequence
pep.cor$seq <- substr(pep.cor$seq, 3, nchar(pep.cor$seq)-2)

# need at least 5 amino acids before and after phospho-site
pep.cor$nchar <- nchar(pep.cor$seq)
pep.cor$firstPhos <- stringi::stri_locate_first_fixed(pep.cor$seq, "*")[,1]
pep.cor$lastPhos <- stringi::stri_locate_last_fixed(pep.cor$seq, "*")[,1]
pep.cor$nToLastPhos <- pep.cor$nchar - pep.cor$lastPhos
pep.cor$MIT <- FALSE
pep.cor[pep.cor$firstPhos > 5 & pep.cor$nToLastPhos > 4,]$MIT <- TRUE


# filter for foreground (Spearman.q <= 0.05) & background (all peptides)
fore <- pep.cor[pep.cor$Spearman.q <= 0.05 & pep.cor$MIT,]
back <- pep.cor[pep.cor$MIT,] # background must include foreground

# extract peptide sequences
fore.seq <- fore$seq
back.seq <- back$seq

# save
write.table(fore.seq, "Phospho_peptide_foreground.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE)
write.table(back.seq, "Phospho_peptide_background.txt", sep="\t", 
            row.names = FALSE, col.names = FALSE)

# run MIT kinase enrichment: https://kinase-library.mit.edu/ea?a=ps
# serine/threonine and tyrosine kinases
# default settings otherwise