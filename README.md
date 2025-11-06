# Multi-omics integration of malignant peripheral nerve sheath tumors (MPNST) identifies potential targets based on chromosome 8q (chr8q) status
Since chr8q gain is associated with high-grade transformation in MPNST and 
inferior overall survival, we integrate multi-omics data to understand drivers 
and potential targets based on chr8q status. To do this, we collected new 
proteomics data and ran correlations between omics expression and chr8q copy
number in MPNST patient-derived xenografts (PDX) in addition to pathway 
analyses, network analyses, drug sensitivity predictions, fluorescent in situ
hybridization, and viability studies.

## Install dependencies
For this analysis there are numerous dependencies, thy are listed in 
the [installDependences.R](./installDependencies.R) script. 

## Initial proteomics processing located in [proteomics](./proteomics)
here we have the proteomics processing files. While many of them are designed
to be run internally at PNNL, they can be used to get a sense of how
the Mass spectrometry reads were mapped to proteins and phosphosites. 

## Chr8q correlation analysis located in [figures](./figures)
This directory contains the tools needed to regenerate the figures for the 
Garana et al. manuscript.

## FAK Analysis located int [FAK Analysis](./FAK analysis)
THis directory contains additional analysis studying the role of FAK in
chromosome 8 mediated activity in MPNST. 

