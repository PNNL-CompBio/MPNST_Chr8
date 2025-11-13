# Multi-omics integration of malignant peripheral nerve sheath tumors (MPNST) identifies potential targets based on chromosome 8q (chr8q) status
Since chr8q gain is associated with high-grade transformation in MPNST and 
inferior overall survival, we integrate multi-omics data to understand drivers 
and potential targets based on chr8q status. To do this, we collected new 
proteomics data and ran correlations between omics expression and chr8q copy
number in MPNST patient-derived xenografts (PDX) in addition to pathway 
analyses, network analyses, drug sensitivity predictions, fluorescent in situ
hybridization, and viability studies.

## Install dependencies
For this analysis there are numerous dependencies which are listed in 
the [installDependences.R](./installDependencies.R) script. 

## Initial proteomics processing located in [proteomics](./proteomics)
These proteomics processing files are used to generate dataframes of relative protein abundance. 
While many of them are designed to be run internally at PNNL, they can be used to get a sense of how
the mass spectrometry reads were mapped to proteins and phosphosites. 

### Create study design tables
0-create_study_design_tables.Rmd

### Initial global proteomics processing
1-process_global_data.Rmd

### Initial phospho proteomics processing
2-process_phospho_data.Rmd

### Prep KSTAR input
2.5-process_KSTAR_input.Rmd
Not currently used but available in case of future need.

### Normalization and batch correction
3-normalize_and_batch_correction.Rmd

### Upload crosstabs to Synapse
4-push_to_synapse.Rmd

## Analysis located in [proteomics](./proteomics)
This directory contains the scripts for multi-omics analysis of chr8q 
gain in MPNST for the Garana et al. manuscript.

### Multi-omics correlations and GSEA
5-panSEA_using_helper.R
Uses functions available in panSEA_helper_20240913.R script.

### Simulation study for GSEA adapted to shuffle tied genes/proteins
GSEA_ties_simulation_20240909_v2.R
Also generates plots for figure S2

## Figure generation located in [figures](./figures)
This directory contains the tools needed to regenerate the figures for the 
Garana et al. manuscript.
### Figure 1: correlations
Figure1_20250410.R

### Figure S1: histograms
FigureS1_20250409.R

### Figure 2: enrichment analyses
- Figure2_20250429.R
- networkAnalysis.R

### Figure S3: drivers of enrichment
- Transcription factor targets: FigureS3_TF_20250418.R
- Kinase substrates: FigureS3_kinase_20250418.R

### Figure 3: drug sensitivity predictions
Figure3_20250513.R

### Drug screen
- drugScreening.R: only used for IC50 t-tests
- fit_curve.py: calculates area under the curve values; not currently used; 
adapted from: https://github.com/PNNL-CompBio/coderdata/blob/main/coderbuild/utils/fit_curve.py

## FAK analysis located in [Chr8_FAK.R](./proteomics/Chr8_FAK.R)
This directory contains additional analysis studying the role of FAK in
chromosome 8 mediated activity in MPNST. This analysis is not included in the 
multi-omics study but rather in a separate chr8 study focused on FAK.
