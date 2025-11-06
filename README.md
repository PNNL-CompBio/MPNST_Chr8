# Multi-omics integration of malignant peripheral nerve sheath tumors (MPNST) identifies potential targets based on chromosome 8q (chr8q) status
Since chr8q gain is associated with high-grade transformation in MPNST and 
inferior overall survival, we integrate multi-omics data to understand drivers 
and potential targets based on chr8q status. To do this, we collected new 
proteomics data and ran correlations between omics expression and chr8q copy
number in MPNST patient-derived xenografts (PDX) in addition to pathway 
analyses, network analyses, drug sensitivity predictions, fluorescent in situ
hybridization, and viability studies.

## Initial proteomics processing located in [proteomics](./proteomics)
NOTE: these processing scripts are designed to be run at PNNL, in the internal
network, so will not be able to be run outside. 
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
This folder also contains specific analyses that were run for the paper.
### Dependency installation
installDependencies.R

### Multi-omics correlations and GSEA
5-panSEA_using_helper.R
Uses functions available in panSEA_helper_20240913.R script.

### Simulation study for GSEA adapted to shuffle tied genes/proteins
GSEA_ties_simulation_20240909_v2.R
Also generates plots for figure S2

