# Proteomics analysis 

This folder contains the anlaysis for the proteomics analysis for this project.

## Proteomics raw data
The first few files require access to PNNL, where the data was processed. The
files are created and then uploaded to synapse.

### Create study design tables
0-create_study_design_tables.Rmd

### Initial global proteomics processing
1-process_global_data.Rmd

### Initial phospho proteomics processing
2-process_phospho_data.Rmd

### Prep KSTAR input
2.5-process_KSTAR_input.Rmd
Not currently used but available in case of future need. 

## Downstream analysis 

This includes normalization and upload. 

TODO: update these files to work with synapse

### Normalization and batch correction
3-normalize_and_batch_correction.Rmd

### Upload crosstabs to Synapse
4-push_to_synapse.Rmd

### Multi-omics correlations and GSEA
5-panSEA_using_helper.R
Uses functions available in panSEA_helper_20240913.R script.

### Simulation study for GSEA adapted to shuffle tied genes/proteins
GSEA_ties_simulation_20240909_v2.R
Also generates plots for figure S2
