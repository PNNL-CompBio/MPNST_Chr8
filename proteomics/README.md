# Proteomics processing

This folder contains the files needed to process the proteomics from the mass
spectrometry instrument.

## On Prem analyses
These analysis run at the PNNL site and will then upload the files to synapse. 
### Create study design tables
0-create_study_design_tables.Rmd 

### Initial global proteomics processing
1-process_global_data.Rmd

### Initial phospho proteomics processing
2-process_phospho_data.Rmd

### Prep KSTAR input
2.5-process_KSTAR_input.Rmd Not currently used but available in case of future need.

### Normalization and batch correction
3-normalize_and_batch_correction.Rmd

### Upload crosstabs to Synapse
4-push_to_synapse.Rmd

## Analyses for paper
Once the initial analyses are run, this code will store all the data on the local 
hard drive and upload the necessary results to Synapse. 

### Multi-omics correlations and GSEA
5-panSEA_using_helper.R Uses functions available in panSEA_helper_20240913.R script.

### Simulation study for GSEA adapted to shuffle tied genes/proteins
GSEA_ties_simulation_20240909_v2.R 
Also generates plots for figure S2
