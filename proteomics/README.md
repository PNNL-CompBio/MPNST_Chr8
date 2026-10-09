# Proteomics processing of MPNST

This folder contains the files needed to process the proteomics from the mass
spectrometry instrument. It cannot be run from outside the PNNL intranet. 

## Create study design
These analysis run at the PNNL site and will then upload the files to synapse. 
### Create study design tables
0-create_study_design_tables.Rmd 

To map data across TMT plexes, we must first create a study design file to
be used by subsequent analysis steps. The study design for this study is
included in the [0-create_study_design_tables.Rmd](0-create_study_design_tables.Rmd) file. 

## Initial processing
The next steps are to align the peptides to known protein sequences and collect
the quantities across proteins and phosphosites. The global data is processed
in the [1-process_global_data.Rmd](2-process_phospho_data.Rmd) markdown and the 
phosphoproteomic data is processed in the 
[2-process_phospho_data.Rmd](2-process_phospho_data.Rmd) markdown.


## Normalization and batch correction
3-normalize_and_batch_correction.Rmd

## Upload crosstabs to Synapse
4-push_to_synapse.Rmd

## Additional analyses

### Multi-omics correlations and GSEA
5-panSEA_using_helper.R Uses functions available in panSEA_helper_20240913.R script.

### Simulation study for GSEA adapted to shuffle tied genes/proteins
GSEA_ties_simulation_20240909_v2.R 
Also generates plots for figure S2
