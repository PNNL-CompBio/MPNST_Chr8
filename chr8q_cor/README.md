# Multi-omics integration of malignant peripheral nerve sheath tumors (MPNST) identifies potential targets based on chromosome 8q (chr8q) status

These scripts are designed to be run after those in the [proteomics](./proteomics)
directory, where the primary proteomics ws done.  

The analysis for this manuscript integrates copy number and RNASeq measurements
of diverse MPNST together with proteomics measurements. The proteomics analysis
is stores in [the proteomics directory](../proteomics) and stores the files
on Synapse at http://synapse.org/chr8_mpnst. Once you request access to the
data files you can reproduce the analysis in the paper below.

# Basic analysis 

We first calculate the chr8q copy numbers and correlations across molecules.
NOTE: always return to this directory as a working directory after every script,
as there are many subdirectories created that can make things very confusing. 

```
#end to end analysis runs and stores lots of files
source('00_computeChr8Expression.R')
source("01_enrichments_copynumber.R")
source('02_enrichments_protRNA.R')
source('03_DMEA_correlations.R')
source('04_networkAnalysis.R') #this requires Cytoscape to be open to run
source("05_tableCompilation.R")

```

This generates the data and tables needed for the figures. 
Below is a breakdown of the scripts, scroll down for figure instructions


### Compute chr8 copy numbers
The first step to analyzing the data is to calculate the copy number of chr8
across the samples. This is currently run by the 
[00_computeChrCorrelations.R](./00_computeChrCorrelations.R)
script.

This script downloads the relevant data to the environment, creates a folder
with today's date and saves the copy number data. It also creates a folder
with Supplemental tables 1 and 2. 

### Copy chromosome enrichments
We then calculate positional enrichment to ensure that chromosome 8 activity
is indeed altered in the samples in 
[01_enrichments_copynumber.R](./01_enrichments_copynumber.R). This will
save files to the working directory with the current date. 

### Calculate mRNA-Prot correlations
The next step is to compute the correlations between copy number
computed in the previous step and mRNA/protein levels to identify functional pathway
enrichment.  This is calculated by
[02_enrichments_protRNA.R](./02_enrichments_protRNA.R). 

### DMEA analysis
The next step is to carry out the drug correlations used for Figure 3, which
are slightly different than the mechanism enrichment analyses. There's a chance
that the DMEA is run more than once when we run [03_DMEA_correlations.R]('./03_DMEA_correlations.R')

### Network analysis

We then used the stored files to build the networks used in Figures 2 and 3. 
This is [04_networkAnalysis.R](./04_networkAnalysis.R). If you want to use
updated data, then you must update the three synapse IDs that pull the most 
recent data. 

### Compile tables
This requires that you run figure 3 code first, but will collect the files 
in the local directory so that they can be supplemental tables
[05_tableCompilation.R](./05_tableCompilation.R)

# Figure panel generation

Now that the calculations have been made we can generate the figure panels in order
through the following scripts. The scripts assume that the above scripts have been run
and uploaded to synapse. However, if you run the following scripts they will refer
to previously run synapse files (currently from the 2026-01-22 run of the code). 

### Figure 1: correlations
To build the figures in a `figures` directory, first start with 
[Figure1_correlations.R](./Figure1_correlations.R). This requires 5 files from synapse
generated from above, so you will need to update it if you re-ran the analysis


### Figure 2: enrichment analyses
To visualize the enrichment analysis run [Figure2_enrichment_centrality.R](./Figure2_enrichment_centrality.R)


### Figure 3: drug sensitivity predictions
Next run [Figure3_drugCors.R](./Figure3_drugCors.R) Creates the figures for Figure 3.

### Figure S1: histograms
FigureS1_histograms.R

### Figure S3, S4: drivers of enrichment
Kinase substrates are plotted here: FigureS3_kinase.R

Transcription factor targets are plotted here but we don't reference in paper: FigureS4_TF.R

## Upload to synapse

There are many files that might be used later that are uploaded to synapse using
the [06_save_to_synapse.R](./03_save_to_synapse.R) script. It is designed to take
the current working directory and run immediately after the first three scripts.
NOTE: you will require upload access to this repository. 

# Older code
The following files were not validated but also run parts of the analysis. 

### Drug screen
- drugScreening.R: only used for IC50 t-tests
- fit_curve.py: calculates area under the curve values; not currently used; 
adapted from: https://github.com/PNNL-CompBio/coderdata/blob/main/coderbuild/utils/fit_curve.py


