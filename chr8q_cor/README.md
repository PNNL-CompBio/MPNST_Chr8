# Multi-omics integration of malignant peripheral nerve sheath tumors (MPNST) identifies potential targets based on chromosome 8q (chr8q) status

These scripts are designed to be run after those in the [proteomics](./proteomics) directory, 
which runs the correlation and enrichment statistics and uploads them to synapse. 

The analysis for this manuscript integrates copy number and RNASeq measurements
of diverse MPNST together with proteomics measurements. The proteomics analysis
is stores in [the proteomics directory](../proteomics) and stores the files
on Synapse at http://synapse.org/chr8_mpnst. Once you request access to the
data files you can reproduce the analysis in the paper below.

# Basic analysis 
First we carry out the analysis and upload to synapse

## Compute chr8 copy numbers
The first step to analyzing the data is to calculate the copy number of chr8
across the samples. This is currently run by the [00_computeChrCorrelations.R](./00_computeChrCorrelations.R)
script.

This script downloads the relevant data to the environment, creates a folder
with todays date and saves the copy number data. 

## Copy chromosome enrichments
We then calculate positional enrichment to ensure that chromosome 8 activity
is indeed altered in the samples in 
[01_enrichments_copynumber.R](./01_enrichments_copynumber.R). This will
save files to the working directory with the current date. 

## Calculate mRNA-Prot correlations
The next step is to compute the correlations between copy number
computed in the previous step and mRNA/protein levels to identify functional pathway
enrichment and drug response measurements.  This is calculalted by
[02_enrichments_protRNA.R](./02_enrichments_protRNA.R). 

## Upload to synapse

There are many files that might be used later that are uploaded to synapse using
the [03_save_to_synapse.R](./03_save_to_synapse.R) script. It is designed to take
the current working directory and run immediately after the first four scripts.

# Figure panel generation

Now that the calculations have been made we can generate the figure panels in order
through the following scripts.

## Figure 1: correlations
Figure1_20250410.R

## Figure S1: histograms
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


