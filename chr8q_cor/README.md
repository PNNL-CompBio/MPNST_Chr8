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

## Compute chr8 copy numbers
The first step to analyzing the data is to calculate the copy number of chr8
across the samples. This is currently run by the 
[00_computeChrCorrelations.R](./00_computeChrCorrelations.R)
script.

This script downloads the relevant data to the environment, creates a folder
with today's date and saves the copy number data. 

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

## Network analysis

We then used the stored files to build the networks used in Figures 2 and 3. 
This is [04_networkAnalysis.R](./04_networkAnalysis.R)

# Figure panel generation

Now that the calculations have been made we can generate the figure panels in order
through the following scripts.

## Figure 1: correlations
To build the figures in a `figures` directory, first start with 
[Figure1_20250410.R](./Figure1_20250410.R)


### Figure 2: enrichment analyses
To visualize the enrichment analysis run [Figure2_20250429.R](./Figure2_2025_0429.R)


### Figure 3: drug sensitivity predictions
Next run [Figure3_20250513.R](./Figure3_20250513.R) which both runs the DMEA 
analysis and stores the files and also creates the figures for Figure 3.

### Compile tables
This requires that you run figure 3 code first, but will collect the files 
in the local directory so that they can be supplemental tables
[05_tableCompilation.R](./05_tableCompilation.R)


# Older code
The following files were not validated but also run parts of the analysis. 

## Figure S1: histograms
FigureS1_20250409.R

### Figure S3: drivers of enrichment
- Transcription factor targets: FigureS3_TF_20250418.R
- Kinase substrates: FigureS3_kinase_20250418.R

### Drug screen
- drugScreening.R: only used for IC50 t-tests
- fit_curve.py: calculates area under the curve values; not currently used; 
adapted from: https://github.com/PNNL-CompBio/coderdata/blob/main/coderbuild/utils/fit_curve.py


