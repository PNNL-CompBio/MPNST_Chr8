# Multi-omics integration of malignant peripheral nerve sheath tumors (MPNST) identifies potential targets based on chromosome 8q (chr8q) status
Since chr8q gain is associated with high-grade transformation in MPNST and 
inferior overall survival, we integrate multi-omics data to understand drivers 
and potential targets based on chr8q status. To do this, we collected new 
proteomics data and ran correlations between omics expression and chr8q copy
number in MPNST patient-derived xenografts (PDX) in addition to pathway 
analyses, network analyses, drug sensitivity predictions, fluorescent in situ
hybridization, and viability studies.

The analysis for this manuscript integrates copy number and RNASeq measurements
of diverse MPNST together with proteomics measurements. The proteomics analysis
is stores in [the proteomics directory](../proteomics) and stores the files
on Synapse at http://synapse.org/chr8_mpnst. Once you request access to the
data files you can reproduce the analysis in the paper below.

## Compute chr8 copy numbers
The first step to analyzing the data is to calculate the copy number of chr8
across the samples. This is currently run by the [5-panSEA_using_helper.R]()
script. Be mindful of running this script as it creates numerous 
sub-directories.

The results of this script include:
- supplemental data tables in the `data` directory.
- copy number-associated analysis in the 
`analysis/Chr8_quant_20250409/positional_medians/20250616/Copy\ Number/`
directory
- RNA-seq related correlations in the
- Protein-related correlations in the 

## Figure panel generation

Figure 1 collects the [chr8 copy number values]() computed by the above step
and uploaded to synapse and plots 

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


