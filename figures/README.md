# Multi-omics integration of malignant peripheral nerve sheath tumors (MPNST) identifies potential targets based on chromosome 8q (chr8q) status
Since chr8q gain is associated with high-grade transformation in MPNST and 
inferior overall survival, we integrate multi-omics data to understand drivers 
and potential targets based on chr8q status. To do this, we collected new 
proteomics data and ran correlations between omics expression and chr8q copy
number in MPNST patient-derived xenografts (PDX) in addition to pathway 
analyses, network analyses, drug sensitivity predictions, fluorescent in situ
hybridization, and viability studies.

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


