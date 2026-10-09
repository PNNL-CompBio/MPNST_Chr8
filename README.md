# Multi-omics integration of malignant peripheral nerve sheath tumors (MPNST) 
Since chr8q gain is associated with high-grade transformation in MPNST and 
inferior overall survival, we integrate multi-omics data to understand drivers 
and potential targets based on chr8q status. To do this, we include two studies
of MPNST omics measruements. 


## Chr8q correlation analysis
The analysis code for the project is included in this repository. It requires
a large amount of proteomic data processing, contained in the 
[proteomics](./proteomics) directory. The code included in this folder is
for information only, as it requires access to the PNNL intranet to operate.
The scripts to run the analysis and  genrerate figures for this paper are included in the 
[chr8q_cor](./chr8q_cor) directory.

The manuscript can be found [on BioRxiv](https://www.biorxiv.org/content/10.64898/2026.01.26.701599v1) and access to the data can be requested on [Synapse](http://synapse.org/mpnst_chr8). 

## FAK Analysis 
We also carry out omics analysis implicating FAK in cancer progression. 
The code to carry out the FAK in the [FAK Analysis](./FAK analysis)
directory. This also requires the dependencies below.

## Install dependencies
The analysis for both manuscripts requires numerous dependencies, they are 
listed in the [installDependences.R](./installDependencies.R) script. 


