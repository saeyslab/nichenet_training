# NicheNet Training

Welcome to the NicheNet training repository. This repository contains the materials and code as presented at the VIB training *Study intercellular communication with NicheNet* as presented on the 31st of May 2024.


## Objectives of the training

This training discusses how to analyze cell-cell communication from scRNA-seq data via the NicheNet analysis framework. The benefits and limitations of NicheNet are highlighted and compared to other approaches. The materials guide the course participant in applying NicheNet to their datasets. Most recent analysis features from NicheNet-v2 are also be covered. Finally, the materials also include how to analyze intercellular communication in multi-sample multi-condition scRNA-seq datasets with MultiNicheNet.

## Required skills

Participants should have some experience with R and know the basics of scRNA-seq data analysis.

## References and resources

The slides are available in the `slides` folder. 

The training materials are based on the following publications:  

1. **NicheNet**:   
The NicheNet paper is available at [https://www.nature.com/articles/s41592-019-0667-57](https://www.nature.com/articles/s41592-019-0667-5). The NicheNet R package is available at [https://github.com/saeyslab/nichenetr](https://github.com/saeyslab/nichenetr) and contains many additional vignettes covering various follow-up analyses, visualizations and customizing your own prior model.  

2. **NicheNet tutorial paper**:  
A tutorial paper on NicheNet that includes code for the case-control example is available as a pre-print at [https://arxiv.org/abs/2404.16358](https://arxiv.org/abs/2404.16358).

3. **MultiNicheNet**:  
The MultiNicheNet preprint is available at  [https://www.biorxiv.org/content/10.1101/2023.06.13.544751v1](https://www.biorxiv.org/content/10.1101/2023.06.13.544751v1). The MultiNicheNet R package with accompanying vignettes is available at [https://github.com/saeyslab/multinichenetr](https://github.com/saeyslab/multinichenetr).


## Installation instructions

NicheNicheNet and MultiNicheNet are implemented as R packages. To run the code in this repository, you will need the following software installed:  
  - *R*: can be downloaded at [https://cran.r-project.org/](https://cran.r-project.org/).  
  - *RStudio*: while not required, this is highly recommended for an easier working environment.
RStudio Desktop can be downloaded at [https://posit.co/products/open-source/rstudio/](https://posit.co/products/open-source/rstudio/).  

To install the required R packages, run the following code in your R console:

```R
if(!requireNamespace("BiocManager", quietly = TRUE)) { 
  install.packages("BiocManager")  
} 

if(!requireNamespace("devtools", quietly = TRUE)) { 
  install.packages("devtools")  
} 

if(!requireNamespace("limma", quietly = TRUE)) { 
  BiocManager::install("limma", update = FALSE) 
} 

if(!requireNamespace("ComplexHeatmap", quietly = TRUE)) { 
  BiocManager::install("ComplexHeatmap", update = FALSE) 
} 

if(!requireNamespace("circlize", quietly = TRUE)) { 
  install.packages("circlize")  
} 

if(!requireNamespace("DiagrammeR", quietly = TRUE)) { 
  install.packages("DiagrammeR")  
} 

if(!requireNamespace("tidyverse", quietly = TRUE)) { 
  install.packages("tidyverse")  
} 

devtools::install_github("saeyslab/nichenetr") 
```  

### Downloading the NicheNet prior model

Running NicheNet always requires general input data (**prior knowledge networks**). These files can be downloaded from [https://zenodo.org/records/7074291](https://zenodo.org/records/7074291).  

## Available case studies

### 1. Case-Control: Studying immune cell interactions in lymph nodes after viral infection

In this example we use mouse NICHE-seq data from [Medaglia et al (2017)](https://www.science.org/doi/10.1126/science.aao4277) to explore intercellular communication in the T cell area in the inguinal lymph node before and 72 hours after lymphocytic choriomeningitis virus (LCMV) infection. Specifically, we will prioritize which ligands can best explain the downstream changes after LCMV infection in CD8 T cells as the receiver population.

The vignette is available as `case_control_example.Rmd`. The code is also available as an R script in `case_control_example.R`.


### 2. Case-Control: Studying TME interactions in the context of anti-PD1 immunotherapy

The vignette for the case control example can be adapted for a scRNAseq data set from breast cancer biopsies of patients receiving anti-PD1 immune-checkpoint blockade therapy. [Bassez et al (2021))](https://www.nature.com/articles/s41591-021-01323-8) collected from each patient one tumor biopsy before anti-PD1 therapy (“pre-treatment”) and one during subsequent surgery (“on-treatment”) A single-cell map of intratumoral changes during anti-PD1 treatment of patients with breast cancer. Based on additional scTCR-seq results, they identified one group of patients with clonotype expansion as response to the therapy (“E”) and one group with only limited or no clonotype expansion (“NE”). Specifically, we want to identify the most important interactions driving pre-therapy differences in CD4 T cells between expander and non-expander patients. The data is available at [https://zenodo.org/records/11394036](https://zenodo.org/records/11394036).


### 3. Cell localization & differentiation: Studying the Kupffer cell niche

This vignettes uses a subset of the mouse liver scRNA-seq data generated in the [Guilliams et al (2022) paper](https://www.sciencedirect.com/science/article/pii/S0092867421014811) and focuses on identifying communication signals that determine Kupffer cell identity. Specifically, we will prioritize which ligands can best explain the differences between Kupffer cells and other liver macrophages.

The vignette is available as `differentiation_example.Rmd`. 
