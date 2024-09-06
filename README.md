# NicheNet Training

Welcome to the NicheNet training repository. This repository contains the materials and code as presented at the VIB training *Study intercellular communication with NicheNet* as presented on the 31st of May 2024.


## Objectives of the training

This training discusses how to analyze cell-cell communication from scRNA-seq data via the NicheNet analysis framework. The benefits and limitations of NicheNet are highlighted and compared to other approaches. The materials guide the course participant in applying NicheNet to their datasets. Most recent analysis features from NicheNet-v2 are also be covered. Finally, the materials also include how to analyze intercellular communication in multi-sample multi-condition scRNA-seq datasets with MultiNicheNet.

### Required skills

Participants should have some experience with R and know the basics of scRNA-seq data analysis.

## Installation

NicheNicheNet and MultiNicheNet are implemented as R packages. To run the code in this repository, you will need the following software installed:  
  - **R**: can be downloaded at [https://cran.r-project.org/](https://cran.r-project.org/).  
  - **RStudio**: while not required, this is highly recommended for an easier working environment.
RStudio Desktop can be downloaded at [https://posit.co/products/open-source/rstudio/](https://posit.co/products/open-source/rstudio/).  

To install the required R packages, run the following code in your R console:

```R
tools <- c("BiocManager", "devtools", "circlize", "DiagrammeR", "tidyverse")

for (tool in tools){
  if(!requireNamespace(tool, quietly = TRUE)) { 
    install.packages(tool)  
  } 
}

if(!requireNamespace("ComplexHeatmap", quietly = TRUE)) { 
  BiocManager::install("ComplexHeatmap", update = FALSE) 
} 

devtools::install_github("saeyslab/nichenetr") 
```  

Alternatively, if you are familiar with Docker, you can download our image at:

```
docker pull csangara/nichenetr:latest
```

This can be run with:
```
docker run -p 8787:8787 -e PASSWORD=bioc csangara/nichenetr:latest 
```
Then, go to `localhost:8787` in your web browser and login with the user **rstudio** and password **bioc**.

## Downloading additional files

To easily access the tutorial files, we recommend you to clone this repository:
```
git clone https://github.com/saeyslab/nichenet_training.git
```
Then open the `nichenet_training.Rproj` in RStudio, which will set `nichenet_training/` as the working directory. When creating a new folder to store the networks and Seurat objects, we assume that you are in this directory.

### NicheNet networks

Running NicheNet always requires general input data (**prior knowledge networks**). These files can be downloaded from [https://zenodo.org/records/7074291](https://zenodo.org/records/7074291).

You can also download them via the command line with: `wget https://zenodo.org/api/records/7074291/files-archive`

Once the files have been downloaded, create a `networks/` folder and extract the files there. You can run the script `check_network_files.R` to see if all required networks are correctly downloaded.

### Case studies

We provide three NicheNet case studies in the `vignettes/` folder that are further explained below. The Seurat objects used in each study can be downloaded at the following links: [[1]](https://zenodo.org/record/3531889/files/seuratObj.rds) [[2]](https://zenodo.org/records/11400203/files/seurat_obj_lite.rds) [[3]](https://zenodo.org/records/5840787/files/seurat_obj_subset_integrated_zonation.rds)

Or via the command line:
```
wget https://zenodo.org/record/3531889/files/seuratObj.rds
wget https://zenodo.org/records/11400203/files/seurat_obj_lite.rds
wget https://zenodo.org/records/5840787/files/seurat_obj_subset_integrated_zonation.rds
```
Create a new `data/` folder and place the objects there.

#### 1. Case-Control: Studying immune cell interactions in lymph nodes after viral infection (`case_control_example`)

In this example we use mouse NICHE-seq data from [Medaglia et al. (2017)](https://www.science.org/doi/10.1126/science.aao4277) to explore intercellular communication in the T cell area in the inguinal lymph node before and 72 hours after lymphocytic choriomeningitis virus (LCMV) infection. Specifically, we will prioritize which ligands can best explain the downstream changes after LCMV infection in CD8 T cells as the receiver population.

#### 2. Case-Control: Studying TME interactions in the context of anti-PD1 immunotherapy (`MultiNicheNet_demo`)

The vignette for the case control example can be adapted for a scRNA-seq data set from breast cancer biopsies of patients receiving anti-PD1 immune-checkpoint blockade therapy. [Bassez et al. (2021)](https://www.nature.com/articles/s41591-021-01323-8) collected from each patient one tumor biopsy before anti-PD1 therapy (“pre-treatment”) and one during subsequent surgery (“on-treatment”) A single-cell map of intratumoral changes during anti-PD1 treatment of patients with breast cancer. Based on additional scTCR-seq results, they identified one group of patients with clonotype expansion as response to the therapy (“E”) and one group with only limited or no clonotype expansion (“NE”). Specifically, we want to identify the most important interactions driving pre-therapy differences in CD4 T cells between expander and non-expander patients.

#### 3. Cell localization & differentiation: Studying the Kupffer cell niche (`differentiation_example`)

This vignettes uses a subset of the mouse liver scRNA-seq data generated by [Guilliams et al. (2022)](https://www.sciencedirect.com/science/article/pii/S0092867421014811) and focuses on identifying communication signals that determine Kupffer cell identity. Specifically, we will prioritize which ligands can best explain the differences between Kupffer cells and other liver macrophages.

## References and resources

The slides are available in the `slides` folder. 

The training materials are based on the following publications:  

- **NicheNet**: The NicheNet paper is available at [https://www.nature.com/articles/s41592-019-0667-57](https://www.nature.com/articles/s41592-019-0667-5). The NicheNet R package is available at [https://github.com/saeyslab/nichenetr](https://github.com/saeyslab/nichenetr) and contains many additional vignettes covering various follow-up analyses, visualizations and customizing your own prior model.  

- **NicheNet tutorial paper**: A tutorial paper on NicheNet that includes code for the case-control example is available as a pre-print at [https://arxiv.org/abs/2404.16358](https://arxiv.org/abs/2404.16358).

- **MultiNicheNet**: The MultiNicheNet preprint is available at  [https://www.biorxiv.org/content/10.1101/2023.06.13.544751v1](https://www.biorxiv.org/content/10.1101/2023.06.13.544751v1). The MultiNicheNet R package with accompanying vignettes is available at [https://github.com/saeyslab/multinichenetr](https://github.com/saeyslab/multinichenetr).

