FROM bioconductor/bioconductor_docker:devel-R-4.4.1

RUN R -e 'devtools::install_github("saeyslab/nichenetr")'
