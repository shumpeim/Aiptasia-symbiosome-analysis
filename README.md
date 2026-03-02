# Aiptasia-symbiosome-analysis
Code for performing analyses and generating figures from Fragpipe output for Aiptasia symbiosome proteomics experiment.

## Requirements
- R 4.3.3
- Bioconductor 3.18

## Setup
1. Install R 4.3.3 from https://cran.r-project.org
2. Open R and run:
   install.packages("BiocManager")
   BiocManager::install(version = "3.18")
   install.packages("renv")
   renv::restore()
4. Open and knit Aiptasia_only_search.Rmd
