# Aiptasia-symbiosome-analysis

## Description
This repository contains R code used to generate all figures (found in the Data_analyses_files folder) and their underlying data (found in the Inputs folder) from Maruyama et al. 2025 https://www.biorxiv.org/content/10.1101/2025.10.09.679812v1.full

The repository also contains several additional control analyses done to assess the contribution of _Breviolum minutum_ peptides and the effect of database size to the results reported in the manuscript. These analyses and conclusions are described in detail in the Ambiguous_peptide_analysis_for_proteomics.md file.

## Requirements for running the "Data_analyses.rmd" code
- R 4.3.3
- Bioconductor 3.18

## Setup
1. Install R 4.3.3 from https://cran.r-project.org
2. Open "Aiptasia Symbiosome analysis.Rproj" in R and run:
   install.packages("BiocManager")
   BiocManager::install(version = "3.18")
   install.packages("renv")
   renv::restore()
4. Open and knit the rmd file.

## To see summary of analysis
Please open "Data_analyses.md" on Github to see output from rmd file.
