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
4. Open and knit the rmd files!

## Description
These analyses include code used to analyze symbiosome proteomics data outputted from Fragpipe.
They include Fragpipe output from:
1. An Aiptasia-only genome search
2. A combined search including Aiptasia and _Breviolum minutum_ genomes

These analyses were normalized and analyzed for differential enrichment using NormalyzerDE.
Plots include volcano plots labeling previously known symbiosome proteins, PCA plots, Cell-compartment GO-term enrichment plots, plots analyzing the enrichment of transmembrane proteins and plots analyzing the proportion of proteins that are transcriptionally upregulated with symbiosis.

Also included are raw data sheets (in the folder "Inputs") and code used to generate every plot.

## To see summary of analysis
Please open "Data_analyses.md" on Github to see output from rmd file.
