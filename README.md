# Aiptasia-symbiosome-analysis

## Description
This repository contains R code used to generate all figures (found in the Data_analyses_files folder) and their underlying data (found in the Inputs folder) from Maruyama et al. 2025 https://www.biorxiv.org/content/10.1101/2025.10.09.679812v1.full

They also contain several analyses used to assess the contribution of _Breviolum minutum_ peptides to the results reported in the manuscript. These analyses and conclusions are described in detail below.

## Additional analyses of proteomics data
Because peptides detected in proteomics may ambiguously map to both Aiptasia and _B. minutum_ proteins, we performed a second Fragpipe analysis that searched the raw spectra against a combined database that included both Aiptasia and _B. minutum_ proteins. We then performed all analyses described in Figure 2 of the manuscript (i.e. normalization, differential enrichment, identification of known symbiosome proteins, transmembrane protein analyses, comparison to RNAseq data, and GO-term enrichment analyses - all reported in the "Data_analyses.md" file). We first examined the peptide output from this combined search to quantify the proportion of peptides that ambiguously mapped to proteins from both organisms. We found that 1.04% of the peptides that mapped to Aiptasia proteins also mapped to _B. minutum_ proteins, suggesting that peptides of ambiguous taxonomic origin should have a minor effect in our downstream analyses. However, 25/200 (12.5%) of the symbiosome-enriched proteins from the Aiptasia-only search and 8/183 (4%) of the symbiosome-enriched proteins from the combined search were not consistently symbiosome-enriched (see protein_comparison_output.csv). 

To test if the ambiguous peptides were causing this incongruency, we removed all ambiguous peptides that mapped to both organisms from the peptide outputs of both the Aiptasia-only search and the combined search. We then aggregated the peptides into proteins using a simplified method by summing the peptide intensities (Carrillo et al. 2010 https://doi.org/10.1093/bioinformatics/btp610). We then analyzed the data through the same analysis pipeline used to analyze the original results (i.e. normalization, differential enrichment, identification of known symbiosome proteins). We found that differential enrichment analysis using these new protein lists that removed ambiguous peptides resulted in the loss of 2/200 (1%) symbiosome-enriched proteins from the Aiptasia-only search and the loss of 4/183 (2.2%) symbiosome-enriched proteins from the combined Aiptasia and _B. minutum_ search. Among the proteins that were lost from the symbiosome-enriched lists in both analyses was Rab11A, a known symbiosome protein (Chen et al. 2005 https://doi.org/10.1016/j.bbrc.2005.10.133). All together, these results suggest that the ambiguous peptides were not driving the incongruence that we observed when comparing the symbiosome-enriched protein lists from the Aiptasia-only search and the combined Aiptasia and _B. minutum search_. Instead, the discrepancy may be caused by the increased stringency in peptide-spectrum match assignment when searching against the larger combined genome database, which expands the search space and shifts the false discovery rate threshold (Lee et al. 2023 https://doi.org/10.1128/msystems.00678-22).

To test the effect of database size on the inconsistencies we observed in symbiosome-enriched proteins from the two search strategies, we performed a third search that scrambled the _B. minutum_ proteins. First, we scrambled each protein sequence in the _B. minutum_ genome file using the shuffleseq command from EMBOSS ver. 2.2.0 (https://www.bioinformatics.nl/cgi-bin/emboss/shuffleseq). This function preserves the amino acid identity of each protein, but randomizes the primary sequence. Next, we combined the Aiptasia genome with the scrambled _B. minutum_ genome and added decoys and contaminants using the “Add decoy” function in Fragpipe-GUI. We then followed the same previously described methods to quantify raw spectra using Fragpipe and analyze the output. The results from this new search showed that only 0.07% of Aiptasia peptides mapped to scrambled _B. minutum_ proteins and differential enrichment analysis found 180 symbiosome-enriched proteins. However, 28/200 (14%) of the symbiosome-enriched proteins from the Aiptasia-only search and 19/183 (10.4%) of the symbiosome-enriched proteins from the combined search were not symbiosome-enriched in this new search that scrambled _B. minutum_ proteins. All together, these results support the hypothesis that the incongruency in symbiosome-enriched proteins between the search strategies was primarily driven by the larger database size, rather than from ambiguous peptides that map to both Aiptasia and _B. minutum_ proteins in the combined search.

Regardless of the incongruencies in the Aiptasia proteins identified by the Aiptasia-only and combined Aiptasia and _B. minutum_ search strategies, the results of follow-up analyses were largely consistent between all of the search results and analyses (i.e., 1. Successful enrichment/identification of known symbiosomal proteins. 2. Membrane transporters, lysosomal proteins, and vesicle trafficking proteins on the symbiosome, 3. Each protein, for localization by immunofluorescence and functional testing with shRNA and CRISPR/Cas9, was enriched on the symbiosome in both cases) (all data and figures are available in this repository in "Data_analyses.md"). 



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
