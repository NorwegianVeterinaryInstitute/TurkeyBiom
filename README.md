# Code repository for the TurkeyBiom project

This repository holds the R code used to analyse and generate figures for the publication "Monensin phase-out in Norwegian turkey production decreases **Bifidobacterium** spp. abundance while enhancing microbial diversity".
All scripts can be found in the `scripts` directory.

- `01_sample_overview.R`: Script for generating figures and tables used to get an overview of the data.
- `02_16S_analysis.R`: Script for analysing the 16S data, used to determine which DNA extraction protocol to use for the shotgun data.
- `03_reads_and_correlation.R`: Basic QC of the metagenomic reads and number of classified/unclassified reads.
- `04_nonpareil_test_curves.R`: Nonpareil curves for the mock and negative control.
- `05_shotgun_mock_analysis.R`: Script for analysing the mock and negative control shotgun data.
- `06_generate_samplesheet.R`: Script for generating the input samplesheet to Taxprofiler.
- `07_generate_phyloseq_object.R`: Script for generating the clean phyloseq object from the kraken-biom file.
- `08_relative_abundance.R`: Script for calculating the relative abundances, running statistics, and generating figures.
- `09_calculate_diversity.R`: Script for calculating alpha- and beta-diversity, including RDA, general statistics and figures.
- `10_DESeq2_analysis.R`: Script for running the DESeq2 differential abundance analysis.
