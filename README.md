# README.md

The directories in this repository contain the scripts used to process and analyze the 16S sequence data for "The role of adaptive immunity as an ecological filter on the gut microbiota in zebrafish"

## Directories and scripts

The scripts should be run in the following order. The input for these scripts is raw, demultiplexed 16S amplicon Illumina sequencing data in FastQ format. Auxilliary scripts called within the main scripts can be found in a "Support_scripts" directory in each main directory below

1. Sequence_processing_and_OTU_picking/
	1. assemble_filter_cat_16S.job
	2. derep_cluster_ID_16S.job
2. R_analysis/
	1. transforming_or_rarefying_counts.R
	2. basic_ecological_analysis.R
	3. neutral_theory_model.R
	4. phyloseq_to_lefse_format.R