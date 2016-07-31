# README.md

The directories in this repository contain the scripts used to process and analyze the 16S sequence data for "The role of adaptive immunity as an ecological filter on the gut microbiota in zebrafish"

## Directories and scripts

The scripts should be run in the following order. The input for these scripts is raw, demultiplexed 16S amplicon Illumina sequencing data in FastQ format. Auxilliary scripts called within the main scripts can be found in a "Support\_scripts" directory in each main directory below.

1. Sequence\_processing\_and\_OTU\_picking/
	1. assemble\_filter\_cat\_16S.job
	2. derep\_cluster\_ID\_16S.job
2. R\_analysis/
	1. transforming\_or\_rarefying\_counts.R
	2. basic\_ecological\_analysis.R
	3. neutral\_theory\_model.R
	4. phyloseq\_to\_lefse\_format.R

## Other required programs

Programs called by the scripts in 1\_Sequence\_processing\_and\_OTU\_picking/

- flash/1.2.7
- fastx_toolkit/0.0.13
- bowtie/2.2.1
- usearch/7.0.1090 -OR- uclust/1.2.22
- mafft (v7.029b)
- fasttree/2.1.4
- RDPTools/140616

## Modular scripts
Some users voiced an interest in being able to run each step of the pipeline individually to see the output. The Modular\_scripts directory contains the exact same pipeline as Sequence\_processing\_and\_OTU\_picking/, but broken up in just such a fashion.

Because many of the steps after the first few scripts take very little time to run, it may not be necessary or even expedient to submit them as jobs. Therefore, all of these scripts can also be run as shell scripts within the raw reads directory (the same directory supplied to the `-d` PBS flag.