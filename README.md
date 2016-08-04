# README.md

The directories in this repository contain a pipeline to process and analyze the 16S sequence data on the ACISS cluster at the University of Oregon. It should hopefully be easy to customize for use on other computers.

## Directories and scripts

The scripts should be run in the following order. The input for these scripts is raw, demultiplexed 16S amplicon Illumina sequencing data in FastQ format. Auxilliary scripts called within the main scripts can be found in "Support\_scripts"

- 1\_assemble\_filter\_cat\_16S.job
- 2\_derep\_cluster\_ID\_16S.job

## Other required programs

Programs called by the scripts:

- flash/1.2.7
- fastx_toolkit/0.0.13
- bowtie/2.2.1
- usearch/7.0.1090 -OR- uclust/1.2.22
- mafft (v7.029b)
- fasttree/2.1.4
- RDPTools/140616

## Modular scripts

Some users voiced an interest in being able to run each step of the pipeline individually to see the output. The Modular\_scripts directory contains the exact same pipeline as the primary two scripts, but broken up in just such a fashion.

Because many of the steps after the first few scripts take very little time to run, it may not be necessary or even expedient to submit them as jobs. Therefore, all of these scripts can also be run as shell scripts within the raw reads directory (the same directory supplied to the `-d` PBS flag.