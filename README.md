# Scripts used in the Lyve-SET paper

Not all scripts were ultimately used.  However, here is a discription of some of the more important scripts used.  Most or all scripts can be run with no options or with `--help` for additional information.

## compareSnps.sh, compareTrees.sh

Compares all pipelines in a given directory.  Looks for the Lyve-SET, kSNP, RealPhy, SNP-Pipeline, and SNVPhyl project directories.

## distanceMatrixToImage.pl

Creates a heatmap image from a 2d distance matrix

## launch_ksnp.sh, launch_realphy.sh, launch_set_alreadyShuffled.sh, launch_snp-pipeline.sh

Runs a given SNP pipeline in a directory.  There must be a "reads" subdirectory where interleaved fastq.gz files are present.  The first parameter must be a reference assembly fasta file.

## pairwiseAlleleDistances.pl

Uncover distances between genomes in allele differences. The input file must be an alleles CSV file from BioNumerics.

## pairwiseScatterplot.R

Create a scatter plot between two or more pipelines

## pipelineSnpsVsPipeline.pl

Compare SNP positions in one pipeline vs another and generate sensitivity and specificity

