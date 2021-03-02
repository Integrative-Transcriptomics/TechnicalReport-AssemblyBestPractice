#!/bin/bash

## Author: Simon Hackl, 28.02.2021
## Project: Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads
## This script is used for sub-sampling (and coverage analysis of these sub-samples) of the CFT073 and RN4220 reference data sets.
## - For sub-sampling the tool rasusa is used:
##   It was decided to generate down-samples with expected 200, 150, 100, 80, 60, 40, 20, 15, 10, 8, 6, 4, 2 and 1 X coverage for the CFT073 and RN4220 reference sets.
##   NOTE: rasusa expects the genome size as a parameter to run, for CFT073 the given length of the reference genome was rounded off to 5.2 million base pairs.
##   As the RN4220 reference genome is highly fragmented it was decided to use 2.8 million base pairs, the median length of 12,401 genome assemblies of S. aureus 
##   (source: https://www.ncbi.nlm.nih.gov/genome/?term=Staphylococcus%20aureus[Organism]&cmd=DetailsSearch), as the genomes length.
## - For the coverage analysis the same procedure as described in S3_readQualityControl.sh is applied (mapping with minimap2 and per position coverage computation with samtools depth).
##   The resulting .cov files will be stored at './results/coverage/<SAMPLE-IDENTIFIER>'.

# Generate directories to store the output files.
mkdir ./results
mkdir ./results/coverage
mkdir ./results/coverage/CFT073
mkdir ./results/coverage/RN4220

# Specify the number of threads to use.
THREADS=6

# The following variables are used to access tools used in this script.
# Change these to point to the executables of the respecitve programs if necessary.
MINIMAP2="minimap2"
SAMTOOLS="samtools"
RASUSA="./tools/rasusa"


for I in 200 150 100 80 60 40 20 15 10 8 6 4 2 1
do
	# Generate sub-sample with I X depth of coverage...
	# ...of the RN4220 (Guppy 3.2.1.0) longreads and map the sub-sampled reads to the reference genome.
	${RASUSA} -c $I -g 2.8MB -i ./data/reads/RN4220/QNFLR049AW~guppy3210-t.fastq -o $DIR/Data/ONT-reads/SAureus-RN4220/QNFLR049AW~guppy3210-t-${I}X.fastq
	${MINIMAP2} -x map-ont -a -t $THREADS -o ./results/coverage/RN4220/QNFLR049AW~guppy3210-t-${I}X.sam ./data/references/RN4220.fasta ./data/reads/RN4220/QNFLR049AW~guppy3210-t-${I}X.fastq
	${SAMTOOLS} sort -O bam -o ./results/coverage/RN4220/QNFLR049AW~guppy3210-t-${I}X.bam ./results/coverage/RN4220/QNFLR049AW~guppy3210-t-${I}X.sam
	${SAMTOOLS} depth -a -o ./results/coverage/RN4220/QNFLR049AW~guppy3210-t-${I}X.cov ./results/coverage/RN4220/QNFLR049AW~guppy3210-t-${I}X.bam
	rm ./results/coverage/RN4220/QNFLR049AW~guppy3210-t-${I}X.sam
	rm ./results/coverage/RN4220/QNFLR049AW~guppy3210-t-${I}X.bam
	# ...of the RN4220 (Guppy 4.0.1.1) longreads and map the sub-sampled reads to the reference genome.
	${RASUSA} -c $I -g 2.8MB -i ./data/reads/RN4220/QNFLR049AW~guppy4011-t.fastq -o $DIR/Data/ONT-reads/SAureus-RN4220/QNFLR049AW~guppy4011-t-${I}X.fastq
	${MINIMAP2} -x map-ont -a -t $THREADS -o ./results/coverage/RN4220/QNFLR049AW~guppy4011-t-${I}X.sam ./data/references/RN4220.fasta ./data/reads/RN4220/QNFLR049AW~guppy4011-t-${I}X.fastq
	${SAMTOOLS} sort -O bam -o ./results/coverage/RN4220/QNFLR049AW~guppy4011-t-${I}X.bam ./results/coverage/RN4220/QNFLR049AW~guppy4011-t-${I}X.sam
	${SAMTOOLS} depth -a -o ./results/coverage/RN4220/QNFLR049AW~guppy4011-t-${I}X.cov ./results/coverage/RN4220/QNFLR049AW~guppy4011-t-${I}X.bam
	rm ./results/coverage/RN4220/QNFLR049AW~guppy4011-t-${I}X.sam
	rm ./results/coverage/RN4220/QNFLR049AW~guppy4011-t-${I}X.bam
	# ...of the CFT073 longreads map the sub-sampled reads to the reference genome.
	${RASUSA} -c $I -g 5.2MB -i ./data/reads/CFT073/SRR8494940-t.fastq -o ./data/reads/CFT073/SRR8494940-t-${I}X.fastq
	${MINIMAP2} -x map-ont -a -t $THREADS -o ./results/coverage/CFT073/SRR8494940-t-${I}X.sam ./data/references/CFT073.fasta ./data/reads/CFT073/SRR8494940-t-${I}X.fastq
	${SAMTOOLS} sort -O bam -o ./results/coverage/CFT073/SRR8494940-t-${I}X.bam ./results/coverage/CFT073/SRR8494940-t-${I}X.sam
	${SAMTOOLS} depth -a -o ./results/coverage/CFT073/SRR8494940-t-${I}X.cov ./results/coverage/CFT073/SRR8494940-t-${I}X.bam
	rm ./results/coverage/CFT073/SRR8494940-t-${I}X.sam
	rm ./results/coverage/CFT073/SRR8494940-t-${I}X.bam
done