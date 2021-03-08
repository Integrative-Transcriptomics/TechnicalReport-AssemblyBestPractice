#!/bin/bash

## Author: Simon Hackl, 28.02.2021
## Project: Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads
## This script is used for assembly polishing with short reads using pilon. It was mainly used to polish consensus assemblies computed with trycycler (manually).
## To use the script, three input parameters have to be specified:
## - The first and second specify the input paired-end short reads (FASTQ files) to be used for polishing. 
## - The third specifies the assembly (FASTA file) to be polished.
## Author: Simon Hackl, 23.02.2021 
## Project: Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads

# Specify the number of threads to use.
THREADS=6

# The NAME variable is used for correct naming.
NAME="${3%.*}"

# Generate directories to store the output files.
mkdir ./temp
mkdir ./temp/pilon-${NAME}

# The following variables are used to access tools used in this script.
# Change these to point to the executables of the respecitve programs if necessary.
PILON="pilon"
BWA="bwa"
SAMTOOLS="samtool"

# The script follows the instruction of the pilon wiki (https://github.com/broadinstitute/pilon/wiki/Requirements-&-Usage).
# First the short reads are mapped to the draft assembly using the burrows-wheeler aligner. The result of the mapping is sorted, converted to BAM format and indexed.
# Pilon uses this indexed BAM and the draft assembly as input. 
${BWA} index ${3}
${BWA} mem -t $THREADS ${3} ${1} ${2} > ./temp/pilon-${NAME}/aln.sam
${SAMTOOLS} sort -O bam -@ $THREADS -o ./temp/pilon-${NAME}/aln.bam ./temp/pilon-${NAME}/aln.sam
${SAMTOOLS} index -b -@ $THREADS ./temp/pilon-${NAME}/aln.bam
${PILON} --genome $3 --frags ./temp/pilon-${NAME}/aln.bam --outdir `dirname ${3}` --threads $THREADS
rm -r ./temp/pilon-${NAME}