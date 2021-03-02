#!/bin/bash

## Author: Simon Hackl, 28.02.2021
## Project: Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads
## This script can be used for assembly polishing with medaka. It was mainly used to polish consensus assemblies computed with trycycler (manually).
## To use the script, two input parameters have to be specified:
## - The first specifies the assembly (FASTA file) to be polished.
## - The second specifies the reference long read set (FASTQ file) to be used for polishing.
## NOTE: The specified long reads should be the same as used for generating the specified assembly.

# Generate directories to store the output files.
mkdir ./results
mkdir ./results/assemblies
mkdir ./temp
mkdir ./temp/medaka-${NAME}

# Specify the number of threads to use.
THREADS=6

# The NAME and SAMPLE variables are used for correct naming.
NAME="${2%.*}"
SAMPLE="${2%-*}"

# The following variables are used to access tools used in this script.
# Change these to point to the executables of the respecitve programs if necessary.
MEDAKA="./tools/medaka"

source ${MEDAKA}/bin/activate
medaka_consensus -i ${1} -d ${2} -o ./temp/medaka-${NAME} -t ${THREADS}
cp ./temp/medaka-${NAME}/consensus.fasta ./results/assemblies/${SAMPLE}/${NAME}_medaka.fasta
rm -r ./temp/medaka-${NAME}