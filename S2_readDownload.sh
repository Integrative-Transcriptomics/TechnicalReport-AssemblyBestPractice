#!/bin/bash
## This script is used to download the CFT073 and MGH78578 public read sets.
## Author: Simon Hackl, 23.02.2021
## Project: Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads

# Generate directories to store the reference FASTA files, genome annotations in GFF3 format and temporary data.
mkdir ./data
mkdir ./data/reads
mkdir ./data/reads/CFT073
mkdir ./data/reads/MGH78578
# Definition of variables referring to the tools used. Change these to point to the executables of the programs if appropriate (i.e. if such are not accessible via path variables).
FASTQ-DUMP="fastq-dump"

# The reads are downloaded with fastq-dump by specifying the respective run identifiers.
${FASTQ-DUMP} --outdir ./data/reads/CFT073/ --split-files SRR8482585
${FASTQ-DUMP} --outdir ./data/reads/CFT073/ SRR8494940
${FASTQ-DUMP} --outdir ./data/reads/MGH78578/ --split-files SRR8482567
${FASTQ-DUMP} --outdir ./data/reads/MGH78578/ SRR8494915
## NOTE: The RN4220 read sets are not publicy available. These were stored in the ./data/reads/RN4420 directory as follows:
## * QNFLR049AW~guppy3210.fastq, longreads basecalled with Guppy version 3.2.1.0
## * QNFLR049AW~guppy4011.fastq, longreads basecalled with Guppy version 4.0.1.1
## * QNFLR056AF_1.fastq, shortreads (first of paired end, forward strand)
## * QNFLR056AF_2.fastq, shorteads (second of paired end, reverse strand)