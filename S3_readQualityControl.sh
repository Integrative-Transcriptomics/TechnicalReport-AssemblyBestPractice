#!/bin/bash

## Author: Simon Hackl, 28.02.2021
## Project: Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads
## This script is used to run quality control and trimming of the read sets:
## - Long reads are be trimmed using Porechop, short reads undergo no trimming.
##   The porechop logs will be stored at './results/porechop/<SAMPLE-IDENTIFIER>'
## - Both, long and short reads, are analyzed using FastQC.
##   The .html and .zip reports will be stored at './results/fastqc/<SAMPLE-IDENTIFIER>'
## - In addition minimap2 is used to map the reads to the respective references.
##   For minimap2 the parameters `-a` and `-x map-ont` (for long reads) or `-x sr` (for short reads) were set,
##   in order to output SAM format files and adjsut some mapping parameters accounting for the read types, respectively.
##   The resulting SAM files are sorted and converted to BAM files using `samtools sort -O bam` and on these
##   `samtools depth -a` (`-a` to write all positions including such with zero coverage) is called.
##   The final .cov files will be stored at './results/coverage/<SAMPLE-IDENTIFIER>'.

# Generate directories to store the output files.
mkdir ./results
mkdir ./results/porechop
mkdir ./results/fastqc
mkdir ./results/coverage
mkdir ./temp
mkdir ./results/porechop/MGH78578
mkdir ./results/fastqc/MGH78578
mkdir ./results/coverage/MGH78578
mkdir ./results/porechop/CFT073
mkdir ./results/fastqc/CFT073
mkdir ./results/coverage/CFT073
mkdir ./results/porechop/RN4220
mkdir ./results/fastqc/RN4220
mkdir ./results/coverage/RN4220

# Specify the number of threads to use.
THREADS=6

# The following variables are used to access tools used in this script.
# Change these to point to the executables of the respecitve programs if necessary.
MINIMAP2="minimap2"
SAMTOOLS="samtools"
PORECHOP="./tools/porechop/porechop-runner.py"
FASTQC="fastqc"

# Run quality control of MGH78578 reference dataset.
${PORECHOP} -i  ./data/reads/MGH78578/SRR8494915.fastq -o ./data/reads/MGH78578/SRR8494915-t.fastq -t $THREADS >> ./results/porechop/MGH78578/porechop-SRR8494915.log
${FASTQC} -d ./temp -o ./results/fastqc/MGH78578/ -t $THREADS ./data/reads/MGH78578/SRR8494915-t.fastq
${MINIMAP2} -x map-ont -a -t $THREADS -o ./results/coverage/MGH78578/SRR8494915-t.sam ./data/references/MGH78578.fasta ./data/reads/MGH78578/SRR8494915-t.fastq
${SAMTOOLS} sort -O bam -o ./results/coverage/MGH78578/SRR8494915-t.bam ./results/coverage/MGH78578/SRR8494915-t.sam
${SAMTOOLS} depth -a -o ./results/coverage/MGH78578/SRR8494915-t.cov ./results/coverage/MGH78578/SRR8494915-t.bam
rm ./results/coverage/MGH78578/SRR8494915-t.sam
rm ./results/coverage/MGH78578/SRR8494915-t.bam
${FASTQC} -d ./temp -o ./results/fastqc/MGH78578/ -t $THREADS ./data/reads/MGH78578/SRR8482567_1.fastq ./data/reads/MGH78578/SRR8482567_2.fastq
${MINIMAP2} -x sr -a -t $THREADS -o ./results/coverage/MGH78578/SRR8482567_1_2.sam ./data/references/MGH78578.fasta ./data/reads/MGH78578/SRR8482567_1.fastq ./data/reads/MGH78578/SRR8482567_2.fastq
${SAMTOOLS} sort -O bam -o ./results/coverage/MGH78578/SRR8482567_1_2.bam ./results/coverage/MGH78578/SRR8482567_1_2.sam
${SAMTOOLS} depth -a -o ./results/coverage/MGH78578/SRR8482567_1_2.cov ./results/coverage/MGH78578/SRR8482567_1_2.bam
rm ./results/coverage/MGH78578/SRR8482567_1_2.sam
rm ./results/coverage/MGH78578/SRR8482567_1_2.bam

# Run quality control of CFT073 reference dataset.mkdir
${PORECHOP} -i  ./data/reads/CFT073/SRR8494940.fastq -o ./data/reads/CFT073/SRR8494940-t.fastq -t $THREADS >> ./results/porechop/CFT073/porechop-SRR8494940.log
${FASTQC} -d ./temp -o ./results/fastqc/CFT073/ -t $THREADS ./data/reads/CFT073/SRR8494940-t.fastq
${MINIMAP2} -x map-ont -a -t $THREADS -o ./results/coverage/CFT073/SRR8494940-t.sam ./data/references/CFT073.fasta ./data/reads/CFT073/SRR8494940-t.fastq
${SAMTOOLS} sort -O bam -o ./results/coverage/CFT073/SRR8494940-t.bam ./results/coverage/CFT073/SRR8494940-t.sam
${SAMTOOLS} depth -a -o ./results/coverage/CFT073/SRR8494940-t.cov ./results/coverage/CFT073/SRR8494940-t.bam
rm ./results/coverage/CFT073/SRR8494940-t.sam
rm ./results/coverage/CFT073/SRR8494940-t.bam
${FASTQC} -d ./temp -o ./results/fastqc/CFT073/ -t $THREADS ./data/reads/CFT073/SRR8482585_1.fastq ./data/reads/CFT073/SRR8482585_2.fastq
${MINIMAP2} -x sr -a -t $THREADS -o ./results/coverage/CFT073/SRR8482585_1_2.sam ./data/references/CFT073.fasta ./data/reads/CFT073/SRR8482585_1.fastq ./data/reads/CFT073/SRR8482585_2.fastq
${SAMTOOLS} sort -O bam -o ./results/coverage/CFT073/SRR8482585_1_2.bam ./results/coverage/CFT073/SRR8482585_1_2.sam
${SAMTOOLS} depth -a -o ./results/coverage/CFT073/SRR8482585_1_2.cov ./results/coverage/CFT073/SRR8482585_1_2.bam
rm ./results/coverage/CFT073/SRR8482585_1_2.sam
rm ./results/coverage/CFT073/SRR8482585_1_2.bam

# Run quality control of RN4220 reference dataset.
${PORECHOP} -i  ./data/reads/RN4220/QNFLR049AW~guppy3210.fastq -o ./data/reads/RN4220/QNFLR049AW~guppy3210-t.fastq -t $THREADS >> ./results/porechop/RN4220/porechop-QNFLR049AW~guppy3210.log
${FASTQC} -d ./temp -o ./results/fastqc/RN4220/ -t $THREADS ./data/reads/RN4220/QNFLR049AW~guppy3210-t.fastq
${MINIMAP2} -x map-ont -a -t $THREADS -o ./results/coverage/RN4220/QNFLR049AW~guppy3210-t.sam ./data/references/RN4220.fasta ./data/reads/RN4220/QNFLR049AW~guppy3210-t.fastq
${SAMTOOLS} sort -O bam -o ./results/coverage/RN4220/QNFLR049AW~guppy3210-t.bam ./results/coverage/RN4220/QNFLR049AW~guppy3210-t.sam
${SAMTOOLS} depth -a -o ./results/coverage/RN4220/QNFLR049AW~guppy3210-t.cov ./results/coverage/RN4220/QNFLR049AW~guppy3210-t.bam
rm ./results/coverage/RN4220/QNFLR049AW~guppy3210-t.sam
rm ./results/coverage/RN4220/QNFLR049AW~guppy3210-t.bam
${PORECHOP} -i  ./data/reads/RN4220/QNFLR049AW~guppy4011.fastq -o ./data/reads/RN4220/QNFLR049AW~guppy4011-t.fastq -t $THREADS >> ./results/porechop/RN4220/porechop-QNFLR049AW~guppy4011.log
${FASTQC} -d ./temp -o ./results/fastqc/RN4220/ -t $THREADS ./data/reads/RN4220/QNFLR049AW~guppy4011-t.fastq
${MINIMAP2} -x map-ont -a -t $THREADS -o ./results/coverage/RN4220/QNFLR049AW~guppy4011-t.sam ./data/references/RN4220.fasta ./data/reads/RN4220/QNFLR049AW~guppy4011-t.fastq
${SAMTOOLS} sort -O bam -o ./results/coverage/RN4220/QNFLR049AW~guppy4011-t.bam ./results/coverage/RN4220/QNFLR049AW~guppy4011-t.sam
${SAMTOOLS} depth -a -o ./results/coverage/RN4220/QNFLR049AW~guppy4011-t.cov ./results/coverage/RN4220/QNFLR049AW~guppy4011-t.bam
rm ./results/coverage/RN4220/QNFLR049AW~guppy4011-t.sam
rm ./results/coverage/RN4220/QNFLR049AW~guppy4011-t.bam
${FASTQC} -d ./temp -o ./results/fastqc/RN4220/ -t $THREADS ./data/reads/RN4220/QNFLR056AF_1.fastq ./data/reads/RN4220/QNFLR056AF_2.fastq
${MINIMAP2} -x sr -a -t $THREADS -o ./results/coverage/RN4220/QNFLR056AF_1_2.sam ./data/references/RN4220.fasta ./data/reads/RN4220/QNFLR056AF_1.fastq ./data/reads/RN4220/QNFLR056AF_2.fastq
${SAMTOOLS} sort -O bam -o ./results/coverage/RN4220/QNFLR056AF_1_2.bam ./results/coverage/RN4220/QNFLR056AF_1_2.sam
${SAMTOOLS} depth -a -o ./results/coverage/RN4220/QNFLR056AF_1_2.cov ./results/coverage/RN4220/QNFLR056AF_1_2.bam
rm ./results/coverage/RN4220/QNFLR056AF_1_2.sam
rm ./results/coverage/RN4220/QNFLR056AF_1_2.bam