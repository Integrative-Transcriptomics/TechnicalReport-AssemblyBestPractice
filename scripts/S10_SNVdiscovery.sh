 #!/bin/bash

## Author: Simon Hackl, 28.02.2021
## Project: Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads
## This script is used to run bcftools for SNV discovery from mapping the reference reads to the reference genome for the CFT073 reference set.
## The SNV discovery is run three times, one time for either of the short and long reads solely and one time for combining both read sets.
## NOTE: It was decided so only discover substitutions and not indels as SNVs for simplicity.
##
## The resulting .vcf files will be stored at './results/SNVdiscovery/CFT073/bcftools/<READ-RUN-IDENTIFIERS>-bcftools.vcf'

# Generate directories to store the output files.
mkdir ./results
mkdir ./results/SNVdiscovery
mkdir ./results/SNVdiscovery/CFT073
mkdir ./results/SNVdiscovery/CFT073/bcftools

# Specify the number of threads to use.
THREADS=6

# The following variables are used to access tools used in this script.
# Change these to point to the executables of the respecitve programs if necessary.
BCFTOOLS="bcftools"
MINIMAP2="minimap2"
SAMTOOLS="samtools"

# Map the long and short reads to the reference genome, the same procedure as in `S3_readQualityControl.sh` is applied.
${MINIMAP2} -x map-ont -a -t $THREADS -o ./temp/SRR8494940.sam ./data/references/CFT073.fasta ./data/reads/CFT073/SRR8494940-t.fastq
${SAMTOOLS} sort -O bam -o ./temp/SRR8494940.bam ./temp/SRR8494940.sam
rm ./temp/SRR8494940.sam
${MINIMAP2} -x sr -a -t $THREADS -o ./temp/SRR8482585.sam ./data/references/CFT073.fasta ./data/reads/CFT073/SRR8482585_1.fastq ./data/reads/CFT073/SRR8482585_2.fastq
${SAMTOOLS} sort -O bam -o ./temp/SRR8482585.bam ./temp/SRR8482585.sam
rm ./temp/SRR8482585.sam

# The reference genome has to be indexed for using bcftools.
${SAMTOOLS} faidx ./data/references/CFT073.fasta

# The mpielup sub-command is used to compute genotype likelihoods from a BAM file. The option --skip-indels is defined to only account for substitutions. 
# The call sub-command creates a VCF file from these genotype likelihoods. The following optional parameters are set:
#	--ploidy 1; defines that only one allele exists in the considered organism, as it is the case for monoploid bacteria (CFT073 reference set with no plasmids).
#	-v; only variant sites will be printed.
#	-m; specifies the usage of the multiallelic-caller (instead of the consensus-caller; the multiallelic-caller is described to overcome 'some limitations' of the consensus-caller, however, these are not further specified, source: http://samtools.github.io/bcftools/bcftools.html#call)
${BCFTOOLS} mpileup --threads ${THREADS} --skip-indels -f ./data/references/CFT073.fasta -O b -o ./temp/SRR8482585-bcftools.bcf ./temp/SRR8482585.bam
${BCFTOOLS} call --threads ${THREADS} --ploidy 1 -m -v -o ./results/SNVdiscovery/CFT073/bcftools/SRR8482585-bcftools.vcf ./temp/SRR8482585-bcftools.bcf
rm ./temp/SRR8482585-bcftools.bcf
${BCFTOOLS} mpileup --threads ${THREADS} --skip-indels -f ./data/references/CFT073.fasta -O b -o ./temp/SRR8494940-bcftools.bcf ./temp/SRR8494940.bam
${BCFTOOLS} call --threads ${THREADS} --ploidy 1 -m -v -o ./results/SNVdiscovery/CFT073/bcftools/SRR8494940-bcftools.vcf ./temp/SRR8494940-bcftools.bcf
rm ./temp/SRR8494940-bcftools.bcf
${BCFTOOLS} mpileup --threads ${THREADS} --skip-indels -f ./data/references/CFT073.fasta -O b -o ./temp/SRR8494940~SRR8482585-bcftools.bcf ./temp/SRR8494940.bam ./temp/SRR8482585.bam
${BCFTOOLS} call --threads ${THREADS}     --ploidy 1 -m -v -o ./results/SNVdiscovery/CFT073/bcftools/SRR8494940~SRR8482585-bcftools.vcf ./temp/SRR8494940~SRR8482585-bcftools.bcf
rm ./temp/SRR8494940~SRR8482585-bcftools.bcf
rm ./temp/SRR8494940.bam
rm ./temp/SRR8482585.bam