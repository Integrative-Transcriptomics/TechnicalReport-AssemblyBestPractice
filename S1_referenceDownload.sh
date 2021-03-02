 #!/bin/bash

## Author: Simon Hackl, 28.02.2021
## Project: Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads
## This script is used to download the CFT073, MGH78578 and RN4220 reference genomes and gene annotations from the NCBI databses.

# Generate directories to store the output files.
mkdir ./data
mkdir ./data/references
mkdir ./data/temp

# `wget` is used on an URL accessing the sequence viewer of the NCBI nuccore datbase with the respective sample IDs. 
# The `sed` command is used to remove empty lines. If necessary, the `cat` command is used to merge multiple filed 
# (i.e. if plasmids are present or the genome is fragmented).

# Download CFT073 reference genome and annotation.
wget -O ./data/references/CFT073.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=NC_004431.1"
sed -i '/^$/d' ./data/references/CFT073.fasta
wget -O ./data/references/CFT073.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NC_004431.1"
sed -i '/^$/d' ./data/references/CFT073.gff3

# Download the MGH78578 reference genome. The reference genome contains 5 plasmids.
wget -O ./data/temp/chromosome.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=NC_009648.1"
wget -O ./data/temp/p1.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=NC_009649.1"
wget -O ./data/temp/p2.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=NC_009650.1"
wget -O ./data/temp/p3.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=NC_009651.1"
wget -O ./data/temp/p4.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=NC_009652.1"
wget -O ./data/temp/p5.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=NC_009653.1"
# The plasmids and the 'chromosome' are copied in one FASTA file using the `cat` command.
cat ./data/temp/* > ./data/references/MGH78578.fasta
sed -i '/^$/d' ./data/references/MGH78578.fasta
rm -r ./data/temp
mkdir ./data/temp

# Download the MGH78578 annotation.
wget -O ./data/temp/chromosome.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NC_009648.1"
wget -O ./data/temp/p1.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff33&id=NC_009649.1"
wget -O ./data/temp/p2.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NC_009650.1"
wget -O ./data/temp/p3.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NC_009651.1"
wget -O ./data/temp/p4.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NC_009652.1"
wget -O ./data/temp/p5.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NC_009653.1"
# All six files are again copied in one GFF3 file.
cat ./data/temp/* > ./data/references/MGH78578.gff3
sed -i '/^$/d' ./data/references/MGH78578.gff3
rm -r ./data/temp
mkdir ./data/temp

# Download the RN4220 reference genome. As the reference genome is fragmented in contigs 1 to 118, a for-loop
# is used to access all 118 contigs. The `sed` command is used to merge all downloaded files.
for i in $(seq -f "%03g" 1 118)
do
	wget -O ./data/temp/NZ_AFGU01000$i.1.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=NZ_AFGU01000$i.1"
done
cat ./data/temp/* > ./data/references/RN4220.fasta
sed -i '/^$/d' ./data/references/RN4220.fasta
rm -r ./data/temp
mkdir ./data/temp

# The same procedure is applied to download the RN4220 genome annotation.
for i in $(seq -f "%03g" 1 118)
do
	wget -O ./data/temp/NZ_AFGU01000$i.1.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NZ_AFGU01000$i.1"
done
cat ./data/temp/* > ./data/references/RN4220.gff3
sed -i '/^$/d' ./data/references/RN4220.gff3
rm -r ./data/temp