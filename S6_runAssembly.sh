#!/bin/bash

## Author: Simon Hackl, 28.02.2021
## Project: Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads
## This is a generalized script to run the assemblers unicycler, flye, canu, raven and haslr on one of the reference data sets or one of its sub samples.
## To use the script, two necessary and one optional input parameters have to be specified.
## - The first defines the assembler to use and should be one of [unicycler~hybrid, haslr, unicycler~longread, flye, canu, raven]
## - The second defines the reference set and should be one of [CFT073, MGH78578, RN4220~guppy3210, RN4220~guppy4011] 
## - The optional third specifies a sub sample coverage and should be one of [200X, 150X, 100X, 80X, 60X, 40X, 20X, 15X, 10X, 8X, 6X, 4X, 2X, 1X]
## E.g. calling `./S6_runAssemblies.sh unicycler~hybrid CFT073 100X` would run unicycler in hybrid assembly mode on the 100X sub sample of the CFT073 reference set.
## NOTE: All assemblers are run with default parameters. If any exceptions from this note are made they are documented in this script.
## 
## The output files of each run can be kept by switching the CLEAN variable below to yes. If this is done, the result files will be stored at
## './results/assembly-out/<SAMPLE-IDENTIFIER>/<ASSEMBLER-IDENTIFIER>-<OPTIONAL-COVERAGE>'.
## The final assembly will be stored at './results/assemblies/<SAMPLE-IDENTIFIER>~<OPTIONAL-COVERAGE>-<ASSEMBLER-IDENTIFIER>'.
## If any separate polishing tools are run, these will be appended to the respective directory and file names with an underscore.

# Generate directories to store the output files.
mkdir ./results
mkdir ./results/assemblies
mkdir ./results/assembly-out
mkdir ./temp

# Specify the number of threads to use.
THREADS=6

# If CLEAN is set to 'true' only the final assembly (FASTA) will be stored.
CLEAN=no

# The following variables are used to access tools used in this script.
# Change these to point to the executables of the respecitve programs if necessary.
UNICYCLER="./tools/unicycler/unicycler-runner.py"
HASLR="./tools/haslr/bin/haslr.py"
FLYE="flye"
RAVEN="raven"
CANU="./tools/canu/bin/canu"
PILON="./tools/pilon/pilon-1.23.jar"
SPADES="./tools/spades/spades.py"
MEDAKA="./tools/medaka"

# The following variables are used only for correct naming of output directories.
ASSEMBLY_TAG=""
SAMPLE_TAG=""

# The following cases define the read sets to use, genome size and naming variables (SAMPLE, *_TAG) depending on the specified sample.
# NOTE: Some assemblers expect the genome size for working, for CFT073 and MGH78578 the given length of the genome was rounded off to 5.2 million and 5.3 million base pairs, respectively.
# As the RN4220 reference genome is highly fragmented it was decided to use 2.8 million base pairs, the median length of 12,401 genome assemblies of S. aureus 
# (source: https://www.ncbi.nlm.nih.gov/genome/?term=Staphylococcus%20aureus[Organism]&cmd=DetailsSearch), as the genome lenght for the RN4220 reference set.
if [[ "$2" == CFT073 ]]
then
	mkdir ./results/assembly-out/CFT073
	SAMPLE=CFT073
	SIZE=5.2
	READS_S1=./data/reads/${SAMPLE}/SRR8482585_1.fastq
	READS_S2=./data/reads/${SAMPLE}/SRR8482585_2.fastq
	if [[ "$3" == *X ]]
	then
		READS_L=./data/reads/${SAMPLE}/SRR8494940-t-${3}.fastq
		ASSEMBLY_TAG=~${3}
		SAMPLE_TAG=-${3}
	else
		READS_L=./data/reads/${SAMPLE}/SRR8494940-t.fastq
	fi
fi
if [[ "$2" == MGH78578 ]]
then
	mkdir ./results/assembly-out/MGH78578
	SAMPLE=MGH78578
	SIZE=5.3
	READS_S1="./data/reads/${SAMPLE}/SRR8482567_1.fastq"
	READS_S2="./data/reads/${SAMPLE}/SRR8482567_2.fastq"
	READS_L="./data/reads/${SAMPLE}/SRR8494915-t.fastq"
fi
if [[ "$2" == RN4220~guppy4011 ]]
then
	mkdir ./results/assembly-out/RN4220~guppy4011
	SAMPLE=RN4220
	SIZE=2.8
	READS_S1=./data/reads/${SAMPLE}/QNFLR056AF_1.fastq
	READS_S2=./data/reads/${SAMPLE}/QNFLR056AF_2.fastq
	if [[ "$3" == *X ]]
	then
		READS_L=./data/reads/${SAMPLE}/QNFLR049AW~guppy4011-t-${3}.fastq
		ASSEMBLY_TAG=~${3}
		SAMPLE_TAG=-${3}
	else
		READS_L=./data/reads/${SAMPLE}/QNFLR049AW~guppy4011-t.fastq
	fi
fi
if [[ "$2" == RN4220~guppy3210 ]]
then
	mkdir ./results/assembly-out/RN4220~guppy3210
	SAMPLE=RN4220
	SIZE=2.8
	READS_S1=./data/reads/${SAMPLE}/QNFLR056AF_1.fastq
	READS_S2=./data/reads/${SAMPLE}/QNFLR056AF_2.fastq
	if [[ "$3" == *X ]]
	then
		READS_L=./data/reads/${SAMPLE}/QNFLR049AW~guppy3210-t-${3}.fastq
		ASSEMBLY_TAG=~${3}
		SAMPLE_TAG=-${3}
	else
		READS_L=./data/reads/${SAMPLE}/QNFLR049AW~guppy3210-t.fastq
	fi
fi

# The following cases will run the specified assembler.
if [[ "$1" == unicycler~hybrid ]]
then
	mkdir ./results/assembly-out/${2}/unicycler~hybrid${SAMPLE_TAG}
	$UNICYCLER -t $THREADS -l $READS_L -1 $READS_S1 -2 $READS_S2 -o ./results/assembly-out/${2}/unicycler~hybrid${SAMPLE_TAG} --pilon_path $PILON --spades_path $SPADES
	cp ./results/assembly-out/${2}/unicycler~hybrid${SAMPLE_TAG}/assembly.fasta ./results/assemblies/${SAMPLE}${ASSEMBLY_TAG}-unicycler~hybrid.fasta
fi
if [[ "$1" == haslr ]]
then
	mkdir ./results/assembly-out/${2}/haslr${SAMPLE_TAG}
	$HASLR -g ${SIZE}m -t $THREADS -x nanopore -l $READS_L -s $READS_S1 $READS_S2 -o ./results/assembly-out/${2}/haslr${SAMPLE_TAG}
	cp ./results/assembly-out/${2}/haslr${SAMPLE_TAG}/asm_contigs_*_sim0.85/asm.final.fa ./results/assemblies/${SAMPLE}${ASSEMBLY_TAG}-haslr.fasta
fi
if [[ "$1" == unicycler~longread ]]
then
	mkdir ./results/assembly-out/${2}/unicycler~longread${SAMPLE_TAG}
	$UNICYCLER -t $THREADS -l $READS_L -o ./results/assembly-out/${2}/unicycler~longread${SAMPLE_TAG} --pilon_path $PILON --spades_path $SPADES
	cp ./results/assembly-out/${2}/unicycler~longread${SAMPLE_TAG}/assembly.fasta ./results/assemblies/${SAMPLE}${ASSEMBLY_TAG}-unicycler~longread.fasta
fi
if [[ "$1" == flye ]]
then
	mkdir ./results/assembly-out/${2}/flye${SAMPLE_TAG}
	if [[ "$2" == RN4220~guppy4011 ]]
	then
		# Flye failed the assembly for the RN4220 longreads basecalled with Guppy version 4.0.1.1 unless the --asm-coverage parameter was set to 50. 
		$FLYE -t $THREADS --nano-raw $READS_L -o ./results/assembly-out/${2}/flye${SAMPLE_TAG} -g ${SIZE}m --asm-coverage 50
	else
		$FLYE -t $THREADS --nano-raw $READS_L -o ./results/assembly-out/${2}/flye${SAMPLE_TAG} -g ${SIZE}m
	fi
	cp ./results/assembly-out/${2}/flye${SAMPLE_TAG}/assembly.fasta ./results/assemblies/${SAMPLE}${ASSEMBLY_TAG}-flye.fasta
fi
if [[ "$1" == raven ]]
then
	mkdir ./results/assembly-out/${2}/raven${SAMPLE_TAG}
	$RAVEN -t $THREADS $READS_L >> ./results/assembly-out/${2}/raven${SAMPLE_TAG}/log.txt
	#Raven prints the assembled contigs to stdout. These are stored in a temp. log.txt file. The output is converted to a proper FASTA file (limiting line length to 80 symbols by using the `fold` command).
	fold -w80 ./results/assembly-out/${2}/raven${SAMPLE_TAG}/log.txt > ./results/assemblies/${SAMPLE}${ASSEMBLY_TAG}-raven.fasta
	rm -r ./results/assembly-out/${2}/raven${SAMPLE_TAG}
fi
if [[ "$1" == canu ]]
then
	mkdir ./results/assembly-out/${2}/canu${SAMPLE_TAG}
	$CANU -p $SAMPLE -d ./results/assembly-out/${2}/canu${SAMPLE_TAG} -nanopore-raw $READS_L genomesize=${SIZE}m -minInputCoverage=0 -stopOnLowCoverage=0
	source ${MEDAKA}/bin/activate
	medaka_consensus -i $READS_L -d ./results/assembly-out/${2}/canu${SAMPLE_TAG}/assembly.fasta -o ./results/assembly-out/${2}/canu${SAMPLE_TAG} -t $THREADS
	cp ./results/assembly-out/${2}/canu${SAMPLE_TAG}/assembly.fasta ./results/assemblies/${SAMPLE}${ASSEMBLY_TAG}-canu.fasta
	cp ./results/assembly-out/${2}/canu${SAMPLE_TAG}/consensus.fasta ./results/assemblies/${SAMPLE}${ASSEMBLY_TAG}-canu_medaka.fasta
fi

# If CLEAN was set to yes, remove the files output by the assemblers and store only the final assembly (FASTA).
if [[ "$CLEAN" == yes ]]
then
	rm -r ./results/assembly-out/${2}/${1}${SAMPLE_TAG}
fi