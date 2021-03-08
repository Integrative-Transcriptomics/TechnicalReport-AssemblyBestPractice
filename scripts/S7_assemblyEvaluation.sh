#!/bin/bash

## Author: Simon Hackl, 28.02.2021
## Project: Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads
## This is a generalized script to run quast for assembly evaluation.
## To use the script, one necessary and one optional input parameter has to be specified.
## - The first defines the reference set for which the assemblies are evaluated, i.e. all computed assemblies of one reference set are evaluated jointly, and should be one of [CFT073, MGH78578, RN4220~guppy3210, RN4220~guppy4011]
## - The optional second specifies a sub sample coverage and should be one of [200X, 150X, 100X, 80X, 60X, 40X, 20X, 15X, 10X, 8X, 6X, 4X, 2X, 1X]
## E.g. calling `./S7_assemblyEvaluation.sh CFT073 100X` would evaluate all assemblies conducted for the CFT073 reference set with the 100X coverage sub sample.
##
## NOTE: After Quast has been installed, it allows to alter the features that are written to the report output files; for this project Quast was installed locally at
## ./tools/quast and the Python script ./tools/quast/quast_libs/reporting.py was edited to yield the following content in the order variable:
#    order = [NAME, CONTIGS__FOR_THRESHOLDS, TOTALLENS__FOR_THRESHOLDS, CONTIGS, LARGCONTIG, TOTALLEN, REFLEN, ESTREFLEN, GC, REFGC,
#             N50, NG50, N75, NG75, L50, LG50, L75, LG75,
#             TOTAL_READS, LEFT_READS, RIGHT_READS,
#             MAPPED_READS_PCNT, REF_MAPPED_READS_PCNT,
#             PROPERLY_PAIRED_READS_PCNT, REF_PROPERLY_PAIRED_READS_PCNT,
#             DEPTH, REF_DEPTH, COVERAGE_1X_THRESHOLD, REF_COVERAGE_1X_THRESHOLD,
#             LARGE_MIS_EXTENSIVE, MISASSEMBL, MIS_RELOCATION, MIS_TRANSLOCATION, MIS_INVERTION, MISCONTIGS, MISCONTIGSBASES,
#             MISLOCAL, MIS_SCAFFOLDS_GAP, MIS_LOCAL_SCAFFOLDS_GAP,
#             STRUCT_VARIATIONS, POTENTIAL_MGE, UNALIGNED_MISASSEMBLED_CTGS,
#             UNALIGNED, UNALIGNEDBASES, MAPPEDGENOME, DUPLICATION_RATIO, AVE_READ_SUPPORT,
#             UNCALLED_PERCENT, SUBSERROR, MISMATCHES, INDELSERROR, INDELS, MIS_SHORT_INDELS, MIS_LONG_INDELS, GENES, OPERONS,
#             BUSCO_COMPLETE, BUSCO_PART,
#             PREDICTED_GENES_UNIQUE, PREDICTED_GENES, RNA_GENES,
#             LARGALIGN, TOTAL_ALIGNED_LEN, NA50, NGA50, NA75, NGA75, LA50, LGA50, LA75, LGA75,
#             KMER_COMPLETENESS, KMER_CORR_LENGTH, KMER_MIS_LENGTH, KMER_MISASSEMBLIES]
## This will write the total number of mismatches and indels and information about misassemblies into the .txt and .tsv reports what is not done by default.
##
## NOTE: In order to run Quast with the parameters specified below (--gene-finding) a GeneMarkS license key has to be acquired at http://topaz.gatech.edu/GeneMark/license_download.cgi
## and placed at the home directory of the user running this script.
##
## The respective Quast reports will be written to the './results/quast/<SAMPLE-IDENTIFIER>~<OPTIONAL-COVERAGE>' directory.

# Specify the number of threads to use.
THREADS=6

# The following variables are used to access tools used in this script.
# Change these to point to the executables of the respecitve programs if necessary.
QUAST="./tools/quast/quast.py"

# If the optional coverage parameter is set, a TAG variable is used for correct naming.
if [[ "$2" == *X ]]
then
	TAG=~${2}
else
	TAG=""
fi

# The SAMPLE parameter is used for pointing to the correct reference identifier, 
# i.e. it cuts the ~guppyX.X.X.X suffix from the RN4220 reference set names.
SAMPLE="${1%~*}"

# If the sample is RN4220, apply the --fragmented parameter. Thereby, Quast tries to filter out misassembly events caused by a fragmented reference genome.
if [[ "$SAMPLE" == "RN4220" ]]
then
	FRAG=--fragmented
else
	FRAG=""
fi

# Generate directories to store the output files.
mkdir ./results
mkdir ./results/quast
mkdir ./results/quast/${1}${TAG}

# The ASSEMBLIES and LABELS variables will store a list of all relevant assembly (FASTA) files and shortened labels used in the quast reports (i.e. only the name of the assembler), respectively.
# NOTE: In order to run this script all assemblies have to be stored as done with the S6_runAssemblies.sh script.
ASSEMBLIES=$(find ./results/assemblies/${1}${TAG}-*.fasta)
ASSEMBLIES=`printf '%s\n' "${ASSEMBLIES[@]}" | paste -sd ' '`
LABELS=$(find ./results/assemblies/${1}${TAG}-*.fasta | sed 's!.*/!!' | sed 's!.*-!!' | sed 's/.fasta/,/g')
LABELS=`printf '%s\n' "${LABELS[@]::-1}" | paste -sd ''`

# Finally, Quast is run on the specified files.
${QUAST} -o ./results/quast/${1}${TAG} -t $THREADS -r ./data/references/${SAMPLE}.fasta -g ./data/references/${SAMPLE}.gff3 --gene-finding --no-icarus --no-plots $FRAG -l $LABELS $ASSEMBLIES
