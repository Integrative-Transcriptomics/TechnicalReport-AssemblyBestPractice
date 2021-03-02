# Researchproject:
### Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads
***
# Project description
The advent of third generation long read sequencing technologies, i.e. Oxford Nanopore Technologies (ONT) devices, greatly expanded the possibilities for de novo assembly in terms of the contiguity of the reconstructed genomes. However, current long read technologies suffer from a not insignificant error rate compared to second generation short read sequencing technologies. Mainly due to new developments in base calling algorithms (the translation of signals from the sequencing device into DNA sequences) and further development of the protein pores used in ONT devices, this error rate is constantly decreasing. Over the last years many new algorithmic approaches have been developed to efficiently use these long reads alone or in combination with short reads, i.e. in a hybrid approach, for de novo assembly or error correction of long reads.

In the course of the project _Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads_ we asked to what extent a quality difference between long read only and hybrid approaches for de novo assembly is noticeable. To answer this, a selection of promising assemblers was chosen based on recently published reviews, a de novo assembly was performed with these and the resulting assemblies were subsequently evaluated according to commonly used criteria.

Another aspect is the cost efficiency of using long reads. Increased read coverage can theoretically correct errors in the long reads, but this comes with an increase of the sequencing cost. We have investigated how different depths of coverage of the long reads affect the quality of the computed assemblies.

Moreover, the applicability of long reads for single nucleotide variation discovery, relative to one reference genome, was studied by using either high qualitative de novo assemblies or the read sets alone.

This repository grants access to the projects report and supplementary data. An extended description of the conducted methods is given and scripts are provided with which the obtained results can be reproduced.
***
# Best practice recommendation
Our findings were summarized in a best practice recommendation.
***
# Extended description of methods and material
### Note on the structure of this repository
Most of the results for this project were obtained by running **Shell** scripts and **Jupyter** (**Python 3**) notebooks. Each of these scripts is provided herein and is named with the prefix _Sx\__ were _x_ is a number indicating the order in which the scripts are ment to be executed (the scripts starting with _Tx\__ were used to ease some manual command calls). Each script comprises a comprehensive documentation on its own and therefore only a top-level description on how and what for the scripts are used is given. In addition manual console commands and tools were run which will be explained in greater detail.

For members of the _Integrated Transcriptomics_ research group a fully executed version of this repository comprising all result files is available at the server _romanov_. 

The most important results, i.e. those which can be shared on github due to their size, are provided in the _./supplementary_files_ directory and were generated mostly manually. By cloning this repository and executing all of the steps described below the directories _./data_ and _./results_ will be generated, featuring the following file catalog:
```
├───data
│   ├───reads 
│   │   └───<SAMPLE-ID>
│   │        // one dir. for each reference set:
│   │        // contains short reads and raw, trimmed
│   │        // and sub sampled long reads.
│   └───references
│        // contains reference genomes (fasta) and
│        // annotations (gff3).
└───results
    ├───assemblies
    │    // contains the final fasta output of each assembler run.
    ├───assembly-out
    │   └───<SAMPLE-ID> 
    │        // one dir. for each reference set (for RN4220
    │        // the basecaller version is appended): stores
    │        // all output files of each assembler run.
    ├───coverage
    │   └───<SAMPLE-ID>
    │        // one dir. for each reference set: contains files (cov)
    │        // storing coverage information for each read set.
    ├───fastqc
    │   └───<SAMPLE-ID>
    │        // one dir. for each reference set: contains reports
    │        // of quality control for each read set (zip and html).
    ├───porechop
    │   └───<SAMPLE-ID>
    │        // one dir. for each reference set: contains files (txt)
    │        // storing information about long read trimming
    ├───quast
    │    // one dir. for each pair of reference set and sub sample
    │    // coverage, contain information about assembly evaluation
    └───SNVdiscovery
        └───CFT073
            ├───bcftools
            ├───dnadiff
            └───mauve
                 // each dir. stores output regarding SNV discovery
                 // with the respective tool
```

### Dependencies and software used
In the following a tabular overview of all software requirements to reproduce the results of this project are given. As some of the tools were accessible via the systems path variables of the system on which the scripts were run, but others were installed locally in a _./tools/_ directory, variables are defined in each script to access the executables of the tools; if necessary these variables can easily be changed in the respective scripts.

| Name | Type | Version | Access | Publication |
|-|-|-|-|-|
| FastQC | read quality control | 0.11.5 | [bioinformatics.babraham.ac.uk/projects/fastqc/](bioinformatics.babraham.ac.uk/projects/fastqc/) | - |
| Porechop | long read trimming | 0.2.4 | [github.com/rrwick/Porechop](github.com/rrwick/Porechop) | - |
| Pilon | short read polishing | 1.23 | [github.com/broadinstitute/pilon](github.com/broadinstitute/pilon) | [doi.org/10.1371/journal.pone.0112963](doi.org/10.1371/journal.pone.0112963) |
| Medaka | long read polishing | 1.2.3 | [github.com/nanoporetech/medaka](github.com/nanoporetech/medaka) | - |
| Unicycler | short-/long read/hybrid assembler | 0.4.9b | [github.com/rrwick/Unicycler](github.com/rrwick/Unicycler) | [doi.org/10.1371/journal.pcbi.1005595](doi.org/10.1371/journal.pcbi.1005595) |
| Flye | long read assembler | 2.7.1-b1590 | [github.com/fenderglass/Flye](github.com/fenderglass/Flye) | [doi.org/10.1038/s41587-019-0072-8](doi.org/10.1038/s41587-019-0072-8) |
| Canu | long read assembler | 2.1.1 | [github.com/marbl/canu](github.com/marbl/canu) | [doi.org/10.1101/gr.215087.116](doi.org/10.1101/gr.215087.116) |
| Raven | long read assembler | 1.2.2 | [github.com/lbcb-sci/raven][github.com/lbcb-sci/raven] | [doi.org/10.1101/2020.08.07.242461](doi.org/10.1101/2020.08.07.242461) |
| Trycycler | long read consensus | 0.3.3 | [github.com/rrwick/Trycycler](github.com/rrwick/Trycycler) | - |
| Haslr | hybrid assembler | 0.8a1 | [github.com/vpc-ccg/haslr](github.com/vpc-ccg/haslr) | [doi.org/10.1101/2020.01.27.921817](doi.org/10.1101/2020.01.27.921817) |
| quast | assembly evaluation | 5.0.2 | [github.com/ablab/quast](github.com/ablab/quast) | [doi.org/10.1093/bioinformatics/bty266](doi.org/10.1093/bioinformatics/bty266) |
| rasusa | read subsampling | 0.3.0 | [github.com/mbhall88/rasusa](github.com/mbhall88/rasusa) | [doi.org/10.5281/zenodo.3731394](doi.org/10.5281/zenodo.3731394) |
| samtools | sequencing data manipulation | 1.11 | [github.com/samtools/samtools](github.com/samtools/samtools) | [doi.org/10.1093/bioinformatics/btp352](doi.org/10.1093/bioinformatics/btp352) |
| bcftools | genomic variant discovery (from reads) | 1.11-21-g0987715 | [github.com/samtools/bcftools](github.com/samtools/bcftools) | [dx.doi.org/10.1093%2Fbioinformatics%2Fbtr509](dx.doi.org/10.1093%2Fbioinformatics%2Fbtr509) |
| DNAdiff | genomic variant discovery (from genome alignment) | 1.3 | [dargithub.com/marbl/MUMmer3](github.com/marbl/MUMmer3) | [doi.org/10.1186/gb-2004-5-2-r12](doi.org/10.1186/gb-2004-5-2-r12) |
| mauve | genomic variant discovery (from genome alignment) | 2.4.0 | [darlinglab.org/mauve/mauve.html](darlinglab.org/mauve/mauve.html) | [doi.org/10.1371/journal.pone.0011147](doi.org/10.1371/journal.pone.0011147) |
| bwa | mapping/alignment | 0.7.17-r1188 | [github.com/lh3/bwa](github.com/lh3/bwa) | [arxiv.org/abs/1303.3997](arxiv.org/abs/1303.3997) |
| minimap2 | mapping/alignment | 2.17-r974-dirty | [github.com/lh3/minimap2](github.com/lh3/minimap2) | [doi.org/10.1093/bioinformatics/bty191](doi.org/10.1093/bioinformatics/bty191) |

### Download of reference genomes and annotations
By running `S1_referenceDownload.sh` the reference genomes and gene annotations of these are downloaded and stored in the _`_./data/references/_`_ directory as _.fasta_ and _.gff3_ files, named after the reference set identifiers. The three references include:
- _Escherichia coli_ strain CFT073, named as CFT073
- _Klebsiella pneumoniae_ strain MGH78578, named as MGH78578 (includes five plasmids)
- _Staphylococcus aureus_ strain RN4220, named as RN4220 (fragmented genome of 118 contigs)

The script basically calls `wget` on URLs referring to the NCBI sequence viewer with the database specified as nuccore, the report specified as fasta or gff3 and the id specified as the NC and NZ identifiers of the respective references. More information on the identifiers and their publishers can be obtained from the project report. If a reference genome consists of more than one file (i.e. due to plasmids or multiple contigs), all the respective files are first stored in a _./temp_ directory and then combined into one file using the `cat` command. Finally `sed -i '/^$/d'` is called on each file to remove empty lines.

### Download of public read data sets
By running `S2_readDownload.sh` the publicy available read sets of the CFT073 and MGH78578 reference sets are downloaded into the _./data/reads/CFT073_ and _./data/reads/MGH78578_ directory, respectively. The script calls `fastq-dump` on the respective SRR identifiers of the read sets. The exact identifiers and origin of the reads are described in the project report. For the short reads the additional parameter `--split-files` was set, as these are paired-end reads, in order to obtain two separate files (containing the forward- and reverse-strand reads). After executing the following files will be accessible
- SRR8494940.fastq (CFT073 long reads)
- SRR8482585_1.fastq and SRR8482585_2.fastq (CFT073 short reads)
- SRR8494915.fastq (MGH78578 long reads)
- SRR8482567_1.fastq and SRR8482567_2.fastq (MGH78578 short reads)

The read sets of the RN4220 reference are not publicy available and are located on a server of the Universität Tübingen. The respective files were manually copied and stored as follows in the project directory _./data/reads/RN4220_
- QNFLR056AF_1.fastq and QNFLR056AF_2.fastq (RN4220 short reads)
- QNFLR049AW~guppy3210.fastq (RN4220 long reads, basecalled with Guppy 3.2.1.0)
- QNFLR049AW~guppy4011.fastq (RN4220 long reads, basecalled with Guppy 4.0.1.1)

### Quality control of the read sets
### Downsampling of the read sets
### Conducting the de novo assemblies 
### Evaluating the de novo assemblies
### Discovery of single nucleotide variations 