# Researchproject:
### Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads
***
# Project description
The advent of third generation long read sequencing technologies, i.e. Oxford Nanopore Technologies (ONT) devices, greatly expanded the possibilities for de novo assembly in terms of the contiguity of the reconstructed genomes. However, current long read technologies suffer from a not insignificant error rate compared to second generation short read sequencing technologies. Mainly due to new developments in base calling algorithms (the translation of signals from the sequencing device into DNA read sequences) and advancements of the protein pores used in ONT devices, this error rate is constantly decreasing. Over the last years many new algorithmic approaches have been developed to efficiently use these long reads alone or in combination with short reads, i.e. in a hybrid approach, for de novo assembly or error correction of long reads.

In the course of the project _Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads_ we asked to what extent a quality difference between long read only and hybrid approaches for de novo assembly is noticeable. To answer this, a selection of promising assemblers was chosen based on recently published reviews, a de novo assembly was performed with these and the resulting assemblies were subsequently evaluated according to conventional criteria.

Another aspect is the cost efficiency of using long reads. Increased read coverage can theoretically be used to correct errors in the long reads, but this comes at an increase of the sequencing cost. We have investigated how different depths of coverage of the long reads affect the quality of assemblies.

Moreover, the applicability of long reads for single nucleotide variation discovery, relative to one reference genome, was studied by using either high qualitative de novo assemblies or the reads alone.

This repository grants access to the projects report and supplementary data. An extended description of the conducted methods is given and scripts are provided with which the obtained results can be reproduced.
***
# Best practice recommendation
Our findings were summarized in a best practice recommendation.
***
# Extended description of methods and material
### Note on the structure of this repository
Most of the results for this project were obtained by running **Shell** scripts and **Jupyter** (**Python 3**) notebooks. Each of these scripts is provided herein and is named with the prefix _Sx\__ were _x_ is a number indicating the order in which the scripts are ment to be executed (the scripts starting with _Tx\__ were used to ease some manual command calls). Each script comprises a documentation on its own and therefore only a top-level description on how and what for the scripts are used is given here. In addition to these scripts, manual console commands and tools were run which will be explained in more detail.

For members of the _Integrated Transcriptomics_ research group a fully executed version of this repository comprising all result files is available at the server _romanov_. 

The most important results, i.e. those which can be shared on github due to their size, are provided in the _./supplementary_files_ directory and were collated manually. By cloning this repository and executing all of the steps described below the directories _./data_ and _./results_ will be generated, featuring the following file catalog:
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
    │        // one dir. for each reference set: contains files (.cov)
    │        // storing coverage information of each read set.
    ├───fastqc
    │   └───<SAMPLE-ID>
    │        // one dir. for each reference set: contains quality
    │        // control reports for each read set (.zip and .html).
    ├───porechop
    │   └───<SAMPLE-ID>
    │        // one dir. for each reference set: contains reports (.txt)
    │        // of long read trimming
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
| FastQC | read quality control | 0.11.5 | [bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | - |
| Porechop | long read trimming | 0.2.4 | [github.com/rrwick/Porechop](https://www.github.com/rrwick/Porechop) | - |
| Pilon | short read polishing | 1.23 | [github.com/broadinstitute/pilon](https://www.github.com/broadinstitute/pilon) | [doi.org/10.1371/journal.pone.0112963](https://www.doi.org/10.1371/journal.pone.0112963) |
| Medaka | long read polishing | 1.2.3 | [github.com/nanoporetech/medaka](https://www.github.com/nanoporetech/medaka) | - |
| Unicycler | short-/long read/hybrid assembler | 0.4.9b | [github.com/rrwick/Unicycler](https://www.github.com/rrwick/Unicycler) | [doi.org/10.1371/journal.pcbi.1005595](https://www.doi.org/10.1371/journal.pcbi.1005595) |
| Flye | long read assembler | 2.7.1-b1590 | [github.com/fenderglass/Flye](https://www.github.com/fenderglass/Flye) | [doi.org/10.1038/s41587-019-0072-8](https://www.doi.org/10.1038/s41587-019-0072-8) |
| Canu | long read assembler | 2.1.1 | [github.com/marbl/canu](https://www.github.com/marbl/canu) | [doi.org/10.1101/gr.215087.116](https://www.doi.org/10.1101/gr.215087.116) |
| Raven | long read assembler | 1.2.2 | [github.com/lbcb-sci/raven](https://www.github.com/lbcb-sci/raven) | [doi.org/10.1101/2020.08.07.242461](https://www.doi.org/10.1101/2020.08.07.242461) |
| Trycycler | long read consensus | 0.3.3 | [github.com/rrwick/Trycycler](https://www.github.com/rrwick/Trycycler) | - |
| Haslr | hybrid assembler | 0.8a1 | [github.com/vpc-ccg/haslr](https://www.github.com/vpc-ccg/haslr) | [doi.org/10.1101/2020.01.27.921817](https://www.doi.org/10.1101/2020.01.27.921817) |
| quast | assembly evaluation | 5.0.2 | [github.com/ablab/quast](https://www.github.com/ablab/quast) | [doi.org/10.1093/bioinformatics/bty266](https://www.doi.org/10.1093/bioinformatics/bty266) |
| rasusa | read subsampling | 0.3.0 | [github.com/mbhall88/rasusa](https://www.github.com/mbhall88/rasusa) | [doi.org/10.5281/zenodo.3731394](https://www.doi.org/10.5281/zenodo.3731394) |
| samtools | sequencing data manipulation | 1.11 | [github.com/samtools/samtools](https://www.github.com/samtools/samtools) | [doi.org/10.1093/bioinformatics/btp352](https://www.doi.org/10.1093/bioinformatics/btp352) |
| bcftools | genomic variant discovery (from reads) | 1.11-21-g0987715 | [github.com/samtools/bcftools](https://www.github.com/samtools/bcftools) | [dx.doi.org/10.1093%2Fbioinformatics%2Fbtr509](https://www.dx.doi.org/10.1093%2Fbioinformatics%2Fbtr509) |
| DNAdiff | genomic variant discovery (from genome alignment) | 1.3 | [dargithub.com/marbl/MUMmer3](https://www.github.com/marbl/MUMmer3) | [doi.org/10.1186/gb-2004-5-2-r12](https://www.doi.org/10.1186/gb-2004-5-2-r12) |
| mauve | genomic variant discovery (from genome alignment) | 2.4.0 | [darlinglab.org/mauve/mauve.html](https://www.darlinglab.org/mauve/mauve.html) | [doi.org/10.1371/journal.pone.0011147](https://www.doi.org/10.1371/journal.pone.0011147) |
| bwa | mapping/alignment | 0.7.17-r1188 | [github.com/lh3/bwa](https://www.github.com/lh3/bwa) | [arxiv.org/abs/1303.3997](https://www.arxiv.org/abs/1303.3997) |
| minimap2 | mapping/alignment | 2.17-r974-dirty | [github.com/lh3/minimap2](https://www.github.com/lh3/minimap2) | [doi.org/10.1093/bioinformatics/bty191](https://www.doi.org/10.1093/bioinformatics/bty191) |

### Download of reference genomes and annotations
The reference genomes and respective gene annotations were downloaded by running (`S1_referenceDownload.sh`)[./scripts/S1_referenceDownload.sh]. Thereby the gathered files were stored in the _./data/references/_ directory as _.fasta_ and _.gff3_ files, named after the reference set identifiers. The three references include:
- _Escherichia coli_ strain CFT073, named as CFT073
- _Klebsiella pneumoniae_ strain MGH78578, named as MGH78578 (includes five plasmids)
- _Staphylococcus aureus_ strain RN4220, named as RN4220 (fragmented genome of 118 contigs)

The tag \<SAMPLE-ID\> will be used in the following to refer to one of CFT073, MGH78578 or RN4220. The script basically calls `wget` on URLs referring to the NCBI sequence viewer with the database specified as nuccore, the report specified as fasta or gff3 and the id specified as the NC and NZ identifiers of the respective references. More information on the identifiers and their publishers can be obtained from the project report. If a reference genome consists of more than one file (i.e. due to plasmids or multiple contigs), all the respective files are first stored in a _./temp_ directory and then combined into one file using the `cat` command. Finally `sed -i '/^$/d'` is called on each file to remove empty lines.

### Download of public read data sets
By running (`S2_readDownload.sh`)[./scripts/S2_readDownload.sh] the publicy available read sets of the CFT073 and MGH78578 sets were downloaded into the _./data/reads/CFT073_ and _./data/reads/MGH78578_ directory, respectively. The script calls `fastq-dump` on the respective SRR identifiers of the read sets. The BioSample accession identifiers of the single read sets can be obtained from the project report. The run identifiers were taken from the respective database entries. For the short reads the additional parameter `--split-files` was set, since these are paired-end reads, in order to obtain two separate files (containing the forward- and reverse-strand reads). After executing the following files will be accessible in the directories described above:
- SRR8494940.fastq (CFT073 long reads)
- SRR8482585_1.fastq and SRR8482585_2.fastq (CFT073 short reads)
- SRR8494915.fastq (MGH78578 long reads)
- SRR8482567_1.fastq and SRR8482567_2.fastq (MGH78578 short reads)

The read sets of the RN4220 reference are not publicy available and are located on a server of the Universität Tübingen. The respective files were manually copied and stored as follows in the project directory _./data/reads/RN4220_:
- QNFLR056AF_1.fastq and QNFLR056AF_2.fastq (RN4220 short reads)
- QNFLR049AW~guppy3210.fastq (RN4220 long reads, basecalled with Guppy 3.2.1.0)
- QNFLR049AW~guppy4011.fastq (RN4220 long reads, basecalled with Guppy 4.0.1.1)

### Quality control of the read sets
A quality control step of all downloaded red sets was executed by running (`S3_readQualityControl.sh`)[./scripts/S3_readQualityControl.sh]. Therewith, the following three steps are executed:
- All long read sets are processed with Porechop for adapter removal. The trimmed long read _.fastq_ files were named with a _-t_ suffix in their filenames and stored in the same directory as the untrimmed reads. The Porechop reports are stored in the _./results/porechop/\<SAMPLE-ID\>/Porechop-\<READ-FILENAME\>.log_ directories and are moreover accessible via (supplementary file 1)[./supplementary_files/F1_porechop.zip].
- To control the quality of the input reads before conducting the assemblies, FastQC was run on all read sets. One _.html_ and one _.zip_ file is generated for each long read and both short read files as _./results/fastqc/\<SAMPLE-ID\>/\<READ-FILENAME\>_fastqc(.zip or .html)_. The reports are in addition accessible via (supplementary file 2)[./supplementary_files/F2_fastqc.zip].
- Lastly, to check for the depth and breadth of coverage across the reference genome the reads were first mapped to the respective genomes using `minimap2` (specifying the `-x map-ont` and `-x sr` parameter for long and short read mapping, respectively, and the `-a` parameter in order to produce output in SAM format). Next, the SAM files were sorted and converted to BAM files using `samtools sort` with the `-O bam` parameter set. The `samtools depth` command was run on all these BAM files, which basically outputs the per position depth of coverage, with the `-a` parameters set to include positions with zero coverage into the report. The resulting files are stored as _./results/coverage/\<SAMPLE-ID\>/\<READ-FILENAME\>.cov (note that one common file for the short reads named with a \_1\_2 suffix is created).

### Sub-sampling of the read sets
To assess the effect of decreasing long read coverage depth on the assembly quality the tool Rasusa was used to generate sub samples with expected coverages of 200X, 150X, 100X, 80X, 60X, 40X, 20X, 15X, 10X, 8X, 6X, 4X, 2X and 1X of the CFT073 and RN4220 long reads. To automate this process the script (`S4_readSubsampling.sh`)[./scripts/S4_readSubsampling.sh] was used. The script conducts the following steps:
- Run Rasusa on the respective long read files with the pre-defined coverage and an estimated genome length. For the genome length of CFT073 the length of the reference genome was rounded of to 5.2 million. For the RN4220 reference genome it was decided to use 2.8 million, the median length of 12,401 genome assemblies of S. aureus  (source: [www.ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov/genome/?term=Staphylococcus%20aureus[Organism]&cmd=DetailsSearch)), as the genomes length due to its fragmentation. Each sub sample will be stored as _./data/reads/\<SAMPLE-ID\>/\<READ-FILENAME\>-\<EXPECTED-COVERAGE\>.fastq_.
- The third step of the section before is executed.

### Analysing the read set coverage and quality
To summarize the per read sequence qualities from the FastQC reports and to analyze the .cov files the Python notebook (`S5_coverageAnalysis.ipynb`)[./scripts/S5_coverageAnalysis.ipynb] was used. The notebook was written with Python version 3.9.1 and documents each of the conducted steps, which basically comprise the following:
- The mean per read Phred quality score from the FastQC report, the mean depth of coverage and percentage of genomic positions with zero coverage as well as the expected percentage of genomic positions with at least coverage depth one according to the Lander-Waterman model are collected and written to (supplementary file 3)[./supplementary_files/F3_readStatistics.csv] (note: the mean per read quality is only reported for non sub-sample reads).
- Each .cov file is used to generate a coverage profile: These coverage profiles split all genomic positions into a set of chunks of fixed length. For each chunk the mean depth of coverage is plotted.

### Running the de novo assemblies
The next step comprised running the de novo assemblies on the different read sets. To simplify this process the script (`S6_runAssembly.sh`)[./scripts/S6_runAssembly.sh] was used. The script requiers the definition of three parameters:
1. the assembler to use, one of unicycler~hybrid, haslr, unicycler~longread, flye, canu or raven
2. the reference set (i.e. the reads) to use, one of CFT073, MGH78578, RN4220~guppy3210 or RN4220~guppy4011
3. (optional) a sub sample coverage, one of 200X, 150X, 100X, 80X, 60X, 40X, 20X, 15X, 10X, 8X, 6X, 4X, 2X or 1X

The script will then run the assembler on the defined reference set and (optional) on an sub sample of the long reads. Note that if MGH78578 and a sub sample coverage is defined, no assembly will be conducted as no sub samples of the MGH78578 reference set were generated. The default parameters of the assemblers were used, with the following exceptions:
- For Flye and the RN4220 reads basecalled with Guppy 4.0.1.1 the parameter `--asm-coverage 50` had to be set in order to the disjointig assembly step being completed. The parameter will reduce the reads used for this initial step to a coverage of 50X.
- For Canu the parameters `-minInputCoverage=0 -stopOnLowCoverage=0` were set in order to prevent the abortion of the assembly process for low coverage sub samples.
- HASLR, Flye and Canu require the specification of an approximate genome length. For these assemblers the values `5.2m`, `5.3m` and `2.8m` are used for the CFT073, MGH78578 and RN4220 reference sets, respectively.
- Raven does not produce any log files but print its output to `stdout`. For this reason only the final assemblies of Raven are stored.

The log files of each assembler are stored at _./results/assembly-out/\<SAMPLE-ID\>/\<ASSEMBLER-ID\>-\<OPTIONAL-COVERAGE\>/_ were \<ASSEMBLER-ID\>, \<SAMPLE-ID\> and \<OPTIONAL-COVERAGE\> are the first, second and third input parameters of the mentioned script, respectively. If no third parameter is given the directory name _...\<ASSEMBER-ID\>_ alone is chosen. The final assemblies (in fasta format) are copied to _./results/assemblies/\<SAMPLE-ID\>~\<OPTIONAL-COVERAGE\>-\<ASSEMBLER-ID\>.fasta_.

The script was run for all available assemblers for the non-sub sample reference sets, all available assemblers for the CFT073 sub samples and Unicycler in hybrid mode for the RN4220 sub samples. After execution the resulting fasta files were checked and empty fasta files were deleted.

### Running Trycycler
Trycycler is a tool to compute consensus sequences from (long read) assemblies and is executed in six steps which are illustrated [here](https://github.com/rrwick/Trycycler/wiki/Illustrated-pipeline-overview) and described in more detail [here](https://github.com/rrwick/Trycycler/wiki/How-to-run-Trycycler). For this project, Trycycler was executed manually. The output directories were chosen in the same manner as for the other assemblies with the assembler ID _"trycycler"_. Trycycler was run on all available long read assemblies (i.e. from Unicycler, Raven, Flye and Canu) of the non sub-sample reference sets and the sub-sample assemblies of the CFT073 reference reads. During the single steps, the following adjustments were made:
- For CFT073:
    - The Canu assembly polished with Medaka could not be circularized during the `reconcile` step. Thus, the unpolished Canu assembly was used.
    - A contig of the Canu assembly in cluster one required additional trimming for circularization during the `reconcile` step. For this reason the parameter `--max_trim_seq 100000` was set for the `reconcile` step.
    - The contigs of  the Canu and Raven assemblies in cluster one showed indels of about 800 base pairs length. The parameter `--max_index_size 800` was set for the `reconcile` step to resolve this issue.
- For CFT073 sub samples:
    - The Canu assembly showed circularization issues for the CFT073 sub samples and was removed before the `reconcile` step.
    - For low coverage sub samples indel errors emerged which prevented circularization. To resolve this issue the parameter `--max_index_size 800` was set for the `reconcile` step.
    - The Raven assembly showed circularization issues for the CFT073 sub sample with 20X coverage and was removed before the `reconcile` step.
- For MGH78578:
    - The Canu assembly polished with Medaka could not be circularized during the `reconcile` step. Thus, the unpolished Canu assembly was used.
    - The third cluster generated by the `cluster` step comprised only one short (82,712 base pairs) contig from the Unicycler assembly and was deleted.
    - The second cluster generated by the `cluster` step comprised one contig of the Canu assembly that was much longer than the other contigs in the cluster (256,899 base pairs versus approx. 186,000 base pairs of the other contigs). The contig was therefore deleted.
    - A contig of the Canu assembly in cluster one required the trimming of 76,343 base pairs for circularization during the `reconcile` step. For this reason the parameter `--max_trim_seq 100000` was set for the `reconcile` step.
- For RN4220 (Guppy version 4.0.1.1):
    - A contig of the Canu assembly in cluster one required the trimming of 69,107 base pairs for circularization during the `reconcile` step. For this reason the parameter `--max_trim_seq 70000` was set for the `reconcile` step.
- For RN4220 (Guppy version 3.2.1.0):
    - The Canu assembly could not be circularized during the `reconcile` step and was removed.

After running the `consensus` step the final assemblies were copied to _./results/assemblies/\<SAMPLE-ID\>~\<OPTIONAL-COVERAGE\>-trycycler.fasta_. These were further processed by polishing with Medaka and Pilon. The resulting polished assemblies were copied to _./results/assemblies/\<SAMPLE-ID\>~\<OPTIONAL-COVERAGE\>-trycycler_medaka.fasta_ and _./results/assemblies/\<SAMPLE-ID\>~\<OPTIONAL-COVERAGE\>-trycycler_medaka_pilon.fasta_. Thereby, ~\<OPTIONAL-COVERAGE\> is only used for the file names if Trycycler was run on assemblies computed for a read sub sample. Both polishing tools were run with the (`T1_runMedaka.sh`)[./scripts/T1_runMedaka.sh] and (`T2_runPilon.sh`)[./scripts/T2_runPilon.sh] scripts.

### Evaluating the de novo assemblies
To evaluate the computed assemblies the tool Quast was used. The script (`S7_assemblyEvaluation.sh`)[./scripts/S7_assemblyEvaluation.sh] was used as a wrapper for Quast. For re-running it requires the definition of two parameters:
1. the reference set for which the assemblies are evaluated, i.e. all computed assemblies of one reference set are evaluated jointly. One of CFT073, MGH78578, RN4220~guppy3210 or RN4220~guppy4011
2. (optional) specifies a sub sample coverage, one of 200X, 150X, 100X, 80X, 60X, 40X, 20X, 15X, 10X, 8X, 6X, 4X, 2X, 1X

The script will basically collect all assemblies matching the specification from the input parameters from the _./results/assemblies/_ directory, the reference genome and the reference genome annotation and run Quast with those. In addition the `--gene-finding`, `--no-icarus` and `--no-plots` parameters are automatically set. The `--gene-finding` parameter specifies that GeneMarkS is run on the input assemblies for gene prediction, the `--no-icarus` parameter specifies that no [Icarus](http://quast.sourceforge.net/icarus.html) plots are generated and the `--no-plots` parameter specifies that the plots generated by Quast are not saved as extra files (they are still included in a _.html_ report). If one of the RN4220 reference sets are specified the additional `--fragmented` parameter is set, which results in Quast trying to detect misassembly events due to a fragmented reference genome (what is clearly the case for the RN4220 reference).

Quast allows to alter the features that are written to the _.txt_ _.tex_ and _.tsv_ reports. For this project information about the total number of mismatches, indels, indel lengths and misassemblies was added. The full feature list is specified in the (`S7_assemblyEvaluation.sh`)[./scripts/S7_assemblyEvaluation.sh] file. General information about adjusting the Quast reports can be found in the (Quast manual)[http://quast.sourceforge.net/docs/manual.html]. In order to run gene prediction with GeneMarkS with Quast (which is specified by setting the `--gene-finding` parameter) a license key for GeneMarkS has to be acquired (here)[http://topaz.gatech.edu/GeneMark/license_download.cgi] and placed in the users home directory. The Quast reports were stored as _./results/quast/\<SAMPLE-ID\>~\<OPTIONAL-COVERAGE\>_. The _.html_ and _.tsv_ reports output by Quast can be accessed via (supplementary file 4)[./supplementary_files/F4_quastReports.zip].

Further analysis of the Quast reports were conducted with the Python notebooks (`S8_assemblyQualityVisualization.ipynb`)[./scripts/S8_assemblyQualityVisualization.ipynb] and (`S9_subsampleAssemblyQualityVisualization`)[./scripts/S9_subsampleAssemblyQualityVisualization.ipynb].

### Discovery of single nucleotide variations
Finally, it was investigated how hybrid and long read only approaches can be used for the discovery of single nucleotide variations. Thereby it was decided to only investigate this for the CFT073 reference set. First, either the short- and long read sets alone or both in combination were used as input for bcftools. bcftools requires the reads to be mappped to a reference genome first. This was conducted again by using minimap2 for mapping and samtools for sorting and conversion into BAM format with the same parameters as described before. These BAM files are used as input for the `bcftools mpielup` and `bcftools call` commands. For the `mpileup` command the parameter `--skip-indels` was set to ignore insertions and deletions and for the `call` command the parameters `--ploidy 1 -m -v` were set to account for a haploid reference genome, use an updated variant caller algorithm and output only variant sites respectively. These steps were executed automatically with the script (`S10_SNVdiscovery.sh`)[./scripts/S10_SNVdiscovery.sh]. The resulting _.vcf_ files were stored at _./results/SNVdiscovery/CFT073/bcftools/_.

Second, the hybrid Unicycler and long read consensus Trycycler (polished with Medaka) assemblies were used as input for DNAdiff and Mauve. For DNAdiff the following commands were run
```
dnadiff -p trycycler_medaka-dnadiff ./data/references/CFT073.fasta ./results/assemblies/CFT073-trycycler_medaka.fasta
dnadiff -p unicycler~hybrid-dnadiff ./data/references/CFT073.fasta ./results/assemblies/CFT073-unicycler~hybrid.fasta
```
and the results were stored at the _./results/SNVdiscovery/CFT073/DNAdiff/_ directory. Mauve was run as GUI program: The two assemblies were aligned with the CFT073 reference genome using progressiveMauve and the SNVs were exported via the Tools \> Export \> Export SNPs menu option and stored at the _./results/SNVdiscovery/CFT073/mauve/_ directory. The results were compared and analyzed with the notebook (`S11_SNVDiscoveryAnalysis.ipynb`)[./scripts/S11_SNVDiscoveryAnalysis.ipynb].

***

If you have any questions about the content of this project, feel free to contact me at simon.hackl@student.uni-tuebingen.de