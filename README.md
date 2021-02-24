# Researchproject:
### Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads
***
***
# Project description
The advent of third generation long read sequencing technologies, i.e. Oxford Nanopore Technologies (ONT) -ION devices, greatly expanded the possibilities for de novo assembly in terms of the contiguity of the reconstructed genomes. However, current long read technologies suffer from a not insignificant error rate compared to second generation short read sequencing technologies. Mainly due to new developments in base calling algorithms (the translation of signals from the sequencing device into DNA sequences) and further development of the protein pores used in ONT devices, this error rate is constantly decreasing. In addition, many new algorithmic approaches have been developed over the last years to efficiently use these long reads alone or in combination with short reads, i.e. in a hybrid approach, for de novo assembly or error correction of long reads.

In the course of the project _Elaboration of a best practice for hybrid and long read de novo assembly of bacterial genomes utilizing Illumina and Oxford Nanopore Technologies reads_ we asked to what extent a quality difference between long read only and hybrid approaches for de novo assembly is noticeable. To answer this, a selection of promising assemblers was chosen based on current benchmark reviews, a de novo assembly was performed with these and the resulting assemblies were subsequently evaluated according to common criteria.

Another aspect is the cost efficiency of using long reads. Increased read coverage can theoretically correct errors in the long reads, but this also increases the sequencing cost. We have investigated how different depths of coverage of the long reads affect the quality of the computed assemblies.

To also include the feasibility of the created de novo assemblies being used for downstream applications we studied the discovery of base substitutions by using the most qualitative de novo assemblies relative to the reference genomes and compared this to the discovery of such substitutions by using the input reads alone.

Our findings were summarized in a best practice recommendation.

This repository grants access to the projects report, supplementary data, the most important result files but especially an extended description of the methods conducted to enable the reproduction of the results.
***
# Best practice recommendation
***
***
# Extended description of methods
### General note on the structure of this repository
Most of the results for this project were obtained by running `bash` scripts and `Jupyter`/`python` notebooks which are supplied in this repository. Some steps were conducted manually and for these the explicit commands will be given below. Not all but the most important results are stored within this repository and can be accessed without re-running the scripts. The structure of the directories reflects the one which is created by running all scripts and should therefore not be changed in order to prevent errors during the execution of the provided scripts. As a naming convention, the prefix _S1_ to _S7_ is used for the main scripts which also give the order in which the scripts are supposed to be run. The _T1_ and _T2_ scripts are helper scripts used during the manual application of Trycycler and a _R_ prefix is used for `Jupyter` notebooks used for pre-processing the generated results file in the respective stages.

### Dependencies and software used
In the following a tabular overview of all software requirements to reproduce the results of this project are given. As some of the tools were accessible via the systems path variables on which the scripts were run but others were installed locally in a _./tools/_ directory, variables are defined in each script to access the executables of each tool; if necessary these variables can easily be changed in the respective scripts.

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
### Download of public read data sets
### Quality control of the read data sets
### Downsampling of the read sets
### Conducting the de novo assemblies
### Evaluating the de novo assemblies
### Discovery of single nucleotide variations 