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
### General note on the file-structure
### Dependencies
### Download of reference genomes and annotations
### Download of public read data sets
### Quality control of the read data sets
### Downsampling of the read sets
### Conducting the de novo assemblies
### Evaluating the de novo assemblies
### Discovery of single nucleotide variations