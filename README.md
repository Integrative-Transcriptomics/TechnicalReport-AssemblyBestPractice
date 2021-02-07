# Researchproject: 
### Comparison of hybrid/longread assembly software and elaboration of a best practice
***
***
# Overview

With the advance of new sequencing technologies, the demands on and possibilities of bioinformatic software for data analysis and genome reconstruction are increasing. Accordingly, a large variety of tools for genome reconstruction has been established over the years.<sup>1</sup>
This research project tries to give an overview of the state of the art software for whole genome assembly of prokaryotes: Based on a literature review a promising pre-selection of assembly software is made and it is investigated which metrics are useful to assess assembly quality. In a subsequent analysis the following questions are dealt with in particular:

- Do some assemblers perfrom significantly better than others?
- What edge cases (regarding read data and the underlying genome) are problematic for which assemblers?
- Can longread only assemblies compete with hybrid assemblies?
- How low can read coverage become until a significant drop in assembly quality is observed? (This also aims to reduce costs for sequencing experiments.)

The final goal of the project is the elaboration of a best practice (possibly comprising a script) based on the findings of the assembly software evaluation. The evaluation covers different project phases which are discussed in more detail below. To keep track of the projects progress the different phases are marked with a checkmar once completed.

- [Literature research and acquisition of reference data](#literature-research-and-acquisition-of-reference-data) :heavy_check_mark:
- [Assembly of reference datasets](#assembly-of-reference-datasets) :heavy_check_mark:
- [Downsampling of reads to decrease coverage depth](#downsampling-of-reads-to-decrease-coverage-depth) :heavy_check_mark:
- [Detection of SNPs](#detection-of-snps) :heavy_check_mark:
- [Conclusion and best practice](#conclusion-and-best-practice) :heavy_check_mark:
- [References](#references)
&nbsp;
***
# Literature research and acquisition of reference data

Based on benchmarks and comparisons of longread as well as hybrid assembly software published from 2019 to 2020<sup>2-5</sup> promising candidates for further testing, based on the conclusions of the respective authors, were selected.

The same articles were examined with regard to the authors' evaluation methods and the accessibility of the reference data sets used.

In addition it was searched for recently developed assemblers, polishing software and evaluation methods which may not have been covered in the benchmarks.

A tabular overview of the software covered by the investigated literature is given below. Entries with a * were not covered in any publications but still considered for further steps.

| Long-read assembler                                                                                         | Hybrid assembler                                   | Assembly polisher<sup>i</sup>                                                                             | Assembly evaluation                           |
|-------------------------------------------------------------------------------------------------------------|----------------------------------------------------|-----------------------------------------------------------------------------------------------|-----------------------------------------------|
| [Unicycler](https://github.com/rrwick/Unicycler)                                                            | [Unicycler](https://github.com/rrwick/Unicycler)   | [Minipolish](https://github.com/rrwick/Minipolish) ([Racon](https://github.com/isovic/racon)) | [Quast](https://github.com/ablab/quast)       |
| [Canu](https://github.com/marbl/canu)                                                                       | [SPAdes-hybrid](https://github.com/ablab/spades)   | [NextPolish](https://github.com/Nextomics/NextPolish)                                         | [ALE](https://github.com/sc932/ALE)           |
| [Flye](https://github.com/fenderglass/Flye)                                                                 | [MaSuRCA](https://github.com/alekseyzimin/masurca) | [Pilon](https://github.com/broadinstitute/pilon)                                              | [MUMmer4](https://github.com/mummer4/mummer)  |
| [Miniasm](https://github.com/lh3/miniasm)/[Minipolish](https://github.com/rrwick/Minipolish)                                                        | [HASLR](https://github.com/vpc-ccg/haslr)*         | [Medaka](https://github.com/nanoporetech/medaka)*                                             | [REAPR](https://www.sanger.ac.uk/tool/reapr/) |
| [NECAT](https://github.com/xiaochuanle/NECAT)                                                               |                                                    |                                                                                               | [GMASS](https://github.com/jkimlab/GMASS)*    |
| [Raven](https://github.com/lbcb-sci/raven)                                                                  |                                                    |                                                                                               |                                               |
| [Redbean](https://github.com/ruanjue/wtdbg2)                                                                |                                                    |                                                                                               |                                               |
| [Shasta](https://github.com/chanzuckerberg/shasta)                                                          |                                                    |                                                                                               |                                               |
| [NextDenovo](https://github.com/Nextomics/NextDenovo)/[NextPolish](https://github.com/Nextomics/NextPolish) |                                                    |                                                                                               |                                               |

<sup>i</sup> <small>Assembly polisher software uses short- or long-reads to refine an assembly and to correct errors.</small> 

Based on the authors conclusions **Unicycler (hybrid and longread)**, **Canu**, **Flye** and **Raven** were selected as promising candidates for further testing. In addition the newly developed (January 2020) hybrid assembler **HASLR** was included, though it was not considered in the benchmarks and the longread consensus tool **Trycycler** was used.

Reference genomes, comprising [*Escherichia coli* strain CFT073](https://www.ncbi.nlm.nih.gov/nuccore/NC_004431.1) and [*Klebsiella pneumoniae* strain MGH78578](https://www.ncbi.nlm.nih.gov/nuccore/NC_009648.1) with 5 plasmids ([pKN3](https://www.ncbi.nlm.nih.gov/nuccore/NC_009649.1), [pKN4](https://www.ncbi.nlm.nih.gov/nuccore/NC_009650.1), [pKN5](https://www.ncbi.nlm.nih.gov/nuccore/NC_009651.1), [pKN7](https://www.ncbi.nlm.nih.gov/nuccore/NC_009652.1), [pKN8](https://www.ncbi.nlm.nih.gov/nuccore/NC_009653.1)), and the [associated raw sequencing data](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA422511) were collected from De Maio et al.<sup>5</sup>. In addition, short- and longreads of [*Staphylococcus aureus* strain RN4220](https://www.ncbi.nlm.nih.gov/nuccore/NZ_AFGU00000000.1) were provided by the group: The longreads are provided in two versions with different basecaller versions being used to call the nucleotides. The reference genome was taken from Nair et al.<sup>6</sup>; the RN4220 reference is fragmented into 179 contigs.

***
# Assembly of reference datasets

## Procedure

Before the assemblies were conducted, possible adapters of the ONT reads were removed using Porechop and the quality of all reads was validated using FastQC. For each of the reference genomes, short- and long-reads were mapped against the reference using minimap2, the depth of coverage was calculated using samtools depth and a depth of coverage plot was created to get insights into which genomic regions are actually supported by the input data.

Next, all of the above listed tools were run with each of the reference datasets. For *S. aureus* two assemblies were run, one with the old and one with the newer reads (regarding the basecaller version):
  - As all tools except canu, which performs an internal read correction, were described as being applicable to raw reads, no pre-assembly longread correction was performed.
  - As canu does not include an automated assembly polishing step Medaka was used to polish the canu assemblies.

An illustrated pipeline of **Trycycler** can be seen [here](https://github.com/rrwick/Trycycler/wiki/Illustrated-pipeline-overview). It comprises a more manual procedure in which different longread assemblies are clustered based on their distance calculated with [**Mash**](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x), sizes and indels.
  - For all trycycler runs as many as possible longread assemblies were considered.
  - In some cases contigs had to be removed due to low similarity to the other input contigs. This will be listed in the results section.

Furthermore, all assemblies were evaluated on a finer level with Quast: Quast provides the most important metrics for the evaluation of assemblies like GC content, length comparisons, indicators for miss assemblies like repetitive regions, local miss assemblies like indels, insertions and gaps as well as statistics for read coverage and structural variations.

## Results

### FastQC readquality

The FastQC reports revealed mean short read qualities of above 28 for all samples and all read positions and long read qualities ranging from 12 to 18 for the E. coli and K. pneumoniae reads, 14 to 24 for the old basecaller version S. aureus reads and 14 to 32 for the new basecaller version S. aureus reads.

### Depth of coverage across references

The following plots depict the average (across a specific resolution, e.g. 10 kBp intervals) depth of coverage of the E. coli and K. pneumoniae references. The left side shows longread coverage while the right side shows short read coverage.

##### Depth of coverage, reference E. coli CFT073
<img align="center" src="images/depthOfCoverage_CFT073_longread.png" width="400" > ![](depthOfCoverage_CFT073_longread.png)
<img align="center" src="images/depthOfCoverage_CFT073_shortread.png" width="400" > ![](depthOfCoverage_CFT073_shortread.png)

##### Depth of coverage, reference K. pneumoniae MGH78578, genome
<img align="center" src="images/depthOfCoverage_MGH78578-chr_longread.png" width="400" > ![](depthOfCoverage_MGH78578-chr_longread.png)
<img align="center" src="images/depthOfCoverage_MGH78578-chr_shortread.png" width="400" > ![](depthOfCoverage_MGH78578-chr_shortread.png)

##### Depth of coverage, reference K. pneumoniae MGH78578, plasmid one
<img align="center" src="images/depthOfCoverage_MGH78578-p1_longread.png" width="400" > ![](depthOfCoverage_MGH78578-p1_longread.png)
<img align="center" src="images/depthOfCoverage_MGH78578-p1_shortread.png" width="400" > ![](depthOfCoverage_MGH78578-p1_shortread.png)

##### Depth of coverage, reference K. pneumoniae MGH78578, plasmid two
<img align="center" src="images/depthOfCoverage_MGH78578-p2_longread.png" width="400" > ![](depthOfCoverage_MGH78578-p2_longread.png)
<img align="center" src="images/depthOfCoverage_MGH78578-p2_shortread.png" width="400" > ![](depthOfCoverage_MGH78578-p2_shortread.png)

##### Depth of coverage, reference K. pneumoniae MGH78578, plasmid three
<img align="center" src="images/depthOfCoverage_MGH78578-p3_longread.png" width="400" > ![](depthOfCoverage_MGH78578-p3_longread.png)
<img align="center" src="images/depthOfCoverage_MGH78578-p3_shortread.png" width="400" > ![](depthOfCoverage_MGH78578-p3_shortread.png)

##### Depth of coverage, reference K. pneumoniae MGH78578, plasmid four
<img align="center" src="images/depthOfCoverage_MGH78578-p4_longread.png" width="400" > ![](depthOfCoverage_MGH78578-p4_longread.png)
<img align="center" src="images/depthOfCoverage_MGH78578-p4_shortread.png" width="400" > ![](depthOfCoverage_MGH78578-p4_shortread.png)

##### Depth of coverage, reference K. pneumoniae MGH78578, plasmid five
<img align="center" src="images/depthOfCoverage_MGH78578-p5_longread.png" width="400" > ![](depthOfCoverage_MGH78578-p5_longread.png)
<img align="center" src="images/depthOfCoverage_MGH78578-p5_shortread.png" width="400" > ![](depthOfCoverage_MGH78578-p5_shortread.png)

As the S. aureus reference genome was highly fragmented it was decided to depict average depth of coverage and length per contig. The contigs were sorted by their length. We see a low correlation of length and coverage.

##### Depth of coverage, reference S. aureus RN4220
<img align="center" src="images/depthOfCoverage_RN4220_shortread.png" width="400" > ![](depthOfCoverage_RN4220_shortread.png)
<img align="center" src="images/depthOfCoverage_RN4220_shortread.png" width="400" > ![](depthOfCoverage_RN4220_longread_old.png)

The following table depicts the percentages (relative to the reference genome length) of zero coverage depth and the mean coverage across the whole reference genome of the different read sets.

| read set                           | genome fraction with 0X coverage | mean coverage across genome |
|------------------------------------|:--------------------------------:|:---------------------------:|
| E. coli CFT073, shortreads         |                          0.046 % |                     93.09 X |
| E. coli CFT073, longreads          |                          0.000 % |                    244.72 X |
| K. pneumoniae MGH78578, shortreads |                         10.924 % |                     69.35 X |
| K. pneumoniae MGH78578, longreads  |                         11.238 % |                    239.32 X |
| S. aureus RN4220, shortreads       |                          0.006 % |                    499.03 X |
| S. aureus RN4220, longreads (old)  |                          0.008 % |                    528.24 X |
| S. aureus RN4220, longreads (new)  |                          0.005 % |                    534.18 X |

### Completion of assemblies

  - Except for HASLR on the E. coli CFT073 dataset all assemblies were able to complete without errors.
  - HASLR did only output a very short assembly of the K. pneumoniae MGH78578 genome.
  - Raven did only output a very short assembly of the S. aureus RN4220 genome for both longread versions.

### Trycycler runs

For all Trycycler runs all of the respective longread assemblies were used as input. In some cases whole assemblies or some contigs had to be deleted due to dissimilarities:
  - The reconcile step showed a circularization issue for the Canu assemblies that were polished with Medaka of the E. coli and K. pneumoniae references. Thus for the two references the original Canu assemblies were used.
  - The original Canu assemblies still showed circularization issues for the mentioned references, which, however, were resolved by increasing the maximal trimming length for circularization to 100,000 bp.
  - One contig of the Unicycler assembly of the K. pneumoniae reference was clustered alone, thus it showed no similarity with the contigs of the other assemblers and was removed.
  - One contig of the Canu assembly of the K. pneumoniae reference showed nearly the double length of the contigs it was clustered with, thus it was removed.
  - The contigs of the Unicycler and Raven assemblies of the E. coli reference exhibit an indel of size close to 800. This issue was resolved by increasing the maximal allowed indel length to the respective value.
  - The very short contig that was produced by the Raven assembler for the S. aureus reference was clustered with one contig of the Unicycler assembly. Thus it was retained.

### Quast results

The full Quast .html reports are deposited on this git. Here, only some of the reported features, considered worth discussing, are presented. All quast commands were run with the *--no-icarus* option and the *--fragmented* option was used for the S. aureus assemblies.

The following table depicts the genomic fraction of the assemblies in percent and the number of different missassemblie events considered by Quast:
  - *relocations* are gaps or overlaps of at least 1kbp of the alignment of the assembly to the reference.
  - *translocations* and *inversions* are alignments on different chromosomes and strands, respectively.

|                          | E. coli CFT073 | K. pneumoniae MGH78578 | S. aureus RN4220 |
|--------------------------|----------------|------------------------|------------------|
| Genome fraction (%)      | > 99.9 %       | > 86.5 %               | > 99.5 %         |
| Number of relocations    | 3-4            | 41 - 53                | 0                |
| Number of translocations | 0              | 0                      | 9-14 (10-13)     |
| Number of inversions     | 0              | 2                      | 0                |

While the number of relocations of the K. pneumoniae and the number of translocations of the S. aureus assemblies can be easily explaiend by the structure of the references (regions of no coverage and multiple contigs), the number of relocations in E. coli assemblies and inversions in K. pneumoniae may be real differences from the reference genome to the genome of the sample that was actually sequenced.

In addition the number of contigs differed much for the hybrid versus longread approach on the K. pneumoniae reference. For all other references all assemblers agreed on 1 to 2 contigs.
  - In K. pneumoniae the Unicycler (hybrid) approach produced 16 contigs versus 2 to 3 contigs for all other assemblers; This may be due to the fact that some regions in reference plasmids have a low short read coverage but no long read coverage. As short reads were only used in the Unicycler (hybrid) approach to construct the architecture of the genome, the plasmids may have been sequenced as a set of very short fragments. This is also supported by the fact that the approach produced 4 contigs larger than 50,000 bp which is more close to the number of contigs produced by the other approaches.

The most significant differences of the assemblies were detected in the number of indels, mismatches and the number of genes that were predicted by GeneMarkS as well as the number of annotated genes that are being completely recovered from the assemblies. These features are depicted below: 

<img align="center" src="images/quast-CFT073.png" width="400" > ![](quast-CFT073.png)
<img align="center" src="images/quast-MGH78578.png" width="400" > ![](quast-MGH78578.png)

<img align="center" src="images/quast-RN4220old.png" width="400" > ![](quast-RN4220old.png)
<img align="center" src="images/quast-RN4220new.png" width="400" > ![](quast-RN4220new.png)

  - From the pictures it is clearly visible that the number of predicted genes exceeds the number of annotated genes for the longread assemblies, except the S. aureus assembly with higher quality longreads. On the example of the CFT073 Canu and Unicycler (hybrid) assemblies it was investigated how many of the predicted genes are located inside an annotated/*true* gene and how the lengths of the annotated/*true* genes are distributed.

<img align="center" src="images/geneOverlaps-Canu-CFT073.png" width="400" > ![](geneOverlaps-Canu-CFT073.png)
<img align="center" src="images/geneLengths-Canu-CFT073.png" width="400" > ![](geneLengths-Canu-CFT073.png)

<img align="center" src="images/geneOverlapgs-UnicHybrid-CFT073.png" width="400" > ![](geneOverlapgs-UnicHybrid-CFT073.png)
<img align="center" src="images/geneLengths-UnicHybrid-CFT073.png" width="400" > ![](geneLengths-UnicHybrid-CFT073.png)

The investigation reveals that indeed the predicted genes are much shorter than the annotated genes and multiple predicted genes are located inside annotated genes for the longread E. coli assembly. It can be suggested, that the low longread quality, which may have also led to the high number of mismatches, disturb the GeneMarkS prediction. However, here it is questionable why Quast reports the annotated genes as being completely recovered.

## Conclusion
  - All longread assemblers (except the one that did not produce correct output) managed to reconstruct the major parts of the reference genomes.
  - Comparing the longread assemblers to each other we see the main differences in the number of indels and mismatches compared to the respective reference. 
Above all we do not see that one specific assembler outperforms all other assemblers in all features or at least in one specific feature across all reference species.
  - The hybrid assembler (Unicycler) outperforms the longread assemblers regarding the number of indels and mismatches.
  - The fact that the Trycycler assemblies that were polished with short reads using Pilon perform close the the Unicycler (hybrid) assembler indicates that the performance difference is due to the quality differences of the reads or the fact that more information may be contained in the shortreads
  - The latter may be rejected as we have the same breadth and high enough depth of coverage for both read types and the S. aureus longreads with higher quality resulted in better results for the longread assemblers. Here, especially the number of indels was greatly decreased.
  
***
# Downsampling of reads to decrease coverage depth

## Procedure

Next the effect of the coverage depth of longreads on the assembly quality was investigated. The tool [**Rasusa**](https://github.com/mbhall88/rasusa) was used to generate subsamples of the longreads with intended 200X, 150X, 100X, 80X, 60X, 40X, 20X, 15X, 10X, 8X, 6X, 4X, 2X and 1X coverage depth. Firtstly, the E. coli CFT073 reference set was investigated, due to its completenes regarding the given reference genome:
- Similar as the procedure applied to the reference read sets, each downsampled read set was mapped against the E. coli CFT073 reference genome using minimap2. Based on this the per base depth of coverage was measured using samtools depth. From this data, the mean depth of coverage (the sum of per base coverage depths divided by the genome length) was calculated. In addition the genomic fraction with zero coverage was determined.
- For each of the subsamples all assemblers (except HASLR) were run. For each coverage depth a Quast report comprising all assemblies was generated.

Secondly, the downsampling was conducted for the S. aureus RN4220 reference set to investigate the effect of the improved basecaller (and thus read quality) on the downsampling.
- Here, only the hybrid approach of Unicycler was used for assembly as it yielded the best results in the previous steps.

TODO: Check actual coverage of S. aureus downsamples and update x-labels on plots.
 
## Results

### Depth of coverage control - E. coli subsamples

| E. coli CFT073 downsample intended coverage | genome fraction with 0X coverage | mean coverage across genome |
|---------------------------------------------|----------------------------------|-----------------------------|
|                                        200X |                          0.000 % |                    187.95 X |
|                                        150X |                          0.000 % |                    140.93 X |
|                                        100X |                          0.000 % |                     93.92 X |
|                                         80X |                          0.000 % |                     75.18 X |
|                                         60X |                          0.000 % |                     56.33 X |
|                                         40X |                          0.001 % |                     37.58 X |
|                                         20X |                          0.001 % |                     18.82 X |
|                                         15X |                          0.002 % |                     13.54 X |
|                                         10X |                          0.010 % |                      9.07 X |
|                                          8X |                          0.233 % |                      7.24 X |
|                                          6X |                          0.672 % |                      5.44 X |
|                                          4X |                          1.844 % |                      3.77 X |
|                                          2X |                         14.345 % |                      1.88 X |
|                                          1X |                         38.647 % |                      0.94 X |

### Quast results - E. coli subsamples

<img align="center" src="images/quast-subsample-CFT073.png" width="800" > ![](quast-subsample-CFT073.png)

### Depth of coverage control - S. aureus subsamples

### Quast results - S. aureus subsamples

<img align="center" src="images/quast-subsample-RN4220.png" width="800" > ![](quast-subsample-RN4220.png)

## Conclusion
  - Down to a coverage of 40 X the assemblers could finish the assemblies correctly and the assembly statistics, regarding the quast report, do not deviate significantly from the full data set.
  - The unicycler hybrid approach was (due to using short reads) able to finish the assembly with all downsamples, but interestingly the number of contigs increased for lower coverages. This implies an important role of longreads for the contiguity of the assemblies.
  - For samples below 40X, Trycycler was unable to complete the reconcile step. In general, the lower the coverage, the less assemblers were able to complete the assembly.
  - However, the significant drops in performance for samples with < 6 X coverage may be mainly due to large fractions of the genome being not covered by at least one read. On the other hand one also observes perfomrance drops for samples with low coverage, but < 1% of the genome being not covered by at leas one read.
  - Finally, the results show, that a coverage of 40 X to 80 X may be sufficient to yield optimal results.
  - Regarding the improved basecaller quality of the S. aureus downsamples it is evident, that the higher quality reads enabled a more contigous assembly even at low coverages below 6X due to reduced number of contigs and a higher genome fraction and discovered genomic features. However, the number of mismatches increased slightly. This may presumably be due to higher fraction of the genome being completely assembled and thus being available for evaluating the number of mismatches at all. 
  
***
# Detection of SNPs

It was assesed in which extend the use of longreads, shortreads or both influences the detection of SNPs relative to the used reference genomes. Thereby, the main question was if using longreads only would result in the same SNP calls as using shortreads only. To validate the derived SNP sets they were compared to those derived from comparing the reference genome with the most qualitative assembly computed using the Unicycler hybrid approach using mauve and MUMmer.  

## Procedure

It was decided to conduct this study only with the E. coli CFT073 reference data set as for this set the most contigous reference genome and an appropriate read coverage was given and furthermore only investigate base changes and to exclude indels (mainly due to simplicity reasons). Firstly, short- and longreads were mapped to the reference genome using minimap2. The resulting .sam files were sorted and converted to .bam files. These were then used as input for bcftools mpielup command (using -O b, --skip-indels and --max-depth 1000 parameters in order to skip indel calling, produce .bcf formatted output and use all available reads per position, respecitvely) in three modes:
  - using the mapping of the longreads only 
  - using the mapping of the shortreads only
  - using the mapping of both long- and shortreads

The obtained .bcf files were further used as input for the bcftools call command (using --ploidy 1 parameter to account for the haploidity of bacteria) resulting in a .vcf file. 

Next, SNP sets were extracted based on the given QUAL column entry, representing the -10 log<sub>10</sub>( probability of no variant ), being any value, higher or equal than 100 or higher or equal than 200. Thus, the two last filterings yield SNPs with error probabilities of 10<sup>-10</sup> and 10<sup>-20</sup> according to bcftools.

Secondly, additional SNP sets were extraced by (i) computing a whole genome alignment of the reference genome and the Unicycler hybrid assembly and exporting SNPs using mauve (http://darlinglab.org/mauve/mauve.html) and (ii) using the reference genome and the Unicycler hybrid assembly as input for the MUMmer wrapper DNAdiff. The second approach was also used by De Maio et al.<sup>5</sup>, from which the reference E. coli CFT073 read set was taken, to detect SNPs relative to the assembly. Thereby, they detected a total number of 141 SNPs (whereby it is not clear which exact pipeline and parameters were used for the Unicycler assembly and which exact positions in the reference genome were considered as SNPs).

All derived SNP sets were comparatively visualized using Venny (https://bioinfogp.cnb.csic.es/tools/venny/).

## Results

##### Comparison of SNP sets detected by mapping different read sets to the reference genome, SNPs with any quality (QUAL) used
<img align="center" src="images/venn-reads-anyQUAL.png" width="500" > ![](venn-reads-anyQUAL.png)

##### Comparison of SNP sets detected using different read sets to the reference genome, SNPs with quality (QUAL) higher or equal to 100 used
<img align="center" src="images/venn-reads-100QUAL.png" width="500" > ![](venn-reads-100QUAL.png)

##### Comparison of SNP sets detected using different read sets to the reference genome, SNPs with quality (QUAL) higher or equal to 200 used
<img align="center" src="images/venn-reads-200QUAL.png" width="500" > ![](venn-reads-200QUAL.png)

##### Comparison of the SNPs detected using any read set (SNPs with any quality (QUAL)) and the SNPs derived by computing a WGA using mauve or mapping by MUMmer, using the assembly of the Unicycler hybrid approach, respectively
<img align="center" src="images/venn-reads-assembly-anyQUAL.png" width="500" > ![](venn-reads-assembly-anyQUAL.png)

##### Comparison of the SNPs detected using any read set (SNPs with quality (QUAL) higher or equal to 200 used) and the SNPs derived by computing a WGA using mauve or mapping by MUMmer, using the assembly of the Unicycler hybrid approach, respectively
<img align="center" src="images/venn-reads-assembly-200QUAL.png" width="500" > ![](venn-reads-assembly-200QUAL.png)


## Conclusion

TODO

***
# Summary and best practice

TODO
***
## References
<sup>1</sup> [Alice Maria Giani, et al. *Long walk to genomics: History and current approaches to genome sequencing and assembly.* Computational and Structural Biotechnology Journal, Volume 18, Pages 9-19, 2020.](https://www.sciencedirect.com/science/article/pii/S2001037019303277)

<sup>2</sup> [Goldstein, S., Beka, L., Graf, J. et al. *Evaluation of strategies for the assembly of diverse bacterial genomes using MinION long-read sequencing.* BMC Genomics, Volume 20, Number 23, 2019.](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5381-7#citeas)

<sup>3</sup> [Chen, Z., Erickson, D.L. & Meng, J. *Benchmarking hybrid assembly approaches for genomic analyses of bacterial pathogens using Illumina and Oxford Nanopore sequencing.* BMC Genomics, Volume 21, Number 631, 2020.](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07041-8)

<sup>4</sup> [Wick RR and Holt KE. *Benchmarking of long-read assemblers for prokaryote whole genome sequencing [version 3; peer review: 4 approved].* F1000Research 2020, 8:2138](https://f1000research.com/articles/8-2138/v3)

<sup>5</sup> [De Maio N, Shaw LP, Hubbard A, et al. *Comparison of long-read sequencing technologies in the hybrid assembly of complex bacterial genomes.* Microb Genom. 5(9), 2019.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6807382/)

<sup>6</sup>[Dhanalakshmi Nair, Guido Memmi, et al. *Whole-Genome Sequencing of Staphylococcus aureus Strain RN4220, a Key Laboratory Strain Used in Virulence Research, Identifies Mutations That Affect Not Only Virulence Factors but Also the Fitness of the Strain* Journal of Bacteriology  193 (9), Apr 2011](https://jb.asm.org/content/193/9/2332)
