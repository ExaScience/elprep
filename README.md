# Overview

elPrep is a high-performance tool for analyzing .sam/.bam files (up to and including variant calling) in sequencing pipelines. 
The key advantage of elPrep is that it only performs a single-pass to process a .sam/.bam 
file, independent of the number of processing steps that need to be applied in a particular pipeline, greatly improving 
runtime performance.

elPrep is designed as an in-memory and multi-threaded application to fully take advantage of the processing power available 
with modern servers. Its software architecture is based on functional programming techniques, which allows easy composition 
of multiple alignment filters and optimizations such as loop merging. To make this possible, elPrep proposes a couple of new algorithms. 
For example, for duplicate marking we devised an algorithm that expresses the operation as a single-pass filter using memoization 
techniques and hierarchical hash tables for efficient parallel synchronisation. For base quality score recalibration (BQSR) 
we designed a parallel range-reduce algorithm. For variant calling, we designed a parallel algorithm that overlaps execution 
as much as possible with other phases in the pipeline.

Our benchmarks show that elPrep executes a 5-step variant calling best practices pipeline (sorting, duplicate marking, base 
quality score recalibration and application, and variant calling) between 6-10 times faster than other tools for whole-exome 
data, and 8-20x faster for whole-genome data.

The main advantage of elPrep is very fast execution times on high-end servers, as is available for example 
through cloud computing or custom server setups. We do not recommend using elPrep on laptops, desktops, 
or low-end servers. Please consult the system requirements below for more details. 

elPrep is being developed at the [ExaScience Life Lab](http://www.exascience.com) at Imec. For questions, use our 
mailing list (below), our github page, or contact us via exascience@imec.be.


Fig. 1 Improvements with elprep 5 wrt runtime, RAM, and disk use for a variant calling best practices pipeline on a 50x 
Platinum NA12878 WGS aligned against hg38. elPrep combines the execution of the 5 pipeline steps for efficient 
parallel execution.

![NA12878 Platinum Genome run](https://github.com/ExaScience/elprep/raw/master/images/wgs-benchmarks-elprep5.jpg)

For more benchmark details, please consult our publication list.

# Advantages

The advantages of elPrep include:

* efficient multi-threaded execution
* operates mainly in-memory, few intermediate files are generated
* 100% equivalent output to results produced by other widely used tools
* compatible with existing tools
* modular, easy to add and remove pipeline steps

# Availability

elPrep is released and distributed under a dual-licensing scheme. It is released as an open-source project under the terms of the GNU Affero General Public License 
version 3 as published by the Free Software Foundation, with Additional Terms. Please see the file LICENSE.txt for a 
copy of the GNU Affero Public License and the Additional Terms. For inquiries about the premium licensing option contact
 us via exascience@imec.be.

We also provide a download of a precompiled binary.

## Binaries

elPrep 5 binaries can be downloaded from this website:

* https://www.imec-int.com/en/expertise/lifesciences/genomics/dna-sequence-analysis-software

## GitHub

The elPrep source code is freely available on GitHub. elPrep is implemented in Go and tested for Linux.

elPrep GitHub URL:

* https://github.com/exascience/elprep

## Dependencies

elPrep works with the .sam, .bam, and .vcf formats as input/output. Previously, there was a dependency on samtools to 
read and write .bam files, but since elPrep4.0, .bam files are directly supported by elPrep, with no need for samtools 
to be present anymore. If you need support for .cram files, consider converting them to/from .bam files before/after 
elPrep using samtools, or other alternatives. There was previously also a dependency on bcftools to read and 
write .vcf.gz and .bcf files, but since elPrep5.0, .vcf.gz are directly supported by elPrep as well. If you need support
 for .bcf files, consider converting them to/from .vcf.gz files before/after elPrep using bcftools, or other alternatives.

elPrep relies on its own .elsites file format for representing known sites during base quality score recalibration. 
Such .elsites files can be generated from .vcf files using the elPrep vcf-to-elsites command, and from .bed files using 
bed-to-elsites. elPrep also uses its ows .elfasta file format for representing references during base quality score 
recalibration and variant calling. They can be generated from .fasta files using the elPrep fasta-to-elfasta command.

There are no dependencies on other tools.

## Building

The following is only relevant information if you wish to build elPrep yourself. It is not necessary to use the 
elPrep binary we provide.

elPrep (since version 3.0) is implemented in Go. Please make sure you have a working installation of 
Go. You can either install Go from the [Go website](https://golang.org). Alternatively, most package managers provide 
options to install Go as a development tool. Check the documentation of the package manager of your Linux distribution 
for details.

First checkout the elPrep sources using the following command:

		go get -u github.com/exascience/elprep

This downloads the elPrep Go source code, and creates the elprep binary in your configured Go home folder, for example
 ~/go/bin/elprep. See the [GOPATH](https://golang.org/cmd/go/#hdr-GOPATH_environment_variable) variable for your Go home
  folder.

Add the binary to your path, for example:
  
          export PATH=$PATH:~/go/bin

## Compatibility

elPrep has been developed for Linux and has not been tested for other operating systems. We have tested elPrep with the 
following Linux distributions:

* Ubuntu 14.04.5 LTS
* Manjaro Linux
* Red Hat Enterprise Linux 6.4 and 6.5
* Amazon Linux 2

# Memory Requirements

## RAM

elPrep is designed to operate in memory, i.e. data is stored in RAM during computation. As long as you do not use the 
in-memory sort, mark duplicates filter, base recalibration, or variant calling, elPrep operates as a streaming tool, and
 peak memory use is limited to a few GB.

elPrep also provides a tool for splitting .sam files per chromosome - or better per groups of chromosomes - and 
guarantees that processing these split files and then merging the results does not lose information when compared to 
processing a .sam file as a whole. Using the split/merge tool greatly reduces the RAM required to process a .sam file, 
but it comes at the cost of an additional processing step.

We recommend the following minimum of RAM when executing memory-intensive operations such as sorting, duplicate marking,
 base quality score recalibration and haplotype caller:

* whole-genome 30x: 128 GB RAM using the elprep split/filter/merge mode (sfm)
* whole-genome 50x: 200 GB RAM using the elprep split/filter/merge mode (sfm)
* whole-exome  30x: 24GB RAM using the elprep split/filter/merge mode (sfm)
* whole-exome  30x: 92GB RAM using the elprep in memory mode (filter)

These numbers are only estimates, and the actual RAM use may look different for your data sets.

## Disk Space

elPrep by default does not write any intermediate files, and therefore does not require additional (peak) disk space 
beyond what is needed for storing the input and output files. If you use the elPrep split and merge tools, elPrep 
requires additional disk space equal to the size of the input file.

# Mailing List and Contact

Use the Google [forum](https://groups.google.com/d/forum/elprep) for discussions. You need a Google account to subscribe
 through the forum URL. You can also subscribe without a Google account by sending an email to 
 elprep+subscribe@googlegroups.com.

You can also contact us via exascience@imec.be directly.

For inquiries about commercial licensing options contact us via exascience@imec.be. 

# Citing elPrep

Please cite the following articles:

Herzeel C, Costanza P, Decap D, Fostier J, Verachtert W (2019) elPrep 4: A multithreaded framework for sequence analysis. PLoS ONE 14(2): e0209523. https://doi.org/10.1371/journal.pone.0209523

Herzeel C, Costanza P, Decap D, Fostier J, Reumers J (2015) elPrep: High-Performance Preparation of Sequence Alignment/Map Files for Variant Calling. PLoS ONE 10(7): e0132868. https://doi.org/10.1371/journal.pone.0132868

Costanza P, Herzeel C, Verachter W (2019) A comparison of three programming languages for a full-fledged next-generation sequencing tool. BMC Bioinformatics 2019 20:301. https://doi.org/10.1186/s12859-019-2903-5

If performance is below your expectations, please contact us first before reporting your results.

# Examples

## Variant calling pipeline for whole-genome data (WGS)

The following elprep command shows a 5-step variant calling best practices pipeline on WGS data: 

    elprep sfm NA12878.input.bam NA12878.output.bam
               --mark-duplicates --mark-optical-duplicates NA12878.output.metrics
               --sorting-order coordinate
               --bqsr NA12878.output.recal --known-sites dbsnp_138.hg38.elsites,Mills_and_1000G_gold_standard.indels.hg38.elsites
               --bqsr-reference hg38.elfasta
               --haplotypecaller NA128787.output.vcf.gz

The command executes a pipeline that consists of 5 steps: sorting, PCR and optical duplicate marking, base quality score 
recalibration and application, and variant calling.

We can break up the command as follows:

 - The sfm subcommand tells elprep to run in sfm (split/filter/merge) mode. This is generally the preferred mode for WGS
  data, unless the data has very low coverage (<= 10x).
    
 - The input file is "NA12878.input.bam".
    
 - Output is written to a file "NA12878.bam" that contains the result of modifying the input bam file by performing 
 duplicate marking, sorting, and base quality score recalibration and application.
    
 - The flags --mark-duplicates and --mark-optical-duplicates instruct elprep to perform PCR and optical duplicate 
 marking respectively. The statistics generated by this are written to a file "NA12878.output.metrics".
    
 - The flag --sorting-order tells elprep to sort the input bam file by coordinate order.
    
 - The flag --bqsr instructs elprep to perform base quality score recalibration. The statistics generated by this are 
 written to a file "NA12878.output.recal". The --bqsr flags also need to know the reference fasta file with which the 
 input bam was created, cf. the "--bqsr-reference hg38.elfasta" option. Note the file extension ".elfasta". elPrep requires 
 converting the fasta file to this format before running the pipeline via the command "elprep fasta-to-elfasta 
 hg38.fasta hg38.elfasta". The --bqsr option also needs to know the known variant sites, passed via the "--known-sites 
 dbsnp_138.hg38.elsites" option. Note the file extension ".elsites". elPrep requires converting vcf files to this format
  before running the pipeline via the command "elprep vcf-to-elsites dbsnp_138.hg38.vcf dbsnp_138.hg38.elsites".
    
 - The flag --haplotypecaller instructs elprep to perform variant calling. It uses the same reference fasta as the one 
 passed for --bqsr (via --bqsr-reference). The result of this step is written are written to a file "NA12878.output.vcf.gz".

For details, consult the manual reference pages.

## Variant calling pipeline for whole-exome data (WES)

The following elprep command shows a 5-step variant calling best practices pipeline on WES data:

    elprep sfm NA12878.input.bam  NA12878.output.bam 
               --mark-duplicates --mark-optical-duplicates NA12878.output.metrics 
               --sorting-order coordinate 
               --bqsr NA12878.output.recal --known-sites dbsnp_137.hg19.elsites,Mills_and_1000G_gold_standard.indels.hg19.elsites 
               --bqsr-reference hg19.elfasta 
               --haplotypecaller NA12878.output.vcf.gz 
               --target-regions nexterarapidcapture_expandedexome_targetedregions.bed 

elPrep uses an internal ".elfasta" format for representing fasta files, which can be created using the 
"elprep fasta-to-eflasta" command before running the pipeline. Similarly, elPrep uses an internal format for representing 
vcf files containing known variant sites (.elsites), which can be created using the command "elprep vcf-to-elsites".

For details, consult the manual reference pages.

# Manual Reference Pages

## Name

### elprep filter - a commandline tool for filtering and updating .sam/.bam files and variant calling

## Synopsis

	elprep filter input.sam output.sam --mark-duplicates --mark-optical-duplicates output.metrics
	                                   --sorting-order coordinate 
	                                   --bqsr output.recal --bqsr-reference hg38.elfasta --known-sites dbsnp_138.hg38.elsites 
	                                   --haplotypecaller output.vcf.gz

	elprep filter input.bam output.bam --mark-duplicates --mark-optical-duplicates output.metrics 
	                                   --sorting-order coordinate 
	                                   --bqsr output.recal --bqsr-reference hg38.elfasta --known-sites dbsnp_138.hg38.elsites 
	                                   --haplotypecaller output.vcf.gz

	elprep filter /dev/stdin /dev/stdout --mark-duplicates --mark-optical-duplicates output.metrics 
	                                     --sorting-order coordinate 
	                                     --bqsr output.recal --bqsr-reference hg38.elfasta --known-sites dbsnp_138.hg38.elsites	
	                                     --haplotypecaller output.vcf.gz

## Description

The elprep filter command requires two arguments: the input file and the output file. The input/output format can be 
.sam or .bam. elPrep determines the format of the input by analyzing the actual contents of the input file. The format 
of the output file is determined by looking at the file extension. elPrep also allows to use /dev/stdin and /dev/stdout 
as respective input or output sources for using Unix pipes. When doing so, elPrep assumes output is in .sam format, 
which can be changed by additional parameters (see below).

The elprep filter command-line tool has three types of command options: filters, which implement actual .sam/.bam 
manipulations, sorting options, and execution-related options, for example for setting the number of threads. For 
optimal performance, issue a single elprep filter call that combines all filters you wish to apply.

The order in which command options are passed is ignored. For optimal performance, elPrep always applies 
filters/operations in the following order:

1. filter-unmapped-reads or filter-unmapped-reads-strict
2. filter-mapping-quality
3. filter-non-exact-mapping-reads or filter-non-exact-mapping-reads-strict
4. filter-non-overlapping-reads
5. clean-sam
6. replace-reference-sequences
7. replace-read-group
8. mark-duplicates
9. mark-optical-duplicates
10. bqsr
11. remove-duplicates
12. remove-optional-fields
13. keep-optional-fields
14. haplotypecaller

Sorting is done after filtering.

Please also see the **elprep sfm** command.

### Unix pipes

elPrep is compatible with Unix pipes and allows using /dev/stdin and /dev/stdout as input or output sources. elPrep 
analyzes the input from /dev/stdin to determine if it is in .sam or .bam format, and assumes that output to /dev/stdout 
is in .sam format, unless specified otherwise (see below).

## Filter Command Options

### --replace-reference-sequences file

This filter is used for replacing the header of a .sam/.bam file by a new header. The new header is passed as a single 
argument following the command option. The format of the new header can either be a .dict file or another .sam/.bam file from which you wish to extract the new header.

All alignments in the input file that do not map to a chromosome that is present in the new header are removed. 
Therefore, there should be some overlap between the old and the new header for this command option to be meaningful. 
The option is typically used to reorder the reference sequence dictionary in the header.

Replacing the header of a .sam/.bam file may destroy the sorting order of the file. In this case, the sorting order in 
the header is set to "unknown" by elPrep in the output file (cf. the 'so' tag).

### --filter-unmapped-reads

Removes all alignments in the input file that are unmapped. An alignment is determined unmapped when bit 0x4 of its FLAG
 is set, conforming to the SAM specification.

### --filter-unmapped-reads-strict

Removes all alignments in the input file that are unmapped. An alignment is determined unmapped when bit 0x4 of its FLAG
 is set, conforming to the SAM specification. Also removes alignments where the mapping position (POS) is 0 or where the 
 reference sequence name (RNAME) is *. Such alignments are considered unmapped by the SAM specification, but some 
 alignment programs may not mark the FLAG of those alignments as unmapped.

### --filter-mapping-quality mapping-quality

Remove all alignments with mapping quality lower than the given mapping quality.

### --filter-non-exact-mapping-reads

Removes all alignments where the mapping is not an exact match with the reference, albeit soft-clipping is allowed. This
 filter checks the CIGAR string and only allows occurences of M and S.

### --filter-non-exact-mapping-reads-strict

Removes all alignments where the mapping is not an exact match with reference or not a unique match. This filter checks 
for each read that the following optional fields are present with the following values: X0=1 (unique mapping), X1=0 
(no suboptimal hit), XM=0 (no mismatch), XO=0 (no gap opening), XG=0 (no gap extension).

### --filter-non-overlapping-reads bed-file

Removes all reads where the mapping positions do not overlap with any region specified in the bed file. Specifically, 
either the start or end of the read's mapping position must be contained in an interval, or the read is removed from the
 output.
 
 This option produces a different result from --target-regions option. For the difference between both options and 
 details on the algorithms, please consult our latest publication.

### --replace-read-group read-group-string

This filter replaces or adds read groups to the alignments in the input file. This command option takes a single 
argument, a string of the form "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" where the names following ID:, PL:, 
PU:, etc. can be any user-chosen name conforming to the SAM specification. See SAM Format Specification Section 1.3 for 
details: The string passed here can be any string conforming to a header line for tag @RG, omitting the tag @RG itself, 
and using whitespace as separators for the line instead of TABs.

### --mark-duplicates

This filter marks the duplicate reads in the input file by setting bit 0x400 of their FLAG conforming to the SAM 
specification. For details on the algorithm and comparison to other tools, please consult our publication list.

### --mark-optical-duplicates file

When the --mark-duplicates filter is passed, one can also pass --mark-optical-duplicates. This option makes sure that 
optical duplicate marking is performed and a metrics file is generated that contains read statistics such as number of 
unmapped reads, secondary reads, duplicate reads, optical duplicates, library size estimate, etc. For details on the 
algorithm and comparison to other tools, please consult our publication list.

The metrics file generated by --mark-optical-duplicates is compatible with MultiQC for visualisation.

### --optical-duplicates-pixel-distance nr

This option allows specifying the pixel distance that is used for optical duplicate marking. This option is only usable 
in conjunction with --mark-optical-duplicates. The default value for the pixel distance is 100. In general, a pixel 
distance of 100 is recommended for data generated using unpatterned flowcells (e.g. HiSeq2500) and a pixel distance of 
2500 is recommended for patterned flowcells (e.g. NovaSeq/HiSeq4000). 

### --remove-duplicates

This filter removes all reads marked as duplicates. Duplicate reads are reads where their FLAG's bit 0x400 is set 
conforming the SAM specification. 

### --remove-optional-fields [all | list]

This filter removes for each alignment either all optional fields or all optional fields specified in the given list. 
The list of optional fields to remove has to be of the form "tag1, tag2, ..." where tag1, tag2, etc are the tags of the 
optional fields that need to be deleted.

### --keep-optional-fields [none | list]

This filter removes for each alignment either none of its optional fields, or all optional fields except those specified
 in the given list. The list of optional fields to keep has to be of the form "tag1, tag2, ..." where tag1, tag2, etc 
 are the tags of the optional fields that need to be kept in the output.

### --clean-sam

This filter fixes alignments in two ways:

* it soft clips alignments that hang off the end of its reference sequence
* it sets the mapping quality to 0 if the alignment is unmapped

This filter is similar to the CleanSam command of Picard.

### --bqsr recal-file

This filter performs base quality score recalibration. The recal-file is used for logging the recalibration tables 
computed during base recalibration. This file is compatible with MultiQC for visualisation.

There are additional elprep options that can be used for configuring the base quality score recalibration:

* --bqsr-reference elfasta (required)
* --known-sites list (optional)
* --quantize-levels nr (optional)
* --sqq list (optional)
* --max-cycle nr (optional)
* --target-regions bed-file (optional)

See detailed descriptions of these options next.

### --bqsr-reference elfasta-file

This option is used to pass a reference file for base quality score recalibration (--bqsr). The reference file must be 
in the .elfasta format, specific to elPrep. 

You can create an .elfasta file from a .fasta file using the elprep command fasta-to-elfasta. For example:

	elprep fasta-to-elfasta ucsc.hg19.fasta ucsc.hg19.elfasta
	

You only need to pass this option once if you are using both the --bqsr and --haplotypecaller options (which both 
require passing a reference file).

### --known-sites list 

This option is used to pass a number of known polymorphic sites that will be excluded during base recalibration (--bqsr)
.  The list is a list of files in the .elsites format, specific to elPrep. For example:

	--known-sites Mills_and_1000G_gold_standard.indels.hg19.elsites,dbsnp_137.hg19.elsites

You can create .elsites files from .vcf or .bed files using the vcf-to-elsites and bed-to-elsites parameters respectively. For example:

	elprep vcf-to-elsites dbsnp_137.hg19.vcf dbsnp_137.hg19.elsites

### --quantize-levels nr

This option is used to specify the number of levels for quantizing quality scores during base quality score recalibration (--bqsr). The default value is 0.

### --sqq list

This option is used to indicate to use static quantized quality scores to a given number of levels during base quality score recalibration (--bqsr). This list should be of the form "[nr, nr, nr]". The default value is [].

### --max-cycle nr

This option is used to specify the maximum cycle value during base quality score recalibration (--bqsr). The default value is 500.

### --target-regions bed-file

This option can be used to restrict which reads the base recalibration operates on by passing a .bed file that lists 
which genomic regions to consider. When this option is used, the reads that fall out of the specified regions are 
removed from the output .bam file. The option is for example used when processing exomes. 

This option produces a different result from --filter-non-overlapping-reads option. For the difference between both 
options and details on the algorithms, please consult our latest publication.

### --bqsr-tablename-prefix prefix

This option can be used to determine the prefix for the table names when logging the recalibration tables. The default
value ensures that the output is compatible with MultiQC. It is normally not necessary to set this option.

### --mark-optical-duplicates-intermediate file

This option is used in the context of filtering files created using the elprep split command. It is used internally by 
the elprep sfm command, but can be used when writing your own split/filter/merge scripts.

This option tells elPrep to perform optical duplicate marking and to write the result to an intermediate metrics file.
The intermediate metrics file generated this way can later be merged with other intermediate metrics files, see the 
merge-optical-duplicates-metrics command.

### --bqsr-tables-only table-file

This option is used in the context of filtering files created using the elprep split command. It is used internally by 
the elprep sfm command, but can be used when writing your own split/filter/merge scripts.

This option tells elPrep to perform base quality score recalibration and to write the result of the recalibration to an 
intermediate table file. This table file will need to be merged with other intermediate recalibration results during the
application of the base quality score recalibration. See the --bqsr-apply option.

### --bqsr-apply path

This option is used when filtering files created by the elprep split command. It is used internally by the elprep sfm 
command, and can be used when writing your own split/filter/merge scripts.

This option is used for applying base quality score recalibration on an input file. It expects a path parameter that 
refers to a directory that contains intermediate recalibration results for multiple files created using the 
--bqsr-tables-only option.

### --haplotypecaller vcf-file

This option performs variant calling for detecting germline SNPS and indels. The vcf-file is used for storing the vcf
 output. This file can be in gzipped format.

There are additional elprep options that can be used for configuring the haplotype variant caller:

* --bqsr-reference elfasta (required)
* --reference-confidence [GVCF | BP_RESOLUTION | NONE] (optional)
* --sample-name name (optional)
* --activity-profile igv-file (optional)
* --assembly-regions igv-file (optional)
* --assembly-region-padding nr (optional)
* --target-regions (optional)

See detailed descriptions of these options next.

### --bqsr-reference elfasta

This option is used to pass a reference file for variant calling (--haplotypecaller). The reference file must be in the 
.elfasta format, specific to elPrep. 

You can create an .elfasta file from a .fasta file using the elprep command fasta-to-elfasta. For example:

	elprep fasta-to-elfasta ucsc.hg19.fasta ucsc.hg19.elfasta

You only need to pass this option once if you are using both the --bqsr and --haplotypecaller options (which both 
require passing a reference file). 

### --reference-confidence [GVCF | BP_RESOLUTION | NONE]

This option is used to set the mode for emitting reference confidence scores when performing variant calling 
(--haplotypecaller). There are three options to choose from:

* GVCF (default): emit the GVCF format, i.e. the reference model is written with condensed non-variant blocks
* BP_RESOLUTION: the reference model is emitted site by site
* NONE: reference confidence calls are not emitted

### --sample-name name

The elPrep haplotypecaller (--haplotypecaller) only works for single samples. In case the input .bam file contains 
multiple samples, the --sample-name option can be used to select the sample reads on which to operate on.

### --activity-profile igv-file

Use this option to output the activity profile calculated by the haplotypecaller to the given file in IGV format.

### --activity-regions igv-file

This option can be used to output the assembly regions calculated by haplotypecaller to the speficied file in IGV format
.

### --assembly-region-padding nr

This option specfies the number of additional bases to include around each assembly region for variant calling.

### --target-regions bed-file

By default, the haplotypecaller scans the full genome for variants. Use this option to restrict the variant caller to 
specific regions by passing a .bed file. It is for example used when processing 
exomes.

You only need to pass this option once if you are using both the --bqsr and --haplotypcaller options.

## Sorting Command Options

### --sorting-order [keep | unknown | unsorted | queryname | coordinate]

This command option determines the order of the alignments in the output file. The command option must be followed by 
one of five possible orders:

1. *keep*: The original order of the input file is preserved in the output file. This is the default setting when the 
--sorting-order option is not passed. Some filters may change the order of the input, in which case elPrep forces a sort
 to recover the order of the input file.
2. *unknown*: The order of the alignments in the output file is undetermined, elPrep performs no sorting of any form. 
The order in the header of the output file will be *unknown*.
3. *unsorted*: The alignments in the output file are unsorted, elPrep performs no sorting of any form. The order in the 
header of the output file will be *unsorted*.
4. *queryname*: The output file is sorted according to the query name. The sort is enforced and guaranteed to be 
executed. If the original input file is already sorted by query name and you wish to avoid a sort with elPrep, use the 
*keep* option instead.
5. *coordinate*: The output file is sorted according to coordinate order. The sort is enforced and guaranteed to be 
executed. If the original input file is already sorted by coordinate order and you wish to avoid a sort with elPrep, 
use the *keep* option instead.

## Execution Command Options

### --nr-of-threads number

This command option sets the number of threads that elPrep uses during execution. The default number of threads is equal
 to the number of cpu threads.

It is normally not necessary to set this option. elPrep by default allocates the optimal number of threads.

### --timed

This command option is used to time the different phases of the execution of the elprep command, e.g. time spent on 
reading from file into memory, filtering, sorting, etc.

It is normally not necessary to set this option. It is only useful to get some details on where execution time is spent.

### --log-path path

This command option is used to specify a path where elPrep can store log files. The default path is the logs folder in 
your home path (~/logs).

## Format conversion tools

elPrep uses internal formats for representing .vcf, .bed, or .fasta files used by specific filter/sfm options. elPrep 
provides commands for creating these files from existing .vcf, .bed or .fasta files.

## Name

### elprep vcf-to-elsites - a commandline tool for converting a .vcf file to an .elsites file

## Synposis

	elprep vcf-to-elsites input.vcf output.elsites --log-path /home/user/logs

## Description

Converts a .vcf file to an .elsites file. Such a file can be passed to the --known-sites suboption of the --bqsr option.

## Options

### --log-path path

Sets the path for writing a log file.

## Name

### elprep bed-to-elsites - a commandline tool for converting a .bed file to an .elsites file

## Synposis

	elprep bed-to-elsites input.bed output.elsites --log-path /home/user/logs
	
## Description

Converts a .bed file to an .elsites file. Such a file can be passed to the --known-sites suboption of the --bqsr option.

## Options

### --log-path path

Sets the path for writing a log file.

## Name

### elprep fasta-to-elfasta - a commandline tool for converting a .fasta file to an .elfasta file

## Synopsis

	elprep fasta-to-elfasta input.fasta output.elfasta --log-path /home/user/logs
	
## Description

Converts a .fasta file to an .elfasta file. The --bqsr-reference suboption of the --bqsr and --haplotypecaller options 
requires an .elfasta file.

## Options

### --log-path path

Sets the path for writing a log file.

## Split and Merge tools

The elprep split command can be used to split up .sam files into smaller files that store the reads "per chromosome," or
 more precisely groups of chromosomes. elPrep determines the chromosomes by analyzing the sequence dictionary in the 
 header of the input file and generates a split file for groups of chromosomes that are roughly equal in size and that 
 stores all read pairs that map to that group of chromosomes. elPrep additionally creates a file for storing the 
 unmapped reads, and in the case of paired-end data, also a file for storing the pairs where reads map to different 
 chromosomes. elPrep also duplicates the latter pairs across chromosome files so that preparation pipelines have access 
 to all information they need to run correctly. Once processed, use the elprep merge command to merge the split files 
 back into a single .sam file.

Splitting the .sam file into smaller files for processing "per chromosome" is useful for reducing the memory pressure as
 these split files are typically significantly smaller than the  input file as a whole. Splitting also makes it possible
  to parallelize the processing of a single .sam file by distributing the different split files across different 
  processing nodes.

We provide an sfm command that executes a pipeline while silently using the elprep filter and split/merge tools. It is 
of course possible to write scripts to combine the filter and split/merge tools yourself.
We provide a recipe for writing your own split/filter/merge scripts on our github wiki.

## Name

### elprep sfm - a commandline tool for filtering and updating .sam/.bam files and variant calling "per chromosome"

## Synopsis

	elprep sfm input.sam output.sam 
	           --mark-duplicates --mark-optical-duplicates output.metrics 
	           --sorting-order coordinate 
	           --bqsr output.recal --bqsr-reference hg38.elfasta --known-sites dbsnp_138.hg38.elsites 
	           --haplotypecaller output.vcf.gz

	elprep sfm input.bam output.bam 
	           --mark-duplicates --mark-optical-duplicates output.metrics 
	           --sorting-order coordinate 
	           --bqsr output.recal --bqsr-reference hg38.elfasta --known-sites dbsnp_138.hg38.elsites 
	           --haplotypecaller output.vcf.gz
	           
	elprep sfm input.bam output.bam 
               --mark-duplicates --mark-optical-duplicates output.metrics 
               --sorting-order coordinate 
               --bqsr output.recal --bqsr-reference hg38.elfasta --known-sites dbsnp_138.hg38.elsites 
               --haplotypecaller output.vcf.gz

## Description

The elprep sfm command is a drop-in replacement for the elprep filter command that minimises the use of RAM. For this, 
it silently calls the elprep split and merge tools to split up the data "per chromosome" for processing, which requires 
less RAM than processing a .sam/.bam file as a whole (see Split and Merge tools).

## Options

The elprep sfm command has the same options as the elprep filter command, with the following additions.

### --intermediate-files-output-type [sam | bam]

This command option sets the format of the split files. By default, elprep uses the same format as the input file for 
the split files. Changing the intermediate file output type may improve either runtime (.sam) or reduce peak disk usage 
(.bam).

### --tmp-path

This command option is used to specify a path where elPrep can store temporary files that are created (and deleted) by 
the split and merge commands that are silently called by the elprep sfm command. The default path is the folder from 
where you call elprep sfm.

### --single-end

Use this command option to indicate the sfm command is processing single-end data. This information is important for the
 split/merge tools to operate correcly. For more details, see the description of the elprep split and elprep merge commands.

### --contig-group-size number

This command option is passed to the split tool.

The elprep split command groups the sequence dictionary entries for deciding how to split up the input data. The goal is
 to end up with groups of sequence dictionary entries (contigs) for which the total length (sum of LN tags) is roughly 
 the same among all groups. By default, the elprep split command identifies the sequence dictionary entry with the 
 largest length (LN) and chooses this as a target size for the groups. 

The --contig-group-size option allows configuring a specific group size. This size may be smaller than some of the 
sequence dictionary entries: elprep split will attempt to create as many groups of contigs of the chosen size, and 
contigs which are "too long" will be put in their own group.

Configuring the contig group size has an impact on how large the split files are that are generated by the elprep split 
command. Consequently, this also impacts how much RAM elprep uses for processing the split files. The default group size
 determines the minimum amount of RAM that is necessary to process a .sam/.bam file without information loss. 

The default value for the --contig-group-size option is 0. For this, elprep split makes groups based on the sequence 
dictionary entry with the biggest length (LN).

Choosing the value 1 for the --contig-group-size tells elprep split to split the data "per chromosome", i.e. a split 
file is created for each entry in the sequence dictionary.

## Name

### elprep split - a commandline tool for splitting .sam/.bam files per chromosome so they can be processed without information loss

## Synopsis

	elprep split [sam-file | /path/to/input/] /path/to/output/ --output-prefix "split-sam-file" --output-type sam 
	--nr-of-threads $threads --single-end

## Description

The elprep split command requires two arguments: 1) the input file or a path to multiple input files and 2) a path to a directory where elPrep can store the split files. The input file(s) can be .sam or .bam. It is also possible to use /dev/stdin as the input for using Unix pipes. There are no structural requirements on the input file(s) for using elprep split. For example, it is not necessary to sort the input file, nor is it necessary to convert to .bam or index the input file.

Warning: If you pass a path to multiple input files to the elprep split command, elprep assumes that they all have the same (or compatible) headers, and just picks the first header it finds as the header for all input files. elprep currently does not make an attempt to resolve potential conflicts between headers, especially with regard to the @SQ, @RG, or @PG header records. We will include proper merging of different SAM/BAM files in a future version of elprep. In the meantime, if you need proper merging of SAM/BAM files, please use samtools merge, Picard MergeSamFiles, or a similar tool. (If such a tool produces SAM file as output, it can be piped into elprep using Unix pipes.)

elPrep creates the output directory denoted by the output path, unless the directory already exists, in which case elPrep may override the existing files in that directory. Please make sure elPrep has the correct permissions for writing that directory.

By default, the elprep split command assumes it is processing pair-end data. The flag --single-end can be used for processing single-end data. The output will look different for paired-end and single-end data.

### Paired-end data (default)

The split command outputs two types of files:

1. a subdirectory "/path/to/output/splits/". The split command groups the entries in the sequence dictionary of the input file and creates a file for each of these groups containing all reads that map to that group.
2. a "/path/to/output/output-prefix-spread.output-type" file containing all reads of which the mate maps to a different entry in the sequence dictionary of the input file.

To process the files created by the elprep split command, one needs to call the elprep filter command for each entry in the path/to/output/splits/ directory as well as the /path/to/output/output-prefix-spread.output-type file. The output files produced this way need to be merged with the elprep merge command. This is implemented by the elprep sfm command.

### Single-end data (--single-end)

The split command groups entries in the sequence dictionary of the input file and creates a file for each of these groups that contain all reads that map to that group, and writes those files to the /path/to/output/ directory.

To process the files created by the elprep split --single-end command, one needs to call the elprep filter command for each entry in the /path/to/output/ directory. The output files produces this way need to be merged with the elprep merge command. This is implemented by the elprep sfm command.

## Options

### --output-prefix name

The split command groups entries in the sequence dictionary. The purpose of this grouping is to create groups of which the lengths of the entries (LN tags) add up to roughly the same size. 

The names of the split files created by elprep split are generated by combing a prefix and a chromosome group name. The --output-prefix option sets that prefix. 

For example, if the prefix is "NA12878", and the sfm command creates N groups for the sequence dictionary of the input file, then the names of the split files will be "NA12878-group1.output-type", "NA12878-group2.output-type", "NA12878-group3.output-type", and so on. A seperate file for the unmapped reads is created, e.g. "NA12878-unmapped.output-type".

If the user does not specify the --output-prefix option, the name of the input file, minus the file extension, is used as a prefix.

### --output-type [sam | bam]

This command option sets the format of the split files. By default, elprep uses the same format as the input file for the split files.

### --nr-of-threads number

This command option sets the number of threads that elPrep uses during execution for parsing/outputting .sam/.bam data. The default number of threads is equal to the number of cpu threads.

It is normally not necessary to set this option. elPrep by default allocates the optimal number of threads.

### --single-end

When this command option is set, the elprep split command will treat the data as single-end data. When the option is not used, the elprep split command will treat the data as paired-end data.

### --log-path path

Sets the path for writing a log file.

### --contig-group-size number

The elprep split command groups the sequence dictionary entries for deciding how to split up the input data. The --contig-group-size options allows configuring a specific group size. See the description of --contig-group-size for the elprep sfm command for more details.

## Name

### elprep merge - a commandline tool for merging .sam/.bam files created by elprep split

## Synopsis

	elprep merge /path/to/input/ sam-output-file --nr-of-threads $threads --single-end

## Description

The elprep merge command requires two arguments: a path to the files that need to be merged, and an output file. Use this command to merge files created with elprep split. The output file can be .sam or .bam. It is also possible to use /dev/stdout as output when using Unix pipes for connecting other tools.

## Options

### --nr-of-threads number

This command option sets the number of threads that elPrep uses during execution for parsing/outputting .sam/.bam data. The default number of threads is equal to the number of cpu threads.

It is normally not necessary to set this option. elPrep by default allocates the optimal number of threads.

### --single-end

This command option tells the elprep merge command to treat the data as single-end data. When this option is not used, elprep merge assumes the data is paired-end, expecting the data is merging to be generated by the elprep split command accordingly.

### --log-path path

Sets the path for writing a log file.

## Name

### elprep merge-optical-duplicate-metrics - a commandline tool for merging intermediate metrics files created by the --mark-optical-duplicates-intermediate option

## Synopsis

	elprep merge-optical-duplicates-metrics input-file output-file metrics-file /path/to/intermediate/metrics --remove-duplicates

## Description

The elprep merge-optical-duplicates-metrics command requires four arguments: 
the names of the original input and output .sam/.bam files for which the metrics are calculated, 
the metrics file to which the merged metrics should be written, and a path to the intermediate metrics files that need 
to be merged (and were generated using --mark-optical-duplicates-intermediate). 

## Options

### --nr-of-threads number

This command option sets the number of threads that elPrep uses during execution for parsing/outputting .sam/.bam data. The default number of threads is equal to the number of cpu threads.

It is normally not necessary to set this option. elPrep by default allocates the optimal number of threads.

## --remove-duplicates

Pass this option if the metrics were generated for a file for which the duplicates were removed. This information will 
be included in the merged metrics file.

# Extending elPrep

If you wish to extend elPrep, for example by adding your own filters, please consult our [API documentation](https://godoc.org/github.com/ExaScience/elprep).

# Acknowledgements

Many thanks for testing, bug reports, or contributions:

Amin Ardeshirdavani

Pierre Bourbon

Benoit Charloteaux

Richard Corbett

Didier Croes

Matthias De Smet

Keith James

Leonor Palmeira

Joke Reumers

Geert Vandeweyer
