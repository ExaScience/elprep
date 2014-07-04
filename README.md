# Overview

elPrep is a high-performance tool for preparing .sam/.bam files for variant calling in sequencing pipelines. It can be used as a drop-in replacement for SAMtools/Picard, and was extensively tested with the GATK best practices pipeline for variant analysis.

elPrep is designed as an in-memory and multi-threaded application to fully take advantage of the processing power available with modern servers. Its software architecture is based on functional programming techniques, which allow to easily compose multiple alignment filters and perform optimizations such as loop merging.

The table below shows the execution time of the preparation phases in a whole-genome-sequencing pipeline for NA12878 on a 4x10-core Intel Xeon E7-4870 server with 512GB RAM. The preparation phases in this example include filtering unmapped reads, sorting for coordinate order, marking duplicates, adding read group information, and reordering the reference sequence dictionary for GATK. 

The table shows that preparation phase with elPrep is 25x faster than using SAMtools/Picard. The output of elPrep is 100% equivalent to output produced by SAMtools/Picard.

	Preparation of NA12878 on one 4x10-core Intel Xeon E7-4870 server
	
							Number of Threads	Execution Time	
		SAMtools/Picard					   40	   29h 44m 55s				
		elPrep							   40	    1h 11m 38s

elPrep is being developed at the [ExaScience Life Lab](http://www.exascience.com), a collaboration between Imec, Intel, and Janssen Pharmaceutica.

# Advantages

The advantages of elPrep include:

* efficient multi-threaded execution
* operates completely in-memory, no intermediate files are generated 
* 100% equivalent output to output produced by Picard-tools and SAMtools for overlapping functionality
* compatible with existing tools such as GATK, SAMtools, and Picard-tools
* extensible implementation

# Availability 

elPrep is released as an open-source project under a BSD 3-Clause License (BSD 2.0). We also provide a download of a precompiled binary to which a specific license applies.

## Binaries

You can download a precompiled binary of elPrep [here](https://github.com/ExaScience/elprep/releases) upon accepting the license agreement. This binary was created using the 64-bit LispWorks 6.1 Professional Edition for Linux. 

## GitHub

The elPrep source code is freely available on GitHub. elPrep is implemented in Common Lisp using the LispWorks compiler for Linux.

	elPrep GitHub clone URL:
	
		https://github.com/ExaScience/elprep.git

## Dependencies

elPrep works with the .sam and .bam formats as input/output. To use .bam SAMtools must be installed:

	
		http://samtools.sourceforge.net

The .cram format will be supported once the successor of samtools-0.1.19 (which supports .cram) is officially released.

The elPrep implementation depends on the Claws and cl-date-time-parser libraries:

	Claws GitHub clone URL:
	
		https://github.com/pcostanza/claws.git
		
	cl-date-time-parser GitHub clone URL:
	
		https://github.com/ExaScience/cl-date-time-parser

The elPrep implementation also depends on the string-case library which is available through the [quicklisp](http://www.quicklisp.org) package manager. 

Installing these libraries is only necessary if you wish to build elPrep yourself. It is not necessary to use the elPrep binary we provide above.

## Building

A build script is provided with the source code if you wish to build elPrep yourself.

	Building elPrep using LispWorks:
	
		lispworks -build save-elprep-script.lisp

Please ensure that the elPrep repository and its dependencies are visible to the asdf path of your LispWorks compiler.

## Compatibility 

The output of elPrep is compatible (as input) with:

* gatk 3.1-1
* gatk 2.7.2
* gatk 1.6
* samtools-0.1.19
* picard-tools-1.113

elPrep may be compatible with other versions of these tools as well, but this has not been tested.

elPrep has been developed for Linux and has not been tested for other operating systems. We have tested elPrep with the following Linux distributions:

* Ubuntu 12.04.3 LTS
* Manjaro Linux
* Red Hat Enterprise Linux 6.4 and 6.5

# Memory Requirements

## RAM

elPrep is designed to operate in-memory, i.e. data is stored in RAM during computation. As long as you do not use the in-memory sort or mark duplicates filter, elPrep operates as a streaming tool and peek memory use is limited to a few GB.

The in-memory sort and mark duplicates filter require keeping the entire input file in memory and therefore use an amount of RAM that is proportional to the size of the input file. As a rule of thumb, elPrep requires 6x times more RAM memory than the size of the input file in .sam format when it is used for sorting or marking duplicates. 

If your machine has less RAM than your input requires this way, you have a couple of options. One solution is to split up your input file per chromosomal region, which is for example also done to run analyses in a distributed setup. A variety of tools such as SAMtools provide this functionality.

Alternatively, you may configure a disk to extend the swap space of your server. Using swap space may have a penalty on execution time for a job compared to running the same job fully in RAM, but it does not change how elPrep is used. If configuring additional swap space is not an option, you may also run elPrep with the gc execution option. This option triggers more aggressive garbage collection during execution and may save peek memory use, but may significantly slow down execution compared to running fully in-memory. See our manual reference pages for more details. 

## Disk Space

elPrep does not write any intermediate files, and therefore does not require additional (peek) disk space beyond what is needed for storing the input and output files.  

# Mailing List

Use the Google [forum](https://groups.google.com/d/forum/elprep) for discussions. You need a Google account to subscribe through the forum URL. You can also subscribe without a Google account by sending an email to elprep+subscribe@googlegroups.com.

# Citing elPrep  

As of today there is no publication on elPrep yet. Please cite the ExaScience Life Lab with a link to the GitHub repository:

		https://github.com/ExaScience/elprep

# Demo

We have a seperate [GitHub repository](https://github.com/ExaScience/elprep-demo) with demo scripts that use elPrep for procesing a subset of NA12878 that maps to chromosome 22.

# Manual Reference Pages

## Name

### elprep - a commandline tool for filtering and updating .sam/.bam files in preparation for variant calling

## Synopsis

	elprep input.sam output.sam --filter-unmapped-reads --replace-refernce-sequences $gatkdict --replace-read-group "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" --mark-duplicates --sorting-order coordinate --nr-of-threads $threads
	
	elprep input.bam output.bam --filter-unmapped-reads --replace-refernce-sequences $gatkdict --replace-read-group "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" --mark-duplicates --sorting-order coordinate --nr-of-threads $threads
	
	elprep input.cram output.cram --filter-unmapped-reads --replace-refernce-sequences $gatkdict --replace-read-group "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" --mark-duplicates --sorting-order coordinate --nr-of-threads $threads
	
	elprep /dev/stdin /dev/stdout --filter-unmapped-reads --replace-refernce-sequences $gatkdict --replace-read-group "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" --mark-duplicates --sorting-order coordinate --nr-of-threads $threads

## Description

The elprep command requires two arguments: the input file and the output file. The input/output format can be .sam or .bam. (.cram will be supported in the near future). To use .bam (or .cram), SAMtools must be installed. elPrep determines the format by looking at the file extension. elPrep also allows to use /dev/stdin and /dev/stdout as respective input or output sources for using unix pipes. When doing so, elPrep assumes the input and output are in .sam format.

The elprep commandline tool has three types of command options: filters, which implement actual .sam/.bam/.cram manipulations, sorting options, and execution-related options, e.g. for setting the number of threads. For optimal performance, issue a single elprep call that combines all filters you wish to apply. 

The order in which command options are passed is ignored. For optimal performance, elPrep always applies filters in the following order:

1. filter-unmapped-reads
2. clean-sam
3. replace-reference-sequences
4. replace-read-group
5. mark-duplicates

Sorting is done after filtering.

## Filter Command Options  

### --replace-reference-sequences file

This filter is used for replacing the header of a .sam/.bam/.cram file by a new header. The new header is passed as a single argument following the command option. The format of the new header can either be a .dict file, e.g. ucsc.hg19.dict from the GATK bundle, or another .sam/.bam/.cram file from which you wish to extract the new header. 

All alignments in the input file that do not map to a chromosome that is present in the new header, are removed. Therefore there must be some overlap between the old and the new header for this command option to be meaningful. The option is typically used to reorder the reference sequence dictionary in the header, e.g. to reflect the order required by GATK.

Replacing the header of a .sam/.bam/.cram file may destroy the sorting order of the file. In this case, the sorting order in the header is set to "unknown" by elPrep in the output file (cf. the so tag). 

### --filter-unmapped-reads [strict]

Removes all alignments in the input file that are unmapped. An alignment is determined unmapped when bit 0x4 of its FLAG is set, conform the SAM specification.

There user may additionally speficy the *strict* option:

* strict: Also removes alignments where the mapping position (POS) is 0 or where the reference sequence name (RNAME) is *. Such alignments are considered unmapped by the SAM speficication, but some alignment programs may not mark the FLAG of those alignments as unmapped. This option is recommended when you are targetting older versions of GATK (cf. GATK 1.6).

### --replace-read-group read-group-string

This filter is replaces or adds read groups to the alignments in the input file. This command option takes a single argument, a string of the form "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" where the names following ID:, PL:, PU:, etc. can be any user-chosen name conform the SAM specification. See SAM Format Specification Section 1.3 for details: The string passed here can be any string conforming to a header line for tag @RG, omitting the tag @RG itself, and using whitespace as separators for the line instead of TABs.

### --mark-duplicates [remove]

This filter marks the duplicate reads in the input file by setting bit 0x400 of their FLAG conform the SAM specification. The algorithm underlying this option is the same as the one used in picard-tools.

The user may additionally pass the *remove* flag so that the duplicate reads are not written to the output file.

### --clean-sam

This filter fixes alignments in two ways:

* it soft clips alignments that hang off the end of its reference sequence
* it sets the mapping quality to 0 if the alignment is unmapped

This filter is similar to the CleanSam command of Picard-tools.

## Sorting Command Options

### --sorting-order [keep | unknown | unsorted | queryname | coordinate]

This command option determines the order of the alignments in the output file. The command option is followed by one of five possible orders:

1. *keep*: The original order of the input file is preserved in the output file. This is the default setting. Some filters may change the order of an input file, in which case elPrep forces a sort to preserve the order of the input file.
2. *unknown*: The order of the alignments in the output file is undetermined, elPrep performs no sorting of any form. The order in the header of the output file will be *unknown*.
3. *unsorted*: The alignments in the output file are unsorted, elPrep performs no sorting of any form. The order in the header of the output file will be *unsorted*.
4. *queryname*: The output file is sorted according to query name. The sort is enforced and guaranteed to be executed. If the original input file is already sorted by query name and you wish to avoid a sort with elPrep, use the *keep* option instead.
5. *coordinate*: The output file is sorted according to coordinate order. The sort is enforced and guaranteed to be executed. If the original input file is already sorted by coordinate order and you wish to avoid a sort with elPrep, use the *keep* option instead.

## Execution Command Options

### --nr-of-threads number

This command option sets the number of threads that elPrep uses during execution. The default number of threads is 1.

### --gc-on [0 | 1 | 2 ]

This option configures garbage collection during the elPrep execution. By default, elPrep is optimised for performance, but not for memory use. Therefore, elPrep may use a large amount of RAM and/or virtual memory / swap space during its execution (see Memory Requirements). If you want to use elPrep in memory-constrained environments this option may be helpful.

Specifically, if you use the mark-duplicates and/or the sorting option in elPrep, memory use will be high even with a high gc-on setting, at the expense of a significantly slower execution time. However, if you do not use the mark-duplicates option, and/or use the sorting-order unsorted option, setting gc-on high can significantly decrease memory use, possibly at the expense of a slower runtime.

There are three options:

* 0: elPrep avoids garbage collection. This is the default setting. Use this option when enough RAM and virtual memory / swap space is available.
* 1: elPrep performs a garbage collect after the filtering and sorting phases. This can reduce peek memory use, but slows down execution somewhat.
* 2: elPrep performs garbage collection interleaved with the execution at regular intervals. This may reduce peek memory significantly, but may slow down execution considerably.

If the option is not passed explicitly, elPrep assumes —gc-on 0 is intended.

<!--### --timed

When this option is passed, elPrep times and prints the execution spent per phase --filtering, sorting, I/O-- to the standard output.
-->

