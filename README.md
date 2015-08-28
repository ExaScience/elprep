# Overview

elPrep is a high-performance tool for preparing .sam/.bam/.cram files for variant calling in sequencing pipelines. It can be used as a drop-in replacement for SAMtools/Picard, and was extensively tested with different pipelines for variant analysis with GATK. The key advantage of elPrep is that it only performs a single-pass to process a .sam/.bam/.cram file, independent of the number of filters that need to be applied in a particular pipeline.

elPrep is designed as an in-memory and multi-threaded application to fully take advantage of the processing power available with modern servers. Its software architecture is based on functional programming techniques, which allows easy composition of multiple alignment filters and optimizations such as loop merging. To make this possible, elPrep proposes a couple of new algorithms. For example, for duplicate marking we devised an algorithm that expresses the operation as a single-pass filter using memoization techniques and hierarchical hash tables for efficient parallel synchronisation.

The main advantage of elPrep is very fast execution times on high-end backend servers, as is available for example through Amazon cloud computing services or custom server setups. We do not recommend using elPrep on laptops, desktops, or low-end servers. Please consult the system requirements below for more details.

The table below shows the execution time of the preparation phases in a whole-genome-sequencing pipeline for NA12878 on a 4x10-core Intel Xeon E7-4870 server with 512GB RAM. The sam file is split up and processed per chromosome. The preparation phases in this example include filtering unmapped reads, sorting for coordinate order, marking duplicates, adding read group information, and reordering the reference sequence dictionary for GATK.

The table shows that the preparation phase with elPrep is 25x faster than with SAMtools and Picard. The output of elPrep is 100% equivalent to output produced by SAMtools and Picard.

	Preparation of NA12878 on one 4x10-core Intel Xeon E7-4870 server

							Number of Threads	Execution Time
		SAMtools/Picard					   40	   29h 44m 55s
		elPrep							   40	    1h 11m 38s

We also excuted the same pipeline on a 2x8-core hyperthreaded Intel Xeon E5-2680 server with 256GB RAM. On this configuration, elPrep is 10.5 times faster than SAMtools and Picard.

	Preparation of NA12878 on one 2x8-core hyperthreaded Intel Xeon E5-2680 server

							Number of Threads	Execution Time
		SAMtools/Picard					   32	   19h 19m 29s
		elPrep							   32	    1h 54m 51s

elPrep is being developed at the [ExaScience Life Lab](http://www.exascience.com), a collaboration between Imec, Intel, and Janssen Pharmaceutica. For questions, use our mailing list (below) or contact [Charlotte Herzeel](https://github.com/caherzee) directly.

# Advantages

The advantages of elPrep include:

* efficient multi-threaded execution
* operates completely in-memory, no intermediate files are generated
* 100% equivalent output to output produced by SAMtools and Picard for overlapping functionality
* compatible with existing tools such as GATK, SAMtools, and Picard
* modular implementation

# Availability

elPrep is released as an open-source project under a BSD 3-Clause License (BSD 2.0). We also provide a download of a precompiled binary to which a specific license applies.

## Binaries

You can download a precompiled binary of elPrep [here](https://github.com/ExaScience/elprep/releases) upon accepting the license agreement. This binary was created using the 64-bit LispWorks 7.0 Professional Edition for Linux.

## GitHub

The elPrep source code is freely available on GitHub. elPrep is implemented in Common Lisp using the LispWorks and SBCL compilers for Linux.

	elPrep GitHub clone URL:

		https://github.com/ExaScience/elprep.git

## Dependencies

elPrep works with the .sam, .bam, and .cram formats as input/output. To use .bam or .cram, SAMtools must be installed in addition to the elPrep binary:


		http://www.htslib.org

The .cram format is only supported for samtools-1.0 and up. It is in principle possible to use another tool for .bam and .cram compression, see our example scripts for more details.

We provide example scripts of how to combine elPrep with GNU parallel, for example for distributed execution. Please install GNU parallel to run these scripts:

		http://www.gnu.org/software/parallel/

## Building

The following is only relevant information if you wish to build elPrep yourself. It is not necessary to use the elPrep binary we provide above.

A build script is provided with the source code if you wish to build elPrep yourself.

	Building elPrep using LispWorks:

		sh make-lispworks-binary.sh

	Building elPrep using SBCL:

		sh make-sbcl-binary.sh

The elPrep implementation depends on the Claws and cl-date-time-parser libraries:

	Claws GitHub clone URL:

		https://github.com/pcostanza/claws.git

	cl-date-time-parser GitHub clone URL:

		https://github.com/ExaScience/cl-date-time-parser

The elPrep implementation also depends on the string-case library which is available through the [quicklisp](http://www.quicklisp.org) package manager.

Please ensure that the elPrep repository and its dependencies are visible to the asdf path of your LispWorks or SBCL compiler.

### SBCL-specific notes

Data sets in sequencing pipelines can become very large. Therefore you may need to tweak the default heap size of SBCL to be able to use an SBCL-based version of elPrep.

Specifically, you may need to increase the value of \*gencgc-card-bytes\* before building SBCL itself. In the source distribution of SBCL, you can edit this variable in the file src/compiler/target/backend-parms.lisp, where target is the platform where you want to run elPrep, so for example in src/compiler/x86-64/backend-parms.lisp.

You may also need to specify a larger dynamic space size when building elPrep. You can do this by building elPrep as follows:

	 	sbcl --dynamic-space-size NNNNNN ...

See the SBCL build script in make-sbcl-binary.sh, which includes an already relatively large setting for the dynamic space size.

## Compatibility

The output of elPrep is compatible (as input) with:

* gatk 3.1-1
* gatk 2.7.2
* gatk 1.6
* samtools-0.1.19
* samtools-1.0, samtools-1.1, samtools-1.2
* picard-tools-1.113

elPrep may be compatible with other versions of these tools as well, but this has not been tested.

elPrep has been developed for Linux and has not been tested for other operating systems. We have tested elPrep with the following Linux distributions:

* Ubuntu 12.04.3 LTS
* Manjaro Linux
* Red Hat Enterprise Linux 6.4 and 6.5

## Docker image

You can download a precompiled Docker containter for elPrep [here](https://hub.docker.com/r/caherzee/elprep/). For building the Docker image yourself, execute the following command in the directory that contains the elPrep Dockerfile:

	docker build -t elprep .

For information on how to use the elPrep Docker container, please consult the documentation on [the elprep page on the docker hub](https://hub.docker.com/r/caherzee/elprep/). 

# Memory Requirements

## RAM

elPrep is designed to operate in memory, i.e. data is stored in RAM during computation. As long as you do not use the in-memory sort or mark duplicates filter, elPrep operates as a streaming tool, and peak memory use is limited to a few GB.

The in-memory sort and mark duplicates filter require keeping the entire input file in memory, and therefore use an amount of RAM that is proportional to the size of the input file. As a rule of thumb, elPrep requires 6x times more RAM memory than the size of the input file in .sam format when it is used for sorting or marking duplicates. 

elPrep provides a tool for splitting .sam files "per chromosome," and guarantees that processing these split files and then merging the results does not lose information when compared to processing a .sam file as a whole. Using the split/merge tool greatly reduces the RAM required to process a .sam file, but it comes at the cost of an additional processing step. For example, for the hg19 reference, the amount of RAM required is roughly a factor 0.6 of the input file in .sam format when using the split and merge tools.

For example, in our experience, for whole-genome sequencing (+500GB .sam or 99GB .bam), a server with at least 256GB and splitting the data per chromosome is required. For exome sequencing (+-40GB .sam or +-8GB .bam), processing the data per chromsome requires a server with 30GB RAM.

If your machine has less RAM than your input requires (after splitting), you have a couple of options. One solution is to (further) split up your input file per genomic region. You can do this by splitting up the reference sequence dictionary into smaller genomic regions. We plan to add support for BED files to facilitate this. (Please contact us if you need this urgently.)

Alternatively, you may configure a disk to extend the swap space of your server. Using swap space may have a penalty on execution time for a job compared to running the same job fully in RAM, but it does not change how elPrep is used. If configuring additional swap space is not an option, you may also run elPrep with the gc execution option. This option triggers more aggressive garbage collection during execution and may save peak memory use, but may significantly slow down execution compared to running fully in memory. See our manual reference pages below for more details.

## Disk Space

elPrep by default does not write any intermediate files, and therefore does not require additional (peak) disk space beyond what is needed for storing the input and output files. If you use the elPrep split and merge tools, elPrep requires additional disk space equal to the size of the input file.

# Mailing List and Contact

Use the Google [forum](https://groups.google.com/d/forum/elprep) for discussions. You need a Google account to subscribe through the forum URL. You can also subscribe without a Google account by sending an email to elprep+subscribe@googlegroups.com.

You can also contact [Charlotte Herzeel](https://github.com/caherzee) directly.

# Citing elPrep

Please cite the following article:

Herzeel C, Costanza P, Decap D, Fostier J, Reumers J (2015) elPrep: High-Performance Preparation of Sequence Alignment/Map Files for Variant Calling. PLoS ONE 10(7): e0132868. doi:10.1371/journal.pone.0132868

Please use the LispWorks version when running and reporting benchmarks, since the SBCL version is in general slower. If performance is below your expectations, please contact us first before reporting your results.

# Demo

We have a seperate [GitHub repository](https://github.com/ExaScience/elprep-demo) with demo scripts that use elPrep for procesing a subset of NA12878 that maps to chromosome 22.

# Manual Reference Pages

## Name

### elprep - a commandline tool for filtering and updating .sam/.bam/.cram files in preparation for variant calling

## Synopsis

	elprep input.sam output.sam --filter-unmapped-reads --replace-reference-sequences $gatkdict --replace-read-group "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" --mark-duplicates --sorting-order coordinate --nr-of-threads $threads

	elprep input.bam output.bam --filter-unmapped-reads --replace-reference-sequences $gatkdict --replace-read-group "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" --mark-duplicates --sorting-order coordinate --nr-of-threads $threads

	elprep input.cram output.cram --filter-unmapped-reads --replace-reference-sequences $gatkdict --replace-read-group "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" --mark-duplicates --sorting-order coordinate --nr-of-threads $threads

	elprep /dev/stdin /dev/stdout --filter-unmapped-reads --replace-reference-sequences $gatkdict --replace-read-group "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" --mark-duplicates --sorting-order coordinate --nr-of-threads $threads

## Description

The elprep command requires two arguments: the input file and the output file. The input/output format can be .sam, .bam. or .cram. To use .bam or .cram, SAMtools must be installed. elPrep determines the format by looking at the file extension. elPrep also allows to use /dev/stdin and /dev/stdout as respective input or output sources for using Unix pipes. When doing so, elPrep assumes the input and output are in .sam format.

The elprep commandline tool has three types of command options: filters, which implement actual .sam/.bam/.cram manipulations, sorting options, and execution-related options, for example for setting the number of threads. For optimal performance, issue a single elprep call that combines all filters you wish to apply.

The order in which command options are passed is ignored. For optimal performance, elPrep always applies filters in the following order:

1. filter-unmapped-reads
2. clean-sam
3. replace-reference-sequences
4. replace-read-group
5. mark-duplicates

Sorting is done after filtering.

### Unix pipes

elPrep is compatible with Unix pipes and allows using /dev/stdin and /dev/stdout as input or output sources. elPrep assumes that input and output on /dev/stdin and /dev/stdout are in .sam format.

### Using .bam or .cram

elPrep by default uses SAMtools for compression and decompression of .bam and .cram files. When the user specifies the input or ouput to be in .bam or .cram, a call to SAMtools is generated to compress or decompress the data. The SAMtools and elPrep commands are connected via Unix pipes to avoid generating intermediate files.

For example, when the user executes the following elPrep command:

	 elprep input.bam output.bam

The actual command that is executed is similar to the following:

	samtools view -h input.bam /dev/stdout | elprep /dev/stdin /dev/stdout | samtools view -b -o output.bam /dev/stdin

elPrep may pass additional parameters to SAMtools which are not shown above, for example for setting the number of compression threads. You are, of course, free to specify the SAMtools commands for .bam/.cram conversion manually in your own scripts. The elPrep download includes the Python scripts elprep\_io\_wrapper.py and elprep.py that illustrate how to do this.

### .cram Conversion Options

When specifying that the output is in .cram format, it is required to specify the reference file to be used for the .cram compression via one of the --reference-t or --reference-T options.

### --reference-T reference-fasta

A fasta format reference file used by SAMtools for .cram compression, optionally compressed with bgzip and indexed by samtools faidx. elPrep uses it to fill in the "-T" option when calling the samtools view command for converting to .cram.

For example, when the user executes the following elPrep command:

	 elprep input.sam output.cram --reference-T ucsc.hg19.fasta

The actual command that is executed is similar to the following:

	 elprep input.sam /dev/stdout | samtools view -C -T ucsc.hg19.fasta -o output.cram /dev/stdin

elPrep may pass additional options to the SAMtools command shown above.

See the [SAMtools manual](http://www.htslib.org/man/samtools/) for more documentation on .cram conversion and the "-T" option to the view command in particular.

### --reference-t reference-fai

A tab-delimited file, where a first column lists the reference names and a second column lists the lengths of those references; for example, a .fai file generated with samtools faidx. elPrep uses it to fill in the "-t" option when calling the samtools view command for converting to .cram.

For example, when the user executes the following elPrep command:

	 elprep input.sam output.cram --reference-t ucsc.hg19.fai

The actual command that is executed is similar to the following:

	 elprep input.sam /dev/stdout | samtools view -C -t ucsc.hg19.fai -o output.cram /dev/stdin

elPrep may pass additional options to the SAMtools command shown above.

See the [SAMtools manual](http://www.htslib.org/man/samtools/) for more documentation on .cram conversion and the "-t" option to the view command in particular.

## Filter Command Options

### --replace-reference-sequences file

This filter is used for replacing the header of a .sam/.bam/.cram file by a new header. The new header is passed as a single argument following the command option. The format of the new header can either be a .dict file, for example ucsc.hg19.dict from the GATK bundle, or another .sam/.bam/.cram file from which you wish to extract the new header.

All alignments in the input file that do not map to a chromosome that is present in the new header are removed. Therefore, there should be some overlap between the old and the new header for this command option to be meaningful. The option is typically used to reorder the reference sequence dictionary in the header, for example to reflect the order required by GATK.

Replacing the header of a .sam/.bam/.cram file may destroy the sorting order of the file. In this case, the sorting order in the header is set to "unknown" by elPrep in the output file (cf. the 'so' tag).

### --filter-unmapped-reads [strict]

Removes all alignments in the input file that are unmapped. An alignment is determined unmapped when bit 0x4 of its FLAG is set, conforming to the SAM specification.

There user may additionally specify the *strict* option:

* strict: Also removes alignments where the mapping position (POS) is 0 or where the reference sequence name (RNAME) is *. Such alignments are considered unmapped by the SAM specification, but some alignment programs may not mark the FLAG of those alignments as unmapped. This option is recommended when you are targeting older versions of GATK (cf. GATK 1.6).

### --replace-read-group read-group-string

This filter is replaces or adds read groups to the alignments in the input file. This command option takes a single argument, a string of the form "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" where the names following ID:, PL:, PU:, etc. can be any user-chosen name conforming to the SAM specification. See SAM Format Specification Section 1.3 for details: The string passed here can be any string conforming to a header line for tag @RG, omitting the tag @RG itself, and using whitespace as separators for the line instead of TABs.

### --mark-duplicates [remove]

This filter marks the duplicate reads in the input file by setting bit 0x400 of their FLAG conforming to the SAM specification. The criteria underlying this option are the same as the ones used in Picard.

The user may additionally pass the *remove* flag so that the duplicate reads are not written to the output file.

### --clean-sam

This filter fixes alignments in two ways:

* it soft clips alignments that hang off the end of its reference sequence
* it sets the mapping quality to 0 if the alignment is unmapped

This filter is similar to the CleanSam command of Picard.

## Sorting Command Options

### --sorting-order [keep | unknown | unsorted | queryname | coordinate]

This command option determines the order of the alignments in the output file. The command option is followed by one of five possible orders:

1. *keep*: The original order of the input file is preserved in the output file. This is the default setting. Some filters may change the order of the input, in which case elPrep forces a sort to recover the order of the input file.
2. *unknown*: The order of the alignments in the output file is undetermined, elPrep performs no sorting of any form. The order in the header of the output file will be *unknown*.
3. *unsorted*: The alignments in the output file are unsorted, elPrep performs no sorting of any form. The order in the header of the output file will be *unsorted*.
4. *queryname*: The output file is sorted according to the query name. The sort is enforced and guaranteed to be executed. If the original input file is already sorted by query name and you wish to avoid a sort with elPrep, use the *keep* option instead.
5. *coordinate*: The output file is sorted according to coordinate order. The sort is enforced and guaranteed to be executed. If the original input file is already sorted by coordinate order and you wish to avoid a sort with elPrep, use the *keep* option instead.

## Execution Command Options

### --nr-of-threads number

This command option sets the number of threads that elPrep uses during execution. The default number of threads is 1.

### --gc-on [0 | 1 | 2 ]

This option configures garbage collection during the elPrep execution. By default, elPrep is optimised for performance, but not for memory use. Therefore, elPrep may use a large amount of RAM and/or virtual memory / swap space during its execution (see Memory Requirements). If you want to use elPrep in memory-constrained environments, this option may be helpful.

Specifically, if you use the mark-duplicates and/or the sorting option in elPrep, memory use will be high even with a high gc-on setting, at the expense of a significantly slower execution time. However, if you do not use the mark-duplicates option, and/or use the sorting-order unsorted option, setting gc-on high can significantly decrease memory use, possibly at the expense of a slower runtime.

There are three options:

* 0: elPrep avoids garbage collection. This is the default setting. Use this option when enough RAM and virtual memory / swap space is available.
* 1: elPrep performs a garbage collect after the filtering and sorting phases. This can reduce peak memory use, but slows down execution somewhat.
* 2: elPrep performs garbage collection interleaved with the execution at regular intervals. This may reduce peak memory significantly, but may slow down execution considerably.

If the option is not passed explicitly, elPrep chooses default settings for gc-on. If elprep is asked to perform duplicate marking (--mark-duplicates) or sorting (--sorting-order or --replace-reference-sequences that triggers a sort), it assumes --gc-on 0 is intented. Otherwise the default is --gc-on 2.

### --split-file

This option tells elPrep it is sorting a file created by the elprep split command. Using this option allows elprep to make use of meta information generated by the elprep split command for optimal parallel sorting. 

## Split and Merge tools

The elprep split command can be used to split up .sam files into smaller files that store the reads "per chromosome". elPrep determines the "chromosomes" by analyzing the sequence dictionary in the header of the input file and generates a split file for each chromosome that stores all read pairs that map to that chromosome. Additonally, elPrep creates a file for storing the unmapped reads, as well as a file for storing the pairs where reads map to different chromosomes. elPrep also duplicates the latter pairs across chromosome files so that preparation pipelines have access to all information they need to run correctly. Once processed, use the elprep merge command to merge the split files back into a single .sam file.

Splitting the .sam file into smaller files for processing "per chromosome" is useful for reducing the memory pressure as these split files are typically significantly smaller than the  input file as a whole. Splitting also makes it possible to parallelize the processing of a single .sam file by distributing the different split files across different processing nodes.

We provide an example python script "elprep-sfm.py" that illustrates how to use the split and merge commands to process a .sam file. We also provide an example python scrip "elprep-sfm-gnupar.py" that uses GNU parallel for optimal execution on multi-socket servers. For a detailed description, see below.

## Name

### elprep split - a commandline tool for splitting .sam/.bam/.cram files per chromosome so they can be processed without information loss

## Synopsis

	elprep split [sam-file | /path/to/input/] /path/to/output/ --output-prefix "split-sam-file" --output-type sam --nr-of-threads $threads

## Description

The elprep split command requires two arguments: 1) the input file or a path to multiple input files and 2) a path to a directory where elPrep can store the split files. The input file(s) can be .sam, .bam, or .cram. It is also possible to use /dev/stdin as the input for using Unix pipes. There are no structural requirements on the input file(s) for using elprep split. For example, it is not necessary to sort the input file, nor is it necessary to convert to .bam or index the input file.

elPrep creates the output directory denoted by the output path, unless the directory already exists, in which case elPrep may override the existing files in that directory. Please make sure elPrep has the correct permissions for writing that directory.

The split command outputs two types of files:

1. a subdirectory "/path/to/output/splits/" containing a file per entry in the sequence dictionary of the input file that contains all reads that map to that entry.
2. a "/path/to/output/output-prefix-spread.output-type" file containing all reads of which the mate maps to a different entry in the sequence dictionary of the input file.

To process the files created by the elprep split command, one needs to call the elprep command for each entry in the path/to/output/splits/ directory as well as the /path/to/output/output-prefix-spread.output-type file. To process the files in the "split" directory one needs to pass the --split-file option to the elprep command. No special options are necessary for processing the reads in the spread reads file (/path/to/output/output-prefix-spread.output-type). Finally, the output files produced this way, need to be merged with the elprep merge command. For example scripts see "elprep-sfm.py" and "elprep-sfm-gnupar.py".

## Options

### --output-prefix name

The names of the split files created by elprep split are generated by combing a prefix and a chromosome name. The --output-prefix option sets that prefix. For example, if the prefix is "NA12878", and sequence dictionary of the input file contains the chromosomes "chr1", "chr2", and "chr3", and so on, then the names of the split files will be "NA12878-chr1.sam", "NA12878-chr2.sam", "NA12878-chr3.sam", and so on.

If the user does not specify the --output-prefix option, the name of the input file, minus the file extension, is used as a prefix.

### --output-type [sam | bam | cram]

This command option sets the format of the split files. By default, elprep uses the same format as the input file for the split files.

### --nr-of-threads number

This command option sets the number of threads that elPrep uses during execution for converting between .bam/.sam formats. The default number of threads is 1. This option is only useful when the input file is in .bam format or when the --output-type of the split files is chosen to be .bam.

## Name

### elprep merge - a commandline tool for merging .sam/.bam/.cram files created by elprep split

## Synopsis

	elprep merge /path/to/input/ sam-output-file --nr-of-threads $threads

## Description

The elprep merge command requires two arguments: a path to the files that need to be merged, and an output file. Use this command to merge files created with elprep split. The output file can be .sam, .bam, or .cram. It is also possible to use /dev/stdout as output when using Unix pipes for connecting other tools.

## Options

### --nr-of-threads number

This command option sets the number of threads that elPrep uses during execution for converting between .bam/.sam formats. The default number of threads is 1. This option is only useful when the files to be marged are in .bam format or when the output-file is in .bam format.

# Python Scripts

The Python scripts have been tested with Python 2.7.3.

## Name

### elprep-sfm.py - a Python script that illustrates the use of elprep split and merge

## Synopsis

	./elprep-sfm.py $input $output --filter-unmapped-reads strict --replace-reference-sequences ucsc.hg19.dict --replace-read-group "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" --mark-duplicates --sorting-order coordinate --nr-of-threads $threads --intermediate-files-output-type [sam | bam | cram]

## Description

A Python script that combines the elprep split, filter, and merge commands. The script may be used as a drop-in replacement for the elprep filter command, except that:

* It first calls elprep split to split up the input file per chromosome. The script creates a temp folder in the current working directory for storing the split files created this way.
* It then calls the elprep filter command for processing the split files one by one.
* Finally, it calls the elprep merge command to create the output file from the split files.

## Special Options

### --intermediates-files-output-type [sam | bam | cram]

The output type that will be used for the intermediate split files. The default output type for intermediate files is .sam.

### Name

### elprep-sfm-gnupar.py - a Python script that illustrates the use of elprep split and merge and GNU parallel for optimal execution on a multi-socket server

## Synopsis

	./elprep-sfm-gnupar.py $input $output --filter-unmapped-reads strict --replace-reference-sequences ucsc.hg19.dict --replace-read-group "ID:group1 LB:lib1 PL:illumina PU:unit1 SM:sample1" --mark-duplicates --sorting-order coordinate --nr-of-threads $threads --nr-of-jobs $jobs --intermediate-files-output-type [sam | bam | cram]

## Description

A Python script that combines the elprep split, filter, and merge commands with GNU parallel. The script may be used as a drop-in replacement for the elprep filter command, except that it adds a --nr-of-jobs option. This option is used to set the number of elprep calls that can be processed in parallel.

The script is structured as follows:

* It first calls elprep split to split up the input file per chromosome. The script creates a temp folder in the current working directory for storing the split files created this way.
* It then calls the elprep filter command for processing the split files in parallel using GNU parallel. There will be [--nr-of-jobs] split files executed in parallel, which each use [--nr-of-threads] threads for their execution.
* It then calls the elprep filter command for processing the spread reads file created by the elprep split command.
* Finally, the script calls the elprep merge command to combine the results of processing the split files and the spread reads file into a single output.

Using this script is useful when using a multi-socket server. For example, when using a server with two sockets and two 12 core processors, running up to two elprep commands in parallel with each 12 threads is optimal from a performance perspective. The command could thus look as follows:

	./elprep-sfm-gnupar.py input.bam output.bam ... --nr-of-threads 12 --nr-of-jobs 2

## Special Options

### --nr-of-jobs number

The number of elprep commands that will be executed in parallel.

### --intermediates-files-output-type [sam | bam | cram]

The output type that will be used for the intermediate split files. The default output type for intermediate files is .sam.

## Name

### elprep\_io\_wrapper.py - a Python script that illustrates how to use piping with elPrep

A Python script that implements the .bam/.cram conversion through piping to SAMtools. The script is loaded by the elprep-sfm.py and elprep.py scripts. You can edit the SAMtools calls to tweak the parameters to fit your needs, or use this script as a starting point for piping to other .bam/.cram conversion tools than SAMtools.

## Name

### elprep.py - a Python script that wraps the elPrep binary

A Python script that wraps the elPrep binary and loads the elprep\_io\_wrapper.py script for .bam/.cram conversion using SAMtools rather than using the internal piping of elPrep. You can use this script as a drop-in replacement for the elprep binary. 

## Name

### elprep\_entrypoint.py - A Python scripts that can be used as an entrypoint for Docker.

A Python scripts that wraps the elprep.py, elprep-sfm.py and elprep-sfm-gnupar.py scripts. It can for example be used as an entrypoint to Docker.

# Extending elPrep

If you wish to extend elPrep, for example by adding your own filters, please consult our [API documentation](http://exascience.github.io/elprep/elprep-package/index.html). 

<!--### --timed

When this option is passed, elPrep times and prints the execution spent per phase --filtering, sorting, I/O-- to the standard output.
-->
