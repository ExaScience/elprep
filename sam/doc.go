// Package sam is a library for parsing and representing SAM files,
// and for efficiently executing sequencing pipelines on
// .sam/.bam/.cram files, taking advantage of modern multi-core
// processors.
//
// Modifications to headers and alignments are expressed as
// filters. The library comes with a number of commonly used
// pre-defined filters, but you can also define and use your own
// filters. A pipeline can be executed with the RunPipeline method of
// the PipelineInput interface, which accepts SAM/BAM/CRAM files as
// input and/or output sources, but can also operate on an in-memory
// representation of such files. PipelineInput and PipelineOutput can
// be implemented to also operate on other input/output sources, such
// as databases.
//
// elPrep provides high-level Filter and AlignmentFilter types that
// operate on SAM file header and alignment structs. elPrep then uses
// the pargo library for expressing pipelines of such filters for
// efficient parallel execution. It is normally not necessary to deal
// with pargo pipelines directly, but you can check the documentation
// at https://godoc.org/github.com/ExaScience/pargo/pipeline for
// details of pargo pipelines if necessary.
package sam
