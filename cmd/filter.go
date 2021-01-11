// elPrep: a high-performance tool for analyzing SAM/BAM files.
// Copyright (c) 2017-2020 imec vzw.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version, and Additional Terms
// (see below).

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Affero General Public License for more details.

// You should have received a copy of the GNU Affero General Public
// License and Additional Terms along with this program. If not, see
// <https://github.com/ExaScience/elprep/blob/master/LICENSE.txt>.

package cmd

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"strings"

	"github.com/exascience/elprep/v5/fasta"

	"github.com/exascience/elprep/v5/internal"

	"github.com/exascience/elprep/v5/bed"
	"github.com/exascience/elprep/v5/filters"
	"github.com/exascience/elprep/v5/sam"
	"github.com/exascience/elprep/v5/utils"
	psync "github.com/exascience/pargo/sync"
)

func parseAndMergeSpreadFile(spreadFile string, reads *sam.Sam) {
	file := sam.Open(spreadFile)
	defer file.Close()
	spreadReads := sam.NewSam()
	file.RunPipeline(spreadReads, []sam.Filter{func(_ *sam.Header) sam.AlignmentFilter {
		contigs, ok := reads.Header.Contigs()
		if !ok {
			log.Panic("Cannot call haplotypes on split file without contig information.")
		}
		contigCheck := make(map[string]bool)
		for _, contig := range contigs {
			contigCheck[contig] = true
		}
		return func(aln *sam.Alignment) bool {
			return contigCheck[aln.RNAME]
		}
	},
		filters.AddREFID,
	}, sam.Keep)
	reads.Alignments = sam.By(sam.CoordinateLess).ParallelStableMerge(reads.Alignments, spreadReads.Alignments)
}

func runRemainingPipeline(
	phase int64,
	fileOut, outFormat string, sortingOrder sam.SortingOrder,
	filteredReads *sam.Sam, filters2 []sam.Filter,
	hc *filters.HaplotypeCaller, sampleName string, vcfOutput string, bedRegions bed.Bed, spreadFile string,
	timed bool, profile string) {
	switch fileOut {
	case "/dev/null", "/dev/zero":
		// hc must be non-nil
		phase++
		timedRun(timed, profile, "Preparing variant calling.", phase, func() {
			if spreadFile != "" {
				parseAndMergeSpreadFile(spreadFile, filteredReads)
			}
			filters2 = append(filters2, filters.FilterReadsBySampleName(&sampleName), filters.HaplotypeCallAln)
			filteredReads.RunPipeline(filteredReads, filters2, sortingOrder)
		})
		filters2 = nil
		go runtime.GC()
		phase++
		timedRun(timed, profile, "Calling variants.", phase, func() {
			vcfPathname := internal.FilepathAbs(vcfOutput)
			internal.MkdirAll(filepath.Dir(vcfPathname), 0700)
			hc.CallVariants(filteredReads, sampleName, vcfPathname, bedRegions)
		})
	default:
		if hc == nil {
			// Write output to file
			phase++
			timedRun(timed, profile, "Write to file.", phase, func() {
				pathname := internal.FilepathAbs(fileOut)
				internal.MkdirAll(filepath.Dir(pathname), 0700)
				output := sam.Create(pathname, outFormat)
				defer output.Close()
				filteredReads.RunPipeline(output, filters2, sortingOrder)
			})
		} else {
			if len(filters2) > 0 {
				phase++
				timedRun(timed, profile, "Preparing output of final SAM file.", phase, func() {
					filteredReads.RunPipeline(filteredReads, filters2, sortingOrder)
				})
				filters2 = nil
				go runtime.GC()
			}
			if spreadFile != "" {
				phase++
				timedRun(timed, profile, "Merging in spread reads.", phase, func() {
					parseAndMergeSpreadFile(spreadFile, filteredReads)
				})
				go runtime.GC()
			}
			keptReads := new(sam.Sam)
			*keptReads = *filteredReads
			phase++
			timedRun(timed, profile, "Write to file.", phase, func() {
				pathname := internal.FilepathAbs(fileOut)
				internal.MkdirAll(filepath.Dir(pathname), 0700)
				output := sam.Create(pathname, outFormat)
				defer output.Close()
				filteredReads.RunPipeline(output, nil, sortingOrder)
			})
			go runtime.GC()
			phase++
			timedRun(timed, profile, "Preparing variant calling.", phase, func() {
				filters3 := []sam.Filter{filters.FilterReadsBySampleName(&sampleName), filters.HaplotypeCallAln}
				keptReads.RunPipeline(keptReads, filters3, sortingOrder)
			})
			go runtime.GC()
			phase++
			timedRun(timed, profile, "Calling variants.", phase, func() {
				vcfPathname := internal.FilepathAbs(vcfOutput)
				internal.MkdirAll(filepath.Dir(vcfPathname), 0700)
				hc.CallVariants(keptReads, sampleName, vcfPathname, bedRegions)
			})
		}
	}
}

// Run the best practices pipeline. Version that uses an intermediate
// slice so that sorting, mark-duplicates, BSQR, and HaplotypeCaller are supported.
func runBestPracticesPipelineIntermediateSam(
	fileIn, fileOut, outFormat string, sortingOrder sam.SortingOrder,
	filters1, filters2 []sam.Filter,
	opticalDuplicateFilter func(reads *sam.Sam),
	baseRecalibrator *filters.BaseRecalibrator, quantizeLevels int, sqqList []uint8, maxCycle int, recalFile string,
	hc *filters.HaplotypeCaller, sampleName, vcfOutput string, bedRegions bed.Bed, spreadFile string,
	timed bool, profile string) {

	// Input and first filters and sorting
	filteredReads := sam.NewSam()
	phase := int64(1)
	timedRun(timed, profile, "Reading SAM into memory and applying filters.", phase, func() {
		input := sam.Open(internal.FilepathAbs(fileIn))
		defer input.Close()
		input.RunPipeline(filteredReads, filters1, sortingOrder)
	})
	filters1 = nil
	go runtime.GC()
	if opticalDuplicateFilter != nil {
		phase++
		timedRun(timed, profile, "Marking optical duplicates.", phase, func() {
			opticalDuplicateFilter(filteredReads)
		})
		opticalDuplicateFilter = nil
		go runtime.GC()
	}
	if sortingOrder != sam.Unsorted {
		sortingOrder = sam.Keep
	}
	// Base recalibration
	if baseRecalibrator != nil {
		// filter out reads not in target regions if bedRegions is passed
		if bedRegions != nil {
			phase++
			targetRegionsFilter := filters.RemoveNonOverlappingReads(bedRegions)
			timedRun(timed, profile, "Filtering out reads not in target regions.", phase, func() {
				filteredReads.RunPipeline(filteredReads, []sam.Filter{targetRegionsFilter}, sortingOrder)
			})
			targetRegionsFilter = nil
			go runtime.GC()
		}
		var baseRecalibratorTables *filters.BaseRecalibratorTables
		phase++
		timedRun(timed, profile, "Base recalibration.", phase, func() {
			baseRecalibratorTables = baseRecalibrator.Recalibrate(filteredReads, maxCycle)
		})
		baseRecalibrator = nil
		go runtime.GC()
		// Finalize BQSR tables + log recal file
		phase++
		timedRun(timed, profile, "Finalize BQSR tables.", phase, func() {
			baseRecalibratorTables.FinalizeBQSRTables()
			pathname := internal.FilepathAbs(recalFile)
			internal.MkdirAll(filepath.Dir(pathname), 0700)
			baseRecalibratorTables.PrintBQSRTables(pathname)
		})
		phase++
		timedRun(timed, profile, "Apply BQSR.", phase, func() {
			n := len(filteredReads.Alignments) / 4096
			if m := 2 * runtime.GOMAXPROCS(0); m > n {
				n = m
			}
			filteredReads.NofBatches(n)
			filteredReads.RunPipeline(filteredReads, []sam.Filter{baseRecalibratorTables.ApplyBQSR(quantizeLevels, sqqList, maxCycle)}, sortingOrder)
		})
		baseRecalibratorTables = nil
		go runtime.GC()
	}
	runRemainingPipeline(
		phase, fileOut, outFormat, sortingOrder,
		filteredReads, filters2,
		hc, sampleName, vcfOutput, bedRegions, spreadFile,
		timed, profile)
}

func runBestPracticesPipelineWithBQSRApplyOnly(
	input, output, outFormat string, sortingOrder sam.SortingOrder,
	filters []sam.Filter, baseRecalibratorTables *filters.BaseRecalibratorTables, recalFile string,
	timed bool, profile string) {
	// Finalize BQSR tables + log recal file
	timedRun(timed, profile, "Finalize BQSR tables", 1, func() {
		baseRecalibratorTables.FinalizeBQSRTables()
		pathname := internal.FilepathAbs(recalFile)
		internal.MkdirAll(filepath.Dir(pathname), 0700)
		baseRecalibratorTables.PrintBQSRTables(pathname)
	})
	baseRecalibratorTables = nil
	go runtime.GC()
	timedRun(timed, profile, "Running pipeline.", 2, func() {
		input := sam.Open(internal.FilepathAbs(input))
		defer input.Close()
		pathname := internal.FilepathAbs(output)
		internal.MkdirAll(filepath.Dir(pathname), 0700)
		output := sam.Create(pathname, outFormat)
		defer output.Close()
		input.RunPipeline(output, filters, sortingOrder)
	})
}

func runBestPracticesPipelineIntermediateSamWithBQSRApplyOnly(
	input, output, outFormat string, sortingOrder sam.SortingOrder,
	filters []sam.Filter, baseRecalibratorTables *filters.BaseRecalibratorTables, recalFile string,
	hc *filters.HaplotypeCaller, sampleName, vcfOutput string, bedRegions bed.Bed, spreadFile string,
	timed bool, profile string) {
	// Finalize BQSR tables + log recal file
	timedRun(timed, profile, "Finalize BQSR tables", 1, func() {
		baseRecalibratorTables.FinalizeBQSRTables()
		pathname := internal.FilepathAbs(recalFile)
		internal.MkdirAll(filepath.Dir(pathname), 0700)
		baseRecalibratorTables.PrintBQSRTables(pathname)
	})
	baseRecalibratorTables = nil
	go runtime.GC()
	filteredReads := sam.NewSam()
	timedRun(timed, profile, "Running pipeline.", 2, func() {
		input := sam.Open(internal.FilepathAbs(input))
		defer input.Close()
		input.RunPipeline(filteredReads, filters, sortingOrder)
	})
	filters = nil
	go runtime.GC()
	runRemainingPipeline(
		3, output, outFormat, sortingOrder,
		filteredReads, nil,
		hc, sampleName, vcfOutput, bedRegions, spreadFile,
		timed, profile)
}

func runBestPracticesPipelineIntermediateSamWithBQSRCalculateTablesOnly(
	fileIn, fileOut, outFormat string, sortingOrder sam.SortingOrder,
	filters1, filters2 []sam.Filter,
	opticalDuplicateFilter func(reads *sam.Sam),
	baseRecalibrator *filters.BaseRecalibrator, tableFile string, maxCycle int,
	bedRegions bed.Bed, timed bool, profile string) {
	// Input and first filters and sorting
	filteredReads := sam.NewSam()
	phase := int64(1)
	timedRun(timed, profile, "Reading SAM into memory and applying filters.", phase, func() {
		input := sam.Open(internal.FilepathAbs(fileIn))
		defer input.Close()
		input.RunPipeline(filteredReads, filters1, sortingOrder)
	})
	filters1 = nil
	go runtime.GC()
	// Marking optical duplicates
	if opticalDuplicateFilter != nil {
		phase++
		timedRun(timed, profile, "Marking optical duplicates.", phase, func() {
			opticalDuplicateFilter(filteredReads)
		})
		opticalDuplicateFilter = nil
		go runtime.GC()
	}
	// Filter out reads based on target regions
	if bedRegions != nil {
		phase++
		targetRegionsFilter := filters.RemoveNonOverlappingReads(bedRegions)
		timedRun(timed, profile, "Filtering out reads not in target regions.", phase, func() {
			filteredReads.RunPipeline(filteredReads, []sam.Filter{targetRegionsFilter}, sortingOrder)
		})
		targetRegionsFilter = nil
		go runtime.GC()
	}
	// Base recalibration
	var baseRecalibratorTables *filters.BaseRecalibratorTables
	phase++
	timedRun(timed, profile, "Base recalibration", phase, func() {
		baseRecalibratorTables = baseRecalibrator.Recalibrate(filteredReads, maxCycle)
	})
	// Print bqsr table to file
	phase++
	timedRun(timed, profile, "Print intermediate BQSR tables", phase, func() {
		pathname := internal.FilepathAbs(tableFile)
		internal.MkdirAll(filepath.Dir(pathname), 0700)
		baseRecalibratorTables.PrintBQSRTablesToIntermediateFile(pathname)
	})
	baseRecalibratorTables = nil
	go runtime.GC()
	// Write output to file
	phase++
	timedRun(timed, profile, "Write to file.", phase, func() {
		pathname := internal.FilepathAbs(fileOut)
		internal.MkdirAll(filepath.Dir(pathname), 0700)
		output := sam.Create(pathname, outFormat)
		defer output.Close()
		if sortingOrder != sam.Unsorted {
			sortingOrder = sam.Keep
		}
		filteredReads.RunPipeline(output, filters2, sortingOrder)
	})
}

// Run the best practices pipeline. Version that doesn't use an
// intermediate slice when neither sorting nor mark-duplicates are
// needed.
func runBestPracticesPipeline(
	fileIn, fileOut, outFormat string, sortingOrder sam.SortingOrder,
	filters []sam.Filter,
	timed bool, profile string) {
	timedRun(timed, profile, "Running pipeline.", 1, func() {
		input := sam.Open(internal.FilepathAbs(fileIn))
		defer input.Close()
		pathname := internal.FilepathAbs(fileOut)
		internal.MkdirAll(filepath.Dir(pathname), 0700)
		output := sam.Create(pathname, outFormat)
		defer output.Close()
		input.RunPipeline(output, filters, sortingOrder)
	})
}

// FilterHelp is the help string for this command.
const FilterHelp = "\nfilter parameters:\n" +
	"elprep filter sam-file sam-output-file\n" +
	"[--output-type [sam | bam]]\n" +
	"[--replace-reference-sequences sam-file]\n" +
	"[--filter-unmapped-reads]\n" +
	"[--filter-unmapped-reads-strict]\n" +
	"[--filter-mapping-quality mapping-quality]\n" +
	"[--filter-non-exact-mapping-reads]\n" +
	"[--filter-non-exact-mapping-reads-strict]\n" +
	"[--filter-non-overlapping-reads bed-file]\n" +
	"[--replace-read-group read-group-string]\n" +
	"[--mark-duplicates]\n" +
	"[--mark-optical-duplicates file]\n" +
	"[--optical-duplicates-pixel-distance nr]\n" +
	"[--remove-duplicates]\n" +
	"[--remove-optional-fields [all | list]]\n" +
	"[--keep-optional-fields [none | list]]\n" +
	"[--sorting-order [keep | unknown | unsorted | queryname | coordinate]]\n" +
	"[--clean-sam]\n" +
	"[--reference elfasta]\n" +
	"[--bqsr recal-file]\n" +
	"[--quantize-levels nr]\n" +
	"[--sqq list]\n" +
	"[--max-cycle nr]\n" +
	"[--known-sites list]\n" +
	"[--haplotypecaller vcf-file]\n" +
	"[--reference-confidence [GVCF | BP_RESOLUTION | NONE]\n" +
	"[--sample-name sample-name]\n" +
	"[--activity-profile igv-file]\n" +
	"[--assembly-regions igv-file]\n" +
	"[--assembly-region-padding nr]\n" +
	"[--target-regions]\n" +
	"[--nr-of-threads nr]\n" +
	"[--timed]\n" +
	"[--log-path path]\n"

// FilterExtendedHelp is the extended help string for this command.
const FilterExtendedHelp = FilterHelp +
	"[--mark-optical-duplicates-intermediate file]\n" +
	"[--bqsr-tables-only table-file]\n" +
	"[--bqsr-apply path]\n" +
	"[--recal-file file]\n" +
	"[--spread-file file]\n" +
	"[--pg-cmd-line cmd]\n"

// Filter implements the elprep filter command.
func Filter() {
	var (
		outputType                                               string
		replaceReferenceSequences                                string
		filterUnmappedReads, filterUnmappedReadsStrict           bool
		filterMappingQuality                                     int
		filterNonExactMappingReads                               bool
		filterNonExactMappingReadsStrict                         bool
		filterNonOverlappingReads                                string
		replaceReadGroup                                         string
		markDuplicates, removeDuplicates                         bool
		markOpticalDuplicates, markOpticalDuplicatesIntermediate string
		opticalDuplicatesPixelDistance                           int
		removeOptionalFields                                     string
		keepOptionalFields                                       string
		sortingOrderString                                       string
		cleanSam                                                 bool
		bqsr                                                     string
		bqsrTablesOnly                                           string
		bqsrApplyFromTables                                      string
		quantizeLevels                                           int
		sqq                                                      string
		maxCycle                                                 int
		knownSites                                               string
		recalFile                                                string
		nrOfThreads                                              int
		timed                                                    bool
		profile                                                  string
		logPath                                                  string
		renameChromosomes                                        bool
		referenceElFasta                                         string
		vcfOutput                                                string
		assemblyRegionPadding                                    int
		referenceConfidence                                      string
		sampleName                                               string
		activityProfile, assemblyRegions                         string
		targetRegions                                            string
		spreadFile                                               string
		pgCmdLine                                                string
		randomSeedFile                                           string
		clearDuplicateFlag                                       bool
	)

	var flags flag.FlagSet

	flags.StringVar(&outputType, "output-type", "", "format of the output file")
	flags.StringVar(&replaceReferenceSequences, "replace-reference-sequences", "", "replace the existing header by a new one")
	flags.BoolVar(&filterUnmappedReads, "filter-unmapped-reads", false, "remove all unmapped alignments")
	flags.BoolVar(&filterUnmappedReadsStrict, "filter-unmapped-reads-strict", false, "remove all unmapped alignments, taking also POS and RNAME into account")
	flags.IntVar(&filterMappingQuality, "filter-mapping-quality", 0, "output only reads that equal or exceed given mapping quality")
	flags.BoolVar(&filterNonExactMappingReads, "filter-non-exact-mapping-reads", false, "output only exact mapping reads (soft-clipping allowed) based on cigar string (only M,S allowed)")
	flags.BoolVar(&filterNonExactMappingReadsStrict, "filter-non-exact-mapping-reads-strict", false, "output only exact mapping reads (soft-clipping allowed) based on optional fields X0=1, X1=0, XM=0, XO=0, XG=0")
	flags.StringVar(&filterNonOverlappingReads, "filter-non-overlapping-reads", "", "output only reads that overlap with the given regions (bed format)")
	flags.StringVar(&replaceReadGroup, "replace-read-group", "", "add or replace alignment read groups")
	flags.BoolVar(&markDuplicates, "mark-duplicates", false, "mark duplicates")
	flags.StringVar(&markOpticalDuplicates, "mark-optical-duplicates", "", "mark optical duplicates")
	flags.StringVar(&markOpticalDuplicatesIntermediate, "mark-optical-duplicates-intermediate", "", "mark optical duplicates intermediate file (only for split files)")
	flags.IntVar(&opticalDuplicatesPixelDistance, "optical-duplicates-pixel-distance", 100, "pixel distance used for optical duplicate marking")
	flags.BoolVar(&removeDuplicates, "remove-duplicates", false, "remove duplicates")
	flags.StringVar(&removeOptionalFields, "remove-optional-fields", "", "remove the given optional fields")
	flags.StringVar(&keepOptionalFields, "keep-optional-fields", "", "remove all except for the given optional fields")
	flags.StringVar(&sortingOrderString, "sorting-order", string(sam.Keep), "determine output order of alignments, one of keep, unknown, unsorted, queryname, or coordinate")
	flags.BoolVar(&cleanSam, "clean-sam", false, "clean the sam file")
	flags.StringVar(&bqsr, "bqsr", "", "base quality score recalibration")
	flags.StringVar(&bqsrTablesOnly, "bqsr-tables-only", "", "base quality score recalibration table calculation (only with split/merge)")
	flags.StringVar(&bqsrApplyFromTables, "bqsr-apply", "", "base quality score recalibration application (only with split/merge)")
	flags.IntVar(&quantizeLevels, "quantize-levels", 0, "number of levels to be used for quantizing recalibrated base qualities (only with --bqsr)")
	flags.StringVar(&sqq, "sqq", "", "levels to be used for statically quantizing recalibrated base qualities (only with --bqsr)")
	flags.IntVar(&maxCycle, "max-cycle", 500, "maximum cycle length (only with --bqsr)")
	flags.StringVar(&knownSites, "known-sites", "", "list of vcf files containing known sites for base recalibration (only with --bqsr)")
	flags.StringVar(&recalFile, "recal-file", "", "log file for recalibration tables (only with --bqsr-apply)")
	flags.IntVar(&nrOfThreads, "nr-of-threads", 0, "number of worker threads")
	flags.BoolVar(&timed, "timed", false, "measure the runtime")
	flags.StringVar(&profile, "profile", "", "write a runtime profile to the specified file(s)")
	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")
	flags.BoolVar(&renameChromosomes, "rename-chromosomes", false, "")
	flags.StringVar(&referenceElFasta, "reference", "", "reference used for base quality score recalibration and haplotypecaller (elfasta format)")
	flags.StringVar(&vcfOutput, "haplotypecaller", "", "invoke the haplotypecaller")
	flags.IntVar(&assemblyRegionPadding, "assembly-region-padding", 100, "padding around assembly regions during variant calling (only with --haplotypecaller)")
	flags.StringVar(&referenceConfidence, "reference-confidence", "GVCF", "reference confidence mode (GVCF, BP_RESOLUTION, NONE) (only with --haplotypecaller)")
	flags.StringVar(&sampleName, "sample-name", "", "sample name to use in haplotypecaller (only with --haplotypecaller)")
	flags.StringVar(&activityProfile, "activity-profile", "", "IGV file to write the activity profile to (only with --haplotypecaller)")
	flags.StringVar(&assemblyRegions, "assembly-regions", "", "IGV file to write the assembly regions to (only with --haplotypecaller)")
	flags.StringVar(&targetRegions, "target-regions", "", "bed file specifying which regions in the genome to look at (only with --bqsr and --haplotypecaller)")
	flags.StringVar(&spreadFile, "spread-file", "", "spread file for calling haplotype caller on a split file (only with --haplotypecaller)")
	flags.StringVar(&pgCmdLine, "pg-cmd-line", "", "program command line to be stored in the header (only for sfm subcommands)")
	if internal.PedanticMode {
		flags.StringVar(&randomSeedFile, "random-seed-file", "", "random seed file (only for sfm subcommands in pedantic mode")
	}
	flags.BoolVar(&clearDuplicateFlag, "clear-duplicate-flag", false, "clear the duplicate flag in every read")

	flags.StringVar(&internal.BQSRTablenamePrefix, "bqsr-tablename-prefix", "GATK", "prefix to be used in BQSR table reports")

	parseFlags(flags, 4, FilterHelp)

	input := getFilename(os.Args[2], FilterHelp)
	output := getFilename(os.Args[3], FilterHelp)

	setLogOutput(logPath)

	// sanity checks

	var sanityChecksFailed bool

	if !checkOutputFormat(outputType) {
		sanityChecksFailed = true
	}

	if !checkExist("", input) {
		sanityChecksFailed = true
	}
	if !checkCreate("", output) {
		sanityChecksFailed = true
	}

	if replaceReferenceSequences != "" && !checkExist("--replace-reference-sequences", replaceReferenceSequences) {
		sanityChecksFailed = true
	}
	if filterNonOverlappingReads != "" && !checkExist("--filter-non-overlapping-reads", filterNonOverlappingReads) {
		sanityChecksFailed = true
	}
	if markOpticalDuplicates != "" && !checkCreate("--mark-optical-duplicates", markOpticalDuplicates) {
		sanityChecksFailed = true
	}
	if markOpticalDuplicatesIntermediate != "" && !checkCreate("--mark-optical-duplicates-intermediate", markOpticalDuplicatesIntermediate) {
		sanityChecksFailed = true
	}
	if bqsr != "" && !checkCreate("--bqsr", bqsr) {
		sanityChecksFailed = true
	}
	if bqsrTablesOnly != "" && !checkCreate("--bqsr-tables-only", bqsrTablesOnly) {
		sanityChecksFailed = true
	}
	if bqsrApplyFromTables != "" && !checkExist("--bqsr-apply", bqsrApplyFromTables) {
		sanityChecksFailed = true
	}
	if referenceElFasta != "" && !checkExist("--reference", referenceElFasta) {
		sanityChecksFailed = true
	}

	var knownSitesList []string

	if knownSites != "" {
		knownSitesList = strings.Split(knownSites, ",")
		for _, k := range knownSitesList {
			if !checkExist("--known-sites", k) {
				sanityChecksFailed = true
			}
		}
	}

	if recalFile != "" && !checkCreate("--recal-file", recalFile) {
		sanityChecksFailed = true
	}

	if profile != "" && !checkCreate("--profile", profile) {
		sanityChecksFailed = true
	}

	sortingOrder := sam.SortingOrder(sortingOrderString)

	switch sortingOrder {
	case sam.Keep, sam.Unknown, sam.Unsorted, sam.Queryname, sam.Coordinate:
	default:
		sanityChecksFailed = true
		log.Println("Error: Invalid sorting-order: ", sortingOrder)
	}

	if (replaceReferenceSequences != "") && (sortingOrder == sam.Keep) {
		log.Println("Warning: Requesting to keep the order of the input file while replacing the reference sequence dictionary may force an additional sorting phase to ensure the original sorting order is respected.")
	}

	if keepOptionalFields != "" && removeOptionalFields != "" {
		sanityChecksFailed = true
		log.Println("Error: Cannot use --keep-optional-fields and --remove-optional-fields in the same filter command.")
	}

	if nrOfThreads < 0 {
		sanityChecksFailed = true
		log.Println("Error: Invalid nr-of-threads: ", nrOfThreads)
	}

	if markOpticalDuplicates != "" && !markDuplicates {
		sanityChecksFailed = true
		log.Println("Error: Cannot use --mark-optical-duplicates without also using --mark-duplicates.")
	}

	if markOpticalDuplicatesIntermediate != "" && !markDuplicates {
		sanityChecksFailed = true
		log.Println("Error: Cannot use --mark-optical-duplicates-intermediate without also using --mark-duplicates.")
	}

	if opticalDuplicatesPixelDistance != 100 && markOpticalDuplicates == "" && markOpticalDuplicatesIntermediate == "" {
		sanityChecksFailed = true
		log.Println("Error: cannot use --optical-duplicates-pixel-distance without also using --mark-optical-duplicates or --mark-optical-duplicates-intermediate")
	}

	if bqsr != "" {
		if bqsrTablesOnly != "" {
			sanityChecksFailed = true
			log.Println("Error: cannot use both --bsqr and --bqsr-tables-only in the same command")
		}
		if bqsrApplyFromTables != "" {
			sanityChecksFailed = true
			log.Println("Error: cannot use both --bqsr and --bqsr-apply in the same command")
		}
		if !checkBQSROptions(referenceElFasta, bqsr, recalFile) {
			sanityChecksFailed = true
		}
	} else if bqsrTablesOnly != "" {
		if bqsrApplyFromTables != "" {
			sanityChecksFailed = true
			log.Println("Error: cannot use both --bqsr-tables-only and --bqsr-apply in the same command")
		}
		if !checkBQSRTablesOnlyOptions(bqsrTablesOnly, referenceElFasta) {
			sanityChecksFailed = true
		}
	} else if bqsrApplyFromTables != "" {
		if !checkBQSRApplyOptions(bqsrApplyFromTables, recalFile) {
			sanityChecksFailed = true
		}
	} else if !checkNonBQSROptions(quantizeLevels, sqq, knownSites) {
		sanityChecksFailed = true
	}

	if vcfOutput == "" {
		switch output {
		case "/dev/null", "/dev/zero":
			sanityChecksFailed = true
			log.Println("Error: Neither writing a SAM output, nor a VCF output.")
		}
	} else if !checkCreate("--haplotypecaller", vcfOutput) || !checkHaplotypecallerOptions(referenceElFasta) {
		sanityChecksFailed = true
	}

	if activityProfile != "" && !checkCreate("--activity-profile", activityProfile) {
		sanityChecksFailed = true
	}

	if assemblyRegions != "" && !checkCreate("--assembly-regions", assemblyRegions) {
		sanityChecksFailed = true
	}

	referenceConfidence = strings.ToUpper(referenceConfidence)
	switch referenceConfidence {
	case "GVCF", "BP_RESOLUTION", "NONE":
	default:
		sanityChecksFailed = true
		log.Println("Error: --reference-confidence must be one of GVCF, BP_RESOLUTION, or NONE")
	}

	if assemblyRegionPadding < 0 {
		sanityChecksFailed = true
		log.Println("Error: --assembly-region-padding must be >= 0")
	} else if assemblyRegionPadding > math.MaxInt32 {
		sanityChecksFailed = true
		log.Println("Error: --assembly-region-padding to large")
	}

	if targetRegions != "" {
		if bqsr == "" && vcfOutput == "" && bqsrTablesOnly == "" {
			sanityChecksFailed = true
			log.Println("Error: --target-regions used without --bqsr or --bqsr-tables-only or --haplotypecaller")
		}
	}

	if targetRegions != "" && filterNonOverlappingReads != "" {
		if targetRegions != filterNonOverlappingReads {
			log.Println("Warning: using --filter-non-overlapping-reads and --target-regions with different .bed files.")
		}
	}

	if filterNonOverlappingReads != "" && vcfOutput == "" && targetRegions == "" {
		log.Println("Warning: using --filter-non-overlapping-reads and --haplotypecaller, but not --target-regions. Haplotypecaller will run for full genome.")
	}

	if spreadFile != "" {
		if vcfOutput == "" {
			sanityChecksFailed = true
			log.Println("Error: Cannot use --spread-file without --haplotypecaller")
		}
		if !checkExist("--spread-file", spreadFile) {
			sanityChecksFailed = true
		}
	}

	if randomSeedFile != "" && !checkExist("--random-seed-file", randomSeedFile) {
		sanityChecksFailed = true
	}

	if sanityChecksFailed {
		fmt.Fprint(os.Stderr, FilterHelp)
		os.Exit(1)
	}

	// building filters and output command line

	var command bytes.Buffer
	fmt.Fprint(&command, os.Args[0], " filter ", input, " ", output)

	if outputType != "" {
		fmt.Fprint(&command, " --output-type ", outputType)
	}

	var filters1, filters2 []sam.Filter

	if filterUnmappedReadsStrict {
		filters1 = append(filters1, filters.RemoveUnmappedReadsStrict)
		fmt.Fprint(&command, " --filter-unmapped-reads-strict")
	} else if filterUnmappedReads {
		filters1 = append(filters1, filters.RemoveUnmappedReads)
		fmt.Fprint(&command, " --filter-unmapped-reads")
	}

	if filterMappingQuality > 0 {
		filterMappingQualityFilter := filters.RemoveMappingQualityLessThan(filterMappingQuality)
		filters1 = append(filters1, filterMappingQualityFilter)
		fmt.Fprint(&command, " --filter-mapping-quality ", filterMappingQuality)
	}

	if filterNonExactMappingReads {
		filters1 = append(filters1, filters.RemoveNonExactMappingReads)
		fmt.Fprint(&command, " --filter-non-exact-mapping-reads")
	}

	if filterNonExactMappingReadsStrict {
		filters1 = append(filters1, filters.RemoveNonExactMappingReadsStrict)
		fmt.Fprint(&command, " --filter-non-exact-mapping-reads-strict")
	}

	var filterNonOverlappingReadsBed bed.Bed

	if filterNonOverlappingReads != "" {
		filterNonOverlappingReadsBed = bed.ParseBed(filterNonOverlappingReads)
		filterNonOverlappingReadsFilter := filters.RemoveNonOverlappingReads(filterNonOverlappingReadsBed)
		filters1 = append(filters1, filterNonOverlappingReadsFilter)
		fmt.Fprint(&command, " --filter-non-overlapping-reads ", filterNonOverlappingReads)
	}

	if clearDuplicateFlag {
		filters1 = append(filters1, filters.ClearDuplicateFlag)
	}

	var targetRegionsBed bed.Bed

	if targetRegions != "" {
		fmt.Fprint(&command, " --target-regions ", targetRegions)
		// parse the bed, used for haplotype-caller and/or bqsr
		targetRegionsBed = bed.ParseBed(targetRegions)
		// filter only executed after mark optical duplicates
	}

	if renameChromosomes {
		filters1 = append(filters1, filters.RenameChromosomes)
		fmt.Fprint(&command, " --rename-chromosomes")
	}

	if cleanSam {
		filters1 = append(filters1, filters.CleanSam)
		fmt.Fprint(&command, " --clean-sam")
	}

	if replaceReferenceSequences != "" {
		filters1 = append(filters1, filters.ReplaceReferenceSequenceDictionaryFromSamFile(replaceReferenceSequences))
		fmt.Fprint(&command, " --replace-reference-sequences ", replaceReferenceSequences)
	}

	if replaceReadGroup != "" {
		filters1 = append(filters1, filters.AddOrReplaceReadGroup(sam.ParseHeaderLineFromString(replaceReadGroup)))
		fmt.Fprint(&command, " --replace-read-group \"", replaceReadGroup, "\"")
	}

	if replaceReferenceSequences != "" || markDuplicates ||
		sortingOrder == sam.Coordinate || sortingOrder == sam.Queryname ||
		vcfOutput != "" {
		filters1 = append(filters1, filters.AddREFID)
	}

	var pairs *psync.Map
	var markDupsFilter sam.Filter

	if markDuplicates {
		alsoOpticals := markOpticalDuplicates != "" || markOpticalDuplicatesIntermediate != ""
		markDupsFilter, _, pairs = filters.MarkDuplicates(alsoOpticals)
		filters1 = append(filters1, markDupsFilter)
		fmt.Fprint(&command, " --mark-duplicates")
	}

	var opticalDuplicatesFilter func(alns *sam.Sam)

	if markOpticalDuplicates != "" {
		opticalDuplicatesFilter = func(alns *sam.Sam) {
			ctr := filters.MarkOpticalDuplicates(alns, pairs, opticalDuplicatesPixelDistance)
			filters.PrintDuplicatesMetrics(markOpticalDuplicates, command.String(), ctr)
		}
		fmt.Fprint(&command, " --mark-optical-duplicates ", markOpticalDuplicates)
		fmt.Fprint(&command, " --optical-duplicates-pixel-distance ", opticalDuplicatesPixelDistance)
	}

	if markOpticalDuplicatesIntermediate != "" {
		opticalDuplicatesFilter = func(alns *sam.Sam) {
			ctr := filters.MarkOpticalDuplicates(alns, pairs, opticalDuplicatesPixelDistance)
			filters.PrintDuplicatesMetricsToIntermediateFile(markOpticalDuplicatesIntermediate, ctr)
		}
		fmt.Fprint(&command, " --mark-optical-duplicates-intermediate ", markOpticalDuplicatesIntermediate)
		fmt.Fprint(&command, " --optical-duplicates-pixel-distance ", opticalDuplicatesPixelDistance)
	}

	if markOpticalDuplicates == "" && markOpticalDuplicatesIntermediate == "" {
		// nil tables for gc
		pairs = nil
	}

	filters1 = append(filters1, filters.RemoveOptionalReads)

	if removeDuplicates {
		filters2 = append(filters2, filters.RemoveDuplicateReads)
		fmt.Fprint(&command, " --remove-duplicates")
	}

	if bqsr != "" {
		// filters created later
		fmt.Fprint(&command, " --bqsr ", bqsr)
		fmt.Fprint(&command, " --quantize-levels ", quantizeLevels)
		fmt.Fprint(&command, " --max-cycle ", maxCycle)
	}

	if bqsrTablesOnly != "" {
		// filters created later
		fmt.Fprint(&command, " --bqsr-tables-only ", bqsrTablesOnly)
		fmt.Fprint(&command, " --max-cycle ", maxCycle)
	}

	if bqsrApplyFromTables != "" {
		// filters created later
		fmt.Fprint(&command, " --bqsr-apply ", bqsrApplyFromTables)
		fmt.Fprint(&command, " --quantize-levels ", quantizeLevels)
		fmt.Fprint(&command, " --max-cycle ", maxCycle)
	}

	if vcfOutput != "" {
		fmt.Fprint(&command, " --haplotypecaller ", vcfOutput)
		fmt.Fprint(&command, " --reference-confidence ", referenceConfidence)
		if sampleName != "" {
			fmt.Fprint(&command, " --sample-name ", sampleName)
		}
		if activityProfile != "" {
			fmt.Fprint(&command, " --activity-profile ", activityProfile)
		}
		if assemblyRegions != "" {
			fmt.Fprint(&command, " --assembly-regions ", assemblyRegions)
		}
		fmt.Fprint(&command, " --assembly-region-padding ", assemblyRegionPadding)
	}

	if bqsr != "" || bqsrTablesOnly != "" || vcfOutput != "" {
		fmt.Fprint(&command, " --reference ", referenceElFasta)
	}

	if spreadFile != "" {
		fmt.Fprint(&command, " --spread-file ", spreadFile)
	}

	var sqqList []uint8

	if (bqsr != "" || (bqsrApplyFromTables != "")) && sqq != "" {
		sqqs := strings.Split(sqq, ",")
		for _, sqq := range sqqs {
			i := internal.ParseUint(strings.TrimSpace(sqq), 10, 32)
			if i > 93 {
				log.Panicf("invalid sqq value %v", i)
			}
			sqqList = append(sqqList, uint8(i))
		}
		fmt.Fprint(&command, " --sqq ", sqq)
	}

	if (bqsr != "" || (bqsrTablesOnly != "")) && knownSites != "" {
		for i, k := range knownSitesList {
			knownSitesList[i] = strings.TrimSpace(k)
		}
		fmt.Fprint(&command, " --known-sites ", knownSites)
	}

	if (bqsr != "" || (bqsrApplyFromTables != "")) && recalFile != "" {
		fmt.Fprint(&command, " --recal-file ", recalFile)
	}

	if removeOptionalFields != "" {
		if removeOptionalFields == "all" {
			filters2 = append(filters2, filters.KeepOptionalFields(nil))
		} else {
			tags := strings.Split(removeOptionalFields, ",")
			for i, tag := range tags {
				tags[i] = strings.TrimSpace(tag)
			}
			filters2 = append(filters2, filters.RemoveOptionalFields(tags))
		}
		fmt.Fprint(&command, " --remove-optional-fields \"", removeOptionalFields, "\"")
	}

	if keepOptionalFields != "" {
		if keepOptionalFields == "none" {
			filters2 = append(filters2, filters.KeepOptionalFields(nil))
		} else {
			tags := strings.Split(keepOptionalFields, ",")
			for i, tag := range tags {
				tags[i] = strings.TrimSpace(tag)
			}
			filters2 = append(filters2, filters.KeepOptionalFields(tags))
		}
		fmt.Fprint(&command, " --keep-optional-fields \"", keepOptionalFields, "\"")
	}

	fmt.Fprint(&command, " --sorting-order ", sortingOrder)

	if nrOfThreads > 0 {
		runtime.GOMAXPROCS(nrOfThreads)
		fmt.Fprint(&command, " --nr-of-threads ", nrOfThreads)
	}

	if timed {
		fmt.Fprint(&command, " --timed")
	}

	if profile != "" {
		fmt.Fprint(&command, " --profile ", profile)
	}

	if logPath != "" {
		fmt.Fprint(&command, " --log-path ", logPath)
	}

	if randomSeedFile != "" {
		fmt.Fprint(&command, " --random-seed-file ", randomSeedFile)
	}

	if pgCmdLine != "" {
		fmt.Fprint(&command, " --pg-cmd-line ", pgCmdLine)
	}

	commandString := command.String()

	pgCmdString := commandString
	if pgCmdLine != "" {
		pgCmdString = pgCmdLine
	}
	filters1 = append([]sam.Filter{filters.AddPGLine(utils.StringMap{
		"ID": utils.ProgramName + " " + utils.ProgramVersion,
		"PN": utils.ProgramName,
		"VN": utils.ProgramVersion,
		"DS": utils.ProgramURL,
		"CL": pgCmdString,
	})}, filters1...)

	// executing command

	log.Println("Executing command:\n", commandString)

	var referenceMap *fasta.MappedFasta
	if referenceElFasta != "" {
		referenceMap = fasta.OpenElfasta(referenceElFasta)
		defer referenceMap.Close()
	}

	if bqsrTablesOnly != "" {
		baseRecalibrator := filters.NewBaseRecalibrator(knownSitesList, referenceMap)
		runBestPracticesPipelineIntermediateSamWithBQSRCalculateTablesOnly(
			input, output, outputType, sortingOrder,
			filters1, filters2, opticalDuplicatesFilter,
			baseRecalibrator, bqsrTablesOnly, maxCycle, targetRegionsBed,
			timed, profile)
		return
	}

	if bqsrApplyFromTables != "" {
		recalFile := internal.FilepathAbs(recalFile)
		baseRecalibratorTables := filters.LoadAndCombineBQSRTables(bqsrApplyFromTables)
		filters2 = append(filters2, baseRecalibratorTables.ApplyBQSR(quantizeLevels, sqqList, maxCycle))
		if vcfOutput != "" {
			var activityProfileOut, assemblyRegionsOut io.WriteCloser
			if activityProfile != "" {
				activityProfileOut = internal.FileCreate(activityProfile)
				defer internal.Close(activityProfileOut)
			}
			if assemblyRegions != "" {
				assemblyRegionsOut = internal.FileCreate(assemblyRegions)
				defer internal.Close(assemblyRegionsOut)
			}
			haplotypeCaller := filters.NewHaplotypeCaller(
				referenceMap, referenceConfidence, int32(assemblyRegionPadding),
				activityProfileOut, assemblyRegionsOut,
				randomSeedFile, pgCmdString,
			)
			defer haplotypeCaller.Close()
			runBestPracticesPipelineIntermediateSamWithBQSRApplyOnly(
				input, output, outputType, sortingOrder,
				append(filters2, filters.AddREFID), baseRecalibratorTables, recalFile,
				haplotypeCaller, sampleName, vcfOutput, targetRegionsBed, spreadFile,
				timed, profile)
			return
		}
		runBestPracticesPipelineWithBQSRApplyOnly(
			input, output, outputType, sortingOrder,
			filters2, baseRecalibratorTables, recalFile,
			timed, profile)
		return
	}

	if markDuplicates ||
		sortingOrder == sam.Coordinate || sortingOrder == sam.Queryname ||
		(replaceReferenceSequences != "" && sortingOrder == sam.Keep) ||
		bqsr != "" || vcfOutput != "" {
		var baseRecalibrator *filters.BaseRecalibrator
		if bqsr != "" {
			baseRecalibrator = filters.NewBaseRecalibrator(knownSitesList, referenceMap)
		}
		var haplotypeCaller *filters.HaplotypeCaller
		if vcfOutput != "" {
			var activityProfileOut, assemblyRegionsOut io.WriteCloser
			if activityProfile != "" {
				activityProfileOut = internal.FileCreate(activityProfile)
				defer internal.Close(activityProfileOut)
			}
			if assemblyRegions != "" {
				assemblyRegionsOut = internal.FileCreate(assemblyRegions)
				defer internal.Close(assemblyRegionsOut)
			}
			haplotypeCaller = filters.NewHaplotypeCaller(
				referenceMap, referenceConfidence, int32(assemblyRegionPadding),
				activityProfileOut, assemblyRegionsOut,
				randomSeedFile, pgCmdString,
			)
			defer haplotypeCaller.Close()
		}
		runBestPracticesPipelineIntermediateSam(
			input, output, outputType, sortingOrder,
			filters1, filters2,
			opticalDuplicatesFilter,
			baseRecalibrator, quantizeLevels, sqqList, maxCycle, internal.FilepathAbs(bqsr),
			haplotypeCaller, sampleName, vcfOutput, targetRegionsBed, spreadFile,
			timed, profile)
		return
	}

	runBestPracticesPipeline(input, output, outputType, sortingOrder, append(filters1, filters2...), timed, profile)
}
