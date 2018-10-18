// elPrep: a high-performance tool for preparing SAM/BAM files.
// Copyright (c) 2017, 2018 imec vzw.

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
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"

	"github.com/exascience/elprep/v4/bed"
	"github.com/exascience/elprep/v4/filters"
	"github.com/exascience/elprep/v4/sam"
	"github.com/exascience/elprep/v4/utils"
	psync "github.com/exascience/pargo/sync"
)

// Run the best practices pipeline. Version that uses an intermediate
// slice so that sorting and mark-duplicates are supported.
func runBestPracticesPipelineIntermediateSam(fileIn, fileOut string, sortingOrder sam.SortingOrder, filters, filters2 []sam.Filter, opticalDuplicateFilter func(reads *sam.Sam) error, deterministic, timed bool, profile string) error {
	filteredReads := sam.NewSam()
	phase := int64(1)
	err := timedRun(timed, profile, "Reading SAM into memory and applying filters.", phase, func() (err error) {
		pathname, err := filepath.Abs(fileIn)
		if err != nil {
			return err
		}
		input, err := sam.Open(pathname)
		if err != nil {
			return err
		}
		defer func() {
			nerr := input.Close()
			if err == nil {
				err = nerr
			}
		}()
		return input.RunPipelineFI(filteredReads, filters, sortingOrder, deterministic)
	})
	if err != nil {
		return err
	}
	go runtime.GC()
	if opticalDuplicateFilter != nil {
		phase++
		err = timedRun(timed, profile, "Marking optical duplicates.", phase, func() (err error) {
			return opticalDuplicateFilter(filteredReads)
		})
		if err != nil {
			return err
		}
		go runtime.GC()
	}
	phase++
	return timedRun(timed, profile, "Write to file.", phase, func() (err error) {
		pathname, err := filepath.Abs(fileOut)
		if err != nil {
			return err
		}
		err = os.MkdirAll(filepath.Dir(pathname), 0700)
		if err != nil {
			return err
		}
		output, err := sam.Create(pathname)
		if err != nil {
			return err
		}
		defer func() {
			nerr := output.Close()
			if err == nil {
				err = nerr
			}
		}()
		if sortingOrder != sam.Unsorted {
			sortingOrder = sam.Keep
		}
		return filteredReads.RunPipeline(output, filters2, sortingOrder)
	})
}

func runBestPracticesPipelineIntermediateSamWithBQSR(fileIn, fileOut string, sortingOrder sam.SortingOrder, filters1, filters2 []sam.Filter, opticalDuplicateFilter func(reads *sam.Sam) error, baseRecalibrator *filters.BaseRecalibrator, quantizeLevels int, sqqList []uint8, recalFile string, deterministic, timed bool, profile string) error {
	// Input and first filters and sorting
	filteredReads := sam.NewSam()
	phase := int64(1)
	err := timedRun(timed, profile, "Reading SAM into memory and applying filters.", phase, func() (err error) {
		pathname, err := filepath.Abs(fileIn)
		if err != nil {
			return err
		}
		input, err := sam.Open(pathname)
		if err != nil {
			return err
		}
		defer func() {
			nerr := input.Close()
			if err == nil {
				err = nerr
			}
		}()
		return input.RunPipelineFI(filteredReads, filters1, sortingOrder, deterministic)
	})
	if err != nil {
		return err
	}
	go runtime.GC()
	if opticalDuplicateFilter != nil {
		phase++
		err = timedRun(timed, profile, "Marking optical duplicates.", phase, func() (err error) {
			return opticalDuplicateFilter(filteredReads)
		})
		if err != nil {
			return err
		}
		go runtime.GC()
	}
	if sortingOrder != sam.Unsorted {
		sortingOrder = sam.Keep
	}
	// Base recalibration
	var baseRecalibratorTables filters.BaseRecalibratorTables
	phase++
	err = timedRun(timed, profile, "Base recalibration", phase, func() error {
		baseRecalibratorTables = baseRecalibrator.Recalibrate(filteredReads)
		return baseRecalibratorTables.Err()
	})
	if err != nil {
		return err
	}
	// Finalize BQSR tables + log recal file
	phase++
	err = timedRun(timed, profile, "Finalize BQSR tables", phase, func() error {
		baseRecalibratorTables.FinalizeBQSRTables()
		pathname, err := filepath.Abs(recalFile)
		if err != nil {
			return err
		}
		err = os.MkdirAll(filepath.Dir(pathname), 0700)
		if err != nil {
			return err
		}
		return baseRecalibratorTables.PrintBQSRTables(pathname)
	})
	if err != nil {
		return err
	}
	go runtime.GC()
	phase++
	err = timedRun(timed, profile, "Apply BQSR", phase, func() error {
		n := len(filteredReads.Alignments) / 4096
		if m := 2 * runtime.GOMAXPROCS(0); m > n {
			n = m
		}
		filteredReads.NofBatches(n)
		return filteredReads.RunPipeline(filteredReads, []sam.Filter{baseRecalibratorTables.ApplyBQSR(quantizeLevels, sqqList)}, sortingOrder)
	})
	if err != nil {
		return err
	}
	go runtime.GC()
	// Write output to file
	phase++
	return timedRun(timed, profile, "Write to file.", phase, func() (err error) {
		pathname, err := filepath.Abs(fileOut)
		if err != nil {
			return err
		}
		err = os.MkdirAll(filepath.Dir(pathname), 0700)
		if err != nil {
			return err
		}
		output, err := sam.Create(pathname)
		if err != nil {
			return err
		}
		defer func() {
			nerr := output.Close()
			if err == nil {
				err = nerr
			}
		}()
		return filteredReads.RunPipeline(output, filters2, sortingOrder)
	})
}

func runBestPracticesPipelineIntermediateSamWithBQSRApplyOnly(input, output string, sortingOrder sam.SortingOrder, filters []sam.Filter, baseRecalibratorTables filters.BaseRecalibratorTables, recalFile string, timed bool, profile string) error {
	// Finalize BQSR tables + log recal file
	err := timedRun(timed, profile, "Finalize BQSR tables", 1, func() error {
		baseRecalibratorTables.FinalizeBQSRTables()
		pathname, err := filepath.Abs(recalFile)
		if err != nil {
			return err
		}
		err = os.MkdirAll(filepath.Dir(pathname), 0700)
		if err != nil {
			return err
		}
		return baseRecalibratorTables.PrintBQSRTables(recalFile)
	})
	if err != nil {
		return err
	}
	return timedRun(timed, profile, "Running pipeline.", 2, func() (err error) {
		pathname, err := filepath.Abs(input)
		if err != nil {
			return err
		}
		input, err := sam.Open(pathname)
		if err != nil {
			return err
		}
		defer func() {
			nerr := input.Close()
			if err == nil {
				err = nerr
			}
		}()
		pathname, err = filepath.Abs(output)
		if err != nil {
			return err
		}
		err = os.MkdirAll(filepath.Dir(pathname), 0700)
		if err != nil {
			return err
		}
		output, err := sam.Create(pathname)
		if err != nil {
			return err
		}
		defer func() {
			nerr := output.Close()
			if err == nil {
				err = nerr
			}
		}()
		return input.RunPipeline(output, filters, sortingOrder)
	})
}

func runBestPracticesPipelineIntermediateSamWithBQSRCalculateTablesOnly(fileIn, fileOut string, sortingOrder sam.SortingOrder, filters1, filters2 []sam.Filter, opticalDuplicateFilter func(reads *sam.Sam) error, baseRecalibrator *filters.BaseRecalibrator, tableFile string, timed bool, profile string) error {
	// Input and first filters and sorting
	filteredReads := sam.NewSam()
	phase := int64(1)
	err := timedRun(timed, profile, "Reading SAM into memory and applying filters.", phase, func() (err error) {
		pathname, err := filepath.Abs(fileIn)
		if err != nil {
			return err
		}
		input, err := sam.Open(pathname)
		if err != nil {
			return err
		}
		defer func() {
			nerr := input.Close()
			if err == nil {
				err = nerr
			}
		}()
		return input.RunPipeline(filteredReads, filters1, sortingOrder)
	})
	if err != nil {
		return err
	}
	go runtime.GC()
	// Marking optical duplicates
	if opticalDuplicateFilter != nil {
		phase++
		err = timedRun(timed, profile, "Marking optical duplicates.", phase, func() (err error) {
			return opticalDuplicateFilter(filteredReads)
		})
		if err != nil {
			return err
		}
		go runtime.GC()
	}
	// Base recalibration
	var baseRecalibratorTables filters.BaseRecalibratorTables
	phase++
	err = timedRun(timed, profile, "Base recalibration", phase, func() error {
		baseRecalibratorTables = baseRecalibrator.Recalibrate(filteredReads)
		return baseRecalibratorTables.Err()
	})
	if err != nil {
		return err
	}
	// Print bqsr table to file
	phase++
	err = timedRun(timed, profile, "Print intermediate BQSR tables", phase, func() error {
		pathname, err := filepath.Abs(tableFile)
		if err != nil {
			return err
		}
		err = os.MkdirAll(filepath.Dir(pathname), 0700)
		if err != nil {
			return err
		}
		return baseRecalibratorTables.PrintBQSRTablesToIntermediateFile(pathname)
	})
	if err != nil {
		return err
	}
	go runtime.GC()
	// Write output to file
	phase++
	return timedRun(timed, profile, "Write to file.", phase, func() (err error) {
		pathname, err := filepath.Abs(fileOut)
		if err != nil {
			return err
		}
		err = os.MkdirAll(filepath.Dir(pathname), 0700)
		if err != nil {
			return err
		}
		output, err := sam.Create(pathname)
		if err != nil {
			return err
		}
		defer func() {
			nerr := output.Close()
			if err == nil {
				err = nerr
			}
		}()
		if sortingOrder != sam.Unsorted {
			sortingOrder = sam.Keep
		}
		return filteredReads.RunPipeline(output, filters2, sortingOrder)
	})
}

// Run the best practices pipeline. Version that doesn't use an
// intermediate slice when neither sorting nor mark-duplicates are
// needed.
func runBestPracticesPipeline(fileIn, fileOut string, sortingOrder sam.SortingOrder, filters []sam.Filter, timed bool, profile string) error {
	return timedRun(timed, profile, "Running pipeline.", 1, func() (err error) {
		pathname, err := filepath.Abs(fileIn)
		if err != nil {
			return err
		}
		input, err := sam.Open(pathname)
		if err != nil {
			return err
		}
		defer func() {
			nerr := input.Close()
			if err == nil {
				err = nerr
			}
		}()
		pathname, err = filepath.Abs(fileOut)
		if err != nil {
			return err
		}
		err = os.MkdirAll(filepath.Dir(pathname), 0700)
		if err != nil {
			return err
		}
		output, err := sam.Create(pathname)
		if err != nil {
			return err
		}
		defer func() {
			nerr := output.Close()
			if err == nil {
				err = nerr
			}
		}()
		return input.RunPipeline(output, filters, sortingOrder)
	})
}

// FilterHelp is the help string for this command.
const FilterHelp = "\nfilter parameters:\n" +
	"elprep filter sam-file sam-output-file\n" +
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
	"[--remove-duplicates]\n" +
	"[--remove-optional-fields [all | list]]\n" +
	"[--keep-optional-fields [none | list]]\n" +
	"[--sorting-order [keep | unknown | unsorted | queryname | coordinate]]\n" +
	"[--clean-sam]\n" +
	"[--bqsr recal-file]\n" +
	"[--bqsr-reference elfasta]\n" +
	"[--quantize-levels nr]\n" +
	"[--sqq list]\n" +
	"[--known-sites list]\n" +
	"[--nr-of-threads nr]\n" +
	"[--timed]\n" +
	"[--log-path path]\n"

	//FilterExtendedHelp is the extended help string for this command.
const FilterExtendedHelp = FilterHelp +
	"[--mark-optical-duplicates-intermediate file]\n" +
	"[--bqsr-tables-only table-file]\n" +
	"[--bqsr-apply path]\n" +
	"[--recal-file file]\n"

// Filter implements the elprep filter command.
func Filter() error {
	var (
		replaceReferenceSequences                                string
		filterUnmappedReads, filterUnmappedReadsStrict           bool
		filterMappingQuality                                     int
		filterNonExactMappingReads                               bool
		filterNonExactMappingReadsStrict                         bool
		filterNonOverlappingReads                                string
		replaceReadGroup                                         string
		markDuplicates, markDuplicatesDet, removeDuplicates      bool
		markOpticalDuplicates, markOpticalDuplicatesIntermediate string
		removeOptionalFields                                     string
		keepOptionalFields                                       string
		sortingOrderString                                       string
		cleanSam                                                 bool
		bqsr                                                     string
		referenceElFasta                                         string
		bqsrTablesOnly                                           string
		bqsrApplyFromTables                                      string
		quantizeLevels                                           int
		sqq                                                      string
		knownSites                                               string
		recalFile                                                string
		deterministic                                            bool
		nrOfThreads                                              int
		timed                                                    bool
		profile                                                  string
		logPath                                                  string
		renameChromosomes                                        bool
	)

	var flags flag.FlagSet

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
	flags.BoolVar(&markDuplicatesDet, "mark-duplicates-deterministic", false, "mark duplicates deterministically")
	flags.BoolVar(&removeDuplicates, "remove-duplicates", false, "remove duplicates")
	flags.StringVar(&removeOptionalFields, "remove-optional-fields", "", "remove the given optional fields")
	flags.StringVar(&keepOptionalFields, "keep-optional-fields", "", "remove all except for the given optional fields")
	flags.StringVar(&sortingOrderString, "sorting-order", string(sam.Keep), "determine output order of alignments, one of keep, unknown, unsorted, queryname, or coordinate")
	flags.BoolVar(&cleanSam, "clean-sam", false, "clean the sam file")
	flags.StringVar(&bqsr, "bqsr", "", "base quality score recalibration")
	flags.StringVar(&bqsrTablesOnly, "bqsr-tables-only", "", "base quality score recalibration table calculation (only with split/merge)")
	flags.StringVar(&bqsrApplyFromTables, "bqsr-apply", "", "base quality score recalibration application (only with split/merge)")
	flags.StringVar(&referenceElFasta, "bqsr-reference", "", "reference used for base quality score recalibration (elfasta format)")
	flags.IntVar(&quantizeLevels, "quantize-levels", 0, "number of levels to be used for quantizing recalibrated base qualities (only with --bqsr)")
	flags.StringVar(&sqq, "sqq", "", "levels to be used for statically quantizing recalibrated base qualities (only with --bqsr)")
	flags.StringVar(&knownSites, "known-sites", "", "list of vcf files containing known sites for base recalibration (only with --bqsr)")
	flags.StringVar(&recalFile, "recal-file", "", "log file for recalibration tables (only with --bqsr-apply)")
	flags.BoolVar(&deterministic, "deterministic", false, "run elprep deterministically")
	flags.IntVar(&nrOfThreads, "nr-of-threads", 0, "number of worker threads")
	flags.BoolVar(&timed, "timed", false, "measure the runtime")
	flags.StringVar(&profile, "profile", "", "write a runtime profile to the specified file(s)")
	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")
	flags.BoolVar(&renameChromosomes, "rename-chromosomes", false, "")

	parseFlags(flags, 4, FilterHelp)

	input := getFilename(os.Args[2], FilterHelp)
	output := getFilename(os.Args[3], FilterHelp)

	setLogOutput(logPath)

	// sanity checks

	var sanityChecksFailed bool

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
	if referenceElFasta != "" && !checkExist("--bqsr-reference", referenceElFasta) {
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

	if markDuplicatesDet {
		log.Println("--mark-duplicates-deterministic is deprecated; just use --mark-duplicates instead.")
		markDuplicates = true
	}

	if markOpticalDuplicates != "" && !markDuplicates {
		sanityChecksFailed = true
		log.Println("Error: Cannot use --mark-optical-duplicates without also using --mark-duplicates.")
	}

	if markOpticalDuplicatesIntermediate != "" && !markDuplicates {
		sanityChecksFailed = true
		log.Println("Error: Cannot use --mark-optical-duplicates-intermediate without also using --mark-duplicates.")
	}

	if deterministic && (bqsrTablesOnly != "" || bqsrApplyFromTables != "") {
		sanityChecksFailed = true
		log.Println("Error: deterministic option is not yet supported for --bqsr-tables-only or --bqsr-apply")
	}

	if !checkBQSROptions(bqsr != "", (bqsrTablesOnly != ""), referenceElFasta, quantizeLevels, sqq, knownSites, bqsr, recalFile) {
		sanityChecksFailed = true
	}
	if !checkBQSRTablesOnlyOptions((bqsrTablesOnly != ""), bqsrTablesOnly, referenceElFasta) {
		sanityChecksFailed = true
	}
	if !checkBQSRApplyOptions((bqsrApplyFromTables != ""), bqsrApplyFromTables, recalFile) {
		sanityChecksFailed = true
	}

	if sanityChecksFailed {
		fmt.Fprint(os.Stderr, FilterHelp)
		os.Exit(1)
	}

	// building filters and output command line

	var command bytes.Buffer
	fmt.Fprint(&command, os.Args[0], " filter ", input, " ", output)

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

	if filterNonOverlappingReads != "" {
		parsedBed, err := bed.ParseBed(filterNonOverlappingReads)
		if err != nil {
			return err
		}
		filterNonOverlappingReadsFilter := filters.RemoveNonOverlappingReads(parsedBed)
		filters1 = append(filters1, filterNonOverlappingReadsFilter)
		fmt.Fprint(&command, " --filter-non-overlapping-reads ", filterNonOverlappingReads)
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
		replaceReferenceSequencesFilter, err := filters.ReplaceReferenceSequenceDictionaryFromSamFile(replaceReferenceSequences)
		if err != nil {
			return err
		}
		filters1 = append(filters1, replaceReferenceSequencesFilter)
		fmt.Fprint(&command, " --replace-reference-sequences ", replaceReferenceSequences)
	}

	if replaceReadGroup != "" {
		record, err := sam.ParseHeaderLineFromString(replaceReadGroup)
		if err != nil {
			return err
		}
		filters1 = append(filters1, filters.AddOrReplaceReadGroup(record))
		fmt.Fprint(&command, " --replace-read-group ", replaceReadGroup)
	}

	if (replaceReferenceSequences != "") || markDuplicates ||
		(sortingOrder == sam.Coordinate) || (sortingOrder == sam.Queryname) {
		filters1 = append(filters1, filters.AddREFID)
	}

	var fragments, pairs *psync.Map
	var markDupsFilter sam.Filter

	if markDuplicates {
		alsoOpticals := markOpticalDuplicates != "" || markOpticalDuplicatesIntermediate != ""
		markDupsFilter, fragments, pairs = filters.MarkDuplicates(alsoOpticals)
		filters1 = append(filters1, markDupsFilter)
		fmt.Fprint(&command, " --mark-duplicates")
	}

	var opticalDuplicatesFilter func(alns *sam.Sam) error

	if markOpticalDuplicates != "" {
		opticalDuplicatesFilter = func(alns *sam.Sam) error {
			ctr := filters.MarkOpticalDuplicates(alns, fragments, pairs, deterministic)
			if err := ctr.Err(); err != nil {
				return err
			}
			return filters.PrintDuplicatesMetrics(input, output, markOpticalDuplicates, removeDuplicates, ctr)
		}
		fmt.Fprint(&command, " --mark-optical-duplicates ", markOpticalDuplicates)
	}

	if markOpticalDuplicatesIntermediate != "" {
		opticalDuplicatesFilter = func(alns *sam.Sam) error {
			ctr := filters.MarkOpticalDuplicates(alns, fragments, pairs, deterministic)
			if err := ctr.Err(); err != nil {
				return err
			}
			return filters.PrintDuplicatesMetricsToIntermediateFile(markOpticalDuplicatesIntermediate, ctr)
		}
		fmt.Fprint(&command, " --mark-optical-duplicates-intermediate ", markOpticalDuplicatesIntermediate)
	}

	if markOpticalDuplicates == "" && markOpticalDuplicatesIntermediate == "" {
		// nil tables for gc
		fragments = nil
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
		fmt.Fprint(&command, " --bqsr-reference ", referenceElFasta)
		fmt.Fprint(&command, " --quantize-levels ", quantizeLevels)
	}

	if bqsrTablesOnly != "" {
		// filters created later
		fmt.Fprint(&command, " --bqsr-tables-only ", bqsrTablesOnly)
		fmt.Fprint(&command, " --bqsr-reference ", referenceElFasta)
		fmt.Fprint(&command, " --quantize-levels ", quantizeLevels)
	}

	if bqsrApplyFromTables != "" {
		// filters created later
		fmt.Fprint(&command, " --bqsr-apply ", bqsrApplyFromTables)
		fmt.Fprint(&command, " --quantize-levels ", quantizeLevels)
	}

	var sqqList []uint8

	if (bqsr != "" || (bqsrTablesOnly != "")) && sqq != "" {
		sqqs := strings.Split(sqq, ",")
		for _, sqq := range sqqs {
			if i, err := strconv.ParseUint(strings.TrimSpace(sqq), 10, 32); err != nil {
				return err
			} else if i > 93 {
				return fmt.Errorf("invalid sqq value %v", i)
			} else {
				sqqList = append(sqqList, uint8(i))
			}
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

	if deterministic {
		fmt.Fprint(&command, " --deterministic")
	}

	// Currently, we need the deterministic option only when optical duplicates are marked.
	// This might change in the future.
	if opticalDuplicatesFilter == nil {
		deterministic = false
	}

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

	commandString := command.String()

	filters1 = append([]sam.Filter{filters.AddPGLine(utils.StringMap{
		"ID": ProgramName + " " + ProgramVersion,
		"PN": ProgramName,
		"VN": ProgramVersion,
		"DS": ProgramURL,
		"CL": commandString,
	})}, filters1...)

	// executing command

	log.Println("Executing command:\n", commandString)

	if bqsr != "" {
		recalFile, err := filepath.Abs(bqsr)
		if err != nil {
			return err
		}
		baseRecalibrator := filters.NewBaseRecalibrator(knownSitesList, referenceElFasta)
		return runBestPracticesPipelineIntermediateSamWithBQSR(input, output, sortingOrder, filters1, filters2, opticalDuplicatesFilter, baseRecalibrator, quantizeLevels, sqqList, recalFile, deterministic, timed, profile)
	}

	if bqsrTablesOnly != "" {
		baseRecalibrator := filters.NewBaseRecalibrator(knownSitesList, referenceElFasta)
		return runBestPracticesPipelineIntermediateSamWithBQSRCalculateTablesOnly(input, output, sortingOrder, filters1, filters2, opticalDuplicatesFilter, baseRecalibrator, bqsrTablesOnly, timed, profile)
	}

	if bqsrApplyFromTables != "" {
		recalFile, err := filepath.Abs(recalFile)
		if err != nil {
			return err
		}
		baseRecalibratorTables, err := filters.LoadAndCombineBQSRTables(bqsrApplyFromTables)
		if err != nil {
			return err
		}
		filters2 = append(filters2, baseRecalibratorTables.ApplyBQSR(quantizeLevels, sqqList))
		return runBestPracticesPipelineIntermediateSamWithBQSRApplyOnly(input, output, sortingOrder, filters2, baseRecalibratorTables, recalFile, timed, profile)
	}

	if markDuplicates ||
		(sortingOrder == sam.Coordinate) || (sortingOrder == sam.Queryname) ||
		((replaceReferenceSequences != "") && (sortingOrder == sam.Keep)) {
		return runBestPracticesPipelineIntermediateSam(input, output, sortingOrder, filters1, filters2, opticalDuplicatesFilter, deterministic, timed, profile)
	}
	return runBestPracticesPipeline(input, output, sortingOrder, append(filters1, filters2...), timed, profile)
}
