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
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strings"

	"github.com/exascience/elprep/v5/internal"
	"github.com/exascience/elprep/v5/sam"
)

// MergeHelp is the help string for this command.
const MergeHelp = "\nmerge parameters:\n" +
	"elprep merge /path/to/input sam-output-file\n" +
	"[--output-type [sam | bam]]\n" +
	"[--single-end]\n" +
	"[--ignore-spread-file] \n" +
	"[--nr-of-threads n]\n" +
	"[--timed]\n" +
	"[--log-path path]\n"

// Merge implements the elprep merge command.
func Merge() {
	var (
		outputType                         string
		profile, logPath                   string
		nrOfThreads                        int
		singleEnd, timed, ignoreSpreadFile bool
	)

	var flags flag.FlagSet

	flags.StringVar(&outputType, "output-type", "", "format of the output file")
	flags.BoolVar(&singleEnd, "single-end", false, "when splitting single-end data")
	flags.BoolVar(&ignoreSpreadFile, "ignore-spread-file", false, "when merging data after applying haplotype caller")
	flags.IntVar(&nrOfThreads, "nr-of-threads", 0, "number of worker threads")
	flags.BoolVar(&timed, "timed", false, "measure the runtime")
	flags.StringVar(&profile, "profile", "", "write a runtime profile to the specified file(s)")
	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")

	parseFlags(flags, 4, MergeHelp)

	input := getFilename(os.Args[2], MergeHelp)
	output := getFilename(os.Args[3], MergeHelp)

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

	if profile != "" && !checkCreate("--profile", profile) {
		sanityChecksFailed = true
	}

	fullInputPath := internal.FilepathAbs(input)
	fullInputPath, filesToMerge := internal.Directory(fullInputPath)
	if len(filesToMerge) < 2 && !singleEnd {
		log.Printf("Given directory %v does not contain /splits/ directory and/or spread reads file. These should have been created by an elprep split invocation.\n", input)
		sanityChecksFailed = true
	}

	if nrOfThreads < 0 {
		sanityChecksFailed = true
		log.Println("Error: Invalid nr-of-threads: ", nrOfThreads)
	}

	if sanityChecksFailed {
		fmt.Fprint(os.Stderr, MergeHelp)
		os.Exit(1)
	}

	firstFile := sam.Open(filepath.Join(fullInputPath, filesToMerge[0]))
	header := firstFile.ParseHeader()
	firstFile.Close()

	var inputPrefix string
	for _, file := range filesToMerge {
		if idx := strings.LastIndex(file, "-unmapped"); idx >= 0 {
			inputPrefix = file[:idx]
			break
		}
	}

	var inputExtension string
	switch ext := filepath.Ext(filesToMerge[0]); ext {
	case sam.SamExt, sam.BamExt:
		inputExtension = ext[1:]
	default:
		inputExtension = "sam"
	}

	// building output command line

	var command bytes.Buffer
	fmt.Fprint(&command, os.Args[0], " merge ", input, " ", output)
	if outputType != "" {
		fmt.Fprint(&command, " --output-type ", outputType)
	}
	if singleEnd {
		fmt.Fprint(&command, " --single-end ")
	}
	if ignoreSpreadFile {
		fmt.Fprint(&command, " --ignore-spread-file ")
	}
	if nrOfThreads > 0 {
		runtime.GOMAXPROCS(nrOfThreads)
		fmt.Fprint(&command, " --nr-of-threads ", nrOfThreads)
	}
	if timed {
		fmt.Fprint(&command, " --timed ")
	}
	if logPath != "" {
		fmt.Fprint(&command, " --log-path ", logPath)
	}

	// executing command

	log.Println("Executing command:\n", command.String())

	output = internal.FilepathAbs(output)

	switch header.HDSO() {
	case sam.Coordinate:
		if singleEnd {
			timedRun(timed, profile, "Merging single-end split files sorted by coordinate.", 1, func() {
				sam.MergeSingleEndFilesSplitPerChromosome(fullInputPath, output, inputPrefix, inputExtension, outputType, header)
			})
			return
		}
		// check if we need to merge in the spread file
		if ignoreSpreadFile {
			timedRun(timed, profile, "Merging paired-end split files sorted by coordinate.", 1, func() {
				sam.MergeSortedFilesSplitPerChromosomeWithoutSpreadFile(fullInputPath, output, inputPrefix, inputExtension, outputType, header)
			})
		} else {
			timedRun(timed, profile, "Merging paired-end split files sorted by coordinate.", 1, func() {
				sam.MergeSortedFilesSplitPerChromosome(fullInputPath, output, inputPrefix, inputExtension, outputType, header)
			})
		}
	case sam.Queryname:
		log.Panic("Merging of files sorted by queryname not yet implemented.")
		return
	default:
		if singleEnd {
			timedRun(timed, profile, "Merging unsorted single-end split files.", 1, func() {
				sam.MergeSingleEndFilesSplitPerChromosome(fullInputPath, output, inputPrefix, inputExtension, outputType, header)
			})
			return
		}
		timedRun(timed, profile, "Merging unsorted paired-end split files.", 1, func() {
			sam.MergeUnsortedFilesSplitPerChromosome(fullInputPath, output, inputPrefix, inputExtension, outputType, header)
		})
		return
	}
}
