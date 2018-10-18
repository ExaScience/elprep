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

	"github.com/exascience/elprep/v4/sam"
)

// SplitHelp is the help string for this command.
const SplitHelp = "\nsplit parameters:\n" +
	"elprep split (sam-file | /path/to/input/) /path/to/output/\n" +
	"[--output-prefix name]\n" +
	"[--output-type [sam | bam]]\n" +
	"[--single-end]\n" +
	"[--nr-of-threads nr]\n" +
	"[--timed]\n" +
	"[--log-path path]\n" +
	"[--contig-group-size nr]\n"

// Split implements the elprep split command.
func Split() error {
	var (
		contigGroupSize                            int
		outputPrefix, outputType, profile, logPath string
		nrOfThreads                                int
		singleEnd, timed                           bool
	)

	var flags flag.FlagSet

	flags.IntVar(&contigGroupSize, "contig-group-size", 0, "maximum sum of reference sequence lengths for creating groups of reference sequences")
	flags.StringVar(&outputPrefix, "output-prefix", "", "prefix for the output files")
	flags.StringVar(&outputType, "output-type", "", "format of the output files")
	flags.BoolVar(&singleEnd, "single-end", false, "when splitting single-end data")
	flags.IntVar(&nrOfThreads, "nr-of-threads", 0, "number of worker threads")
	flags.BoolVar(&timed, "timed", false, "measure the runtime")
	flags.StringVar(&profile, "profile", "", "write a runtime profile to the specified file(s)")
	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")

	parseFlags(flags, 4, SplitHelp)

	input := getFilename(os.Args[2], SplitHelp)
	output := getFilename(os.Args[3], SplitHelp)

	ext := filepath.Ext(input)
	if outputPrefix == "" {
		base := filepath.Base(input)
		outputPrefix = base[:len(base)-len(ext)]
	}
	if outputType == "" {
		switch ext {
		case sam.SamExt, sam.BamExt:
			outputType = ext[1:]
		default:
			outputType = "sam"
		}
	}

	setLogOutput(logPath)

	// sanity checks

	var sanityChecksFailed bool

	if !checkExist("", input) {
		sanityChecksFailed = true
	}

	if filepath.Dir(output) != filepath.Clean(output) {
		log.Printf("Given output path is not a path: %v.\n", output)
		sanityChecksFailed = true
	}

	if profile != "" && !checkCreate("--profile", profile) {
		sanityChecksFailed = true
	}

	if nrOfThreads < 0 {
		sanityChecksFailed = true
		log.Println("Error: Invalid nr-of-threads: ", nrOfThreads)
	}

	if sanityChecksFailed {
		fmt.Fprint(os.Stderr, SplitHelp)
		os.Exit(1)
	}

	// building output command line

	var command bytes.Buffer
	fmt.Fprint(&command, os.Args[0], " split ", input, " ", output)
	fmt.Fprint(&command, " --output-prefix ", outputPrefix)
	fmt.Fprint(&command, " --output-type ", outputType)
	if singleEnd {
		fmt.Fprint(&command, " --single-end ")
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

	fullInput, err := filepath.Abs(input)
	if err != nil {
		return err
	}

	fullOutput, err := filepath.Abs(output)
	if err != nil {
		return err
	}

	err = os.MkdirAll(output, 0700)
	if err != nil {
		return err
	}

	if singleEnd {
		err := timedRun(timed, profile, "Splitting single-end files.", 1, func() (err error) {
			return sam.SplitSingleEndFilePerChromosome(fullInput, fullOutput, outputPrefix, outputType, contigGroupSize)
		})
		return err
	}
	err = timedRun(timed, profile, "Splitting paired-end files.", 1, func() (err error) {
		return sam.SplitFilePerChromosome(fullInput, fullOutput, outputPrefix, outputType, contigGroupSize)
	})
	return err
}
