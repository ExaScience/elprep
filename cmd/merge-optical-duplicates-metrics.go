// elPrep: a high-performance tool for preparing SAM/BAM files.
// Copyright (c) 2017-2019 imec vzw.

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

	"github.com/exascience/elprep/v4/filters"
)

// MergeOpticalDuplicatesMetricsHelp is the help string for this command.
const MergeOpticalDuplicatesMetricsHelp = "\nmerge-optical-duplicates-metrics parameters:\n" +
	"elprep merge-optical-duplicates-metrics sam-input-file sam-output-file metrics-file /path/to/intermediate/metrics\n" +
	"[--remove-duplicates]\n" +
	"[--nr-of-threads nr]\n" +
	"[--timed]\n" +
	"[--log-path path]\n"

// Merge implements the elprep merge command.
func MergeOpticalDuplicatesMetrics() error {
	var (
		profile, logPath        string
		nrOfThreads             int
		timed, removeDuplicates bool
	)

	var flags flag.FlagSet

	flags.IntVar(&nrOfThreads, "nr-of-threads", 0, "number of worker threads")
	flags.BoolVar(&timed, "timed", false, "measure the runtime")
	flags.BoolVar(&removeDuplicates, "remove-duplicates", false, "use when duplicates were removed during duplicate marking")
	flags.StringVar(&profile, "profile", "", "write a runtime profile to the specified file(s)")
	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")

	parseFlags(flags, 6, MergeOpticalDuplicatesMetricsHelp)

	input := getFilename(os.Args[2], MergeOpticalDuplicatesMetricsHelp)
	output := getFilename(os.Args[3], MergeOpticalDuplicatesMetricsHelp)
	metrics := getFilename(os.Args[4], MergeOpticalDuplicatesMetricsHelp)
	intermediateMetrics := getFilename(os.Args[5], MergeOpticalDuplicatesMetricsHelp)

	setLogOutput(logPath)

	// sanity checks

	var sanityChecksFailed bool

	if !checkExist("", input) {
		log.Println("Warning: Input file does not exist: ", input)
	}

	if !checkExist("", intermediateMetrics) {
		sanityChecksFailed = true
	}

	if profile != "" && !checkCreate("--profile", profile) {
		sanityChecksFailed = true
	}

	metricsDir, err := filepath.Abs(intermediateMetrics)
	if err != nil {
		return err
	}

	if nrOfThreads < 0 {
		sanityChecksFailed = true
		log.Println("Error: Invalid nr-of-threads: ", nrOfThreads)
	}

	if sanityChecksFailed {
		fmt.Fprint(os.Stderr, MergeOpticalDuplicatesMetricsHelp)
		os.Exit(1)
	}

	// building output command line

	var command bytes.Buffer
	fmt.Fprint(&command, os.Args[0], " merge-optical-duplicates-metrics ", input, " ", output, " ", metrics, " ", intermediateMetrics)
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
	if removeDuplicates {
		fmt.Fprint(&command, " --remove-duplicates")
	}

	// executing command

	log.Println("Executing command:\n", command.String())

	var ctr filters.DuplicatesCtrMap

	// merge intermediate metrics files
	err = timedRun(timed, profile, "Loading and combining duplicate metrics.", 1, func() error {
		ctr = filters.LoadAndCombineDuplicateMetrics(metricsDir)
		return ctr.Err()
	})
	if err != nil {
		return err
	}
	return timedRun(timed, profile, "Printing comdined duplicate metrics.", 2, func() error {
		return filters.PrintDuplicatesMetrics(input, output, metrics, removeDuplicates, ctr)
	})
}
