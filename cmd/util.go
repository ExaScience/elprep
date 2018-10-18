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
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"
)

const (
	// ProgramName is "elprep"
	ProgramName = "elprep"
	// ProgramVersion is the version of the elprep binary
	ProgramVersion = "4.0.0"
	// ProgramURL is the repository for the elprep source code
	ProgramURL = "http://github.com/exascience/elprep"
)

// ProgramMessage is the first line printed when the elprep binary is
// called.
var ProgramMessage string

func init() {
	ProgramMessage = fmt.Sprintln("\n", ProgramName, "version", ProgramVersion, "compiled with", runtime.Version(), "- see", ProgramURL, "for more information.")
}

// HelpMessage is printed to show the --help and --help-extended flags
const HelpMessage = "Print command details:\n" +
	"[--help]\n" +
	"[--help-extended]\n"

func getFilename(s, help string) string {
	switch s {
	case "-h", "--h", "-help", "--help":
		fmt.Fprint(os.Stderr, help)
		os.Exit(0)
	default:
		if strings.HasPrefix(s, "-") || strings.HasPrefix(s, "--") {
			log.Println("Filename(s) in command line missing.")
			fmt.Fprint(os.Stderr, help)
			os.Exit(1)
		}
	}
	return s
}

func parseFlags(flags flag.FlagSet, requiredArgs int, help string) {
	if len(os.Args) < requiredArgs {
		fmt.Fprintln(os.Stderr, "Incorrect number of parameters.")
		fmt.Fprint(os.Stderr, help)
		os.Exit(1)
	}
	flags.SetOutput(ioutil.Discard)
	if err := flags.Parse(os.Args[requiredArgs:]); err != nil {
		x := 0
		if err != flag.ErrHelp {
			fmt.Fprintln(os.Stderr, err)
			x = 1
		}
		fmt.Fprint(os.Stderr, help)
		os.Exit(x)
	}
	if flags.NArg() > 0 {
		fmt.Fprintln(os.Stderr, "Cannot parse remaining parameters:", flags.Args())
		fmt.Fprint(os.Stderr, help)
		os.Exit(1)
	}
}

func logCheckFile(parameter, format string, v ...interface{}) {
	if parameter != "" {
		log.Printf(format+" for command line parameter %v.\n", append(v, parameter)...)
	} else {
		log.Printf(format+".\n", v...)
	}
}

func checkExist(parameter, filename string) bool {
	if len(filename) == 0 {
		logCheckFile(parameter, "Error: Missing filename")
		return false
	}
	if filename[0] == '-' {
		logCheckFile(parameter, "Error: Missing filename before %v", filename)
		return false
	}
	if _, err := os.Stat(filename); err == nil {
		return true
	} else if os.IsNotExist(err) {
		logCheckFile(parameter, "Error: File %v does not exist", filename)
		return false
	} else if os.IsPermission(err) {
		logCheckFile(parameter, "Error: No permission to read file %v", filename)
		return false
	} else {
		logCheckFile(parameter, "Error %v when trying to access file %v", err, filename)
		return false
	}
}

func checkCreate(parameter, filename string) bool {
	if len(filename) == 0 {
		logCheckFile(parameter, "Error: Missing filename")
		return false
	}
	if filename[0] == '-' {
		logCheckFile(parameter, "Error: Missing filename before %v", filename)
		return false
	}
	if _, err := os.Stat(filename); err == nil {
		// Assume that the file has been written by previous elPrep runs, and can be overwritten.
		return true
	}
	err := os.MkdirAll(filepath.Dir(filename), 0700)
	if err == nil {
		err = ioutil.WriteFile(filename, nil, 0666)
	}
	if err != nil {
		if os.IsPermission(err) {
			logCheckFile(parameter, "Error: No permission to create file %v", filename)
		} else {
			logCheckFile(parameter, "Error %v when trying to create file %v", err, filename)
		}
		return false
	}
	_ = os.Remove(filename)
	return true
}

func checkBQSROptions(bqsr, bqsrTablesOnly bool, elFasta string, quantizationLevel int, sqq string, knownSites, bqsrRecalFile, recalFile string) bool {
	if bqsr {
		if bqsrRecalFile == "" {
			log.Println("Error: Attempt to calculate base recalibration without specifying a log file for the recalibration tables.")
			return false
		}
		if elFasta == "" {
			log.Println("Error: Attempt to calculate base recalibration without specifying a reference file. Please add --bqsr-reference option to your call.")
			return false
		}
		if recalFile != "" {
			log.Println("Warning: The --recal-file option is set with using --bqsr file. The parameter is ignored and file is used instead.")
		}
	} else {
		if !bqsrTablesOnly {
			if quantizationLevel != 0 {
				log.Println("Warning: The --quantization-level optional flag is set without using --bqsr. This parameter is ignored because base recalibration is not requested.")
			}
			if sqq != "" {
				log.Println("Warning: The --sqq optional flag is set without using --bqsr. This parameter is ignored because base recalibration is not requested.")
			}
			if knownSites != "" {
				log.Println("Warning: The --known-sites optional flag is set without using --bqsr. This parameter is ignored because base recalibration is not requested.")
			}
		}
	}
	return true
}

func checkBQSRTablesOnlyOptions(bqsr bool, tableFile, elFasta string) bool {
	if bqsr {
		if tableFile == "" {
			log.Println("Error: Attempt to calculate base recalibration tables without specifying a table file.")
			return false
		}
		if elFasta == "" {
			log.Println("Error: Attempt to calculate base recalibration tables without specifying a reference file. Please add the --bqsr-reference option to your call.")
			return false
		}
	}
	return true
}

func checkBQSRApplyOptions(bqsr bool, path, recalFile string) bool {
	if bqsr {
		if path == "" {
			log.Println("Error: Attempt to apply base recalibration without specifying a path to the tables calculated with --bqsr-tables-only.")
			return false
		}
		if recalFile == "" {
			log.Println("Error: Attempt to apply base recalibration without specifying a log file for the recalibration tables. Please add the --recal-file option to your call.")
			return false
		}
	}
	return true
}

func createLogFilename() string {
	t := time.Now()
	zone, _ := t.Zone()
	return fmt.Sprintf("logs/elprep/elprep-%d-%02d-%02d-%02d-%02d-%02d-%v.log", t.Year(), t.Month(), t.Day(), t.Hour(), t.Minute(), t.Second(), zone)
}

var logOutput io.Writer

func setLogOutput(path string) {
	logPath := createLogFilename()
	var fullPath string
	if path == "" {
		fullPath = filepath.Join(os.Getenv("HOME"), logPath)
	} else {
		fullPath = filepath.Join(path, logPath)
	}
	if err := os.MkdirAll(filepath.Dir(fullPath), 0700); err != nil {
		log.Fatal(err, ", while trying to create directories for log file ", fullPath)
	}
	f, err := os.Create(fullPath)
	if err != nil {
		log.Fatal(err, ", while trying to create log file ", fullPath)
	}
	fmt.Fprintln(f, ProgramMessage)
	logOutput = io.MultiWriter(f, os.Stderr)
	log.SetOutput(logOutput)
	log.Println("Created log file at", fullPath)
	log.Println("Command line:", os.Args)
}

func timedRun(timed bool, profile, msg string, phase int64, f func() error) error {
	if profile != "" {
		filename := profile + strconv.FormatInt(phase, 10) + ".prof"
		file, err := os.Create(filename)
		if err != nil {
			return fmt.Errorf("%v, while creating file %v for a CPU profile", err, filename)
		}
		if err = pprof.StartCPUProfile(file); err != nil {
			return fmt.Errorf("%v, while starting a CPU profile", err)
		}
		defer pprof.StopCPUProfile()
	}
	if timed {
		log.Println(msg)
		start := time.Now()
		defer func() {
			end := time.Now()
			log.Println("Elapsed time: ", end.Sub(start))
		}()
	}
	return f()
}
