package cmd

import (
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"
	"time"
)

const (
	// ProgramName is "elprep"
	ProgramName = "elprep"
	// ProgramVersion is the version of the elprep binary
	ProgramVersion = "3.06"
	// ProgramURL is the repository for the elprep source code
	ProgramURL = "http://github.com/exascience/elprep"
)

// ProgramMessage is the first line printed when the elprep binary is
// called.
var ProgramMessage string

func init() {
	ProgramMessage = fmt.Sprintln(ProgramName, "version", ProgramVersion, "- see", ProgramURL, "for more information.")
}

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

func checkCramOutputOptions(extension, fai, fasta string) (string, string, bool) {
	switch extension {
	case "cram", ".cram":
		if (fai == "") && (fasta == "") {
			log.Println("Error: Attempt to output to cram without specifying a reference file. Please add --reference-t or --reference-T to your call.")
			return "", "", false
		}
		if (fai != "") && (fasta != "") {
			log.Println("Error: Attempt to output to cram while specifying both a .fai and a .fasta file as a reference. Please add only one of --reference-t or --reference-T to your call, but not both.")
			return "", "", false
		}
		return fai, fasta, true
	}
	if fai != "" {
		log.Println("Warning: The .fai file specified by the --reference-t parameter is ignored because the output is not to cram.")
	}
	if fasta != "" {
		log.Println("Warning. The .fasta file specified by the --reference-T parameter is ignored because the output is not to cram.")
	}
	return "", "", true
}

func createLogFilename() string {
	t := time.Now()
	zone, _ := t.Zone()
	return fmt.Sprintf("logs/elprep/elprep-%d-%02d-%02d-%02d-%02d-%02d-%v.log", t.Year(), t.Month(), t.Day(), t.Hour(), t.Minute(), t.Second(), zone)
}

func setLogOutput() {
	logPath := createLogFilename()
	fullPath := filepath.Join(os.Getenv("HOME"), logPath)
	if err := os.MkdirAll(filepath.Dir(fullPath), 0700); err != nil {
		log.Fatal(err, ", while trying to create directories for log file ", fullPath)
	}
	f, err := os.Create(fullPath)
	if err != nil {
		log.Fatal(err, ", while trying to create log file ", fullPath)
	}
	fmt.Fprintln(f, ProgramMessage)
	log.SetOutput(io.MultiWriter(f, os.Stderr))
	log.Println("Created log file at", fullPath)
	log.Println("Command line:", os.Args)
}
