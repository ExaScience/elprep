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

	"github.com/exascience/elprep/internal"
	"github.com/exascience/elprep/sam"
)

// MergeHelp is the help string for this command.
const MergeHelp = "Merge parameters:\n" +
	"elprep merge /path/to/input sam-output-file\n" +
	"[--single-end]\n" +
	"[--nr-of-threads n]\n" +
	"[--reference-t fai-file]\n" +
	"[--reference-T fasta-file]\n"

// Merge implements the elprep merge command.
func Merge() error {
	var (
		referenceFai, referenceFasta string
		nrOfThreads                  int
		singleEnd                    bool
	)

	var flags flag.FlagSet

	flags.BoolVar(&singleEnd, "single-end", false, "when splitting single-end data")
	flags.IntVar(&nrOfThreads, "nr-of-threads", 0, "number of worker threads")
	flags.StringVar(&referenceFai, "reference-t", "", "specify a .fai file for cram output")
	flags.StringVar(&referenceFasta, "reference-T", "", "specify a .fasta file for cram output")

	if len(os.Args) < 4 {
		fmt.Fprintln(os.Stderr, "Incorrect number of parameters.")
		fmt.Fprint(os.Stderr, MergeHelp)
		os.Exit(1)
	}

	input := getFilename(os.Args[2], MergeHelp)
	output := getFilename(os.Args[3], MergeHelp)

	if err := flags.Parse(os.Args[4:]); err != nil {
		x := 0
		if err != flag.ErrHelp {
			fmt.Fprintln(os.Stderr, err.Error())
			x = 1
		}
		fmt.Fprint(os.Stderr, MergeHelp)
		os.Exit(x)
	}

	setLogOutput()

	// sanity checks

	var sanityChecksFailed bool

	referenceFai, referenceFasta, success := checkCramOutputOptions(filepath.Ext(output), referenceFai, referenceFasta)
	sanityChecksFailed = !success

	fullInputPath, err := internal.FullPathname(input)
	if err != nil {
		return err
	}
	filesToMerge, err := internal.Directory(fullInputPath)
	if err != nil {
		log.Printf("Given directory %v causes error %v.\n", input, err)
		sanityChecksFailed = true
	} else if len(filesToMerge) < 2 && !singleEnd {
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

	firstFile, err := sam.Open(filepath.Join(fullInputPath, filesToMerge[0]), true)
	if err != nil {
		return err
	}
	header, _, err := sam.ParseHeader(firstFile.Reader)
	if err != nil {
		return err
	}
	err = firstFile.Close()
	if err != nil {
		return err
	}

	var inputPrefix string
	for _, file := range filesToMerge {
		if idx := strings.LastIndex(file, "-unmapped"); idx >= 0 {
			inputPrefix = file[:idx]
			break
		}
	}

	var inputExtension string
	switch ext := filepath.Ext(filesToMerge[0]); ext {
	case sam.SamExt, sam.BamExt, sam.CramExt:
		inputExtension = ext[1:]
	default:
		inputExtension = "sam"
	}

	// building output command line

	var command bytes.Buffer
	fmt.Fprint(&command, os.Args[0], " merge ", input, " ", output)
	if singleEnd {
		fmt.Fprint(&command, " --single-end ")
	}
	if nrOfThreads > 0 {
		runtime.GOMAXPROCS(nrOfThreads)
		fmt.Fprint(&command, " --nr-of-threads ", nrOfThreads)
	}
	if referenceFai != "" {
		fmt.Fprint(&command, " --reference-t ", referenceFai)
	}
	if referenceFasta != "" {
		fmt.Fprint(&command, " --reference-T ", referenceFasta)
	}

	// executing command

	log.Println("Executing command:\n", command.String())

	output, err = internal.FullPathname(output)
	if err != nil {
		return err
	}

	switch header.HDSO() {
	case sam.Coordinate:
		if singleEnd {
			return sam.MergeSingleEndFilesSplitPerChromosome(fullInputPath, output, referenceFai, referenceFasta, inputPrefix, inputExtension, header)
		}
		return sam.MergeSortedFilesSplitPerChromosome(fullInputPath, output, referenceFai, referenceFasta, inputPrefix, inputExtension, header)
	case sam.Queryname:
		log.Fatal("Merging of files sorted by queryname not yet implemented.")
		panic("Unreachable code.")
	default:
		if singleEnd {
			return sam.MergeSingleEndFilesSplitPerChromosome(fullInputPath, output, referenceFai, referenceFasta, inputPrefix, inputExtension, header)
		}
		return sam.MergeUnsortedFilesSplitPerChromosome(fullInputPath, output, referenceFai, referenceFasta, inputPrefix, inputExtension, header)
	}
}
