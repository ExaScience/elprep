package cmd

import (
	"bytes"
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"runtime"

	"github.com/exascience/elprep/internal"
	"github.com/exascience/elprep/sam"
)

// SplitHelp is the help string for this command.
const SplitHelp = "Split parameters:\n" +
	"elprep split (sam-file | /path/to/input/) /path/to/output/\n" +
	"[--output-prefix name]\n" +
	"[--output-type [sam | bam | cram]]\n" +
	"[--single-end]\n" +
	"[--nr-of-threads nr]\n" +
	"[--reference-t fai-file]\n" +
	"[--reference-T fasta-file]\n"

// Split implements the elprep split command.
func Split() error {
	var (
		outputPrefix, outputType, referenceFai, referenceFasta string
		nrOfThreads                                            int
		singleEnd                                              bool
	)

	var flags flag.FlagSet

	flags.StringVar(&outputPrefix, "output-prefix", "", "prefix for the output files")
	flags.StringVar(&outputType, "output-type", "", "format of the output files")
	flags.BoolVar(&singleEnd, "single-end", false, "when splitting single-end data")
	flags.IntVar(&nrOfThreads, "nr-of-threads", 0, "number of worker threads")
	flags.StringVar(&referenceFai, "reference-t", "", "specify a .fai file for cram output")
	flags.StringVar(&referenceFasta, "reference-T", "", "specify a .fasta file for cram output")

	if len(os.Args) < 4 {
		fmt.Fprintln(os.Stderr, "Incorrect number of parameters.")
		fmt.Fprint(os.Stderr, SplitHelp)
		os.Exit(1)
	}

	input := getFilename(os.Args[2], SplitHelp)
	output := getFilename(os.Args[3], SplitHelp)

	if err := flags.Parse(os.Args[4:]); err != nil {
		x := 0
		if err != flag.ErrHelp {
			fmt.Fprintln(os.Stderr, err)
			x = 1
		}
		fmt.Fprint(os.Stderr, SplitHelp)
		os.Exit(x)
	}

	ext := filepath.Ext(input)
	if outputPrefix == "" {
		base := filepath.Base(input)
		outputPrefix = base[:len(base)-len(ext)]
	}
	if outputType == "" {
		switch ext {
		case sam.SamExt, sam.BamExt, sam.CramExt:
			outputType = ext[1:]
		default:
			outputType = "sam"
		}
	}

	setLogOutput()

	// sanity checks

	var sanityChecksFailed bool

	referenceFai, referenceFasta, success := checkCramOutputOptions(outputType, referenceFai, referenceFasta)
	sanityChecksFailed = !success

	if filepath.Dir(output) != filepath.Clean(output) {
		log.Printf("Given output path is not a path: %v.\n", output)
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
	if referenceFai != "" {
		fmt.Fprint(&command, " --reference-t ", referenceFai)
	}
	if referenceFasta != "" {
		fmt.Fprint(&command, " --reference-T ", referenceFasta)
	}

	// executing command

	log.Println("Executing command:\n", command.String())

	fullInput, err := internal.FullPathname(input)
	if err != nil {
		return err
	}

	fullOutput, err := internal.FullPathname(output)
	if err != nil {
		return err
	}

	err = os.MkdirAll(output, 0700)
	if err != nil {
		return err
	}

	if singleEnd {
		return sam.SplitSingleEndFilePerChromosome(fullInput, fullOutput, outputPrefix, outputType, referenceFai, referenceFasta)
	}
	return sam.SplitFilePerChromosome(fullInput, fullOutput, outputPrefix, outputType, referenceFai, referenceFasta)
}
