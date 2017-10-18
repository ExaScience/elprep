package cmd

import (
	"bytes"
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"runtime/pprof"
	"strconv"
	"time"

	"github.com/exascience/elprep/bed"
	"github.com/exascience/elprep/internal"
	"github.com/exascience/elprep/sam"
	"github.com/exascience/elprep/utils"
)

func timedRun(timed bool, profile, msg string, phase int64, f func() error) error {
	if profile != "" {
		file, _ := os.Create(profile + strconv.FormatInt(phase, 10) + ".prof")
		pprof.StartCPUProfile(file)
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

/*
Run the best practices pipeline. Version that uses an intermediate
slice so that sorting and mark-duplicates are supported.
*/
func runBestPracticesPipelineIntermediateSam(fileIn, fileOut, fai, fasta, sortingOrder string, filters, filters2 []sam.Filter, timed bool, profile string) error {
	filteredReads := sam.NewSam()
	err := timedRun(timed, profile, "Reading SAM into memory and applying filters.", 1, func() (err error) {
		pathname, err := internal.FullPathname(fileIn)
		if err != nil {
			return err
		}
		input, err := sam.Open(pathname, false)
		if err != nil {
			return err
		}
		defer func() {
			nerr := input.Close()
			if err == nil {
				err = nerr
			}
		}()
		return input.SamReader().RunPipeline(filteredReads, filters, sortingOrder)
	})
	if err != nil {
		return err
	}
	return timedRun(timed, profile, "Write to file.", 2, func() (err error) {
		pathname, err := internal.FullPathname(fileOut)
		if err != nil {
			return err
		}
		err = os.MkdirAll(filepath.Dir(pathname), 0700)
		if err != nil {
			return err
		}
		output, err := sam.Create(pathname, fai, fasta)
		if err != nil {
			return err
		}
		defer func() {
			nerr := output.Close()
			if err == nil {
				err = nerr
			}
		}()
		if sortingOrder != "unsorted" {
			sortingOrder = "keep"
		}
		return filteredReads.RunPipeline(output.SamWriter(), filters2, sortingOrder)
	})
}

/*
Run the best practices pipeline. Version that doesn't use an
intermediate slice when neither sorting nor mark-duplicates are
needed.
*/
func runBestPracticesPipeline(fileIn, fileOut, fai, fasta, sortingOrder string, filters []sam.Filter, timed bool, profile string) error {
	return timedRun(timed, profile, "Running pipeline.", 1, func() (err error) {
		pathname, err := internal.FullPathname(fileIn)
		if err != nil {
			return err
		}
		input, err := sam.Open(pathname, false)
		if err != nil {
			return err
		}
		defer func() {
			nerr := input.Close()
			if err == nil {
				err = nerr
			}
		}()
		pathname, err = internal.FullPathname(fileOut)
		if err != nil {
			return err
		}
		output, err := sam.Create(pathname, fai, fasta)
		if err != nil {
			return err
		}
		defer func() {
			nerr := output.Close()
			if err == nil {
				err = nerr
			}
		}()
		return input.SamReader().RunPipeline(output.SamWriter(), filters, sortingOrder)
	})
}

const FilterHelp = "Filter parameters:\n" +
	"elprep filter sam-file sam-output-file\n" +
	"[--replace-reference-sequences sam-file]\n" +
	"[--filter-unmapped-reads]\n" +
	"[--filter-unmapped-reads-strict]\n" +
	"[--filter-non-exact-mapping-reads]\n" +
	"[--filter-non-exact-mapping-reads-strict]\n" +
	"[--filter-non-overlapping-reads bed-file]\n" +
	"[--replace-read-group read-group-string]\n" +
	"[--mark-duplicates]\n" +
	"[--remove-duplicates]\n" +
	"[--sorting-order [keep | unknown | unsorted | queryname | coordinate]]\n" +
	"[--clean-sam]\n" +
	"[--nr-of-threads nr]\n" +
	"[--timed]\n" +
	"[--reference-t fai-file]\n" +
	"[--reference-T fasta-file]\n"

/*
Filter implements the elprep filter command.
*/
func Filter() error {
	var (
		replaceReferenceSequences                                     string
		filterUnmappedReads, filterUnmappedReadsStrict                bool
		filterNonExactMappingReads                                    bool
		filterNonExactMappingReadsStrict                              bool
		filterNonOverlappingReads                                     string
		replaceReadGroup                                              string
		markDuplicates, markDuplicatesDeterministic, removeDuplicates bool
		sortingOrder                                                  string
		cleanSam                                                      bool
		nrOfThreads                                                   int
		timed                                                         bool
		profile                                                       string
		reference_t, reference_T                                      string
		renameChromosomes                                             bool
	)

	var flags flag.FlagSet

	flags.StringVar(&replaceReferenceSequences, "replace-reference-sequences", "", "replace the existing header by a new one")
	flags.BoolVar(&filterUnmappedReads, "filter-unmapped-reads", false, "remove all unmapped alignments")
	flags.BoolVar(&filterUnmappedReadsStrict, "filter-unmapped-reads-strict", false, "remove all unmapped alignments, taking also POS and RNAME into account")
	flags.BoolVar(&filterNonExactMappingReads, "filter-non-exact-mapping-reads", false, "output only exact mapping reads (soft-clipping allowed) based on cigar string (only M,S allowed)")
	flags.BoolVar(&filterNonExactMappingReadsStrict, "filter-non-exact-mapping-reads-strict", false, "output only exact mapping reads (soft-clipping allowed) based on optional fields X0=1, X1=0, XM=0, XO=0, XG=0")
	flags.StringVar(&filterNonOverlappingReads, "filter-non-overlapping-reads", "", "output only reads that overlap with the given regions (bed format)")
	flags.StringVar(&replaceReadGroup, "replace-read-group", "", "add or replace alignment read groups")
	flags.BoolVar(&markDuplicates, "mark-duplicates", false, "mark duplicates")
	flags.BoolVar(&markDuplicatesDeterministic, "mark-duplicates-deterministic", false, "mark duplicates deterministically")
	flags.BoolVar(&removeDuplicates, "remove-duplicates", false, "remove duplicates")
	flags.StringVar(&sortingOrder, "sorting-order", "keep", "determine output order of alignments, one of keep, unknown, unsorted, queryname, or coordinate")
	flags.BoolVar(&cleanSam, "clean-sam", false, "clean the sam file")
	flags.IntVar(&nrOfThreads, "nr-of-threads", 0, "number of worker threads")
	flags.BoolVar(&timed, "timed", false, "measure the runtime")
	flags.StringVar(&profile, "profile", "", "write a runtime profile to the specified file(s)")
	flags.StringVar(&reference_t, "reference-t", "", "specify a .fai file for cram output")
	flags.StringVar(&reference_T, "reference-T", "", "specify a .fasta file for cram output")
	flags.BoolVar(&renameChromosomes, "rename-chromosomes", false, "")

	if len(os.Args) < 4 {
		fmt.Fprint(os.Stderr, FilterHelp)
		os.Exit(1)
	}

	input := getFilename(os.Args[2], FilterHelp)
	output := getFilename(os.Args[3], FilterHelp)

	if err := flags.Parse(os.Args[4:]); err != nil {
		x := 0
		if err != flag.ErrHelp {
			fmt.Fprintln(os.Stderr, err.Error())
			x = 1
		}
		fmt.Fprint(os.Stderr, FilterHelp)
		os.Exit(x)
	}

	setLogOutput()

	// sanity checks

	sanityChecksFailed := false

	reference_t, reference_T, success := checkCramOutputOptions(filepath.Ext(output), reference_t, reference_T)
	sanityChecksFailed = !success

	switch sortingOrder {
	case "keep", "unknown", "unsorted", "queryname", "coordinate":
	default:
		sanityChecksFailed = true
		log.Println("Error: Invalid sorting-order: ", sortingOrder)
	}

	if (replaceReferenceSequences != "") && (sortingOrder == "keep") {
		log.Println("Warning: Requesting to keep the order of the input file while replacing the reference sequence dictionary may force an additional sorting phase to ensure the original sorting order is respected.")
	}

	if nrOfThreads < 0 {
		sanityChecksFailed = true
		log.Println("Error: Invalid nr-of-threads: ", nrOfThreads)
	}

	if sanityChecksFailed {
		fmt.Fprint(os.Stderr, FilterHelp)
		os.Exit(1)
	}

	// building filters and output command line

	var command bytes.Buffer
	fmt.Fprint(&command, os.Args[0], " filter ", input, " ", output)

	if reference_t != "" {
		fmt.Fprint(&command, " --reference-t ", reference_t)
	}

	if reference_T != "" {
		fmt.Fprint(&command, " --reference-T ", reference_T)
	}

	var filters, filters2 []sam.Filter

	if filterUnmappedReadsStrict {
		filters = append(filters, sam.FilterUnmappedReadsStrict)
		fmt.Fprint(&command, " --filter-unmapped-reads-strict")
	} else if filterUnmappedReads {
		filters = append(filters, sam.FilterUnmappedReads)
		fmt.Fprint(&command, " --filter-unmapped-reads")
	}

	if filterNonExactMappingReads {
		filters = append(filters, sam.FilterNonExactMappingReads)
		fmt.Fprint(&command, " --filter-non-exact-mapping-reads")
	}

	if filterNonExactMappingReadsStrict {
		filters = append(filters, sam.FilterNonExactMappingReadsStrict)
		fmt.Fprint(&command, " --filter-non-exact-mapping-reads-strict")
	}

	if filterNonOverlappingReads != "" {
		parsedBed, err := bed.ParseBed(filterNonOverlappingReads)
		if err != nil {
			return err
		}
		filterNonOverlappingReadsFilter := sam.FilterNonOverlappingReads(parsedBed)
		filters = append(filters, filterNonOverlappingReadsFilter)
		fmt.Fprint(&command, " --filter-non-overlapping-reads ", filterNonOverlappingReads)
	}

	if renameChromosomes {
		filters = append(filters, sam.RenameChromosomes)
		fmt.Fprint(&command, " --rename-chromosomes")
	}

	if cleanSam {
		filters = append(filters, sam.CleanSam)
		fmt.Fprint(&command, " --clean-sam")
	}

	if replaceReferenceSequences != "" {
		replaceReferenceSequencesFilter, err := sam.ReplaceReferenceSequenceDictionaryFromSamFile(replaceReferenceSequences)
		if err != nil {
			return err
		}
		filters = append(filters, replaceReferenceSequencesFilter)
		fmt.Fprint(&command, " --replace-reference-sequences ", replaceReferenceSequences)
	}

	if replaceReadGroup != "" {
		record, err := sam.ParseHeaderLineFromString(replaceReadGroup)
		if err != nil {
			return err
		}
		filters = append(filters, sam.AddOrReplaceReadGroup(record))
		fmt.Fprint(&command, " --replace-read-group ", replaceReadGroup)
	}

	if (replaceReferenceSequences != "") ||
		markDuplicatesDeterministic || markDuplicates ||
		(sortingOrder == "coordinate") || (sortingOrder == "queryname") {
		filters = append(filters, sam.AddREFID)
	}

	if markDuplicatesDeterministic {
		filters = append(filters, sam.MarkDuplicates(true))
		fmt.Fprint(&command, " --mark-duplicates-deterministic")
	} else if markDuplicates {
		filters = append(filters, sam.MarkDuplicates(false))
		fmt.Fprint(&command, " --mark-duplicates")
	}

	filters = append(filters, sam.FilterOptionalReads)

	if removeDuplicates {
		filters2 = append(filters2, sam.FilterDuplicateReads)
		fmt.Fprint(&command, " --remove-duplicates")
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

	commandString := command.String()

	filters = append([]sam.Filter{sam.AddPGLine(utils.StringMap{
		"ID": ProgramName + " " + ProgramVersion,
		"PN": ProgramName,
		"VN": ProgramVersion,
		"DS": ProgramURL,
		"CL": commandString,
	})}, filters...)

	// executing command

	log.Println("Executing command:\n", commandString)

	if markDuplicatesDeterministic || markDuplicates ||
		(sortingOrder == "coordinate") || (sortingOrder == "queryname") ||
		((replaceReferenceSequences != "") && (sortingOrder == "keep")) {
		return runBestPracticesPipelineIntermediateSam(input, output, reference_t, reference_T, sortingOrder, filters, filters2, timed, profile)
	} else {
		return runBestPracticesPipeline(input, output, reference_t, reference_T, sortingOrder, filters, timed, profile)
	}
}
