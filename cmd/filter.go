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
	"strings"
	"time"

	"github.com/exascience/elprep/bed"
	"github.com/exascience/elprep/filters"
	"github.com/exascience/elprep/internal"
	"github.com/exascience/elprep/sam"
	"github.com/exascience/elprep/utils"
)

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

// Run the best practices pipeline. Version that uses an intermediate
// slice so that sorting and mark-duplicates are supported.
func runBestPracticesPipelineIntermediateSam(fileIn, fileOut, fai, fasta string, sortingOrder sam.SortingOrder, filters, filters2 []sam.Filter, timed bool, profile string) error {
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
		if sortingOrder != sam.Unsorted {
			sortingOrder = sam.Keep
		}
		return filteredReads.RunPipeline(output.SamWriter(), filters2, sortingOrder)
	})
}

// Run the best practices pipeline. Version that doesn't use an
// intermediate slice when neither sorting nor mark-duplicates are
// needed.
func runBestPracticesPipeline(fileIn, fileOut, fai, fasta string, sortingOrder sam.SortingOrder, filters []sam.Filter, timed bool, profile string) error {
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

// FilterHelp is the help string for this command.
const FilterHelp = "Filter parameters:\n" +
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
	"[--remove-duplicates]\n" +
	"[--remove-optional-fields [all | list]]\n" +
	"[--keep-optional-fields [none | list]]\n" +
	"[--sorting-order [keep | unknown | unsorted | queryname | coordinate]]\n" +
	"[--clean-sam]\n" +
	"[--nr-of-threads nr]\n" +
	"[--timed]\n" +
	"[--reference-t fai-file]\n" +
	"[--reference-T fasta-file]\n"

// Filter implements the elprep filter command.
func Filter() error {
	var (
		replaceReferenceSequences                                     string
		filterUnmappedReads, filterUnmappedReadsStrict                bool
		filterMappingQuality                                          int
		filterNonExactMappingReads                                    bool
		filterNonExactMappingReadsStrict                              bool
		filterNonOverlappingReads                                     string
		replaceReadGroup                                              string
		markDuplicates, markDuplicatesDeterministic, removeDuplicates bool
		removeOptionalFields                                          string
		keepOptionalFields                                            string
		sortingOrderString                                            string
		cleanSam                                                      bool
		nrOfThreads                                                   int
		timed                                                         bool
		profile                                                       string
		referenceFai, referenceFasta                                  string
		renameChromosomes                                             bool
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
	flags.BoolVar(&markDuplicatesDeterministic, "mark-duplicates-deterministic", false, "mark duplicates deterministically")
	flags.BoolVar(&removeDuplicates, "remove-duplicates", false, "remove duplicates")
	flags.StringVar(&removeOptionalFields, "remove-optional-fields", "", "remove the given optional fields")
	flags.StringVar(&keepOptionalFields, "keep-optional-fields", "", "remove all except for the given optional fields")
	flags.StringVar(&sortingOrderString, "sorting-order", string(sam.Keep), "determine output order of alignments, one of keep, unknown, unsorted, queryname, or coordinate")
	flags.BoolVar(&cleanSam, "clean-sam", false, "clean the sam file")
	flags.IntVar(&nrOfThreads, "nr-of-threads", 0, "number of worker threads")
	flags.BoolVar(&timed, "timed", false, "measure the runtime")
	flags.StringVar(&profile, "profile", "", "write a runtime profile to the specified file(s)")
	flags.StringVar(&referenceFai, "reference-t", "", "specify a .fai file for cram output")
	flags.StringVar(&referenceFasta, "reference-T", "", "specify a .fasta file for cram output")
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
			fmt.Fprintln(os.Stderr, err)
			x = 1
		}
		fmt.Fprint(os.Stderr, FilterHelp)
		os.Exit(x)
	}

	setLogOutput()

	// sanity checks

	var sanityChecksFailed bool

	referenceFai, referenceFasta, success := checkCramOutputOptions(filepath.Ext(output), referenceFai, referenceFasta)
	sanityChecksFailed = !success

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

	if sanityChecksFailed {
		fmt.Fprint(os.Stderr, FilterHelp)
		os.Exit(1)
	}

	// building filters and output command line

	var command bytes.Buffer
	fmt.Fprint(&command, os.Args[0], " filter ", input, " ", output)

	if referenceFai != "" {
		fmt.Fprint(&command, " --reference-t ", referenceFai)
	}

	if referenceFasta != "" {
		fmt.Fprint(&command, " --reference-T ", referenceFasta)
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

	if (replaceReferenceSequences != "") ||
		markDuplicatesDeterministic || markDuplicates ||
		(sortingOrder == sam.Coordinate) || (sortingOrder == sam.Queryname) {
		filters1 = append(filters1, filters.AddREFID)
	}

	if markDuplicatesDeterministic {
		filters1 = append(filters1, filters.MarkDuplicates(true))
		fmt.Fprint(&command, " --mark-duplicates-deterministic")
	} else if markDuplicates {
		filters1 = append(filters1, filters.MarkDuplicates(false))
		fmt.Fprint(&command, " --mark-duplicates")
	}

	filters1 = append(filters1, filters.RemoveOptionalReads)

	if removeDuplicates {
		filters2 = append(filters2, filters.RemoveDuplicateReads)
		fmt.Fprint(&command, " --remove-duplicates")
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

	if markDuplicatesDeterministic || markDuplicates ||
		(sortingOrder == sam.Coordinate) || (sortingOrder == sam.Queryname) ||
		((replaceReferenceSequences != "") && (sortingOrder == sam.Keep)) {
		return runBestPracticesPipelineIntermediateSam(input, output, referenceFai, referenceFasta, sortingOrder, filters1, filters2, timed, profile)
	}
	return runBestPracticesPipeline(input, output, referenceFai, referenceFasta, sortingOrder, append(filters1, filters2...), timed, profile)
}
