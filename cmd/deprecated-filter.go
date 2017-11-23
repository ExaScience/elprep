package cmd

import (
	"bytes"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"

	"github.com/exascience/elprep/sam"
	"github.com/exascience/elprep/utils"
)

type commandLine []string

func (cmd *commandLine) pop() (string, bool) {
	if line := (*[]string)(cmd); len(*line) > 0 {
		entry := (*line)[0]
		*cmd = (*line)[1:]
		return entry, true
	}
	return "", false
}

func (cmd commandLine) peek() string {
	if line := []string(cmd); len(line) > 0 {
		return line[0]
	}
	return ""
}

func appendIf(slice []sam.Filter, filters ...sam.Filter) []sam.Filter {
	for _, filter := range filters {
		if filter != nil {
			slice = append(slice, filter)
		}
	}
	return slice
}

// DeprecatedFilterHelp is the help string for this command.
const DeprecatedFilterHelp = "Filter parameters: (deprecated, please use the filter command instead)\n" +
	"elprep sam-file sam-output-file\n" +
	"[--replace-reference-sequences sam-file]\n" +
	"[--filter-unmapped-reads [strict]]\n" +
	"[--replace-read-group read-group-string]\n" +
	"[--mark-duplicates [remove]]\n" +
	"[--sorting-order [keep | unknown | unsorted | queryname | coordinate]]\n" +
	"[--clean-sam]\n" +
	"[--nr-of-threads nr]\n" +
	"[--timed]\n" +
	"[--reference-t fai-file]\n" +
	"[--reference-T fasta-file]\n"

/*
DeprecatedFilter parses the command line for elprep filter in the
style of previous elprep versions (1.x and 2.x), for backwards
compatibility. This command is deprecated and will be removed at a
later stage.
*/
func DeprecatedFilter() error {
	setLogOutput()
	log.Println("Warning: Calling elprep without a command to invoke the filter functionality is depecratead. Please use the filter command instead.")
	cmdLine := commandLine(os.Args[1:])
	sortingOrder := sam.Keep
	nrOfThreads := 0
	var markDuplicates, markDuplicatesDeterministic, timed bool
	var replaceRefSeqDictFilter,
		removeUnmappedReadsFilter,
		replaceReadGroupFilter,
		markDuplicatesFilter,
		removeDuplicatesFilter,
		cleanSamFilter,
		renameChromosomesFilter sam.Filter
	var refSeqDict, filterUnmappedArg, readGroupString, referenceFai, referenceFasta, profile string
	var filenames []string
	for entry, found := cmdLine.pop(); found; entry, found = cmdLine.pop() {
		switch entry {
		case "-h":
			fmt.Fprint(os.Stderr, DeprecatedFilterHelp)
			os.Exit(0)
		case "--replace-reference-sequences":
			refSeqDict, _ = cmdLine.pop()
			r, err := sam.ReplaceReferenceSequenceDictionaryFromSamFile(refSeqDict)
			if err != nil {
				return err
			}
			replaceRefSeqDictFilter = r
		case "--filter-unmapped-reads":
			if cmdLine.peek() == "strict" {
				filterUnmappedArg, _ = cmdLine.pop()
				removeUnmappedReadsFilter = sam.FilterUnmappedReadsStrict
			} else {
				removeUnmappedReadsFilter = sam.FilterUnmappedReads
			}
		case "--replace-read-group":
			readGroupString, _ = cmdLine.pop()
			record, err := sam.ParseHeaderLineFromString(readGroupString)
			if err != nil {
				return err
			}
			replaceReadGroupFilter = sam.AddOrReplaceReadGroup(record)
		case "--mark-duplicates":
			markDuplicates = true
			for {
				if next := cmdLine.peek(); next == "remove" {
					cmdLine.pop()
					removeDuplicatesFilter = sam.FilterDuplicateReads
				} else if next == "deterministic" {
					cmdLine.pop()
					markDuplicatesDeterministic = true
				} else {
					break
				}
			}
			markDuplicatesFilter = sam.MarkDuplicates(markDuplicatesDeterministic)
		case "--sorting-order":
			if so := cmdLine.peek(); (so == "") || strings.Contains(so, "--") {
				sortingOrder = sam.Keep
			} else {
				cmdLine.pop()
				sortingOrder = so
				switch sortingOrder {
				case sam.Keep, sam.Unknown, sam.Unsorted, sam.Queryname, sam.Coordinate:
				default:
					log.Printf("Invalid sorting-order %v.\n", sortingOrder)
					fmt.Fprint(os.Stderr, DeprecatedFilterHelp)
					os.Exit(1)
				}
			}
		case "--clean-sam":
			cleanSamFilter = sam.CleanSam
		case "--nr-of-threads":
			sval, _ := cmdLine.pop()
			if val, err := strconv.ParseInt(sval, 10, 64); err != nil {
				return err
			} else if val < 0 {
				log.Printf("Invalid nr-of-threads option %v.\n", val)
			} else {
				nrOfThreads = int(val)
				runtime.GOMAXPROCS(nrOfThreads)
			}
		case "--gc-on":
			if lvl, _ := cmdLine.pop(); (lvl == "") || strings.Contains(lvl, "--") {
				// ignored
			} else if val, err := strconv.ParseInt(lvl, 10, 8); err != nil {
				return err
			} else {
				switch val {
				case 0, 1, 2:
					log.Println("Warning: The gc-on option is not supported anymore in this version of elPrep.")
				default:
					log.Printf("Invalid gc-on option %v.\n", lvl)
					fmt.Fprint(os.Stderr, DeprecatedFilterHelp)
					os.Exit(1)
				}
			}
		case "--timed":
			timed = true
		case "--profile":
			profile, _ = cmdLine.pop()
		case "--reference-t":
			if ref := cmdLine.peek(); (ref == "") || strings.Contains(ref, "--") {
				log.Println("Please provide reference file with --reference-t.")
				fmt.Fprint(os.Stderr, DeprecatedFilterHelp)
				os.Exit(1)
			} else {
				referenceFai, _ = cmdLine.pop()
			}
		case "--reference-T":
			if ref := cmdLine.peek(); (ref == "") || strings.Contains(ref, "--") {
				log.Println("Please provide reference file with --reference-T.")
				fmt.Fprint(os.Stderr, DeprecatedFilterHelp)
				os.Exit(1)
			} else {
				referenceFasta, _ = cmdLine.pop()
			}
		case "--rename-chromosomes":
			renameChromosomesFilter = sam.RenameChromosomes
		case "--split-file":
			log.Println("Warning: The split-file option is not necessary anymore in this version of elPrep.")
		default:
			filenames = append(filenames, entry)
		}
	}
	if len(filenames) != 2 {
		log.Printf("Incorrect number of parameters: Expected 2, got %v.\n", filenames)
		fmt.Fprint(os.Stderr, DeprecatedFilterHelp)
		os.Exit(1)
	}
	if (replaceRefSeqDictFilter != nil) && (sortingOrder == sam.Keep) {
		log.Print("Warning: Requesting to keep the order of the input file while replacing the reference sequence dictionary may force an additional sorting phase to ensure the original sorting order is respected.")
	}
	if (filepath.Ext(filenames[1]) == sam.CramExt) && (referenceFai == "") && (referenceFasta == "") {
		log.Println("Error: Attempting to output to cram without specifying a reference file. Please add --reference-t or --reference-T to your call.")
		fmt.Fprint(os.Stderr, DeprecatedFilterHelp)
		os.Exit(1)
	}
	var s bytes.Buffer
	fmt.Fprint(&s, os.Args[0], " filter ", filenames[0], " ", filenames[1])
	if referenceFai != "" {
		fmt.Fprint(&s, " --reference-t ", referenceFai)
	}
	if referenceFasta != "" {
		fmt.Fprint(&s, " --reference-T ", referenceFasta)
	}
	if removeUnmappedReadsFilter != nil {
		switch filterUnmappedArg {
		case "strict":
			fmt.Fprint(&s, " --filter-unmapped-reads-strict")
		case "":
			fmt.Fprint(&s, " --filter-unmapped-reads")
		default:
			log.Fatalf("Invalid argument %v to --filter-unmapped-reads.", filterUnmappedArg)
		}
	}
	if renameChromosomesFilter != nil {
		fmt.Fprint(&s, " --rename-chromosomes")
	}
	if cleanSamFilter != nil {
		fmt.Fprint(&s, " --clean-sam")
	}
	if replaceRefSeqDictFilter != nil {
		fmt.Fprint(&s, " --replace-reference-sequences ", refSeqDict)
	}
	if replaceReadGroupFilter != nil {
		fmt.Fprint(&s, " --replace-read-group ", readGroupString)
	}
	if markDuplicatesFilter != nil {
		if markDuplicatesDeterministic {
			fmt.Fprint(&s, " --mark-duplicates-deterministic")
		} else {
			fmt.Fprint(&s, " --mark-duplicates")
		}
		if removeDuplicatesFilter != nil {
			fmt.Fprint(&s, " --remove-duplicates")
		}
	}
	fmt.Fprint(&s, " --sorting-order ", sortingOrder)
	if nrOfThreads > 0 {
		fmt.Fprint(&s, " --nr-of-threads ", nrOfThreads)
	}
	if timed {
		fmt.Fprint(&s, " --timed")
	}
	if profile != "" {
		fmt.Fprint(&s, " --profile ", profile)
	}
	cmdString := s.String()
	filters := appendIf([]sam.Filter{sam.AddPGLine(utils.StringMap{
		"ID": ProgramName + " " + ProgramVersion,
		"PN": ProgramName,
		"VN": ProgramVersion,
		"DS": ProgramURL,
		"CL": cmdString,
	})}, removeUnmappedReadsFilter, renameChromosomesFilter, cleanSamFilter, replaceRefSeqDictFilter, replaceReadGroupFilter)
	if (replaceRefSeqDictFilter != nil) || (markDuplicatesFilter != nil) ||
		(sortingOrder == sam.Coordinate) || (sortingOrder == sam.Queryname) {
		filters = append(filters, sam.AddREFID)
	}
	filters = appendIf(filters, markDuplicatesFilter)
	filters = append(filters, sam.FilterOptionalReads)
	var filters2 []sam.Filter
	if removeDuplicatesFilter != nil {
		filters2 = append(filters2, removeDuplicatesFilter)
	}
	log.Println("Executing command:\n", cmdString)
	if markDuplicates || (sortingOrder == sam.Coordinate) || (sortingOrder == sam.Queryname) ||
		((replaceRefSeqDictFilter != nil) && (sortingOrder == sam.Keep)) {
		return runBestPracticesPipelineIntermediateSam(filenames[0], filenames[1], referenceFai, referenceFasta, sortingOrder, filters, filters2, timed, profile)
	}
	return runBestPracticesPipeline(filenames[0], filenames[1], referenceFai, referenceFasta, sortingOrder, filters, timed, profile)
}
