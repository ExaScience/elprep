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
	"os/exec"
	"path"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"time"

	"github.com/exascience/elprep/v4/filters"
	"github.com/exascience/elprep/v4/internal"

	"github.com/exascience/elprep/v4/sam"
)

// SfmHelp is the help string for this command.
const SfmHelp = "\nsfm parameters:\n" +
	"elprep sfm sam-file sam-output-file\n" +
	"[--replace-reference-sequences sam-file]\n" +
	"[--filter-unmapped-reads]\n" +
	"[--filter-unmapped-reads-strict]\n" +
	"[--filter-mapping-quality mapping-quality]\n" +
	"[--filter-non-exact-mapping-reads]\n" +
	"[--filter-non-exact-mapping-reads-strict]\n" +
	"[--filter-non-overlapping-reads bed-file]\n" +
	"[--replace-read-group read-group-string]\n" +
	"[--mark-duplicates]\n" +
	"[--mark-optical-duplicates file]\n" +
	"[--optical-duplicates-pixel-distance nr]\n" +
	"[--remove-duplicates]\n" +
	"[--remove-optional-fields [all | list]]\n" +
	"[--keep-optional-fields [none | list]]\n" +
	"[--sorting-order [keep | unknown | unsorted | queryname | coordinate]]\n" +
	"[--clean-sam]\n" +
	"[--bqsr]\n" +
	"[--bqsr-reference elfasta]\n" +
	"[--quantize-levels nr]\n" +
	"[--sqq list]\n" +
	"[--known-sites list]\n" +
	"[--nr-of-threads nr]\n" +
	"[--timed]\n" +
	"[--log-path path]\n" +
	"[--intermediate-files-output-prefix name]\n" +
	"[--intermediate-files-output-type [sam | bam]]\n" +
	"[--single-end]\n" +
	"[--contig-group-size nr]\n"

	// CombinedSfmFilterHelp is a help string that combines the help strings for the filter and sfm commands.
const CombinedSfmFilterHelp = "filter/sfm parameters:\n" +
	"elprep [filter | sfm] sam-file sam-output-file\n" +
	"[--replace-reference-sequences sam-file]\n" +
	"[--filter-unmapped-reads]\n" +
	"[--filter-unmapped-reads-strict]\n" +
	"[--filter-mapping-quality mapping-quality]\n" +
	"[--filter-non-exact-mapping-reads]\n" +
	"[--filter-non-exact-mapping-reads-strict]\n" +
	"[--filter-non-overlapping-reads bed-file]\n" +
	"[--replace-read-group read-group-string]\n" +
	"[--mark-duplicates]\n" +
	"[--mark-optical-duplicates file]\n" +
	"[--optical-duplicates-pixel-distance nr]\n" +
	"[--remove-duplicates]\n" +
	"[--remove-optional-fields [all | list]]\n" +
	"[--keep-optional-fields [none | list]]\n" +
	"[--sorting-order [keep | unknown | unsorted | queryname | coordinate]]\n" +
	"[--clean-sam]\n" +
	"[--bqsr recal-file]\n" +
	"[--bqsr-reference elfasta]\n" +
	"[--quantize-levels nr]\n" +
	"[--sqq list]\n" +
	"[--known-sites list]\n" +
	"[--nr-of-threads nr]\n" +
	"[--timed]\n" +
	"[--log-path path]\n" +
	"[--intermediate-files-output-prefix name] (sfm only)\n" +
	"[--intermediate-files-output-type [sam | bam]] (sfm only)\n" +
	"[--single-end] (sfm only)\n" +
	"[--contig-group-size nr] (sfm only)\n"

// Sfm implements the elprep sfm command.
func Sfm() error {
	var (
		replaceReferenceSequences                           string
		filterUnmappedReads, filterUnmappedReadsStrict      bool
		filterMappingQuality                                int
		filterNonExactMappingReads                          bool
		filterNonExactMappingReadsStrict                    bool
		filterNonOverlappingReads                           string
		replaceReadGroup                                    string
		markDuplicates, markDuplicatesDet, removeDuplicates bool
		markOpticalDuplicates                               string
		opticalDuplicatesPixelDistance                      int
		removeOptionalFields                                string
		keepOptionalFields                                  string
		sortingOrderString                                  string
		cleanSam                                            bool
		bqsr                                                string
		referenceElFasta                                    string
		quantizeLevels                                      int
		sqq                                                 string
		knownSites                                          string
		deterministic                                       bool
		nrOfThreads                                         int
		timed                                               bool
		profile                                             string
		logPath                                             string
		renameChromosomes                                   bool
		outputPrefix, outputType                            string
		contigGroupSize                                     int
		singleEnd                                           bool
	)

	var flags flag.FlagSet

	// filter flags
	flags.StringVar(&replaceReferenceSequences, "replace-reference-sequences", "", "replace the existing header by a new one")
	flags.BoolVar(&filterUnmappedReads, "filter-unmapped-reads", false, "remove all unmapped alignments")
	flags.BoolVar(&filterUnmappedReadsStrict, "filter-unmapped-reads-strict", false, "remove all unmapped alignments, taking also POS and RNAME into account")
	flags.IntVar(&filterMappingQuality, "filter-mapping-quality", 0, "output only reads that equal or exceed given mapping quality")
	flags.BoolVar(&filterNonExactMappingReads, "filter-non-exact-mapping-reads", false, "output only exact mapping reads (soft-clipping allowed) based on cigar string (only M,S allowed)")
	flags.BoolVar(&filterNonExactMappingReadsStrict, "filter-non-exact-mapping-reads-strict", false, "output only exact mapping reads (soft-clipping allowed) based on optional fields X0=1, X1=0, XM=0, XO=0, XG=0")
	flags.StringVar(&filterNonOverlappingReads, "filter-non-overlapping-reads", "", "output only reads that overlap with the given regions (bed format)")
	flags.StringVar(&replaceReadGroup, "replace-read-group", "", "add or replace alignment read groups")
	flags.BoolVar(&markDuplicates, "mark-duplicates", false, "mark duplicates")
	flags.BoolVar(&markDuplicatesDet, "mark-duplicates-deterministic", false, "mark duplicates deterministically")
	flags.BoolVar(&removeDuplicates, "remove-duplicates", false, "remove duplicates")
	flags.StringVar(&markOpticalDuplicates, "mark-optical-duplicates", "", "mark optical duplicates")
	flags.IntVar(&opticalDuplicatesPixelDistance, "optical-duplicates-pixel-distance", 100, "pixel distance used for optical duplicate marking")
	flags.StringVar(&removeOptionalFields, "remove-optional-fields", "", "remove the given optional fields")
	flags.StringVar(&keepOptionalFields, "keep-optional-fields", "", "remove all except for the given optional fields")
	flags.StringVar(&sortingOrderString, "sorting-order", string(sam.Keep), "determine output order of alignments, one of keep, unknown, unsorted, queryname, or coordinate")
	flags.BoolVar(&cleanSam, "clean-sam", false, "clean the sam file")
	flags.StringVar(&bqsr, "bqsr", "", "base quality score recalibration")
	flags.StringVar(&referenceElFasta, "bqsr-reference", "", "reference used for base quality score recalibration (elfasta format)")
	flags.IntVar(&quantizeLevels, "quantize-levels", 0, "number of levels to be used for quantizing recalibrated base qualities (only with --bqsr)")
	flags.StringVar(&sqq, "sqq", "", "levels to be used for statically quantizing recalibrated base qualities (only with --bqsr)")
	flags.StringVar(&knownSites, "known-sites", "", "list of vcf files containing known sites for base recalibration (only with --bqsr)")
	flags.BoolVar(&deterministic, "deterministic", false, "run elprep deterministically (currently not supported in sfm mode)")
	flags.IntVar(&nrOfThreads, "nr-of-threads", 0, "number of worker threads")
	flags.BoolVar(&timed, "timed", false, "measure the runtime")
	flags.StringVar(&profile, "profile", "", "write a runtime profile to the specified file(s)")
	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")
	flags.BoolVar(&renameChromosomes, "rename-chromosomes", false, "")
	//split/merge flags
	flags.StringVar(&outputPrefix, "intermediate-files-output-prefix", "", "prefix for the output files")
	flags.StringVar(&outputType, "intermediate-files-output-type", "", "format of the output files")
	flags.IntVar(&contigGroupSize, "contig-group-size", 0, "maximum sum of reference sequence lengths for creating groups of reference sequences")
	flags.BoolVar(&singleEnd, "single-end", false, "when splitting single-end data")

	parseFlags(flags, 4, SfmHelp)

	input := getFilename(os.Args[2], SfmHelp)
	output := getFilename(os.Args[3], SfmHelp)

	setLogOutput(logPath)

	// sanity checks

	var sanityChecksFailed bool

	if !checkExist("", input) {
		sanityChecksFailed = true
	}
	if !checkCreate("", output) {
		sanityChecksFailed = true
	}

	if replaceReferenceSequences != "" && !checkExist("--replace-reference-sequences", replaceReferenceSequences) {
		sanityChecksFailed = true
	}
	if filterNonOverlappingReads != "" && !checkExist("--filter-non-overlapping-reads", filterNonOverlappingReads) {
		sanityChecksFailed = true
	}
	if markOpticalDuplicates != "" && !checkCreate("--mark-optical-duplicates", markOpticalDuplicates) {
		sanityChecksFailed = true
	}
	if bqsr != "" && !checkCreate("--bqsr", bqsr) {
		sanityChecksFailed = true
	}
	if referenceElFasta != "" && !checkExist("--bqsr-reference", referenceElFasta) {
		sanityChecksFailed = true
	}

	if knownSites != "" {
		knownSitesList := strings.Split(knownSites, ",")
		for _, k := range knownSitesList {
			if !checkExist("--known-sites", k) {
				sanityChecksFailed = true
			}
		}
	}

	if profile != "" && !checkCreate("--profile", profile) {
		sanityChecksFailed = true
	}

	if deterministic {
		sanityChecksFailed = true
		log.Println("Error: deterministic mode currently not supported in sfm")
	}

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

	if markDuplicatesDet {
		log.Println("--mark-duplicates-deterministic is deprecated; just use --mark-duplicates instead.")
		markDuplicates = true
	}

	if !checkBQSROptions(bqsr != "", false, referenceElFasta, quantizeLevels, sqq, knownSites, bqsr, "") {
		sanityChecksFailed = true
	}

	if markOpticalDuplicates != "" && !markDuplicates {
		sanityChecksFailed = true
		log.Println("Error: Cannot use --mark-optical-duplicates without also using --mark-duplicates.")
	}

	if sanityChecksFailed {
		fmt.Fprint(os.Stderr, SfmHelp)
		os.Exit(1)
	}

	// building commandline arguments and output command line

	var splitArgs, mergeArgs, filterArgs, filterArgs2 []string
	var command bytes.Buffer
	fmt.Fprint(&command, os.Args[0], " sfm ", input, " ", output)

	if filterUnmappedReadsStrict {
		fmt.Fprint(&command, " --filter-unmapped-reads-strict")
		filterArgs = append(filterArgs, "--filter-unmapped-reads-strict")
	} else if filterUnmappedReads {
		fmt.Fprint(&command, " --filter-unmapped-reads")
		filterArgs = append(filterArgs, "--filter-unmapped-reads")
	}

	if filterMappingQuality > 0 {
		fmt.Fprint(&command, " --filter-mapping-quality ", filterMappingQuality)
		filterArgs = append(filterArgs, "--filter-mapping-quality", strconv.Itoa(filterMappingQuality))
	}

	if filterNonExactMappingReads {
		fmt.Fprint(&command, " --filter-non-exact-mapping-reads")
		filterArgs = append(filterArgs, "--filter-non-exact-mapping-reads")
	}

	if filterNonExactMappingReadsStrict {
		fmt.Fprint(&command, " --filter-non-exact-mapping-reads-strict")
		filterArgs = append(filterArgs, "--filter-non-exact-mapping-reads-strict")
	}

	if filterNonOverlappingReads != "" {
		fmt.Fprint(&command, " --filter-non-overlapping-reads ", filterNonOverlappingReads)
		filterArgs = append(filterArgs, "--filter-non-overlapping-reads", filterNonOverlappingReads)
	}

	if renameChromosomes {
		fmt.Fprint(&command, " --rename-chromosomes")
		filterArgs = append(filterArgs, "--rename-chromosomes")
	}

	if cleanSam {
		fmt.Fprint(&command, " --clean-sam")
		filterArgs = append(filterArgs, "--clean-sam")
	}

	if replaceReferenceSequences != "" {
		fmt.Fprint(&command, " --replace-reference-sequences ", replaceReferenceSequences)
		filterArgs = append(filterArgs, "--replace-reference-sequences", replaceReferenceSequences)
	}

	if replaceReadGroup != "" {
		fmt.Fprint(&command, " --replace-read-group ", replaceReadGroup)
		filterArgs = append(filterArgs, "--replace-read-group", replaceReadGroup)
	}

	if markDuplicates {
		fmt.Fprint(&command, " --mark-duplicates")
		filterArgs = append(filterArgs, "--mark-duplicates")
	}

	if markOpticalDuplicates != "" {
		fmt.Fprint(&command, " --mark-optical-duplicates ", markOpticalDuplicates)
		fmt.Fprint(&command, " --optical-duplicates-pixel-distance ", opticalDuplicatesPixelDistance)
	}

	if removeDuplicates {
		fmt.Fprint(&command, " --remove-duplicates")
		filterArgs = append(filterArgs, "--remove-duplicates")
	}

	if bqsr != "" {
		fmt.Fprint(&command, " --bqsr ", bqsr)
		fmt.Fprint(&command, " --bqsr-reference ", referenceElFasta)
		fmt.Fprint(&command, " --quantize-levels ", quantizeLevels)
		filterArgs = append(filterArgs, "--bqsr-reference", referenceElFasta, "--quantize-levels", strconv.Itoa(quantizeLevels))
	}

	if bqsr != "" && sqq != "" {
		fmt.Fprint(&command, " --sqq ", sqq)
		filterArgs2 = append(filterArgs2, "--sqq", sqq)
	}

	if bqsr != "" && knownSites != "" {
		fmt.Fprint(&command, " --known-sites ", knownSites)
		filterArgs = append(filterArgs, "--known-sites", knownSites)
	}

	if removeOptionalFields != "" {
		fmt.Fprint(&command, " --remove-optional-fields \"", removeOptionalFields, "\"")
		if bqsr != "" {
			filterArgs2 = append(filterArgs2, "--remove-optional-fields", removeOptionalFields)
		} else {
			filterArgs = append(filterArgs, "--remove-optional-fields", removeOptionalFields)
		}
	}

	if keepOptionalFields != "" {
		fmt.Fprint(&command, " --keep-optional-fields \"", keepOptionalFields, "\"")
		if bqsr != "" {
			filterArgs2 = append(filterArgs2, "--keep-optional-fields", keepOptionalFields)
		} else {
			filterArgs = append(filterArgs, "--keep-optional-fields", keepOptionalFields)
		}
	}

	fmt.Fprint(&command, " --sorting-order ", sortingOrder)
	filterArgs = append(filterArgs, "--sorting-order", sortingOrderString)

	if nrOfThreads > 0 {
		runtime.GOMAXPROCS(nrOfThreads)
		fmt.Fprint(&command, " --nr-of-threads ", nrOfThreads)
		filterArgs = append(filterArgs, "--nr-of-threads", strconv.Itoa(nrOfThreads))
		splitArgs = append(splitArgs, "--nr-of-threads", strconv.Itoa(nrOfThreads))
		mergeArgs = append(mergeArgs, "--nr-of-threads", strconv.Itoa(nrOfThreads))
	}

	if timed {
		fmt.Fprint(&command, " --timed")
		filterArgs = append(filterArgs, "--timed")
		filterArgs2 = append(filterArgs2, "--timed")
		mergeArgs = append(mergeArgs, "--timed")
		splitArgs = append(splitArgs, "--timed")
	}

	if profile != "" {
		fmt.Fprint(&command, " --profile ", profile)
		filterArgs = append(filterArgs, "--profile", profile)
		filterArgs2 = append(filterArgs2, "--profile", profile)
		splitArgs = append(splitArgs, "--profile", profile)
		mergeArgs = append(mergeArgs, "--profile", profile)
	}

	if logPath != "" {
		fmt.Fprint(&command, " --log-path ", logPath)
		filterArgs = append(filterArgs, "--log-path", logPath)
		filterArgs2 = append(filterArgs2, "--log-path", logPath)
		splitArgs = append(splitArgs, "--log-path", logPath)
		mergeArgs = append(mergeArgs, "--log-path", logPath)
	}

	ext := filepath.Ext(input)
	if outputPrefix == "" {
		base := filepath.Base(input)
		outputPrefix = base[:len(base)-len(ext)]
	}
	fmt.Fprint(&command, " --intermediate-files-output-prefix ", outputPrefix)
	splitArgs = append(splitArgs, "--output-prefix", outputPrefix)

	if outputType == "" {
		outputType = ext[1:]
	}
	fmt.Fprint(&command, " --intermediate-files-output-type ", outputType)
	splitArgs = append(splitArgs, "--output-type", outputType)

	if contigGroupSize > 0 {
		contigGroupSizeString := strconv.FormatInt(int64(contigGroupSize), 10)
		splitArgs = append(splitArgs, "--contig-group-size", contigGroupSizeString)
	}

	if singleEnd {
		splitArgs = append(splitArgs, "--single-end")
		mergeArgs = append(mergeArgs, "--single-end")
	}

	commandString := command.String()

	// executing command

	log.Println("Executing command:\n", commandString)

	// split command
	timeStamp := time.Now().Format(time.RFC3339)
	splitsName := fmt.Sprintf("elprep-splits-%s", timeStamp)
	splitsDir, err := filepath.Abs(splitsName)
	if err != nil {
		return err
	}
	splitsDir = splitsDir + string(filepath.Separator)
	splitOpt := []string{"split", os.Args[2], splitsDir}
	splitArgs = append(splitOpt, splitArgs...)
	log.Println("Splitting...")
	splitCmd := exec.Command(os.Args[0], splitArgs...)
	splitCmd.Stderr = os.Stderr
	err = splitCmd.Run()
	if err != nil {
		return err
	}

	// set up directory for metrics
	metricsName := fmt.Sprintf("elprep-metrics-%s", timeStamp)
	metricsDir := ""
	if markOpticalDuplicates != "" {
		metricsDir, err = filepath.Abs(metricsName)
		if err != nil {
			return err
		}
		metricsDir = metricsDir + string(filepath.Separator)
		err = os.MkdirAll(filepath.Dir(metricsDir), 0700)
		if err != nil {
			return err
		}
	}

	// filter commands
	mergeName := fmt.Sprintf("elprep-splits-processed-%s", timeStamp) + string(filepath.Separator)
	mergeDir, err := filepath.Abs(mergeName)
	if err != nil {
		return err
	}
	mergeDir = mergeDir + string(filepath.Separator)
	splitFilesDir := splitsDir
	if !singleEnd {
		splitFilesDir = path.Join(splitsDir, "splits") + string(filepath.Separator)
	}
	files, err := internal.Directory(splitFilesDir)
	if err != nil {
		return err
	}
	log.Println("Filtering...")
	if bqsr != "" {
		// phase 1: Recalibration
		tabsDir, err := filepath.Abs("elprep-tabs-" + timeStamp)
		if err != nil {
			return err
		}
		tabsDir = tabsDir + string(filepath.Separator)
		err = os.MkdirAll(filepath.Dir(tabsDir), 0700)
		if err != nil {
			return err
		}
		for _, fileName := range files {
			inFile := path.Join(splitFilesDir, fileName)
			outFile := path.Join(mergeDir, fileName)
			tableName := fmt.Sprintf("%s.elrecal", fileName)
			tableFile := path.Join(tabsDir, tableName)
			fileFilterArgs := append(filterArgs, "--bqsr-tables-only", tableFile)
			filterOpt := []string{"filter", inFile, outFile}
			fileFilterArgs = append(filterOpt, fileFilterArgs...)
			if markOpticalDuplicates != "" {
				metricsFile := path.Join(metricsDir, fileName)
				fileFilterArgs = append(fileFilterArgs, []string{"--mark-optical-duplicates-intermediate", metricsFile}...)
				opticalDuplicatesPixelDistanceString := strconv.FormatInt(int64(opticalDuplicatesPixelDistance), 10)
				fileFilterArgs = append(fileFilterArgs, []string{"--optical-duplicates-pixel-distance", opticalDuplicatesPixelDistanceString}...)

			}
			filterCmd := exec.Command(os.Args[0], fileFilterArgs...)
			filterCmd.Stderr = os.Stderr
			err = filterCmd.Run()
			if err != nil {
				return err
			}
		}
		// spread file
		if !singleEnd {
			spreadFileName := fmt.Sprintf("%s-spread.%s", outputPrefix, outputType)
			spreadInFile := path.Join(splitsDir, spreadFileName)
			spreadOutFile := path.Join(mergeDir, spreadFileName)
			spreadTableFileName := fmt.Sprintf("%s-spread.elrecal", outputPrefix)
			spreadTableFile := path.Join(tabsDir, spreadTableFileName)
			fileFilterArgs := append(filterArgs, "--bqsr-tables-only", spreadTableFile)
			filterOpt := []string{"filter", spreadInFile, spreadOutFile}
			fileFilterArgs = append(filterOpt, fileFilterArgs...)
			if markOpticalDuplicates != "" {
				metricsFile := path.Join(metricsDir, spreadFileName)
				fileFilterArgs = append(fileFilterArgs, []string{"--mark-optical-duplicates-intermediate", metricsFile}...)
				opticalDuplicatesPixelDistanceString := strconv.FormatInt(int64(opticalDuplicatesPixelDistance), 10)
				fileFilterArgs = append(fileFilterArgs, []string{"--optical-duplicates-pixel-distance", opticalDuplicatesPixelDistanceString}...)
			}
			filterCmd := exec.Command(os.Args[0], fileFilterArgs...)
			filterCmd.Stderr = os.Stderr
			err = filterCmd.Run()
			if err != nil {
				return err
			}
		}
	} else {
		// no bqsr
		for _, fileName := range files {
			inFile := path.Join(splitFilesDir, fileName)
			outFile := path.Join(mergeDir, fileName)
			filterOpt := []string{"filter", inFile, outFile}
			fileFilterArgs := append(filterOpt, filterArgs...)
			if markOpticalDuplicates != "" {
				metricsFile := path.Join(metricsDir, fileName)
				fileFilterArgs = append(fileFilterArgs, []string{"--mark-optical-duplicates-intermediate", metricsFile}...)
			}
			filterCmd := exec.Command(os.Args[0], fileFilterArgs...)
			filterCmd.Stderr = os.Stderr
			err = filterCmd.Run()
			if err != nil {
				return err
			}
		}
		// spread file
		if !singleEnd {
			spreadFileName := fmt.Sprintf("%s-spread.%s", outputPrefix, outputType)
			spreadInFile := path.Join(splitsDir, spreadFileName)
			spreadOutFile := path.Join(mergeDir, spreadFileName)
			filterOpt := []string{"filter", spreadInFile, spreadOutFile}
			fileFilterArgs := append(filterOpt, filterArgs...)
			if markOpticalDuplicates != "" {
				metricsFile := path.Join(metricsDir, spreadFileName)
				fileFilterArgs = append(fileFilterArgs, []string{"--mark-optical-duplicates-intermediate", metricsFile}...)
				opticalDuplicatesPixelDistanceString := strconv.FormatInt(int64(opticalDuplicatesPixelDistance), 10)
				fileFilterArgs = append(fileFilterArgs, []string{"--optical-duplicates-pixel-distance", opticalDuplicatesPixelDistanceString}...)
			}
			filterCmd := exec.Command(os.Args[0], fileFilterArgs...)
			filterCmd.Stderr = os.Stderr
			err = filterCmd.Run()
			if err != nil {
				return err
			}
		}
	}
	err = os.RemoveAll(splitsDir)
	if err != nil {
		return err
	}
	log.Println("Merging...")
	// merge
	if bqsr != "" {
		mergeOpt := []string{"merge", mergeDir, "/dev/stdout"}
		mergeArgs = append(mergeOpt, mergeArgs...)
		mergeCmd := exec.Command(os.Args[0], mergeArgs...)
		outPipe, err := mergeCmd.StdoutPipe()
		if err != nil {
			return err
		}
		mergeCmd.Stderr = os.Stderr
		// phase 2: apply bqsr
		log.Println("Filtering...")
		filterOpt2 := []string{"filter", "/dev/stdin", output}
		filterArgs2 = append(filterOpt2, filterArgs2...)
		tabsDir, err := filepath.Abs("elprep-tabs-" + timeStamp)
		if err != nil {
			return err
		}
		tabsDir = tabsDir + string(filepath.Separator)
		filterArgs2 = append(filterArgs2, "--bqsr-apply", tabsDir, "--recal-file", bqsr)
		applyBqsrCommand := exec.Command(os.Args[0], filterArgs2...)
		applyBqsrCommand.Stdin = outPipe
		applyBqsrCommand.Stderr = os.Stderr
		err = mergeCmd.Start()
		if err != nil {
			return err
		}
		err = applyBqsrCommand.Start()
		if err != nil {
			return err
		}
		err = mergeCmd.Wait()
		if err != nil {
			return err
		}
		err = applyBqsrCommand.Wait()
		if err != nil {
			return err
		}
		err = os.RemoveAll(tabsDir)
		if err != nil {
			return err
		}
	} else {
		mergeOpt := []string{"merge", mergeDir, output}
		mergeArgs = append(mergeOpt, mergeArgs...)
		mergeCmd := exec.Command(os.Args[0], mergeArgs...)
		mergeCmd.Stderr = os.Stderr
		err = mergeCmd.Run()
		if err != nil {
			return err
		}
	}
	err = os.RemoveAll(mergeDir)
	if err != nil {
		return err
	}
	// merge intermediate metrics files
	if markOpticalDuplicates != "" {
		ctr := filters.LoadAndCombineDuplicateMetrics(metricsDir)
		if err := ctr.Err(); err != nil {
			return err
		}
		if err = filters.PrintDuplicatesMetrics(input, output, markOpticalDuplicates, removeDuplicates, ctr); err != nil {
			return err
		}
		if err = os.RemoveAll(metricsDir); err != nil {
			return err
		}
	}
	return nil
}
