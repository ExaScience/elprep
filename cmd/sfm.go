// elPrep: a high-performance tool for analyzing SAM/BAM files.
// Copyright (c) 2017-2020 imec vzw.

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
	"io/ioutil"
	"log"
	"math"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"

	"github.com/exascience/elprep/v5/filters"
	"github.com/exascience/elprep/v5/internal"
	"github.com/google/uuid"

	"github.com/exascience/elprep/v5/sam"
)

// SfmHelp is the help string for this command.
const SfmHelp = "\nsfm parameters:\n" +
	"elprep sfm sam-file sam-output-file\n" +
	"[--output-type [sam | bam]]\n" +
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
	"[--reference elfasta]\n" +
	"[--quantize-levels nr]\n" +
	"[--sqq list]\n" +
	"[--max-cycle nr]\n" +
	"[--known-sites list]\n" +
	"[--haplotypecaller vcf-file]\n" +
	"[--reference-confidence [GVCF | BP_RESOLUTION | NONE]\n" +
	"[--sample-name sample-name]\n" +
	"[--activity-profile igv-file]\n" +
	"[--assembly-regions igv-file]\n" +
	"[--assembly-region-padding nr]\n" +
	"[--target-regions bed-file]\n" +
	"[--nr-of-threads nr]\n" +
	"[--timed]\n" +
	"[--log-path path]\n" +
	"[--intermediate-files-output-prefix name]\n" +
	"[--intermediate-files-output-type [sam | bam]]\n" +
	"[--tmp-path path]\n" +
	"[--single-end]\n" +
	"[--contig-group-size nr]\n"

// CombinedSfmFilterHelp is a help string that combines the help strings for the filter and sfm commands.
const CombinedSfmFilterHelp = "filter/sfm parameters:\n" +
	"elprep [filter | sfm] sam-file sam-output-file\n" +
	"[--output-type [sam | bam]]\n" +
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
	"[--reference elfasta]\n" +
	"[--quantize-levels nr]\n" +
	"[--sqq list]\n" +
	"[--max-cycle nr]\n" +
	"[--known-sites list]\n" +
	"[--haplotypecaller vcf-file]\n" +
	"[--reference-confidence [GVCF | BP_RESOLUTION | NONE]\n" +
	"[--sample-name sample-name]\n" +
	"[--activity-profile igv-file]\n" +
	"[--assembly-regions igv-file]\n" +
	"[--assembly-region-padding nr]\n" +
	"[--target-regions bed-file]\n" +
	"[--nr-of-threads nr]\n" +
	"[--timed]\n" +
	"[--log-path path]\n" +
	"[--intermediate-files-output-prefix name] (sfm only)\n" +
	"[--intermediate-files-output-type [sam | bam]] (sfm only)\n" +
	"[--tmp-path path]\n" +
	"[--single-end] (sfm only)\n" +
	"[--contig-group-size nr] (sfm only)\n"

// Sfm implements the elprep sfm command.
func Sfm() {
	var (
		outputType                                       string
		replaceReferenceSequences                        string
		filterUnmappedReads, filterUnmappedReadsStrict   bool
		filterMappingQuality                             int
		filterNonExactMappingReads                       bool
		filterNonExactMappingReadsStrict                 bool
		filterNonOverlappingReads                        string
		replaceReadGroup                                 string
		markDuplicates, removeDuplicates                 bool
		markOpticalDuplicates                            string
		opticalDuplicatesPixelDistance                   int
		removeOptionalFields                             string
		keepOptionalFields                               string
		sortingOrderString                               string
		cleanSam                                         bool
		bqsr                                             string
		referenceElFasta                                 string
		quantizeLevels                                   int
		sqq                                              string
		maxCycle                                         int
		knownSites                                       string
		vcfOutput                                        string
		assemblyRegionPadding                            int
		referenceConfidence                              string
		sampleName                                       string
		activityProfile, assemblyRegions                 string
		targetRegions                                    string
		nrOfThreads                                      int
		timed                                            bool
		profile                                          string
		logPath, tmpPath                                 string
		renameChromosomes                                bool
		intermediateOutputPrefix, intermediateOutputType string
		contigGroupSize                                  int
		singleEnd                                        bool
		clearDuplicateFlag                               bool
	)

	var flags flag.FlagSet

	// filter flags
	flags.StringVar(&outputType, "output-type", "", "format of the final output file")
	flags.StringVar(&replaceReferenceSequences, "replace-reference-sequences", "", "replace the existing header by a new one")
	flags.BoolVar(&filterUnmappedReads, "filter-unmapped-reads", false, "remove all unmapped alignments")
	flags.BoolVar(&filterUnmappedReadsStrict, "filter-unmapped-reads-strict", false, "remove all unmapped alignments, taking also POS and RNAME into account")
	flags.IntVar(&filterMappingQuality, "filter-mapping-quality", 0, "output only reads that equal or exceed given mapping quality")
	flags.BoolVar(&filterNonExactMappingReads, "filter-non-exact-mapping-reads", false, "output only exact mapping reads (soft-clipping allowed) based on cigar string (only M,S allowed)")
	flags.BoolVar(&filterNonExactMappingReadsStrict, "filter-non-exact-mapping-reads-strict", false, "output only exact mapping reads (soft-clipping allowed) based on optional fields X0=1, X1=0, XM=0, XO=0, XG=0")
	flags.StringVar(&filterNonOverlappingReads, "filter-non-overlapping-reads", "", "output only reads that overlap with the given regions (bed format)")
	flags.StringVar(&replaceReadGroup, "replace-read-group", "", "add or replace alignment read groups")
	flags.BoolVar(&markDuplicates, "mark-duplicates", false, "mark duplicates")
	flags.BoolVar(&removeDuplicates, "remove-duplicates", false, "remove duplicates")
	flags.StringVar(&markOpticalDuplicates, "mark-optical-duplicates", "", "mark optical duplicates")
	flags.IntVar(&opticalDuplicatesPixelDistance, "optical-duplicates-pixel-distance", 100, "pixel distance used for optical duplicate marking")
	flags.StringVar(&removeOptionalFields, "remove-optional-fields", "", "remove the given optional fields")
	flags.StringVar(&keepOptionalFields, "keep-optional-fields", "", "remove all except for the given optional fields")
	flags.StringVar(&sortingOrderString, "sorting-order", string(sam.Keep), "determine output order of alignments, one of keep, unknown, unsorted, queryname, or coordinate")
	flags.BoolVar(&cleanSam, "clean-sam", false, "clean the sam file")
	flags.StringVar(&bqsr, "bqsr", "", "base quality score recalibration")
	flags.StringVar(&referenceElFasta, "reference", "", "reference used for base quality score recalibration and haplotypecaller (elfasta format)")
	flags.IntVar(&quantizeLevels, "quantize-levels", 0, "number of levels to be used for quantizing recalibrated base qualities (only with --bqsr)")
	flags.StringVar(&sqq, "sqq", "", "levels to be used for statically quantizing recalibrated base qualities (only with --bqsr)")
	flags.IntVar(&maxCycle, "max-cycle", 500, "maximum cycle length (only with --bqsr)")
	flags.StringVar(&knownSites, "known-sites", "", "list of vcf files containing known sites for base recalibration (only with --bqsr)")
	flags.StringVar(&vcfOutput, "haplotypecaller", "", "invoke the haplotypecaller")
	flags.IntVar(&assemblyRegionPadding, "assembly-region-padding", 100, "padding around assembly regions during variant calling (only with --haplotypecaller)")
	flags.StringVar(&referenceConfidence, "reference-confidence", "GVCF", "reference confidence mode (GVCF, BP_RESOLUTION, NONE) (only with --haplotypecaller)")
	flags.StringVar(&sampleName, "sample-name", "", "sample name to use in haplotypecaller (only with --haplotypecaller)")
	flags.StringVar(&activityProfile, "activity-profile", "", "IGV file to write the activity profile to (only with --haplotypecaller)")
	flags.StringVar(&assemblyRegions, "assembly-regions", "", "IGV file to write the assembly regions to (only with --haplotypecaller)")
	flags.StringVar(&targetRegions, "target-regions", "", "bed file specifying which regions in the genome to look at (only with --bqsr and --haplotypecaller)")
	flags.IntVar(&nrOfThreads, "nr-of-threads", 0, "number of worker threads")
	flags.BoolVar(&timed, "timed", false, "measure the runtime")
	flags.StringVar(&profile, "profile", "", "write a runtime profile to the specified file(s)")
	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")
	flags.StringVar(&tmpPath, "tmp-path", "", "write split files to a specified directory")
	flags.BoolVar(&renameChromosomes, "rename-chromosomes", false, "")
	//split/merge flags
	flags.StringVar(&intermediateOutputPrefix, "intermediate-files-output-prefix", "", "prefix for the intermediate output files")
	flags.StringVar(&intermediateOutputType, "intermediate-files-output-type", "", "format of the intermediate output files")
	flags.IntVar(&contigGroupSize, "contig-group-size", 0, "maximum sum of reference sequence lengths for creating groups of reference sequences")
	flags.BoolVar(&singleEnd, "single-end", false, "when splitting single-end data")
	flags.BoolVar(&clearDuplicateFlag, "clear-duplicate-flag", false, "clear the duplicate flag in every read")

	flags.StringVar(&internal.BQSRTablenamePrefix, "bqsr-tablename-prefix", "GATK", "prefix to be used in BQSR table reports")

	parseFlags(flags, 4, SfmHelp)

	input := getFilename(os.Args[2], SfmHelp)
	output := getFilename(os.Args[3], SfmHelp)

	setLogOutput(logPath)

	// sanity checks

	var sanityChecksFailed bool

	if !checkOutputFormat(outputType) {
		sanityChecksFailed = true
	}

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
	if referenceElFasta != "" && !checkExist("--reference", referenceElFasta) {
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

	if vcfOutput == "" {
		switch output {
		case "/dev/null", "/dev/zero":
			sanityChecksFailed = true
			log.Println("Error: Neither writing a SAM output, nor a VCF output.")
		}
	} else if !checkCreate("--haplotypecaller", vcfOutput) {
		sanityChecksFailed = true
	}

	if activityProfile != "" && !checkCreate("--activity-profile", activityProfile) {
		sanityChecksFailed = true
	}

	if assemblyRegions != "" && !checkCreate("--assembly-regions", assemblyRegions) {
		sanityChecksFailed = true
	}

	referenceConfidence = strings.ToUpper(referenceConfidence)
	switch referenceConfidence {
	case "GVCF", "BP_RESOLUTION", "NONE":
	default:
		sanityChecksFailed = true
		log.Println("Error: --reference-confidence must be one of GVCF, BP_RESOLUTION, or NONE")
	}

	if assemblyRegionPadding < 0 {
		sanityChecksFailed = true
		log.Println("Error: --assembly-region-padding must be >= 0")
	} else if assemblyRegionPadding > math.MaxInt32 {
		sanityChecksFailed = true
		log.Println("Error: --assembly-region-padding to large")
	}

	if targetRegions != "" {
		if bqsr == "" && vcfOutput == "" {
			sanityChecksFailed = true
			log.Println("Error: --target-regions used without --bqsr or --haplotypecaller")
		}
	}

	if targetRegions != "" && filterNonOverlappingReads != "" {
		if targetRegions != filterNonOverlappingReads {
			log.Println("Warning: using --filter-non-overlapping-reads and --target-regions with different .bed files.")
		}
	}

	if filterNonOverlappingReads != "" && vcfOutput == "" && targetRegions == "" {
		log.Println("Warning: using --filter-non-overlapping-reads and --haplotypecaller, but not --target-regions. Haplotypecaller will run for full genome.")
	}

	if profile != "" && !checkCreate("--profile", profile) {
		sanityChecksFailed = true
	}

	if singleEnd && vcfOutput != "" {
		sanityChecksFailed = true
		log.Println("Error: Haplotypecaller (sfm) not supported on single end reads.")
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

	if bqsr != "" {
		if !checkBQSROptions(referenceElFasta, bqsr, "") {
			sanityChecksFailed = true
		}
	} else if !checkNonBQSROptions(quantizeLevels, sqq, knownSites) {
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

	if outputType != "" {
		fmt.Fprint(&command, " --output-type ", outputType)
		if bqsr != "" {
			filterArgs2 = append(filterArgs2, "--output-type", outputType)
		} else {
			mergeArgs = append(mergeArgs, "--output-type", outputType)
		}
	}

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

	if clearDuplicateFlag {
		fmt.Fprint(&command, " --clear-duplicate-flag")
		filterArgs = append(filterArgs, "--clear-duplicate-flag")
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
		fmt.Fprint(&command, " --reference ", referenceElFasta)
		fmt.Fprint(&command, " --quantize-levels ", quantizeLevels)
		fmt.Fprint(&command, " --max-cycle ", maxCycle)
		maxCycleString := strconv.Itoa(maxCycle)
		filterArgs = append(filterArgs, "--reference", referenceElFasta, "--max-cycle", maxCycleString)
		filterArgs2 = append(filterArgs2, "--quantize-levels", strconv.Itoa(quantizeLevels), "--max-cycle", maxCycleString)
	}

	if bqsr != "" && sqq != "" {
		fmt.Fprint(&command, " --sqq ", sqq)
		filterArgs2 = append(filterArgs2, "--sqq", sqq)
	}

	if bqsr != "" && knownSites != "" {
		fmt.Fprint(&command, " --known-sites ", knownSites)
		filterArgs = append(filterArgs, "--known-sites", knownSites)
	}

	if targetRegions != "" {
		fmt.Fprint(&command, " --target-regions ", targetRegions)
		// parameter added later to filter args because spread file is treated special
	}

	if vcfOutput != "" {
		fmt.Fprint(&command, " --haplotypecaller ", vcfOutput)
		if bqsr != "" {
			filterArgs2 = append(filterArgs2, "--reference", referenceElFasta)
			if sampleName != "" {
				filterArgs2 = append(filterArgs2, "--sample-name", sampleName)
			}
			if activityProfile != "" {
				filterArgs2 = append(filterArgs2, "--activity-profile", activityProfile)
			}
			if assemblyRegions != "" {
				filterArgs2 = append(filterArgs2, "--assembly-regions", assemblyRegions)
			}
			if referenceConfidence != "" {
				filterArgs2 = append(filterArgs2, "--reference-confidence", referenceConfidence)
			}
			// --target-regions parameter added later because spread file is treated special
			filterArgs2 = append(filterArgs2, "--assembly-region-padding", strconv.Itoa(assemblyRegionPadding))
		} else {
			filterArgs = append(filterArgs, "--reference", referenceElFasta)
			if sampleName != "" {
				filterArgs = append(filterArgs, "--sample-name", sampleName)
			}
			if activityProfile != "" {
				filterArgs = append(filterArgs, "--activity-profile", activityProfile)
			}
			if assemblyRegions != "" {
				filterArgs = append(filterArgs, "--assembly-regions", assemblyRegions)
			}
			if referenceConfidence != "" {
				filterArgs = append(filterArgs, "--reference-confidence", referenceConfidence)
			}
			filterArgs = append(filterArgs, "--assembly-region-padding", strconv.Itoa(assemblyRegionPadding))
		}
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

	if tmpPath != "" {
		fmt.Fprint(&command, " --tmp-path ", tmpPath)
		if err := os.MkdirAll(tmpPath, 0700); err != nil {
			log.Panic(err, ", while trying to create directories for split files ", tmpPath)
		}
	}

	ext := filepath.Ext(input)
	if intermediateOutputPrefix == "" {
		base := filepath.Base(input)
		intermediateOutputPrefix = base[:len(base)-len(ext)]
	}
	fmt.Fprint(&command, " --intermediate-files-output-prefix ", intermediateOutputPrefix)
	splitArgs = append(splitArgs, "--output-prefix", intermediateOutputPrefix)

	if intermediateOutputType == "" {
		if ext == "" {
			intermediateOutputType = "sam"
		} else {
			intermediateOutputType = ext[1:]
		}
	}
	fmt.Fprint(&command, " --intermediate-files-output-type ", intermediateOutputType)
	splitArgs = append(splitArgs, "--output-type", intermediateOutputType)

	if contigGroupSize > 0 {
		contigGroupSizeString := strconv.FormatInt(int64(contigGroupSize), 10)
		splitArgs = append(splitArgs, "--contig-group-size", contigGroupSizeString)
	}

	if singleEnd {
		splitArgs = append(splitArgs, "--single-end")
		mergeArgs = append(mergeArgs, "--single-end")
	}

	commandString := command.String()

	filterArgs = append(filterArgs, "--pg-cmd-line", commandString)

	// executing command

	log.Println("Executing command:\n", commandString)

	// split command
	uuidStamp := uuid.New().String()
	splitsName := filepath.Join(tmpPath, fmt.Sprintf("elprep-splits-%s", uuidStamp))
	splitsDir := internal.FilepathAbs(splitsName) + string(filepath.Separator)
	splitOpt := []string{"split", input, splitsDir}
	splitArgs = append(splitOpt, splitArgs...)
	log.Println("Splitting...")
	splitCmd := exec.Command(os.Args[0], splitArgs...)
	if input == "/dev/stdin" {
		splitCmd.Stdin = os.Stdin
	}
	splitCmd.Stderr = os.Stderr
	internal.RunCmd(splitCmd)

	// set up directory for metrics
	metricsDir := ""
	if markOpticalDuplicates != "" {
		metricsName := filepath.Join(tmpPath, fmt.Sprintf("elprep-metrics-%s", uuidStamp))
		metricsDir = internal.FilepathAbs(metricsName) + string(filepath.Separator)
		internal.MkdirAll(filepath.Dir(metricsDir), 0700)
	}

	// Set up directory for vcf files
	vcfDir := ""
	randomSeedFile := ""
	if vcfOutput != "" {
		vcfName := filepath.Join(tmpPath, fmt.Sprintf("elprep-vcfs-%s", uuidStamp))
		vcfDir = internal.FilepathAbs(vcfName) + string(filepath.Separator)
		internal.MkdirAll(filepath.Dir(vcfDir), 0700)
		if internal.PedanticMode {
			if rsFile, err := ioutil.TempFile("", "random-seed"); err != nil {
				log.Panicf("Cannot created temporary file for random seeds: %v", err)
			} else {
				randomSeedFile = rsFile.Name()
				if _, err := fmt.Fprintln(rsFile, "init"); err != nil {
					log.Panic(err)
				}
				internal.Close(rsFile)
			}
		}
	}

	// filter commands
	mergeName := filepath.Join(tmpPath, fmt.Sprintf("elprep-splits-processed-%s", uuidStamp))
	mergeDir := internal.FilepathAbs(mergeName) + string(filepath.Separator)
	splitFilesDir := splitsDir
	if !singleEnd {
		splitFilesDir = path.Join(splitsDir, "splits")
	}
	splitFilesDir, files := internal.Directory(splitFilesDir)
	spreadOutFile := ""
	if bqsr != "" {
		log.Println("Filtering (phase 1)...")
		// phase 1: Recalibration
		tabsName := filepath.Join(tmpPath, fmt.Sprintf("elprep-tabs-%s", uuidStamp)) + string(filepath.Separator)
		tabsDir := internal.FilepathAbs(tabsName) + string(filepath.Separator)
		internal.MkdirAll(filepath.Dir(tabsDir), 0700)
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
				fileFilterArgs = append(fileFilterArgs, "--mark-optical-duplicates-intermediate", metricsFile)
				opticalDuplicatesPixelDistanceString := strconv.FormatInt(int64(opticalDuplicatesPixelDistance), 10)
				fileFilterArgs = append(fileFilterArgs, "--optical-duplicates-pixel-distance", opticalDuplicatesPixelDistanceString)

			}
			if targetRegions != "" {
				fileFilterArgs = append(fileFilterArgs, "--target-regions", targetRegions)
			}
			filterCmd := exec.Command(os.Args[0], fileFilterArgs...)
			filterCmd.Stderr = os.Stderr
			internal.RunCmd(filterCmd)
		}
		// spread file
		if !singleEnd {
			spreadFileName := fmt.Sprintf("%s-spread.%s", intermediateOutputPrefix, intermediateOutputType)
			spreadInFile := path.Join(splitsDir, spreadFileName)
			spreadOutFile = path.Join(mergeDir, spreadFileName)
			spreadTableFileName := fmt.Sprintf("%s-spread.elrecal", intermediateOutputPrefix)
			spreadTableFile := path.Join(tabsDir, spreadTableFileName)
			fileFilterArgs := append(filterArgs, "--bqsr-tables-only", spreadTableFile)
			filterOpt := []string{"filter", spreadInFile, spreadOutFile}
			fileFilterArgs = append(filterOpt, fileFilterArgs...)
			if markOpticalDuplicates != "" {
				metricsFile := path.Join(metricsDir, spreadFileName)
				fileFilterArgs = append(fileFilterArgs, "--mark-optical-duplicates-intermediate", metricsFile)
				opticalDuplicatesPixelDistanceString := strconv.FormatInt(int64(opticalDuplicatesPixelDistance), 10)
				fileFilterArgs = append(fileFilterArgs, "--optical-duplicates-pixel-distance", opticalDuplicatesPixelDistanceString)
			}
			if targetRegions != "" {
				fileFilterArgs = append(fileFilterArgs, "--target-regions", targetRegions)
			}
			filterCmd := exec.Command(os.Args[0], fileFilterArgs...)
			filterCmd.Stderr = os.Stderr
			internal.RunCmd(filterCmd)
		}
	} else {
		// no bqsr
		// first spread file
		if vcfDir != "" {
			log.Println("Filtering (phase 1) and variant calling...")
		} else {
			log.Println("Filtering (phase 1)...")
		}
		if !singleEnd {
			spreadFileName := fmt.Sprintf("%s-spread.%s", intermediateOutputPrefix, intermediateOutputType)
			spreadInFile := path.Join(splitsDir, spreadFileName)
			spreadOutFile = path.Join(mergeDir, spreadFileName)
			filterOpt := []string{"filter", spreadInFile, spreadOutFile}
			fileFilterArgs := append(filterOpt, filterArgs...)
			if markOpticalDuplicates != "" {
				metricsFile := path.Join(metricsDir, spreadFileName)
				fileFilterArgs = append(fileFilterArgs, "--mark-optical-duplicates-intermediate", metricsFile)
				opticalDuplicatesPixelDistanceString := strconv.FormatInt(int64(opticalDuplicatesPixelDistance), 10)
				fileFilterArgs = append(fileFilterArgs, "--optical-duplicates-pixel-distance", opticalDuplicatesPixelDistanceString)
			}
			filterCmd := exec.Command(os.Args[0], fileFilterArgs...)
			filterCmd.Stderr = os.Stderr
			internal.RunCmd(filterCmd)
		}
		// rest of the files
		sort.Strings(files)
		for _, fileName := range files {
			inFile := path.Join(splitFilesDir, fileName)
			outFile := path.Join(mergeDir, fileName)
			filterOpt := []string{"filter", inFile, outFile}
			fileFilterArgs := append(filterOpt, filterArgs...)
			if markOpticalDuplicates != "" {
				metricsFile := path.Join(metricsDir, fileName)
				fileFilterArgs = append(fileFilterArgs, "--mark-optical-duplicates-intermediate", metricsFile)
				opticalDuplicatesPixelDistanceString := strconv.FormatInt(int64(opticalDuplicatesPixelDistance), 10)
				fileFilterArgs = append(fileFilterArgs, "--optical-duplicates-pixel-distance", opticalDuplicatesPixelDistanceString)
			}
			unmappedFileName := fmt.Sprintf("%s-unmapped.%s", intermediateOutputPrefix, intermediateOutputType)
			if vcfDir != "" && fileName != unmappedFileName { // skip unmapped reads for variant calling
				vcfFile := path.Join(vcfDir, fmt.Sprintf("%s.vcf.gz", fileName))
				// point output to intermediate vcfFile
				fileFilterArgs = append(fileFilterArgs, "--haplotypecaller", vcfFile)
				// add spread file
				if spreadOutFile != "" {
					fileFilterArgs = append(fileFilterArgs, "--spread-file", spreadOutFile)
				}
				if targetRegions != "" {
					fileFilterArgs = append(fileFilterArgs, "--target-regions", targetRegions)
				}
				if randomSeedFile != "" {
					fileFilterArgs = append(fileFilterArgs, "--random-seed-file", randomSeedFile)
				}
			}
			filterCmd := exec.Command(os.Args[0], fileFilterArgs...)
			filterCmd.Stderr = os.Stderr
			internal.RunCmd(filterCmd)
		}
	}
	if err := os.RemoveAll(splitsDir); err != nil {
		log.Panic(err)
	}
	// merge and filtering phase 2 and variant calling
	if bqsr != "" && vcfOutput == "" {
		log.Println("Filtering (phase 2) and merging sam/bam...")
		mergeOpt := []string{"merge", mergeDir, "/dev/stdout"}
		mergeArgs = append(mergeOpt, mergeArgs...)
		mergeCmd := exec.Command(os.Args[0], mergeArgs...)
		outPipe, err := mergeCmd.StdoutPipe()
		if err != nil {
			log.Panic(err)
		}
		mergeCmd.Stderr = os.Stderr
		// phase 2: apply bqsr
		filterOpt2 := []string{"filter", "/dev/stdin", output}
		filterArgs2 = append(filterOpt2, filterArgs2...)
		tabsName := filepath.Join(tmpPath, fmt.Sprintf("elprep-tabs-%s", uuidStamp)) + string(filepath.Separator)
		tabsDir := internal.FilepathAbs(tabsName) + string(filepath.Separator)
		filterArgs2 = append(filterArgs2, "--bqsr-apply", tabsDir, "--recal-file", bqsr, "--bqsr-tablename-prefix", internal.BQSRTablenamePrefix)
		applyBqsrCommand := exec.Command(os.Args[0], filterArgs2...)
		applyBqsrCommand.Stdin = outPipe
		if output == "/dev/stdout" {
			applyBqsrCommand.Stdout = os.Stdout
		}
		applyBqsrCommand.Stderr = os.Stderr
		if err := mergeCmd.Start(); err != nil {
			log.Panic(err)
		}
		if err := applyBqsrCommand.Start(); err != nil {
			log.Panic(err)
		}
		if err := mergeCmd.Wait(); err != nil {
			log.Panic(err)
		}
		if err := applyBqsrCommand.Wait(); err != nil {
			log.Panic(err)
		}
		if err := os.RemoveAll(tabsDir); err != nil {
			log.Panic(err)
		}
	} else if bqsr != "" && vcfOutput != "" {
		// apply bqsr and haplotype caller
		// set up directory for metrics
		log.Println("Filtering (phase 2) and variant calling...")
		mergeName := filepath.Join(tmpPath, fmt.Sprintf("elprep-merge2-%s", uuidStamp))
		newMergeDir := internal.FilepathAbs(mergeName) + string(filepath.Separator)
		internal.MkdirAll(filepath.Dir(newMergeDir), 0700)
		_, files := internal.Directory(mergeDir)
		spreadFileName := fmt.Sprintf("%s-spread.%s", intermediateOutputPrefix, intermediateOutputType)
		spreadInFile := path.Join(mergeDir, spreadFileName)
		unmappedFileName := fmt.Sprintf("%s-unmapped.%s", intermediateOutputPrefix, intermediateOutputType)
		tabsName := filepath.Join(tmpPath, fmt.Sprintf("elprep-tabs-%s", uuidStamp)) + string(filepath.Separator)
		tabsDir := internal.FilepathAbs(tabsName) + string(filepath.Separator)
		// apply bqsr on the spread file so we are sure the processed spread file is used by haplotypecaller
		spreadOutFile := ""
		if spreadInFile != "" {
			spreadOutFile = path.Join(mergeDir, fmt.Sprintf("apply-bqsred-%s", spreadFileName))
			filterOpt2 := []string{"filter", spreadInFile, spreadOutFile}
			filterArgs3 := append(filterArgs2, "--bqsr-apply", tabsDir, "--recal-file", bqsr, "--bqsr-tablename-prefix", internal.BQSRTablenamePrefix)
			filterArgs3 = append(filterOpt2, filterArgs3...)
			filterCmd := exec.Command(os.Args[0], filterArgs3...)
			filterCmd.Stderr = os.Stderr
			internal.RunCmd(filterCmd)
		}
		sort.Strings(files)
		for _, fileName := range files {
			// skip spread file
			if fileName != spreadFileName {
				// apply bqsr and haplotype filter command for each split file
				inFile := path.Join(mergeDir, fileName)
				outFile := path.Join(newMergeDir, fileName)
				vcfOutFile := path.Join(vcfDir, fmt.Sprintf("%s.vcf.gz", fileName))
				filterOpt2 := []string{"filter", inFile, outFile}
				filterArgs3 := append(filterArgs2, "--bqsr-apply", tabsDir, "--recal-file", bqsr, "--bqsr-tablename-prefix", internal.BQSRTablenamePrefix)
				// skip unmapped reads for variant calling
				if fileName != unmappedFileName {
					filterArgs3 = append(filterArgs3, "--haplotypecaller", vcfOutFile, "--pg-cmd-line", commandString)
					if randomSeedFile != "" {
						filterArgs3 = append(filterArgs3, "--random-seed-file", randomSeedFile)
					}
				}
				if targetRegions != "" && fileName != unmappedFileName {
					filterArgs3 = append(filterArgs3, "--target-regions", targetRegions)
				}
				filterArgs3 = append(filterOpt2, filterArgs3...)
				if spreadOutFile != "" && fileName != unmappedFileName { // skipped unmapped reads for variant calling
					filterArgs3 = append(filterArgs3, "--spread-file", spreadOutFile)
				}
				filterCmd := exec.Command(os.Args[0], filterArgs3...)
				filterCmd.Stderr = os.Stderr
				internal.RunCmd(filterCmd)
			}
		}
		// remove intermediate recal tables
		if err := os.RemoveAll(tabsDir); err != nil {
			log.Panic(err)
		}
		log.Println("Merging sam/bam...")
		// merge
		mergeOpt := []string{"merge", newMergeDir, output}
		mergeArgs = append(mergeArgs, "--ignore-spread-file")
		mergeArgs = append(mergeOpt, mergeArgs...)
		mergeCmd := exec.Command(os.Args[0], mergeArgs...)
		if output == "/dev/stdout" {
			mergeCmd.Stdout = os.Stdout
		}
		mergeCmd.Stderr = os.Stderr
		internal.RunCmd(mergeCmd)
		if err := os.RemoveAll(newMergeDir); err != nil {
			log.Panic(err)
		}
	} else {
		log.Println("Merging sam/bam...")
		mergeOpt := []string{"merge", mergeDir, output}
		mergeArgs = append(mergeOpt, mergeArgs...)
		mergeCmd := exec.Command(os.Args[0], mergeArgs...)
		if output == "/dev/stdout" {
			mergeCmd.Stdout = os.Stdout
		}
		mergeCmd.Stderr = os.Stderr
		internal.RunCmd(mergeCmd)
	}
	if err := os.RemoveAll(mergeDir); err != nil {
		log.Panic(err)
	}
	// merge intermediate metrics files
	if markOpticalDuplicates != "" {
		log.Println("Merging metrics...")
		ctr := filters.LoadAndCombineDuplicateMetrics(metricsDir)
		filters.PrintDuplicatesMetrics(markOpticalDuplicates, commandString, ctr)
		if err := os.RemoveAll(metricsDir); err != nil {
			log.Panic(err)
		}
	}
	// merge intermediate vcf files
	if vcfOutput != "" {
		log.Println("Merge vcf...")
		filters.CombineVcfOutputs(vcfDir, vcfOutput)
		if err := os.RemoveAll(vcfDir); err != nil {
			log.Panic(err)
		}
	}
}
