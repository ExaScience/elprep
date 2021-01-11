// elPrep: a high-performance tool for analyzing SAM/BAM files.
// Copyright (c) 2020 imec vzw.

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

package filters

import (
	"fmt"
	"io"
	"log"
	"math"
	"runtime"
	"sort"
	"sync"
	"time"

	"github.com/exascience/elprep/v5/internal"
	"github.com/exascience/pargo/pipeline"

	"github.com/exascience/elprep/v5/vcf"

	"github.com/exascience/elprep/v5/intervals"

	"github.com/exascience/elprep/v5/bed"

	"github.com/exascience/elprep/v5/utils"

	"github.com/exascience/pargo/parallel"

	"github.com/exascience/elprep/v5/fasta"

	"github.com/exascience/elprep/v5/sam"
)

var (
	// DP represenst depth in VCF files.
	DP = utils.Intern("DP")

	// ExcessHet represenst Excess Heterozygosity in VCF files.
	ExcessHet = utils.Intern("ExcessHet")
	// MLEAC represenst the maximum likelihood expectation for the allele count in VCF files.
	MLEAC = utils.Intern("MLEAC")
	// MLEAF represents the maximum likelihood expectation for the allele frequency in VCF files.
	MLEAF = utils.Intern("MLEAF")
	// LowQual represents low quality variants in VCF files.
	LowQual = utils.Intern("LowQual")
	// Raw_MQandDP represents raw squared mapping quality and depth in VCF files.
	RAW_MQandDP = utils.Intern("RAW_MQandDP")
	// MQ represents root mean square of mapping quality in VCF files
	MQ = utils.Intern("MQ")
	// BaseQRankSum represents rank sum test of alt vs ref base qualities in VCF files.
	BaseQRankSum = utils.Intern("BaseQRankSum")
	// MQRankSum represents rank sum test of alt vs ref read mapping qualities in VCF files.
	MQRankSum = utils.Intern("MQRankSum")
	// ReadPosRankSum represents rank sum test of alt vs ref read position bias in VCF files.
	ReadPosRankSum = utils.Intern("ReadPosRankSum")

	// AC represents allele count in genotypes for each ALT allele in VCF files.
	AC = utils.Intern("AC")
	// AF represents allele frequency for each ALT allele in VCF files.
	AF = utils.Intern("AF")
	// AN represents total number of alleles in called genotypes in VCF files.
	AN = utils.Intern("AN")

	// AD represents allelic depths in VCF files.
	AD = utils.Intern("AD")
	// GQ represents genotype quality in VCF files.
	GQ = utils.Intern("GQ")
	// MIN_DP represents minimum depths in GVCF blocks in VCF files.
	MIN_DP = utils.Intern("MIN_DP")
	// PL represents likelihoods for genotypes in VCF files.
	PL = utils.Intern("PL")
	// SB represents strand bias in VCF files.
	SB = utils.Intern("SB")

	// FS represents phred-scaled p-value using Fisher's exact test to detect strand bias in VCF files.
	FS = utils.Intern("FS")
	// QD represents variant confidence/quality by depth in VCF files.
	QD = utils.Intern("QD")
	// SOR represents symmetric odds ratio of 2x2 contingency table to detect strand bias in VCF files.
	SOR = utils.Intern("SOR")

	noVariationFormatForGvcf       = []utils.Symbol{vcf.GT, DP, GQ, MIN_DP, PL}
	noVariationFormatForNonGvcf    = []utils.Symbol{vcf.GT, AD, DP, GQ, PL}
	noCallGT                       = []int32{-1, -1}
	noVariationGT                  = []int32{0, 0}
	cigarConsumesReferenceBasesOrS = map[byte]int32{'M': 1, 'D': 1, 'N': 1, '=': 1, 'X': 1, 'S': 1}
)

const (
	nonRef = "<NON_REF>"

	log10InformativeThreshold = 0.2
)

/*
	note: high quality soft clips are treated weirdly
	there was a discussion to remove them completely,
	but it was decided to ignore this for now
	we may have to check that later again
*/

type confidenceMode uint8

const (
	none confidenceMode = iota
	bpResolution
	gvcf
)

type (
	// an assemblyRegion is an either active or inactive short stretch of reads
	assemblyRegion struct {
		// the contig that this region falls in
		contig string
		// the reference
		reference []byte
		// the reads that cover this region
		alns []*sam.Alignment
		// the activity states for each reference location that this region covers
		supportingStates []float64
		// start and end reference locations
		start, end, extension, contigLength int32
		// active or not?
		isActive bool
		// a side channel for handling deletions in calculateGenotypes
		deletions sideChannel
	}

	haplotype struct {
		isRef    bool
		bases    string
		location int32
		events   []*vcf.Variant
		cigar    []sam.CigarOperation
		score    float64
	}

	// parameters that influence how assembly regions are determined
	HaplotypeCaller struct {
		reference                                        *fasta.MappedFasta
		activityProfile, assemblyRegions                 io.Writer
		activeProbThreshold                              float64
		maxProbPropagationDistance                       int32
		minRegionSize, maxRegionSize, padding            int32
		minBaseQual                                      byte
		confidenceMode                                   confidenceMode
		refPseudocount, snpPseudocount, indelPseudocount float64
		log10Priors                                      [3]float64
		log10ACgt0Prior                                  float64
		standardConfidenceForCalling                     float64
		standardConfidenceForCallingByMin10              float64 // standard confidence for calling is divided by -10 almost every time it is accessed
		standardConfidenceForActivityByMin10             float64
		random                                           *internal.Rand
		maxReadsPerAlignmentStart                        int32
		indelSizeToEliminateInRefModel                   int32
		useSoftClippedBases                              bool
		gqBands                                          []int
		commandLine                                      string
		randomSeedFile                                   string
	}
)

func (region *assemblyRegion) paddedStart() int32 {
	return maxInt32(1, region.start-region.extension)
}

func (region *assemblyRegion) paddedEnd() int32 {
	return minInt32(region.contigLength, region.end+region.extension)
}

func makeRandom(randomSeedFile string) *internal.Rand {
	if randomSeedFile != "" {
		rsf := internal.FileOpen(randomSeedFile)
		var content string
		fmt.Fscanln(rsf, &content)
		internal.Close(rsf)
		if content != "init" {
			seed := internal.ParseInt(content, 10, 64)
			return internal.NewRandReflect(seed)
		}
	}
	return internal.NewRand(47382911)
}

// NewHaplotypeCaller creates an object that contains the relevent
// parameters for the haplotypecaller.
func NewHaplotypeCaller(
	reference *fasta.MappedFasta,
	referenceConfidence string,
	assemblyRegionPadding int32,
	activityProfile, assemblyRegions io.Writer,
	randomSeedFile string,
	commandLine string) *HaplotypeCaller {

	var confidenceMode confidenceMode
	var standardConfidenceForCalling float64

	switch referenceConfidence {
	case "NONE":
		confidenceMode = none
		standardConfidenceForCalling = 30.0 // todo: command line parameter; the default depends on confidenceMode
	case "GVCF":
		confidenceMode = gvcf
		standardConfidenceForCalling = -0.0 // todo: command line parameter; the default depends on confidenceMode
	case "BP_RESOLUTION":
		confidenceMode = bpResolution
		standardConfidenceForCalling = -0.0 // todo: command line parameter; the default depends on confidenceMode
	default:
		log.Panicf("invalid reference confidence mode %v", referenceConfidence)
	}

	standardConfidenceForActivity := math.Min(4, standardConfidenceForCalling)

	var log10Priors [3]float64
	heterozygosity := 0.001        // todo: command line parameter
	indelHeterozygosity := 1.25e-4 // todo: command line parameter
	heterozygosityStddev := 0.01   // todo: command line parameter

	log10Heterozygosity := log10(heterozygosity)
	log10Priors[1] = log10Heterozygosity - log10(1)
	log10Priors[2] = log10Heterozygosity - log10(2)
	log10Sum := approximateLog10SumLog10(log10Priors[1], log10Priors[2])
	if log10Sum >= 0 {
		log.Panicf("heterozygosity %v is too large for total ploidy 2", heterozygosity)
	}
	log10Priors[0] = log10OneMinusPow10(log10Sum)

	refPseudocount := heterozygosity / math.Pow(heterozygosityStddev, 2)
	snpPseudocount := heterozygosity * refPseudocount
	indelPseudocount := indelHeterozygosity * refPseudocount

	return &HaplotypeCaller{
		reference:                            reference,
		activityProfile:                      activityProfile,
		assemblyRegions:                      assemblyRegions,
		activeProbThreshold:                  0.002,                              // todo: command line parameter
		maxProbPropagationDistance:           50 + int32(len(gaussianKernel)>>1), // todo: command line parameter
		minRegionSize:                        50,                                 // todo: command line parameter
		maxRegionSize:                        300,                                // todo: command line parameter
		padding:                              assemblyRegionPadding,
		minBaseQual:                          10,             // todo: commond line parameter
		confidenceMode:                       confidenceMode, // todo: command line parameter
		refPseudocount:                       refPseudocount,
		snpPseudocount:                       snpPseudocount,
		indelPseudocount:                     indelPseudocount,
		log10Priors:                          log10Priors,
		log10ACgt0Prior:                      approximateLog10SumLog10(log10Priors[1], log10Priors[2]),
		standardConfidenceForCalling:         standardConfidenceForCalling, // todo: command line parameter
		standardConfidenceForCallingByMin10:  standardConfidenceForCalling / -10,
		standardConfidenceForActivityByMin10: standardConfidenceForActivity / -10,
		random:                               makeRandom(randomSeedFile),
		maxReadsPerAlignmentStart:            50,   // todo: command line parameter
		useSoftClippedBases:                  true, // todo: command line parameter
		indelSizeToEliminateInRefModel:       10,   // todo: command line parameter
		gqBands: []int{ // todo: command line parameter (must be sorted, only unique entries, and if 0 and 100 not provided, start with 0 and end in 100)
			0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
			10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
			20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
			30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
			40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
			50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
			60, 70, 80, 90, 99, 100,
		},
		commandLine:    commandLine,
		randomSeedFile: randomSeedFile,
	}
}

// apply a Gaussian curve to the activity state at a reference location,
// and then store the values in the states vector
func (hc *HaplotypeCaller) bandPassProcessState(states []float64, pos int32, state float64) {
	filterSize := int32(len(gaussianKernel) >> 1)
	bandStart := -filterSize
	if pos+bandStart < 0 {
		bandStart = -pos
	}
	bandEnd := filterSize
	if pos+bandEnd >= int32(len(states)) {
		bandEnd = int32(len(states)) - 1 - pos
	}
	for i := bandStart; i <= bandEnd; i++ {
		states[pos+i] += state * gaussianKernel[i+filterSize]
	}
}

// before calling bandPassProcessStates, check for high-quality soft clips mean
// if this is above a certain threshold, then instead of calling bandPassProcessStates just once
// call it 2*int(hqSoftClipsMean)+1 times [!!!]
// (this behavior doesn't make sense, but it's not yet clear how to address it as of mid-April 2019)
func (hc *HaplotypeCaller) processState(states []float64, pos int32, state float64, hqSoftClipsMean float64) {
	if state > 0 {
		if hqSoftClipsMean > averageHighQualitySoftClipsBasesThreshold {
			numHQClips := int32(hqSoftClipsMean)
			if hc.maxProbPropagationDistance < numHQClips {
				numHQClips = hc.maxProbPropagationDistance
			}
			/*
				for i := -numHQClips; i <= numHQClips; i++ {
					hc.bandPassProcessState(states, pos, state) // I believe this should be pos+i
				}
			*/
			// optimized version that implements the effective semantics
			hc.bandPassProcessState(states, pos, state*float64(2*int(numHQClips)+1))
		} else {
			hc.bandPassProcessState(states, pos, state)
		}
	} else if state < 0 {
		log.Panic("Need to cover the case in processState in HaplotypeCaller where state is less than zero.")
	}
}

// compute the assembly regions for the given reads and the given activiy states vector.
func (hc *HaplotypeCaller) computeAssemblyRegions(region computeRegion, states []float64) (regions []*assemblyRegion) {
	start, stop := region.start, region.stop
	for size := stop - start; size > 0; size = stop - start {
		isActiveRegion := states[start-region.start] > hc.activeProbThreshold
		maxRegionEnd := start + hc.maxRegionSize
		var maxActivityBoundary int32
		if size < hc.maxRegionSize {
			maxActivityBoundary = stop
		} else {
			maxActivityBoundary = maxRegionEnd
		}
		endOfActiveRegion := start + 1
		for ; endOfActiveRegion < maxActivityBoundary; endOfActiveRegion++ {
			if states[endOfActiveRegion-region.start] > hc.activeProbThreshold != isActiveRegion {
				break
			}
		}
		if isActiveRegion && endOfActiveRegion == maxRegionEnd {
			minI := endOfActiveRegion - 1
			minP := math.MaxFloat64
			top := minI
			if top == stop-1 {
				top--
			}
			bottom := start + hc.minRegionSize - 1 // assumption: minRegionSize > 1
			for i := top; i >= bottom; i-- {
				if cur := states[i-region.start]; cur < minP && cur <= states[i+1-region.start] && cur < states[i-1-region.start] {
					minI = i
					minP = cur
				}
			}
			endOfActiveRegion = minI + 1
		}
		var supportingStates []float64
		if hc.activityProfile != nil {
			supportingStates = states[start-region.start : endOfActiveRegion-region.start]
		}
		regions = append(regions, &assemblyRegion{
			contig:           region.contigName,
			reference:        region.reference,
			supportingStates: supportingStates,
			start:            start + 1,
			end:              endOfActiveRegion,
			extension:        hc.padding,
			contigLength:     region.contigLength,
			isActive:         isActiveRegion,
		})
		start = endOfActiveRegion
	}
	return
}

// fill in the reads that overlap this assembly region.
// maxReferenceLength is the maximum reference length for the reads in this contig.
// this is used to ensure that we don't have to look at all the reads from this contig,
// which speeds up this part considerably.
func (hc *HaplotypeCaller) fillRegionAlns(region *assemblyRegion, alns []*sam.Alignment, maxReferenceLength int32) {
	region.alns, _ = alnSlice(alns, region.paddedStart(), region.paddedEnd(), maxReferenceLength)
}

// activityInfo records the activity states and high-quality soft clips means
// for a (subslice of a) reference region
type activityInfo struct {
	start                     int32
	isActive, hqSoftClipsMean []float64
}

type computeRegion struct {
	contigName                                    string
	reference                                     []byte
	alns                                          []*sam.Alignment
	start, stop, maxReferenceLength, contigLength int32
}

func (region computeRegion) paddedStart(hc *HaplotypeCaller) int32 {
	return maxInt32(0, region.start-hc.padding)
}

func (region computeRegion) paddedStop(hc *HaplotypeCaller) int32 {
	return minInt32(region.contigLength, region.stop+hc.padding)
}

func (hc *HaplotypeCaller) finalizeAssemblyRegion(region *assemblyRegion) {
	paddedStart := region.paddedStart()
	paddedEnd := region.paddedEnd()
	aln := new(sam.Alignment)
	for i := 0; i < len(region.alns); {
		*aln = *region.alns[i]
		hardClipLowQualEnds(aln, hc.minBaseQual-1)
		if hc.useSoftClippedBases {
			if wellDefined, _ := hasWellDefinedFragmentSize(aln); wellDefined {
				revertSoftClippedBases(aln)
			} else {
				hardClipSoftClippedBases(aln)
			}
		} else {
			hardClipSoftClippedBases(aln)
		}
		if !isStrictUnmapped(aln) {
			hardClipAdaptorSequence(aln)
		}
		if aln.SEQ.Len() > 0 && sam.ReadLengthFromCigar(aln.CIGAR) > 0 {
			hardClipToRegion(aln, int(paddedStart), int(paddedEnd))
			if readOverlapsRegion(aln, paddedStart, paddedEnd) {
				region.alns[i] = aln
				i++
				aln = new(sam.Alignment)
				continue
			}
		}
		region.alns = append(region.alns[:i], region.alns[i+1:]...)
	}

	sam.By(sam.CoordinateLess).ParallelStableSort(region.alns)
	forEachReadPair(region.alns, cleanOverlappingReadPair)
}

func (hc *HaplotypeCaller) writeVcfHeader(input *sam.Sam, sampleName string, vcfFile *vcf.OutputFile) {
	vcfHeader := vcf.NewHeader()
	vcfHeader.Meta["elPrepCommandLine"] = []interface{}{
		&vcf.MetaInformation{
			ID: utils.Intern(utils.ProgramName),
			Fields: utils.StringMap{
				"CommandLine": hc.commandLine,
				"Version":     utils.ProgramVersion,
				"URL":         utils.ProgramURL,
				"Date":        time.Now().Format("Mon Jan 02 15:04:05 MST 2006"),
			},
		},
	}
	vcfHeader.Infos = []*vcf.FormatInformation{
		{
			ID:          BaseQRankSum,
			Description: "Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities",
			Number:      1,
			Type:        vcf.Float,
		},
		{
			ID:          DP,
			Description: "Approximate read depth; some reads may have been filtered",
			Number:      1,
			Type:        vcf.Integer,
		},
		{
			ID:          utils.Intern("DS"),
			Description: "Were any of the samples downsampled?",
			Number:      0,
			Type:        vcf.Flag,
		},
		{
			ID:          ExcessHet,
			Description: "Phred-scaled p-value for exact test of excess heterozygosity",
			Number:      1,
			Type:        vcf.Float,
		},
		{
			ID:          utils.Intern("InbreedingCoeff"),
			Description: "Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation",
			Number:      1,
			Type:        vcf.Float,
		},
		{
			ID:          MLEAC,
			Description: "Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed",
			Number:      vcf.NumberA,
			Type:        vcf.Integer,
		},
		{
			ID:          MLEAF,
			Description: "Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed",
			Number:      vcf.NumberA,
			Type:        vcf.Float,
		},
		{
			ID:          MQRankSum,
			Description: "Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities",
			Number:      1,
			Type:        vcf.Float,
		},
		{
			ID:          ReadPosRankSum,
			Description: "Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias",
			Number:      1,
			Type:        vcf.Float,
		},
	}
	vcfHeader.Formats = []*vcf.FormatInformation{
		{
			ID:          AD,
			Description: "Allelic depths for the ref and alt alleles in the order listed",
			Number:      vcf.NumberR,
			Type:        vcf.Integer,
		},
		{
			ID:          DP,
			Description: "Approximate read depth (reads with MQ=255 or with bad mates are filtered)",
			Number:      1,
			Type:        vcf.Integer,
		},
		{
			ID:          GQ,
			Description: "Genotype Quality",
			Number:      1,
			Type:        vcf.Integer,
		},
		{
			ID:          vcf.GT,
			Description: "Genotype",
			Number:      1,
			Type:        vcf.String,
		},
		{
			ID:          PL,
			Description: "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification",
			Number:      vcf.NumberG,
			Type:        vcf.Integer,
		},
	}
	vcfHeader.Meta["FILTER"] = []interface{}{
		&vcf.MetaInformation{
			ID:          LowQual,
			Description: "Low quality",
		},
	}
	for _, sq := range input.Header.SQ {
		contig := sq["SN"]
		vcfHeader.Meta["contig"] = append(vcfHeader.Meta["contig"], interface{}(
			&vcf.MetaInformation{
				ID:     utils.Intern(contig),
				Fields: utils.StringMap{"length": fmt.Sprint(len(hc.reference.Seq(contig)))},
			}))
	}
	vcfHeader.Meta["source"] = []interface{}{"HaplotypeCaller"}
	if hc.confidenceMode == none {
		vcfHeader.Infos = append(vcfHeader.Infos,
			&vcf.FormatInformation{
				ID:          AC,
				Description: "Allele count in genotypes, for each ALT allele, in the same order as listed",
				Number:      vcf.NumberA,
				Type:        vcf.Integer,
			},
			&vcf.FormatInformation{
				ID:          AF,
				Description: "Allele Frequency, for each ALT allele, in the same order as listed",
				Number:      vcf.NumberA,
				Type:        vcf.Float,
			},
			&vcf.FormatInformation{
				ID:          AN,
				Description: "Total number of alleles in called genotypes",
				Number:      1,
				Type:        vcf.Integer,
			},
			&vcf.FormatInformation{
				ID:          FS,
				Description: "Phred-scaled p-value using Fisher's exact test to detect strand bias",
				Number:      1,
				Type:        vcf.Float,
			},
			&vcf.FormatInformation{
				ID:          MQ,
				Description: "RMS Mapping Quality",
				Number:      1,
				Type:        vcf.Float,
			},
			&vcf.FormatInformation{
				ID:          QD,
				Description: "Variant Confidence/Quality by Depth",
				Number:      1,
				Type:        vcf.Float,
			},
			&vcf.FormatInformation{
				ID:          SOR,
				Description: "Symmetric Odds Ratio of 2x2 contingency table to detect strand bias",
				Number:      1,
				Type:        vcf.Float,
			},
		)
	} else {
		vcfHeader.Infos = append(vcfHeader.Infos, &vcf.FormatInformation{
			ID:          RAW_MQandDP,
			Description: "Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.",
			Number:      2,
			Type:        vcf.Integer,
		})
		vcfHeader.Formats = append(vcfHeader.Formats,
			&vcf.FormatInformation{
				ID:          PGT,
				Description: "Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another",
				Number:      1,
				Type:        vcf.String,
			},
			&vcf.FormatInformation{
				ID:          PID,
				Description: "Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group",
				Number:      1,
				Type:        vcf.String,
			},
			&vcf.FormatInformation{
				ID:          PS,
				Description: "Phasing set (typically the position of the first variant in the set)",
				Number:      1,
				Type:        vcf.Integer,
			},
			&vcf.FormatInformation{
				ID:          SB,
				Description: "Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.",
				Number:      4,
				Type:        vcf.Integer,
			},
		)
		vcfHeader.Meta["ALT"] = []interface{}{
			&vcf.MetaInformation{
				ID:          utils.Intern("NON_REF"),
				Description: "Represents any possible alternative allele at this location",
			},
		}
		if hc.confidenceMode == gvcf {
			vcfHeader.Infos = append(vcfHeader.Infos, &vcf.FormatInformation{
				ID:          vcf.END,
				Description: "Stop position of the interval",
				Number:      1,
				Type:        vcf.Integer,
			})
			vcfHeader.Formats = append(vcfHeader.Formats,
				&vcf.FormatInformation{
					ID:          MIN_DP,
					Description: "Minimum DP observed within the GVCF block",
					Number:      1,
					Type:        vcf.Integer,
				})
			for i := 1; i < len(hc.gqBands); i++ {
				min := hc.gqBands[i-1]
				max := hc.gqBands[i]
				vcfHeader.Meta[fmt.Sprintf("GVCFBlock%d-%d", min, max)] = []interface{}{fmt.Sprintf("minGQ=%d(inclusive),maxGQ=%d(exclusive)", min, max)}
			}
		}
	}
	sort.Slice(vcfHeader.Formats, func(i, j int) bool {
		return *vcfHeader.Formats[i].ID < *vcfHeader.Formats[j].ID
	})
	sort.Slice(vcfHeader.Infos, func(i, j int) bool {
		return *vcfHeader.Infos[i].ID < *vcfHeader.Infos[j].ID
	})
	vcfHeader.Columns = append(vcfHeader.Columns, "FORMAT", sampleName)
	vcfHeader.Format(vcfFile.Writer)
}

// CallVariants is the main entry point into the haplotypecaller.
func (hc *HaplotypeCaller) CallVariants(input *sam.Sam, sampleName, vcfOutput string, bedRegions bed.Bed) {
	var validContigs map[string]bool
	if headerContigs, ok := input.Header.Contigs(); ok {
		validContigs = make(map[string]bool)
		for _, contig := range headerContigs {
			validContigs[contig] = true
		}
	}

	var processedBedRegions utils.SmallMap
	if bedRegions != nil {
		processedBedRegions = intervals.FromBedOrdered(bedRegions)
		parallel.Range(0, len(processedBedRegions), 0, func(low, high int) {
			for i := low; i < high; i++ {
				processedBedRegions[i].Value = intervals.ParallelFlatten(processedBedRegions[i].Value.([]intervals.Interval))
			}
		})
	}

	vcfFile := vcf.Create(vcfOutput, "", -1)
	defer vcfFile.Close()

	if hc.assemblyRegions != nil {
		fmt.Fprintln(hc.assemblyRegions, "#track graphType=line")
		fmt.Fprintln(hc.assemblyRegions, "Chromosome\tStart\tEnd\tFeature\tAssemblyRegions")
	}
	if hc.activityProfile != nil {
		fmt.Fprintln(hc.activityProfile, "#track graphType=line")
		fmt.Fprintln(hc.activityProfile, "Chromosome\tStart\tEnd\tFeature\tActivityProfile")
	}

	hc.writeVcfHeader(input, sampleName, vcfFile)

	alnsByContig := make(map[string][]*sam.Alignment)
	{
		contigBoundaries := make([]int, len(input.Header.SQ)+1)
		parallel.Range(1, len(contigBoundaries), len(contigBoundaries), func(low, high int) {
			for i := int32(low); i < int32(high); i++ {
				contigBoundaries[i] = sort.Search(len(input.Alignments), func(index int) bool {
					refid := input.Alignments[index].REFID()
					return refid == -1 || refid >= i
				})
			}
		})
		if validContigs != nil {
			for index, record := range input.Header.SQ {
				contigName := record["SN"]
				if validContigs[contigName] {
					alnsByContig[contigName] = input.Alignments[contigBoundaries[index]:contigBoundaries[index+1]]
				}
			}
		} else {
			for index, record := range input.Header.SQ {
				alnsByContig[record["SN"]] = input.Alignments[contigBoundaries[index]:contigBoundaries[index+1]]
			}
		}
		input.Alignments = nil
	}
	computeRegionChannel := make(chan computeRegion, runtime.GOMAXPROCS(0))
	if len(processedBedRegions) > 0 {
		go func() {
			defer close(computeRegionChannel)
			for _, entry := range processedBedRegions {
				contigName := *entry.Key
				if validContigs != nil && !validContigs[contigName] {
					continue
				}
				contigRegions := entry.Value.([]intervals.Interval)
				contigIndex := utils.Find(input.Header.SQ, func(record utils.StringMap) bool {
					return record["SN"] == contigName
				})
				if contigIndex < 0 {
					log.Panicf("contig %v from BED file not found in SAM file header", contigName)
				}
				contigLength := sam.SQLN(input.Header.SQ[contigIndex])
				contigAlns := hc.downsample(alnsByContig[contigName])
				delete(alnsByContig, contigName)
				region := computeRegion{
					contigName:         contigName,
					reference:          hc.reference.Seq(contigName),
					maxReferenceLength: hc.maxReferenceLength(contigAlns),
					contigLength:       contigLength,
				}
				for _, bedRegion := range contigRegions {
					region.start, region.stop = bedRegion.Start, bedRegion.End
					slicedAlns, firstContigAlnsIndex := alnSlice(contigAlns, region.paddedStart(hc), region.paddedStop(hc), region.maxReferenceLength)
					region.alns = slicedAlns
					computeRegionChannel <- region
					for i := 0; i < firstContigAlnsIndex; i++ {
						contigAlns[i] = nil
					}
					contigAlns = contigAlns[firstContigAlnsIndex:]
				}
			}
		}()
	} else {
		go func() {
			defer close(computeRegionChannel)
			for _, contig := range input.Header.SQ {
				contigName, ok := contig["SN"]
				if !ok {
					log.Panic("contig name not found in sequence dictionary")
				}
				if validContigs != nil && !validContigs[contigName] {
					continue
				}
				contigLength := sam.SQLN(contig)
				contigAlns := hc.downsample(alnsByContig[contigName])
				delete(alnsByContig, contigName)
				maxReferenceLength := hc.maxReferenceLength(contigAlns)
				computeRegionChannel <- computeRegion{
					contigName:         contigName,
					reference:          hc.reference.Seq(contigName),
					alns:               contigAlns,
					start:              0,
					stop:               contigLength,
					maxReferenceLength: maxReferenceLength,
					contigLength:       contigLength,
				}
			}
		}()
	}

	var assemblyRegionPipeline pipeline.Pipeline
	assemblyRegionPipeline.Source(pipeline.NewSingletonChan(computeRegionChannel))
	assemblyRegionPipeline.SetVariableBatchSize(1, 1)
	var previousRegion *assemblyRegion
	assemblyRegionChannel := make(chan *assemblyRegion, runtime.GOMAXPROCS(0))
	assemblyRegionPipeline.Add(
		pipeline.LimitedPar(runtime.GOMAXPROCS(0), pipeline.Receive(func(_ int, data interface{}) interface{} {
			region := data.(computeRegion)
			regionSize := region.stop - region.start
			activityStates := make([]float64, regionSize)
			var q pipeline.Pipeline
			start := region.start
			q.Source(pipeline.NewFunc(-1, func(size int) (interface{}, int, error) {
				if start >= region.stop {
					return nil, 0, nil
				}
				stop := start + int32(size)
				if stop >= region.stop {
					stop = region.stop
					size = int(stop - start)
				}
				result := [2]int32{start, stop}
				start = stop
				return result, size, nil
			}))
			q.SetVariableBatchSize(512, 512)
			q.Add(
				pipeline.LimitedPar(0, pipeline.Receive(func(_ int, data interface{}) interface{} {
					loc := data.([2]int32)
					startLoc, stopLoc := loc[0]+1, loc[1]+1
					activity := activityInfo{
						start:           loc[0],
						isActive:        make([]float64, stopLoc-startLoc),
						hqSoftClipsMean: make([]float64, stopLoc-startLoc),
					}
					startPOS := startLoc - region.maxReferenceLength
					startAln := sort.Search(len(region.alns), func(index int) bool {
						return region.alns[index].POS >= startPOS
					})
					stopAln := sort.Search(len(region.alns), func(index int) bool {
						return region.alns[index].POS > stopLoc
					})
					forEachPileup(region.alns[startAln:stopAln], startLoc, stopLoc, func(p *pileup) {
						isActiveProb, hqSoftClipsMean := hc.isActive(p, fasta.ToUpperAndN(region.reference[p.location-1]))
						activity.isActive[p.location-startLoc] = isActiveProb
						activity.hqSoftClipsMean[p.location-startLoc] = hqSoftClipsMean
					})
					return activity
				})),
				pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
					activity := data.(activityInfo)
					for i := int32(0); i < int32(len(activity.isActive)); i++ {
						hc.processState(activityStates, activity.start-region.start+i, activity.isActive[i], activity.hqSoftClipsMean[i])
					}
					return nil
				})),
			)
			internal.RunPipeline(&q)
			assemblyRegions := hc.computeAssemblyRegions(region, activityStates)
			parallel.Range(0, len(assemblyRegions), 0, func(low, high int) {
				for i := low; i < high; i++ {
					hc.fillRegionAlns(assemblyRegions[i], region.alns, region.maxReferenceLength)
				}
			})
			return assemblyRegions
		})),
		pipeline.StrictOrd(pipeline.ReceiveAndFinalize(func(_ int, data interface{}) interface{} {
			// set up side channel for deletions to be used in calculateGenotypes
			for _, region := range data.([]*assemblyRegion) {
				if previousRegion == nil || previousRegion.contig != region.contig {
					region.deletions.makeInitial()
				} else {
					region.deletions.linkFrom(previousRegion.deletions)
				}
				previousRegion = region
				assemblyRegionChannel <- region
			}
			return data
		}, func() {
			close(assemblyRegionChannel)
		})),
	)
	var waitForProfiles sync.WaitGroup
	if hc.activityProfile != nil || hc.assemblyRegions != nil {
		profilesChannel := make(chan []*assemblyRegion, 2*runtime.GOMAXPROCS(0))
		waitForProfiles.Add(1)
		go func() {
			defer waitForProfiles.Done()
			for data := range profilesChannel {
				printAssemblyRegions(hc.assemblyRegions, hc.activityProfile, data)
			}
		}()
		assemblyRegionPipeline.Add(pipeline.StrictOrd(pipeline.ReceiveAndFinalize(func(_ int, data interface{}) interface{} {
			profilesChannel <- data.([]*assemblyRegion)
			return data
		}, func() {
			close(profilesChannel)
		})))
	}

	var waitForAssemblyRegions sync.WaitGroup
	waitForAssemblyRegions.Add(1)
	go func() {
		defer waitForAssemblyRegions.Done()
		assemblyRegionPipeline.Run()
	}()

	var variantCallPipeline pipeline.Pipeline
	variantCallPipeline.Source(pipeline.NewSingletonChan(assemblyRegionChannel))
	variantCallPipeline.SetVariableBatchSize(1, 1)
	var finalVariant *vcf.Variant
	var finalVariantOK bool

	variantCallPipeline.Add(
		pipeline.LimitedPar(runtime.GOMAXPROCS(0), pipeline.Receive(func(_ int, data interface{}) interface{} {
			region := data.(*assemblyRegion)
			return hc.callRegion(input.Header, region).makeVariantBlock(region.contig)
		})),
		pipeline.StrictOrd(hc.makeVariantCombinerFilter(&finalVariant, &finalVariantOK)),
		pipeline.LimitedPar(runtime.GOMAXPROCS(0), pipeline.Receive(func(_ int, data interface{}) interface{} {
			variants := data.(fullVariants)
			records := make([][]byte, 0, len(variants))
			var buf []byte
			for _, variant := range variants {
				buf = variant.Format(buf)
				records = append(records, append([]byte(nil), buf...))
				buf = buf[:0]
			}
			return records
		})),
		pipeline.StrictOrd(pipeline.ReceiveAndFinalize(func(_ int, data interface{}) interface{} {
			records := data.([][]byte)
			for _, record := range records {
				_, _ = vcfFile.Write(record)
			}
			return nil
		}, func() {
			if finalVariantOK {
				if _, err := vcfFile.Write(finalVariant.Format(nil)); err != nil {
					assemblyRegionPipeline.SetErr(err)
				}
			}
		})),
	)
	internal.RunPipeline(&variantCallPipeline)
	waitForProfiles.Wait()
	waitForAssemblyRegions.Wait()
	if err := assemblyRegionPipeline.Err(); err != nil {
		log.Panic(err)
	}
}
