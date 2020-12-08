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
	"path"
	"sort"

	"github.com/exascience/elprep/v5/internal"
	"github.com/exascience/elprep/v5/vcf"

	"github.com/exascience/elprep/v5/sam"
	"github.com/exascience/pargo/parallel"
)

func isGoodCigar(cigars []sam.CigarOperation) bool {
	if len(cigars) == 0 {
		return false
	}
	for _, op := range cigars {
		if op.Length == 0 {
			return false
		}
	}
	index := 0
	switch cigars[index].Operation {
	case 'H':
		if index++; index == len(cigars) {
			return false
		}
		switch cigars[index].Operation {
		case 'S':
			if index++; index == len(cigars) {
				return false
			}
		}
	case 'S', 'P':
		if index++; index == len(cigars) {
			return false
		}
	}
	switch cigars[index].Operation {
	case 'M', '=', 'X', 'N':
		index++
	case 'I':
		if index++; index < len(cigars) {
			switch cigars[index].Operation {
			case 'I', 'D', 'S', 'H':
				return false
			}
		}
	default:
		return false
	}
	for index < len(cigars) {
		switch op := cigars[index].Operation; op {
		case 'M', '=', 'X', 'N':
			index++
		case 'I', 'D':
			if index++; index < len(cigars) {
				switch cigars[index].Operation {
				case 'I', 'D', 'S', 'H':
					return false
				}
			} else if op == 'D' {
				return false
			}
		case 'P':
			if index++; index < len(cigars) {
				switch cigars[index].Operation {
				case 'P', 'S', 'H':
					return false
				}
			} else {
				return false
			}
		case 'S':
			if index++; index < len(cigars) {
				if cigars[index].Operation != 'H' {
					return false
				}
				index++
			}
			return index == len(cigars)
		case 'H':
			return index+1 == len(cigars)
		default:
			return false
		}
	}
	return true
}

// HaplotypeCallAln filters out the reads that
// the haplotypecaller cannot process.
func HaplotypeCallAln(hdr *sam.Header) sam.AlignmentFilter {
	return func(aln *sam.Alignment) bool {
		// this is not the same as recalibrateAln
		// we could generalize the two functions, but meh
		_, found := aln.TAGS.Get(sr)
		if found {
			return false
		}
		if aln.FlagSome(sam.Secondary | sam.Duplicate | sam.QCFailed) {
			return false
		}
		refLength := sam.ReferenceLengthFromCigar(aln.CIGAR)
		if refLength == 0 {
			return false
		}
		result := !isStrictUnmapped(aln) &&
			aln.POS > 0 && refLength >= 0 &&
			alignmentAgreesWithHeader(hdr, aln) &&
			aln.SEQ.Len() == int(sam.ReadLengthFromCigar(aln.CIGAR)) &&
			aln.MAPQ >= 20 &&
			aln.MAPQ != 255 &&
			aln.RG() != nil &&
			aln.SEQ.Len() == len(aln.QUAL) &&
			aln.SEQ.Len() > 0 &&
			isGoodCigar(aln.CIGAR) &&
			!cigarContainsN(aln.CIGAR)
		aln.TAGS = nil
		return result
	}
}

// running average is used to compute high-quality soft clips means
type runningAverage struct {
	observations float64
	mean         float64
}

func (avg runningAverage) add(observation float64) runningAverage {
	var newAvg runningAverage
	newAvg.observations = avg.observations + 1
	newAvg.mean = avg.mean + (observation-avg.mean)/newAvg.observations
	return newAvg
}

const (
	log10One      = 0.0
	log10Ploidy   = 0.3010299956639812
	log10OneThird = -0.47712125471966244
)

const (
	jacobianLogTableMaxTolerance = 8
	jacobianLogTableStep         = 0.0001
	jacobianLogTableInvStep      = 1 / jacobianLogTableStep
)

// some math stuff
// we use tables computed by gatk to ensure correctness
func jacobianLog(difference float64) float64 {
	return jacobianLogTable[int(math.Round(difference*jacobianLogTableInvStep))]
}

func approximateLog10SumLog10(a, b float64) float64 {
	if a > b {
		a, b = b, a
	}
	if math.IsInf(a, -1) {
		return b
	}
	if diff := b - a; diff < jacobianLogTableMaxTolerance {
		return b + jacobianLog(diff)
	}
	return b
}

func alnSlice(alns []*sam.Alignment, regionStart, regionEnd, maxReferenceLength int32) (result []*sam.Alignment, firstIndex int) {
	lowestPOS := regionStart - maxReferenceLength + 1
	j := sort.Search(len(alns), func(index int) bool {
		return alns[index].POS >= lowestPOS
	})
	for i := j; i < len(alns); i++ {
		aln := alns[i]
		if aln.POS > regionEnd {
			break
		}
		if aln.End() >= regionStart {
			if len(result) == 0 {
				firstIndex = i
			}
			result = append(result, aln)
		}
	}
	return
}

// some math stuff
var (
	ln10              = math.Log(10)
	log1mexpThreshold = math.Log(0.5)
)

func log1mexp(a float64) float64 {
	if a > 0 {
		return math.NaN()
	}
	if a == 0 {
		return math.Inf(-1)
	}
	if a < log1mexpThreshold {
		return math.Log1p(-math.Exp(a))
	}
	return math.Log(-math.Expm1(a))
}

func log10OneMinusPow10(a float64) float64 {
	if a > 0 {
		return math.NaN()
	}
	if a == 0 {
		return math.Inf(-1)
	}
	b := a * ln10
	return log1mexp(b) / ln10
}

func (hc *HaplotypeCaller) downsample(alns []*sam.Alignment) (result []*sam.Alignment) {
	if hc.maxReadsPerAlignmentStart < 1 {
		return alns
	}
	pos := int32(1)
	total := int32(0)
	var cur []*sam.Alignment
	var insert int
	for _, aln := range alns {
		if aln.POS == pos {
			total++
			if total <= hc.maxReadsPerAlignmentStart {
				cur = append(cur, aln)
			} else {
				slot := hc.random.Int31n(total)
				if slot < hc.maxReadsPerAlignmentStart {
					cur[slot] = aln
				}
			}
		} else {
			copy(alns[insert:], cur)
			insert += len(cur)
			pos = aln.POS
			total = 1
			cur = cur[:0]
			cur = append(cur, aln)
		}
	}
	copy(alns[insert:], cur)
	newEnd := insert + len(cur)
	for i := newEnd; i < len(alns); i++ {
		alns[i] = nil
	}
	return alns[:newEnd]
}

func (hc *HaplotypeCaller) maxReferenceLength(alns []*sam.Alignment) int32 {
	return int32(parallel.RangeReduceInt(0, len(alns), 0,
		func(low, high int) int {
			var maxReferenceLength int32
			for i := low; i < high; i++ {
				if length := sam.ReferenceLengthFromCigar(alns[i].CIGAR); length > maxReferenceLength {
					maxReferenceLength = length
				}
			}
			return int(maxReferenceLength)
		},
		func(x, y int) int {
			if x > y {
				return x
			}
			return y
		}))
}

// FilterReadsBySampleName filters out reads that do not belong to exactly one sample.
// If *sampleName != "", ensure only reads are passed through that match this sample name.
// If *sampleName == "", check that the header has only one (or no) sample name, and don't filter reads.
// In the latter case, after use of this filter in a pipeline, *sampleName will contain the unique sample
// name detected from the header, if any.
func FilterReadsBySampleName(sampleName *string) sam.Filter {
	return func(hdr *sam.Header) sam.AlignmentFilter {
		if *sampleName != "" {
			var validRGs []string
			for _, rg := range hdr.RG {
				if rg["SM"] == *sampleName {
					id, ok := rg["ID"]
					if !ok {
						log.Panicf("Unexpected read group without an ID for sample name %v.", *sampleName)
					}
					validRGs = append(validRGs, id)
				}
			}
			if len(validRGs) == 0 {
				log.Panicf("No read group available with requested sample name %v.", *sampleName)
			}
			return func(aln *sam.Alignment) bool {
				rg := aln.RG()
				if rg == nil {
					return false
				}
				rgn := rg.(string)
				for _, validRG := range validRGs {
					if rgn == validRG {
						return true
					}
				}
				return false
			}
		}
		if len(hdr.RG) > 0 {
			sm, ok := hdr.RG[0]["SM"]
			if ok {
				*sampleName = sm
				for i := 1; i < len(hdr.RG); i++ {
					sm, ok := hdr.RG[i]["SM"]
					if !ok {
						log.Panicf("Unexpected read group %v without a sample name; first sample name is %v.", hdr.RG[i]["ID"], *sampleName)
					}
					if sm != *sampleName {
						log.Panicf("Multiple sample names present. Ensure to request a particular sample name.")
					}
				}
			} else {
				for i := 1; i < len(hdr.RG); i++ {
					sm, ok := hdr.RG[i]["SM"]
					if ok {
						log.Panicf("Unexpected read group %v without a sample name; some sample name is %v.", hdr.RG[0]["ID"], sm)
					}
				}
			}
		}
		return nil
	}
}

func readOverlapsRegion(aln *sam.Alignment, regionStart, regionEnd int32) bool {
	if aln.SEQ.Len() == 0 {
		return false
	}
	alnStart, alnEnd := aln.POS, aln.End()
	if alnStart > alnEnd {
		return false
	}
	return alnStart <= regionEnd && regionStart <= alnEnd
}

func forEachReadPair(alns []*sam.Alignment, f func(aln1, aln2 *sam.Alignment)) {
	m := make(map[string]*sam.Alignment)
	for _, aln2 := range alns {
		if !aln2.IsMultiple() || aln2.IsNextUnmapped() || aln2.PNEXT == 0 || aln2.PNEXT > aln2.End() {
			continue
		}
		if aln1, ok := m[aln2.QNAME]; ok {
			delete(m, aln2.QNAME)
			f(aln1, aln2)
		} else {
			m[aln2.QNAME] = aln2
		}
	}
}

const pcrSnvErrorRate = 1e-4

var (
	pcrSnvErrorQual     = math.Round(-10 * log10(pcrSnvErrorRate))
	halfPcrSnvErrorQual = uint8(int(pcrSnvErrorQual) / 2)
)

func cleanOverlappingReadPair(aln1, aln2 *sam.Alignment) {
	if aln1.RNAME != aln2.RNAME {
		return
	}
	soft1 := softStart(aln1)
	soft2 := softStart(aln2)
	if soft1 >= soft2 {
		soft1, soft2 = soft2, soft1
		aln1, aln2 = aln2, aln1
	}
	if aln1.End() < aln2.POS {
		return
	}
	readBases, fallsInsideOrJustBeforeDeletionOrSkippedRegion := computeReadCoordinateForReferenceCoordinate(aln1.CIGAR, soft1, int(aln2.POS))
	if readBases == -1 {
		log.Panic("Unexpected coordinate.")
	}
	if fallsInsideOrJustBeforeDeletionOrSkippedRegion {
		readBases++
	}
	if nofOverlappingBases := minInt(aln1.SEQ.Len()-readBases, aln2.SEQ.Len()); nofOverlappingBases > 0 {
		q1 := make([]byte, len(aln1.QUAL))
		copy(q1, aln1.QUAL)
		q2 := make([]byte, len(aln2.QUAL))
		copy(q2, aln2.QUAL)
		for index2 := 0; index2 < nofOverlappingBases; index2++ {
			index1 := readBases + index2
			base1 := aln1.SEQ.Base(index1)
			base2 := aln2.SEQ.Base(index2)
			if base1 == base2 {
				q1[index1] = minUint8(q1[index1], halfPcrSnvErrorQual)
				q2[index2] = minUint8(q2[index2], halfPcrSnvErrorQual)
			} else {
				q1[index1] = 0
				q2[index2] = 0
			}
		}
		aln1.QUAL = q1
		aln2.QUAL = q2
	}
}

const (
	readLengthFilterThreshold  = 10
	readQualityFilterThreshold = 20
)

func filterNonPassingReads(region *assemblyRegion) (removed []*sam.Alignment) {
	i := 0
	l := len(region.alns)
	for j := 0; j < l; j++ {
		aln := region.alns[i]
		if aln.SEQ.Len() < readLengthFilterThreshold ||
			aln.MAPQ < readQualityFilterThreshold ||
			(aln.IsMultiple() && !aln.IsNextUnmapped() && aln.RNEXT != "=" && aln.RNEXT != aln.RNAME) {
			region.alns = append(region.alns[:i], region.alns[i+1:]...)
			region.alns = append(region.alns, aln)
		} else {
			i++
		}
	}
	removed = region.alns[i:]
	region.alns = region.alns[:i]
	return
}

// generate the activity and region files
func printAssemblyRegions(regionFile, activityFile io.Writer, regions []*assemblyRegion) {
	for _, region := range regions {
		if regionFile != nil {
			fmt.Fprintf(regionFile, "%v\t%v\t%v\tend-marker\t0.00000\n", region.contig, region.start-1, region.start)
			var activity float64
			if region.isActive {
				activity = 1.0
			} else {
				activity = -1.0
			}
			fmt.Fprintf(regionFile, "%v\t%v\t%v\tsize=%v\t%.5f\n", region.contig, region.start-1, region.end, region.end-region.start+1, activity)
		}
		if activityFile != nil {
			for index, state := range region.supportingStates {
				if state > 1.0 {
					state = 1.0
				}
				fmt.Fprintf(activityFile, "%v\t%v\t%v\tstate\t%.5f\n", region.contig, int(region.start)-1+index, int(region.start)+index, state)
			}
			region.supportingStates = nil // we don't need this anymore after printing
		}
	}
}

// CombineVcfOutputs combines multiple VCF files.
// This is used in the sfm command.
func CombineVcfOutputs(vcfPath, vcfOutput string) {
	vcfPath, files := internal.Directory(vcfPath)
	sort.Strings(files)
	input := vcf.Open(path.Join(vcfPath, files[0]))
	header, _ := vcf.ParseHeader(input.Reader)
	output := vcf.Create(vcfOutput, "", -1)
	defer output.Close()
	header.Format(output.Writer)
	internal.Copy(output.Writer, input.Reader)
	input.Close()
	for _, filename := range files[1:] {
		nextInput := vcf.Open(path.Join(vcfPath, filename))
		vcf.SkipHeader(nextInput.Reader)
		internal.Copy(output.Writer, nextInput.Reader)
		nextInput.Close()
	}
}
