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
	"math"
	"sort"

	"github.com/exascience/pargo/parallel"

	"github.com/exascience/elprep/v5/sam"

	"github.com/exascience/elprep/v5/vcf"
)

type trimmingResult struct {
	needsTrimming                      bool
	leftFlankStart, leftFlankEnd       int32
	rightFlankStart, rightFlankEnd     int32
	extendedSpanStart, extendedSpanEnd int32
	callableSpanStart, callableSpanEnd int32
}

func (hc *HaplotypeCaller) trim(region *assemblyRegion, variationEvents map[int32]*vcf.Variant) trimmingResult {
	if len(variationEvents) == 0 {
		return trimmingResult{
			needsTrimming:  false,
			leftFlankStart: region.start,
			leftFlankEnd:   region.end,
		}
	}

	variantSpanStart := int32(math.MaxInt32)
	variantSpanEnd := int32(math.MinInt32)
	var withinActiveRegion int
	foundNonSnp := false
	for _, vc := range variationEvents {
		if vc.Pos <= region.end {
			if end := vc.End(); end >= region.start {
				if !foundNonSnp {
					if len(vc.Ref) != 1 {
						foundNonSnp = true
					} else {
						for _, a := range vc.Alt {
							if len(a) != 1 {
								foundNonSnp = true
								break
							}
						}
					}
				}
				if vc.Pos < variantSpanStart {
					variantSpanStart = vc.Pos
				}
				if end > variantSpanEnd {
					variantSpanEnd = end
				}
				withinActiveRegion++
			}
		}
	}

	if withinActiveRegion == 0 {
		return trimmingResult{
			needsTrimming:  false,
			leftFlankStart: region.start,
			leftFlankEnd:   region.end,
		}
	}

	var padding int32
	if foundNonSnp {
		padding = 150
	} else {
		padding = 20
	}

	maximumSpanStart := maxInt32(region.start-25, 1)
	maximumSpanEnd := minInt32(region.end+25, region.contigLength)
	idealSpanStart := maxInt32(variantSpanStart-padding, 1)
	idealSpanEnd := minInt32(variantSpanEnd+padding, region.contigLength)
	finalSpanStart := minInt32(maxInt32(maximumSpanStart, idealSpanStart), variantSpanStart)
	finalSpanEnd := maxInt32(minInt32(maximumSpanEnd, idealSpanEnd), variantSpanEnd)

	var callableSpanStart, callableSpanEnd int32
	if hc.confidenceMode == none {
		callableSpanStart = variantSpanStart
		callableSpanEnd = variantSpanEnd
	} else {
		callableSpanStart = maxInt32(variantSpanStart, region.start)
		callableSpanEnd = minInt32(variantSpanEnd, region.end)
	}

	var leftFlankStart, leftFlankEnd, rightFlankStart, rightFlankEnd int32

	if region.start < callableSpanStart {
		leftFlankStart = region.start
		leftFlankEnd = callableSpanStart - 1
	}
	if region.end > callableSpanEnd {
		rightFlankStart = callableSpanEnd + 1
		rightFlankEnd = region.end
	}

	return trimmingResult{
		needsTrimming:     true,
		leftFlankStart:    leftFlankStart,
		leftFlankEnd:      leftFlankEnd,
		rightFlankStart:   rightFlankStart,
		rightFlankEnd:     rightFlankEnd,
		extendedSpanStart: finalSpanStart,
		extendedSpanEnd:   finalSpanEnd,
		callableSpanStart: variantSpanStart,
		callableSpanEnd:   variantSpanEnd,
	}
}

func trimRegion(region *assemblyRegion, spanStart, spanEnd, extendedSpanStart, extendedSpanEnd int32) *assemblyRegion {
	subActiveStart := maxInt32(region.start, spanStart)
	subActiveEnd := minInt32(region.end, spanEnd)
	requiredOnRight := maxInt32(extendedSpanEnd-subActiveEnd, 0)
	requiredOnLeft := maxInt32(subActiveStart-extendedSpanStart, 0)
	requiredExtension := minInt32(maxInt32(requiredOnLeft, requiredOnRight), region.extension)
	result := assemblyRegion{
		contig:       region.contig,
		reference:    region.reference,
		start:        subActiveStart,
		end:          subActiveEnd,
		extension:    requiredExtension,
		contigLength: region.contigLength,
		isActive:     region.isActive,
		deletions:    region.deletions,
	}
	extendedLocStart := result.paddedStart()
	extendedLocEnd := result.paddedEnd()
	trimmedReads := make([]*sam.Alignment, 0, len(region.alns))
	newAln := new(sam.Alignment)
	for _, aln := range region.alns {
		*newAln = *aln
		hardClipToRegion(newAln, int(extendedLocStart), int(extendedLocEnd))
		if readOverlapsRegion(newAln, extendedLocStart, extendedLocEnd) {
			trimmedReads = append(trimmedReads, newAln)
			newAln = new(sam.Alignment)
		}
	}
	sam.By(sam.CoordinateLess).ParallelStableSort(trimmedReads)
	result.alns = trimmedReads
	return &result
}

func trimRegion1(region *assemblyRegion, spanStart, spanEnd, extension int32) *assemblyRegion {
	extendedStart := maxInt32(1, spanStart-extension)
	extendedEnd := minInt32(spanEnd+extension, region.contigLength)
	return trimRegion(region, spanStart, spanEnd, extendedStart, extendedEnd)
}

func addCigarElement(newCigar []sam.CigarOperation, pos, start, end int32, ce sam.CigarOperation) ([]sam.CigarOperation, int32) {
	length := minInt32(pos+ce.Length-1, end) - maxInt32(pos, start) + 1
	if length > 0 {
		newCigar = append(newCigar, sam.CigarOperation{length, ce.Operation})
	}
	return newCigar, pos + ce.Length
}

func (h *haplotype) trim(spanStart, spanEnd int32) *haplotype {
	newStart := spanStart - h.location
	newEnd := spanEnd - h.location

	refPos := int32(0)
	basesPos := int32(0)
	basesStart := int32(-1)
	basesStop := int32(-1)

getBasesCoveringRefInterval:
	for _, ce := range h.cigar {
		switch ce.Operation {
		case 'I':
			basesPos += ce.Length
		case 'M', 'X', '=':
			if newStart >= refPos && newStart < refPos+ce.Length {
				basesStart = basesPos + newStart - refPos
			}
			if newEnd >= refPos && newEnd < refPos+ce.Length {
				basesStop = basesPos + newEnd - refPos
				break getBasesCoveringRefInterval
			}
			refPos += ce.Length
			basesPos += ce.Length
		case 'D':
			if (newStart >= refPos && newStart < refPos+ce.Length) ||
				(newEnd >= refPos && newEnd < refPos+ce.Length) {
				return nil
			}
			refPos += ce.Length
		}
	}
	newBases := h.bases[basesStart : basesStop+1]

	newCigar := make([]sam.CigarOperation, 0, len(h.cigar))
	pos := int32(0)
	for _, ce := range h.cigar {
		if pos > newEnd {
			break
		}
		switch ce.Operation {
		case 'M', 'X', '=', 'D':
			newCigar, pos = addCigarElement(newCigar, pos, newStart, newEnd, ce)
		case 'S', 'I':
			if pos >= newStart {
				newCigar = append(newCigar, ce)
			}
		}
	}
	if ce := newCigar[0]; ce.Operation == 'I' || ce.Operation == 'D' {
		return nil
	}
	if ce := newCigar[len(newCigar)-1]; ce.Operation == 'I' || ce.Operation == 'D' {
		return nil
	}
	for i := 1; i < len(newCigar); i++ {
		if newCigar[i-1].Operation == newCigar[i].Operation {
			newCigar[i-1].Length += newCigar[i].Length
			newCigar = append(newCigar[:i], newCigar[i+1:]...)
		} else {
			i++
		}
	}
	return &haplotype{
		isRef:    h.isRef,
		bases:    newBases,
		location: spanStart,
		cigar:    newCigar,
		score:    h.score,
	}
}

func (hc *HaplotypeCaller) callRegion(hdr *sam.Header, region *assemblyRegion) variantSlice {
	variants := hc.makeVariantSlice()

	if !region.isActive || len(region.alns) == 0 {
		region.deletions.noDeletions()
		if hc.confidenceMode == none {
			return variants
		}
		hc.finalizeAssemblyRegion(region)
		return hc.referenceModelForNoVariation(variants, region)
	}
	hc.finalizeAssemblyRegion(region)

	haplotypes := hc.assembleReads(region)

	variationEvents := parallel.RangeReduce(0, len(haplotypes), 0, func(low, high int) interface{} {
		variationEvents := make(map[int32]*vcf.Variant)
		for i := low; i < high; i++ {
			events := makeEventMap(fmt.Sprintf("HC%d", i), region.contig, haplotypes[i], region.reference, nil)
			haplotypes[i].events = events
			for _, vc := range events {
				variationEvents[vc.Pos] = vc
			}
		}
		return variationEvents
	}, func(left, right interface{}) interface{} {
		l := left.(map[int32]*vcf.Variant)
		r := right.(map[int32]*vcf.Variant)
		// note: don't swap l, r based on map sizes, because order matters!
		for key, value := range l {
			r[key] = value
		}
		return r
	}).(map[int32]*vcf.Variant)

	trimmingResult := hc.trim(region, variationEvents)

	if !trimmingResult.needsTrimming {
		region.deletions.noDeletions()
		if hc.confidenceMode == none {
			return variants
		}
		return hc.referenceModelForNoVariation(variants, region)
	}

	var regionForGenotyping *assemblyRegion
	if hc.confidenceMode == none {
		regionForGenotyping = trimRegion(region,
			trimmingResult.extendedSpanStart, trimmingResult.extendedSpanEnd,
			trimmingResult.extendedSpanStart, trimmingResult.extendedSpanEnd,
		)
	} else {
		regionForGenotyping = trimRegion(region,
			trimmingResult.callableSpanStart, trimmingResult.callableSpanEnd,
			trimmingResult.extendedSpanStart, trimmingResult.extendedSpanEnd,
		)
	}
	trimmedHaplotypes := make([]*haplotype, 0, len(haplotypes))
trimHaplotypes:
	for _, h := range haplotypes {
		if trimmedHaplotype := h.trim(regionForGenotyping.paddedStart(), regionForGenotyping.paddedEnd()); trimmedHaplotype != nil {
			for i, th := range trimmedHaplotypes {
				if trimmedHaplotype.bases == th.bases {
					if trimmedHaplotype.isRef {
						trimmedHaplotypes[i] = trimmedHaplotype
					}
					continue trimHaplotypes
				}
			}
			trimmedHaplotypes = append(trimmedHaplotypes, trimmedHaplotype)
		}
	}
	sort.SliceStable(trimmedHaplotypes, func(i, j int) bool {
		bi := trimmedHaplotypes[i].bases
		bj := trimmedHaplotypes[j].bases
		if len(bi) < len(bj) {
			return true
		} else if len(bi) > len(bj) {
			return false
		} else {
			return bi < bj
		}
	})
	haplotypes = trimmedHaplotypes
	refHaplotype := -1
	variationPresent := false
	for i, h := range trimmedHaplotypes {
		if h.isRef {
			refHaplotype = i
			if variationPresent {
				break
			}
		} else {
			variationPresent = true
			if refHaplotype >= 0 {
				break
			}
		}
	}

	if !variationPresent {
		region.deletions.noDeletions()
		if hc.confidenceMode == none {
			return variants
		}
		return hc.referenceModelForNoVariation(variants, region)
	}

	for i := 0; i < len(regionForGenotyping.alns); {
		if aln := regionForGenotyping.alns[i]; aln.SEQ.Len() < 10 {
			regionForGenotyping.alns = append(regionForGenotyping.alns[:i], regionForGenotyping.alns[i+1:]...)
		} else {
			i++
		}
	}

	filteredReads := filterNonPassingReads(regionForGenotyping)

	if len(regionForGenotyping.alns) == 0 {
		region.deletions.noDeletions()
		if hc.confidenceMode == none {
			return variants
		}
		return hc.referenceModelForNoVariation(variants, region)
	}

	readLikelihoods := computeReadLikelihoods(haplotypes, regionForGenotyping.alns)
	realignReadsToTheirBestHaplotype(readLikelihoods, haplotypes)
	returnCalls, _ := hc.assignGenotypeLikelihoods(hdr, regionForGenotyping, filteredReads, haplotypes, readLikelihoods)

	if len(returnCalls) == 0 {
		if hc.confidenceMode == none {
			return variants
		}
		return hc.referenceModelForNoVariation(variants, region)
	}

	if hc.confidenceMode == none {
		parallel.Range(0, len(returnCalls), 0, func(low, high int) {
			for i := low; i < high; i++ {
				call := returnCalls[i]
				if value, ok := call.Info.Get(RAW_MQandDP); ok {
					raw := value.([]interface{})
					sum := raw[0].(int)
					depth := raw[1].(int)
					rms := math.Sqrt(float64(sum) / float64(depth))
					call.Info.Delete(RAW_MQandDP)
					call.Info.Set(MQ, formatf(rms, 2))
				}
				computeGenotypeFormat(call)
			}
		})
		return variants.addFullVariants(returnCalls)
	}

	parallel.Range(0, len(returnCalls), 0, func(low, high int) {
		for i := low; i < high; i++ {
			computeGenotypeFormat(returnCalls[i])
		}
	})

	if trimmingResult.leftFlankEnd >= trimmingResult.leftFlankStart {
		variants = hc.referenceModelForNoVariation(variants, trimRegion1(region, trimmingResult.leftFlankStart, trimmingResult.leftFlankEnd, region.extension))
	}

	variants = hc.calculateRefConfidence(variants, regionForGenotyping, readLikelihoods.alns, returnCalls)

	if trimmingResult.rightFlankEnd >= trimmingResult.rightFlankStart {
		variants = hc.referenceModelForNoVariation(variants, trimRegion1(region, trimmingResult.rightFlankStart, trimmingResult.rightFlankEnd, region.extension))
	}

	return variants
}
