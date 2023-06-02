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
	"math"

	"github.com/exascience/elprep/v5/utils/nibbles"

	"github.com/exascience/elprep/v5/vcf"

	"github.com/exascience/elprep/v5/fasta"

	"github.com/bits-and-blooms/bitset"

	"github.com/exascience/elprep/v5/sam"
)

func getBasesAndBaseQualitiesAlignedOneToOne(aln *sam.Alignment) (sam.Sequence, []byte) {
	bases := aln.SEQ
	baseQualities := aln.QUAL
	for _, element := range aln.CIGAR {
		if op := element.Operation; op == 'I' || op == 'D' {
			var nofRefBases int
			for _, element := range aln.CIGAR {
				nofRefBases += int(cigarConsumesReferenceBasesOrS[element.Operation]) * int(element.Length)
			}
			paddedBases := nibbles.Make2(0, nofRefBases)
			paddedBaseQualities := make([]byte, 0, nofRefBases)
			var pos int32
			for _, element := range aln.CIGAR {
				if operatorConsumesReadBases.Contains(element.Operation) {
					end := pos + element.Length
					if operatorConsumesReferenceBases.Contains(element.Operation) {
						paddedBases = paddedBases.AppendSlice(nibbles.Nibbles(bases.Slice(int(pos), int(end))))
						paddedBaseQualities = append(paddedBaseQualities, baseQualities[pos:end]...)
					}
					pos = end
				} else if operatorConsumesReferenceBases.Contains(element.Operation) {
					for j := int32(0); j < element.Length; j++ {
						paddedBases = paddedBases.Append('-')
						paddedBaseQualities = append(paddedBaseQualities, 0)
					}
				}
			}
			return sam.Sequence(paddedBases), paddedBaseQualities
		}
	}
	return bases, baseQualities
}

const (
	aMask = 1 << iota
	cMask
	gMask
	tMask
	rMask = aMask | gMask
	yMask = cMask | tMask
	sMask = cMask | gMask
	wMask = aMask | tMask
	kMask = gMask | tMask
	mMask = aMask | cMask
	bMask = cMask | gMask | tMask
	dMask = aMask | gMask | tMask
	hMask = aMask | cMask | tMask
	vMask = aMask | cMask | gMask
	nMask = aMask | cMask | gMask | tMask
	xMask = 0
	uMask = tMask
)

var baseToMask = map[byte]int{
	'a': aMask, 'A': aMask,
	'c': cMask, 'C': cMask,
	'g': gMask, 'G': gMask,
	't': tMask, 'T': tMask,
	'r': rMask, 'R': rMask,
	'y': yMask, 'Y': yMask,
	's': sMask, 'S': sMask,
	'w': wMask, 'W': wMask,
	'k': kMask, 'K': kMask,
	'm': mMask, 'M': mMask,
	'b': bMask, 'B': bMask,
	'd': dMask, 'D': dMask,
	'h': hMask, 'H': hMask,
	'v': vMask, 'V': vMask,
	'n': nMask, 'N': nMask,
	'x': xMask, 'X': xMask,
	'u': uMask, 'U': uMask,
}

func nucleotideIntersect(a, b byte) bool {
	return (baseToMask[a] & baseToMask[b]) != 0
}

func calculateBaselineMismatchQualities(readBases sam.Sequence, readQualities []byte, readStart int32, ref []byte, refIndex, regionPaddedEnd int32) []int32 {
	n := minInt32(int32(readBases.Len())-readStart, regionPaddedEnd-refIndex)
	results := make([]int32, n)
	var sum int32
	for i := n - 1; i >= 0; i-- {
		readBase := readBases.Base(int(readStart + i))
		refBase := ref[refIndex+i]
		if !nucleotideIntersect(readBase, refBase) && readBase != '-' {
			sum += int32(readQualities[readStart+i])
		}
		results[i] = sum
	}
	return results
}

func traverseEndOfReadForIndelMismatches(
	informativeBases *bitset.BitSet,
	readStart int32,
	readBases sam.Sequence,
	readQualities []byte,
	lastReadBaseToMarkAsIndelRelevent int32,
	secondaryReadBreakPosition int32,
	refIndex int32,
	regionPaddedEnd int32,
	ref []byte,
	baselineMismatchSums []int32,
	insertionLength, deletionLength int32,
) {
	globalMismatchCostForReadAlignedToReference := baselineMismatchSums[0]
	var baseQualitySum int32
	nofBasesToDirectlyCompare := minInt32(
		int32(readBases.Len())-readStart-insertionLength,
		regionPaddedEnd-refIndex-deletionLength,
	)
	for readOffset, refOffset := nofBasesToDirectlyCompare+insertionLength-1, nofBasesToDirectlyCompare+deletionLength-1; readOffset >= 0 && refOffset >= 0; readOffset, refOffset = readOffset-1, refOffset-1 {
		readBase := readBases.Base(int(readStart + readOffset))
		refBase := ref[refIndex+refOffset]
		if !nucleotideIntersect(readBase, refBase) && readBase != '-' {
			baseQualitySum += int32(readQualities[readStart+readOffset])
			if baseQualitySum > globalMismatchCostForReadAlignedToReference {
				break
			}
		}
		siteOfRealComparisonPoint := minInt32(readOffset, refOffset)
		if readBases.Base(int(readStart+siteOfRealComparisonPoint)) != '-' &&
			readStart+siteOfRealComparisonPoint < lastReadBaseToMarkAsIndelRelevent &&
			readStart+siteOfRealComparisonPoint <= secondaryReadBreakPosition &&
			baselineMismatchSums[siteOfRealComparisonPoint] >= baseQualitySum {
			informativeBases.Set(uint(readStart + siteOfRealComparisonPoint))
		}
	}
}

func (hc *HaplotypeCaller) readHasNoPlausibleIndelsOfMaxIndelSize(cache map[*sam.Alignment]*bitset.BitSet, aln *sam.Alignment, readStart int32, ref []byte, refIndex, regionPaddedEnd int32) bool {
	if cached, ok := cache[aln]; ok {
		return cached.Test(uint(readStart))
	}
	readLength := aln.SEQ.Len()
	informativeBases := bitset.New(uint(readLength))
	if int32(readLength)-readStart < hc.indelSizeToEliminateInRefModel ||
		regionPaddedEnd-refIndex < hc.indelSizeToEliminateInRefModel {
		cache[aln] = informativeBases
		return false
	}
	// todo: see https://github.com/broadinstitute/gatk/issues/5646 for possible changes in the future
	// note: somebody actually seems to be working on this, se wo have to pay attention
	secondaryReadBreakPosition := int32(readLength) - hc.indelSizeToEliminateInRefModel
	readBases, readQualities := getBasesAndBaseQualitiesAlignedOneToOne(aln)
	if int32(readBases.Len())-readStart <= hc.indelSizeToEliminateInRefModel {
		cache[aln] = informativeBases
		return false
	}
	lastReadBaseToMarkAsIndelRelevant := regionPaddedEnd - refIndex + readStart + 1
	referenceWasShorter := int32(readBases.Len()) >= lastReadBaseToMarkAsIndelRelevant
	if !referenceWasShorter {
		lastReadBaseToMarkAsIndelRelevant = int32(readBases.Len()) - hc.indelSizeToEliminateInRefModel
	}
	baselineMismatchSums := calculateBaselineMismatchQualities(readBases, readQualities, readStart, ref, refIndex, regionPaddedEnd)
	for indelSize := int32(1); indelSize <= hc.indelSizeToEliminateInRefModel; indelSize++ {
		traverseEndOfReadForIndelMismatches(
			informativeBases,
			readStart,
			readBases,
			readQualities,
			lastReadBaseToMarkAsIndelRelevant,
			secondaryReadBreakPosition,
			refIndex,
			regionPaddedEnd,
			ref,
			baselineMismatchSums,
			0, indelSize,
		)
		traverseEndOfReadForIndelMismatches(
			informativeBases,
			readStart,
			readBases,
			readQualities,
			lastReadBaseToMarkAsIndelRelevant,
			secondaryReadBreakPosition,
			refIndex,
			regionPaddedEnd,
			ref,
			baselineMismatchSums,
			indelSize, 0,
		)
	}
	if lastReadBaseToMarkAsIndelRelevant <= secondaryReadBreakPosition {
		for i := uint(0); i < uint(lastReadBaseToMarkAsIndelRelevant); i++ {
			informativeBases.Flip(i)
		}
		if referenceWasShorter {
			informativeBases.Clear(uint(lastReadBaseToMarkAsIndelRelevant - 1))
		}
	} else {
		for i, j := uint(0), uint(secondaryReadBreakPosition+1); i < j; i++ {
			informativeBases.Flip(i)
		}
	}
	cache[aln] = informativeBases
	return informativeBases.Test(uint(readStart))
}

const maxIndelInformativeReads = 40

var indelPLs [maxIndelInformativeReads + 1][3]float64

func init() {
	log10_0 := log10(0)
	log10_1 := log10(1)
	log10_2 := log10(2)
	denominator := -log10_2
	indelErrorRate := -4.5
	indelQual := byte(math.Round(indelErrorRate * -10))
	noIndelLikeliHood := qualToProbLog10[indelQual]
	indelLikeliHood := float64(indelQual) / -10
	indelPLs[1] = [3]float64{
		noIndelLikeliHood,
		approximateLog10SumLog10(noIndelLikeliHood+log10_1, indelLikeliHood+log10_1) + denominator,
		approximateLog10SumLog10(noIndelLikeliHood+log10_0, indelLikeliHood+log10_2) + denominator,
	}
	for i := 2; i <= maxIndelInformativeReads; i++ {
		indelPLs[i][0] = indelPLs[i-1][0] + indelPLs[1][0]
		indelPLs[i][1] = indelPLs[i-1][1] + indelPLs[1][1]
		indelPLs[i][2] = indelPLs[i-1][2] + indelPLs[1][2]
	}
}

func normalizeFromLog10(likelihoods [3]float64) [3]float64 {
	maxValue := math.Max(likelihoods[0], math.Max(likelihoods[1], likelihoods[2]))
	normalized := [3]float64{
		math.Pow(10, likelihoods[0]-maxValue),
		math.Pow(10, likelihoods[1]-maxValue),
		math.Pow(10, likelihoods[2]-maxValue),
	}
	sum := normalized[0] + normalized[1] + normalized[2]
	normalized[0] /= sum
	normalized[1] /= sum
	normalized[2] /= sum
	return normalized
}

func getGQLog10FromLikelihoods(likelihoods [3]float64) float64 {
	qual := likelihoods[0] - math.Max(likelihoods[1], likelihoods[2])
	if qual < 0 {
		normalized := normalizeFromLog10(likelihoods)
		return log10(1 - normalized[0])
	} else {
		return -1 * qual
	}
}

func computeGQ(pls [3]int) int {
	// Bubble sort first pass
	if pls[0] > pls[1] {
		pls[0], pls[1] = pls[1], pls[0]
	}
	if pls[1] > pls[2] {
		pls[1], pls[2] = pls[2], pls[1]
	}
	// Bubble sort second pass
	if pls[0] > pls[1] {
		pls[0], pls[1] = pls[1], pls[0]
	}
	if pls[1] > pls[2] {
		pls[1], pls[2] = pls[2], pls[1]
	}
	// return second smallest minus smallest element
	return pls[1] - pls[0]
}

func (hc *HaplotypeCaller) calculateRefConfidence(variants variantSlice, region *assemblyRegion, alns []*sam.Alignment, calls []*vcf.Variant) variantSlice {
	ref := region.reference
	cache := make(map[*sam.Alignment]*bitset.BitSet, len(alns))
	paddedEnd := region.paddedEnd()
	sam.By(sam.CoordinateLess).ParallelStableSort(alns)
	forEachPileupIncludingEmpty(alns, region.start, region.end+1, func(p *pileup) {
		var overlappingSite *vcf.Variant
		for _, vc := range calls {
			if p.location >= vc.Pos && p.location <= vc.End() && (overlappingSite == nil || vc.Pos > overlappingSite.Pos) {
				overlappingSite = vc
			}
		}
		if overlappingSite != nil && overlappingSite.Pos == p.location {
			variants = variants.addFullVariant(overlappingSite)
			return // next pileup
		}
		// cf. ReferenceConfidenceModel::makeReferenceConfidenceVariantContext
		refConfidence := calculateGenotypeLikelihoodsOfRefVsAny(p, fasta.ToUpperAndN(ref[p.location-1]), 6, isAltAfterAssembly)
		// cf. ReferenceConfidenceModel::doIndelRefConfCalc
		for i := 1; i < len(refConfidence.genotypeLikelihoods); i++ {
			refConfidence.genotypeLikelihoods[i] = math.Min(refConfidence.genotypeLikelihoods[0], refConfidence.genotypeLikelihoods[i])
		}
		// cf. ReferenceConfidenceModel::calcNReadsWithNoPlausibleIndelsReads
		var nInformative int
		for _, element := range p.filteredElements {
			cigarElement := element.aln.CIGAR[element.cigarIndex]
			if cigarElement.Operation == 'D' {
				continue
			}
			if element.cigarOffset == cigarElement.Length-1 {
				nextCigarElement, ok := element.nextOp()
				if op := nextCigarElement.Operation; ok && (op == 'D' || op == 'I') {
					continue
				}
			}
			offset := cigarConsumesReferenceBasesOrS[cigarElement.Operation] * element.cigarOffset
			for _, el := range element.aln.CIGAR[:element.cigarIndex] {
				offset += cigarConsumesReferenceBasesOrS[el.Operation] * el.Length
			}
			if hc.readHasNoPlausibleIndelsOfMaxIndelSize(cache, element.aln, offset, ref, p.location-1, paddedEnd) {
				if nInformative++; nInformative > maxIndelInformativeReads {
					nInformative = maxIndelInformativeReads
					break
				}
			}
		}
		indelGLs := indelPLs[nInformative]
		gqSnpGLs := getGQLog10FromLikelihoods(refConfidence.genotypeLikelihoods)
		gqIndelGLs := getGQLog10FromLikelihoods(indelGLs)
		var leastConfidenceGLs [3]float64
		if gqIndelGLs > gqSnpGLs {
			leastConfidenceGLs = indelGLs
		} else {
			leastConfidenceGLs = refConfidence.genotypeLikelihoods
		}
		adjust := math.Max(leastConfidenceGLs[0], math.Max(leastConfidenceGLs[1], leastConfidenceGLs[2]))
		leastConfidencePLs := [3]int{
			int(math.Round(math.Min(-10*(leastConfidenceGLs[0]-adjust), math.MaxInt32))),
			int(math.Round(math.Min(-10*(leastConfidenceGLs[1]-adjust), math.MaxInt32))),
			int(math.Round(math.Min(-10*(leastConfidenceGLs[2]-adjust), math.MaxInt32))),
		}
		gq := minInt(computeGQ(leastConfidencePLs), maxGenotypeQual)
		variants = variants.addReference(region.contig, ref, singleVariant{
			location: p.location,
			dp:       refConfidence.refDepth + refConfidence.nonRefDepth,
			ad:       [2]int{refConfidence.refDepth, refConfidence.nonRefDepth},
			pls:      leastConfidencePLs,
			gq:       gq,
		})
	})
	return variants
}

func (hc *HaplotypeCaller) referenceModelForNoVariation(variants variantSlice, region *assemblyRegion) variantSlice {
	filterNonPassingReads(region)
	return hc.calculateRefConfidence(variants, region, region.alns, nil)
}
