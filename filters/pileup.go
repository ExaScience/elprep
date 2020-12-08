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
	"log"
	"math"

	"github.com/exascience/elprep/v5/sam"
)

const (
	highQualitySoftClipThreshold              = 28
	averageHighQualitySoftClipsBasesThreshold = 6
)

func countHighQualitySoftClips(aln *sam.Alignment) (result int32) {
	var alignPos int32
	for _, op := range aln.CIGAR {
		switch op.Operation {
		case 'S':
			for i := int32(0); i < op.Length; i++ {
				if aln.QUAL[alignPos+i] > highQualitySoftClipThreshold {
					result++
				}
			}
			fallthrough
		case 'M', 'I', '=', 'X':
			alignPos += op.Length
		}
	}
	return
}

/*
Pileup handling:
For each position in the reference, the algorithm bulids a pileup.
This means that bases that are only in the reads (clips, paddings, inserts)
are not considered.
*/

// A pileupElement represent one base in one read at a particular position
type pileupElement struct {
	// the read
	aln *sam.Alignment
	// the index within the read
	baseIndex int32
	// the cigar element on that position
	cigarIndex int32
	// the offset within that cigar element
	cigarOffset int32
	// the number of high-quality soft clips in the cigar string
	// (total number, not just the current position)
	nofHighQualitySoftClips int32
}

// the base in the read for the pileup element
func (p *pileupElement) base() byte {
	return p.aln.SEQ.Base(int(p.baseIndex))
}

// the quality of the base in the read for the pileup element
func (p *pileupElement) qual() byte {
	return p.aln.QUAL[p.baseIndex]
}

// the cigar operation that covers the base in the read for the pileup element
func (p *pileupElement) op() sam.CigarOperation {
	return p.aln.CIGAR[p.cigarIndex]
}

// the cigar operation that precedes the current cigar operation of the pileup element, if any
func (p *pileupElement) prevOp() (sam.CigarOperation, bool) {
	cigar := p.aln.CIGAR
	if p.cigarOffset > 0 {
		return cigar[p.cigarIndex], true
	}
	if p.cigarIndex > 0 {
		return cigar[p.cigarIndex-1], true
	}
	return sam.CigarOperation{}, false
}

// the cigar operation that precedes the current cigar operation of the pileup element, if any
// however, this method oly returns operations that match bases in the reference
// (so no clips, paddings, or inserts)
func (p *pileupElement) prevOnGenomeOp() (sam.CigarOperation, bool) {
	cigar := p.aln.CIGAR
	if p.cigarOffset > 0 {
		switch op := cigar[p.cigarIndex]; op.Operation {
		case 'M', '=', 'X', 'D':
			return op, true
		}
	}
	for i := p.cigarIndex - 1; i >= 0; i-- {
		switch op := cigar[i]; op.Operation {
		case 'M', '=', 'X', 'D':
			return op, true
		}
	}
	return sam.CigarOperation{}, false
}

// the cigar operation that succeeds the current cigar operation of the pileup element, if any
func (p *pileupElement) nextOp() (sam.CigarOperation, bool) {
	cigar := p.aln.CIGAR
	op := cigar[p.cigarIndex]
	if p.cigarOffset+1 < op.Length {
		return op, true
	}
	if nextCigarIndex := int(p.cigarIndex) + 1; nextCigarIndex < len(cigar) {
		return cigar[nextCigarIndex], true
	}
	return sam.CigarOperation{}, false
}

// the cigar operation that succeeds the current cigar operation of the pileup element, if any
// however, this method oly returns operations that match bases in the reference
// (so no clips, paddings, or inserts)
func (p *pileupElement) nextOnGenomeOp() (sam.CigarOperation, bool) {
	cigar := p.aln.CIGAR
	op := cigar[p.cigarIndex]
	if p.cigarOffset+1 < op.Length {
		switch op.Operation {
		case 'M', '=', 'X', 'D':
			return op, true
		}
	}
	for i := int(p.cigarIndex) + 1; i < len(cigar); i++ {
		switch op := cigar[i]; op.Operation {
		case 'M', '=', 'X', 'D':
			return op, true
		}
	}
	return sam.CigarOperation{}, false
}

// the first pileup element in the next cigar operation following the current pileup element
// however, ignore clips, paddings, and inserts
func (p *pileupElement) firstPileupElementInNextCigarOperation() bool {
	cigar := p.aln.CIGAR
	cigarLen := int32(len(cigar))
	for p.cigarIndex++; p.cigarIndex < cigarLen; p.cigarIndex++ {
		switch op := cigar[p.cigarIndex]; op.Operation {
		case 'H', 'P':
		case 'I', 'S':
			p.baseIndex += op.Length
		case 'D', 'N':
			p.cigarOffset = 0
			return true
		case 'M', '=', 'X':
			p.baseIndex++
			p.cigarOffset = 0
			return true
		default:
			log.Panicf("Invalid cigar operation %c for read %v.", op.Operation, p.aln.QNAME)
		}
	}
	return false
}

// first pileup element in the read
// however, ignore clips, paddings, and inserts
func firstPileupElement(aln *sam.Alignment) (*pileupElement, bool) {
	p := pileupElement{
		aln:                     aln,
		baseIndex:               -1,
		cigarIndex:              -1,
		cigarOffset:             0,
		nofHighQualitySoftClips: countHighQualitySoftClips(aln),
	}
	return &p, p.firstPileupElementInNextCigarOperation()
}

// next pileup element after the current one
// ignores clips, paddings, and inserts
func (p *pileupElement) nextPileupElement() bool {
	op := p.aln.CIGAR[p.cigarIndex]
	if p.cigarOffset++; p.cigarOffset < op.Length {
		switch op.Operation {
		case 'M', '=', 'X':
			p.baseIndex++
		}
		return true
	}
	return p.firstPileupElementInNextCigarOperation()
}

func firstPileupElementAtOrAbove(aln *sam.Alignment, location int32) (*pileupElement, int32) {
	p := pileupElement{
		aln:                     aln,
		baseIndex:               0,
		cigarIndex:              -1,
		cigarOffset:             0,
		nofHighQualitySoftClips: countHighQualitySoftClips(aln),
	}
	cigar := aln.CIGAR
	loc := aln.POS
	for i, op := range cigar {
		switch op.Operation {
		case 'H', 'P':
		case 'I', 'S':
			p.baseIndex += op.Length
		case 'D', 'N':
			if loc+op.Length > location {
				delta := maxInt32(location-loc, 0)
				p.baseIndex--
				p.cigarIndex = int32(i)
				p.cigarOffset = delta
				return &p, maxInt32(location, loc)
			}
			loc += op.Length
		case 'M', '=', 'X':
			if loc+op.Length > location {
				delta := maxInt32(location-loc, 0)
				p.baseIndex += delta
				p.cigarIndex = int32(i)
				p.cigarOffset = delta
				return &p, maxInt32(location, loc)
			}
			p.baseIndex += op.Length
			loc += op.Length
		default:
			log.Panicf("Invalid cigar operation %c for read %v.", op.Operation, aln.QNAME)
		}
	}
	return nil, -1
}

// a pileup contains all pileup elements at a particular reference location
type pileup struct {
	location                      int32
	allElements, filteredElements []pileupElement
}

// filter out the pileup elements that are only in the adaptor boundaries of the corresponding reads
// (we need to keep the reads in the original pileup because we may need them for later pileups)
func (p *pileup) filterAdaptors() {
	p.filteredElements = p.filteredElements[:0]
	for _, element := range p.allElements {
		if element.aln.TLEN > 100 {
			p.filteredElements = append(p.filteredElements, element)
			continue
		}
		boundary, _, wellDefinedSize := computeAdaptorBoundary(element.aln)
		if !wellDefinedSize {
			p.filteredElements = append(p.filteredElements, element)
			continue
		}
		if element.aln.IsReversed() {
			if int(p.location) > boundary {
				p.filteredElements = append(p.filteredElements, element)
				continue
			}
		} else {
			if int(p.location) < boundary {
				p.filteredElements = append(p.filteredElements, element)
				continue
			}
		}
	}
}

func (p *pileup) advanceToFirstElementAbove(alns []*sam.Alignment, low, high int32) (firstIndexAbove int, firstElementAbove *pileupElement, firstLocAbove int32) {
	firstIndexAbove = len(alns)
	firstLocAbove = high
	for i, aln := range alns {
		element, loc := firstPileupElementAtOrAbove(aln, low)
		if loc < 0 {
			continue
		} else if loc == low {
			p.allElements = append(p.allElements, *element)
		} else { // loc > low
			firstIndexAbove = i
			firstElementAbove = element
			firstLocAbove = loc
			return
		}
	}
	return
}

func (p *pileup) nextPileup() {
	for i := 0; i < len(p.allElements); {
		if p.allElements[i].nextPileupElement() {
			i++
		} else {
			p.allElements = append(p.allElements[:i], p.allElements[i+1:]...)
		}
	}
}

func forEachPileup(alns []*sam.Alignment, low, high int32, f func(p *pileup)) {
	if 1 >= high || low >= high {
		return
	}
	p := &pileup{location: low}
	firstIndexAbove, firstElementAbove, firstLocAbove := p.advanceToFirstElementAbove(alns, low, high)
	for ; len(p.allElements) > 0 && p.location < firstLocAbove; p.nextPileup() {
		if p.filterAdaptors(); len(p.filteredElements) > 0 {
			f(p)
		}
		if p.location++; p.location >= high {
			return
		}
	}
	if firstLocAbove >= high {
		return
	}
	p.location = firstLocAbove
	if firstElementAbove != nil {
		p.allElements = append(p.allElements, *firstElementAbove)
		firstElementAbove = nil
	}
	for i := firstIndexAbove + 1; i < len(alns); i++ {
		aln := alns[i]
		element, ok := firstPileupElement(aln)
		if !ok {
			continue
		}
		for ; len(p.allElements) > 0 && p.location < aln.POS; p.nextPileup() {
			if p.filterAdaptors(); len(p.filteredElements) > 0 {
				f(p)
			}
			if p.location++; p.location >= high {
				return
			}
		}
		if aln.POS >= high {
			return
		}
		p.location = aln.POS
		p.allElements = append(p.allElements, *element)
	}
	for ; len(p.allElements) > 0; p.nextPileup() {
		if p.filterAdaptors(); len(p.filteredElements) > 0 {
			f(p)
		}
		if p.location++; p.location >= high {
			return
		}
	}
}

func forEachPileupIncludingEmpty(alns []*sam.Alignment, low, high int32, f func(p *pileup)) {
	if 1 >= high || low >= high {
		return
	}
	p := &pileup{location: low}
	firstIndexAbove, firstElementAbove, firstLocAbove := p.advanceToFirstElementAbove(alns, low, high)
	for ; p.location < firstLocAbove; p.nextPileup() {
		p.filterAdaptors()
		f(p)
		if p.location++; p.location >= high {
			return
		}
	}
	if firstElementAbove != nil {
		p.allElements = append(p.allElements, *firstElementAbove)
		firstElementAbove = nil
	}
	for i := firstIndexAbove + 1; i < len(alns); i++ {
		aln := alns[i]
		element, ok := firstPileupElement(aln)
		if !ok {
			continue
		}
		for ; p.location < aln.POS; p.nextPileup() {
			p.filterAdaptors()
			f(p)
			if p.location++; p.location >= high {
				return
			}
		}
		p.allElements = append(p.allElements, *element)
	}
	for ; len(p.allElements) > 0; p.nextPileup() {
		p.filterAdaptors()
		f(p)
		if p.location++; p.location >= high {
			return
		}
	}
	p.filteredElements = nil
	for ; p.location < high; p.location++ {
		f(p)
	}
}

// is at least one of the cigar operations left or right to this pileup element a soft clip?
func isNextToSoftClip(element pileupElement) bool {
	if op, ok := element.prevOp(); ok {
		if op.Operation == 'S' {
			return true
		}
	}
	if op, ok := element.nextOp(); ok {
		if op.Operation == 'S' {
			return true
		}
	}
	return false
}

func isAltBeforeAssembly(element pileupElement, refBase byte) bool {
	if element.base() != refBase {
		return true
	}
	if element.op().Operation == 'D' {
		return true
	}
	if op, ok := element.prevOp(); ok {
		switch op.Operation {
		case 'I', 'S':
			return true
		}
	}
	if op, ok := element.nextOp(); ok {
		switch op.Operation {
		case 'I', 'S':
			return true
		}
	}
	if op, ok := element.prevOnGenomeOp(); ok {
		if op.Operation == 'D' {
			return true
		}
	}
	if op, ok := element.nextOnGenomeOp(); ok {
		if op.Operation == 'D' {
			return true
		}
	}
	return false
}

func isAltAfterAssembly(element pileupElement, refBase byte) bool {
	if element.base() != refBase {
		return true
	}
	if element.op().Operation == 'D' {
		return true
	}
	return false
}

type referenceConfidence struct {
	refDepth, nonRefDepth int
	genotypeLikelihoods   [3]float64
	hqSoftClipsMean       float64
}

// this is used to calculate the activity for a particular pileup
func calculateGenotypeLikelihoodsOfRefVsAny(p *pileup, refBase byte, minBaseQual byte, isAlt func(pileupElement, byte) bool) (result referenceConfidence) {
	var readCount float64
	var hqSoftClips runningAverage

	for _, element := range p.filteredElements {
		var qual byte
		if element.op().Operation == 'D' {
			qual = 30
		} else if qual = element.qual(); qual <= minBaseQual {
			continue
		}
		readCount++
		isAlt := isAlt(element, refBase)
		refLikelihood := qualToProbLog10[qual]
		nonRefLikelihood := float64(qual)/-10.0 + log10OneThird
		if isAlt {
			refLikelihood, nonRefLikelihood = nonRefLikelihood, refLikelihood
			result.nonRefDepth++
		} else {
			result.refDepth++
		}
		result.genotypeLikelihoods[0] += refLikelihood + log10Ploidy
		result.genotypeLikelihoods[1] += approximateLog10SumLog10(
			refLikelihood+log10One,
			nonRefLikelihood+log10One,
		)
		result.genotypeLikelihoods[2] += nonRefLikelihood + log10Ploidy
		if isAlt && isNextToSoftClip(element) {
			hqSoftClips = hqSoftClips.add(float64(element.nofHighQualitySoftClips))
		}
	}

	denominator := readCount * log10Ploidy
	result.genotypeLikelihoods[0] -= denominator
	result.genotypeLikelihoods[1] -= denominator
	result.genotypeLikelihoods[2] -= denominator
	result.hqSoftClipsMean = hqSoftClips.mean
	return
}

// determine if a pileup is active with regard to the base at the reference location
// also compute the high-quality soft clips mean for this location
func (hc *HaplotypeCaller) isActive(p *pileup, refBase byte) (isActiveProb float64, hqSoftClipsMean float64) {
	if len(p.filteredElements) == 0 {
		return
	}

	refConfidence := calculateGenotypeLikelihoodsOfRefVsAny(p, refBase, hc.minBaseQual, isAltBeforeAssembly)
	genotypeLikelihoods := refConfidence.genotypeLikelihoods
	hqSoftClipsMean = refConfidence.hqSoftClipsMean

	adjust := genotypeLikelihoods[0]
	if l := genotypeLikelihoods[1]; l > adjust {
		adjust = l
	}
	if l := genotypeLikelihoods[2]; l > adjust {
		adjust = l
	}
	for i := range genotypeLikelihoods {
		if adjusted := -10 * (genotypeLikelihoods[i] - adjust); adjusted > math.MaxInt32 {
			genotypeLikelihoods[i] = math.MaxInt32 / -10.0
		} else {
			genotypeLikelihoods[i] = math.Round(adjusted) / -10.0
		}
	}

	log10ACeq0Likelihood := genotypeLikelihoods[0]
	log10ACeq0Prior := hc.log10Priors[0]
	log10ACeq0Posterior := log10ACeq0Likelihood + log10ACeq0Prior

	for AC := 1; AC < len(hc.log10Priors); AC++ {
		if hc.log10Priors[AC]+genotypeLikelihoods[AC] > log10ACeq0Posterior {
			log10ACgt0Likelihood := approximateLog10SumLog10(genotypeLikelihoods[1], genotypeLikelihoods[2])
			log10ACgt0Posterior := log10ACgt0Likelihood + hc.log10ACgt0Prior
			log10PosteriorNormalizationConstant := approximateLog10SumLog10(log10ACeq0Posterior, log10ACgt0Posterior)
			normalizedLog10ACeq0Posterior := log10ACeq0Posterior - log10PosteriorNormalizationConstant
			if normalizedLog10ACeq0Posterior < hc.standardConfidenceForActivityByMin10 {
				isActiveProb = 1.0 - math.Pow(10.0, normalizedLog10ACeq0Posterior)
			}
			return
		}
	}

	return
}
