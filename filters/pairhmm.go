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
	"strings"
	"sync"

	"github.com/exascience/elprep/v5/sam"
	"github.com/exascience/pargo/parallel"
)

type float64Matrix struct {
	cols  int
	array []float64
}

func (m *float64Matrix) ensureSize(rows, cols int) {
	m.cols = cols
	totalSize := rows * cols
	if totalSize <= cap(m.array) {
		m.array = m.array[:totalSize]
		for i := range m.array {
			m.array[i] = 0
		}
	} else {
		m.array = make([]float64, totalSize)
	}
}

func (m *float64Matrix) rowView(row int) []float64 {
	offset := row * m.cols
	return m.array[offset : offset+m.cols]
}

type pairHMMMatrices struct {
	match, insertion, deletion float64Matrix
}

var pairHMMMatricesPool = sync.Pool{New: func() interface{} { return new(pairHMMMatrices) }}

func getPairHMMMatrices() *pairHMMMatrices {
	return pairHMMMatricesPool.Get().(*pairHMMMatrices)
}

func putPairHMMMatrices(p *pairHMMMatrices) {
	pairHMMMatricesPool.Put(p)
}

func (p *pairHMMMatrices) ensureSize(readBases, alleleBases int) {
	parallel.Do(
		func() { p.match.ensureSize(readBases, alleleBases) },
		func() { p.insertion.ensureSize(readBases, alleleBases) },
		func() { p.deletion.ensureSize(readBases, alleleBases) },
	)
}

func modifiedQuality(aln *sam.Alignment, index int) byte {
	qual := aln.QUAL[index]
	if qual > aln.MAPQ {
		qual = aln.MAPQ
	}
	if qual < 18 {
		return 6
	}
	return qual
}

func findNumberOfForwardRepetitions(repeatUnit, testString string) (nofRepeats int) {
	repeatLength := len(repeatUnit)
	for len(testString) >= repeatLength && strings.HasPrefix(testString, repeatUnit) {
		nofRepeats++
		testString = testString[repeatLength:]
	}
	return nofRepeats
}

func findNumberOfBackwardRepetitions(repeatUnit, testString string) (nofRepeats int) {
	repeatLength := len(repeatUnit)
	for len(testString) >= repeatLength && strings.HasSuffix(testString, repeatUnit) {
		nofRepeats++
		testString = testString[:len(testString)-repeatLength]
	}
	return nofRepeats
}

func findTandemRepeatUnits(readBases string, offset int) (bestRepeatUnit string, maxRL int) {
	offset1 := offset + 1
	var maxBW int
	bestBWRepeatUnit := readBases[offset:offset1]
	bwTestString := readBases[:offset1]
	for str := 1; str <= 8; str++ {
		repeatOffset := offset1 - str
		if repeatOffset < 0 {
			break
		}
		repeatUnit := readBases[repeatOffset:offset1]
		maxBW = findNumberOfBackwardRepetitions(repeatUnit, bwTestString)
		if maxBW > 1 {
			bestBWRepeatUnit = repeatUnit
			break
		}
	}
	bestRepeatUnit = bestBWRepeatUnit
	maxRL = maxBW

	if offset1 < len(readBases) {
		var maxFW int
		bestFWRepeatUnit := readBases[offset1 : offset1+1]
		fwTestString := readBases[offset1:]
		for str := 1; str <= 8; str++ {
			repeatOffset := offset1 + str
			if repeatOffset > len(readBases) {
				break
			}
			repeatUnit := readBases[offset1:repeatOffset]
			maxFW = findNumberOfForwardRepetitions(repeatUnit, fwTestString)
			if maxFW > 1 {
				bestFWRepeatUnit = repeatUnit
				break
			}
		}
		if bestFWRepeatUnit != bestBWRepeatUnit {
			testString := readBases[:offset1]
			maxBW = findNumberOfBackwardRepetitions(bestFWRepeatUnit, testString)
		}
		maxRL = maxFW + maxBW
		bestRepeatUnit = bestFWRepeatUnit
	}

	if maxRL > 20 {
		maxRL = 20
	}
	return bestRepeatUnit, maxRL
}

func matchProbs(alnBases string, index int) (matchToMatch, matchToIndel float64) {
	var repeatLength int
	if index == len(alnBases)-1 {
		repeatLength = 21
	} else {
		_, repeatLength = findTandemRepeatUnits(alnBases, index)
	}
	return matchToMatchProb[repeatLength], matchToIndelProb[repeatLength]
}

var (
	initialCondition      = math.Pow(2, 1020)
	initialConditionLog10 = log10(initialCondition)
	indelToIndel          = qualityToErrorProbability(10)
	indelToMatch          = 1 - indelToIndel
)

const globalReadMismappingRate = 45 / -10.0

type readLikelihoods struct {
	alns   []*sam.Alignment
	values map[*haplotype][]float64
}

// note: this may set some entries of alns to nil
// this is ok, because we deal with it properly afterwards
func computeReadLikelihoods(haplotypes []*haplotype, alns []*sam.Alignment) readLikelihoods {
	var maxReadLength, maxHaplotypeLength int
	parallel.Do(
		func() {
			maxReadLength = parallel.RangeReduceInt(0, len(alns), 0, func(low, high int) int {
				var max int
				for i := low; i < high; i++ {
					if l := alns[i].SEQ.Len(); l > max {
						max = l
					}
				}
				return max
			}, maxInt)
		},
		func() {
			maxHaplotypeLength = parallel.RangeReduceInt(0, len(haplotypes), 0, func(low, high int) int {
				var max int
				for i := low; i < high; i++ {
					if l := len(haplotypes[i].bases); l > max {
						max = l
					}
				}
				return max
			}, maxInt)
		},
	)
	result := readLikelihoods{
		alns:   alns,
		values: make(map[*haplotype][]float64, len(haplotypes)),
	}
	for _, haplotype := range haplotypes {
		result.values[haplotype] = make([]float64, len(alns))
	}

	parallel.Range(0, len(result.alns), len(result.alns), func(low, high int) {
		for readIndex := low; readIndex < high; readIndex++ {
			aln := result.alns[readIndex]
			alnBases := aln.SEQ.AsString()
			matchProbCache := make([][2]float64, len(aln.QUAL))
			parallel.Range(0, len(aln.QUAL), 0, func(low, high int) {
				for i := low; i < high; i++ {
					matchToMatch, matchToIndel := matchProbs(alnBases, i)
					matchProbCache[i] = [2]float64{matchToMatch, matchToIndel}
				}
			})
			parallel.Range(0, len(haplotypes), len(haplotypes), func(low, high int) {
				p := getPairHMMMatrices()
				defer putPairHMMMatrices(p)

				p.ensureSize(maxReadLength+1, maxHaplotypeLength+1)

				for haplotypeIndex := low; haplotypeIndex < high; haplotypeIndex++ {
					haplotype := haplotypes[haplotypeIndex]
					initialValue := initialCondition / float64(len(haplotype.bases))
					pDeletion0 := p.deletion.rowView(0)
					for j := 0; j <= maxHaplotypeLength; j++ {
						pDeletion0[j] = initialValue
					}
					for i := range aln.QUAL {
						x := alnBases[i]
						qual := modifiedQuality(aln, i)
						matchPrior := 1 - qualToErrorProb[qual]
						nonMatchPrior := qualToErrorProb[qual] / 3
						// matchToMatch, matchToIndel := matchProbs(alnBases, i)
						cachedMatchProbs := matchProbCache[i]
						matchToMatch, matchToIndel := cachedMatchProbs[0], cachedMatchProbs[1]

						// note: it's important to get the row views for performance
						pMatchI := p.match.rowView(i)
						pMatchI1 := p.match.rowView(i + 1)
						pInsertionI := p.insertion.rowView(i)
						pInsertionI1 := p.insertion.rowView(i + 1)
						pDeletionI := p.deletion.rowView(i)
						pDeletionI1 := p.deletion.rowView(i + 1)

						for j := 0; j < len(haplotype.bases); j++ {
							y := haplotype.bases[j]
							var prior float64
							if x == y || x == 'N' || y == 'N' {
								prior = matchPrior
							} else {
								prior = nonMatchPrior
							}
							pMatchI1[j+1] = prior * (pMatchI[j]*matchToMatch +
								pInsertionI[j]*indelToMatch +
								pDeletionI[j]*indelToMatch)
							pInsertionI1[j+1] = pMatchI[j+1]*matchToIndel + pInsertionI[j+1]*indelToIndel
							pDeletionI1[j+1] = pMatchI1[j]*matchToIndel + pDeletionI1[j]*indelToIndel
						}
					}
					var sum float64
					pMatchEnd := p.match.rowView(len(aln.QUAL))
					pInsertionEnd := p.insertion.rowView(len(aln.QUAL))
					for j := 1; j <= len(haplotype.bases); j++ {
						sum += pMatchEnd[j] + pInsertionEnd[j]
					}
					result.values[haplotype][readIndex] = log10(sum) - initialConditionLog10
				}
			})
		}
	})

	if len(haplotypes) > 1 {
		for r := range result.alns {
			bestLikelihood := math.Inf(-1)
			for _, haplotype := range haplotypes {
				if !haplotype.isRef {
					if likelihood := result.values[haplotype][r]; likelihood > bestLikelihood {
						bestLikelihood = likelihood
					}
				}
			}
			if !math.IsInf(bestLikelihood, -1) {
				worstLikelihoodCap := bestLikelihood + globalReadMismappingRate
				for _, haplotype := range haplotypes {
					if l := result.values[haplotype]; l[r] < worstLikelihoodCap {
						l[r] = worstLikelihoodCap
					}
				}
			}
		}
	}

checkPoorlyModeledReads:
	for i := 0; i < len(result.alns); {
		maxErrorsForReads := math.Min(2, math.Ceil(float64(len(result.alns[i].QUAL))*0.02))
		log10MaxLikelihoodForTrueAllele := maxErrorsForReads * -4.0
		for _, haplotype := range haplotypes {
			if result.values[haplotype][i] >= log10MaxLikelihoodForTrueAllele {
				i++
				continue checkPoorlyModeledReads
			}
		}
		result.alns = append(result.alns[:i], result.alns[i+1:]...)
		for h, l := range result.values {
			result.values[h] = append(l[:i], l[i+1:]...)
		}
	}

	return result
}
