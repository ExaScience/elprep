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
	"log"
	"math"
	"sort"
	"strings"

	"github.com/exascience/elprep/v5/utils"

	"github.com/exascience/elprep/v5/fasta"

	"github.com/exascience/elprep/v5/sam"
	"github.com/exascience/elprep/v5/vcf"
	"github.com/exascience/pargo/parallel"
)

var singleAlleles = map[byte]string{
	'A': "A",
	'C': "C",
	'G': "G",
	'T': "T",
	'N': "N",
}

func makeBlock(vc1, vc2 *vcf.Variant) {
	if len(vc1.Ref) == 1 && len(vc1.Alt[0]) == 1 {
		if vc1.Ref == vc2.Ref {
			vc1.Alt[0] += vc2.Alt[0][1:]
		} else {
			vc1.Ref = vc2.Ref
			vc1.SetEnd(vc2.End())
		}
	} else {
		var insertion, deletion *vcf.Variant
		if alt := vc1.Alt[0]; len(vc1.Ref) == 1 && len(alt) > 1 && vc1.Ref[0] == alt[0] {
			insertion, deletion = vc1, vc2
		} else {
			insertion, deletion = vc2, vc1
		}
		vc1.Ref = deletion.Ref
		vc1.Alt[0] = insertion.Alt[0]
		vc1.SetEnd(deletion.End())
	}
}

func addEvent(m []*vcf.Variant, vc *vcf.Variant) []*vcf.Variant {
	index := sort.Search(len(m), func(i int) bool {
		return m[i].Pos >= vc.Pos
	})
	if index == len(m) {
		m = append(m, vc)
		if vc.Pos < m[0].Pos {
			copy(m[1:], m)
			m[0] = vc
		}
		return m
	}
	if prev := m[index]; prev.Pos == vc.Pos {
		makeBlock(prev, vc)
		return m
	}
	m = append(m, vc)
	copy(m[index+1:], m[index:])
	m[index] = vc
	return m
}

func makeEventMap(source, contig string, h *haplotype, ref []byte, startPosKeySet map[int32]bool) (result []*vcf.Variant) {
	refPos := h.location
	if refPos < 1 {
		return result
	}

	alignment := h.bases
	var alignmentPos int32

	for cigarIndex, ce := range h.cigar {
		switch ce.Operation {
		case 'I':
			if refPos > 1 {
				refByte := fasta.ToUpperAndN(ref[refPos-2])
				if refAllele, ok := singleAlleles[refByte]; ok && cigarIndex > 0 && cigarIndex < len(h.cigar)-1 {
					insertionBases := alignment[alignmentPos : alignmentPos+ce.Length]
					allSimple := true
					for i := 0; i < len(insertionBases); i++ {
						if _, ok := simpleBaseToBaseIndexTable[insertionBases[i]]; !ok {
							allSimple = false
							break
						}
					}
					if allSimple {
						insertionStart := refPos - 1
						vc := &vcf.Variant{
							Source: source,
							Chrom:  contig,
							Pos:    insertionStart,
							Ref:    refAllele,
							Alt:    []string{string(append([]byte(refAllele), insertionBases...))},
						}
						result = addEvent(result, vc)
						if startPosKeySet != nil {
							startPosKeySet[insertionStart] = true
						}
					}
				}
			}
			alignmentPos += ce.Length
		case 'S':
			alignmentPos += ce.Length
		case 'D':
			if refPos > 1 {
				refByte := fasta.ToUpperAndN(ref[refPos-2])
				if refAllele, ok := singleAlleles[refByte]; ok {
					deletionBases := ref[refPos-2 : refPos-1+ce.Length]
					allSimple := true
					for _, b := range deletionBases {
						if _, ok := simpleBaseToBaseIndexTable[b]; !ok {
							allSimple = false
							break
						}
					}
					if allSimple {
						deletionStart := refPos - 1
						var deletions strings.Builder
						for _, b := range deletionBases {
							deletions.WriteByte(fasta.ToUpperAndN(b))
						}
						vc := &vcf.Variant{
							Source: source,
							Chrom:  contig,
							Pos:    deletionStart,
							Ref:    deletions.String(),
							Alt:    []string{refAllele},
						}
						vc.SetEnd(deletionStart + ce.Length)
						result = addEvent(result, vc)
						if startPosKeySet != nil {
							startPosKeySet[deletionStart] = true
						}
					}
				}
			}
			refPos += ce.Length
		case 'M', '=', 'X':
			for offset := int32(0); offset < ce.Length; offset++ {
				refByte := fasta.ToUpperAndN(ref[refPos-1+offset])
				altByte := alignment[alignmentPos+offset]
				if refByte != altByte {
					if _, ok := simpleBaseToBaseIndexTable[refByte]; ok {
						if _, ok := simpleBaseToBaseIndexTable[altByte]; ok {
							matchStart := refPos + offset
							vc := &vcf.Variant{
								Source: source,
								Chrom:  contig,
								Pos:    matchStart,
								Ref:    singleAlleles[refByte],
								Alt:    []string{singleAlleles[altByte]},
							}
							result = addEvent(result, vc)
							if startPosKeySet != nil {
								startPosKeySet[matchStart] = true
							}
						}
					}
				}
			}
			refPos += ce.Length
			alignmentPos += ce.Length
		default:
			log.Panicf("Unknown cigar operation %v.", ce.Operation)
		}
	}

	return result
}

func getOverlappingEvents(loc int32, haplotypes []*haplotype) map[*haplotype][]*vcf.Variant {
	result := make(map[*haplotype][]*vcf.Variant, len(haplotypes))
	for _, h := range haplotypes {
		var overlaps []*vcf.Variant
		for _, variant := range h.events {
			if variant.Pos <= loc && variant.End() >= loc {
				overlaps = append(overlaps, variant)
			}
		}
		result[h] = overlaps
	}
	return result
}

type locationAndAlleles struct {
	pos  int32
	ref  string
	alts string
}

func makeLocationAndAlleles(v *vcf.Variant) locationAndAlleles {
	var alts strings.Builder
	for _, alt := range v.Alt {
		alts.WriteString(alt)
		alts.WriteByte(',')
	}
	return locationAndAlleles{
		pos:  v.Pos,
		ref:  v.Ref,
		alts: alts.String(),
	}
}

var spanDel = []string{"*"}

func computeActiveVariantContextsWithSpanDelsReplaced(loc int32, haplotypes []*haplotype, overlaps map[*haplotype][]*vcf.Variant, ref []byte) (result []*vcf.Variant) {
	uniqueLocationsAndAlleles := make(map[locationAndAlleles]bool)
	refAllele := singleAlleles[fasta.ToUpperAndN(ref[loc-1])]
	var replacement *vcf.Variant
	for _, h := range haplotypes {
		for _, variant := range overlaps[h] {
			locationAndAlleles := makeLocationAndAlleles(variant)
			if !uniqueLocationsAndAlleles[locationAndAlleles] {
				uniqueLocationsAndAlleles[locationAndAlleles] = true
				if variant.Pos != loc {
					if replacement == nil {
						replacement = &vcf.Variant{
							Chrom: variant.Chrom,
							Pos:   loc,
							Ref:   refAllele,
							Alt:   spanDel,
						}
						replacement.SetEnd(loc)
					}
					variant = replacement
				}
				result = append(result, variant)
			}
		}
	}
	return result
}

func sortBySources(events []*vcf.Variant) {
	sources := make([]string, 0, len(events))

collectSources:
	for _, e := range events {
		for _, s := range sources {
			if s == e.Source {
				continue collectSources
			}
		}
		sources = append(sources, e.Source)
	}

	sort.SliceStable(events, func(i, j int) bool {
		si, sj := events[i].Source, events[j].Source
		for _, s := range sources {
			if sj == s {
				return false
			}
			if si == s {
				return true
			}
		}
		return false
	})
}

func variantSize(vc *vcf.Variant) int32 {
	return vc.End() - vc.Start() + 1
}

func addAllele(alleles []string, allele string) []string {
	for _, a := range alleles {
		if a == allele {
			return alleles
		}
	}
	return append(alleles, allele)
}

func isSymbolicAllele(a string) bool {
	if len(a) <= 1 {
		return false
	}
	if first := a[0]; first == '<' || first == '.' {
		return true
	}
	if last := a[len(a)-1]; last == '>' || last == '.' {
		return true
	}
	for i := 0; i < len(a); i++ {
		if c := a[i]; c == '[' || c == ']' {
			return true
		}
	}
	return false
}

var emptyId = []string{"."}

type alleleMap struct {
	alleles    []string // alleles[0] is ref, all other are alt
	haplotypes map[string][]*haplotype
}

func newAlleleMap(refAllele string) *alleleMap {
	return &alleleMap{
		alleles:    []string{refAllele},
		haplotypes: map[string][]*haplotype{refAllele: nil},
	}
}

func (m *alleleMap) addAllele(allele string) {
	m.alleles = append(m.alleles, allele)
	m.haplotypes[allele] = nil
}

func (m *alleleMap) remove(allele string) {
	for i, a := range m.alleles {
		if a == allele {
			m.alleles = append(m.alleles[:i], m.alleles[i+1:]...)
			break
		}
	}
	delete(m.haplotypes, allele)
}

func (m *alleleMap) maybeAdd(allele string, h *haplotype) {
	if hs, ok := m.haplotypes[allele]; ok {
		m.haplotypes[allele] = append(hs, h)
	}
}

func (m *alleleMap) add(allele string, h *haplotype) {
	if hs, ok := m.haplotypes[allele]; ok {
		m.haplotypes[allele] = append(hs, h)
	} else {
		m.alleles = append(m.alleles, allele)
		m.haplotypes[allele] = []*haplotype{h}
	}
}

const maxAcceptableAlleleCount = 44

type alleleScoredByHaplotype struct {
	allele                     string
	isRef                      bool
	bestScore, secondBestScore float64
}

type readAlleleLikelihoods struct {
	alleles []string
	alns    []*sam.Alignment
	values  map[string][]float64
}

func marginalize(likelihoods readLikelihoods, alleleMapper *alleleMap, start, stop int32) (result readAlleleLikelihoods) {
	result = readAlleleLikelihoods{
		alleles: alleleMapper.alleles,
		alns:    make([]*sam.Alignment, 0, len(likelihoods.alns)),
		values:  make(map[string][]float64),
	}
	readsToKeep := make([]int, 0, len(likelihoods.alns))
	for r, aln := range likelihoods.alns {
		if rstart, rend := aln.POS, aln.End(); (start >= rstart && start <= rend) || (stop >= rstart && stop <= rend) || (rstart >= start && rend <= stop) {
			readsToKeep = append(readsToKeep, r)
			result.alns = append(result.alns, aln)
		}
	}
	for _, allele := range alleleMapper.alleles {
		newValues := make([]float64, len(readsToKeep))
		for i := range newValues {
			newValues[i] = math.Inf(-1)
		}
		for _, h := range alleleMapper.haplotypes[allele] {
			oldValues := likelihoods.values[h]
			for newIndex, oldIndex := range readsToKeep {
				if v := oldValues[oldIndex]; v > newValues[newIndex] {
					newValues[newIndex] = v
				}
			}
		}
		result.values[allele] = newValues
	}
	return
}

func (likelihoods *readAlleleLikelihoods) updateNonRef(refAllele string, altAlleles []string) {
	nonRefAlleleIndex := len(altAlleles) - 1
	qualifiedAlleleLikelihoods := make([]float64, 0, len(altAlleles))
	var nonRefLikelihoods []float64
	if values, ok := likelihoods.values[nonRef]; ok {
		nonRefLikelihoods = values
	} else {
		nonRefLikelihoods = make([]float64, len(likelihoods.alns))
		negInf := math.Inf(-1)
		for i := range nonRefLikelihoods {
			nonRefLikelihoods[i] = negInf
		}
		likelihoods.values[nonRef] = nonRefLikelihoods
	}
	for aln := range likelihoods.alns {
		bestLikelihood := math.Inf(-1)
		for _, allele := range likelihoods.alleles {
			if likelihood := likelihoods.values[allele][aln]; likelihood > bestLikelihood {
				bestLikelihood = likelihood
			}
		}
		qualifiedAlleleLikelihoods = qualifiedAlleleLikelihoods[:0]
		if alleleLikelihood := likelihoods.values[refAllele][aln]; !math.IsNaN(alleleLikelihood) && alleleLikelihood < bestLikelihood {
			qualifiedAlleleLikelihoods = append(qualifiedAlleleLikelihoods, alleleLikelihood)
		}
		for i := 0; i < nonRefAlleleIndex; i++ {
			if alleleLikelihood := likelihoods.values[altAlleles[i]][aln]; !math.IsNaN(alleleLikelihood) && alleleLikelihood < bestLikelihood {
				qualifiedAlleleLikelihoods = append(qualifiedAlleleLikelihoods, alleleLikelihood)
			}
		}
		var median float64
		switch len(qualifiedAlleleLikelihoods) {
		case 0:
			if len(altAlleles) <= 1 {
				nonRefLikelihoods[aln] = math.NaN()
			} else {
				nonRefLikelihoods[aln] = bestLikelihood
			}
			continue
		case 1:
			median = qualifiedAlleleLikelihoods[0]
		default:
			sort.Float64s(qualifiedAlleleLikelihoods)
			half := (len(qualifiedAlleleLikelihoods) + 1) / 2
			lo := qualifiedAlleleLikelihoods[half-1]
			hi := qualifiedAlleleLikelihoods[half]
			// the following is the same as (hi + lo)/2 for even lengths, except for NaN, +Inf, -Inf
			// the following is the same as lo for odd lengths, except for NaN, +Inf, -Inf
			median = lo + (float64(1-len(qualifiedAlleleLikelihoods)%2)/2)*(hi-lo)
		}
		if math.IsNaN(median) && nonRefAlleleIndex > 1 {
			median = bestLikelihood
		}
		nonRefLikelihoods[aln] = median
	}
}

// note: see the VCF specification to understand how the order of genotypes is determined
// (which leads to this strange nesting of two loops)
func forEachAltGenotype(ref string, alt []string, refAndOneAltComponent, oneAltComponent func(index int, allele string), twoAltComponents func(index int, allele0, allele1 string)) {
	index := 1
	for j, a := range alt {
		refAndOneAltComponent(index, a)
		index++
		for i := 0; i < j; i++ {
			twoAltComponents(index, alt[i], a)
			index++
		}
		oneAltComponent(index, a)
		index++
	}
}

func computeSingleComponentGenotypeLikelihood(likelihoods []float64) float64 {
	var gl float64
	for _, l := range likelihoods {
		gl += l + log10Ploidy
	}
	return gl
}

func computeTwoComponentGenotypeLikelihood(likelihoods0, likelihoods1 []float64) float64 {
	var gl float64
	for r := range likelihoods0 {
		gl += approximateLog10SumLog10(likelihoods0[r], likelihoods1[r])
	}
	return gl
}

func findBestAlleles(nofAlleles int, gls []float64) (allele1, allele2, bestIndex int) {
	k := 0
	maxgl := math.Inf(-1)
	for j := 0; j < nofAlleles; j++ {
		for i := 0; i <= j; i, k = i+1, k+1 {
			if gl := gls[k]; gl > maxgl {
				maxgl, allele1, allele2, bestIndex = gl, i, j, k
			}
		}
	}
	return
}

func contains(s string, l []string) bool {
	for _, ls := range l {
		if s == ls {
			return true
		}
	}
	return false
}

func subsetAlleles(vc *vcf.Variant, gls []float64, allelesSubset []string) (plsSubset []interface{}, glsSubset []float64) {
	var glsSub [(8*8 + 8) / 2]float64
	maxGL := gls[0]
	glsSub[0] = maxGL
	rgindex := 1
	refAndOrOneComponent := func(index int, alt string) {
		if contains(alt, allelesSubset) {
			gl := gls[index]
			if gl > maxGL {
				maxGL = gl
			}
			glsSub[rgindex] = gl
			rgindex++
		}
	}
	forEachAltGenotype(vc.Ref, vc.Alt,
		refAndOrOneComponent,
		refAndOrOneComponent,
		func(index int, alt1, alt2 string) {
			if contains(alt1, allelesSubset) && contains(alt2, allelesSubset) {
				gl := gls[index]
				if gl > maxGL {
					maxGL = gl
				}
				glsSub[rgindex] = gl
				rgindex++
			}
		})
	var sum float64
	for i := 0; i < rgindex; i++ {
		e := glsSub[i] - maxGL
		sum += e
		glsSub[i] = e
	}
	if sum < -0.1 {
		// informative
		plsSubset = make([]interface{}, rgindex)
		for i := 0; i < rgindex; i++ {
			if adjusted := -10 * glsSub[i]; adjusted > math.MaxInt32 {
				plsSubset[i] = int(math.MaxInt32)
			} else {
				plsSubset[i] = int(math.Round(adjusted))
			}
		}
		return plsSubset, glsSub[:rgindex]
	}
	return nil, nil
}

type alleleFrequency struct {
	log10Posteriors [2]float64
	log10pRef       map[string]float64
	alleleCounts    []interface{}
}

var nonInformativeAlleleFrequency = alleleFrequency{
	log10pRef: make(map[string]float64),
}

const alleleFrequencyDummyPrior = -1e-10

func init() {
	afLog10Posteriors := [2]float64{alleleFrequencyDummyPrior, math.Inf(-1) + alleleFrequencyDummyPrior}
	log10Sum := log10SumLog10(afLog10Posteriors[0], afLog10Posteriors[1])
	afLog10Posteriors[0] -= log10Sum
	afLog10Posteriors[1] -= log10Sum
	nonInformativeAlleleFrequency.log10Posteriors = afLog10Posteriors
}

var log10TwoComponentCombinationCount = log10Gamma(3) - log10Gamma(2) - log10Gamma(2)

func log10NormalizedGenotypePosteriors(log10Posteriors []float64, vc *vcf.Variant, reducedAlt []string, reducedGls [(8*8 + 8) / 2]float64, log10AlleleFrequencies map[string]float64) {
	log10AlleleFrequenciesForRef := log10AlleleFrequencies[vc.Ref]
	maxValue := reducedGls[0] + 2*log10AlleleFrequenciesForRef
	refWithAlt := log10TwoComponentCombinationCount + log10AlleleFrequenciesForRef
	log10Posteriors[0] = maxValue
	var maxIndex int
	forEachAltGenotype(vc.Ref, reducedAlt,
		func(index int, alt string) {
			value := refWithAlt + reducedGls[index] + log10AlleleFrequencies[alt]
			if value > maxValue {
				maxValue = value
				maxIndex = index
			}
			log10Posteriors[index] = value
		},
		func(index int, alt string) {
			value := reducedGls[index] + 2*log10AlleleFrequencies[alt]
			if value > maxValue {
				maxValue = value
				maxIndex = index
			}
			log10Posteriors[index] = value
		}, func(index int, alt1, alt2 string) {
			value := log10TwoComponentCombinationCount + reducedGls[index] +
				log10AlleleFrequencies[alt1] + log10AlleleFrequencies[alt2]
			if value > maxValue {
				maxValue = value
				maxIndex = index
			}
			log10Posteriors[index] = value
		})
	var log10Sum float64
	if math.IsInf(maxValue, -1) {
		log10Sum = maxValue
	} else {
		sum := 1.0
		for i, value := range log10Posteriors {
			if i == maxIndex || math.IsInf(value, -1) {
				continue
			}
			sum += math.Pow(10, value-maxValue)
		}
		log10Sum = maxValue + log10(sum)
	}
	for i := range log10Posteriors {
		log10Posteriors[i] -= log10Sum
	}
}

func log10SumLog10(a, b float64) float64 {
	if a > b {
		return a + log10(1+math.Pow(10, b-a))
	}
	return b + log10(1+math.Pow(10, a-b))
}

func log10SumLog10Slice(values []float64) float64 {
	if len(values) == 0 {
		return math.Inf(-1)
	}
	maxValue := values[0]
	var maxIndex int
	for i := 1; i < len(values); i++ {
		if v := values[i]; v > maxValue {
			maxValue = v
			maxIndex = i
		}
	}
	if math.IsInf(maxValue, -1) {
		return maxValue
	}
	sum := 1.0
	for i, v := range values {
		if i == maxIndex || math.IsInf(v, -1) {
			continue
		}
		sum += math.Pow(10, v-maxValue)
	}
	return maxValue + log10(sum)
}

func (hc *HaplotypeCaller) computeLog10NonRef(vc *vcf.Variant, reducedAlt []string, reducedPls []interface{}) alleleFrequency {
	if len(reducedPls) == 0 {
		// non informative case
		af := nonInformativeAlleleFrequency
		af.alleleCounts = make([]interface{}, len(reducedAlt))
		for i := range reducedAlt {
			af.alleleCounts[i] = int(0)
		}
		return af
	}

	numAlleles := len(reducedAlt) + 1
	priorPseudocounts := make([]float64, 0, numAlleles)
	priorPseudocounts = append(priorPseudocounts, hc.refPseudocount)
	for _, a := range reducedAlt {
		if len(a) <= 1 || isSymbolicAllele(a) {
			priorPseudocounts = append(priorPseudocounts, hc.indelPseudocount)
		} else {
			priorPseudocounts = append(priorPseudocounts, hc.snpPseudocount)
		}
	}

	var reducedGls [(8*8 + 8) / 2]float64

	for i := range reducedPls {
		reducedGls[i] = float64(reducedPls[i].(int)) / -10
	}

	alleleCounts := make(map[string]float64)
	newAlleleCounts := make(map[string]float64)
	log10AlleleFrequencies := make(map[string]float64)
	flatLog10AlleleFrequency := -log10(float64(numAlleles))
	vcRef := vc.Ref
	log10AlleleFrequencies[vcRef] = flatLog10AlleleFrequency
	for _, a := range reducedAlt {
		log10AlleleFrequencies[a] = flatLog10AlleleFrequency
	}
	log10Posteriors := make([]float64, len(reducedPls))
	posteriorPseudocounts := make([]float64, numAlleles)
	for {
		newAlleleCounts[vcRef] = math.Inf(-1)
		for _, a := range reducedAlt {
			newAlleleCounts[a] = math.Inf(-1)
		}
		log10NormalizedGenotypePosteriors(log10Posteriors, vc, reducedAlt, reducedGls, log10AlleleFrequencies)
		newAlleleCounts[vcRef] = log10SumLog10(newAlleleCounts[vcRef], log10Posteriors[0]+log10Ploidy)
		forEachAltGenotype(vcRef, reducedAlt,
			func(index int, alt string) {
				inc := log10Posteriors[index] + log10One
				newAlleleCounts[vcRef] = log10SumLog10(newAlleleCounts[vcRef], inc)
				newAlleleCounts[alt] = log10SumLog10(newAlleleCounts[alt], inc)
			}, func(index int, alt string) {
				newAlleleCounts[alt] = log10SumLog10(newAlleleCounts[alt], log10Posteriors[index]+log10Ploidy)
			}, func(index int, alt1, alt2 string) {
				inc := log10Posteriors[index] + log10One
				newAlleleCounts[alt1] = log10SumLog10(newAlleleCounts[alt1], inc)
				newAlleleCounts[alt2] = log10SumLog10(newAlleleCounts[alt2], inc)
			})
		for a, x := range newAlleleCounts {
			newAlleleCounts[a] = math.Pow(10, x)
		}

		nac := newAlleleCounts[vcRef]
		alleleCountsMaximumDifference := math.Abs(alleleCounts[vcRef] - nac)
		sum := priorPseudocounts[0] + nac
		posteriorPseudocounts[0] = sum
		for i, a := range reducedAlt {
			nac := newAlleleCounts[a]
			if value := math.Abs(alleleCounts[a] - nac); value > alleleCountsMaximumDifference {
				alleleCountsMaximumDifference = value
			}
			value := priorPseudocounts[i+1] + nac
			sum += value
			posteriorPseudocounts[i+1] = value
		}
		alleleCounts, newAlleleCounts = newAlleleCounts, alleleCounts

		log10AlleleFrequencies[vcRef] = log10(posteriorPseudocounts[0] / sum)
		for i, a := range reducedAlt {
			log10AlleleFrequencies[a] = log10(posteriorPseudocounts[i+1] / sum)
		}

		if alleleCountsMaximumDifference <= 0.1 {
			break
		}
	}

	var log10PNoVariant float64

	log10NormalizedGenotypePosteriors(log10Posteriors, vc, reducedAlt, reducedGls, log10AlleleFrequencies)
	var nonVariantPosteriors [2]float64
	nonVariantPosteriors[0] = log10Posteriors[0]
	nvpLength := 1
	forEachAltGenotype(vcRef, reducedAlt, func(index int, alt string) {
		if alt == "*" {
			nonVariantPosteriors[nvpLength] = log10Posteriors[index]
			nvpLength++
		}
	}, func(_ int, _ string) {}, func(_ int, _, _ string) {})
	if nvpLength == 1 {
		log10PNoVariant = log10Posteriors[0]
	} else {
		log10PNoVariant = math.Min(0, log10SumLog10(nonVariantPosteriors[0], nonVariantPosteriors[1]))
	}

	log10POfZeroCountsByAllele := make(map[string]float64)

	if numAlleles == 2 {
		log10POfZeroCountsByAllele[reducedAlt[0]] = log10PNoVariant
	} else {
		log10AbsentPosteriors := make(map[string][]float64)
		posterior := log10Posteriors[0]
		for _, a := range reducedAlt {
			log10AbsentPosteriors[a] = append(log10AbsentPosteriors[a], posterior)
		}
		refAndOrOneComponent := func(index int, alt string) {
			posterior := log10Posteriors[index]
			for _, a := range reducedAlt {
				if a != alt {
					log10AbsentPosteriors[a] = append(log10AbsentPosteriors[a], posterior)
				}
			}
		}
		forEachAltGenotype(vcRef, reducedAlt,
			refAndOrOneComponent,
			refAndOrOneComponent,
			func(index int, alt1, alt2 string) {
				posterior := log10Posteriors[index]
				for _, a := range reducedAlt {
					if a != alt1 && a != alt2 {
						log10AbsentPosteriors[a] = append(log10AbsentPosteriors[a], posterior)
					}
				}
			})
		for _, a := range reducedAlt {
			log10POfZeroCountsByAllele[a] = math.Min(0, log10SumLog10Slice(log10AbsentPosteriors[a]))
		}
	}

	log10PosteriorYesNo := [2]float64{log10PNoVariant, log10OneMinusPow10(log10PNoVariant)}
	afLog10Posteriors := [2]float64{
		log10PosteriorYesNo[0] + alleleFrequencyDummyPrior,
		log10PosteriorYesNo[1] + alleleFrequencyDummyPrior,
	}
	log10Sum := log10SumLog10(afLog10Posteriors[0], afLog10Posteriors[1])
	afLog10Posteriors[0] -= log10Sum
	afLog10Posteriors[1] -= log10Sum

	intAlleleCounts := make([]interface{}, len(reducedAlt))
	for i, a := range reducedAlt {
		intAlleleCounts[i] = int(int32(math.Round(alleleCounts[a])))
	}

	return alleleFrequency{
		log10Posteriors: afLog10Posteriors,
		log10pRef:       log10POfZeroCountsByAllele,
		alleleCounts:    intAlleleCounts,
	}
}

func isVcCoveredByDeletion(deletions *deletionsHandler, vc *vcf.Variant) bool {
	deletions.wg.Wait()
	for i := 0; i < len(deletions.slice); {
		if d := deletions.slice[i]; d.end < vc.Pos {
			deletions.slice = append(deletions.slice[:i], deletions.slice[i+1:]...)
		} else if d.start != vc.Pos {
			return true
		} else {
			i++
		}
	}
	return false
}

func (hc *HaplotypeCaller) computeOutputAlleles(merged *vcf.Variant, reducedAlt []string, af alleleFrequency, deletions *deletionsHandler) (outputAlleles []string, mleCounts []interface{}, siteIsMonomorphic bool) {
	siteIsMonomorphic = true
	referenceAlleleSize := int32(len(merged.Ref))
	if len(reducedAlt) == 1 && reducedAlt[0] == nonRef {
		isPlausible := af.log10pRef[nonRef]+1.0e-10 < hc.standardConfidenceForCallingByMin10
		siteIsMonomorphic = siteIsMonomorphic && !isPlausible
		outputAlleles = reducedAlt
		mleCounts = af.alleleCounts
		if referenceAlleleSize > 0 {
			deletions.wg.Wait()
			deletions.slice = append(deletions.slice, deletion{merged.Pos, merged.Pos + referenceAlleleSize})
		}
	} else {
		outputAlleles = make([]string, 0, len(reducedAlt))
		mleCounts = make([]interface{}, 0, len(reducedAlt))
		for i, a := range reducedAlt {
			isPlausible := af.log10pRef[a]+1.0e-10 < hc.standardConfidenceForCallingByMin10
			isSpuriousSpanningDeletion := a == "*" && !isVcCoveredByDeletion(deletions, merged)
			toOutput := (hc.confidenceMode != none || isPlausible || a == nonRef) && !isSpuriousSpanningDeletion
			siteIsMonomorphic = siteIsMonomorphic && !(isPlausible && !isSpuriousSpanningDeletion)
			if toOutput {
				outputAlleles = append(outputAlleles, a)
				mleCounts = append(mleCounts, af.alleleCounts[i])
				deletionSize := referenceAlleleSize
				if !isSymbolicAllele(a) {
					deletionSize -= int32(len(a))
				}
				if deletionSize > 0 {
					deletions.wg.Wait()
					deletions.slice = append(deletions.slice, deletion{merged.Pos, merged.Pos + deletionSize})
				}
			}
		}
	}
	return outputAlleles, mleCounts, siteIsMonomorphic
}

func (hc *HaplotypeCaller) calculateGenotypes(variant *vcf.Variant, pls []interface{}, gls []float64, deletions *deletionsHandler) (*vcf.Variant, []float64) {
	if len(variant.Alt) > 49 {
		return nil, nil
	}

	reducedAlt := variant.Alt
	reducedPls := pls
	var nofNonProperAlts int
	nofReducedAlt := len(reducedAlt)
	hasNonRef := hc.confidenceMode != none // <=> contains(nonRef, variant.Alt)
	if hasNonRef {
		nofNonProperAlts = 1
		nofReducedAlt--
	}
	if nofReducedAlt > 6 {
		a1, a2, _ := findBestAlleles(len(reducedAlt)+1, gls)
		var bestAlt, oneIndex int
		if a1 > 0 && a1 < nofReducedAlt+1 {
			bestAlt = 1
			oneIndex = a1
		}
		if a2 > 0 && a2 < nofReducedAlt+1 && a1 != a2 {
			bestAlt++
			oneIndex = a2
		}
		reducedAlt = make([]string, 6+nofNonProperAlts)
		copy(reducedAlt[:6], variant.Alt)
		switch bestAlt {
		case 1:
			if oneIndex--; oneIndex > 5 {
				reducedAlt[5] = variant.Alt[oneIndex]
			}
		case 2:
			if a1, a2 = a1-1, a2-1; a1 > 4 {
				reducedAlt[4] = variant.Alt[a1]
				reducedAlt[5] = variant.Alt[a2]
			} else if a2 > 5 {
				reducedAlt[5] = variant.Alt[a2]
			}
		}
		if hasNonRef {
			reducedAlt[6] = nonRef
		}
		reducedPls, _ = subsetAlleles(variant, gls, reducedAlt)
	}
	af := hc.computeLog10NonRef(variant, reducedAlt, reducedPls)
	outputAlleles, mleCounts, siteIsMonomorphic := hc.computeOutputAlleles(variant, reducedAlt, af, deletions)
	if len(outputAlleles) == 1 && outputAlleles[0] == "*" {
		return nil, nil
	}
	if siteIsMonomorphic && (len(outputAlleles) == 0 || outputAlleles[0] != nonRef) {
		return nil, nil
	}
	var log10Confidence float64
	if hc.confidenceMode != none || !siteIsMonomorphic {
		log10Confidence = af.log10Posteriors[0]
	} else {
		log10Confidence = af.log10Posteriors[1]
	}
	if log10Confidence == 0 && math.Signbit(log10Confidence) {
		log10Confidence = 0
	}
	phredScaledConfidence := -10 * log10Confidence
	if phredScaledConfidence == 0 && math.Signbit(phredScaledConfidence) {
		phredScaledConfidence = 0
	}
	var filter []utils.Symbol
	if !(phredScaledConfidence >= hc.standardConfidenceForCalling) {
		if len(outputAlleles) == 0 || outputAlleles[0] != nonRef {
			return nil, nil
		}
		filter = []utils.Symbol{LowQual}
	}
	call := &vcf.Variant{
		Source: "HC_call",
		Chrom:  variant.Chrom,
		Pos:    variant.Pos,
		Ref:    variant.Ref,
		Alt:    outputAlleles,
		Qual:   phredScaledConfidence,
		Filter: filter,
		Info:   variant.Info,
	}
	var gt vcf.Genotype
	if len(outputAlleles) == 0 {
		pls = nil
		gls = nil
		call.Alt = nil
		gt.GT = noVariationGT
	} else {
		pls, gls = subsetAlleles(variant, gls, outputAlleles)
		if pls == nil {
			gt.GT = noCallGT
		} else {
			a1, a2, bestgl := findBestAlleles(len(outputAlleles)+1, gls)
			gt.GT = []int32{int32(a1), int32(a2)}
			gt.Data.Set(PL, pls)
			qual := math.Inf(-1)
			for i, g := range gls {
				if i != bestgl && g >= qual {
					qual = g
				}
			}
			qual = gls[bestgl] - qual
			var log10PError float64
			if qual < 0 {
				maxValue := gls[bestgl]
				var sum float64
				for i := range gls {
					v := math.Pow(10, gls[i]-maxValue)
					gls[i] = v
					sum += v
				}
				log10PError = log10(1 - gls[bestgl]/sum)
			} else {
				log10PError = -qual
			}
			gt.Data.Set(GQ, minInt(int(math.Round(log10PError*-10)), maxGenotypeQual))
		}
	}
	if len(mleCounts) > 0 {
		call.Info.Set(MLEAC, mleCounts)
		mleFrequencies := make([]interface{}, len(mleCounts))
		var idiv int
		for _, a := range gt.GT {
			if a != -1 {
				idiv++
			}
		}
		if idiv == 0 {
			for i := range mleCounts {
				mleFrequencies[i] = math.NaN()
			}
		} else {
			fdiv := float64(idiv)
			for i, c := range mleCounts {
				mleFrequencies[i] = math.Min(1, float64(c.(int))/fdiv)
			}
		}
		call.Info.Set(MLEAF, mleFrequencies)
	}
	call.GenotypeData = []vcf.Genotype{gt}

	return call, gls
}

type (
	rank struct {
		value float64
		rank  float32
		isAlt bool
	}

	rankSumTest struct {
		alts, refs []rank
	}
)

func calcAlignmentByteArrayOffset(cigar []sam.CigarOperation, offset int) int {
	var pos, alignmentPos int
	for _, ce := range cigar {
		l := int(ce.Length)
		switch ce.Operation {
		case 'I', 'S':
			pos += l
			if pos >= offset {
				return alignmentPos
			}
		case 'D':
			alignmentPos += l
		case 'M', '=', 'X':
			if pos+l-1 >= offset {
				return alignmentPos + offset - pos
			} else {
				pos += l
				alignmentPos += l
			}
		}
	}
	return alignmentPos
}

func computeDiploidGenotypeCounts(vc *vcf.Variant, gls []float64) (refCount, hetCount, homCount int) {
	if _, ok := vc.GenotypeData[0].Data.Get(PL); !ok {
		return
	}

	var (
		idxAA = 0
		idxAB = 1
		idxBB = 2
	)

	gt := vc.GenotypeData[0].GT

	if len(vc.Alt) != 1 {
		if gt[0] != gt[1] && gt[0] != 0 && gt[1] != 0 { // HET
			return 0, 0, 1
		}

		if gt[1] != 0 {
			a2 := int(gt[1])
			idxAB = (a2 * (a2 + 1) / 2) + 0
			idxBB = idxAB + a2
		} else if gt[0] != 0 {
			a2 := int(gt[0])
			idxAB = (a2 * (a2 + 1) / 2) + 0
			idxBB = idxAB + a2
		}
	}

	log10Sum := log10SumLog10Slice(gls)

	refCount = int(math.Round(math.Pow(10, gls[idxAA]-log10Sum)))
	hetCount = int(math.Round(math.Pow(10, gls[idxAB]-log10Sum)))
	homCount = int(math.Round(math.Pow(10, gls[idxBB]-log10Sum)))
	return
}

const minNeededValue = 1.0e-16

func exactTest(hetCount, refCount, homCount int) float64 {
	var obsHomR, obsHomC int
	if refCount < homCount {
		obsHomR, obsHomC = refCount, homCount
	} else {
		obsHomR, obsHomC = homCount, refCount
	}
	rareCopies := 2*obsHomR + hetCount

	if rareCopies <= 1 {
		return 0.5
	}

	N := hetCount + obsHomC + obsHomR

	probs := make([]float64, rareCopies+1)
	mid := rareCopies * (2*N - rareCopies) / (2*N - 1)
	if mid%2 != rareCopies%2 {
		mid++
	}
	probs[mid] = 1
	sum := 1.0
	curHets := mid
	curHomR := (rareCopies - mid) / 2
	curHomC := N - curHets - curHomR
	for curHets >= 2 {
		potentialProb := probs[curHets] * float64(curHets*(curHets-1)) / float64(4*(curHomR+1)*(curHomC+1))
		if potentialProb < minNeededValue {
			break
		}
		probs[curHets-2] = potentialProb
		sum += potentialProb
		curHets -= 2
		curHomR++
		curHomC++
	}

	curHets = mid
	curHomR = (rareCopies - mid) / 2
	curHomC = N - curHets - curHomR
	for curHets <= rareCopies-2 {
		potentialProb := probs[curHets] * 4 * float64(curHomR*curHomC) / float64((curHets+2)*(curHets+1))
		if potentialProb < minNeededValue {
			break
		}
		probs[curHets+2] = potentialProb
		sum += potentialProb
		curHets += 2
		curHomR--
		curHomC--
	}

	rightPval := probs[hetCount] / (2 * sum)
	if hetCount == rareCopies {
		return rightPval
	}
	var probSum float64
	for i := hetCount + 1; i < len(probs); i++ {
		probSum += probs[i]
	}
	return rightPval + probSum/sum
}

var phredScaledMinPValue = -10 * log10(minNeededValue)

func calculateEH(vc *vcf.Variant, gls []float64) (phredPval float64) {
	refCount, hetCount, homCount := computeDiploidGenotypeCounts(vc, gls)
	pval := exactTest(hetCount, refCount, homCount)
	if pval < 10e-60 {
		return phredScaledMinPValue
	}
	return -10 * log10(pval)
}

var sqrt2 = math.Sqrt(2)

func erfInv(x float64) float64 {
	w := -math.Log((1 - x) * (1 + x))
	var p float64
	if w < 6.25 {
		w -= 3.125
		p = -3.6444120640178196996e-21
		p = -1.685059138182016589e-19 + p*w
		p = 1.2858480715256400167e-18 + p*w
		p = 1.115787767802518096e-17 + p*w
		p = -1.333171662854620906e-16 + p*w
		p = 2.0972767875968561637e-17 + p*w
		p = 6.6376381343583238325e-15 + p*w
		p = -4.0545662729752068639e-14 + p*w
		p = -8.1519341976054721522e-14 + p*w
		p = 2.6335093153082322977e-12 + p*w
		p = -1.2975133253453532498e-11 + p*w
		p = -5.4154120542946279317e-11 + p*w
		p = 1.051212273321532285e-09 + p*w
		p = -4.1126339803469836976e-09 + p*w
		p = -2.9070369957882005086e-08 + p*w
		p = 4.2347877827932403518e-07 + p*w
		p = -1.3654692000834678645e-06 + p*w
		p = -1.3882523362786468719e-05 + p*w
		p = 0.0001867342080340571352 + p*w
		p = -0.00074070253416626697512 + p*w
		p = -0.0060336708714301490533 + p*w
		p = 0.24015818242558961693 + p*w
		p = 1.6536545626831027356 + p*w
	} else if w < 16.0 {
		w = math.Sqrt(w) - 3.25
		p = 2.2137376921775787049e-09
		p = 9.0756561938885390979e-08 + p*w
		p = -2.7517406297064545428e-07 + p*w
		p = 1.8239629214389227755e-08 + p*w
		p = 1.5027403968909827627e-06 + p*w
		p = -4.013867526981545969e-06 + p*w
		p = 2.9234449089955446044e-06 + p*w
		p = 1.2475304481671778723e-05 + p*w
		p = -4.7318229009055733981e-05 + p*w
		p = 6.8284851459573175448e-05 + p*w
		p = 2.4031110387097893999e-05 + p*w
		p = -0.0003550375203628474796 + p*w
		p = 0.00095328937973738049703 + p*w
		p = -0.0016882755560235047313 + p*w
		p = 0.0024914420961078508066 + p*w
		p = -0.0037512085075692412107 + p*w
		p = 0.005370914553590063617 + p*w
		p = 1.0052589676941592334 + p*w
		p = 3.0838856104922207635 + p*w
	} else if !(math.IsInf(w, 1) || math.IsInf(w, -1)) {
		w = math.Sqrt(w) - 5.0
		p = -2.7109920616438573243e-11
		p = -2.5556418169965252055e-10 + p*w
		p = 1.5076572693500548083e-09 + p*w
		p = -3.7894654401267369937e-09 + p*w
		p = 7.6157012080783393804e-09 + p*w
		p = -1.4960026627149240478e-08 + p*w
		p = 2.9147953450901080826e-08 + p*w
		p = -6.7711997758452339498e-08 + p*w
		p = 2.2900482228026654717e-07 + p*w
		p = -9.9298272942317002539e-07 + p*w
		p = 4.5260625972231537039e-06 + p*w
		p = -1.9681778105531670567e-05 + p*w
		p = 7.5995277030017761139e-05 + p*w
		p = -0.00021503011930044477347 + p*w
		p = -0.00013871931833623122026 + p*w
		p = 1.0103004648645343977 + p*w
		p = 4.8499064014085844221 + p*w
	} else {
		p = math.Inf(1)
	}

	return p * x
}

func (test rankSumTest) mannWithneyU() (float64, bool) {
	n1 := len(test.alts)
	n2 := len(test.refs)

	if n1 == 0 || n2 == 0 {
		return 0, false
	}

	ranks := append(test.alts, test.refs...)

	sort.SliceStable(ranks, func(i, j int) bool {
		return ranks[i].value < ranks[j].value
	})

	for i := range ranks {
		ranks[i].rank = float32(i + 1)
	}

	var nties float64

	for i := 0; i < len(ranks); {
		rank := ranks[i].rank
		count := 1
		for j := i + 1; j < len(ranks) && ranks[j].value == ranks[i].value; j++ {
			rank += ranks[j].rank
			count++
		}
		if count > 1 {
			rank /= float32(count)
			for j, k := i, i+count; j < k; j++ {
				ranks[j].rank = rank
			}
			if count != len(ranks) {
				fcount := float64(count)
				nties += math.Pow(fcount, 3) - fcount
			}
		}
		i += count
	}

	var r float32
	for _, rank := range ranks {
		if rank.isAlt {
			r += rank.rank
		}
	}

	u := float64(r - float32((n1*(n1+1))/2))

	var z float64

	if n1 >= 10 || n2 >= 10 {
		m := float64(n1*n2) / 2
		correction := -0.5
		if nties == 0 {
			correction = 0
		}
		sigma := math.Sqrt((float64(n1*n2) / 12) * (float64(n1+n2+1) - nties/float64((n1+n2)*(n1+n2-1))))
		z = (u - m - correction) / sigma
	} else {
		newUDelta := float64((n1 * (n1 + 1)) / 2)
		histogram := make(map[int]int)
		var totalSum int

		permutation := make([]int, n1+n2)
		for i := n1; i < len(permutation); i++ {
			permutation[i] = 1
		}
		for {
			var newU float64
			for i, grouping := range permutation {
				if grouping == 0 {
					newU += float64(ranks[i].rank)
				}
			}
			newU -= newUDelta
			histogram[int(math.Round(2*newU))]++
			totalSum++

			// compute next permutation
			k := -1
			for i := len(permutation) - 2; i >= 0; i-- {
				if permutation[i] < permutation[i+1] {
					k = i
					break
				}
			}
			if k == -1 {
				break
			}
			l := -1
			for i, k1 := len(permutation)-1, k+1; i >= k1; i-- {
				if permutation[k] < permutation[i] {
					l = i
					break
				}
			}
			permutation[k], permutation[l] = permutation[l], permutation[k]
			for begin, end := k+1, len(permutation)-1; begin < end; begin, end = begin+1, end-1 {
				permutation[begin], permutation[end] = permutation[end], permutation[begin]
			}
		}
		u2 := int(math.Round(2 * u))
		sumOfAllSmallerBins := float64(histogram[u2]) / 2
		histogramKeys := make([]int, 0, len(histogram))
		for key := range histogram {
			if key < u2 {
				histogramKeys = append(histogramKeys, key)
			}
		}
		sort.Ints(histogramKeys)
		for _, key := range histogramKeys {
			sumOfAllSmallerBins += float64(histogram[key])
		}
		p := sumOfAllSmallerBins / float64(totalSum)
		z = sqrt2 * erfInv(2*p-1)
	}

	return z, !math.IsNaN(z)
}

func getDeviancePart(x, mu float64) (ret float64) {
	if d, t := x-mu, x+mu; math.Abs(d) < 0.1*t {
		v := d / t
		s1 := v * d
		s := math.NaN()
		ej := 2 * x * v
		v *= v
		for j := 1; s1 != s; j++ {
			s = s1
			ej *= v
			s1 += ej / float64(j*2+1)
		}
		ret = s1
	} else {
		ret = x*math.Log(x/mu) + mu - x
	}
	return
}

const relErr = 1 - 10e-7

var (
	halfLog2Pi = 0.5 * math.Log(6.283185307179586)

	exactStirlingErrors = [...]float64{
		0, 0.15342640972002736, 0.08106146679532726, 0.05481412105191765, 0.0413406959554093,
		0.03316287351993629, 0.02767792568499834, 0.023746163656297496, 0.020790672103765093,
		0.018488450532673187, 0.016644691189821193, 0.015134973221917378, 0.013876128823070748,
		0.012810465242920227, 0.01189670994589177, 0.011104559758206917, 0.010411265261972096,
		0.009799416126158804, 0.009255462182712733, 0.008768700134139386, 0.00833056343336287,
		0.00793411456431402, 0.007573675487951841, 0.007244554301320383, 0.00694284010720953,
		0.006665247032707682, 0.006408994188004207, 0.006171712263039458, 0.0059513701127588475,
		0.0057462165130101155, 0.005554733551962801,
	}
)

func getStirlingError(z float64) (ret float64) {
	if z < 15 {
		z2 := 2 * z
		if math.Floor(z2) == z2 {
			ret = exactStirlingErrors[int(z2)]
		} else {
			lg, _ := math.Lgamma(z + 1)
			ret = lg - (z+0.5)*math.Log(z) + z - halfLog2Pi
		}
	} else {
		z2 := z * z
		ret = (0.08333333333333333 - (0.002777777777777778-(7.936507936507937e-4-(5.952380952380953e-4-8.417508417508417e-4/z2)/z2)/z2)/z2) / z
	}
	return ret
}

func logBinomialProbability(x, n int, p, q float64) (ret float64) {
	fn := float64(n)
	if x == 0 {
		if p < 0.1 {
			ret = -getDeviancePart(fn, fn*q) - fn*p
		} else {
			ret = fn * math.Log(q)
		}
	} else if x == n {
		if q < 0.1 {
			ret = -getDeviancePart(fn, fn*p) - fn*q
		} else {
			ret = fn * math.Log(p)
		}
	} else {
		fx := float64(x)
		fnx := float64(n - x)
		ret = getStirlingError(fn) - getStirlingError(fx) - getStirlingError(fnx) -
			getDeviancePart(fx, fn*p) - getDeviancePart(fnx, fn*q)
		f := 6.283185307179586 * fx * fnx / fn
		ret += -0.5 * math.Log(f)
	}
	return
}

type hypergeometricDistribution struct {
	populationSize, numberOfSuccesses, sampleSize int
	lowerDomain, upperDomain                      int
	p, q, p3                                      float64
}

func makeHypergeometricDistribution(populationSize, numberOfSuccesses, sampleSize int) hypergeometricDistribution {
	fPopSize := float64(populationSize)
	p := float64(sampleSize) / fPopSize
	q := float64(populationSize-sampleSize) / fPopSize
	return hypergeometricDistribution{
		populationSize:    populationSize,
		numberOfSuccesses: numberOfSuccesses,
		sampleSize:        sampleSize,
		lowerDomain:       maxInt(0, numberOfSuccesses-(populationSize-sampleSize)),
		upperDomain:       minInt(sampleSize, numberOfSuccesses),
		p:                 p,
		q:                 q,
		p3:                logBinomialProbability(sampleSize, populationSize, p, q),
	}
}

func (dist hypergeometricDistribution) logProbability(x int) (ret float64) {
	if x >= dist.lowerDomain && x <= dist.upperDomain {
		p1 := logBinomialProbability(x, dist.numberOfSuccesses, dist.p, dist.q)
		p2 := logBinomialProbability(dist.sampleSize-x, dist.populationSize-dist.numberOfSuccesses, dist.p, dist.q)
		ret = p1 + p2 - dist.p3
	} else {
		ret = math.Inf(-1)
	}
	return
}

var minLog10ScaledQual = math.Log10(math.SmallestNonzeroFloat64)

func (hc *HaplotypeCaller) annotateCall(call *vcf.Variant, likelihoods readAlleleLikelihoods, gls []float64) {
	// Genotype annotations: DepthPerAlleleBySample, DepthPerSampleHC, StrandBiasBySample
	// InfoFieldAnnotations: Coverage, RMSMappingQuality, ExcessHet,
	// (BaseQualityRankSumTest, MappingQualityRankSumTest, ReadPosRankSumTest) extends RankSumTest
	// ChromosomeCounts, FisherStrand, StrandOddsRatio, QualityByDepth

	callGT := call.GenotypeData[0].GT

	if hc.confidenceMode == none {
		var an int // AN / ChromosomeCounts
		for _, g := range callGT {
			if g >= 0 {
				an++
			}
		}
		if an > 0 {
			var ac, af []interface{} // AC, AF / ChromosomeCounts
			for i := 1; i <= len(call.Alt); i++ {
				var iac int
				for _, g := range callGT {
					if int(g) == i {
						iac++
					}
				}
				ac = append(ac, iac)
				if an == 0 {
					af = append(af, 0.0)
				} else {
					af = append(af, float64(iac)/float64(an))
				}
			}
			call.Info.Set(AN, an)
			call.Info.Set(AC, ac)
			call.Info.Set(AF, af)
		}
	}

	alleleCounts := make(map[string]int) // AD / DepthPerAlleleBySample
	var depth int                        // DP / DepthPerSampleHC
	var contingencyTable [4]int          // SB / StrandBiasBySample
	var squareSum, numReadsUsed int      // RAW_MQandDP / RMSMappingQuality

	var (
		baseQuality    rankSumTest // BaseQRankSum
		mappingQuality rankSumTest // MQRankSum
		readPosition   rankSumTest // ReadPosRankSum
	)

	for r, aln := range likelihoods.alns {
		if mq := aln.MAPQ; mq != 255 {
			squareSum += int(mq) * int(mq)
			numReadsUsed++
		}

		bestAllele, bestLikelihood := call.Ref, likelihoods.values[call.Ref][r]
		secondBestAllele, secondBestLikelihood := "", math.Inf(-1)
		for _, a := range call.Alt {
			if likelihood := likelihoods.values[a][r]; likelihood > bestLikelihood {
				secondBestAllele, secondBestLikelihood = bestAllele, bestLikelihood
				bestAllele, bestLikelihood = a, likelihood
			} else if likelihood > secondBestLikelihood {
				secondBestAllele, secondBestLikelihood = a, likelihood
			}
		}
		if bestLikelihood-secondBestLikelihood < log10InformativeThreshold {
			if bestAllele != call.Ref {
				if likelihood := likelihoods.values[call.Ref][r]; bestLikelihood-likelihood <= log10InformativeThreshold {
					secondBestAllele, secondBestLikelihood = bestAllele, bestLikelihood
					bestAllele, bestLikelihood = call.Ref, likelihood
				}
			}
		}
		if secondBestAllele != "" && bestLikelihood-secondBestLikelihood > log10InformativeThreshold { // informative
			depth++
			alleleCounts[bestAllele]++
		}

		bestAllele, bestLikelihood = "", math.Inf(-1)
		secondBestAllele, secondBestLikelihood = "", math.Inf(-1)
		for _, a := range likelihoods.alleles {
			if likelihood := likelihoods.values[a][r]; likelihood > bestLikelihood {
				secondBestAllele, secondBestLikelihood = bestAllele, bestLikelihood
				bestAllele, bestLikelihood = a, likelihood
			} else if likelihood > secondBestLikelihood {
				secondBestAllele, secondBestLikelihood = a, likelihood
			}
		}
		if bestLikelihood-secondBestLikelihood < log10InformativeThreshold {
			if bestAllele != call.Ref {
				if likelihood := likelihoods.values[call.Ref][r]; bestLikelihood-likelihood <= log10InformativeThreshold {
					secondBestAllele, secondBestLikelihood = bestAllele, bestLikelihood
					bestAllele, bestLikelihood = call.Ref, likelihood
				}
			}
		}
		if secondBestAllele != "" && bestLikelihood-secondBestLikelihood > log10InformativeThreshold { // informative
			if bestAllele == call.Ref {
				contingencyTable[(aln.FLAG&sam.Reversed)>>4]++
			} else if contains(bestAllele, call.Alt) {
				contingencyTable[2+((aln.FLAG&sam.Reversed)>>4)]++
			}
			if aln.MAPQ != 0 && aln.MAPQ != 255 {
				// BaseQualityRankSumTest, MappingQualityRankSumTest
				var isRef, isAlt bool
				if isRef = bestAllele == call.Ref; !isRef {
					isAlt = contains(bestAllele, call.Alt)
				}
				if isRef || isAlt {
					softStart := softStart(aln)
					leftmostSafePos := maxInt(softStart, int(call.Pos))
					readCoord, ok := getReadCoordinateForReferenceCoordinate(aln.CIGAR, softStart, leftmostSafePos, right)
					if ok {
						baseQual := float64(aln.QUAL[readCoord])
						mappingQual := float64(aln.MAPQ)
						if isRef {
							baseQuality.refs = append(baseQuality.refs, rank{value: baseQual, isAlt: false})
							mappingQuality.refs = append(mappingQuality.refs, rank{value: mappingQual, isAlt: false})
						} else {
							baseQuality.alts = append(baseQuality.alts, rank{value: baseQual, isAlt: true})
							mappingQuality.alts = append(mappingQuality.alts, rank{value: mappingQual, isAlt: true})
						}
					} else {
						log.Panic("getReadCoordinateForReferenceCoordinate failed in annotateCall")
					}
					if softEnd(aln) >= int(call.Pos) {
						// ReadPosRankSumTest
						if pos := int(call.Pos); pos != leftmostSafePos {
							readCoord, ok = getReadCoordinateForReferenceCoordinate(aln.CIGAR, softStart, pos, right)
						}
						if ok && !isInsideDeletion(aln.CIGAR, readCoord) {
							var leadingHardClips, trailingHardClips int
							if ce := aln.CIGAR[0]; ce.Operation == 'H' {
								leadingHardClips = int(ce.Length)
							}
							if ce := aln.CIGAR[len(aln.CIGAR)-1]; ce.Operation == 'H' {
								trailingHardClips = int(ce.Length)
							}
							readPos := leadingHardClips + calcAlignmentByteArrayOffset(aln.CIGAR, readCoord)
							nofAlignedBases := nofAlignedBasesWithSoftClips(aln.CIGAR)
							numOriginalBases := nofAlignedBases + leadingHardClips + trailingHardClips
							if readPos > numOriginalBases/2 {
								readPos = numOriginalBases - (readPos + 1)
							}
							rp := float64(readPos)
							if isRef {
								readPosition.refs = append(readPosition.refs, rank{value: rp, isAlt: false})
							} else {
								readPosition.alts = append(readPosition.alts, rank{value: rp, isAlt: true})
							}
						}
					}
				}
			}
		}
	}
	if dp := len(likelihoods.alns); dp > 0 {
		call.Info.Set(DP, dp)
		call.Info.Set(RAW_MQandDP, []interface{}{squareSum, numReadsUsed})
	}
	if zscore, ok := baseQuality.mannWithneyU(); ok {
		call.Info.Set(BaseQRankSum, formatf(zscore, 3))
	}
	if zscore, ok := mappingQuality.mannWithneyU(); ok {
		call.Info.Set(MQRankSum, formatf(zscore, 3))
	}
	if zscore, ok := readPosition.mannWithneyU(); ok {
		call.Info.Set(ReadPosRankSum, formatf(zscore, 3))
	}
	for _, g := range callGT {
		if g >= 0 {
			call.Info.Set(ExcessHet, formatf(calculateEH(call, gls), 4))
			ad := make([]interface{}, len(call.Alt)+1)
			ad[0] = alleleCounts[call.Ref]
			for i, a := range call.Alt {
				ad[i+1] = alleleCounts[a]
			}
			call.GenotypeData[0].Data.Set(AD, ad)
			call.GenotypeData[0].Data.Set(DP, depth)
			if hc.confidenceMode != none {
				call.GenotypeData[0].Data.Set(SB, []interface{}{
					contingencyTable[0],
					contingencyTable[1],
					contingencyTable[2],
					contingencyTable[3],
				})
			}
			break
		}
	}
	if hc.confidenceMode == none {
		// SOR / StrandOddsRatio
		t00 := float64(contingencyTable[0]) + 1
		t01 := float64(contingencyTable[1]) + 1
		t10 := float64(contingencyTable[2]) + 1
		t11 := float64(contingencyTable[3]) + 1
		ratio := (t00/t01)*(t11/t10) + (t01/t00)*(t10/t11)
		if t00 > t01 {
			t00, t01 = t01, t00
		}
		refRatio := t00 / t01
		if t10 > t11 {
			t10, t11 = t11, t10
		}
		altRatio := t10 / t11
		sor := math.Log(ratio) + math.Log(refRatio) - math.Log(altRatio)
		call.Info.Set(SOR, formatf(sor, 3))

		// FS / FisherStrand
		if sum := contingencyTable[0] + contingencyTable[1] + contingencyTable[2] + contingencyTable[3]; sum <= 2 {
			call.Info.Set(FS, "0.000")
		} else {
			if sum > 400 {
				normFactor := float64(sum) / 200
				contingencyTable[0] = int(float64(contingencyTable[0]) / normFactor)
				contingencyTable[1] = int(float64(contingencyTable[1]) / normFactor)
				contingencyTable[2] = int(float64(contingencyTable[2]) / normFactor)
				contingencyTable[3] = int(float64(contingencyTable[3]) / normFactor)
			}
			m := contingencyTable[0] + contingencyTable[1]
			n := contingencyTable[2] + contingencyTable[3]
			k := contingencyTable[0] + contingencyTable[2]
			lo := maxInt(0, k-n)
			hi := minInt(k, m)
			var pValue float64
			if hi > lo {
				dist := makeHypergeometricDistribution(m+n, m, k)
				logds := make([]float64, 0, hi+1-lo)
				for i := lo; i <= hi; i++ {
					logds = append(logds, dist.logProbability(i))
				}
				threshold := logds[contingencyTable[0]-lo] * relErr
				for i := 0; i < len(logds); {
					if d := logds[i]; d <= threshold {
						logds[i] = d * math.Log10E
						i++
					} else {
						logds = append(logds[:i], logds[i+1:]...)
					}
				}
				pValue = math.Abs(-10 * math.Max(math.Log10(math.Max(math.Min(math.Pow(10, log10SumLog10Slice(logds)), 1), 1e-320)), minLog10ScaledQual))
			}
			call.Info.Set(FS, formatf(pValue, 3))
		}

		// QD / QualityByDepth
		for _, g := range callGT {
			if g > 0 {
				var qdDepth int
				if depth == 0 {
					qdDepth = len(likelihoods.alns)
				} else {
					qdDepth = depth
				}
				if qdDepth > 0 {
					qd := call.Qual.(float64) / float64(qdDepth)
					qd = fixTooHighQD(qd)
					call.Info.Set(QD, formatf(qd, 2))
				}
				break
			}
		}
	}
}

func isBiallelic(call *vcf.Variant) bool {
	if len(call.Alt) == 1 {
		return true
	}
	if len(call.Alt) == 2 {
		return contains(nonRef, call.Alt)
	}
	return false
}

func constructHaplotypeMapping(calls []*vcf.Variant, calledHaplotypes map[*haplotype]bool) map[*vcf.Variant]map[*haplotype]bool {
	haplotypeMap := make(map[*vcf.Variant]map[*haplotype]bool)
	for _, call := range calls {
		if !isBiallelic(call) {
			haplotypeMap[call] = make(map[*haplotype]bool)
			continue
		}
		alt := call.Alt[0]
		hapsWithAllele := make(map[*haplotype]bool)
		if alt == "*" {
			for h := range calledHaplotypes {
				for _, vc := range h.events {
					if vc.Pos < call.Pos && vc.End() >= call.Pos {
						hapsWithAllele[h] = true
					}
				}
			}
		} else {
			for h := range calledHaplotypes {
				for _, vc := range h.events {
					if vc.Pos == call.Pos && contains(alt, vc.Alt) {
						hapsWithAllele[h] = true
					}
				}
			}
		}
		haplotypeMap[call] = hapsWithAllele
	}
	return haplotypeMap
}

const (
	phase01 = "0|1"
	phase10 = "1|0"
)

type phaseSetID struct {
	id    int
	phase string
}

func containsAll(m1, m2 map[*haplotype]bool) bool {
	for k2 := range m2 {
		if !m1[k2] {
			return false
		}
	}
	return true
}

func intersectionIsEmpty(m1, m2 map[*haplotype]bool) bool {
	if len(m2) < len(m1) {
		m1, m2 = m2, m1
	}
	for k := range m1 {
		if m2[k] {
			return false
		}
	}
	return true
}

func constructPhaseSetMapping(calls []*vcf.Variant, haplotypeMap map[*vcf.Variant]map[*haplotype]bool, totalHaps int) (phaseSetMapping map[*vcf.Variant]phaseSetID, uniqueCounter int) {
	nofCalls := len(calls)
	phaseSetMapping = make(map[*vcf.Variant]phaseSetID)

	for i := 0; i < nofCalls-1; i++ {
		call := calls[i]
		haplotypesWithCall := haplotypeMap[call]
		if len(haplotypesWithCall) == 0 {
			continue
		}

		callIsOnAllHaps := len(haplotypesWithCall) == totalHaps

		for j := i + 1; j < nofCalls; j++ {
			comp := calls[j]
			haplotypesWithComp := haplotypeMap[comp]
			if len(haplotypesWithComp) == 0 {
				continue
			}

			compIsOnAllHaps := len(haplotypesWithComp) == totalHaps
			if callIsOnAllHaps || compIsOnAllHaps || (len(haplotypesWithCall) == len(haplotypesWithComp) && containsAll(haplotypesWithCall, haplotypesWithComp)) {
				if _, ok := phaseSetMapping[call]; !ok {
					if _, ok := phaseSetMapping[comp]; ok {
						return nil, 0
					}
					phaseSetMapping[call] = phaseSetID{uniqueCounter, phase01}
					phaseSetMapping[comp] = phaseSetID{uniqueCounter, phase01}
					uniqueCounter++
				} else if _, ok := phaseSetMapping[comp]; !ok {
					phaseSetMapping[comp] = phaseSetMapping[call]
				}
			} else if len(haplotypesWithCall)+len(haplotypesWithComp) == totalHaps {
				if intersectionIsEmpty(haplotypesWithCall, haplotypesWithComp) {
					if _, ok := phaseSetMapping[call]; !ok {
						if _, ok := phaseSetMapping[comp]; ok {
							for key := range phaseSetMapping {
								delete(phaseSetMapping, key)
							}
							return nil, 0
						}
						phaseSetMapping[call] = phaseSetID{uniqueCounter, phase01}
						phaseSetMapping[comp] = phaseSetID{uniqueCounter, phase10}
						uniqueCounter++
					} else if _, ok := phaseSetMapping[comp]; !ok {
						callPhase := phaseSetMapping[call]
						if callPhase.phase == phase01 {
							phaseSetMapping[comp] = phaseSetID{callPhase.id, phase10}
						} else {
							phaseSetMapping[comp] = phaseSetID{callPhase.id, phase01}
						}
					}
				}
			}
		}
	}
	return
}

var (
	// PID is the HaplotypeCaller phasing ID.
	PID = utils.Intern("PID")
	// PGT is the HaplotypeCaller phasing genotype.
	PGT = utils.Intern("PGT")
	// PS is the HaplotypeCaller phasing set ID.
	PS = utils.Intern("PS")
)

func phaseVC(vc *vcf.Variant, ID string, phaseGT string, phaseSetID int) {
	gt := vc.GenotypeData[0]
	if phaseGT == phase10 && gt.GT[0] >= 0 && gt.GT[1] > 0 && gt.GT[0] != gt.GT[1] {
		vc.GenotypeData[0].GT[0], vc.GenotypeData[0].GT[1] = gt.GT[1], gt.GT[0]
	}
	vc.GenotypeData[0].Phased = true
	vc.GenotypeData[0].Data.Set(PID, ID)
	vc.GenotypeData[0].Data.Set(PGT, phaseGT)
	vc.GenotypeData[0].Data.Set(PS, phaseSetID)
}

func constructPhaseGroups(calls []*vcf.Variant, phaseSetMapping map[*vcf.Variant]phaseSetID, uniqueCounterEndValue int) {
	for count := 0; count < uniqueCounterEndValue; count++ {
		var firstIndex int
		var firstCall *vcf.Variant
		var callPhase phaseSetID
		for firstIndex, firstCall = range calls {
			var ok bool
			callPhase, ok = phaseSetMapping[firstCall]
			if ok && callPhase.id == count {
				break
			}
		}
		uniqueID := fmt.Sprintf("%d_%s_%s", firstCall.Pos, firstCall.Ref, firstCall.Alt[0])
		phaseSetID := int(firstCall.Pos)
		phaseVC(firstCall, uniqueID, callPhase.phase, phaseSetID)
		for index := firstIndex + 1; index < len(calls); index++ {
			call := calls[index]
			var ok bool
			if callPhase, ok = phaseSetMapping[call]; ok && callPhase.id == count {
				phaseVC(call, uniqueID, callPhase.phase, phaseSetID)
			}
		}
	}
}

func decomposeHaplotypesIntoVariants(haplotypes []*haplotype, region *assemblyRegion) []int32 {
	startPosKeySet := parallel.RangeReduce(0, len(haplotypes), 0, func(low, high int) interface{} {
		startPosKeySet := make(map[int32]bool)
		for i := low; i < high; i++ {
			haplotypes[i].events = makeEventMap(fmt.Sprintf("HC%d", i), region.contig, haplotypes[i], region.reference, startPosKeySet)
		}
		return startPosKeySet
	}, func(x, y interface{}) interface{} {
		sx := x.(map[int32]bool)
		sy := y.(map[int32]bool)
		if len(sy) > len(sx) {
			sx, sy = sy, sx
		}
		for p := range sy {
			sx[p] = true
		}
		return sx
	}).(map[int32]bool)
	startPositions := make([]int32, 0, len(startPosKeySet))
	for p := range startPosKeySet {
		startPositions = append(startPositions, p)
	}
	sort.Slice(startPositions, func(i, j int) bool {
		return startPositions[i] < startPositions[j]
	})
	return startPositions
}

func makeMergedVariant(events []*vcf.Variant) *vcf.Variant {
	sortBySources(events)
	first := events[0]
	name := first.Source
	refAllele := first.Ref
	for i := 1; i < len(events); i++ {
		ref := events[i].Ref
		if len(ref) > len(refAllele) {
			refAllele = ref
		}
	}
	var altAlleles []string
	longestVC := first
	longestSize := variantSize(longestVC)
	for _, vc := range events {
		if size := variantSize(vc); size > longestSize {
			longestVC = vc
			longestSize = size
		}
		if refAllele == vc.Ref {
			for _, a := range vc.Alt {
				altAlleles = addAllele(altAlleles, a)
			}
		} else {
			extraBases := refAllele[len(vc.Ref):]
			for _, a := range vc.Alt {
				if a == "*" {
					altAlleles = addAllele(altAlleles, a)
				} else if !isSymbolicAllele(a) {
					altAlleles = addAllele(altAlleles, a+extraBases)
				}
			}
		}
	}
	merged := &vcf.Variant{
		Source: name,
		Chrom:  longestVC.Chrom,
		Pos:    longestVC.Pos,
		ID:     emptyId,
		Ref:    refAllele,
		Alt:    altAlleles,
	}
	merged.SetEnd(longestVC.End())
	return merged
}

func createAlleleMapper(variant *vcf.Variant, haplotypes []*haplotype, overlaps map[*haplotype][]*vcf.Variant, loc int32) *alleleMap {
	vRef := variant.Ref
	alleleMapper := newAlleleMap(vRef)
	// the next few lines establish an order of keys in the alleleMapper
	for _, a := range variant.Alt {
		if !isSymbolicAllele(a) {
			alleleMapper.addAllele(a)
		}
	}
	for _, h := range haplotypes {
		spanningEvents := overlaps[h]
		if len(spanningEvents) == 0 {
			alleleMapper.haplotypes[vRef] = append(alleleMapper.haplotypes[vRef], h)
			continue
		}
		for _, spanningEvent := range spanningEvents {
			if spanningEvent.Pos == loc {
				if firstAlt := spanningEvent.Alt[0]; firstAlt == "*" {
					alleleMapper.maybeAdd("*", h)
				} else {
					alleleMapper.maybeAdd(firstAlt+vRef[len(spanningEvent.Ref):], h)
				}
			} else {
				alleleMapper.add("*", h)
				break
			}
		}
	}
	return alleleMapper
}

func reduceAltAlleles(variant *vcf.Variant, alleleMapper *alleleMap) {
	scored := make([]alleleScoredByHaplotype, 0, len(alleleMapper.alleles))
	bestScore := math.Inf(-1)
	secondBestScore := math.Inf(-1)
	for _, h := range alleleMapper.haplotypes[variant.Ref] {
		if h.score > bestScore {
			bestScore, secondBestScore = h.score, bestScore
		} else if h.score > secondBestScore {
			secondBestScore = h.score
		}
	}
	scored = append(scored, alleleScoredByHaplotype{
		allele:          variant.Ref,
		isRef:           true,
		bestScore:       bestScore,
		secondBestScore: secondBestScore,
	})
	for i := 1; i < len(alleleMapper.alleles); i++ {
		allele := alleleMapper.alleles[i]
		bestScore := math.Inf(-1)
		secondBestScore := math.Inf(-1)
		for _, h := range alleleMapper.haplotypes[allele] {
			if h.score > bestScore {
				bestScore, secondBestScore = h.score, bestScore
			} else if h.score > secondBestScore {
				secondBestScore = h.score
			}
		}
		scored = append(scored, alleleScoredByHaplotype{
			allele:          allele,
			isRef:           false,
			bestScore:       bestScore,
			secondBestScore: secondBestScore,
		})
	}
	sort.SliceStable(scored, func(i, j int) bool {
		l, r := scored[i], scored[j]
		if l.isRef && !r.isRef {
			return true
		}
		if !l.isRef && r.isRef {
			return false
		}
		if l.bestScore > r.bestScore {
			return true
		}
		if l.bestScore < r.bestScore {
			return false
		}
		if l.secondBestScore > r.secondBestScore {
			return true
		}
		if l.secondBestScore < r.secondBestScore {
			return false
		}
		return l.allele < r.allele
	})
	allelesToRemove := make(map[string]bool, len(scored)-maxAcceptableAlleleCount)
	for _, toRemove := range scored[maxAcceptableAlleleCount:] {
		alleleMapper.remove(toRemove.allele)
		allelesToRemove[toRemove.allele] = true
	}
	for i := 0; i < len(variant.Alt); {
		if allelesToRemove[variant.Alt[i]] {
			variant.Alt = append(variant.Alt[:i], variant.Alt[i+1:]...)
		} else {
			i++
		}
	}
}

func calculateGenotypeLikelihoods(variant *vcf.Variant, likelihoods readAlleleLikelihoods) (gls []float64, pls []interface{}) {
	denominator := float64(len(likelihoods.alns)) * log10Ploidy
	nofAlleles := len(variant.Alt) + 1
	gls = make([]float64, (nofAlleles*nofAlleles+nofAlleles)/2)
	likelihoodForRef := likelihoods.values[variant.Ref]
	maxGL := computeSingleComponentGenotypeLikelihood(likelihoodForRef) - denominator
	gls[0] = maxGL
	forEachAltGenotype(variant.Ref, variant.Alt,
		func(index int, alt string) {
			gl := computeTwoComponentGenotypeLikelihood(likelihoodForRef, likelihoods.values[alt]) - denominator
			if gl > maxGL {
				maxGL = gl
			}
			gls[index] = gl
		}, func(index int, alt string) {
			gl := computeSingleComponentGenotypeLikelihood(likelihoods.values[alt]) - denominator
			if gl > maxGL {
				maxGL = gl
			}
			gls[index] = gl
		}, func(index int, alt1, alt2 string) {
			gl := computeTwoComponentGenotypeLikelihood(likelihoods.values[alt1], likelihoods.values[alt2]) - denominator
			if gl > maxGL {
				maxGL = gl
			}
			gls[index] = gl
		})
	pls = make([]interface{}, len(gls))
	for i := range gls {
		var pl int
		var gl float64
		if adjusted := -10 * (gls[i] - maxGL); adjusted > math.MaxInt32 {
			pl = math.MaxInt32
			gl = float64(math.MaxInt32) / -10
		} else {
			r := math.Round(adjusted)
			pl = int(r)
			gl = r / -10
		}
		pls[i] = pl
		gls[i] = gl
	}
	return gls, pls
}

func (hc *HaplotypeCaller) computeReadAlleleLikelihoodsForAnnotation(call *vcf.Variant, alleleLikelihoods readAlleleLikelihoods, likelihoods readLikelihoods, alleleMapper *alleleMap, filteredAlns []*sam.Alignment) readAlleleLikelihoods {
	// todo: in prepareReadAlleleLikelihoods, check for sample contamination
	cstart := call.Pos
	cend := call.End()
	if true { // todo: !isSampleContaminationPresent
		for i := 0; i < len(alleleLikelihoods.alns); {
			aln := alleleLikelihoods.alns[i]
			// the following overlap test is convoluted, but this is necessary because rend may not always be >= rstart
			if rstart, rend := aln.POS, aln.End(); (cstart >= rstart && cstart <= rend) || (cend >= rstart && cend <= rend) || (rstart >= cstart && rend <= cend) {
				i++
			} else {
				alleleLikelihoods.alns = append(alleleLikelihoods.alns[:i], alleleLikelihoods.alns[i+1:]...)
				for _, a := range alleleLikelihoods.alleles {
					values := alleleLikelihoods.values[a]
					alleleLikelihoods.values[a] = append(values[:i], values[i+1:]...)
				}
			}
		}
		if len(alleleLikelihoods.alleles) != len(call.Alt)+1 {
			alleleLikelihoods.updateNonRef(call.Ref, call.Alt)
		}
	} else {
		// we need to do this again because of potential sample contamination!
		alleleLikelihoods = marginalize(likelihoods, alleleMapper, cstart, cend)
		if hc.confidenceMode != none {
			alleleLikelihoods.alleles = append(alleleLikelihoods.alleles, nonRef)
			alleleLikelihoods.updateNonRef(call.Ref, call.Alt)
		}
	}
	for _, aln := range filteredAlns {
		if rstart, rend := aln.POS, aln.End(); rstart <= cend && cstart <= rend {
			alleleLikelihoods.alns = append(alleleLikelihoods.alns, aln)
			for _, a := range alleleLikelihoods.alleles {
				alleleLikelihoods.values[a] = append(alleleLikelihoods.values[a], 0)
			}
		}
	}
	return alleleLikelihoods
}

func reverseTrimAlleles(call *vcf.Variant) {
	trim := len(call.Ref) - 1
	if trim < 1 {
		return
	}
	for _, a := range call.Alt {
		if isSymbolicAllele(a) {
			continue
		}
		for i := 0; i <= trim; i++ {
			if i == len(a) {
				trim = i - 1
				break
			}
			if a[len(a)-i-1] != call.Ref[len(call.Ref)-i-1] {
				trim = i
				break
			}
		}
		if trim < 1 {
			return
		}
	}
	call.Ref = call.Ref[:len(call.Ref)-trim]
	for i, a := range call.Alt {
		if !isSymbolicAllele(a) {
			call.Alt[i] = a[:len(a)-trim]
		}
	}
}

func computeGenotypeFormat(call *vcf.Variant) {
	sort.Slice(call.Info, func(i, j int) bool {
		return *call.Info[i].Key < *call.Info[j].Key
	})
	sort.Slice(call.GenotypeData[0].Data, func(i, j int) bool {
		return *call.GenotypeData[0].Data[i].Key < *call.GenotypeData[0].Data[j].Key
	})
	genotypeKeys := make([]utils.Symbol, 0, len(call.GenotypeData[0].Data)+1)
	genotypeKeys = append(genotypeKeys, vcf.GT)
	for _, entry := range call.GenotypeData[0].Data {
		genotypeKeys = append(genotypeKeys, entry.Key)
	}
	call.GenotypeFormat = genotypeKeys
}

func (hc *HaplotypeCaller) assignGenotypeLikelihoods(hdr *sam.Header, region *assemblyRegion, filteredAlns []*sam.Alignment, haplotypes []*haplotype, likelihoods readLikelihoods) (returnCalls []*vcf.Variant, calledHaplotypes map[*haplotype]bool) {
	startPositions := decomposeHaplotypesIntoVariants(haplotypes, region)

	calledHaplotypes = make(map[*haplotype]bool)
	deletions := region.deletions.handleDeletions()

	var containsCalls bool

	for _, loc := range startPositions {
		if loc < region.start || loc > region.end {
			continue
		}

		overlaps := getOverlappingEvents(loc, haplotypes)
		events := computeActiveVariantContextsWithSpanDelsReplaced(loc, haplotypes, overlaps, region.reference)

		if len(events) == 0 {
			continue
		}

		merged := makeMergedVariant(events)
		alleleMapper := createAlleleMapper(merged, haplotypes, overlaps, loc)
		if len(alleleMapper.alleles) > maxAcceptableAlleleCount {
			reduceAltAlleles(merged, alleleMapper)
		}

		alleleLikelihoods := marginalize(likelihoods, alleleMapper, maxInt32(merged.Pos-2, 1), minInt32(merged.End()+2, region.contigLength))

		// todo: handle --comtamination-fraction-to-filter parameter here

		if hc.confidenceMode != none {
			merged.Alt = append(merged.Alt, nonRef)
			alleleLikelihoods.alleles = append(alleleLikelihoods.alleles, nonRef)
			alleleLikelihoods.updateNonRef(alleleLikelihoods.alleles[0], alleleLikelihoods.alleles[1:])
		}

		gls, pls := calculateGenotypeLikelihoods(merged, alleleLikelihoods)
		if call, gls := hc.calculateGenotypes(merged, pls, gls, deletions); call != nil {
			if !containsCalls {
				for _, g := range call.GenotypeData[0].GT {
					if g >= 0 {
						containsCalls = true
						break
					}
				}
			}
			alleleLikelihoods = hc.computeReadAlleleLikelihoodsForAnnotation(call, alleleLikelihoods, likelihoods, alleleMapper, filteredAlns)
			hc.annotateCall(call, alleleLikelihoods, gls)
			if len(call.Alt) > 0 && len(call.Alt) != len(merged.Alt) {
				reverseTrimAlleles(call)
			}
			returnCalls = append(returnCalls, call)
			for _, h := range alleleMapper.haplotypes[call.Ref] {
				calledHaplotypes[h] = true
			}
			for _, a := range call.Alt {
				for _, h := range alleleMapper.haplotypes[a] {
					calledHaplotypes[h] = true
				}
			}
		}
	}

	deletions.close()

	if hc.confidenceMode != none {
		if !containsCalls {
			return nil, nil
		}
		haplotypeMap := constructHaplotypeMapping(returnCalls, calledHaplotypes)
		phaseSetMapping, uniqueCounterEndValue := constructPhaseSetMapping(returnCalls, haplotypeMap, len(calledHaplotypes)-1)
		constructPhaseGroups(returnCalls, phaseSetMapping, uniqueCounterEndValue)
	}

	return returnCalls, calledHaplotypes
}
